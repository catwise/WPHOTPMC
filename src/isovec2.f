 	SUBROUTINE isovec2 (nx,ny,array,R,x0,y0,rmax,
     1      ratio,phimax,iso,nx_iso,ny_iso,zlim)
	
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        real*4 rad(100),rval(100),ddx(100),ddy(100)
	integer*2 iso(100,100),ibf(100,100)
        real*4 array(nx,ny)
        logical debug

	debug = .false.

c	write (6,*) 'ISOVEC ',zlim

c flush iso
        do j=1,100
        do i=1,100
                iso(i,j) = 0
        enddo
        enddo


	idex = 0
	nz = 0
	rvec = R
	ivec = nint(rvec)
	i0 = nint(x0)
	j0 = nint(y0)
	ih = i0 + ivec
	il = i0 - ivec
	jl = j0 - ivec
	jh = j0 + ivec

	idex0 = int(R)
        jdex0 = int(R)

	iso(idex0,jdex0) = 1

	if (il.lt.1) il=1
	if (jl.lt.1) jl=1
	if (ih.gt.nx) ih=nx
	if (jh.gt.ny) jh=ny


ccccccccccccccccccccccccccccccccccccccccccccccc
	del = 1.0
	theta = 1.0 - del

	do 1000 k=1,360
	  theta = theta + del
	
	  if (theta.le.90.0) goto 10
	  if (theta.le.180.0) goto 11
	  if (theta.le.270.0) goto 12
	  if (theta.gt.360.) goto 1000
	  goto 13

ccccccc  QUAD 1
 10	 i=1
 
c	write (6,*) 'QUAD 1 begin ',theta,i0,j0
	 n = 0

	 if (theta.lt.45.0) then
		NNL = i0
		NNH = ih
		IDN = 1
		deta = theta 
	 else
		NNL = j0
		NNH = jh
		IDN = 1
		deta =  90.0 - theta 
	 endif

	    do 50 M = NNL,NNH, IDN
		if (theta.lt.45.0) then
		  i = M
		  dx = i - i0
		  dy = dx * tan(deta/57.2957795)
		  j = nint(dy) + j0
		else
		   j = M
		  dy = j - j0
		  dx = dy * tan(deta/57.2957795)
		  i = nint(dx) + i0
		endif


		dr = sqrt ((dx**2) + (dy**2) )
                if (dr.gt.rvec) goto 50

		val = array(i,j)
                if (val.lt.-99.) goto 50

		n = n + 1
                rad(n) = dr
                rval(n) = val
		ddx (n) = dx
		ddy (n) = dy

		if (val.lt.zlim-0.01) goto 50

 50	 continue

c now find n-sigma perimeter

	 idex1 = 0
	 idex2 = 0

	 do L=1,n
		val = rval(L)
		if (val.eq.zlim) then
			idex1 = L
			idex2 = 0
			goto 51
		endif

		if (val.lt.zlim) then
			idex1 = L-1
			idex2 = L
			goto 51
		endif

	enddo

 51     if (idex1.eq.0) then
		idex = 0
		jdex = 0
		drval = 0.
	else if (idex2.eq.0) then
		dx = ddx(idex1)
		dy = ddy (idex1)
		drval = rad (idex1)
                idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)
	else
		dr1 = rad (idex1)
		dr2 = rad (idex2)
		rv1 = rval (idex1)
		rv2 = rval (idex2)     

		call interp (dr1,dr2,rv1,rv2,zlim,drval)


		if (theta.lt.45.0) then

		  dy =  drval * sin(deta/57.2957795)
		  dx =  drval * cos (deta/57.2957795)
		
		else
		  
	  	  dx =  drval * sin(deta/57.2957795)
                  dy =  drval * cos (deta/57.2957795)
                
		endif

		idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)

c		write (6,*) 'dr1,dr2,rv1,rv2,zlim,drval',dr1,dr2,rv1,rv2,zlim,drval
c		write (6,*) 'dx,dy,idex,jdex ',dx,dy,idex,jdex	
	endif	


	if (drval.gt.0) then

		do L=1,n
                  rr = rad(L)
		  if (rr.lt.drval) then

		if (theta.lt.45.0) then

                  dy =  rr * sin(deta/57.2957795)
                  dx =  rr * cos (deta/57.2957795)

                else

                  dx =  rr * sin(deta/57.2957795)
                  dy =  rr * cos (deta/57.2957795)

                endif

                ii = idex0 + nint(dx)
                jj = jdex0 + nint(dy)

		if (iso (ii,jj).le.1) iso (ii,jj) = 1

		endif
		
		enddo

	endif

	if ( (idex.gt.0).and.(jdex.gt.0) ) then
                iso (idex,jdex) = 100
        endif


	goto 1000


ccccccc  QUAD 2
 11	 i=1
 
c	write (6,*) 'QUAD 2 begin ',theta,i0,j0
	 n = 0

	 if (theta.gt.135.0) then
		NNL = i0
		NNH = il
		IDN = -1
		deta = 180.0 - theta 
	 else
		NNL = j0
		NNH = jh
		IDN = 1
		deta =  theta - 90.0 
	 endif

	    do 60 M = NNL,NNH, IDN
		if (theta.gt.135.0) then
		  i = M
		  dx = i - i0
		  dy = -dx * tan(deta/57.2957795)
		  j = nint(dy) + j0
		else
		   j = M
		  dy = j - j0
		  dx = -dy * tan(deta/57.2957795)
		  i = nint(dx) + i0
		endif


		dr = sqrt ((dx**2) + (dy**2) )
                if (dr.gt.rvec) goto 60

		val = array(i,j)
                if (val.lt.-99.) goto 60

		n = n + 1
                rad(n) = dr
                rval(n) = val
		ddx (n) = dx
		ddy (n) = dy

		if (val.lt.zlim-0.01) goto 60

 60	 continue

c now find n-sigma perimeter

	 idex1 = 0
	 idex2 = 0
	 do L=1,n
		val = rval(L)
		if (val.eq.zlim) then
			idex1 = L
			idex2 = 0
			goto 61
		endif

		if (val.lt.zlim) then
			idex1 = L-1
			idex2 = L
			goto 61
		endif

	enddo

 61     if (idex1.eq.0) then
		idex = 0
		jdex = 0
	else if (idex2.eq.0) then
		dx = ddx(idex1)
		dy = ddy (idex1)

                idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)
	else
		dr1 = rad (idex1)
		dr2 = rad (idex2)
		rv1 = rval (idex1)
		rv2 = rval (idex2)     

		call interp (dr1,dr2,rv1,rv2,zlim,drval)


		if (theta.gt.135.0) then

		  dy =  drval * sin(deta/57.2957795)
		  dx =  -drval * cos (deta/57.2957795)
		
		else
		  
	  	  dx =  -drval * sin(deta/57.2957795)
                  dy =  drval * cos (deta/57.2957795)
                
		endif

		idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)

	endif	

	        if (drval.gt.0) then

                do L=1,n
                  rr = rad(L)
                  if (rr.lt.drval) then


                if (theta.gt.135.0) then

                  dy =  rr * sin(deta/57.2957795)
                  dx =  -rr * cos (deta/57.2957795)

                else
                 
                  dx =  -rr * sin(deta/57.2957795)
                  dy =  rr * cos (deta/57.2957795)

                endif

                ii = idex0 + nint(dx)
                jj = jdex0 + nint(dy)

		if (iso (ii,jj).le.1) iso (ii,jj) = 1

		endif

		enddo

		endif


	if ( (idex.gt.0).and.(jdex.gt.0) ) then
		iso (idex,jdex) = 100
	endif


	goto 1000


ccccccc  QUAD 3
 12	 i=1
 
c	write (6,*) 'QUAD 3 begin ',theta,i0,j0
	 n = 0

	 if (theta.lt.225.0) then
		NNL = i0
		NNH = il
		IDN = -1
		deta = theta - 180.0 
	 else
		NNL = j0
		NNH = jl
		IDN = -1
		deta = 270.0 - theta 
	 endif

	    do 70 M = NNL,NNH, IDN
		if (theta.lt.225.0) then
		  i = M
		  dx = i - i0
		  dy = dx * tan(deta/57.2957795)
		  j = nint(dy) + j0
		else
		   j = M
		  dy = j - j0
		  dx = dy * tan(deta/57.2957795)
		  i = nint(dx) + i0
		endif


		dr = sqrt ((dx**2) + (dy**2) )
                if (dr.gt.rvec) goto 70

		val = array(i,j)
                if (val.lt.-99.) goto 70

		n = n + 1
                rad(n) = dr
                rval(n) = val
		ddx (n) = dx
		ddy (n) = dy

		if (val.lt.zlim-0.01) goto 70

 70	 continue

c now find n-sigma perimeter

	 idex1 = 0
	 idex2 = 0
	 do L=1,n
		val = rval(L)
		if (val.eq.zlim) then
			idex1 = L
			idex2 = 0
			goto 71
		endif

		if (val.lt.zlim) then
			idex1 = L-1
			idex2 = L
			goto 71
		endif

	enddo

 71     if (idex1.eq.0) then
		idex = 0
		jdex = 0
	else if (idex2.eq.0) then
		dx = ddx(idex1)
		dy = ddy (idex1)

                idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)
	else
		dr1 = rad (idex1)
		dr2 = rad (idex2)
		rv1 = rval (idex1)
		rv2 = rval (idex2)     

		call interp (dr1,dr2,rv1,rv2,zlim,drval)


		if (theta.lt.225.0) then

		  dy = -drval * sin(deta/57.2957795)
		  dx = -drval * cos (deta/57.2957795)
		
		else
		  
	  	  dx = -drval * sin(deta/57.2957795)
                  dy = -drval * cos (deta/57.2957795)
                
		endif

		idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)

c		write (6,*) 'dr1,dr2,rv1,rv2,zlim,drval',dr1,dr2,rv1,rv2,zlim,drval
c		write (6,*) 'dx,dy,idex,jdex ',dx,dy,idex,jdex	
	endif	


	        if (drval.gt.0) then

                do L=1,n
                  rr = rad(L)
                  if (rr.lt.drval) then


		if (theta.lt.225.0) then

                  dy = -rr * sin(deta/57.2957795)
                  dx = -rr * cos (deta/57.2957795)

                else
                 
                  dx = -rr * sin(deta/57.2957795)
                  dy = -rr * cos (deta/57.2957795)

                endif

                ii = idex0 + nint(dx)
                jj = jdex0 + nint(dy)

		if (iso (ii,jj).le.1) iso(ii,jj) = 1

		endif
		enddo
		endif

	if ( (idex.gt.0).and.(jdex.gt.0) ) then
		iso (idex,jdex) = 100
	endif


	goto 1000


ccccccc  QUAD 4
 13	 i=1
 
c	write (6,*) 'QUAD 4 begin ',theta,i0,j0
	 n = 0

	 if (theta.gt.315.0) then
		NNL = i0
		NNH = ih
		IDN = 1
		deta = 360 - theta 
	 else
		NNL = j0
		NNH = jl
		IDN = -1
		deta =  theta  - 270.0
	 endif

	    do 80 M = NNL,NNH, IDN
		if (theta.gt.315.0) then
		  i = M
		  dx = i - i0
		  dy = -dx * tan(deta/57.2957795)
		  j = nint(dy) + j0
		else
		   j = M
		  dy = j - j0
		  dx = -dy * tan(deta/57.2957795)
		  i = nint(dx) + i0
		endif


		dr = sqrt ((dx**2) + (dy**2) )
                if (dr.gt.rvec) goto 80

		val = array(i,j)
                if (val.lt.-99.) goto 80

		n = n + 1
                rad(n) = dr
                rval(n) = val
		ddx (n) = dx
		ddy (n) = dy

		if (val.lt.zlim-0.01) goto 80

 80	 continue

c now find n-sigma perimeter

	 idex1 = 0
	 idex2 = 0
	 do L=1,n
		val = rval(L)
c		write (6,*) L,ddx(l),ddy(l),rad(L),val
		if (val.eq.zlim) then
			idex1 = L
			idex2 = 0
			goto 81
		endif

		if (val.lt.zlim) then
			idex1 = L-1
			idex2 = L
			goto 81
		endif

	enddo

 81     if (idex1.eq.0) then
		idex = 0
		jdex = 0
	else if (idex2.eq.0) then
		dx = ddx(idex1)
		dy = ddy (idex1)

                idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)
	else
		dr1 = rad (idex1)
		dr2 = rad (idex2)
		rv1 = rval (idex1)
		rv2 = rval (idex2)     

		call interp (dr1,dr2,rv1,rv2,zlim,drval)


		if (theta.gt.315.0) then

		  dy = -drval * sin(deta/57.2957795)
		  dx =  drval * cos (deta/57.2957795)
		
		else
		  
	  	  dx =  drval * sin(deta/57.2957795)
                  dy = -drval * cos (deta/57.2957795)
                
		endif

		idex = idex0 + nint(dx)
                jdex = jdex0 + nint(dy)

c		write (6,*) 'dr1,dr2,rv1,rv2,zlim,drval',dr1,dr2,rv1,rv2,zlim,drval
c		write (6,*) 'dx,dy,idex,jdex ',dx,dy,idex,jdex	
	endif	

        if (drval.gt.0) then

                do L=1,n
                  rr = rad(L)
                  if (rr.lt.drval) then

		if (theta.gt.315.0) then

                  dy = -rr * sin(deta/57.2957795)
                  dx =  rr * cos (deta/57.2957795)

                else

                  dx =  rr * sin(deta/57.2957795)
                  dy = -rr * cos (deta/57.2957795)

                endif

                ii = idex0 + nint(dx)
                jj = jdex0 + nint(dy)

		if (iso (ii,jj).le.1) iso(ii,jj) = 1
		
		endif
		enddo
		endif


 
	if ( (idex.gt.0).and.(jdex.gt.0) ) then
		iso (idex,jdex) = 100
	endif


	goto 1000

 1000	continue


c com ellipse parameters

        rmin = 99.
        rmax = 0.

	dmaxx = 0
	dmaxy = 0

        idel = 0
        jdel = 0
	IRR = nint(R)

        nj=0
        do j=1,(IRR*2)-1
		jj = (j - jdex0) + j0
        do i=1,(IRR*2)-1
		ii = (i - idex0) + i0

                ival = iso(i,j)
                dx = (i*1.) - idex0
                dy = (j*1.) - jdex0
                dr = sqrt( (dx**2) + (dy**2))

                if (ival.eq.100) then
                        if (dr.lt.rmin) then
                                rmin=dr
                                dminx = dx
                                dminy = dy
                        endif

                        nj=nj+1
                        if (dr.gt.rmax) then
                                rmax=dr
                                dmaxx = dx
                                dmaxy = dy
                        endif

                        if (abs(dx).gt.idel) idel = abs(dx)
                        if (abs(dy).gt.jdel) jdel = abs(dy)

                endif
        enddo
        enddo


        if (rmax.ge.1) then
c due tp pixelation, rmin uncertainty of 1
                ratio = (rmin+1.) / rmax
                if (ratio.gt.1.) ratio=1.0
        else
                ratio = 1.
        endif


        if (dmaxx.eq.0.) then
                phimax = 90.
        else
          phimax = 57.2957795 * atan (dmaxy / dmaxx)
        endif


c	write (6,*) 'rmin,rmax,phimax',rmin,rmax,phimax

cc now collapse isomask to smaller box

        itemp = 2 * (idel / 2)
        idel = itemp + 5

        itemp = 2 * (jdel / 2)
        jdel = itemp + 5

        il2 = idex0 - idel
        ih2 = idex0 + idel
        jl2 = jdex0 - jdel
        jh2 = jdex0 + jdel


        jdex=0
        do j=jl2,jh2
           idex=0
           jdex=jdex+1
        do i=il2,ih2
           idex=idex+1
           ibf(idex,jdex) = iso(i,j)
        enddo
        enddo

        do j=1,jdex
        do i=1,idex
                iso(i,j)=ibf(i,j)
        enddo
        enddo

        nx_iso=idex
        ny_iso=jdex

	 return
	end

