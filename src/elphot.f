	subroutine elphot (J,nf,ib,nnx,nny,nx,ny,sarray,MASK,iMask,Unc,
     1     rmax,zmed,csig,
     1     xc,yc,rat,phimax,ntot,ztot,xint, flag, iflag, SUM2err, fbits )

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (pi=3.1415926536, valmin = -500., 
     1             rtrunc = 0.5, rtruncmax = 1800.)
	real*4  sarray(nnx,nny,nf,4), Unc(nnx,nny,nf,4)
	integer flag, iflag, iMask(nnx,nny,nf,4),fatal,ibit,BitSet,fbits(32)
	integer*2 MASK(nnx,nny,4)

	logical the_end

c	fatal = 2**2 + 2**18            ! bits 2 and 18 are regarded as fatal

	fatal = 0
        do I = 1,32
                ibit = fbits(I)
                if (ibit.ge.0) then
                  fatal = fatal + (2**ibit)
		else
                        !quit
                        goto 74
                endif
        enddo
 74     I = 0

	if (rmax.lt.rtrunc) rmax = rtrunc
	if (rmax.gt.rtruncmax) rmax=rtruncmax

	the_end = .false.

	SUM2err = 0.
	ntot = 0
	xint = 0.
	ztot = 0.
	iflag = 0

        jstrt=yc-rmax-2
	jstrt = max(jstrt,1)
    	jend=yc+rmax+2
	jend=min(jend,ny)
    	istrt=xc-rmax-2
	istrt = max(istrt,1)
    	iend=xc+rmax+2
	iend=min(iend,nx)

	do 100 JJ=jstrt,jend
	      	dy=(jj*1.)-yc
	do 101 II=istrt,iend
		dx=(ii*1.)-xc

                call w_ell (rat,phimax,dx,dy,dr)

		if ( (dr.ge.0.).and.(dr.le.rmax+1.20) ) then
c inside  ellipse !!!!

			if ((jj.lt.1).or.(jj.gt.ny)) then
				the_end=.true.
                                goto 100
                        endif

			if ((ii.lt.1).or.(ii.gt.nx)) then
				the_end=.true.
				goto 101
			endif

			del = abs(dr-rmax)
			if (del.le.1.25) then
                                !  check frac pixel
		call fracpix (xc,yc,ii,jj,rmax,rat,
     1                  phimax,farea)

                        else
                                  farea = 1.0
                        endif

			z = sarray (ii,jj,J,ib) - zmed  ! subtract the background
			if (MASK  (ii,jj,ib).lt.0) then
				flag = 1  ! star contamination
			endif

ccccc Check the iMask for bad pixels, etc
cccccc   fatal == bits 2 or 18

			if (farea.gt.0.) then
			  iMaskval = iMask(ii,jj,J,ib)
			  if (iMaskval.gt.0) iflag = 3

			  icheck = 0
			  icheck = iand(iMaskval, fatal)
			  if (icheck.gt.0) iflag = 2   !   bit 2 or 18 set

			endif

		        if (z.gt.Valmin) then

		       		xint = xint + (z*farea)   ! integrated flux
				ztot = ztot + farea       ! total area
		       		ntot=ntot+nint(farea)

				errval = Unc (ii,jj,J,ib)
				SUM2err = SUM2err + (errval**2)

			endif

		  endif

 101		continue
 100	  continue

	  ntot = nint(ztot)

	return
	end

c  (x0,y0) = coordinates of source 
c  i,j = pixel of interest
c  R is target radius
c  rat & phimax the elliptical params
c  farea is fractional area of pixel

c modified (sped up) by T. Conrow

	subroutine fracpix (x0,y0,i,j,R,rat,phimax,farea)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (eps = 0.001)

	if (R.le.1.) then
	  block = 15.0
        else if (R.le.3.2) then
          block = 10.
	else
	   block = 5.
	endif

	area = block**2
	nb = block
	MID = nint(block/2.0)

	R2 = (R+eps)**2

	nsum = 0

	if(rat .lt. 0.99) then

	  do l=1,nb
		DY = ((L - MID) / block) + j - Y0
          do k=1,nb
                DX = ((k - MID) / block) + i - X0
		call w_ell2 (rat,phimax,dx,dy,dr2)
		if (R2.ge.DR2) nsum=nsum+1
	  enddo
	  enddo

	else

	  do l=1,nb
		DY2 = ( ((L - MID) / block) + j - Y0 )**2
	  do k=1,nb
                DX2 = ( ((k - MID) / block) + i - X0 )**2
		if (R2.ge.DX2+DY2) nsum=nsum+1
	  enddo
	  enddo

	endif

	farea = (nsum*1.0) / area

	return
	end


        subroutine w_ell2 (e,phimax,dx,dy,a2)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        parameter (pi = 3.14159265)

        d=-phimax*pi/180.
        r2 = (dx**2) + (dy**2)

        if(dx.ne.0) then
          t=atan(dy/dx)
        else
          t=pi/2.
        endif

        ct = cos(t)
        st = sin(t)
        cd = cos(d)
        sd = sin(d)

        c1= (((e * ct)**2 ) + (st**2) ) * (cd**2)
        c2 = (1. - (e**2)) * 2. * st * ct * sd * cd
        c3 = (((e * st)**2 ) + (ct**2) ) * (sd**2)

        a2 = (r2 / (e**2)) * (c1 + c2 + c3)

        return
        end


        subroutine w_ell (e,phimax,dx,dy,a)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        parameter (pi = 3.14159265)

        d=-phimax*pi/180.
        r2 = (dx**2) + (dy**2)

        if(dx.ne.0) then
          t=atan(dy/dx)
        else
          t=pi/2.
        endif

        ct = cos(t)
        st = sin(t)
        cd = cos(d)
        sd = sin(d)

        c1= (((e * ct)**2 ) + (st**2) ) * (cd**2)
        c2 = (1. - (e**2)) * 2. * st * ct * sd * cd
        c3 = (((e * st)**2 ) + (ct**2) ) * (sd**2)

        a = sqrt( (r2 / (e**2)) * (c1 + c2 + c3) )

        return
        end

