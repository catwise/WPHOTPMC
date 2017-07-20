      SUBROUTINE MOMENT(DATA,N,AVE,SDEV)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION DATA(N)
      real*8 s,var,p,t

      IF(N.eq.1) then
		AVE = DATA(1)
		SDEV = 0.
		goto 99
      else IF(N.Lt.1)then
c		write (6,*) 'error N < 1'
		goto 99
      endif

      S=0.d0
      DO 11 J=1,N
        S=S+(DATA(J)*1.d0)
11    CONTINUE
	
      AVE=sngl(S)/N
      VAR=0.d0
      DO 12 J=1,N
        S=(DATA(J)-AVE)*1.d0
	P=S*S
        VAR=VAR+P
12    CONTINUE
      VAR=VAR/(N-1)
      t = var ** 0.5
      SDEV=sngl(t)
     
 99	i=1 
      RETURN
      END


	SUBROUTINE ModMOMENT(DATA,N,SDEV)

        implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION DATA(N)
      real*8 s,var,p,t

      Ave = 0.  ! this is the mode
      SDEV = 0.

      IF(N.eq.1) then
                SDEV = 0.
                goto 99
      else IF(N.Lt.1)then
                goto 99
      endif

      S=0.d0
      VAR=0.d0

      DO 12 J=1,N
        S=(DATA(J)-AVE)*1.d0
        P=S*S
        VAR=VAR+P
12    CONTINUE
      VAR=VAR/(N-1)
      t = var ** 0.5
      SDEV=sngl(t)

 99     i=1
      RETURN
      END


      SUBROUTINE MDIAN1(X,N,XMED)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION X(N)
      
      N2=N/2
      CALL SORT(N,X)
      
      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
      ELSE
        XMED=X(N2+1)
      ENDIF


      RETURN
      END

	SUBROUTINE MDIAN2(X,Y,N,XMED, YMED)

        implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION X(N), Y(N)

      N2=N/2
      CALL SORT2(N,X, Y)

      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
	XMED2=0.5*(Y(N2)+Y(N2+1))
      ELSE
        XMED=X(N2+1)
	YMED=Y(N2+1)
      ENDIF


      RETURN
      END


	SUBROUTINE MDIAN3(X,Y,Z,N,XMED, YMED, Zmed)

        implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION X(N), Y(N), Z(N)

	if (n.gt.0) then
	call  MOMENT(X,N,XMED,SDEV)
	call  MOMENT(Y,N,YMED,SDEV)
	call  MOMENT(Z,N,ZMED,SDEV)
	return

	else


      N2=N/2
      CALL SORT3(N,X, Y, Z)

      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
        YMED=0.5*(Y(N2)+Y(N2+1))
	ZMED=0.5*(Z(N2)+Z(N2+1))
      ELSE
        XMED=X(N2+1)
        YMED=Y(N2+1)
	ZMED=Z(N2+1)
      ENDIF

	endif

      RETURN
      END


      SUBROUTINE IRASTAT(X,N,XQ1,XQ2,XQ3)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION X(N)

      CALL SORT(N,X)

c first compute 50% value

      N2=nint(N/2.0)
      XQ1 = X(N2)

c now compute 68%
      N2=nint(N * 0.68)
      XQ2 = X(N2)

c now compute 87%
      N2=nint(N * 0.87)
      XQ3 = X(N2)
  

      RETURN
      END


      SUBROUTINE QUART3(X,N,XQ)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION X(N)

      z = (3.*N)/4.
      N2= nint(z)

      CALL SORT(N,X)

      XQ = X(N2+1)

      RETURN
      END

      SUBROUTINE QUART(X,N,XQ)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION X(N)

      z = (1.*N)/4.
      N2= nint(z)

      CALL SORT(N,X)

      XQ = X(N2+1)

      RETURN
      END


	subroutine stats (nx,ny,iskip,il,ih,jl,jh,
     1    zl0,zh0,arr,nit,zmean,zstd,zmed)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (nz = 1000 * 1000,converge = 1.e-4)

	dimension arr(nx,ny)
	real*4 z(nz)

	nmin = 10 ** 2

c	write (6,*) 'here',il,ih,jl,jh,nit

c iterate

c	write (6,*) 'iter   mean       sigma       median'
c	write (6,*) il,ih,jl,jh
c	write (6,*) zl0,zh0,nit
	zmean=0.
	zstd=0.
	zmed=0.

	zmean_old = 0.
	do j=1,nit
		if (j.eq.1) then
			zl = zl0
			zh = zh0
		else
			fact = 3. 
			zl = (-fact * zstd) + zmed
			zh =  (fact * zstd) + zmed
		endif

		n = 0
		jd = 0
		do 57  k = jl,jh
			jd = jd + 1
			if (jd.eq.iskip) then
				jd=0
				goto 57
			endif	
			id=0
		do l = il,ih
			id = id + 1
			val = arr(l,k)
			if ( (val.ge.zl).and.(val.le.zh)) then
				if (id.lt.iskip) then
				 n = n + 1
				 if (n.gt.nz) n=nz
				 z(n) = val
				else
				 id = 0
				endif
			endif

		enddo
 57		continue


c		write (6,*)j,n,zmean,zstd,zmed

		if (n.gt.nmin) then
		  call MOMENT(z,n,zmean,zstd)
		  call MDIAN1(z,n,zmed)
		else
			zmean = -99.
			zmed = -99.
			zstd = -99.
			goto 99
		endif

c		write (6,*) zmean,zstd,zmed
c		write (6,*) 'res ',n,zl,zh,zmean,zstd,zmed

		ratio = (zmean - zmean_old) / zmean
c		if (abs(ratio).le.converge) goto 99
		zmean_old = zmean

  	enddo
  99 	i=1

	return
	end


	subroutine fast_stats (nx,ny,include,il,ih,jl,jh,
     1    zl0,zh0,arr,nit,zmean,zstd,zmed)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (nz = 500 * 1000)
	dimension arr(nx,ny)
	real*4 z(nz)

c iterate

c	write (6,*) 'iter   mean       sigma       median'
c	write (6,*) il,ih,jl,jh
c	write (6,*) zl0,zh0

	zmean_old = 0.
	do j=1,nit
		if (j.eq.1) then
			zl = zl0
			zh = zh0
		else
			fact = 3. 
			zl = (-fact * zstd) + zmed
			zh =  (fact * zstd) + zmed
		endif

		n = 0
		jd = 0
		do 577  k = jl,jh
			jd = jd + 1
			if (jd.lt.include) goto 577
			jd=0
			id=0
		do 57 l = il,ih
			id = id + 1
			if (id.lt.include) goto 57
			id = 0
			val = arr(l,k)
			if ((val.ge.zl).and.(val.le.zh)) then
				n = n + 1
				if (n.gt.nz) n=nz
				z(n) = val
			endif
 57		continue
 577		continue

c		write (6,*)n,zmean,zstd,zmed

		if (n.gt.12) then
		  call MOMENT(z,n,zmean,zstd)
c		  call MDIAN1(z,n,zmed)
		  zzlow = zmean - 10.
		  zzhigh = zmean + 10.
		  width = 0.01
	 	  z16 = 0.
		  z84 = 0.
		  call bin_stats (zzlow,zzhigh,width,z,N,zmed,hsig, z16, z84)
		else
			zmean = -99.
			zmed = -99.
			zstd = -99.
			goto 99
		endif

c		write (6,*) zmean,zstd,zmed
c		write (6,*) n,zl,zh,zmean,zstd,zmed

c		ratio = (zmean - zmean_old) / zmean
c		if (abs(ratio).le.converge) goto 99
c		zmean_old = zmean

  	enddo
 99	i=1
	return
	end


	subroutine bin_stats (blow,bhigh,bwidth,data,N,zmedian,hsig, z16,z84)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (max = 25000, zfloor = -947.)
	real*4 data(N), data2(max), modsig
	integer count (max)


	ntry = 0
	z50 = 0.
	z16 = 0.
	z84 = 0.
	hsig = 0.
	znorm = 0.
	zmode = 0.
	modsig = 0.

 47	do j=1,max
		count(j) = 0
	enddo

	ntot = 0
	iff = 0
	modedex = 0
	modemax = 0
	do 50 j=1,N
		val = data(j)

		if (val.eq.zfloor) iff=1
		if (val.lt.blow) goto 50
		if (val.gt.bhigh) goto 50

		dval = val - blow
		idex = 1 + int(dval / bwidth)	

c	write (6,*) val,blow,bhigh, idex

		if (idex.gt.max) goto 50

		ntot = ntot + 1
		count(idex) = count(idex) + 1

		if (count(idex).gt.modemax) then
			modemax = count(idex)
		 	modedex = idex
		endif

 50	continue

c	write (6,*) 'binstats ',ntot,bwidth

	zmode = -9999.
	if (modedex.gt.0) then
	  zmode =  blow + (modedex * bwidth) - (bwidth / 2.0)
	endif

c	write (6,*) 'zmode = ',bwidth,zmode
	zhigh = zmode + (bwidth /2.)

	m = 0
	do  j=1,N
                val = data(j)
		if ( (val.ge.-500.).and.(val.le.zhigh)) then
		   m=m+1
		   data2(m) = zmode - val
c		   m = m + 1
c		   data2(m) = val - zmode	
		endif
	enddo

	zAVE = 0.
	zSDEV = 0.
	if (m.gt.0) then
	  call MOMENT(data2,m,zAVE,zSDEV)
	  call ModMOMENT(data2,m,zSDEV)
	  modsig = zSDEV
	endif

c	write (6,*) 'zmode stats ',m,modsig

cc now consider the neighbors, do a intensity-weighted moment
	if (modedex.gt.2) then
	 iC1 = count(modedex-2) 
	 iC2 = count(modedex-1)
	 iC3 = count(modedex+1)
         iC4 = count(modedex+2)

c	 write (6,*) iC1,iC2,modemax,iC3,iC4
	 itot = iC1 + iC2 + modemax + iC3 + iC4
	 t1 = -2. * bwidth * iC1 
	 t2 = -1. * bwidth * iC2 
	 t3 =  1. * bwidth * iC3 
         t4 =  2. * bwidth * iC4 
	 zmoment = (t1 + t2 + t3 + t4) / (itot * 1.)
c	 write (6,*) 'moment: ',zmoment
	 zmode = zmode + zmoment
c	 write (6,*) 'new zmode = ',zmode
	endif


	if (modsig.gt.0.) then
		zmedian = zmode
		hsig = modsig
		return  
	endif

	

	if (ntot.lt.9) then
		ntry = ntry + 1
		if (ntry.le.5) then
c	write (6,*) 'nstry ',ntry,bwidth
			bwidth = bwidth * 1.5
			goto 47
		endif

		zmedian = -9999.
		if (iff.eq.1) zmedian = zfloor
		goto 66
	endif

	i50 = nint( (ntot*1.) / 2.0)
	i84 = nint( 0.84 * ntot * 1. )
	i16 = nint( 0.16 * ntot * 1. )
	znorm = sqrt (ntot*1.)


	isum = 0
        isum_prev = 0
	zmedian = -9999.
	z84  = -9999.
	z16 = -9999.


	do j=1,max

		isum = isum + count(j)

		vv =  blow + (j * bwidth) - (bwidth / 2.0)

		if (isum.eq.i50) zmedian = vv
		if (isum.eq.i16) z16 = vv
		if (isum.eq.i84) z84 = vv

c	if (count(j).gt.0) write (120,*) j,vv,vv-zmode,count(j)

		if ((isum.gt.i50).and.(zmedian.lt.-999)) then
 		  ! interpolate
		  idel = isum - isum_prev
		  ilow = i50 - isum_prev
		  ratio = (ilow*1.) / (idel*1.)
		  del = bwidth * ratio
		  zmedian = (del + blow) + ((j-1) * bwidth) - (bwidth / 2.0)
		endif

		if ((isum.gt.i16).and.(z16.lt.-999)) then
                  ! interpolate
                  idel = isum - isum_prev
                  ilow = i16 - isum_prev
                  ratio = (ilow*1.) / (idel*1.)
                  del = bwidth * ratio
                  z16 = (del + blow) + ((j-1) * bwidth) - (bwidth / 2.0)
                endif

		if ((isum.gt.i84).and.(z84.lt.-999)) then
                  ! interpolate
                  idel = isum - isum_prev
                  ilow = i84 - isum_prev
                  ratio = (ilow*1.) / (idel*1.)
                  del = bwidth * ratio
                  z84 = (del + blow) + ((j-1) * bwidth) - (bwidth / 2.0)
                endif

		isum_prev = isum
	enddo
 51 	i=1
 66	i=1

	hsig = (z84 - z16)/2.0   ! rms of the pixel distribution
c	hsig = hsig / znorm  ! rms in the mean (median)

c	write (88,*) '# ',bwidth,zmedian
c	close (88)

c	close (120)

	return
	end			



c  weighted mean
      SUBROUTINE WMOMENT(DATA,w,N,AVE,SDEV,sdev2)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION DATA(N),w(N)
      real*8 s,var,var2,p,t,t2

      IF(N.eq.1) then
                AVE = DATA(1)
                SDEV = 0.
                goto 99
      else IF(N.Lt.1)then
c                write (6,*) 'error N < 1'
                goto 99
      endif

	wt = 0.
	wt2 = 0.
	S=0.d0
	do j=1,N
	  wt = wt + w(j)
	  wt2 = wt2 + (w(j)**2)

	  S=S+(DATA(J)*w(j)*1.d0)

	enddo


	AVE=sngl(S / wt*1.d0 )

       VAR=0.d0
	sum = 0.  
       DO 12 J=1,N
        S=(DATA(J)-ave)*1.d0
        P=S*S
        sum=sum+(P*w(j))
12    CONTINUE

	Var = Sum / (wt*1.d0)  !  biased weighted sample variance

	Var2 = Sum * (wt / ((wt**2) - wt2) * 1.d0)

      t = var ** 0.5
      t2 = var2 ** 0.5

      SDEV=sngl(t)
      sdev2 = sngl(t2)

 99     i=1
      RETURN
      END


      SUBROUTINE RMS(DATA,N,AVE)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION DATA(N)
      real*8 s,var,p,t

      IF(N.eq.1) then
                AVE = DATA(1)
                goto 99
      else IF(N.Lt.1)then
c                write (6,*) 'error N < 1'
                goto 99
      endif

      S=0.d0
      DO 11 J=1,N
        S=S+(DATA(J)**2.d0)
11    CONTINUE

      VAR=S/(N*1.d0)
      t = var ** 0.5
      ave=sngl(t)

 99     i=1
      RETURN
      END



c weighted median

	subroutine minabs( n, f, g, K, w, alpha)

       implicit integer (i-n)
       implicit real*4 (a-h)
       implicit real*4 (o-z)

	integer i,n,low,large,ml,mh,mlt,mht,itry,j
	real*4 gn,gp,t,alpha,er,gnt,gpt,gplx,gmix,small,grad
        real*4    w(n),f(n), g(n)
	integer k(n)
	logical debug

c	debug = .true.
	debug = .false.    !  JWF B21116

	small = 0.
	do i=1,n
		k(i) = i
	enddo

	low = 1
	large = n

	if (debug) write (6,*) 'Weighted Median'
	if (debug) write (6,*) low,large
	ml = n
	mh = 1
	gn = 0.
	gp = 0.
        t = 0
        grad = 0
        gmix = 0
	do itry = 1,n 


	 if (debug) write (6,*) (large-low)/3+itry, large-low+1

	if (debug) write (6,*)mod( (large-low)/3+itry, large-low+1)
	if (debug) write (6,*) low,low+mod( (large-low)/3+itry, large-low+1) 

		j = k( low+mod( (large-low)/3+itry, large-low+1) )

	 if (debug)  write (6,*) 'j = ',j
	if (j.gt.n) goto 47


		if( abs(w(j)) .ne. 0.) then

			t = f(j)/(w(j))
			f(j) = w(j)*t

			do i = low,large 
			  j = k(i)
			  er = f(j)-w(j)*t
			  g(j) = 0.
			  if    ( er.gt.small) then
					g(j) = -w(j)
			  else if ( er.lt.-small) then
					g(j) = +w(j)
			  else
					g(j) = 0.
			  endif

			enddo

	if (debug) write (6,*) 'call split'

			call split(low,large,k,g,mlt,mht)
	if (debug) write (6,*) 'done',mlt,mht

			gnt = gn
			do i = low,mlt
			   gnt = gnt+g(k(i))
			enddo
			gpt = gp
			do i = mht,large
				gpt = gpt+g(k(i))
			enddo
			gplx = 0.
			gmix = 0.
			do i = mlt,mht 
			     j = k(i)
			     if (w(j).lt.0.) then
				gplx = gplx-w(j)
				gmix = gmix+w(j)
			     endif
			if (w(j).gt.0.) then
				gplx = gplx+w(j)
				gmix = gmix-w(j)
			endif

			enddo

			grad = gnt+gpt
			if ((grad+gplx)*(grad+gmix).lt.0.) then
				goto 47
			endif

			if (grad.ge.0.) low = mht+1
			if (grad.le.0.) large = mlt-1
			if (low.gt.large) goto 47
			if (grad.ge.0.) gn = gnt+gmix
			if (grad.le.0.) gp = gpt+gplx
			if (grad+gplx.eq.0.) ml = mlt
	endif

 47	if (grad+gmix.eq.0.) mh = mht

	enddo

	alpha = - t

	return
	end


	subroutine split(low,large,k,g,ml,mh)

c	given g(k(i)),i=low,large
c	then rearrange k(i),i=low,large and find ml and mh so that
c	(g(k(i)),i=low,(ml-1)) < 0. and
c	(g(k(i)),i=ml,mh)=0. and
c	(g(k(i)),i=(mh+1),large) > 0.

	implicit integer (i-n)
       implicit real*4 (a-h)
       implicit real*4 (o-z)

	real*4    g(1)
	integer k(1)
	integer low,large,ml,mh,keep,i,ii

	ml = low
	mh = large
	do 500 kkk=1,99999

	  ml = ml-1

	  do 501 kkkk=1,9999
		ml = ml+1
		if (g(k(ml)).ge.0.) goto 502
 501	  continue

 502	  do 503 kkkk=1,99999
		mh = mh+1

		do 504 kkkkk=1,99999
			mh = mh-1
			if (g(k(mh)).le.0.) goto 505
 504		continue

 505		keep = k(mh)
		k(mh) = k(ml)
		k(ml) = keep

		if (g(k(ml)).ne.g(k(mh))) goto 500 ! Break out of enclosing "repeat{"
			
		do i = ml,mh 
			ii = i
			if (g(k(i)).ne.0.0) goto 30
		enddo

	 	return	! Break out of two enclosing "repeat{"'s

 30  		keep = k(mh)

		k(mh) = k(ii)
		k(ii) = keep

 503		continue

 500	continue

	return
	end



	subroutine bincounts (blow,bhigh,bwidth,data,N,max,count)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (zfloor = -947.)
	real*4 data(N)
	integer count (max)

	do j=1,max
		count(j) = 0
	enddo

	ntot = 0
	iff = 0
	do 50 j=1,N
		val = data(j)
		if (val.eq.zfloor) iff=1
		if (val.lt.blow) goto 50
		if (val.gt.bhigh) goto 50

		dval = val - blow
		idex = 1 + int(dval / bwidth)	
		if (idex.gt.max) goto 50

		ntot = ntot + 1
		count(idex) = count(idex) + 1

c	if (idex.eq.24) write (67,*) val

 50	continue


	return
	end			


	subroutine iSTATS (data,ztmp,n,zlow,zhigh,iter,ave,sdev,xmed,m)

        implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 data(n),ztmp(n)

	ave = -9.
	xmed = -9.
	sdev = 0.
	m=0

	do K = 1,iter

c	write (6,*) K,ave,sdev,xmed

	   m=0
	   do J = 1,n
		if ((data(J).gt.zlow).and.(data(J).le.zhigh)) then
		   m=m+1
		   ztmp (m) = data(J)
		endif
	   enddo

	   if (m.lt.5) then
		return
	   endif

	   call MOMENT(ztmp,m,AVE,SDEV)
	   call MDIAN1(ztmp,m,XMED)
 
	   zlow = xmed - (3. * sdev)
	   zhigh = xmed + (3. * sdev)

	enddo


	return
	end

