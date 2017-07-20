	subroutine twod_ellip2 (nx,ny,iso,bamin,rmin,thetamin,coadd)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (mmax=100,rsemi_min = 2.,maxper = 1000,lsize = 100*100)

	integer*2 iso(100,100)
	real*4 perimeter(maxper,2),stats(maxper)
	real*4 stats2(maxper)
	character*(*) coadd
	real*4 larray (mmax*mmax)
	real*4 bog(mmax,mmax),min_ba
	character*72 ff
	logical debug

	debug = .false.
c	debug = .true.

	do j=1,maxper
		perimeter(j,1) = 0
		perimeter(j,2) = 0
		stats(j) = 0
	enddo

c clean up mask
        do    j=2,ny-1
        do 90 i=2,nx-1
                ival = iso(i,j)
                if (ival.le.1) goto 90

		n = 0

		iv = iso(i-1,j)
		if (iv.eq.1) n=n+1
		iv = iso(i+1,j)
                if (iv.eq.1) n=n+1
		iv = iso(i,j-1)
                if (iv.eq.1) n=n+1
		iv = iso(i,j+1)
                if (iv.eq.1) n=n+1

		if (n.lt.1) iso(i,j) = 0

 90	continue
	enddo


cccc

	if (debug) then
	do j=1,mmax
         do i=1,mmax
                bog(i,j) = 0.
         enddo
         enddo


	jj= 0
        do j=ny+1,ny*2
		jj=jj+1
        do i=1,nx
		bog(i,j) = iso(i,jj)
	enddo
	enddo

	endif

	idex0 = (nx / 2) + 1
	jdex0 = (ny / 2) + 1

	nper = 0
	n=0
        do    j=2,ny-1
        do 50 i=2,nx-1
                ival = iso(i,j)
                if (ival.le.1) goto 50

c ok, we have a perimeter point
		if (nper.eq.maxper) goto 50

		dx = i - idex0
		dy = j - jdex0 
		nper = nper + 1
		perimeter(nper,1) = dx
		perimeter(nper,2) = dy

		dr = sqrt (  (dx**2) + (dy**2) )
c		n=n+1
c		if (n.le.500) z(n) = dr
 50	continue
	enddo


	if (nper.le.5) goto 47

	if (debug) then
  	 write (6,*) nper,' perimeter pixels'

	 do jj=1,nper
		write (30,*) perimeter(jj,1),perimeter(jj,2)
	 enddo

	endif

c  axil ratio

	min_ba = 0.2   ! minimum b/a that we will solve for
	dba = 0.01     ! step size
	ndba = 1 + ((1.0 - min_ba) / dba)
	
	ba = min_ba - dba

	rmin = 0.
	drsigmin = 999.
	thetamin = 0.
	fracmin = 1.0
	bamin = 0.

	do L = 1,ndba
		ba = ba + dba

c	write (6,*) ba

		dtheta = 1.0
		ntheta = (180.0/dtheta) 	
		theta = -90 - dtheta

		do M = 1,ntheta	
		  theta = theta + dtheta
	
		  do k = 1,nper
		    dx = perimeter(k,1)
		    dy = perimeter(k,2)
		
		    call w_ell (ba,theta,dx,dy,dr)
		    stats(k) = dr
		  enddo

		  call MOMENT (stats,nper,drmean,drsig)

		  nper2 = 0
		  drlow = drmean - (3.0*drsig)
		  drh =   drmean + (3.0*drsig)

		  do LL = 1,nper
			val = stats(LL)
			if ( (val.gt.drlow).and.(val.lt.drh)) then
				nper2 = nper2 + 1
				stats2(nper2) = val
			endif
		  enddo

		  if (nper2.gt.3) then
			call MOMENT (stats2,nper2,drmean,drsig)
		  else
			drsig = 999.
		  endif

c percentage of mean
		  if (drmean.gt.0.) then
		    frac = drsig / drmean
		  else
		    frac = 1.0
		  endif

		  if (frac.lt.fracmin) then
			rmin = drmean
			drsigmin = drsig
			thetamin = theta
			fracmin = frac
			bamin = ba
		  endif
		enddo

c	write (6,*) bamin,thetamin
	
	enddo

	jj= jdex0 + (ny*2)
	do j=(ny*2)+1 , ny*3
	do i=1,nx
		dx = i - idex0
		dy = j - jj
		call w_ell (bamin,thetamin,dx,dy,dr)

		dd = abs(dr - rmin)
		if ((dr.le.rmin).and.(debug)) bog(i,j) = 100
	enddo
	enddo	

	if (debug) then
	 ic = 0
	 do jj=1,ny*3
         do ii=1,nx
             ic = ic + 1
             larray(ic) = bog (ii,jj)
         enddo
         enddo
	 ff = 'peri.fits'
	 write (6,*) ' '
	 write (6,*) 'isophot mask & best fit ellipse mask: peri.fits'
c	call wimage (nx,ny*3,larray,lsize,coadd,ff,iref,jref,ra,dec)

	 nnn = ny*3
	 call wimage  (nx,nnn,lsize,larray,coadd,ff)

	endif

 47	ii=1

c	if (bamin.gt.0.85) thetamin=0.
c	if (rmin.le.3.0) thetamin=0.

	return
	end

