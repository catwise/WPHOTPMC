	SUBROUTINE twod_ellip_new (nx,ny,iso,nmax,bamin,
     1     rmin,thetamin,fracmin,coadd)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (mmax=750,rsemi_min = 4.5,maxper = 2000)
	integer*2 iso(nmax,nmax)
	real*4 perimeter(maxper,2),stats(maxper)
	real*4 stats2(maxper)
	character*(*) coadd
c	real*4 larray (mmax*mmax)
	real*4 bog(mmax,mmax),min_ba
	character*72 ff
	logical debug

	debug = .false.

	do j=1,maxper
		perimeter(j,1) = 0
		perimeter(j,2) = 0
		stats(j) = 0
	enddo
c
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
c	write (6,*) 'center of iso ',nx,ny,idex0,jdex0

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
c		perimeter(nper,1) = dx-0.5
c		perimeter(nper,2) = dy-0.5

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
		write (300,*) perimeter(jj,1),perimeter(jj,2)
	 enddo

	endif

c  axil ratio

	min_ba = 0.1   ! minimum b/a that we will solve for
	dba = 0.020     ! step size
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
c		  drlow = Rfid - (3.0*drsig)
c		  drh = Rfid + (1.5*drsig)

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
	enddo


	if (debug) then

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

c	 ic = 0
cc        do jj=1,ny*2
c	 do jj=1,ny*3
c         do ii=1,nx
c             ic = ic + 1
c             larray(ic) = bog (ii,jj)
c         enddo
c         enddo
c	 ff = 'peri'
c	 write (6,*) ' '
c	 write (6,*) 'isophot mask & best fit ellipse mask: peri.fits'
c	 call writeimage (nx,ny*3,larray,ic,coadd,ff,ii,jj,crv1,crv2)

	endif

 47	ii=1

	return
	end

