	subroutine GALcphot (Jsrc,ib,ncoaddsize,nx,ny,sarray,Unc, Cov, CoMSK,iCOADD,
     1     Nap,R, rat, phimax, zmed,csig,
     1     xc,yc,ntot,ztot,xint, flag, iflag, SUM2err, fbits, Lconf)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (pi=3.1415926536, valmin = -500., valmax = 1.e7) 
	real*4  sarray(ncoaddsize,ncoaddsize,4), Unc(ncoaddsize,ncoaddsize,4)
	real*4  Cov(ncoaddsize,ncoaddsize,4)
	real*4 R, xint,ztot, SUM2err, Lconf
	integer flag, fatal,ibit,BitSet, fbits(32)
	integer iflag,ntot
c	integer MSK(ncoaddsize,ncoaddsize,4)
	integer*2 CoMSK(ncoaddsize,ncoaddsize,4)
	integer iCOADD(ncoaddsize,ncoaddsize,4), imask
	integer trigger, satrigger


c	fatal = 2**2 + 2**18            ! bits 2 and 18 are regarded as fatal
	fatal = 0
        do I = 1,32
                ibit = fbits(I)
                if (ibit.gt.0) then
                  fatal = fatal + (2**ibit)
		else
                        !quit
                        goto 74
                endif
        enddo
 74     I = 0


	flag = 0
	iflag = 0
	ntot = 0
	xint = 0.
	ztot = 0.
	SUM2err = 0.
        trigger = 0
	satrigger = 0

	i0 = nint(xc)
	j0 = nint(yc)

	Rmax = R
	Rap = R
	Imax = nint(Rmax)

        jstrt=j0-Imax-2
	jstrt = max(jstrt,1)
    	jend=j0+Imax+2
	jend=min(jend,ny)
    	istrt=i0-Imax-2
	istrt = max(istrt,1)
    	iend=i0+Imax+2
	iend=min(iend,nx)

c	write (6,*) istrt,iend,jstrt,jend

	do 100 JJ=jstrt,jend
	      	dy=(jj*1.)-yc
	do 101 II=istrt,iend
		dx=(ii*1.)-xc

		call w_ell (rat,phimax,dx,dy,dr)
		
		if (dr.gt.Rmax+1.20) goto 101  ! end

c	write (6,*) ii,jj,ib

		z = sarray (ii,jj,ib) - zmed  ! subtract the background

		if (dr.le.R+1.20)  then  ! inside aperture or just outside boundary

		del = dr - R
		if (abs(del).le.1.25) then   ! near edge
                        !  check frac pixel
		    call fracpix (xc,yc,ii,jj,R,rat,phimax,farea)
                else
                    farea = 1.0
                endif

ccccc Check the iMask for bad pixels, etc
cccccc   fatal == bits 2 or 18

		if (del.le.0) then  ! within  aperture

c			   MSK(ii,jj,ib) = 0   ! temporary -- disable the MSK since it is broken

		    imask = iCOADD  (ii,jj,ib)
		    if ( (imask.ne.0).and.(imask.ne.Jsrc) ) then
			flag = 1  ! star contamination
		    endif


		    if ((Cov (ii,jj,ib).lt.1.1).and.(trigger.eq.0)) then
                         if (Cov (ii,jj,ib).le.0.) then
                                 iflag = 2 + iflag   ! poor coverage
                                 trigger = 1
                         else if (CoMSK(ii,jj,ib).gt.0) then
                           iflag = 2 +  iflag  !  bad pixels with low coverage
                           trigger = 1
                         endif
                    endif

                    if ((CoMSK(ii,jj,ib).eq.100).and.(satrigger.eq.0)) then
                         iflag= 16 + iflag  !  saturation
                         Satrigger = 1
                    endif




                endif


		errval = Unc (ii,jj,ib)

c  RSS the confusion noise into the UNC value

		if ((z.gt.Valmin).and.(z.lt.Valmax)) then

		  	xint = xint + (z*farea)   ! integrated flux
			ztot = ztot + farea       ! total area
		       	ntot=ntot+nint(farea)
			
			if ((errval.ge.0.).and.(errval.lt.Valmax)) then

			  SUM2err = SUM2err + (errval**2) + (Lconf**2)  ! this is the poisson error due to the source itself

			else
				
			  iflag = max(2,iflag)    ! bad news w/ UNC

			endif

		else
			iflag =  max(2,iflag)     ! bad news
		endif

		endif


 101		continue
 100	  continue

	  ntot = nint(ztot)

c	write (6,*) R, rat, phimax
c	write (6,*) ztot, xint

	return
	end

