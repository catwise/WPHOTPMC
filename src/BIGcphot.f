	subroutine BIGcphot (Jsrc,ib,ncoaddsize,nx,ny,sarray,Unc,Cov, CoMSK, iCOADD,
     1     Nap,R,zmed,csig,
     1     xc,yc,ntot,ztot,xint, flag, iflag, SUM2err, fbits, Lconf, cov_ave, cov_std)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (pi=3.1415926536, valmin = -500., valmax = 1.e7) 
	real*4  sarray(ncoaddsize,ncoaddsize,4), Unc(ncoaddsize,ncoaddsize,4)
	real*4  cov(ncoaddsize,ncoaddsize,4)
	real*4 R(19), xint(19),ztot(19), SUM2err(19), Lconf
	real*4 zcov (9999)
	integer flag(19), fatal,ibit,BitSet, fbits(32)
	integer iflag(19),ntot(19)
	integer*2 CoMSK(ncoaddsize,ncoaddsize,4)
	integer iCOADD(ncoaddsize,ncoaddsize,4), imask
	integer trigger(19), satrigger(19)

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

	do Map = 1,nap	
	   flag(Map) = 0
	   iflag(Map) = 0
	   ntot(Map) = 0
	   xint(Map) = 0.
	   ztot(Map) = 0.
	   SUM2err(Map) = 0.
	   trigger(Map) = 0
	  satrigger(Map) = 0
	enddo

	iflag = 0

	i0 = nint(xc)
	j0 = nint(yc)

	Rmax = 0.
	do K=1,nap
	  Rmax = max (Rmax, R(K) )
	enddo

	Imax = nint(Rmax)

        jstrt=j0-Imax-2
	jstrt = max(jstrt,1)
    	jend=j0+Imax+2
	jend=min(jend,ny)
    	istrt=i0-Imax-2
	istrt = max(istrt,1)
    	iend=i0+Imax+2
	iend=min(iend,nx)

	rat = 1.
	phimax = 0.

	sum_cov = 0.
	nsum_cov = 0

c	if ( (Jsrc.eq.13).and.(ib.eq.1)) write (6,*) 'DEBUG'

c	write (6,*) istrt,iend,jstrt,jend, zmed
c	if (ib.eq.1) write (6,*) 'BIGCphot ',R(1)

	do 100 JJ=jstrt,jend
	      	dy=(jj*1.)-yc
	do 101 II=istrt,iend
		dx=(ii*1.)-xc

		dr = sqrt ( (dx**2) + (dy**2) )
		if (dr.gt.Rmax+1.20) goto 101  ! end

c	write (6,*) ii,jj,ib

		z = sarray (ii,jj,ib) - zmed  ! subtract the background

		do Map = 1, nap

			Rap = R(Map)

			if (dr.le.Rap+1.20)  then  ! inside aperture or just outside boundary


			! coverage for the standard aperture
			  if (MAP.eq.1) then
				nsum_cov = nsum_cov + 1
				if (nsum_cov.lt.99999) zcov(nsum_cov)= Cov (ii,jj,ib)  ! THJ 07Jan2018
			  endif

			  del = dr - Rap
			  if (abs(del).le.1.25) then   ! near edge
                            !  check frac pixel
			    call fracpix (xc,yc,ii,jj,Rap,rat,phimax,farea)
                          else
                                  farea = 1.0
                          endif

ccccc Check the iMask for bad pixels, etc
cccccc   fatal == bits 2 or 18

cvalue    Condition
c-----    ------------------------------------------------------
c  0      nominal
c  1      source confusion
c  2      bad or fatal pixels:  bit 2 or 18
c  4      non-zero bit flag tripped (other than 2 or 18)
c  8      corruption
c 16      saturation
c 32      upper limit

			if (del.le.0) then  ! within circ aperture

c			   MSK(ii,jj,ib) = 0   ! temporary -- disable the MSK since it is broken

			    imask = iCOADD  (ii,jj,ib)
			    if ( (imask.ne.0).and.(imask.ne.Jsrc) ) then
				flag(MAP) = 1  ! star contamination
			    endif


c	if (ib.eq.4) write (6,*) ii,jj,Cov (ii,jj,ib),CoMSK(ii,jj,ib)

c	if ( (Jsrc.eq.13).and.(ib.eq.1).and.(MAP.eq.1) ) write (6,*) ii,jj,Cov (ii,jj,ib),CoMSK(ii,jj,ib),iflag(MAP)

			    if ((Cov (ii,jj,ib).lt.1.1).and.(trigger(map).eq.0)) then
				if (Cov (ii,jj,ib).le.0.) then
					iflag(MAP) = 2 + iflag(MAP)    ! poor coverage
					trigger(map) = 1
				else if (CoMSK(ii,jj,ib).gt.0) then 
                                  iflag(MAP) = 2 +  iflag(MAP)  !  bad pixels with low coverage
				  trigger(map) = 1
                                endif
			    endif

			    if ((CoMSK(ii,jj,ib).eq.100).and.(satrigger(map).eq.0)) then 
				iflag(MAP) = 16 + iflag(MAP)  !  saturation
				Satrigger(map) = 1
			    endif

c			   if (MSK(ii,jj,ib).gt.0) iflag(Map) = 3
c			   icheck = 0
c                           icheck = iand(MSK(ii,jj,ib), fatal)
c                           if (icheck.gt.0) then
c				iflag(Map) = 2   !   fatal
c			   endif
c

                        endif


			  errval = Unc (ii,jj,ib)

c  RSS the confusion noise into the UNC value

		          if ((z.gt.Valmin).and.(z.lt.Valmax)) then

		       		xint(Map) = xint(Map) + (z*farea)   ! integrated flux
				ztot(Map) = ztot(Map) + farea       ! total area
		       		ntot(Map)=ntot(Map)+nint(farea)
			
				if ((errval.ge.0.).and.(errval.lt.Valmax)) then


				  SUM2err(Map) = SUM2err(Map) + (errval**2) + (Lconf**2)  ! this is the poisson error due to the source itself

c	write (6,*) 'errval = ',ii,jj,ib,errval, SUM2err(Map)

				else
				
				  iflag(Map) = max(2,iflag(Map))   ! bad news w/ UNC

				endif

			  else
				iflag(Map) = max(2,iflag(Map))    ! bad news
			  endif

		      endif

		enddo  ! MAP



 101		continue
 100	  continue
 

 102	  do Map = 1, nap
		ntot(Map) = nint(ztot(Map))
	  enddo

c	write (6,*) 'down here ',nsum_cov
	if (nsum_cov.gt.1.) then

	   call MOMENT(zcov,nsum_cov,cov_ave, cov_std)

	else
		cov_ave = 0.
		cov_std = 0.
	endif

c	write (6,*) 'SUM2err(1) = ',SUM2err(1)


	return
	end

