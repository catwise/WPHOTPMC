	subroutine cphot (Jsrc,J,nf,nfpix,ib,nnx,nny,nx,ny,sarray,MASK,iMask,Unc,
     1     Nap,R,zmed,csig,
     1     xc,yc,ntot,ztot,xint, flag, iflag, SUM2err, SPIT, MIPS, fbits )

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (pi=3.1415926536, valmin = -500., valmax = 1.e7) 
	real*4  sarray(nnx,nny,nfpix), Unc(nnx,nny,nfpix)
	real*4 R(19), xint(19),ztot(19), SUM2err(19)
	integer flag(19), iMask(nnx,nny,nfpix),fatal,ibit,BitSet, fbits(32)
	integer iflag(19),ntot(19)
c	integer MASK(nnx,nny,nf,4), IImask
	integer MASK(1,1,1) ! Mask being dummied out since it's never assigned -- TPC
        integer IImask
	logical SPIT, MIPS

c	fatal = 2**2 + 2**18            ! bits 2 and 18 are regarded as fatal
c	if (SPIT) then
c		fatal = (2**1) + (2**8) + (2**9) + (2**14)
c		! MIPS >>   fatal  = (2**10) + (2**11)
c		if ((MIPS).and.(ib.eq.4)) fatal  = (2**1) + (2**8) + (2**9) + (2**14) + (2**10) + (2**11)
c	endif

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


	do 100 JJ=jstrt,jend
	      	dy=(jj*1.)-yc
	do 101 II=istrt,iend
		dx=(ii*1.)-xc

		dr = sqrt ( (dx**2) + (dy**2) )
		if (dr.gt.Rmax+1.20) goto 101  ! end

		z = sarray (ii,jj,J) - zmed  ! subtract the background

		do Map = 1, nap

			Rap = R(Map)

			if (dr.le.Rap+1.20)  then  ! inside aperture or just outside boundary

			  del = dr - Rap
			  if (abs(del).le.1.25) then   ! near edge
                            !  check frac pixel
			    call fracpix (xc,yc,ii,jj,Rap,rat,phimax,farea)
                          else
                                  farea = 1.0
                          endif

ccccc Check the iMask for bad pixels, etc
cccccc   fatal == bits 2 or 18

			  icheck = 0
		  	  if (del.le.0) then  ! within circ aperture

c			    IImask = MASK  (ii,jj,J) ! Not used
			    IIMask = 0
			    if ( (IImask.ne.0).and.(IImask.ne.Jsrc) ) then
				flag(MAP) = 1  ! star contamination
			    endif

			    iMaskval = iMask(ii,jj,J)
			    if (iMaskval.gt.0) iflag(Map) = 3

			    icheck = iand(iMaskval, fatal)
			    if (icheck.gt.0) iflag(Map) = 2   !   bit 2 or 18 set

			  endif

			  errval = Unc (ii,jj,J)

		          if (icheck.eq.0.and.(z.gt.Valmin).and.(z.lt.Valmax).and.(errval.ge.0.).and.(errval.lt.Valmax)) then

		       		xint(Map) = xint(Map) + (z*farea)   ! integrated flux
				ztot(Map) = ztot(Map) + farea       ! total area
		       		ntot(Map) = ntot(Map)+nint(farea)

				if ((errval.ge.0.).and.(errval.lt.Valmax)) then
					SUM2err(Map) = SUM2err(Map) + (errval**2)  ! this is the poisson error due to the source itself
				else
					 iflag(Map) = 2   ! bad UNC
				endif

			  else
				iflag(Map) = 2   ! bad news
			  endif

		      endif

		enddo  ! MAP



 101		continue
 100	  continue
 

 102	  do Map = 1, nap
		ntot(Map) = nint(ztot(Map))
	  enddo

	return
	end

