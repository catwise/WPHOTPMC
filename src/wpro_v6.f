c
c Counters and dimensions:
c
c nnx/y = allocated size of pixel arrays
c nx/yfp = mappsf dimension = frame dimensions
c nx/y = pixel array area loaded with data
c nx/yorig = frame dimensions
c nf = frame count (for array allocations except pixels)
c nfpix = pixel array allocation
c npsize = allocated size of each PSF dimaension
c npnmax = number of PSFs per band
c nall = allocated source count for results for whole tile
c nmax/nsrc = no. active sources for this region
c nsrc_adb = total count of active deblends
c nsrc_out = total output source count
c
c "Mask" is iMask in caller
c
	subroutine wpro_v6 (fbits,nnx,nny,nxfp,nyfp,nx,ny, nxorig, nyorig,
     1      Array,Unc,Mask,pix_order,nf,nfpix,
     1	    wflag,psfs,psfuncs,MapPSFs,npsize,npnmax,npsf,pxcent,pycent,ppix,
     1	    fwhm,nmax,nsrc,IDlist,HDR,RAlist,DEClist,Xpos,Ypos,xTRANS,yTRANS,JD,
     1	    Rstann,Rstwid,LBACK, LSIG, Lconf, MAGSTD, eMAGSTD,FLGSTD, nall,
     1	    Rsat, SatNum,zero, adb_nmax,
     1	    Table, framefluxes,frameuncs,framechis,nsrc_out, nsrc_adb, istat, SPIT, MIPS,
     +      nf0,basename,order)

!-----------------------------------------------------------------------------
! Do profile-fitting photometry on WISE images.  Results are output in the
! form of an array named "Table", arranged as follows:
!	Table(n,1)	=	RA [deg] (single precision version)
!	Table(n,2)	=	Dec [deg] (single precision version)
!	Table(n,3)	=	uncertainty in RA [arcsec]
!	Table(n,4)	=	uncertainty in Dec [arcsec]
!	Table(n,5)	=	cross term (sigxy) in RA,Dec [arcsec]
!	Table(n,6:9)	=	fluxes in bands 1,...4 [dn]
!	Table(n,10:13)	=	uncert. in fluxes in bands 1,...4 [dn]
!	Table(n,14:17)	=	reduced chi squared in bands 1,...4
!	Table(n,18)	=	overall reduced chi squared
!	Table(n,19)	=	blend number
!	Table(n,20)	=	number of actively deblended components
!	Table(n,21:24)	=	fraction of pixels affected by tempcal
!				masking in bands 1,...4
!	Table(n,25:28)	=	fraction of saturated pixels
!				in bands 1,...4                          non-PM obs epoch
!	Table(n,55)	=	nIters (no. of iterations in fit)
!	Table(n,56)	=	nSteps (no. of steps on last iteration)
!       ---------------------------------------------------------------------------------
! 	Table(n,29)	=	Proper motion in RA [arcsec/yr]          PM parameters
! 	Table(n,30)	=	Proper motion in Dec [arcsec/yr]
!	Table(n,31)	=	Uncertainty in RA proper motion ["/yr]
!	Table(n,32)	=	Uncertainty in Dec proper motion ["/yr]
!	Table(n,33)	=	int(JD0)  [mean obs epoch] ! JWF B21225 from here down
!	Table(n,34)	=	frac(JD0) [mean obs epoch]
!	Table(n,35)	=	uncertainty in RA [arcsec]
!	Table(n,36)	=	uncertainty in Dec [arcsec]
!     	Table(n,37)	=	cross term (sigxy) in RA,Dec [arcsec]
!	Table(n,38)	=	int(RA*100)   [RA_pm]
!	Table(n,39)	=	frac(RA*100)  [RA_pm]
!	Table(n,40)	=	int(Dec*100)  [Dec_pm]
!	Table(n,41)	=	frac(Dec*100) [Dec_pm]
!	Table(n,42:45)	=	fluxes in bands 1,...4 [dn]
!	Table(n,46:49)	=	uncert. in fluxes in bands 1,...4 [dn]
!	Table(n,50:53)	=	reduced chi squared in bands 1,...4
!	Table(n,54)	=	overall reduced chi squared
!	Table(n,57)	=	nIters (no. of iterations in fit)
!	Table(n,58)	=	nSteps (no. of steps on last iteration)
!	Table(n,59)	=	pmqual (coded for nblend0,swapped,radial
!                   position difference nonPM-MeanObsEpoch PM)
!   Table(n,60:63)  = fitradii for bands 1..4 [asec]
!   Table(n,64:65)  = P vector components
!
!	where n is the source number (n = 1,...nsrc)
!
!	Note: The RA uncertainty is expressed as a great-circle angle.
!
! Double precision versions of the estimated source positions are output in
! arrays RAlist(n), DEClist(n), where n = 1,...nsrc.  In doing so, the input
! versions of those arrays are overwritten.
!-----------------------------------------------------------------------------

        implicit real*4 (a-h,o-z)
	implicit integer (i-n)

        integer*4 nf0                ! JWF B30427
        integer*2 order(nf0)         ! JWF B30503
        character*200 basename(nf0)  ! JWF B30427

! Set some parameter values.
	include 'wpro.inc'
	integer, parameter :: maxblend = 100
	integer, parameter :: maxblend2 = 200
	integer, parameter :: maxbands = 4
	real, parameter :: dtor = 0.0174533
	real, parameter :: tle = 1.08574	! 2.5*log10(e)

! Declare the input arrays.
   	real(8) RAlist(nmax), DEClist(nmax),JD(nf,4)
   	real(4) Array(nnx,nny,nfpix),Unc(nnx,nny,nfpix)
	real(4) psfs(npsize,npsize,npnmax,4)
	real(4) psfuncs(npsize,npsize,npnmax,4)
	real(4) pxcent(4),pycent(4),ppix(4),fwhm(4)
	real(4) Xpos(nmax,nf,4), Ypos(nmax,nf,4), Rstann(4), Rstwid(4)
	real(4) LBACK(nmax,nf,4), LSIG(nmax,nf,4), Lconf(nmax,nf,4)
	real(4) MAGSTD(nmax,4), eMAGSTD(nmax,4), Rsat(nall,4), zero(4)
c	real(4) Table(nmax,32)
c       real(4) Table(nmax,59)  ! JWF B30221
        real(4) Table(nmax,65)  ! JWF B31209
	real(4) framefluxes(nmax,nf,4),frameuncs(nmax,nf,4),framechis(nmax,nf,4)
	integer FLGSTD(nmax,4), adb_nmax, IDlist(nmax)
	integer Mask(nnx,nny,nfpix)
        integer fatal, trouble, ntroublebits, troublebits(6), saturation, fbits(32)
	integer satbits(9),satpix(9),ignorebits(5),ignorepix, IDpri
        integer*4 pix_order(nf,4)
   	integer(2) wflag(nf,4), MapPSFs(nxfp,nyfp,4)
	integer(2) xTRANS(nf,4),yTRANS(nf,4), SatNum(nall,4)
	character(150000) HDR(nf,4), Hfits

! Declare the internal arrays.
   	real(8), allocatable :: RAlistnew(:), DEClistnew(:)
	real(8), allocatable :: x8set(:,:,:),y8set(:,:,:)
c	real(4), allocatable :: p(:), grad(:), pmost(:)   ! JWF B30207
	real(4), allocatable :: p(:), pmost(:)            ! JWF B30207
	real(4), allocatable :: p_PM(:), p_preADB(:)      ! JWF B30207
	real(4), allocatable :: Xposnew(:,:,:),Yposnew(:,:,:)
	real(4), allocatable :: LBACKnew(:,:,:),LSIGnew(:,:,:)
	real(4), allocatable :: Rsatnew(:,:)
c	integer(4), allocatable :: npassprev(:),offset(:,:,:)              ! JWF B31206
	integer(4), allocatable :: npassprev(:),offset(:,:,:),IDlistnew(:) ! JWF B31206
	logical(1), allocatable :: onframe(:,:,:)
	real(8) ra8, dec8, x8, y8, x8pri, y8pri, ra8off, dec8off
	real(8) RApri,Decpri,JD0
	real(8) racan(maxblend),deccan(maxblend),racannew,deccannew
	real(4) pxc(maxbands),pyc(maxbands)
	real(4) xpixscale(maxframes),ypixscale(maxframes),tframe(maxframes)
	real(4) fitradii0(maxbands), fitradii(maxbands), fitradii_std(maxbands)
	real(4) pscale(maxbands),psamp(maxbands),sigfac(maxbands)
c       maxpix is set to 200000
	real(4) dvalue(maxpix,maxbands),varpsf(maxpix,maxbands),sigset(maxpix)
	real(4) sigdvalue(maxpix,maxbands),sigdvalue0(maxpix,maxbands)
	real(4) flux(maxbands,maxblend),sigflux(maxbands,maxblend)
        real*8  flux_preADB(maxbands,maxblend)  !  JWF B30226
	real(4) flux_nonPM(maxbands,maxblend),sigflux_nonPM(maxbands,maxblend)! JWF B30208
	real(4) finalflux(maxbands,maxblend),flux0(maxbands,maxblend)
	real(4) frameflux(maxframes,maxbands,maxblend)
	real(4) framesig(maxframes,maxbands,maxblend)
	real(4) framechi(maxframes,maxbands)
	real(4) rchisqb(maxbands),ss(maxblend),backbuf(100000)
	real(4) sigx(maxblend),sigy(maxblend),sigxy(maxblend)
	real(4) apflux(maxbands,maxblend), sigapflux(maxbands)
	real(4) wtapflux(maxbands),p_final(maxblend2)
	real(4) flux_init(maxbands,maxblend),rchisqb_init(maxbands)
	real(4) fsat(maxbands),ftroubles(maxbands),freject(maxbands)
	real(4) troublerpix, troublersq
	integer WCS(maxframesin,maxbands), offscl
	integer ifq(maxframes),ibq(maxframes),kpn(maxframes)
	integer ivalue(maxpix,maxbands), jvalue(maxpix,maxbands)
	integer ifvalue(maxpix,maxbands),primask(maxpix,maxbands)
	integer neighbors(maxblend),npsf(maxbands),nnpsf(maxbands)
	integer icomp(maxblend),iband(maxbands),nvp(maxbands)
   	integer nxorig(*), nyorig(*), nx(*),ny(*),nnx,nny, nmax, nout, pcount(maxbands)
	integer nvalues(maxbands),npvalues(maxbands),npall(maxbands)
	integer isatsampmin(maxbands),satsamples(maxbands)
	integer(2) nlayers(maxbands)
	logical(4) debug,newlayer(maxbands),posconstrain, SPIT, MIPS
	logical(4) fullset,ok,goodcomp(maxblend),goodcomp_init(maxblend)
        logical*4 goodcomp_preADB(maxblend)
	logical(4) it_is_saturated, newpix

	common /funccom/ racan,deccan,nvalues,ivalue,jvalue,dvalue,
     *	    sigdvalue,ifvalue,WCS,ifq,ibq,primask,nframes,
     *	    kpn,nnpsf,pxc,pyc,psamp,nblend,goodcomp,icp,
     *	    xpixscale,ypixscale,flux,rchisq,rchisqb,tframe
	common /poscom/ posconstrain
c
        logical*4 ok2, DoThisOne(maxblend), NotRunaway           !  JWF B30208
        real*8 sumJD, sumWgtJD, WgtJD, Dtmp1, Dtmp2              !  JWF B30207
        real*8 sumLconf(4), sumvarpsf(4)                         !  JWF B30419
        real*4 covRpmR, covDpmD, covpmRD, dist, cd, pmra, pmdec  !  JWF B21221
        real*4 blendist, blendist2(4), medsig(4), Skew, hival    !  JWF B30419/B306046
        integer*4 iter, nsteps, nsteps2, iter2,                  !  JWF B30129
     +            m, m0, mf, kk, jj, nFusable, nblend_preADB,    !  JWF B30208
     +            icp_preADB, icomp_preADB(maxblend),            !  JWF B30319
     +            nblend_nonPM, kTarget, nUsed, nitmmm           !  JWF B30325/B30606
c
        include 'jwfcom.f'             ! JWF B21219
c
	data ignorebits /10, 4, 21, 27, 28/ ! If any of these bits is set, then
					   ! ignore that pixel for the purpose
					   ! of finding earliest saturation
	data satbits /10, 11, 12, 13, 14, 15, 16, 17, 18/
					   ! These bits indicate saturation at
					   ! various samples up the ramp
c
        Real*4       RBlank           ! JWF B21120
        Integer*4    IBlank           ! JWF B21120
        Data         IBlank/-1/       ! JWF B21120
        equivalence (Rblank, IBlank)  ! JWF B21120

        integer*4 kStat, nblend0, ntry                              ! JWF B30312/B30609
        logical*4 GotTarget, Swapped                                ! JWF B30312
        logical*4 GotAllTargets, GotTargetK(10)
        data      GotTarget, GotTargetK, GotAllTargets/12*.false./  ! JWF B30412
        data      hival/1.0e20/                                     ! JWF B30606
        real*4    dRA
        integer*4 nDoF, nDoF_pm, nWarnPMNaN, NdistWarn  ! JWF B30505
        real*8    vsum, racan0, deccan0                 ! JWF B30501
        data kTarget/0/, nWarnPMNaN/0/, NdistWarn/0/    ! JWF B30505
        real*4    dammmmean, dammmmed, dammmmode, dammmsigma, BGsig
c
	integer nmmmwarn, maxmmmwarn
	data nmmmwarn/0/, maxmmmwarn/50/
c
        character*200 FSimgNam, FSimgNam0                           ! JWF B30824
        character*20  NumStr                                        ! JWF B30824
        character*2   TrgStr                                        ! JWF B30824
        character*1   BndStr                                        ! JWF B30824
        integer*4 naxis, naxes(2), status, blocksize, outunit,      ! JWF B30824
     +            group, fpixel, nelements, kpix, bitpix, npix      ! JWF B30824
        logical*4 simple, extend, zexist                            ! JWF B30824
   	real*4,   allocatable :: pstamp1(:), pstamp2(:)             ! JWF B30824
       logical IzBad
c
c
c
c
c
        if (GotAllTargets .and. QuitEarly) return
c
c	JD0 = 55197.d0	! Modified Julian Date at 2010 Jan 1 0h UT
c       JD0 = 55400.d0	! Modified Julian Date at 2010 Jul 23 0h UT, JWF B21105
 	debug = .false.
c	debug = .true.

	print *, '$Id: wpro_v6.f 13438 2014-02-27 18:10:00Z tim $'
c	print *,'WPRO Version 6;  2012-Mar-22'
	print *,'wpro_v6 for WSDS v7;  2013-Dec-05'      !  JWF B31205
        print *,'maxsteps =        ',maxsteps            !  JWF B30219
        print *,'MeanObsEpochType =',MeanObsEpochType    !  JWF B30219

	istat = 0

! Set templates for interpreting pixel masks.

! count glitch pixels or persistent latent pixels detected in dynasky
	ntroublebits = 2
	troublebits(1) = 21
	troublebits(2) = 28
	trouble = 0
	do i=1,ntroublebits
	  trouble = trouble + 2**troublebits(i)
	end do
	troublerpix = 1.25 ! trouble bit search radius, frame pixels

	saturation = 0                  ! bits 10-18 => saturation
	do i = 10,18
	    saturation = saturation + 2**i
	enddo
	    

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

	satpix = 2**satbits
	ignorepix = sum(2**ignorebits)

! Set the fitting radii.
	if (minval(fwhm) > 0.) then
	    fitradii0 = 1.25*fwhm		! arcsec
	else
	    print *,'WPRO: FWHM values not set'
	    istat = 1
	    return
	endif

! Set passive deblending parameters.
	blendist = 2.*maxval(fwhm)	! critical separation for passive
     +               *RblendFac         ! deblending [arcsec]
        blendist2(1) = 2.0*RBlendFac*fwhm(1)
        blendist2(2) = 2.0*RBlendFac*fwhm(2)
        blendist2(3) = 2.0*RBlendFac*fwhm(3)
        blendist2(4) = 2.0*RBlendFac*fwhm(4)

	pratio	= 0.1			! critical flux ratio for passive
					! deblending
c	nbmax = 3			! maximum number of passively deblended
					! components (now in jwfcom)
	pdb_dchi = 0.5			! minimum acceptable reduction in the
					! reduced chi squared resulting from
					! adding a new passively-deblended
					! source component

! Set active deblending parameters.
!	adb_nmax			= maximum allowable number of
!					  components which can be added during 
!					  active deblending. It is set in the
!					  calling program (WPHot).  A zero
!					  value turns off active deblending.
	adb_chiThresh = 1.5		! threshold value of the reduced chi
					! squared for active deblending to be
					! attempted
	adb_dchi = 0.25	        	! minimum acceptable reduction in the
					! reduced chi squared due to active 
					! deblending
	adb_chimax = 3.			! maximum acceptable value of reduced
					! chi squared after deblending
	adb_pksep = 1.6*minval(fwhm)	! minimum separation for which detector
					! can distinguish sources [arcsec]
	adb_pksep = 2.2*minval(fwhm)	! minimum separation for which detector
					! can distinguish sources [arcsec]
	adb_minsep = 0.65*minval(fwhm)	! lower limit for separation of
					! deblended sources [arcsec]

        nsentinel = 0                   ! Count out of range pixels

	if (nf > maxframesin) then
	    print *,'WPRO: Need to increase maxframesin'
	    istat = 1
	    return
	endif

! Normalize the PSFs.
	do ib = 1,maxbands
	  if (any(wflag(1:nf,ib)==1)) then
	    do k = 1,npnmax
		psfsum = sum(psfs(1:npsf(ib),1:npsf(ib),k,ib))
		do j = 1,npsize
		do i = 1,npsize
		    psfs(i,j,k,ib) = psfs(i,j,k,ib)/psfsum
		    psfuncs(i,j,k,ib) = psfuncs(i,j,k,ib)/psfsum
		enddo
		enddo
	    enddo
	  endif
	enddo

! Set up various quantities to be passed in common to FUNC.
	nnpsf = npsf
	pxc = pxcent
	pyc = pycent
	psamp = ppix

c	allocate (RAlistnew(nmax))                 ! JWF B30312
c	allocate (DEClistnew(nmax))                ! JWF B30312
c	allocate (Xposnew(nmax,nf,4))              ! JWF B30312
c	allocate (Yposnew(nmax,nf,4))              ! JWF B30312
c	allocate (LBACKnew(nmax,nf,4))             ! JWF B30312
c	allocate (LSIGnew(nmax,nf,4))              ! JWF B30312
c	allocate (Rsatnew(nall,4))                 ! JWF B30312
c	allocate (npassprev(nmax))                 ! JWF B30312
c	allocate (x8set(maxblend,nf,maxbands))     ! JWF B30312
c	allocate (y8set(maxblend,nf,maxbands))     ! JWF B30312
c	allocate (offset(maxblend,nf,maxbands))    ! JWF B30312
c

c	print *,'=========== allocating with nmax,nall,nf,maxblend,maxbands = ',nmax,nall,nf,maxblend,maxbands !!! TPC

	allocate (RAlistnew(nmax),    DEClistnew(nmax),    Xposnew(nmax,nf,4),
     +            Yposnew(nmax,nf,4), LBACKnew(nmax,nf,4), LSIGnew(nmax,nf,4),
     +	          Rsatnew(nall,4),    npassprev(nmax),     x8set(maxblend,nf,maxbands),
     +            y8set(maxblend,nf,maxbands), offset(maxblend,nf,maxbands),
     +            IDlistnew(nmax), stat = kStat)           ! JWF B30312/B31205
        if (kStat .ne. 0) then
          print *,'ERROR (wpro_v6): cannot allocate arrays; kStat =',kStat
          istat = kStat
          return
        end if
c
	npassprev = 0
 	RAlistnew  = 0.0d0               ! JWF B30312
 	DEClistnew = 0.0d0               ! JWF B30312
 	Xposnew    = 0.0                 ! JWF B30312
 	Yposnew    = 0.0                 ! JWF B30312
 	LBACKnew   = 0.0                 ! JWF B30312
 	LSIGnew    = 0.0                 ! JWF B30312
 	Rsatnew    = 0.0                 ! JWF B30312
 	x8set      = 0.0d0               ! JWF B30312
 	y8set      = 0.0d0               ! JWF B30312
 	offset     = 0                   ! JWF B30312

! Initialize coordinate transformation routine.
	do ifs = 1,nf
	do ib = 1,maxbands
	  WCS(ifs,ib) =  -1
	  if (wflag(ifs,ib) == 1) then

	    Hfits = HDR(ifs,ib)
	    call wcsinit(Hfits,iwcs)
	    if ( iwcs < 0) then
		print *,'*** WCSinit ERROR',iwcs
		call exit(1)
	    endif
	    WCS(ifs,ib) = iwcs

	   endif

	enddo
	enddo

! Find out which candidates are on which frames.
c	print *,'    allocating with nmax,nf,maxbands = ',nmax,nf,maxbands !!! TPC
	allocate (onframe(nmax,nf,maxbands))
	do n = 1,nsrc
	    do ifs = 1,nf
	    do ib = 1,maxbands
	      itrans = Xpos(n,ifs,ib) - xTRANS(ifs,ib)
	      jtrans = Ypos(n,ifs,ib) - yTRANS(ifs,ib)
	      if(itrans > 0 .and. itrans <= nx(ib) .and.
     *	         jtrans > 0 .and. jtrans <= ny(ib) .and.
     *		 Xpos(n,ifs,ib) > 0. .and. Xpos(n,ifs,ib) <= nxorig(ib) .and.
     *		 Ypos(n,ifs,ib) > 0. .and. Ypos(n,ifs,ib) <= nyorig(ib)) then
		onframe(n,ifs,ib) = .true.
	      else
		onframe(n,ifs,ib) = .false.
	      endif
	    enddo
	    enddo
	enddo

! Get mean values of sigapflux for each band.
	sigapflux = 0.
	wtapflux = 0.
        do n = 1,nsrc
	    do ib = 1,4
		if (eMAGSTD(n,ib) < 9. .and. eMAGSTD(n,ib) > 0) then
		    apflx = 10.**(-0.4*(MAGSTD(n,ib) - zero(ib)))
		    sigapflux(ib) = sigapflux(ib) + eMAGSTD(n,ib)*apflx/tle
		    wtapflux(ib) = wtapflux(ib) + 1.
		endif
	    enddo
	enddo
	sigapflux = sigapflux/wtapflux
					
! Process each candidate in detection list.
	ihalf = nsrc / 2
        iq = nsrc / 4
        i3q = (nsrc * 3 ) / 4
	where(Lconf < 0.) Lconf = 0.
	framefluxes = -1.e35
	frameuncs = -1.e35
	framechis = -1.e35
	nsrc_adb  = 0
	m         = 0
        mf        = 0
	SatNum    = 0

        do n = 1,nsrc				! <<<<<<< Begin source loop
	  if (any(onframe(n,1:nf,1:maxbands))) then
             Swapped = .false.
c	     if (n.eq.ihalf) write (6,*) '50% complete'
c            if (n.eq.iq) write (6,*) '25% complete'
c            if (n.eq.i3q) write (6,*) '75% complete'
c            print *,'(317) Doing source #', n  !  JWF B21109 dbg

! Set fitting radii.
	    do ib = 1,4
		fitradii(ib) = min(2.*Rsat(n,ib), (npsf(ib)*ppix(ib)*3600.)/2.)
	    enddo
	    where(fitradii < fitradii0) fitradii = fitradii0
c	    sigfac = sqrt((fitradii**2 - Rsat(n,1:4)**2))/fitradii0     ! JWF B30930
       	    sigfac = sqrt(abs(fitradii**2 - Rsat(n,1:4)**2))/fitradii0  ! JWF B30930
	    where(sigfac < 1.) sigfac = 1.

! Establish focal-plane sampling intervals at primary candidate location.
! Also, find out which PSFs to use for this source.
	    cd = cos(DEClist(n)*dtor)

	    do ib = 1,maxbands
		pscale(ib) = 0.
		pcount(ib) = 0
	    enddo

	    iframe = 0
	    do ifs = 1,nf
	    do ib = 1,maxbands
	      if (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
		iframe = iframe + 1
		if (iframe > maxframes) then
		    print *,'WPRO: Need to increase maxframes'
		    istat = 1
		    return
		endif
		ifq(iframe) = ifs
		ibq(iframe) = ib
		x8 = Xpos(n,ifs,ib) * 1.d0
		y8 = Ypos(n,ifs,ib) * 1.d0
		iwcs = WCS(ifs,ib)
		call pix2wcs (iwcs, x8, y8, RA8, DEC8)
		call pix2wcs (iwcs, x8+1.d0, y8, RA8off, DEC8off)
                dRA = RA8off-RA8
                if (abs(dRA) .gt. 180.0) dRA = dRA - sign(360.0,dRA)
		xpixscale(iframe) =
     *		    sqrt((dRA*cd)**2 + (DEC8off-DEC8)**2)
		call pix2wcs (iwcs, x8, y8+1.d0, RA8off, DEC8off)
                dRA = RA8off-RA8
                if (abs(dRA) .gt. 180.0) dRA = dRA - sign(360.0,dRA)
		ypixscale(iframe) =
     *		    sqrt((dRA*cd)**2 + (DEC8off-DEC8)**2)
		pscale(ib) = pscale(ib) +
     *		    (xpixscale(iframe) + ypixscale(iframe))/2.
		pcount(ib) = pcount(ib) + 1

! psfset gets indexed by frame;  psfs are indexed by ID# and band.
! The ID# is obtained by looking up MapPSFs for the Xpos,Ypos.
		ifp = min(max(nint(x8), 1), nxfp)
		jfp = min(max(nint(y8), 1), nyfp)
		kpn(iframe) = MapPSFs(ifp,jfp,ib)
	      endif
	    enddo
	    enddo
	    nframes = iframe
c
	    where(pcount /= 0) pscale = 3600.*pscale/pcount

! Look for any immediate neighbors. ! JWF B30220: the primary candidate is neighbor #1
            if (BlendType .eq. 1) then
              call find_neighbors(n,RAlist,DEClist,MAGSTD,eMAGSTD,onframe,nmax,
     *		nf,maxbands,nsrc,blendist,pratio,nbmax,nblend,neighbors)
            else
              call find_neighbors_perband(n,RAlist,DEClist,MAGSTD,eMAGSTD,onframe,nmax,
     *		nf,maxbands,nsrc,blendist2,pratio,nbmax,nblend,neighbors)
            endif
	    if (any(Rsat(n,1:4) /= 0.)) nblend = 1
	    nblend0 = nblend
	    icp = 1		! component number for primary

            do 5 k = 1, 10
              GotTarget = .false.
              do 3 kk = 1, nblend
                GotTarget = GotTarget .or.
     +                    ((dabs(RAlist(neighbors(kk)) -TargetRA(k))  .lt. 1.0d-5)      ! JWF dbg
     +               .and. (dabs(Declist(neighbors(kk))-TargetDec(k)) .lt. 1.0d-5))     ! JWF dbg
                if (GotTarget) kTarget = k
3             continue
              if (GotTarget) then
                sumLconf  = 0.0d0
                sumvarpsf = 0.0d0
                if   ((dabs(RAlist(n)-TargetRA(k))   .lt. 1.0d-5)          ! JWF dbg
     +          .and. (dabs(Declist(n)-TargetDec(k)) .lt. 1.0d-5)) then    ! JWF dbg
                  print *,'Target',k,'identified as input source no.', n   ! JWF dbg
                  GotTargetK(k) = .true.
                else
                  print *,'Target',k,
     1              'identified in blend group of input source no.',n   ! JWF dbg
                  print *,'List numbers in blend:',neighbors(1:nblend)
                end if
                if (BlendType .eq. 1) then
                  print *,'blendist =',blendist
                else
                  print *,'blendist2 =',blendist2
                end if
                GotAllTargets = .true.
                do 4 kk = 1, 10
                  if ((TargetRA(kk) .le. 360.0d0)
     +                 .and. .not.GotTargetK(kk)) then
                    GotAllTargets = .false.
                    go to 6
                  end if
4               continue
                go to 6
              end if
5           continue
6           continue
            if (GotTarget) then
              print *,'wpro_v6(561): fitradii:', fitradii
              print *,'wpro_v6(562): sigfac:', sigfac
            end if

! Find blend components on all frames.
	    do nn = 1,nblend
	 	racan(nn) = RAlist(neighbors(nn)) * 1.d0
	  	deccan(nn) = DEClist(neighbors(nn)) * 1.d0
              do ifs = 1,nf
	  	do ib = 1,maxbands
	 	    if (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
			iwcs = WCS(ifs,ib)
			offscl = -1
			call wcs2pix(iwcs,racan(nn),deccan(nn),x8,y8,offscl)
			x8set(nn,ifs,ib) = x8
			y8set(nn,ifs,ib) = y8
			offset(nn,ifs,ib) = offscl
	 	    endif
	  	enddo
	  	enddo
	    enddo

! Gather data for this candidate.
	    nvalues = 0
	    npvalues = 0
	    npall = 0
	    nlayers = 0
	    fsat = 0.		! fraction of saturated pixels
	    freject = 0.	! fraction of pixels excluded from solution
	    primask = 0.
	    varpsf = 0.

	    do nn = 1,nblend
              if (GotTarget) then
                print *,'nn, neighbors(nn):',nn, neighbors(nn) ! JWF B30430
                print *,'MAGSTD(neighbors(nn),1:4), eMAGSTD(neighbors(nn),1:4):',
     +                   MAGSTD(neighbors(nn),1:4), eMAGSTD(neighbors(nn),1:4)
              end if
	      racan(nn) = RAlist(neighbors(nn)) * 1.d0
	      deccan(nn) = DEClist(neighbors(nn)) * 1.d0
	      do ib = 1,maxbands
		if (eMAGSTD(neighbors(nn),ib) > 0. .and.
     *		    eMAGSTD(neighbors(nn),ib) < 9. .and.
     *		    MAGSTD(neighbors(nn),ib) > -50. .and.
     *		    MAGSTD(neighbors(nn),ib) < 50.) then
		  apflux(ib,nn) = 10.**(-0.4*(MAGSTD(neighbors(nn),ib)
     *		    - zero(ib)))
		else
		  apflux(ib,nn) = 0.
		endif
		flux(ib,nn) = 0.
		sigflux(ib,nn) = 0.
	      enddo
            enddo                               ! JWF B30502

	    do nn = 1,nblend                    ! JWF B30502
              if (GotTarget) then
                print *,'apflux(1:4,nn):',apflux(1:4,nn)  ! JWF B30430
              end if
	      iframe = 0
	      if (nn==1) isatsampmin = 10
	      do ifs = 1,nf
	 	newlayer = .true.
	      do ib = 1,maxbands
	       if (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
                ifs_pix = pix_order(ifs,ib)
		if(ifs_pix .lt. 1 .or. ifs_pix .gt. nfpix) then
                  print *, '***ERROR: Pix_order sync error: ifs,ib,ifs_pix = ',ifs,ib,ifs_pix
		  print *,pix_order
                  call exit(99)
                endif
		iframe = iframe + 1
c		tframe(iframe) = JD(ifs,ib) - JD0   ! don't have JD0 yet
		iwcs = WCS(ifs,ib)
		offscl = -1
		call wcs2pix(iwcs,racan(nn),deccan(nn),x8,y8,offscl)
		ifp = min(max(nint(x8), 1), nxfp)
		jfp = min(max(nint(y8), 1), nyfp)

		k = MapPSFs(ifp,jfp,ib)
		xprat = xpixscale(iframe)/psamp(ib)
		yprat = ypixscale(iframe)/psamp(ib)
		rn = xprat*yprat

		i0 = nint(Xpos(neighbors(nn),ifs,ib))
		j0 = nint(Ypos(neighbors(nn),ifs,ib))
		fitpix = fitradii(ib)/pscale(ib)
	  	ifit = nint(fitpix + 0.5)
		ilo = max(i0-ifit,1)
		ihi = min(i0+ifit,nxorig(ib))
		jlo = max(j0-ifit,1)
		jhi = min(j0+ifit,nyorig(ib))
		rfitsq = fitpix**2
		back = LBACK(neighbors(nn),ifs,ib)
                if (GotTarget) print *,'wpro_v6(624): local background in band',ib,
     +           'for blend component',nn,'=',back ! JWF B30430
		do j = jlo,jhi
		do i = ilo,ihi
		    call checkpix(i,j,iframe,ib,nn,nvalues,ivalue,jvalue,
     *			ifvalue,maxpix,maxbands,newpix)
		    dx = i - Xpos(neighbors(nn),ifs,ib)
		    dy = j - Ypos(neighbors(nn),ifs,ib)
		    rsq = dx**2 + dy**2
		    itrans = i - xTRANS(ifs,ib)
		    jtrans = j - yTRANS(ifs,ib)

		    if(rsq <= rfitsq .and. itrans >= 1 .and. itrans <= nx(ib)
     *			.and. jtrans >= 1 .and. jtrans <= ny(ib) .and. newpix)
     *			 then

		      if(Array(itrans,jtrans,ifs_pix) .lt. -1e20) then
! Sentinel pixel encountered
		        nsentinel = nsentinel + 1
! Mask out from here on utilizing the outlier rejection bit
			Mask(itrans,jtrans,ifs_pix) = ior(Mask(itrans,jtrans,ifs_pix),2**27)
		        if(nsentinel .le. 10) then
		          print *,'===WARNING !!! Sentinel pixel encountered. band=',ib,
     1			   ', frame#=',ifs,', ifs_pix=',ifs_pix,', i,j=',i,j,
     2			   ', itrans,jtrans=',itrans,jtrans
		        endif
		      endif

		      npall(ib) = npall(ib) + 1	! total pixels in fitting region
		      it_is_saturated = .false.
		      if(iand(Mask(itrans,jtrans,ifs_pix),saturation)/=0)
     *			    it_is_saturated = .true.
		      if (it_is_saturated) fsat(ib) = fsat(ib) + 1.

! If the pixel is saturated due to large source flux, then find out the first
! sample up the ramp in which saturation occurred.
		      if(nn==1) then

			if(iand(Mask(itrans,jtrans,ifs_pix), ignorepix)==0
     *			    .and. it_is_saturated) then
			    isatsamp = 10
			    isamp = 0
			    do while(isamp < 9 .and. isatsamp==10)
				isamp = isamp + 1
				if(iand(Mask(itrans,jtrans,ifs_pix),
     *				    satpix(isamp)) /= 0) isatsamp = isamp
			    enddo
			    if (isatsamp<isatsampmin(ib))
     *				isatsampmin(ib) = isatsamp
			endif
		      endif ! (nn==1)

		      if(iand(Mask(itrans,jtrans,ifs_pix), fatal)==0 .and.
     *			    back > -900. .and. nvalues(ib)<maxpix) then
			if(nn==1) then
			    npvalues(ib) = npvalues(ib) + 1
			    if (newlayer(ib)) then
				nlayers(ib) = nlayers(ib) + 1
				newlayer(ib) = .false.
                                if (GotTarget) print *,'band',ib,'uses ', ! JWF B30429
     +                             basename(order(ifs))(1:9)              ! JWF B30429
			    endif
			endif
			nvalues(ib) = nvalues(ib) + 1
			if (nvalues(ib)==maxpix)
     *			    print *,'WARNING!!!  WPRO: Need to increase maxpix'
			ivalue(nvalues(ib),ib) = i
			jvalue(nvalues(ib),ib) = j
			dvalue(nvalues(ib),ib) = Array(itrans,jtrans,ifs_pix) - back
			sigdvalue0(nvalues(ib),ib) = Unc(itrans,jtrans,ifs_pix)
c                       if (GotTarget) then                                       ! JWF B30426
c                         print *,'(itrans,jtrans,ifs,ib):', itrans,jtrans,ifs,ib ! JWF B30426
c                         print *,'(Array,Unc,Mask,back):',                       ! JWF B30426
c    +                              Array(itrans,jtrans,ifs_pix),                 ! JWF B30426
c    +                                Unc(itrans,jtrans,ifs_pix),                 ! JWF B30426
c    +                               Mask(itrans,jtrans,ifs_pix), back            ! JWF B30426
c                       end if                                                    ! JWF B30426
			vsum = 0.0d0
			do nnn = 1,nblend
			  ipp = nint(pxc(ib) + (i - x8set(nnn,ifs,ib))*xprat)
			  jpp = nint(pyc(ib) + (j - y8set(nnn,ifs,ib))*yprat)
			  if (ipp>0 .and. ipp<=npsf(ib) .and. jpp>0 .and.
     *			   jpp<=npsf(ib) .and. offset(nnn,ifs,ib)>=0) vsum =
     *			    vsum + (rn*psfuncs(ipp,jpp,k,ib)*apflux(ib,nnn))**2
			enddo
			varpsf(nvalues(ib),ib) = vsum +
     *			    Lconf(neighbors(nn),ifs,ib)**2
			ifvalue(nvalues(ib),ib) = iframe
			if (nn == 1) primask(nvalues(ib),ib) = 1
		      else
			freject(ib) = freject(ib) + 1.
		      endif ! (iand(Mask(itrans,jtrans,ifs_pix), fatal)==0 .and.
        		    !	  back > -900. .and. nvalues(ib)<maxpix)
		    endif ! (rsq <= rfitsq .and. itrans >= 1 .and. itrans <= nx(ib)
     		          !  .and. jtrans >= 1 .and. jtrans <= ny(ib) .and. newpix)
		enddo ! i = ilo,ihi
		enddo ! j = jlo,jhi
	       endif
	      enddo ! ib = 1,maxbands
	      enddo ! do ifs = 1,nf

	    enddo ! nn = 1,nblend

	    where(npall /= 0) fsat = fsat/npall
	    where(npall /= 0) freject = freject/npall
	    satsamples = 0

	    do ib = 1,maxbands
		if (isatsampmin(ib) < 10) satsamples(ib) = isatsampmin(ib)

! Set Poisson noise to be constant over fitting region, and equal to median
! value.
		do nv = 1,nvalues(ib)
		    sigset(nv) = sigdvalue0(nv,ib)
		enddo
		call medsort(sigset,nvalues(ib),sigmed,q1,q2)
		do nv = 1,nvalues(ib)
		    sigdvalue(nv,ib) = sqrt(sigmed**2 + varpsf(nv,ib))
		enddo
	    enddo

! Calculate maximum likelihood nonPM solution.
	    nparms = 2*nblend
	    allocate (p(nparms))
	    allocate (pmost(nparms))
	    goodcomp = .false.
	    do nn = 1,nblend
		goodcomp(nn) = .true.
	    enddo
	    p = 0.                          ! position OFFSETS start out zero

! Do the minimization.
c           ftol = 1.e-3     ! JWF B30312
	    ftol = ftol_npm  ! JWF B30312
	    fret = xpixscale(1)/100.
c	    posconstrain = .true.         ! JWF B30308
	    posconstrain = poscon0        ! JWF B30308; namelist control, default T
	    flux = 0.
            if (GotTarget) then ! JWF dbg
             print *,'about to call most_wise for source #',n,' ...'  ! JWF dbg
             print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
c            print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg  [not set yet]
             print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
             print *,'p:',p(1:2*nblend) ! JWF dbg
             print *,'racan:', racan(1:nblend) ! JWF dbg
             print *,'deccan:',deccan(1:nblend) ! JWF dbg
             sigset = varpsf(1:nvalues(1),1)
             call tjsort(nvalues(1),sigset)
             print *,'varpsf(1:',nvalues(1),'1):', sigset(1:nvalues(1))
             sigset = varpsf(1:nvalues(2),2)
             call tjsort(nvalues(2),sigset)
             print *,'varpsf(1:',nvalues(2),'2):', sigset(1:nvalues(2))
             sigset = varpsf(1:nvalues(3),3)
             call tjsort(nvalues(3),sigset)
             print *,'varpsf(1:',nvalues(3),'3):', sigset(1:nvalues(3))
             sigset = varpsf(1:nvalues(4),4)
             call tjsort(nvalues(4),sigset)
             print *,'varpsf(1:',nvalues(4),'4):', sigset(1:nvalues(4))
            end if ! JWF dbg
	    call most_wise(p,nparms,psfs,npsize,npnmax,ftol,iter,fret,nsteps,.false.) ! JWF B30206
	    fret = func(p,psfs,npsize,npnmax,nDoF,GotTarget)
	    pmost = p
            if (GotTarget) then ! JWF dbg
             print *,'back from most_wise'  !  JWF dbg
             print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
c            print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg  [not set yet]
             print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
             print *,'p:',p(1:2*nblend) ! JWF dbg
             print *,'racan:', racan(1:nblend) ! JWF dbg
             print *,'deccan:',deccan(1:nblend) ! JWF dbg
            end if

! Check for runaway fit.   [JWF: this is where the RA:Dec solutions are checked]
	    fitmaxsq = (maxval(PMfac*fitradii)/3600.)**2
            NotRunaway = .true.                               ! JWF B30302
	    do nn = 1,nblend
	      if (p(nn)**2+p(nblend+nn)**2 >= fitmaxsq) then  ! JWF B21207
                if (WarnRunaway) then                                ! JWF B30329
                 print *,'Post most_wise: discarding blend component',
     +                 nn,'distance:',sqrt(p(nn)**2+p(nblend+nn)**2) ! JWF B21207
                 print *,'max distance allowed:',sqrt(fitmaxsq)      ! JWF B30329
                 print *,'mdet source no.:',n                        ! JWF B30328
                 print *,'mdet RA, Dec:',racan(nn),deccan(nn)        ! JWF B30329
                end if
                if (nn==1) then                                     ! JWF B30302
                  nBadBlend = nBadBlend + 1                         ! JWF B21207
                  NotRunaway = .false.                              ! JWF B30302
                end if                                              ! JWF B30302
		do ib = 1,maxbands
		  flux(ib,nn) = 0.
	        enddo
	        goodcomp(nn) = .false.
	      endif
	      icomp(nn) = nn
	      ss(nn) = sum((flux(1:maxbands,nn)/sigapflux(1:maxbands))**2)
	    enddo

	    nbands = 0
	    do ib = 1,maxbands
		if(nvalues(ib) > 8) then
		    nbands = nbands + 1
		    iband(nbands) = ib
		endif
	    enddo

            if (GotTarget) then                                   ! JWF dbg
              print *,'nvalues:', nvalues                         ! JWF dbg
              print *,'nlayers:', nlayers                         ! JWF dbg
              print *,'nbands, iband:', nbands, iband(1:nbands)   ! JWF dbg
            end if                                                ! JWF dbg

            if (TossZeroFlux) then
              if (.not.(any(flux(iband(1:nbands),1) /= 0.))) then    ! JWF B30207: Ken
                deallocate(p)                                        ! says toss entire
                deallocate(pmost)                                    ! solution if the
                if (NotRunaway) nAllZero(1) = nAllZero(1) + 1        ! primary component
                if (WarnZeroFlux) then                               ! has all zero fluxes
                  print *,'Post most_wise: '
     +                  //'Discarding all-zero primary component for '
     +                  //'source',n,'; RA:Dec =',RAlist(n),Declist(n)
                  print *,'nbands:', nbands
                  print *,'nvalues:', nvalues
                  print *,'flux(iband(1:nbands),icomp(1)):', flux(iband(1:nbands),icomp(1))
                end if
                go to 2000
              end if
            end if
c
	    if (nblend > 1) then
! Check for redundant components.  Start by ordering the components in
! decreasing strength.
		nswap = 0
		do nn = 1,nblend-1
		    do nnp = nn+1,nblend
			if(ss(nnp) > ss(nn)) then
			    icompnn = icomp(nn)
			    ssnn = ss(nn)

			    icomp(nn) = icomp(nnp)
			    ss(nn) = ss(nnp)

			    icomp(nnp) = icompnn
			    ss(nnp) = ssnn
			    nswap = nswap + 1
			endif
		    enddo
		enddo

                if (GotTarget) then ! JWF dbg
                 print *,'about to do hypothesis testing for source #',n,' ...'  ! JWF dbg
                 print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
                 print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
                 print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
                 print *,'p:',p(1:2*nblend) ! JWF dbg
                 print *,'racan:', racan(1:nblend) ! JWF dbg
                 print *,'deccan:',deccan(1:nblend) ! JWF dbg
                 print *,'nvalues:', nvalues           ! JWF B30509
                end if ! JWF dbg

! Now do some hypothesis testing.
                icp = icomp(1)       ! JWF B30110 [copied here from below at Ken's suggestion]
		flux_init =flux
		rchisqb_init = rchisqb
		rchisq_init = rchisq
		goodcomp_init = goodcomp
		fullset = .true.
		npassive = 0
		do while(fullset .and. npassive < nblend-1)
		    npassive = npassive + 1
		    goodcomp = .false.
		    do nn = 1,npassive
			goodcomp(icomp(nn)) = .true.
		    enddo
                    if (GotTarget) then
                      print *,'about to call ppack; pmost =',pmost
                      print *,'nblend, goodcomp:',nblend,goodcomp(1:nblend)
                      print *,'p:',p
                    end if
		    call ppack(pmost,nblend,goodcomp,p)	! p now has good
							! components only
                    if (GotTarget) then
                      print *,'back from ppack; pmost =',pmost
                      print *,'nblend, goodcomp:',nblend,goodcomp(1:nblend)
                      print *,'p:',p
                    end if
		    fret = func(p,psfs,npsize,npnmax,nDoF,.false.)
		    if(rchisq <= rchisq_init + pdb_dchi) fullset=.false.
		enddo
		if (fullset) then
		    flux = flux_init
		    rchisqb = rchisqb_init
		    rchisq = rchisq_init
		    npassive = nblend
		    goodcomp = goodcomp_init
		endif
	    else
		npassive = 1
	    endif

            if (GotTarget) then
              print *,'finished hypothesis testing for source #',n,' ...'  ! JWF dbg
              print *,'about to call ppack again; pmost =',pmost
              print *,'nblend, goodcomp:',nblend,goodcomp(1:nblend)
              print *,'p:',p
              print *,'nvalues:', nvalues           ! JWF B30509
            end if

  	    call ppack(pmost,nblend,goodcomp,p)

            if (GotTarget) then ! JWF dbg
             print *,'back from ppack for source #',n,' ...'  ! JWF dbg
             print *,'pmost:',pmost
             print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
             print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
             print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
             print *,'p:',p(1:2*nblend) ! JWF dbg
             print *,'racan:', racan(1:nblend) ! JWF dbg
             print *,'deccan:',deccan(1:nblend) ! JWF dbg
             print *,'nvalues:', nvalues           ! JWF B30509
            end if ! JWF dbg
c
            if (allocated(p_preADB)) then                          ! JWF B30329: this is
             if (WarnAlloc8err) then                               ! now known not to
              print *,'WARNING: p_preADB already allocated!'       ! be a problem, but
              print *,'         n:',n                              ! maybe we want to know
              print *,'         nblend_preADB:',nblend_preADB      ! JWF B30402
              if (nblend_preADB .gt. 0) then                       ! JWF B30402
               print *,'       p_preADB:',p_preADB(1:2*nblend_preADB)     ! JWF B30402
               print *,'    flux_preADB:',flux_preADB(1:4,1:nblend_preADB)! JWF B30402
               print *,'goodcomp_preADB:',goodcomp_preADB(1:nblend_preADB)! JWF B30402
               print *,'   icomp_preADB:',icomp_preADB(1:nblend_preADB)   ! JWF B30402
               print *,'     icp_preADB:',icp_preADB               ! JWF B30402
              end if                                               ! JWF B30402
              print *,'       de-allocating for re-allocation'     ! JWF B30329
              print *,'       non-PM values for current source:'   ! JWF B30402
              print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'icomp:',(icomp(nn),nn = 1, nblend)          ! JWF dbg
              print *,'flux:',(flux(1:4,nn), nn = 1, nblend)       ! JWF dbg
              print *,'p:',p(1:2*nblend)                           ! JWF dbg
              print *,'racan:', racan(1:nblend)                    ! JWF dbg
              print *,'deccan:',deccan(1:nblend)                   ! JWF dbg
             end if
             deallocate(p_preADB)                                  ! JWF B30329
            end if
            allocate (p_preADB(2*nblend)) !  JWF B30205
            p_preADB = p                  !  JWF B30209 - save for possible use
            nblend_preADB = nblend        !               as init for PM solution
            goodcomp_preADB = goodcomp    !  ditto
            flux_preADB = flux            !  also this, &c.
            icp_preADB = icp              !  JWF B30319
            icomp_preADB = icomp          !  JWF B30319
c
! Do active deblending if necessary.
	    nactive = 0
	    if (adb_nmax /= 0 .and. maxval(Rsat(n,1:4)) < adb_pksep) then

! Start by recalculating rchisq using a larger fitting radius.
	     fitradii_std = fitradii
	     fitradii = 2.2*fwhm	! arcsec
	     where(fitradii < fitradii_std) fitradii = fitradii_std
	     flux0 = flux
	     primask = 0.
	     varpsf = 0.
	     nvalues = 0
	     do nn = 1,nblend0
	      iframe = 0
	      do ifs = 1,nf
	      do ib = 1,maxbands
	       if (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
                ifs_pix = pix_order(ifs,ib)
		iframe = iframe + 1
c		tframe(iframe) = JD(ifs,ib) - JD0    ! don't have JD0 yet
		iwcs = WCS(ifs,ib)
		offscl = -1
		call wcs2pix(iwcs,racan(nn),deccan(nn),x8,y8,offscl)
		ifp = min(max(nint(x8), 1), nxfp)
		jfp = min(max(nint(y8), 1), nyfp)
		k = MapPSFs(ifp,jfp,ib)
		xprat = xpixscale(iframe)/psamp(ib)
		yprat = ypixscale(iframe)/psamp(ib)
		rn = xprat*yprat
		i0 = nint(Xpos(neighbors(nn),ifs,ib))
		j0 = nint(Ypos(neighbors(nn),ifs,ib))
		fitpix = fitradii(ib)/pscale(ib)
	  	ifit = nint(fitpix + 0.5)
		ilo = max(i0-ifit,1)
		ihi = min(i0+ifit,nxorig(ib))
		jlo = max(j0-ifit,1)
		jhi = min(j0+ifit,nyorig(ib))
		rfitsq = fitpix**2
		back = LBACK(neighbors(nn),ifs,ib)
                if (GotTarget) print *,'wpro_v6(1025): local background in band',
     +              ib,'=',back ! JWF B30430
		do j = jlo,jhi
		do i = ilo,ihi
		    call checkpix(i,j,iframe,ib,nn,nvalues,ivalue,jvalue,
     *			ifvalue,maxpix,maxbands,newpix)
		    rsq = (i - Xpos(neighbors(nn),ifs,ib))**2 +
     *			(j - Ypos(neighbors(nn),ifs,ib))**2
		    itrans = i - xTRANS(ifs,ib)
		    jtrans = j - yTRANS(ifs,ib)
		    if(rsq <= rfitsq .and. itrans >= 1 .and. itrans <= nx(ib)
     *			.and. jtrans >= 1 .and. jtrans <= ny(ib) .and. newpix)
     *			then
		      if(iand(Mask(itrans,jtrans,ifs_pix), fatal)==0 .and.
     *			    back > -900. .and. nvalues(ib)<maxpix) then
			nvalues(ib) = nvalues(ib) + 1
			ivalue(nvalues(ib),ib) = i
			jvalue(nvalues(ib),ib) = j
			dvalue(nvalues(ib),ib) = Array(itrans,jtrans,ifs_pix) - back
			sigdvalue0(nvalues(ib),ib) = Unc(itrans,jtrans,ifs_pix)
			vsum = 0.0d0
			do nnn = 1,nblend0
			 if (goodcomp(nnn)) then
			  ipp = nint(pxc(ib) + (i - x8set(nnn,ifs,ib))*xprat)
			  jpp = nint(pyc(ib) + (j - y8set(nnn,ifs,ib))*yprat)
			  if (ipp>0 .and. ipp<=npsf(ib) .and. jpp>0 .and.
     *			   jpp<=npsf(ib) .and. offset(nnn,ifs,ib)>=0) vsum =
     *			    vsum + (rn*psfuncs(ipp,jpp,k,ib)*flux(ib,nnn))**2
			  endif
			enddo
			varpsf(nvalues(ib),ib) = vsum +
     *			    Lconf(neighbors(nn),ifs,ib)**2
			ifvalue(nvalues(ib),ib) = iframe
			if (nn == 1) primask(nvalues(ib),ib) = 1
		      endif
		    endif
		enddo
		enddo
	       endif
	      enddo
	      enddo
	     enddo
	     do ib = 1,maxbands
		do nv = 1,nvalues(ib)
		    sigset(nv) = sigdvalue0(nv,ib)
		enddo
		call medsort(sigset,nvalues(ib),sigmed,q1,q2)
		do nv = 1,nvalues(ib)
		    sigdvalue(nv,ib) = sqrt(sigmed**2 + varpsf(nv,ib))
		enddo
	     enddo
	     fret = func(p,psfs,npsize,npnmax,nDoF,.false.)
             if (GotTarget) then ! JWF dbg
               print *,'after recalculating rchisq using a larger fitting radius:'
               print *,'nvalues:', nvalues           ! JWF B30509
               sigset = dvalue(1:nvalues(1),1)
               call tjsort(nvalues(1),sigset)
               print *,'dvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
               sigset = dvalue(1:nvalues(2),2)
               call tjsort(nvalues(2),sigset)
               print *,'dvalue(1:',nvalues(2),'2):', sigset(1:nvalues(2))
               sigset = dvalue(1:nvalues(3),3)
               call tjsort(nvalues(3),sigset)
               print *,'dvalue(1:',nvalues(3),'3):', sigset(1:nvalues(3))
               sigset = dvalue(1:nvalues(4),4)
               call tjsort(nvalues(4),sigset)
               print *,'dvalue(1:',nvalues(4),'4):', sigset(1:nvalues(4))
               sigset = sigdvalue(1:nvalues(1),1)
               call tjsort(nvalues(1),sigset)
               print *,'sigdvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
               sigset = sigdvalue(1:nvalues(2),2)
               call tjsort(nvalues(2),sigset)
               print *,'sigdvalue(1:',nvalues(2),'2):', sigset(1:nvalues(2))
               sigset = sigdvalue(1:nvalues(3),3)
               call tjsort(nvalues(3),sigset)
               print *,'sigdvalue(1:',nvalues(3),'3):', sigset(1:nvalues(3))
               sigset = sigdvalue(1:nvalues(4),4)
               call tjsort(nvalues(4),sigset)
               print *,'sigdvalue(1:',nvalues(4),'4):', sigset(1:nvalues(4))
               sigset = varpsf(1:nvalues(1),1)
               call tjsort(nvalues(1),sigset)
               print *,'varpsf(1:',nvalues(1),'1):', sigset(1:nvalues(1))
               sigset = varpsf(1:nvalues(2),2)
               call tjsort(nvalues(2),sigset)
               print *,'varpsf(1:',nvalues(2),'2):', sigset(1:nvalues(2))
               sigset = varpsf(1:nvalues(3),3)
               call tjsort(nvalues(3),sigset)
               print *,'varpsf(1:',nvalues(3),'3):', sigset(1:nvalues(3))
               sigset = varpsf(1:nvalues(4),4)
               call tjsort(nvalues(4),sigset)
               print *,'varpsf(1:',nvalues(4),'4):', sigset(1:nvalues(4))
             end if
	     if (rchisq > adb_chiThresh) then
              if (GotTarget) then ! JWF dbg
               print *,'(813) about to call act_wise.....'   ! JWF dbg
               print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
               print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
               print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
               print *,'sigflux:',(sigflux(1:4,nn), nn = 1, nblend) ! JWF dbg
               print *,'racan:', racan(1:nblend) ! JWF dbg
               print *,'deccan:',deccan(1:nblend) ! JWF dbg
               print *,'p:', p(1:2*nblend)
               print *,'nvalues:', nvalues           ! JWF B30509
              end if
		call act_wise(adb_nmax, adb_dchi, adb_chimax, adb_pksep,
     *		    adb_minsep, apflux, psfs,npsize,npnmax, varpsf,nlayers,
     *		    n,onframe,nmax,nf,maxbands,wflag,nf,p, p_final,nactive,
     *              iter2,nsteps2,nDoF,.false.) ! JWF B30117
              if (GotTarget) then ! JWF dbg
               print *,'(827) back from act_wise; nactive =',nactive   ! JWF dbg
               print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
               print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
               print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
               print *,'sigflux:',(sigflux(1:4,nn), nn = 1, nblend) ! JWF dbg
               print *,'racan:', racan(1:nblend) ! JWF dbg
               print *,'deccan:',deccan(1:nblend) ! JWF dbg
               print *,'p:', p(1:2*nblend)
               print *,'p_final:', p_final(1:2*nblend)
               print *,'nvalues:', nvalues           ! JWF B30509
              end if
		if (nactive /= 0) then
                    nsteps = nsteps2
                    iter   = iter2
		    deallocate(p)
		    mblend = ntrue(goodcomp,nblend)
		    allocate (p(2*mblend+2))
		    p = p_final(1:2*mblend+2)
		else
		    fitradii = fitradii_std
		    flux = flux0
		endif
	     else
		fitradii = fitradii_std
		flux = flux0
	     endif
	    endif
	    nsrc_adb = nsrc_adb + nactive    ! [JWF: end of active deblend loop]


!---------------------------------------------------- ADDED 12/15/09
! Re-evaluate pixel uncertainties in light of the new flux and position
! estimates.

	    mblend = ntrue(goodcomp,nblend)
	    nc = 0
	    ncp = 0
	    ftroubles = 0.	! fraction of pixels affected by tempcal masking

	    do nn = 1,nblend
	      if (goodcomp(nn)) then
	        nc = nc + 1
		if (nn==icomp(1)) ncp = nc
	 	racannew = racan(nn) + p(nc)/cos(deccan(nn)*dtor)
	  	deccannew = deccan(nn) + p(mblend+nc)
	  	do ifs = 1,nf
	  	do ib = 1,maxbands
	 	    if (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
			iwcs = WCS(ifs,ib)
			offscl = -1
			call wcs2pix(iwcs,racannew,deccannew,x8,y8,offscl)
			x8set(nn,ifs,ib) = x8
			y8set(nn,ifs,ib) = y8
			offset(nn,ifs,ib) = offscl
	 	    endif
	  	enddo ! ib = 1,maxbands
	  	enddo ! ifs = 1,nf
	      endif
	    enddo ! nn = 1,nblend
	    icp = icomp(1)
            ok2 = .true.
	    if (ncp /= 0) then
		RApri = racan(icp) + p(ncp)/cos(deccan(icp)*dtor)
		Decpri = deccan(icp) + p(mblend+ncp)
c           else                             !  JWF B21130
c             ok2 = .false.                  !  JWF B21207
c             print *,'WARNING (wpro_v6[808]): ncp never set, RApri & Decpri undefined!' ! JWF B21130
	    endif
	    primask = 0
	    varpsf = 0.
	    nvalues = 0
	    nvp = 0
	    do nn = 1,nblend0
	      iframe = 0
	      do ifs = 1,nf
	      do ib = 1,maxbands
c		print *,'----------- src ',n,',   comp ',nn, ',   band ',ib,', frame ', ifr, ' ... onframe = ',onframe(n,ifs,ib) !!! TPC
	       if (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
                ifs_pix = pix_order(ifs,ib)
		iframe = iframe + 1
c		tframe(iframe) = JD(ifs,ib) - JD0    ! don't have JD0 yet
		iwcs = WCS(ifs,ib)
		offscl = -1
		call wcs2pix(iwcs,racan(nn),deccan(nn),x8,y8,offscl)
		if (ncp /= 0) then
		    iprioffscl = -1
		    call wcs2pix(iwcs,RApri,Decpri,x8pri,y8pri,iprioffscl)
		endif
		ifp = min(max(nint(x8), 1), nxfp)
		jfp = min(max(nint(y8), 1), nyfp)
		k = MapPSFs(ifp,jfp,ib)
		xprat = xpixscale(iframe)/psamp(ib)
		yprat = ypixscale(iframe)/psamp(ib)
		rn = xprat*yprat
		i0 = nint(Xpos(neighbors(nn),ifs,ib))
		j0 = nint(Ypos(neighbors(nn),ifs,ib))

! Recalculate local background.
		rinsq = Rstann(ib)**2
		routsq = (Rstann(ib) + Rstwid(ib))**2
		irin = nint(Rstann(ib))
		irout = nint(Rstann(ib) + Rstwid(ib))
		ilo = max(i0-irout,1)
		ihi = min(i0+irout,nxorig(ib))
		jlo = max(j0-irout,1)
		jhi = min(j0+irout,nyorig(ib))
                nblank = 0
		kb = 0

                do j = jlo,jhi
		do i = ilo,ihi
		    rsq = (i - Xpos(neighbors(nn),ifs,ib))**2 +
     *			(j - Ypos(neighbors(nn),ifs,ib))**2
		    itrans = i - xTRANS(ifs,ib)
		    jtrans = j - yTRANS(ifs,ib)
		    if (rsq >= rinsq .and. rsq < routsq         ! JWF B21109
     *		    	.and. itrans >= 1 .and. itrans <= nx(ib)
     *			.and. jtrans >= 1 .and. jtrans <= ny(ib)) then
     		      if (iand(Mask(itrans,jtrans,ifs_pix), fatal)==0) then
			kb = kb+1
			backbuf(kb) = Array(itrans,jtrans,ifs_pix)
			if(backbuf(kb) .eq. -9991) nblank = nblank + 1 ! Debug aid TPC
                      end if
		    endif
		enddo ! i = ilo,ihi
		enddo ! j = jlo,jhi

		if (kb /= 0) then
c
c================= code seriously modified by TPC & JWF on B30606 ==============

                  call sort(kb,backbuf)

                  jj = BGtype(ib)
		  if (kb .lt. 2) jj = 1 ! punt
		  ntry = 0

 7		  continue

		  ntry = ntry + 1
                  !!! Debug TPC
                  if(ntry .gt. 1 .and. dammmsigma .le. -98 .and. nmmmwarn .lt. maxmmmwarn) then
                    print *, '=== Warning: Identical values in frame data.'
                    print *, '     Src,blendsrc,comp,band=',n,neighbors(nn),nn,ib
                    print *, '     kb,nblank=',kb,nblank
                    print *, '     x,y=',Xpos(neighbors(nn),ifs,ib),Ypos(neighbors(nn),ifs,ib)
                    print *, '     ifr,frame=',ifs,basename(order(ifs))(1:9)
                    print *, '     i/jlo/hi=',ilo,ihi,jlo,jhi
                    print *, '     irin,irout=',irin,irout
		    !print *, '     pixels=',backbuf(1:kb)
		    nmmmwarn = nmmmwarn + 1
		  endif


c                  print *, '    kb,backbuf(1),backbuf(kb),jj,ntry = ',kb,backbuf(1),backbuf(kb),jj,ntry  !!! TPC

                  Select Case (jj)
c
                  Case(1)               ! original mode estimator; the
                                        ! negative kb below signals medsort
                                        ! that it doesn't have to sort
                    call medsort(backbuf,-kb,bmed,q1,q2)
		    bsig = bmed - q1
                    back = estmode(backbuf,kb,bmed,bsig)
                    nUsed = kb
                    Skew  = -99.999
                    BGsig = bsig/sqrt(max(float(Nused-1),1.0))
c
                  Case(2)                           ! mmm mode
                    call mmm(backbuf,kb,hival, dammmmean, dammmmed, back,
     +                       dammmsigma, Skew, nUsed, nitmmm)
                    if (dammmsigma .lt. 0.0) then
                      jj = 1
                      go to 7
                    end if
                    BGsig = dammmsigma/sqrt(max(float(Nused-1),1.0))
c
                  Case(3)                           ! mmm median
                    call mmm(backbuf,kb,hival, dammmmean, back, dammmmode,
     +                       dammmsigma, Skew, nUsed, nitmmm)
                    if (dammmsigma .lt. 0.0) then
                      jj = 1
                      go to 7
                    end if
                    BGsig = dammmsigma/sqrt(max(float(Nused-1),1.0))
c
                  Case(4)                           ! mmm mean
                    call mmm(backbuf,kb,hival, back, dammmmed, dammmmode,
     +                       dammmsigma, Skew, nUsed, nitmmm)
                    if (dammmsigma .lt. 0.0) then
                      jj = 1
                      go to 7
                    end if
                    BGsig = dammmsigma/sqrt(max(float(Nused-1),1.0))
c
                  Case(5)                           ! recursive median distance
                                                    ! rejection & trimmed mean
                    call rmdr(backbuf,kb,BGtrimFrac,SkThresh,back,nUsed,Skew,bsig)
                    BGsig = dammmsigma/sqrt(max(float(Nused-1),1.0))
c
                  End Select
c
c		    if(ib .eq. 4) then
c		      write(112,*) neighbors(nn),nn,ifs,basename(order(ifs))(1:9),ib,
c     1                             Xpos(neighbors(nn),ifs,ib),Ypos(neighbors(nn),ifs,ib),
c     2                             bsig,q1,q2,bmed,back,irin,irout,kb,
c     3                             rmdrback,nrmdr,rmdrskew,
c     4                             dammmback,dammmmean,dammmmed,dammmskew,dammmsigma,nmmm,nUsed,
c     5                             ' -- '
c     6                             ,backbuf(1:kb),' ---- '
c                    endif
c
c=============== end of code seriously modified by TPC & JWF on B30606 =========
c
		else
		  back = LBACK(neighbors(nn),ifs,ib)
                  nUsed = 0
                  Skew  = -99.999
		endif
                if (GotTarget) write(6,'(a,i2,a,f8.2,a,i5,a,f7.3,a,f8.2)')
     +             'wpro_v6(1260): local background in band',
     +              ib,'=',back,'; nUsed =',nUsed,'; Skew =',Skew,
     +              '; BGsig =',BGsig                                ! JWF B30611
		fitpix = fitradii(ib)/pscale(ib)
	  	ifit = nint(fitpix + 0.5)
		ilo = max(i0-ifit,1)
		ihi = min(i0+ifit,nxorig(ib))
		jlo = max(j0-ifit,1)
		jhi = min(j0+ifit,nyorig(ib))
		rfitsq = fitpix**2
		troublersq = troublerpix**2
		do j = jlo,jhi
		do i = ilo,ihi
		    call checkpix(i,j,iframe,ib,nn,nvalues,ivalue,jvalue,
     *			ifvalue,maxpix,maxbands,newpix)
		    dx = i - Xpos(neighbors(nn),ifs,ib)
		    dy = j - Ypos(neighbors(nn),ifs,ib)
		    rsq = dx**2 + dy**2
		    itrans = i - xTRANS(ifs,ib)
		    jtrans = j - yTRANS(ifs,ib)
		    if(rsq <= rfitsq .and. itrans >= 1 .and. itrans <= nx(ib)
     *			.and. jtrans >= 1 .and. jtrans <= ny(ib) .and. newpix)
     *			then

! Count troublesome pixels even if they're masked out in later processing
		      if(back .gt. -900. .and. nvalues(ib)<maxpix .and.
     *                   goodcomp(nn) .and. nn .eq. icomp(1) .and. ncp .ne. 0 .and.
     *                   iand(Mask(itrans,jtrans,ifs_pix), trouble) /= 0) then
                        ! is primary component and one or more trouble bits set
			
			! Get pixel distance from refined frame pixel position
			dxt = i - x8set(nn,ifs,ib)
			dyt = j - y8set(nn,ifs,ib)
			rsqt = dxt**2 + dyt**2

                        if(rsqt .le. troublersq) then ! inside trouble radius
		          do ii=1,ntroublebits ! step through trouble bits
			    if(iand(Mask(itrans,jtrans,ifs_pix), 2**troublebits(ii)) /= 0)
     *                      then ! This trouble bit is set
			      if(mod(int(ftroubles(ib)/10**(ii-1)),10) .lt. 9) then
			        ftroubles(ib) = ftroubles(ib) + 10**(ii-1)
			      endif ! Counter not saturated
			      if(abs(dxt) .le. 0.5 .and. abs(dyt) .le. 0.5) then
                                ! Tougher test; looking for contamination of the wpro
                                ! position center pixel exactly
			        if(mod(int(ftroubles(ib)/10**(ii+2-1)),10) .lt. 9) then
				  ftroubles(ib) = ftroubles(ib) + 10**(ii+2-1)
			        endif ! Counter not saturated
			      endif ! On center pixel
			    endif ! This trouble bit is set
			  end do ! step through trouble bits
			endif ! inside trouble radius

		      endif ! is primary component and one or more trouble bits set

		      if(iand(Mask(itrans,jtrans,ifs_pix), fatal)==0 .and.
     *			    back > -900. .and. nvalues(ib)<maxpix) then
			nvalues(ib) = nvalues(ib) + 1
			ivalue(nvalues(ib),ib) = i
			jvalue(nvalues(ib),ib) = j
			dvalue(nvalues(ib),ib) = Array(itrans,jtrans,ifs_pix) - back
			sigdvalue0(nvalues(ib),ib) = Unc(itrans,jtrans,ifs_pix)
			vsum = 0.0d0
			do nnn = 1,nblend
			 if (goodcomp(nnn)) then
			  ipp = nint(pxc(ib) + (i - x8set(nnn,ifs,ib))*xprat)
			  jpp = nint(pyc(ib) + (j - y8set(nnn,ifs,ib))*yprat)
			  if (ipp>0 .and. ipp<=npsf(ib) .and. jpp>0 .and.
     *			   jpp<=npsf(ib) .and. offset(nnn,ifs,ib)>=0) vsum =
     *			    vsum + (rn*psfuncs(ipp,jpp,k,ib)*flux(ib,nnn))**2
			 endif
			enddo  ! nnn = 1,nblend
			varpsf(nvalues(ib),ib) = vsum +
     *			    Lconf(neighbors(nn),ifs,ib)**2
                        if (GotTarget) then
                           sumvarpsf(ib) = sumvarpsf(ib) + varpsf(nvalues(ib),ib)
                           print *,'(nn,neighbors(nn),ifs,ib,i,j):',
     +                     nn,neighbors(nn),ifs,ib,i,j,'; Lconf(neighbors(nn),ifs,ib):',
     +                     Lconf(neighbors(nn),ifs,ib)
                           sumLconf(ib) = sumLconf(ib) + Lconf(neighbors(nn),ifs,ib)
                        end if
			ifvalue(nvalues(ib),ib) = iframe
			if (nn==icomp(1)) primask(nvalues(ib),ib) = 1
			if (ncp /= 0 .and. iprioffscl >= 0) then
			    rsqpri = (i - x8pri)**2 + (j - y8pri)**2
			    if (rsqpri <= rfitsq) nvp(ib) = nvp(ib) + 1
			endif
		      endif
		    endif
		enddo ! i = ilo,ihi
		enddo ! j = jlo,jhi
	       endif ! onframe
	      enddo ! ib = 1,maxbands
	      enddo ! ifs = 1,nf
	    enddo     ! nn = 1,nblend0
	    do ib = 1,maxbands
		do nv = 1,nvalues(ib)
		    sigset(nv) = sigdvalue0(nv,ib)
		enddo
		call medsort(sigset,nvalues(ib),sigmed,q1,q2)
                medsig(ib) = sigmed                ! JWF B30419
		do nv = 1,nvalues(ib)
		    sigdvalue(nv,ib) = sqrt(sigmed**2 + varpsf(nv,ib))
		enddo
	    enddo

!---------------------------------------------------- ADDED 12/15/09

	    posconstrain = .false.
	    fret = func(p,psfs,npsize,npnmax,nDoF,GotTarget)
	    call framefunc(p,psfs,npsize,npnmax,n,onframe,nmax,nf,maxbands,
     *		wflag,nf,frameflux,framechi,GotTarget)
	    mblend = ntrue(goodcomp,nblend)

            if (GotTarget) then ! JWF dbg
             print *,'after recalculating local background for source #',n,' ...'  ! JWF dbg
             print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
             print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
             print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
             print *,'p:',p(1:2*nblend) ! JWF dbg
             print *,'racan:', racan(1:nblend) ! JWF dbg
             print *,'deccan:',deccan(1:nblend) ! JWF dbg
             print *,'nvalues:', nvalues           ! JWF B30509
             sigset = dvalue(1:nvalues(1),1)
             call tjsort(nvalues(1),sigset)
             print *,'dvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
             sigset = dvalue(1:nvalues(2),2)
             call tjsort(nvalues(2),sigset)
             print *,'dvalue(1:',nvalues(2),'2):', sigset(1:nvalues(2))
             sigset = dvalue(1:nvalues(3),3)
             call tjsort(nvalues(3),sigset)
             print *,'dvalue(1:',nvalues(3),'3):', sigset(1:nvalues(3))
             sigset = dvalue(1:nvalues(4),4)
             call tjsort(nvalues(4),sigset)
             print *,'dvalue(1:',nvalues(4),'4):', sigset(1:nvalues(4))
             sigset = sigdvalue(1:nvalues(1),1)
             call tjsort(nvalues(1),sigset)
             print *,'sigdvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
             sigset = sigdvalue(1:nvalues(2),2)
             call tjsort(nvalues(2),sigset)
             print *,'sigdvalue(1:',nvalues(2),'2):', sigset(1:nvalues(2))
             sigset = sigdvalue(1:nvalues(3),3)
             call tjsort(nvalues(3),sigset)
             print *,'sigdvalue(1:',nvalues(3),'3):', sigset(1:nvalues(3))
             sigset = sigdvalue(1:nvalues(4),4)
             call tjsort(nvalues(4),sigset)
             print *,'sigdvalue(1:',nvalues(4),'4):', sigset(1:nvalues(4))
             sigset = sigdvalue0(1:nvalues(1),1)
             call tjsort(nvalues(1),sigset)
             print *,'sigdvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
             sigset = sigdvalue0(1:nvalues(2),2)
             call tjsort(nvalues(2),sigset)
             print *,'sigdvalue0(1:',nvalues(2),'2):', sigset(1:nvalues(2))
             sigset = sigdvalue0(1:nvalues(3),3)
             call tjsort(nvalues(3),sigset)
             print *,'sigdvalue0(1:',nvalues(3),'3):', sigset(1:nvalues(3))
             sigset = sigdvalue0(1:nvalues(4),4)
             call tjsort(nvalues(4),sigset)
             print *,'sigdvalue0(1:',nvalues(4),'4):', sigset(1:nvalues(4))
             sigset = varpsf(1:nvalues(1),1)
             call tjsort(nvalues(1),sigset)
             print *,'varpsf(1:',nvalues(1),'1):', sigset(1:nvalues(1))
             sigset = varpsf(1:nvalues(2),2)
             call tjsort(nvalues(2),sigset)
             print *,'varpsf(1:',nvalues(2),'2):', sigset(1:nvalues(2))
             sigset = varpsf(1:nvalues(3),3)
             call tjsort(nvalues(3),sigset)
             print *,'varpsf(1:',nvalues(3),'3):', sigset(1:nvalues(3))
             sigset = varpsf(1:nvalues(4),4)
             call tjsort(nvalues(4),sigset)
             print *,'varpsf(1:',nvalues(4),'4):', sigset(1:nvalues(4))
            end if ! JWF dbg

! Final check for runaway fit.     [JWF: again only for RA:Dec]
	    nc = 0
	    ncp = 0
            if (.not.ok2) then
              print *,'WARNING (after act_wise): "(n==icp)" test never T on source #', n
              nOK2b = nOK2b + 1
            end if
            NotRunaway = .true.                                        ! JWF B30302
	    do nn = 1,nblend
	      if (goodcomp(nn)) then
		nc = nc + 1
		if (nn==icomp(1)) ncp = nc
		if (p(nc)**2+p(mblend+nc)**2 >= fitmaxsq) then         ! JWF B21207
                 if (WarnRunaway) then                                 ! JWF B30329
                   print *,'Post act_wise: discarding blend component',
     +                   nc,'distance:',sqrt(p(nc)**2+p(mblend+nc)**2) ! JWF B21207
                   print *,'max distance allowed:',sqrt(fitmaxsq)      ! JWF B30329
                   print *,'mdet source no.:',n                        ! JWF B30328
                   print *,'mdet RA, Dec:',racan(nn),deccan(nn)        ! JWF B30329
                  end if
                  if (nn==1) then                                      ! JWF B30302
                    nBadBlend = nBadBlend + 1                          ! JWF B21207
                    NotRunaway = .false.                               ! JWF B30302
                  end if                                               ! JWF B30302
	      	  do ib = 1,maxbands
		    flux(ib,nn) = 0.
		  enddo
		endif
	      endif
	    enddo

	    nbands = 0
	    do ib = 1,maxbands
		if(nvalues(ib) > 8) then
		    nbands = nbands + 1
		    iband(nbands) = ib
		endif
	    enddo

            if (TossZeroFlux) then
              if (.not.(any(flux(iband(1:nbands),icomp(1)) /= 0.))) then ! JWF B30207: Ken
                deallocate(p)                                            ! says toss entire
                deallocate(pmost)                                        ! solution if the
                if (NotRunaway) nAllZero(2) = nAllZero(2) + 1            ! primary component
                if (WarnZeroFlux) then                                   ! has all zero fluxes
                  print *,'Post act_wise: '
     +                  //'Discarding all-zero primary component for '
     +                  //'source',n,'; RA:Dec =',RAlist(n),Declist(n)
                  print *,'nbands:', nbands
                  print *,'nvalues:', nvalues
                  print *,'flux(iband(1:nbands),icomp(1)):', flux(iband(1:nbands),icomp(1))
                end if
                deallocate(p_preADB)
 	        go to 2000
              end if
            end if
c
! Update the passive blend list.
	    do nn = 1,npassive
     		npassprev(neighbors(icomp(nn))) =
     *		    max(npassprev(neighbors(icomp(nn))),npassive)
	    enddo

c	    nbands = 0
c	    do ib = 1,maxbands
c		if(nvalues(ib) > 8) then
c		    nbands = nbands + 1
c		    iband(nbands) = ib
c		endif
c	    enddo

	    nactiveg = 0
	    if (nactive > 0) then
	      do nn = nblend-nactive+1,nblend
		if (any(flux(iband(1:nbands),nn) > 0.)) nactiveg = nactiveg + 1
	      enddo
	      if (nactiveg > 0 .and.
     *		.not.any(flux(iband(1:nbands),icomp(1)) /= 0.))
     *		nactiveg = nactiveg - 1
	    endif

	    if ((ncp /= 0 .or. nactiveg /= 0) .and. nbands /= 0) then
	     finalflux = flux

             if (GotTarget) then ! JWF dbg
              print *,'after updating the passive blend list for source #',n,' ...'  ! JWF dbg
              print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
              print *,'p:',p(1:2*nblend) ! JWF dbg
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
              print *,'nactive, nactiveg:', nactive, nactiveg
              print *,'nvalues:', nvalues           ! JWF B30509
             end if ! JWF dbg


! Do the uncertainty calculation, based on full error model including
! variable Poisson noise.    !!!!!! JWF B30610: include background uncertainty??
	     do ib = 1,maxbands
	     do nv = 1,nvalues(ib)
		sigdvalue(nv,ib) = sqrt(sigdvalue0(nv,ib)**2 + varpsf(nv,ib))
	     enddo
	     enddo
             if (GotTarget) then
               print *,'sumvarpsf:', sumvarpsf(1:4)
               print *,'sumLconf: ', sumLconf(1:4)
               print *,'medsig: ', medsig
               print *,'nvalues:', nvalues             ! JWF B30509
               sigset = sigdvalue(1:nvalues(1),1)
               call tjsort(nvalues(1),sigset)
               print *,'sigdvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
               sigset = sigdvalue(1:nvalues(2),2)
               call tjsort(nvalues(2),sigset)
               print *,'sigdvalue(1:',nvalues(2),'2):', sigset(1:nvalues(2))
               sigset = sigdvalue(1:nvalues(3),3)
               call tjsort(nvalues(3),sigset)
               print *,'sigdvalue(1:',nvalues(3),'3):', sigset(1:nvalues(3))
               sigset = sigdvalue(1:nvalues(4),4)
               call tjsort(nvalues(4),sigset)
               print *,'sigdvalue(1:',nvalues(4),'4):', sigset(1:nvalues(4))
               sigset = sigdvalue0(1:nvalues(1),1)
               call tjsort(nvalues(1),sigset)
               print *,'sigdvalue(1:',nvalues(1),'1):', sigset(1:nvalues(1))
               sigset = sigdvalue0(1:nvalues(2),2)
               call tjsort(nvalues(2),sigset)
               print *,'sigdvalue0(1:',nvalues(2),'2):', sigset(1:nvalues(2))
               sigset = sigdvalue0(1:nvalues(3),3)
               call tjsort(nvalues(3),sigset)
               print *,'sigdvalue0(1:',nvalues(3),'3):', sigset(1:nvalues(3))
               sigset = sigdvalue0(1:nvalues(4),4)
               call tjsort(nvalues(4),sigset)
               print *,'sigdvalue0(1:',nvalues(4),'4):', sigset(1:nvalues(4))
               sigset = varpsf(1:nvalues(1),1)
               call tjsort(nvalues(1),sigset)
               print *,'varpsf(1:',nvalues(1),'1):', sigset(1:nvalues(1))
               sigset = varpsf(1:nvalues(2),2)
               call tjsort(nvalues(2),sigset)
               print *,'varpsf(1:',nvalues(2),'2):', sigset(1:nvalues(2))
               sigset = varpsf(1:nvalues(3),3)
               call tjsort(nvalues(3),sigset)
               print *,'varpsf(1:',nvalues(3),'3):', sigset(1:nvalues(3))
               sigset = varpsf(1:nvalues(4),4)
               call tjsort(nvalues(4),sigset)
               print *,'varpsf(1:',nvalues(4),'4):', sigset(1:nvalues(4))
               print *,'unsorted sigdvalue(1:',nvalues(1),'1):', sigdvalue(1:nvalues(1),1)
               print *,'unsorted sigdvalue(1:',nvalues(2),'2):', sigdvalue(1:nvalues(2),2)
               print *,'unsorted sigdvalue(1:',nvalues(3),'3):', sigdvalue(1:nvalues(3),3)
               print *,'unsorted sigdvalue(1:',nvalues(4),'4):', sigdvalue(1:nvalues(4),4)
               print *,'wpro_v6 about to call unspew'  ! JWF B30509
             end if

 	     call unspew(p,psfs,npsize,npnmax,varpsf,nlayers,n,onframe,
     *		nmax,nf,maxbands,wflag,nf,framesig,sigx,sigy,sigxy,sigflux,
     +          ok,GotTarget)       ! JWF B30510
c            print *,'back from unspew'  !  JWF dbg
	     fret = func(p,psfs,npsize,npnmax,nDoF,GotTarget)
	     flux = finalflux

             if (GotTarget) then ! JWF dbg
              print *,'after calling unspew and func for source #',n,' ...'  ! JWF dbg
              print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'sigflux:',(sigflux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
              print *,'p:',p(1:2*nblend) ! JWF dbg
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
              print *,'nactive, nactiveg:', nactive, nactiveg
             end if ! JWF dbg

!-----------------------------------------------------------------------
! For debugging purposes.

	     if (debug) then
	      icp = icomp(1)
	      if (ncp /= 0) then
                print *,'======================================================='
                print *,'(990) accessing p at ncp =',ncp ! JWF dbg
                print *,'(991) accessing p at mblend+ncp =',mblend+ncp ! JWF dbg
	  	write(6,'("Source#",i5,"  dRA =",f6.2," +/-",f6.2,'
     *		    //'"  dDec =",f6.2," +/-",f6.2," arcsec")')
     *		    n,3600.*p(ncp),sigx(icp),3600.*p(mblend+ncp),sigy(icp)
	  	if (ok) write(6,'(14x,"sigxy =",f6.2," arcsec")') sigxy(icp)
		print *,'        RA,Dec:',
     *		    racan(icp)+p(ncp)/cos(deccan(icp)*dtor),
     *		    deccan(icp)+p(mblend+ncp)
	  	print *,'          flux:',flux(1:4,icp)
	  	print *,'   (cf: apflux:',apflux(1:4,icp),')'
	  	print *,'       sigflux:',sigflux(1:4,icp)
	      else
		write(6,'("Source#",i5)') n
	      endif
	      if(nactive /= 0) then
		print *,'   Actively deblended components:'
		nc = mblend-nactive
		do ic = nblend-nactive+1,nblend
		 if (goodcomp(ic)) then
		  nc = nc + 1
                print *,'(1015) accessing p at nc =',nc ! JWF dbg
                print *,'(1016) accessing p at mblend+nc =',mblend+nc ! JWF dbg
	  	  write(6,'(8x,"Comp#",i1,"  dRA =",f6.2," +/-",f6.2,'
     *		    //'"  dDec =",f6.2," +/-",f6.2," arcsec")')
     *		    ic-nblend+1,3600.*p(nc),sigx(ic),
     *		    3600.*p(mblend+nc),sigy(ic)
	  	  print *,'          flux:',flux(1:4,ic)
	  	  print *,'       sigflux:',sigflux(1:4,ic)
		 endif
		enddo
	       endif
	       if (nf > 1) then
		print *,'   Flux estimates from individual frames:'
		  do ifs = 1,nf
		    if (any(onframe(n,ifs,1:4) .and. wflag(ifs,1:4)==1)) then
		      print *,'          ',ifs,': ',frameflux(ifs,1:4,icp)
		      print *,'             +/- ',framesig(ifs,1:4,icp)
		      print *,'           rchisq: ',framechi(ifs,1:4)
		    endif
		  enddo
	       endif
	       print *,'          rchisq =',rchisq
	       print *,'  Individual rch =',rchisqb
	       print *,'  function value =',fret,', reached in',iter,
     *		' iterations'
	       print *,'          #data values/band:',nvalues,
     *		'  ( total =',sum(nvalues),')'
	       print *,'          #primary pix/band:',nvp
	       print *,'          nblend =',max(npassprev(n),npassive)+nactiveg,
     *		'  cf: #components in this fit =',mblend
	       write(6,'(11x,"fraction of saturated pixels:"'//
     *		',4f7.4)') fsat
	       if (any(fsat>0)) write(6,
     *		'(11x,"sample #s of earliest saturation:",4i3)') satsamples
	       write(6,'(11x,"fraction of pixels with glitch masking:"'//
     *		',4f7.4)') ftroubles
	       write(6,'(8x,"fraction of pixels rejected for any reason:"'//
     *		',4f7.4)') freject
                print *,'======================================================='
	     endif ! (debug)
!-----------------------------------------------------------------------

! If the primary component of the blend is not the original candidate, then
! move the original candidate down the list so that it eventually gets
! processed.
	     npri = neighbors(icomp(1))
             IDpri = IDlist(npri)
	     if (npri /= n) then
                nBlendSwap = nBlendSwap + 1
                If (GotTarget .or. (nBlendSwap .le. MaxWarnBlendSwap)) then
                  print *,'wpro_v6(1173): blend swap for source',n
                  print *,'               new primary is source',npri
                  print *,'wpro_v6(1174): mblend:',mblend
                  print *,'wpro_v6(1175): icomp:', icomp(1:mblend)
                  print *,'wpro_v6: old primary mdet position:',
     +                     racan(1), deccan(1)
                  print *,'wpro_v6: new primary mdet position:',
     +                     racan(icomp(1)), deccan(icomp(1))
                  if (GotTarget) then
                    GotTargetK(kTarget) = .true.
                    GotAllTargets = .true.
                    do 8 kk = 1, 10
                      if ((TargetRA(kk) .le. 360.0d0)
     +                      .and. .not.GotTargetK(kk)) then
                         GotAllTargets = .false.
                      end if
8                   continue
                  end if
                end if

                Swapped    = .true.
	       	satsamples = 0
		RAlist(npri) = RAlist(n)
		DEClist(npri) = DEClist(n)
                IDlist(npri) = IDlist(n)
		do ib = 1,4
		    MAGSTD(npri,ib) = MAGSTD(n,ib)
		    eMAGSTD(npri,ib) = eMAGSTD(n,ib)
		    Rsat(npri,ib) = Rsat(n,ib)
	          do ifs = 1,nf
		    Xpos(npri,ifs,ib) = Xpos(n,ifs,ib)
		    Ypos(npri,ifs,ib) = Ypos(n,ifs,ib)
		    LBACK(npri,ifs,ib) = LBACK(n,ifs,ib)
		    LSIG(npri,ifs,ib) = LSIG(n,ifs,ib)
		    onframe(npri,ifs,ib) = onframe(n,ifs,ib)
		  enddo
		enddo
	     endif
             if (GotTarget) then ! JWF dbg
              print *,'after possible blend swap for source #',n,' ...'  ! JWF dbg
              print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'sigflux:',(sigflux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
              print *,'p:',p(1:2*nblend) ! JWF dbg
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
              print *,'nactive, nactiveg:', nactive, nactiveg
              print *,'starting Table-load loop for non-PM'
             end if ! JWF dbg

! Add the results to output table and update frame-dependent information.
             m0 = m
	     nc = 0
             racan0  = racan(icomp(1))  ! JWF B30502
             deccan0 = deccan(icomp(1)) ! JWF B30502
	     do nn = 1,nblend		! Begin loop over components. --------------------------------------------------------
              if (GotTarget) print *,'checking component no.',nn,'for non-PM Table load'                                     !
              DoThisOne(nn) = .false.                                                                                        !
	      if (goodcomp(nn)) then   !--------------------------------------------------------------------------------     !
	       nc = nc + 1                                                                                             !     !
               if (.not.(any(flux(iband(1:nbands),nn) /= 0.))) nAllZero(3) = nAllZero(3) + 1                           !     !
	       if ((nn==icomp(1) .or. (nn>nblend-nactive)) .and.                                                       !     !
     *	       any(flux(iband(1:nbands),nn) /= 0.)) then     !---------------------------------------------------      !     !
		m = m+1                                                                                         !      !     !
                mf = m                                                                                          !      !     !
                DoThisOne(nn) = .true.                                                                          !      !     !
		RAlistnew(m)  = racan(nn) + p(nc)/cos(deccan(nn)*dtor) ! NOTE: check for 0-360 error            !      !     !
		DEClistnew(m) = deccan(nn) + p(mblend+nc)              ! and past-the-pole error TBD            !      !     !
                IDlistnew(m)  = IDpri                                                                           !      !     !
                if (DumPvec) write (76,'(3F10.5,2F12.7,2I7,I4)') 3600.0*p(nc), 3600.0*p(mblend+nc),             !      !     !
     +             3600.0*sqrt(p(nc)**2+p(mblend+nc)**2), RAlistnew(m), DEClistnew(m), iter, nsteps, nblend0    !      !     !
		Table(m,1)  = RAlistnew(m)                                                                      !      !     !
		Table(m,2)  = DEClistnew(m)                                                                     !      !     !
		Table(m,3)  = sigx(nn)                                                                          !      !     !
		Table(m,4)  = sigy(nn)                                                                          !      !     !
                Table(m,64) = 3600.0*p(nc)                                                                      !      !     !
                if (Table(m,64) .gt.  99.99999) Table(m,64) =  99.99999                                         !      !     !
                if (Table(m,64) .lt. -99.99999) Table(m,64) = -99.99999                                         !      !     !
                Table(m,65) = 3600.0*p(mblend+nc)                                                               !      !     !
                if (Table(m,65) .gt.  99.99999) Table(m,65) =  99.99999                                         !      !     !
                if (Table(m,65) .lt. -99.99999) Table(m,65) = -99.99999                                         !      !     !
		if (ok) Table(m,5) = sigxy(nn)                                                                  !      !     !
                if (GotTarget) then                                                                             !      !     !
                  print *,'Loading non-PM Table slot',m,'for source no.',n                                      !      !     !
                  print *,'nn, nc:',nn,nc                                                                       !      !     !
                  print *,'RA, Dec:',RAlistnew(m),DEClistnew(m)                                                 !      !     !
                  print *,'sigra, sigdec, sigradec:', Table(m,3:5)                                              !      !     !
                  print *,'sigfac:',sigfac                                                                      !      !     !
                end if                                                                                          !      !     !
		do i = 1,nbands         !-----------------------------------------------------------------      !      !     !
		  if (rchisqb(iband(i)) /= 0. .and. nvp(iband(i)) /= 0) then  !-----------------------   !      !      !     !
		    Table(m,5+iband(i)) = flux(iband(i),nn)                                          !   !      !      !     !
		    Table(m,9+iband(i)) = sigflux(iband(i),nn)*sigfac(iband(i))                      !   !      !      !     !
		    Table(m,13+iband(i)) = rchisqb(iband(i))                                         !   !      !      !     !
                    if (GotTarget) then                                                              !   !      !      !     !
                      print *,'band',iband(i),':'                                                    !   !      !      !     !
                      print *,'flux, sigflux, rchisqb:',Table(m,5+iband(i)),                         !   !      !      !     !
     +                         Table(m,9+iband(i)),Table(m,13+iband(i))                              !   !      !      !     !
                    end if                                                                           !   !      !      !     !
		    Table(m,20+iband(i)) = ftroubles(iband(i))                                        !   !      !      !     !
		    Table(m,24+iband(i)) = fsat(iband(i))                                            !   !      !      !     !
		    Table(m,59+iband(i)) = fitradii(iband(i))                                        !   !      !      !     !
		    if (Table(m,59+iband(i)) .gt. 999.99) Table(m,59+iband(i)) = 999.99              !   !      !      !     !
		    SatNum(m,iband(i)) = satsamples(iband(i))                                        !   !      !      !     !
		    do ifs = 1,nf   !-----------------------------------------------------------     !   !      !      !     !
			if(onframe(n,ifs,iband(i)) .and.                                       !     !   !      !      !     !
     *					wflag(ifs,iband(i))==1) then   !------------------     !     !   !      !      !     !
			    ff = frameflux(ifs,iband(i),nn)                              !     !     !   !      !      !     !
			    fs = framesig(ifs,iband(i),nn)                               !     !     !   !      !      !     !
                            if (fs .ge. WarnBigFrameSig) then                            !     !     !   !      !      !     !
                              print *,'WARNING: framesig out of limits'                  !     !     !   !      !      !     !
                              print *,'         ifs,iband(i),nn:',ifs,iband(i),nn        !     !     !   !      !      !     !
                              print *,'         neighbors(nn):',neighbors(nn)            !     !     !   !      !      !     !
                              print *,'         frameflux,framesig:', ff,fs              !     !     !   !      !      !     !
                              print *,'         racan(nn),deccan(nn)',                   !     !     !   !      !      !     !
     +                                          racan(nn),deccan(nn)                     !     !     !   !      !      !     !
                              print *,'         RAlistnew,Declistnew:',                  !     !     !   !      !      !     !
     +                                          RAlistnew(m),Declistnew(m)               !     !     !   !      !      !     !
                            end if                                                       !     !     !   !      !      !     !
			    fc = framechi(ifs,iband(i))                                  !     !     !   !      !      !     !
			    if (ff /= 0. .and. fs /= 0.) then   !----------------        !     !     !   !      !      !     !
				framefluxes(m,ifs,iband(i)) = ff                !        !     !     !   !      !      !     !
				frameuncs(m,ifs,iband(i)) = fs                  !        !     !     !   !      !      !     !
				framechis(m,ifs,iband(i)) = fc                  !        !     !     !   !      !      !     !
				if (GotTarget) then                             !        !     !     !   !      !      !     !
                                  print *,'frame',ifs,'flux,unc,chi2,S/N:',     !        !     !     !   !      !      !     !
     +                                     ff,fs,fc,ff/fs                       !        !     !     !   !      !      !     !
                                end if                                          !        !     !     !   !      !      !     !
			    endif ! (ff /= 0. .and. fs /= 0.)  !-----------------        !     !     !   !      !      !     !
			endif ! (onframe(n,ifs,iband(i))) &c. !---------------------------     !     !   !      !      !     !
		    enddo  ! ifs = 1,nf    !----------------------------------------------------     !   !      !      !     !
		  endif ! (rchisqb(iband(i)) /= 0. .and. nvp(iband(i)) /= 0) !------------------------   !      !      !     !
		enddo ! i = 1,nbands   !------------------------------------------------------------------      !      !     !
                if (UseNonPMposns) then                                                                         !      !     !
                  racan(nn)  = RAlistnew(m)      ! JWF B30220 TBD  set up for PM to start with these            !      !     !
                  deccan(nn) = Declistnew(m)     ! JWF B30220 TBD         and zero p values                     !      !     !
                end if                                                                                          !      !     !
		Table(m,18) = rchisq                                                                            !      !     !
		Table(m,19) = max(npassprev(neighbors(icomp(1))),npassive) +                                    !      !     !
     *		    nactiveg                                                                                    !      !     !
		Table(m,20) = nactiveg                                                                          !      !     !
                if (SingleFrame) then                      ! JWF B31206                                         !      !     !
                  Table(m,59) = nblend0                                                                         !      !     !
                  if (Swapped) Table(m,59) = -Table(m,59)                                                       !      !     !
                end if                                                                                          !      !     !
c                                                                                                               !      !     !
c========================= code added by JWF on B21120 ============================                             !      !     !
                Table(m,55) = iter                                                                              !      !     !
                Table(m,56) = nsteps                                                                            !      !     !
               end if ! ((nn==icomp(1) .or. (nn>nblend-nactive)) .and.any(flux(iband(1:nbands),nn) /= 0.))-------      !     !
              end if ! (goodcomp(nn)) ----------------------------------------------------------------------------------     !
             end do ! nn = 1,nblend  (End of loop over components for nonPM)--------------------------------------------------
             flux_nonPM    = flux            ! JWF B30208 - save for frame subtraction
             sigflux_nonPM = sigflux         ! JWF B30208 - save for frame subtraction
             nblend_nonPM  = nblend          ! JWF B30325 - save for frame subtraction
             if (GotTarget) then
               print *,'finished Table-load for non-PM; nblend,mblend:',nblend, mblend
               print *,'flux_nonPM =',flux_nonPM(1:4,1:nblend)
               print *,'sigflux_nonPM =',sigflux_nonPM(1:4,1:nblend)
               print *,'nblend_nonPM =',nblend_nonPM
             end if
c
c   Compute the mean observational epoch      ! JWF B30207
c
             jj = MeanObsEpochType
             if (SingleFrame) jj = 1          ! JWF B31122
             if (jj .eq. 0) then
               JD0 = MJD0
               go to 70
             end if
9            sumJD = 0.0d0
             Select Case (jj)
c
             Case(1)                           ! unweighted
              do 10 kk = 1, nframes
                sumJD = sumJD + JD(ifq(kk),ibq(kk))
10            continue
              JD0 = sumJD/dfloat(nframes)
c
             Case(2)                           ! flux weighting
              sumWgtJD = 0.0d0
              nFusable = 0
              do 20 kk = 1, nframes
                WgtJD = frameflux(ifq(kk),ibq(kk),icomp(1))
                if (WgtJD .gt. 0.0) then
                  sumJD    = sumJD + JD(ifq(kk),ibq(kk))*WgtJD
                  sumWgtJD = sumWgtJD + WgtJD
                  nFusable = nFusable + 1
                end if
20            continue
              if (nFusable .gt. nframes/2) then
                JD0 = sumJD/sumWgtJD
              else
                nWgtFail = nWgtFail + 1
                jj = 1
                go to 9
              end if
c
             Case(3)                           ! S/N weighting
              sumWgtJD = 0.0d0
              nFusable = 0
              do 30 kk = 1, nframes
                Dtmp1 = frameflux(ifq(kk),ibq(kk),icomp(1))
c                print *,'wpro_v6(1294): kk, Dtmp1:', kk, Dtmp1 ! JWF dbg
                if (Dtmp1 .le. 0.0) go to 30
                Dtmp2 = framesig(ifq(kk),ibq(kk),icomp(1))
c                print *,'wpro_v6(1297): kk, Dtmp2:', kk, Dtmp2 ! JWF dbg
                if (Dtmp2 .le. 0.0) go to 30
                WgtJD    = Dtmp1/Dtmp2
                sumJD    = sumJD + JD(ifq(kk),ibq(kk))*WgtJD
                sumWgtJD = sumWgtJD + WgtJD
                nFusable = nFusable + 1
c                print *,'wpro_v6(1303): kk,sumWgtJD,nFusable:',kk,sumWgtJD,nFusable ! JWF dbg
c                print *,'Dtmp1,Dtmp2,WgtJD:',Dtmp1,Dtmp2,WgtJD ! JWF dbg
30            continue
              if (nFusable .gt. nframes/2) then
                JD0 = sumJD/sumWgtJD
c                print *,'wpro_v6(1308): sumJD,sumWgtJD,JD0:',sumJD,sumWgtJD,JD0 ! JWF dbg
              else
                nWgtFail = nWgtFail + 1
                jj = 1
                go to 9
              end if
c
             Case(4)                           ! inverse variance weighting
              sumWgtJD = 0.0d0
              nFusable = 0
              do 40 kk = 1, nframes
                Dtmp1 = framesig(ifq(kk),ibq(kk),icomp(1))
                if (Dtmp1 .le. 0.0) go to 40
                WgtJD    = 1.0d0/Dtmp1**2
                sumJD    = sumJD + JD(ifq(kk),ibq(kk))*WgtJD
                sumWgtJD = sumWgtJD + WgtJD
                nFusable = nFusable + 1
40            continue
              if (nFusable .gt. nframes/2) then
                JD0 = sumJD/sumWgtJD
              else
                nWgtFail = nWgtFail + 1
                jj = 1
                go to 9
              end if
c
             Case(5)                           ! inverse chi-square weighting
              sumWgtJD = 0.0d0
              nFusable = 0
              do 50 kk = 1, nframes
                Dtmp1 = framechi(ifq(kk),ibq(kk))
                if (Dtmp1 .le. 0.0) go to 50
                WgtJD    = 1.0d0/Dtmp1
                sumJD    = sumJD + JD(ifq(kk),ibq(kk))*WgtJD
                sumWgtJD = sumWgtJD + WgtJD
                nFusable = nFusable + 1
50            continue
              if (nFusable .gt. nframes/2) then
                JD0 = sumJD/sumWgtJD
              else
                nWgtFail = nWgtFail + 1
                jj = 1
                go to 9
              end if
c
             End Select
c
70           do 80 kk = 1, nframes
               tframe(kk) = JD(ifq(kk),ibq(kk)) - JD0
80           continue
c
             if (SingleFrame) then            ! JWF B31122
               MJD0 = JD0                     ! JWF B31122
               go to 999                      ! JWF B31122
             end if                           ! JWF B31122
c
! Calculate maximum likelihood PM solution.
             if (GotTarget) then
              print *,'wpro_v6(1361): mblend =',mblend,'; nblend =',nblend,'; nblend_preADB =',nblend_preADB ! JWF dbg B30216
              print *,'wpro_v6(1362): p(1:2*mblend) =',p(1:2*mblend)           ! JWF dbg B30216
              print *,'wpro_v6(1363): p_preADB(1:2*nblend_preADB) =',p_preADB(1:2*nblend_preADB)  ! JWF dbg B30216
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
             end if
             if (PMinitADB) then
              nparms = 2*mblend + 2	! the extra two parameters represent
c					! the proper motion of the primary,
c					! dRAcosDec/dt, dDec/dt in deg/day
c                                       ! CORRECTION (JWF, B21127): deg/yr;
c                                       ! values in p array must be divided
c                                       ! by 365.25 when multiplying by time
c                                       ! interval in units of Julian Days
              allocate (p_PM(nparms))
              p_PM = 0.0
              do 100 k = 1, 2*mblend
                p_PM(k) = p(k)
100           continue
              if (GotTarget) then
                print *,'PMinitADB: p_PM =',p_PM(1:2*mblend+2)  ! JWF dbg
              end if
             else
              nblend   = nblend_preADB
              goodcomp = goodcomp_preADB   ! JWF B20226: func_pm and dfunc_pm need this
              flux     = flux_preADB       ! JWF B30226
              icomp    = icomp_preADB
              icp      = icp_preADB
              mblend   = ntrue(goodcomp,nblend)
              nparms   = 2*mblend + 2
              allocate (p_PM(nparms))
              p_PM = 0.0
              do 110 k = 1, 2*mblend
                p_PM(k) = p_preADB(k)
110           continue
              if (GotTarget) then
                print *,'.not.PMinitADB: p_PM =',p_PM(1:2*mblend+2)  ! JWF dbg
              end if
             end if
             p_PM(1)        = 0.0  ! start PM with racan & deccan updated during
             p_PM(mblend+1) = 0.0  ! non-PM Table loading instead of p values
c
             deallocate(p)
             deallocate(pmost)
             deallocate(p_preADB)
	     allocate (p(nparms))
	     allocate (pmost(nparms))
             p = p_PM
             deallocate(p_PM)
             if (GotTarget) then
              print *,'wpro_v6(1409): mblend =',mblend,'; nparms =',nparms ! JWF dbg
              print *,'wpro_v6(1410): p(1:2*mblend+2) =',p(1:2*mblend+2) ! JWF dbg
              print *,'p(1:2*mblend+2) =',p(1:2*mblend+2) ! JWF dbg
              print *,'goodcomp:',goodcomp(1:nblend)      ! JWF dbg
              print *,'icomp:', icomp(1:nblend)           ! JWF dbg
              print *,'flux:', flux(1:4,1:nblend)         ! JWF dbg
              print *,'racan:', racan(1:nblend)           ! JWF dbg
              print *,'deccan:',deccan(1:nblend)          ! JWF dbg
             end if
c
	     ftol = ftol_pm    ! JWF B30226
	     fret = xpixscale(1)/PMstepFac ! JWF B30524
c	     posconstrain = .true.         ! JWF B30311
	     posconstrain = poscon0        ! JWF B30311
             if (InitPMflux0) flux = 0.0   ! JWF B30312
             if (InitPMposn0) p    = 0.0   ! JWF B30312
             if (GotTarget) then
              print *,'wpro_v6: about to call most_wise_pm....'            ! JWF dbg
              print *,'wpro_v6: using ftol =',ftol
              print *,'p(1:2*mblend+2) =',p(1:2*mblend+2) ! JWF dbg
              print *,'flux:', flux(1:4,1:nblend)  ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
             end if
c
 	     call most_wise_pm(p,nparms,psfs,npsize,npnmax,ftol,iter,fret,
     +                         maxsteps,nsteps,ok2,GotTarget) ! JWF B30117
             if (GotTarget) then
              print *,'wpro_v6: back from most_wise_pm'   ! JWF dbg
              print *,'wpro_v6(1423): mblend =',mblend,'; nparms =',nparms ! JWF dbg
              print *,'wpro_v6(1424): p(1:2*mblend+2) =',p(1:2*mblend+2) ! JWF dbg
              print *,'flux:', flux(1:4,1:nblend)  ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
             end if
             if (.not.ok2) Nok2b = Nok2b + 1
             if (IzBad(p(2*mblend+1)) .or. IzBad(p(2*mblend+2))) then
               nPMNaN = nPMNaN + 1
               if (nWarnPMNaN .lt. MaxWarnPMNaN)then                    ! JWF B30502
                 nWarnPMNaN = nWarnPMNaN + 1                            ! JWF B30502
                 print *,'WARNING(wpro_v6): NaN PM solution on source', ! JWF B30502
     +                    n,'; mdet RA, Dec:',racan0,deccan0            ! JWF B30502
               end if                                                   ! JWF B30502
             end if                                                     ! JWF B30502
	     posconstrain = .false.        ! JWF B30311
	     fret = func_pm(p,psfs,npsize,npnmax,nDoF_pm,GotTarget)
	     pmost = p
             finalflux = flux

             if (GotTarget) then ! JWF dbg
              print *,'about to call unspew_pm for source #',n,' ...'  ! JWF dbg
              print *,'mblend =',mblend                                ! JWF dbg
              print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
              print *,'p:',p(1:2*mblend+2) ! JWF dbg
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
              print *,'nactive, nactiveg:', nactive, nactiveg
             end if ! JWF dbg
c
             call unspew_pm(p,psfs,npsize,npnmax,varpsf,nlayers,n,onframe,
     *		nmax,nf,maxbands,wflag,nf,framesig,sigx,sigy,sigxy,sigflux,
     *		sigpmr,sigpmd,covRpmR,covDpmD,covpmRD,ok,GotTarget)   ! JWF B21221
c            print *,'back from unspew_pm'  !  JWF dbg
	     fret = func_pm(p,psfs,npsize,npnmax,nDoF_pm,GotTarget)
	     flux = finalflux

             if (GotTarget) then ! JWF dbg
              print *,'after calling unspew_pm and func_pm for source #',n,' ...'  ! JWF dbg
              print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'icomp:',(icomp(nn),nn = 1, nblend) ! JWF dbg
              print *,'flux:',(flux(1:4,nn), nn = 1, nblend) ! JWF dbg
              print *,'finalflux:', finalflux(1:4,1:nblend)  ! JWF dbg
              print *,'p:',p(1:2*mblend+2) ! JWF dbg
              print *,'racan:', racan(1:nblend) ! JWF dbg
              print *,'deccan:',deccan(1:nblend) ! JWF dbg
              print *,'nactive, nactiveg:', nactive, nactiveg
             end if ! JWF dbg
c
! Add the PM results to output table and update frame-dependent information.
             m = m0
	     nc = 0
999	     do nn = 1, nblend_nonPM  		! Begin loop over components.
              if (GotTarget) print *,'checking component no.',nn,'for PM Table load'
	      if (goodcomp(nn) .or. (nn .gt. nblend)) then
               nc = nc + 1
	       if (((nn==icomp(1)) .or. (DoThisOne(nn)) .or. (nn .gt. nblend))
     +            .and. (any(flux_nonPM(iband(1:nbands),nn) /= 0.))) then
                if (SingleFrame) go to 1000       ! JWF B31206
		m = m+1
                if (GotTarget) then
                  print *,'loading PM data into Table slot',m
                  print *,'nn, nc:',nn,nc
                  if (nn .gt. nblend) then
                    print *,'RA, Dec:',RAlistnew(m),DEClistnew(m)
                    print *,'This is an active deblend, jumping to frame subtraction'
                  else
                    print *,'RA, Dec:',racan(nn)+p(nc)/cos(deccan(nn)*dtor),deccan(nn)+p(mblend+nc)
                  end if
                end if
 		do i = 1,nbands
		  if (rchisqb(iband(i)) /= 0. .and. nvp(iband(i)) /= 0) then
                    if (tdbg) then
                      Table(m,41+iband(i)) = nvalues(iband(i))
                      Table(m,45+iband(i)) = nlayers(iband(i))
                    else
                      Table(m,41+iband(i)) = flux(iband(i),nn)
                      Table(m,45+iband(i)) = sigflux(iband(i),nn)*sigfac(iband(i))
                    end if
		    Table(m,49+iband(i)) = rchisqb(iband(i))
                    if (GotTarget) then
                      print *,'band',iband(i),':'
                      print *,'flux, sigflux, rchisqb:',
     +                 Table(m,41+iband(i)),Table(m,45+iband(i)),Table(m,49+iband(i))
                    end if
		  endif
		enddo ! i = 1,nbands
                if (tdbg) then
                  Table(m,52) = nDoF
                  Table(m,53) = nDoF_pm
                end if
		Table(m,54) = rchisq
c
c========================= code added by JWF on B21120 ============================
                pmra  = 0.0
                pmdec = 0.0
                if (nn .ne. icomp(1)) then      ! this is an active deblend
                  if (GotTarget) print *,'(1149) took active-deblend branch, setting PM null' ! JWF dbg
                  Table(m,29) = RBlank
                  Table(m,30) = RBlank
                  Table(m,31) = RBlank
                  Table(m,32) = RBlank
                  Table(m,35) = sigx(nn)
                  Table(m,36) = sigy(nn)
                  Table(m,59) = 0.0
                  go to 1000
                else                             ! primary component
c========================= end of code added by JWF on B21120 =====================
c
                  if (GotTarget) print *,'(1149) took primary-component branch' ! JWF dbg
                  Table(m,29) = 3600.*p(2*mblend+1)	! RA PM  JWF B21121
		  Table(m,30) = 3600.*p(2*mblend+2)	! Dec PM  JWF B21121
c
                  cd  = cos(deccan(nn)*dtor)
                  dRA = abs(racan(nn)+p(nc)/cd-RAlistnew(m))
                  if (dRA .ge. 180.0) dRA = dRA - sign(360.0,dRA)
                  dist = 3600.0*sqrt((dRA*cd)**2 + (deccan(nn)+p(mblend+nc)-Declistnew(m))**2)
                  if (GenJD0posns) then   ! JWF B30221: NOTE works best if ireg=jreg=1
                    pmra  = Table(m,29)
                    pmdec = Table(m,30)
                    if (dist  .gt.  99.999) dist  =  99.999
                    if (pmra  .lt. -99.999) pmra  = -99.999
                    if (pmra  .gt.  99.999) pmra  =  99.999
                    if (pmdec .lt. -99.999) pmdec = -99.999
                    if (pmdec .gt.  99.999) pmdec =  99.999
                    write(75,'(4F12.7,F10.4,F11.4,F13.4,F7.3,2F8.3)')
     +                RAlistnew(m), DEClistnew(m),
     +                racan(nn)+p(nc)/cd,
     +                deccan(nn)+p(mblend+nc), sigx(nn), sigy(nn), sigxy(nn),
     +                dist, pmra, pmdec
                  end if
c
                  nPMsolns    = nPMsolns + 1
		  Table(m,31) = sigpmr			! RA PM unc
		  Table(m,32) = sigpmd			! Dec PM unc
                  pmra  = p(2*mblend+1)/365.25d0     ! JWF B21127  [deg/day]
                  pmdec = p(2*mblend+2)/365.25d0     ! JWF B21127  [deg/day]
                  Table(m,35) = sqrt(sigx(nn)**2 + ((MJD0-JD0)*sigpmr/365.25)**2
     +                                           + 2.0*(MJD0-JD0)*covRpmR)
                  Table(m,36) = sqrt(sigy(nn)**2 + ((MJD0-JD0)*sigpmd/365.25)**2
     +                                           + 2.0*(MJD0-JD0)*covDpmD)
c
                  if (DumPvec) write (77,'(6F10.5,2F12.7,2I10,I4)')
     +               3600.0*p(nc), 3600.0*p(mblend+nc),
     +               3600.0*sqrt(p(nc)**2+p(mblend+nc)**2),
     +               Table(m,29), Table(m,30), sqrt(Table(m,29)**2+Table(m,30)**2),
     +               racan(nn)+p(nc)/cd,
     +               deccan(nn)+p(mblend+nc), iter, nsteps, nblend0
c
                end if ! (nn .ne. icomp(1))
                if (GotTarget) then                        ! JWF dbg
                  print *,'RA and Dec for non-PM:', RAlistnew(m), DEClistnew(m) ! JWF dbg
                  print *,'RA and Dec for obs-epoch:', racan(nn)+p(nc)/cd,
     +                                                 deccan(nn)+p(mblend+nc) ! JWF dbg
                  print *,'Mean Obs Epoch:', JD0 ! JWF dbg
                end if                                            ! JWF dbg
                Table(m,33) = int(JD0)                  ! JWF B30107
                Table(m,34) = JD0 - int(JD0)            ! JWF B30503
                Table(m,37) = 0.0
   	        if (ok) Table(m,37) = sigxy(nn)  ! JWF B21222 [this is done above,
c                                                  assume should be here too];
c                                                  NOTE: covpmRD negligible
                if (GotTarget) then                            ! JWF dbg
                  print *,'sigRA, SigDec, SigRADec for obs-epoch:',
     +                     sigx(nn),sigy(nn),sigxy(nn)         ! JWF dbg
                end if                                         ! JWF dbg
c   NOTE: position here is at MeanObsEpoch; transform to RefEpoch is in writeTABLE
                Table(m,38) = int(100.0*(racan(nn) + p(nc)/cd))
                Table(m,39) =     100.0*(racan(nn) + p(nc)/cd)
     +                      - int(100.0*(racan(nn) + p(nc)/cd))
                Table(m,40) = int(100.0*(deccan(nn) + p(mblend+nc)))
                Table(m,41) =     100.0*(deccan(nn) + p(mblend+nc))
     +                      - int(100.0*(deccan(nn) + p(mblend+nc)))
c
                Table(m,57) = iter
                Table(m,58) = nsteps
                if (dist .gt. 9.99) then
                  if (NdistWarn .lt. 10) print *,'WARNING(wpro_v6): dist =',dist
                  NdistWarn = NdistWarn + 1
                  if (NdistWarn .eq. 10) print *,'wpro_v6: last warning of this type'
                  dist = 9.99
                end if

                Table(m,59) = nblend0 + 0.001*dist
                if (Swapped) Table(m,59) = -Table(m,59)
c
c w?mpro_pm, w?sigmpro_pm, w?snr_pm will be done in writeTABLE
c
c====== resume with flux subtraction from frames ========================================
1000            iframe = 0
		do ifs = 1,nf
		do ib = 1,maxbands
	 	  if(onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1) then
                    ifs_pix = pix_order(ifs,ib)
		    iframe = iframe + 1
	    	    iwcs = WCS(ifs,ib)
        	    offscl = -1
c        	    call wcs2pix(iwcs, RAlistnew(m), DEClistnew(m),  ! JWF B21127
c    *			x8, y8, offscl)                              ! JWF B21127
c========================= code added by JWF on B21127 ============================
                    if (PMsubtract .and. (nn==icomp(1))) then        ! JWF B30126
                      tframe(iframe) = JD(ifs,ib) - JD0
                      dRA  = pmra*tframe(iframe)/dcos(dtor*DEClistnew(m))
                      dDec = pmdec*tframe(iframe)
                    else                                             ! JWF B30126
                      dRA  = 0.0                                     ! JWF B30126
                      dDec = 0.0                                     ! JWF B30126
                    end if                                           ! JWF B30126
                    call wcs2pix(iwcs,
     +                           RAlistnew(m)+dRA,
     +                           DEClistnew(m)+dDec,
     +                           x8, y8, offscl)
c     if (n .eq. 100) print *,'(1199): dRA,dDec = ',dRA, dDec,
c    + ' on frame',ifs  ! JWF dbg
c========================= end of code added by JWF on B21127 =====================
        	    if (offscl >= 0) then
			Xposnew(m,ifs,ib) = x8
			Yposnew(m,ifs,ib) = y8
	    	    endif
                    if (GotTarget) then
                      print *,'frame subtract: m,ifs,ib =',m,ifs,ib
                      print *,'Xposnew(m,ifs,ib),Yposnew(m,ifs,ib):',
     +                         Xposnew(m,ifs,ib),Yposnew(m,ifs,ib)
                    end if
		    LBACKnew(m,ifs,ib) = LBACK(n,ifs,ib)
		    LSIGnew(m,ifs,ib) = LSIG(n,ifs,ib)
		    Rsatnew(m,ib) = Rsat(n,ib)
!
! Subtract estimated contribution of source from observed image and update
! uncertainty map.
		    if(flux_nonPM(ib,nn)/sigflux_nonPM(ib,nn) > 2. .and. nvp(ib) /= 0) then
                      if (GotTarget) then
                        print *,'subtracting flux',flux_nonPM(ib,nn),'in band',ib,
     +                          'for source nn =',nn,'in frame ifs =',ifs
                        if (WrtFSimg) then                                   ! JWF B30824
                          write(BndStr,'(I1)') ib                            ! JWF B30824
                          write(TrgStr,'(I2.2)') kTarget                     ! JWF B30824
                          write(NumStr,*) ifs                                ! JWF B30824
                          call trimblnk(NumStr)                              ! JWF B30824
                          FSimgNam0 = 'FSimage_'//TrgStr//'_W'//BndStr//'_fr'! JWF B30824
     +                                         //NumStr(1:lnblnk(NumStr))    ! JWF B30824
                          blocksize=1                                        ! JWF B30824
                          simple=.true.                                      ! JWF B30824
                          bitpix=-32                                         ! JWF B30824
                          naxis=2                                            ! JWF B30824
                          extend=.false.                                     ! JWF B30824
                          group=1                                            ! JWF B30824
                          fpixel=1                                           ! JWF B30824
                          kpix = 0                                           ! JWF B30824
                          npix = 0                                           ! JWF B30824
                        end if                                               ! JWF B30824
                      end if
		      ifp = min(max(nint(x8), 1), nxfp)
		      jfp = min(max(nint(y8), 1), nyfp)
		      k = MapPSFs(ifp,jfp,ib)
	    	      xprat = xpixscale(iframe)/psamp(ib)
	    	      yprat = ypixscale(iframe)/psamp(ib)
		      ipsfspan = nint(npsf(ib)/xprat)
		      jpsfspan = nint(npsf(ib)/yprat)
		      ipos = nint(x8)
		      jpos = nint(y8)
		      ilo = max(ipos-ipsfspan/2, 1)
		      ihi = min(ipos+ipsfspan/2, nxorig(ib))
		      jlo = max(jpos-jpsfspan/2, 1)
		      jhi = min(jpos+jpsfspan/2, nyorig(ib))
                      if (GotTarget .and. WrtFSimg) then                     ! JWF B30824
                        naxes(1) = ihi - ilo + 1                             ! JWF B30824
                        naxes(2) = jhi - jlo + 1                             ! JWF B30824
                        print *,'jlo,jhi,ilo,ihi:',jlo,jhi,ilo,ihi           ! JWF B30829
                        print *,'ipsfspan, jpsfspan:',ipsfspan,jpsfspan      ! JWF B30829
                        print *,'npsf(ib),xprat,yprat:',npsf(ib),xprat,yprat ! JWF B30829
                        print *,'FSimg NAXIS1&2 =',naxes                     ! JWF B30824
                        nelements=naxes(1)*naxes(2)                          ! JWF B30824
                        allocate(pstamp1(nelements))                         ! JWF B30824
                        allocate(pstamp2(nelements))                         ! JWF B30824
                        pstamp1 = 0.0                                        ! JWF B30824
                        pstamp2 = 0.0                                        ! JWF B30824
                      end if                                                 ! JWF B30824
		      do j = jlo,jhi
		      do i = ilo,ihi
                        if (GotTarget .and. WrtFSimg) kpix = kpix + 1        ! JWF B30824
			itrans = i - xTRANS(ifs,ib)
			jtrans = j - yTRANS(ifs,ib)
			if (itrans >= 1 .and. itrans <= nx(ib) .and.
     *			    jtrans >= 1 .and. jtrans <= ny(ib)) then
			  xpsample = pxcent(ib) + (i-x8)*xprat
			  ypsample = pycent(ib) + (j-y8)*yprat
			  pval =xprat*yprat*bilint(psfs,npsize,npsize,
     *			    npnmax,4,npsf(ib),npsf(ib),
     *			    kpn(iframe),ib,xpsample,ypsample)
                          if (GotTarget .and. WrtFSimg) then                 ! JWF B30824
                            pstamp1(kpix) = Array(itrans,jtrans,ifs_pix)     ! JWF B30824
                            npix = npix + 1                                  ! JWF B30824
                            if (pstamp1(kpix) .le. -900) pstamp1(kpix) = 0.0 ! JWF B30829
                          end if                                             ! JWF B30824
                          Array(itrans,jtrans,ifs_pix) =
     *			    Array(itrans,jtrans,ifs_pix) - flux_nonPM(ib,nn)*pval
                          if (GotTarget .and. WrtFSimg) then                 ! JWF B30824
                            pstamp2(kpix) = Array(itrans,jtrans,ifs_pix)     ! JWF B30824
                            if (pstamp2(kpix) .le. -900) pstamp2(kpix) = 0.0 ! JWF B30829
                          end if                                             ! JWF B30829
			  ipp = nint(xpsample)
			  jpp = nint(ypsample)
			  rfitsq = (fitradii(ib)/pscale(ib))**2
			  if ((i-x8)**2+(j-y8)**2 <= rfitsq .and.
     *			    ipp>0 .and. ipp<=npsf(ib) .and.
     *			    jpp>0 .and. jpp<=npsf(ib)) then
			      sigpval = xprat*yprat*psfuncs(ipp,jpp,k,ib)
			      Unc(itrans,jtrans,ifs_pix) =
     *				sqrt(Unc(itrans,jtrans,ifs_pix)**2
     *				+ (flux_nonPM(ib,nn)*sigpval)**2)
			  endif
			endif
		      enddo
		      enddo
c============================== code added by JWF B30824 ===================================
                      if (GotTarget .and. WrtFSimg .and. (npix .gt. 0)) then
                        status=0
                        call ftgiou(outunit,status)
                        kpix = 0
1100                    kpix = kpix + 1
                        if (kpix .gt. 999) then
                          print *,'ERROR (wpro_v6): cannot construct unique FITS name'
                          print *,'                 turning off WrtFSimg'
                          WrtFSimg = .false.
                          go to 1200
                        end if
                        write(NumStr,*) kpix
                        call trimblnk(NumStr)
                        FSimgNam = FSimgNam0(1:lnblnk(FSimgNam0))//'_a_'
     +                           //NumStr(1:lnblnk(NumStr))//'.fits'
                        print *,'FS image name (before flux subtraction): ',FSimgNam
                        call Access(FSimgNam(1:LNBlnk(FSimgNam)),zexist,.false.)
                        if (zexist) go to 1100
                        status=0
                        call ftinit(outunit,FSimgNam,blocksize,status)
                        status=0
                        call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)
                        status=0
                        call ftppre(outunit,group,fpixel,nelements,pstamp1,status)
                        status=0
                        call ftclos(outunit, status)
                        deallocate(pstamp1)
c
                        FSimgNam = FSimgNam0(1:lnblnk(FSimgNam0))//'_b_'
     +                           //NumStr(1:lnblnk(NumStr))//'.fits'
                        print *,'FS image name (after flux subtraction): ',FSimgNam
                        status=0
                        call ftinit(outunit,FSimgNam,blocksize,status)
                        status=0
                        call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)
                        status=0
                        call ftppre(outunit,group,fpixel,nelements,pstamp2,status)
                        status=0
                        call ftclos(outunit, status)
                        status=0
                        call ftfiou(outunit, status)
                        deallocate(pstamp2)
                      end if
                      if (GotTarget .and. WrtFSimg .and. (npix .eq. 0)) then
                        print *,'No flux subtraction done after all'
                        deallocate(pstamp1)
                        deallocate(pstamp2)
                      end if
1200                  continue
c=========================== end of code added by JWF B30824 ===============================
		    endif ! (flux_nonPM(ib,nn)/sigflux_nonPM(ib,nn) > 2. &c.)
		  endif ! (onframe(n,ifs,ib) .and. wflag(ifs,ib) == 1)
		enddo ! ib = 1,maxbands
		enddo ! ifs = 1,nf
               endif ! ((nn==icomp(1)) .or. (DoThisOne(nn)))
	      endif ! (goodcomp(nn) .or. (nn .gt. nblend))
	     enddo	! end loop over components
            else
              nTossedItAll = nTossedItAll + 1                                   ! JWF B30402
              if (WarnAlloc8err) then                                           ! JWF B30402
                print *,'WARNING: we threw away everything for mdet source',n   ! JWF B30402
                print *,'         non-PM values for current source:'            ! JWF B30402
                print *,'nblend, goodcomp:',nblend,(goodcomp(nn),nn = 1, nblend)! JWF dbg
                print *,'icomp:',(icomp(nn),nn = 1, nblend)                     ! JWF dbg
                print *,'flux:',(flux(1:4,nn), nn = 1, nblend)                  ! JWF dbg
                print *,'p:',p(1:2*nblend)                                      ! JWF dbg
                print *,'racan:', racan(1:nblend)                               ! JWF dbg
                print *,'deccan:',deccan(1:nblend)                              ! JWF dbg
              end if                                                            ! JWF B30402
	    endif ! ((ncp /= 0 .or. nactiveg /= 0) .and. nbands /= 0)
            deallocate(p)
	    deallocate(pmost)
            if (m .ne. mf) then
              print *,'ERROR: sync lost between nonPM and PM output!'
              print *,'       n =',n,'; m =',m,'; mf =',mf,'; resetting m to',mf
              print *,'       mdet RA, Dec:',racan(n),deccan(n)        ! JWF B30329
              m = mf
            end if
	  endif   ! if (any(onframe(n,1:nf,1:maxbands)))
2000      continue                                   ! jump here if bailing out on primary
          if (GotAllTargets .and. QuitEarly) go to 2001  ! JWF dbg
	enddo					! <<<<<<< End source loop
2001    continue  ! JWF dbg
	nsrc_out = m
        Nwpros   = Nwpros + mf

        if(nsentinel .gt. 0) print *,'===WARNING: Encountered ',nsentinel,' sentinel pixels.'

! Close coordinate transformation routine.
	do ifs = 1,nf
        do ib = 1,maxbands
	    iwcs =  WCS(ifs,ib)
	    if (wflag(ifs,ib) == 1) then
	      call wcsclose(iwcs)
	    endif
        enddo
        enddo

! Replace the input RAlist, DEClist, Xpos, Ypos, LBACK, and LSIG with the
! new values.
	RAlist = RAlistnew
	DEClist = DEClistnew
        IDlist = IDlistnew
	Xpos = Xposnew
	Ypos = Yposnew
	LBACK = LBACKnew
	LSIG = LSIGnew
	Rsat = Rsatnew

	deallocate(RAlistnew)
	deallocate(DEClistnew)
	deallocate(Xposnew)
	deallocate(Yposnew)
	deallocate(LBACKnew)
	deallocate(LSIGnew)
	deallocate(Rsatnew)
	deallocate(onframe)
	deallocate(npassprev)
	deallocate(x8set)
	deallocate(y8set)
	deallocate(offset)
	return

	end


	integer(4) function ntrue(x,n)

! Count the number of true values in x.

	integer(4), allocatable :: markers(:)
	logical(4) x(*)
	integer(4) n

	allocate (markers(n))
	markers = 1

	ntrue = sum(markers, mask=x(1:n))

	deallocate(markers)
	return

	end


	subroutine checkpix(i,j,iframe,ib,nn,nvalues,ivalue,jvalue,ifvalue,
     *	    maxpix,maxbands,newpix)

! Check to see if this pixel is already in the list.

        implicit real*4 (a-h,o-z)
	implicit integer (i-n)
	integer nvalues(maxbands),ivalue(maxpix,maxbands),
     *	    jvalue(maxpix,maxbands),ifvalue(maxpix,maxbands)
	logical(4) newpix

	newpix = .true.

	if (nn /= 1) then
	    do nv = 1,nvalues(ib)
		if (i==ivalue(nv,ib) .and. j==jvalue(nv,ib) .and.
     *		    iframe==ifvalue(nv,ib)) then
		  newpix = .false.
		  return
		endif
	    enddo
	endif

	return

	end



	subroutine ppack(pmost,nblend,goodcomp,p)

! Rearrange the parameter vector to exclude bad components.

	implicit real*4 (a-h,o-z)
	implicit integer (i-n)
	real(4) pmost(*),p(*)
	logical(4) goodcomp(*)

! Unpack the parameter vector.
c	pmr = p(2*nblend+1)      ! JWF B30208
c	pmd = p(2*nblend+2)      ! JWF B30208
	mblend = ntrue(goodcomp,nblend)
	nc = 0

	do n = 1,nblend
	    if (goodcomp(n)) then
		nc = nc + 1
		p(nc) = pmost(n)
		p(mblend+nc) = pmost(nblend+n)
	    endif
	enddo

c	p(2*mblend+1) = pmr      ! JWF B30208
c	p(2*mblend+2) = pmd      ! JWF B30208
	return

	end
c
c======================== following code added B30824 by JWF ===================
c
      subroutine trimblnk(string)
      character*(*) string
      integer*4 lnblnk
c
      if (lnblnk(string) .eq. 0) return
c
10    if (string(1:1) .ne. ' ') return
      string = string(2:lnblnk(string))
      go to 10
c
      end
