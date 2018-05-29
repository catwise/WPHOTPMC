c-------------------------------------------------------------------------------------
c Wise Photometry System
c   T. Jarrett & K. Marsh
c  
c Version 1.0: Feb 23, 2008
c Version 2.0; April, 2008   modified to handle multi-frame
c  May 29.2008
c  June 20.2008, looking for memory leaks
c  Aug 24, 2008;  updated WPRO functionality
c  July 2009 -- significant alterations to improve memory usage
c  PROPER MOTION version:  2012 Feb 28, K. Marsh
c  PROPER MOTION version:  2013 Mar 06, J. Fowler
c  UnWISE (coadds) Processing, 1st version:  2017 March 29, THJ 
c  UnWISE (coadds) Processing, 2nd version:  2017 May 03 ; THJ
c  STD and PSF scaling ; 2018, Jan/Feb  THJ
c  W1 cryo scaling for PSFunc and UNC, JWF B80316
c vsn 4.4  B80408: added STDfloor &c.
c vsn 4.4  B80410: added file names to std uncertainty proc status msgs
c vsn 4.4  B80411: ID unc frames by center pixel values
c vsn 4.4  B80412: ID unc frames by pixel sums
c vsn 4.4  B80504: proc. last std frames; added W2 pre-hibernation scale
c vsn 4.4  B80509: same as B80504 but compiled with Tom's gfortran options
c vsn 4.5  B80524: call exit(64) if an img frame is missing
c vsn 4.5  B80529: added more call exit(64) statements for more errors
c 
c-------------------------------------------------------------------------------------

      program WPHot

      implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      character*150000, allocatable :: HDR0(:,:), HDR(:,:)
      character*150000 Hfits, HeadNull
      character*3000 string

      character*200 fname,s0,namlis,finf,arg(99),fout
      character*200 ofile,ifile,command2,hfile,metaf
      character*200 outfile , Hnull
      character*200 psfdir,calbname,qdir
      character*200 Pnames (100,4), Punc(100,4)
      character*200 hname,mdexname,coname, mdex, IOCname
      character*200 coaddfits(4), icoaddfits(4),svbfits (4)
      character*200 COUfits(4),cocovfits(4), comskfits(4), xscf, xscl, pspit

      character*200, allocatable :: msk(:,:),uncn(:,:),fram(:,:)
      character*200, allocatable :: frames(:), path(:), basename(:)

      character*200 comment, mipspsf, mipsunc, gztmp
      character*50 command,ffram, ftab, imID, DATE

      character*10 intword,uncword,covword,mdexword,level,clevel
c      character*10 version ! Now handled by 'vsn' below
      character*15 What
      character*1 cc,TYPE

c Main MEP table variables
      character*200 meproot                                       ! CJG B30228
      character*9, allocatable :: namtable(:,:)                   ! CJG B30212
      real*8, allocatable :: jdtable(:,:)                         ! CJG B30212
      real*4, allocatable :: meptable (:,:,:)                     ! CJG B30228
      character*9, allocatable :: namtable_tmp(:,:)               ! TPC
      real*8, allocatable :: jdtable_tmp(:,:)                     ! TPC
      real*4, allocatable :: meptable_tmp (:,:,:)                 ! TPC
        integer nfmep, nfmep_old, ngood, nmep, nxy                  ! TPC

      real*8, allocatable :: RAlist_0(:), DEClist_0(:)
c      real*8, allocatable :: RAlist (:), DEClist (:)              ! JWF B21221
      real*8, allocatable :: RAlist (:), DEClist (:), deltaJD(:)  ! JWF B21221
      real*8, allocatable :: RAsublist (:), DECsublist (:)
      real*8, allocatable :: tmpRAlist (:), tmpDEClist (:), wproJD(:), varJD(:), JD(:,:)
        real*4, allocatable :: Rsat(:,:),Rsatsub(:,:)
      real*4, allocatable :: Table_0 (:,:),Table (:,:),XSCprox (:)
      real*4, allocatable :: Tarray (:,:,:), XSCshape(:,:)
      real*4, allocatable :: Larray (:), SNRlist (:,:,:)
        integer*4, allocatable :: iLarray(:), iTarray(:,:)
      real*4, allocatable :: Array (:,:,:), Unc (:,:,:)
        integer, allocatable :: Mask (:,:,:), iMask(:,:,:)
      real*4, allocatable :: Xpos (:,:,:), Ypos (:,:,:)
      real*4, allocatable :: LBACK (:,:,:), LSIG (:,:,:), Lconf(:,:,:)
      real*4, allocatable :: GlobalStats (:,:)
c      real*4, allocatable :: AppCorr (:,:)
      real*4, allocatable :: tmpMAG (:,:), etmpMAG (:,:), tmpFLUX(:,:), etmpFLUX(:,:)
      real*4, allocatable :: varFLUX(:),eVarFLUX(:), VarRchi2(:)
      real*4, allocatable :: ztmp (:), eztmp (:), peztmp(:)
      real*4, allocatable :: PSFs (:,:,:,:), PSFuncs(:,:,:,:)
      real*4, allocatable :: MAGSTD (:,:), eMAGSTD (:,:)
      real*4, allocatable :: framefluxes(:,:,:),frameuncs(:,:,:)
      real*4, allocatable :: framechis(:,:,:)

      integer*2, allocatable :: xTRANS (:,:), yTRANS (:,:), SatNum(:,:) , SatNumSub(:,:) , nDF (:,:)
        ! this set of indices determine the order of the frame file read
      integer*2, allocatable :: tnew(:),torder (:),order (:), preorder (:,:), fractive(:,:)
        ! New index to map from order() array to per-band offset into the pixel arrays - TPC
      integer*4, allocatable :: pix_order(:,:)
        integer*4                 npixfr_band(4), nactive_band(4)

      real*4  zero(4),z1,z2,z3,z4
        real*4  Rcirc (19,4), Rapc(19), sky, sdev
      real*4  mfluxtmp(19), mefluxtmp(19), mzmag(19), mzerr(19)
      real*4  Rcoadd(4),cmag(19,4), cemag(19,4)
      real*4  cozero(5), IRACapcor(5), cov_ave(5), cov_std(5)
      real*4  STDscale (4)    !  thj 22Jan2018
      real*4  PSFuncscale (4)    !  thj 22Feb2018
      real*4  W1cryoPSFuncscale, PrehibSTDScale(4)  !  JWF B80316
      real*4  STDbias(4)                            !  JWF B80402

      real*4, allocatable :: COADD(:,:,:), COUNC(:,:,:),SVB(:,:,:), COCOV (:,:,:)
      integer*2, allocatable :: COmsk(:,:,:)
      real*4, allocatable :: WPROmagR (:,:,:), medsky(:,:), medskysig(:,:), formalsig(:,:)
      real*4, allocatable ::  slope (:,:), eslope (:,:), slopechi2 (:,:)

      real*4  SVBvalcenter(4)
      real*4  crval1,crval2,
     1     cdelt1,cdelt2,crot,crpix1,crpix2,scrot(4)
c      real*4  scrval1(4),scrval2(4),
c     1     scdelt1(4),scdelt2(4),scrot(4),scrpix1(4),scrpix2(4)

      real*4  cscale (4), pscale (4), Rstap (4), RstapF (4), Rstann (4), Rstwid (4)
      real*4  metval, adb_alloscale
      real*4 second(2),second2(2)
      real*4  pcdelt1,pcdelt2,pcrota,pcrpix1,pcrpix2,zz,pcrval1,pcrval2,pcrot
      real*4 pxcent(4), pycent(4),ppix(4),fwhm(4)
      real*4 SingleArray (19,4), eSingleArray (19,4)
      real*4 BookeepApCorr (4)

      real*4  fluxcon(5), funfactor(5),  Fcorr (4), BGmax (4)
      real*4  IRACphotcorr (256,256,4)
      real*4 CmosaicCorr (4), galmag(4), galerr(4), galba (4), galpa (4), galRsemi (4)

      real*4  CoaddXY (2,4), ChiFac   ! chi2 threshold for chkvar
      real*4 dks                                            ! CJG B30320

      real*8 ra0, dec0, xpix,ypix, tJD, sumJD, MJDhibern8, MJDendCryoPSF

      integer NmerSTD(19,4), SingleFlg (19,4), coFLAG(19,4), NbannAve(4)
      integer N_M_wpro, M_M_wpro, M_M_var
      integer nx,ny,nsx(4),nsy(4),lsize,  WCS0, offscl, cband0, WCSc (4)
      integer ncx(4), ncy(4)
      integer nsubx(4),nsuby(4), i0_reg(4), j0_reg(4)
      integer calgridX,calgridY,nsrc, ntmp (4)

      integer*4 mfiles                                       ! CJG B30226
      integer*4 kStat                                        ! JWF B30404
      real*4    SatMag(2), SatFlux(2)                        ! JWF B30614
      real*8    SatMJD                                       ! JWF B30614
      logical*4 NotSat                                       ! JWF B30614
      integer*4 Access                                       ! JWF B60711
      real*8    PixSum                                       ! JWF B80412
      integer   doreg, nbuf_band(4) ! TPC
      real*4    buf_as

        integer*4, allocatable ::  Rho12(:), Rho23(:), Rho34(:), P12(:), P23(:), P34(:)
c      real*4, allocatable ::  Fmin(:,:), Fmax(:,:)
c      real*4, allocatable ::  DeltaMag (:,:)
        real*4, allocatable ::  mLogQ(:,:)
      real*4, allocatable ::  ks (:,:)                       ! CJG B30319
      real*8, allocatable ::  AVE_mJD(:,:),mJD_min(:,:),mJD_max(:,:)
      real*8, allocatable ::  AVE_mJD_var(:,:),mJD_min_var(:,:),mJD_max_var(:,:)

      integer*2, allocatable :: wflag(:,:), wflagactive(:,:), Mdupe(:,:)
      integer*2, allocatable :: MapPSFs(:,:,:)
      integer*2, allocatable :: N_Mstat(:,:), M_Mstat(:,:)
      integer, allocatable :: tmpFLG(:,:), WCS(:), FLGSTD (:,:), IDlist(:)
        integer, allocatable :: mdetIDlist(:), tmpIDlist(:), IDlist_0(:)
      integer, allocatable :: iCOADD (:,:,:)
      integer*2 Kdex, fmode
      integer nsr , nmax,  nfi, nf, nwpro,SingleFlag, mergetype
      integer npsf(4), npn(4), adb_nmax
      integer miflag(19), fbits (32),  galflag (4), fatal
      integer region_order(225), region_check(225) ! TPC

      real*4 mmLogQ

      logical debug,smode,verbose,unitest,SPIT, MIPS, NEP, SEP,erase,zexist, cexist (4), sexist (4), docentroid
      logical doweight, doSVB, pointless, IzBad, Circinus, doioc, gotscale (4), dosolo
      logical mepwrite                              ! CJG B30228
      logical*4 NamWrt           ! JWF B30304
      logical*4 sort_subframes ! TPC
      data      NamWrt,sort_subframes/.true.,.false./   ! JWF B30304
      logical STDfloor(4)               ! JWF B80408
      data    STDfloor/4*.false./       ! JWF B80408
      data    PixSum/0.0d0/             ! JWF B80412

      data fluxcon/0.1088,0.1388,0.5952,0.2021,0.0454/   ! Spitzer conversion from MJy/st to Dn/s

c      data Fcorr /31., 34., 40., 115./ ! WISE, old jarrett values based on WAPP repeatability
c      data Fcorr /32., 37., 58., 221./ ! WISE, values based on FMasci theoretical values

c       data Fcorr /11., 10., 15., 13./   ! correlation factor for IRAC (old WAPP-calibrated values)
c       data Fcorr/11.4, 11.9, 17.3, 25.66/  ! correlation factor for IRAC based on Frank's calculations

c      data Fcorr /11., 10., 13., 8./   ! correlation factor for IRAC,,MIPS-24
c       data Fcorr /11., 13., 11., 20./   ! correlation factor for IRAC (DeepCarol) OLD
c       data Fcorr /22.1, 23.3, 33.5, 49.5/ ! correlation factor for IRAC (DeepCarol)
c
c
c      data mosaicCorr/0.183, 0.231, 0.583, 0.371/   !  mag aper corrections for WISE mosaics
c
c       real*4 fDup                                   ! JWF B21206
c
        include 'jwfcom.f'                            ! JWF B21219
        data      Nok2a,Nok2b/2*0/, nPMsolns/0/       ! JWF B30227
        data      maxsteps/250/, nBlendSwap/0/        ! JWF B30221
        data      PMsubtract/.false./, Nwpros/0/      ! JWF B30126
        data      PMinitADB/.false./                  ! JWF B30207
        data      MeanObsEpochType/2/, nWgtFail/0/    ! JWF B30208
        data      TossZeroFlux/.false./               ! JWF B30221
        data      GenJD0posns/.false./                ! JWF B30221
        data      InitPMflux0,InitPMposn0/2*.false./  ! JWF B30312
        data      ftol_npm,ftol_pm/2*1.0e-3/          ! JWF B30312
        data      nSingMatFrmFunc/0/, RblendFac/1.0/  ! JWF B30402
        data      WarnAlloc8err/.false./              ! JWF B30402
        data      nTossedItAll/0/, BlendType/2/       ! JWF B30404
        data      mLogQfac/1.0/, tdbg/.false./        ! JWF B30419
        data      BGtype/1,3,1,3/, BGtrimFrac/0.20/   ! JWF B30604     !  TJ fixed the W3 slot to mirror W1 ; 14Dec2017
        data      SkThresh/9.9e25/                    ! JWF B30604
        data      SrcSubSNR/2.,2.,2.,2./              ! TPC B30605
        data      WarnBigFrameSig/9.9e25/             ! JWF B30711
        data      SingleFrame,postcryo/2*.false./     ! JWF B31121
        data      STDscale,PSFuncscale/8*1.0/         ! JWF B80226
        data      W1cryoPSFuncscale/1.0/              ! JWF B80316
        data      PrehibSTDScale/4*1.0/               ! JWF B80316
        data      STDbias/4*0.0/                      ! JWF B80402
c
      namelist/WPHpars/nx,ny,nxAWAIC,nyAWAIC,
     1     zero, cozero, CmosaicCorr, fbits,
     1     Naper, Rapc, IRACapcor, Rannulus, Rwidth,
     1     edgebuf, fwhm, Fcorr, BGmax, ChiFac,
     1     adb_nmax, adb_alloscale,
     1     ireg, jreg, nbuf,
     1     STDscale,PSFuncscale, STDbias, STDfloor,       ! THJ  22Jan/Feb2018, JWF B80316-B80408: these new parameters scale 
     +     W1cryoPSFuncscale,PrehibSTDScale,              !      the STD (unc) and PSFunc images, respectively
     +     MJDendCryoPSF, MJDhibern8,                     ! JWF B80322
     1     doSVB, hname,  Hnull, MJD0, PMfac, minPMsig,   ! JWF B21207
     +     maxsteps, GenJD0posns, TossZeroFlux,           ! JWF B30117
     +     PMsubtract, PMinitADB, MeanObsEpochType,       ! JWF B30208
     +     TargetRA, TargetDec, QuitEarly, WarnZeroFlux,  ! JWF B30227
     +     NamWrt, poscon0, InitPMflux0, InitPMposn0,     ! JWF B30312
     +     ftol_npm, ftol_pm, UseNonPMposns,              ! JWF B30312
     +     max_mep_lines, WarnRunaway,WarnAlloc8err,      ! CJG B30219
     +     RblendFac, BlendType, mLogQfac, tdbg,          ! JWF B30404
     +     MaxWarnPMNaN, pminFac,                         ! JWF B30507
     +     region_order,sort_subframes,                   ! TPC
     +     MaxWarnBlendSwap, CorrPMerr,                   ! JWF B30520/B30521
     +     PMstepFac, DumPvec, nbmax, BGtype, BGtrimFrac, ! JWF B30524/B30529/B30604
     +     SkThresh, SrcSubSNR, SatMag, SatMJD,           ! JWF/TPC B30605
     +     WarnBigFrameSig, WrtFSimg, SingleFrame,        ! JWF B30711/B30824/B31203
     +     postcryo                                       ! JWF B31203

      data max_mep_lines/999999/                        ! CJG B30228
c
      data MJD0/56700.d0/, nBadBlend/0/, nPMNaN/0/,     ! JWF B30130 ;  THJ 24Feb2018
     +     PMfac/1.0/, minPMsig/0.01/, nAllZero/3*0/,   ! JWF B30130,B80311
     +     TargetRA,TargetDec/20*9.9d9/,                ! JWF B30227
     +     QuitEarly/.false./, WarnZeroFlux/.true./,    ! JWF B30227
     +     poscon0/.true./, UseNonPMposns/.true./       ! JWF B30308
     +     WarnRunaway/.false./, ChiFac/3.0/,           ! JWF B30308
     +     MaxWarnPMNaN/0/, pminFac/1.0e-4/,            ! JWF B30507
     +     MaxWarnBlendSwap/0/, CorrPMerr/.false./,     ! JWF B30520/B30521
     +     PMstepFac/100.0/, DumPvec/.false./, nbmax/3/ ! JWF B30524/B30529
     +     nFallbackMMM/0/, SatMag/8.0,7.0/,            ! JWF B30606/B30614
     +     SatMJD/55468.7775/, WrtFSimg/.false./,       ! JWF B30614/B30824
     +     MJDhibern8/55593.96d0/,                      ! JWF B80322
     +     MJDendCryoPSF/55480.0d0/                     ! JWF B80322
c
      integer*4 N_adb, nBorderViolators, nbordrej, nduprej, nnanposrej
      data N_adb/0/, nBorderViolators/0/, nbordrej/0/, nduprej/0/, nnanposrej/0/
c
      character*11 vsn                ! JWF B30507
      character*8  cdate, ctime       ! JWF B21109
      integer*4    jdate(3),jtime(3)  ! JWF B21109
      integer*4    IArgc,nCWchk       ! JWF B80404
      data         vsn/'4.5  B80529'/ ! JWF
      common /vdt/ cdate,ctime,vsn    ! JWF B30507
      logical findpeak                ! JWF B60714
      logical DidCryo                 ! JWF B80307
      data    DidCryo/.false./        ! JWF B80307
      integer*4, allocatable :: nBand(:) ! JWF B80318
      real*8,    allocatable ::   MJD(:) ! JWF B80318

c
      call signon('WPHotpmc')          ! JWF B30507
      if (IArgc() .eq. 0) then
        print *,'For instructions on running this program: ask Tom Jarrett'
        call signoff('WPHotpmc')
        stop
      end if
c
c     print *, 'Revision = $Id: WPHotpm.f 13022 2013-11-09 00:13:05Z tim $'
c
      iverify = 0  !  this is flipped to one by writeTABLE when an output file is written

      dd9 = 0.
        call  ETime (Second2, dd9)
        dd9 = second2(1)  ! User time


      doSVB = .false.
      SPIT = .false.
      NEP = .false.
      SEP = .false.
      MIPS = .false.
      circinus = .false.

      gotscale = .false.

      docentroid = .false.
c      docentroid = .true.

      doweight = .true.
c        doweight = .false.

      unitest = .false.
       debug = .false.
c      debug = .true.
      verbose = .false.
      mepwrite = .false.                                          ! CJG B30213
      mfiles = 0

      istat = 0

c      write (6,'(a,a)') 'WPHOT version ',version
c      write (6,*) ' '


      if (SPIT) then
            write (6,*) '*********** SPITZER MODE Invoked ********************************************'
            if (MIPS) then
                   write (6,*) '   &&&&&&&&&  MIPS mode invoked  &&&&&&&&&'
c                   mipspsf = 'mips_24_5000c_256x256.fits'
c                   mipsunc = 'mips_unc_24_5000c_256x256.fits'
                   mipspsf = 'mips_24_256x256.fits'
                  mipsunc = 'mips_unc_24_256x256.fits'

                  if (NEP) then
c                         mipspsf = '/wise/jarrett/cal/psf/IRAC/coadd/nep_MIPS/IRACPSF-w5-psf-wpro-01x01-01x01.fits'
c                          mipsunc = '/wise/jarrett/cal/psf/IRAC/coadd/nep_MIPS/IRACPSF-w5-psfunc-wpro-01x01-01x01.fits'
                          mipspsf = '/wise/jarrett/NEP/psf/MIPS_NEP_PSF.fits'
                          mipsunc = '/wise/jarrett/NEP/psf/MIPS_NEP_PSFunc.fits'
                        else
                          mipspsf = '/wise/jarrett/cal/psf/IRAC/coadd/sep_MIPS/IRACPSF-w5-psf-wpro-01x01-01x01.fits'
                          mipsunc = '/wise/jarrett/cal/psf/IRAC/coadd/sep_MIPS/IRACPSF-w5-psfunc-wpro-01x01-01x01.fits'
                        endif

                        if (Circinus) then
                            mipspsf = '/wise/jarrett/circinus/combine/psf/MIPS24-Circinus-psf-wpro-01x01-01x01.fits'
                            mipsunc = '/wise/jarrett/circinus/combine/psf/MIPS24-Circinus-psfunc-wpro-01x01-01x01.fits'
                        endif

            endif
      else
            MIPS =.false.
      endif


c read the command arguments
      
      narg = 0
      do j=1,99
        L = 0
        call getarg (j,s0)
        L = numchar (s0)
        if (L.le.0) then
            ! we are done reading
            goto 111
        endif
        narg=narg+1
          if (narg .gt. 99) go to 112
        arg(narg) = s0
      enddo

 111      if (verbose) write (6,*) narg,'  arguments loaded'

c      if (narg.ge.99) then! JWF B21223
112     if (narg.gt.99) then
            write (6,*) 'ERROR -- max arguments reached; exiting ...'
            call exit(9)
      endif


      call parg (narg,arg,level,ofile,intword,uncword,covword,mdexword,
     1    ifile,namlis,metaf,smode,verbose, pointless, unitest,
     1     qdir,psfdir,calbname,calgridX,calgridY,
     1     mdexname,coname, clevel, imID, xscf, xscl, IOCname, detsnr,
     1     meproot,doreg)                                             ! CJG B30228, TPC

      LM = numchar(meproot)
      if(LM .gt. 1) mepwrite = .true.      ! CJG B30213


c      zexist = .false.                                   ! JWF B60711
c      inquire (directory = qdir,exist=zexist)            ! JWF B60711
      zexist = Access(qdir(1:LNBlnk(qdir)),' ') .eq. 0   ! JWF B60711

      if (.not.zexist) then
            write (6,*) 'ERROR -- QA directory does not exist; exiting ...'
                call exit(9)
        endif

c      zexist = .false.                                      ! JWF B60711
c      inquire (directory = psfdir,exist=zexist)             ! JWF B60711
      zexist = Access(psfdir(1:LNBlnk(psfdir)),' ') .eq. 0  ! JWF B60711

        if (.not.zexist) then
                write (6,*) 'ERROR -- PSF directory does not exist; exiting ...'
                call exit(9)
        endif


      if (unitest) write (6,*) 'Executing UNIT-TEST mode'


      if (calgridX.gt.0) then

       if (debug) then
       write (6,'(a)') psfdir
       write (6,'(a)') calbname
       write (6,*) calgridX,calgridY
       endif

      else

            write (6,*) '** Error:  calgridX <= 0'
            call exit(5)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c namelist

      BGmax = 0.
      do j=1,32
        fbits(j) = -1
      enddo
      do j=1,225
        region_order(j) = j
        region_check(j) = 0
      enddo

        open (unit=10,file=namlis)
        read (10,WPHpars)
        close (10)
        if (NamWrt) write(6, WPHpars)  ! JWF B30304

        if (PMinitADB) print *,                                      ! JWF B30305
     +     'WARNING: PMinitADB selected; this is totally untested!'  ! JWF B30305

      if(doreg .lt. 0 .or. doreg .gt. ireg*jreg) then
        write (6,*) '***ERROR -- Bad -doreg parameter setting; ',
     1                'doreg,ireg,jreg=',doreg,ireg,jreg,'. ',
     1                'Exiting ...'
        call exit(64)
      endif

        if(ireg.lt.1 .or. ireg.gt.15 .or. jreg.lt.1 .or. jreg.gt.15 .or. ireg.ne.jreg) then
          print *,'***ERROR: Ireg or jreg out of range (1-15) or not equal:',ireg,jreg
          call exit(64)
        endif
        nreg = 0
        do j=1,ireg*jreg
        kreg = region_order(j)
          if(kreg .le. 0) exit
          if(kreg .gt. ireg*jreg) then
            print *,'***ERROR: Requested region ',kreg,' is out of range.'
            call exit(64)
          endif
        region_check(kreg) = region_check(kreg) + 1
          if(region_check(kreg) .gt. 1) then
            print *,'***ERROR: Region ',kreg,' specified more than once.'
            call exit(64)
          endif
          nreg = nreg + 1
        enddo
        if(doreg .eq. 0) print *, 'Processing', nreg, 'out of', ireg*jreg, 'regions:', (region_order(i),i=1,nreg)
        if(doreg .gt. 0) print *, 'Processing ONLY region', doreg
c       if(nreg .lt. ireg*jreg .or. doreg) print *,'!!! Not all regions being processed.'   ! JWF B60711
        if((nreg .lt. ireg*jreg) .or. (doreg .ne. 0)) print *,'!!! Not all regions being processed.' ! JWF B60711

        if ((MeanObsEpochType .lt. 0) .or. (MeanObsEpochType .gt. 5))
     +  then
          print *,'WARNING: illegal MeanObsEpochType:', MeanObsEpochType
          MeanObsEpochType = 2
          print *,'         reset to 2'
        end if

        if ((BlendType .lt. 1) .or. (BlendType .gt. 2)) then
          print *,'***ERROR: illegal BlendType specification:', BlendType
          call exit(64)
        end if

c
        if (GenJD0posns) then
          open(75, file = 'JD0_posns_'//ofile)
          write(75,'(a)') '|   ra_npm  |   dec_npm |   ra_JD0  |'
     +       //'   dec_JD0 |sigra_JD0|sigdec_JD0|sigradec_JD0|'
     +       //' dist |  PMRA | PMDec |'
          write(75,'(a)') '|     r     |      r    |     r     |'
     +       //'      r    |     r   |     r    |      r     |'
     +       //'   r  |   r   |   r   |'
          write(75,'(a)') '|    deg    |     deg   |    deg    |'
     +       //'     deg   |   asec  |   asec   |    asec    |'
     +       //' asec |asecpyr|asecpyr|'
          write(75,'(a)') '|    null   |    null   |    null   |'
     +       //'    null   |   null  |   null   |    null    |'
     +       //' null |  null |  null |'
        end if
        if (DumPvec) then
          open(76, file = 'Pvec_'//ofile)
          write (76,'(a)') '|   p_x   |   p_y   |  p_rad  |     ra    |    dec    |nIters|nSteps|nb0|'
          write (76,'(a)') '|    r    |    r    |    r    |     r     |     r     |   i  |   i  | i |'
          write (76,'(a)') '|   asec  |   asec  |   asec  |    deg    |    deg    |   -  |   -  | - |'
          write (76,'(a)') '|   null  |   null  |   null  |    null   |   null    | null | null |nul|'
          open(77, file = 'Pvec_pm_'//ofile)
          write (77,'(a)')
     +   '|   p_x   |   p_y   |  p_rad  |   pm_x  |   pm_y  |  pm_rad |'
     +              //'     ra    |    dec    |nIters_pm|nSteps_pm|nb0|'
          write (77,'(a)')
     +   '|    r    |    r    |    r    |    r    |    r    |    r    |'
     +              //'     r     |     r     |    i    |    i    | i |'
          write (77,'(a)')
     +   '|   asec  |   asec  |   asec  | asecpyr | asecpyr | asecpyr |'
     +              //'    deg    |    deg    |    -    |    -    | - |'
          write (77,'(a)')
     +   '|   null  |   null  |   null  |   null  |   null  |   null  |'
     +              //'    null   |   null    |   null  |   null  |nul|'
        end if

        SatFlux(1) = 10.0**(0.4*(zero(1)-SatMag(1)))      ! JWF B30614
        SatFlux(2) = 10.0**(0.4*(zero(2)-SatMag(2)))      ! JWF B30614
c       print *,'SatFlux:', SatFlux                       ! JWF B30614

      ncoaddsize = nxAWAIC

      if (SPIT) then  ! mosaic zero points for pscale = 1 arcsec

      else

        do ib=1,4
            cozero(ib) = zero(ib)
         enddo

      endif

        Raper = Rapc (1)   ! standard aperture in arcsec

c       do j=1,Naper
c        write (6,*) Rapc(J)
c       enddo

        if (verbose) write (6,*) 'namelist loaded'


ccc  TJ  07-Sep-2017
	fatal = 0
        do I = 1,32
                ibit = fbits(I)
                if (ibit.ge.0) then
                  fatal = fatal + (2**ibit)
                else
                        !quit
                        goto 741
                endif
        enddo
 741     I = 0


cccccccccccccccccccccccccccccccccccccccccccccccccccc

       call countfile (ifile,nfi)
c       write (6,*) 'frames = ',nfi

        allocate (wflag(nfi,4))
        if (.not.allocated(wflag)) then
          write (6,*) 'allocation failed for wflag'
          istat = 5
          call exit(istat)
        end if


        allocate (path(nfi))   ! number of frames
        if (.not.allocated(path)) then
          write (6,*) 'allocation failed for path'
          istat = 5
          call exit(istat)
        end if

        allocate (basename(nfi))
        if (.not.allocated(basename)) then
          write (6,*) 'allocation failed for basename'
          istat = 5
          call exit(istat)
        end if

cccccccccccccccccccccccccccccccccccccccccccc
c load the input list of frames

      if (verbose) write (6,*) 'load input list of frames '
      call pinput (ifile,nfi,nf,path,basename,wflag,smode)

      if(nf .ne. nfi) then
          print *,'***ERROR: countfile and pinput disagree about the frame count: ',
     1            nfi, ' vs. ', nf
          call exit(64)
        endif        

c       don't need nfi below here anymore - TPC

c      write (6,*) 'testdone',nf  ! TEMP

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c now deal with the PSFS

         allocate (MapPSFs(nx,ny,4))                            ! ADDED
         if (.not.allocated(MapPSFs)) then                      ! ADDED
           write (6,*) 'allocation failed for MapPSFs'          ! ADDED
           istat = 5                                            ! ADDED
         end if                                                 ! ADDED

      if (verbose) write (6,*) 'load PSFs'

         npn = 0
         call  PSFnames (psfdir,calbname,calgridX,calgridY,Pnames,Punc,
     1      npn,MapPSFs,nx,ny,MIPS,mipspsf,mipsunc, nf, wflag)



         nnx=0
         nny=0
         npsize = -1
         do ib=1,4

c          if ( wflag ( 1,ib ) .eq.1 ) then
         if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

           do j=1,npn(ib)                                         ! MODIFIED/MOVED
           fname = Pnames(j,ib)

           call headpar(nnx,nny,fname,pcrval1,
     1         pcrval2,pcdelt1,pcdelt2,pcrota,pcrpix1,pcrpix2,
     1          xlo,xhi,ylo,yhi,istat)

            if (istat /= 0) call exit(istat)
            if (max(nnx,nny) > npsize) npsize = max(nnx,nny)
           enddo

          endif

         enddo

         if (verbose) write (6,*) 'npsize = ',npsize

         npnmax = maxval(npn)

        if (verbose) write (6,*) 'npnmax = ',npnmax,npsize

         allocate (PSFs(npsize,npsize,npnmax,4))
         if (.not.allocated(PSFs)) then
           write (6,*) 'allocation failed for PSFs'
           istat = 5
          end if

          allocate (PSFuncs(npsize,npsize,npnmax,4))
         if (.not.allocated(PSFuncs)) then
           write (6,*) 'allocation failed for PSFuncs'
           istat = 5
         end if


ccccccccccccccccccccccccccccccccccccccccccccccc
      if ((nf.eq.1).or.(smode)) then
            ! do not divide the coadd region
                ireg = 1
                jreg = 1
c      else if ((ireg.eq.1).and.(jreg.eq.1)) then
c            ! force 2x2, seems to run faster
c            ireg = 2
c            jreg = 2
      endif


      iblow = 1
      ibhigh = 4
      ib0 = 1
      doIOC = .false.

      if (smode) then
         write (6,*) '%%Invoking single-band mode'
       write (6,*) '%%band engaged: ',(wflag(1,ib),ib=1,4)

         do ib=1,4
            if (wflag(1,ib).eq.1) then
             ib0 = ib
           iblow = ib
           ibhigh = ib
            endif
         enddo

       LI = numchar(IOCname)
       if (LI.gt.1) then
            doIOC = .true.
       endif
       

        endif

      if (verbose) write (6,*) 'frames = ',nf

      if (nf.lt.1) then
         write (6,*) 'ERROR -- number of frames is too small ; exiting '
         call exit(64)
      endif

      do j=1,nf
        L1 = numchar(path(j))
        L2 = numchar(basename(j))
        if (verbose) write (6,'(i5,1x,a,1x,a,1x,4i3)') j,path(j)(1:L1),basename(j)(1:l2),(wflag(j,ib),ib=1,4)
      enddo


        L2 = numchar(ofile)

      outfile  = ofile(1:L2) 
c        metaf =  ofile(1:L2-4) // '_WPHOT-meta.tbl'

      L = numchar (outfile)
      if (verbose) write (6,'(a,a)') 'output file = ',outfile (1:L)

        imeta = 99

      L = numchar (metaf)
      if (verbose) write (6,'(a,a)') 'output meta  = ',metaf (1:L)

        call initmeta (metaf,imeta, level, ifile, namlis, nf, wflag,
     1    psfdir,calbname,calgridX,calgridY, vsn, STDscale,PSFuncscale,
     +    W1cryoPSFuncscale,PrehibSTDScale,STDbias)   !  THJ 22Jan2018, 22Feb2018, JWF B80316,B80402

c  write to meta
      What = "nf"
      TYPE = "i"
      metval = nf*1.
      comment = 'number of frames'
      iband = 0
        if (smode) iband = ib0

      call MetWr (imeta, iband, What, TYPE, metval, comment)

      do ib=1,4
        What = "npn"
        TYPE = "i"
        metval = npn(ib)*1.
        comment = 'PSFs'
        iband = ib
        if (smode) then
          if (ib.eq.ib0) call MetWr (imeta, iband, What, TYPE, metval, comment)
        else
          call MetWr (imeta, iband, What, TYPE, metval, comment)
        endif
      enddo

      allocate (fram(nf,4))
        if (.not.allocated(fram)) then
          write (6,*) 'allocation failed for fram'
          istat = 5
        end if

      allocate (order(nf))
        if (.not.allocated(order)) then
          write (6,*) 'allocation failed for order'
          istat = 5
        end if
      allocate (preorder(nf,4))
        if (.not.allocated(preorder)) then
          write (6,*) 'allocation failed for preorder'
          istat = 5
        end if
      allocate (fractive(nf,4))
        if (.not.allocated(fractive)) then
          write (6,*) 'allocation failed for fractive'
          istat = 5
        end if
      allocate (torder(nf))
        if (.not.allocated(torder)) then
          write (6,*) 'allocation failed for torder'
          istat = 5
        end if
      allocate (tnew(nf))
        if (.not.allocated(tnew)) then
          write (6,*) 'allocation failed for tnew'
          istat = 5
        end if


      allocate (msk(nf,4))
        if (.not.allocated(msk)) then
          write (6,*) 'allocation failed for msk'
          istat = 5
        call exit(istat)
        end if

      allocate (uncn(nf,4))
        if (.not.allocated(uncn)) then
          write (6,*) 'allocation failed for uncn'
          istat = 5
        call exit(istat)
        end if


      allocate (HDR0(nf,4))
        if (.not.allocated(HDR0)) then
          write (6,*) 'allocation failed for HDR0'
          istat = 5
        call exit(istat)
        end if

c      allocate (mdex(nf))
c        if (.not.allocated(mdex)) then
c          write (6,*) 'allocation failed for mdex'
c        istat = 5
c        call exit(istat)
c        end if


      allocate (WCS(nf))
        if (.not.allocated(WCS)) then
          write (6,*) 'allocation failed for WCS'
          istat = 5
        call exit(istat)
        end if



      Lextra = nx * ny * 4 * 4
      Lsize = Lextra
      if (SPIT) Lextra = nxAWAIC*nyAWAIC
      allocate (Larray(Lextra))
        if (.not.allocated(Larray)) then
          write (6,*) 'allocation failed for Larray'
          istat = 5
        call exit(istat)
        end if
      allocate (iLarray(Lextra))
        if (.not.allocated(iLarray)) then
          write (6,*) 'allocation failed for iLarray'
          istat = 5
        call exit(istat)
        end if

      
      allocate (Mdupe(nxAWAIC,nyAWAIC))
        if (.not.allocated(Mdupe)) then
          write (6,*) 'allocation failed for Mdupe'
          istat = 5
          call exit(istat)
        end if

      do j=1,nyAWAIC
      do i=1,nxAWAIC
        Mdupe(i,j) = 0
      enddo
      enddo



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c construct table names

      nmax = 0
      NSUMmax = 0

      J=1
        L1 = numchar(path(j))
        L2 = numchar(basename(j))
      L3 = numchar(mdexword)
      L4 = numchar (level)

      j=1
      mdex = path(j)(1:L1) // '/' //
     1          basename(j)(1:L2) //
     1          '-' // mdexword (1:L3) //
     1          '-' // level (1:L4) //
     1          '.tbl'

c      if ((nf.gt.1).or.(SPIT)) then  !  multi-frame
            L2 = numchar(mdexname)
            mdex =  mdexname(1:L2)
c      endif


        if (smode) then

         do ib=1,4
            write (cc,'(i1)') ib

            if (wflag(j,ib).eq.1) then
c           mdex = path(j)(1:L1) // '/' //
c     1          basename(j)(1:L2) //
c     1          '-w' // cc // 
c     1          '-' // mdexword (1:L3) //
c     1          '-' // level (1:L4) //
c     1          '.tbl'
          endif            

         enddo

        endif


        L = numchar(mdex)       
        
          if (verbose) write (6,'(a)') mdex(1:L)

        call countfile (mdex,nsr)
        if (nsr.eq.0) then
            write (6,'(a)') mdex(1:L)
            write (6,*) 'caution -- src table file is empty'
c            call exit (9)
        endif

        nmax = max (nmax, nsr )

        NSUMmax = max(NSUMmax, nsr)


      if (verbose) write (6,*) 'nmax sources = ',nmax
      if (verbose) write (6,*) 'total sources = ',NSUMmax



      What = "doSVB"
      TYPE = "i"
      metval = 0.
      if (doSVB) metval = 1.
      comment = 'Apply SVB subtraction from each frame?  0 = false, 1 = true'
      iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)


      What = "adb_nmax"
      TYPE = "i"
        metval = adb_nmax
      comment = 'active deblend throttle: maximum allowable number of deblend components'
      iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "adb_alloscale"
        TYPE = "r"
        metval = adb_alloscale
        comment = 'active deblend allocation scaling factor'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)


c scale the maximum number of sources to allow WPRO to actively deblend

      NSUMmax = nint(NSUMmax * adb_alloscale )

        What = "NSUMmax"
        TYPE = "i"
        metval = NSUMmax*1.
        comment = 'number of sources allocated'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)



c allocate to allow the total number of sources for all frames

      allocate (RAlist(NSUMmax))
        if (.not.allocated(RAlist)) then
          write (6,*) 'allocation failed for RAlist'
          istat = 5
        end if

        allocate (DEClist(NSUMmax))
        if (.not.allocated(DEClist)) then
          write (6,*) 'allocation failed for DEClist'
          istat = 5
        end if

      allocate (deltaJD(NSUMmax))
        if (.not.allocated(deltaJD)) then
          write (6,*) 'allocation failed for deltaJD'
          istat = 5
        end if

        allocate (mdetIDlist(NSUMmax))
        if (.not.allocated(mdetIDlist)) then
          write (6,*) 'allocation failed for mdetIDlist'
          istat = 5
        end if

      allocate (tmpRAlist(NSUMmax))
        if (.not.allocated(tmpRAlist)) then
          write (6,*) 'allocation failed for tmpRAlist'
          istat = 5
        end if

        allocate (tmpDEClist(NSUMmax))
        if (.not.allocated(tmpDEClist)) then
          write (6,*) 'allocation failed for tmpDEClist'
          istat = 5
        end if

        allocate (tmpIDlist(NSUMmax))
        if (.not.allocated(tmpIDlist)) then
          write (6,*) 'allocation failed for tmpIDlist'
          istat = 5
        end if


      allocate (Rsat(NSUMmax,4))
        if (.not.allocated(Rsat)) then
          write (6,*) 'allocation failed for Rsat'
          istat = 5
        end if

      allocate (Rsatsub(NSUMmax,4))
        if (.not.allocated(Rsatsub)) then
          write (6,*) 'allocation failed for Rsatsub'
          istat = 5
        end if

      allocate (SatNum(NSUMmax,4))                    ! ADDED FEB 1
        if (.not.allocated(SatNum)) then                ! ADDED FEB 1
          write (6,*) 'allocation failed for SatNum'    ! ADDED FEB 1
          istat = 5                                     ! ADDED FEB 1
        end if  
      SatNum = 0

         allocate (SatNumSub(NSUMmax,4))                   
        if (.not.allocated(SatNumSub)) then               
          write (6,*) 'allocation failed for SatNumSub'   
          istat = 5                                     
        end if


c      allocate (Table_0(NSUMmax,32))
c      allocate (Table_0(NSUMmax,59))       ! JWF B30221
      allocate (Table_0(NSUMmax,65))       ! JWF B31209
        if (.not.allocated(Table_0)) then
          write (6,*) 'allocation failed for table_0'
          istat = 5
          call exit (istat)
        end if
      Table_0 = -99.99   ! initialize matrix with null values  !  TJ 23Aug2016
      ntable_0 = 0

      allocate (RAlist_0(NSUMmax))
        if (.not.allocated(RAlist_0)) then
          write (6,*) 'allocation failed for ralist_0'
          istat = 5
          call exit (istat)
        end if

      allocate (DEClist_0(NSUMmax))
        if (.not.allocated(DEClist_0)) then
          write (6,*) 'allocation failed for declist_0'
          istat = 5
          call exit (istat)
        end if

      allocate (IDlist_0(NSUMmax))
        if (.not.allocated(IDlist_0)) then
          write (6,*) 'allocation failed for idlist_0'
          istat = 5
          call exit (istat)
        end if

      allocate (N_Mstat(NSUMmax,4))
      if (.not.allocated(N_Mstat)) then
          write (6,*) 'allocation failed for N_Mstat'
          istat = 5
          call exit (istat)
        end if

      allocate (M_Mstat(NSUMmax,4))
        if (.not.allocated(M_Mstat)) then
          write (6,*) 'allocation failed for M_Mstat'
          istat = 5
          call exit (istat)
        end if

      allocate (WPROmagR(NSUMmax,4,3))
      if (.not.allocated(WPROmagR)) then
          write (6,*) 'allocation failed for WPROmagR'
          istat = 5
          call exit (istat)
        end if

      allocate (mLogQ(NSUMmax,4))
        if (.not.allocated(mLogQ)) then
          write (6,*) 'allocation failed for mLogQ'
          istat = 5
          call exit (istat)
        end if
c      allocate (Fmin(NSUMmax,4))
c        if (.not.allocated(Fmin)) then
c          write (6,*) 'allocation failed for Fmin'
c          istat = 5
c          call exit (istat)
c        end if
c      allocate (Fmax(NSUMmax,4))
c        if (.not.allocated(Fmax)) then
c          write (6,*) 'allocation failed for Fmax'
c          istat = 5
c          call exit (istat)
c        end if

c      allocate (DeltaMag(NSUMmax,4))
c      if (.not.allocated(DeltaMag)) then
c          write (6,*) 'allocation failed for DeltaMag'
c          istat = 5
c          call exit (istat)
c        end if

      allocate (KS(NSUMmax,4))
      if (.not.allocated(KS)) then
          write (6,*) 'allocation failed for KS'
          istat = 5
          call exit (istat)
        end if

      allocate (nDF(NSUMmax,4))
      if (.not.allocated(nDF)) then
          write (6,*) 'allocation failed for nDF'
          istat = 5
          call exit (istat)
        end if


      allocate (Rho12(NSUMmax))
        if (.not.allocated(Rho12)) then
          write (6,*) 'allocation failed for Rho12'
          istat = 5
          call exit (istat)
        end if
      allocate (Rho23(NSUMmax))
        if (.not.allocated(Rho23)) then
          write (6,*) 'allocation failed for Rho23'
          istat = 5
          call exit (istat)
        end if
       allocate (Rho34(NSUMmax))
        if (.not.allocated(Rho34)) then
          write (6,*) 'allocation failed for Rho34'
          istat = 5
          call exit (istat)
        end if
      allocate (P12(NSUMmax))
        if (.not.allocated(P12)) then
          write (6,*) 'allocation failed for P12'
          istat = 5
          call exit (istat)
        end if
      allocate (P23(NSUMmax))
        if (.not.allocated(P23)) then
          write (6,*) 'allocation failed for P23'
          istat = 5
          call exit (istat)
        end if
      allocate (P34(NSUMmax))
        if (.not.allocated(P23)) then
          write (6,*) 'allocation failed for P34'
          istat = 5
          call exit (istat)
        end if

      allocate (AVE_mJD(NSUMmax,4))
        if (.not.allocated(AVE_mJD)) then
          write (6,*) 'allocation failed for AVE_mJD'
          istat = 5
          call exit (istat)
        end if

      allocate (mJD_min(NSUMmax,4))
        if (.not.allocated(mJD_min)) then
          write (6,*) 'allocation failed for mJD_min'
          istat = 5
          call exit (istat)
        end if

      allocate (mJD_max(NSUMmax,4))
        if (.not.allocated(mJD_max)) then
          write (6,*) 'allocation failed for mJD_max'
          istat = 5
          call exit (istat)
        end if

      allocate (AVE_mJD_var(NSUMmax,4))
        if (.not.allocated(AVE_mJD_var)) then
          write (6,*) 'allocation failed for AVE_mJD_var'
          istat = 5
          call exit (istat)
        end if

      allocate (mJD_min_var(NSUMmax,4))
        if (.not.allocated(mJD_min_var)) then
          write (6,*) 'allocation failed for mJD_min_var'
          istat = 5
          call exit (istat)
        end if

      allocate (mJD_max_var(NSUMmax,4))
        if (.not.allocated(mJD_max_var)) then
          write (6,*) 'allocation failed for mJD_max_var'
          istat = 5
          call exit (istat)
        end if


c      allocate (slope(NSUMmax,4))
c        if (.not.allocated(slope)) then
c          write (6,*) 'allocation failed for slope'
c          istat = 5
c          call exit (istat)
c        end if

c      allocate (eslope(NSUMmax,4))
c        if (.not.allocated(eslope)) then
c          write (6,*) 'allocation failed for eslope'
c          istat = 5
c          call exit (istat)
c        end if

c      allocate (slopechi2(NSUMmax,4))
c        if (.not.allocated(slopechi2)) then
c          write (6,*) 'allocation failed for slopechi2'
c          istat = 5
c          call exit (istat)
c        end if


      allocate (medsky(NSUMmax,4))
        if (.not.allocated(medsky)) then
          write (6,*) 'allocation failed for medsky'
          istat = 5
          call exit (istat)
        end if
      allocate (medskysig(NSUMmax,4))
        if (.not.allocated(medskysig)) then
          write (6,*) 'allocation failed for medskysig'
          istat = 5
          call exit (istat)
        end if
      allocate (formalsig(NSUMmax,4))
        if (.not.allocated(formalsig)) then
          write (6,*) 'allocation failed for formalsig'
          istat = 5
          call exit (istat)
        end if

      allocate (tmpFLUX(nf,19))
        if (.not.allocated(tmpFLUX)) then
           write (6,*) 'allocation failed for tmpFLUX'
          istat = 5
        end if
      allocate (etmpFLUX(nf,19))
        if (.not.allocated(etmpFLUX)) then
           write (6,*) 'allocation failed for etmpFLUX'
          istat = 5
        end if

      mergetype = 1  !  inverse variance weighting
      if (doweight) then
            mergetype = 1  !  inverse variance weighting
      else
            mergetype = 0  !  straight averaging
      endif

      What = "mergetype"
        TYPE = "i"
        metval = mergetype
        comment = 'weighting scheme:  0 = none, 1 = inverse variance'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)

      allocate (tmpMAG(nf,19))
        if (.not.allocated(tmpMAG)) then
         write (6,*) 'allocation failed for tmpMAG'
        istat = 5 
        end if

      allocate (etmpMAG(nf,19))
        if (.not.allocated(etmpMAG)) then
           write (6,*) 'allocation failed for etmpMAG'
          istat = 5
        end if

      allocate (tmpFLG(nf,19))
        if (.not.allocated(tmpFLG)) then
           write (6,*) 'allocation failed for tmpFLG'
          istat = 5
        end if



      allocate (ztmp(nf))
        if (.not.allocated(ztmp)) then
           write (6,*) 'allocation failed for ztmp'
          istat = 5
        end if
      allocate (eztmp(nf))
        if (.not.allocated(eztmp)) then
           write (6,*) 'allocation failed for eztmp'
          istat = 5
        end if
      allocate (peztmp(nf))
        if (.not.allocated(peztmp)) then
           write (6,*) 'allocation failed for peztmp'
          istat = 5
        end if
      


      allocate (GlobalStats(nf,2))
      if (.not.allocated(GlobalStats)) then
        write (6,*) 'allocation failed for GlobalStats'
         istat = 5
        end if

c===================== beginning of code added by CJG B30213 =========================

c allocate an multiple epoch table

      if(mepwrite)then

c === TPC
c          NF/5 = approx. 2*avg. depth for a given source, usually.
c          Most of the time this should be more than enough.
c          If it needs to be bigger it will be reallocated later.

           nfmep = nf/5

           if(nfmep .lt. 10) then
               if(nf .le. 10) nfmep = nf
               if(nf .gt. 10) nfmep = 10 
           end if

         write (6,*) 'Initial allocation for MEP tables = ',nsummax, nfmep
c === TPC

         allocate (mepTable(nsummax,20,nfmep))
         if (.not.allocated(mepTable)) then
            write (6,*) 'allocation failed for mepTable'
            istat = 5
            call exit (istat)
         endif
         
         mepTable = -99.99  ! TJ 23aug2016  
         
         if (debug) write (6,*) 'allocate name table ',nsummax, nfmep
         
         allocate (namTable(nsummax,nfmep))
         if (.not.allocated(namTable)) then
            write (6,*) 'allocation failed for namTable'
            istat = 5
            call exit (istat)
         endif

           namtable = '99999z999'

         if (debug) write (6,*) 'allocate jdtable table ',nsummax, nfmep

         allocate (jdTable(nsummax,nfmep))
         if (.not.allocated(jdTable)) then
            write (6,*) 'allocation failed for jdTable'
            istat = 5
            call exit (istat)
         endif

         jdtable = -99.d19

        else
c         allocate (mepTable(nmax,20,nfmep), namTable(nmax,nfmep), jdTable(nmax,nfmep),
c    +              stat = kStat)           ! JWF B30404
          allocate (mepTable(1,1,1), namTable(1,1), jdTable(1,1),
     +              stat = kStat)           ! JWF B30404
          if (kStat .ne. 0) then
            print *,'***ERROR: cannot allocate mep arrays; kStat =',kStat
            istat = 64
          end if
      endif

c===================== end of code added by CJG B30213 =====================


cccccccccccc
c quit if any allocation failures

      if (istat.gt.0) then
            write (6,*) 'quitting ...'
            call exit (istat)
      endif

c      if (SPIT) then
ccccc  get null fits header
c        write (6,'(a,a)') 'Hnull = ',Hnull(1:72)
c       call  readFhead(Hnull,HeadNull)
c      endif


cccccccccccccccccccccccccccccccccccccccccccc
c  load table of sources

      What = "ndet"
        TYPE = "i"
        comment = 'number of detections'

      NsrcAll = 0

      write (What,'(a)') 'ndet'

      nsrc = 0
      call loadtable (NSUMmax,mdex,nmax,RAlist,DEClist,mdetIDlist,
     +                  Rsat,J,nf,nsrc,NsrcAll,unitest)
      if (verbose) write (6,*) nsrc,' RA/DEC positions loaded'

      Rsatsub = 0.0

      if ( nsrc . le. 0 ) then
            write (6,*) nsrc,' RA/DEC positions loaded'
            write (6,*) 'caution -- src table empty'
      endif

       metval = nsrc* 1.
       iband = 0
       if (smode) iband = ib0
       call MetWr (imeta, iband, What, TYPE, metval, comment)


ccccccccccccccccccccc

      What = "NsrcAll"
        TYPE = "i"
        metval = NsrcAll*1.
        comment = 'number of sources loaded'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c load in the PSFs

      if (verbose) write (6,*) 'loading PSFs'

      do ib=1,4

        do j=1,npn(ib)                           ! MODIFIED/MOVED

c       if ( wflag ( 1,ib ) .eq.1 ) then
       if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

          fname = Pnames(j,ib)
          MM = numchar(fname)
c          if (verbose) write (6,'(a)') fname(1:MM)
          write (6,'(a)') fname(1:MM)
          call readimage
     1          (nnx,nny,lsize,Larray,iLarray,fname,pcrval1,pcrval2,
     1          pcdelt1,ppix(ib),pcrot,pxcent(ib),pycent(ib), tJD)

         if (debug) write (6,'(2i4,1x,a,2i5,f12.7,2f7.2)') ib,j,fname(1:MM),nnx,nny,ppix(ib),pxcent(ib),pycent(ib)
          if (verbose) write (6,*) ib,'  PSF loaded with size & pix scale: ',nnx,nny,ppix(ib)

           npsf(ib) = nnx

           ic = 0
           do jj=1,nny
           do ii=1,nnx
                ic=ic+1
                PSFs(ii,jj,j,ib) = Larray(ic)
           enddo
           enddo

        fname = Punc(j,ib)
          MM = numchar(fname)
c         if (verbose) write (6,'(a)') fname(1:MM)

        call readimage
     1     (nnx,nny,lsize,Larray,iLarray,fname,
     1     zz,zz,
     1     zz,zz,zz,zz,zz,tJD)

c         if (verbose) write (6,'(2i4,1x,a,2i5)') ib,j,fname(1:MM),nnx,nny
c          if (verbose) write (6,*) ib,'  PSFunc loaded with size: ',nnx,nny


         if (nnx.ne.npsf(ib)) then
            write (6,*) 'ERROR -- PSF and PSFunc have mismatched dimensions'
            call exit (64)
         endif

         ic = 0
           do jj=1,nny
           do ii=1,nnx
                ic=ic+1
      
            if (SPIT) then   ! SEP

               if (NEP) then
                 if (ib.eq.1) Larray(ic) = Larray(ic) * 1.0
                     if (ib.eq.2) Larray(ic) = Larray(ic) * 0.8
                     if (ib.eq.3) Larray(ic) = Larray(ic) * 0.7
                     if ((ib.eq.4).and.(.not.MIPS)) Larray(ic) = Larray(ic) * 0.7

               else if (SEP) then
                 if (ib.eq.1) Larray(ic) = Larray(ic) * 1.05
                 if (ib.eq.2) Larray(ic) = Larray(ic) * 1.25
                 if (ib.eq.3) Larray(ic) = Larray(ic) * 1.4
                 if ((ib.eq.4).and.(.not.MIPS)) Larray(ic) = Larray(ic) * 1.5

               else
                 if (ib.eq.1) Larray(ic) = Larray(ic) * 0.25
                     if (ib.eq.2) Larray(ic) = Larray(ic) * 0.25
                     if (ib.eq.3) Larray(ic) = Larray(ic) * 0.5
                     if ((ib.eq.4).and.(.not.MIPS)) Larray(ic) = Larray(ic) * 0.5
               endif

            
            endif

                PSFuncs(ii,jj,j,ib) = Larray(ic)  * PSFuncscale(ib)

           enddo
           enddo

        endif

         enddo

         enddo


ccccccccccccccccccccccc  load the aperture corrections
c         Ncorrdim = calgridX * calgridY

c         allocate (AppCorr (Ncorrdim,4))   !  grid location, band

c         if (.not.allocated(AppCorr)) then
c          write (6,*) 'allocation failed for AppCorr'
c          istat = 5
c          call exit (istat)
c         end if

c         call LoadApCorr (Raper,apcorrdir,calbname,calgridX,calgridY,Ncorrdim,AppCorr,napcorr,smode,verbose)

            
c       if (SPIT) then
c      ! aperture corrections for R=3 pixels, Annulus = 10-20 pixels
c            AppCorr (1,1) = -0.1153
c            AppCorr (1,2) = -0.1162
c            AppCorr (1,3) = -0.12788
c            AppCorr (1,4) = -0.2141
c
c            if (MIPS) then
c                  AppCorr (1,3) = -0.2141   ! IRAC-4
c                  AppCorr (1,4) = -0.1658       ! MIPS-24, R = 13"
c                   AppCorr (1,4) = -0.54239      ! MIPS-24, R=6.6"
c            endif
c
c       endif



ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c namelist logistics

      allocate (Tarray(nx,ny,2))
        if (.not.allocated(Tarray)) then
          write (6,*) 'allocation failed for Tarray'
          istat = 5
        call exit(istat)
        end if
      allocate (iTarray(nx,ny))
        if (.not.allocated(iTarray)) then
          write (6,*) 'allocation failed for iTarray'
          istat = 5
        call exit(istat)
        end if


cccccccccccccccc
ccc  Get the coadd headers

      if (verbose) write (6,*) 'load the coadd FITS headers'
      
      do ib=1,4
                cexist(ib) = .false.
                sexist(ib) = .false.
      enddo

        LC = numchar(coname)
        LL = numchar (clevel)

      cband0 = 0

        do ib=1,4

            cexist(ib) = .false.
            sexist(ib) = .false.
            ncx(ib) = 0
                ncy(ib) = 0
            WCSc(ib) = -1
            scrot(ib) = 0  ! TEMPORARY -- get this value from the header

c            if (wflag(1,ib) .eq. 1) then
            if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

              write (cc,'(i1)') ib
              if ((MIPS).and.(ib.eq.4)) then
                        cc = '5'
               endif

              coaddfits(ib) = coname(1:LC) // '-w' // cc // '-img-' // clevel(1:LL) // '.fits'  ! unwise

              L = numchar(coaddfits(ib))
              if (verbose) write (6,'(a,a)')'coadd ', coaddfits(ib)(1:L)

              zexist = .false.
                  erase = .false.
c                 call access(coaddfits(ib),zexist,erase)                              ! JWF B60711
                  zexist = Access(coaddfits(ib)(1:LNBlnk(coaddfits(ib))),' ') .eq. 0   ! JWF B60711
                  if (.not.zexist) then
                   ! see if the file is compressed
                    gztmp = coaddfits(ib)(1:L) // '.gz'
                    zexist = .false.
                    erase = .false.
c                   call access(gztmp,zexist,erase)                      ! JWF B60711
					zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
                    if (zexist) then
                        coaddfits(ib) = gztmp
                    endif
                  endif

ccccccccccc
            zexist = .false.
                erase = .false.
c               call access(coaddfits(ib),zexist,erase)                              ! JWF B60711 
                zexist = Access(coaddfits(ib)(1:LNBlnk(coaddfits(ib))),' ') .eq. 0   ! JWF B60711

            if (zexist) then

               if(cband0 .lt. 1) cband0 = ib 

               cexist(ib) = .true.
               call  readFhead(coaddfits(ib),Hfits)


               command = 'NAXIS1'
               val = 0.
                 call keyhead (Hfits,command,s0)
                 L = numchar(s0)
                 if (L.gt.0) then
                     read (s0,*) val
                 endif
                 ncx(ib) = int(val)

      
               command = 'NAXIS2'
               val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                   endif
                   ncy(ib) = int(val)


               command = 'CDELT2'
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) cd2
               else ! look for the CD matrix
                        call CDmatrix (Hfits, cd1, cd2)
                   endif


               cscale(ib) = abs(cd2*3600.)
               if (verbose) write (6,*) 'coadd pixel scale: ',ib,cscale(ib)

               command = 'CRVAL1'
               val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                   endif
                   crval1  = val

               command = 'CRVAL2'
               val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                   endif
                   crval2  = val


               command = 'CROTA2'
               val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                   endif
                   scrot(ib)  = val

            
              if (.not.SPIT) then
               command = 'MAGZP'
               val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                     zero(ib)  = val
                 cozero(ib) = zero(ib)
               endif

               DATE = 'null'
               command = 'DATE'
               call keyhead (Hfits,command,DATE)
               L = numchar(DATE)
               if (L.le.0) then
                  DATE = 'null'
               endif

              endif


c      write (6,*) crval1,crval2

                  Ra0 = crval1 * 1.d0
                  Dec0 = crval2 * 1.d0
                  call WhatPos (1,WCSc(ib),ncx(ib),ncy(ib), cscale(ib), Hfits, ra0, dec0, x0, y0, igo, edgebuf)


              if (verbose) write (6,'(i2,2i6,2f8.2,a,f8.4)') ib,ncx(ib),ncy(ib),x0, y0,' coadd pixel scale: ',cscale(ib)


            else
               write (6,*) 'ERROR -- coadd does not exist:'
                   L = numchar(coaddfits(ib))
                   write (6,'(a)') coaddfits(ib)(1:L)
               write (6,*) '   exiting'
                   call exit(64)

              endif
         endif

      enddo ! ib=1,4
      
      print *,'Selected band ',cband0,' as coadd representative band.'

        if(cband0 .lt. 1) then
           write(6,*) '***ERROR: No valid Coadd found for assigning cband0.'
           call exit(64)
        endif

      do ib=1,4

        if (ib.eq.1) What = "zero"
        if (ib.eq.2) What = "zero"
        if (ib.eq.3) What = "zero"
        if (ib.eq.4) What = "zero"

c        What = "zero"

        TYPE = "r"
        metval = zero(ib)
        iband = ib
        if (ib.eq.1) comment = 'zero point magnitude, w1'
        if (ib.eq.2) comment = 'zero point magnitude, w2'
        if (ib.eq.3) comment = 'zero point magnitude, w3'
        if (ib.eq.4) comment = 'zero point magnitude, w4'


        if (smode) then
                if (ib.eq.ib0) then
                  call MetWr (imeta, iband, What, TYPE, metval, comment)
                endif
        else
                call MetWr (imeta, iband, What, TYPE, metval, comment)

        endif

      if (ib.eq.1) What = "cozero" 
        if (ib.eq.2) What = "cozero"
        if (ib.eq.3) What = "cozero"
        if (ib.eq.4) What = "cozero"

        metval = cozero(ib)
        if (ib.eq.1) comment = 'coadd zero point magnitude, w1'
        if (ib.eq.2) comment = 'coadd zero point magnitude, w2'
        if (ib.eq.3) comment = 'coadd zero point magnitude, w3'
        if (ib.eq.4) comment = 'coadd zero point magnitude, w4'

        if (SPIT) comment = 'IRAC coadd zero point magnitude'

        if (MIPS) then
          if (ib.eq.3) metval = cozero(4)
          if (ib.eq.4) metval = cozero(5)
          if (ib.eq.4) comment = 'MIPS coadd zero point magnitude'
        endif

        if (smode) then
                if (ib.eq.ib0) then
                  call MetWr (imeta, iband, What, TYPE, metval, comment)
                endif
        else
                call MetWr (imeta, iband, What, TYPE, metval, comment)

        endif

        enddo


        zmin = 9999.
        do ib=1,4
                zmin = min (zero(ib),zmin)
        enddo

        if (zmin.le.0.) then
                write (6,*) 'ERROR -- zero mags incorrectly set ',zero
                call exit(64)
        endif

        if ((nx.le.0).or.(ny.le.0)) then
                write (6,*) 'ERROR -- initial array size too small ',nx,ny
                call exit(64)
        endif


cccccccccccccccc
c  now frame names and  get the pixel scale

      if (verbose) write (6,*) 'load frame headers into memory'

      LL = numchar(level)

        do j=1,nf

          L1 = numchar(path(j))
          L2 = numchar(basename(j))
          L3 = numchar(intword)

        WCS(J) = -1

          do 6000 ib=1,4

        ib0=ib
         if (MIPS) then
            if (ib.eq.3) then
                  ib0 = 4
            else if (ib.eq.4) then
                        ib0 = 5
            endif
         endif

           write (cc,'(i1)') ib0
         fram(j,ib) = 'null'
         msk (j,ib) = 'null'
         uncn (j,ib) = 'null'

         preorder (j,ib) = 0

           if (wflag(j,ib).eq.1) then

                fram(j,ib) = path(j)(1:L1) // '/' //
     1          basename(j)(1:L2) // '-w' // cc //
     1          '-img-' // level(1:LL) // '.fits'  ! unwise
            L = numchar(fram(j,ib))
                
            msk(j,ib) = path(j)(1:L1) // '/' //
     1          basename(j)(1:L2) // '-w' // cc //
c     1          '-invvar-' // level(1:LL) // '.fits'       !  TJ 07Sep2017
     1          '-std-' // level(1:LL) // '.fits'       !  TJ 07Sep2017
            L = numchar(msk(j,ib))

              uncn(j,ib) = path(j)(1:L1) // '/' //
     1          basename(j)(1:L2) // '-w' // cc //
     1          '-std-' // level(1:LL) // '.fits'    ! unwise
                L = numchar(uncn(j,ib))


            zexist = .false.
            erase = .false.
c           call access(fram(j,ib),zexist,erase)
			zexist = Access(fram(j,ib)(1:LNBlnk(fram(j,ib))),' ') .eq. 0   ! JWF B60711

            if (.not.zexist) then
              ! see if the file is compressed
              gztmp = fram(j,ib)(1:L) // '.gz'
              zexist = .false.
              erase = .false. 
c              call access(gztmp,zexist,erase)
			zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711

              if (zexist) then
                  fram(j,ib) = gztmp
                  L = L + 3
              endif

            endif


c load the FITS header into memory

                if (verbose) write (6,'(a,i1,1x,i5,1x,a,a)') 
     1               'ib,j=',ib,j,'frame=',fram(j,ib)(1:L)

                if (zexist) then
                  call  readFhead(fram(j,ib),Hfits)
                  MM = numchar(Hfits)
                  HDR0 (j,ib) = Hfits
                else
                   if (SPIT) then
                     Hfits = HeadNull
                     MM = numchar(Hfits)
                     HDR0 (j,ib) = Hfits
                   else
                     LLL = numchar (fram(j,ib))
                     write (6,*) '===WARNING -- INT frame does not exist: ',fram(j,ib)(1:LLL)
                 fram(j,ib) = 'null'
                 msk (j,ib) = 'null'
                     uncn (j,ib) = 'null'
                 wflag(j,ib) = 0
                 HDR0 (j,ib) = 'null'
                 call exit(64)         ! JWF B80524
                 goto 6000
c                     write (6,*) '  exiting '
c                     call exit(9)
                  endif

                endif


            zexist = .false.
                erase = .false.
c                call access(uncn(j,ib),zexist,erase)
			zexist = Access(uncn(j,ib)(1:LNBlnk(uncn(j,ib))),' ') .eq. 0   ! JWF B60711
            if (.not.zexist) then
                  ! see if the file is compressed
              gztmp = uncn(j,ib)(1:L) // '.gz'
                  zexist = .false.
                  erase = .false.
c                  call access(gztmp,zexist,erase)
 		    	  zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
              if (zexist) then
                        uncn(j,ib) = gztmp
                  endif

                endif
            if(.not.zexist) then
                     LLL = numchar (fram(j,ib))
                     write (6,*) '***ERROR -- UNC frame does not exist: ',fram(j,ib)(1:LLL)
                 call exit(64)
            endif


            zexist = .false.
                erase = .false.
c                call access(msk(j,ib),zexist,erase)
 		    	 zexist = Access(msk(j,ib)(1:LNBlnk(msk(j,ib))),' ') .eq. 0   ! JWF B60711
                if (.not.zexist) then
                  ! see if the file is compressed
		  L = numchar(msk(j,ib))  ! TJ 07Sep  ! TJ 07Sep
                  gztmp = msk(j,ib)(1:L) // '.gz'
                  zexist = .false.
                  erase = .false.
c                 call access(gztmp,zexist,erase)
 		    	  zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
                  if (zexist) then
                        msk(j,ib) = gztmp
                  endif
                endif
            if(.not.zexist) then
                     LLL = numchar (msk(j,ib))
                     write (6,*) '***WARNING -- MSK frame does not exist: ',msk(j,ib)(1:LLL)
		     msk(j,ib) = "null"  
c                 call exit(9)
                  call exit(64)
            endif


            if ((zexist).and.(ireg.gt.1)) then   ! dead-zone mitigation
              
                command = 'WCROTA2'
                   val = 0.
               ybuf = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                 if (abs(val).gt.0.) then 
                  
                  if (abs(Dec0).gt.85.) then   !  the WCROTA2 is becoming chaotic at the equ poles
                   val = 45.   !  assume the worst case and bump up the buffer
                  endif


                  call polygrid(ireg,val,ybuf)
c                  write (6,*) 'rotation angle, buffer = ',val,ybuf

                  if (nint(ybuf).gt.nbuf) then
                  nbuf   = nint(ybuf)  ! bump up the grid-to-grid overlap to remove geometric dead zones
                  endif
                 endif
               endif

            endif


            if ( (zexist).and.(.not.gotscale(ib))) then
c            if (j.eq.1) then

               gotscale(ib) = .true.

               command = 'NAXIS1'
                   val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0) 
                   if (L.gt.0) then
                     read (s0,*) val
                   endif
                   nsx(ib) = int(val)


                   command = 'NAXIS2'
                   val = 0.
                   call keyhead (Hfits,command,s0)
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) val
                   endif
                   nsy(ib) = int(val)

                   command = 'CDELT2'
                   call keyhead (Hfits,command,s0) 
                   L = numchar(s0)
                   if (L.gt.0) then
                     read (s0,*) cd2
                   else ! look for the CD matrix
                        call CDmatrix (Hfits, cd1, cd2)
                   endif


              pscale (ib) = abs(cd2 * 3600.)
              if (verbose) write (6,*) 'frame pixel scale: ',ib,pscale (ib)
            endif

           endif ! wflag(j,ib)

 6000       enddo      ! ib Band

        enddo ! j frame

      TYPE = "i" 
        metval = nbuf
        write (comment,'(a,i3)') 'buffer size in frame pixels'
      write (What,'(a)') "nbuf"
        iband = 0
        if (smode) iband = ib
        call MetWr (imeta, iband, What, TYPE, metval, comment)

      write (6,*) 'nbuf = ',nbuf,' w1-3 frame pixels'

      ! Make buffer uniform accross bands (esp. w4) -- TPC

      buf_as = nbuf*pscale(1)*3600

      write (6,*) 'buf = ',buf_as,' arc-seconds'

        do ib=1, 4
        nbuf_band(ib) = nint(buf_as/3600/pscale(ib))
      enddo

      write (6,*) 'nbuf per band = ',nbuf_band,' pixels'

c      call exit(0)  ! TEMP

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Divide the WISE footprint

      if (verbose) write (6,*) 'Divisions: ',ireg,jreg

      nallocate = 0
      nsubxX = 0
        nsubyY = 0
      iquad_frame = 0
        jquad_frame = 0
      ra_reg = 0.
        dec_reg = 0.

      dosolo = .true.   !  this means only one grid is happening
      if (ireg*jreg.gt.1) dosolo = .false.

        do 8000 kreg_ix=1,ireg*jreg

c        Do regions in the specified order. usually kreg = kreg_ix
         kreg = region_order(kreg_ix)
c        If this spot in the region order array is zero, terminate region processing
       if(kreg .le. 0) then
         print *
         print *,'!!! Terminating region processing.'
           goto 8010
         endif

c        If doreg > 0, a single region has been requested -- TPC
         if(doreg .gt. 0) then
           if(kreg .ne. doreg) then
           print *
           print *,'!!! Skipping region ',kreg
             goto 8000
           endif
         endif

       if (verbose) then
         write (6,*) ' '
         write (6,*) 'Region: ',kreg,' of ',ireg*jreg,' ;  frames = ',nf
       endif

c      write (6,*) 'sitting at kreg loop; hit enter'
c      read (5,'(a)') s0


c      command = "free -ot"
c        s0 = ""
c        call unix (command,s0)

       igot = 0
c       nsubxX = 0
c       nsubyY = 0

       dd1 = 0.
         call  ETime (Second, dd1)
         dd1 = second(1)  ! User time

       nactive_band(1:4) = 0
       fractive(:,1:4) = 0

       do 8001 ib=1,4

        ohyea = 0
        ntmp (ib) = 0


c        if (wflag(1,ib).eq.1) then
        if (any(wflag(1:nf,ib)==1)) then

         ntmp (ib) = nf

         do Jfr = 1,nf
          if (WCS(Jfr).ge.0) call wcsclose(WCS(Jfr))
          WCS(Jfr) = -1
           enddo

           nbuf_coadd = nint(buf_as / 3600 / cscale(ib))

           call divider (kreg,ncx(ib),ncy(ib),ireg,jreg, i0, j0, isize, jsize)
c           write (6,*) 'Zcoadd division ',ib, kreg, i0, j0, isize, jsize ! TEMP


           il = i0 - (isize/2.0)
           il=max(il,1)
           ih = (il + isize) - 1
           jl = j0 - (jsize/2.0)
           jl=max(jl,1)
           jh = (jl + jsize) - 1

! this is the hard limit used to select extracted sources; see Table_0 biz


           il_0 = il - 3
           ih_0 = ih + 3
           jl_0 = jl - 3
           jh_0 = jh + 3

           il_0=max(il_0,1)
           jl_0=max(jl_0,1)
           ih_0=min(ih_0,ncx(ib))
           jh_0=min(jh_0,ncy(ib))

      if (verbose) then
      write (6,*) 'limits'
      write (6,*) ' selecting extracted sources'
      write (6,*) il_0,ih_0,jl_0,jh_0
      endif

! add a buffer zone to deal with edge-effects

c           il = il - (nbuf_coadd/2)
c           ih = ih + (nbuf_coadd/2)
c           jl = jl - (nbuf_coadd/2)
c           jh = jh + (nbuf_coadd/2)

         move = nbuf_coadd
         il = il - move
           ih = ih + move
           jl = jl - move
           jh = jh + move

         il=max(il,1)
         jl=max(jl,1)
         ih=min(ih,ncx(ib))
         jh=min(jh,ncy(ib))

      if (verbose) then
        write (6,*) 'selecting input sources'
        write (6,*) il,ih,jl,jh
      endif

      
c      write (6,*) (cscale(ib) / pscale(ib))*((ih-il)+1),(cscale(ib) / pscale(ib))*((jh-jl)+1)
c      write (6,*) isize * cscale(ib) / pscale(ib), jsize * cscale(ib) / pscale(ib)

c          write (6,*) il,ih,jl,jh

c           Expand region size to include a buffer.
c           This usage of nbuf_band makes the w4 region the same size (in true angle)
c           as the other bands. --- TPC      
          nsubx(ib) = nint(isize * cscale(ib) / pscale(ib)) + (nbuf_band(ib)*2)
            nsuby(ib) = nint(jsize * cscale(ib) / pscale(ib)) + (nbuf_band(ib)*2)

          nsubx(ib) = min (nsubx(ib), nsx(ib))
          nsuby(ib) = min (nsuby(ib), nsy(ib))

cc check band 4 for out of bounds
ccc         This appears to make the max. frame region size that of the first
ccc         region in the first band processed. This will also be the size
ccc         allocated for frame pixels in all bands (nsubxx, nsubyy). --- TPC
          if (ib.gt.1) then
            if (nsubxX.gt.0) then
                  nsubx(ib) = min (nsubx(ib), nsubxX)
            endif
            if (nsubyY.gt.0) then
                        nsuby(ib) = min (nsuby(ib), nsubyY)
                endif
          endif



          i0_reg(ib) = nsubx(ib)/2.
            j0_reg(ib) = nsuby(ib)/2.

c            if (verbose) write (6,*) 'frame sub size = ',nsubx(ib),nsuby(ib)
          if (verbose) write (6,'(a,i3,i3,2i6,2i6,2i5)') 'coadd division (reg, band, i0,j0, Xcosize,Ycosize, nsubx,nsuby): '
     1      ,kreg,ib,i0, j0, isize, jsize,nsubx(ib),nsuby(ib)


c         if (igot.eq.0) then !  get the sublist of sources

c region center
            xpix = i0 * 1.d0
            ypix = j0 * 1.d0

            call pix2wcs (WCSc(ib), xpix, ypix, RA0, DEC0)
          if (debug) write (6,'(a,2f10.5)') 'Region center coords: ',RA0, DEC0
         
          ra_reg = RA0 * 1.0
          dec_reg = DEC0 * 1.0


        dd90 = 0.
         call  ETime (Second, dd90)
         dd90 = second(1)  ! User time

c      dosolo = .true.  ! TEMP  !!!

      if (dosolo) then  !  if only one grid, then don't bother with re-ordering
       nactive = nf
       do Jfr=1,nf
         order (Jfr) = Jfr
         if(wflag(jfr,ib) .eq. 1) then
             nactive_band(ib) = nactive_band(ib) + 1
             fractive(jfr,ib) = 1
           endif
       enddo
       goto 524
      endif

c          if ((ib.eq.1).and.(verbose)) write (6,*) 'BEGIN frame Re-Order book-keeping  ; nf = ',nf

        nactive=0  ! Nactive as used here is per-band; below it will be reused for all bands
          do Jfr=1,nf
           if (fram(Jfr,ib)(1:4).ne.'null')  then
            Hfits = HDR0 (Jfr, ib)
            if (WCS(Jfr).ge.0) call wcsclose(WCS(Jfr))
            WCS(Jfr) = -1
            call wcsinit(Hfits,WCS(Jfr))

            call wcs2pix(WCS(Jfr), ra0, dec0, xpix, ypix, offscl)

            xquad_frame = xpix * 1.0
            yquad_frame = ypix * 1.0

c            write (6,*) 'grid center pix (frame space):',Jfr,xquad_frame, yquad_frame

            iquad_frame = nint(xquad_frame)
            jquad_frame = nint(yquad_frame)

            if (nsy(ib).eq.nsuby(ib)) jquad_frame = j0_reg(ib)
            if (nsx(ib).eq.nsubx(ib)) iquad_frame = i0_reg(ib)

             igg = 0
             do j=1,nsy(ib)
                 jdel_yz = j - jquad_frame
                 jdex = j0_reg(ib) + jdel_yz
               do 520 i=1,nsx(ib)
                 idel_xz = i - iquad_frame
                 idex = i0_reg(ib) + idel_xz

                 if (idex.lt.1) goto 520
                 if (jdex.lt.1) goto 520
                 if (idex.gt.nsubx(ib)) goto 520
                 if (jdex.gt.nsuby(ib)) goto 520

                  igg=1
              goto 521
 520           continue
             enddo
 521         if (igg.eq.1) then
c            write (6,*) '** active frame: ',Jfr
            nactive = nactive + 1
            torder (nactive) = Jfr
                fractive(jfr,ib) = 1
            order (nactive) = 0
            tnew (nactive) = 0
           endif

           endif
       enddo                  ! Jfr

       write (6,*) 'total number of active frames for band ',ib,': ',nactive,'  out of: ',nf

       nactive_band(ib) = nactive

c   zero out the frames that are not active
        do jj = nactive+1,nf
            torder (jj) = 0
            order (jj) = 0
            tnew (jj) = 0
        enddo

c        do jord = 1,nactive
c            write (6,*) jord,torder(jord)
c        enddo


ccccc now see which active frames are new compared to *preorder*
ccccc  but read the preorder list according to what was last read-in -- i.e., reverse the order

      nord = 0
      nnew = 0

      do 522 kord = nf,1,-1
        if (preorder(kord,ib).eq.0) goto 522

        do jord = 1,nactive
            if (torder(jord).eq.preorder(kord,ib)) then
              ! active frame in previous list
              nord = nord + 1
              order (nord) = torder(jord)
              torder(jord) = 0
              goto 522
            endif
        enddo
 522      continue

c      write (6,*) ' '  ! --- TPC added debug vvvvvvvvvvvvv
c        write (6,*) 'Loading ',nord,' old frames in band ',ib,': ',order(1:nord)
c        do jj=1,nord
c                write (6,*) order (jj)
c        enddo


 ! now see what is new

      do jord = 1,nactive
            if (torder(jord).gt.0) then
              nnew = nnew+1
                  tnew (nnew) = torder(jord)
              nord = nord + 1
              order (nord) = torder(jord)
            endif
      enddo

c       write (6,*) ' '  ! --- TPC added debug vvvvvvvvvvvv
c         write (6,*) 'Loading ',nnew,' new frames in band ',ib,': ',tnew(1:nnew)
c         do jj=1,nnew
c                write (6,*) tnew (jj)
c         enddo
         

c  now fill the preorder
       do jord = 1, nord
            preorder (jord,ib) = order (jord)
       enddo
       do jord = nord+1,nf
            preorder (jord,ib) = 0
       enddo


       write (6,*) ' ' ! --- TPC added debug
         write (6,*) 'New pre-order counts for band ',ib,
     1               '; nactive, nord (should agree) = ',nactive,nord
c       if ((ib.ne.0).and.(verbose))  then
c         write (6,*) 'Frame Read Order:',ib
c           do jj=1,nactive
c                write (6,'(2i5)') preorder (jj,ib), order (jj)
cc             write (6,'(i5)') order (jj)
c           enddo
c           write (6,*) 'END Re-Order book-keeping, now apply '
c        endif


  ! TEMP
c        if  (kreg.eq.11) call exit(0)
c        goto 8000


        ntmp (ib) = nactive

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Target of 'if(dosolo)'

 524         if (igot.eq.0) then  ! first band with frames: get the sublist of sources

          print *,'Final boundaries for input src selection: ',
     1              'band = ',ib,', il,ih,jl,jh = ', il,ih,jl,jh

          igot = 1
          nsubsrc = 0
          do Jsrc = 1,nsrc

              Ra0 = RAlist(Jsrc)
                  Dec0 = Declist(Jsrc)

                  call WhatPos (1,WCSc(ib),ncx(ib),ncy(ib), cscale(ib), Hfits, ra0, dec0, x0, y0, igo, edgebuf)

c      if ( (RA0.gt.180.1435).and.(RA0.lt.180.1437).and.(Dec0.gt.0.995).and.(Dec0.lt.0.996)) then
c            write (6,'(a,2f10.5,i5,2f9.1,i3,4i6)') '******got it ',RA0,Dec0,kreg,x0,y0,igo,il,ih,jl,jh
c            write (49,'(a,2f10.5,i5,2f9.1,6i5)') 'got it ',RA0,Dec0,kreg,x0,y0,igo,il,ih,jl,jh
c      endif

              ix0 = nint(x0)
              iy0 = nint(y0)

              if ( (ix0.ge.il).and.(ix0.le.ih).and.(iy0.ge.jl).and.(iy0.le.jh) ) then
                  nsubsrc = nsubsrc + 1
                  tmpRAlist (nsubsrc) = Ra0
                  tmpDeclist (nsubsrc) = Dec0
                  tmpIDlist (nsubsrc) = mdetIDlist(Jsrc)

                  do IK = 1,4
                    Rsatsub (nsubsrc,IK) = Rsat(Jsrc,IK)
                  enddo

c                  write (66,'(2f10.5, 2i6)') Ra0,dec0,ix0,iy0  ! TEMP
                  
              endif
          enddo

          if ((kreg_ix.eq.1).and.(verbose)) write (6,*) nsubsrc,'  sources loaded for 1st sub-region ', kreg

          if (kreg.lt.10) then
            write (What,'(a,i1)') "nsubsrc_",kreg
          else if (kreg.lt.100) then
            write (What,'(a,i2)') "nsubsrc_",kreg
          else
            write (What,'(a,i3)') "nsubsrc_",kreg
          endif

             TYPE = "i"
             metval = nsubsrc
                  write (comment,'(a,i2)') 'source loaded, within region: ',kreg
             iband = 0
           if (smode) iband = ib
             call MetWr (imeta, iband, What, TYPE, metval, comment)


           nsubsrcA = nint(nsubsrc * adb_alloscale )  ! allow room for active deblending

          if (nallocate.eq.0) then
                nsubxX = nsubx(ib)
                nsubyY = nsuby(ib)
            endif

        endif ! igot

       endif  !  any(wflag(ib))

         if (verbose) write (6,*) 'frame sub size = ',ib,nsubx(ib),nsuby(ib),nsubxX, nsubyY

 8001    continue   ! ib


ccc run through each order(band) and create one list

      if (.not.dosolo) then

      nord = 0
      ib1 = 1

      write (6,*) 'Per band frame counts = ',ntmp

      write (6,*) 'Frame index lists per band, before merging bands:'
      write (6,*) '  band 1 frame indicies = ',(preorder(i,1),i=1,ntmp(1))
      write (6,*) '  band 2 frame indicies = ',(preorder(i,2),i=1,ntmp(2))
      write (6,*) '  band 3 frame indicies = ',(preorder(i,3),i=1,ntmp(3))
      write (6,*) '  band 4 frame indicies = ',(preorder(i,4),i=1,ntmp(4))

      do k=1,ntmp(ib1)
            ival0 = preorder (k,ib1)
            nord=nord + 1
            torder (nord) = ival0

            idex = 0
            ib2 = 2
            do k2=1,ntmp(ib2)
              if (preorder (k2,ib2).eq.ival0) then
                  idex = k2
                  goto 321
              endif
            enddo

 321            if (idex.gt.0) then  ! W2 common value with W1
                preorder (idex,ib2) = 0
            endif      

            idex = 0
                ib3 = 3
                do k3=1,ntmp(ib3)
                  if (preorder (k3,ib3).eq.ival0) then
                        idex = k3
                        goto 322
                  endif
                enddo

 322            if (idex.gt.0) then  ! W3 common value with W1
                    preorder (idex,ib3) = 0
                endif

            idex = 0
                ib4 = 4
                do k4=1,ntmp(ib4)
                  if (preorder (k4,ib4).eq.ival0) then
                        idex = k4
                        goto 323
                  endif
                enddo

 323            if (idex.gt.0) then  ! W4 common value with W1
                    preorder (idex,ib4) = 0
                endif
      enddo
                  

      ib2 = 2
        do k=1,ntmp(ib2)
                ival0 = preorder (k,ib2)

            if (ival0.gt.0) then
                 nord=nord + 1
                 torder (nord) = ival0

             print *, 'In band 2 but not in band 1: #',nord,' = ',ival0

                 idex = 0
                 ib3 = 3
                 do k3=1,ntmp(ib3)
                   if (preorder (k3,ib3).eq.ival0) then
                         idex = k3
                         goto 325
                   endif
                 enddo

 325             if (idex.gt.0) then  ! W3 common value with W2
                    preorder (idex,ib3) = 0
                 endif

                 idex = 0
                 ib4 = 4
                 do k4=1,ntmp(ib4)
                  if (preorder (k4,ib4).eq.ival0) then
                        idex = k4
                        goto 326
                  endif
                 enddo

 326             if (idex.gt.0) then  ! W4 common value with W2
                    preorder (idex,ib4) = 0
                 endif

            endif

        enddo

      ib3 = 3
        do k=1,ntmp(ib3)
                ival0 = preorder (k,ib3)

                if (ival0.gt.0) then
                 nord=nord + 1
                 torder (nord) = ival0

             print *, 'In band 3 but not in bands 1,2: #',
     1                    nord,' = ',ival0

                 idex = 0
                 ib4 = 4
                 do k4=1,ntmp(ib4)
                  if (preorder (k4,ib4).eq.ival0) then
                        idex = k4
                        goto 327
                  endif
                 enddo

 327             if (idex.gt.0) then  ! W4 common value with W3
                    preorder (idex,ib4) = 0
                 endif

                endif

        enddo

      ib4 = 4
        do k=1,ntmp(ib4)
                ival0 = preorder (k,ib4)

                if (ival0.gt.0) then
                 nord=nord + 1
                 torder (nord) = ival0

             print *, 'In band 4 but not in bands 1,2,3: #',
     1                    nord,' = ',ival0

                endif

        enddo

      nactive = nord

c      write (6,*) ntmp
      write (6,*) 'allocating ',nactive,'  frames'

      if(sort_subframes) then
        print *,'Sorting frames into ascending index order ...'
        call sorti2(nactive,torder)
      endif

      write (6,*) 'this is the list used for all bands: '
      do k=1,nactive
         order (k) = torder (k)
         write (6,'(a,i5,a,i5,a,a)') 'k=',k,', order=',order(k),
     1           ', fid=',basename(order(k))(1:9)
         do ib = 1,4
           preorder (k,ib) = torder (k)
         enddo
      enddo

      do k=nactive+1,nf
          do ib = 1,4
            preorder(k,ib) = 0
          enddo
      enddo

      endif ! .not.dosolo

c         if  (kreg.eq.7) call exit(0)
c         goto 8000

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc allocate & load the frames


c         if (nallocate.eq.0) then
c                  nsubxX = nsubx(ib)
c                  nsubyY = nsuby(ib)
c         endif

          if (nallocate.gt.0) then

                if (verbose) write (6,*) 'de-allocating'

              deallocate (pix_order)
                deallocate (Xpos)
                deallocate (Ypos)
                deallocate (SNRlist)
                deallocate (LBACK)
                deallocate (LSIG)
                deallocate (Lconf)
                deallocate (MAGSTD)
                deallocate (eMAGSTD)
                deallocate (framefluxes)
                deallocate (frameuncs)
            deallocate (framechis)
                deallocate (FLGSTD)
                deallocate (RAsublist)
                deallocate (DECsublist)
                deallocate (IDlist)
                deallocate (Table)

            deallocate (Array)
              deallocate (Unc)
              deallocate (iMask)
              deallocate (Mask)
            deallocate (wflagactive)
            deallocate (HDR)

            deallocate (xTRANS)
            deallocate (yTRANS)
            deallocate (varFLUX)
            deallocate (eVarFLUX)
            deallocate (wproJD)
            deallocate (varJD)
              deallocate (JD)
            deallocate (VarRchi2)

            endif

          if (verbose) write (6,*) 'allocate mini-arrays with size: ',nsubxX,nsubyY

          allocate (pix_order(nactive,4))
          if (.not.allocated(pix_order)) then
            write (6,*) 'allocation failed for pix_order'
            istat = 5
            call exit(istat)
          end if

          allocate (wflagactive(nactive,4))
            if (.not.allocated(wflagactive)) then
                write (6,*) 'allocation failed for wflagactive'
                istat = 5
                call exit(istat)
            end if

      
          allocate (HDR(nactive,4))
            if (.not.allocated(HDR)) then
              write (6,*) 'allocation failed for HDR'
              istat = 5
              call exit(istat)
            end if

          nactive_pix = sum(nactive_band(1:4))
            ! Add room for sentinels
            nactive_pix = nactive_pix + 5
          print *,'Pixel arrays being allocated with length = ',nactive_pix,
     1              ' = sum(',nactive_band,') + 5'      

          allocate (Array(nsubxX,nsubyY,nactive_pix))
            if (.not.allocated(Array)) then
              write (6,*) 'allocation failed for array'
              istat = 5
            call exit(istat)
            end if

          Array = -1.0e32

          allocate (Unc(nsubxX,nsubyY,nactive_pix))
            if (.not.allocated(Unc)) then
              write (6,*) 'allocation failed for unc'
              istat = 5
            call exit(istat)
            end if

          allocate (iMask(nsubxX,nsubyY,nactive_pix))   ! this is the MSK image
            if (.not.allocated(iMask)) then
              write (6,*) 'allocation failed for iMask'
              istat = 5
            call exit(istat)
            end if

!!! Mask is never assigned a value other than zero, so I'm
!!! taking it out of action -- TPC
c          allocate (Mask(nsubxX,nsubyY,nactive,4))
          allocate (Mask(1,1,1))
            if (.not.allocated(Mask)) then
              write (6,*) 'allocation failed for Mask'
              istat = 5
            call exit(istat)
            end if

          UNC = 0.
          iMASK = 0
          Mask = 0


         allocate (RAsublist(nsubsrcA))
           if (.not.allocated(RAsublist)) then
             write (6,*) 'allocation failed for RAsublist'
             istat = 5
           end if

           allocate (DECsublist(nsubsrcA))
           if (.not.allocated(DECsublist)) then
             write (6,*) 'allocation failed for DECsublist'
             istat = 5
           end if

           allocate (IDlist(nsubsrcA))
           if (.not.allocated(IDlist)) then
             write (6,*) 'allocation failed for IDlist'
             istat = 5
           end if


         allocate (Xpos(nsubsrcA,nactive,4))
           if (.not.allocated(Xpos)) then
            write (6,*) 'allocation failed for Xpos'
            istat = 5
           end if

      ! NEW
         Xpos = 0.

           allocate (Ypos(nsubsrcA,nactive,4))
           if (.not.allocated(Ypos)) then
             write (6,*) 'allocation failed for Ypos'
             istat = 5
           end if

      ! NEW
         Ypos = 0.

         allocate (SNRlist(nsubsrcA,nactive,4))
           if (.not.allocated(SNRlist)) then
             write (6,*) 'allocation failed for SNRlist'
             istat = 5
           end if

      ! NEW
           SNRlist = 0.
      

           allocate (LBACK(nsubsrcA,nactive,4))
           if (.not.allocated(LBACK)) then
             write (6,*) 'allocation failed for LBACK'
             istat = 5
           end if
         LBACK = -999.99

           allocate (LSIG(nsubsrcA,nactive,4))
           if (.not.allocated(LSIG)) then
             write (6,*) 'allocation failed for LSIG'
             istat = 5
           end if
         LSIG = 0.

           allocate (Lconf(nsubsrcA,nactive,4))
           if (.not.allocated(Lconf)) then
             write (6,*) 'allocation failed for Lconf'
             istat = 5
           end if
         Lconf = 0.

         allocate (MAGSTD(nsubsrcA,4))
           if (.not.allocated(MAGSTD)) then
             write (6,*) 'allocation failed for MAGSTD'
             istat = 5
           end if

        
           allocate (eMAGSTD(nsubsrcA,4))
           if (.not.allocated(eMAGSTD)) then
             write (6,*) 'allocation failed for eMAGSTD'
             istat = 5
           end if

           allocate (framefluxes(nsubsrcA,nactive,4))
           if (.not.allocated(framefluxes)) then
             write (6,*) 'allocation failed for framefluxes'
             istat = 5
           end if
           allocate (frameuncs(nsubsrcA,nactive,4))
           if (.not.allocated(frameuncs)) then
             write (6,*) 'allocation failed for frameuncs'
             istat = 5
           end if
         allocate (framechis(nsubsrcA,nactive,4))
           if (.not.allocated(framechis)) then
             write (6,*) 'allocation failed for framechis'
             istat = 5                                          
           end if                                         
         framechis = -1.e35
           framefluxes = -1.e35
           frameuncs = -1.e35

           allocate (FLGSTD(nsubsrcA,4))
           if (.not.allocated(FLGSTD)) then
             write (6,*) 'allocation failed for FLGSTD'
             istat = 5
           end if


       allocate (xTRANS(nactive,4))
         if (.not.allocated(xTRANS)) then
              write (6,*) 'allocation failed for xTRANS'
              istat = 5
         end if

         allocate (yTRANS(nactive,4))
         if (.not.allocated(yTRANS)) then
             write (6,*) 'allocation failed for yTRANS'
             istat = 5
         end if 


      allocate (varFLUX(nactive))
        if (.not.allocated(varFLUX)) then
           write (6,*) 'allocation failed for varFLUX'
           istat = 5
        end if
        allocate (eVarFLUX(nactive))
        if (.not.allocated(eVarFLUX)) then
           write (6,*) 'allocation failed for eVarFLUX'
           istat = 5
        end if

      allocate (VarRchi2(nactive))
        if (.not.allocated(VarRchi2)) then
           write (6,*) 'allocation failed for VarRchi2'
           istat = 5
        end if

        allocate (wproJD(nactive))
        if (.not.allocated(wproJD)) then
           write (6,*) 'allocation failed for wproJD'
           istat = 5
        end if
        allocate (varJD(nactive))
        if (.not.allocated(varJD)) then
           write (6,*) 'allocation failed for varJD'
           istat = 5
        end if
      allocate (JD(nactive,4))
        if (.not.allocated(JD)) then
          write (6,*) 'allocation failed for JD'
          istat = 5
        end if

      allocate (nBand(nactive_pix))                   ! JWF B80318
      if (.not.allocated(nBand)) then
        write (6,*) 'allocation failed for nBand'
        istat = 5
        call exit(istat)
      end if
      nBand = 0
  
      allocate (MJD(nactive_pix))                     ! JWF B80318
      if (.not.allocated(MJD)) then
        write (6,*) 'allocation failed for MJD'
        istat = 5
        call exit(istat)
      end if
      MJD = 1.0d0


      nallocate = nallocate + 1

      do Jsrc = 1, nsubsrc
        RAsublist(Jsrc) = tmpRAlist (Jsrc)
        DECsublist(Jsrc) = tmpDEClist (Jsrc)
c        IDlist(Jsrc) = Jsrc
          IDlist(Jsrc) = tmpIDlist(Jsrc)
      enddo

c         endif  ! igot ;  arrays allocated; src sublist made


ccc  load the pixel arrays

      if (verbose) write (6,*) 'loading mini-frames'

      wflagactive = 0

! New indexing of frames !!! TPC
        npixfr = 0
      npixfr_band(:) = 0
      pix_order(:,:) = -999999999

      do 8002 ib=1,4

	icountbad = 0   ! TJ 07Sep2017

        do Jord=1,nactive
            xTRANS(Jord, ib) = 0
                yTRANS(Jord, ib) = 0
            JD(Jord,ib) = 0.
        enddo


        if (any(wflag(1:nf,ib)==1)) then

          ! Insert a sentinel layer to help ensure we're indexing correctly.
          ! Leave at the initialized value of -1e32
        npixfr = npixfr + 1

        do Jord = 1,nactive
          Jfr = order (Jord)

          wflagactive (Jord,ib) = 0
	  icountbad = 0   ! TJ 07Sep2017


            ! Consider only frame previously recognized as active for this band/region
          if ((fram(Jfr,ib)(1:4).ne.'null') .and. (Jfr.gt.0) .and.
     1          fractive(jfr,ib).gt.0)  then

c            write (6,'(2i5,2x,a)') jfr,ib,fram(Jfr,ib)(1:72)

           wflagactive (Jord,ib) = wflag (Jfr,ib)

             ! Load based on band-dependent sequence - TPC
           npixfr_band(ib) = npixfr_band(ib) + 1
           npixfr = npixfr + 1
             ! Map from Jord=1:nactive (and order(jord)=1:nfr)) to position in pixel arrays
             pix_order(jord,ib) = npixfr


           do NUT = 1,3 ! TJ 21July2017  , don't need to set the mask, assume all zero's;   using INVVAR for mask

              if (NUT.eq.1) s0 = fram(Jfr,ib)
              if (NUT.eq.2) s0 = uncn(Jfr,ib)
              if (NUT.eq.3) s0 = msk(Jfr,ib)
              print *,'CatWISE - NUT, Jord, Jfr, ib:',NUT, Jord, Jfr, ib
              print *,'CatWISE - s0: ',s0(1:lnblnk(s0))

            if (NUT.eq.1) then
              nsx(ib) = 0
                nsy(ib) = 0

              xTRANS(Jord, ib) = 0
              yTRANS(Jord, ib) = 0
            endif


	    if (s0(1:4).eq.'null') then  ! TJ 07Sep2017
		Larray = 0.
	    else
              call readimage
     1             (nsizex,nsizey,lsize,Larray,iLarray,s0,
     1                  crval1,crval2,cdelt1,cdelt2,crot,crpix1,crpix2, tJD)
	    endif

            if (NUT.eq.1) then
            nsx(ib) = nsizex
            nsy(ib) = nsizey

            if (JD(Jord,ib).le.0.)  then
              JD (Jord,ib) = tJD
c             write (6,*)  Jord, Jfr, ib, JD (Jord,ib),'  ',s0(1:132)
            endif

            endif


            ic = 0
            do jj=1,nsizey
            do ii=1,nsizex
                  ic=ic+1
                  if(NUT .le. 2) then
                    Tarray(ii,jj, NUT) = Larray(ic)
                  else
                    ! Integer mask array
c                    iTarray(ii,jj) = iLarray(ic)
c	  	     iTarray(ii,jj) = int(Larray(ic))  ! unwise

		    !TJ   STD mask biz
		    iTarray(ii,jj) = 0  ! nominal
		    if (Larray(ic).le.0.) then
				iTarray(ii,jj) = 2  !  this should kill the pixel
		    endif

                  endif

                  if ((SPIT).and.(NUT.le.2))  then 
                   ! Spitzer only
                   !  convert SB to DN/s and multiply by the Loc-Dep Photom correction
                    if ( (MIPS).and.(ib.eq.4) ) then
                      Tarray(ii,jj, NUT) = Tarray(ii,jj, NUT) / fluxcon(5)
                    else if ( (MIPS).and.(ib.eq.3) ) then
                      Tarray(ii,jj, NUT) =  Tarray(ii,jj, NUT) / fluxcon(4)
                    else
                      Tarray(ii,jj, NUT) = Tarray(ii,jj, NUT) / fluxcon(ib)
                    endif
                  endif
                enddo
            enddo


c            write (6,*) ib,'  image loaded with size: ',nsizex,nsizey,nsx(ib),nsy(ib)

            if (NUT.eq.1) then

            Hfits = HDR0 (Jfr, ib)
            HDR(Jord,ib) = Hfits
            
            if (WCS(Jord).ge.0) call wcsclose(WCS(Jord))
            WCS(Jord) = -1
            call wcsinit(Hfits,WCS(Jord))


cccc  load the position arrays
            do  Jsrc = 1,nsubsrc

              Ra0 = RAsublist(Jsrc) * 1.d0
              Dec0 = Decsublist(Jsrc) * 1.d0

              ijo = 0
              Hfits = HDR0(Jfr,ib)

              call WhatPos (Jord,WCS(Jord),nsx(ib),nsy(ib),pscale(ib), Hfits, ra0, dec0, x0, y0, ijo, edgebuf) ! frame positions

              if ((ib.gt.1).and.(ijo.eq.1)) then
                if ( any ( wflag ( 1:nf,1 ) ==1 ) ) then
                  if ( (XPos (Jsrc,Jord,1).eq.0.).and.(YPos (Jsrc,Jord,1).eq.0.)) then
c                  ijo = 0    !  THJ 20Dec2017, this line kills Option 1b
                  endif
                endif
              endif

              if (ijo.eq.1) then

                XPos (Jsrc,Jord,ib) = x0
                YPos (Jsrc,Jord,ib) = y0

              else
                XPos (Jsrc,Jord,ib) = 0.
                YPos (Jsrc,Jord,ib) = 0.
              endif

            enddo            ! sources (jsrc)

            if (debug) write (6,*)  '%%%%%    X-Y position matrix loaded'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccq
c region center in frame space

c               debug = .true.

            ra0 = ra_reg * 1.d0
            dec0 = dec_reg * 1.d0

            call wcs2pix(WCS(Jord), ra0, dec0, xpix, ypix, offscl)

            xquad_frame = xpix * 1.0
            yquad_frame = ypix * 1.0

            if (debug) then
              write (6,*) 'quad center pix coords in frame space: ',Jord,xquad_frame, yquad_frame 
            endif

            iquad_frame = nint(xquad_frame)
            jquad_frame = nint(yquad_frame)
             
            if (nsy(ib).eq.nsuby(ib)) jquad_frame = j0_reg(ib)
            if (nsx(ib).eq.nsubx(ib)) iquad_frame = i0_reg(ib)

         
            if (debug) then
              write (6,*) ' '
              write (6,*)
     1              ' actual transformation ',
     1                  '(add these numbers to the new position xy): '
              write (6,*) iquad_frame - i0_reg(ib)
              write (6,*) jquad_frame - j0_reg(ib)
            endif

             xTRANS(Jord, ib) = iquad_frame - i0_reg(ib)
               yTRANS(Jord, ib) = jquad_frame - j0_reg(ib)


            endif ! nut=1

cc here we actually load the sub-images Array, UNC, iMask

            PixSum = 0.0d0
            do j=1,nsy(ib)
            jdel_yz = j - jquad_frame
            jdex = j0_reg(ib) + jdel_yz
            do 120 i=1,nsx(ib)
              idel_xz = i - iquad_frame
              idex = i0_reg(ib) + idel_xz
        
              if (idex.lt.1) goto 120
              if (jdex.lt.1) goto 120
              if (idex.gt.nsubx(ib)) goto 120
              if (jdex.gt.nsuby(ib)) goto 120

c                 if (NUT.eq.1) write (6,*) i,j,idex,jdex

              if (NUT.eq.1) Array(idex,jdex,npixfr) = Tarray(i,j, NUT)
              if (NUT.eq.2) then
                UNC(idex,jdex,npixfr) = Tarray(i,j, NUT)  * STDscale(ib)  ! THJ 22Jan2018 -- scale up the UNCs
                if (UNC(idex,jdex,npixfr) .gt. 0.0) then
                  if (STDfloor(ib)) then
                    if (UNC(idex,jdex,npixfr) .lt. STDbias(ib))
     +                  UNC(idex,jdex,npixfr)   =  STDbias(ib)
                  else
                    UNC(idex,jdex,npixfr) = UNC(idex,jdex,npixfr) + STDbias(ib)
                  end if
                end if
                PixSum = PixSum + UNC(idex,jdex,npixfr)
              end if

	          if (NUT.eq.3) then
	            iMASK(idex,jdex,npixfr) = iTarray(i,j)
		        if (iTarray(i,j).gt.0) icountbad=icountbad+1
              endif

 120        continue      ! i
            enddo            ! j
            
            if (NUT .eq. 2) then
              MJD(npixfr)   = tJD
              nBand(npixfr) = ib
              nCWchk        = npixfr
              s0            = uncn(Jfr,ib)
              print *
              print *,'CatWISE: uncertainty image ',s0(1:lnblnk(s0))
              print *,'CatWISE: npixfr, MJD(npixfr), nBand(npixfr):',
     +                 npixfr, MJD(npixfr), nBand(npixfr)
              print *,'CatWISE PixSum:', PixSum
              print *
            end if

          enddo            ! NUT = 1-3


	  write (6,*) '**** Badpix from STD mask:  ',ib,Jfr,icountbad
	  icountbad = 0


c vvvvvvvvvvvvvvv debug vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
          if (debug) then
            if ((ib.eq.1).and.(Jord.eq.92)) then
            ic = 0
            do j=1,nsuby(ib)
              do i=1,nsubx(ib)
                ic=ic+1
                larray(ic) = Array(i,j,pix_order(jord,ib))
                if (larray(ic).gt.0.) ohyea = 1
              enddo
            enddo           
            if (Jord.le.9) then
              write (fout,'(a,i1,a)') 'temp',Jord,'.fits'
            else if (Jord.le.99) then
              write (fout,'(a,i2,a)') 'temp',Jord,'.fits'
            else
              write (fout,'(a,i3,a)') 'temp',Jord,'.fits'
            endif
            ohyea = 1
            write (6,'(a)') fout(1:50)
            call wfits (nsubx(ib),nsuby(ib),lsize,larray,fout)
c               call exit(0)   
            endif
            if ((Jord.eq.1).and.(ib.eq.1)) then
c               if (ib.eq.9) then
            ic=0
            do j=1,nsy(ib)
              do i=1,nsx(ib)
                ic=ic+1
                larray(ic) = Tarray(i,j,1)
              enddo
            enddo
            fout = 'parent.fits'
            write (6,'(a)') fout(1:50)
            call wfits (nsx(ib),nsy(ib),lsize,larray,fout)
            do  Jsrc = 1,nsubsrc
              xx = XPos (Jsrc,Jord,ib) - xTRANS(Jord, ib)
              yy = YPos (Jsrc,Jord,ib) - yTRANS(Jord, ib)
c                write (67,'(2f10.5,4f9.2)') RAsublist(Jsrc),Decsublist(Jsrc),XPos (Jsrc,Jord,ib), YPos (Jsrc,Jord,ib), xx, yy
              write (22,'(2f10.5,4f9.2)') RAsublist(Jsrc),Decsublist(Jsrc),XPos (Jsrc,Jord,ib), YPos (Jsrc,Jord,ib), xx, yy
            enddo
c              call exit(0)
            endif            ! if ((Jord.eq.1).and.(ib.eq.1)) then
          endif            !!! debug
c ^^^^^^^^ debug ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        endif                  ! if ((fram(Jfr,ib)(1:4).ne.'null').and.(Jfr.gt.0))  then


      enddo ! Jord, ending loop in which we have ordered the read according to what was read in the previous grid

      endif                  ! any wflag

c	   write (6,*) '**** Badpix from STD mask:  ',ib,icountbad  ! TJ 07Sep


 8002      continue            ! ib


      call  ETime (Second, dd2)
        dd2 = second(1)
        dd = dd2 - dd1

      if(verbose) then
        write (6,*) 'total number of frames examined', nf
        write (6,*) 'total number of active/allocated frames ',nactive
        write (6,'(a,f10.3,a)') 'frame load DTime ',dd,' sec'
        print *
          print *, 'Pixel frame count = ',npixfr, '(no. allocated = ',nactive_pix,')'
          print *, 'Pixel frame counts per band = ',npixfr_band
        print *, 'Mapping from overall frame sequence to location in pixel arrays (fr=frpix):'
        nperrow = 15
          do ib=1,4
            print *, '  Band = ',ib
          nrows = nactive/nperrow
          if(mod(nactive,nperrow) .ne. 0) nrows = nrows + 1
            do i=1,nrows
              jbot = (i-1)*nperrow+1
              jtop = min(i*nperrow,nactive)
            nrow = jtop - (i-1)*nperrow
c              write(6,'(a,<nrow>(i4,a,i5,a,a9,a))') '    ',               ! JWF B60711; format repeats should
              write(6,'(a,(i4,a,i5,a,a9,a))') '    ',                      !             work without this
     1                (j,'=',pix_order(j,ib),':',basename(order(j)),', ', j=jbot,jtop)
            end do
          end do
      endif

        if(npixfr .gt. nactive_pix) then
          print *,'***ERROR: Read pixels for ',npixfr,' bandframes. This exceeds the ',
     1            nactive_pix,' bandframes allocated.'
          call exit(64)
        endif

!!  at this point, all of the frames are loaded; source positions are
!!  loaded: RAsublist, DECsublist, XPos, YPos
c
c----------------------- CatWISE-specific code changes - JWF B80316 --------------------
c                                                   Assumes only 1 PSF per band, 641x641
c                                                           Assumes images are 2048x2048
      print *,'Start of CatWISE-specific cryo check; no. files:', nCWchk
      do 12345 Jord = 1, nCWchk                     
        print *,'Checking frame JORD =', Jord,
     +          ' for cryo status; MJD =', MJD(Jord)
        if ((MJD(Jord) .gt. 2.0d0) .and.
     +      (MJD(Jord) .lt. MJDhibern8)) then
          if ((MJD(Jord) .lt. MJDendCryoPSF) .and. ! Peter's end-of-cryo-PSF epoch
     +        (nBand(Jord) .eq. 1) .and. .not.DidCryo) then    
            print *,'W1 cryo detected; Jord, MJD(Jord):', Jord, MJD(Jord)
            print *,'W1 PSFunc rescaled by', W1cryoPSFuncscale
            do 1234 jj = 1, 641
              do 1233 ii = 1, 641
                PSFuncs(ii,jj,1,1) = PSFuncs(ii,jj,1,1)*W1cryoPSFuncscale
1233          continue
1234        continue
            DidCryo = .true.
          end if
          PixSum = 0.0d0
          do 5432 j = 1, 2048
            do 5431 i = 1, 2048
              PixSum = PixSum + UNC(i,j,Jord)
              UNC(i,j,Jord) = UNC(i,j,Jord)*PrehibSTDScale(nBand(Jord))
5431        continue
5432      continue
          print *,'Pre-hibernation std image for frame no.',Jord,
     +            ' band', nBand(Jord),
     +            ' rescaled by', PrehibSTDScale(nBand(Jord))
          print *,'Previous PixSum =', PixSum
        end if
12345 continue
      print *,'End of CatWISE-specific cryo check'
c
c-----------------End of CatWISE-specific code changes - JWF B80316 --------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c we need to measure the local background and compute the flux before we call WPRO
c set the standard aperture parameters in  pixels
c NOTE -- this is now calibrated for COADD pixels

      if (verbose) then
        write (6,*) ' '
        write (6,*) 'initialize aperture sizes for coadd pixels'
      endif

      do ib=1,4

         Rstap (ib) = 0.
         RstapF (ib) = 0.
         Rstann (ib) = 0.
         Rstwid (ib) = 0.

c       if (wflag(1,ib) .eq. 1) then
       if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

         RR = Rapc(2)
         do Map = 2, Naper  ! these are the circular apertures
            RR = Rapc(Map)

c            if (wflag(1,ib) .eq. 1) then
              if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then
              Rcirc (Map,ib) = RR / cscale(ib)   ! units = pixels
              if ((ib.eq.4).and.(.not.SPIT)) Rcirc (Map,ib) = Rcirc (Map,ib) * 2.  !  increase the radius due to the bigger beam
              if ((MIPS).and.(ib.eq.4))  Rcirc (Map,ib) = Rcirc (Map,ib) * 2. 
            endif

            if (kreg_ix.le.1 .or. doreg.gt.0) then    ! TPC
             write (What,'(a,i1)') "Rap_",Map-1
               TYPE = "r"
               metval = Rcirc (Map,ib) * cscale(ib)
               comment = 'circular aperture radius in arcsec'
               iband = ib
               call MetWr (imeta, iband, What, TYPE, metval, comment)
            endif


         enddo


      ! standard aperture  & annulus
          Rstap (ib) = Raper / cscale(ib)   ! units = coadd pixels
          RstapF (ib) = Raper / pscale(ib)   ! units = frame pixels
          Rstann (ib) = Rannulus / pscale(ib)   ! units = frame pixels
          Rstwid (ib) = Rwidth / pscale(ib)   ! units = frame pixels

          if ((ib.eq.4).and.(.not.SPIT)) then ! 
            Rstap (ib) = Raper * 2. / cscale(ib)  !  increase the radius due to the bigger beam
            RstapF (ib) = Raper * 2. / pscale(ib)
          endif

          if ((MIPS).and.(ib.eq.4))  Rstap (ib) = Raper * 1.8045/ cscale(ib)  !   MIPS

          Rcirc (1,ib) = Rstap (ib)     ! standard aperture in coadd pixels

          if (kreg_ix.le.1 .or. doreg.gt.0) then ! TPC
              write (What,'(a,i1)') "Stap"
                TYPE = "r"
                metval = Raper
            if (ib.eq.4) metval = Raper * 2. 
                comment = 'standard aperture radius in arcsec'
                iband = ib
                call MetWr (imeta, iband, What, TYPE, metval, comment)
          endif
          
        endif

        if (SPIT) then
          TYPE = "r"
          metval = IRACapcor(ib)
          comment = 'IRAC aperture correction'
          if (MIPS) then
                if (ib.eq.4) comment = 'MIPS aperture correction'
          endif
          iband = ib
          call MetWr (imeta, iband, What, TYPE, metval, comment)

        CmosaicCorr(ib) = IRACapcor(ib)

         endif

      enddo


       if (verbose) then 
        write (6,*) '%%%%%%%%%%%%%%%%%%%%%%'

       write (6,*) '------------'
         write (6,*) 'Running Pre-WPRO'
         write (6,*) '------------'

      endif

         call  ETime (Second, dd100)
          dd100 = second(1)
          dd = dd100 - dd90

         write (6,'(a,f10.3,a)')' Ttest time ',dd,' sec'
c      call exit(0)

c  photometry; call the wrapper  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      dd1 = 0.
      call  ETime (Second, dd1)
      dd1 = second(1)  ! User time

      fmode = 0   ! fast mode; use to compute quick standard ap photometry


      if (verbose) write (6,*) 'call Phot_Wrapper '
        if (verbose) write (6,*) ' computing local backgrounds and simple aperture photometry for '
        if (verbose) write (6,*) nsubsrc,' sources and ',nactive,' frames'

c  notes: nsubsrcA

       call Phot_Wrapper (fmode,wflagactive,nsubsrcA,nsubsrc,
     1       nsubxX,nsubyY,nsubx,nsuby, nactive, nactive_pix,
     1       Array,MASK,iMask, Unc, pix_order, SNRlist,
     1       Xpos,Ypos,zero, xTRANS, yTRANS,
     1       RstapF, Rstann, Rstwid, BGmax,
     1       LBACK, LSIG, Lconf, MAGSTD, eMAGSTD, FLGSTD, NSUMmax,Rsatsub, pscale,
     1       NtotSrc, SPIT, MIPS, fbits, NbannAve, IOCname, doIOC)


c      write (6,*) NbannAve(2)
c      write (6,*) 'here ',LBACK (1616,18,1),LBACK (40282,18,1)

      call  ETime (Second, dd2)
      dd2 = second(1)

c      write (6,*) second(1),second(2),dd2

      dd = dd2 - dd1

       if (verbose) write (6,*) 'total number of sources processed for region ',kreg,
     1                           ' = ',nsubsrc
c      if (verbose) rtime = (NtotSrc*1.) / dd
      if (verbose) write (6,'(a,f10.3,a)')' Phot_Wrapper DTime ',dd,' sec'
c      if (verbose) write (6,*) 'sources per second ',rtime


ccccccccc
c  allocate the big output table for Wpro

      NmaxSrc = nsubsrcA

      if (debug) write (6,*) 'allocate big table ',NmaxSrc

c      allocate (Table(NmaxSrc,32))
c      allocate (Table(NmaxSrc,59))   ! JWF B30221
      allocate (Table(NmaxSrc,65))   ! JWF B31209
        if (.not.allocated(Table)) then
          write (6,*) 'allocation failed for table'
          istat = 5
        call exit (istat)
        end if

      Table = -99.99   ! initialize matrix with null values  !  TJ 23Aug2016

c      ! TEMP
c      do Jsrc = 1,nsubsrc
c            Ra0 = RAsublist(Jsrc) * 1.d0
c            Dec0 = Decsublist(Jsrc) * 1.d0
c            write (67,'(2f10.5)') Ra0,Dec0
c      enddo
c      call exit(0)

      if (debug) then

        do Jsrc = 1,nsubsrc

          Ra0 = RAsublist(Jsrc) * 1.d0
          Dec0 = Decsublist(Jsrc) * 1.d0

c                 if ( (RA0.gt.180.1435).and.(RA0.lt.180.1437).and.(Dec0.gt.0.995).and.(Dec0.lt.0.996)) then
          do ib=1,1
            do kk=1,nactive
            xx = Xpos (Jsrc,kk,ib)- xTRANS(kk, ib)
            yy = Ypos (Jsrc,kk,ib)- yTRANS(kk, ib)
            
c               if (MAGSTD(Jsrc,ib).lt.20.) then
            write (29,'(3i6,4f8.1,2f12.3,f9.3)') Jsrc,ib,kk,xx,yy,Xpos (Jsrc,kk,ib),Ypos (Jsrc,kk,ib),
     1            LBACK (Jsrc,kk,ib),LSIG(Jsrc,kk,ib),MAGSTD(Jsrc,ib)
c            endif

            enddo
          enddo

c            endif

        enddo

        goto 887

        iunit=59
        do ib=1,4
          iunit=iunit+1
          write (iunit,*) 'id, ra, dec, x, y, LB, sig, SNRpeak, mag, dmag, flag'

c      if (wflag(1,ib).eq.1) then
          if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

            do Jsrc = 1,nsubsrc

      
            xx = Xpos (Jsrc,5,ib)- xTRANS(5, ib)
            yy = Ypos (Jsrc,5,ib)- yTRANS(5, ib)

            if ((xx.gt.0.).and.(yy.gt.0.)) then
              
              write (iunit,66) RAsublist(Jsrc), DECsublist(Jsrc),
     1                     ( Xpos (Jsrc,Jfr,ib),Ypos (Jsrc,Jfr,ib),
     1                 xx, yy,
     1                     LBACK (Jsrc,Jfr,ib),LSIG(Jsrc,Jfr,ib),
     1                 SNRlist(Jsrc,Jfr,ib),
     1                 MAGSTD(Jsrc,ib),eMAGSTD(Jsrc,ib),FLGSTD(Jsrc,ib),Jfr=5,5)
 66              format (2f11.5, 10 (4f9.2, f9.4,f8.4, f8.1, f8.3,f6.3,i3)  )

            endif
            
            enddo ! jsrc
            close (iunit)

          endif ! any wflag

        enddo ! ib=1,4

      endif ! debug



cccccccccccccccccccccccccccccccccccccccccccccc  WPRO

 887    dd1 = 0.
        call  ETime (Second, dd1)
        dd1 = second(1)


      Nwpro = nsubsrc

      Nwpro_out = 0
      Nwpro_adb = 0

      if (debug) then
      write (6,*) nsubxX,nsubyY
      write (6,*) nsubx
      write (6,*) nsuby
      write (6,*) xTRANS(1, 1),yTRANS(1, 1)
      write (6,*) nx,ny
       write (6,*) nsubsrcA, nsubsrc
      write (6,*) nsx
      write (6,*) nsy

      do Jsrc = 1,nsubsrc
            write (6,'(i6,2f9.2,3f10.3,2f8.2)') Jsrc,Xpos (Jsrc,1,1),Ypos (Jsrc,1,1),
     1       LBACK (Jsrc,1,1), LSIG (Jsrc,1,1), Lconf (Jsrc,1,1), MAGSTD(Jsrc,1), eMAGSTD(Jsrc,1)
      enddo

      endif

      if (verbose) write (6,*) 'call WPRO'

        SatNumSub = 0


      call wpro_v6 (fbits,nsubxX,nsubyY,nx,ny,nsubx,nsuby,nsx,nsy,
     1      Array,Unc,iMask,pix_order,nactive,nactive_pix,wflagactive,
     1      PSFs,PSFuncs,MapPSFs,npsize,npnmax,npsf,pxcent,pycent,ppix,
     1      fwhm,nsubsrcA, nsubsrc, IDlist, HDR, RAsublist,DECsublist,Xpos,Ypos, xTRANS, yTRANS, JD, Rstann, Rstwid,
     1      LBACK,LSIG, Lconf, MAGSTD, eMAGSTD,FLGSTD, NSUMmax,Rsatsub, SatNumSub, zero, adb_nmax, Table,
     1      framefluxes,frameuncs,framechis, Nwpro_out, Nwpro_adb, istatus, SPIT, MIPS,
     +      nf, basename, order)
c
        N_adb = N_adb + Nwpro_adb
c
c      do j=1,nwpro_out
c         write (75,'(2i6,2f12.6,2f13.4,7f12.4)') kreg,j,Table(j,1:2),
c    *            Table(j,6:7),Table(j,10:11),Table(j,18),Table(j,29:32)
c      enddo

c      do kk=1,4
c      write (6,*) (SatNumSub(kk,ib),ib=1,4)
c      enddo


          if (istatus /= 0) then      
             write (6,*) 'ERROR -- WPRO fails; status non-zero  ',istatus
             call exit(64)
        endif

c      write (6,*) ' '
c      write (6,*) table (Jsrc,1),table(Jsrc,2)
c      write (6,*) table (Jsrc,8),table(Jsrc,12)
c      write (6,*) table (Jsrc,16), table(Jsrc,18),table(Jsrc,19), table(Jsrc,20)
c      write (6,*) table (Jsrc,25),table (Jsrc,26),table (Jsrc,27),table (Jsrc,29)
c      call exit(0)

c       write (6,*) table (1,6),table (1,7),table (1,8),table (1,9)
c       write (6,*) table (1,10),table (1,11),table (1,12),table (1,13)
c       write (6,*) table (1,25),table (1,26),table (1,27),table (1,28)

c      write (6,*) ' '
c      write (6,*) table (2,6),table (2,7),table (2,8),table (2,9)
c       write (6,*) table (2,10),table (2,11),table (2,12),table (2,13)


!-----------------------------------------------------------------------------
! Do profile-fitting photometry on WISE images.  Results are output in the
! form of an array named "Table", arranged as follows:
!               Table(n,1)      =       RA [deg] (single precision version)
!               Table(n,2)      =       Dec [deg] (single precision version)
!               Table(n,3)      =       uncertainty in RA [arcsec]
!               Table(n,4)      =       uncertainty in Dec [arcsec]
!               Table(n,5)      =       cross term (sigxy) in RA,Dec [arcsec]
!               Table(n,6:9)    =       fluxes in bands 1,...4 [dn]
!               Table(n,10:13)  =       uncert. in fluxes in bands 1,...4 [dn]
!               Table(n,14:17)  =       reduced chi squared in bands 1,...4
!               Table(n,18)     =       overall reduced chi squared
!               Table(n,19)     =       blend number
!               Table(n,20)     =       number of actively deblended components
!               Table(n,21:24)  =       fraction of pixels affected by latents
!                                       in bands 1,...4
!               Table(n,25:28)  =       fraction of saturated pixels
!                                       in bands 1,...4
!             Table(n,29)      =      Proper motion in RA [arcsec/yr]
!             Table(n,30)      =      Proper motion in Dec [arcsec/yr]
!            Table(n,31)      =      Uncertainty in RA proper motion ["/yr]
!            Table(n,32)      =      Uncertainty in Dec proper motion ["/yr]
!            Table(n,33)      =      RA proper motion to MJD0 from JD0 [RA deg] ! JWF B21222
!            Table(n,34)      =      Dec proper motion to MJD0 from JD0 [Dec deg] ! JWF B21222
!  (see wpro_v6 for full list)

      write (6,*) 

 4777      if (verbose) write (6,*) 'wpro complete'


      call  ETime (Second, dd2)
      dd2 = second(1)

        dd=dd2-dd1

       if (verbose)  write (6,*)
     +      'total number of sources processed by WPRO:',Nwpro_out
c      if (verbose)  rtime = nactive * (Nwpro_out*1.) / dd
       if (verbose)  write (6,'(a,f10.3,a)')' DTime ',dd,' sec'
c      if (verbose)  write (6,*) 'sources per second / frame ',rtime


      if (debug) then


      do j=1,nwpro_out

        if (IzBad(Table (j,3))) Table (j,3) = 999.
        if (IzBad(Table (j,4))) Table (j,4) = 999.
        Table (j,3) = min (Table (j,3), 999.)
        Table (j,4) = min (Table (j,4), 999.)

        Ra0 = table (j,1)
          Dec0 = table (j,2)

        if ( (RA0.gt.180.733).and.(RA0.lt.180.74).and.(Dec0.gt.-0.94).and.(Dec0.lt.-0.92)) then
            write (6,'(i6,2f10.5,20f9.3)') J,ra0,dec0, (LBACK (J,1,ib),ib=1,4)

          write (6,*)  Table (J,5+4), Table (J,9+4), Table (J,13+4)
          call exit(0)
        endif


c         if ( (RA0.gt.180.1435).and.(RA0.lt.180.1437).and.(Dec0.gt.0.995).and.(Dec0.lt.0.996)) then


c          write (57,'(i6,4f12.5,2f12.1,5f10.1,f5.0)') j,Table (j,1),Table (j,2),  Table (j,3),Table (j,4),
c     1       (Table (j,ii),ii=6,6), (Table (j,ii),ii=10,10)
c     1         ,(Table (j,ii),ii=14,17),
c     1     (Table (j,ii),ii=18,19)

c           write (58,*) Table (j,3),Table (j,4)

c         endif

      enddo

      endif ! debug


       if (nwpro_out.le.0) then
c                 write (6,*) 'ERROR -- WPRO fails, no sources extracted'
             write (6,*) 'WARNING -- no sources extracted ', Nwpro_out
                 call exit(64)
          endif


c  write to meta
      if (kreg.lt.10) then
            write (What,'(a,i1)') "Nwpro_out_",kreg
      else if (kreg.lt.100) then
                write (What,'(a,i2)') "Nwpro_out_",kreg
      else 
            write (What,'(a,i3)') "Nwpro_out_",kreg
      endif


        TYPE = "i"
        metval = Nwpro_out*1.
      write (comment,'(a,i2)') 'number of WPRO sources extracted, within region: ',kreg
      iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)

      if (kreg.lt.10) then
                write (What,'(a,i1)') "Nwpro_adb_",kreg
        else if (kreg.lt.100) then
                write (What,'(a,i2)') "Nwpro_adb_",kreg
        else
                write (What,'(a,i3)') "Nwpro_adb_",kreg
        endif
        TYPE = "i"
        metval = Nwpro_adb*1.
        comment = 'number of WPRO actively deblended sources'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)

      Nwpro =  Nwpro_out

c       Dup-checking search radius
      FW = 6.                  ! W1 asecs
      if (SPIT) FW = 2.
      del_pix = nint((FW / 2.5) / cscale(cband0)) ! pixels
      dmin2 = del_pix**2
      idel_pix = nint(del_pix)

        !! TPC debug
      write (6,*) 'Border and dup checking of ',Nwpro_out,' sources for region ',kreg,' (#',kreg_ix,'):'
      write (6,*) '  xmin/maxpix=',il_0,ih_0, '; ymin/maxpix=',jl_0,jh_0
      write (6,*) '  nx/yCoadd=',nxAWAIC,nyAWAIC
      write (6,*) '  searchRpix=',del_pix,', as int=',idel_pix

      do 600 Jsrc=1,Nwpro_out

        Ra0 = RAsublist(Jsrc)
        Dec0 = DECsublist(Jsrc)

        if(Ra0.ne.Ra0 .or. Dec0.ne.Dec0) then
          ! NaN position detected. Reject and keep count.
            RAsublist(Jsrc) = -9.
            DECsublist(Jsrc) = 0.
          nnanposrej = nnanposrej + 1
            goto 600
        endif

c         print *,'Doing source',jsrc,'; RA, Dec:',ra0,dec0 ! JWF dbg

        if (clevel(1:2).eq.'1b') then
            call WhatPos(1,WCSc(4),ncx(4),ncy(4), cscale(4), Hfits, ra0, dec0, x0, y0, igo, edgebuf)
        else
           call WhatPos(1,WCSc(cband0),ncx(cband0),ncy(cband0),cscale(cband0),
     +                     Hfits, ra0, dec0, x0, y0, igo, edgebuf)
        endif

        ico = nint(x0)
          jco = nint(y0)
c         print *,'(ico,jco):', ico,jco ! JWF dbg

c apply the hard limits
c         if ((ico.le.1).or.(jco.le.1).or.(ico.ge.nxAWAIC).or.(jco.ge.nyAWAIC)) then


c      if (Jsrc.eq.1) then
c        write (6,*) x0,y0
c        write (6,*) ico,jco
c        call exit (0)
c      endif

        if ((ico.le.il_0).or.(jco.le.jl_0).or.(ico.ge.nxAWAIC).or.(jco.ge.nyAWAIC).or.
     1       (ico.ge.ih_0).or.(jco.ge.jh_0)) then
                if (verbose .and. (ireg*jreg .eq. 1)) then
                  print *,
     +              'Removing off-border source with (RA,Dec,ix,iy):',
     +            RAsublist(Jsrc), DECsublist(Jsrc), ico, jco
                  nBorderViolators = nBorderViolators + 1
                end if
            if(mod(nbordrej,500) .eq. 0) then
              print *,'Border rejection (#',nbordrej,') of source #',Jsrc,'/',Nwpro_out,':'
              print *,'  RA,Dec,ico,jco = ',RAsublist(Jsrc), DECsublist(Jsrc), ico, jco
              print *,'  Border xlo,hi, ylo,hi = ',il_0,ih_0,jl_0,jh_0
            endif
            RAsublist(Jsrc) = -9.
            DECsublist(Jsrc) = 0.
            nbordrej = nbordrej + 1
            goto 600
        endif


c check for Dupes
c see if any sources within 1/3 of FWHM

        if ((nf.gt.1).and.(ireg*jreg.gt.1).and.doreg.eq.0) then  ! only do this for multiple frames & divisions

         il = ico - idel_pix
         il = max (il,1)
         il = min (il,nxAWAIC)

         ih = ico + idel_pix
         ih = max (ih,1)
         ih = min (ih,nxAWAIC)

         jl = jco - idel_pix
         jl = max (jl,1)
           jl = min (jl,nyAWAIC)

           jh = jco + idel_pix
         jh = max (jh,1)
           jh = min (jh,nyAWAIC)

         do jj = jl,jh
           dy2 = (jj-jco)**2
           do ii= il,ih
             dx2 = (ii-ico)**2
             dr2 = dy2+dx2

c             if ((Mdupe (ii,jj).ge.1).and.(dr2.le.dmin2)) then
               ! Reject dups if the matching source is not from the current region
             if (Mdupe(ii,jj).gt.0 .and. dr2.le.dmin2 .and. Mdupe(ii,jj).ne.kreg) then ! TPC
                 ! dupe
             RAsublist(Jsrc) = -9. ! kill the source
             DECsublist(Jsrc) = 0.
             nduprej = nduprej + 1
c               write (86,'(a,2i4,2f11.5,2i6,1x,2i6,f9.2)') 'dupe ',ireg,jreg,ra0,dec0,ico,jco,ii,jj,sqrt(dr2) ! TEMP
                 ! Exit now that it's rejected -- TPC
c             goto 610 
             goto 600
             endif
           enddo
         enddo

 610         continue

c         Mdupe (ico,jco) = Mdupe (ico,jco) + 1
           Mdupe (ico,jco) = kreg   ! TPC

       endif                  ! (nf.gt.1).and.(ireg*jreg.gt.1).and.doreg.eq.0

       if (RAsublist(Jsrc).ge.0.) then

          ntable_0 = ntable_0 + 1

          RAlist_0 (ntable_0) = RAsublist(Jsrc)
          DEClist_0 (ntable_0) = DECsublist(Jsrc)
          IDlist_0 (ntable_0) = IDlist(Jsrc)
c           print *,'keeping this one' ! JWF dbg
c          do L=1,30
c          do L=1,32    !   JWF B21015
c          do L=1,59    !   JWF B30221
          do L=1,65    !   JWF B30221
            Table_0 (ntable_0,L) = Table (Jsrc, L)
          enddo

          do ib=1,4
              SatNum (ntable_0,ib) = SatNumSub (Jsrc,ib)
          enddo

c===================== start of code added by CJG B30213 =========================

          if(mepwrite) then

c =============== TPC
c              Count just good frames
             ngood = 0
             do Jfr = 1, nactive
                  nxy = 0
                  do ib=1, 4
                    if(xpos(jsrc,jfr,ib).gt.0 .and. ypos(jsrc,jfr,ib).gt.0) then
                        nxy = nxy + 1
                    end if
                  end do
                  if(nxy .gt. 0) ngood = ngood + 1
               end do

               if(ngood+1 .gt. nfmep) then ! TPC
             nfmep_old = nfmep
                 nfmep = int(ngood*1.5 + 1)
             if(nfmep .gt. nf) nfmep = nf
             print *,'%%% Reallocating MEP tables from ', nfmep_old,
     +                   'to ', nfmep
             allocate (mepTable_tmp(nsummax,20,nfmep),
     +                     namTable_tmp(nsummax,nfmep),
     +                     jdTable_tmp(nsummax,nfmep),
     +                     stat = istat)
                 if (istat .ne. 0) then
                   print *,'ERROR: cannot allocate temp MEP arrays; ',
     +                     'Stat =',iStat,', nfmep = ',nfmep
                   call exit(64)
                 end if
             meptable_tmp = -99
             namtable_tmp = '99999z999'
             jdtable_tmp = -99.d19
             meptable_tmp(:,:,1:nfmep_old) = meptable
             namtable_tmp(:,  1:nfmep_old) = namtable
             jdtable_tmp (:,  1:nfmep_old) = jdtable
c                Mov_alloc is a built-in function. It ...
c                . Deallocates then re-allocates meptable, namtable, jdtable.
c                . Assigns *_tmp space to meptable, ...
c                . Deallocates *_tmp, ...
             call move_alloc(meptable_tmp,meptable)
             call move_alloc(namtable_tmp,namtable)
             call move_alloc(jdtable_tmp,jdtable)
               end if
c ============== TPC

               nmep = 0
             do 901 Jfr = 1, nactive ! frames

c                 Skip bad frames to save memory
                  nxy = 0                    ! TPC  vvvvvvvvvvv
                  do ib=1, 4
                    if(xpos(jsrc,jfr,ib).gt.0 .and. ypos(jsrc,jfr,ib).gt.0) then
                        nxy = nxy + 1
                    end if
                  end do
                  if(nxy .eq. 0) goto 901    ! TPC ^^^^^^^^^^^^

                  nmep = nmep + 1
              
              namtable(ntable_0,nmep) = basename(order(jfr))(1:9)
              jdtable(ntable_0,nmep) = JD(jfr,1)
              
              do ib=1,4
                 ibloc=1+(ib-1)*2
                 meptable(ntable_0,ibloc,nmep) = xpos(jsrc,jfr,ib)
                 meptable(ntable_0,ibloc+1,nmep) = ypos(jsrc,jfr,ib)
              enddo
              
              do ib=1,4
                 ibloc=9+(ib-1)*3
                 meptable(ntable_0,ibloc,nmep) = framefluxes(jsrc,jfr,ib)
                 meptable(ntable_0,ibloc+1,nmep) = frameuncs(jsrc,jfr,ib)
                 meptable(ntable_0,ibloc+2,nmep) = framechis(jsrc,jfr,ib)
              enddo
              
 901             enddo

          endif  ! mepwrite
c===================== end of code added by CJG B30213 =========================

        else
! rejected for some reason elsewhere; dups and out-of-region sources already skipped
          print *,'!!! Source ',jsrc,' in region ',kreg,' rejected for mysterious reasons!' ! TPC
          goto 600
        endif  ! rasublist(jsrc) >= 0

      ib=1

c  frame-fluxes and N/M statistics ;  N of M;  N_M
c      write (6,*) 'variability stats'

       do 701 ib = 1,4     ! band
c          print *,'doing band',ib ! JWF dbg
           medsky (ntable_0,ib) = -999.999
           medskysig(ntable_0,ib) = 0.
           formalsig(ntable_0,ib) = 0.
c================== JWF B31203: jump var stuff this if SingleFrame = T =================
             if (SingleFrame) go to 697
           N_Mstat (ntable_0,ib) = 0
           M_Mstat (ntable_0,ib) = 0
c           slope (ntable_0,ib) = 0.
c           eslope (ntable_0,ib) = 0.
c           slopechi2 (ntable_0,ib) = 0.

           mLogQ (ntable_0,ib) = -1.0
c           DeltaMag (ntable_0,ib) = 99.
           ks (ntable_0,ib) = 99.
c           Fmin (ntable_0,ib) = -99999.999
c           Fmax (ntable_0,ib) = -99999.999
           Rho12(ntable_0) = -999
           Rho23(ntable_0) = -999
           Rho34(ntable_0) = -999
           P12(ntable_0) = -999
           P23(ntable_0) = -999
           P34(ntable_0) = -999

           AVE_mJD(ntable_0,ib)  = -99.
             mJD_min(ntable_0,ib) = -99.
             mJD_max(ntable_0,ib) = -99.

c      write (6,*) 'debug ',ib

c           if (wflag(1,ib).eq.0) goto 701

           M_M_wpro = 0
           N_M_wpro = 0
           M_M_var  = 0 ! For use in the variability calculations and array indicies
697           nmer = 0
           do 702 Jfr = 1, nactive   ! frames
c               print *,'doing frame',jfr ! JWF dbg
c               if (wflagactive(Jfr,ib).eq.0) print *,'tossing on wflagactive(Jfr,ib).eq.0' ! JWF dbg
            if (wflagactive(Jfr,ib).eq.0) goto 702

            igo = 0
            xx =  XPos (Jsrc,Jfr,ib)
            yy =  YPos (Jsrc,Jfr,ib)
            call WhatLoc (nsx(ib),nsy(ib),pscale(ib), xx, yy, igo, edgebuf)   ! check for edge violation

c      if (jsrc.eq.1797) then
c      write (86,*) jsrc,ib,jfr,xx,yy,igo
c      print *,' jsrc,ib,jfr,xx,yy,igo:',jsrc,ib,jfr,xx,yy,igo ! JWF dbg
c      endif

            if (igo.eq.0) goto 702


                xx =  XPos (Jsrc,Jfr,ib)-xTRANS(Jfr, ib)  ! x frame position
                yy =  YPos (Jsrc,Jfr,ib)-yTRANS(Jfr, ib)  ! y frame position
            if (xx.lt.1) igo = 0
            if (yy.lt.1) igo = 0
            if (xx.gt.nsubx(ib)) igo = 0
            if (yy.gt.nsuby(ib)) igo = 0
            if ((XPos (Jsrc,Jfr,ib).eq.0.).and.(YPos (Jsrc,Jfr,ib).eq.0.)) igo = 0

c                call WhatLoc (nsubx(ib),nsuby(ib),pscale(ib), xx, yy, igo, edgebuf)   ! check for edge violation
c            if ((XPos (Jsrc,Jfr,ib).eq.0.).and.(YPos (Jsrc,Jfr,ib).eq.0.)) igo = 0


            if (igo.eq.0) goto 702

c                 print *,'WPHotpm: (Jsrc,Jfr,ib) =',Jsrc,Jfr,ib                ! JWF B30308
c                 print *,'framefluxes(Jsrc,Jfr,ib) =',framefluxes(Jsrc,Jfr,ib) ! JWF B30308
c                 print *,'  frameuncs(Jsrc,Jfr,ib) =',frameuncs(Jsrc,Jfr,ib)   ! JWF B30308
c                 print *,'  framechis(Jsrc,Jfr,ib) =',framechis(Jsrc,Jfr,ib)   ! JWF B30308

                if (SingleFrame) go to 698

                if ((framefluxes (Jsrc,Jfr,ib).gt.-999.).and.
     1              (frameuncs(Jsrc,Jfr,ib).gt.0.)) then
                                                         ! JWF B30614: Reject W1 & W2 saturated
c N/M counts for output
              M_M_wpro = M_M_wpro + 1
                  WproJD(M_M_wpro) = JD(Jfr,ib)
              SNR = framefluxes(Jsrc,Jfr,ib) / frameuncs(Jsrc,Jfr,ib)
              if (SNR.ge.3.) N_M_wpro = N_M_wpro + 1

c Save variability-related data. Filter out Post-cryo saturated frames.
                   NotSat = .true.                       ! fluxes in post-cryo
                   if ((ib .lt. 3) .and. (JD(Jfr,ib) .gt. SatMJD))
     +                  NotSat = framefluxes(Jsrc,Jfr,ib) .lt. SatFlux(ib)
                   if (NotSat) then
                        M_M_var = M_M_var + 1
                        VarFLUX(M_M_var) = framefluxes(Jsrc,Jfr,ib)
                        eVarFLUX(M_M_var) = frameuncs(Jsrc,Jfr,ib)
                        varJD(M_M_var) = JD(Jfr,ib)
                        VarRchi2(M_M_var) = framechis(Jsrc,Jfr,ib)
                   end if
                end if

698            if ( LBACK (Jsrc,Jfr,ib).gt.-999.) then

                  nmer=nmer+1
                        ztmp(nmer) = LBACK(Jsrc,Jfr,ib)
                        eztmp(nmer)  = LSIG(Jsrc,Jfr,ib)
                        peztmp(nmer) = Lconf(Jsrc,Jfr,ib)

c      if (ib.eq.1) write (6,*) 'sky ',LBACK(Jsrc,Jfr,ib),LSIG(Jsrc,Jfr,ib),Lconf(Jsrc,Jfr,ib)

            endif


 702           continue  ! Jfr

             if (SingleFrame) go to 699

           N_Mstat(ntable_0,ib) = N_M_wpro
             M_Mstat(ntable_0,ib) = M_M_wpro


c       if (ib.eq.1) write (6,*) 'here ',N_M_wpro,M_M_wpro  !TEMP

           if (.not.SPIT) then

c      write (6,*) ' '
c            write (6,*) M_M_var
c             do jj=1,M_M_var
c             write (6,'(2i4,3e13.5)') ib,jj,VarFLUX(jj),eVarFLUX(jj),VarRchi2(jj)  ! TEMP
c             enddo

            mmLogQ = -1.0
            nnDF = 0
            zFmin = -9999.999
            zFmax = -9999.999
            mLogQ(ntable_0,ib) =  mmLogQ
c            Fmin (ntable_0,ib) = zFmin
c            Fmax (ntable_0,ib) = zFmax
c            DeltaMag (ntable_0,ib) = 99.
            nDF (ntable_0,ib) = 0
c            dDeltaMag = 99.
c            ChiFac = 3.0  ! TEMP , this should be in the namelist

            if (M_M_var.gt.0.) then

            call ChkVar(VarFLUX, eVarFLUX,VarRchi2, varJD,  M_M_var,
     1                      ib,ChiFac, mmLogQ, dks, nnDF)

            endif

            if (nnDF.gt.0) then
              nDF (ntable_0,ib) = nnDF

              mLogQ(ntable_0,ib) =  mmLogQ
c              Fmin (ntable_0,ib) = zFmin
c              Fmax (ntable_0,ib) = zFmax

c              if ((dDeltaMag.gt.0.).and.(dDeltaMag.lt.90.)) then
c                  DeltaMag (ntable_0,ib) = dDeltaMag
c              endif

              ks(ntable_0,ib) = dks                         ! CJG B30319
c              write(6,*)'ks,dks: ',ks(ntable_0,ib),dks

c              if ( (zFmin.gt.0.).and.(zFmax.gt.0.) ) then
c                  DeltaMag (ntable_0,ib) = 2.5*log10(zFmax/zFmin)
c              endif

            endif

c             write (6,*) 'mLogQ, Fmin, Fmax, nDF: ',ib,mLogQ(ntable_0,ib),Fmin(ntable_0,ib),nDF

c       write (6,*) M_M_var
c      do kk=1,M_M_var
c        write (6,*) kk,varJD(kk),VarFLUX(kk)
c      enddo
             call varworks (ntable_0, Nactive,ib, M_M_var,VarFLUX,eVarFLUX,zero(ib),mergetype,varJD,
     1    NSUMmax,WPROmagR)

c Computation of JD min/max/mean for all frames used by wpro and for variability
            sumJD = 0.
            sumJDvar = 0
            if (M_M_Wpro.gt.0) then

               mJD_min(ntable_0,ib) = WproJD(1)
               mJD_max(ntable_0,ib) = WproJD(1)

               if (M_M_var.gt.0) then
                 mJD_min_var(ntable_0,ib) = varJD(1)
                 mJD_max_var(ntable_0,ib) = varJD(1)
               endif

               do kk=1,M_M_Wpro
                  sumJD =  sumJD + WproJD(kk)
                  if (WproJD(kk) .lt. mJD_min(ntable_0,ib)) mJD_min(ntable_0,ib) = WproJD(kk)
                  if (WproJD(kk) .gt. mJD_max(ntable_0,ib)) mJD_max(ntable_0,ib) = WproJD(kk)
                  if(kk .le. M_M_var) then
                    sumJDvar =  sumJDvar + varJD(kk)
                    if (varJD(kk) .lt. mJD_min_var(ntable_0,ib)) mJD_min_var(ntable_0,ib) = varJD(kk)
                    if (varJD(kk) .gt. mJD_max_var(ntable_0,ib)) mJD_max_var(ntable_0,ib) = varJD(kk)
                        endif
               enddo

               AVE_mJD(ntable_0,ib) = sumJD / (M_M_wpro * 1.d0)
               if(M_M_var .gt. 0) AVE_mJD_var(ntable_0,ib) = sumJDvar / (M_M_var * 1.d0)

c !!! The *mjd_var arrays aren't currently used. Eventually they should be passed to writeTABLE.

             endif


c            write (6,*) ib,AVE_mJD(ntable_0,ib)
c            write (6,*) mJD_min(ntable_0,ib),mJD_max(ntable_0,ib)


c              do jj=1,M_M_var
c                 write (6,'(i4,2e13.5)') jj,VarFLUX(jj),eVarFLUX(jj)
c               enddo
c      write (6,*) ' '
           endif


! SKY STATS
699            if (nmer.gt.1) then

            call MDIAN3(ztmp,eztmp,peztmp,Nmer,medsky(ntable_0,ib), medskysig(ntable_0,ib), formalsig(ntable_0,ib))

            else if (nmer.eq.1) then

                medsky(ntable_0,ib) = ztmp(1)
                medskysig(ntable_0,ib) = eztmp(1)
            formalsig(ntable_0,ib) = peztmp(1)

            endif

 701        continue  ! ib

          if (SingleFrame) go to 600

cccc correlation coeffificients;  NOTE:  the P values below are *really* Q values in the output -- see the SIS
      Rho12(ntable_0) = -999
      Rho23(ntable_0) = -999
      Rho12(ntable_0) = -999
      P12(ntable_0)  = -999
      P23(ntable_0)  = -999
      P34(ntable_0)  = -999
        n12 = 0
        n23 = 0
        n34 = 0
      call  GetCorrs(iRho12, iRho23, iRho34, iP12, iP23, iP34, n12,n23, n34)

        if (n12.gt.0) then
         Rho12(ntable_0) = iRho12
         P12(ntable_0) = iP12
      endif

      if (n23.gt.0) then
         Rho23(ntable_0) = iRho23
         P23(ntable_0) = iP23
      endif
      
      if (n34.gt.0) then
         Rho34(ntable_0) = iRho34
           P34(ntable_0) = iP34
      endif

c      write (6,*) 'Rho12, Rho23, Rho34, P12, P23, P34: ', Rho12, Rho23, Rho34, P12, P23, P34
c      call exit(0)



cccc cccc cccc cccc cccc cccc cccc cccc cccc


c      write (52,'(i6,4(2i4,1x),4(3f8.2,1x))') ntable_0,(N_Mstat(ntable_0,ib),M_Mstat(ntable_0,ib),ib=1,4),
c     1    (WPROmagR (ntable_0,ib,1),WPROmagR (ntable_0,1,2),WPROmagR (ntable_0,ib,3), ib=1,4)

c      write (53,'(i6,4(i3,f9.3,f9.3,f9.3))') ntable_0,(M_Mstat(ntable_0,ib),
c     1    medsky(ntable_0,ib),medskysig(ntable_0,ib),formalsig(ntable_0,ib),ib=1,4)

c      write (6,'(i6,4(i3,f9.3,f9.3,f9.3))') ntable_0,(M_Mstat(ntable_0,ib),
c     1    medsky(ntable_0,ib),medskysig(ntable_0,ib),formalsig(ntable_0,ib),ib=1,4)



 600      continue   ! Jsrc

c      call exit(0)

 8000      continue                ! kreg_ix    ;   end of focal plane division

 8010      continue            ! Early exit from the region loop

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      deallocate(MapPSFs)                                 ! ADDED
        if (istatus /= 0) call exit(istatus)
c      goto 886

      What = "ntable_0"
        TYPE = "i"
        metval = ntable_0*1.
        comment = 'total number of extracted sources'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)


      if (ireg*jreg.gt.1) pointless = .false.

c      write (6,*) 'here'
c      write (6,*) nsubxX,nsubyY
      
      call QA_stats (nsubxX,nsubyY,nsubx,nsuby,Array,pix_order,Larray,Lsize,nf,
     1     nfpix,NSUMmax,Table_0, zero, wflag, ntable_0,
     1     crval1,crval2,qdir,
     1     imeta, What, comment, Type,
     1     basename, coname, smode, pointless, fram, SPIT,MIPS,fluxcon, detsnr)


c reset the RAlist,DEClist arrays;   single-precision
      do j=1,ntable_0
            RAlist(J) = RAlist_0(J)
            DEClist(J) = DEClist_0(J)
      enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c TEMP
      if (debug) then
      ic=0
            do j=1,nyAWAIC
            do i=1,nxAWAIC
             ic=ic+1
             larray(ic) = Mdupe (i,j)
            enddo
            enddo

            fout = 'dupe.fits'
            write (6,'(a)') fout(1:50)
            call wfits (nxAWAIC,nyAWAIC,lsize,larray,fout)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c   COADD time -- initialization jazz

 886      if (verbose) write (6,*) 'deallocating'
c      if (nf.gt.1) then                ! JWF B21114 commented out for test
c       deallocate (Xpos)               ! JWF B21114 commented out for test
c        deallocate (Ypos)               ! JWF B21114 commented out for test
c      endif                            ! JWF B21114 commented out for test

        deallocate (SNRlist)
        deallocate (LBACK)
        deallocate (LSIG)
        deallocate (Lconf)
        deallocate (MAGSTD)
        deallocate (eMAGSTD)
        deallocate (framefluxes)
        deallocate (frameuncs)
        deallocate (FLGSTD)
        deallocate (RAsublist)
        deallocate (DECsublist)
        deallocate (IDlist)
c        deallocate (Table)

      deallocate (Array)
        deallocate (Unc)
        deallocate (iMask)
        deallocate (Mask)
      deallocate (HDR)


c      command = "free -ot"
c        s0 = ""
c        call unix (command,s0)
c      write (6,*) 'done'


      if (verbose) write (6,*) 'allocate coadd arrays'

      allocate (COADD(nxAWAIC,nyAWAIC,4))
        if (.not.allocated(COADD)) then
          write (6,*) 'allocation failed for COADD'
          istat = 5
          call exit(istat)
        end if

      allocate (COUNC(nxAWAIC,nyAWAIC,4))
        if (.not.allocated(COUNC)) then
          write (6,*) 'allocation failed for COUNC'
          istat = 5
          call exit(istat)
        end if

      allocate (iCOADD(nxAWAIC,nyAWAIC,4))
        if (.not.allocated(iCOADD)) then
          write (6,*) 'allocation failed for iCOADD'
          istat = 5
          call exit(istat)
        end if

      allocate (COCOV(nxAWAIC,nyAWAIC,4))
        if (.not.allocated(COCOV)) then
          write (6,*) 'allocation failed for COCOV'
          istat = 5
          call exit(istat)
        end if

      allocate (COmsk(nxAWAIC,nyAWAIC,4))
        if (.not.allocated(COmsk)) then
          write (6,*) 'allocation failed for COmsk'
          istat = 5
          call exit(istat)
        end if



ccccc load COADD arrays

      do ib=1,4
                cexist(ib) = .false.
                sexist(ib) = .false.
        enddo

        LC = numchar(coname)
        LL = numchar (clevel)

        do ib=1,4

                cexist(ib) = .false.
                sexist(ib) = .false.

c                if (wflag(1,ib) .eq. 1) then
            if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

                  write (cc,'(i1)') ib
                  if ((MIPS).and.(ib.eq.4)) then
                        cc = '5'
                  endif

                  coaddfits(ib) = coname(1:LC) // '-w' // cc // '-img-' // clevel(1:LL) // '.fits'    ! unwise
              COUfits(ib)   = coname(1:LC) // '-w' // cc // '-std-' // clevel(1:LL) // '.fits'    ! unwise
c              COCOVfits(ib) = coname(1:LC) // '-w' // cc // '-invvar-' // clevel(1:LL) // '.fits'    ! unwise
	      COCOVfits(ib) = coname(1:LC) // '-w' // cc // '-n-' // clevel(1:LL) // '.fits'    ! unwise
              COmskfits(ib) = coname(1:LC) // '-w' // cc // '-msk-' // clevel(1:LL) // '.fits'    ! unwise

              if (clevel(1:LL).eq.'1b') then  ! using 1b frames for "coadd" measurements; msk = cmsk
                  COmskfits(ib) = coname(1:LC) // '-w' // cc // '-msk-' // clevel(1:LL) // '.fits'
              endif

	write (6,'(a,a)') 'coadd unc ',COUfits(1)(1:99)
	write (6,'(a,a)') 'coadd cov  ',COCOVfits(1)(1:99)
	write (6,'(a,a)') 'coadd msk ',COmskfits(1)(1:99)


               L = numchar(coaddfits(ib))
                  zexist = .false.
                  erase = .false.
c                 call access(coaddfits(ib),zexist,erase)
 		    	  zexist = Access(coaddfits(ib)(1:LNBlnk(coaddfits(ib))),' ') .eq. 0   ! JWF B60711
                  if (.not.zexist) then
                   ! see if the file is compressed
                    gztmp = coaddfits(ib)(1:L) // '.gz'
                    zexist = .false.
                    erase = .false.
c                   call access(gztmp,zexist,erase)
  		    	    zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
                    if (zexist) then
                        coaddfits(ib) = gztmp
                    endif
                  endif

ccccccccccc
                zexist = .false.
                erase = .false.
c               call access(coaddfits(ib),zexist,erase)
 		    	zexist = Access(coaddfits(ib)(1:LNBlnk(coaddfits(ib))),' ') .eq. 0   ! JWF B60711

                if (zexist) then
                   cexist(ib) = .true.
               call readimage
     1     (nnnx,nnny,lsize,Larray,iLarray,coaddfits(ib),
     1     zzz,zzz,zzz,zzz,zzz,zzz,zzz, tJD)

            if (verbose) write (6,*) cexist(ib),' ',coaddfits(ib)(1:99)
            if (verbose) write (6,*) ib,'  coadd image loaded with size: ',ncx(ib),ncy(ib)
            if (verbose) write (6,*) 'coadd pixel scale: ',cscale(ib)

             if ((ncx(ib).gt.ncoaddsize).or.(ncx(ib).gt.ncoaddsize)) then
                        write (6,*) 'ERROR -- coadd exceed ',
     +                    'allocations; max = ncoaddsize x ncoaddsize'
                        print *,'ncoaddsize =',ncoaddsize,'; squared =',ncoaddsize**2
                        call exit(64)
                  endif

            ic = 0
                do jj=1,ncy(ib)
                do ii=1,ncx(ib)
                        ic=ic+1
                        coadd(ii,jj,ib) = Larray(ic)
                  COmsk(ii,jj,ib) = 0

                  if (SPIT)  then   !  convert SB to DN/s and multiply by the Loc-Dep Photom correction
                          if ( (MIPS).and.(ib.eq.4) ) then
                           coadd(ii,jj,ib) = coadd(ii,jj,ib) / fluxcon(5)
                          else if ( (MIPS).and.(ib.eq.3) ) then
                           coadd(ii,jj,ib) = coadd(ii,jj,ib) / fluxcon(4)
                          else
                          coadd(ii,jj,ib) =  coadd(ii,jj,ib) / fluxcon(ib)
                          endif
                      endif

                  if (Larray(ic).eq.-9991) then
                                coadd(ii,jj,ib) = -9991
                        endif

                  enddo
                  enddo

            else
             write (6,*) 'ERROR -- coadd does not exist:'
             L = numchar(coaddfits(ib))
            write (6,'(a)') coaddfits(ib)(1:L)
                 call exit(64)

            endif


c ancillary coadd mosaics
            L = numchar(COmskfits(ib))
                zexist = .false.
                erase = .false.
c               call access(COmskfits(ib),zexist,erase)
 		    	zexist = Access(COmskfits(ib)(1:LNBlnk(COmskfits(ib))),' ') .eq. 0   ! JWF B60711
                if (.not.zexist) then
                  ! see if the file is compressed
                  gztmp = COmskfits(ib)(1:L) // '.gz'
                  zexist = .false.
                  erase = .false.
c                 call access(gztmp,zexist,erase)
   		    	  zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
                  if (zexist) then
                        COmskfits(ib) = gztmp
                  endif
                endif

            L = numchar(COUfits(ib))
                zexist = .false.
                erase = .false.
c               call access(COUfits(ib),zexist,erase)
 		    	zexist = Access(COUfits(ib)(1:LNBlnk(COUfits(ib))),' ') .eq. 0   ! JWF B60711
                if (.not.zexist) then
                  ! see if the file is compressed
                  gztmp = COUfits(ib)(1:L) // '.gz'
                  zexist = .false.
                  erase = .false.
c                 call access(gztmp,zexist,erase)
   		    	  zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
                  if (zexist) then
                        COUfits(ib) = gztmp
                  endif
                endif


            L = numchar(COCOVfits(ib))
                zexist = .false.
                erase = .false.
c               call access(COCOVfits(ib),zexist,erase)
   		    	zexist = Access(COCOVfits(ib)(1:LNBlnk(COCOVfits(ib))),' ') .eq. 0   ! JWF B60711
                if (.not.zexist) then
                  ! see if the file is compressed
                  gztmp = COCOVfits(ib)(1:L) // '.gz'
                  zexist = .false.
                  erase = .false.
c                 call access(gztmp,zexist,erase)
   		    	  zexist = Access(gztmp(1:LNBlnk(gztmp)),' ') .eq. 0   ! JWF B60711
                  if (zexist) then
                        COCOVfits(ib) = gztmp
                  endif
                endif

            

c COUNC
            zexist = .false.
                erase = .false.
c               call access(coufits(ib),zexist,erase)
   		    	zexist = Access(coufits(ib)(1:LNBlnk(coufits(ib))),' ') .eq. 0   ! JWF B60711

                if (zexist) then
               call readimage
     1     (nnnx,nnny,lsize,Larray,iLarray,coufits(ib),
     1     zzz,zzz,zzz,zzz,zzz,zzz,zzz, tJD)

                if (verbose) write (6,*) coufits(ib)(1:99)
                if (verbose) write (6,*) ib,'  coadd unc image loaded '

                 if ((nnnx.gt.ncoaddsize).or.(nnny.gt.ncoaddsize)) then
                        write (6,*) 'ERROR -- coadd unc exceed ',
     +                    'allocations; max = ncoaddsize x ncoaddsize'
                        print *,'ncoaddsize =',ncoaddsize,'; squared =',ncoaddsize**2
                        call exit(64)
                  endif

                  ic = 0
                  do jj=1,ncy(ib)
                  do ii=1,ncx(ib)
                        ic=ic+1
                        counc(ii,jj,ib) = Larray(ic)

             if (SPIT)  then   !  convert SB to DN/s and multiply by the Loc-Dep Photom correction
                   if ( (MIPS).and.(ib.eq.4) ) then
                        COUNC(ii,jj,ib) =  COUNC(ii,jj,ib) / fluxcon(5)
                   else if ( (MIPS).and.(ib.eq.3) ) then
                        COUNC(ii,jj,ib) =   COUNC(ii,jj,ib) / fluxcon(4)
                   else
                       COUNC(ii,jj,ib) =  COUNC(ii,jj,ib) / fluxcon(ib)
                   endif

            else
                  counc(ii,jj,ib) = Larray(ic)

                endif

                  enddo
                  enddo


            endif


ccccc   Coadd Coverage

            if (clevel(1:2).eq.'1b') then  ! using 1b frames for "coadd" measurements;  force COV = 1

              do jj=1,ncy(ib)
                    do ii=1,ncx(ib)
                        cocov(ii,jj,ib) = 1.
                    enddo
                  enddo

            endif


                zexist = .false.
                erase = .false.
c               call access(cocovfits(ib),zexist,erase)
   		    	zexist = Access(cocovfits(ib)(1:LNBlnk(cocovfits(ib))),' ') .eq. 0   ! JWF B60711

                if (zexist) then
               call readimage
     1     (nnnx,nnny,lsize,Larray,iLarray,cocovfits(ib),
     1     zzz,zzz,zzz,zzz,zzz,zzz,zzz, tJD)

                if (verbose)  write (6,*) cocovfits(ib)(1:99)
                if (verbose) write (6,*) ib,'  coadd cov image loaded'

                  if ((nnnx.gt.ncoaddsize).or.(nnny.gt.ncoaddsize)) then
                        write (6,*) 'ERROR -- coadd cov exceed ',
     +                    'allocations; max = ncoaddsize x ncoaddsize'
                        print *,'ncoaddsize =',ncoaddsize,'; squared =',ncoaddsize**2
                        call exit(64)
                  endif

                  ic = 0
                  do jj=1,ncy(ib)
                  do ii=1,ncx(ib)
                        ic=ic+1
c                        cocov(ii,jj,ib) = Larray(ic)
			cocov(ii,jj,ib) = iLarray(ic) * 1.  ! unwise
	
	       ! unwise  (this is no longer needed;  not using the invvar image !!
c                        if (clevel(1:1).eq.'u') then
c                               if (ib.eq.1)  cocov(ii,jj,ib) = cocov(ii,jj,ib) * 300.   ! INVVAR scaling to get roughly the coverage number
c                               if (ib.eq.2)  cocov(ii,jj,ib) = cocov(ii,jj,ib) * 3480.   ! INVVAR scaling to get roughly the coverage number
c                        endif

                  enddo
                  enddo


                endif


ccccc   Coadd mask

                zexist = .false.
                erase = .false.
c               call access(comskfits(ib),zexist,erase)
   		    	zexist = Access(comskfits(ib)(1:LNBlnk(comskfits(ib))),' ') .eq. 0   ! JWF B60711

                if (zexist) then
                   call readimage
     1     (nnnx,nnny,lsize,Larray,iLarray,comskfits(ib),
     1     zzz,zzz,zzz,zzz,zzz,zzz,zzz, tJD)

c                if (verbose)  write (6,*) comskfits(ib)(1:99)
                if (verbose) write (6,*) ib,'  coadd mask image loaded'

                  if ((nnnx.gt.ncoaddsize).or.(nnny.gt.ncoaddsize)) then
                        write (6,*) 'ERROR -- coadd msk exceed ',
     +                    'allocations; max = ncoaddsize x ncoaddsize'
                        print *,'ncoaddsize =',ncoaddsize,'; squared =',ncoaddsize**2
                        call exit(64)
                  endif

                  ic = 0
                  do jj=1,ncy(ib)
                  do ii=1,ncx(ib)
                        ic=ic+1
c                        comsk(ii,jj,ib) = iLarray(ic)
			comsk(ii,jj,ib) = int(Larray(ic))  ! unwise
                  enddo
                  enddo


                endif




        endif

      enddo


      NsrcAll = ntable_0

      if (verbose)  write (6,*) 'call CoaddStarMask '

      do ib=1,4

c       if (wflag(1,ib).gt.0) then
       if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

      call CoaddStarMask  (ib,NSUMmax,NsrcAll,
     1       ncoaddsize,ncx, ncy, WCSc,fwhm,cscale,
     1       nf,iCOADD, Rstann, Rstwid, RAlist, Declist, SPIT)

      if (ib.eq.-1) then

        ic=0
        do jj=1,ncy ( ib )
        do ii=1,ncx ( ib )
        ic=ic+1
        larray(ic)=iCOADD(ii,jj,ib)*1.
        enddo
        enddo

        fout = 'icoadd.fits'
        call  wimage (ncx ( ib ),ncy ( ib ),lextra,larray,coaddfits(ib),fout)


      endif

       endif

      enddo


cccccccccc
c XSC
c
cccc XSClist


      write (6,*) 'allocate XSCprox with ',ntable_0
        allocate (XSCprox(ntable_0))
        if (.not.allocated(XSCprox)) then
          write (6,*) 'allocation failed for XSCprox'
          istat = 5
        end if


        ngal = 0

        zexist = .false.
        erase = .false.
c       call access(xscf,zexist,erase)
   		zexist = Access(xscf(1:LNBlnk(xscf)),' ') .eq. 0   ! JWF B60711
        if ((zexist).and.(.not.smode))  then

      write (6,*) 'allocate XSCshape with ',ntable_0
         allocate (XSCshape(ntable_0,5))
         if (.not.allocated(XSCshape)) then
          write (6,*) 'allocation failed for XSCshape'
          istat = 5
         end if

          zexist = .false.
          erase = .false.
c         call access(xscl,zexist,erase)
    	  zexist = Access(xscl(1:LNBlnk(xscl)),' ') .eq. 0   ! JWF B60711
          if (zexist) then

            call loadXSC (nf,wflag,xscf, xscl, WCSc, ncx,ncy, ngal, XSCshape,
     1       ntable_0, NSUMmax, RAlist, Declist,  XSCprox, SPIT, MIPS)

          endif

        endif

      ngal = 0

        zexist = .false.
        erase = .false.
c       call access(xscf,zexist,erase)
   	    zexist = Access(xscf(1:LNBlnk(xscf)),' ') .eq. 0   ! JWF B60711
        if (zexist) then

      if (allocated(XSCshape)) deallocate (XSCshape)

        allocate (XSCshape(NsrcAll,5))
          if (.not.allocated(XSCshape)) then
           write (6,*) 'allocation failed for XSCshape'
           istat = 5
          end if
      
        deallocate (XSCprox)
          allocate (XSCprox(NsrcAll))
          if (.not.allocated(XSCprox)) then
            write (6,*) 'allocation failed for XSCprox'
            istat = 5
          end if

          zexist = .false.
          erase = .false.
c         call access(xscl,zexist,erase)
   	      zexist = Access(xscl(1:LNBlnk(xscl)),' ') .eq. 0   ! JWF B60711
          if ((zexist).and.(.not.smode))  then

            call loadXSC (nf,wflag,xscf, xscl, WCSc, ncx,ncy, ngal, XSCshape,
     1       NsrcAll, NSUMmax, RAlist, Declist,  XSCprox, SPIT, MIPS)

          if (verbose) write (6,*) 'number of XSC galaxies loaded ',ngal

          endif

        endif

      What = "ngal"
        TYPE = "i"
        metval = ngal*1.
        comment = 'number of XSC sources matched'
        iband = 0
        if (smode) iband = ib0
        call MetWr (imeta, iband, What, TYPE, metval, comment)

cccc
cccc  COADD measurements

      Table = -99.99  ! TJ 23aug2016

      dd1 = 0.
      call  ETime (Second, dd1)
      dd1 = second(1)

      ih = NsrcAll/2
      i3q =(3* NsrcAll)/4
      iq = NsrcAll/4

      if (verbose) write (6,*) 'WAPPCO:   run through WPRO extractions ',NsrcAll   ! WAPPCO  Wappco  WappCo  wappco

      NaperC = Naper

c      command = "free -ot"
c        s0 = ""
c        call unix (command,s0)
c        write (6,*) 'done'

c      write(6,*)'ks(2,1): ',ks(2,1)

      do 700  Jsrc = 1,NsrcAll    ! run through each WPRO source

        if ((Jsrc.eq.iq).and.(verbose)) write (6,*) '25% WAPPco complete'
        if ((Jsrc.eq.ih).and.(verbose)) write (6,*) '50% WAPPco complete'
        if ((Jsrc.eq.i3q).and.(verbose)) write (6,*) '75% WAPPco complete'

c        ra0 = Table_0(Jsrc,1) * 1.d0
c        dec0 = Table_0(Jsrc,2) * 1.d0
       ra0 = RAlist_0(Jsrc)
       dec0 = DEClist_0(Jsrc)
c       nblend = int(Table_0(Jsrc,19))    ! JWF B30424
       nblend = nint(Table_0(Jsrc,19))   ! JWF B30424

c      write (6,*) 'DEBUG-1 ',Jsrc,ra0,dec0, nblend

        chx = 0.
          chy = 0.

        do 592 ib=1,4

              Rcoadd(ib) = 0.
              cov_ave(ib) = 0.
                  cov_std(ib) = 0.

              CoaddXY (1,ib) = 0.
              CoaddXY (2,ib) = 0.

              do KK=1,NaperC
                  cmag (KK,ib) = 99.
                  cemag (KK,ib) = 0.
                  coflag (KK,ib) = 0
              enddo

c              if (cexist(ib).and.(wflag(1,ib).gt.0)) then
              if ( cexist(ib) .and. ( any ( wflag ( 1:nf,ib ) ==1 )) ) then

c      write (6,*) ib,WCSc ( ib ), ncx ( ib ) ,ncy ( ib ), cscale(ib)

                  call WhatPos ( 1,WCSc ( ib ) ,ncx ( ib ) ,ncy ( ib ) , cscale(ib), '0', ra0, dec0, x0, y0, igo, edgebuf )
      
                    ipix = nint(x0)
                  jpix = nint(y0)

                  CoaddXY (1,ib) = x0
                  CoaddXY (2,ib) = y0

                  if (ipix.lt.1) goto 592
                  if (jpix.lt.1) goto 592
                  if (ipix.gt.ncx(ib)) goto 592
                        if (jpix.gt.ncy(ib)) goto 592

c  compute the centroid using
c                  if ((nblend.eq.1).and.(SPIT)) then
                  
                  if ((ib.eq.4).and.(docentroid)) then  ! to deal with lousy W4 astrometry


c      if (Jsrc.lt.5) goto 700
                    if (ib.le.4) then
                      xmom = x0
                      ymom = y0
c      if ((Jsrc.eq.5).and.(ib.eq.3)) write (6,*) x0,y0

                      do KITE = 1,3
                       if (KITE.eq.1) then
                        findpeak = .true.
                        ibox = 9
                        else if (KITE.le.3) then
                         findpeak = .false.
                         ibox = 9
                       endif
                       call PosCoMoment (ib,ncoaddsize,ncx(ib) ,ncy(ib),Coadd,xmom,ymom,ibox,cx,cy,findpeak)
                       xmom = cx
                       ymom = cy

c            if ((Jsrc.eq.5).and.(ib.eq.3)) write (6,*) kite,xmom,ymom
                     enddo

                     chx = xmom-x0
                     chy = ymom-y0
                    endif
                  else
                     chx = 0.
                     chy = 0.
                  endif

                  x0 = x0 + chx
                  y0 = y0 + chy
                  ipix = nint(x0)
                        jpix = nint(y0)
                  CoaddXY (1,ib) = x0
                        CoaddXY (2,ib) = y0


c      if (Jsrc.eq.5) then
c            write (6,*) ib,x0,y0,chx,chy
c            if (ib.eq.4) call exit(0)
c      endif


                  Rl = 0.
                  Rh = 0.
                  call CoaddAnnulus (Jsrc,ib,ncoaddsize,ncx(ib) ,ncy(ib),Coadd,Counc,iCOADD,
     1                  x0,y0, Rstann, Rstwid, nannC, Csky, Csdev, 
     1                  Nbann, Rl, Rh, pscale(ib), cscale(ib), Fcorr(ib), BGmax(ib), SKY_pop)

c      !TEMPORARY
c                  Csky = Csky * 1.

c      write (6,*) x0,y0,Coadd(3601,2709,4),Counc(3601,2709,4)
c      write (6,*) Rstann
c      write (6,*) Rstwid
c      write (6,*) Rl, Rh
c       write (6,*)  Nbann,nannC, Csky, Csdev



                  T2 = (Csdev**2) - (SKY_pop**2)

                        if (T2.ge.0.) then
                                confusion_noise = sqrt (T2)
                        else
                                confusion_noise = 0.
                        endif

c                  write (54,'(i6,i3,2f8.1,i6,i6,3f12.4)') Jsrc,ib,x0,y0,nannC,Nbann,Csky, Csdev,SKY_pop

                        Rcoadd(ib) = Rcirc(1,ib)*cscale(ib)   ! arcseconds

                        if ((Csky.lt.-500.).or.(Csdev.le.0.)) goto 592   ! give up


                  call CoaddPhot (Jsrc,ib,ncoaddsize,ncx, ncy,
     1                  Coadd, Counc, CoCov, CoMsk, iCOADD, cscale(ib),
     1                 x0,y0, NaperC, Rcirc , nannC, Csky, Csdev, Nbann, cozero,
     1                 mfluxtmp, mefluxtmp, mzmag, mzerr, miflag, SPIT, MIPS, fbits, 
     1                 Fcorr(ib), confusion_noise, cov_ave(ib), cov_std(ib))


                  if (debug) then
c                  do KK=1,NaperC
		   do KK=2,3    ! THJ 14Dec2017

		     if ((Jsrc.gt.17000).and.(Jsrc.lt.19000)) then
                    write (6,*) Jsrc,ib,KK,Rcirc(KK,ib), Rcirc(KK,ib)*cscale(ib),        
     1                     mfluxtmp(kk),mefluxtmp(kk),mzmag(kk),mzerr(kk), miflag(kk)
		     endif

                  enddo
                  endif

c      if (ib.eq.4) then
c            kk=1
c            write (6,*) mzmag(KK),CmosaicCorr(ib),mzerr(KK),miflag(KK)
c            call exit(0)
c      endif


                  do KK=1,NaperC

                     if (SPIT) then

                      Rcoadd(ib) = 3.   ! arcsec
                      cmag(KK,ib) = mzmag(KK)
                      if (KK.eq.1) cmag(KK,ib) = mzmag(KK) - IRACapcor(ib)   ! with the IRAC aperture correction

                      if (MIPS) then
                       if (ib.eq.4) then 
                                 Rcoadd(ib) = 9.  !  arcsec
                             endif
                            endif
                  
                    else
                        cmag(KK,ib) = mzmag(KK)
                        if (KK.eq.1) cmag(KK,ib) = mzmag(KK) - CmosaicCorr(ib)   ! with the aperture correction
                    endif

                    cemag(KK,ib) = mzerr(KK)
                          coflag (KK,ib) = miflag(KK)

!!!!!!!!!!!!!!!! special case, when nf = 1 (single frame coverage), check that the coverage metric is unity

                    if ((nf.eq.1).and.(cov_ave(ib).lt.2.).and.(.not.SPIT).and.(coflag (KK,ib).ne.32)) then
                        idex = nint(x0)
                        jdex = nint(y0)
                        ! check the core and adjacent pixels
                        do JJJ = jdex-1,jdex+1
                        do III = idex-1,idex+1
                          if (CoCov(iii,jjj,ib).lt.0.9) coflag (KK,ib) = max (2,coflag (KK,ib))
                          if (CoMsk(iii,jjj,ib).gt.0) coflag (KK,ib) = max (2,coflag (KK,ib))
                        enddo
                        enddo
                              
                    endif

                  enddo

      

ccccccccccccc  XSC ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c XSCshape (nsrc,1) = rk20fe   # corrected
c XSCshape (nsrc,2) = ba       # corrected
c XSCshape (nsrc,3) = pa
c XSCshape (nsrc,4) = rk20fe   # corrected for wise W4
c XSCshape (nsrc,5) = ba       # corrected for wise W4

                galflag(ib) = 0
                galmag(ib) = 99.
                galerr(ib) = 0.
                 galba (ib) = 0.
                galpa (ib) = 0. 
                galRsemi(ib) = 0.

               if (ngal.gt.0) then

                 if ((XSCprox(Jsrc).ge.0.).and.(XSCprox(Jsrc).le.2.)) then  ! galaxy core
                  r20 = XSCshape (Jsrc, 1)
                  ba = XSCshape (Jsrc, 2)
                  pa = XSCshape (Jsrc, 3)

c      if (ib.eq.1) write (86,*) 'DEB ',x0,y0, r20,ba,pa

                  rot_angle = pa + scrot(ib)  ! take into account the rotation of the image (if any)

                  call dcon (rot_angle,pa_math)

c      if (ib.eq.1) write (86,*)  'DEB ',scrot(ib),pa_math


                  if (ib.eq.4) then
                    r20 = XSCshape (Jsrc, 4)
                    ba = XSCshape (Jsrc, 5)
                  endif

                  if (SPIT) then
                    r20 = XSCshape (Jsrc, 1)
                    ba = XSCshape (Jsrc, 2)

c                    ! minimum radius should be the IRAC standard aperture == 12.22 arcsec
c                    r20 = max (r20, 12.22)

                  endif

                  if ((MIPS).and.(ib.eq.4)) then
                      r20 = XSCshape (Jsrc, 4)
                      ba = XSCshape (Jsrc, 5)

c                      r20 = max (r20, 12.22)

                  endif


                  pscale_factor = pscale(ib) / cscale(ib)
                        RLannulus = Rstann(ib)*pscale_factor  ! coadd inner annulus radius in pixels

                  r20pix = r20 /  cscale(ib)   ! semi-major axis in pixels
                  RR = r20pix * 1.2

                  RLannulus = max (RLannulus, RR)
c limit the background annulus to a max radius of 3 arcmin
                  RmaxAnn = 3.0 * 60. / cscale(ib)

                  RLannulus = min (RLannulus, RmaxAnn)

                  r20pix = min (r20pix, RLannulus)

                  galRsemi (ib) = r20pix * cscale(ib)  ! arcseconds
                  galba (ib) = ba
                  galpa (ib) = pa

c            write (6,*) ib, x0,y0
c              write (6,*) r20, ba, pa

c                  call CoaddAnnulus (Jsrc,ib,ncoaddsize,ncx(ib),ncy(ib),Coadd, Counc,iCOADD, ! JWF B21114
c    1                  x0,y0, Rstann, Rstwid, nannC2, Csky2, Csdev2,                              ! JWF B21114
c    1                  Nbann2, Rlannulus, Rhannulus, pscale(ib), cscale(ib), Fcorr(ib), SKY_pop2) ! JWF B21114
c=================================== mod by JWF B21114 ========================================================
                  call CoaddAnnulus (Jsrc,ib,ncoaddsize,ncx(ib),ncy(ib),
     +                  Coadd, Counc,iCOADD,
     1                  x0,y0, Rstann, Rstwid, nannC2, Csky2, Csdev2,
     1                  Nbann2, Rlannulus, Rhannulus, pscale(ib), cscale(ib),
     +                  Fcorr(ib), BGmax(ib), SKY_pop2)
c=================================== end of mod by JWF B21114 =================================================


c      write (6,*) RLannulus, Rhannulus
c      write (6,*) Nbann2,nannC2, Csky2, Csdev2

                  galmag(ib) = 99.
                  galerr(ib) = 0.

c      if (ib.eq.1) write (86,*)  'DEB ',r20pix, ba, pa_math

                  call CoaddPhotGAL (Jsrc,ib,ncoaddsize,ncx, ncy,
     1                  Coadd, Counc,  CoCov, CoMsk,  iCOADD,pscale(ib), cscale(ib),
     1                 x0,y0, NaperC, r20pix, ba, pa_math, nannC2, Csky2, Csdev2, Nbann2, cozero,
     1                 galmag(ib), galerr(ib), galflag(ib), SPIT, MIPS, fbits,
     1                 Fcorr(ib), confusion_noise)

c                  write (6,*) r20pix, ba, pa_math
c                  write (6,*) nannC2, Csky2, Csdev2
c                  write (6,*) galmag(ib), galerr(ib), galflag(ib)
c                  write (6,*) cozero(ib)
c                  write (6,*) ' '



                endif
                endif


                endif

 592      continue  ! ib

c      write(6,*)'ks(3,1): ',ks(3,1)

ccccccccccccccc  WRITE IT OUT  ******************************************************************

c                ngal = 0

                iunit=75
            munit=85
                call writeTABLE (iunit,munit,NSUMmax,Table_0,     ! CJG B30228
     1           NsrcAll, zero,wflag, nf, nfmep,             ! TPC
     1           meptable, namtable, jdtable,                     ! CJG B30213  
     1           ireg,jreg,Jsrc,
     1           nsubsrcA,Xpos, YPos, RAlist, DEClist, IDList_0, SatNUM, iverify,
     1           NaperC,
     1           N_Mstat,M_Mstat,
     1           medsky,medskysig, formalsig,
     1           cscale, Rcirc,
     1           hname, outfile,meproot,mepwrite,Fcorr,              ! CJG B30228
     1           Rcoadd,cmag,cemag, coflag, imID, cov_ave, cov_std,
     1           WPROmagR, CmosaicCorr,
     1           CoaddXY, mLogQ, ks, nDF,  Rho12, Rho23, Rho34, P12, P23, P34,
     1           AVE_mJD,mJD_min,mJD_max,
     1           ngal,XSCprox,BGmax,
     1           galRsemi,galba,galpa,galmag,galerr,galflag, smode, verbose, SPIT,
     1           mfiles)                                            ! CJG B30228

 700      continue  ! Jsrc

      close (iunit)
c===================== start of code added by CJG B30228 =======================

      if(mepwrite) then
          What = "mfiles"
          TYPE = "i"
          metval = mfiles*1.0
          comment = 'number of MEP files.'
          iband = 0
          if (smode) iband = ib0
          call MetWr (imeta, iband, What, TYPE, metval, comment)
        write(6,*)'deallocating meptable, jdtable, and namtable.'
        deallocate (mepTable)      ! CJG B30228
        deallocate (jdTable)      ! CJG B30228
        deallocate (namTable)      ! CJG B30228
      end if

c===================== end of code added by CJG B30228 =======================

      call  ETime (Second, dd2)
      dd2 = second(1)

        dd=dd2-dd1

c       if (verbose) write (6,*) 'total number of sources processed ',NsrcAll
c        rtime = nf * (NsrcAll*1.) / dd
        if (verbose) write (6,'(a,f10.3,a)')'  DTime ',dd,' sec'
c        if (verbose) write (6,*) 'sources per second / frame ',rtime


      goto 847

cccc
c statistics on the WPRO-WAPPco mag

        L = numchar (coname)
      LQ = numchar (qdir)     

! strip the path off the name
        idex = 0
        do k=L,1,-1
          if (coname(k:k).eq.'/') then
                idex = k
                goto 950
          endif
        enddo

 950    if (idex.gt.0) then
                coname = coname(idex+1:L)
                L = numchar (coname)
        endif


      do ib=1,4

c       if ( wflag ( 1,ib ) .eq.1 ) then
        if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then

        LD = numchar (DATE)

        m=0
        do Jsrc=1,NsrcAll
            if (Table(Jsrc,ib).gt.-5.) then
                  m=m+1
                  XSCprox(m) = Table(Jsrc,ib)
c                  if (ib.eq.1) write (86,*) m,XSCprox(m)
            endif
        enddo
      
        zlow = -1.
        zhigh = 1.
        iter = 7
        call iSTATS (XSCprox,larray,m,zlow,zhigh,iter,ave,sdev,xmed,mv)

        if (ib.eq.1) then
                s0 = qdir(1:LQ) // '/' // coname(1:L) // '-w1-mdexqamdiff.tbl'
          else if (ib.eq.2) then
                s0 = qdir(1:LQ) // '/' // coname(1:L) // '-w2-mdexqamdiff.tbl'
          else if (ib.eq.3) then
                s0 = qdir(1:LQ) // '/' // coname(1:L) // '-w3-mdexqamdiff.tbl'
          else if (ib.eq.4) then
                s0 = qdir(1:LQ) // '/' // coname(1:L) // '-w4-mdexqamdiff.tbl'
          endif

        open (unit=32,file=s0)
          write (32,'(a)') '\ global WPRO-WAPPco statistics'
        write (32,'(a,i5)') '\nf = ',nf
          write (32,'(a)') '| N   | median |  ave   |  sdev  | nf  |      DATE             |'
        write (32,'(a)') '| --  | mag    |  mag   |  mag   | --  |      --               |'
        write (32,'(a)') '| --  | null   |  null  |  null  | --  |      --               |'

        if ((mv.le.0).or.(ave.lt.-5.)) then
          write (32,'(i6,a,a,a,i6,1x,a)') mv,'  null   ','  null   ','  null   ',nf,DATE(1:LD)
        else
            write (32,'(i6,3f9.3,i6,1x,a)') mv,Xmed,Ave,sdev,nf,DATE(1:LD)
        endif
          close (32)
          write (6,'(a,i3,i6,3f9.3)') 'global WPRO-WAPPco stats (n,med,ave,sig): ',ib,mv,Xmed,Ave,sdev

       endif

      enddo


 847      if (iverify.eq.0) then

        ! do this step when absolutely no sources were extracted
        L = numchar(outfile)
        if (verbose) write (6,*) 'WARNING -- no sources extracted'
        if (verbose) write (6,'(a,1x,a)') 'write table: ',outfile(1:L)
        open (unit=iunit,file=outfile)
        write (iunit,'(a,i7)') '\Nsrc = 0'

        Lmax = 9999
          if (smode) Lmax = 1630
        open (unit=11,file=hname)
          do jj=1,999
           read (11,'(a)',end=910) string
           L = numchar(string)
           L = min (L,Lmax)
           write (iunit,'(a)') string(1:L)
          enddo
 910    close (11)
      close (iunit)

      endif


      close (75)

      call  ETime (Second2, dd10)
        dd10 = second2(1)

        dd=dd10-dd9

        if (Nwpros-nBorderViolators .ne. 0) then
          fDup = float(NsrcAll)/float(Nwpros-nBorderViolators)
        else
          fDup = 1.0
        end if
        print *,'Total no. of input sources:          ', nsr
        print *,'Total no. of output sources:         ', NsrcAll
        print *,'No. of wpro_v6 sources:              ', Nwpros
        if (ireg*jreg .eq. 1)
     +    print *,'No. of deleted off-tile sources:     ', nBorderViolators                 ! JWF B30228
        print *,'No. of extractions rejected as dups: ', nduprej                          ! TPC
        print *,'No. of extractions outside region:   ', nbordrej                         ! TPC
        print *,'No. of extractions w/ NaN positions: ', nnanposrej                       ! TPC
        print *,'Duplication factor:           ', 1.0/fDup                                ! JWF B30228
        print *,'No. of PM solutions:                 ', NInt(fDup*nPMsolns)              ! JWF B30228
        print *,'No. of active deblends:              ', NInt(fDup*N_adb)                 ! JWF B30228
        print *,'No. of runaway/bad fits:             ', NInt(fDup*nBadBlend)             ! JWF B30228
        if (TossZeroFlux) then                                                            ! JWF B30404
          print *,'No. of all-flux-zero rejects(1):     ', NInt(fDup*nAllZero(1))         ! JWF B30228
          print *,'No. of all-flux-zero rejects(2):     ', NInt(fDup*nAllZero(2))         ! JWF B30228
          print *,'No. of all-flux-zero rejects(3):     ', NInt(fDup*nAllZero(3))         ! JWF B30228
        else                                                                              ! JWF B30404
          print *,'No. of all-flux-zero rejects:        ', NInt(fDup*nAllZero(3))         ! JWF B30404
        end if                                                                            ! JWF B30404
        print *,'No. of NaN PM solutions:             ', NInt(fDup*nPMNaN)                ! JWF B30228
        if (Nok2b .gt. 0)                                                                 ! JWF B21211
     +   print *,'No. "(n==icp)" test failures:        ', NInt(fDup*Nok2b)                ! JWF B30228
        if (nWgtFail .gt. 0)                                                              ! JWF B30207
     +   print *,'No. of MeanObsEpoch Weight Fallbacks:', NInt(fDup*nWgtFail)             ! JWF B30228
        if (nBlendSwap .gt. -1)                                                           ! JWF B30221
     +   print *,'No. of Blend Swaps:                  ', NInt(fDup*nBlendSwap)           ! JWF B30228
        if (nSingMatFrmFunc .gt. 0)                                                       ! JWF B30221
     +   print *,'No. of framefunc singular matrices:  ', NInt(fDup*nSingMatFrmFunc)      ! JWF B30228
        if (nTossedItAll .gt. 0)                                                          ! JWF B30402
     +   print *,'No. of complete solutions tossed:    ', NInt(fDup*nTossedItAll)         ! JWF B30402
        if (nFallbackMMM .gt. 0)                                                          ! JWF B30606
     +   print *,'No. of mmm fallbacks:                ', NInt(fDup*nFallbackMMM)         ! JWF B30606
        print *,'WPHotpm terminating; istat =         ',istat   ! JWF B21204
        write (6,'(a,f10.3,a)')' total WPHot DTime ',dd,' sec'

        call signoff('WPHotpmc')        ! JWF B30507

      call exit (istat)
      end
      
      include 'includes_wphot.f'
        
