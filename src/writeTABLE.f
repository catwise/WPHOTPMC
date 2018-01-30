c modify to write nulls == Null

	subroutine writeTABLE (iunit,munit,                          ! CJG B30228
     1    nmax,Table, nout, zero, wflag, nf, nfmep,             ! TPC
     1    meptable, namtable, jdtable,                               ! CJG B30213
     1    ireg,jreg,Jsrc,
     1    npos,XPos, YPos, RAlist, DEClist, IDlist, SatNUM, iverify,
     1    Naper,
     1    N_Mstat, M_Mstat,
     1    medsky,medskysig, formalsig,
     1    cscale, R,
     1    hname,outf,meproot,mepwrite,Fcorr,                         ! CJG B30228
     1    Rcoadd,cmag,cemag,coflag, imID, cov_ave, cov_std,
     1    WPROmagR, coApCorr,
     1    CoaddXY,mLogQ, KS, nDF, Rho12, Rho23, Rho34, P12, P23, P34,
     1    AVE_mJD,mJD_min,mJD_max,
     1    ngal,XSCprox,BGmax,
     1    galRsemi,galba,galpa,galmag,galerr,galflag, smode, verbose, SPIT,
     1    mfiles)                                             ! CJG B30228

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*8 RAlist(nmax),DEClist(nmax), dvalue
	real*8 AVE_mJD(nmax,4),mJD_min(nmax,4),mJD_max(nmax,4)
	real*4 XSCprox(nout), galRsemi(4),galba(4),galpa(4),galmag(4),galerr(4)
	real*4 XPOS (npos, nf, 4), YPOS (npos, nf, 4)
	real*4 xfpos(4), yfpos(4)                                   ! CJG B30213
c   	real*4 table (nmax,25),  R(19,4), cscale(4), Rap(19,4)
c   	real*4 table (nmax,59),  R(19,4), cscale(4), Rap(19,4)  !  JWF B30221
   	real*4 table (nmax,65),  R(19,4), cscale(4), Rap(19,4)  !  JWF B31211

	real*4 MAGSTD(19,4), eMAGSTD(19,4), zero (4)
	real*4 MAGWPRO (4), eMAGWPRO (4), rchi(4), BGmax(4)
	real*4 fsignal(4), ffnoise(4)                               ! CJG B30213
	real*4 formalsig(nmax,4), medsky(nmax,4),medskysig(nmax,4), SNR(4)
	real*4 Rcoadd(4),cmag(19,4),cemag(19,4), coApCorr(4)
	real*4 mLogQ(nmax,4),WPROmagR(nmax,4,3), cov_ave(5), cov_std(5)
c	real*4 slope(nmax,4), eslope(nmax,4), slopechi2(nmax,4)
	integer*4 Rho12(nmax), Rho23(nmax), Rho34(nmax), P12(nmax),P23(nmax), P34(nmax)
        integer*4 IDlist(nmax)                                      ! JWF B31206
c	real*4 Fmin(nmax,4), Fmax(nmax,4)
	real*4 KS(nmax,4)       !  CJG B30319
        real*4 value, dtor      !  JWF B30107
        real*8 dRA, dDec, JD0   !  JWF B30107
        data dtor/0.0174533/    !  JWF B30107
        real*4 snr_pm(4), mpro_pm(4), sigmpro_pm(4), rchi_pm(4)  ! JWF B30211

	real*4 mx(4),my(4),mmag(4),msig(4),mchi(4),mflux(4),munc(4)  ! CJG B30219

	real*8 jdtable(nmax,nfmep)                                  ! CJG B30213
	real*4 meptable(nmax,20,nfmep)                              ! CJG B30213
	character*9 namtable(nmax,nfmep)                            ! TPC
	logical mepwrite, wrframe(nfmep)                            ! CJG B30213

c	integer*4 mwrote,msources,mfiles                            ! CJG B30225
	integer*4 mfiles                                            ! CJG B30225
        integer*4 dummy_count                                       ! TPC

        Real*4       realnan                    ! TPC
        Integer*4    irealnan                   ! TPC
        Data         irealnan/-1/               ! TPC  
        equivalence (realnan, irealnan)         ! TPC

	real*4 Fcorr(4), CoaddXY(2,4)
	integer kill(4), galflag (4)
   	integer nout,nmax,NumMer(19,4), FLG(19,4)
	integer coflag(19,4)
	integer point,anyxy,nwframes                                 ! CJG B30228
	integer*2 temp(4), wflag (nf,4), N_Mstat(nmax,4), M_Mstat(nmax,4), SatNum(nmax,4), nDF(nmax,4)

	character*2 sNaper
	character*(*) outf,hname, imID
	character*(*) meproot                                       ! CJG B30228
	character*3 NaN                                             ! CJG B30228
	character*50 formula, srcID
	character*24 MEPsrcID     ! TPC
	character*25 s0,s1,format
	character*3000 header, BigString
        character*13 str13      !  JWF B30107
        character*12 str12      !  JWF B31206
        character*10 str10      !  JWF B21012
        character*9  str9       !  JWF B21012
        character*7  str7       !  JWF B30129
        character*3  str3       !  JWF B31207
        character*420 AllNull   !  JWF B30221 - put in DATA to make static
        character*283 AllNullpc !  JWF B31206 - post-cryo version
	character*9 framename
        data   AllNull/' for active deblends not kept in PM solution'/ ! JWF B30211
        data AllNullpc/' for active deblends not kept in PM solution'/ ! JWF B31206
	logical debug, IzBad, writePOS, smode, verbose, SPIT
c
        include 'jwfcom.f'             ! JWF B21219
        integer*4 Itemp, numbands      ! JWF B30129/B31205
c
	dummy_count = 501938894 ! TPC - dummy count for top of MEP files
	writePOS = .false.
	debug = .false.
c
        if (postcryo) then             ! JWF B31205
          numbands = 2                 ! JWF B31205
        else                           ! JWF B31205
          numbands = 4                 ! JWF B31205
        end if                         ! JWF B31205
c
c this sets the header length
	Lmax = 9999
	if (smode) Lmax = 1630


c	write (6,*) 'WRITE OUT DEBUG ',SingleFrame,smode,postcryo

c==========================start of code added by CJG B30213============================

c open/close multi-epoch output files...

	
	if(mepwrite) then 

	   if(next_file) then 

	   
	      close(unit=munit)
	      open(unit=munit,file=meproot,status='old', 
     1             form='unformatted', recl=1, access='DIRECT')
	      write(munit,rec=2)msources
c	      write(88,*)msources
	      close(unit=munit)

	      write(6,*)'setting counters to zero for next mep file.'
	      mfiles=mfiles+1
	      msources=0
	      mwrote=0
	      next_file=.false.


	      if(jsrc.le.nout)then

		 L = numchar(meproot)
		 point=L-1
		 write(6,'(''Closing mepfile '',a, '' with source '',i7)')
     1	             meproot(1:L),jsrc-1
		 do i=1,L
		    if(meproot(i:i).eq.'.') point=i
		 enddo
		 read(meproot(point-3:point-1),'(i3)')ithfile
		 ithfile=ithfile+1
		 meproot(point-3:point-1) = '000'
		 if (ithfile.lt.10) then
		    write(meproot(point-1:point-1),'(i1)')ithfile
		 elseif (ithfile.ge.10.and.ithfile.lt.100) then
		    write(meproot(point-2:point-1),'(i2)')ithfile
		 elseif (ithfile.ge.100) then 
		    write(meproot(point-3:point-1),'(i3)')ithfile
		 endif
		 
		 if(verbose) write(6,'(''Opening mepfile '',a,'' with source '',i7)')
     1	             meproot(1:L),jsrc
		 
		 open (unit=munit,file=meproot,status='unknown',
     1	             form='unformatted',access='SEQUENTIAL',recl=33)	
		 write(unit=munit) dummy_count   ! TPC
c		 write(88,*)dummy_count
	      endif

	   endif
	endif
	   
c==========================end of code added by CJG B30213============================



c  header biz

	do ib=1,4
	 temp(ib) = 0
	 if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then
		temp(ib) = 1
	 endif

	 do Map = 1, Naper
          Rap (Map,ib) = cscale(ib) * R(Map,ib)
          enddo
	enddo

	if (Jsrc.eq.1) then
c
          AllNull   = '     null         null        null      null      null'
     +          //'       null       null      null     null     null     null'
     +          //'     null     null     null      null         null'
     +          //'          null        null          null        null'
     +          //'          null        null         null       null'
     +          //'         null     null        null         null     null'
     +          //'         null        null      null        null      null'
     +          //'         null     null    null        null'
c
          AllNullpc = '     null          null        null      null     null'
     +          //'        null       null      null     null     null     null'
     +          //'     null      null         null          null        null'
     +          //'         null       null         null     null        null'
     +          //'         null       null    null      null      null'

	  iverify = 1
	  L = numchar(outf)
	  if (verbose) write (6,'(a,1x,a)') 'write table: ',outf(1:L)
c	  if (verbose) write (6,'(a,1x,a)') 'write meptable: ',meproot  ! CJG B30213

	  if(mepwrite) then	! CJG B30213
             if (verbose) write (6,'(a,1x,a)') 'write meptable: ',meproot  ! JWF B31122
	     open (unit=munit,file=meproot,status='unknown', ! CJG B30228
     1          form='unformatted',access='SEQUENTIAL',recl=33) ! CJG B30228
	     write(unit=munit) dummy_count   ! TPC
c	     write(88,*)dummy_count
	     write(6,*)'setting counters to zero at first open.' ! CJG B30228
	     mwrote=0		! CJG B30228
	     msources=0		! CJG B30228
	     next_file=.false.	! CJG B30228

	  endif			! CJG B30228

	  open (unit=iunit,file=outf)

	  write (iunit,'(a,i7)') '\Nsrc = ',nout

	  write (iunit,'(a,i5)') '\ frames engaged: ',nf
	  write (iunit,'(a,4i3)') '\ bands engaged: ',(temp(ib),ib=1,4)
	  write (iunit,'(a,4f7.3)') '\ zero mags(band): ',(zero(ib),ib=1,4)

	  do ib0=1,numbands

	  if (temp(ib0).eq.1) then


	  write (iunit,'(a,i2,a,f7.2,a,2x,f7.2,a)') '\ band = ',ib0,'  standard Rap(band) = ',
     1     Rap(1,ib0),' arcsec , ',R(1,ib0),' pix'

	   endif

	  enddo

	  do ib0=1,numbands
	   if (debug) junit = 100 + ib0
	   if (temp(ib0).eq.1) then

	  if (debug) then
	    write (junit,'(a)') '|   ra    |    dec  | WPRO  | ePRO  |  rchi | nap |  WAPP |  eWAPP| flg| fcorr |'
	  endif

c	  write (6,*) Naper,Naper-1

	  if (Naper-1.lt.10) then
	    write (sNaper,'(i1)') Naper-1
	  else
	    write (sNaper,'(i2)') Naper-1
	  endif

	MMM = numchar(sNaper)
c	write (6,*) sNaper(1:MMM)

	  formula = '(a,i2,a,' // sNaper(1:MMM) // 'f7.2,a,2x,' // sNaper(1:MMM) // 'f7.2,a)'
c	  formula = '(a,5f7.2,a,2x,5f7.2,a)'
          LF = numchar(formula)

c	write (6,'(a)') formula(1:lf)

	   write (iunit,formula(1:lf)) '\ band = ',ib0,'  circ apertures Rap = ',
     1      (Rap(Map,ib0),Map=2,Naper),' arcsec , ',(R(Map,ib0),Map=2,Naper),' coadd pix'

	  endif

	  enddo

          write (iunit,'(a,f12.6)') '\ MJD0 = ',MJD0   !  JWF B21206

	  open (unit=11,file=hname)
	  do jj=1,999
	   read (11,'(a)',end=910) header
	   L = numchar(header)
	   L = min (L,Lmax)
	   write (iunit,'(a)') header(1:L)
	  enddo
 910	close (11)

	endif



cccccccccccccccccccccccccccccccccccccc


	BigString = ''
	LBig = 0

	write(srcID,'(a,a,a,i6.6)') ' ',imID(1:lnblnk(imID)), '-',Jsrc
	Lsrc = numchar(srcID)


! 'source id jazz ',Jsrc

	L1 = LBIG + 1
        L2 = L1 + (Lsrc-1)
        LBIG = L2
        BigString (L1:L2) = srcID(1:Lsrc)

        Lgap = 25 - Lsrc

	s0 = '                         '
        L1 = LBIG + 1
        L2 = L1 + (Lgap-1)
        LBIG = L2
	L = 1 + (L2-l1)
        BigString (L1:L2) = s0(1:L)

	ivalue = Jsrc
        format = 'i7'
        call int_stringform (ivalue, format, s0, L)
	L1 = LBIG + 1
	L2 = L1 + (L-1)
	LBIG = L2
	BigString (L1:L2) = s0(1:L)


c   coordinates

	dvalue = RAlist(Jsrc)   !  RA

	!  check that RA is not out of bounds; if so, wrap about 360 degrees
	if (dvalue.gt.360.d0) then
		dvalue = dvalue - 360.
	else if (dvalue.lt.0.) then
		dvalue = 360.d0 + dvalue
	endif

	format = 'f12.7'
	call double_stringform (dvalue, format, s0, L)
	L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

 	dvalue = Declist(Jsrc)  ! DEC
	!  check that DEC is not out of bounds; if so, wrap about +-90 degrees
        if (dvalue.gt.90.d0) then
                dvalue = 180.d0 - dvalue    !  JWF B21015
        else if (dvalue.lt.-90.d0) then
                dvalue = -(180. + dvalue)
        endif

	call double_stringform (dvalue, format, s0, L)
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

c==========================start of code added by CJG B30221============================

c       determine how many and which frames to write out...

	if(mepwrite)then

	   nwframes=0
	   do ii=1,nfmep
	      anyxy=0
	      wrframe(ii)=.false.
	      do ib=1,4          ! JWF B31205: should this 4 be changed to numbands?

		 ibloc=1+(ib-1)*2
		 valuex = meptable(jsrc,ibloc,ii)
		 valuey = meptable(jsrc,ibloc+1,ii)
		 if (valuex.gt.0 .and. valuey.gt.0) anyxy=anyxy+1

	      enddo

	      if(anyxy.gt.0) then
		 nwframes=nwframes+1
		 wrframe(ii)=.true.
	      endif
c       write(6,*)jsrc,ii,value,anyxy,nwframes,wrframe(ii)
	   enddo


c       write source identifier line to meptable              ! CJG B30213

	   MEPsrcID = srcID    ! TPC
	   write(munit)MEPsrcID,jsrc,ralist(jsrc),declist(jsrc),nwframes ! CJG, TPC
c	   write(88,*)MEPsrcID,jsrc,ralist(jsrc),declist(jsrc),nwframes ! CJG, TPC
	   mwrote=mwrote+1
	   msources=msources+1

c       now go through all epochs

	   nmep = 0

	   do ii=1,nfmep

	      if(wrframe(ii)) then

	         nmep = nmep + 1

		 do ib=1,4       ! JWF B31205: should this 4 be changed to numbands?

		    kill(ib) = 0
		 
		    MAGWPRO (ib) = 99. ! null value
		    eMAGWPRO (ib) = 0. ! null value
		    rchi(ib) = 0.
		    fsignal(ib) = 0.
		    ffnoise(ib) = 0.
		    
		    ibloc=1+(ib-1)*2
		    xfpos(ib) = meptable(jsrc,ibloc,nmep)
		    yfpos(ib) = meptable(jsrc,ibloc+1,nmep)
		    
		    ibloc=9+(ib-1)*3
		    fsignal(ib)=meptable(jsrc,ibloc,nmep)
		    ffnoise(ib)=meptable(jsrc,ibloc+1,nmep)
		    rchi(ib)=meptable(jsrc,ibloc+2,nmep)
		    
		    if (ffnoise(ib).le.0.)  then 
		       kill(ib) = 1
		    else
		       Thresh = 2. * ffnoise(ib) ! upper limit, if needed
		       if (fsignal(ib).lt.Thresh) then ! compute upper limit
			  if (fsignal(ib).lt.0.) then
			     upper_lim = Thresh
			  else
			     upper_lim = fsignal(ib) + thresh
			  endif
			  fsignal(ib) = upper_lim ! set the signal to the upper limit value
			  ffnoise(ib) = 0. ! null the flux noise for upper limits
		       endif
		    endif
		 
		    if (fsignal(ib).gt.0.) then
		    
		       MAGWPRO (ib) = zero(ib) - (2.5*log10(fsignal(ib)))
		       if (ffnoise(ib).gt.0.) then
			  eMAGWPRO (ib) = (1.0857 * ffnoise(ib))  / fsignal(ib)
		       else	! upper limit
			  eMAGWPRO (ib) = 9.99 ! indicates upper limit
		       endif
		       
		    endif
		    
		    if (MAGWPRO (ib).lt.-9.99) MAGWPRO (ib) = -9.99
		    if (MAGWPRO (ib).gt.30.) MAGWPRO (ib) = 99.

c                   Reset to original values before the upper limit computation -- TPC
		    fsignal(ib)=meptable(jsrc,ibloc,nmep)  ! TPC
		    ffnoise(ib)=meptable(jsrc,ibloc+1,nmep)
		    
		 enddo ! Band
	      
		 framename = namtable(Jsrc,nmep)
		 Lname = numchar(framename)
	      
		 
c       append frame/band specific frame positions....
		 
		 do ib=1,4          ! JWF B31205: should this 4 be changed to numbands?
		    
		    inull = 0
		    value = xfpos(ib)
		    mx(ib)=value
		    value = min (value, 9999.999)
		    if (value.le.0.) then
		       inull=1
		    endif
		    format = 'f9.3'
		    if (IzBad(value)) inull = 1
		 
		    if (inull.eq.1) then
		       mx(ib)=realnan                            ! CJG B30219
		    endif
		    
		    value = yfpos(ib)
		    my(ib)=value
		    value = min (value, 9999.999)
		    if (value.le.0.) then
		       inull=1
		    endif
		    if (IzBad(value)) inull = 1
		 
		    if (inull.eq.1) then
		       my(ib)=realnan                            ! CJG B30219
		    endif
		 enddo ! band
	      
	      
c       append magnitudes, errors, and rchis

c       write(6,*)'magwpro, ib, kill(ib)'
	      
		 do ib=1,4          ! JWF B31205: should this 4 be changed to numbands?
		 
		    inull=0
		    value = MAGWPRO(ib)
		    mmag(ib)=value
c       write(6,*)value, ib, kill(ib)
		    if (value.gt.30.) inull=1
		    if (value.lt.-30.) inull=1
		    if (IzBad(value)) inull = 1
		    
		    if ((kill(ib).eq.1).or.(inull.eq.1)) then
		       mmag(ib) = realnan                       ! CJG B30219
		       rchi(ib) = 0.
		    endif
		    
		    value = eMAGWPRO(ib)
		    msig(ib)=value
c       write(6,*)'emagwpro: ',value
		    if (value.le.0.) inull=1
		    if (value.gt.5.) inull=1
		    if (IzBad(value)) inull = 1
		    
		    if ((kill(ib).eq.1).or.(inull.eq.1)) then
		       msig(ib) = realnan                       ! CJG B30219
		    endif
		    
		    inull = 0
		    value = rchi(ib)
		    mchi(ib)=value
c       write(6,*)'rchi(i): ',value
		    if (value.le.0.) inull = 1
		    if (value.gt.99999.) inull = 1
		    if (IzBad(value)) inull = 1
		    if (inull.eq.1) value = 0.
		    if ((kill(ib).eq.1).or.(inull.eq.1)) then
		       mchi(ib) = realnan                      ! CJG B30219
		    endif
		 
		 enddo ! band
	      
c       append frame fluxes and uncertainties...
		 
		 do ib=1,4      ! JWF B31205: should this 4 be changed to numbands?
c       write(6,*)'fsignal(ib), ffnoise(ib), ib: ',fsignal(ib),
c	1	      ffnoise(ib), ib

		    mflux(ib)=fsignal(ib)
		    
		    if ((kill(ib).eq.1).or.(temp(ib).eq.0)) then
		       mflux(ib) = realnan                     ! CJG B30219
		    endif
		    
		    munc(ib)=ffnoise(ib)
		    if (kill(ib).eq.1) then
		       munc(ib) = realnan                      ! CJG B30219
		    endif
		 enddo


		 write(munit)(mx(im),my(im),mmag(im),msig(im),mchi(im),
     1                 mflux(im),munc(im),im=1,4),jdtable(jsrc,nmep),
     1                 framename(1:9)     ! JWF B31205: should this 4 be changed to numbands?

c		 write(88,*)(mx(im),my(im),mmag(im),msig(im),mchi(im),
c     1                 mflux(im),munc(im),im=1,4),jdtable(jsrc,nmep),
c     1                 framename(1:9)    ! JWF B31205: should this 4 be changed to numbands?


		 mwrote=mwrote+1
		 if(mwrote.ge.max_mep_lines)next_file=.true.
c		 write(6,*)jsrc,msources,ii,nmep,mwrote,next_file
		 

	      endif ! wrframe(ii)
	      
	   enddo ! ii=1,nfmep

	   if(jsrc.eq.nout) then 
	      close(unit=munit)
	      open(unit=munit,file=meproot,status='old', 
     1          form='unformatted', recl=1, access='DIRECT')
	      write(munit,rec=2)msources
	      close(unit=munit)
	      mfiles=mfiles+1
	   endif



	endif   ! mepwrite

c========================end of code added by CJG B30213======================
	

c regular mdex output	   


	do ii=3,5

	  value = Table (Jsrc,ii)

	  inull = 0

	  if ((value.gt.99.999).or.(value.lt.-99.999)) then
		  inull = 1
		  value = 99.999
	  endif

	  if (IzBad(value)) inull = 1

	  format = 'f9.4'
	  call real_stringform (value, format, s0, L)

	  if (inull.eq.1) then
		call nullswap (s0,s1)
		s0 = s1(1:L)
	  endif

          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	enddo


c  pixel location
ccccc   only write out pixels when they are meaningful

	if ((ireg*jreg.eq.1).and.(nf.eq.1)) then
		writePOS = .true.
	endif



	do ib=1,numbands

	  value = 0.
	  inull = 0
	  if ((temp(ib).eq.1).and.(writePOS)) then
		inull = 0
c          else
c		inull = 1
c		value = 0.
	  endif


	  if (writePOS) then
		 value = XPos (Jsrc,1,ib)
	  else
		value = CoaddXY (1,ib)
	  endif

	  value = min (value, 9999.999)
	  if (value.le.0.) inull=1
	  format = 'f9.3'
	  if (IzBad(value)) inull = 1


          call real_stringform (value, format, s0, L)
	  if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	  if (writePOS) then
		value = YPos (Jsrc,1,ib)
	  else
                value = CoaddXY (2,ib)
          endif

	  value = min (value, 9999.999)
	  if (value.le.0.) inull=1
	  if (IzBad(value)) inull = 1

	  call real_stringform (value, format, s0, L)
          if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	enddo


c    median sky & sig
	do ib=1,numbands

	 kill(ib) = 0

c	write (6,*) ib,Jsrc,medsky(Jsrc,ib),medskysig(Jsrc,ib)

	  if (temp(ib).eq.1) then
                inull = 0
		
		if ((medsky(Jsrc,ib).lt.-999.).or.(medsky(Jsrc,ib).gt.BGmax(ib))) then
			inull=1
			kill(ib) = 1
		endif

		if ((medskysig(Jsrc,ib).le.0.).or.(medskysig(Jsrc,ib).gt.999.)) then
			inull=1
			kill(ib) = 1
                endif

		medsky(Jsrc,ib) = min (medsky(Jsrc,ib), BGmax(ib))
		medsky(Jsrc,ib) = max (medsky(Jsrc,ib), -99.)
		medskysig(Jsrc,ib) = min (medskysig(Jsrc,ib), 999.)
		medskysig(Jsrc,ib) = max (medskysig(Jsrc,ib), 0.)


		if ((formalsig(Jsrc,ib).lt.-999.).or.(formalsig(Jsrc,ib).gt.999.)) then
                        inull=1
                        kill(ib) = 1
                endif

		formalsig(Jsrc,ib) =  min (formalsig(Jsrc,ib), 999.)
                formalsig(Jsrc,ib) = max (formalsig(Jsrc,ib), 0.)


          else
                inull = 1
		kill(ib) = 1
          endif



	  value = medsky(Jsrc,ib)
	  format = 'f10.3'
          call real_stringform (value, format, s0, L)
	  if (IzBad(value)) inull = 1

          if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	  value = medskysig(Jsrc,ib)
	   format = 'f8.3'
	  call real_stringform (value, format, s0, L)
	  if (IzBad(value)) inull = 1

          if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	  value = formalsig(Jsrc,ib)
          format = 'f8.3'
          call real_stringform (value, format, s0, L)
          if (IzBad(value)) inull = 1

          if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	enddo

cccccccccccccccccccccccccccc                   ! JWF B31206
cc WPRO fitting radii
c
        do 500 ib = 1, numbands
          write(Str7,'(F7.2)') Table(Jsrc,59+ib)
          BigString = BigString(1:lnblnk(BigString))//Str7
500     continue
        LBIG = lnblnk(BigString)
c
cccccccccccccccccccccccccccc
cc  WPRO fluxes, SNRs, mags
	
	do 600 ib=1,numbands

	  MAGWPRO (ib) = 99.  ! null value
	  eMAGWPRO (ib) = 0.  ! null value
	  SNR(ib) = 0.  ! null value
	  rchi(ib) = 0.
	  kill(ib) = 0


	  if ((temp(ib).eq.0).or.(kill(ib).eq.1)) goto 600

	   signal = Table (Jsrc,5+ib)
	   fnoise = Table (Jsrc,9+ib)
	   rchi(ib) = Table (Jsrc,13+ib)

c	write (6,*) 'here ',ib,signal, fnoise, rchi(ib), kill (ib)

c	if (jsrc.eq.374) write (6,*) 'here ',ib,signal,fnoise,rchi(ib)

	   if (SPIT) then
		if (ib.le.3) then
		   fnoise = Table (Jsrc,9+ib) * sqrt(Fcorr(ib))   ! bump up the wpro noise to account for coadd correlated pixels
		else 
		   fnoise = Table (Jsrc,9+ib) * sqrt(Fcorr(ib))
		endif
	   endif

c	write (6,*) ib,signal,fnoise
c	   if ((fnoise.le.0.).or.(fnoise.gt.1.e8)) then ! problem with flux noise;  kill the source
	   if (fnoise.le.0.)  then ! problem with flux noise;  kill the source
		kill(ib) = 1
                goto 600
           endif

c	write (6,*) 'kill value ',ib,kill(ib)

	   SNR(ib) = 0.
	   if ((signal.gt.0.).and.(fnoise.gt.0.)) then
	     SNR(ib) = signal / fnoise
	     SNR (ib) = min (SNR(ib),9999.)

c		if (SNR(ib).le.0.05) then
c			write (6,'(i4,f12.4,f12.4,f12.4)')  ib,signal,fnoise,SNR(ib)
c		endif

	   endif

	   Thresh = 2. * fnoise  ! upper limit, if needed

	   if (signal.lt.Thresh) then  ! compute upper limit

	      if (signal.lt.0.) then
		 upper_lim = Thresh
	      else
		 upper_lim = signal + thresh
	      endif
		
	      signal = upper_lim   ! set the signal to the upper limit value
	      fnoise = 0.  ! null the flux noise for upper limits
	   endif

	   if (signal.le.0.) goto 600   ! all else fails, give up on source flux

           MAGWPRO (ib) = zero(ib) - (2.5*log10(signal))

c	if (jsrc.eq.374) write (6,*) 'there ',ib,zero(ib),signal,MAGWPRO (ib)


	   if (fnoise.gt.0.) then

	     eMAGWPRO (ib) = (1.0857 * fnoise)  / signal 
c	     eMAGWPRO (ib) = max (eMAGWPRO (ib), 0.0)
c	     eMAGWPRO (ib) = min (eMAGWPRO (ib), 9.99)

	   else  ! upper limit

	     eMAGWPRO (ib) = 9.99  ! indicates upper limit

	   endif

	  if (MAGWPRO (ib).lt.-9.99) MAGWPRO (ib) = -9.99
c	  if (eMAGWPRO (ib).lt.-9.99) eMAGWPRO (ib) = 0.0
	  if (MAGWPRO (ib).gt.30.) MAGWPRO (ib) = 99.

 600	continue


! SNR
	do ib=1,numbands
	
	 inull = 0
	 if (SNR(ib).le.0.) then
		inull=1
		SNR(ib) = 0.
	 endif

	 value = SNR(ib)
          format = 'f7.1'
          call real_stringform (value, format, s0, L)
	  if (IzBad(value)) inull = 1

          if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	enddo



! raw Fluxes (WPRO)

	format = '1pe12.4'
	do ib=1,numbands
	   call real_stringform (Table(Jsrc,ib+5), format, s0, L)

	   if ((kill(ib).eq.1).or.(temp(ib).eq.0)) then
		call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
	   L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   s0 = '            '
	   L1 = LBIG + 1
	   L2 = L1 + 1
	   LBIG = L2
           BigString (L1:L2) = s0(1:L)


	   fnoise = Table (Jsrc,9+ib)
	   if (SPIT) then
		fnoise = Table (Jsrc,9+ib) * sqrt(Fcorr(ib)) ! bump up the wpro noise to account for coadd correlated pixels
           endif

	   call real_stringform (fnoise, format, s0, L)
	   if (kill(ib).eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
	   L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)
	enddo


! WPRO mags

	do ib=1,numbands

	   inull=0
	   value = MAGWPRO(ib)
	   if (value.gt.30.) inull=1
	   if (value.lt.-30.) inull=1
	   if (IzBad(value)) inull = 1
	   if (temp(ib).eq.0) inull = 1

	   format = 'f7.3'
	   call real_stringform (value, format, s0, L)
           if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
		rchi(ib) = 0.
           endif
	   L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   s0 = '            '
	   L1 = LBIG + 1
           L2 = L1 + 2
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   value = eMAGWPRO(ib)
	   if (value.le.0.) inull=1
	   if (value.gt.5.) inull=1
	   if (IzBad(value)) inull = 1

	   call real_stringform (value,format, s0, L)
           if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
           L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   inull = 0
	   value = rchi(ib)
	   if (value.le.0.) inull = 1
	   if (IzBad(value)) inull = 1
	   format = '1pe11.3'
	   if (inull.eq.1) value = 0.
	   call real_stringform (value,format, s0, L)
           if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
           L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	enddo



! cummulative rCHI

	inull = 0
	value = Table(Jsrc,18)
	if (value.le.0.) inull=1
	if (value.gt.99999.) inull=1
	if (IzBad(value)) inull = 1
	if (inull.eq.1) value = 0.

	format = '1pe11.3'
	call real_stringform (value,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

! nblend, active deblend

	do ii=19,20
	 ivalue = int(Table(Jsrc,ii))
	 format = 'i4'
	 call int_stringform (ivalue,format, s0, L)
	 L1 = LBIG + 1
         L2 = L1 + (L-1)
         LBIG = L2
         BigString (L1:L2) = s0(1:L)
	enddo


	
! latent fraction (aka troublesome pixels/bits)
c	do ii=21,24                      ! JWF B31207
	do ii=21,20+numbands             ! JWF B31207
		inull=0
		value = Table(Jsrc,ii)
		ivalue = value
		if (IzBad(value)) inull = 1

		if ((value.lt.0.).or.(value.gt.9999.00)) then
			inull=1
			value = 0.
		endif
		format ='i8.4'   ! [TPC] modified to emphasize non-fractional values for this flag
		call int_stringform (ivalue,format, s0, L)
		if (inull.eq.1) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
        	endif

		L1 = LBIG + 1
        	L2 = L1 + (L-1)
        	LBIG = L2
        	BigString (L1:L2) = s0(1:L)
	enddo

! saturation fraction
!               Table(n,25:28)  =       fraction of saturated pixels

c	do ii=25,28                      ! JWF B31207
	do ii=25,24+numbands             ! JWF B31207
                inull=0
                value = Table(Jsrc,ii)
                if (IzBad(value)) inull = 1

                if ((value.lt.0.).or.(value.gt.1.00)) then
                        inull=1
                        value = 0.
                endif
                format ='f8.5'
                call real_stringform (value,format, s0, L)
                if (inull.eq.1) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif

                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)
        enddo

ccccccccc
c  saturation level flag

        write (s0,'(a,4i1)') '   ',(SatNum (Jsrc,ib),ib=1,4)
        L = 7
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 5'

!  COADD standard aperture
!  MAGWPRO (ib) = zero(ib) - (2.5*log10(signal))

	do MAP = 1, Naper
        do ib=1,numbands

	  zmag = cmag(Map,ib)
	  flux_error  = cemag(Map,ib) / 1.0857


	  if (flux_error.gt.0.) then

		S_N = 1./flux_error

		if ((S_N.lt.2.0).or.(coflag (Map,ib).eq.32)) then  ! derive upper limit


		  Tthresh =  medskysig(Jsrc,ib) * 2.

		  if (coflag (Map,ib).eq.32) then
			if (zmag.lt.30.) then
			 signal = 10**((zero(ib)-zmag)/2.5)
			 thresh = signal
			 upper_lim = thresh
			else
			  signal = 0.
                          upper_lim = Tthresh
			endif

	          else if (zmag.lt.30.) then

		    signal = 10**((zero(ib)-zmag)/2.5)
		    znoise = signal * flux_error
		    thresh = 2. * znoise
		    upper_lim = signal + thresh

		  else

		    signal = Tthresh
		    thresh = Tthresh
		    upper_lim = Tthresh
		    if (flux_error.lt.5.) then
		      znoise = signal * flux_error
		      thresh = 2. * znoise
		      upper_lim = thresh
		    endif

		  endif

		  zmag = zero(ib) - (2.5*log10(upper_lim))

		   coflag (Map,ib) = 32  ! upper limit
                   cmag(Map,ib) = zmag
                   cemag(Map,ib) = 0.


		endif  ! upper limit condition
	    endif

	enddo
	enddo

c	write (6,*) 'DEBUG 6'
c standard ap
	Map = 1
	do ib=1,numbands

	  inull = 0
	  value = cmag(Map,ib)
          if (value.gt.30.) inull=1
          if (value.lt.-10.) inull = 1
	  if (IzBad(value)) inull = 1


c	  ivalue = NumMer(Map,ib)
	  ivalue  = M_Mstat(Jsrc,ib)

c	write (6,*) 'DEBUG 6b'

	  format = 'i4'

c	write (6,*) 'DEBUG 6c'

	  value = cmag(Map,ib)
          if (value.gt.30.) inull=1
          if (value.lt.-10.) inull = 1
          if (IzBad(value)) inull = 1
	  format = 'f7.3'
	  call real_stringform (value,format, s0, L)
	  if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 6d'

	  value = cemag(Map,ib)
          if (value.gt.10.) inull=1
          if ((value.le.0.).or.(coflag(Map,ib).eq.32)) inull = 1
	  if (IzBad(value)) inull = 1

	  call real_stringform (value,format, s0, L)
          if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 6e'

	  ivalue = coflag(Map,ib)
	  if (ivalue.eq.32) inull = 0
	  format = 'i6'
          call int_stringform (ivalue,format, s0, L)
	  if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 6f'

	  value = cov_ave(ib)
          if (value.gt.99999.) inull=1   !  thj 07Jan2018
          if (value.le.0.) inull = 1
          if (IzBad(value)) inull = 1

          format = 'f8.3'
	  if (value.gt.999.) format = 'f8.2'
	  if (value.gt.9999.) format = 'f8.1'
          call real_stringform (value,format, s0, L)
          if (inull.eq.1) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 6g'

	  value = coApCorr(ib)
	  format = 'f7.3'
          call real_stringform (value,format, s0, L)
	  if (IzBad(value)) inull = 1

          if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 6h'


	if (.not.SPIT) then
cc compare Standard Ap with WPRO
c cmag(Map,ib)  ceag(Map,ib) coflag(Map,ib)
c MAGWPRO (ib) eMAGWPRO (ib)   Table(n,25:28) fraction of saturated pixels

c	 value = -9.
c	 zTable (Jsrc, ib ) = value
	 igo = 0

c	write (6,*) 'DEBUG 6i'

c	 if (Table(Jsrc,19).le.1.) igo=igo+1    ! nblend
c	 if ((MAGWPRO (ib) .lt. 25.) .and. (cmag(Map,ib) .lt. 25.) ) igo=igo+1    ! mags
c	 if (Table (Jsrc,24+ib).le.0.) igo = igo + 1  ! Saturation
c	 if ( (eMAGWPRO (ib) .gt. 0.0) .and. (eMAGWPRO (ib) .le. 0.025) )  igo=igo+1    ! WPRO SNR
c	 if ( (cemag(Map,ib) .gt. 0.0) .and. (cemag(Map,ib) .le. 0.025) )  igo=igo+1    ! WAPPco SNR
c	 if ( coflag (Map,ib) .le. 0 ) igo=igo+1    ! WAPPco flag
c
c	 if (igo.eq.6) then
c		value = MAGWPRO (ib) - cmag(Map,ib)
c		zTable (Jsrc, ib ) = value
c	 endif

	endif

c	if ((ib.eq.1).and.(igo.eq.6)) write (85,*) Jsrc,MAGWPRO (ib),cmag(Map,ib), value

	enddo


c	write (6,'(a)') BigString (1:LBIG)

!  set of circular apertures

c	write (6,*) 'DEBUG 7'

	do MAP = 2, Naper

	do ib=1,numbands

	  inull = 0
	  value = cmag(Map,ib)
          if (value.gt.30.) inull=1
          if (value.lt.-10.) inull = 1
	  if (IzBad(value)) inull = 1

	   s0 = '            '
           L1 = LBIG + 1
           L2 = L1 + 2
           LBIG = L2
           BigString (L1:L2) = s0(1:L)


	  format = 'f7.3'
	  call real_stringform (value,format, s0, L)
	  if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	   s0 = '            '
           L1 = LBIG + 1
           L2 = L1 + 2
           LBIG = L2
           BigString (L1:L2) = s0(1:L)


	  value = cemag(Map,ib)
          if (value.gt.10.) inull=1
	  if ((value.le.0.).or.(coflag(Map,ib).eq.32)) inull = 1
	  if (IzBad(value)) inull = 1

	  call real_stringform (value,format, s0, L)
          if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	   s0 = '            '
           L1 = LBIG + 1
           L2 = L1 + 2
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	  ivalue = coflag(Map,ib)
	  if (ivalue.eq.32) inull = 0
	  format = 'i5'
          call int_stringform (ivalue,format, s0, L)
	  if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)


	enddo

	enddo  ! MAP


	if (smode) goto 90  ! don't write the N/M stats or the galaxy biz

        if (SingleFrame) go to 70

! N_M stats

c	write (6,*) 'DEBUG 8'

        do ib=1,numbands
		inull= 0

		ivalue = N_Mstat(Jsrc,ib)
		if ( (ivalue.lt.0).or.(ivalue.gt.9999) ) inull=1
		format = 'i7'
		call int_stringform (ivalue,format, s0, L)
		if ((kill(ib).eq.1).or.(inull.eq.1)) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
          	endif
          	L1 = LBIG + 1
          	L2 = L1 + (L-1)
          	LBIG = L2
          	BigString (L1:L2) = s0(1:L)

		ivalue = M_Mstat(Jsrc,ib)
		if ( (ivalue.lt.0).or.(ivalue.gt.9999) ) inull=1
                format = 'i6'
                call int_stringform (ivalue,format, s0, L)
                if ((kill(ib).eq.1).or.(inull.eq.1)) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

c WPRO source repeatability

		value = WPROmagR (Jsrc,ib,1)
		if (value.gt.30.) inull=1
                if (value.le.-10.) inull = 1
                if (IzBad(value)) inull = 1

		format = 'f8.3'
                call real_stringform (value,format, s0, L)
		if (inull.eq.1) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		value = WPROmagR (Jsrc,ib,2) ! delta-mag
                if (value.gt.20.) inull=1
                if (value.le.0.) inull = 1
                if (IzBad(value)) inull = 1
                format = 'f8.3'
                call real_stringform (value,format, s0, L)
		if (inull.eq.1) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		value = WPROmagR (Jsrc,ib,3) ! delta-mag / sqrt (n)
                if (value.gt.10.) inull=1
                if (value.le.0.) inull = 1
                if (IzBad(value)) inull = 1
                format = 'f8.3'
                call real_stringform (value,format, s0, L)
		if (inull.eq.1) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

ccccc new variability stuff
		inull = 0
		value = KS(Jsrc,ib)  !  Stetson K index
		if ((value.lt.0.).or.(value.ge.99.)) inull = 1
		if (IzBad(value)) inull = 1
c		write(6,*)'ks,inull: ',value,inull
		format = 'f12.5'
                call real_stringform (value,format, s0, L)
                if (inull.eq.1) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
		L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		inull = 0
		ivalue = nDF(Jsrc,ib)
                if ((ivalue.lt.0).or.(ivalue.gt.9999)) inull = 1
                format = 'i6'
                call int_stringform (ivalue,format, s0, L)
                if (inull.eq.1) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		value = mLogQ(Jsrc,ib)
		if ((value.lt.0).or.(value.gt.99.99)) inull = 1
                format = 'f6.2'
                call real_stringform (value,format, s0, L)
                if (inull.eq.1) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

cc   min/max & mean mJD
	        dvalue = mJD_min(Jsrc,ib)
                format = 'f18.8'
                call double_stringform (dvalue, format, s0, L)
                if (dvalue .le. 0.d0) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		dvalue = mJD_max(Jsrc,ib)
                format = 'f18.8'
                call double_stringform (dvalue, format, s0, L)
                if (dvalue .le. 0.d0) then
                  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		dvalue = AVE_mJD(Jsrc,ib)
c	write (6,*) dvalue
		format = 'f18.8'
        	call double_stringform (dvalue, format, s0, L)
		if (dvalue .le. 0.d0) then
		  call nullswap (s0,s1)
                  s0 = s1(1:L)
                endif
        	L1 = LBIG + 1
        	L2 = L1 + (L-1)
        	LBIG = L2
        	BigString (L1:L2) = s0(1:L)

	enddo

c variability correlation coefficients
	inull = 0
        ivalue = Rho12(Jsrc)  
        if ((ivalue.lt.-99).or.(ivalue.gt.99)) inull = 1
        format = 'i6'
        call int_stringform (ivalue,format, s0, L)
        if (inull.eq.1) then
          call nullswap (s0,s1)
          s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)
        if (postcryo) go to 50  ! JWF B40107

	inull = 0
        ivalue = Rho23(Jsrc)
        if ((ivalue.lt.-99).or.(ivalue.gt.99)) inull = 1
	call int_stringform (ivalue,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

	inull = 0
        ivalue = Rho34(Jsrc)
        if ((ivalue.lt.-99).or.(ivalue.gt.99)) inull = 1
        call int_stringform (ivalue,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

50	inull = 0
        ivalue = P12(Jsrc)   !  NOTE:  this is really the Q value
        if ((ivalue.lt.0).or.(ivalue.gt.9)) inull = 1
        call int_stringform (ivalue,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)
        if (postcryo) go to 70  ! JWF B40107

	inull = 0
        ivalue = P23(Jsrc) !  NOTE:  this is really the Q value
        if ((ivalue.lt.0).or.(ivalue.gt.9)) inull = 1
        call int_stringform (ivalue,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

	inull = 0
        ivalue = P34(Jsrc) !  NOTE:  this is really the Q value
        if ((ivalue.lt.0).or.(ivalue.gt.99)) inull = 1
        call int_stringform (ivalue,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

c	write (6,*) 'DEBUG 9'

c Galaxy photometry
c XSCprox(nout), galRsemi(4),galba(4),galpa(4),galmag(4),galerr(4)

c	write (6,*) 'here'

70	inull = 0
	if (ngal.le.0) then
	 	inull = 1
		value = -1.
	else
		value = XSCprox(Jsrc)
	endif

	format = 'f8.2'
	call real_stringform (value,format, s0, L)
        if ((inull.eq.1).or.(value.lt.0.)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
		inull = 1
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)

	do ib=1,numbands
                value = galRsemi(ib)
                format = 'f8.2'
                call real_stringform (value,format, s0, L)
		
   		if ((inull.eq.1).or.(value.le.0.)) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		value = galba(ib)
                format = 'f7.3'
                call real_stringform (value,format, s0, L)

                if ((inull.eq.1).or.(value.le.0.)) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		value = galpa(ib)
                format = 'f7.1'
                call real_stringform (value,format, s0, L)

                if ((inull.eq.1).or.(galba(ib).le.0.)) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)


		value = galmag(ib)
                format = 'f8.3'
                call real_stringform (value,format, s0, L)

                if ((inull.eq.1).or.(value.gt.30.)) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)

		value = galerr(ib)
                format = 'f8.3'
                call real_stringform (value,format, s0, L)

                if ((inull.eq.1).or.(galmag(ib).gt.30.)) then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)


		ivalue = galflag(ib)
                format = 'i7'
                call int_stringform (ivalue,format, s0, L)
                if ((inull.eq.1).or.(galmag(ib).gt.30.))  then
                        call nullswap (s0,s1)
                        s0 = s1(1:L)
                endif
                L1 = LBIG + 1
                L2 = L1 + (L-1)
                LBIG = L2
                BigString (L1:L2) = s0(1:L)


	enddo

ccccccccccc

c	write (6,*) 'DEBUG 10'

c================== start of code added by JWF B30129 =====================
c
        Itemp = NInt(Table(Jsrc,55))             ! nIters (nonPM)
        if (Itemp .gt. 999999) Itemp = 999999
        write(Str7,'(I7)') Itemp
        BigString = BigString(1:lnblnk(BigString))//Str7
        Itemp = NInt(Table(Jsrc,56))             ! nSteps (nonPM)
        if (Itemp .gt. 999999) Itemp = 999999
        write(Str7,'(I7)') Itemp
        BigString = BigString(1:lnblnk(BigString))//Str7
c                                                           ! MDET info
        write(Str7,'(I7)') IDlist(Jsrc)          ! mdetID
        BigString = BigString(1:lnblnk(BigString))//Str7
        write(Str10,'(F10.5)') Table(Jsrc,64)    ! p(nc)
        BigString = BigString(1:lnblnk(BigString))//Str10
        write(Str10,'(F10.5)') Table(Jsrc,65)    ! p(mblend+nc)
        BigString = BigString(1:lnblnk(BigString))//Str10
        LBIG = lnblnk(BigString)
c
        if (SingleFrame) go to 80                ! JWF B31122
c                                                           ! MeanObsMJD

c	write (6,*) 'DEBUG in WRITETABLE MeanObsMJD'
c	write (6,*) ' Jsrc ',Jsrc,Table(Jsrc,33),Table(Jsrc,34)

        if (Table(Jsrc,33) .gt. 0.0) then
          JD0 = dble(Table(Jsrc,33)) + dble(Table(Jsrc,34)) !  JWF B30107
          write (str13,'(f13.6)') JD0   !  JWF B30107
        else
          if (postcryo) then
            BigString = BigString(1:lnblnk(BigString))//AllNullpc
          else
            BigString = BigString(1:lnblnk(BigString))//AllNull
          end if
          go to 90
        end if
        BigString = BigString(1:lnblnk(BigString))//str13   ! JWF B30107
        LBig = lnblnk(BigString)
c
        if (IzBad(Table(Jsrc,29))) then                    ! ra_pm
          dRA = 0.0
        else
c                                 1314900.0 = 3600*365.25
          value  = Table(Jsrc,29)/1314900.0     !  PMRA [deg/day] JWF B30122
          dRA    = value*(MJD0-JD0)/cos(DEClist(Jsrc)*dtor) ! JWF B30107 [obs-epoch Dec good enough]
        end if
c       dvalue = RAlist(Jsrc) + dRA             ! JWF B30107 [add offset to ref epoch]
        dvalue = 0.01d0*(dble(Table(Jsrc,38))+dble(Table(Jsrc,39))) + dRA

	!  check that RA is not out of bounds; if so, wrap about 360 degrees
	if (dvalue.gt.360.d0) then
		dvalue = dvalue - 360.
	else if (dvalue.lt.0.) then
		dvalue = 360.d0 + dvalue
	endif

	format = 'f12.7'
	call double_stringform (dvalue, format, s0, L)
	L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)
c
        if (IzBad(Table(Jsrc,30))) then                    ! dec_pm
          dDec = 0.0
        else
          value  = Table(Jsrc,30)/1314900.0   !  PMDec [deg/day] JWF B30122
          dDec   = value*(MJD0-JD0)
        end if
c	dvalue = Declist(Jsrc) + dDec           ! JWF B30107 [add offset to ref epoch]
        dvalue = 0.01d0*(dble(Table(Jsrc,40))+dble(Table(Jsrc,41))) + dDec
	!  check that DEC is not out of bounds; if so, wrap about +-90 degrees
        if (dvalue.gt.90.d0) then
                dvalue = 180.d0 - dvalue    !  JWF B21015
        else if (dvalue.lt.-90.d0) then
                dvalue = -(180. + dvalue)
        endif

	call double_stringform (dvalue, format, s0, L)
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)
c                                           ! sigra_pm, sigdec_pm, sigradec_pm
	do ii=35, 37

	  value = Table (Jsrc,ii)

	  inull = 0

	  if ((value.gt.99.999).or.(value.lt.-99.999)) then
		  inull = 1
		  value = 99.999
	  endif

	  if (IzBad(value)) inull = 1

          if (ii .eq. 35) then
            format = 'f9.4'
          else if (ii .eq. 36) then
            format = 'f10.4'
          else
            format = 'f12.4'
          end if

	  call real_stringform (value, format, s0, L)

	  if (inull.eq.1) then
		call nullswap (s0,s1)
		s0 = s1(1:L)
	  endif

          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)

	enddo
c                                           ! PMRA, PMDec, sigPMRA, sigPMDec
c       if (IzBad(Table(Jsrc,29))) then          ! JWF B30502
        if   (IzBad(Table(Jsrc,29)) .or. IzBad(Table(Jsrc,30))
     +   .or. IzBad(Table(Jsrc,31)) .or. IzBad(Table(Jsrc,32))) then
c                       |  PMRA   |  PMDec  | sigPMRA|sigPMDec|
            BigString = BigString(1:lnblnk(BigString))
     +               //'      null      null     null     null'
        else
          value = Table(Jsrc,29)
          if ((abs(value) .ge. 0.01) .and. (abs(value) .le. 999.9999))
     +    then
            write(Str10,'(F10.4)') value
          else
            write(Str10,'(1pE10.2)') value
          end if
          BigString = BigString(1:lnblnk(BigString))//Str10
c
          value = Table(Jsrc,30)
          if ((abs(value) .ge. 0.01) .and. (abs(value) .le. 999.9999))
     +    then
            write(Str10,'(F10.4)') value
          else
            write(Str10,'(1pE10.2)') value
          end if
          BigString = BigString(1:lnblnk(BigString))//Str10
c
          value = Table(Jsrc,31)
          if (value .gt.  999.9999) value =  999.9999
          write(Str9,'(F9.4)') value
          BigString = BigString(1:lnblnk(BigString))//Str9
c
          value = Table(Jsrc,32)
          if (value .gt.  999.9999) value =  999.9999
          write(Str9,'(F9.4)') value
          BigString = BigString(1:lnblnk(BigString))//Str9
        end if
        LBIG = lnblnk(BigString)
c                                                ! w?snr_pm, w?flux_pm, w?mpro_pm,
	do 700 ib=1,numbands                     ! ?sigmpro_pm, w?rchi2_pm

	  mpro_pm(ib)    = 99. ! null value
	  sigmpro_pm(ib) = 0.  ! null value
	  snr_pm(ib)     = 0.  ! null value
          rchi_pm(ib)    = 0.
          kill(ib)       = 0

	  if ((temp(ib).eq.0).or.(kill(ib).eq.1)) goto 700

          signal = Table (Jsrc,41+ib)
	  fnoise = Table (Jsrc,45+ib)
	  rchi_pm(ib) = Table (Jsrc,49+ib)

	  if (fnoise.le.0.)  then ! problem with flux noise;  kill the source
	    kill(ib) = 1
            goto 700
          endif

	  if ((signal.gt.0.).and.(fnoise.gt.0.)) then
	    snr_pm(ib) = signal / fnoise
	    snr_pm(ib) = min (snr_pm(ib),9999.)
	  endif

	  Thresh = 2. * fnoise  ! upper limit, if needed

	  if (signal.lt.Thresh) then  ! compute upper limit
	     if (signal.lt.0.) then
		 upper_lim = Thresh
	     else
		 upper_lim = signal + thresh
	     endif
	     signal = upper_lim   ! set the signal to the upper limit value
	     fnoise = 0.  ! null the flux noise for upper limits
	   endif

	   if (signal.le.0.) goto 700   ! all else fails, give up on source flux

           mpro_pm(ib) = zero(ib) - (2.5*log10(signal))

	   if (fnoise.gt.0.) then
	     sigmpro_pm(ib) = (1.0857 * fnoise)  / signal
	   else  ! upper limit
	     sigmpro_pm(ib) = 9.99  ! indicates upper limit
	   endif

	  if (mpro_pm(ib).lt.-9.99) mpro_pm(ib) = -9.99
	  if (mpro_pm(ib).gt.30.) mpro_pm(ib) = 99.

 700	continue
c
	do ib=1,numbands                         ! w?snr_pm
	 inull = 0
	 if (snr_pm(ib).le.0.) then
	    inull=1
	    snr_pm(ib) = 0.
	 endif

	 value = snr_pm(ib)
          format = 'f9.1'
          call real_stringform (value, format, s0, L)
	  if (IzBad(value)) inull = 1

          if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
          endif
          L1 = LBIG + 1
          L2 = L1 + (L-1)
          LBIG = L2
          BigString (L1:L2) = s0(1:L)
	enddo
c                                                ! w?flux_pm
	format = '1pe12.4'
	do ib=1,numbands
	   call real_stringform (Table(Jsrc,ib+41), format, s0, L)

	   if ((kill(ib).eq.1).or.(temp(ib).eq.0)) then
		call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
	   L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   s0 = '            '
	   L1 = LBIG + 1
	   L2 = L1 + 1
	   LBIG = L2
           BigString (L1:L2) = s0(1:L)


	   fnoise = Table (Jsrc,45+ib)

	   call real_stringform (fnoise, format, s0, L)
	   if (kill(ib).eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
	   L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)
	enddo
c
c                                                ! w?mpro_pm, w?sigmpro_pm, w?rchi2_pm
c
	do ib=1,numbands
	   inull=0
	   value = mpro_pm(ib)
	   if (value.gt.30.) inull=1
	   if (value.lt.-30.) inull=1
	   if (IzBad(value)) inull = 1
	   if (temp(ib).eq.0) inull = 1

	   format = 'f7.3'
	   call real_stringform (value, format, s0, L)
           if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
		rchi(ib) = 0.
           endif
c	   L1 = LBIG + 1   ! JWF B30213
	   L1 = LBIG + 4   ! JWF B30213
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   s0 = '            '
c	   L1 = LBIG + 1   ! JWF B30213
	   L1 = LBIG + 4   ! JWF B30213
           L2 = L1 + 2
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   value = sigmpro_pm(ib)
	   if (value.le.0.) inull=1
	   if (value.gt.5.) inull=1
	   if (IzBad(value)) inull = 1

	   call real_stringform (value,format, s0, L)
           if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
           L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)

	   inull = 0
	   value = rchi_pm(ib)
	   if (value.le.0.) inull = 1
	   if (IzBad(value)) inull = 1
	   format = '1pe11.3'
	   if (inull.eq.1) value = 0.
	   call real_stringform (value,format, s0, L)
           if ((kill(ib).eq.1).or.(inull.eq.1)) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
           endif
           L1 = LBIG + 1
           L2 = L1 + (L-1)
           LBIG = L2
           BigString (L1:L2) = s0(1:L)
	enddo
c
c                                                ! rchi2_pm
c
	inull = 0
	value = Table(Jsrc,54)
	if (value.le.0.) inull=1
	if (value.gt.99999.) inull=1
	if (IzBad(value)) inull = 1
	if (inull.eq.1) value = 0.

	format = '1pe11.3'
	call real_stringform (value,format, s0, L)
        if (inull.eq.1) then
                call nullswap (s0,s1)
                s0 = s1(1:L)
        endif
        L1 = LBIG + 1
        L2 = L1 + (L-1)
        LBIG = L2
        BigString (L1:L2) = s0(1:L)
c
80      Itemp = abs(Table(Jsrc,59))              ! pmcode/bc
        if (IzBad(Table(Jsrc,59))) then                        ! JWF B30505
          print *,'writeTABLE(WARNING): pmcode NaN/Inf for output source',Jsrc
          print *,'nulling out pmcode...'
          Itemp = 0
        end if
        if (Itemp .eq. 0) then
          Str7 = '   null'
        else
          value = 1.0e5*(abs(Table(Jsrc,59)) - float(Itemp))
          if (value .gt. 999.0) value = 999.0
          if (Itemp .gt. 9) Itemp = 9
          if (Table(Jsrc,59) .lt. 0.0) then
            if (SingleFrame) then
              write(Str3,'(I2,''Y'')') Itemp
              BigString = BigString(1:lnblnk(BigString))//Str3
              go to 90
            else
              write(Str7,'(I3,''Y'',I3.3)') Itemp, NInt(value)
              BigString = BigString(1:lnblnk(BigString))//Str7
            end if
          else
            if (SingleFrame) then
              write(Str3,'(I2,''N'')') Itemp
              BigString = BigString(1:lnblnk(BigString))//Str3
              go to 90
            else
              write(Str7,'(I3,''N'',I3.3)') Itemp, NInt(value)
              BigString = BigString(1:lnblnk(BigString))//Str7
            end if
          end if
        end if
c
        Itemp = NInt(Table(Jsrc,57))             ! nIters_pm
        if (Itemp .gt. 999999) Itemp = 999999
        write(Str10,'(I10)') Itemp
        BigString = BigString(1:lnblnk(BigString))//Str10
        Itemp = NInt(Table(Jsrc,58))             ! nSteps_pm
        if (Itemp .gt. 999999) Itemp = 999999
        write(Str10,'(I10)') Itemp
        BigString = BigString(1:lnblnk(BigString))//Str10
c
c==================== end of code added by JWF B30129 =====================
c
90      LBig = lnblnk(BigString)
c
   	write (iunit,'(a)') BigString (1:LBig)

c
c	if (Jsrc.gt.99) call exit(0)
c
	return
	end
c
c==========================================================================
c
         Function IzBad(R4)

         Logical   IzBad
         Real*4    R4, R4a
         Integer*4 I4, MaskExp, IAnd

         Equivalence (I4, R4a)
    
         Data MaskExp/2139095040/
         R4a = R4
         IzBad = IAnd(I4,MaskExp)  .eq. MaskExp

         Return
         End


