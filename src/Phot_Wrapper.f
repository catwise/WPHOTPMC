	subroutine Phot_Wrapper (fmode,wflag,NSUMmax,nsrc,
     1       nnx,nny,nx,ny,nf,nfpix,
     1       Array,MASK, iMAsk, Unc, pix_order, SNRlist,
     1       Xpos,Ypos,zero, xTRANS, yTRANS,
     1       Rstap, Rstann, Rstwid, BGmax,
     1       LBACK, LSIG, Lconf, MAGSTD, eMAGSTD, FLGSTD, Nsat, Rsat, pscale,
     1       NtotSrc,SPIT, MIPS, fbits, NbannAve, IOCname, doIOC)


	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*200 IOCname, fname
	character*1 s0
        integer*4 pix_order(nf,4)
	real*4 Array (nnx,nny,nfpix), Unc(nnx,nny,nfpix)
        integer*4 iMask (nnx,nny,nfpix)
	real*4 Xpos (NSUMmax,nf,4),Ypos (NSUMmax,nf,4), zero(4)
	real*4 Rstap(4), Rstann(4), Rstwid(4),Raper(19,4)
	real*4 SNRlist(NSUMmax,nf,4)
	real*4 LBACK(NSUMmax,nf,4), LSIG(NSUMmax,nf,4), Lconf(NSUMmax,nf,4)
        real*4 eMAGSTD(NSUMmax,4), MAGSTD(NSUMmax,4)
	real*4 fluxtmp(19), efluxtmp(19), zmag(19), zerr(19)
	real*4 BGmax(4), Rsat (Nsat,4), pscale(4), pct
	integer nx(4),ny(4), flag, nsrc, FLGSTD (NSUMmax,4), ngood_band(4), ngood, ntot(19)
	integer iflag(19),fbits(32)
c       integer MASK(nnx,nny,nf,4)
        integer MASK(1,1,1) ! Mask being dummied out since it's never assigned anything -- TPC
	integer NbannAve(4), nave(4)
	integer*2 wflag(nf,4), fmode
	integer*2 xTRANS (nf,4), yTRANS (nf,4)
	logical SPIT, MIPS, doIOC, doConf

	real*4, allocatable :: flux(:,:), sigflux(:,:)
	real*4 fluxmed, q1, q2, sigfluxmed, wtdflux(4), wtdsigflux(4), fluxwt(4)
	integer*4 stat,nflux(4),used(4)

	write (6,*) 'executing PHOT_WRAPPER ',nsrc

	doConf = .true.  !  compute and pass the confusion noise
        pct = 0.1 ! Percentile trim to use for trimmed flux mean
	ngood = 0
	ngood_band(1:4) = 0

	if (doIOC) then  ! should be single frame, single band

		L = numchar (IOCname)

		 do ib = 1,4


c                  if (wflag(1,ib) .gt. 0) then
		  if ( any ( wflag ( 1:nfi,ib ) ==1 ) ) then
		

			write (s0,'(i1)') ib

			fname = IOCname(1:L) // '-w' // s0 // '-mdexAphot.tbl'
			M = numchar(fname)
	
			write (6,'(a)') fname(1:M)

			open (unit=81,file=fname)
			write (81,'(a,i3)') '\band = ',ib
			write (81,'(a,f8.2,a)') '\Rap = ',Rstap(ib),' pix'
			write (81,'(a)') '|    x    |     y   |     sky |   skysig|  amag |amagsig| flg|'
			write (81,'(a)') '|   pix   |    pix  |      dn |      dn |  mag  |   mag | -- |'
		  endif

	  	 enddo


	endif

        ! allocate temporary storage for er-frame fluxes for later median'ing
	allocate(flux(nf,4), sigflux(nf,4),   stat=stat)
	if(stat .ne. 0) then
	  print *,'ERROR (phot_wrapper): cannot allocate flux/sigflux array; Stat =',stat
	  call exit(64)
	endif

	do ib=1,4
	 NbannAve(ib) = 0
	 nave (ib) = 0
	enddo

	NtotSrc = 0
	nap = 1  ! small aperture
	

c	write (6,*) NsrcAll,nf,wflag(1,1),wflag(2,1),wflag(3,1)

        LBACK(:,:,:) = -999.99
        LSIG(:,:,:) = 0.
        SNRlist(:,:,:) = 0.
        Lconf(:,:,:) = 0.


	do 500 Jsrc = 1, nsrc

c	   if (Jsrc.lt.200) goto 500
	 MAGSTD (Jsrc,1:4) = 99.
	 eMAGSTD (Jsrc,1:4) = 0.0
         FLGSTD (Jsrc,1:4) = 0

c	write (6,*) 'DEBUG  aperures',Rstap(1:nap)

c        Temps for computing median flux, and from that, e/magstd
         flux(:,:) = 0
         sigflux(:,:) = 0
	 nflux(:) = 0
	 fluxwt(:) = 0
	 wtdflux(:) = 0
	 wtdsigflux(:) = 0

	 do Jfr = 1, nf

	 do 300 ib = 1,4

	  if (wflag(Jfr,ib) .eq. 0) goto  300
	  
	  do Map=1,nap
		Raper(Map,ib) = Rstap(ib) 
	  enddo

c	write (6,*) 'DEBUG-x',ib

c	write (6,*) 'DEBUG-y',ib

c	write (6,*) 'DEBUG ',Jsrc,Jfr,ib
	if (XPos (Jsrc,Jfr,ib).le.0.) goto 300
        if (YPos (Jsrc,Jfr,ib).le.0.) goto 300

c	write (6,*) 'here we go ',Jfr,ib
c	write (6,*)  XPos (Jsrc,Jfr,ib),xTRANS(Jfr, ib)
c	write (6,*)  YPos (Jsrc,Jfr,ib),yTRANS(Jfr, ib)

	  xx = XPos (Jsrc,Jfr,ib) - (xTRANS(Jfr, ib)*1.)
          yy = YPos (Jsrc,Jfr,ib) - (yTRANS(Jfr, ib)*1.)

	  Rcrater = Rsat (Jsrc,ib) / pscale (ib)

c	write (6,*) 'Rcrater = ',ib,Rsat (Jsrc,ib),Rcrater

c	! TEMP
c	if ( (ib.eq.1).and.(Jsrc.eq.565) ) 
c     1      write (64,'(2I4,4f9.1)') Jsrc,Jfr,XPos (Jsrc,Jfr,ib), YPos (Jsrc,Jfr,ib), xx, yy

c	if ( (ib.eq.1).and.(Jsrc.eq.576) )
c     1      write (64,'(2I4,4f9.1)') Jsrc,Jfr,XPos (Jsrc,Jfr,ib), YPos (Jsrc,Jfr,ib), xx, yy


          if (xx.le.0) goto 300
          if (yy.le.0) goto 300
	  if (xx.gt.nnx) goto 300
	  if (yy.gt.nny) goto 300

	  NtotSrc = NtotSrc + 1

	    mode = 0
	    if (fmode.eq.1) mode = 1 !  use the MASK

c	if ((Jsrc.eq.94).and.(ib.eq.1).and.(Jfr.eq.18)) write (6,*) Jsrc,Jfr,ib,xx,yy

c	write (6,*) 'call AnnulusStats'
	    Nbann = 0
	    call AnnulusStats (Jsrc,pix_order(Jfr,ib),nf,nfpix,ib,nnx,nny,nx(ib),ny(ib),
     1            Array,MASK,iMask, Unc,
     1            xx,yy, Rstann, Rstwid, nann, sky, sdev, Nbann, Rl, Rh, mode, SPIT, fbits, zero,
     1            SKYrms, Hsig, RMSformal, BGmax(ib), Rcrater)
	    if (Nbann.gt.0) then
		nave(ib) = nave(ib)+1
		NbannAve(ib) = NbannAve(ib) + Nbann
	    endif
c
c	write (6,*) 'DEBUG sky,sdev ',sky,sdev


	!  RMSformal is the skysig derived from AWAIC UNC images
	
	    T2 = (sdev**2) - (RMSformal**2)

	    if ((T2.ge.0.).and.(doConf)) then
	      confusion_noise = sqrt (T2)
	    else
	      confusion_noise = 0.
	    endif

	    LBACK (Jsrc,Jfr,ib) = sky
            LSIG (Jsrc,Jfr,ib) = sdev
	    Lconf (Jsrc,Jfr,ib) = confusion_noise

	    if (sky.gt.-500.) then

	     ibox = 3
	     call find_peak (pix_order(Jfr,ib),ib,nf,nfpix,nnx,nny,nx(ib),ny(ib),
     1          Array,ibox,xx,yy,ip,jp,valmax,ibog)
	     SNRlist (Jsrc,Jfr,ib) = (valmax-sky) / sdev


c	write (6,*) 'DEBUG circphot '
              call CircPhot (Jsrc,pix_order(Jfr,ib),nf,nfpix,ib,nnx,nny,nx,ny,Array,MASK, iMask, Unc,
     1            xx,yy, nap, Raper , nann, sky, sdev, Nbann, zero,
     1            fluxtmp, efluxtmp, zmag, zerr, iflag, ntot, SPIT, MIPS, fbits )


	      Map = 1

c	 if (ib.eq.3) write (6,*) xx,yy,zero(ib),fluxtmp(1), efluxtmp(1), zmag(1), zerr(1)
c	if (jsrc.gt.20) return

c	write (6,*) 'DEBUG ',xx,yy,zero(ib),fluxtmp(1), efluxtmp(1), zmag(1), zerr(1)

c	      if ((zmag(MAP).gt.0.).and.(zmag(MAP).lt.25.)) then
              if(efluxtmp(map) .gt. 0 .and. fluxtmp(map) .gt. -500 .and. sky .gt. -500) then

		nflux(ib) = nflux(ib) + 1

		flux(nflux(ib),ib) = fluxtmp(map)
		sigflux(nflux(ib),ib) = efluxtmp(map)

		fluxwt(ib) = fluxwt(ib) + ntot(map)
		wtdflux(ib) = wtdflux(ib) + fluxtmp(map)*ntot(map)
		wtdsigflux(ib) = wtdsigflux(ib) + (fluxtmp(map)**2)*ntot(map)

c		MAGSTD (Jsrc,ib) = zmag(Map)
c		eMAGSTD (Jsrc,ib) = zerr(Map)

	      else  ! negative flux;  use upperlimit

c		if (MAGSTD (Jsrc,ib).gt.25.) then
c
c		  thresh = 2. * sdev
c		  MAGSTD (Jsrc,ib) = zero(ib) - (2.5*log10(thresh))
c		  eMAGSTD (Jsrc,ib) = 0.50 ! 2*sigma
c		  FLGSTD (Jsrc,IB) =  32 ! upper limit
c
c		endif

	     endif

	    else
c		FLGSTD (Jsrc,IB) =  8  !  sky corruption
            endif


c	MAGSTD (Jsrc,ib) = 25.  ! TEMP
c	eMAGSTD (Jsrc,ib) = 0.5

c	   write (71,'(i3,2f10.2,i6,2f8.2,2f7.3,i4)') ib,xx,yy,nann, sky, sdev, 
c     1       MAGSTD (K,J,ib), eMAGSTD (K,J,ib), FLGSTD (K,J,IB)


	    if (doIOC) then  ! should be single frame, single band
		write (81,'(2f10.2,2f10.2,2f8.3,i5)')  xx,yy,sky, sdev,MAGSTD (Jsrc,ib), eMAGSTD (Jsrc,ib), 
     1           FLGSTD (Jsrc,IB)
	    endif

 300	    continue  ! ib

	  enddo  ! frames

c	write (6,*) 'DEBUG ',Jsrc,xx,yy,nflux(1:2)
c	 write (6,*) 'DEBUG ',flux(1:nflux(1),1)
c	 write (6,*) 'DEBUG ',sigflux(1:nflux(1),1)

c         compute mags from median fluxes or trimmed average fluxes
	  igotit = 0
	  do ib=1,4
	    if(nflux(ib) .gt. 0) then
	      call trimsort(flux(1:nflux(ib),ib), nflux(ib), fluxmed, q1, q2, pct, fluxtrim, sigfluxtrim)
	      if(nflux(ib) .gt. 7) then
                ! Use +/-1 sig points in the flux distribution
		sigfluxmed = (q2-q1)/2/sqrt(1.0*nflux(ib))
              else
		! For low coverage, use the median sigflux
	        call medsort(sigflux(1:nflux(ib),ib), nflux(ib), sigfluxmed, q1, q2)
		sigfluxmed = sigfluxmed/sqrt(1.0*nflux(ib)) ! Very hand-wavey
              endif

c	write (6,*) 'DEBUG tricky ',ib,nflux(ib),fluxmed, fluxtrim, sigfluxtrim, sigfluxmed

	      if(fluxwt(ib) .gt. 50 .and. wtdflux(ib) .gt. 0) then
		wtdflux(ib) = wtdflux(ib)/fluxwt(ib)
		wtdsigflux(ib) = sqrt(max(0.0,(wtdsigflux(ib)/fluxwt(ib) - wtdflux(ib)**2)/fluxwt(ib)))
	      else 
		wtdflux(ib) = 99
                wtdsigflux(ib) = 0
	      end if
	      used(ib) = 0
	      fluxest = -99
	      sigest = 0
	      !if(wtdflux(ib) .gt. 0 .and. wtdsigflux(ib) .gt. 0) then
	      !  fluxest = wtdflux(ib)
	      !  sigest = wtdsigflux(ib)
	      !  used(ib) = 1
	      if(fluxtrim .gt. 0 .and. sigfluxtrim .gt. 0) then
                 ! Prefer the trimmed ave. if available
                fluxest = fluxtrim
                sigest = sigfluxtrim
		used(ib) = 2
              else  ! For single-frame use in frame pipeline or for unwise
                fluxest = fluxmed
                sigest = sigfluxmed
                used(ib) = 3
              endif
	      !!!print *, 'ib,jsrc,n,fluxtrim,sigfluxtrim,fluxmed,sigfluxmed = ',ib,jsrc,nflux(ib),fluxtrim,sigfluxtrim,fluxmed,sigfluxmed
	      if(fluxest .gt. 0 .and. sigest .gt. 0 .and. sigest .le. 0.5*fluxest) then
		! Compute mag and emag
	        MAGSTD(Jsrc,ib) = zero(ib) - (2.5*log10(fluxest))
                eMAGSTD(Jsrc,ib) = 1.0857 * sigest/fluxest
                ngood_band(ib) = ngood_band(ib) + 1
		if(igotit .eq. 0) ngood = ngood + 1
		igotit = 1
              else
		! Illegal value(s)
                MAGSTD(Jsrc,ib) = 99
                eMAGSTD(Jsrc,ib) = 0
                FLGSTD(Jsrc,ib) = 8 ! no value
	      endif
            else
	      ! Insufficient good measurements
              MAGSTD(Jsrc,ib) = 99
              eMAGSTD(Jsrc,ib) = 0
              FLGSTD(Jsrc,ib) = 8 ! no value
	    end if
          end do

	  if(mod(jsrc-1,250) .eq. 0) then
            print *,'--- #',Jsrc
            print *, 'MagSTD:  ',MAGSTD(Jsrc,:) 
            print *, 'eMagSTD: ',eMAGSTD(Jsrc,:)
            print *, 'WT:      ',fluxwt
            print *, 'N:       ',nflux
            print *, 'Used:    ',used
          endif

 500	continue ! jsrc

	print *,'Recovered mags for ',ngood,' sources.'
	print *,'Count per band = ', ngood_band

	! average Nbann
	
	do ib=1,4
	  if (nave(ib).gt.0) then
		NbannAve(ib) = nint(1.*NbannAve(ib) / (nave(ib)*1.))
	  endif
	enddo

	deallocate(flux, sigflux)

	return
	end
