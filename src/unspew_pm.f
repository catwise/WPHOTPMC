c	subroutine unspew(p,psfs,npsize,npnmax,varpsf,nlayers,isn,onframe,
c    *	    nmax,nf,mxbands,wflag,nfi,framesig,sigx,sigy,sigxy,sigflux,
c    *	    sigpmr,sigpmd,ok)
	subroutine unspew_pm(p,psfs,npsize,npnmax,varpsf,nlayers,isn,onframe,
     *	    nmax,nf,mxbands,wflag,nfi,framesig,sigx,sigy,sigxy,sigflux,
     *	    sigpmr,sigpmd,covRpmR,covDpmD,covpmRD,ok,debug)          ! JWF B30412
 
! UNcertainties of Simultaneous Parameter Estimates for Wise. The output
! array, gamma, contains the a posteriori covariance of the parameter estimates.
 
	implicit real*4 (a-h,o-z)
	implicit integer (i-n)
	include 'wpro.inc'
	integer, parameter :: maxblend = 100
	integer, parameter :: maxpsize =  256
	integer, parameter :: maxbands = 4
	integer, parameter :: maxsegs = 121
	real, parameter :: dtor = 0.0174533
 
	real(4), allocatable :: gamma(:,:)
	real(4), allocatable :: xi(:), eta(:)
	real(4), allocatable :: xsource(:,:), ysource(:,:)
	real(4), allocatable :: xsourcep(:,:), ysourcep(:,:)
	real(4), allocatable :: xsourceq(:,:), ysourceq(:,:)
	real(4), allocatable :: G(:,:)
	real(4), allocatable :: sumxx(:),sumyy(:),sumxy(:)
	real(4), allocatable :: sumxxt(:),sumyyt(:),sumxyt(:)
	real(4), allocatable :: sumxxtt(:),sumyytt(:),sumxytt(:)
        real*4   sumHxx, sumHyy, sumHxy, sumHxxt, sumHyyt,
     +           sumHxyt, sumHxxtt, sumHyytt, sumHxytt
 
	real(8) ra8, dec8, x8, y8
	real(8) racan(maxblend),deccan(maxblend)
	real(4) p(*)
	real(4) psfs(npsize,npsize,npnmax,4)
	real(4) sigx(maxblend),sigy(maxblend),sigxy(maxblend)
	real(4) dvalue(maxpix,maxbands), sigdvalue(maxpix,maxbands)
	real(4) msigdvalue(maxpix,maxbands)
	real(4) varpsf(maxpix,maxbands),rchisqb(maxbands)
	real(4) pxc(maxbands),pyc(maxbands),psamp(maxbands)
	real(4) pint(maxpix,maxbands,maxblend)
	real(4) pintx(maxpix,maxbands,maxblend)
	real(4) pinty(maxpix,maxbands,maxblend)
	real(4) xpixscale(maxframes),ypixscale(maxframes),tframe(maxframes)
	real(4) flux(maxbands,maxblend),sigflux(maxbands,maxblend)
	real(4) framesig(maxframes,maxbands,maxblend)
	integer ivalue(maxpix,maxbands), jvalue(maxpix,maxbands)
	integer ifvalue(maxpix,maxbands),primask(maxpix,maxbands)
	integer ifset(maxpix,maxbands)
	integer ib(maxbands),nvalues(maxbands)
	integer WCS(maxframesin,maxbands), offscl, nnpsf(maxbands)
	integer ifq(maxframes),ibq(maxframes),kpn(maxframes)
c	integer(2) nlayers(maxframes),wflag(nfi,4)         ! JWF B60714
	integer(2) nlayers(maxbands),wflag(nfi,4)          ! JWF B60714
	logical(4) ok,goodcomp(maxblend),goodcompz(maxblend)
	logical(1) onframe(nmax,nf,mxbands)
 
c       integer*4 ndump    !  JWF B21010
        logical   IzBad    !  JWF B21010
        logical*4 debug    !  JWF B30412
        real*4     dra, ddec(maxbands)                ! JWF B30513
        real*4    covRpmR, covDpmD, covpmRD           ! JWF B21221
c       data      ndump/0/ !  JWF B21010
 
	common /funccom/ racan,deccan,nvalues,ivalue,jvalue,dvalue,
     *	    sigdvalue,ifvalue,WCS,ifq,ibq,primask,nframes,kpn,
     *	    nnpsf,pxc,pyc,psamp,nblend,goodcomp,icp,
     *	    xpixscale,ypixscale,flux,rchisq,rchisqb,tframe
c
        include 'jwfcom.f'             ! JWF B21219
c
	allocate(xi(nblend))
	allocate(eta(nblend))
	allocate(xsource(nblend,nframes))
	allocate(ysource(nblend,nframes))
	allocate(xsourcep(nblend,nframes))
	allocate(ysourcep(nblend,nframes))
	allocate(xsourceq(nblend,nframes))
	allocate(ysourceq(nblend,nframes))
 
	mblend = ntrue(goodcomp,nblend)
 	nc = 0
	do n = 1,nblend
	    if (goodcomp(n)) then
		nc = nc + 1
		xi(n) = p(nc)			! offsets east and north from
		eta(n) = p(mblend+nc)		! candidate location [deg]
	    endif
	enddo
	pmr = p(2*mblend+1)/365.25	! deg/day
	pmd = p(2*mblend+2)/365.25	! deg/day
 
! Which bands are present?
	nbands = 0	
	do i = 1,maxbands
	    if(nvalues(i) > nblend+2) then
		nbands = nbands + 1
		ib(nbands) = i
	    endif
	enddo
 
! Make sure all components have nonzero flux in at least one band.
	goodcompz = goodcomp
	do n = 1,nblend
	  if (goodcomp(n)) then
	    if (.not.any(flux(ib(1:nbands),n) /= 0.)) goodcompz(n) = .false.
	  endif
	enddo
 
! Allocate space for the Fisher information matrix.
	mblend = ntrue(goodcompz,nblend)
	ng0 = (nbands + 2)*mblend
	ng = ng0 + 2
	allocate (gamma(ng,ng))

	step = 0.1		! step size for gradient calculation [pixels]

! Find pixel locations of source at offset (xi,eta).
	do iframe = 1,nframes
	    iwcs = WCS(ifq(iframe),ibq(iframe))
	    ddec(ibq(iframe)) = step*xpixscale(iframe)
	    dra = ddec(ibq(iframe))/cos(dtor*deccan(1))

	    do n = 1,nblend
        	ra8 = racan(n) + xi(n)/cos(dtor*deccan(n))
        	dec8 = deccan(n) + eta(n)
		if (n==icp) then
		    ra8 = ra8 + pmr*tframe(iframe)/cos(dtor*deccan(n))
		    dec8 = dec8 + pmd*tframe(iframe)
		endif
        if (debug) then
          print *,'unspew_pm(128): doing blend component no.',n
          print *,'ra8,dec8:', ra8, dec8
          print *,'pmr, pmd, iframe, tframe(iframe):',pmr, pmd, iframe, tframe(iframe)
        end if
        offscl = -1
        call wcs2pix(iwcs, ra8, dec8, x8, y8, offscl)
        if (offscl .ne. 0) then
		  xsource(n,iframe) = 1.
		  ysource(n,iframe) = 1.
        else
		  xsource(n,iframe) = x8
		  ysource(n,iframe) = y8
		endif
        if (debug) print *,'iframe,x8,y8:',iframe,x8,y8 ! JWF dbg
        offscl = -1
        call wcs2pix(iwcs, ra8+dra, dec8, x8, y8, offscl)
        if (offscl .ne. 0) then
		  xsourcep(n,iframe) = 1.
		  ysourcep(n,iframe) = 1.
        else
		  xsourcep(n,iframe) = x8
		  ysourcep(n,iframe) = y8
		endif
        offscl = -1
        call wcs2pix(iwcs, ra8, dec8+ddec(ibq(iframe)), x8, y8, offscl)
        if (offscl .ne. 0) then
		  xsourceq(n,iframe) = 1.
		  ysourceq(n,iframe) = 1.
        else
		  xsourceq(n,iframe) = x8
		  ysourceq(n,iframe) = y8
		endif
	    enddo
	enddo

! Interpolate the PSF and its gradient.
	msigdvalue = sigdvalue
	do i = 1,nbands
	    npsf = nnpsf(ib(i))
 
	    do k = 1,nvalues(ib(i))
		if (nlayers(ib(i)) > 1) msigdvalue(k,ib(i)) = sqrt(
     *		    sigdvalue(k,ib(i))**2 + (nlayers(ib(i))-1) *
     *		    varpsf(k,ib(i)))
		iframe = ifvalue(k,ib(i))
		ifset(k,i) = ifq(iframe)
		pmin = pminFac*psfs(nint(pxc(ib(i))),nint(pyc(ib(i))), ! JWF B30507
     *		    kpn(iframe),ib(i))
	        xprat = xpixscale(iframe)/psamp(ib(i))
	        yprat = ypixscale(iframe)/psamp(ib(i))
	        rn = xprat*yprat
		ik = ivalue(k,ib(i))
		jk = jvalue(k,ib(i))
	    	do n = 1,nblend
		    xx = pxc(ib(i)) + (ik-xsource(n,iframe))*xprat
		    yy = pyc(ib(i)) + (jk-ysource(n,iframe))*yprat
		    pint(k,i,n) =rn*max(bilint(psfs,npsize,npsize,
     *			npnmax,4,npsf,npsf,kpn(iframe),ib(i),xx,yy),
     *			pmin)
		    xx = pxc(ib(i)) + (ik-xsourcep(n,iframe))*xprat
		    yy = pyc(ib(i)) + (jk-ysourcep(n,iframe))*yprat
		    pintx(k,i,n) = (rn*max(bilint(psfs,npsize,npsize,
     *			npnmax,4,npsf,npsf,kpn(iframe),ib(i),xx,yy),
     *			pmin) - pint(k,i,n))/ddec(ibq(iframe))
		    xx = pxc(ib(i)) + (ik-xsourceq(n,iframe))*xprat
		    yy = pyc(ib(i)) + (jk-ysourceq(n,iframe))*yprat
		    pinty(k,i,n) = (rn*max(bilint(psfs,npsize,npsize,
     *			npnmax,4,npsf,npsf,kpn(iframe),ib(i),xx,yy),
     *			pmin) - pint(k,i,n))/ddec(ibq(iframe))
		enddo
	    enddo
	enddo
 
	allocate (G(ng,ng))
        if (debug) print *,'unspew_pm(201): G allocated with ng =', ng !  JWF dbg
	allocate (sumxx(nbands))
	allocate (sumyy(nbands))
	allocate (sumxy(nbands))
	allocate (sumxxt(nbands))
	allocate (sumyyt(nbands))
	allocate (sumxyt(nbands))
	allocate (sumxxtt(nbands))
	allocate (sumyytt(nbands))
	allocate (sumxytt(nbands))
	G = 0.
	gamma = 0.
 
! Start by calculating the flux sigmas for individual frames for the benefit
! of N out of M statistics.
c                                 But only upon request!   JWF B30518

        if (debug) then                                !   JWF B30518
	do ifs = 1,nf		! begin frameset loop
	 framesig(ifs,1:4,1:nblend) = 0.
	 gamma = 0.
	 if (any(onframe(isn,ifs,1:4) .and. wflag(ifs,1:4)==1)) then
	  nc = 0
	  do n = 1,nblend
	  if (goodcompz(n)) then
	   nc = nc+1
	   mc = 0
	   do m = 1,nblend
	   if (goodcompz(m)) then
	    mc = mc+1
	    do i = 1,nbands
		sumHxx = 0.
		sumHyy = 0.
		sumHxy = 0.
		sumHxxt = 0.
		sumHyyt = 0.
		sumHxyt = 0.
		sumHxxtt = 0.
		sumHyytt = 0.
		sumHxytt = 0.
		do k = 1,nvalues(ib(i))
		  if(ifset(k,i)==ifs) then
		    sumHxx = sumHxx + pintx(k,i,m)*pintx(k,i,n)/
     *			sigdvalue(k,ib(i))**2       ! BUGFIX 3/18/13
		    sumHyy = sumHyy + pinty(k,i,m)*pinty(k,i,n)/
     *			sigdvalue(k,ib(i))**2       ! BUGFIX 3/18/13
		    sumHxy = sumHxy + pintx(k,i,m)*pinty(k,i,n)/
     *			sigdvalue(k,ib(i))**2       ! BUGFIX 3/18/13
		    if (m==icp .and. n==icp) then
			tk = tframe(ifvalue(k,ib(i)))
			sumHxxt = sumHxxt + tk*pintx(k,i,m)*pintx(k,i,n)/
     *			    sigdvalue(k,ib(i))**2   ! BUGFIX 3/18/13
			sumHyyt = sumHyyt + tk*pinty(k,i,m)*pinty(k,i,n)/
     *			    sigdvalue(k,ib(i))**2   ! BUGFIX 3/18/13
			sumHxyt = sumHxyt + tk*pintx(k,i,m)*pinty(k,i,n)/
     *			    sigdvalue(k,ib(i))**2   ! BUGFIX 3/18/13
			sumHxxtt = sumHxxtt + tk**2 * pintx(k,i,m)*pintx(k,i,n)/
     *			    sigdvalue(k,ib(i))**2   ! BUGFIX 3/18/13
			sumHyytt = sumHyytt + tk**2 * pinty(k,i,m)*pinty(k,i,n)/
     *			    sigdvalue(k,ib(i))**2   ! BUGFIX 3/18/13
			sumHxytt = sumHxytt + tk**2 * pintx(k,i,m)*pinty(k,i,n)/
     *			    sigdvalue(k,ib(i))**2   ! BUGFIX 3/18/13
		    endif
		  endif
		enddo
		sumxx(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxx
		sumyy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyy
		sumxy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxy
		if (m==icp .and. n==icp) then
		    sumxxt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxxt
		    sumyyt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyyt
		    sumxyt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxyt
		    sumxxtt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxxtt
		    sumyytt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyytt
		    sumxytt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxytt
		endif
	    enddo
c           if (debug) print *,'unspew_pm(275): mc, nc =', mc, nc  ! JWF dbg
	    G(mc,nc) = sum(sumxx)
	    G(nc,mc) = G(mc,nc)
            if (debug) print *,'unspew_pm(278): mblend+mc, mblend+nc =', mblend+mc,mblend+nc ! JWF dbg
	    G(mblend+mc,mblend+nc) = sum(sumyy)
	    G(mblend+nc,mblend+mc) = G(mblend+mc,mblend+nc)
	    G(mc,mblend+nc) = sum(sumxy)
	    G(mblend+nc,mc) = G(mc,mblend+nc)
	    do i = 1,nbands
		jj = 2*mblend + (mc-1)*nbands + i
		kk = 2*mblend + (nc-1)*nbands + i
		sumjk = 0.
		sumjn = 0.
		sumjnn = 0.
		do k = 1,nvalues(ib(i))
		  if(ifset(k,i)==ifs) then
		    sumjk = sumjk +
     *			pint(k,i,m)*pint(k,i,n)/sigdvalue(k,ib(i))**2 !BUGFIX 3/18/13
		    sumjn = sumjn +
     *			pintx(k,i,n)*pint(k,i,m)/sigdvalue(k,ib(i))**2 !BUGFIX 3/18/13
		    sumjnn = sumjnn +
     *			pinty(k,i,n)*pint(k,i,m)/sigdvalue(k,ib(i))**2 !BUGFIX 3/18/13
		  endif
		enddo
c               print *,'unspew(299): jj, kk =', jj, kk  ! JWF dbg
		G(jj,kk) = sumjk
		G(kk,jj) = G(jj,kk)
		G(jj,nc) = -flux(ib(i),n)*sumjn
		G(nc,jj) = G(jj,nc)
		G(jj,mblend+nc) = -flux(ib(i),n)*sumjnn
		G(mblend+nc,jj) = G(jj,mblend+nc)
	    enddo  ! i = 1, nbands
	   endif  !  (goodcompz(m))
	   enddo  ! m = 1,nblend
	  endif  ! (goodcompz(n))
	  enddo  ! n = 1,nblend

c         print *,'unspew(311): ng0, ng0+1, ng0+2 =', ng0, ng0+1, ng0+2  ! JWF dbg

	  G(ng0+1, ng0+1) = sum(sumxxtt)
	  G(ng0+2, ng0+2) = sum(sumyytt)
	  G(ng0+1, ng0+2) = sum(sumxytt)
	  G(ng0+2, ng0+1) = G(ng0+1, ng0+2)
	  G(1, ng0+1) = sum(sumxxt)
	  G(ng0+1, 1) = G(1, ng0+1)
	  G(1, ng0+2) = sum(sumxyt)
	  G(ng0+2, 1) = G(1, ng0+2)
c         print *,'unspew(320): mblend+1, ng0+1 =', mblend+1, ng0+1  ! JWF dbg
	  G(mblend+1, ng0+1) = G(1, ng0+2)
	  G(ng0+1, mblend+1) = G(mblend+1, ng0+1)
	  G(mblend+1, ng0+2) = sum(sumyyt)
c	  G(ng0, mblend+1) = G(mblend+1, ng0+2)  ! commented out by JWF B21130
	  G(ng0+2, mblend+1) = G(mblend+1, ng0+2)  ! JWF B21130 [guessing this is what was intended]
	  call matinv(G,gamma,ng,ok)
	  do n = 1,ng
	    if((.not.ok .or. gamma(n,n) <= 0. .or. isnan(gamma(n,n)))
     *		.and. G(n,n) /= 0.) gamma(n,n) = 1./G(n,n)
	  enddo
	  nc = 0
	  do nn = 1,nblend
	   if (goodcompz(nn)) then
	    nc = nc + 1
	    do i = 1,nbands
		n = 2*mblend + (nc-1)*nbands + i
		framesig(ifs,ib(i),nn) = sqrt(max(gamma(n,n),0.))
	    enddo
            write (6,'(a,i2,i4,1p4e11.3)') 'unspew_pm(345):',
     +             nn, ifs, (framesig(ifs,kk,nn), kk = 1, 4)
	   endif
           write (6,'(a,i2,i4,1p4e11.3)') 'unspew_pm(349):', nn, ifs,
     +            3600.*sqrt(max(gamma(nc,nc),0.)),                 ! sigra_pm
     +            3600.*sqrt(max(gamma(mblend+nc,mblend+nc),0.)),   ! sigdec_pm
     +            3600.*365.25*sqrt(max(gamma(ng0+1,ng0+1),0.0)),   ! sigPMRA
     +            3600.*365.25*sqrt(max(gamma(ng0+2,ng0+2),0.0))    ! sigPMDec
	  enddo
	 endif
	enddo		! end frameset loop
        end if          !  debug

! Now calculate matrix elements for full multiframe solution.
	nc = 0
	do n = 1,nblend
	 if (goodcompz(n)) then
	  nc = nc+1
	  mc = 0
	  do m = 1,nblend
	   if (goodcompz(m)) then
	    mc = mc+1
 
! Positional variables.
	    do i = 1,nbands
		sumHxx = 0.
		sumHyy = 0.
		sumHxy = 0.
		sumHxxt = 0.
		sumHyyt = 0.
		sumHxyt = 0.
		sumHxxtt = 0.
		sumHyytt = 0.
		sumHxytt = 0.
		do k = 1,nvalues(ib(i))
                    sumHxx = sumHxx + pintx(k,i,m)*pintx(k,i,n)/
     *			msigdvalue(k,ib(i))**2
		    sumHyy = sumHyy + pinty(k,i,m)*pinty(k,i,n)/
     *			msigdvalue(k,ib(i))**2
		    sumHxy = sumHxy + pintx(k,i,m)*pinty(k,i,n)/
     *			msigdvalue(k,ib(i))**2
		    if (m==icp .and. n==icp) then
                        tk = tframe(ifvalue(k,ib(i)))
                        if (CorrPMerr) then
                          sumHxxt  = sumHxxt + tk*pintx(k,i,m)*pintx(k,i,n)/
     *		       	             msigdvalue(k,ib(i))**2
			  sumHyyt  = sumHyyt + tk*pinty(k,i,m)*pinty(k,i,n)/
     *			             msigdvalue(k,ib(i))**2
			  sumHxyt  = sumHxyt + tk*pintx(k,i,m)*pinty(k,i,n)/
     *			             msigdvalue(k,ib(i))**2
			  sumHxxtt = sumHxxtt + tk**2 * pintx(k,i,m)*pintx(k,i,n)/
     *			             msigdvalue(k,ib(i))**2
			  sumHyytt = sumHyytt + tk**2 * pinty(k,i,m)*pinty(k,i,n)/
     *			             msigdvalue(k,ib(i))**2
			  sumHxytt = sumHxytt + tk**2 * pintx(k,i,m)*pinty(k,i,n)/
     *			             msigdvalue(k,ib(i))**2
                        else
                          sumHxxt  = sumHxxt + tk*pintx(k,i,m)*pintx(k,i,n)/
     *		       	             sigdvalue(k,ib(i))**2
			  sumHyyt  = sumHyyt + tk*pinty(k,i,m)*pinty(k,i,n)/
     *			             sigdvalue(k,ib(i))**2
			  sumHxyt  = sumHxyt + tk*pintx(k,i,m)*pinty(k,i,n)/
     *			             sigdvalue(k,ib(i))**2
			  sumHxxtt = sumHxxtt + tk**2 * pintx(k,i,m)*pintx(k,i,n)/
     *			             sigdvalue(k,ib(i))**2
			  sumHyytt = sumHyytt + tk**2 * pinty(k,i,m)*pinty(k,i,n)/
     *			             sigdvalue(k,ib(i))**2
			  sumHxytt = sumHxytt + tk**2 * pintx(k,i,m)*pinty(k,i,n)/
     *			             sigdvalue(k,ib(i))**2
                        end if
		    endif
		enddo
		sumxx(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxx
		sumyy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyy
		sumxy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxy
		if (m==icp .and. n==icp) then
		    sumxxt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxxt
		    sumyyt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyyt
		    sumxyt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxyt
		    sumxxtt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxxtt
		    sumyytt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyytt
		    sumxytt(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxytt
		endif
	    enddo
	    G(mc,nc) = sum(sumxx)
	    G(nc,mc) = G(mc,nc)
	    G(mblend+mc,mblend+nc) = sum(sumyy)
	    G(mblend+nc,mblend+mc) = G(mblend+mc,mblend+nc)
	    G(mc,mblend+nc) = sum(sumxy)
	    G(mblend+nc,mc) = G(mc,mblend+nc)
 
! Terms involving fluxes.
	    do i = 1,nbands
		jj = 2*mblend + (mc-1)*nbands + i
		kk = 2*mblend + (nc-1)*nbands + i
		sumjk = 0.
		sumjn = 0.
		sumjnn = 0.
		do k = 1,nvalues(ib(i))
		    sumjk = sumjk +
     *			pint(k,i,m)*pint(k,i,n)/msigdvalue(k,ib(i))**2
		    sumjn = sumjn +
     *			pintx(k,i,n)*pint(k,i,m)/msigdvalue(k,ib(i))**2
		    sumjnn = sumjnn +
     *			pinty(k,i,n)*pint(k,i,m)/msigdvalue(k,ib(i))**2
		enddo
		G(jj,kk) = sumjk
		G(kk,jj) = G(jj,kk)
		G(jj,nc) = -flux(ib(i),n)*sumjn
		G(nc,jj) = G(jj,nc)
		G(jj,mblend+nc) = -flux(ib(i),n)*sumjnn
		G(mblend+nc,jj) = G(jj,mblend+nc)
	    enddo
 
	   endif
	  enddo
	 endif
	enddo
 
! Now add in the elements corresponding to proper motion.
	G(ng0+1, ng0+1) = sum(sumxxtt)
	G(ng0+2, ng0+2) = sum(sumyytt)
	G(ng0+1, ng0+2) = sum(sumxytt)
	G(ng0+2, ng0+1) = G(ng0+1, ng0+2)
	G(1, ng0+1) = sum(sumxxt)
	G(ng0+1, 1) = G(1, ng0+1)
	G(1, ng0+2) = sum(sumxyt)
	G(ng0+2, 1) = G(1, ng0+2)
	G(mblend+1, ng0+1) = G(1, ng0+2)
	G(ng0+1, mblend+1) = G(mblend+1, ng0+1)
	G(mblend+1, ng0+2) = sum(sumyyt)
	G(ng0+2, mblend+1) = G(mblend+1, ng0+2)
 
! Invert the matrix
	call matinv(G,gamma,ng,ok)
	do n = 1,ng
	    if((.not.ok .or. gamma(n,n) <= 0. .or. isnan(gamma(n,n)))
     *		.and. G(n,n) /= 0.) gamma(n,n) = 1./G(n,n)
	enddo
 
! Extract uncertainty values from the inverse matrix.
	nc = 0
	do nn = 1,nblend
	  if (goodcompz(nn)) then
	    nc = nc + 1
	    sigx(nn) = 3600.*sqrt(max(gamma(nc,nc),0.))
	    sigy(nn) = 3600.*sqrt(max(gamma(mblend+nc,mblend+nc),0.))
	    if (ok) then
		sigxy(nn) = 3600.*sqrt(abs(gamma(nc,mblend+nc)))
		sigxy(nn) = min(sigxy(nn), sqrt(sigx(nn)*sigy(nn)/2.))
		if (gamma(nc,mblend+nc) < 0) sigxy(nn) = -sigxy(nn)
	    endif
	    do i = 1,nbands
		j = 2*mblend + (nc-1)*nbands + i
		sigflux(ib(i),nn) = sqrt(max(gamma(j,j),0.))
	    enddo
	  endif
	enddo
c       if (ndump .lt. 100) print *,'unspew(474): gamma(ng0+1,ng0+1) =',gamma(ng0+1,ng0+1) ! JWF dbg
c       if (ndump .lt. 100) print *,'unspew(475): gamma(ng0+2,ng0+2) =',gamma(ng0+2,ng0+2) ! JWF dbg
c       if (ndump .lt. 100) print *,'unspew(476): minPMsig =', minPMsig  ! JWF dbg
c	sigpmr = 3600.*365.25*sqrt(max(gamma(ng0+1,ng0+1),0.)) ! arcsec/yr; JWF B21207
c	sigpmd = 3600.*365.25*sqrt(max(gamma(ng0+2,ng0+2),0.)) ! arcsec/yr; JWF B21207
        if (IzBad(gamma(ng0+1,ng0+1))) then        ! JWF B21207
          sigpmr = 999.9                           ! JWF B21207
        else                                       ! JWF B21207
          sigpmr = 3600.*365.25*sqrt(max(gamma(ng0+1,ng0+1),0.0))   ! arcsec/yr; JWF B21207
c         if (sigpmr .lt. minPMsig) sigpmr = minPMsig               ! JWF B21210
          sigpmr = sqrt(sigpmr**2 + minPMsig**2)                    ! JWF B30313
        end if                                     ! JWF B21207
        if (IzBad(gamma(ng0+2,ng0+2))) then        ! JWF B21207
          sigpmd = 999.9                           ! JWF B21207
        else                                       ! JWF B21207
	  sigpmd = 3600.*365.25*sqrt(max(gamma(ng0+2,ng0+2),0.0))   ! arcsec/yr; JWF B21207
c         if (sigpmd .lt. minPMsig) sigpmd = minPMsig               ! JWF B21210
          sigpmd = sqrt(sigpmd**2 + minPMsig**2)                    ! JWF B30313
        end if
c================== start of code added by JWF B21221 =====================
        if (IzBad(gamma(1,ng-1))) then              ! RA/PMRA error covariance
         covRpmR = 0.0
        else
          covRpmR = gamma(1,ng-1)
        end if
        if (IzBad(gamma(2,ng))) then                ! Dec/PMDec error covariance
         covDpmD = 0.0
        else
          covDpmD = gamma(2,ng)
        end if
        if (IzBad(gamma(ng-1,ng-1))) then           ! PMRA/PMDec error covariance
         covpmRD = 0.0
        else
          covpmRD = gamma(ng-1,ng-1)
        end if
c==================== end of code added by JWF B21221 =====================
c
c================== start of code added by JWF B21010 =====================
c       go to 200
        if (debug)  then
c         ndump = ndump + 1
c         if (ndump .eq. 1) open (40, file = 'ErrorCovariances.txt')
c         write (40,'(''Source #'',I5,''; ng ='',I4)') ndump, ng
          write (6, '(''; ng ='',I4)') ng
          write(6,'(a)') 'Error Covariance Matrix'
          do 120 nn = 1, ng
            write (6,'(1p99e10.2)') (gamma(j,nn), j = 1, ng)
120       continue
          write (6,'(''-----------------------------------'')')
          write(6,'(a)') 'Correlation Matrix'
          do 140 nn = 1, ng
            do 130 j = 1, ng
              G(j,nn) = gamma(j,nn)/sqrt(gamma(j,j)*gamma(nn,nn))
130         continue
140       continue
c         write (40,'(''Correlation:'')')
          do 150 nn = 1, ng
c           write (40,'(1p99e10.2)') (G(j,nn), j = 1, ng)
            write (6, '(1p99e10.2)') (G(j,nn), j = 1, ng)
150       continue
c         write (40,'(''-----------------------------------'')')
          write (6,'(''-----------------------------------'')')
        end if
200     continue
c==================== end of code added by JWF B21010 =====================
 
	deallocate(xi)
	deallocate(eta)
	deallocate(xsource)
	deallocate(ysource)
	deallocate(xsourcep)
	deallocate(ysourcep)
	deallocate(xsourceq)
	deallocate(ysourceq)
	deallocate(G)
	deallocate(gamma)
	deallocate(sumxx)
	deallocate(sumyy)
	deallocate(sumxy)
	deallocate(sumxxt)
	deallocate(sumyyt)
	deallocate(sumxyt)
	deallocate(sumxxtt)
	deallocate(sumyytt)
	deallocate(sumxytt)
	return
 
	end
