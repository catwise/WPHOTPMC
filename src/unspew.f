	subroutine unspew(p,psfs,npsize,npnmax,varpsf,nlayers,isn,onframe,
     *	    nmax,nf,mxbands,wflag,nfi,framesig,sigx,sigy,sigxy,sigflux,ok,debug)
 
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
c	integer(2) bset(maxframes),nlayers(maxframes),wflag(nfi,4)  ! JWF B30207
c	integer(2) nlayers(maxframes),wflag(nfi,4)                  ! JWF B30207, B60714
	integer(2) nlayers(maxbands),wflag(nfi,4)                   ! JWF B60714
        logical*4  debug                                            ! JWF B30510
        real*4     dra, ddec(maxbands)                              ! JWF B30513
	logical(4) ok,goodcomp(maxblend),goodcompz(maxblend)
	logical(1) onframe(nmax,nf,mxbands)
 
	common /funccom/ racan,deccan,nvalues,ivalue,jvalue,dvalue,
     *	    sigdvalue,ifvalue,WCS,ifq,ibq,primask,nframes,kpn,
     *	    nnpsf,pxc,pyc,psamp,nblend,goodcomp,icp,
     *	    xpixscale,ypixscale,flux,rchisq,rchisqb,tframe
c
        include 'jwfcom.f'                       ! JWF B30507
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
	ng = (nbands + 2)*mblend
	allocate (gamma(ng,ng))
 
	step = 0.1		! step size for gradient calculation [pixels]

        if (debug) then
          print *,'unspew: ================================================' ! JWF dbg
          print *,'   RA,Dec:', racan(1), deccan(1) ! JWF dbg
          print *,'   mblend, nblend, nbands:', mblend, nblend, nbands
          print *,'   goodcomp: ', goodcomp(1:nblend)
          print *,'   goodcompz:', goodcompz(1:nblend)
        end if
 
! Find pixel locations of source at offset (xi,eta).
	do iframe = 1,nframes
	    iwcs = WCS(ifq(iframe),ibq(iframe))
	    ddec(ibq(iframe)) = step*xpixscale(iframe)
	    dra = ddec(ibq(iframe))/cos(dtor*deccan(1))
            if (debug) then
              write(6,'(a,3i5,1p2e11.3)') 'unspew(114):',
     +              ibq(iframe), iframe, ifq(iframe),
     +              dra, ddec(ibq(iframe))
              write(6,'(a,3i5,1p2e11.3)') 'unspew(117):',
     +              iframe, ibq(iframe), ifq(iframe),
     +              xpixscale(iframe), ypixscale(iframe)
              write(6,'(a,3i5,1p2e11.3)') 'unspew(120):',
     +              iframe, ibq(iframe), ifq(iframe),
     +              xpixscale(ifq(iframe)), ypixscale(ifq(iframe))
            end if

	    do n = 1,nblend
        	ra8 = racan(n) + xi(n)/cos(dtor*deccan(n))
        	dec8 = deccan(n) + eta(n)
        	offscl = -1
        	call wcs2pix(iwcs, ra8, dec8, x8, y8, offscl)
        	if (offscl .ne. 0) then
		    xsource(n,iframe) = 1.
		    ysource(n,iframe) = 1.
        	else
		    xsource(n,iframe) = x8
		    ysource(n,iframe) = y8
		endif
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
                if (debug) then                         ! JWF B30512
                  write(6,'(a,i2,2I4,1p4e11.3)') 'unspew(165):',
     +            ib(i), iframe, ifq(iframe), pmin, xprat, yprat, rn
                end if
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
		enddo ! n = 1,nblend
                if (debug) then                            ! JWF B30512
                  write (6,'(a,I2,2I5,1p3e11.3)') 'unspew(183):',
     +             ib(i), ivalue(k,ib(i)), jvalue(k,ib(i)),
     +             pintx(k,i,1), pinty(k,i,1), pint(k,i,1)
                end if
	    enddo  ! k = 1,nvalues(ib(i))
	enddo ! i = 1, nbands
 
	allocate (G(ng,ng))
	allocate (sumxx(nbands))
	allocate (sumyy(nbands))
	allocate (sumxy(nbands))
	G = 0.
	gamma = 0.

! Start by calculating the flux sigmas for individual frames for the benefit
! of N out of M statistics.
	do ifs = 1,nf		! begin frameset loop
	 framesig(ifs,1:4,1:nblend) = 0.
	 gamma = 0.
	 if (any(onframe(isn,ifs,1:4) .and. wflag(ifs,1:4)==1)) then
          if (debug) print *,'unspew(219): doing frame',ifs,'ng =',ng     ! JWF B30702
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
		do k = 1,nvalues(ib(i))
		  if(ifset(k,i)==ifs) then
		    sumHxx = sumHxx + pintx(k,i,m)*pintx(k,i,n)/
     *			sigdvalue(k,ib(i))**2        ! BUGFIX 3/18/13
		    sumHyy = sumHyy + pinty(k,i,m)*pinty(k,i,n)/
     *			sigdvalue(k,ib(i))**2        ! BUGFIX 3/18/13
		    sumHxy = sumHxy + pintx(k,i,m)*pinty(k,i,n)/
     *			sigdvalue(k,ib(i))**2        ! BUGFIX 3/18/13
		  endif
		enddo
		sumxx(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxx
		sumyy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyy
		sumxy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxy
	    enddo
	    G(mc,nc) = sum(sumxx)
	    G(nc,mc) = G(mc,nc)
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
     *			pint(k,i,m)*pint(k,i,n)/sigdvalue(k,ib(i))**2  ! BUGFIX 3/18/13
		    sumjn = sumjn +
     *			pintx(k,i,n)*pint(k,i,m)/sigdvalue(k,ib(i))**2 ! BUGFIX 3/18/13
		    sumjnn = sumjnn +
     *			pinty(k,i,n)*pint(k,i,m)/sigdvalue(k,ib(i))**2 ! BUGFIX 3/18/13
		  endif
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
	  call matinv(G,gamma,ng,ok)
          if (debug) then
            print *,'unspew(281): G matrix:'                       ! JWF B30702
            do kk = 1, ng
              write(6,'(1pE13.5)') (G(jj,kk), jj = 1, ng)
            end do
            if (.not.ok) print *,'unspew(282): matrix is singular' ! JWF B30702
            print *,'unspew(286): gamma matrix:'                   ! JWF B30702
            do kk = 1, ng
              write(6,'(1pE13.5)') (gamma(jj,kk), jj = 1, ng)
            end do
          end if ! debug
	  do n = 1,ng
	    if((.not.ok .or. gamma(n,n) <= 0. .or. isnan(gamma(n,n)))
     *		.and. G(n,n) /= 0.) gamma(n,n) = 1./G(n,n)
            if (debug) then
              if((.not.ok .or. gamma(n,n) <= 0. .or. isnan(gamma(n,n)))
     *	          .and. G(n,n) /= 0.)
     +           print *,'n =',n,'; gamma(n,n) = 1./G(n,n) =',gamma(n,n)
            end if
	  enddo
	  nc = 0
	  do nn = 1,nblend
	   if (goodcompz(nn)) then
	    nc = nc + 1
	    do i = 1,nbands
		n = 2*mblend + (nc-1)*nbands + i
		framesig(ifs,ib(i),nn) = sqrt(max(gamma(n,n),0.))
                if (debug) print *,'unspew(295): framesig(ifs,ib(i),nn) =',framesig(ifs,ib(i),nn) !JWF B30702
	    enddo
	   endif  !  (goodcompz(nn))
	  enddo  !  nn = 1,nblend
	 endif  !  (any(onframe(isn,ifs,1:4) .and. wflag(ifs,1:4)==1))
	enddo		! end frameset loop, ifs = 1,nf

! Now calculate matrix elements.
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
                if (debug .and. (nc .eq. 1) .and. (mc .eq. 1)) then
                  print *,'band',ib(1),
     +                    'sumHxx, sumHyy, sumHxy summation terms and factors'
                end if !  debug
		do k = 1,nvalues(ib(i))
		    sumHxx = sumHxx + pintx(k,i,m)*pintx(k,i,n)/
     *			msigdvalue(k,ib(i))**2
		    sumHyy = sumHyy + pinty(k,i,m)*pinty(k,i,n)/
     *			msigdvalue(k,ib(i))**2
		    sumHxy = sumHxy + pintx(k,i,m)*pinty(k,i,n)/
     *			msigdvalue(k,ib(i))**2
                    if (debug .and. (nc .eq. 1) .and. (mc .eq. 1)) then
                      print *,'k, i, m, n:',k, i, m, n
                      write(6,'(a,1p5e11.3)')
     +                'msigdvalue,pintx(k,i,m),pintx(k,i,n),'
     +                         //'pinty(k,i,m),pinty(k,i,n):',
     +                 msigdvalue(k,ib(i)),pintx(k,i,m),pintx(k,i,n),
     +                                     pinty(k,i,m),pinty(k,i,n)
c
c                     print *,'k, i, m, n, msigdvalue(k,ib(i)):',
c    +                         k,i,m,n,msigdvalue(k,ib(i))
c                     print *,'sumHxx,pintx(k,i,m),pintx(k,i,n):',
c    +                         sumHxx,pintx(k,i,m),pintx(k,i,n)
c                     print *,'sumHyy,pinty(k,i,m),pinty(k,i,n):',
c    +                         sumHyy,pinty(k,i,m),pintx(k,i,n)
c                     print *,'sumHxy,pintx(k,i,m),pinty(k,i,n):',
c    +                         sumHxy,pintx(k,i,m),pinty(k,i,n)
                    end if  !  debug
		enddo  !  k = 1,nvalues(ib(i))
		sumxx(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxx
		sumyy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHyy
		sumxy(i) = flux(ib(i),m)*flux(ib(i),n)*sumHxy
                if (debug .and. (nc .eq. 1) .and. (mc .eq. 1)) then
                  print *,'flux(ib(i),m),flux(ib(i),n):',flux(ib(i),m),flux(ib(i),n)
                  print *,'sumxx, sumyy, sumxy, i:',sumxx(i), sumyy(i), sumxy(i), i
                end if  !  debug
	    enddo  !  i = 1,nbands
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
		enddo  !  k = 1,nvalues(ib(i))
		G(jj,kk) = sumjk
		G(kk,jj) = G(jj,kk)
		G(jj,nc) = -flux(ib(i),n)*sumjn
		G(nc,jj) = G(jj,nc)
		G(jj,mblend+nc) = -flux(ib(i),n)*sumjnn
		G(mblend+nc,jj) = G(jj,mblend+nc)
	    enddo  !  i = 1,nbands

	   endif  !  (goodcompz(m))
	  enddo  !  m = 1,nblend
	 endif  !  (goodcompz(n))
	enddo  !  n = 1,nblend
c
c================== start of code added by JWF B30510 =====================
        if (debug)  then
          print *,'unspew(339) -------------------------------'
          write (6, '(''; ng ='',I4)') ng
          write(6,'(a)') 'Fisher Information Matrix'
          do 20 nn = 1, ng
            write (6,'(1p99e10.2)') (G(j,nn), j = 1, ng)
20        continue
          write (6,'(''-----------------------------------'')')
        end if
c==================== end of code added by JWF B30510 =====================
c
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
c
c================== start of code added by JWF B30510 =====================
        if (debug)  then
          print *,'unspew(364) -------------------------------'
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
          do 150 nn = 1, ng
            write (6, '(1p99e10.2)') (G(j,nn), j = 1, ng)
150       continue
          write (6,'(''-----------------------------------'')')
        end if
c==================== end of code added by JWF B30510 =====================
c
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
	return
 
	end
