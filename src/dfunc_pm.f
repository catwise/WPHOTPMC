c	subroutine dfunc(p,psfs,npsize,npnmax,grad,fp)       ! JWF B21207
	subroutine dfunc_pm(p,psfs,npsize,npnmax,grad,fp,ok2)   ! JWF B21207
 
! Calculate the gradient of the chi squared function for a blend of NBLEND
! components located at positions given by the vector P.
 
	implicit real*4 (a-h,o-z)
	implicit integer (i-n)
	include 'wpro.inc'
	integer, parameter :: maxblend = 100
	integer, parameter :: maxblend2 = 200
	integer, parameter :: maxbands = 4
	real, parameter :: dtor = 0.0174533
 
	real(4), allocatable :: pint(:,:)
	real(4), allocatable :: xi(:), eta(:), model(:)
	real(4), allocatable :: xsource(:,:), ysource(:,:)
	real(4), allocatable :: A(:,:), b(:)                                  ! JWF B30506
	real(4), allocatable :: AA(:,:), Ai(:,:)                              ! JWF B30506
	real(4), allocatable :: sumb(:), sumbp(:), fluxb(:)
	integer, allocatable :: npointsb(:)
 
c	real(4) model(maxpix)                                                 ! JWF B30506
c	real(4) xsource(maxblend,maxframes), ysource(maxblend,maxframes)
c	real(4) xi(maxblend), eta(maxblend)
c	real(4) A(maxblend,maxblend), b(maxblend)                             ! JWF B30501
c	real(4) AA(maxblend,maxblend), Ai(maxblend,maxblend)                  ! JWF B30501
c	real(4) sumb(maxbands), sumbp(maxbands), fluxb(maxblend)              ! JWF B30506
c	integer npointsb(maxbands)
 
	real(8) ra8, dec8, x8, y8
	real(8) racan(maxblend),deccan(maxblend)
	real(4) p(*),grad(*)
	real(4) dvalue(maxpix,maxbands),sigdvalue(maxpix,maxbands)
	real(4) tvalue(maxpix,maxbands)
	real(4) pxc(maxbands),pyc(maxbands),psamp(maxbands)
	real(4) psfs(npsize,npsize,npnmax,4)
	real(4) xpixscale(maxframes),ypixscale(maxframes),tframe(maxframes)
	real(4) flux(maxbands,maxblend),rchisqb(maxbands)
	real(4) xsourcep(maxblend,maxframes)
	real(4) ysourcep(maxblend,maxframes)
	real(4) xsourceq(maxblend,maxframes)
	real(4) ysourceq(maxblend,maxframes)
	real(4) pintx(maxpix,maxblend),pinty(maxpix,maxblend)
	real(4) gsum(maxbands,maxblend2)
	integer ivalue(maxpix,maxbands),jvalue(maxpix,maxbands)
	integer ifvalue(maxpix,maxbands),primask(maxpix,maxbands)
	integer ib(maxbands),nvalues(maxbands)
	integer WCS(maxframesin,maxbands), offscl, nnpsf(maxbands)
	integer ifq(maxframes),ibq(maxframes)
	integer kpn(maxframes)
	logical(4) ok,goodcomp(maxblend)
        integer*4 never     !  JWF dbg
        data      never/0/  !  JWF dbg
        logical*4 ok2       !  JWF B21207
 
	common /funccom/ racan,deccan,nvalues,ivalue,jvalue,dvalue,
     *	    sigdvalue,ifvalue,WCS,ifq,ibq,primask,nframes,
     *	    kpn,nnpsf,pxc,pyc,psamp,nblend,goodcomp,icp,
     *	    xpixscale,ypixscale,flux,rchisq,rchisqb,tframe
c
        include 'jwfcom.f'                       ! JWF B30507
c
	allocate(xi(nblend))
	allocate(eta(nblend))
	allocate(xsource(nblend,nframes))
	allocate(ysource(nblend,nframes))
 
        ok2 = .true.  ! JWF B21207
 
! Unpack the parameter vector.
	mblend = ntrue(goodcomp,nblend)
	xi = 0.
	eta = 0.
	nc = 0
	do n = 1,nblend
	    if (goodcomp(n)) then
		nc = nc + 1
		xi(n) = p(nc)			! offsets east and north from
		eta(n) = p(nc+mblend)		! candidate location [deg]
	    endif
	enddo
	yr = 365.25
	pmr = p(2*mblend+1)/yr	! deg/day
	pmd = p(2*mblend+2)/yr	! deg/day
 
	step = 0.1				! step size for gradient
						! calculation [pixels]
 
! Which bands are present?
	nbands = 0	
	do i = 1,maxbands
	    if(nvalues(i) > nblend+2) then
		nbands = nbands + 1
		ib(nbands) = i
	    endif
	enddo
 
! Find pixel locations of source at offset (xi,eta).
	do iframe = 1,nframes
	    iwcs = WCS(ifq(iframe),ibq(iframe))
	    ddec = step*xpixscale(iframe)
	    dra = ddec/cos(dtor*deccan(1))
 
	    do n = 1,nblend
        	ra8 = racan(n) + xi(n)/cos(dtor*deccan(n))
        	dec8 = deccan(n) + eta(n)
		if (n==icp) then
		    ra8 = ra8 + pmr*tframe(iframe)/cos(dtor*deccan(n))
		    dec8 = dec8 + pmd*tframe(iframe)
		endif
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
        	call wcs2pix(iwcs, ra8, dec8+ddec, x8, y8, offscl)
        	if (offscl .ne. 0) then
		    xsourceq(n,iframe) = 1.
		    ysourceq(n,iframe) = 1.
        	else
		    xsourceq(n,iframe) = x8
		    ysourceq(n,iframe) = y8
		endif
	    enddo
	enddo
 
! Estimate flux at each band.
	allocate (A(mblend,mblend))
	allocate (b(mblend))
	allocate (AA(mblend,mblend))
	allocate (Ai(mblend,mblend))
	allocate (fluxb(mblend))
	allocate (sumb(nbands))
	allocate (sumbp(nbands))
	allocate (npointsb(nbands))
 
	npointsb = 0
 
	do i = 1,nbands
	    nv = nvalues(ib(i))
 
! Interpolate the PSF.
	    allocate (pint(nv,nblend))
	    allocate (model(nv))
 
	    do k = 1,nv
		iframe = ifvalue(k,ib(i))
		tvalue(k,ib(i)) = tframe(iframe)
		pmin = pminFac*psfs(nint(pxc(ib(i))),nint(pyc(ib(i))),  ! JWF B30507
     *		    kpn(iframe),ib(i))
	        xprat = xpixscale(iframe)/psamp(ib(i))
	        yprat = ypixscale(iframe)/psamp(ib(i))
		dtheta = step*xpixscale(iframe)
 
	    	do n = 1,nblend
		    xpsample = pxc(ib(i)) +
     *			(ivalue(k,ib(i)) - xsource(n,iframe))*xprat
		    ypsample = pyc(ib(i)) +
     *			(jvalue(k,ib(i)) - ysource(n,iframe))*yprat
		    pint(k,n) = xprat*yprat*max(bilint(psfs,npsize,
     *			npsize,npnmax,4,nnpsf(ib(i)),nnpsf(ib(i)),
     *			kpn(iframe),ib(i),xpsample,ypsample), pmin)
		    xpsample = pxc(ib(i)) +
     *			(ivalue(k,ib(i)) - xsourcep(n,iframe))*xprat
		    ypsample = pyc(ib(i)) +
     *			(jvalue(k,ib(i)) - ysourcep(n,iframe))*yprat
		    pintx(k,n) =(xprat*yprat*max(bilint(psfs,npsize,
     *			npsize,npnmax,4,nnpsf(ib(i)),nnpsf(ib(i)),
     *			kpn(iframe),ib(i),xpsample,ypsample),pmin)
     *			- pint(k,n))/dtheta
		    xpsample = pxc(ib(i)) +
     *			(ivalue(k,ib(i)) - xsourceq(n,iframe))*xprat
		    ypsample = pyc(ib(i)) +
     *			(jvalue(k,ib(i)) - ysourceq(n,iframe))*yprat
		    pinty(k,n) = (xprat*yprat*max(bilint(psfs,npsize,
     *			npsize,npnmax,4,nnpsf(ib(i)),nnpsf(ib(i)),
     *			kpn(iframe),ib(i),xpsample,ypsample),pmin)
     *			- pint(k,n))/dtheta
		enddo
	    enddo
 
! Calculate the matrix elements.
	    nc = 0
	    do n = 1,nblend
	      if (goodcomp(n)) then
		nc = nc + 1
		b(nc) = 0.0
		do k = 1,nv
		    b(nc) = b(nc) + dvalue(k,ib(i)) *
     *			pint(k,n)/sigdvalue(k,ib(i))**2
		enddo
		mc = 0
		do m = 1,nblend
		  if (goodcomp(m)) then
		    mc = mc + 1
		    A(mc,nc) = 0.0
		    do k = 1,nv
			A(mc,nc) = A(mc,nc) +
     *			    pint(k,m)*pint(k,n)/sigdvalue(k,ib(i))**2
		    enddo
		  endif
		enddo
	      endif
	    enddo
 
	    AA = matmul(transpose(A), A)
	    call matinv(AA,Ai,mblend,ok)
 
	    if (ok) then
		fluxb = matmul(matmul(Ai,transpose(A)), b)
	    else
		print *,'DFUNC: singular matrix'
		fluxb = 0.0
	    endif
 
	    nc = 0
	    do n = 1,nblend
	      if (goodcomp(n)) then
		nc = nc + 1
		fluxb(nc) = max(fluxb(nc), 0.0)
		flux(ib(i),n) = fluxb(nc)
	      endif
	    enddo
 
! Calculate residuals.
	    do k = 1,nv
		model(k) = 0.0
	    enddo
	    nc = 0
	    do n = 1,nblend
	      if (goodcomp(n)) then
		nc = nc + 1
		do k = 1,nv
		    model(k) = model(k) + fluxb(nc)*pint(k,n)
		enddo
	      endif
	    enddo
 
	    sumb(i) = sum(((dvalue(1:nv, ib(i)) - model(1:nv))/
     *		sigdvalue(1:nv,ib(i)))**2)
	    sumbp(i) = sum((primask(1:nv,ib(i)) *
     *		(dvalue(1:nv,ib(i)) - model(1:nv))/
     *		sigdvalue(1:nv,ib(i)))**2)
 
	    nc = 0
            ncp = -9 !  JWF B21126
c           print *,'dfunc(245): icp, nblend =', icp, nblend  !  JWF dbg
c           print *,'dfunc(246): goodcomp =', goodcomp(1:nblend)  !  JWF dbg
	    do n = 1,nblend
	      if (goodcomp(n)) then
		nc = nc + 1
		if (n==icp) ncp = nc
		gsum(i,nc) = -2.*fluxb(nc)*sum(pintx(1:nv,n)*
     *		    (dvalue(1:nv,ib(i)) - model(1:nv))/
     *		    sigdvalue(1:nv,ib(i))**2)
		gsum(i,mblend+nc) = -2.*fluxb(nc)*sum(pinty(1:nv,n)*
     *		    (dvalue(1:nv,ib(i)) - model(1:nv))/
     *		    sigdvalue(1:nv,ib(i))**2)
	      endif
	    enddo
            if (ncp .gt. 0) then             ! JWF B21127
 	      gsum(i,2*mblend+1) = -2.*fluxb(ncp)*sum(pintx(1:nv,icp)*
     *		    (dvalue(1:nv,ib(i)) - model(1:nv))*tvalue(1:nv,ib(i))/
     *		    sigdvalue(1:nv,ib(i))**2)
	      gsum(i,2*mblend+2) = -2.*fluxb(ncp)*sum(pinty(1:nv,icp)*
     *		    (dvalue(1:nv,ib(i)) - model(1:nv))*tvalue(1:nv,ib(i))/
     *		    sigdvalue(1:nv,ib(i))**2)
            else                             ! JWF B21127
              gsum(i,2*mblend+1) = 0.0       ! JWF B21127
              gsum(i,2*mblend+2) = 0.0       ! JWF B21127
              ok2 = .false.                  ! JWF B21207
c             print *,'ERROR (dfunc): "(n==icp)" test never T'  ! JWF B21127
c             call exit(64)                  ! JWF B21127
              never = never + 1              ! JWF dbg
c             print *,'(dfunc): "(n==icp)" test never T, occurence #',never ! JWF dbg
            end if                           ! JWF B21127
c
	    npointsb(i) = sum(primask(1:nv,ib(i)))
 
	    deallocate(pint)
	    deallocate(model)
	enddo
 
	rchisq = sum(sumbp)/max((sum(npointsb) - (nbands+2)), 1)
	do i = 1,nbands
	    rchisqb(ib(i)) = sumbp(i)/max((npointsb(i) - 3), 1)
	enddo
	sigpos = maxval(xpixscale)/2.
c	sigpm = sigpos/10.	! NEW  [may or may not do this]
c	dt = maxval(tframe(1:nframes)) - minval(tframe(1:nframes)) ! JWF B21229 [Ken says not needed]
	fp = sum(sumb) + sum((xi/sigpos)**2 + (eta/sigpos)**2)
c    *     + (pmr*yr/sigpm)**2 + (pmd*yr/sigpm)**2	! NEW  [may or may not do this]
	nc = 0
 
	do n = 1,nblend
	  if (goodcomp(n)) then
	    nc = nc + 1
	    grad(nc) = sum(gsum(1:nbands,nc)) + 2.*xi(nc)/sigpos**2
	    grad(mblend+nc) = sum(gsum(1:nbands,mblend+nc)) +
     *		2.*eta(nc)/sigpos**2
	  endif
	enddo
	grad(2*mblend+1) = sum(gsum(1:nbands,2*mblend+1))
	grad(2*mblend+2) = sum(gsum(1:nbands,2*mblend+2))
c	grad(2*mblend+1)=sum(gsum(1:nbands,2*mblend+1))+2.*pmr*yr/sigpm**2 ! NEW [may or may not do this]
c	grad(2*mblend+2)=sum(gsum(1:nbands,2*mblend+2))+2.*pmd*yr/sigpm**2 ! NEW [may or may not do this]
 
	if (isnan(fp)) then
	    fp = 0.
	    nc = 0
	    do n = 1,nblend
	      if (goodcomp(n)) then
		nc = nc + 1
		grad(nc) = 2.*xi(nc)/sigpos**2
		grad(mblend+nc) = 2.*eta(nc)/sigpos**2
	      endif
	    enddo
 	    grad(2*mblend+1:2*mblend+2) = 0.
c	    grad(2*mblend+1) = 2.*pmr*yr/sigpm**2	! NEW [may or may not do this]
c	    grad(2*mblend+2) = 2.*pmd*yr/sigpm**2	! NEW [may or may not do this]
	endif
 
	deallocate(A)
	deallocate(b)
	deallocate(AA)
	deallocate(Ai)
	deallocate(xi)
	deallocate(eta)
	deallocate(xsource)
	deallocate(ysource)
	deallocate(fluxb)
	deallocate(sumb)
	deallocate(sumbp)
	deallocate(npointsb)
 
	return
 
	end
