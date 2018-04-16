	real*4 function func(p,psfs,npsize,npnmax,nDoF,dbg)
 
! Calculate the reduced chi squared for blend of NBLEND components located
! at positions given by the vector P.
 
	implicit real*4 (a-h,o-z)
	implicit integer (i-n)
	include 'wpro.inc'
	integer, parameter :: maxblend = 100
	integer, parameter :: maxpsize =  256
	integer, parameter :: maxbands = 4
	integer, parameter :: maxsegs = 121
	real, parameter :: dtor = 0.0174533
 
	real(4), allocatable :: pint(:,:)
 	real(4), allocatable :: A(:,:), b(:)                                  ! JWF B30506
 	real(4), allocatable :: AA(:,:), Ai(:,:)                              ! JWF B30506
	real(4), allocatable :: xi(:), eta(:), model(:)
	real(4), allocatable :: xsource(:,:), ysource(:,:)
	real(4), allocatable :: sumb(:), sumbp(:), fluxb(:)
	integer, allocatable :: npointsb(:)
 
c	real(4) model(maxpix)                                                 ! JWF B30506
c	real(4) xsource(maxblend,maxframes), ysource(maxblend,maxframes)
c	real(4) xi(maxblend), eta(maxblend)
c       real(4) A(maxblend,maxblend), b(maxblend)                             ! JWF B30501
c	real(4) AA(maxblend,maxblend), Ai(maxblend,maxblend)                  ! JWF B30501
c	real(4) sumb(maxbands), sumbp(maxbands), fluxb(maxblend)              ! JWF B30506
c	integer npointsb(maxbands)
 
	real(8) ra8, dec8, x8, y8
	real(8) racan(maxblend),deccan(maxblend)
	real(4) p(*)
	real(4) psfs(npsize,npsize,npnmax,4)
	real(4) dvalue(maxpix,maxbands),sigdvalue(maxpix,maxbands)
	real(4) pxc(maxbands),pyc(maxbands),psamp(maxbands)
	real(4) xpixscale(maxframes),ypixscale(maxframes),tframe(maxframes)
	real(4) flux(maxbands,maxblend),rchisqb(maxbands)
	integer ivalue(maxpix,maxbands),jvalue(maxpix,maxbands)
	integer ifvalue(maxpix,maxbands),primask(maxpix,maxbands)
	integer ib(maxbands),nvalues(maxbands)
	integer WCS(maxframesin,maxbands), offscl, nnpsf(maxbands)
	integer ifq(maxframes),ibq(maxframes),kpn(maxframes)
	logical(4) posconstrain,ok,goodcomp(maxblend)
        logical*4 dbg  ! JWF B30315
        integer*4 nDoF ! JWF B30419
 
	common /funccom/ racan,deccan,nvalues,ivalue,jvalue,dvalue,
     *	    sigdvalue,ifvalue,WCS,ifq,ibq,primask,nframes,kpn,
     *	    nnpsf,pxc,pyc,psamp,nblend,goodcomp,icp,
     *	    xpixscale,ypixscale,flux,rchisq,rchisqb,tframe
	common /poscom/ posconstrain
c
        include 'jwfcom.f'                       ! JWF B30507
c
	allocate(xi(nblend))
	allocate(eta(nblend))
	allocate(xsource(nblend,nframes))
	allocate(ysource(nblend,nframes))
 
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
 
        if (dbg) then
          print *,'============================='
          print *,'func(63): nblend, mblend =',nblend, mblend ! JWF dbg
        end if
 
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
                if (dbg) then
                  print *,'func(89): n, iframe =',n,iframe ! JWF dbg
                  print *,'          ra8, dec8 =',ra8,dec8  ! JWF dbg
                  print *,'          xsource, ysource:', xsource(n,iframe),ysource(n,iframe) ! JWF dbg
                end if
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
 
	    allocate (pint(nv,nblend))
	    allocate (model(nv))
 
! Interpolate the PSF.
 
	    do k = 1,nv
		iframe = ifvalue(k,ib(i))
		pmin = pminFac*psfs(nint(pxc(ib(i))),nint(pyc(ib(i))),  ! JWF B30507
     *		    kpn(iframe),ib(i))
	        xprat = xpixscale(iframe)/psamp(ib(i))
	        yprat = ypixscale(iframe)/psamp(ib(i))
 
	    	do n = 1,nblend
		    xpsample = pxc(ib(i)) +
     *			(ivalue(k,ib(i)) - xsource(n,iframe))*xprat
		    ypsample = pyc(ib(i)) +
     *			(jvalue(k,ib(i)) - ysource(n,iframe))*yprat
		    pint(k,n) = xprat*yprat*max(bilint(psfs,npsize,
     *			npsize,npnmax,4,nnpsf(ib(i)),nnpsf(ib(i)),
     *			kpn(iframe),ib(i),xpsample,ypsample), pmin)
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
		print *,'FUNC: singular matrix'
		fluxb = 0.0
	    endif
 
	    nc = 0
	    do n = 1,nblend
		if (goodcomp(n)) then
		    nc = nc + 1
		    if(posconstrain) fluxb(nc) = max(fluxb(nc), 0.0)
		    flux(ib(i),n) = fluxb(nc)
	 	endif
                if (dbg) then
                  print *,'func(179): n, ib(i), flux(ib(i),n):', ! JWF dbg
     +            n,ib(i),flux(ib(i),n) 
                end if
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
	    sumbp(i) = sum((primask(1:nv,ib(i))*(dvalue(1:nv,ib(i))
     *		- model(1:nv))/sigdvalue(1:nv,ib(i)))**2)
	    npointsb(i) = sum(primask(1:nv,ib(i)))
 
	    deallocate(pint)
	    deallocate(model)
	enddo
 
        nDoF = max((sum(npointsb) - (nbands+2)), 1)
	rchisq = sum(sumbp)/nDoF
	do i = 1,nbands
	    rchisqb(ib(i)) = sumbp(i)/max((npointsb(i) - 3), 1)
	enddo
	sigpos = maxval(xpixscale)/2.
	func = sum(sumb) + sum((xi/sigpos)**2 + (eta/sigpos)**2)
 
	if (isnan(func)) func = 0.
 
	deallocate(xi)
	deallocate(eta)
	deallocate(xsource)
	deallocate(ysource)
 
	deallocate(A)
	deallocate(b)
	deallocate(AA)
	deallocate(Ai)
 
	deallocate(fluxb)
	deallocate(sumb)
	deallocate(sumbp)
	deallocate(npointsb)
 
	return
 
	end
