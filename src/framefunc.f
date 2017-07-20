	subroutine framefunc(p,psfs,npsize,npnmax,isn,onframe,nmax,nf,
     *	    mxbands,wflag,nfi,frameflux,framechi,dbg)
 
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
 
	real(4), allocatable :: xi(:), eta(:), model(:)
	real(4), allocatable :: xsource(:,:), ysource(:,:)
	real(4), allocatable :: pint(:,:)
	real(4), allocatable :: A(:,:), b(:),Af(:,:,:),bf(:,:)
	real(4), allocatable :: AA(:,:), Ai(:,:), fluxb(:)
	real(4), allocatable :: sumb(:), sumbp(:)
	integer, allocatable :: npointsb(:)
 
	real(8) ra8, dec8, x8, y8
	real(8) racan(maxblend),deccan(maxblend)
	real(4) p(*)
	real(4) psfs(npsize,npsize,npnmax,4)
	real(4) dvalue(maxpix,maxbands),sigdvalue(maxpix,maxbands)
	real(4) pxc(maxbands),pyc(maxbands),psamp(maxbands)
	real(4) xpixscale(maxframes),ypixscale(maxframes),tframe(maxframes)
	real(4) frameflux(maxframes,maxbands,maxblend),rchisqb(maxbands)
	real(4) framechi(maxframes,maxbands)
	integer ivalue(maxpix,maxbands),jvalue(maxpix,maxbands)
	integer ifvalue(maxpix,maxbands),primask(maxpix,maxbands),ifset(maxpix)
	integer ib(maxbands),nvalues(maxbands)
	integer WCS(maxframesin,maxbands), offscl, nnpsf(maxbands)
	integer ifq(maxframes),ibq(maxframes),kpn(maxframes)
   	integer(2) wflag(nfi,4)
	logical(4) posconstrain,ok,goodcomp(maxblend)
	logical(1) onframe(nmax,nf,mxbands)
        logical*4  dbg  ! JWF B30702
      real(4) flux(maxbands,maxblend)                ! JWF B60714

 
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
            if (dbg) print *,'framefunc(110): nv =',nv,'in band',ib(i) ! JWF B30702
 
! Interpolate the PSF.
	    allocate (pint(nv,nblend))
 
	    do k = 1,nv
		iframe = ifvalue(k,ib(i))
		ifset(k) = ifq(iframe)
		pmin = pminFac*psfs(nint(pxc(ib(i))),nint(pyc(ib(i))), ! JWF B30507
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
	    allocate (Af(nf,mblend,mblend))
	    allocate (bf(nf,mblend))
 
! Calculate the matrix elements.
	    nc = 0
	    do n = 1,nblend
	      if (goodcomp(n)) then
		nc = nc + 1
		do ifs = 1,nf
		    bf(ifs,nc) = 0.
		enddo
		do k = 1,nv
		    bf(ifset(k),nc) = bf(ifset(k),nc) + dvalue(k,ib(i)) *
     *			pint(k,n)/sigdvalue(k,ib(i))**2
		enddo
		mc = 0
		do m = 1,nblend
	 	  if (goodcomp(m)) then
		    mc = mc + 1
		    do ifs = 1,nf
		        Af(ifs,mc,nc) = 0.
		    enddo
		    do k = 1,nv
			Af(ifset(k),mc,nc) = Af(ifset(k),mc,nc) +
     *			    pint(k,m)*pint(k,n)/sigdvalue(k,ib(i))**2
		    enddo
		  endif
		enddo
	      endif
	    enddo
 
	    do ifs = 1,nf
	      if (any(onframe(isn,ifs,1:4) .and. wflag(ifs,1:4)==1)) then
		nc = 0
		do n = 1,nblend
		  if (goodcomp(n)) then
		    nc = nc + 1
		    b(nc) = bf(ifs,nc)
		    mc = 0
		    do m = 1,nblend
		      if (goodcomp(m)) then
			mc = mc + 1
			A(mc,nc) = Af(ifs,mc,nc)
		      endif
		    enddo
		  endif
		enddo
                if (dbg) then                                      ! JWF B30702
                  print *,'framefunc(181): mblend =',mblend
                  print *,'framefunc(182): A matrix for frame',ifs
                  print *,'framefunc(183):',A
                end if
		AA = matmul(transpose(A), A)
		call matinv(AA,Ai,mblend,ok)
 
		if (ok) then
		    fluxb = matmul(matmul(Ai,transpose(A)), b)
                    if (dbg) print *,'framefunc(190): good fluxb =',fluxb
		else
                    if (dbg) print *,'framefunc(192): bad fluxb = 0'
		    fluxb = 0.
                    nSingMatFrmFunc = nSingMatFrmFunc + 1                           ! JWF B30319
c                   if (nSingMatFrmFunc .le. 10)                                    ! JWF B30319
c    +              print *,'FRAMEFUNC: singular matrix'                            ! JWF B30319
c                   if (nSingMatFrmFunc .eq. 10)  print *,'FRAMEFUNC: last warning' ! JWF B30319
c    +              //' for singular matrix; check processing summary'              ! JWF B30319
		endif
 
		nc = 0
		do n = 1,nblend
		  if (goodcomp(n)) then
		    nc = nc + 1
		    if(posconstrain) fluxb(nc) = max(fluxb(nc), 0.)
		    frameflux(ifs,ib(i),n) = fluxb(nc)
		  endif
		enddo
 
! Calculate residuals.
		allocate (model(nv))
		do k = 1,nv
		    model(k) = 0.
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
		sumbp(i) = sum((primask(1:nv,ib(i))*(dvalue(1:nv,ib(i))
     *		    - model(1:nv))/sigdvalue(1:nv,ib(i)))**2)
	        npointsb(i) = sum(primask(1:nv,ib(i)))
		framechi(ifs,ib(i)) = sumbp(i)/max((npointsb(i) - 3), 1)
		deallocate(model)
	      endif
	    enddo	! ifs loop
 
	    deallocate(Af)
	    deallocate(bf)
	    deallocate(pint)
	enddo		! band loop
 
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
