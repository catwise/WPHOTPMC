      subroutine act_wise(adb_nmax, adb_dchi, adb_chimax, adb_pksep,
     *          adb_minsep, apflux, psfs, npsize, npnmax,varpsf,nlayers,numsrc,
     *          onframe,nmax,nf,mxbands,wflag,nfi,p, p_final, nactive,
     *      iter, nsteps, nDoF, dbg) ! JWF B30117

! ACTive deblending for WISE.
!
! Input parameters:
!      adb_nmax      = maximum number of additional components
!                    (set to zero to turn off active deblending)
!      adb_dchi       = minimum acceptable reduction in the reduced chi
!                    squared due to active deblending
!      adb_chimax      = maximum acceptable value of reduced chi squared
!                    after deblending
!      adb_pksep      = minimum separation for which detector can distinguish
!                    sources ["]
!      adb_minsep      = lower limit for separation of deblended sources ["]
!      apflux            = set of aperture fluxes
!      p            = initial set of position parameters
!
! Output parameters:
!      p_final            = final set of position parameters
!      nactive            = number of new components added by active deblending

        implicit real*4 (a-h,o-z)
      implicit integer (i-n)

      include 'wpro.inc'
      integer, parameter :: maxblend = 100
      integer, parameter :: maxblend2 = 200
      integer, parameter :: maxbands = 4
      real, parameter :: dtor = 0.0174533

      real(8) racan(maxblend),deccan(maxblend)
      real(4) p(*),flux(maxbands,maxblend),snr(maxbands)
      real(4) dvalue(maxpix,maxbands), sigdvalue(maxpix,maxbands)
      real(4) apflux(maxbands,maxblend), rchisqb(maxbands)
      real(4) psfs(npsize,npsize,npnmax,4)
      real(4) pxc(maxbands),pyc(maxbands),psamp(maxbands)
      real(4) xpixscale(maxframes),ypixscale(maxframes),tframe(maxframes)
      real(4) p_init(maxblend2),p_final(maxblend2)
      real(4) pinitial(maxblend2),pcurr(maxblend2),pnew(maxblend2)
      real(4) flux_init(maxbands,maxblend),rchisqb_init(maxbands)
       real(4) fluxinitial(maxbands,maxblend)
      real(4) x(maxblend),y(maxblend)
      real(4) varpsf(maxpix,maxbands),framesig(maxframes,maxbands,maxblend)
      real(4) sigx(maxblend),sigy(maxblend),sigxy(maxblend)
      real(4) sigflux(maxbands,maxblend)
      integer nvalues(maxbands),nnpsf(maxbands)
      integer ivalue(maxpix,maxbands), jvalue(maxpix,maxbands)
      integer WCS(maxframesin,maxbands), offscl
      integer ifq(maxframes),ibq(maxframes),kpn(maxframes)
      integer ifvalue(maxpix,maxbands),primask(maxpix,maxbands)
      integer adb_nmax
c      integer(2) bset(maxframes),nlayers(maxbands),wflag(nfi,4) ! JWF B30207
      integer(2) nlayers(maxbands),wflag(nfi,4)                 ! JWF B30207
      logical(4) goodcomp(maxblend),goodcomp_init(maxblend),ok
      logical(1) onframe(nmax,nf,mxbands),finished, outofrange
        integer*4 iter, nsteps, nDoF                              ! JWF B30207
        logical*4 dbg  ! JWF B30315
        real*4    dRA  ! JWF B30327

      common /funccom/ racan,deccan,nvalues,ivalue,jvalue,dvalue,
     *          sigdvalue,ifvalue,WCS,ifq,ibq,primask,nframes,
     *          kpn,nnpsf,pxc,pyc,psamp,nblend,goodcomp,icp,
     *          xpixscale,ypixscale,flux,rchisq,rchisqb,tframe

! Check that arrays are large enough.
      if (nblend+adb_nmax > maxblend) then
          print *,'Cannot do active deblending: ',
     *            'Maximum number of components exceeded'
          nactive = 0
          return
      endif

! Save the initial solution.
      mblend = ntrue(goodcomp,nblend)
      nc = 0
      do n = 1,nblend
        if (goodcomp(n)) then
          nc = nc + 1
          p_init(nc) = p(nc)
          p_init(mblend+nc) = p(mblend+nc)
          do ib = 1,maxbands
            flux_init(ib,n) = flux(ib,n)
          enddo
        endif
      enddo
      do ib = 1,maxbands
          rchisqb_init(ib) = rchisqb(ib)
      enddo
      rchisq_init = rchisq
      nblend_init = nblend
      goodcomp_init = goodcomp

! Add additional component and do flux-only solution at a grid of component
! positions.
      do iadd = 1,adb_nmax
          racan(nblend+iadd) = racan(1)
          deccan(nblend+iadd) = deccan(1)
      enddo
      pcurr = p_init
      phi = func(pcurr,psfs,npsize,npnmax,nDoF,dbg)
      
      ds = adb_minsep                  ! step size for coarse scan [arcsec]
      iscan = nint(2.*adb_pksep/ds)      ! halfwidth of coarse scan [samples]
      nblendmax = nblend + adb_nmax
      finished = .false.

      do while (nblend < nblendmax .and. .not.finished)
          pinitial = pcurr
          fluxinitial = flux
          rchisqinitial = rchisq
          mblend = ntrue(goodcomp,nblend)
          nc = 0
          do n = 1,nblend
            if (goodcomp(n)) then
            nc = nc + 1
            x(n) = pcurr(nc)
            y(n) = pcurr(mblend+nc)
            endif
          enddo

          nblend = nblend + 1
          goodcomp(nblend) = .true.
          nc = 0
          do n = 1,nblend-1
            if (goodcomp(n)) then
            nc = nc + 1
            pcurr(nc) = x(n)
            pcurr(mblend+nc) = y(n)
            endif
          enddo
          rchisqmin = 1.e35

          do ioff = -iscan,iscan 
          do joff = -iscan,iscan
            irsq = ioff**2 + joff**2
            if (irsq /= 0 .and. irsq <= iscan**2) then
            pcurr(nc+1) = x(1) + ioff*ds/3600.
            pcurr(mblend+nc+1) = y(1) + joff*ds/3600.
            phi = func(pcurr,psfs,npsize,npnmax,nDoF,dbg)
            if (rchisq < rchisqmin) then
                pnew = pcurr
                rchisqmin = rchisq
            endif
            endif
          enddo
          enddo

          pcurr = pnew
          mblend = ntrue(goodcomp,nblend)
          np = 2*mblend

! Calculate full maximum likelihood solution based on the model with the
! extra component.
          ftol = 1.e-3
          phimin = xpixscale(1)/8.
c           call most_wise(pcurr,np,psfs,npsize,npnmax,ftol,iter,phimin)        ! JWF B30207
          call most_wise(pcurr,np,psfs,npsize,npnmax,ftol,iter,phimin,nsteps,dbg) ! JWF B30207

! Check that new component is in range (not closer to any other component than
! adb_minsep and not further away from primary than adb_pksep).
          outofrange = .false.
          nc = 0
          do n = 1,nblend
            if (goodcomp(n)) then
            nc = nc + 1
            x(n) = pcurr(nc)
            y(n) = pcurr(mblend+nc)
            endif
          enddo
          distmin = 1.e35
          cd = cos(deccan(1)*dtor)
          nc = 0
          do n = 1,nblend-1
            if (goodcomp(n)) then
            nc = nc + 1
                dRA = abs(racan(1)+x(nblend)-(racan(n)+x(n)))
                if (dRA .ge. 180.0) dRA = dRA - sign(360.0,dRA)
            dist = sqrt((dRA*cd)**2 + (deccan(1)+y(nblend)-(deccan(n)+y(n)))**2)
            if (dist < distmin) distmin = dist
            endif
          enddo
          if (distmin < adb_minsep/3600.)  outofrange = .true.
          pridist = sqrt((x(nblend)-x(1))**2 + (y(nblend)-y(1))**2)
          if (pridist > adb_pksep/3600.)  outofrange = .true.

          if (rchisq>(rchisqinitial-adb_dchi) .or. outofrange) then
            pcurr = pinitial
            flux = fluxinitial
            rchisq = rchisqinitial
            goodcomp(nblend) = .false.
            nblend = nblend - 1
            finished = .true.  
          endif
      enddo

! Check that the primary and extra component didn't get switched.
      if (goodcomp(1)) then
          if(pcurr(1)**2+pcurr(1+mblend)**2 > 
     *            pcurr(mblend)**2+pcurr(2*mblend)**2) then
                xadb = pcurr(mblend)
                yadb = pcurr(2*mblend)
                pcurr(mblend) = pcurr(1)
                pcurr(2*mblend) = pcurr(1+mblend)
                pcurr(1) = xadb
                pcurr(1+mblend) = yadb
          endif
      endif
            
      phi = func(pcurr,psfs,npsize,npnmax,nDoF,dbg)
      rchisqfinal = rchisq
      nactive = nblend - (nblendmax - adb_nmax)
      p_final = pcurr
      call unspew(p_final,psfs,npsize,npnmax,varpsf,nlayers,numsrc,onframe,
     *     nmax,nf,maxbands,wflag,nfi,framesig,sigx,sigy,sigxy,sigflux,ok,dbg)

! Check the quality of the solution.
      nul = 0
      nzs = 0
      nls = 0
      where(sigflux <= 0.) sigflux = 1.e10
      do n = 1,nblend
          if (goodcomp(n)) then
            snr = flux(1:maxbands,n)/sigflux(1:maxbands,n)
            if (.not.any(snr > 2.)) nul = nul + 1
            if (sigx(n)==0. .or. sigy(n)==0. .or. sigxy(n)==0.) nzs = nzs+1
            if (sigx(n) > adb_pksep .or. sigy(n) > adb_pksep) nls = nls+1
          endif
      enddo

! Reject the solution if the reduced chi squared exceeds the threshold, or
! if it failed the quality checks.
      if (rchisqfinal > adb_chimax .or. nul>0 .or. nzs>0 .or. nls>0) then
          p_final = p_init
          flux = flux_init
          rchisqb = rchisqb_init
          rchisq = rchisq_init
          nblend = nblend_init
          goodcomp = goodcomp_init
          nactive = 0
          return
      endif

      return

      end
