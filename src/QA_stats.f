      subroutine QA_stats (nnx,nny,nx,ny,Array,pix_order,Larray,Lsize,nf,
     1     nfpix,nmax,Table, zero,wflag, nsrc,
     1     crval1,crval2,qdir,
     1     imeta, What, comment, Type,
     1     basename, coname, smode, 
     1     pointless, fram, SPIT,MIPS, fluxcon, detsnr)

!-----------------------------------------------------------------------------
! Do profile-fitting photometry on WISE images.  Results are output in the
! form of an array named "Table", arranged as follows:
!       Table(n,1)  =   RA [deg] (single precision version)
!       Table(n,2)  =   Dec [deg] (single precision version)
!       Table(n,3)  =   uncertainty in RA [arcsec]
!       Table(n,4)  =   uncertainty in Dec [arcsec]
!       Table(n,5)  =   cross term (sigxy) in RA,Dec [arcsec]
!       Table(n,6:9)    =   fluxes in bands 1,...4 [dn]
!       Table(n,10:13)  =   uncert. in fluxes in bands 1,...4 [dn]
!       Table(n,14:17)  =   reduced chi squared in bands 1,...4
!       Table(n,18) =   overall reduced chi squared
!       Table(n,19) =   blend number
!       Table(n,20) =   number of actively deblended components
!
!   where n is the source number (n = 1,...nsrc)
!
!   Note: The RA uncertainty is expressed as a great-circle angle.
!
! Double precision versions of the estimated source positions are output in
! arrays RAlist(n), DEClist(n), where n = 1,...nsrc.  In doing so, the input
! versions of those arrays are overwritten.
!-----------------------------------------------------------------------------

      implicit real*4 (a-h,o-z)
      implicit integer (i-n)

      parameter (mx = 99999, ncx=50,ncy = 50)

      character*(*) What, comment, Type,  qdir, basename(nf), fram(nf,4), coname
      character*500 soutf, poutf
      character*72  rname, pname,fout
      character*3 Sfr
      real*4 Array (nnx,nny,nfpix), Larray(Lsize), chimap (ncx,ncy), chiv (9999)
c   real(4) Table(nmax,25), metval, zero(4), zmag(mx), zv(9999)
c   real(4) Table(nmax,59), metval, zero(4), zmag(mx), zv(9999)     ! JWF B30129
      real(4) Table(nmax,65), metval, zero(4), zmag(mx), zv(9999)   ! JWF B31211
      real(4) R4Blank                                               ! JWF B60711
      integer(4) I4Blank                                            ! JWF B60711
      equivalence      (R4blank, I4blank)                           ! JWF B60711
      data I4blank/-1/                                              ! JWF B60711
c     real(4) crval1(4),crval2(4)                                   ! JWF B60714
      real(4) crval1,crval2                                         ! JWF B60714
      real*4  fluxcon(5)
      integer*2 wflag(nf,4)
      integer*4 pix_order(nf,4)
      integer igo(4), counts (mx), ipix (9999),jpix(9999)
      integer n(4), nonly (4)
      integer nx(4),ny(4)

      integer DETgo(4)

      logical smode, pointless, SPIT, MIPS, dochimap

      dochimap = .false.

c   write (6,*) 'QA Statistics'
      write (6,*) 'DETSNR = ',detsnr

c various permutations of source counts

      LQ = numchar (qdir)

      if (SPIT) goto 476

      nbad = 0

      nblend_1 = 0
      nblend_2 = 0
      nblend_3 = 0
      nblend_4 = 0
      nblend_5 = 0

c call wfits (nsubx(ib),nsuby(ib),lsize,larray,fout)


      n12 = 0
      n13 = 0
      n14 = 0
      n23 = 0
      n24 = 0
      n34 = 0
      n123 = 0
      n124 = 0
      n234 = 0
      n1234 = 0

      NdetOnly_1 = 0
      NdetOnly_2 = 0
      NdetOnly_3 = 0
      NdetOnly_4 = 0

      NonlyW1 = 0
      NonlyW2 = 0
      NonlyW3 = 0
      NonlyW4 = 0

      NonlyW12 = 0
      NonlyW13 = 0
      NonlyW14 = 0
      NonlyW23 = 0
      NonlyW24 = 0
            NonlyW34 = 0

      NonlyW123 = 0
      NonlyW124 = 0
      NonlyW134 = 0
      NonlyW234 = 0

      NallW1234 = 0

      do ib=1,4
            n(ib) = 0
            nonly (ib) = 0
      enddo

!                     Table(n,6:9)      =         fluxes in bands 1,...4 [dn]
!                     Table(n,10:13)  =         uncert. in fluxes in bands 1,...4 [dn]
!                     Table(n,18)       =         overall reduced chi squared
!                     Table(n,19)       =         blend number


      Jfr = 1
      nrchi_0_0p5 = 0
      nrchi_0p5_2 = 0
      nrchi_2_5 = 0
      nrchi_5_15 = 0
      nrchi_gt15 = 0

      rmin = 999.
            rmax = 0.
            dmin = 999.
            dmax = -999.

      do Jsrc = 1, nsrc

            ra0  = Table(Jsrc,1)
            dec0  = Table(Jsrc,2)
            call  Cel2Ec(RA0, Dec0, eLon, eLat)
            if (elon.lt.rmin) rmin = elon
                  if (elon.gt.rmax) rmax = elon
                  if (elat.lt.dmin) dmin = elat
                  if (elat.gt.dmax) dmax = elat


            rchi = Table(Jsrc,18)
            nb = Table(Jsrc,19)

            if (rchi.le.0.5) then
            nrchi_0_0p5 = nrchi_0_0p5 + 1
            else if (rchi.le.2.) then
            nrchi_0p5_2 = nrchi_0p5_2 + 1
            else if (rchi.le.5.) then
                        nrchi_2_5 = nrchi_2_5 + 1
            else if (rchi.le.15.) then
                        nrchi_5_15 = nrchi_5_15 + 1
            else if (rchi.gt.15.) then
                        nrchi_gt15 = nrchi_gt15 + 1
            endif


            if (nb.eq.1) then
            nblend_1 = nblend_1 + 1
            else if (nb.eq.2) then
                        nblend_2 = nblend_2 + 1
            else if (nb.eq.3) then
                        nblend_3 = nblend_3 + 1
            else if (nb.eq.4) then
                        nblend_4 = nblend_4 + 1
            else if (nb.ge.5) then
                        nblend_5 = nblend_5 + 1
            endif




            do ib=1,4

            igo(ib) = 0
            DETgo(ib) = 0

c         if (wflag(1,ib).eq.1) then
                  if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then


             Value = Table(Jsrc,5+ib)   !  flux
             Unc = Table(Jsrc,9+ib)   !  flux
                          if (Unc.le.0.) then

                  if (Unc.eq.0.) then
                    write (6,*) 'ERROR -- encountered zero WPRO flux uncertainty'
                    if (nbad.eq.0) write (6,*) 'WARNING  -- encountered zero WPRO flux uncertainty'
                    nbad=nbad+1
                                SNR=0.
                                                       else
                                    if (UNC.gt.-50.) then
                                      write (6,*) 'ERROR -- encountered negative WPRO flux uncertainty: ',ib,Unc
                                      nbad = nbad + 1
                                      SNR  = 0.
                                    endif
                              endif
             else
                  SNR = Value / Unc
             endif

             if (SNR.ge.3.0) igo(ib) = 1
             if (SNR.ge.detsnr) DETgo(ib) = 1


            endif

            enddo   ! ib


            do  ib0=1,4

              if (igo(ib0).gt.0) n(ib0)=n(ib0)+1

              if (igo(ib0).gt.0) then

             ncheck = 0
             do ib=1,4
               if (ib.ne.ib0) then
                  if (igo(ib).gt.0) then
                     ncheck = ncheck + 1
                  endif
               endif
            enddo

            if (ncheck.eq.0) then
                  nonly(ib0) = nonly(ib0) + 1
            endif

              endif
                           enddo


c 12
       if ((igo(1).gt.0).and.(igo(2).gt.0)) n12 = n12 + 1
       if ((igo(1).gt.0).and.(igo(3).gt.0)) n13 = n13 + 1
       if ((igo(1).gt.0).and.(igo(4).gt.0)) n14 = n14 + 1

       if ((igo(2).gt.0).and.(igo(3).gt.0)) n23 = n23 + 1
       if ((igo(2).gt.0).and.(igo(4).gt.0)) n24 = n24 + 1

       if ((igo(3).gt.0).and.(igo(4).gt.0)) n34 = n34 + 1

       if ((igo(1).gt.0).and.(igo(2).gt.0).and.(igo(3).gt.0)) n123 = n123 + 1
       if ((igo(1).gt.0).and.(igo(2).gt.0).and.(igo(4).gt.0)) n124 = n124 + 1
       if ((igo(2).gt.0).and.(igo(3).gt.0).and.(igo(4).gt.0)) n234 = n234 + 1

       if ((igo(1).gt.0).and.(igo(2).gt.0).and.(igo(3).gt.0).and.(igo(4).gt.0)) n1234 = n1234 + 1


c DETSNR stats
       nT = 0
       do ib=1,4
            if (DETgo(ib).gt.0) nT = nT +1
       enddo

       if (nT.eq.1) then
            NdetOnly_1 = NdetOnly_1 + 1
       else if (nT.eq.2) then
                        NdetOnly_2 = NdetOnly_2 + 1
       else if (nT.eq.3) then
                        NdetOnly_3 = NdetOnly_3 + 1
       else if (nT.eq.4) then
                        NdetOnly_4 = NdetOnly_4 + 1
       endif

       if (nT.eq.1) then
            if (DETgo(1).gt.0) NonlyW1 = NonlyW1 + 1
            if (DETgo(2).gt.0) NonlyW2 = NonlyW2 + 1
            if (DETgo(3).gt.0) NonlyW3 = NonlyW3 + 1
            if (DETgo(4).gt.0) NonlyW4 = NonlyW4 + 1
       endif

      if (nT.eq.2) then
                  if ((DETgo(1).gt.0).and.(DETgo(2).gt.0)) NonlyW12 = NonlyW12 + 1 
            if ((DETgo(1).gt.0).and.(DETgo(3).gt.0)) NonlyW13 = NonlyW13 + 1
            if ((DETgo(1).gt.0).and.(DETgo(4).gt.0)) NonlyW14 = NonlyW14 + 1
            if ((DETgo(2).gt.0).and.(DETgo(3).gt.0)) NonlyW23 = NonlyW23 + 1
                        if ((DETgo(2).gt.0).and.(DETgo(4).gt.0)) NonlyW24 = NonlyW24 + 1
                        if ((DETgo(3).gt.0).and.(DETgo(4).gt.0)) NonlyW34 = NonlyW34 + 1
      endif

      if (nT.eq.3) then
                        if ((DETgo(1).gt.0).and.(DETgo(2).gt.0).and.(DETgo(3).gt.0)) NonlyW123 = NonlyW123 + 1
            if ((DETgo(1).gt.0).and.(DETgo(2).gt.0).and.(DETgo(4).gt.0)) NonlyW124 = NonlyW124 + 1
            if ((DETgo(1).gt.0).and.(DETgo(3).gt.0).and.(DETgo(4).gt.0)) NonlyW134 = NonlyW134 + 1
            if ((DETgo(2).gt.0).and.(DETgo(3).gt.0).and.(DETgo(4).gt.0)) NonlyW234 = NonlyW234 + 1
      endif

      if (nT.eq.4) then
            NallW1234 = NallW1234 + 1
      endif
                           enddo ! Jsrc


      What = "NPRO_nbad"
            TYPE = "i"
            metval = nbad*1.
            comment = 'number of bad WPRO (UNC <= 0.0) values'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)


      do ib=1,4

      What = "n3sig"
            TYPE = "i"
            metval = n(ib)*1.
            comment = 'number of WPRO sources with SNR > 3 '
            iband = ib
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      enddo


      do ib=1,4

            What = "n3sig_only"
            TYPE = "i"
            metval = nonly(ib)*1.
            comment = 'number of solo WPRO sources with SNR > 3 '
            iband = ib
            call MetWr (imeta, iband, What, TYPE, metval, comment)

            enddo



      What = "n12"
            TYPE = "i"
            metval = n12*1.
            comment = 'number of W1 and W2 WPRO sources with SNR > 3 '
            iband = 12
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "n13"
            TYPE = "i"
            metval = n13*1.
            comment = 'number of W1 and W3 WPRO sources with SNR > 3 '
            iband = 13
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "n14"
            TYPE = "i"
            metval = n14*1.
            comment = 'number of W1 and W4 WPRO sources with SNR > 3 '
            iband = 14
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)


      What = "n23"
            TYPE = "i"
            metval = n23*1.
            comment = 'number of W2 and W3 WPRO sources with SNR > 3 '
            iband = 23
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "n24"
            TYPE = "i"
            metval = n24*1.
            comment = 'number of W2 and W4 WPRO sources with SNR > 3 '
            iband = 24
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)


      What = "n34"
            TYPE = "i"
            metval = n34*1.
            comment = 'number of W3 and W4 WPRO sources with SNR > 3 '
            iband = 34
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "n123"
            TYPE = "i"
            metval = n123*1.
            comment = 'number of W1, W2 and W3 WPRO sources with SNR > 3 '
            iband = 123
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "n124"
            TYPE = "i"
            metval = n124*1.
            comment = 'number of W1, W2 and W4 WPRO sources with SNR > 3 '
            iband = 124
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)


      What = "n234"
            TYPE = "i"
            metval = n234*1.
            comment = 'number of W2, W3 and W4 WPRO sources with SNR > 3 '
            iband = 234
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)


      What = "n1234"
            TYPE = "i"
            metval = n1234*1.
            comment = 'number of W1, W2, W3 and W4 WPRO sources with SNR > 3 '
            iband = 1234
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

cccc

      What = "detsnr"
      TYPE = "r"
      metval = detsnr
      comment = 'Source detection threshold SNR'
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nsrc"
      TYPE = "i"
      metval = nsrc * 1.
      comment = 'total number of sources'
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NdetOnly_1"
            TYPE = "i"
            metval = NdetOnly_1 * 1.
            comment = 'total number of sources detected in only one band, SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)
              What = "NdetOnly_2"
            TYPE = "i"
            metval = NdetOnly_2 * 1.
            comment = 'total number of sources detected in only two bands, SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NdetOnly_3"
            TYPE = "i"
            metval = NdetOnly_3 * 1.
            comment = 'total number of sources detected in only three bands, SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NdetOnly_4"
            TYPE = "i"
            metval = NdetOnly_4 * 1.
            comment = 'total number of sources detected in four bands, SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW1"
      TYPE = "i"
            metval = NonlyW1 * 1.
            comment = 'total number of sources detected in only W1,  SNR > DETSNR'
            iband = 1
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW2"
            TYPE = "i"
            metval = NonlyW2 * 1.
            comment = 'total number of sources detected in only W2,  SNR > DETSNR'
            iband = 2
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW3"
            TYPE = "i"
            metval = NonlyW3 * 1.
            comment = 'total number of sources detected in only W3,  SNR > DETSNR'
            iband = 3
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW4"
            TYPE = "i"
            metval = NonlyW4 * 1.
            comment = 'total number of sources detected in only W4,  SNR > DETSNR'
            iband = 4
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW12"
            TYPE = "i"
            metval = NonlyW12 * 1.
            comment = 'total number of sources detected in only W1 & W2,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW13"
            TYPE = "i"
            metval = NonlyW13 * 1.
            comment = 'total number of sources detected in only W1 & W3,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW14"
            TYPE = "i"
            metval = NonlyW14 * 1.
            comment = 'total number of sources detected in only W1 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW23"
            TYPE = "i"
            metval = NonlyW23 * 1.
            comment = 'total number of sources detected in only W2 & W3,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW24"
            TYPE = "i"
            metval = NonlyW24 * 1.
            comment = 'total number of sources detected in only W2 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW34"
            TYPE = "i"
            metval = NonlyW34 * 1.
            comment = 'total number of sources detected in only W3 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW123"
            TYPE = "i"
            metval = NonlyW123 * 1.
            comment = 'total number of sources detected in only W1, W2 & W3,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW124"
            TYPE = "i"
            metval = NonlyW124 * 1.
            comment = 'total number of sources detected in only W1, W2 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW134"
            TYPE = "i"
            metval = NonlyW134 * 1.
            comment = 'total number of sources detected in only W1, W3 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "NonlyW234"
            TYPE = "i"
            metval = NonlyW234 * 1.
            comment = 'total number of sources detected in only W2, W3 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

c   What = "NonlyW1234"
      What = "NallW1234"            ! JWF B21113
            TYPE = "i"
c         metval = NonlyW1234 * 1.
            metval = NallW1234 * 1.   ! JWF B21113
            comment = 'total number of sources detected in W1, W2, W3 & W4,  SNR > DETSNR'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)



      What = "nrchi_0_05"
      TYPE = "i"
            metval = nrchi_0_0p5*1.
            comment = 'number of WPRO sources with 0 < rChi^2 <= 0.5'
      iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nrchi_05_2"
            TYPE = "i"
            metval = nrchi_0p5_2*1.
            comment = 'number of WPRO sources with 0.5 < rChi^2 <= 2.0'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nrchi_2_5"
            TYPE = "i"
            metval = nrchi_2_5*1.
            comment = 'number of WPRO sources with 2 < rChi^2 <= 5.0'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nrchi_5_15"
            TYPE = "i"
            metval = nrchi_5_15*1.
            comment = 'number of WPRO sources with 5 < rChi^2 <= 15.0'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

       What = "nrchi_15"
            TYPE = "i"
            metval = nrchi_gt15*1.
            comment = 'number of WPRO sources with rChi^2 > 15.0'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nblend_1"
            TYPE = "i"
            metval = nblend_1*1.
      comment = 'number of WPRO sources with nblend=1'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nblend_2"
            TYPE = "i"
            metval = nblend_2*2.
            comment = 'number of WPRO sources with nblend=2'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nblend_3"
            TYPE = "i"
            metval = nblend_3*3.
            comment = 'number of WPRO sources with nblend=3'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nblend_4"
            TYPE = "i"
            metval = nblend_4*4.
            comment = 'number of WPRO sources with nblend=4'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)

      What = "nblend_5"
            TYPE = "i"
            metval = nblend_5*5.
            comment = 'number of WPRO sources with nblend>=5'
            iband = 0
            call MetWr (imeta, iband, What, TYPE, metval, comment)


c### Chi^2 map

      rname = coname
            L = numchar (rname)

! strip the path off the name
            idex = 0
            do k=L,1,-1
              if (rname(k:k).eq.'/') then
                        idex = k
                        goto 950
              endif
            enddo

 950      if (idex.gt.0) then
                        rname = coname(idex+1:L)
                        L = numchar (rname)
            endif

c         write (6,'(a)') rname(1:L)

      rmid = (rmax+rmin)/2.
            dmid = (dmax+dmin)/2.

c            write (6,*) rmin,rmax
c            write (6,*) dmin,dmax
c            write (6,*) rmid,dmid

            dhalf = (dmax+dmin) / 2.
            dd=dhalf/ 57.2957795

      xscale = (ncx*1.) / ((rmax-rmin)*cos(dd))
            yscale = (ncy*1.) / (dmax-dmin)

            cdelt1 = cos(dd) * (rmax-rmin) / (ncx*1.)
            cdelt2 = (dmax-dmin) / (ncy*1.)

      do IB = 1, 4

      nv = 0
      do Jsrc = 1, nsrc

c       if (wflag(1,ib).eq.1) then
        if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then


                         Value = Table(Jsrc,5+ib)   !  flux
                         Unc = Table(Jsrc,9+ib)   !  flux
             zchi = Table(Jsrc,ib+13)
!                     Table(n,14:17)  =         reduced chi squared in bands 1,...4

             SNR = 0.
             if ((Unc.gt.0.).and.(zchi.gt.0.).and.(zchi.lt.999.)) then
              SNR = Value / Unc
             endif

            if ((SNR.gt.10.).and.(SNR.lt.500.).and.(nv.lt.9999)) then

              nv=nv+1
              ra0 = Table(Jsrc,1)
              dec0 = Table(Jsrc,2)
              call  Cel2Ec(RA0, Dec0, eLon, eLat)
              dd = elat/ 57.2957795
                          dra = (elon - rmid) * cos(dd)
                          xpix = (ncx/2.) + (xscale * dra)
                          ipix (nv) = nint(xpix)
                          ddec = elat - dmid
                          ypix = (ncy/2.) + (yscale * ddec)
                          jpix (nv) = nint(ypix)
              chiv(nv) = zchi

            endif

         endif
      enddo

c   write (6,*) nv

      if (nv.gt.9.) then

       do j=1,ncy
             do i=1,ncx
                        chimap(i,j) = -9.

              mv=0
              do K = 1,nv
                        if ( (ipix(K).eq.i).and.(jpix(K).eq.j).and.(mv.lt.9999)) then
                          mv=mv+1
                          zv(mv) = chiv(K)
                        endif
              enddo

              if (mv.gt.2) then
                        call MDIAN1(zv,mv,XMED)
              else if (mv.gt.0) then
                        call MOMENT(zv,mv,AVE,SDEV)
                        xmed = ave
              else
                        xmed = -9.
              endif

              chimap(i,j) = xmed

              enddo
              enddo

       mv=0

             ic = 0
             do j= 1,ncy
             do i=1,ncx
              ic=ic+1
              larray(ic) = chimap(i,j)

              if ((chimap(i,j).gt.0.).and.(mv.lt.9999)) then
                        mv=mv+1
                        zv(mv) = chimap(i,j)
              endif
             enddo
             enddo


             if ((nf.gt.1).and.(dochimap)) then
        if (ib.eq.1) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w1-mdexqachi.fits'
              else if (ib.eq.2) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w2-mdexqachi.fits'
              else if (ib.eq.3) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w3-mdexqachi.fits'
              else if (ib.eq.4) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w4-mdexqachi.fits'
              endif 

c   write (6,'(a)') soutf

       call wfits_chimap (ncx,ncy,lsize,larray,soutf, rmid,dmid,cdelt1,cdelt2)

      endif


             call MOMENT(zv,mv,AVE,SDEV)
       xmed = ave
             if (mv.gt.2) then
                        call MDIAN1(zv,mv,XMED)
             endif

       if (ib.eq.1) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w1-mdexqachi.tbl'
              else if (ib.eq.2) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w2-mdexqachi.tbl'
              else if (ib.eq.3) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w3-mdexqachi.tbl'
              else if (ib.eq.4) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w4-mdexqachi.tbl'
              endif

        open (unit=32,file=soutf)
        write (32,'(a)') '\ global chi^2 statistics'
        write (32,'(a)') '| N   | median |  ave   |  sdev  |'
        write (32,'(a)') '| i   |   r      |   r      |      r   |'
              write (32,'(a)') '|       |            |            |            |'
        write (32,'(a)') '|       | null   |   null |  null  |'

c       if ((mv.le.0).or.(median.lt.-5.)) then
        if ((mv.le.0).or.(Xmed.lt.-5.)) then       !  JWF B21114
            write (32,'(i6,a)')  mv,'  null        null       null'
        else
            write (32,'(i6,3f9.3)') mv,Xmed,Ave,sdev
        endif
        close (32)
              write (6,*) 'global chi2 stats (n,med,ave,sig): ',ib,mv,Xmed,Ave,sdev


      endif

      enddo  ! IB


ccccccccccccccccccccccccccccccccccccc
c   now compute the starcounts
c
c zmag = zero(ib) - (2.5*log10(xint))
c                     Table(n,6:9)      =         fluxes in bands 1,...4 [dn]


             LQ = numchar (qdir)

      Area = 1.0  ! assume one sq. degree

      do 750 ib=1,4

c              if (wflag(1,ib).eq.0) goto 750
       if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then


c       ra0 = crval1(ib)                         ! JWF B60714: CatWISE
c       dec0 = crval2(ib)                        !             doesn't need
        ra0 = crval1                             !             band-dependent     
        dec0 = crval2                            !             frame centers
        call ast_precess (ra0, dec0, 2000., raout, decout, 1950.)
              call asteqtogal (raout,decout,glon,glat)


        if (ib.eq.1) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w1-mdexqasc.tbl'
              else if (ib.eq.2) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w2-mdexqasc.tbl'
              else if (ib.eq.3) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w3-mdexqasc.tbl'
              else if (ib.eq.4) then
                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w4-mdexqasc.tbl'
              endif


c   write (6,*) 'this is diff counts'
      write (6,'(a)') soutf(1:72)
        open (unit=66,file=soutf)
        write (66,'(a)') '\ Differential source counts in Log (counts / mag / area)'
        write (66,'(a)') '\ Area is assumed to be 1 sq. degree'
        write (66,'(a,1x,i3)') '\ band = ',ib
        write (66,'(a,1x,2f11.5)') '\ image center galactic coordinates = ',glon,glat
        write (66,'(a)') '|  bin |      N |   scount|      slow |      shigh|'
        write (66,'(a)') '|   r  |      i |       r   |       r   |        r  |'
        write (66,'(a)') '|        |        |             |             |             |'
        write (66,'(a)') '|        | null |  null   |   null  |      null |'



        ns = 0
              do Jsrc = 1, nsrc
            flux = Table(Jsrc,5+ib)
            if (flux.gt.0.) then
             ns=ns+1
             if (ns.gt.mx) ns = mx
             zmag(ns) = zero(ib) - (2.5*log10(flux))
            endif
        enddo


        blow = 0.0
        bhigh = 25.
        bwidth = 0.5

        call bincounts (blow,bhigh,bwidth,zmag,Ns,mx,counts)


c       sum = 0.

        do 751 K=1,mx
            zbin = blow + (K * bwidth) - (bwidth / 2.0)
            if (zbin.lt.blow) goto 751
            if (zbin.gt.bhigh) goto 751

            star = 0.
            del_high = 0.
            del_low = 0.

            if (counts(K).gt.0) then
              T = counts(K)*1.
              zcount = T / (bwidth * Area)
              star = log10 (zcount)

              scount = sqrt ( T )

              shigh = (T + scount) / (bwidth * Area)
              star_high = log10 (shigh)
              del_high = star_high - star

              slow = (T - scount) / (bwidth * Area)
              slow = max (0.1, slow)
              star_low = log10 (slow)
              del_low = star - star_low


            endif

c   sum = sum + counts(K)
            write (66,'(f7.3,i7,3f10.4)') zbin,counts(K),star,del_low,del_high

 751        continue

       close (66)

             endif ! wflag conditional

 750      continue  ! ib



cccccccccccccccccccccccccc  Pointless Images

 476      if (.not.pointless) goto 477

      write (6,*) 'generating pointless image'

      do Jfr = 1,nf

              if (Jfr.lt.10) then
                        write (Sfr,'(a,i1)') '00',Jfr
              else if (Jfr.lt.100) then
                        write (Sfr,'(a,i2)') '0',Jfr
              else if (Jfr.lt.1000) then
                        write (Sfr,'(i3)') Jfr
              endif

              rname = basename(Jfr)
              L = numchar(rname)

              do ib=1,4

                  if (ib.eq.1) then
                         poutf = qdir(1:LQ) // '/' // rname(1:L) // '-w1-pointless.fits'
                  else if (ib.eq.2) then
             poutf = qdir(1:LQ) // '/' // rname(1:L) // '-w2-pointless.fits'
                  else if (ib.eq.3) then
                         poutf = qdir(1:LQ) // '/' // rname(1:L) // '-w3-pointless.fits'
                  else if (ib.eq.4) then
                         poutf = qdir(1:LQ) // '/' // rname(1:L) // '-w4-pointless.fits'
                  endif

c                  if (wflag(1,ib).eq.1) then
            if (wflag (jfr,ib) == 1) then


c   write (6,*) nx(ib), ny(ib), Jfr, ib

              ic=0
                    do J=1,ny(ib)
                    do I=1,nx(ib)
            ic=ic+1
                        Larray(ic) = Array (I,J,pix_order(Jfr,ib))
            if (SPIT) then
                  if ((MIPS).and.(ib.eq.4)) then
                    Larray(ic) = Larray(ic) * fluxcon (5)
                  else
                    Larray(ic) = Larray(ic) * fluxcon (ib)
                  endif
            endif
c         if (Larray(ic).lt.-999.) Larray(ic) = 1. / 0.  ! create NaN      ! JWF B60711
            if (Larray(ic).lt.-999.) Larray(ic) = R4Blank  ! create NaN      ! JWF B60711
              enddo
              enddo

              write (6,'(a)') poutf(1:200)

              call  wimage (nx(ib),ny(ib),lsize,larray,fram(Jfr,ib),poutf)

            endif ! wflag

               enddo  ! ib

            enddo  ! Jfr



cccccccccccccccccccccccccccccccc
c  residual image histogram  --- ONLY do this for single band mode

 477      if ((.not.smode).or.(SPIT) ) return

      blow = -999.
            bhigh = 7999.
            bwidth = 1.

      do Jfr = 1,nf
        if (Jfr.lt.10) then
            write (Sfr,'(a,i1)') '00',Jfr
        else if (Jfr.lt.100) then
                        write (Sfr,'(a,i2)') '0',Jfr
        else if (Jfr.lt.1000) then
                        write (Sfr,'(i3)') Jfr
        endif

       rname = basename(Jfr)
       L = numchar(rname)
                     do ib=1,4

        if (ib.eq.1) then
c                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '.' // Sfr // '.res-W1.tbl'
             soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w1-mdexqares.tbl'
              else if (ib.eq.2) then
c                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '.' // Sfr //  '.res-W2.tbl'
            soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w2-mdexqares.tbl'
              else if (ib.eq.3) then
c                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '.' // Sfr //  '.res-W3.tbl'
            soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w3-mdexqares.tbl'
              else if (ib.eq.4) then
c                        soutf = qdir(1:LQ) // '/' // rname(1:L) // '.' // Sfr //  '.res-W4.tbl'
            soutf = qdir(1:LQ) // '/' // rname(1:L) // '-w4-mdexqares.tbl'
              endif

c   write (6,*) 'qa test'
c   write (6,'(a)') soutf(1:72)
c       if (wflag(1,ib).eq.1) then
        if ( any ( wflag ( 1:nf,ib ) ==1 ) ) then


         open (unit=66,file=soutf)
               write (66,'(a)') '\ Histogram of residual (WPRO) image'
         write (66,'(a,1x,i4)') '\ frame index ',Jfr
         write (66,'(a,1x,i3)') '\ band = ',ib
         write (66,'(a,1x,f9.2)') '\ bin width = ',bwidth
               write (66,'(a)') '| bin      |       N |       LogN|       low |       high|'
               write (66,'(a)') '|  r       |       i |       r   |       r   |        r  |'

         nh = 0
         do J=1,ny(ib)
         do I=1,nx(ib)
            val = Array (I,J,pix_order(Jfr,ib))
            if (val.gt.-500.) then
              nh=nh+1
              Larray(nh) = val
            endif
         enddo
         enddo
                     do K=1,mx
            counts(K) = 0
         enddo

         if (nh.gt.0) then
                  call bincounts (blow,bhigh,bwidth,Larray,Nh,mx,counts)
         endif

         on = 0

         do 761 K=1,mx
                        zbin = blow + (K * bwidth) - (bwidth / 2.0)
                        if (zbin.lt.blow) goto 761
                        if (zbin.gt.bhigh) goto 761

                        star = 0.
                        del_high = 0.
                        del_low = 0.

                        if (counts(K).gt.0) then
                          T = counts(K)*1.
                          zcount = T / (bwidth * Area)
                          star = log10 (zcount)

                          scount = sqrt ( T )

                          shigh = (T + scount) / (bwidth * Area)
                          star_high = log10 (shigh)
                          del_high = star_high - star

                          slow = (T - scount) / (bwidth * Area)
                          slow = max (0.1, slow)
                          star_low = log10 (slow)
                          del_low = star - star_low

              on = 1
                        endif

                        if (counts(K).gt.0) then
                  write (66,'(f9.1,i8,3f10.4)') zbin,counts(K),star,del_low,del_high
            endif

 761        continue

             close (66)

        endif

       enddo  ! ib

      enddo  ! Jfr



      return
      end
