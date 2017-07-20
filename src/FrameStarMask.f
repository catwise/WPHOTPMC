      subroutine FrameStarMask (ib,NSUMmax,NsrcAll, Xpos, Ypos, nnx,nny,
     1       nx, ny, nf, nsrc, nmax, MASK, Rstann, Rstwid, SNRlist, Table)

      implicit integer (i-n)
            implicit real*4 (a-h)
            implicit real*4 (o-z)

      parameter (FWHM_pix = 2.0, scale = 1.5)
c     parameter (posuncLIM = 2.0)  ! position uncertatin in arcsec      ! JWF B60711

c   real*4 Table(Nmax,25),posuncLIM
c   real*4 Table(Nmax,59),posuncLIM      ! JWF B30221
      real*4 Table(Nmax,65),posuncLIM      ! JWF B31211
      parameter (posuncLIM = 2.0)  ! position uncertainty in arcsec     ! JWF B60711
      real*4 Xpos(NSUMmax,nf,4), Ypos(NSUMmax,nf,4), SNRlist (NSUMmax,nf,4)
      real*4 scalefact (4), Rstann(4), Rstwid(4)
      integer MASK(nnx,nny, nf, 4), imask
      integer nx(4),ny(4),nsrc
      logical debug

      data scalefact/1.1,1.1,1.0,1.0/      !  scale the FWHM by this number to get the band value
c   debug = .true.
      debug = .false.

cc  initialize mask

      imask = 0

      Rmask = FWHM_pix * scalefact (ib) * scale   ! masking radius in pixels
            IRmask = nint(Rmask)

      Rl = Rstann(ib) - 1.  !  2x times the frame value
            Rh = Rl + Rstwid(ib)  !  2x times the frame value
      Rh = Rh + 1.

c initialize the mask
       do Jfr = 1, nf

       do jj=1,ny(ib)
       do ii=1,nx(ib)
        Mask(ii,jj,Jfr,ib) = imask
       enddo
       enddo

      enddo

      pscale = 2.75
      if (ib.eq.4) pscale = pscale * 2.

      pscale = 1.

      do 100 Jsrc = 1,NsrcAll

        xposunc = Table (Jsrc, 3) / pscale  !  arcsec
        yposunc = Table (Jsrc, 4) / pscale  !  arcsec

        if (xposunc.gt.posuncLIM) goto 100  ! don;t include sourcses with lousy positions
        if (yposunc.gt.posuncLIM) goto 100  ! don;t include sourcses with lousy positions

        do 101 Jfr = 1, nf

        ix = nint(Xpos (Jsrc,Jfr,ib))
        if (ix.lt.1) goto 101
              if (ix.gt.nx(ib)) goto 101

        jy = nint(Ypos (Jsrc,Jfr,ib))
        if (jy.lt.1) goto 101
              if (jy.gt.ny(ib)) goto 101

        SNR = SNRlist (Jsrc,Jfr,ib)
        if (SNR.lt.3.0) goto 101  ! avoid low SNR things


        xx = Xpos (Jsrc,Jfr,ib)
        yy = Ypos (Jsrc,Jfr,ib)

c   if (ib.eq.4) write (111,'(2f9.1,2f8.1,2x,4f8.1)') xx,yy,xposunc,yposunc,(SNRlist (Jsrc,Jfr,ib0),ib0=1,4)

c found source; mask from image

        iil = (ix - IRmask) - 1
        iil = max (iil,1)
        iih = (ix + IRmask) + 1
        iih = min (iih, nx(ib))
        jjl = (jy - IRmask) - 1
        jjl = max (1, jjl)
              jjh = (jy + IRmask) + 1
        jjh = min (jjh, ny(ib))

        n=n+1

        do   jj=jjl,jjh
              dy = (jj*1.) - yy
              dy2 = dy**2
        do   ii=iil,iih
              dx = (ii*1.) - xx
              dx2 = dx**2
              dr = sqrt ( dx2 + dy2 )

              if (dr.le.Rmask) then

            if (MASK (ii,jj,Jfr,ib).eq.0) then
              MASK (ii,jj,Jfr,ib) = Jsrc
            else
              MASK (ii,jj,Jfr,ib) = abs(MASK (ii,jj,Jfr,ib)) + Jsrc
              MASK (ii,jj,Jfr,ib) = -MASK (ii,jj,Jfr,ib)
            endif

              endif

        enddo
        enddo

 101       continue ! Jframes

 100      continue  !  nsrc


      return
      end
