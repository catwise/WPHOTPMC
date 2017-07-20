	subroutine TripNearby (Jsrc0,x0,y0,ib,J,nnx,nny,nx,ny,nmax,nf, nfi, NSUMmax, NsrcAll,
     1       Xpos, Ypos, nsrc, MASK, Rstann, Rstwid, SNRlist)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

c	parameter (FWHM_pix = 2.0, scale = 2.75)   
	parameter (FWHM_pix = 2.0, scale = 1.5)

	real*4 Xpos(NSUMmax,nf,4), Ypos(NSUMmax,nf,4), SNRlist (NSUMmax,nf,4)
	real*4 scalefact (4), Rstann(4), Rstwid(4)
	integer*2 MASK(nnx,nny,4), imask
	integer nx(4),ny(4),nsrc
	logical debug

	data scalefact/1.1,1.1,1.0,1.0/    !  scale the FWHM by this number to get the band value
c	debug = .true.
	debug = .false.

cc  initialize mask

	imask = -947


	 Rmask = FWHM_pix * scalefact (ib) * scale   ! masking radius in pixels
c        Rmask = FWHM_pix *  scale   ! masking radius in pixels
         IRmask = nint(Rmask)

	Rl = Rstann(ib) - 1.
        Rh = Rl + Rstwid(ib)
	Rh = Rh + 1.

	i0 = nint(x0)
        j0 = nint(y0)

c	write (6,*) ib,i0,j0,rh

	il = int (x0 - Rh)
	ih = int (x0 + Rh)
	jl = int (y0 - Rh)
        jh = int (y0 + Rh)
	
	il = max (il,1)
	ih = min (ih,nx(ib))
	jl = max (jl,1)
        jh = min (jh,ny(ib))

c	write (6,*) il,ih,jl,jh, nsrc

c initialize the mask
	 do jj=jl,jh
	 do ii=il,ih
	  Mask(ii,jj,ib) = 1
	 enddo
	 enddo

c mask stars; avoid target star

	do 100 Jsrc = 1,NsrcAll

	  if (Jsrc.eq.Jsrc0) goto 100

	  ix = nint(Xpos (Jsrc,J,ib))
	  if (ix.lt.il) goto 100
          if (ix.gt.ih) goto 100

	  jy = nint(Ypos (Jsrc,J,ib))
	  if (jy.lt.jl) goto 100
          if (jy.gt.jh) goto 100

	  idelx = abs(ix - i0)
	  jdely = abs(jy - j0)

	  if ( (idelx.le.2).and.(jdely.le.2) ) goto 100

	  SNR = SNRlist (Jsrc,J,ib)

	  if (SNR.lt.3.0) goto 100  ! avoid low SNR things


c found source; mask from image

c	write (6,*) ix,jy

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
	      dy = (jj*1.) - Ypos (Jsrc,J,ib)
	      dy2 = dy**2
	  do   ii=iil,iih
	      dx = (ii*1.) - Xpos (Jsrc,J,ib)
	      dx2 = dx**2
	      dr = sqrt ( dx2 + dy2 )

	      if (dr.le.Rmask) then

		MASK (ii,jj,ib) = iMask

	      endif

	  enddo
	  enddo

 100	continue  !  nsrc


	return
	end
