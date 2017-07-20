	subroutine CoaddNearby (Jsrc0,x0,y0,ib,NSUMmax,NsrcAll,
     1     nx, ny, WCSvb,
     1     nf,nsrc, MASK, Rstann, Rstwid, SNRlist, RAlist, Declist)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (FWHM_pix = 4.0, scale = 1.5)

	real*8 ra0,dec0, RAlist(NSUMmax), Declist(NSUMmax)
	real*4 SNRlist (NSUMmax,nf,4)
	real*4 scalefact (4), Rstann(4), Rstwid(4)
	integer*2 MASK(5401,5401,4), imask
	integer nx(4),ny(4),nsrc, WCSvb (4)
	logical debug

	data scalefact/1.1,1.1,1.0,1.0/    !  scale the FWHM by this number to get the band value
c	debug = .true.
	debug = .false.

cc  initialize mask

	imask = -947


	 Rmask = FWHM_pix * scalefact (ib) * scale   ! masking radius in pixels
         IRmask = nint(Rmask)

	Rl = (Rstann(ib)*2.) - 1.  !  2x times the frame value
        Rh = Rl + (Rstwid(ib)*2.)  !  2x times the frame value
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
	ppp = 1.375
	edgebuf = 0.

	do 100 Jsrc = 1,NsrcAll

	  if (Jsrc.eq.Jsrc0) goto 100

	  Ra0 = RAlist(Jsrc) * 1.d0
          Dec0 = Declist(Jsrc) * 1.d0
	  igo = 1
	  call WhatPos ( 1,WCSvb(ib) ,nx ( ib ) ,ny ( ib ) , ppp, '0', ra0, dec0, xx, yy, igo, edgebuf )

c	 write (69,*) ra0,dec0,xx,yy

	  ix = nint(xx)
	  if (ix.lt.il) goto 100
	  if (ix.gt.ih) goto 100

          jy = nint(yy)
	  if (jy.lt.jl) goto 100
          if (jy.gt.jh) goto 100


	  idelx = abs(ix - i0)
	  jdely = abs(jy - j0)

	  if ( (idelx.le.2).and.(jdely.le.2) ) goto 100

	  SNR = SNRlist (Jsrc,1,ib)

	  if (SNR.lt.3.0) goto 100  ! avoid low SNR things


c	write (6,*) xx,yy,SNR, IRmask

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

		MASK (ii,jj,ib) = iMask

	      endif

	  enddo
	  enddo

 100	continue  !  nsrc


	return
	end
