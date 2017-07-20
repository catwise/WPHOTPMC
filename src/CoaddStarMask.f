	subroutine CoaddStarMask (ib,NSUMmax,NsrcAll,
     1     ncoaddsize,nx, ny, WCSvb,fwhm,cscale,
     1     nf, MASK, Rstann, Rstwid, RAlist, Declist, SPIT)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (FWHM_pix = 4.0, scale = 1.5)

	real*8 ra0,dec0, RAlist(NSUMmax), Declist(NSUMmax)
c	real*4 SNRlist (NSUMmax,nf,4)
	real*4 scalefact (4), Rstann(4), Rstwid(4)
	real*4 fwhm(4),cscale(4)
	integer MASK(ncoaddsize,ncoaddsize,4), imask
	integer nx(4),ny(4),nsrc, WCSvb (4)
	logical debug, SPIT

	data scalefact/1.1,1.1,1.0,1.0/    !  scale the FWHM by this number to get the band value
c 	debug = .true.
	debug = .false.


c	write (6,*) 'inside ',ib

cc  initialize mask
c	FW = FWHM_pix
	FW = fwhm (ib) / cscale (ib)

c	write (6,*)  fwhm (ib), cscale (ib), FW
c	call exit(0)

c	if (SPIT) FW = 2.5

	imask = 0

 	Rmask = FW * scalefact (ib) * scale   ! masking radius in pixels
        IRmask = nint(Rmask)

	Rl = (Rstann(ib)*2.) - 1.  !  2x times the frame value
        Rh = Rl + (Rstwid(ib)*2.)  !  2x times the frame value
	Rh = Rh + 1.

c initialize the mask
	 do jj=1,ny(ib)
	 do ii=1,nx(ib)
	  Mask(ii,jj,ib) = imask
	 enddo
	 enddo

	ppp = 1.375
	edgebuf = 0.

c	write (6,*) 'A'

	do 100 Jsrc = 1,NsrcAll

	  Ra0 = RAlist(Jsrc) * 1.d0
          Dec0 = Declist(Jsrc) * 1.d0
	  igo = 1

c	write (6,*) Jsrc,nx ( ib ) ,ny ( ib )
c	write (6,'(2f10.5)') Ra0,Dec0
	
	  call WhatPos ( 1,WCSvb(ib) ,nx ( ib ) ,ny ( ib ) , ppp, '0', ra0, dec0, xx, yy, igo, edgebuf )

c	write (6,*) ra0,dec0,xx,yy
c	 write (69,*) ra0,dec0,xx,yy

	  ix = nint(xx)
	  if (ix.lt.1) goto 100
	  if (ix.gt.nx(ib)) goto 100

          jy = nint(yy)
	  if (jy.lt.1) goto 100
          if (jy.gt.ny(ib)) goto 100

c	  SNR = SNRlist (Jsrc,1,ib)
c	  if (SNR.lt.3.0) goto 100  ! avoid low SNR things

c	write (6,*) Jsrc, ix,jy, IRmask

c found source; mask from image

	  iil = (ix - IRmask) - 1
	  iil = max (iil,1)
	  iih = (ix + IRmask) + 1
	  iih = min (iih, nx(ib))
	  jjl = (jy - IRmask) - 1
	  jjl = max (1, jjl)
          jjh = (jy + IRmask) + 1
	  jjh = min (jjh, ny(ib))

c	write (6,*) jjl,jjh, iil,iih

c	  n=n+1     ! JWF B21114 [never initialized, never used]

	  do   jj=jjl,jjh
	      dy = (jj*1.) - yy
	      dy2 = dy**2
	  do   ii=iil,iih
	      dx = (ii*1.) - xx
	      dx2 = dx**2
	      dr = sqrt ( dx2 + dy2 )

	      if (dr.le.Rmask) then

		if (MASK (ii,jj,ib).eq.0) then
		  MASK (ii,jj,ib) = Jsrc
		else
		  MASK (ii,jj,ib) = abs(MASK (ii,jj,ib)) + Jsrc
		  MASK (ii,jj,ib) = -MASK (ii,jj,ib)
		endif

c	write (6,*) Jsrc, ii,jj,MASK (ii,jj,ib)

	      endif

	  enddo
	  enddo

 100	continue  !  nsrc

c	write (6,*) 'done'

	return
	end
