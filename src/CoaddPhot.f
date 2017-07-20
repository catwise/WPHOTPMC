	subroutine CoaddPhot (Jsrc,ib,ncoaddsize,nx,ny,
     1            Array, Unc, Cov, CoMSK, iCOADD,cscale,
     1            x0,y0, Nap, Raper, nann, sky, sdev,  Nbann,zero,
     1            xint,efluxdn,zmag, zerr, flag, SPIT, MIPS, fbits, Fcorr, Lconf, cov_ave, cov_std)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 Array (ncoaddsize,ncoaddsize,4), Unc(ncoaddsize,ncoaddsize,4) 
	real*4 Cov (ncoaddsize,ncoaddsize,4)
	real*4  Rap(19), Raper(19,4), sky, sdev, zero(5), Lconf
	real*4 x0,y0, ztot(19), cscale
	real*4  xint(19),efluxdn(19),zmag(19), zerr(19), SUM2err(19), Fcorr
	integer nx(4),ny(4), fbits(32)
	integer ntot(19),iflag(19),flag(19)
	integer  iCOADD (ncoaddsize,ncoaddsize,4)
c	integer  MSK(ncoaddsize,ncoaddsize,4)
	integer*2 CoMSK(ncoaddsize,ncoaddsize,4)
	logical SPIT, MIPS

	thresh = 2. * sdev  ! for upper limits  ! now includes the Fcorr scaling

c	if (ib.eq.4) write (6,*) sdev,thresh


	NAPC = Nap

	do Map=1,NAPC
	 flag(Map) = 0
	 Rap(Map) = Raper(Map,ib) 

c	write (6,*) 'radius in pixels ',ib,Rap(Map)

	 if (SPIT) then
                Rap(1) = 2. / cscale ! arcsec to pixels
                if (MIPS) then
                        if (ib.eq.4) Rap(1) = 9. / cscale ! arcsec to pixels
                endif
         endif


	enddo

c	write (6,*) Rap(1), x0,y0, ib, sky, sdev


	call  BIGcphot (Jsrc,ib,ncoaddsize,nx(ib),ny(ib),
     1     Array,Unc, Cov, CoMSK, iCOADD,
     1     NapC, Rap,sky, sdev,
     1     x0,y0,ntot,ztot,xint, flag, iflag, SUM2err, fbits, Lconf, cov_ave, cov_std)

c	write (6,*) ib, zero(ib), xint(1), flag(1), iflag(1)

c	 if (ib.eq.2) then
c                if ((x0.gt.2330.).and.(x0.lt.2350.).and.(y0.gt.2830.).and.(y0.lt.2850.)) then
c                  write (6,*) x0,y0, zero(ib)
c		write (6,*) Array (2347, 2849,ib), Array (2346, 2849,ib), Array (2348, 2849,ib)
c                write (6,*) sky, sdev, napc
c                  write (6,*) rap(1), ntot(1), xint(1), zero(ib) - (2.5*log10(xint(1)))
c        call exit(0)
c                endif
c
c        endif


c	if (ib.eq.4) then
c	  zmag(1) = zero(ib) - (2.5*log10(xint(1)))
c	  write (6,*) ib, zero(ib), xint(1), zmag(1)
c	endif


cvalue    Condition
c-----    ------------------------------------------------------
c  0      nominal
c  1      source confusion
c  2      bad or fatal pixels:  bit 2 or 18
c  4      non-zero bit flag tripped (other than 2 or 18)
c  8      corruption
c 16      saturation
c 32      upper limit

	do Map=1,NAPC

c	if ((MAP.eq.1).and.(ib.eq.1)) write (6,*) ib, x0,y0,Rap(Map),ntot(Map),xint(Map),SUM2err(Map),sky, sdev
c	if (ib.eq.4)  write (6,*) x0,y0,Rap(Map),ntot(Map),ztot(Map),xint(Map),SUM2err(Map)

	if (flag(Map).lt.0) then
		flag(Map) = 1   !  star contamination
	endif

c	if ((ib.eq.1).and.(map.eq.1).and.(Jsrc.eq.13)) then
c		write (6,*) 'debug ',iflag(1),flag(1)
c	endif

	if (iflag(Map).gt.0) then
		flag(Map) = iflag (Map) + flag(Map)
	endif

	zmag(Map) = 99.
	zerr(Map) = 0.
	if (xint(Map).gt.0.) then

	  zmag(Map) = zero(ib) - (2.5*log10(xint(Map)))
	  if (zmag(Map).gt.99.) then
		zmag(Map) = 99.
		flag(Map) = 8 + flag(Map) 
	  endif
	else

	  flag(Map) =  8 + flag(Map)  ! negative flux or corruption

	  ! get an upper limit
	  ncoadd = 1
	  call photom_error (ntot(Map),xint(Map),nann,sdev,Nbann,zerr(Map),efluxdn(Map),ncoadd, SUM2err(Map), Fcorr)
	  flag(Map) =  32  ! upper limit

	  if (efluxdn(Map).gt.0.) then
		 zmag(Map) = zero(ib) - (2.5*log10(2. * efluxdn(Map)))
		 zerr(Map) = 0.
	  endif

c	write (6,*) 'neg flux'
c	write (6,*) x0,y0,Rap(Map),ntot(Map),ztot(Map),xint(Map),SUM2err(Map), zerr(Map),efluxdn(Map)
c	call exit(0)

	endif

	 ncoadd = 1
	 call photom_error (ntot(Map),xint(Map),nann,sdev,Nbann,zerr(Map),efluxdn(Map),ncoadd, SUM2err(Map), Fcorr)

c	write (6,*) 'error done ',ntot(Map),xint(Map),zmag(Map),zerr(Map),efluxdn(Map),sdev
c	if (MAP.eq.2) call exit(0)


c	if (Jsrc.gt.2400) then
c		if (ib.eq.1) then
c		  write (6,*) 'error done ',ntot(Map),xint(Map),zmag(Map),zerr(Map),efluxdn(Map),sdev,nann
c	call exit(0)
c		endif
c
c	endif

	enddo  ! map


c	write (6,66) x0,y0,Rap(ib),ztot,ntot,xint,efluxdn,zmag, zerr, flag
c66	format (2f8.2,f7.2,f10.3,i5,2f10.1,1x,2f8.3, i5)
c	if ((zmag.eq.99.).and.(iflag.eq.1)) then
c		write (6,*) ztot,xint,zmag,flag
c	endif
c	if (ib.eq.1) call exit (0)



	return
	end
