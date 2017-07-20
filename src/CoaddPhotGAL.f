	subroutine CoaddPhotGAL (Jsrc,ib,ncoaddsize,nx,ny,
     1            Array, Unc,  Cov, CoMSK,  iCOADD,pscale,cscale,
     1            x0,y0, Nap, Raper, ba, pa_math,  nann, sky, sdev,  Nbann,zero,
     1            zmag, zerr, flag, SPIT, MIPS, fbits, Fcorr, Lconf)


	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 Array (ncoaddsize,ncoaddsize,4), Unc(ncoaddsize,ncoaddsize,4) 
	real*4 Cov (ncoaddsize,ncoaddsize,4)
	real*4  Rap, Raper, sky, sdev, zero(5), Lconf
	real*4 x0,y0, ztot, pscale,cscale
	real*4  xint,efluxdn,zmag, zerr, SUM2err, Fcorr
	integer nx(4),ny(4), fbits(32)
	integer ntot,iflag,flag
	integer  iCOADD (ncoaddsize,ncoaddsize,4)
	integer*2 CoMSK(ncoaddsize,ncoaddsize,4)
c	integer  MSK(ncoaddsize,ncoaddsize,4)
	logical SPIT, MIPS

	thresh = 2. * sdev  ! for upper limits  ! now includes the Fcorr scaling

c	if (ib.eq.4) write (6,*) sdev,thresh

	pscale_factor = pscale / cscale

	NAPC = Nap

	flag = 0


	call  GALcphot (Jsrc,ib,ncoaddsize,nx(ib),ny(ib),
     1     Array,Unc, Cov, CoMSK,iCOADD,
     1     NapC, Raper, ba, pa_math, sky, sdev,
     1     x0,y0,ntot,ztot,xint, flag, iflag, SUM2err, fbits, Lconf)


	if (flag.lt.0) then
		flag = 1   !  star contamination
	endif

	if (iflag.gt.0) then
		flag = iflag + flag
	endif

	zmag = 99.
	zerr = 0.
	if (xint.gt.0.) then
	  zmag = zero(ib) - (2.5*log10(xint))
	  if (zmag.gt.99.) then
		zmag = 99.
c		flag = 5
		flag = 8 + flag
	  endif
	else
		flag =  8 + flag  ! negative flux or corruption
	endif

	ncoadd = 1
	call photom_error (ntot,xint,nann,sdev,Nbann,zerr,efluxdn,ncoadd, SUM2err, Fcorr)


	return
	end
