	subroutine CircPhot (Jsrc,J,nf,nfpix,ib,nnx,nny,nx,ny,
     1            Array,MASK,iMask, Unc,
     1            x0,y0, Nap, Raper, nann, sky, sdev, Nbann, zero,
     1            xint,efluxdn,zmag, zerr, flag, ntot, SPIT, MIPS, fbits)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 Array (nnx,nny,nfpix), Unc(nnx,nny,nfpix) 
	real*4  Rap(19), Raper(19,4), sky, sdev, zero(4)
	real*4 x0,y0, ztot(19)
	real*4  xint(19),efluxdn(19),zmag(19), zerr(19), SUM2err(19)
	integer nx(4),ny(4), iMask(nnx,nny,nfpix)
	integer ntot(19),iflag(19),flag(19), fbits(32)
c	integer MASK(nnx,nny,nf,4)
	integer MASK(1,1,1) ! Mask being dummied out since it's never assigned anything - TPC
	logical SPIT,MIPS


	thresh = 2. * sdev  ! for upper limits

	do Map=1,Nap
	 flag(Map) = 0
	 Rap(Map) = Raper(Map,ib)

	enddo

	call  cphot (Jsrc,J,nf,nfpix,ib,nnx,nny,nx(ib),ny(ib),
     1     Array,MASK,iMask,Unc,
     1     Nap, Rap,sky, sdev,
     1     x0,y0,ntot,ztot,xint, flag, iflag, SUM2err, SPIT, MIPS, fbits)

c	write (6,*) 'cphot'
	


	do Map=1,nap

c	  write (6,*) x0,y0,Rap(Map),ntot(Map),ztot(Map),xint(Map),SUM2err(Map)

	if (flag(Map).lt.0) then
		flag(Map) = 1   !  star contamination
	endif

	if (iflag(Map).gt.0) then
		flag(Map) = iflag (Map)
	endif

	zmag(Map) = 99.
	zerr(Map) = 0.
	if (xint(Map).gt.0.) then
	  zmag(Map) = zero(ib) - (2.5*log10(xint(Map)))
	  if (zmag(Map).gt.99.) then
		zmag(Map) = 99.
		flag(Map) = 5
	  endif
	else
c	  flag(Map) = 5   ! negative flux or corruption

c	write (6,*) 'neg flux'
c	write (6,*) x0,y0,Rap(Map),ntot(Map),ztot(Map),xint(Map),SUM2err(Map)

	endif

	 ncoadd = 1
	 Fcorr = 1.
	 call photom_error (ntot(Map),xint(Map),nann,sdev,Nbann,zerr(Map),efluxdn(Map),ncoadd, SUM2err(Map), Fcorr)
c	write (6,*) 'error done ',ntot(Map),xint(Map),zmag(Map),zerr(Map),efluxdn(Map)
c	if (MAP.eq.2) call exit(0)

c	write (6,*) sdev,Nbann,zerr(Map),efluxdn(Map)
c	if (ib.eq.4) call exit(0)

	enddo  ! map


c	write (6,66) x0,y0,Rap(ib),ztot,ntot,xint,efluxdn,zmag, zerr, flag
c66	format (2f8.2,f7.2,f10.3,i5,2f10.1,1x,2f8.3, i5)
c	if ((zmag.eq.99.).and.(iflag.eq.1)) then
c		write (6,*) ztot,xint,zmag,flag
c	endif
c	if (ib.eq.1) call exit (0)



	return
	end
