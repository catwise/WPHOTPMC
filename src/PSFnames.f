	subroutine PSFnames (psfdir,calbname,calgridX,calgridY, Pnames,
     *	    Punc,npn,MapPSFs,nx,ny, MIPS, mipspsf, mipsunc, nfi, wflag)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*(*) psfdir,calbname, Pnames(100,4), Punc(100,4), mipspsf, mipsunc
	character*500 string,fname
	character*25 s0,gridsz,gloc
	character*2 band
	character*1 c1,c2

	integer calgridX,calgridY,npn(4)
	integer(2) MapPSFs(nx,ny,4), wflag(nfi,4)

	logical MIPS

	MapPSFs = 1
	write (c1,'(i1)') calgridX
	write (c2,'(i1)') calgridY

	gridsz = '0' // c1 // 'x0' // c2

	write (6,'(a)') gridsz

	L = numchar (psfdir)
	M = numchar (calbname)

c construct names

	ibmax = 4
	if (MIPS) then
		ibmax = 3
	endif


	do ib=1,ibmax


c	 if ( wflag ( 1,ib ) .eq.1 ) then
	 if ( any ( wflag ( 1:nfi,ib ) ==1 ) ) then


	 ib0 = ib
	 if ((MIPS).and.(ib.eq.3)) then
		ib0 = 4
	 endif

	 write (band,'(a,i1)') 'w',ib0

	 ndim = 0

	 do j=1,calgridY
	   write (c2,'(i1)') j
	 do i=1,calgridX
	   ndim = ndim + 1

	   write (c1,'(i1)') i

	   gloc = '0' // c1 // 'x0' // c2

	   Pnames (ndim,ib) = psfdir(1:L) // '/' // calbname(1:M) //
     1      '-' // band //  '-psf-wpro-' // gridsz(1:5) //
     1      '-' // gloc(1:5) // '.fits'

	   Punc(ndim,ib) = psfdir(1:L) // '/' // calbname(1:M) //
     1      '-' // band //  '-psfunc-wpro-' // gridsz(1:5) //
     1      '-' // gloc(1:5) // '.fits'

c	   write (6,'(a)') Pnames (ndim,ib)
c	   write (6,'(i4,2x,a)') ndim,Punc(ndim,ib)(1:92)

	   call headpar(nxp,nyp,Pnames(ndim,ib),crval1,crval2,cdelt1,
     *		cdelt2,crot,crpix1,crpix2,xlo,xhi,ylo,yhi,istat)
	   ilo = max(nint(xlo),1)
	   ihi = min(nint(xhi),nx)
	   jlo = max(nint(ylo),1)
	   jhi = min(nint(yhi),ny)

	   do jj = jlo,jhi
	   do ii = ilo,ihi
		MapPSFs(ii,jj,ib) = ndim
	   enddo
	   enddo

	 enddo
	 enddo

	npn(ib) = ndim

	endif

	enddo ! ib

	if (MIPS) then

	 M = numchar(mipspsf)
	 ndim = 1
	 ib = 4
	 Pnames (ndim,ib) = psfdir(1:L) // '/' // mipspsf(1:M) 

	 Pnames (ndim,ib) = mipspsf(1:M)

	 M = numchar(mipsunc)
	 Punc(ndim,ib) = psfdir(1:L) // '/' // mipsunc(1:M)

	 Punc(ndim,ib) = mipsunc(1:M)

	 npn(ib) = ndim

	endif


	return
	end
