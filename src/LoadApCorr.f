	subroutine LoadApCorr (Rstap,apcorrdir,calbname,calgridX,calgridY,Ncorrdim,AppCorr,napcorr,smode,verb)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*(*) apcorrdir,calbname
	character*500 string,fname
	character*25 s0,gridsz,gloc
	character*2 band
	character*1 c1,c2

c	real*4 AppCorr(999,2,Ncorrdim,4)
	real*4 AppCorr(Ncorrdim,4), R(4)
	integer calgridX,calgridY,Ncorrdim

	logical verb,smode

	do ib=1,3
	 R(ib) = Rstap  ! standard aperture
	enddo
	R(4) = Rstap * 2.  !  W4 has double the pixel scale


	write (c1,'(i1)') calgridX
	write (c2,'(i1)') calgridY

	gridsz = '0' // c1 // 'x0' // c2

c	if (verb) write (6,'(a)') gridsz

	L = numchar (apcorrdir)
	M = numchar (calbname)

c construct names

	do ib=1,4

	write (band,'(a,i1)') 'w',ib

	ndim = 0

	do j=1,calgridY
	   write (c2,'(i1)') j
	do i=1,calgridX
	   ndim = ndim + 1

	   write (c1,'(i1)') i

	   gloc = '0' // c1 // 'x0' // c2

	   fname = apcorrdir(1:L) // '/' // calbname(1:M) //
     1      '-' // band //  '-apcorr-' // gridsz(1:5) //
     1      '-' // gloc(1:5) // '.tbl'

c	   if (verb) write (6,'(a)') fname(1:92)

	   Rclose = 9999.
	   correction_mag = 0.

	   open (unit=11,file=fname)
	   napcorr = 0
	   do 100 K=1,999
	    read (11,'(a)',end=99) string
	    if (string(1:1).eq.'|') goto 100
	    if (string(1:1).eq.'\') goto 100
	    if (string(1:1).eq.'/') goto 100

	    napcorr = napcorr + 1

	    read (string,*) rap,rarc,napp,fapcorr,zmag, delta_mag

	    Rdiff = abs(rarc - R(ib))
	    if (Rdiff.lt.Rclose) then
		 Rclose = Rdiff
		  correction_mag = delta_mag
	    endif
		
 100      continue
 99       close (11)

	  AppCorr(ndim,ib) = correction_mag
c	  if (smode) AppCorr(ndim,ib) = 0.0

c	    AppCorr(napcorr,1,ndim,ib) = rarc  ! radius in arcsec
c	    AppCorr(napcorr,2,ndim,ib) = delta_mag   ! mag  correction (add this to the measured mag to get the corrected mag)

c 100	  continue
c99	  close (11)


	  if (napcorr.eq.0) then
		write (6,*) 'ERROR -- napcorr = ',napcorr
		call exit (9)
	  endif

c	  write (6,'(a,i2,a,i3,a,a,a,i2,a)') 'band=',ib,' # = ',napcorr,
c     1      '  apercorr points loaded for grid loc: ',
c     1      gloc(1:5),' (',ndim,')'

	enddo
	enddo


	enddo

	return
	end

