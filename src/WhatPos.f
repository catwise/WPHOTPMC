	subroutine WhatPos (Jfr,WCS,nx,ny, pscale, Hfits, Ra8, Dec8, x0, y0, igo, barc)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*8 ra8,dec8,x8,y8
	integer WCS, offscl
	character*(*) Hfits


c  RA,DEC >> x,y

c	write (6,*) 'here WhatPos'

	
	if (WCS.lt.0) then

 47	  call wcsinit(Hfits,WCS)

c	write (6,*) 'ere ',WCS


   	  if ( WCS .lt. 0) then
               print *,'*** WCSinit ERROR ',WCS
              call exit(1)
          endif

	endif


        offscl = -1
        call wcs2pix(WCS, ra8, dec8, x8, y8, offscl)


        if(offscl .lt. 0) then
                 print *,'*** WCSPlay ERROR ',offscl,': Can''t convert.'
                 call exit(1)
         endif

	x0 = x8 * 1.0
	y0 = y8 * 1.0


c   	call wcsclose(WCS)

	buffer = barc / pscale  ! pixels
c	buffer = 4  ! pixels
	
	xhigh = (nx*1.) - buffer
	yhigh = (ny*1.) - buffer

c	if (Jfr.eq.13) then
c		write (6,*) nx,ny,buffer,x0,y0,xhigh,yhigh
c	endif
		

	if (x0.lt.buffer) return
	if (x0.gt.xhigh) return
	if (y0.lt.buffer) return
        if (y0.gt.yhigh) return

	igo = 1

	return
	end
