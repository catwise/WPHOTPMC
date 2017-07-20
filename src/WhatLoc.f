	subroutine WhatLoc (nx,ny, pscale, x0, y0, igo, barc)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	igo = 0

c	barc = 15. ! arcsc

	buffer = barc / pscale  ! pixels
	
	xhigh = (nx*1.) - buffer
	yhigh = (ny*1.) - buffer

	if (x0.lt.buffer) return
	if (x0.gt.xhigh) return
	if (y0.lt.buffer) return
        if (y0.gt.yhigh) return

	igo = 1

	return
	end
