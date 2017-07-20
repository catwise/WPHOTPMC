	subroutine interp (xl,xh,yl,yh,yval,xval)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


	xval = 0.	
	xt = xh - xl
	yt = yh - yl

	if (yt.eq.0.) goto 99

	dy = yval - yl

	ratio = dy / yt

	xval = (ratio * xt) + xl

 99	i=1

	return
	end

	subroutine interpx (xl,xh,yl,yh,xval,yval)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


	yval = 0.	
	xt = xh - xl
	yt = yh - yl

	if (xt.eq.0.) goto 99

	dx = xval - xl

	ratio = dx / xt

	yval = (ratio * yt) + yl

 99	i=1

	return
	end
