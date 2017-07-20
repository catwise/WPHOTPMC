	subroutine centroid (ib,Jfr,nf,nnx,nny,nx,ny,array,x,y,ibox,cx,cy)
c  centroid determination;
	
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 array(nnx,nny,nf,4)
c define aperature size

	iflag = 0

	nbox = ibox 
	ip = nint(x)
	jp = nint(y)
	
	rmax = nbox/2.

 47	jmax = jp + (nbox/2)
	jmin = jp - (nbox/2)

	imax = ip + (nbox/2)
	imin = ip - (nbox/2)
	
	if (jmax.gt.ny) iflag = 1
	if (imax.gt.nx) iflag = 1
	if (jmin.lt.1) iflag = 1
	if (imin.lt.1) iflag = 1

	if (iflag.eq.1) then
		return
	endif

	sumx = 0.
	sumy = 0.
	sum2 = 0.

	do j=jmin,jmax
	    dy = (j - y)*1.
	do i = imin,imax
	    dx = (i - x)*1.

		val = array(i,j,Jfr,ib) 
		if (val.gt.0.) then
c			sumx = sumx + (1. * i * val)
c			sumy = sumy + (1. * j * val)				
			
			sumx = sumx + (dx * val)
			sumy = sumy + (dy * val)
			sum2 = sum2 + val
		
c		else
c			! bad pixels; give up
c			cx = x
c			cy = y
c			return
		endif
	enddo
	enddo

	if (sum2.le.0.) then
		return
	else
c		cx = sumx / sum2
c		cy = sumy /sum2

		dcx = sumx / sum2
		dcy = sumy /sum2

		cx = x + dcx
		cy = y + dcy

c	if (cx.lt.-500.) write (6,*) x,y,dcx, dcy, cx,cy

		if ((abs(dcx).gt.rmax).or.(abs(dcy).gt.rmax)) then
                        cx = x
                        cy = y
                endif


c		dx = abs(cx-x)
c		dy = abs(cy-y)

c		if ((dx.gt.rmax).or.(dy.gt.rmax)) then
c			cx = x
c			cy = y
c		endif


	endif

	return
	end
      
