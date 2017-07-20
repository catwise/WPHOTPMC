	subroutine PosCoMoment (ib,nco,nx,ny,array,x,y,ibox,cx,cy,findpeak)
c  centroid determination;
	
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 array(nco,nco,4)
	logical findpeak

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


	if (findpeak) then

	 ipeak = ip
	 jpeak = jp
	 pval = 0.
	 do j=jmin,jmax
         do i = imin,imax
	   val = array(i,j,ib)
	   if (val.ge.pval) then
		ipeak = i
		jpeak = j
		pval = val
	   endif
	 enddo
	 enddo

	 ip = ipeak
	 jp = jpeak

	 x = ip * 1.
	 y = jp * 1.
	 findpeak = .false.
	 goto 47


	endif





	sumx = 0.
	sumy = 0.
	sum2 = 0.

c	write (6,*) jmin,jmax, imin,imax

	do j=jmin,jmax
c	    dy = (j - y)*1.
	    dy = (j - jp) *1.
	do i = imin,imax
c	    dx = (i - x)*1.
	    dx = (i - ip) * 1.

		val = array(i,j,ib) 
		if (val.gt.0.) then
c			sumx = sumx + (1. * i * val)
c			sumy = sumy + (1. * j * val)				
			
			sumx = sumx + (dx * val)
			sumy = sumy + (dy * val)
			sum2 = sum2 + val
c	if (ib.eq.3) write (6,*) i,j,dx,dy,sumx,sumy
		
		else
			! bad pixels; give up

			cx = x
			cy = y
			return
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
c	write (6,*) 'hit the max ',dcx,dcy,rmax
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

c	shiftx = dcx
c	shifty = dcy

c	write (6,*) 'shiftx,shifty',shiftx,shifty

	return
	end
      
