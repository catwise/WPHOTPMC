cccccccccccccccccccccccc
c image moment

	subroutine rot_moment (nx,ny,array,x0,y0,
     1    radius,iorder,zminor,zmajor,anglemin)


	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        real*4 array(nx,ny)

	il = nint(x0 - (radius*sqrt(2.)))
        ih = nint(x0 + (radius*sqrt(2.)))
        jl = nint(y0 - (radius*sqrt(2.)))
        jh = nint(y0 + (radius*sqrt(2.)))

        if (il.lt.1) il=1
        if (jl.lt.1) jl=1
        if (ih.gt.nx) ih=nx
        if (jh.gt.ny) jh=ny

	dangle = 2.0
	if (anglemin.eq.0.) then

	  angle = 0. - dangle
	  nstep= (90./dangle) + 1

	else
	  angle = anglemin-dangle
	  nstep=1
	endif

	anglemin = 0.
        zminor=0.
        zmajor=0.
        ratmin = 9999.
	
	do 400 KKK = 1,nstep

	angle=angle+dangle
c	if (angle.gt.90.) goto 400

	sum = 0.
	sumx = 0.0
	tot = 0.0
	sumy = 0.
        n = 0
        do j = jl,jh
          dy = (j*1.) - y0
        do 50 i = il,ih
          dx = (1.*i) - x0

c rotate axis
	  zang = angle
	  call rot_axis (dx,dy,zang,dxP,dyP)
	  dr = ( (dxP**2) + (dyP**2) ) ** 0.5
 
c         dr = ( (dx**2) + (dy**2) ) ** 0.5

          if (dr.gt.radius) goto 50

          val = array(i,j)
          if (val.lt.-100) goto 50
c	  if (val.eq.0.) goto 50

	  termx = (dxP ** iorder) *  val
	  termy = (dyP ** iorder) *  val

	   sumx = sumx + termx
	   sumy = sumy + termy
	n=n+1
c	write (6,'(i5,2f7.1,4f8.2)') n,dx,dy,termx,termy,sumx,sumy
	   tot = tot + val

 50     continue
        enddo

c	write (6,*) sumx,sumy

	if (tot.ne.0.) then
		zmomx = sumx / tot
		zmomy = sumy / tot

		if (zmomx.gt.zmomy) then
			ratty = (zmomy/zmomx) ** (1. / iorder)
			Angmaj = angle
		else
			ratty = (zmomx/zmomy) ** (1. / iorder)
			Angmaj = angle+90.0
		endif

	else
		zmomx = 0.
		zmomy = 0.
		ratty = 999.
		Angmaj = 0
	endif

		if (ratty.lt.ratmin) then
			anglemin = Angmaj
			zminor = min(zmomx,zmomy)
			zmajor = max(zmomx,zmomy)
			ratmin = ratty
		endif


c	if (iorder.eq.2) then
c	  zz = angle + 90.0
c          if (zz.gt.90.) zz = zz - 180.
c	  if (zz.eq.90.0) zz=-90.0
c	  write (62,*) zz,zminor,zmajor,ratty,anglemin
c	endif

c	write (6,*) zminor,zmax,anglemin

 400	continue
 	

	return
	end

c    rotate axis ANGLE

	subroutine rot_axis (x0,y0,angle,x,y)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 x0,y0,angle,x,y


	if (angle.eq.0.) then
		x = x0
		y = y0
		goto 47
	else if (abs(angle).eq.90.) then
		x = y0
		y = -x0
		goto 47
	else if (abs(angle).eq.180.) then
		x = -x0
		y = -y0
		goto 47
	else if (abs(angle).eq.270.) then
		x = y0
		y = x0
		goto 47
	endif

	t1 = y0 / sin(angle/57.2957795)
	t2 = x0 / cos(angle/57.2957795)
	t3 = tan(angle/57.2957795)
	t4 = 1./tan(angle/57.2957795)


c	write (6,*) angle, t1,t2,t3,t4

	term1 = t1 - t2
	term2 = t3 + t4

	Y = term1 / term2


	t1 = x0 / cos(angle/57.2957795)
	t2 = Y * tan(angle/57.2957795) 

	X = t1 + t2


 47	return
	end

