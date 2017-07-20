cccccccccccccccccccccccc
c image moment

	subroutine imoment (nx,ny,nsx,nsy,array,x0,y0,
     1    radius,iorder,zmomx,zmomy)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        dimension array(nx,ny)

	il = nint(x0 - radius)
        ih = nint(x0 + radius)
        jl = nint(y0 - radius)
        jh = nint(y0 + radius)

        if (il.lt.1) il=1
        if (jl.lt.1) jl=1
        if (ih.gt.nsx) ih=nsx
        if (jh.gt.nsy) jh=nsy

c	write (6,*) il,ih,jl,jh

	sumx = 0.0
	tot = 0.0
	sumy = 0.
        n = 0
        do j = jl,jh
          dy = (j*1.) - y0
        do 50 i = il,ih
          dx = (1.*i) - x0
          dr = ( (dx**2) + (dy**2) ) ** 0.5

          if (dr.gt.radius) goto 50

          val = array(i,j)
          if (val.lt.-100) goto 50
c	  if (val.eq.0.) goto 50

	  termx = (dx ** iorder) *  val
	  termy = (dy ** iorder) *  val

	   sumx = sumx + termx
	   sumy = sumy + termy
	n=n+1
	   tot = tot + val

 50     continue
        enddo


c	if (tot.gt.0.) then
	if (n.gt.1) then
		zmomx = sumx / tot
		zmomy = sumy / tot
	else
		zmomx = 0.
		zmomy = 0.
	endif

	return
	end

	
