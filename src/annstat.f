

	subroutine annstat (nx,ny,array,x0,y0,
     1    rinner,router,ave,zmed,z_68,z_87,sdev)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	dimension array(nx,ny)
	real*4 z(5000)

c annlus: find median (or quartile) value

	il = nint(x0) - router
	ih = nint(x0) + router
	jl = nint(y0) - router
	jh = nint(y0) + router

	if (il.lt.1) il=1
	if (jl.lt.1) jl=1
	if (ih.gt.nx) ih=nx
	if (jh.gt.ny) jh=ny


        zl = -50.
        zh = 5000.

	n = 0
	do  j = jl,jh
	  dy = (j*1.) - y0
	do 50 i = il,ih
	  dx = (1.*i) - x0
	  dr = ( (dx**2) + (dy**2) ) ** 0.5
 
	  val = array(i,j) 

	  if ((dr.lt.rinner).or.(dr.gt.router)) goto 50
	  if ((val.lt.zl).or.(val.gt.zh)) goto 50
	  if (n.gt.5000) goto 50

	  n = n + 1
	  z(n) = val

 50  	continue
	enddo
	if (n.gt.10) then
	  call IRASTAT (z,n,zmed,z_68,z_87)
	  call MOMENT(z,N,AVE,SDEV)
	else	
		zmed = 0.
		z_68 = 0.
		z_87 = 0.
		ave = 0.
	endif


  99 	i=1


	return
	end


