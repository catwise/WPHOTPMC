	subroutine mannuli (nx,ny,array,x0,y0,rinner,
     1     router,zmed,zmean,zstd,floor)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	dimension array(nx,ny)
	real*4 z(5000)
	logical doquart

	doquart = .false.
	if (zmed.lt.0.) doquart = .true.

	zstd = zstd * 5.

c annlus: find median (or quartile) value

	il = nint(x0) - router
	ih = nint(x0) + router
	jl = nint(y0) - router
	jh = nint(y0) + router

	if (il.lt.1) il=1
	if (jl.lt.1) jl=1
	if (ih.gt.nx) ih=nx
	if (jh.gt.ny) jh=ny


	do jit=1,4
                if ((jit.eq.1).and.(zmed.le.0.)) then
                        zl = -50.
                        zh = 5000.
                else
                        fact = 3.
                        zl = (-fact * zstd) + zmed
                        zh =  (fact * zstd) + zmed
                endif

		n = 0
		do j = jl,jh
		  dy = (j*1.) - y0
		do 50 i = il,ih
		  dx = (1.*i) - x0
		  dr = ( (dx**2) + (dy**2) ) ** 0.5
  
		  val = array(i,j) 

		  if ((dr.lt.rinner).or.(dr.gt.router)) goto 50
		  if (val.lt.floor) goto 50
		  if (val.lt.zl) goto 50
		  if (val.gt.zh) goto 50

		  n = n + 1
		  if (n.gt.5000) goto 50
		  z(n) = val

 50	  	continue
		enddo


		if (n.gt.2) then
		  call MOMENT(z,n,zmean,zstd)
		  call MDIAN1(z,n,zmed)
		  call QUART (z,n,zquart)
		else if (n.ge.1) then
			zmean = z(1)
			zmed = z(1)
			zstd = 0.
			zquart = zmed
			goto 99
		else	
			zmean = -99.
			zmed = -99.
			zstd = -99.
			zquart = -99.
			goto 99
		endif


  	enddo
  99 	i=1

	if (doquart) zmed = zquart

	return
	end


