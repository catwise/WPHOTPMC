	subroutine smoment (nx,ny,array,x0,y0,
     1     radius,zmoment,rat,phimax)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	dimension array(nx,ny)
	logical doell

	doell=.true.
	if (rat.ge.0.95) doell=.false.

c annlus: find median value

	il = nint(x0 - radius)
	ih = nint(x0 + radius)
	jl = nint(y0 - radius)
	jh = nint(y0 + radius)

	if (il.lt.1) il=1
	if (jl.lt.1) jl=1
	if (ih.gt.nx) ih=nx
	if (jh.gt.ny) jh=ny

	sum = 0.0
	tot=0.
	n = 0
	do j = jl,jh
	  dy = (j*1.) - y0
	do 50 i = il,ih
	  dx = (1.*i) - x0

	  if (doell) then
                        call w_ell (rat,phimax,dx,dy,dr)
          else
                        dr = ( (dx**2) + (dy**2) ) ** 0.5
          endif

	  if (dr.gt.radius) goto 50
 
	  val = array(i,j) 
	  if (val.lt.-100) goto 50

	  term = dr *  val
	  sum = sum + term
	  tot=tot+val

 50  	continue
	enddo


c normalize by integrated flux

	if (tot.gt.0.) then
	  zmoment = sum / tot
	  zmoment = min (zmoment,80.)
	else
	  zmoment=0.
	endif

	return
	end

