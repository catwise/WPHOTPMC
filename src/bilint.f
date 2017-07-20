	real*4 function bilint(a,nxdim,nydim,maxsegs,maxbands,nx,ny,k,ib,x,y)

! Linearly interpolate array A at pixel location (x,y) which may be fractional.

	implicit real*4(a-h,o-z)
	implicit integer(i-n)
	real*4 a(nxdim,nydim,maxsegs,maxbands),x,y

	ix = ifix(x*1.0)
	iy = ifix(y*1.0)

	if (ix < 1 .or. iy < 1 .or. ix >= nx .or. iy >= ny) then
	    bilint = 0.
	    return
	endif

	aiy = a(ix,iy,k,ib) + (x - ix)*(a(ix+1,iy,k,ib)-a(ix,iy,k,ib)) 
	aiyp = a(ix,iy+1,k,ib) + (x - ix)*(a(ix+1,iy+1,k,ib)-a(ix,iy+1,k,ib))
	bilint = aiy + (y - iy)*(aiyp - aiy)
	return

	end
