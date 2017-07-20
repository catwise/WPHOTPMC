	subroutine find_peak (J,ib,nf,nfpix,nnx,nny,nx,ny,array,ibox,x0,y0,ip,jp,valmax,ibog)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 Array (nnx,nny,nfpix)

	ibog = 0
	ip = nint(x0)
	jp = nint(y0)

	jmax = nint(y0) + nint(ibox/2.)
	jmin = nint(y0) - nint(ibox/2.)
	imax = nint (x0) + nint(ibox/2.)
	imin = nint (x0) - nint(ibox/2.)

	if (jmax.gt.ny) jmax=ny
	if (imax.gt.nx) imax=nx
	if (jmin.lt.1) jmin=1
	if (imin.lt.1) imin=1

	valmax = -900.
	do    jj = jmin,jmax
	do 50 ii = imin,imax
		val = array(ii,jj,J)
		if (val.lt.-99.) ibog=1
		if (val.gt.valmax) then
			valmax = val
			ip = ii
			jp = jj
		endif
 50	continue
	enddo

	return
	end
