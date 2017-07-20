	subroutine GetApCorr (x,y,ib,nx,ny,calgridX,calgridY,Ncorrdim,AppCorr,ndim,correction_mag)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 AppCorr(Ncorrdim,4)
	integer calgridX,calgridY,Ncorrdim


	xgrid = (nx*1.) / (calgridX*1.)
	ygrid = (ny*1.) / (calgridY*1.)

	igrid = int(xgrid)
	jgrid = int(ygrid)

	idex = 1 + int (x / xgrid)
	jdex = 1 + int (y / ygrid)

	idex = max (idex,1)
	idex = min (idex, calgridX)
	jdex = max (jdex,1)
        jdex = min (jdex, calgridY)


	ndim = 0

        do j=1,calgridY
        do i=1,calgridX
	  ndim = ndim + 1
	  if ((i.eq.idex).and.(j.eq.jdex)) then
		correction_mag =  AppCorr(ndim,ib) 
		return
	  endif
	enddo
	enddo


	write (6,*) 'ERROR -- aperture correction mag not found'
	call exit(9)

	return
	end

