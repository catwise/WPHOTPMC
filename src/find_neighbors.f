	subroutine find_neighbors(n,RAlist,DEClist,MAGSTD,eMAGSTD,onframe,
     *	    nmax,nf,maxbands,nsrc,blendist,pratio,nbmax,nblend,neighbors)

! Find neigbhors of nth candidate source.

	implicit real*4 (a-h,o-z)
	implicit integer (i-n)

	real, parameter :: dtor = 0.0174533

   	real(8) RAlist(*), DEClist(*), dra
	real(4) MAGSTD(nmax,4),eMAGSTD(nmax,4),dmag(4)
	integer neighbors(*)
	logical(1) onframe(nmax,nf,maxbands)

	rcritsq = (blendist/3600.)**2
	dmagmax = -2.5*alog10(pratio*1.0)
	cd = cos(DEClist(n)*dtor)
	nblend = 1
	neighbors(nblend) = n
	if (maxval(eMAGSTD(n,1:4)) == 0) return ! 0 means no value
	i = n

	do while (i < nsrc .and. nblend < nbmax)
	    i = i+1
	    if (any(onframe(i,1:nf,1:maxbands))) then
	        dra = abs(RAlist(i)-RAlist(n))
		if(dra .ge. 180)  dra = dra - 360.0
		rsq = (dra*cd)**2 +
     *		    (DEClist(i)-DEClist(n))**2
		dmag = MAGSTD(i,1:4) - MAGSTD(n,1:4)
		where(eMAGSTD(i,1:4) == 0) dmag = 1.e35

		if (maxval(eMAGSTD(i,1:4)) > 0) then
		  if (rsq <= rcritsq .and. minval(dmag) <= dmagmax) then
		    nblend = nblend + 1
		    neighbors(nblend) = i
		  endif
		endif
	    endif
	enddo

	return

	end

c
c Alternative version ...
c

	subroutine find_neighbors_perband(
     *      n,RAlist,DEClist,MAGSTD,eMAGSTD,onframe,
     *	    nmax,nf,maxbands,nsrc,blendist,pratio,nbmax,nblend,neighbors
     *     )

! Find neigbhors of nth candidate source as above, but using a smarter search
! that depends on band and mag availability

	implicit real*4 (a-h,o-z)
	implicit integer (i-n)

	real, parameter :: dtor = 0.0174533
        real, parameter :: eMagBad = 0

   	real(8) RAlist(*), DEClist(*), dra
	real(4) MAGSTD(nmax,4),eMAGSTD(nmax,4),dmag(4),blendist(4)
	integer neighbors(*)
	logical(1) onframe(nmax,nf,maxbands)

	dmagmax = -2.5*alog10(pratio*1.0)
	cd = cos(DEClist(n)*dtor)
	nblend = 1
	neighbors(nblend) = n
	if (maxval(eMAGSTD(n,1:4)) == eMagBad) return
	i = n

	do while (i < nsrc .and. nblend < nbmax)
	    i = i+1
	    if (any(onframe(i,1:nf,1:maxbands))) then
	        dra = abs(RAlist(i)-RAlist(n))
		if(dra .ge. 180)  dra = dra - 360.0
		rsq = (dra*cd*3600)**2 +
     *		      ((DEClist(i)-DEClist(n))*3600)**2
		dmag = MAGSTD(i,1:4) - MAGSTD(n,1:4)
		where(eMAGSTD(i,1:4) == eMagBad) dmag = 1.e35

		if (maxval(eMAGSTD(i,1:4)) > eMagBad .and.
     *              minval(dmag) <= dmagmax) then
c                 find the largest search radius for bands which both
c                 have non-null mags
		  maxblendist = maxval(blendist,mask=eMAGSTD(i,1:4) > eMagBad
     *                                               .and.
     *                                               eMAGSTD(n,1:4) > eMagBad
     *                                )
c                 maxblendist will be -(max. real*4 value) if the mask is
c                 all .false.
		  if (maxblendist .gt. 0 .and. rsq <= maxblendist**2) then
		    nblend = nblend + 1
		    neighbors(nblend) = i
		  endif
		endif
	    endif
	enddo

	return

	end

