c    coordinate conversion routines for ao's & coadded survey data images
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	--  output from test runs are included below code  --
c
c	for questions contact: internet-	gene@ipac.caltech.edu
c				span   -	romeo::"gene%ipac" 
c				bitnet -	gene%ipac@Hamlet.Bitnet
c
c				or (818) 584-2932  
c					gene kopan  ipac
c
c -------------------------------------------------------------------------
c	test driver
c
	subroutine find (rag,decg,icrp,jcrp,crota2,cra,cdec,
     1    nz,ny,ra,dec,pz,py)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	double precision te2g(3,3),ppy,ppz,ddry,ddrz,
     1    rra,ddec,rrag,ddecg,ccrota2
	real*4 icrp,jcrp

	dry = cdec * 60.
	drz = cra * 60.

        ppy = py * 1.d0
        ppz = pz * 1.d0
        ddry = dry * 1.d0
        ddrz = drz * 1.d0
        RRA = ra*1.0d0
        DDEC = dec *1.d0
        rrag = rag*1.d0
        ddecg = decg*1.d0
        ccrota2 = crota2 * 1.d0

	call cte2g( rrag,ddecg,ccrota2, te2g )
c	call Single_CTE2G(rag,decg,CROTA2,TE2G)

	call skytog(icrp,jcrp,ppy,ppz,ddry,ddrz,te2g,rra,ddec)
	py = sngl (ppy)
	pz = sngl (ppz)

c	call Singel_SKYTOG(PY,PZ,NY,NZ,DRY,DRZ,TE2G,ra,dec)
c	write (6,*) 'fin ',pz,py


	return
	end

c -----------------------------------------------------------------------
      SUBROUTINE Singel_GTOSKY(PY,PZ,NY,NZ,DRY,DRZ,TE2G,ALP,DEL)
C
C     GTOSKY - CONVERT PIXEL (PY,PZ) TO RA & DEC (ALP,DEL)
C
C     PY     =  Y PIXEL
C     PZ     =  Z PIXEL
C     NY     =  NUMBER OF CELLS IN Y DIRECTION
C     NZ     =  NUMBER OF CELLS IN Z DIRECTION
C     DRY    =  CELL SIZE IN Y DIRECTION (ARCMIN)
C     DRZ    =  CELL SIZE IN Z DIRECTION (ARCMIN)
C     TE2G   =  TRANSFORMATION FROM EME TO GRID COORDINATES
C     ALP    =  RIGHT ASCENSION OF PIXEL PY,PZ (DEG, EME50)
C     DEL    =  DECLINATION     OF PIXEL PY,PZ (DEG, EME50)
C
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      INTEGER NY,NZ,NYC,NZC
      REAL*4    ALP,DEL,PY,PZ
      DOUBLE PRECISION  TE2G(3,3),VE(3),VT(3)
      DATA RTD/57.29577951/
C                        COMPUTE GRID VECTOR
      NYC = NY/2 +1
      NZC = NZ/2 +1
      VT(2) = (PY-NYC)*DRY/(60.*RTD)
      VT(3) = (PZ-NZC)*DRZ/(60.*RTD)
      VT2   = 1.-VT(2)*VT(2)-VT(3)*VT(3)
      VT(1) = SQRT( VT2 )
C                         COMPUTE EME VECTOR AND RA,DEC
      VE(1) = TE2G(1,1)*VT(1)+TE2G(2,1)*VT(2)+TE2G(3,1)*VT(3)
      VE(2) = TE2G(1,2)*VT(1)+TE2G(2,2)*VT(2)+TE2G(3,2)*VT(3)
      VE(3) = TE2G(1,3)*VT(1)+TE2G(2,3)*VT(2)+TE2G(3,3)*VT(3)
      CDEL =  VE(2)*VE(2) + VE(3)*VE(3)
      CDEL = SQRT( CDEL )
      SDEL =  VE(1)
      DEL = RTD*ATAN2(  SDEL, CDEL  )
      SALP = -VE(2)
      CALP =  VE(3)
      ALP = RTD*ATAN2(  SALP, CALP  )
      IF( ALP .LT. 0. ) ALP = ALP + 360.
C
      RETURN
      END
C
      
      SUBROUTINE Singel_SKYTOG(PY,PZ,NY,NZ,DRY,DRZ,TE2G,ALP,DEL)
C
C      SKYTOG - CONVERT RA & DEC (ALP,DEL) TO PIXEL (PY,PZ)
C
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      INTEGER NY,NZ,NYC,NZC
      REAL*4    PY,PZ,ALP,DEL
      DOUBLE PRECISION  TE2G(3,3),VE(3),VT(3)
      DATA RTD/57.29577951/
C                           COMPUTE EME VECTOR
C
      NYC = NY/2 +1
      NZC = NZ/2 +1
      SDEL = SIN( sngl(DEL/ RTD*1.d0 ))
      CDEL = COS( sngl(DEL/ RTD*1.d0 ))
      SALP = SIN( sngl(ALP/ RTD*1.d0 ))
      CALP = COS( sngl(ALP/ RTD*1.d0 ))
      VE(1) = +SDEL
      VE(2) = -SALP*CDEL
      VE(3) = +CALP*CDEL
C                           COMPUTE GRID VECTOR AND PY,PZ
      VT(1) = TE2G(1,1)*VE(1)+TE2G(1,2)*VE(2)+TE2G(1,3)*VE(3)
      VT(2) = TE2G(2,1)*VE(1)+TE2G(2,2)*VE(2)+TE2G(2,3)*VE(3)
      VT(3) = TE2G(3,1)*VE(1)+TE2G(3,2)*VE(2)+TE2G(3,3)*VE(3)
      PY  = 60.*RTD*VT(2)/DRY + NYC
      PZ  = 60.*RTD*VT(3)/DRZ + NZC
C
      RETURN
      END
C
      SUBROUTINE Single_CTE2G(ALPG,DELG,CROTA2,TE2G)
C
C    COMPUTE EME50 TO GRID TRANSFORMATION
C
C     TE2G   =  TRANSFORMATION FROM EME TO GRID COORDINATES
C     ALPG   =  RIGHT ASCENSION OF GRID CENTER (DEG, EME50)
C     DELG   =  DECLINATION     OF GRID CENTER (DEG, EME50)
C     CROTA2 =  GRID ROTATION ANGLE  - FITS -  (DEG, EME50)
C
C	NOTE: TWSG IS CLOCKWISE TWIST WRT SOUTH
C
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DOUBLE PRECISION  TE2G(3,3)
      DATA RTD/57.29577951/
C
      TWSG = CROTA2 - 180.0
      SALP = SIN( ALPG/RTD )
      CALP = COS( ALPG/RTD )
      SDEL = SIN( DELG/RTD )
      CDEL = COS( DELG/RTD )
      STWS = SIN( TWSG/RTD )
      CTWS = COS( TWSG/RTD )
      TE2G(1,1) = +SDEL
      TE2G(1,2) = -SALP*CDEL
      TE2G(1,3) = +CALP*CDEL
      TE2G(2,1) = -CDEL*CTWS
      TE2G(2,2) = -CALP*STWS -SALP*SDEL*CTWS
      TE2G(2,3) = -SALP*STWS +CALP*SDEL*CTWS
      TE2G(3,1) = +CDEL*STWS
      TE2G(3,2) = -CALP*CTWS +SALP*SDEL*STWS
      TE2G(3,3) = -SALP*CTWS -CALP*SDEL*STWS
C
      RETURN
      END

	subroutine astgaltoeq (lii,bii,ra,dec)
c AST_GALTOEQ -- Convert galactic coordinates (1950) to equatorial coordinates.

	real*8 lii,bii,ra,dec, LP, BP, LO, BO, GEPOCH
	real*8 ao, ap, a1, b1, a2, b2, degtorad, radtodeg


c Definition of system: Longtitude of pole,Latitude of pole,Longitude of origin
c			Latitude of origin, Epoch of definition

	LP= 	123.00d0		
	BP=	27.40d0		
	LO=	97.7422d0
	BO=	-60.1810d0	
	GEPOCH=	1950.0d0

c   g187.1632+41.5656
c   d135.35+36.50


c  double	lii		# Galactic longitude (degrees)
c  double	bii		# Galactic latitude (degrees)
c  double	ra		# Right ascension (hours)
c  double	dec		# Declination (degrees)
c  double	epoch		# Epoch of coordinates


	ao = DEGTORAD (LO)
	bo = DEGTORAD (BO)
	ap = DEGTORAD (LP)
	bp = DEGTORAD (BP)
	a1 = DEGTORAD (lii)
	b1 = DEGTORAD (bii)

	call ast_coord (ao, bo, ap, bp, a1, b1, a2, b2)


	a2 = mod (24.0d0 + RADTODEG(a2) / 15.0d0, 24.0d0)
	b2 = RADTODEG (b2)

c	# Precess the coordinates
c	call ast_precess (a2, b2, GEPOCH, ra, dec, epoch)

	ra = a2 * 15.
	dec = b2


	return
	end


c AST_COORD -- Convert spherical coordinates to new system.
c
c This procedure converts the longitude-latitude coordinates (a1, b1)
c of a point on a sphere into corresponding coordinates (a2, b2) in a
c different coordinate system that is specified by the coordinates of its
c origin (ao, bo).  The range of a2 will be from -pi to pi.

	subroutine ast_coord (ao, bo, ap, bp, a1, b1, a2, b2)

	real*8 ao, bo,ap, bp,a1, b1,a2, b2
	real*8 sao, cao, sbo, cbo, sbp, cbp, x, y, z, xp, yp, zp, temp

c	ao, bo		# Origin of new coordinates (radians)
c	ap, bp		# Pole of new coordinates (radians)
c	a1, b1		# Coordinates to be converted (radians)
c	a2, b2		# Converted coordinates (radians)


	x = cos (a1) * cos (b1)
	y = sin (a1) * cos (b1)
	z = sin (b1)
	xp = cos (ap) * cos (bp)
	yp = sin (ap) * cos (bp)
	zp = sin (bp)

c Rotate the origin about z.
	sao = sin (ao)
	cao = cos (ao)
	sbo = sin (bo)
	cbo = cos (bo)
	temp = -xp * sao + yp * cao
	xp = xp * cao + yp * sao
	yp = temp
	temp = -x * sao + y * cao
	x = x * cao + y * sao
	y = temp

c Rotate the origin about y.
	temp = -xp * sbo + zp * cbo
	xp = xp * cbo + zp * sbo
	zp = temp
	temp = -x * sbo + z * cbo
	x = x * cbo + z * sbo
	z = temp

c Rotate pole around x.
	sbp = zp
	cbp = yp
	temp = y * cbp + z * sbp
	y = y * sbp - z * cbp
	z = temp

c Final angular coordinates.
	a2 = atan2 (y, x)
	b2 = asin (z)

	return
	end

	real*8 function DEGTORAD (x)
	real*8 x

	DEGTORAD = x / 57.29577951d0

	return
	end

	real*8 function RADTODEG (x)
	real*8 x

	RADTODEG = x * 57.29577951d0

	return
	end


	subroutine coord_pos (x,y,ra,dec,icrp,jcrp,
     1    rot,cra,cdec,nx,ny,raf,decf)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 icrp,jcrp

c *** note :  changed -cra to cra

	if ((raf.eq.0.).and.(decf.eq.0)) then
	  call find2 (ra,dec,icrp,jcrp,rot,cra,cdec,nx,ny,raf,decf,x,y)
	else
	  call find (ra,dec,icrp,jcrp,rot,cra,cdec,nx,ny,raf,decf,x,y)
	endif

c	  write (6,12) raf,decf,x,y
c 12	  format (2f10.4,2x,2i7)


	return	
	end



c    coordinate conversion routines for ao's & coadded survey data images
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	--  output from test runs are included below code  --
c
c	for questions contact: internet-	gene@ipac.caltech.edu
c				span   -	romeo::"gene%ipac" 
c				bitnet -	gene%ipac@Hamlet.Bitnet
c
c				or (818) 584-2932  
c					gene kopan  ipac
c
c -------------------------------------------------------------------------
c	test driver
c
	subroutine find2(rag,decg,icrp,jcrp,
     1    crota2,cra,cdec,nz,ny,ra,dec,pz,py)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	double precision te2g(3,3),ppy,ppz,ddry,
     1    ddrz,rra,ddec,rrag,ddecg,ccrota2
	real*4 icrp,jcrp

	dry = cdec * 60.
	drz = cra * 60.

	ppy = py * 1.d0
	ppz = pz * 1.d0
	ddry = dry * 1.d0
	ddrz = drz * 1.d0
	RRA = ra*1.0d0
	DDEC = dec *1.d0
	rrag = rag*1.d0
	ddecg = decg*1.d0
	ccrota2 = crota2 * 1.d0

	call cte2g( rrag,ddecg,ccrota2,te2g )

	call gtosky( icrp,jcrp,ppy,ppz,ny,nz,ddry,ddrz,te2g,rra,ddec)
	
	ra = sngl (rra)
	dec = sngl (ddec)

c     SUBROUTINE GTOSKY(icrp,jcrp,PY,PZ,NY,NZ,DRY,DRZ,TE2G,ALP,DEL)
C
c    INTEGER NY,NZ,NYC,NZC
c     real*4 icrp,jcrp
c     DOUBLE PRECISION  ALP,DEL,PY,PZ,DRY,DRZ,RTD
c     DOUBLE PRECISION  DATAN2,VT2,SDEL,CDEL,SALP,CALP
c     DOUBLE PRECISION  TE2G(3,3),VE(3),VT(3)

		 
	return
	end


	subroutine asteqtogal (ras,decs,glong,glat)
c -- Convert  equatorial coordinates (1950)  to galactic coordinates 
c	ras,decs in degreees

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 ras,decs,glong,glat
	real*8 lii,bii,ra,dec
	real*8 x,y,z,r,temp,drar,cosdecg,sindecg,cosdecr, degtorad, radtodeg
	real*8 LONGNCP,RAGPOLE,DECGPOLE,GEPOCH

c Definition of system: Longtitude of pole,Latitude of pole,Longitude of origin
c			Latitude of origin, Epoch of definition

	LONGNCP = 123.00d0
	RAGPOLE  = 192.25d0      
	DECGPOLE  = 27.4d0       
	GEPOCH    =  1950.0d0    

	ra = ras / 15.d0
	dec = decs * 1.d0

c Precompute the necessary constants.
	drar = DEGTORAD (15.0d0 * ra - RAGPOLE)
	cosdecg = cos (DEGTORAD (DECGPOLE))
	sindecg = sin (DEGTORAD(DECGPOLE))
	cosdecr = cos (DEGTORAD (dec))

c Compute the tansformation equations
	x = cosdecr * cos (drar)
	y =  cosdecr * sin (drar)
	z = sin (DEGTORAD (dec))
	temp = z * cosdecg - x * sindecg
	z = z * sindecg + x * cosdecg
	x = temp
        r = sqrt (x * x + y * y)
c 
c Compute lii and bii and convert to degrees.
	if (r .lt. 0.0001) then
            lii = 0.0d0
        else
            lii = DEGTORAD (LONGNCP)  + atan2 (-y, x)
	endif

        if (lii .lt. 0.0d0) lii = lii + (2. * 3.141593)

        bii = 1.d0 * atan2 (z, r)
        lii = RADTODEG (lii)
        bii = RADTODEG (1.d0 * bii)

c   g187.1632+41.5656
c   d135.35+36.50

	glong = sngl(lii)
	glat = sngl(bii)

	return	
	end


c AST_PRECESS -- Precess coordinates from epoch1 to epoch2.
c
c The method used here is based on the new IAU system described in the
c supplement to the 1984 Astronomical Almanac.  The precession is
c done in two steps; precess epoch1 to the standard epoch J2000.0 and then
c precess from the standard epoch to epoch2.  The precession between
c any two dates is done this way because the rotation matrix coefficients
c are given relative to the standard epoch.

	subroutine ast_precess (ra, dec, ep1, raout, decout, ep2)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 ra, dec, ep1, raout, decout, ep2
	real*8 ra1, dec1, epoch1, ra2, dec2, epoch2
	real*8 r0(3), r1(3), p(3, 3)
	real*8  DEGTORAD, RADTODEG 

c double	r0[3], r1[3], p[3, 3]
c bool	fp_equald()

	ra1 = ra/15.d0
	dec1 = dec*1.d0
	epoch1 = ep1*1.d0
	epoch2  = ep2*1.d0

c	# If the input epoch is 0 or undefined then assume the input epoch
c	# is the same as the output epoch.  If the two epochs are the same
c	# then return the coordinates from epoch1.

	if ((epoch1.eq.0.) .or. (epoch1.eq.epoch2)) then
	    ra2 = ra1
	    dec2 = dec1
   	    goto 99
	endif


c	# Rectangular equitorial coordinates (direction cosines).
	ra2 = DEGTORAD (ra1 * 15.)
	dec2 = DEGTORAD (dec1)

	r0(1) = cos (ra2) * cos (dec2)
	r0(2) = sin (ra2) * cos (dec2)
	r0(3) = sin (dec2)

c	# If epoch1 is not the standard epoch then precess to the standard
c	# epoch.

	if (epoch1 .ne. 2000.) then
	    call ast_rotmatrix (epoch1, p)

c	    # Note that we multiply by the inverse of p which is the
c	    # transpose of p.

	    r1(1) = p(1, 1) * r0(1) + p(1, 2) * r0(2) + p(1, 3) * r0(3)
	    r1(2) = p(2, 1) * r0(1) + p(2, 2) * r0(2) + p(2, 3) * r0(3)
	    r1(3) = p(3, 1) * r0(1) + p(3, 2) * r0(2) + p(3, 3) * r0(3)
	    r0(1) = r1(1)
	    r0(2) = r1(2)
	    r0(3) = r1(3)

	endif


c	# If epoch2 is not the standard epoch then precess from the standard
c	# epoch to the desired epoch.

	if (epoch2 .ne. 2000.) then

	    call ast_rotmatrix (epoch2, p)
	    r1(1) = p(1, 1) * r0(1) + p(2, 1) * r0(2) + p(3, 1) * r0(3)
	    r1(2) = p(1, 2) * r0(1) + p(2, 2) * r0(2) + p(3, 2) * r0(3)
	    r1(3) = p(1, 3) * r0(1) + p(2, 3) * r0(2) + p(3, 3) * r0(3)

	    r0(1) = r1(1)
            r0(2) = r1(2)
            r0(3) = r1(3)
 
	endif


c	# Convert from radians to hours and degrees.
	ra2 = RADTODEG (atan2 (r0(2), r0(1)) / 15.d0) 
	dec2 = RADTODEG (1.d0 * asin (r0(3)))
	if (ra2 .lt. 0.) ra2 = ra2 + 24
 99	i=1

	ra2 = ra2 * 15.d0
	raout = sngl(ra2) 
	decout = sngl(dec2)

	return
	end





c ROTMATRIX -- Compute the precession rotation matrix from the standard epoch
c J2000.0 to the specified epoch.

	subroutine ast_rotmatrix (epoch, p)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*8 epoch,p(3,3)
	real*8 t, a, b, c, ca, cb, cc, sa, sb, sc
	real*8  DEGTORAD, ast_julday


c double	ast_julday()

c	# The rotation matrix coefficients are polynomials in time measured
c	# in Julian centuries from the standard epoch.  The coefficients are
c	# in degrees.

	t = (ast_julday (epoch) - 2451545.0d0) / 36525d0

	a = t * (0.6406161d0 + t * (0.0000839d0 + t * 0.0000050d0))
	b = t * (0.6406161d0 + t * (0.0003041d0 + t * 0.0000051d0))
	c = t * (0.5567530d0 - t * (0.0001185d0 + t * 0.0000116d0))

c	# Compute the cosines and sines once for efficiency.
	ca = cos (DEGTORAD (a))
	sa = sin (DEGTORAD (a))
	cb = cos (DEGTORAD (b))
	sb = sin (DEGTORAD (b))
	cc = cos (DEGTORAD (c))
	sc = sin (DEGTORAD (c))

c	# Compute the rotation matrix from the sines and cosines.
	p(1, 1) = ca * cb * cc - sa * sb
	p(2, 1) = -sa * cb * cc - ca * sb
	p(3, 1) = -cb * sc
	p(1, 2) = ca * sb * cc + sa * cb
	p(2, 2) = -sa * sb * cc + ca * cb
	p(3, 2) = -sb * sc
	p(1, 3) = ca * sc
	p(2, 3) = -sa * sc
	p(3, 3) = cc

	return
	end


c AST_JULDAY -- Convert epoch to Julian day.

	real*8 function ast_julday (epoch)
	real*8 epoch,jd, J2000, JD2000 , JYEAR

c define  J2000           2000.0D0                # J2000
c define  JD2000          2451545.0D0             # J2000 Julian Date
c define  JYEAR           365.25D0                # Julian year

	J2000 = 2000.0D0 
	JD2000 = 2451545.0D0 
	JYEAR  = 365.25D0  


        jd = JD2000 + (epoch - J2000) * JYEAR
	ast_julday = jd

	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc
c double precision routines



c -----------------------------------------------------------------------
      SUBROUTINE GTOSKY(icrp,jcrp,PY,PZ,NY,NZ,DRY,DRZ,TE2G,ALP,DEL)

C
C     GTOSKY - CONVERT PIXEL (PY,PZ) TO RA & DEC (ALP,DEL)
C
C     PY     =  Y PIXEL
C     PZ     =  Z PIXEL
C     NY     =  NUMBER OF CELLS IN Y DIRECTION
C     NZ     =  NUMBER OF CELLS IN Z DIRECTION
C     DRY    =  CELL SIZE IN Y DIRECTION (ARCMIN)
C     DRZ    =  CELL SIZE IN Z DIRECTION (ARCMIN)
C     TE2G   =  TRANSFORMATION FROM EME TO GRID COORDINATES
C     ALP    =  RIGHT ASCENSION OF PIXEL PY,PZ (DEG, EME50)
C     DEL    =  DECLINATION     OF PIXEL PY,PZ (DEG, EME50)
C
      INTEGER NY,NZ
      real*4 icrp,jcrp
      DOUBLE PRECISION fyc,fzc
      DOUBLE PRECISION  ALP,DEL,PY,PZ,DRY,DRZ,RTD
      DOUBLE PRECISION  DATAN2,VT2,SDEL,CDEL,SALP,CALP
      DOUBLE PRECISION  TE2G(3,3),VE(3),VT(3)
      DATA RTD/57.29577951D0/
C                        COMPUTE GRID VECTOR

C  --- better to use FZC=CRPIX1 and FYC=CRPIX2 from fits header

      FYC = jcrp*1.d0
      FZC = icrp*1.d0
      VT(2) = (PY-FYC)*DRY/(60.*RTD)
      VT(3) = (PZ-FZC)*DRZ/(60.*RTD)
      VT2   = 1.-VT(2)*VT(2)-VT(3)*VT(3)
      VT(1) = DSQRT( VT2 )
C                         COMPUTE EME VECTOR AND RA,DEC
      VE(1) = TE2G(1,1)*VT(1)+TE2G(2,1)*VT(2)+TE2G(3,1)*VT(3)
      VE(2) = TE2G(1,2)*VT(1)+TE2G(2,2)*VT(2)+TE2G(3,2)*VT(3)
      VE(3) = TE2G(1,3)*VT(1)+TE2G(2,3)*VT(2)+TE2G(3,3)*VT(3)
      CDEL =  VE(2)*VE(2) + VE(3)*VE(3)
      CDEL =  DSQRT( CDEL )
      SDEL =  VE(1)
      DEL = RTD*DATAN2(  SDEL, CDEL  )
      SALP = -VE(2)
      CALP =  VE(3)
      ALP = RTD*DATAN2(  SALP, CALP  )
      IF( ALP .LT. 0. ) ALP = ALP + 360.
C
      RETURN
      END
Ccccccccccc
ccccccccccc

      SUBROUTINE SKYTOG(icrp,jcrp,PY,PZ,DRY,DRZ,TE2G,ALP,DEL)

C
C      SKYTOG - CONVERT RA & DEC (ALP,DEL) TO PIXEL (PY,PZ)
C
      real*4 icrp,jcrp
      DOUBLE PRECISION   ALP,DEL,PY,PZ,RTD,FYC,FZC
      DOUBLE PRECISION   SDEL,CDEL,SALP,CALP,DSIN,DCOS    
      DOUBLE PRECISION  TE2G(3,3),VE(3),VT(3),DRY,DRZ
      DATA RTD/57.29577951/
C                           COMPUTE EME VECTOR
C
C  --- better to use FZC=CRPIX1 and FYC=CRPIX2 from fits header

      FYC = jcrp*1.d0
      FZC = icrp*1.d0
      SDEL = DSIN( DEL/RTD )
      CDEL = DCOS( DEL/RTD )
      SALP = DSIN( ALP/RTD )
      CALP = DCOS( ALP/RTD )
      VE(1) = +SDEL
      VE(2) = -SALP*CDEL
      VE(3) = +CALP*CDEL
C                           COMPUTE GRID VECTOR AND PY,PZ
      VT(1) = TE2G(1,1)*VE(1)+TE2G(1,2)*VE(2)+TE2G(1,3)*VE(3)
      VT(2) = TE2G(2,1)*VE(1)+TE2G(2,2)*VE(2)+TE2G(2,3)*VE(3)
      VT(3) = TE2G(3,1)*VE(1)+TE2G(3,2)*VE(2)+TE2G(3,3)*VE(3)
      PY  = 60.d0*RTD*VT(2)/DRY + FYC
      PZ  = 60.d0*RTD*VT(3)/DRZ + FZC
C
      RETURN
      END
C
C
      SUBROUTINE CTE2G(ALPG,DELG,CROTA2,TE2G)
C
C    COMPUTE EME50 TO GRID TRANSFORMATION
C
C     TE2G   =  TRANSFORMATION FROM EME TO GRID COORDINATES
C     ALPG   =  RIGHT ASCENSION OF GRID CENTER (DEG, EME50)
C     DELG   =  DECLINATION     OF GRID CENTER (DEG, EME50)
C     CROTA2 =  GRID ROTATION ANGLE  - FITS -  (DEG, EME50)
C
C	NOTE: TWSG IS RIGHT-HANDED TWIST WRT SOUTH
C
      DOUBLE PRECISION  ALPG,DELG,CROTA2,RTD,TWSG,DSIN,DCOS
      DOUBLE PRECISION  SALP,CALP,SDEL,CDEL,STWS,CTWS
      DOUBLE PRECISION  TE2G(3,3)
      DATA RTD/57.29577951/
C
      TWSG = 180.0 + CROTA2 
      SALP = DSIN( ALPG/RTD )
      CALP = DCOS( ALPG/RTD )
      SDEL = DSIN( DELG/RTD )
      CDEL = DCOS( DELG/RTD )
      STWS = DSIN( TWSG/RTD )
      CTWS = DCOS( TWSG/RTD )
      TE2G(1,1) = +SDEL
      TE2G(1,2) = -SALP*CDEL
      TE2G(1,3) = +CALP*CDEL
      TE2G(2,1) = -CDEL*CTWS
      TE2G(2,2) = -CALP*STWS -SALP*SDEL*CTWS
      TE2G(2,3) = -SALP*STWS +CALP*SDEL*CTWS
      TE2G(3,1) = +CDEL*STWS
      TE2G(3,2) = -CALP*CTWS +SALP*SDEL*STWS
      TE2G(3,3) = -SALP*CTWS -CALP*SDEL*STWS
C
      RETURN
      END


cccc ecliptic coordinates; written by JWF

      subroutine Cel2Ec(RA, Dec, Long, Lat)
c
      real*4 RA, Dec, Long, Lat, SOb, Cob, X, Y, Z, d2r,
     +       cRA, cDec, sRA, sDec, X2, Y2
c
      data d2r/1.745329252e-2/, cOb, sOb/0.91748206, 0.39777716/
c
c-----------------------------------------------------------------------
c
      cRA   = cos(d2r*RA)
      cDec  = cos(d2r*Dec)
      sRA   = sin(d2r*RA)
      sDec  = sin(d2r*Dec)
c
      X =  sDec
      Y = -cDec*sRA
      Z =  cDec*cRA
c
      X2 =  X*Cob + Y*Sob
      Y2 = -X*Sob + Y*Cob
c     Z2 =  Z
c
      Lat  = asin(X2)/d2r
      Long = atan2(-Y2,Z)/d2r
      if (Long .lt. 0.0) Long = Long + 360.0
c
      return
      end

