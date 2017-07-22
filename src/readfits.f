      subroutine readimage 
     1     (nx,ny,lsize,array,iarray,fname, 
     1     crval1,crval2,
     1     cdelt1,cdelt2,crot,crpix1,crpix2, tJD)

C  Read a FITS image and determine the minimum and maximum pixel value.
C  Rather than reading the entire image in
C  at once (which could require a very large array), the image is read
C  in pieces, 100 pixels at a time.  

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      integer status,unit,readwrite,blocksize,naxes(2),nfound
      integer group,firstpix,nbuffer,npixels,i
      real*8 cd1a,cd1b,cd2a,cd2b,tDB,rat,angle, tJD, tJDmin, tJDmax
      real*4 datamin,datamax,nullval,array(lsize)
      integer*4 bitpix,iarray(lsize)
      real*4 crval1,crval2,cdelt1,cdelt2,crot,crpix1,crpix2, zJD
      logical anynull,debug
      character*(*) fname
      character*80 comment
	
	debug = .true.
	debug = .false.

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

	L = numchar(fname)
	if (debug) write (6,'(a,a)') 'reading ',fname(1:L)
	if (debug) write (6,*) status

C  Open the FITS file 
      readwrite=0
      call ftopen(unit,fname,readwrite,blocksize,status)
      if (status .ne. 0)then
          print *,'*** READIMAGE failed to open ',fname(1:L)
          call exit (9)
       end if

	if (debug)  write (6,*) status

C  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

	  if (debug)  write (6,*) status


C  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'*** READIMAGE failed to read the NAXISn keywords: ',fname(1:L)
          call exit (9)
       end if

	nx = naxes(1)
        ny = naxes(2)

c	write (6,*) nx,ny
	if (debug)  write (6,*) status,nx,ny


cccc

		!!  using MIN/MAX to compute JDOBS
	tJD = 0.
	tJDmin = 0.
	tJDmax = 0.

	call ftgkyd(unit, 'MJDMIN', tJDmin, comment, status)  ! min Julian Date
	call ftgkyd(unit, 'MJDMAX', tJDmax, comment, status)  ! max Julian Date

	tJD = (tJDmax+tJDmin)/2.  ! TJ 2017


c        call ftgkyd(unit, 'MJD_OBS', tJD, comment, status)  ! Julian Date

c	write (6,'(a,2f13.6)') 'JD ',tJD,gg+zz

        if (status.gt.0) then
                tJD = 0.
                status=0
        endif

	bitpix = 0
	call ftgkyj(unit, 'BITPIX', bitpix, comment, status)
	if (status.ne.0) then
          print *,'*** READIMAGE failed to read BITPIX: ',fname(1:L)
          call exit (9)
        end if

	crval1 = 0.
	call ftgkye(unit, 'CRVAL1', crval1, comment, status)
	if (status.gt.0) then
		crval1 = 0.
		status=0
	endif

	crval2 = 0.
        call ftgkye(unit, 'CRVAL2', crval2, comment, status)
        if (status.gt.0) then
                crval1 = 0.
                status=0
        endif
	
	call ftgkye(unit, 'CRPIX1', crpix1, comment, status)
        if (status.gt.0) status=0
	
	call ftgkye(unit, 'CRPIX2', crpix2, comment, status)
        if (status.gt.0) status=0

	cdelt1 = 0.
        cdelt2 = 0.
	call ftgkye(unit, 'CDELT1', cdelt1, comment, status)
        if (status.gt.0) status=0
        
        call ftgkye(unit, 'CDELT2', cdelt2, comment, status)
        if (status.gt.0) status=0

	call ftgkye(unit, 'CROTA2', crot, comment, status)
        if (status.gt.0) then
		crot = 0.
		status=0
	endif

c  CD matrix

	 if (debug)  write (6,*) 'here'

	if ((abs(cdelt1).le.1.e-5).and.(abs(cdelt2).le.1.e-5)) then

	  cd1_1=0.
	  cd1_2=0.
	  cd2_1=0.
	  cd2_2=0.

	 call ftgkye(unit, 'CD1_1', cd1_1, comment, status)
	 status=0

	 call ftgkye(unit, 'CD1_2', cd1_2, comment, status)
         status=0

	 call ftgkye(unit, 'CD2_1', cd2_1, comment, status)
         status=0
	
	 call ftgkye(unit, 'CD2_2', cd2_2, comment, status)
         status=0


         if (cd2_2.ne.0.) then
                  rat = CD1_2 / CD2_2
                  angle = -datan (rat) * 57.2957795d0

		  tdb = angle/57.2957795d0

	  	  cd2a = CD2_2 / dcos(tdb)
	  	  cd2b = -CD1_2 / dsin(tdb)

		if ((abs(CD2_2).gt.abs(CD1_2))) then
                      cdelt2 = cd2a*1.
                else
                      cdelt2 = cd2b*1.
                endif

	  	  cd1a = -CD1_1 / dcos(tdb)
	  	  cd1b = -CD2_1 / dsin(tdb)

		if ((abs(CD1_1).gt.abs(CD2_1))) then
                     cdelt1 = cd1a*1.
                else
                     cdelt1 = cd1b*1.
                endif

                  crot = angle*1.

	else if (cd1_1.ne.0.) then

		  	rat = CD2_1 / CD1_1
			angle = datan (rat) * 57.2957795d0
			tdb = angle/57.2957795d0

			cd2a = CD2_2 / dcos(tdb)
			cd2b = -CD1_2 / dsin(tdb)

			if ((abs(CD2_2).gt.abs(CD1_2))) then
			  cdelt2 = cd2a*1.
		 	else
			  cdelt2 = cd2b*1.
			endif

			cd1a = -CD1_1 / dcos(tdb)
			cd1b = -CD2_1 / dsin(tdb)
			if ((abs(CD1_1).gt.abs(CD2_1))) then
			 cdelt1 = cd1a*1.	
		        else
			 cdelt1 = cd1b*1.
			endif

			crot = angle*1.

	  endif

        endif


	 if (debug)  write (6,*) 'here2 ',unit,status

C  Initialize variables
      npixels=naxes(1)*naxes(2)
	if (debug)  write (6,*) npixels

      status = 0
      group=1
      firstpix=1
      nullval=-9991
      datamin=1.0E30
      datamax=-1.0E30
      nbuffer=npixels

      if(bitpix .lt. 0) then
        ! Floating point data
        call ftgpve(unit,group,firstpix,nbuffer,nullval,
     1              array,anynull,status)
      else
        call ftgpvj(unit,group,firstpix,nbuffer,nullval,
     1              iarray,anynull,status)
      endif

	if (status.ne.0) then
          print *,'*** READIMAGE failed to read ',nbuffer,' pixels: ',fname(1:L)
          call exit (9)
        end if


      if (debug)  write (6,*) status


C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

       return
      end

