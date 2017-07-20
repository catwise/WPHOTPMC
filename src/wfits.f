C *************************************************************************
      subroutine wfits (nx,ny,lsize,larray,fout)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

C  Create a FITS primary array containing a 2-D image

      integer status,readwrite,inunit,outunit,blocksize,
     1   bitpix,naxis,naxes(2),nkeys,nspace
      integer i,j,group,fpixel,nelements
      real*4 larray(lsize)
      real*4 ra0,dec0,crpix10,crpix20,cd01,cd02,rot0
      character*80 record,comment
      character*(*) fout
      logical simple,extend,zexist,erase

C  The STATUS parameter must be initialized before using FITSIO.  A
C  positive value of STATUS is returned whenever a serious error occurs.
C  FITSIO uses an `inherited status' convention, which means that if a
C  subroutine is called with a positive input value of STATUS, then the
C  subroutine will exit immediately, preserving the status value. For 
C  simplicity, this program only checks the status value at the end of 
C  the program, but it is usually better practice to check the status 
C  value more frequently.

      status=0
	zexist = .false.
        erase = .true.
        call deletefile(fout,status,zexist, erase)




C  Get  unused Logical Unit Numbers to use to open the FITS files.

        call ftgiou(outunit,status)

	status=0

C  Create the new empty FITS file.  The blocksize parameter is a
C  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      call ftinit(outunit,fout,blocksize,status)


C  This do-loop of calls to FTGREC and FTPREC copies all the keywords from
C  the input to the output FITS file.  Notice that the specified number
C  of rows in the output table, as given by the NAXIS2 keyword, will be
C  incorrect.  This value will be modified later after it is known how many
C  rows will be in the table, so it does not matter how many rows are specified
C  initially.


C  Initialize parameters about the FITS image.
C  BITPIX = 16 means that the image pixels will consist of 16-bit
C  integers.  The size of the image is given by the NAXES values. 
C  The EXTEND = TRUE parameter indicates that the FITS file
C  may contain extensions following the primary array.
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.false.


c	call FTDKEY(outunit,'NAXIS1',status)
c         write (6,*) 'write status ',status
c        status=0
c        call ftpkyj(outunit,'NAXIS1',nx,'array size',status)
c         write (6,*) 'write status ',status
c        status=0
c        call FTDKEY(outunit,'NAXIS2',status)
c         write (6,*) 'write status ',status
c        status=0
c        call ftpkyj(outunit,'NAXIS2',ny,'array size',status)
c        write (6,*) 'write status ',status


C  Write the required header keywords to the file
       call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)

C  Write the array to the FITS file.
C  The last letter of the subroutine name defines the datatype of the
C  array argument; in this case the 'J' indicates that the array has an
C  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
C  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
C  total number of pixels.  GROUP is seldom used parameter that should
C  almost always be set = 1.
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftppre(outunit,group,fpixel,nelements,larray,status)



c	write (6,*) 'write status ',status

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
c	call ftclos(inunit, status)
c      call ftfiou(inunit, status)
      call ftclos(outunit, status)
      call ftfiou(outunit, status)

      return
      end




C *************************************************************************
      subroutine wfits_chimap (nx,ny,lsize,larray,fout, rmid,dmid,xscale,yscale)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

C  Create a FITS primary array containing a 2-D image

      integer status,readwrite,inunit,outunit,blocksize,
     1   bitpix,naxis,naxes(2),nkeys,nspace
      integer i,j,group,fpixel,nelements
      real*4 larray(lsize)
      real*4 ra0,dec0,crpix10,crpix20,cd01,cd02,rot0
      character*80 record,comment
      character*(*) fout
      logical simple,extend,zexist,erase

C  The STATUS parameter must be initialized before using FITSIO.  A
C  positive value of STATUS is returned whenever a serious error occurs.
C  FITSIO uses an `inherited status' convention, which means that if a
C  subroutine is called with a positive input value of STATUS, then the
C  subroutine will exit immediately, preserving the status value. For 
C  simplicity, this program only checks the status value at the end of 
C  the program, but it is usually better practice to check the status 
C  value more frequently.

      status=0
	zexist = .false.
        erase = .true.
        call deletefile(fout,status,zexist, erase)




C  Get  unused Logical Unit Numbers to use to open the FITS files.

        call ftgiou(outunit,status)

	status=0

C  Create the new empty FITS file.  The blocksize parameter is a
C  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      call ftinit(outunit,fout,blocksize,status)


C  This do-loop of calls to FTGREC and FTPREC copies all the keywords from
C  the input to the output FITS file.  Notice that the specified number
C  of rows in the output table, as given by the NAXIS2 keyword, will be
C  incorrect.  This value will be modified later after it is known how many
C  rows will be in the table, so it does not matter how many rows are specified
C  initially.


C  Initialize parameters about the FITS image.
C  BITPIX = 16 means that the image pixels will consist of 16-bit
C  integers.  The size of the image is given by the NAXES values. 
C  The EXTEND = TRUE parameter indicates that the FITS file
C  may contain extensions following the primary array.
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.false.


c	call FTDKEY(outunit,'NAXIS1',status)
c         write (6,*) 'write status ',status
c        status=0
c        call ftpkyj(outunit,'NAXIS1',nx,'array size',status)
c         write (6,*) 'write status ',status
c        status=0
c        call FTDKEY(outunit,'NAXIS2',status)
c         write (6,*) 'write status ',status
c        status=0
c        call ftpkyj(outunit,'NAXIS2',ny,'array size',status)
c        write (6,*) 'write status ',status


C  Write the required header keywords to the file
       call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)


c rmid,dmid,xscale,yscale
        crpix1 = nx/2.
        crpix2 = ny/2.

        status=0
        call FTDKEY(outunit,'CTYPE1',status)
        status=0
        call ftpkys(outunit, 'CTYPE1','ELON-TAN',
     1       'Orthographic Ecliptic Projection',status)

        call FTDKEY(outunit,'CTYPE2',status)
        status=0
        call ftpkys(outunit, 'CTYPE2','ELAT-TAN',
     1       'Orthographic Ecliptic Projection',status)

        call FTDKEY(outunit,'CRVAL1',status)
        status=0
        call ftpkyE(outunit,'CRVAL1',rmid,10,'ELON (deg)',status)

        call FTDKEY(outunit,'CRVAL2',status)
        status=0
        call ftpkye(outunit,'CRVAL2',dmid,10,'ELAT (deg)',status)

        call FTDKEY(outunit,'CRPIX1',status)
        status=0
        call ftpkye(outunit,'CRPIX1',crpix1,8,'pixel ref',status)

        call FTDKEY(outunit,'CRPIX2',status)
        status=0
        call ftpkye(outunit,'CRPIX2',crpix2,8,'pixel ref',status)

        call FTDKEY(outunit,'CDELT1',status)
        status=0
        call ftpkye(outunit,'CDELT1',xscale,10,'pixel scale',status)

        call FTDKEY(outunit,'CDELT2',status)
        status=0
        call ftpkye(outunit,'CDELT2',yscale,10,'pixel scale',status)

        call FTDKEY(outunit,'CROTA2',status)
        status=0
        call ftpkye(outunit,'CROTA2',0.,8,'rotation',status)

        call FTDKEY(outunit,'EQUINOX',status)
        status=0
	tt = 2000.
        call ftpkye(outunit,'EQUINOX',tt,8,'equinox',status)



C  Write the array to the FITS file.
C  The last letter of the subroutine name defines the datatype of the
C  array argument; in this case the 'J' indicates that the array has an
C  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
C  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
C  total number of pixels.  GROUP is seldom used parameter that should
C  almost always be set = 1.
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftppre(outunit,group,fpixel,nelements,larray,status)

	status=0
        call ftpkys(outunit,'COMMENT',
     1     ' Chi-square map','created by T.H. Jarrett,
     1      IPAC/Caltech',status)

c	write (6,*) 'write status ',status

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
c	call ftclos(inunit, status)
c      call ftfiou(inunit, status)
      call ftclos(outunit, status)
      call ftfiou(outunit, status)

      return
      end


