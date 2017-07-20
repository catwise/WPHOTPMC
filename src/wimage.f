C *************************************************************************
	subroutine wimage (nx,ny,lsize,larray,fin,fout)


C  Create a FITS primary array containing a 2-D image

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      integer status,readwrite,inunit,outunit,blocksize,
     1   bitpix,naxis,naxes(2),nkeys,nspace
      integer i,j,group,fpixel,nelements,ncombine
      real*4 larray(lsize)
      real*4 ra0,dec0,crpix10,crpix20,cd01,cd02,rot0
	real*4 sumut
      character*80 record,comment,s0
      character*(*) fin,fout
      logical simple,extend, zexist, erase


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
c	call deletefile(fout,status)

C  Get  unused Logical Unit Numbers to use to open the FITS files.
        call ftgiou(inunit,status)
        call ftgiou(outunit,status)

c	write (6,*) 'test0 ',status,outunit

C  The input FITS file is opened with READONLY access, and the output
C  FITS file is opened with WRITE access.
	readwrite=0
	blocksize=1
        call ftopen(inunit,fin,readwrite,blocksize,status)

c	write (6,*) 'test1 ',status

	status=0
c       readwrite=1
cblocksize=1
c       call ftopen(outunit,fout,readwrite,blocksize,status)


C  Create the new empty FITS file.  The blocksize parameter is a
C  historical artifact and the value is ignored by FITSIO.

	write (6,'(a)') fout(1:50)
      blocksize=1
      call ftinit(outunit,fout,blocksize,status)

c	write (6,*) 'test2 ',status

C  This do-loop of calls to FTGREC and FTPREC copies all the keywords from
C  the input to the output FITS file.  Notice that the specified number
C  of rows in the output table, as given by the NAXIS2 keyword, will be
C  incorrect.  This value will be modified later after it is known how many
C  rows will be in the table, so it does not matter how many rows are specified
C  initially.

C  Find the number of keywords in the input table header.
      call ftghsp(inunit,nkeys,nspace,status)

      do i=1,nkeys
          call ftgrec(inunit,i,record,status)
          call ftprec(outunit,record,status)
      end do

c	write (6,*) 'testH ',status


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

C  Write the required header keywords to the file
c      call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)

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

c	write (6,*) 'test3 ',status

C  Write another optional keyword to the header
C  The keyword record will look like this in the FITS file:
C
C  EXPOSURE=                 1500 / Total Exposure Time
C

c	call FTDKEY(outunit,'NAXIS1',status)
c        status=0
c	call ftpkyj(outunit,'NAXIS1',nx,'array size',status)
c	status=0
c	call FTDKEY(outunit,'NAXIS2',status)
c        status=0
c        call ftpkyj(outunit,'NAXIS2',ny,'array size',status)

c	write (6,*) 'write status ',status

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
	call ftclos(inunit, status)
      call ftfiou(inunit, status)
      call ftclos(outunit, status)
      call ftfiou(outunit, status)

      return
      end


C *************************************************************************
      subroutine deletefile(filename,status, zexist, erase)

C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename
	logical zexist, erase

	call access(filename,zexist,erase)

      
      return
      end


