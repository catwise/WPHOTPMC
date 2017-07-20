	subroutine readFhead(fname,Hdr)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

C  Print out all the header keywords in all extensions of a FITS file

	character*(*) fname,Hdr
	real*4 cdelt1, cdelt2
      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      integer L,numchar,nl,nh
      character record*80
      character*80 comment


C  The STATUS parameter must always be initialized.
      status=0
C     open the FITS file, with read-only access.  The returned BLOCKSIZE
C     parameter is obsolete and should be ignored.

	call ftgiou(unit,status)

	L = numchar(fname)
c	write (6,'(a)') ' '
c	write (6,'(a)') fname(1:L)

	status = 0
      readwrite=0
      call ftopen(unit,fname(1:L),readwrite,blocksize,status)

	if (status.ne.0) then
		write (6,*) 'problem reading header ',status
		write (6,'(a)') fname(1:L)
	endif


c	cdelt1 = 0.
c        cdelt2 = 0.
c        call ftgkye(unit, 'CDELT1', cdelt1, comment, status)
c        if (status.gt.0) status=0
c
c        call ftgkye(unit, 'CDELT2', cdelt2, comment, status)
c        if (status.gt.0) status=0

cwrite (6,*) status

C  The FTGHSP subroutine returns the number of existing keywords in the
C  current header data unit (CHDU), not counting the required END keyword,
      call ftghsp(unit,nkeys,nspace,status)

c	write (6,*) status,nkeys
	Hdr = ''

C  Read each 80-character keyword record, and print it out.
	nl = 1
      do i = 1, nkeys
          call ftgrec(unit,i,record,status)
c         write (6,'(a)') record

	  nh = nl + 79

	  if (nh.gt.150000) goto 47

c	write (6,*) nl,nh
	  Hdr (nl:nh) = record
	  nl = nh + 1

      end do

 47 	call ftclos(unit, status)
        call ftfiou(unit, status)


	return
	end

