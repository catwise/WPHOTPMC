c parse the input file to WPHOT

	subroutine pinput (infile,nin,nf,path,basename,wflag,smode)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*(*) infile,path(nin),basename(nin)
	character*72 s0
	character*500 string

	integer*2 wflag (nin,4)
	integer nin, nf

	logical smode

	open (unit=11,file=infile)
	
	nf = 0
	do 100 j=1,99999
	  read (11,'(a)',end=99) string

	  if (string(1:1).eq.'\') goto 100           ! ' (for syntaxhighlighting
	  if (string(1:1).eq.'|') goto 100

	  L = numchar (string)

	  if (L.lt.10) then
              print *,'===WARNING: frame info file has a short line: ', j
              goto 100
          endif

	  nf = nf + 1
c	  if (nf.ge.15) goto 99  ! TEMP

	  I = 1
	  call sfields (string,I,path(nf))
	  I = 2
	  call sfields (string,I,basename(nf))

	  I = 3
	  call sfields (string,I,s0)
	  read (s0,*) wflag(nf,1)
	  I = I + 1
	  call sfields (string,I,s0)
          read (s0,*) wflag(nf,2)
	  I = I + 1
	  call sfields (string,I,s0)
          read (s0,*) wflag(nf,3)
	  I = I + 1
	  call sfields (string,I,s0)
          read (s0,*) wflag(nf,4)

 100	continue
 99	close (11)

	if (smode) then  ! demand wflag conformance

	Jfr = 1
	ns = 0
	do ib=1,4
	  ns = ns + wflag(Jfr,ib)
	enddo

	if (ns.gt.1) then
		write (6,*) 'wflag = ',(wflag(Jfr,ib),ib=1,4)
		write (6,*) 'ERROR -- ifile inconsistent with single-band mode; exiting'
		call exit (9)
	endif

	endif


	return
	end
	  



