	subroutine parg (narg,arg,level,ofile,intword,uncword,covword,mdexword,
     1     ifile,namlist,metaf, singlemode,verbose,pointless,unitest,
     1     qdir,psfdir,calbname,calgridX,calgridY,
     1     mdexname, coname, clev, imID, xscf, xscl, IOCname, detsnr,
     1     meproot, doreg )                              ! CJG B30228, TPC

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*(*) arg(narg),ofile,intword,uncword,covword,mdexword
	character*(*) ifile,namlist,level,metaf
	character*(*) psfdir,calbname
	character*(*) qdir, mdexname, coname, clev, imID
	character*(*) xscf, xscl, IOCname
	character*(*) meproot                    ! CJG B30228

	character*25 key,s0
	integer nkey,narg,numchar,L
	integer calgridX,calgridY
	logical singlemode,verbose,unitest, pointless
	integer doreg ! TPC

	key = '-level'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          level = arg(nkey+1)(1:L)
	else
		write (6,*) 'ERROR -- level not set ; exiting ...'
		call exit(9)
        endif


	key = '-ofile'

	call findkey (narg,arg,key,nkey)
	if (nkey.gt.0) then
	  L = numchar(arg(nkey+1))
	  ofile = arg(nkey+1)(1:L)
	else
                write (6,*) 'ERROR -- ofile not set ; exiting ...'
                call exit(9)
	endif

	key = '-meta'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          metaf = arg(nkey+1)(1:L)
	else
                write (6,*) 'ERROR -- meta not set ; exiting ...'
                call exit(9)
        endif


c	key = '-opath'
c        call findkey (narg,arg,key,nkey)
c        if (nkey.gt.0) then
c          L = numchar(arg(nkey+1))
c          opath = arg(nkey+1)(1:L)
c        endif

	key = '-int'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          intword = arg(nkey+1)(1:L)
        endif

	key = '-unc'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          uncword = arg(nkey+1)(1:L)
        endif

	key = '-cov'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          covword = arg(nkey+1)(1:L)
        endif

        key = '-mdex'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          mdexword = arg(nkey+1)(1:L)
        endif


	key = '-ifile'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          ifile = arg(nkey+1)(1:L)
	else
                write (6,*) 'ERROR -- ifile not set ; exiting ...'
                call exit(9)
        endif


	key = '-namlis'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          namlist = arg(nkey+1)(1:L)
	else
                write (6,*) 'ERROR -- namelist not set ; exiting ...'
                call exit(9)

        endif

	singlemode = .false.
	nkey = 0
	key = '-single'
	call findkey (narg,arg,key,nkey)

	if (nkey.gt.0) singlemode = .true.

	key = '-qadir'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          qdir = arg(nkey+1)(1:L)
        else
                write (6,*) 'ERROR -- qdir not set ; exiting ...'
                call exit(9)

        endif


c	key = '-calapcorrdir'
c        call findkey (narg,arg,key,nkey)
c        if (nkey.gt.0) then
c          L = numchar(arg(nkey+1))
c         apcorrdir = arg(nkey+1)(1:L)
c	else
c                write (6,*) 'ERROR -- calapcorrdir not set ; exiting ...'
c                call exit(9)
c
c        endif

	key = '-calpsfdir'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
         psfdir = arg(nkey+1)(1:L)
	else
                write (6,*) 'ERROR -- calpsfdir not set ; exiting ...'
                call exit(9)

        endif


	key = '-calbname'

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
         calbname = arg(nkey+1)(1:L)
	else
                write (6,*) 'ERROR -- calbname not set ; exiting ...'
                call exit(9)

        endif

	key = '-calgridszX'
	calgridX = 0

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
         s0 = arg(nkey+1)(1:L)
	 read (s0,*) calgridX
	else
                write (6,*) 'ERROR -- X gridsize not set ; exiting ...'
                call exit(9)

        endif

	key = '-calgridszY'
	calgridY = 0

        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
         s0 = arg(nkey+1)(1:L)
         read (s0,*) calgridY
	else
                write (6,*) 'ERROR -- Y gridsize not set ; exiting ...'
                call exit(9)

        endif


	key = '-mdettab'
	mdexname = 'null'
	call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          mdexname = arg(nkey+1)(1:L)
	endif

	key = '-coadd'
        coname = 'null'
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          coname = arg(nkey+1)(1:L)
        endif

	key = '-clevel'
        clev = 'null'
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          clev = arg(nkey+1)(1:L)
        endif

	key = '-imageid'
        imID = ''
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          imID = arg(nkey+1)(1:L)
        endif


	key = '-xscfile'
        xscf = ''
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          xscf = arg(nkey+1)(1:L)
        endif

	key = '-xsclookup'
        xscl = ''
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          xscl = arg(nkey+1)(1:L)
        endif

	nkey = 0
        key = '-pointless'
	pointless = .false.
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0)  pointless = .true.

	detsnr = 20.0  ! default value
	nkey = 0
        key = '-detsnr'
        call findkey (narg,arg,key,nkey)
	if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
	  s0 = arg(nkey+1)(1:L)
          read (s0,*) detsnr 
        endif


	key = '-ioc'
        IOCname = ''
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
          IOCname = arg(nkey+1)(1:L)
        endif


c==========================start of code added by CJG B30213===========

	key = '-meproot'
	meproot = ''
	call findkey(narg,arg,key,nkey)
	if(nkey.gt.0) then 
	   L = numchar(arg(nkey+1))
	   meproot = arg(nkey+1)(1:L)
	endif

c==========================end of code added by CJG B30213===========

c       !!!! Added by TPC
	doreg = 0  ! default value
	nkey = 0
        key = '-doreg'
        call findkey (narg,arg,key,nkey)
	if (nkey.gt.0) then
          L = numchar(arg(nkey+1))
	  s0 = arg(nkey+1)(1:L)
          read (s0,*) doreg
        endif

        nkey = 0
        key = '-v'
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) verbose = .true.

	nkey = 0
        key = '-u'
        call findkey (narg,arg,key,nkey)
        if (nkey.gt.0) unitest = .true.



	return
	end


	subroutine findkey (narg,arg,key,nkey)
	character*(*) arg(narg),key
	character*200 string
	integer nkey,nk,n,j,narg,numchar

	nk = numchar(key)
	nkey = 0

	do j=1,narg
	  string = arg(j)
	  n = numchar(string)

	  if (string(1:n).eq.key(1:nk)) then
		nkey = j
		return
	  endif
	enddo

	return
	end

	
