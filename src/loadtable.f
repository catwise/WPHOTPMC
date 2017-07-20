	subroutine loadtable (NSUMmax,mdex,nsr,RAlist,DEClist,IDlist,
     +                        Rsat,J,nf,nsrc,NsrcAll,unitest)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*(*) mdex
	character*500 string,header
	character*25 s0

	real*8 RAlist(NSUMmax),DEClist(NSUMmax),ra,dec
	real*4 Xpos(nsr,4),Ypos(nsr,4), Rsat(NSUMmax,4)
	integer nsrc,nsr,IDlist(NSUMmax)
	integer is
	logical unitest


c	write (6,*) 'inside loadtable'
c	write (6,*) NSUMmax,nsr

	open (unit=11,file=mdex)

	header = ''
	ihead = 0
	do 110 K=1,nsr*99
	  read (11,'(a)',end=19) string
          if (string(1:1).eq.'\') goto 110
          if (string(1:1).eq.'|') then 
            if (ihead.eq.0) then
                ihead = 1
                header = string
                goto 19
            endif
	  endif
 110	continue
 19	close (11)


c got the header; get locations of keywords
	L = numchar (header)
	do k=1,L
	  if (header(k:k).eq.'|') header(k:k)=' '
	enddo


	s0 = 'Src'
        call header_parse (header,s0,isrc)


	s0 = 'RA'
        call header_parse (header,s0,ira)
	s0 = 'Dec'
        call header_parse (header,s0,idec)

	s0 = 'SNR'
        call header_parse (header,s0,isnr)

	s0 = 'Rsat1'
        call header_parse (header,s0,irsat1)
	s0 = 'Rsat2'
        call header_parse (header,s0,irsat2)
	s0 = 'Rsat3'
        call header_parse (header,s0,irsat3)
	s0 = 'Rsat4'
        call header_parse (header,s0,irsat4)

c	write (6,*) irsat1,irsat2,irsat3,irsat4
c	call exit(0)


c	s0 = 'x_1'
c        call header_parse (header,s0,ix1)
c	s0 = 'y_1'
c        call header_parse (header,s0,iy1)

c	s0 = 'x_2'
c        call header_parse (header,s0,ix2)
c        s0 = 'y_2'
c        call header_parse (header,s0,iy2)

c	s0 = 'x_3'
c        call header_parse (header,s0,ix3)
c        s0 = 'y_3'
c        call header_parse (header,s0,iy3)

c	s0 = 'x_4'
c        call header_parse (header,s0,ix4)
c        s0 = 'y_4'
c        call header_parse (header,s0,iy4)


	open (unit=11,file=mdex)

	nsrc = 0

	do 100 K=1,nsr*99

	  read (11,'(a)',end=99) string

	  if (string(1:1).eq.'\') goto 100
	  if (string(1:1).eq.'|') goto 100



	  nsrc=nsrc+1
	  NsrcAll = NsrcAll + 1

c	  read (string,*) is,ra,dec,snr

	  call sfields (string,isrc,s0)
	  read (s0,*) is

	  call sfields (string,ira,s0)
          read (s0,*) ra
	  call sfields (string,idec,s0)
          read (s0,*) dec
c	  call sfields (string,isnr,s0)
c          read (s0,*) snr

c	  RAlist (nsrc(J), J) = ra
c	  DEClist (nsrc(J), J) = dec
c	  IDlist (nsrc(J), J) = is

	  RAlist (NsrcAll) = ra
	  DEClist (NsrcAll) = dec
 	  IDlist (NsrcAll) = is

	  do ib=1,4
		Rsat(NsrcAll,ib) = 0.
	  enddo
		
	  if (irsat1.gt.0) then
	    s0 = '0.'
	    call sfields (string,irsat1,s0)
	    read (s0,*) Rsat(NsrcAll,1)
	  endif

	  if (irsat2.gt.0) then
	    s0 = '0.'
	    call sfields (string,irsat2,s0)
            read (s0,*) Rsat(NsrcAll,2)
	  endif

	  if (irsat3.gt.0) then
	    s0 = '0.'
	    call sfields (string,irsat3,s0)
            read (s0,*) Rsat(NsrcAll,3)
	  endif

	  if (irsat4.gt.0) then
	    s0 = '0.'
	    call sfields (string,irsat4,s0)
            read (s0,*) Rsat(NsrcAll,4)
	  endif


c	  do ib = 1,4
	
c	  table (nsrc,1,ib) = is * 1.
c	  table (nsrc,2,ib) = ra
c          table (nsrc,3,ib) = dec
c	  table (nsrc,4,ib) = snr
c
c	  enddo

c	  ib = 1
c	  I = ix1
c	  call sfields (string,I,s0)
c	  if (s0(1:1).eq.'-') s0 = '0.0'
c	  read (s0,*) table (nsrc,5,ib)    !  x

c	  I = iy1
c	  call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,6,ib)    !  y

c	  Xpos(nsrc,ib) = table (nsrc,5,ib) 
c          Ypos(nsrc,ib) = table (nsrc,6,ib)

c	  ib = 2
c          I = ix2
c          call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,5,ib)    !  x

c          I = iy2
c          call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,6,ib)    !  y

c          Xpos(nsrc,ib) = table (nsrc,5,ib)
c          Ypos(nsrc,ib) = table (nsrc,6,ib)


c	  ib = 3
c          I = ix3
c          call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,5,ib)    !  x

c          I = iy3
c          call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,6,ib)    !  y

c          Xpos(nsrc,ib) = table (nsrc,5,ib)
c          Ypos(nsrc,ib) = table (nsrc,6,ib)


c	  ib = 4
c          I = ix4
c          call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,5,ib)    !  x

c          I = iy4
c          call sfields (string,I,s0)
c          if (s0(1:1).eq.'-') s0 = '0.0'
c          read (s0,*) table (nsrc,6,ib)    !  y

c          Xpos(nsrc,ib) = table (nsrc,5,ib)
c          Ypos(nsrc,ib) = table (nsrc,6,ib)


	   if (unitest) then
		if (NsrcAll .gt. 99 ) goto 99
	   endif

 100	continue
 99	close (11)


	return
	end

