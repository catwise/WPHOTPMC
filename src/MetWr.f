	subroutine MetWr (imeta, ib, What, TYPE, metval, comment)
	character*(*) TYPE, What,comment
	character*200 string
	character*6 root
	character*19 s0
	real*4 metval
	integer numchar, LW, imetval, LC, ib, imeta, npad
	integer L, ns, j, Li

	logical doint

	doint = .false.
	if (TYPE(1:1).eq.'i') then
		imetval = int (metval)
		doint = .true.
	endif

	LW = numchar(What)
	LC = numchar(comment)

	root = 'WPHot:'

	string = ' ' // root // What(1:LW)
	L = numchar (string)

	ns = 26
	do j=L+1,ns
	  string(j:j) = ' '
	enddo

	write (s0,'(i4)') ib

	ns = ns + 1
	string (ns:ns+3) = s0(1:4)
	ns = ns + 3

	do j=ns+1,ns+6
          string(j:j) = ' '
        enddo
	ns = ns + 6

	string(ns:ns) = TYPE(1:1)

	do j=ns+1,ns+2
          string(j:j) = ' '
        enddo
        ns = ns + 2


	if (doint) then
		do j=ns+1,ns+10
		 string(j:j) = ' '
		enddo
		ns = ns + 10
		write (s0,'(i7)')  imetval
		string (ns+1:ns+7) = s0(1:7)
		ns = ns + 7

	else
		do j=ns+1,ns+4
                 string(j:j) = ' '
                enddo
                ns = ns + 4
		write (s0,'(e13.5)')  metval
		string (ns+1:ns+13) = s0(1:13)
                ns = ns + 13
	endif

	
	do j=ns+1,ns+3
	 string(j:j) = ' '
	enddo
	ns = ns + 3

	string (ns+1:ns+LC) = comment (1:LC)
	ns = ns + LC

	write (imeta,'(a)') string(1:ns)

	return
	end
