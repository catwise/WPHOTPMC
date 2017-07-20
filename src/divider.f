	subroutine divider (kreg,nx,ny,inum,jnum, i0, j0, isize, jsize)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	isize = nint(nx / inum * 1.) 
	jsize = nint(ny / jnum * 1.) 

	ihalf = nint(isize / 2.)
	jhalf = nint(jsize / 2.)

	n = 0
	j0 = -jhalf

	icount = 0
	do L=1,jnum
		j0 = j0 + (jhalf*2)
		jl = j0 - jhalf + 1
		jh = j0 + jhalf + 1
		jl = max (1,jl)
		jh = min (jh,ny)
		i0 = -ihalf

		ibeg = 1
                iend = inum
                iord = 1

		if (icount.eq.1) then
		 ibeg = inum
		 iend = 1
		 iord = -1
		 icount = 0

		else
			icount = 1
		endif
	

c	write (6,*) 'cheer ',kreg,L,ibeg,iend,iord

	do k=ibeg,iend,iord

		i0 =  ihalf * ( (2*k) - 1)
c		i0 = i0 + (ihalf*2)
		il = i0 - ihalf + 1
		ih = i0 + ihalf + 1
                il = max (1,il)
                ih = min (ih,nx)
		
c		write (6,*) k,i0,j0

		n=n+1
		if (n.eq.kreg) then

c	write (6,*) 'cheer got ',k,i0,j0
			return
		endif
	enddo
	enddo

	return
	end	
