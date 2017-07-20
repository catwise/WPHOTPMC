	subroutine keyhead (Hfits,key,s0)
	
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        character*(*) Hfits,key,s0
        character*80 line

        s0 = ''
        L = numchar (Hfits)

        nlines = nint(l / 80.)
c        write (6,*) 'nlines = ',nlines

        L = numchar (key)

c       write (6,'(a)') key(1:L)

        iend = 0
        do J=1,nlines
          idex = iend + 1
          iend = idex + 79

	  if (iend.gt.150000) return

c	write (6,*) J,idex,iend

          line = Hfits (idex:iend)

c	write (6,'(a)') line(1:80)

          I = 1
          call sfields (line,I,s0)
          M = numchar (s0)

c	write (6,'(a)') s0(1:M)

          if (s0(1:M).eq.key(1:L)) then
                I = 3
                call sfields (line,I,s0)
                return
          endif

	  s0 = ''

        enddo


        return
        end

