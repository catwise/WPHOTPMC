c
c---- NUMCHAR ----------------------------------------------------
c
        INTEGER FUNCTION NUMCHAR(CSTRING)

C       This function determines the length of the character string
C       in cstring.

C       Author: Richard J. Stover
C       Added strong typing: D. Van Buren, T. Jarrett

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

        CHARACTER*(*) CSTRING
c       Integer*4     I
        Integer(4)    I                     !NRG B60616

        IF(CSTRING .EQ. ' ') THEN
          NUMCHAR = 0
        ELSE
        DO 8701 I=LEN(CSTRING),1,-1

        IF (CSTRING(I:I) .NE. ' ' .AND. ICHAR(CSTRING(I:I)) .NE. 0)
     &        GOTO 50
        IF (ICHAR(CSTRING(I:I)) .EQ. 0) CSTRING(I:I) = ' '
 8701      CONTINUE
 50        NUMCHAR=I
         END IF

         RETURN
         END

      integer function fitchar(cstring,chr)
c       written by T. jarrett

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      CHARACTER*(*) CSTRING,chr

      L = numchar (chr)
      ilen = numchar (cstring)

      IF(CSTRING .EQ. ' ') THEN
            fitchar = 0
      else
            do 100 i=1,ilen-(L-1)
                  if (CSTRING(I:I+(L-1)).eq.chr(1:L)) then
                        fitchar = i
                        goto 50
                  endif

 100            continue
            fitchar = 0
      endif       

 50      return
      end

      subroutine sclean (string)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)


      CHARACTER*(*) STRING
      character*500 tmp

      L = numchar (string)

      idex = 1
      do 100 i=1,L
            if (string(i:i).ne.' ') then
                  idex = i
                  goto 50
            endif
 100      continue

 50      tmp(1:L) = string(1:L)
      string = ' '
      string(1:L-(idex-1)) = tmp(idex:L)

      return
      end

c     real*8 function rstringcon (char)
      real(8) function rstringcon (char)          !NRG BG0616

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) char

      read (char,*) rstringcon

      return
      end

c     real*4 function rstrcon  (char)
      real(4) function rstrcon  (char)            !NRG B60616
      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) char
            read (char,*) rstrcon
      return
      end

      function istringcon (char)
      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) char

      read (char,*) istringcon

      return
      end

      subroutine access(char,zexist,erase)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) char
      logical erase,zexist
      character*30 command

      zexist=.false.

      inquire (file=char,exist=zexist)
      if ((zexist).and.(erase)) then
            command = '/bin/rm'
            call unix (command,char)
      endif

      return
      end

      subroutine unix(command,char)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) char,command
      character*500 result

      result = ''

      ic = numchar (command)
c      ichr = len (char)
      ic2 = numchar (char)

      do k=1,ic
      	result (k:k) = command(k:k)
      enddo

      result (ic+1:ic+1)=' '
      do k=1,ic2
      	result(k+ic+1:k+ic+1) = char(k:k)
      enddo

c      ii = len(result)
      ii = numchar (result)
c      write (6,'(a)') result(1:ii)
      i = system (result(1:ii))
      
      if (i .ne. 0) print *,'stringstuff WARNING: system code =',i,   ! JWF B60619
     +                    ' on command ',result(1:ii)                 ! JWF B60619
      
      return
      end


c read string fields
      subroutine sfields (string,nf,sout)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) string,sout
      character*1 ch
      integer fields(5000,2)

c assumptions: 1st field = first non-blank string

      sout = ' '

      L = NUMCHAR (string)
      n = 0

      imax = 0
      do 100 j=1,L

            if (j.le.imax) goto 100
            ch = string(j:j)
            if (ch.ne.' ') then
                  imin = j
                  do i=j+1,L
                        ch = string(i:i)
                        if (ch.eq.' ') then
                              imax = i-1
                              n = n + 1
                              fields (n,1) = imin
                              fields (n,2) = imax
                              if (n.eq.nf) goto 99
                              GOTO 100
                        endif
                  enddo
                  n = n + 1
                  fields (n,1) = imin
                  fields (n,2) = L
                  imax = L
                  if (n.eq.nf) goto 99
            endif
 100      continue

      
 99      if (n.lt.nf) then
            sout = ' '
      else
            i1 = fields (nf,1)
            i2 = fields (nf,2)
            sout = string(i1:i2)
      endif


      return
      end



       subroutine header_parse (head,keyword,idex)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

        character*(*) head,keyword
c       character*25 s0,result,upcase                   ! JWF B60703
        character*50 s0,result

c      result = upcase (keyword)                        ! JWF B60703
c      keyword = result                                 ! JWF B60703

c        call upcase(keyword)                           ! JWF B60703
         call upcase(keyword)                           ! JWF B60703

        if (keyword(1:4).eq.'NULL') then
                idex=0
                return
        endif

c       L = numchar (head)                                              ! JWF B60519
        LK = numchar(keyword)


        imax = 150000
        idex = 0

        do 50 I=1,imax
                call sfields (head,I,s0)
c               call upcase(s0)                            ! JWF B60703
                call upcase(s0)                            ! JWF B60703
c               result = upcase(s0)                        ! JWF B60703
c               s0 = result                                ! JWF B60703
                LL = numchar(s0)
                if (LL.eq.0) goto 47
                if (LL.ne.LK) goto 50

                if (s0(1:LL).eq.keyword(1:LK)) then
                        idex=I
                        goto 47
                endif
 50     continue

 47     return
        end

c
c Stolen from JWF, inserted by TPC
c
c-----------------------------------------------------------------------------
c---------------------------------- JWF B60703 -------------------------------
c----------------un-commented out original subroutine upcase------------------
c-------------------commented out functione upcase----------------------------
c-----------------------------------------------------------------------------
       subroutine upcase(field)
       character*25 field
       character*1  tmpchar
       integer*4    i,k,lnblnk
       byte         tmpbyte
       equivalence (tmpbyte,tmpchar)
       k = lnblnk(field)
       do 10 i = 1, k
         tmpchar = field(i:i)
         if ((tmpbyte .gt. 96) .and. (tmpbyte .lt. 123)) then
           tmpbyte = tmpbyte - 32
           field(i:i) = tmpchar
         end if
10      continue
        return
        end


c     function upcase(string) result(upper)
c     character(len=*), intent(in) :: string
c     character(len=len(string)) :: upper
c      integer :: j, L                                                     ! JWF B60619
c     integer :: j                                                        ! JWF B60619
c
c     do j = 1,len(string)
c       if(string(j:j) >= "a" .and. string(j:j) <= "z") then
c            upper(j:j) = achar(iachar(string(j:j)) - 32)
c       else
c            upper(j:j) = string(j:j)
c       end if
c     end do
c     end function upcase

      subroutine stringswap (string,s1,s2)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) string,s1,s2
      character*20 s0

      L = numchar(string)

      L1 = numchar(s1)
      L2 = numchar(s2)

c SNR1 >> w1snr  (4 goes to 5 chars)

      index = 0
      do 500 J = 0,L-L1
         jl = J+1
         jh = J+L1 
         s0 = string(jl:jh)
         
         if (s0(1:L1).eq.s1(1:l1)) then
            ! match

            ! check neighbors
            if ((string(jl-1:jl-1).ne.' ').and.(string(jl-1:jl-1).ne.'|')) goto 500
            if ((string(jh+1:jh+1).ne.' ').and.(string(jh+1:jh+1).ne.'|')) goto 500

            il = J+1
            ih = J+L1

            if (L2.eq.L1) then

                string(il:ih) = s2(1:l2)
               write (6,*) 'ok'
               return

            else if (L1.gt.L2) then
                  write (6,'(a,a)') 'beware, string to small: ',s2(1:l2)
                  return
            else if (L2.gt.L1) then
               idiff = L2-L1
               il2 = il - idiff
            
               s0 = string(il2:ih)
               call checkstring (s0,L2,igo)

               if (igo.eq.0) then
                  write (6,'(a)')  string(il2:ih)
                  write (6,'(a)') s2(1:l2)

                  string(il2:ih) = s2(1:l2)
                  write (6,*) 'ok'
                  return

               else
                  write (6,'(a,a)')  'trouble: ',s0(1:l2)

            ! pad
                   string = string(1:il-1) // s2(1:l2) // string(ih+1:L)

                  write (6,*) 'Add padding to format statement:  ',idiff
                  write (6,*) 'ok'
                   return
            
               endif

               return
            endif

         endif
 500      continue

      return
      end

      subroutine checkstring (s0,L2,igo)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

        character*(*) s0

      igo = 0
      do K=1,L2
         if (s0(K:K).eq.'|') igo = 1
      enddo

      return
      end



c read string fields
      subroutine sfields2 (string,nf,sout, i1,i2)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) string,sout
      character*1 ch
      integer fields(500,2)

c assumptions: 1st field = first non-blank string

      sout = ' '
      i1 = 0
      i2 = 0

      L = NUMCHAR (string)
      n = 0

      imax = 0
      do 100 j=1,L

            if (j.le.imax) goto 100
            ch = string(j:j)
            if (ch.ne.' ') then
                  imin = j
                  do i=j+1,L
                        ch = string(i:i)
                        if (ch.eq.' ') then
                              imax = i-1
                              n = n + 1
                              fields (n,1) = imin
                              fields (n,2) = imax
                              if (n.eq.nf) goto 99
                              GOTO 100
                        endif
                  enddo
                  n = n + 1
                  fields (n,1) = imin
                  fields (n,2) = L
                  imax = L
                  if (n.eq.nf) goto 99
            endif
 100      continue

      
 99      if (n.lt.nf) then
            sout = ' '
      else
            i1 = fields (nf,1)
            i2 = fields (nf,2)
            sout = string(i1:i2)
      endif


      return
      end


      subroutine nullswap (string,newstring)

      implicit integer (i-n)
        implicit real(4) (a-h)
        implicit real(4) (o-z)

      character*(*) string,newstring
      character*4 null

      null = 'null'

      L = numchar (string)
      newstring = string(1:L)

      if (L.lt.4) return   ! string too small to stuff NULL into

      Ndiff = L - 4

      do K=1,Ndiff
            newstring(K:K) = ' '
      enddo

      newstring(L-3:L) = null

      return
      end

      subroutine real_stringform (value, form, string, L)

      real(4) value
      character*(*) form, string
      character*25 format
      integer L, numchar
      logical IzBad

      if (IzBad(value)) value = 0.

      L = numchar(form)
      format = '(' // form(1:L) // ')'
      write (string,format) value
      L = numchar (string)

      return
      end


      subroutine double_stringform (dvalue, form, string, L)

        real(8) dvalue  
      real(4) value
        character*(*) form, string
        character*25 format
        integer L, numchar
        logical IzBad

      value = dvalue * 1.
        if (IzBad(value)) dvalue = 0.d0

        L = numchar(form)
        format = '(' // form(1:L) // ')'
        write (string,format) dvalue
        L = numchar (string)

        return
        end


      subroutine int_stringform (value, form, string, L)
      real(4) rvalue
      integer value
      character*(*) form, string
        character*25 format
        integer L, numchar
      logical IzBad

      rvalue = value * 1.
        if (IzBad(rvalue)) value = 0

        L = numchar(form)
        format = '(' // form(1:L) // ')'
        write (string,format) value
        L = numchar (string)

        return
        end
