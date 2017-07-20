      Subroutine MatInv(A,Ainv,N,OK)
c-----------------------------------------------------------------------
c     Invert a matrix A(N,N) via LU decomposition;
c     uses work spaces T(N,N) and Indx(N);
c     returns inverse in Ainv(N,N)
c-----------------------------------------------------------------------
c
      Integer*4 N, I, J, Indx(N)
      Real*4    A(N,N), Ainv(N,N), T(N,N), D
      Logical*4 OK
c
c-----------------------------------------------------------------------
c
      Do 20 J = 1, N
        Indx(J) = J
        Do 10 I = 1, N
          T(I,J) = A(I,J)
10      Continue
20    Continue
c
      Call ludcmp(T,N,N,Indx,D,OK)
      If (.not.OK) Return
c
      Do 40 J = 1, N
        Do 30 I = 1, N
          If (I .eq. J) then
            Ainv(I,J) = 1.0
          Else
            Ainv(I,J) = 0.0
          End If
30      Continue
40    Continue
c
      Do 50 J = 1, N
        call lubksb(T,N,N,Indx,Ainv(1,J))
50    Continue
c
      Return
      End
c
c=======================================================================
c
      SUBROUTINE ludcmp(a,n,np,indx,d,OK)
      INTEGER*4 n,np,indx(n)
      REAL*4 d,a(np,np),TINY
      PARAMETER (TINY=1.0e-20)
      INTEGER*4 i,imax,j,k
      REAL*4 aamax,dum,vv(np)
      REAL*8 sum
      logical*4 OK
c
      OK = .True.
      d=1.
      imax = 1                 ! make compiler happy
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
          OK = .False.
          Return
        end if
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*4 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
