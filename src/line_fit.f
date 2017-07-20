	subroutine line_fit (nf,N,zmed,JD8,JD,flux,eflux,slope, eslope, chi2, boffset)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (nterms = 2)

	real*4 flux(nf), eflux(nf), zmed, chi2, chisq, alamda, JD(999)
	real*4 A(nterms), covar (nterms,nterms), alpha (nterms,nterms)
	real*8 JD8(nf), x8, t8
        integer lista(nterms)

c	write (6,*) N, zmed, JD(1)

	x8 = 1.d32
	do j=1,n
		if (JD8(j).lt.x8) then
			x8 = JD8(J)
		endif
	enddo

	
	do j=1,n
		t8 = JD8(j) - x8
		JD(j) = t8 * 1.0
c	write (6,*) JD8(j),x8,JD(j)
	enddo


	call SORT3(N,JD, flux, eflux)

c	x = JD(1)
c	do j=1,n
c	   JD(j) = JD(j) - x
c	enddo


c	write (6,*) N, zmed, JD(1), JD(2)

c	write (6,*) (JD(kk),kk=1,N)
c	write (6,*) (flux(kk),kk=1,N)
c	write (6,*) (eflux(kk),kk=1,N)


	A(1) = zmed
	A(2) = 0.00

	lista(1) = 1
	lista(2) = 2
	alamda = -1.

	CHISQ1=1.E9
	termold = 1.e19
	chisq = 0.


	        alamda=-1  
	        call mrqmin(JD, flux, eflux, N,A,lista,Nterms,covar,alpha,  
     *       Nterms,chisq,alamda)  
	        k=1  
	        itst=0  
c 1           write(*,'(/1x,a,i2,t18,a,f10.4,t43,a,e9.2)') 'Iteration #',k,  
c     *       'Chi-squared:',chisq,'ALAMDA:',alamda  
c	        write(*,'(1x,t5,a,t13,a,t21,a,t29,a,t37,a,t45,a)') 'A(1)','A(2)'
c	        write(*,'(1x,6f10.4)') (a(i),i=1,nterms)  
 1	        k=k+1  
	        ochisq=chisq  
	        call mrqmin(JD, flux, eflux,N,A,lista,Nterms,covar,alpha,
     *       Nterms,chisq,alamda)
	        if (chisq.gt.ochisq) then  
	          itst=0  
	        else if (abs(ochisq-chisq).lt.0.1) then  
	          itst=itst+1  
	        endif  
	        if (itst.lt.5) then  
	          goto 1  
	        endif  
	        alamda=0.0  
	        call mrqmin(JD, flux, eflux,N,A,lista,Nterms,covar,alpha,
     *       Nterms,chisq,alamda)
c	        write(*,*) 'Uncertainties:'  
c	        write(*,'(1x,6f8.4/)') (sqrt(covar(i,i)),i=1,nterms)  

	slope = A(2)
	eslope = sqrt(covar(2,2))
	chi2 = chisq
	boffset = A(1)

	return
	END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC     FUNCTION galaxy profile    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        FUNCTION fun_poly (X,A,N)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


        DIMENSION A(N)

	sum = 0.
	do K = 1, N
		sum = sum + (A(K) * (x**(K-1)) )
	enddo
	fun_poly = sum

c  fun_poly = A(1) + (A(2) * x) + (A(3) * x * x) + (A(4) * x * x * x)

        return
        end


        subroutine funcpoly (x,A,y,dyda,N)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        dimension A(N),dyda(N)

c	write (6,*) n,x, A(1), A(2)

	Y = 0.
        do K = 1, N
                Y = Y + (A(K) * (x**(K-1)) )
		if (K.eq.1) dyda(1) = 1.
		if (K.eq.2) dyda(2) = x
		if (K.eq.3) dyda(3) = x ** 2
		if (K.eq.4) dyda(4) = x ** 3
        enddo

c	write (6,*) 'x,y = ',x, Y,  dyda(1), dyda(2)
c	call exit(0)


c        Y = A(1) + (A(2) * x) + (A(3) * x * x) + (A(4) * x * x * x)

        return
        end




      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *alamda)
	
	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL*4 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,m,mfit
      REAL*4 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      j=0
      do 14 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END


      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL*4 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      PARAMETER (MMAX=20)
      INTEGER mfit,i,j,k,l,m
      REAL*4 dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcpoly(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END


      SUBROUTINE gaussj(a,n,np,b,m,mp)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


      INTEGER m,mp,n,np,NMAX
      REAL*4 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*4 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
c                 write (6,*) 'line_fit error: singular matrix in gaussj'
		 return
c		call exit(9)
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) then
c		write (6,*) 'line_fit error: singular matrix in gaussj'
		return
c                call exit(9)
	endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END


      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


      INTEGER ma,mfit,npc,ia(ma)
      REAL*4 covar(npc,npc)
      INTEGER i,j,k
      REAL*4 swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END


