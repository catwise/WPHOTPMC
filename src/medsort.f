      subroutine medsort(a,n,amed,q1,q2)

! Calculate median after sorting with HEAPSORT algorithm.
! Also calculate the 15.87% and 84.13% quantiles (q1 and q2).

      implicit real(4) (a-h,o-z)
      implicit integer(i-n)
      real(4), allocatable :: as(:)
      real(4) a(*)
        logical sortit

      sortit = .true.
      if(n .lt. 0) then
          sortit = .false.
        n = -n
        endif

      if (n==1) amed = a(1)
      if (n <= 1) return
      allocate (as(n))
      as = a(1:n)
      if(sortit) call sort(n,as)

      if (mod(n,2) == 1) then
          amed = as(n/2+1)
      else
          amed = (as(n/2) + as(n/2+1))/2.
      endif

      q1 = as(max(nint(0.1587*n),1))
      q2 = as(nint(0.8413*n))
      deallocate(as)
      return

      end

      subroutine trimsort(a,n,amed,q1,q2,pct,trimmean,trimsig)

! Calculate median and trimmed average after sorting with HEAPSORT algorithm.
! Also calculate the 15.87% and 84.13% quantiles (q1 and q2).

      implicit real(4) (a-h,o-z)
      implicit integer(i-n)
      real(4), allocatable :: as(:)
      real(4) a(*),amed,q1,q2,pct,trimmean,trimsig
        logical sortit
	  integer(4) ibot, itop                            ! JWF B60713
      

      trimmean = -999
      trimsig = -999

      sortit = .true.
      if(n .lt. 0) then
          sortit = .false.
        n = -n
        endif

      if (n==1) amed = a(1)
      if (n <= 1) return

      allocate (as(n))
      as = a(1:n)

      if(sortit) call sort(n,as)

      if (mod(n,2) == 1) then
          amed = as(n/2+1)
      else
          amed = (as(n/2) + as(n/2+1))/2.
      endif

      q1 = as(max(nint(0.1587*n),1))
      q2 = as(nint(0.8413*n))

      if(pct<=0) then
          !!!print *, '!!!! trimsort: Pct=',pct
          deallocate(as)
          return
        endif

        ! Trimmed average

        bot = int(pct*n)
        top = int((1-pct)*n) + 1
		ibot = bot                                      ! JWF B60713
		itop = top                                      ! JWF B60713

        !!!print *,'--- trimdort: pct,n,bot,top = ',pct,n,bot,top

      if(bot .lt. 1 .or. top .gt. n .or. top-bot .lt. 4) then      
          deallocate(as)
          return ! Not enough points for trim
        endif

        sum = 0
        sumsq = 0
      m = 0
c       do i=bot+1,top-1                             ! JWF B60713
        do i=ibot+1,itop-1                           ! JWF B60713
        m = m + 1
          sum = sum + as(i)
          sumsq = sumsq + as(i)**2
        end do

        trimmean = sum/m
        trimsig = sqrt(max(sumsq/m - trimmean**2, 0.0)*m/(m-1.))/sqrt(1.*m)

      deallocate(as)
      return

      end
