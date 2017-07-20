	subroutine SigClipWmean (n,f,w, f2, w2, AVE,SDEV,sdev2, N_M_2)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (Nmin = 5)

	real*4 f(n), w(n), f2(n), w2(n)
	logical debug

	debug = .false.

   
	N_M_2 = n

c step 0; demand list to be >=Nmin sources

	if (n.lt.Nmin) then

		! do not sig-clip;  instead, compute straight Weighted mean
		
		call WMOMENT(f,w,n,AVE,SDEV,sdev2)
		return

	endif
	
c step 1, sort list

	call SORT2(N,f,w)


c now f should be sorted

	if (debug) then

	do j=1,n
	  write (6,*) j,f(j),w(j)
	enddo

	endif

c  step 2, compute plain median

	N2=N/2

      IF(2*N2.EQ.N)THEN
        XMED=0.5*(f(N2)+f(N2+1))
      ELSE
        XMED=f(N2+1)
      ENDIF

	if (debug) write (6,*) 'median = ',xmed

c  step 3, drop max value;  if large sample, drop the low also
c	if (n.gt.8) then
c	  jl = 2
c	  jh = n-1
c	else
c	  jl = 1
c	  jh = n-1
c	endif

c step 3, drop min/max for large sample
	jl = 1
	jh = n
	if (n.gt.8) then   ! demand 9 samples to min/max reject
         jl = 2
         jh = n-1
	endif

	m=0
	do j=jl,jh
	  m=m+1
	  f2(m)=f(j)
	  w2(m)=w(j)
	enddo

c  step 4, compute weighted mean & sig

	call WMOMENT(f2,w2,m,AVE,SDEV,sdev2)

	if (debug) write (6,*) 'weighted mean = ',m,ave,sdev2

c using sigma, compute low/high thresholds

	do loop = 1, 3

	zhigh = xmed + (3. * sdev2)
	zlow = xmed - (3. * sdev2)

	
	if (debug) write (6,*) 'thresholds ',zlow,zhigh

c eliminate all sources  beyond limits

	m=0
	do 50 j=1,n
	  val = f(j)

	  if (val.lt.zlow) goto 50
	  if (val.gt.zhigh) goto 50
 
	  m=m+1
	  f2(m) = val
	  w2(m) = w(j)

	if (debug) write (6,*) m,val,w(j)

 50	continue


        if (m.lt.3) then
                ! do not sig-clip;  instead, compute Weighted mean

                call WMOMENT(f,w,n,AVE,SDEV,sdev2)
                return
        endif


 	N2=m/2

      IF(2*N2.EQ.m)THEN
        XMED=0.5*(f2(N2)+f2(N2+1))
      ELSE
        XMED=f2(N2+1)
      ENDIF

c now compute weighted moment

 47	call WMOMENT(f2,w2,m,AVE,SDEV,sdev2)

	if (debug) then
	write (6,*) 'median ',xmed
	write (6,*) 'weighted mean = ',m,ave,sdev2
	endif

	N_M_2 = m

	enddo  ! loop

	return
	end
