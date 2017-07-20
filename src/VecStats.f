!   ave == trimmed average using low-side sigma
!   sdev = low-side sigma
!   hsig = histogram sigma (average of low and high-side sigma)
!   zmed == histogram median

	subroutine VecStats (n,vec,m,vout,niter,vlowmin,vhighmax,
     1     ave,sdev,zmed, SKYrms, hsig, doCOADD)   ! where SKYrms is the sky sigma based on the UNC

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 vec(n),vout(n)
	logical doCOADD

c	write (88,*) 'inside vecstats ',SKYrms, vlow,vhigh

	nsigma = 5  ! changed from a value of "5" to avoid more outliers

c	doCOADD = .false.
c	if (vhigh.eq.4999.) doCOADD = .true.

	if (doCOADD) then
		nsigma = 5
		niter = 4
	endif

	hsig = 0.
	sig_low = 0.
	ave = -9999.
        sdev = SKYrms
        zmed = -9999.
	vlow = vlowmin
	vhigh = vhighmax

	do 1000 Jit = 1,niter

	  if (Jit.eq.niter) nsigma = 2

	  m = 0
	  do 50  J = 1, n
	   if (vec(J).lt.vlow) goto 50
	   if (vec(J).ge.vhigh) goto 50
	   m=m+1
	   vout(m)=vec(J)
c	if (jit.eq.1) write (121,*) vec(j), vlow,vigh
 50	  continue

	  if ((Jit.eq.1).and.(m.gt.8)) then
		 call MOMENT(vout,m,AVE,SDEV)  ! mean sky, outliers and all
		 vlow = ave - (sdev * 10.)
		 vhigh = ave + (sdev * 10.)

		 vlow = max (vlow,vlowmin)
		 vhigh = min (vhigh,vhighmax)

	  endif

c	write (88,*) 'kere ',vlow,vhigh,m
c	write (6,*) 'kere ',vlow,vhigh,m,  SKYrms   ! TEMP

	if (m.gt.8) then

c	   if ((doCOADD).and.(Jit.gt.1)) then
c		 call MOMENT(vout,m,AVE,SDEV)
c	   endif

c	   call MOMENT(vout,m,AVE,SDEV)
c	   call MDIAN1(vout,m,zMED)


	   vwidth = sdev / 10.
	   if (doCOADD) vwidth = vwidth / 2.
	   if (Jit.eq.1) vwidth = SKYrms / 3.   ! use the original sigma_pop
	   
           vwidth = max (vwidth,0.0001)
	   vwidth = min (vwidth,SKYrms*5.)


c	 write (88,*) '#  SKYrms = ',SKYrms,'  vwidth = ',vwidth

c	write (6,*) vlow,vhigh,vwidth

           call bin_stats (vlow,vhigh,vwidth,vout,m,zMED, hsig, z16, z84)
	   sdev = zMED - z16   ! this is the low-side sigma
	!  now using the mode and it's sigma ; Jan 19, 2011
	   if (hsig.gt.0.) then
		sdev = hsig
	   endif

c	write (6,*) 'hsig ',hsig,sdev,z84-zmed ! TEMP
c	if (Jit.eq.niter) then
c		write (6,*) hsig,sdev,nsigma
c		call exit(0) ! temporary
c	endif

	else
	   goto 1001
	endif

	vlow = zmed - (sdev * nsigma * 1.)
	vhigh = zmed + (sdev * nsigma * 1.)
	vlow = max (vlow,vlowmin)
        vhigh = min (vhigh,vhighmax)

c	write (6,*) 'deer ',jit, m, ave, sdev, zmed, vlow, vhigh    ! TEMP

 1000	continue


cccc we should now have the mode and the mode's sigma
	if (hsig.gt.0.) then
	   ave = zmed	
	   return
	endif	   
!! otherwise go with the trimmed histogram average

	m = 0
        do 60  J = 1, n
           if (vec(J).lt.vlow) goto 60
           if (vec(J).ge.vhigh) goto 60
           m=m+1
           vout(m)=vec(J)
 60     continue


	! trimmed average
          call MOMENT(vout,m,AVE,stdev)

c	write (6,*) 'trimmed average: ',m,AVE,stdev

c	write (121,'(a,6f10.4)') '# ',ave,stdev,vlow,vhigh, zMED - z16, z84-zmed
c	close (121)

	! TEMPORARY
c	ave = ave * 0.94

c	write (6,*) nsigma
c	write (6,*) zmed,hsig
c	write (6,*) ave, sdev
c	call exit(0)

 1001	return

	end
