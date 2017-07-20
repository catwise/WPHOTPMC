	function estmode(a,n,amed,sigma)

! Estimate mode of the distribution of N values of an array A.

	implicit real(4) (a-h,o-z)
	implicit integer(i-n)
	integer, allocatable :: h(:)
	real(4) a(*)

	frsig = 0.04		! precision, expressed as a fraction of sigma

	maxsigs = 100		! maximum number of sigmas to search on
				! either side of median

	if (n <= 0) return
	binsize = sigma*frsig
	maxdevs = nint(maxsigs/frsig)
	nbins = 2*maxdevs + 1
	allocate (h(nbins))
	h = 0

	do i = 1,n
	    ibin = maxdevs + nint((a(i) - amed)/binsize)
	    if (ibin > 0 .and. ibin <= nbins) h(ibin) = h(ibin) + 1
	enddo

	npeak = -2**30
	do ibin = 1,nbins
	    if (h(ibin) > npeak) then
		npeak = h(ibin)
		ibinmode = ibin
	    endif
	enddo

	estmode = amed + (ibinmode - maxdevs)*binsize
	if (npeak < 6) estmode = amed
	deallocate(h)
	return

	end
c
c=======================================================================
c
	subroutine rmdr(a,kb,BGtrimFrac,SkThresh, back,nUsed,Skew,Sigma)

! Estimate mode of the distribution of kb values of an array backbuf
! using rmdr trimmed average.
! a() is assumed to be sorted.

      real*4    a(kb), back, Avg, Skew, VarSkew, Sigma,
     +          Median, BGtrimFrac, SkThresh, m2, m3, mu1, mu2, mu3
      real*8    Sum, SumSq, SumCb
      integer*4 kb, NRMDR, K1, K2, nUsed, Nmed, K
c
c-----------------------------------------------------------------------
c
      Skew  = 0.0
      nUsed = kb
      if (nUsed .eq. 1) then
        back  = a(1)
        Sigma = back
        return
      end if
c
      K1 = 1
      K2 = kb
c                                        ! This probably can't happen
10    if (nUsed .eq. 1) go to 30         ! because we trap it above & below
      NRMDR = BGtrimFrac*nUsed
      if (NRMDR .gt. nUsed-1) NRMDR = nUsed-1
      if (NRMDR .lt. 1)       NRMDR = 1  ! Must be safe, since nUsed > 1, and
c                                        ! we don't want an infinite loop here
      Do 20 K = 1, NRMDR
        NMed = (K1 + K2 - 1)/2
        If (Mod(nUsed,2) .eq. 0) Then
          Median = 0.5*(a(NMed) + a(NMed+1))
        Else
          Median = a(NMed)
        End If
        If ((a(K2) - Median) .gt. (Median - a(K1))) Then
          K2 = K2 - 1
        Else
          K1 = K1 + 1
        End If
        nUsed = nUsed - 1
20    Continue
c
      if (nUsed .eq. 1) then             ! if nUsed = 1, then
        back  = a(K1)                    ! K1 = K2 at this point
        Skew  = 0.0
        Sigma = back
        return
      end if
c
30    Sum   = 0.0d0
      SumSq = 0.0d0
      SumCb = 0.0d0
c
      Do 40 K = K1, K2
        Sum   = Sum   + a(K)
        SumSq = SumSq + a(K)**2
        SumCb = SumCb + a(K)**3
40    Continue
c
      mu1 = Sum/nUsed
      m2  = SumSq/nUsed
      m3  = SumCb/nUsed
c
      mu2 = abs(m2 - mu1**2)
      mu3 = m3 - 3.0*m2*mu1 + 2.0*mu1**3
c
      Sigma = sqrt(abs(mu2))
      if (mu2 .gt. 0) then
        Skew = mu3/sqrt(mu2)**3
      else
        Skew = 0.0
      end if
c
      if (nUsed .eq. 1) go to 50
      if (nUsed .lt. 18) then
        VarSkew = 0.25           ! clip prior sigma(skewness) at 0.5
      else
        VarSkew = 6.0*float(nUsed-2)/(float(nUsed+1)*float(nUsed+3)) ! Gaussian pop.
      end if
c
      if (Skew**2/VarSkew .gt. SkThresh) go to 10
c
50    back = mu1
c
      return
c
      end



c=======================================================================
c
c               Official DAO version:  1988 July 1
c
c This version of MMM (modified by PBS 1984.IV.10ff) assumes that
c the sky brightnesses in the one-dimensional array SKY are already
c sorted on entering this routine, and that pixels outside the "bad"
c limits have already been eliminated.
c
c This particular version of MMM also takes cognizance of the fact that,
c pixels falling below the LOBAD threshold already having been 
c eliminated, the contaminated sky pixels values overwhelmingly display
c POSITIVE departures from the true value.
c
c If for some reason it is impossible to obtain the mode of the sky
c distribution, this will be flagged by setting SIGMA = -1.0.
c
c Arguments
c
c     SKY (INPUT) is a real vector containing actual sorted sky values.
c    NSKY (INPUT) is the number of defined elements in SKY.
c  SKYMOD (OUTPUT) is the estimated mode of the sky values.
c   SIGMA (OUTPUT) is the computed standard deviation of the peak in
c         the sky histogram.
c    SKEW (OUTPUT) is the computed skewness of the peak in the sky
c         histogram.
c
c=======================================================================
c
      subroutine mmm(sky, nsky, hibad, skymn, skymed, skymod, sigma, skew,
     +               nused, niter)
      implicit none
      integer*4 nsky
c
      real*4 sky(nsky)
      double precision dsqrt, dble
c
      real*4 alog10, amin1, amax1
      double precision sum, sumsq
      real*4 cut, cut1, cut2, delta, skymid, skymed, skymn, skymod
      real*4 sigma, skew, r, sign, x, hibad, center, side, f1p0001
      integer*4 i, minimm, maximm, niter, istep, maxit, minsky, jstep, nused
      logical redo
      data f1p0001/1.0001/
      integer nwarn,maxwarn
      data nwarn/0/,maxwarn/50/
c
c-----------------------------------------------------------------------
c
c SECTION 1
c
      data maxit / 30 /
      data minsky / 20 /

      if (nsky .le. 0) goto 9900
c
c SKYMID is the median value for the whole ensemble of sky pixels.
c Notice that if NSKY is odd, then (NSKY+1)/2 and (NSKY/2)+1 will be the
c same number, whereas if NSKY is even, they will be different numbers.
c This same trick will be used again later.
c
c Initialize the variables for accumulating the mean and standard
c deviation, and initialize the rejection limits.
c

      skymid = 0.5 * (sky((nsky + 1) / 2) + sky((nsky / 2) + 1))

      sum = 0.d0
      sumsq = 0.d0
c
c For the first pass we will consider only pixels in a symmetric 
c interval of brightness values about the median value.  This exploits
c the assumption that all the bad pixels are already rejected from the
c lower end of the brightness range.
c

      cut1 = min(skymid - sky(1),sky(nsky) - skymid,hibad - skymid)

      !!!print *, '--',sky(1),skymid,sky(nsky),cut1

      cut2 = skymid + cut1
c

      cut1 = skymid - cut1

      minimm = 0
      do 1010 i = 1, nsky
      if (sky(i) .lt. cut1) then
      minimm = i
      goto 1010
      end if
      if (sky(i) .gt. cut2) goto 1020
      delta = sky(i) - skymid
      sum = sum + delta
      sumsq = sumsq + (delta ** 2)
      maximm = i
c
c Henceforth in this subroutine, MINIMM will point to the highest value
c rejected at the lower end of the vector, and MAXIMM will point to the
c highest value accepted at the upper end of the vector.
c MAXIMM-MINIMM is the number of pixels within the acceptance range.
c
c Compute mean and sigma (from the first pass).
c

 1010 continue

 1020 continue
      !!!print *, '--- ',cut1,cut2,minimm,maximm,nsky,sky(1),sky(nsky),skymid
      skymed = 0.5 * (sky(((minimm + maximm) + 1) / 2) + sky(((minimm + 
     &maximm) / 2) + 1))

      skymn = sum / dble(maximm - minimm)
      sigma = sqrt((sumsq / dble(maximm - minimm)) - (skymn ** 2))
c
c The middle sky value, SKYMID, was subtracted off up above and added 
c back in down here to reduce the truncation error in the computation 
c of SIGMA.
c Note that this definition of SIGMA is incorrect by a factor of
c SQRT [NSKY/(NSKY-1.)], but for all but pathological cases (where none
c of this can be trusted anyway), it's close enough.
c

      skymn = skymn + skymid

      skymod = skymn
c
c If the mean is less than the mode, that means the contamination is
c slight, and the mean value is what we really want.  Note that this
c introduces a slight bias toward underestimating the sky when
c the scatter in the sky is caused by random fluctuations rather than
c by contamination, but I think this bias is negligible compared to the
c problem of contamination.
c
c-----------------------------------------------------------------------
c
c SECTION 2
c
c Rejection and recomputation loop:
c

      if (skymed .lt. skymn) skymod = (3. * skymed) - (2. * skymn)

      niter = 0
 2000 niter = niter + 1
c
c Compute Chauvenet rejection criterion.
c         

      if ((niter .gt. maxit) .or. ((maximm - minimm) .lt. minsky)) goto
     &9900

      r = log10(float(maximm - minimm))
c
c Compute rejection limits (symmetric about the current mode).
c

      r = max(2.,(((- (.1042 * r)) + 1.1695) * r) + .8895)

      cut = (r * sigma) + (0.5 * abs(skymn - skymod))
      cut = max(1.5,cut)
      cut1 = skymod - cut
c
c Recompute mean and sigma by adding and/or subtracting sky values
c at both ends of the interval of acceptable values.
c
c At each end of the interval, ISTEP will show the direction we have to 
c step through the vector to go from the old partition to the new one.
c Pixels are added or subtracted depending upon whether the limit is 
c moving toward or away from the mode.
c

      cut2 = skymod + cut
c
c Is CUT1 above or below the minimum currently-accepted value?
c

      redo = .false.

      istep = int(sign(f1p0001,cut1 - sky(minimm + 1)))
c
c If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = +1, 
c then we know that at least one pixel must be deleted at the low end.
c

      jstep = (istep + 1) / 2

      if (istep .gt. 0) goto 2120
c
c Quit when SKY(MINIMM) < CUT1 <= SKY(MINIMM+1)
c

 2100 if ((istep .lt. 0) .and. (minimm .le. 0)) goto 2150
c
c If ISTEP is positive, subtract out the sky value at MINIMM+1; if 
c ISTEP is negative, add in the sky value at MINIMM.
c

      if ((sky(minimm) .le. cut1) .and. (sky(minimm + 1) .ge. cut1)) 
     &goto 2150

 2120 continue

      if(minimm + jstep .lt. 1 .or. minimm + jstep .gt. nsky) then
!!!!! Debug TPC
	if(nwarn .lt. maxwarn) then
          print *,'*** MMM range error.'
          print *,'    minimm,maximm,jstep,nsky=',minimm,maximm,jstep,nsky
          print *,'    istep,cut,cut1,cut2=',istep,cut,cut1,cut2
          print *,'    niter,maxit,minsky=',niter,maxit,minsky
          print *,'    skymed,skymod,skymn,r,sigma=',skymed,skymod,skymn
          print *,'    r,sigma=',r,sigma
          print *,'    sky=',sky(1:nsky)
	  nwarn = nwarn + 1
        endif
	sigma = -99.0
	skew = 0.0
	nused = maximm - minimm
	return
      endif

      delta = sky(minimm + jstep) - skymid
      sum = sum - (real(istep) * delta)
      sumsq = sumsq - (real(istep) * (delta ** 2))
      minimm = minimm + istep
c A change has occured

      redo = .true.
c

      goto 2100
c
c Is CUT2 above or below the current maximum?
c

 2150 continue

      istep = int(sign(f1p0001,cut2 - sky(maximm)))
c
c If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = -1, 
c then we know that we must subtract at least one pixel from the high 
c end.
c

      jstep = (istep + 1) / 2

      if (istep .lt. 0) goto 2220
c
c Quit when SKY(MAXIMM) <= CUT2 < SKY(MAXIMM+1)
c

 2200 if ((istep .gt. 0) .and. (maximm .ge. nsky)) goto 2250
c
c If ISTEP is positive, add in the sky value at MAXIMM+1; if ISTEP is 
c negative, subtract off the sky value at MAXIMM.
c

      if ((sky(maximm) .le. cut2) .and. (sky(maximm + 1) .ge. cut2)) 
     &goto 2250

 2220 delta = sky(maximm + jstep) - skymid
      sum = sum + (real(istep) * delta)
      sumsq = sumsq + (real(istep) * (delta ** 2))
      maximm = maximm + istep
c A change has occured

      redo = .true.
c

      goto 2200
c
c Compute mean and sigma (from this pass).
c

 2250 continue

      skymn = sum / dble(maximm - minimm)
      sigma = sqrt((sumsq / dble(maximm - minimm)) - (skymn ** 2))
c
c Obtain the median.  To first approximation, the median would be the
c value of the sky in the middle pixel in the sorted data (if the
c total number is odd) or the mean of the two pixels straddling
c the middle (if the total number of pixels is even).
c
c     SKYMED=0.5*(SKY((MINIMM+MAXIMM+1)/2)+SKY((MINIMM+MAXIMM)/2+1))
c
c However, this is not good enough.  If you look at the estimator for
c the mode, you will note that a tiny change in the list of sky pixels,
c just sufficient to alter the median value of the sky brightness by
c one unit, will change the estimator of the mode by three units.  We
c really want something more robust than this.  As a first attempt
c at a more robust median estimator, I propose to estimate the median
c of the distribution by the mean of the central five percent of sky
c values.  This involves considerable care to make sure you get
c a perfectly symmetric sample of pixels about the median, whether
c there is an even or an odd number of pixels within the acceptance
c interval.
c

      skymn = skymn + skymid

      skymed = 0.0
      x = 0.0
      center = real((minimm + 1) + maximm) / 2.
c

      side = (real(nint(0.05 * real(maximm - minimm))) / 2.) + 0.25

      do 2310 i = nint(center - side), nint(center + side)
      skymed = skymed + sky(i)
c

 2310 x = x + 1.

      skymed = skymed / x
      skymod = skymn
c
c If the mean is less than the mode, that means the contamination is
c slight, and the mean value is what we really want.  Note that this
c introduces a slight bias toward underestimating the sky when
c the scatter in the sky is caused by random fluctuations rather than
c by contamination, but I think this bias is negligible compared to the
c problem of contamination.
c
c If the limits have not yet stopped moving, try again.
c

      if (skymed .lt. skymn) skymod = (3. * skymed) - (2. * skymn)
c
c-----------------------------------------------------------------------
c
c Normal return.
c

      if (redo) goto 2000

      skew = (skymn - skymod) / max(1.,sigma)
      nused = maximm - minimm
c
c-----------------------------------------------------------------------
c
c An error condition has been detected.
c

      return 

 9900 sigma = -1.0
      skew = 0.0
      nused = maximm - minimm
c

      return 
      end
