	subroutine HistStats (n,vec,m,vout,niter,
     1     vlow,vhigh,vwidth,
     1     z50,hsig, zmode)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (max = 50000, nsigma = 3)

	real*4 vec(n),vout(n)
	integer count (max)
	logical debug

	debug = .true.
	debug = .false.

	z16 = 0
	z84 = 0

	it = 0
 47	ztot = (vhigh - vlow) / vwidth
	it = it + 1

	if (nint(ztot).ge.max) then
c	  write (74,*) '** warning -- HistStat bin width too small ',it
	  vwidth = vwidth * 2.
	  goto 47
	endif

	ave = -9999.
        sdev = -9999.
        zmed = -9999.

	do 1000 Jit = 1,niter

	do j=1,max
                count(j) = 0
        enddo

	m = 0
	do 50  J = 1, n
	   if (vec(J).lt.vlow) goto 50
	   if (vec(J).ge.vhigh) goto 50

	   m=m+1
	   vout(m)=vec(J)

	  dval = vec(J) - vlow
          idex = 1 + int(dval / vwidth)
          if (idex.ge.max) goto 50

           count(idex) = count(idex) + 1
	
 50	continue

	zmax = vlow + (max * vwidth) - (vwidth / 2.0)

	if (debug) write (6,*) 'histo pass ',jit,m,zmax


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if (m.lt.3) then
		z84 = -999.
                z50 = -999.
                z16 = -999.
                return
        endif


c mode

	if (Jit.eq.niter) then
	  countmax = 0.
	  do j=1,max
	   zme = vlow + (j * vwidth) - (vwidth / 2.0)
	   if (debug) write (114,*) zme,count(j)

	  if (count(j).gt.countmax) then
		countmax = count(j)
		zmode = zme
	  endif
	  enddo
	endif

	

c 84% quartile
	i84 = nint( 0.84 * m * 1. )
        isum = 0
	isum_prev = 0
        do j=1,max

                isum = isum + count(j)
                if (isum.eq.i84) then
                        z84 = vlow + (j * vwidth) - (vwidth / 2.0)
                        goto 53
                endif

                if (isum.gt.i84) then
c                 interpolate
                  idel = isum - isum_prev
                  ilow = i84 - isum_prev
                  ratio = (ilow*1.) / (idel*1.)
                  del = vwidth * ratio

                  z84 = (del + vlow) + ((j-1) * vwidth) - (vwidth / 2.0)
                  goto 53
                endif

                isum_prev = isum
        enddo
 53     i=1

	if (debug) write (6,*) 'hist 84% ',z84
	
c 50% (median) quartile

	i50 = nint( 0.50 * m * 1. )
        isum = 0
        do j=1,max

                isum = isum + count(j)
                if (isum.eq.i50) then
                        z50 = vlow + (j * vwidth) - (vwidth / 2.0)
                        goto 54
                endif

                if (isum.gt.i50) then
c                 interpolate

                  idel = isum - isum_prev
                  ilow = i50 - isum_prev
                  ratio = (ilow*1.) / (idel*1.)
                  del = vwidth * ratio

                  z50 = (del + vlow) + ((j-1) * vwidth) - (vwidth / 2.0)
                  goto 54
                endif

                isum_prev = isum
        enddo
 54     i=1

	if (debug) write (6,*) 'hist 50% ',z50

c 16%  quartile

        i16 = nint( 0.16 * m * 1. )
        isum = 0
        do j=1,max

                isum = isum + count(j)
                if (isum.eq.i16) then
                        z16 = vlow + (j * vwidth) - (vwidth / 2.0)
                        goto 55
                endif

                if (isum.gt.i16) then
c                 interpolate
                  idel = isum - isum_prev
                  ilow = i16 - isum_prev
                  ratio = (ilow*1.) / (idel*1.)
                  del = vwidth * ratio

                  z16 = (del + vlow) + ((j-1) * vwidth) - (vwidth / 2.0)
                  goto 55
                endif

                isum_prev = isum
        enddo
 55     i=1

	if (debug) write (6,*) 'hist 16% ',z16

	
c done with histogram; compute width == sigma

	hsig = (z84 - z16)/2.0

	if (debug ) write (6,*) Jit, z16, z50, z84, hsig

	 vlow  = z50 - (hsig * nsigma * 1.)
         vhigh = z50 + (hsig * nsigma * 1.)

 1000	continue

 1001	return

	end

