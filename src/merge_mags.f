c this subroutine takes as input the a set of mags and combines them into a merged mag; inverse weighting is applied
c Wmag = weighted mean mag with Sigma Clipping to remove outliers
c Wsig1 = delta_mag derived from "unbiased weighted sample variance"
c Wsig2 = delta_mag derived from "variance of the weighted mean"
c Wsig3 = delta_mag derived from "standard error of the mean"
  
	subroutine MergeMags (n,mag,dmag,flg,zero,wflag,
     1    Wmag, Wsig1, Wsig2, Wsig3, nS, flgmerge)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (m=9999, scale = 10.)
	real*4 mag(n),dmag(n),zero
	real*4 flux (m),dflux(m),w(m),SNR
	real*4 flux2(m), w2(m)
	real*4 Wmag, Wsig1, Wsig2, Wsig3
	integer wflag, flgmerge, flg(n)


	Wmag = 99.
	Wsig1 = 0.
	Wsig2 = 0.
	Wsig3 = 0.

	flgmerge = 0
c convert mags to fluxes

	F0 = zero

	wt = 0.
	nS = 0

	do 500 j=1,n
		if (mag(j).gt.50.) goto 500
		if (dmag(j).le.0.) goto 500

		flgmerge = max (flgmerge , flg(j) )

		nS = nS + 1

		flux(ns) = 10**((zero-mag(j))/2.5)  / scale
c		flux(nS) = F0 * (10**(-mag(j)/2.5) )    ! Jansky

		if (dmag(j).le.0.) then
		  SNR = 100.
		else
		  SNR = 1.0857 / dmag(j)
		  SNR = min (SNR, 100.)
		  SNR = max (SNR, 0.1)
		endif

		dflux(nS) = flux(nS) / SNR

		w(nS) = 1. / (dflux(nS) ** 2)
		if (wflag.eq.0) w(nS) = 1.

		wt = wt + w(nS)

c	write (6,*) mag(j),flux(j),dflux(j),w(j)

 500	continue

	if (nS.le.0) then
		! not enough sources to measure
		return
	endif

	wtvar = 1 ! If wt <= 0
	if (wflag.gt.0) then
	 if (wt.gt.0.) then
	   wtvar = 1. / sqrt (wt)  ! variance of the weighted mean
	 endif
	endif

c	call MOMENT(flux,N,AVE,SDEV)
c	SNR = ave / sdev  ! using the unbiased weighted sample variance
c        zdmag = 1.0857 / SNR
c	write (6,*) 'no-weighting, delta_mag = ',zdmag


	if (nS.eq.1) then

		Ave = flux(1)
		zmag = zero - (2.5 * log10(ave*scale))
		Wmag = zmag
        	Wsig1 = 0.
        	Wsig2 = 0.
        	Wsig3 = 0.

	else

c	call WMOMENT(flux,w,nS,AVE,SDEV,sdev2)
	call  SigClipWmean (nS,flux,w, flux2, w2, AVE,SDEV,sdev2,n_m_2)


c	write (6,*) AVE,SDEV,sdev2

	SNR = ave / sdev2  ! using the unbiased weighted sample variance

	zdmag = 1.0857 / SNR 

	zmag = zero - (2.5 * log10(ave*scale))
c	zmag = -2.5 * log10((ave) / F0)

	if (wflag.gt.0) then
	  SNR2 = ave / wtvar ! using the variance of the weighted mean
	  zdmag2 = 1.0857 / SNR2 ! using the variance of the weighted mean
	else
	  zdmag2 = zdmag
	endif

        sdev3 = sdev2 / sqrt (n*1.)  ! standard dev. of the mean
        SNR3 = ave/sdev3
        zdmag3 = 1.0857 / SNR3

c	write (6,'(i5,4f7.3)') n,zmag,zdmag,zdmag2,zdmag3

	Wmag = zmag
	Wsig1 = zdmag
	Wsig2 = zdmag2
	Wsig3 = zdmag3

	endif

	Wmag = min (Wmag, 99.)
	Wsig1 = min (Wsig1, 9.99)
	Wsig2 = min (Wsig2, 9.99)
	Wsig3 = min (Wsig3, 9.99)
	Wsig1 = max (Wsig1, 0.0)
	Wsig2 = max (Wsig2, 0.0)
	Wsig3 = max (Wsig3, 0.0)

	return
	end
