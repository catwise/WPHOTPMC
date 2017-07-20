c this subroutine takes as input the a set of mags and combines them into a merged mag; inverse weighting is applied
c Wflx = weighted mean mag with Sigma Clipping to remove outliers
c Wsig1 = delta_mag derived from "unbiased weighted sample variance"
c Wsig2 = delta_mag derived from "variance of the weighted mean"
c Wsig3 = delta_mag derived from "standard error of the mean"
  
cMergeWPRO  (Nf,ib, M_M_wpro,Wproflux,eWproflux,zero(ib),mergetype)

	subroutine MergeWPRO (nf,ib,n,flux,dflux,zero,wflag,
     1    Wflx, Wsig1, Wsig2, Wsig3)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (m=999, scale = 1.)
	real*4 flux(nf),dflux(nf),zero
	real*4 z(m),dz(m),w(m),SNR
	real*4 flux2(m), w2(m)
	real*4 Wflx, Wsig1, Wsig2, Wsig3
	integer wflag, flgmerge

	Wflx = 0.
	Wsig1 = 0.
	Wsig2 = 0.
	Wsig3 = 0.

	flgmerge = 0

	F0 = zero

	wtsnr = 0.
	wt = 0.
	nS = 0
	N_M = 0   !  demand the source to have at least SNR>3
	M_M = 0   !  N of M statistics

	do 500 j=1,n

		error = dflux(j)
		if (error.le.0.) goto 500
		if (error.gt.1.e7) goto 500

		nS = nS + 1

		z(ns) = flux(j)
		dz (ns) = dflux(j)

		SNR = flux(j) / dflux(j)

		SNR = min (SNR, 200.)  ! takes care of ridiculously bright things
		SNR = max (SNR, 1.)    ! takes care of negative flux (assume 2sigma)

		w(nS) =  SNR  ** 2  ! this is the weight based on the S/N
		wtsnr = wtsnr + w(nS)
		w(nS) =  1 . / (error**2)   ! this is the weight based on the flux uncertainty

		if (wflag.eq.0) w(nS) = 1.

		wt = wt + w(nS)

 500	continue

	if (nS.le.0) then
		! not enough sources to measure
		return
	endif

	wtvar = 0.
	if (wflag.gt.0) then
	 if (wt.gt.0.) then
	   wtvar = 1. / sqrt (wtsnr) ! variance of the weighted mean using SNR
	 endif
	endif

	if (nS.eq.1) then

		Ave = z(1)
		zmag = 99.
		if (ave.gt.0.) then
		  zmag = zero - (2.5 * log10(ave*scale))
		endif
		Wflx = z(1)
        	Wsig1 = 0.
        	Wsig2 = 0.
        	Wsig3 = 0.

	else

	call  SigClipWmean (nS,z,w, flux2, w2, AVE,SDEV,sdev2, N_M_2)

c	write (6,*) nS,AVE,SDEV,sdev2

	SNR = ave / sdev2  ! using the unbiased weighted sample variance

	if (wflag.gt.0) then
	  SNR2 = 1./wtvar
	else
	  SNR2 = SNR
	endif

        sdev3 = sdev2 / sqrt (n*1.)  ! standard dev. of the mean
        SNR3 = ave/sdev3

	Wflx = Ave
	Wsig1 = SNR   !  zdmag
	Wsig2 = SNR2  !  zdmag2
	Wsig3 = SNR3  !  zdmag3

	endif

	return
	end
