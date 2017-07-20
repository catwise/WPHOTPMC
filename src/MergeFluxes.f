c this subroutine takes as input the a set of mags and combines them into a merged mag; inverse weighting is applied
c Wflx = weighted mean mag with Sigma Clipping to remove outliers
c Wsig1 = delta_mag derived from "unbiased weighted sample variance"
c Wsig2 = delta_mag derived from "variance of the weighted mean"
c Wsig3 = delta_mag derived from "standard error of the mean"
  
	subroutine MergeFluxes (nf,ib,Jsrc,Map,n,flux,dflux,flg,zero,wflag,
     1    Wflx, Wsig1, Wsig2, Wsig3, nS, flgmerge, N_M, M_M, N_M_2)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (m=999, scale = 1.)
	real*4 flux(nf,19),dflux(nf,19),zero
	real*4 z(m),dz(m),w(m),SNR
	real*4 flux2(m), w2(m)
	real*4 Wflx, Wsig1, Wsig2, Wsig3
	integer wflag, flgmerge, flg(nf,19)
	integer N_M, M_M, N_M_2

c	write (6,*) 'inside'
c	do j=1,n
c		write (6,*) j,flux(j,Map),dflux(j,Map)
c	enddo
c	write (6,*) '***'
	 

	Wflx = 0.
	Wsig1 = 0.
	Wsig2 = 0.
	Wsig3 = 0.

	flgmerge = 0
c convert mags to fluxes

	F0 = zero

	wtsnr = 0.
	wt = 0.
	nS = 0
	N_M = 0   !  demand the source to have at least SNR>3
	M_M = 0   !  N of M statistics

	do 500 j=1,n

		error = dflux(j,Map)
		if (error.le.0.) goto 500
		if (error.gt.1.e7) goto 500

		flgmerge = max (flgmerge , flg(j,Map) )
		
		nS = nS + 1

		z(ns) = flux(j,Map)
		dz (ns) = dflux(j,Map)

		SNR = flux(j,Map) / dflux(j,Map)

		if (SNR .ge. 3.0)  N_M = N_M + 1    ! for N out of M statistics

		SNR = min (SNR, 200.)  ! takes care of ridiculously bright things
		SNR = max (SNR, 1.)    ! takes care of negative flux (assume 2sigma)

		w(nS) =  SNR  ** 2  ! this is the weight based on the S/N
		wtsnr = wtsnr + w(nS)

		w(nS) =  1 . / (error**2)   ! this is the weight based on the flux uncertainty

		if (wflag.eq.0) w(nS) = 1.

		wt = wt + w(nS)


 500	continue

	M_M = nS
	N_M_2 = 0

	if (nS.le.0) then
		! not enough sources to measure
		return
	endif

	wtvar = 0.
	if (wflag.gt.0) then
	 if (wt.gt.0.) then
c	   wtvar = 1. / sqrt (wt)  ! variance of the weighted mean
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

c	zdmag = 1.0857 / SNR 
c	zmag = zero - (2.5 * log10(ave*scale))
c	zmag = -2.5 * log10((ave) / F0)

c	if ((ib.eq.3).and.(Jsrc.eq.1000)) then
c	  do jj=1,nS
c	  write (6,*) z(jj),dz(jj),w(jj)
c	  enddo
c	  write (6,*) ave,sdev,sdev2
c	  write (6,*) SNR,1./wtvar 
c	endif


	if (wflag.gt.0) then
	  SNR2 = 1./wtvar
c	  SNR2 = ave / wtvar ! using the variance of the weighted mean
c	  zdmag2 = 1.0857 / SNR2 ! using the variance of the weighted mean
	else
	  SNR2 = SNR
c	  zdmag2 = zdmag
	endif

        sdev3 = sdev2 / sqrt (n*1.)  ! standard dev. of the mean
        SNR3 = ave/sdev3
c        zdmag3 = 1.0857 / SNR3

c	write (6,'(i5,4f7.3)') n,zmag,zdmag,zdmag2,zdmag3

	Wflx = Ave
	Wsig1 = SNR   !  zdmag
	Wsig2 = SNR2  !  zdmag2
	Wsig3 = SNR3  !  zdmag3

	endif

c	Wflx = min (Wflx, 99.)
c	Wsig1 = min (Wsig1, 9.99)
c	Wsig2 = min (Wsig2, 9.99)
c	Wsig3 = min (Wsig3, 9.99)
c	Wsig1 = max (Wsig1, 0.0)
c	Wsig2 = max (Wsig2, 0.0)
c	Wsig3 = max (Wsig3, 0.0)

	return
	end
