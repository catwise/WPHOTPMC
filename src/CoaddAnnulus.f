	subroutine CoaddAnnulus (Jsrc,ib,ncoaddsize,nx,ny,Array,Unc,iCOADD,
     1            x0,y0,  Rstann, Rstwid, nann, sky, sdev, 
     1            Nb, Rl, Rh, pscale,cscale, Fcorr, BGmax, SKY_pop )

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

c readfits.f returns pixel values of -9991 for Nans, so we need to exclude those
	parameter (nmax = 9999, Valmax = 1.e5, Valmin = -9900.0)

	real*4 Array (ncoaddsize,ncoaddsize,4), Unc (ncoaddsize,ncoaddsize,4)
	real*4 Rstann(4), Rstwid(4)
	real*4 x0,y0, pscale,cscale
	real*4 z (nmax), z2 (nmax), Fcorr, BGmax
	integer iCOADD(ncoaddsize,ncoaddsize,4), imask
	logical debug, doCOADD
c
c vvvvvvv stuff for MMM
        real*4 hival, dammmmean, dammmmed, dammmmode, dammmsigma, skew
        integer*4 nused, nitmmm
c
        include 'jwfcom.f'             ! TPC, to get bgtype(4)
c
	integer nmmmwarn, maxmmmwarn
	data nmmmwarn/0/, maxmmmwarn/50/
c

c       sky_pop = 9.9  ! JWF B21114 dbg
c	debug = .true.
	debug = .false.

	doCOADD = .true.


c	write (6,*) 'CoaddAnnulus',ncoaddsize,Array(3601,2709,4),Unc(3601,2709,4)
c       print *,'CoaddAnnulus: Jsrc =',jsrc   ! JWF B21112 dbg


cc   annuluar background

	mode = 0   ! do not use iCOADD mask
c	mode = 1   ! use mask

cccccccccccccccc

	pscale_factor = pscale / cscale

	sky = -999.
	sdev = 0.

	if (Rl.le.0.) then
	  Rl = Rstann(ib)*pscale_factor  !  convert frame value to coadd value
	endif

	Rh = Rl + (Rstwid(ib)*pscale_factor)  !  convert frame value to coadd value
	


	il = int(x0 - Rh)
        ih = nint(x0 + Rh)
        jl = int(y0 - Rh)
        jh = nint(y0 + Rh)

c	write (6,*) ib,Rl,Rh ,'  pixels '
c	write (6,*) 'annulus ',il,ih,jl,jh


 722	nann = 0
        do jj=jl,jh
           dy2 = (jj-y0)**2
        do 890 ii=il,ih
           dx2 = (ii-x0)**2.
           dr = sqrt (dx2+dy2)
           if (dr.lt.rl) goto 890
           if (dr.ge.rh) goto 890

           if (ii.lt.1) goto 890
           if (jj.lt.1) goto 890
           if (ii.gt.nx) goto 890
           if (jj.gt.ny) goto 890
	   if (nann.ge.nmax) goto 890

           val = Array(ii,jj,ib)
	   vunc = Unc (ii,jj,ib)
	   imask = 1

c	write (6,*) ii,jj,dr,Val,vunc

	   if (vunc.gt.Valmax .or. vunc .lt. 0) goto 890
	   if(val .lt. valmin) goto 890

	   imask = iCOADD(ii,jj,ib)
	   if (mode.eq.0) imask=0

c	if (ib.eq.4)  write (6,*) ii,jj,ib,imask,Val,vunc

           if ((imask.eq.0).or.(imask.eq.Jsrc)) then
            nann = nann + 1
            z (nann) = val
            z2 (nann) = vunc
           endif

 890    continue
        enddo

c	write (6,*) nann
c       print *,'nann =',nann  ! JWF B21112 dbg
c	if (ib.eq.4) call exit(0)

	if (nann.lt.29) then
		if (mode.eq.1) then ! try again
			mode = 0
			goto 722
		endif
		return
	endif


ccccccccccccccccccccccccccccccccccccccccc  Compute the RMS using the UNC image

	sum2 = 0.
	Nb = 0
	do jj=1,nann
		if (z2(jj).gt.0.) then
		  sum2 = sum2 + (z2(jj)**2)
		  Nb=Nb+1
		endif
	enddo
c       print *,'Jsrc, sum2, Nb:', jsrc, sum2, nb  ! JWF B21112 dbg

 	SKY_pop = sqrt (sum2 / (Nb**1)   )    ! this is the standard deviation in the population
c	SKY_pop = 1.0          !  JWF B21114 dbg
c       print *,'SKY_pop:', sky_pop  !  JWF B21113 dbg
	SKYrms = sqrt ( sum2 / (Nb**2)   )    ! this is the standard deviation in the mean

c	SKY_pop = sqrt(Fcorr) * SKY_pop
c	SKYrms = sqrt(Fcorr) * SKYrms
c	sdev = SKYrms  ! use the Unc image for the uncertainty


c	write (6,*) nb,skyrms
c	call exit(0)


c inflate if using a median

c	SKYrms = 1.57 * SKYrms 



cccccccc   compute the population statistics; use the results as input to histogram statistics
cccccccc    this will be the median sky

	nit = 2
c	if (mode.eq.0) nit=3

	!  N-sigma rejection, median background

        zlow = -99.
	zhigh = BGmax
        if (zhigh.le.0.) zhigh = 4999.*3.  ! 

	!   ave == trimmed average using low-side sigma
	!   sdev = low-side sigma
	!   hsig = histogram sigma (average of low and high-side sigma)
	!   zmed == histogram median

        call VecStats (nann,z,m,z2,nit,zlow,zhigh,
     1      ave,sdev,xmed, SKY_pop, hsig, doCoadd)  ! sdev is the low-side sigma; ave is the trimmed average using low-side sigma


	! Conventional background
c	sky = xmed
	sdev = hsig   ! use the histogram sigma;  this is the uncertainty per pixel

!!! new as of Feb 08, 2010
	sky = ave   ! trimmed average using low-side sigma; really the mode now
	! sdev is the low-side sigma
!! new as of Jan 19, 2011 -- this is the mode and its sigma


c ================ !!! TPC
c         New BG estimation using MMM
	  if(BGType(ib) .ne. 1) then
	    call sort(nann,z)
	    hival = 1.e20
	    call mmm(z,nann,hival, dammmmean, dammmmed, dammmmode,
     1               dammmsigma, skew, nused, nitmmm)
	    if (dammmsigma .gt. 0.0) then
	      if(BGType(ib) .eq. 2) then
		sky = dammmmode
              else if(BGType(ib) .eq. 3) then
                sky = dammmmed
              else if(BGType(ib) .eq. 4) then
                sky = dammmmean
	      endif
            else
                !!! error in mmm, continue to use sky = mode
	      if(dammmsigma .lt. -98 .and. nmmmwarn .lt. maxmmmwarn) then
		print *,'=== Warning (coaddannulus): MMM failed.'
		print *,'    x,y,band=',x0,y0,ib
		nmmmwarn = nmmmwarn + 1
	      end if
	    endif
	  else ! bgtype(ib) == 1
              !!! continue to use sky = mode
	  endif
c         else bgtype = 1 =>  retain current definition of 'sky' = mode
c ================ !!! TPC



c	write (113,*) sky, sky-sdev,sky+sdev


	if (debug) then

       write (6,'(a,2f8.3,i7,2f11.5,a,f10.4)')
     1     'VecStat annular sky (R-,R+,N,ave,sky,rms): ',
     1      rl,rh,m,ave,sky,' +-',sdev

	endif


	return
	end


