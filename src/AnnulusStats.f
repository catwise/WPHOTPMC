c (MSK is AKA iMask) -- TPC
	subroutine AnnulusStats (Jsrc,J,nf,nfpix,ib,nnx,nny,nx,ny,Array,MASK,MSK,Unc,
     1            x0,y0,  Rstann, Rstwid, nann, sky, sdev, Nbann, Rl, Rh, mode, SPIT, fbits, zero,
     1            skyRMS, hsig, SKY_pop, BGmax, Rcrater)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

c readfits.f returns pixel values of -9991 for Nans, so we need to exclude those
	parameter (nmax = 9999, Valmax = 1.e5, Valmin = -9900.0, badpixlim = -1.0e+20)

c	integer MASK (nnx,nny,nf,4)
	integer MASK (1,1,1) ! Mask being dummied out since it's never assigned anything -- TPC
	integer imask

	real*4 Array (nnx,nny,nfpix), Unc (nnx,nny,nfpix), Rstann(4), Rstwid(4)
	integer fatal,ibit,BitSet, fbits(32), MSK(nnx,nny,nfpix)
	real*4 x0,y0, zero(4), BGmax
	real*4 z (nmax), z2 (nmax), zpure (nmax)
	logical debug, SPIT, doCOADD
c
c vvvvvvv stuff for MMM
        real*4 hival, dammmmean, dammmmed, dammmmode, dammmsigma, skew
        integer*4 nused, nitmmm
c
        include 'jwfcom.f'             ! TPC, to get bgtype(4)
c

	debug = .true.
	debug = .false.

	doCOADD = .false.
c	doCOADD = .true.  ! Temporary

	mode = 0 ! Don't check mask() for star contamination, because mask() hasn't been set

	fatal = 0
        do I = 1,32
                ibit = fbits(I)
                if (ibit.ge.0) then
                  fatal = fatal + (2**ibit)
                else
                        !quit
                        goto 74
                endif
        enddo
 74     I = 0

cc   annuluar background

	sky = -999.99
	sdev = 0.

	Rl = Rstann(ib)
	Rh = Rl + Rstwid(ib)
	RRmax = nx / 2.

cccc  check Crater radius   !  Jan 25, 2011
cc  when Rcrater > Rl/4, use a set large annulus of:  Rl*10

	if (Rcrater.gt.Rl/4.) then
	  Rinner = Rl*10.

c	  Crater_Scale = 10.
c	Crater_Scale = 32.  ! TEMP
c	  Rinner = (Crater_Scale * Rcrater) + Rl

	  Rinner = min (RRmax-1., Rinner)
	  call NannMax (rinner, router)
	  Rl = Rinner
	  Rh = router
c	  write (6,*) 'here; new Rinner/Router: ',Rl,Rh
	endif

	Rl2 = Rl**2
	Rh2 = Rh**2

c	if (Jsrc.lt.5000) return  ! temporary !!
c	write (6,*) 'here jsrc,j,nf,ib,x0,y0 ',Jsrc,J,nf,ib,x0,y0
c	write (6,*) rl,rh

	if (debug) then
	write (13,*) x0,y0, rl,rh
	write (13,*) x0,y0+rl
	write (13,*) x0,y0+rh
	  write (13,*) 'mark ',x0,y0,' el red ',Rl*2.,Rl*2.,' 0'
	  write (13,*) 'mark ',x0,y0,' el red ',Rh*2.,Rh*2.,' 0'
	endif

	il = int(x0 - Rh)
        ih = nint(x0 + Rh)
        jl = int(y0 - Rh)
        jh = nint(y0 + Rh)

	nann = 0
        do jj=jl,jh
           dy2 = (jj-y0)**2
        do 890 ii=il,ih
           dx2 = (ii-x0)**2.
	   dr2 = dx2+dy2
           if (dr2.lt.rl2) goto 890
           if (dr2.ge.rh2) goto 890

           if (ii.lt.1) goto 890
           if (jj.lt.1) goto 890
           if (ii.gt.nx) goto 890
           if (jj.gt.ny) goto 890
	   if (nann.ge.nmax) goto 890

           val = Array(ii,jj,J)
	   vunc = Unc (ii,jj,J)
	   
c	   call pure (ib,zero(ib),val,vpure)  ! this the pure poisson value

	   if (val.le.badpixlim) goto 890  ! Sub-region is partially off the frame
	   if (val.lt.Valmin) goto 890  ! Nan in the image
	   if (vunc.gt.Valmax .or. vunc .lt. 0) goto 890


c check the MSK
	   iMaskval = MSK(ii,jj,J)
           icheck = 0
           icheck = iand(iMaskval, fatal)
           if (icheck.gt.0) then
		 goto 890  ! bad pixel
	   endif


c star contamination
	   imask = 0 ! Belt and suspenders; never look at mask()
c	   imask = MASK(ii,jj,J) ! Not fully allocated or used any more

	   if (mode.eq.0) imask=0 ! Always this

	   if ((imask.eq.0).or.(imask.eq.Jsrc)) then

c	write (6,*) nann
            nann = nann + 1

            z (nann) = val
	    z2 (nann) = vunc
c   	    zpure (nann) = vpure

c	write (6,*) ib,zero(ib),val,vpure

c	write (14,*) nann,val,vunc,nmax  ! temporary

	   endif


 890    continue
        enddo


c	if (ib.eq.3) write (6,*)'DEB ',nann

	if (nann.lt.29) return


cccccccccc
c	call MDIAN1(zpure,Nann,SigPure)   ! this is temporary

ccccccccccccccccccccccccccccccccccccccccc  Compute the RMS using the UNC image

	sum2 = 0.
	Nb = 0
	do jj=1,nann
		if (z2(jj).gt.0.) then
		  sum2 = sum2 + (z2(jj)**2)
		  Nb=Nb+1
		endif
	enddo

	Nbann = Nb
	SKY_pop = sqrt ( sum2 / (Nb**1)   )

	SKYrms = sqrt ( sum2 / (Nb**2)   )   ! sigma in the mean
	sdev = SKY_pop  ! use the Unc image for the uncertainty

c	if (ib.eq.3) write (6,*)  '# sigma_pop, sigma_mean ',SKY_pop,SKYrms

c	write (6,*)  '# sigma_pop, sigma_mean ',SKY_pop,SKYrms  ! TEMP

c	write (6,*) 'sigpure = ',SigPure
c	write (6,*) 'Sky_pop = ',SKY_pop
c	write (6,*) 'Skyrms = ',SKYrms
c	call exit (0)

c	write (88,*) '# sigma_pop, sigma_mean ',SKY_pop,SKYrms
c	write (88,*) '# ',nb,skyrms
c	call exit(0)


c inflate if using a median

c	SKYrms = 1.57 * SKYrms 



cccccccc   compute the population statistics; use the results as input to histogram statistics
cccccccc    this will be the median sky

	nit = 3

	!  N-sigma rejection, median background
	!   ave == trimmed average using low-side sigma
	!   sdev = low-side sigma
	!   hsig = histogram sigma (average of low and high-side sigma)
	!   zmed == histogram median

c update, now using MODE and it's sigma
	! hsig is the MODE sigma
	! zmed is the histogram mode

        zlow = -99.
        zhigh = BGmax
	
	
c	 if (ib.eq.3) write (6,*) zhigh

	if (zhigh.le.0.) BGmax = 5000.*3.
        call VecStats (nann,z,m,z2,nit,zlow,zhigh,
     1      ave,sdev,xmed, SKY_pop, hsig, doCOADD)  ! sdev is the low-side sigma; ave = trimmed average

c	write (6,*) 'VECSTAT ',m,ave,sdev,xmed, hsig  ! TEMP

c	if (ib.eq.3) write (6,*) 'VECSTAT ',m,ave,sdev,xmed, hsig

	! Conventional
c	sky = xmed  ! median sky value
c	sdev = SKY_pop  ! use the Unc image for the uncertainty

	if ( (hsig.gt.0.).and.(hsig.lt.999.) ) then
		sdev = hsig  ! use the annulus histogram RMS for the uncertainty (sigma per pix)
	endif

	!!! new as of Feb 08, 2010
	  if (sdev.gt.0.) then
           sky = ave   ! trimmed average using low-side sigma <<< actually the mode
           ! sdev is the low-side sigma or Hsig

	   !  as of Jan 19, 2011, this is now the mode and its sigma

	  else
		sky = xmed  ! median sky value
		sdev = SKY_pop  ! use the Unc image for the uncertainty
	  endif

c ================ !!! TPC
c         New BG estimation using MMM
	  if(BGType(ib) .ne. 1) then

c	    print *,'--- ib=',ib,', j=',j,', il,ih,jl,jh=',il,ih,jl,jh,', nnx,nny=',nnx,nny,
c     1              ', nx,ny=',nx,ny !!! ,', z(1:',nann,') = ', z(1:nann) !!! TPC
c	    il1 = max(1,il)
c	    ih1 = min(nx,ih)
c	    jl1 = max(1,jl)
c	    jh1 = min(ny,jh)
c	    do jj=jl1,jh1
c	      print *,jj,':',Array(il1:ih1,jj,j)
c	    end do

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
	    endif
	  else ! bgtype(ib) == 1
              !!! continue to use sky = mode
	  endif
c         else bgtype = 1 =>  retain current definition of 'sky' = mode
c ================ !!! TPC

c	if (ib.eq.3) write (6,*) 'sky, sdev ',sky,sdev

	if (.not.SPIT) return

	return



c !!!!!!!!!!!!!!!!!!!!
c Nothing after here matters (except the nannmax subroutine)
c !!!!!!!!!!!!!!!!!!!!


	if (debug) then

       write (6,'(a,2f8.3,i7,2f11.5,a,f10.4)')
     1     'VecStat annular sky (R-,R+,N,ave,sky,rms): ',
     1      rl,rh,m,ave,sky,' +-',sdev


          blow = 0.
          bhigh = 999.
          bwidth = 2.5
          call binit (blow,bhigh,bwidth,z2,m)

	   write (15,*) sky-sdev,sky+sdev,'  VecStat'

	endif

c	if (mode.eq.0) return

	nit = 2
	vwidth = SKYrms / 10.
	zlow =  xmed - (10. * SKYrms)
	zhigh = xmed + (10. * SKYrms)

	if (debug) write (6,*) 'histogram limits ',zlow,zhigh

	zmode = 0.
	call HistStats (nann,z,m,z2,nit,zlow,zhigh,
     1     vwidth, z50,hsig, zmode)

	if (debug) then
	write (6,'(a,2f8.3,i7,f10.4,a,2f10.4)') 
     1    'HistStat annular sky histo (R-,R+,N,sky,hsig, mode): ',
     1     rl,rh,m,z50,' +-',hsig, zmode

c	  write (15,*) z50-hsig,z50+hsig,'  HistStat'

	endif

	if ( (hsig.gt.0.).and.(hsig.lt.999.) ) then
	  sky = z50
	  sdev = hsig
	  nann = m
	else
		sdev = SKYrms
	endif

c	if (ib.eq.1) then
c		write (77,*) skyRMS, sdev, sdev/skyRMS
c	endif

	! use SKYrms for the sky uncertainty

	 sdev = SKYrms

c	if (zmode.gt.0.) then
c		sky = zmode
c	endif

	return
	end


	subroutine binit (blow,bhigh,bwidth,data,N)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (max = 9000, zfloor = -947.)
	real*4 data(N)
	integer count (max)

	do j=1,max
		count(j) = 0
	enddo

	ntot = 0
	iff = 0
	do 50 j=1,N
		val = data(j)
		if (val.eq.zfloor) iff=1
		if (val.lt.blow) goto 50
		if (val.gt.bhigh) goto 50

		dval = val - blow
		idex = 1 + int(dval / bwidth)	
		if (idex.gt.max) goto 50

		ntot = ntot + 1
		count(idex) = count(idex) + 1

 50	continue

	isum = 0
	do j=1,max
	   zbin = blow + (j * bwidth) - (bwidth / 2.0)
	   write (14,*) zbin,count(j)
	enddo
 51 	i=1



 66	i=1

	return
	end			



	subroutine NannMax (rinner, router)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (nmax = 5000)  ! maximum number of pixels in the annulus

	pi = 3.14159


	router = sqrt((rinner**2) + ((nmax * 1.)/pi))

	return
	end
