       subroutine photom_error (npix,sum,nsky,csig,Nbann,zerr,noise_dn,ncoadd, SUM2err, Fcorr)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

c gain is to convert DN into electrons
        parameter (gain = 10.)   
        real*4 noise_dn, Error_median, Pi
	logical upperlim

	upperlim = .false.

	if (sum.lt.0.) upperlim = .true.

	flux = sum

	if (upperlim) then
		flux = 2. * csig   ! 2 sigma upperlim
	endif


	goto 47  ! jump to second method

c Poisson component of the source; sum is in dn units
c   gain is electrons per dn

	T0 = flux / (gain * ncoadd)
	
c integrated noise due to sky noise

	T1 = npix * (csig**2)

c uncertainty due to using a finite annulus

	T2 = ((npix * csig)**2 ) / (nsky*1.)

cc  add in quadrature

c	write (6,*) npix,nsky
c	write (6,*) T0, T1, T2

	noise_dn = sqrt (T0 + T1 + T2)

        flux_error = noise_dn / flux  ! flux (dn)

c	write (6,*) noise_dn, flux, flux_error

	zerr = 1.0857 * flux_error   ! dmag

	zerr = min (zerr, 9.999)
	zerr = max (zerr, 0.000)


ccccccccccccccccccccccccccccccccccc  Alternative method

 47	T1 = Fcorr * SUM2err  ! poisson error and aperture error

	PI = 3.14159
	Error_median = PI / 2. ! error due to using a median background instead of the mean
	T2 = Fcorr * (Error_median/Nbann*1.) * ((npix * csig)**2)   ! annulus error 

	T3 = (0.005 * flux) ** 2  ! add a 1% flat noise to account for the uncertainty in bright sources

c	write (6,*) 'error ',flux, npix, csig, Nbann
c	write (6,*) T1, T2, T3

	noise_dn = sqrt (T1 + T2 + T3)

	flux_error  = noise_dn / flux  ! flux (dn)        N/S

c	write (6,*) flux_error  

c	write (6,*) noise_dn, flux, flux_error

	zerr = 1.0857 * flux_error   ! dmag

        zerr = min (zerr, 9.999)
        zerr = max (zerr, 0.000)

        return
	end
