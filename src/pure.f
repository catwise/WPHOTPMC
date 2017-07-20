	subroutine  pure (ib,zero,fdn,sigdn)

        implicit real*4 (a-h,o-z)
	implicit integer (i-n)

	real*4 fdn,sigdn
	real*4 zero, f0(4)
	real*4 sig_0(4), D(4)

	data f0/311.,168.,30.,7.67/  ! these are the nominal zero points

	data sig_0 /12.6, 10, 19.7, 38.8/
	data D/101.6, 101.1, 1110, 2168/

c	write (6,*) fdn

c convert fdn to uJy

	zmag = zero - (2.5*log10(fdn))

	fluxJy = f0(ib) * 10**(-zmag/2.5)

	fluxuJY = fluxJy * 1.e6

	sig_pixel = sig_0(ib) * sqrt( 1 + (fluxuJY / D(ib)) )  ! uJy

	sigJy = sig_pixel / 1.e6

	zmagsig = -2.5 * log10( sigJy / f0(ib) )

	sigDN = 10**((zero - zmagsig)/2.5) 

c	write (6,*) fdn,fluxuJY,sig_pixel,sigDN

	return
	end
