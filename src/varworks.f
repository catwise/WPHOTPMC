	subroutine varworks (Jsrc,nf,ib,M_M_wpro,Wproflux,eWproflux,zero,mergetype,WproJD,
     1      NSUMmax,WPROmagR)
        implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        real*4 Wproflux(nf),eWproflux(nf),zero, WPROmagR(NSUMmax,4,3)
        real*4 Wflx, Wsig1, Wsig2, Wsig3
	real*4  JD(999)
	real*8  WproJD(nf), x8, dx8, xmax8
        integer  mergetype

	if ((Jsrc.eq.-1660).and.(ib.eq.1)) then
                write (23,*) 'here ',zero
        do kk=1,M_M_wpro

		zz = zero - (2.5*log10(Wproflux(kk)))

                write (23,'(i5,f13.6,f12.3,f8.3)') kk,WproJD(kk),Wproflux(kk),zz
        enddo
	write (23,*) ' '

c        call exit(0)
        endif


	call MergeWPRO  (Nf,ib, M_M_wpro,Wproflux,eWproflux,zero,mergetype,
     1    Wflx, WSNR1, WSNR2, WSNR3)

	WPROmagR (Jsrc,ib,1) = 99.99
        WPROmagR (Jsrc,ib,2) = 0.
        WPROmagR (Jsrc,ib,3) = 0.



        if (Wflx.gt.0.) then

                  WPROmagR (Jsrc,ib,1) = zero - (2.5*log10(Wflx))

                  if (WSNR1.gt.0.) WPROmagR (Jsrc,ib,2) = 1.0857 / WSNR1   ! dmag w/ unbiased weighted mean
                  if (WSNR3.gt.0.) WPROmagR (Jsrc,ib,3) = 1.0857 / WSNR3 ! dmag w/ "standard error of the mean"

                  WPROmagR (Jsrc,ib,2) = min (WPROmagR (Jsrc,ib,2), 99.)
                  WPROmagR (Jsrc,ib,3) = min (WPROmagR (Jsrc,ib,3), 99.)

        else
                  WPROmagR (Jsrc,ib,1) = 99.99
                  WPROmagR (Jsrc,ib,2) = 0.
                  WPROmagR (Jsrc,ib,3) = 0.
        endif


	return
	end
