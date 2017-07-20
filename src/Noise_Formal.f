	subroutine noise_formal (Jsrc,J,nf,ib,nnx,nny,nx,ny,Array,MSK,Unc,
     1            x0,y0,  Rstann, Rstwid, fbits, 
     1            sky_pop)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	parameter (Valmax = 1.e5)

	real*4 Array (nnx,nny,nf,4), Unc (nnx,nny,nf,4), Rstann(4), Rstwid(4)
	integer fatal,ibit,BitSet, fbits(32), MSK(nnx,nny,nf,4)
	real*4 x0,y0

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

	sky = -999.
	sdev = 0.

	Rl = Rstann(ib)
	Rh = Rl + Rstwid(ib)

	il = int(x0 - Rh)
        ih = nint(x0 + Rh)
        jl = int(y0 - Rh)
        jh = nint(y0 + Rh)

	sum2 = 0.
        Nb = 0

	nann = 0
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

           val = Array(ii,jj,J,ib)
	   vunc = Unc (ii,jj,J,ib)
	   
	   if (vunc.gt.Valmax) goto 890

c check the MSK
	   iMaskval = MSK(ii,jj,J,ib)
           icheck = 0
           icheck = iand(iMaskval, fatal)
           if (icheck.gt.0) then
		 goto 890  ! bad pixel
	   endif


	   if (vunc.gt.0.) then
             nann = nann + 1
	     sum2 = sum2 + (vunc**2)
	   endif

 890    continue
        enddo

	if (nann.lt.5) return


	SKYrms = sqrt ( sum2 / (nann**2)   )   ! sigma in the mean
	SKY_pop = sqrt ( sum2 / (Nann * 1.) )  ! population sigma


c	write (6,*) ib,nann,SKYrms
c	if (ib.eq.4) call exit(0)


	return
	end

