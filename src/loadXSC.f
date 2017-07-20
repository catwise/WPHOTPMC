c  these 2mass values are "corrected" for the angular resolution of the beam (see below)
c XSCshape (nsrc,1) = rk20fe   # corrected
c XSCshape (nsrc,2) = ba       # corrected
c XSCshape (nsrc,3) = pa      
c XSCshape (nsrc,4) = rk20fe   # corrected for wise W4 
c XSCshape (nsrc,5) = ba       # corrected for wise W4

	subroutine loadXSC (nfi, wflag,xscf, xscl, WCS, nx,ny, nmatch, XSCshape, 
     1     nstar, NSUMmax, RAlist, Declist,  XSCprox, SPIT, MIPS)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	character*(*) xscf, xscl
	character*500 string,header
	character*25 s0
	integer WCS(4)
	integer nx(4),ny(4), nsrc
	integer*2 wflag(nfi,4)
	
	real*8 ra8,dec8,x8,y8, RAlist(NSUMmax), Declist(NSUMmax)
	real*4 ra,dec, ba,pa, XSCprox(nstar)
	real*4 XSCpos (999,5)
	real*4 XSCshape (nstar,5)
	real*4 Rout(5),BAout(5)
	real*4 BAwin(9), Rwin(9)

	logical SPIT, MIPS, dowrite


	 nsrc = 0

	dowrite = .false.
c	dowrite = .true.

	nmax = 999

	i2mass = 2
	itype=4  !  3 = irac; 4 = wise; 5  = w4
        if (SPIT) itype=3

c	write (6,*) 'inside loadXSC',nx(1),ny(1)

	open (unit=11,file=xscf)

	header = ''
	ihead = 0
	do 110 K=1,999
	  read (11,'(a)',end=19) string
          if (string(1:1).eq.'\') goto 110
          if (string(1:1).eq.'|') then 
            if (ihead.eq.0) then
                ihead = 1
                header = string
                goto 19
            endif
	  endif
 110	continue
 19	close (11)


c got the header; get locations of keywords
	L = numchar (header)
	do k=1,L
	  if (header(k:k).eq.'|') header(k:k)=' '
	enddo


	s0 = 'ra'
        call header_parse (header,s0,ira)
	s0 = 'dec'
        call header_parse (header,s0,idec)

	s0 = 'r_k20fe'
        call header_parse (header,s0,irk)
	s0 = 'sup_ba'
        call header_parse (header,s0,iba)
	s0 = 'sup_phi'
        call header_parse (header,s0,ipa)


	open (unit=11,file=xscf)

	nsrc = 0

	do 100 K=1,1000

	  read (11,'(a)',end=99) string

	  if (string(1:1).eq.'\') goto 100
	  if (string(1:1).eq.'|') goto 100

	  call sfields (string,ira,s0)
          read (s0,*) ra
	  call sfields (string,idec,s0)
          read (s0,*) dec


	  call sfields (string,irk,s0)
	  if (s0(1:4).eq.'null') goto 100  ! do not use NULL radii
          read (s0,*) rk20fe
	  call sfields (string,iba,s0)
	  if (s0(1:4).eq.'null') s0 = "1."
          read (s0,*) ba
	  call sfields (string,ipa,s0)
	  if (s0(1:4).eq.'null') s0 = "0."
          read (s0,*) pa


	  ra8 = ra * 1.d0
	  dec8 = dec * 1.d0

	  do ib=1,4
	   x0 = 0.
	   y0 = 0.

c	   if (wflag(1,ib) .eq. 1) then
	   if ( any ( wflag ( 1:nfi,ib ) ==1 ) ) then

	    offscl = -1
            call wcs2pix(WCS(ib), ra8, dec8, x8, y8, offscl)

	    goto 47
	   endif
	  enddo

 47       x0 = x8 * 1.0
          y0 = y8 * 1.0

   	  if ( (x0.gt.0.) .and. (x0.le.nx(ib)) ) then
	     if ( (y0.gt.0.) .and. (y0.le.ny(ib)) ) then
	       nsrc = nsrc + 1
	       XSCpos (nsrc,1) = ra
               XSCpos (nsrc,2) = dec
	       XSCpos (nsrc,3) = rk20fe
	       XSCpos (nsrc,4) = ba
	       XSCpos (nsrc,5) = pa

	    endif
	  endif


	if (nsrc.eq.nmax) goto 99  ! quit out, max soruces loaded

 100	continue
 99	close (11)

c	write (6,*) nsrc,'  XSC galaxies loaded'


c now  run through source list
c	write (6,*) nstar

	if (dowrite)
     1  write (66,'(a)') '|    ra   |   dec   |rprox  |r_k20fe |sup_ba  |sup_phi |   rW4  |  baW4  |'

	nmatch = 0

c	write (6,*) 'here ',nstar

	do 5000 J = 1, nstar 

	   ra0 = ralist (J) * 1.0
	   dec0 = declist (J) * 1.0
	   dd = dec0/ 57.2957795

c	write (6,'(2f10.5)') ra0,dec0

	   XSCprox(J) = -999.
	   XSCshape (J,1) = 0.
	   XSCshape (J,2) = 0.
	   XSCshape (J,3) = 0.
	   XSCshape (J,4) = 0.
	   XSCshape (J,5) = 0.

	   idex = 0
	   rmin = 9999.
	   do K = 1, nsrc
	     ra =  XSCpos (K,1) 
             dec = XSCpos (K,2) 
	     rk20fe = XSCpos (K,3)
	     ba = XSCpos (K,4)
	     pa = XSCpos (K,5)

	     rcheck = rk20fe * 1.1  ! add 10%

	     call dcon (pa,pa_math)


	     delpos = 9999.
	     call check (ra0,dec0,dd,ra,dec,dra,ddec,delpos)

	     batmp = 1.0
	     call w_ell (batmp,pa_math,-dra,ddec,rmajor)


	     if (rmajor.le.rcheck) then
		if (rmajor.lt.rmin) then
			idex = K
			rmin = rmajor
		endif
	     endif	

	   enddo

	   if (idex.gt.0.) then
		XSCprox(J) = rmin
		rk20fe = XSCpos (idex,3)
                ba = XSCpos (idex,4)
                pa = XSCpos (idex,5)
		nmatch = nmatch + 1
	   else
		goto 5000
	   endif

	   open (unit=21,file=xscl)
	   Dmin = 1.

           do 300 jj=1,99999

            read (21,'(a)',end=22) string
            if (string(1:1).eq.'|') goto 300
            if (string(1:1).eq.'\') goto 300

            read (string,*) Rout(1), BAout(1),
     1       (Rout(ib0),BAout(ib0),ib0=2,5)

             BAfrac = BAout(i2mass) / ba
             if (BAfrac.gt.1.) BAfrac = 1./BAfrac

             Rfrac = Rout(i2mass) / rk20fe
             if (Rfrac.gt.1.) Rfrac = 1./Rfrac

             RQ = sqrt ( (BAfrac**2) + (Rfrac ** 2) )

             DD = 1.41 - RQ
	     if (DD.lt.Dmin) then
                Dmin = DD
                do K=1,5
                   Rwin (K) = Rout(K)
                   BAwin(K) = BAout(K)
                enddo
             endif

 300       continue
 22        close (21)

c	write (6,*) rk20fe, ba
c	write (6,*) Rwin(i2mass), BAwin(i2mass) , Rwin (itype), BAwin(itype)
	

	if (MIPS) then

		R_ratio = Rwin (itype+1) / Rwin(i2mass)   ! these should be MIPS 
        	BA_ratio = BAwin(itype+1) / BAwin(i2mass) ! these should be MIPS

		XSCshape (J,4) = rk20fe * R_ratio
        	XSCshape (J,5) = ba * BA_ratio
        	XSCshape (J,5) = min (XSCshape (J,5), 1.0)

	else if (itype.eq.4) then  !  W4

		R_ratio = Rwin (itype+1) / Rwin(i2mass)   ! these should be W4 
                BA_ratio = BAwin(itype+1) / BAwin(i2mass) ! these should be W4 

                XSCshape (J,4) = rk20fe * R_ratio
                XSCshape (J,5) = ba * BA_ratio
                XSCshape (J,5) = min (XSCshape (J,5), 1.0)

	else
		XSCshape (J,4) = 0.
		XSCshape (J,5) = 0.
	endif

	R_ratio = Rwin (itype) / Rwin(i2mass)
        BA_ratio = BAwin(itype) / BAwin(i2mass)

        XSCshape (J,1) = rk20fe * R_ratio
        XSCshape (J,2) = ba * BA_ratio
        XSCshape (J,2) = min (XSCshape (J,2), 1.0)

	XSCshape (J,3) = pa

c	write (6,*) XSCshape (J,1), XSCshape (J,2)
c	write (6,*) XSCshape (J,4), XSCshape (J,5)
c	write (6,*) ' '

c	if ((idex.gt.0).and.(dowrite)) then
c	   write (66,'(2f10.5,1x,f7.2,5f9.2)') ra0,dec0,XSCprox(J),
c     1 (XSCshape (J,kk),kk=1,5)
c	endif

 5000	continue
	if (dowrite) close (66)
	   
	return
	end

