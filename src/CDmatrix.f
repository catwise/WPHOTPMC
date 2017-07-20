	subroutine CDmatrix (Hfits, cdelt1, cdelt2)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*8 cd1a,cd1b,cd2a,cd2b,tDB,rat,angle

	character*(*) Hfits
	character*50 s0,key


c	write (6,*) 'CDMARIX'

	val = 0.
	key = "CD1_1"
	call keyhead (Hfits,key,s0)
	L = numchar(s0)
	if (L.gt.0) then
                 read (s0,*) val
	endif
	cd1_1 = val


	val = 0.
	key = "CD1_2"
        call keyhead (Hfits,key,s0)
        L = numchar(s0)
        if (L.gt.0) then
                 read (s0,*) val
        endif
        cd1_2 = val

	val = 0.
	key = "CD2_1"
        call keyhead (Hfits,key,s0)
        L = numchar(s0)
        if (L.gt.0) then
                 read (s0,*) val
        endif
        cd2_1 = val

	val = 0.
	key = "CD2_2"
        call keyhead (Hfits,key,s0)
        L = numchar(s0)
        if (L.gt.0) then
                 read (s0,*) val
        endif
        cd2_2 = val

c	write (6,*) cd1_1,cd1_2,cd2_1,cd2_2


         if (cd2_2.ne.0.) then
                  rat = CD1_2 / CD2_2

                  angle = -datan (rat) * 57.2957795d0

		  tdb = angle/57.2957795d0

	  	  cd2a = CD2_2 / dcos(tdb)
	  	  cd2b = -CD1_2 / dsin(tdb)

		if ((abs(CD2_2).gt.abs(CD1_2))) then
                      cdelt2 = cd2a*1.
                else
                      cdelt2 = cd2b*1.
                endif

	  	  cd1a = -CD1_1 / dcos(tdb)
	  	  cd1b = -CD2_1 / dsin(tdb)

		if ((abs(CD1_1).gt.abs(CD2_1))) then
                     cdelt1 = cd1a*1.
                else
                     cdelt1 = cd1b*1.
                endif

                  crot = angle*1.

	else if (cd1_1.ne.0.) then

		  	rat = CD2_1 / CD1_1
			angle = datan (rat) * 57.2957795d0
			tdb = angle/57.2957795d0

			cd2a = CD2_2 / dcos(tdb)
			cd2b = -CD1_2 / dsin(tdb)

			if ((abs(CD2_2).gt.abs(CD1_2))) then
			  cdelt2 = cd2a*1.
		 	else
			  cdelt2 = cd2b*1.
			endif

			cd1a = -CD1_1 / dcos(tdb)
			cd1b = -CD2_1 / dsin(tdb)
			if ((abs(CD1_1).gt.abs(CD2_1))) then
			 cdelt1 = cd1a*1.	
		        else
			 cdelt1 = cd1b*1.
			endif

			crot = angle*1.

        endif


      return
      end

