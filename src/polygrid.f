	subroutine polygrid(igrid,rotA2,y)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

	real*4 A(3),x,y

ccc    rotA2 is the rotational angle

	rot = 0.0  

	if (abs(rotA2).le.45.) then
            rot = -1. * abs(rotA2)

        else if ((abs(rotA2).gt.45.).and.(abs(rotA2).lt.90.)) then
                rot = -1. * abs(rotA2)

        else if ((abs(rotA2).eq.90.).or.(abs(rotA2).eq.180.).or.(abs(rotA2).eq.270.)) then
                rot = 0.
        else if ((abs(rotA2).gt.90.).and.(abs(rotA2).lt.180.)) then
                rot = 90. - abs(rotA2)
        else if ((abs(rotA2).gt.180.).and.(abs(rotA2).lt.270.)) then
                rot = 180. - abs(rotA2)
        else if ((abs(rotA2).gt.270.).and.(abs(rotA2).lt.360.)) then
                rot = 270. - abs(rotA2)
        else if (abs(rotA2).gt.360.) then
                rot = 360. - abs(rotA2)
        endif


c	write (6,*) rot   !  TEMP

	if (igrid.eq.3) then
	    A(1) = -0.6181749
	    A(2) = -6.090606
	    A(3) =  -0.06515151
	else if (igrid.eq.5) then
	     A(1) = -0.6181749
	     A(2) =  -3.686515 
	     A(3) =  -0.03984848
	else if (igrid.eq.7) then
	     A(1) = -0.6636426  
	     A(2) = -2.586819
	     A(3) = -0.02742425
	else if (igrid.eq.9) then
	     A(1) = -0.499993
             A(2) = -2.033484
             A(3) = -0.02196969
	else
	     y = 0.
	     return
	endif


! y is the buffer size "nbuf"

	y = 0.
        do j=1,3
          y = y + (A(j) * rot**(j-1))
        enddo

	y = min (y,150.)  ! absolute max buffer size

	return
	end
