	subroutine countfile (in,nfi)	
	character*(*) in
	character*512 cc
	integer j, nfi

	open (unit=11,file=in)
	nfi = 0
	do j=1,9999999
	 read (11,'(a)',end=99) cc

	  L = numchar(cc)
	  if (L.lt.10) then
              print *,'===WARNING: frame info file has a short line: ', j
              goto 100
          endif
	 if ( (cc(1:1).ne.'|').and.(cc(1:1).ne.'\') ) then  ! modified Sep 12, 2010
	   nfi=nfi+1
	 endif
 100	enddo
 99	close (11)


	return
	end
