	subroutine getFITSkey(Hdr,key,idex)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

C  Print out all the header keywords in all extensions of a FITS file

	character*(*) Hdr, key


	call header_parse (Hdr,key,idex)



	return
	end

