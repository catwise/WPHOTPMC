	subroutine check (ra0,dec0,dd,ra,dec,delra,deldec,delpos)

	implicit integer (i-n)
        implicit real*4 (a-h)

        real*4 ra0,dec0,dd,ra,dec

        delra =  abs(ra - ra0 )
        if (delra .ge. 180.0) delra = delra - sign(360.0,delra)
        delra = 3600. * delra  * cos(dd)
        deldec = 3600. * (dec - dec0 )
        delpos = sqrt( delra**2 + deldec**2)  ! arcsec

        return
        end

