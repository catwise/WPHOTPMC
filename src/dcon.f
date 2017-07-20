	subroutine dcon (pa,angle)
c convert from standard p.a. to mathematical angle (north of east)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

        phimax = pa
        if ((pa.ge.0.).and.(pa.le.90.)) then
                phimax = 90. + pa
        else if (pa.gt.90.) then
                phimax = pa - 90.
        else if (pa.lt.0.) then
                phimax = pa + 90.
        endif
        angle = phimax

        return
	end

	subroutine acon (angle,pa)
c convert to standard p.a. (E of N)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)


        zz = angle + 90.0
        if (zz.gt.90.) zz = zz - 180.
        pa = zz

        return
        end



