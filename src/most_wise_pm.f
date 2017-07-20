c       subroutine most_wise(p,N,psfs,npsize,npnmax,ftol,iter,fret)     ! JWF B21207
        subroutine most_wise_pm(p,N,psfs,npsize,npnmax,ftol,iter,fret,
     +                          maxsteps,nsteps_tot,ok2,dbg) ! JWF B30117
!
! Method Of STeepest descent for WISE.  Minimize function FUNC with
! respect to vector P of length N.  On input, fret should be set to the 
! initial step size. On output, it is the minimum function value.
!
	implicit real*4(a-h,o-z)
	implicit integer(i-n)
c	integer, parameter :: nmax=200, itmax=100, maxsteps=100  ! JWF B30117
	integer, parameter :: nmax=200, itmax=100                ! JWF B30117
        integer*4 maxsteps, nsteps, iter, nsteps_tot, nDoF_pm    ! JWF B30419
	real, parameter :: eps=1.e-10, shrink=0.5
        real(4) p(N), g(nmax), pnew(nmax)
	logical calculate_gradient, finished
        logical*4 ok2, dbg  ! JWF B21207
		real(4) psfs(npsize,npsize,npnmax,4)


        if(fret <= 0.) then
            print *,'MOST_WISE - initial step size not set'

            return
        end if

        step = fret
        if (dbg) print *,'most_wise_pm(26): step, 3600*step:', step,3600*step
        iter = 1
        alpha = 0
	finished = .false.
        nsteps_tot = 0
	calculate_gradient = .true.

	do while (.not.finished)
	    if (calculate_gradient) then
c               print *,'most_wise(30): calling dfunc_pm'    ! JWF dbg
c		call dfunc(p,psfs,npsize,npnmax,g,fp)        ! JWF B21207
		call dfunc_pm(p,psfs,npsize,npnmax,g,fp,ok2) ! JWF B21207
c               print *,'most_wise(32): back from dfunc_pm'  ! JWF dbg
		gradmag = sqrt(sum(g(1:N)**2))

                if (dbg) print *,'most_wise_pm(40): iter, gradmag:',iter,gradmag

		if (gradmag == 0.) then
		    fret = fp
		    return
		endif

		alpha = 1./gradmag
		calculate_gradient = .false.
	    endif

	    stepal = step*alpha
	    fpold = 1.e35
	    do i = 1,N
		pnew(i) = p(i)
	    enddo
	    fpnew = fp    !  fp = chi-square from the dfunc_pm call
	    nsteps = 0

	    do while (fpnew < fpold .and. nsteps < maxsteps)
		nsteps = nsteps + 1
                nsteps_tot = nsteps_tot + 1
		fpold = fpnew
                if (dbg) print *,'most_wise_pm(64): iter, nsteps:',iter,nsteps
		do i = 1,N
		    pnew(i) = pnew(i) - stepal*g(i)
                    if (dbg) then
                      print *,'most_wise_pm(66): i,3600*stepal*g(i):',
     +                       i,3600*stepal*g(i)
                      print *,'most_wise_pm(67): g(i)/gradmag:',g(i)/gradmag
                    end if
		enddo
                if (dbg) then
                      print *,'most_wise_pm(69): pnew:',pnew(1:N)
                      print *,'most_wise_pm(70): pmra, pmdec:',
     +                        3600*pnew(N-1),3600*pnew(N)
                end if
		fpnew = func_pm(pnew,psfs,npsize,npnmax,nDoF_pm,dbg)
                if (dbg) print *,'most_wise_pm(70): fpold,fpnew:',fpold,fpnew
	    enddo

	    do i = 1,N
		pnew(i) = pnew(i) + stepal*g(i)
	    enddo
	    fret = fpold
	    iter = iter + 1
	    delchi = 2.*abs(fret-fp)
	    chimin = ftol*(abs(fret)+abs(fp)+eps)
            if (dbg) then
              print *,'most_wise_pm(86): delchi, chimin:', delchi, chimin
              print *,'most_wise_pm(87): fret,fp:', fret,fp
            end if
            if(fret < fp) then
                if (dbg) then
                  print *,'most_wise_pm(91): p:', p
                end if
		p = pnew(1:N)
		if(delchi <= chimin .or. iter==itmax) then
                    if (dbg) print *,'most_wise_pm(96): delchi <= chimin .or. iter==itmax'
		    finished=.true.
		else
		    fp = fret
		    calculate_gradient = .true.
                    if (dbg) print *,'most_wise_pm(101): calculate_gradient = .true.'
		endif
            else
		if(delchi <= 0.01*chimin .or. iter == itmax) then
                    if (dbg) print *,'most_wise_pm(105): delchi <= 0.01*chimin .or. iter == itmax'
		    fret = fp
		    finished = .true.
		else
		    step = step*shrink
                    if (dbg) print *,'most_wise_pm(110): new step:',step
		endif
	    end if
	enddo

	return

        end
