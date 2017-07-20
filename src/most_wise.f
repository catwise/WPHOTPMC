c       subroutine most_wise(p,N,psfs,npsize,npnmax,ftol,iter,fret) ! JWF B30206
        subroutine most_wise(p,N,psfs,npsize,npnmax,ftol,iter,fret,nsteps_tot,dbg) ! JWF B30206
!
! Method Of STeepest descent for WISE.  Minimize function FUNC with 
! respect to vector P of length N.  On input, fret should be set to the 
! initial step size. On output, it is the minimum function value.
!
      implicit real*4(a-h,o-z)
      implicit integer(i-n)
      integer, parameter :: nmax=200, itmax=100, maxsteps=100
      real, parameter :: eps=1.e-10, shrink=0.5
        real(4) p(N), g(nmax), pnew(nmax)
      logical calculate_gradient, finished
        integer*4 iter, nsteps, nsteps_tot, nDoF
        logical*4 dbg
      real(4) PSFs(npsize,npsize,npnmax,4)

c      write (6,*) 'calling most_wise ',fret

        if(fret <= 0.) then
            print *,'MOST_WISE - initial step size not set'

            return
        end if

        step = fret
        iter = 1
        alpha = 0
      finished = .false.
      calculate_gradient = .true.
        nsteps_tot = 0

      do while (.not.finished)
          if (calculate_gradient) then
            call dfunc(p,psfs,npsize,npnmax,g,fp)
            gradmag = sqrt(sum(g(1:N)**2))

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
          fpnew = fp
          nsteps = 0

          do while (fpnew < fpold .and. nsteps < maxsteps)
            nsteps = nsteps + 1
                nsteps_tot = nsteps_tot + 1
            fpold = fpnew
            do i = 1,N
                pnew(i) = pnew(i) - stepal*g(i)
            enddo
            fpnew = func(pnew,psfs,npsize,npnmax,nDoF,dbg)
          enddo

          do i = 1,N
            pnew(i) = pnew(i) + stepal*g(i)
          enddo
          fret = fpold
          iter = iter + 1
          delchi = 2.*abs(fret-fp)
          chimin = ftol*(abs(fret)+abs(fp)+eps)

            if(fret < fp) then
            p = pnew(1:N)
            if(delchi <= chimin .or. iter==itmax) then
                finished=.true.       
            else
                fp = fret
                calculate_gradient = .true.
            endif
            else
            if(delchi <= 0.01*chimin .or. iter == itmax) then
                fret = fp
                finished = .true.
            else
                step = step*shrink
            endif
          end if
      enddo

      return

        end
