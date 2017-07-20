      subroutine ChkVar(flux, fluxsig, fluxchi, fluxMJD, N,
     +                  Band, ChiFac, mLogQ, dKS, nDF)
c-----------------------------------------------------------------------
c
c  Compute the variability indicator mLogQ for the flux array
c
c  Version 1.0  B10419: initial version
c          1.1  B10608: added photometric chi-square filtering
c          1.2  B10609: blocked null measurements from fluxchi median
c          1.3  B10614: added data dump option via ChiFac < 0
c          1.4  B10616: switched table-file header I/O to Tom's method
c          1.5  B10624: added fluxMJD for correlating fluxes across band
c          1.6  B10713: use rchi2 filtering in cross-band correlation
c          1.7  B10718: added test on N>1 for correlation
c
c  Input
c         flux:    Array of flux values (real*4 array)
c         fluxsig: 1-sigma uncertainties of flux (real*4 array)
c         fluxchi: template-fit chi-squares (real*4 array)
c         fluxMJD: MJD array for fluxes (real*8 array)
c         N:       Number of fluxes and uncertainties (integer*4)
c         Band:    Band number (integer*4)
c         ChiFac : Chi-square threshold = ChiFac*Median(fluxchi)
c
c  Note: for any 0 < J <= N, flux(J), fluxsig(J), fluxchi(J), and
c        fluxMJD(J) all refer to the same frame; ChiFac = 0
c        turns off rchi2 filter & skips sorting for median;
c        ChiFac < 0 turns on data dump, e.g., ChiFac = -K-cf,
c        where K >= 100, means dump first K sources, use ChiFac = |cf|
c
c  Output
c         mLogQ:   -mLogQfac*log(Q), Q = 1-P(chi-square), 0-99.99,
c                  (real*4); mLogQfac scales to a different log base
c         DeltaMag: Full-width magnitude variation - deprecated - CJG
c         nDF:      No. of degrees of freedom in chi-square (integer*4)
c         dKS:       Stetson K index (ApJ 735, 68, 2011)  - CJG
c
c-----------------------------------------------------------------------
c
                     Integer*4  MaxFlux, MaxBuf
                     Parameter (MaxFlux = 1000, MaxBuf = 4*MaxFlux)
c
      character*9 NumStr, CFstr
      character*1 BndStr
      integer*4   N, Band, nBuf(4), k, m, nDF, nSrc, nDump,
     +            Iunit
      real*4      flux(N), fluxsig(N), FBuf(MaxFlux,4), Fmin, Fmax,
     +            AvgFlux, F(N), Wgt(N), AvgF(4), fluxchi(N), r4tmp,
     +            ChiMax, DeltaMag, ChiFac, chiflux(N), CF, mLogQ
      real*4      sdf,sf(N),factn,mag(N),emag(N),avgmag,dks      ! CJG B30319
      real*8      summ,sumn,sumd,sumdk                           ! CJG B30319
      real*8      ChiSq, SumF, sumWgt, Qchisq, Rtmp, fluxMJD(N),
     +            MJDBuf(MaxFlux,4)
      logical*4   FluxOK, OKflux(MaxFlux,4), DoFchisq, DoDump
c
      data      nBuf/4*0/, FBuf/MaxBuf*0.0/, AvgF/4*0.0/, nSrc/0/,
     +          DoDump/.false./, nDump/0/
c
      common/ckvarcom/MJDBuf,nBuf,FBuf,AvgF,OKflux
      include 'jwfcom.f'
c
c-----------------------------------------------------------------------
c
      nBuf(Band) =  0
      SumF       =  0.0d0
      SumWgt     =  0.0d0
      SumFu      =  0.0d0
      Summ       =  0.0d0
      Sumn       =  0.0d0
      nDF        =  0
      Fmax       = -9.9e9
      Fmin       =  9.9e9
      avgmag     =  0.0
      if (N .lt. 2) go to 110
c
      CF       = ChiFac                ! don't modify input parameter
      DoFchisq = (ChiFac .ne. 0.0)
      DoDump   = (ChiFac .lt. 0.0)
      if (DoDump) then
        CF    = abs(ChiFac)
        m     = CF/100.0
        nDump = 100*m
        CF    = CF - nDump
      end if
c
      m = 0
      do 50 k = 1, N                   ! load multi-band arrays
        if (k .gt. MaxFlux) go to 60
        nBuf(Band)              = nBuf(Band) + 1
        FBuf(nBuf(Band),Band)   = flux(k)
        MJDBuf(nBuf(Band),Band) = fluxMJD(k)
        OKflux(nBuf(Band),Band) = FluxOK(flux(k), fluxsig(k))
        if (DoFchisq .and. OKflux(nBuf(Band),Band)) then
          m = m + 1
          chiflux(m) = fluxchi(k)
        end if
50    continue
60    if (DoFchisq .and. (m .lt. 2)) go to 110
c
      if (DoFchisq) then               ! get chi-square threshold
        call TJsort(m,chiflux)
        if (mod(m, 2) .eq. 1) then
          ChiMax = CF*chiflux(m/2)
        else
          ChiMax = CF*(chiflux(m/2)+chiflux(m/2+1))/2.0
        end if
      else
        ChiMax = 9.9e25
      end if
c
      do 100 k = 1, N
        if (OKflux(k,Band) .and.
     +     (fluxchi(k) .lt. ChiMax)) then
          nDF      = nDF + 1
          Wgt(nDF) = 1.0d0/fluxsig(k)**2
          F(nDF)   = flux(k)
          sf(ndf)  = fluxsig(k)
          SumF     = SumF   + Wgt(nDF)*flux(k)
          SumWgt   = SumWgt + Wgt(nDF)
          mag(ndf) = -2.5*log10(flux(k))                      ! CJG  B30319
          emag(ndf) = (1.0857 * fluxsig(k)) / flux(k)         ! CJG  B30319
          Summ    = Summ   + mag(ndf)                         ! CJG  B30319
          Sumn    = Sumn   + 1.d0                             ! CJG  B30319
          if (flux(k) .gt. Fmax) Fmax = flux(k)
          if (flux(k) .lt. Fmin) Fmin = flux(k)
        else
          OKflux(k,Band) = .false.
        end if
100   continue
c
110   if (nDF .lt. 2) then             ! we'll need to subtract 1 below
        nDF      = 0
        mLogQ    = -1.0
        DeltaMag = 99.0
        go to 300
      end if
c
      AvgFlux    = SumF/sumWgt
      AvgF(Band) = AvgFlux
      DeltaMag   = 2.5*log10(Fmax/Fmin)
      if (DeltaMag .gt. 98.0) DeltaMag = 98.0
c
      ChiSq = 0.0d0
      do 200 k = 1, nDF
        ChiSq = ChiSq + Wgt(k)*(F(k) - AvgFlux)**2
200   continue
c
      nDF = nDF - 1                    ! lost one DF in average
      r4tmp = Qchisq(Chisq,nDF)
      if (r4tmp .gt. 0.0) then
        mLogQ = -mLogQfac*log(r4tmp)
      else
        mLogQ = 99.99
      end if
c     if (mLogQ .lt. 0) mLogQ = 0  ! leave negative, will be null in mdex table
      if (mLogQ .gt. 99.99) mLogQ = 99.99
c
c==========================start of code added by CJG B30319=================
c
c      write(6,*)'ndf: ',ndf
      if(ndf.gt.2) then
         sdf=ndf
         factn=sqrt(sdf/(sdf-1))
         avgmag=summ/sumn
         sumd = 0.d0
         sumd2 = 0.d0
         do 250 k = 1, ndf
            deltap = factn * (mag(k) - avgmag) / emag(k)
            sumd = sumd + abs(deltap)
            sumd2 = sumd2 + deltap**2
 250     continue
         
         dks = (1.0 / sqrt(sdf) ) * sumd / sqrt(sumd2)
         
      else
         
         dks = -99.0

      endif
      
         
c      write(6,*)'dks: ',avgmag,sdf,deltap,sumd,sumd2,dks
c
c==========================end of code added by CJG B30319=================
c
c
300   if (DoDump) go to 400
      return
c
400   if (Band .eq. 1) nSrc = nSrc + 1
      if (nSrc .ge. nDump) return
      write (NumStr,'(I9)') nSrc
405   if (NumStr(1:1) .eq. ' ') then
        NumStr = NumStr(2:9)
        if (lnblnk(NumStr) .ne. 0) go to 405
        NumStr = 'xxxxx'               ! this should never happen
      end if
      write (CFstr,'(F9.1)') CF
410   if (CFstr(1:1) .eq. ' ') then
        CFstr = CFstr(2:9)
        if (lnblnk(CFstr) .ne. 0) go to 410
        CFstr = NumStr                 ! this should never happen
      end if
      BndStr = achar(48+Band)
      Iunit = 73
      open (unit=Iunit,
     +  file='/wise/data2/jwf/chkvar/w'
     +                      //BndStr//'/chkvar_w'//BndStr//'_cf'
     +                      //CFstr(1:lnblnk(CFstr))//'_'
     +                      //NumStr(1:lnblnk(NumStr))//'.tbl')
      write(Iunit,'(A,I5)')   '\Ndet =       ', N
      write(Iunit,'(A,F5.2)') '\mLogQ =      ', mLogQ
      write(Iunit,'(A,I5)')   '\Ndf =        ', nDF
      if (nDF .gt. 0) then
        Rtmp = Qchisq(Chisq,nDF)
      else
        Rtmp  = 0.5
        ChiSq = 0.0
        Fmin  = 0.0
        Fmax  = 0.0
      end if
      write(Iunit,'(A,1pE10.3)') '\ChiSq =      ', ChiSq
      write(Iunit,'(A,1pE10.3)') '\QChiSq =     ', Rtmp
      write(Iunit,'(A,F10.3)')   '\Log10Q =     ', dlog10(Rtmp)
      write(Iunit,'(A,F10.3)')   '\ChiFac =     ', CF
      write(Iunit,'(A,F10.3)')   '\ChiMax =     ', ChiMax
      write(Iunit,'(A,F10.3)')   '\ChiMed =     ', ChiMax/CF
      write(Iunit,'(A,F17.3)')   '\Fmin =', Fmin
      write(Iunit,'(A,F17.3)')   '\Fmax =', Fmax
      write(Iunit,'(A,F17.3)')   '\Favg =', AvgFlux
      write(Iunit,'(A,F10.3)')   '\DeltaMag =   ', DeltaMag
      write(Iunit,
     + '(''|  N  |U|   Flux   | SigFlux | ChiFlux|      MJD     |'')')
      write(Iunit,
     + '(''| int |c|   real   |  real   |  real  |     real     |'')')
      write(Iunit,
     + '(''|     | |    DN    |   DN    |        |      MJD     |'')')
      write(Iunit,
     + '(''| null|n|   null   |  null   |  null  |     null     |'')')
      do 500 k = 1, N
        if (OKflux(k,Band)) then
          write(Iunit,'(I6,'' Y'',F11.2,F10.2,F9.2,F15.8)')
     +          k, flux(k), fluxsig(k), fluxchi(k), fluxMJD(k)
        else
          write(Iunit,'(I6,'' N'',F11.2,F10.2,F9.2,F15.8)')
     +          k, flux(k), fluxsig(k), fluxchi(k), fluxMJD(k)
        end if
500   continue
      close(Iunit)
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine GetCorrs(Rho12,Rho23,Rho34, Q12,Q23,Q34, N12,N23,N34)
c-----------------------------------------------------------------------
c
c  Compute the band-to-band sample correlation coefficients and their
c  probabilities for adjacent bands, using flux data buffered by ChkVar
c
c  Output
c         Rho12:  W1W2 correlation coefficient (percentage, integer*4)
c         Rho23:  W2W3 correlation coefficient (percentage, integer*4)
c         Rho34:  W3W4 correlation coefficient (percentage, integer*4)
c         Q12:    -log10(1-P(Rho12)) given no real correlation
c                 (0-9, integer*4)
c         Q23:    -log10(1-P(Rho23)) given no real correlation
c                 (0-9, integer*4)
c         Q34:    -log10(1-P(Rho34)) given no real correlation
c                 (0-9, integer*4)
c         N12:    No. of flux sets in W1 and W2 which were usable
c                 for the Rho12 correlation calculation (integer*4)
c         N23:    No. of flux sets in W2 and W3 which were usable
c                 for the Rho23 correlation calculation (integer*4)
c         N34:    No. of flux sets in W3 and W4 which were usable
c                 for the Rho34 correlation calculation (integer*4)
c
c                 Note: -99 <= Rho12 <= 99, 0 <= Q12 <= 9
c                       -99 <= Rho23 <= 99, 0 <= Q23 <= 9
c                       -99 <= Rho34 <= 99, 0 <= Q34 <= 9
c
c-----------------------------------------------------------------------
c
                     Integer*4  MaxFlux, MaxBuf
                     Parameter (MaxFlux = 1000, MaxBuf = 4*MaxFlux)
c
      integer*4 nBuf(4), Rho12, Rho23, Rho34, Q12, Q23, Q34, k1, k2, k3,
     +          k4, N12, N23, N34
      real*4    FBuf(MaxFlux,4),AvgF(4),Rho
      real*8    Sum11, Sum12, Sum22, Sum23, Sum33, Sum34, Sum44, Qcorr,
     +          MJDBuf(MaxFlux,4), dMJD
      logical*4 OKflux(MaxFlux,4)
      data      dMJD/5.0d-5/           ! 4.32 sec
c
      common/ckvarcom/MJDBuf,nBuf,FBuf,AvgF,OKflux
c
      Rho12 = 0
      Rho23 = 0
      Rho34 = 0
      Q12   = 0
      Q23   = 0
      Q34   = 0
      N12   = 0
      N23   = 0
      N34   = 0
c
      if ((nBuf(1) .gt. 1) .and. (nBuf(2) .gt. 1)) then
        Sum11 = 0.0d0
        Sum12 = 0.0d0
        Sum22 = 0.0d0
        do 120 k2 = 1, nBuf(2)
          if (.not.OKflux(k2,2)) go to 120
          do 110 k1 = 1, nBuf(1)
            if (dabs(MJDBuf(k1,1)-MJDBuf(k2,2)) .gt. dMJD)  go to 110
            if (.not.OKflux(k1,1)) go to 120
            N12   = N12   + 1
            Sum11 = Sum11 + (FBuf(k1,1) - AvgF(1))**2
            Sum12 = Sum12 + (FBuf(k1,1) - AvgF(1))
     +                     *(FBuf(k2,2) - AvgF(2))
            Sum22 = Sum22 + (FBuf(k2,2) - AvgF(2))**2
            go to 120
110       continue
120     continue
        if ((Sum11 .gt. 0.0d0) .and. (Sum22 .gt. 0.0d0)
     +                         .and. (N12 .gt. 1)) then
          Rho   = Sum12/dsqrt(Sum11*Sum22)
          Rho12 = NInt(100.0*Rho)
          if (Rho12 .lt. -99) Rho12 = -99
          if (Rho12 .gt.  99) Rho12 =  99
          Q12 = NInt(-dlog10(Qcorr(Rho, N12)))
          if (Q12 .lt. 0) Q12 = 0
          if (Q12 .gt. 9) Q12 = 9
        else
          Rho12 = 0
          Q12   = 0
          N12   = 0                    ! tag as null
        end if
      end if
c
      if ((nBuf(2) .gt. 1) .and. (nBuf(3) .gt. 1)) then
        Sum22 = 0.0d0
        Sum23 = 0.0d0
        Sum33 = 0.0d0
        do 220 k3 = 1, nBuf(3)
          if (.not.OKflux(k3,3)) go to 220
          do 210 k2 = 1, nBuf(2)
            if (dabs(MJDBuf(k2,2)-MJDBuf(k3,3)) .gt. dMJD)  go to 210
            if (.not.OKflux(k2,2)) go to 220
            N23   = N23   + 1
            Sum22 = Sum22 + (FBuf(k2,2) - AvgF(2))**2
            Sum23 = Sum23 + (FBuf(k2,2) - AvgF(2))
     +                     *(FBuf(k3,3) - AvgF(3))
            Sum33 = Sum33 + (FBuf(k3,3) - AvgF(3))**2
            go to 220
210       continue
220     continue
        if ((Sum22 .gt. 0.0d0) .and. (Sum33 .gt. 0.0d0)
     +                         .and. (N23 .gt. 1)) then
          Rho   = Sum23/dsqrt(Sum22*Sum33)
          Rho23 = NInt(100.0*Rho)
          if (Rho23 .lt. -99) Rho23 = -99
          if (Rho23 .gt.  99) Rho23 =  99
          Q23 = NInt(-dlog10(Qcorr(Rho, N23)))
          if (Q23 .lt. 0) Q23 = 0
          if (Q23 .gt. 9) Q23 = 9
        else
          Rho23 = 0
          Q23   = 0
          N23   = 0                    ! tag as null
        end if
      end if
c
      if ((nBuf(3) .gt. 1) .and. (nBuf(4) .gt. 1)) then
        Sum33 = 0.0d0
        Sum34 = 0.0d0
        Sum44 = 0.0d0
        do 320 k4 = 1, nBuf(4)
          if (.not.OKflux(k4,4)) go to 320
          do 310 k3 = 1, nBuf(3)
            if (dabs(MJDBuf(k3,3)-MJDBuf(k4,4)) .gt. dMJD)  go to 310
            if (.not.OKflux(k3,3)) go to 320
            N34   = N34   + 1
            Sum33 = Sum33 + (FBuf(k3,3) - AvgF(3))**2
            Sum34 = Sum34 + (FBuf(k3,3) - AvgF(3))
     +                     *(FBuf(k4,4) - AvgF(4))
            Sum44 = Sum44 + (FBuf(k4,4) - AvgF(4))**2
            go to 320
310       continue
320     continue
        if ((Sum33 .gt. 0.0d0) .and. (Sum44 .gt. 0.0d0)
     +                         .and. (N34 .gt. 1)) then
          Rho   = Sum34/dsqrt(Sum33*Sum44)
          Rho34 = NInt(100.0*Rho)
          if (Rho34 .lt. -99) Rho34 = -99
          if (Rho34 .gt.  99) Rho34 =  99
          Q34 = NInt(-dlog10(Qcorr(Rho, N34)))
          if (Q34 .lt. 0) Q34 = 0
          if (Q34 .gt. 9) Q34 = 9
        else
          Rho34 = 0
          Q34   = 0
          N34   = 0                    ! tag as null
        end if
      end if
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
1000  nBuf(1) = 0
      nBuf(2) = 0
      nBuf(3) = 0
      nBuf(4) = 0
c
      return
      end
c
c=======================================================================
c
      function FluxOK(Flux, SigF)
c
      logical*4 FluxOK
      real*4    Flux, SigF, FluxMin, SigFMax
c     data      FluxMin/-998.0/, SigFMax/1.0e7/
      data      FluxMin/1.0/, SigFMax/1.0e7/
c
c Now for the input fluxes and their associated error, you will
c need to test that the flux error is positive (i.e., not equal to zero
c or a negative value) and that it does not exceed some huge value.
c This is my test condition:
c         error = flux_unc(j)
c         if (error.le.0.) REJECT
c         if (error.gt.1.e7) REJECT
c
      if (Flux .lt. FluxMin) go to 30
      if (SigF .le. 0.0)     go to 30
      if (SigF .gt. SigFMax) go to 30
c
      FluxOK = .true.
      return
c
30    FluxOK = .false.
      return
      end
c
c=======================================================================
c
      function Qcorr(Rho, N)
c
      integer*4 N
      real*4    Rho
      real*8    Qcorr, Z, Const, Q
      data      Const/0.7071067812d0/
c
      if (abs(Rho) .ge. 1.0) then
        Qcorr = 1.0d-25
        return
      end if
c
      if (N .lt. 4) then
        Qcorr = 0.5d0                  ! 50-50, no significance
        return
      end if
c
      Z = 0.5d0*dlog((1.0d0 + dble(Rho))/(1.0d0 - dble(Rho)))
      Q = derfc(Const*dabs(Z)*dsqrt(dfloat(N)-3.0d0))
      if (Q .lt. 1.0d-25) Q = 1.0d-25
      Qcorr = Q
c
      return
      end
c
c
c=======================================================================
c
      function Qchisq(ChiSq,Ndf)
c
      integer*4 Ndf
      real*8    Qchisq, GammQ, ChiSq, DegFrdm
c
      DegFrdm = Ndf
      Qchisq = Gammq(0.5d0*DegFrdm,0.5d0*ChiSq)
c
      return
      end
c
c=======================================================================
c
      FUNCTION gammq(a,x)
      real*8   a,gammq,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if((x .lt. 0.0d0) .or. (a .le. 0.0d0)) then
        print *,'ERROR: bad arguments in gammq'
        print *,'       x = ',x
        print *,'       a = ',a
        call exit(64)
      end if
      if(x .lt. a+1.0d0)then
        call gser(gamser,a,x,gln)
        gammq=1.0d0 - gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.0d-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x .le. 0.0d0) then
c       if(x .lt. 0.0d0) pause 'x < 0 in gser'
        gamser=0.0d0
        return
      endif
      ap=a
      sum=1.0d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.0d0
        del=del*x/ap
        sum=sum+del
        if(abs(del) .lt. abs(sum)*EPS) goto 1
11    continue
      print *,'WARNING: a too large, ITMAX too small in gser'
      print *,'         a = ',a
      print *,'         x = ',x
      print *,'     ITMAX = ',ITMAX
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.0d-7,FPMIN=1.0d-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.0d0-a
      c=1.0d0/FPMIN
      d=1.0d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2. 0d0
        d=an*d+b
        if(abs(d) .lt. FPMIN) d=FPMIN
        c=b+an/c
        if(abs(c) .lt. FPMIN) c=FPMIN
        d=1.0d0/d
        del=d*c
        h=h*del
        if(abs(del-1.0d0) .lt. EPS) goto 1
11    continue
      print *,'WARNING: a too large, ITMAX too small in gcf'
      print *,'         a = ',a
      print *,'         x = ',x
      print *,'     ITMAX = ',ITMAX
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
c-----------------------------------------------------------------------
      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.0d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C=======================================================================
c                                  Sort for real*4 array
c                                  from Numerical Recipes via T. Jarrett
      SUBROUTINE TJSORT(N,RA)
c
      Integer*4 N,L,IR,J,I
      Real*4 RA(N),RRA
c
      if (n .lt. 2) return
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
c
c=======================================================================
c
      subroutine SignOn(pgmnam)
c
c *** signon- routine which provides sign-on and sign-off messages
c             (orig by John Fowler- mod by Howard McCallon-041214-SIRTF)
c
c     inputs:  pgmnam = program name                                 [call arg]
c
c     outputs: message to stdout
c
      character*(*) pgmnam
      character vsn*11,cdate*8,ctime*8,Fmt*11,FLen*4
      integer*4 onoff,jdate(3),jtime(3),lnblnk
      real*4    dummyt,second(2),etime
c
      common /vdt/ cdate,ctime,vsn
c##
      onoff = 1
c
c         i. obtain date
c
100   cdate = '00-00-00'
c     call idate(jdate(1),jdate(2),jdate(3))    ! Linux call
      call idate(jdate)                         ! JWF B60713 - gfortran call
c
      jdate(3) = mod(jdate(3), 100)
      write(cdate(1:2), '(i2)') jdate(2)
      write(cdate(4:5), '(i2)') jdate(1)
      write(cdate(7:8), '(i2)') jdate(3)
c
      if(cdate(4:4) .eq. ' ') cdate(4:4) = '0'
      if(cdate(7:7) .eq. ' ') cdate(7:7) = '0'
c
c         ii. obtain time
c
      ctime = '00:00:00'
      call itime(jtime)
      write(ctime(1:2), '(i2)') jtime(1)
      write(ctime(4:5), '(i2)') jtime(2)
      write(ctime(7:8), '(i2)') jtime(3)
c
      if(ctime(4:4) .eq. ' ') ctime(4:4) = '0'
      if(ctime(7:7) .eq. ' ') ctime(7:7) = '0'
c
c         iii. set up format for pgmnam
c
      write(Flen,'(I4)') lnblnk(pgmnam)
      Fmt = '(A'//Flen//'$)'
c
c         iv. write out results
c
      write(*,Fmt) pgmnam
      if(onoff .eq. 1) then                      ! sign on
        write(*,301) vsn,cdate,ctime
      else                                       ! sign off
        dummyt = etime(second)
        write(*,302) vsn,cdate,ctime,second
      endif
  301 format(' version: ',a11,' - execution begun on ',a8,' at ',a8)
  302 format(' version: ',a11,' - execution ended on ',a8,' at ',a8
     *    /1x,f9.2,' cpu seconds used;',f8.2,' system seconds used.')
c
      return
c
      entry SignOff(pgmnam)
      OnOff = 2
      go to 100
c
      end
