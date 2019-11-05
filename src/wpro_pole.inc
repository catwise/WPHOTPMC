! wpro.inc:   Set some dimensioning parameters for wpro.f
c
c B91103:
c Here are the numbers of “frames” (actually epochs) for the ecliptic
c poles in NEO3/NEO5:
c
c NEP (2709p666) 149/221
c SEP (0890m667) 149/221
c
cPeter

c	  integer, parameter :: maxframes = 224  ! maximum coverage x 4   ; 08Jan2018  THJ
c     integer, parameter :: maxframes = 320  ! maximum coverage x 4   !  B80622 JWF
      integer, parameter :: maxframes = 480  ! maximum coverage x 4   !  B91103 JWF

c	  integer, parameter :: maxframesin = 112  ! maximum # input frames   ; 08Jan2018  THJ
c     integer, parameter :: maxframesin = 160  ! maximum # input frames   ! B80622 JWF
c     integer, parameter :: maxframesin = 240  ! maximum # input frames   ! B91103 JWF
      integer, parameter :: maxframesin = 320  ! maximum # input frames   ! B91104 JWF
                           ! in any one band

c     integer, parameter :: maxpix = 50000 	   ! maximum number of pixels
      integer, parameter :: maxpix = 72000 	   ! maximum number of pixels !  B80622 JWF
      integer, parameter :: maxpix = 10800 	   ! maximum number of pixels !  B91103 JWF
						   ! per band in the solution 
						   ! for any one source
