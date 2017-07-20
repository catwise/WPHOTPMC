      SUBROUTINE SORTA(R,N)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION R(500,3)

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN

          L=L-1
          RRA=R(L,3)
	  RRB=R(L,2)
	  RRC=R(L,1)

        ELSE

          RRA=R(IR,3)
	  RRB=R(IR,2)
	  RRC=R(IR,1)

          R(IR,3)=R(1,3)
	  R(IR,2)=R(1,2)
	  R(IR,1)=R(1,1)

          IR=IR-1
          IF(IR.EQ.1)THEN

            R(1,3)=RRA
	    R(1,2)=RRB
	    R(1,1)=RRC
            RETURN

          ENDIF
        ENDIF

        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(R(J,3).LT.R(J+1,3))J=J+1
          ENDIF
          IF(RRA.LT.R(J,3))THEN
            R(I,3)=R(J,3)
	    R(I,2)=R(J,2)
	    R(I,1)=R(J,1)

            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF

        GO TO 20
        ENDIF

        R(I,3)=RRA
	R(I,2)=RRB
	R(I,1)=RRC

      GO TO 10
      END


      SUBROUTINE SORT3(N,RA,RB, RC)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION RA(N),RB(N), RC(N)

	if (n.lt.2) goto 99

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
	  RRC=RC(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
	  RRC=RC(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
	  RC(IR)=RC(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
	    RC(1)=RRC
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
            RB(I)=RB(J)
	    RC(I)=RC(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
	RC(I)=RRC
      GO TO 10

 99	i=1

      END


      SUBROUTINE SORT2(N,RA,RB)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION RA(N),RB(N)

	if (n.lt.2) goto 99

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
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
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10

 99	i=1

      END

      SUBROUTINE SORT(N,RA)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      DIMENSION RA(N)

	if (n.lt.2) return

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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Integer sort added by TPC
      SUBROUTINE SORTI4(N,RA)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)
      INTEGER*4 RA,RRA
      DIMENSION RA(N)

      if (n.lt.2) return

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


      SUBROUTINE SORTI2(N,RA)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)
      INTEGER*2 RA,RRA
      DIMENSION RA(N)
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

ccccccccccccccc

      SUBROUTINE SORTN(R,N,M)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      real*4 buf (50), R(15000,8)

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN

          L=L-1
c          RRA=R(L,3)
c	  RRB=R(L,2)
c	  RRC=R(L,1)

	  do ii=1,M
		buf(ii) = R(L,ii)
	  enddo
	  RRA=R(L,3)

        ELSE

	  do ii=1,M
                buf(ii) = R(IR,ii)
          enddo
          RRA=R(IR,3)

	  do ii=1,M
                R(IR,ii) = R(1,ii)
          enddo

          IR=IR-1
          IF(IR.EQ.1)THEN

	    do ii=1,M
                R(1,ii) = buf(ii)
	    enddo

            RETURN

          ENDIF
        ENDIF

        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(R(J,3).LT.R(J+1,3))J=J+1
          ENDIF
          IF(RRA.LT.R(J,3))THEN
	    do ii=1,M
                R(i,ii) = R(J,ii)
            enddo

            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF

        GO TO 20
        ENDIF

	do ii=1,M
                R(i,ii) = buf(ii)
        enddo


      GO TO 10
      END

cccccccccccccc

      SUBROUTINE SORTVB(R,N,M)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      real*4 buf (50), R(9,5000)

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN

          L=L-1
	  do ii=1,M
		buf(ii) = R(ii,L)
	  enddo
	  RRA=R(3,L)

        ELSE

	  do ii=1,M
                buf(ii) = R(ii,IR)
          enddo
          RRA=R(3,IR)

	  do ii=1,M
                R(ii,IR) = R(ii,1)
          enddo

          IR=IR-1
          IF(IR.EQ.1)THEN

	    do ii=1,M
                R(ii,1) = buf(ii)
	    enddo

            RETURN

          ENDIF
        ENDIF

        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(R(3,J).LT.R(3,J+1))J=J+1
          ENDIF
          IF(RRA.LT.R(3,J))THEN
	    do ii=1,M
                R(ii,i) = R(ii,J)
            enddo

            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF

        GO TO 20
        ENDIF

	do ii=1,M
                R(ii,i) = buf(ii)
        enddo


      GO TO 10
      END

cccccccccccccc

      SUBROUTINE SORT3CH(R,N,M,ID)

	implicit integer (i-n)
        implicit real*4 (a-h)
        implicit real*4 (o-z)

      real*4 buf (50), R(15,15000)

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN

          L=L-1
	  do ii=1,M
		buf(ii) = R(ii,L)
	  enddo
	  RRA=R(ID,L)

        ELSE

	  do ii=1,M
                buf(ii) = R(ii,IR)
          enddo
          RRA=R(ID,IR)

	  do ii=1,M
                R(ii,IR) = R(ii,1)
          enddo

          IR=IR-1
          IF(IR.EQ.1)THEN

	    do ii=1,M
                R(ii,1) = buf(ii)
	    enddo

            RETURN

          ENDIF
        ENDIF

        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(R(ID,J).LT.R(ID,J+1))J=J+1
          ENDIF
          IF(RRA.LT.R(ID,J))THEN
	    do ii=1,M
                R(ii,i) = R(ii,J)
            enddo

            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF

        GO TO 20
        ENDIF

	do ii=1,M
                R(ii,i) = buf(ii)
        enddo


      GO TO 10
      END

