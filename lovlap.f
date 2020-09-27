      SUBROUTINE LOVLAP (STRAD,A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      DIMENSION A(MXUSER),B(MXUSER)
      DIMENSION FACT(25)
      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX
      DIMENSION BINCOE(7,7)
      LOGICAL*1 JGO
      DATA BINCOE/7*1.D0,   0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     1 2*0.D0,1.D0,3.D0,6.D0,10.D0,15.D0,   3*0.D0,1.D0,4.D0,10.D0,
     2 20.D0,   4*0.D0,1.D0,5.D0,15.D0,   5*0.D0,1.D0,6.D0,   6*0.D0,
     3 1.D0/
      DATA JGO/.FALSE./
      SAVE FACT,JGO
C
C    SUBROUTINE TO CALCULATE OVERLAP INTEGRALS IN A LOCAL
C    COORDINATE SYSTEM.
C
C    INTEGRALS ARE CALCULATED BY TRANSFORMATION TO ELLIPSOIDAL
C    COORDINATES AND THEREBY EXPRESSED IN TERMS OF C-FUNCTIONS.
C    SEE J.C.P.,24,201. ORIGINALLY WRITTEN BY R.M.STEVENS.
C
C
C    GENERATE FACTORIALS ONLY ONCE.
C
      IF(.NOT.JGO) Then
        JGO=.TRUE.
        FACT(1)=1.D0
        DO I=2,25
          FACT(I)=FACT(I-1)*DFLOAT(I-1)
        End Do
      End If
C
      M2=M1
      STRAD=0.D0
      RHOA=R*SK1
      RHOB=R*SK2
      TERMA=0.5D0**(L1+L2+1) * DSQRT(
     & DFLOAT((L1+L1+1)*(L2+L2+1))
     & *FACT(L1-M1+1)*FACT(L2-M1+1)/
     & (FACT(N1+N1+1)*FACT(N2+N2+1)*FACT(L1+M1+1)*FACT(L2+M1+1))
     & *RHOA**(N1+N1+1)*RHOB**(N2+N2+1))

      JEND=1+((L1-M1)/2)
      KEND=1+((L2-M2)/2)
      IEB=M1+1

      DO J=1,JEND
        JU=J-1
        IAB=N1-L1+JU+JU+1
        ICB=L1-M1-JU-JU+1
        CON1=FACT(L1+L1-JU-JU+1)/(FACT(L1-M1-JU-JU+1)
     &   *FACT(JU+1)*FACT(L1-JU+1))
        DO K=1,KEND
          KU=K-1
          CON12=CON1*FACT(L2+L2-KU-KU+1)/(FACT(L2-M2-KU-KU+1)
     &     *FACT(KU+1)*FACT(L2-KU+1))
          IEV=JU+KU+L2
          IF(2*(IEV/2).NE.IEV) CON12=-CON12
          IBB=N2-L2+KU+KU+1
          IDB=L2-M2-KU-KU+1
          VALUE=0.D0
          DO I6=1,IEB
            DO I5=1,IEB
              VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)
              IEV=I5+I6
              IF(2*(IEV/2).NE.IEV) VALUE1=-VALUE1
              DO I4=1,IDB
                VALUE1=-VALUE1
                VALUE2=BINCOE(IDB,I4)*VALUE1
                DO I3=1,ICB
                  VALUE3=BINCOE(ICB,I3)*VALUE2
                  DO I2=1,IBB
                    VALUE3=-VALUE3
                    VALUE4=BINCOE(IBB,I2)*VALUE3
                    DO I1=1,IAB
                      TERM=VALUE4*BINCOE(IAB,I1)
                      IR=IEB+IEB+ICB+IDB-I3-I4-I6-I6        +I2-1+I1
                      IP=IEB+IEB+ICB+IDB-I3-I4-I5-I5+IAB+IBB-I2+1-I1
                      VALUE=VALUE+A(IP)*B(IR)*TERM
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          STRAD=STRAD+VALUE*CON12
        End Do
      End Do

      STRAD=STRAD*TERMA
      RETURN
      END
