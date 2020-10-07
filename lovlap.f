      SUBROUTINE LOVLAP (STRAD,A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MAXBC=7)
      PARAMETER (MAXFAC=24)
      DIMENSION A(*),B(*)
C    Factorials.
      DIMENSION FACT(0:MAXFAC)
C
      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX
      DIMENSION BINCOE(MAXBC,MAXBC)
C
      LOGICAL*1 JGO
      DATA JGO/.FALSE./
      SAVE FACT,JGO
      SAVE BINCOE
C
C    SUBROUTINE TO CALCULATE OVERLAP INTEGRALS IN A LOCAL
C    COORDINATE SYSTEM.
C
C    INTEGRALS ARE CALCULATED BY TRANSFORMATION TO ELLIPSOIDAL
C    COORDINATES AND THEREBY EXPRESSED IN TERMS OF C-FUNCTIONS.
C    SEE J.C.P.,24,201. ORIGINALLY WRITTEN BY R.M.STEVENS.
C
C    (C-functions can be expressed in terms of the A,B functions)
C
C    Please see J.C.P.,24,201. Ruedenberg Roothaan and Jaunzemis
C
C    Study of Two-Center Integrals Useful in Calculation on Molecular
C    Structure. III. A Unified Treatment of the Hybrid, Coulomb, and
C    One-Electron Integrals
C
C    Chapter 2.b. One-Electron Integrals Expressible in Terms of
C    C-Functions
C
C    1) nuclear-attraction integrals (I. Eq.6) and 2) the kinetic-energy
C    integrals (I. Eq.5) can be expressed in terms of overlap integrals
C
C      (https://aip.scitation.org/doi/abs/10.1063/1.1742457)
C
C
C    GENERATE FACTORIALS ONLY ONCE.
C
      IF(.NOT.JGO) Then
        JGO=.TRUE.
        FACT(0)=1.D0
        DO I=1,MAXFAC
          FACT(I)=FACT(I-1)*DFLOAT(I)
        End Do
C
C    Generate Binomial-Coefients.
C
        DO I=1,MAXBC
          BINCOE(I,1)=1.D0
          BINCOE(I,I)=1.D0
        End Do
        DO I=1,MAXBC-1
          I1=I-1
          DO J=1,I1
            BINCOE(I+1,J+1)=BINCOE(I,J+1)+BINCOE(I,J)
          End Do
        End Do
      End If
C
      M2=M1
      STRAD=0.D0
      RHOA=R*SK1
      RHOB=R*SK2
      TERMA=0.5D0**(L1+L2+1) * DSQRT(
     & DFLOAT((L1+L1+1)*(L2+L2+1))
     & *FACT(L1-M1)*FACT(L2-M1)/
     & (FACT(N1+N1)*FACT(N2+N2)*FACT(L1+M1)*FACT(L2+M1))
     & *RHOA**(N1+N1+1)*RHOB**(N2+N2+1))

      JEND=1+((L1-M1)/2)
      KEND=1+((L2-M2)/2)
      IEB=M1+1

      DO J=1,JEND
        JU=J-1
        IAB=N1-L1+JU+JU+1
        ICB=L1-M1-JU-JU+1
        CON1=FACT(L1+L1-JU-JU)/(FACT(L1-M1-JU-JU)
     &   *FACT(JU)*FACT(L1-JU))
        DO K=1,KEND
          KU=K-1
          CON12=CON1*FACT(L2+L2-KU-KU)/(FACT(L2-M2-KU-KU)
     &     *FACT(KU)*FACT(L2-KU))
          IEV=JU+KU+L2
          IF(Mod(IEV,2).NE.0) CON12=-CON12
          IBB=N2-L2+KU+KU+1
          IDB=L2-M2-KU-KU+1
          VALUE=0.D0
          DO I6=1,IEB
            DO I5=1,IEB
              VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)
              IEV=I5+I6
              IF(Mod(IEV,2).NE.0) VALUE1=-VALUE1
              DO I4=1,IDB
                VALUE1=-VALUE1
                VALUE2=BINCOE(IDB,I4)*VALUE1
                DO I3=1,ICB
                  VALUE3=BINCOE(ICB,I3)*VALUE2
                  DO I2=1,IBB
                    VALUE3=-VALUE3
                    VALUE4=BINCOE(IBB,I2)*VALUE3
                    III=IEB+IEB+ICB+IDB-I3-I4
                    IR =III-I6-I6        +I2-1
                    IP =III-I5-I5+IAB+IBB-I2+1
                    DO I1=1,IAB
                      TERM=VALUE4*BINCOE(IAB,I1)
                      IR=IR+1
                      IP=IP-1
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
