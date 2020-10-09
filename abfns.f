      SUBROUTINE ABFNS(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      Parameter (SMALL=1.0D-7)
      DIMENSION A(*),B(*)
      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL
C
C    SUBROUTINE FOR CALCULATING A AND B FUNCTIONS FOR USE IN LOVLAP.
C
      J=MAXCAL+1
      RHO1=0.5D0*(SK1+SK2)*RR
      RHO2=0.5D0*(SK1-SK2)*RR
      IF((DABS(RHO1).GT.165.D0).OR.(DABS(RHO2).GT.165.D0)) Then
C
C    IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.
C
        DO I=1,J
          A(I)=0.D0
          B(I)=0.D0
        End Do
        Return
      End If
      C=DEXP(-RHO1)
      A(1)=C/RHO1
      DO I=2,J
        A(I)=(DFLOAT(I-1)*A(I-1)+C)/RHO1
      End Do
      IX=J
      IR=DABS(2.D0*RHO2)
      IS=MIN0(IR+1,19)
      IF(RHO2.EQ.0) Then
        DO I=1,IX,2
          B(I)=2.D0/DFLOAT(I)
          B(I+1)=0.D0
        End Do
        Return
      End If

      D=DEXP(RHO2)
      H=1.D0/D
C
C    USE THE DSINH ROUTINE INSTEAD OF SUMMING THE INFINITE SERIES.
C
      !R=2.D0*DSINH(RHO2)
      R=D-H
C
C    AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
C    RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
C
C    pt >> k case. (RHO2 >> K): Upward recursion.
C    pt << k case. (RHO2 << K): Downward recursion.
C
C    Please See JCP.,24,201
C    (See also IJQC 67,199-204 (1998) p202)
C
C    IJQC 63 843(1997) p849.
C    * Upward recursion: from Eq. 50a
C      B(K) = [K*B(K-1) + (-1)**K*EXP(RHO)-EXP(-RHO)]/RHO
C      * Note: 0-based index
C      B(K) = [(K-1)*B(K-1) + (-1)**(K-1)*EXP(RHO)-EXP(-RHO)]/RHO
C      * Note: 1-based index
C    * Downward recursion:
C      B(K-1) = [B(K)*RHO - (-1)**(K-1)*EXP(RHO)+EXP(-RHO)]/(K-1)
C
      B(1)=R/RHO2
      If (IX.LE.RHO2) Then
C
C    Upward recursion case:
C
        DO K=2,IX
          B(K)=((-1)**(K-1)*D-H+DFLOAT(K-1)*B(K-1))/RHO2
        End Do
        Return
      End If

C
C    Using the following procedure to get the B_k first and
C    make B_k-1 series with the Downward recursion
C    for pt << k case. (RHO2 << IX)
C
C    For larger RHO2(beta), Use downward recursion.
C    a) first, calc B(K), then get B(K-1) recursivly.
C
C    The Top B-Function is obtained by summation
C    of the infinite series.
C
      Call GetBk(RHO2,IX,BK,SMALL)
      B(IX)=BK
C
C    After the top B(K) obtained,
C    get the last B(K-1) series by downward recursion.
C
      DO K=IX-1,2,-1
        BK=(BK*RHO2-(-1)**(K)*D+H)/DFLOAT(K)
        B(K)=BK
      End Do
      RETURN
      END

      Subroutine GetBk(RHO,K,BK,SMALL)
      Implicit Real*8 (A-H,O-Z)
C
C    The B-Function can be obtained by summation
C    of the infinite series.
C
C    Almost all B_k obtained within ~20 cycles.
C
      IF(MOD(K,2).EQ.0) Then
        TR=RHO
        BK=-2.D0*TR/DFLOAT(K+1)
        DO J=1,100
          TR=TR*RHO**2/DFLOAT((2*J)*(2*J+1))
          IF(J.GT.5.AND.DABS(TR/BK).LE.SMALL) GoTo 51
          BK=BK-2.D0*TR/DFLOAT(K+1+2*J)
        End Do
      Else
        TR=1.D0
        BK=2.D0*TR/DFLOAT(K)
        DO J=1,100
          TR=TR*RHO**2/DFLOAT((2*J)*(2*J-1))
          IF(J.GT.5.AND.DABS(TR/BK).LE.SMALL) GoTo 51
          BK=BK+2.D0*TR/DFLOAT(K+2*J)
        End Do
      End If
 51   Continue
      Return
      End
