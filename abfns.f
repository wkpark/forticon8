      SUBROUTINE ABFNS(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MXUSER=230)
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
        DO I=1,MXUSER
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
      R=2.D0*DSINH(RHO2)
C
C    AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
C    RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
C
      B(1)=R/RHO2
      DO I=2,IX,IS
        IF(IR.NE.0) Then
C
C    MODIFICATION TO AVOID EXCEEDING STORAGE LIMITS.
C    D. WALLACE 04/14/71
C
          DO K=I,IX
            !IF(((-1)**K).LT.0) Then
            IF(MOD(K,2).NE.0) Then
              B(K)=( R  +DFLOAT(K-1)*B(K-1))/RHO2
            Else
              B(K)=(-D-H+DFLOAT(K-1)*B(K-1))/RHO2
            End If
          End Do
        End If
        IN=I+IS-1
C
C    AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE
C    NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
C    OF THE INFINITE SERIES.
C
        IF(IN.GT.IX) Return
        !IF((-1)**IN.GT.0) Then
        IF(MOD(IN,2).EQ.0) Then
          TR=RHO2
          B(IN)=-2.D0*TR/DFLOAT(IN+1)
          DO J=1,100
            TR=TR*RHO2**2/DFLOAT((2*J)*(2*J+1))
            IF(DABS(TR/B(IN)).LE.SMALL) GoTo 51
            B(IN)=B(IN)-2.D0*TR/DFLOAT(IN+1+2*J)
          End Do
        Else
          TR=1.D0
          B(IN)=2.D0*TR/DFLOAT(IN)
          DO J=1,100
            TR=TR*RHO2**2/DFLOAT((2*J)*(2*J-1))
            IF(DABS(TR/B(IN)).LE.SMALL) GoTo 51
            B(IN)=B(IN)+2.D0*TR/DFLOAT(IN+2*J)
          End Do
        End If
 51     CONTINUE
      End Do
      RETURN
      END
