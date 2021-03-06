C
      PROGRAM PREHUCKEL
      PARAMETER (MAXAT=500)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PROGRAM TO CONVERT THE FILES PRODUCED BY FORTICON8 (filename.TMP)
C AND (filename.EXP) TO A SET OF INPUT FILES FOR THE PSI/88 PROGRAM.
C
C ALL ROUTINES OF PSI/88 APPLY EXCEPT FOR PSI1. THE PSI1 PART OF
C THE PROGRAM WAS RECONSTRUCTED FROM PSI/77 AND HAS BEEN PLACED IN
C THE PSI/88 DIRECTORY (PSI1EHT.FOR). RUNNING PSI1 OF PSI/88 WILL
C CAUSE THE PROGRAM TO BOMB SINCE IT IS NOT CONFIGURED FOR EHMO
C CALCULATIONS.
C                                        JJN   9-1-90
C
C HACKED OUT OF PREPLOT (FOR MOPAC) - DAN SEVERANCE, PURDUE - 9/88
C
C FILES USED:
C
C FILE 18 - FORTICON8 USER-DEFINED PARAMETER FILE (VIA LOCAL
C           MODIFICATION) THE FILE IS NAMED filename.EXP
C FILE 13 - FORTICON8 GRAPH FILE (VIA LOCAL MODIFICATION)
C           THE FILE IS NAMED filename.TMP
C FILE 8  - OUTPUT PSI1EHT   INPUT FILE
C FILE 9  - OUTPUT PSICON    INPUT FILE
C FILE 10 - OUTPUT PSI2      INPUT FILE
C
      PARAMETER (MAXVAL=1000)
      PARAMETER (MAXMO=1000)
      PARAMETER (MAXUSR=230)
      COMMON/OTFILE/VECS(MAXMO,MAXVAL),XYZ(MAXAT,3),NAT(MAXAT),
     . ATS(70),SYMWRT(MAXUSR),ATS2(70)
      Character*2 SYMWRT,ATS,ATS2,SYMBL
      COMMON/INFILE/NUMAT,NORBS,NELECS,HOMO,LUMO,ICHG,HOMTWO,NUSER,
     . HOMO2,SYMBL(MAXUSR),VELEC(MAXUSR),NS(MAXUSR),NP(MAXUSR),
     . ND(MAXUSR)
      Common/INFILE2/EXPS(MAXUSR),EXPP(MAXUSR),EXPD(MAXUSR),
     . EXPD2(MAXUSR),C1(MAXUSR),C2(MAXUSR)
      INTEGER*4 NUMAT,NORBS,NELECS,HOMO,NAT,LUMO,ICHG,VELEC,NUSER,HOMO2
      CHARACTER*72 SUBTIT,TITLE
      Character*80 AMO
      Integer AtoI
      DATA ATS/'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',
     . 'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     . 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
     . 'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
     . 'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os',
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra'/
      DATA ATS2/' K','CA','SC','TI',' V','CR','MN','FE','CO','NI',
     . 'CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR',' Y','ZR',
     . 'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE',
     . ' I','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD',
     . 'TB','DY','HO','ER','TM','YB','LU','HF','TA',' W','RE','OS',
     . 'IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA'/
C
C READ FROM DISK THE FOLLOWING DATA FOR GRAPHICS CALCULATION, IN ORDER:
C
C   1.   TITLE
C   2.   NUMBER OF ATOMS, ORBITAL, ELECTRONS, CHARGE
C   3.   ALL ATOMIC NUMBERS AND ATOMIC COORDINATES
C   4.   EIGENVECTORS
C   5.   IF IT EXISTS, READ THE USER-DEFINED PARAMETERS FROM 18.
C
      READ (13,201) TITLE
201   FORMAT (A80)
C
      READ (13,100) NUMAT,NORBS,NELECS,ICHG,NUSER
100   FORMAT (I3,2I4,2I3)   
C
      DO 1001 I=1,NUMAT
      READ (13,101) NAT(I),(XYZ(I,J), J=1,3)
101   FORMAT (I2,3F12.6)
1001  CONTINUE
C
      READ (13,103) ((VECS(I,J),J=1,NORBS), I=1,NORBS)
103   FORMAT (8F10.6)
C
C   CALCULATE THE NUMBER OF THE HOMO.
C
      TOTAL=NORBS
      HOMTWO=NELECS
      HOMTWO=(HOMTWO/2.0)
      HOMO=NINT(HOMTWO)
      HOMO2=(TOTAL-HOMO+1)
C
C
      LUMO=HOMO2-1
20    WRITE(*,*)
      WRITE (*,'(A,I4,A,I4)')' The HOMO is MO number ',HOMO2,
     *  ' The LUMO is MO number ',LUMO
      WRITE(*,*)
      WRITE (*,'(A,$)') ' WHICH MO DO YOU WISH TO PLOT? '
      READ (*,'(A5)') AMO
      Call UpCase(AMO)
      IMO = -1
      If (AMO.EQ.'HOMO') IMO = HOMO2
      If (AMO.EQ.'LUMO') IMO = LUMO
      If (IMO.EQ.-1) IMO = AtoI(AMO)
      Write(*, '('' Selected MO is '',I5)') IMO
      IF (IMO.LT.0.OR.IMO.GT.NORBS) GO TO 20
C
C WRITE PSI1EHT INPUT FILE
C
      WRITE (8,'(A)') 'EHT   AUTO 0    11'
      WRITE (8,108) ICHG,ICHG,ICHG,ICHG
108   FORMAT (4I2)
C
C   READ THE USER-DEFINED PARAMETERS, IF ANY.
C
      IF (NUSER.EQ.0) GO TO 123
      DO 767 J=1,NUSER
      READ (18,768) SYMBL(J),VELEC(J),NS(J),EXPS(J),NP(J),EXPP(J),
     . ND(J),EXPD(J),C1(J),EXPD2(J),C2(J)
768   FORMAT (A2,2I3,F6.3,I3,F6.3,I3,F6.3,F6.4,F6.3,F6.4)
767   CONTINUE
C
C
123   DO 1002 I=1,NUMAT
      WRITE (8,106) NAT(I),(XYZ(I,J), J=1,3)
      IF (NUSER.EQ.0) GO TO 1002
      IF (NAT(I).LE.18) GO TO 1002
      DO 867 K=1,NUSER
      DO 967 L=1,70
      IF (SYMBL(K).NE.ATS2(L)) GO TO 967
      WRITE (8,555) SYMBL(K),VELEC(K),NS(K),EXPS(K),NP(K),EXPP(K),
     . ND(K),EXPD(K),C1(K),EXPD2(K),C2(K)
555   FORMAT (A2,2I3,F6.3,I3,F6.3,I3,F6.3,F6.4,F6.3,F6.4)
967   CONTINUE
      REWIND 18
867   CONTINUE
1002  CONTINUE
106   FORMAT (I2,8X,3F10.6)
      WRITE (8,'(A)') '99'
      WRITE (8,30) (VECS(IMO,J), J=1,NORBS)
30    FORMAT (8F10.6)
      WRITE (8,'(A)') '0101 1.0'
      WRITE (8,'(A)') '01010000'
      WRITE (8,'(A)') '0.005'
C
C WRITE PSICON/88 INPUT FILE
C
      WRITE (9,'(A)') 'EHT'
      WRITE (9,'(A)') '01010001'
      WRITE (9,'(A)') '  0.075000'
C
C WRITE PSI2/88 INPUT FILE
C
      WRITE (10,'(A)') TITLE
      SUBTIT=' '
      IF (IMO.EQ.(HOMO2+2)) SUBTIT = 'HOMO-2'
      IF (IMO.EQ.(HOMO2+1)) SUBTIT = 'HOMO-1'
      IF (IMO.EQ.HOMO2)     SUBTIT = 'HOMO'
      IF (IMO.EQ.LUMO)     SUBTIT = 'LUMO'
      IF (IMO.EQ.(LUMO-1)) SUBTIT = 'LUMO+1'
      IF (IMO.EQ.(LUMO-2)) SUBTIT = 'LUMO+2'
      WRITE (10,'(A)') SUBTIT
      WRITE (10,'(A/A)') '01','00'
      WRITE (10,'(A)') ' '
      DO 1003 I=1,NUMAT
      WRITE (10,60) NAT(I),(XYZ(I,J),J=1,3)
60    FORMAT (I2,8X,3F10.6)
1003  CONTINUE
      WRITE (10,'(A)') '99'
      X = 10.0
      SCALE = 0.7
      WRITE (10,'(4F10.6)') X,X,X,SCALE
      DO 1011 I=1,NUMAT
      SYMWRT(I)='XX'
1011  CONTINUE
      DO 1005 I=1,NUMAT
      K=0
      DO 1004 J=19,89
      K=K+1
      IF (NAT(I).EQ.J) SYMWRT(I)=ATS(K)
1004  CONTINUE
1005  CONTINUE
      L=1
      DO 1006 I=1,NUMAT
      IF (SYMWRT(I).EQ.'XX') GO TO 1006
      L=L+1
      WRITE (10,1007) SYMWRT(I),L
1007  FORMAT (A2,I2)
1006  CONTINUE
      WRITE (10,'(A)') '02'
      STOP
      END

      Integer Function AtoI(A)
      Character*80 A
      Integer I, J, K
      Integer Zero, Nine
      Logical IsNumber
      Zero = IChar('0')
      Nine = IChar('9')
      IsNumber = .FALSE.
      J = 0
      Do I = 1, 80
          K = IChar(A(I:I))
          If (K.GE.Zero.AND.K.LE.Nine) Then
              J = J * 10 + K - Zero
              IsNumber = .True.
          Else If (IsNumber) Then
              Goto 100
          End If
      End Do

 100  AtoI = J
      Return
      End

      Subroutine UpCase(KeyWrd)
      Character*80 KeyWrd
      ICapA = IChar('A')
      ILowA = IChar('a')
      ILowZ = IChar('z')
      Do I=1, 80
         ILine = IChar(KeyWrd(I:I))
         If (ILine.GE.ILowA.AND.ILine.LE.ILowZ) Then
            KeyWrd(I:I) = CHAR(ILine + ICapA - ILowA)
         EndIf
      End Do
      Return
      End
C
C
