C FORTICON8 (VAX VERSION)  QCPE 517                                         
C
C  THE FOLLOWING CODE IS NOT THE ORIGINAL SOURCE CODE. THIS CODE 
C  INCLUDES MODIFICATIONS TO WRITE OUT PERTINENT INFORMATION TO A
C  TEMPORARY FILE (FOR013) TO BE USED IN MO PLOTTING ROUTINES.
C
C  ALL MODIFICATIONS ARE COMMENTED TO PROVIDE EASE IN IDENTIFYING
C  CHANGES.      JJN   8-7-90
C
C********************************************************************** FORT0001
C                                                                     * FORT0002
C    PROGRAM FORTICON8  COMPLETE FORTRAN VERSION OF ICON8             * FORT0003
C                                                                     * FORT0004
C    THE FOLLOWING SUBPROGRAMS, WHICH EXIST AS ASSEMBLER ROUTINES     * FORT0005
C    IN ICON8, ARE TRANSLATED INTO FORTRAN0  MATRIX, ABFNS,           * FORT0006
C    LOVLAP, GRMSCH, TRNFRM, DSUM, ROTATE, DOT, VECSUM, REDUCE,       * FORT0007
C    FULCHM, AND REDCHM.  FORTICON INLUDES THESE AS WELL AS ALL       * FORT0008
C    THE FORTRAN SUBPROGRAMS OF ICON8.                                * FORT0009
C                                                                     * FORT0010
C********************************************************************** FORT0011
C                                                                       FORT0012
           IMPLICIT REAL*8(A-H,O-Z)                                     FORT0013
C                                                                       FORT0014
C     PROGRAM ICON FOR PERFORMING EXTENDED HUCKEL CALCULATIONS          FORT0015
C     WITH OR WITHOUT CHARGE ITERATION.                                 FORT0016
C     ** QCPE VERSION **                                                FORT0017
C                                                                       FORT0018
C     ** SAMPLE DECK **                                                 FORT0019
C ....0....1....0....2....0....3....0....4....0....5....0....6....0     FORT0020
C                                                                       FORT0021
C ETHYLENE                                                              FORT0022
C   4  2        2                                                       FORT0023
C  0.92665        1.205          0.0                                    FORT0024
C -0.92665        1.205          0.0                                    FORT0025
C  0.92665       -1.205          0.0                                    FORT0026
C -0.92665       -1.205          0.0                                    FORT0027
C  0.0            0.67           0.0                                    FORT0028
C  0.0           -0.67           0.0                                    FORT0029
C  C C                                                                  FORT0030
C                                                                       FORT0031
C
C    REVISED TO ALLOW MORE FLEXIBILITY IN NUMBER OF ATOMS
C    POSSIBLE. MAXATM (MAXIMUM NUMBER OF ATOMS) SET AT 500.
C    MAXIMUM NUMBER OF HEAVY ATOMS SET AT 250. MAXIMUM NUMBER
C    OF USER-DEFINED ATOMS SET TO 230.
C    JJN   8-28-90
C
	   PARAMETER (MAXATM=500)
	   PARAMETER (BB=250)
           PARAMETER (MXUSER=230)
           PARAMETER (MXUSR2=231)
C
C   THE ABOVE PARAMETERS ARE DEFINED AS SUCH TO ELIMINATE THE NEED
C   TO SEARCH THROUGH THE CODE FOR EVERY OCCURRENCE OF THE NUMBERS
C   (AS I HAD TO DO WITH THE ORIGINAL CODE). THE PARAMETER STATEMENTS
C   ARE DEFINED IN EACH SUBROUTINE REQUIRING THEM SO IF THESE NEED TO
C   BE CHANGED AGAIN, JUST FIND THE PARAMETER STATEMENTS. ALSO, THERE
C   ARE A FEW PLACES WHERE THE USE OF THE PARAMETERS IS NOT ALLOWED AND
C   THE ACTUAL NUMBER IS USED. MOST NOTABLY IN DATA STATEMENTS. A QUICK
C   SEARCH THROUGH FOR THE VALUES OF MXUSER AND MXUSR2 WILL BE
C   NECESSARY.       JJN  9-8-90
C
           COMMON/TITLE/AB(10)                                          FORT0032
           COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,        FORT0033
     .     IPRINT,IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                   FORT0034
           LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                       FORT0035
           COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)             FORT0036
           LOGICAL*1 PRT,PUN                                            FORT0037
           INTEGER*2 IOVPOP,IENRGY                                      FORT0038
           COMMON/ATOM/AC(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),
     .     ND(BB),EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),
     .     COULS(BB),COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)
           INTEGER*2 AC,SYMBOL,VELEC                                    FORT0042
           COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,  FORT0043
     .     PRTCYC,ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)                  FORT0044
           REAL*8 LAMPRI                                                FORT0045
           INTEGER*4 PRTCYC                                             FORT0046
           LOGICAL*1 PARTIT,PRINTX,ITABLE                               FORT0047
           COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5), FORT0048
     .     BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),     FORT0049
     .     AD3(5),BD3(5),CD3(5)                                         FORT0050
           COMMON/STARS/STAR,STAR2                                      FORT0051
           INTEGER*2 STAR,STAR2                                         FORT0052
C                                                                       FORT0053
C     CALL INPUT TO READ IN AND PRINT OUT INPUT DATA.                   FORT0054
C                                                                       FORT0055
1000       CALL INPUT(NATOM,NDIM,NTYPE)                                 FORT0056
           IF(IPRINT.LT.-3) GO TO 2000                                  FORT0057
C                                                                       FORT0058
C     CALCULATE MATRIX DIMENSIONS.                                      FORT0059
C                                                                       FORT0060
           NC=NDIM*(NDIM+1)/2                                           FORT0061
           NHS=NC+NC-NDIM                                               FORT0062
           NSS=NHS                                                      FORT0063
C                                                                       FORT0064
C     IF NO WAVEFUNCTIONS NEEDED, SET NHS=0. THIS HAS THE EFFECT        FORT0065
C     OF EQUIVALENCING THE H AND S MATRICES.                            FORT0066
C                                                                       FORT0067
           ONEMAT=.FALSE.                                               FORT0068
           IF(.NOT.ITERAT.AND.IPRINT.LT.-1) ONEMAT=.TRUE.               FORT0069
           IF(IPUNCH.NE.0) GO TO 600                                    FORT0070
           DO 400 I=6,20                                                FORT0071
           IF(PUN(I)) GO TO 600                                         FORT0072
400        CONTINUE                                                     FORT0073
           GO TO 500                                                    FORT0074
600        ONEMAT=.FALSE.                                               FORT0075
500        IF(ONEMAT) NHS=0                                             FORT0076
           IF(METH.LT.3) NTYPE=1                                        FORT0077
           NMD=NTYPE*NTYPE                                              FORT0078
           NCL=NDIM                                                     FORT0079
           IF(NDIM.LT.10) NCL=10                                        FORT0080
           NOCC=(NDIM+1)/2                                              FORT0081
           NHDG=1                                                       FORT0082
           IF(METH.GT.2.AND.L5) NHDG=NDIM                               FORT0083
C                                                                       FORT0084
C     CALL MATRIX TO ALLOCATE SPACE FOR MATRICES.                       FORT0085
C     ORDER OF MATRICES0 H S MAD C SP PD MAXS MAXP MAXD                 FORT0086
C                        COUL0 SORB IOCC HDG                            FORT0087
C                                                                       FORT0088
200        CALL MATRIX(13,NHS,NSS,NMD,NC,NDIM,NDIM,NDIM,NDIM,2*NDIM,    FORT0089
     .     NCL,NDIM,NOCC,NHDG,     NDIM,NDIM,NC,NATOM,NTYPE,NHDG,NC,NC, FORT0090
     .     NC,NC,NC,NC,NDIM)                                            FORT0091
2000       CONTINUE                                                     FORT0092
           GO TO 1000                                                   FORT0093
           END                                                          FORT0094
           BLOCK DATA                                                   FORT0095
C                                                                       FORT0096
C     INITIALIZATION OF INTERNAL ATOMIC DATA. THERE ARE PROVISIONS      FORT0097
C     FOR 20 USER DEFINED ATOMS, 15 INTERNALLY DEFINED ATOMS, AND       FORT0098
C     SPACE FOR 5 MORE TO BE USED EITHER WAY.                           FORT0099
C                                                                       FORT0100
C
C     THIS HAS BEEN MODIFIED TO ALLOW 230 USER DEFINED ATOMS, 15
C     INTERNALLY DEFINED ATOMS AND SPACE FOR 5 MORE TO BE USED EITHER
C     WAY (TOTAL=250).   JJN  8-28-90
C
	   PARAMETER (MAXATM=500)
	   PARAMETER (BB=250)           
           PARAMETER (MXUSER=230)
           PARAMETER (MXUSR2=231)
	   IMPLICIT REAL*8(A-H,O-Z)                                
           COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),
     .     ND(BB),EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),
     .     C2(BB),COULS(BB),COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),
     .     Z(MAXATM)
           INTEGER*2 KEY,SYMBOL,VELEC                              
           COMMON/STARS/STAR,STAR2                                  
           INTEGER*2 STAR,STAR2                                     
           COMMON/START/NUSER                                        
           DATA SYMBOL/230*'**',' C',' N',' O',' F','SI',' P',' S',
     .     'CL','LI','BE',' B','NA','MG','AL',' H',5*'  '/    
           DATA VELEC/230*0,4,5,6,7,4,5,6,7,1,2,3,1,2,3,1,5*0/  
           DATA NS/230*0,4*2,4*3,3*2,3*3,1,5*0/                      
           DATA EXPS/230*0.0D0,1.625D0,1.950D0,2.275D0,2.425D0,1.383D0,
     .     1.6D0,1.817D0,2.033D0,0.65D0,0.975D0,1.3D0,0.733D0,0.95D0,   FORT0114
     .     1.167D0,1.3D0,5*0.0D0/                                    
           DATA COULS/230*0.0D0,-21.4D0,-26.0D0,-32.3D0,-40.0D0,-17.3D0,
     .     -18.6D0,-20.0D0,-30.0D0,-5.4D0,-10.0D0,-15.2D0,-5.1D0,-9.0D0,
     .     -12.3D0,-13.6D0,5*0.0D0/                               
           DATA NP/230*0,4*2,4*3,3*2,3*3,6*0/                       
           DATA EXPP/230*0.0D0,1.625D0,1.950D0,2.275D0,2.425D0,1.383D0,
     .     1.6D0,1.817D0,2.033D0,0.65D0,0.975D0,1.3D0,0.733D0,0.95D0,   FORT0121
     .     1.167D0,6*0.0D0/                                         
           DATA COULP/230*0.0D0,-11.4D0,-13.4D0,-14.8D0,-18.1D0,-9.2D0,
     .     -14.0D0,-13.3D0,-15.0D0,-3.5D0,-6.0D0,-8.5D0,-3.0D0,-4.5D0,  FORT0124
     .     -6.5D0,6*0.0D0/                                        
           DATA ND/234*0,4*3,12*0/                                  
           DATA EXPD/234*0.0D0,1.383D0,1.4D0,1.5D0,2.033D0,12*0.0D0/ 
           DATA COULD/234*0.0D0,-6.0D0,-7.0D0,-8.0D0,-9.0D0,12*0.0D0/
           DATA C1/250*0.0D0/                                          
           DATA C2/250*0.0D0/                                        
           DATA EXPD2/250*0.0D0/                                 
           DATA STAR,STAR2 / ' *','**'/                                 FORT0132
C
C    PROVISION FOR 230 USER DEFINED ATOMS.   JJN 8-28-90
C
           DATA NUSER/231/                                      
C
           END                                                          FORT0134
           SUBROUTINE INPUT(NATOM,NDIM,NTYPE)                           FORT0135
C                                                                       FORT0136
C     SUBROUTINE FOR READING IN AND PRINTING OUT INPUT DATA.            FORT0137
C                                                                       FORT0138
	   PARAMETER (MAXATM=500)
	   PARAMETER (BB=250)
           PARAMETER (MXUSER=230)
           PARAMETER (MXUSR2=231)
           IMPLICIT REAL*8(A-H,O-Z)                                     FORT0139
           COMMON/TITLE/AB(10)                                          FORT0140
           COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,        FORT0141
     .     IPRINT,IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                   FORT0142
           LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                       FORT0143
           COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)             FORT0144
           LOGICAL*1 PRT,PUN                                            FORT0145
           INTEGER*2 IOVPOP,IENRGY                                      FORT0146
           COMMON/ATOM/AC(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),
     .     ND(BB),EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),
     .     C2(BB),COULS(BB),COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),
     .     Z(MAXATM),SYMBL(MAXATM),ATS(89)
           INTEGER*2 AC,SYMBOL,SYMBL,VELEC,ATS                        
           COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,  FORT0151
     .     PRTCYC,ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)                  FORT0152
           REAL*8 LAMPRI                                                FORT0153
           INTEGER*4 PRTCYC                                             FORT0154
           LOGICAL*1 PARTIT,PRINTX,ITABLE                               FORT0155
C                                                                       FORT0156
C     SINCE THE ITERNAL ATOMIC PARAMETERS ( EXPS, EXPP, ETC. ) ARE      FORT0157
C     NOT USED WHEN DOING CHARGE ITERATION ( METH >1 ) THE SPACE        FORT0158
C     ALLOCATED TO THEM CAN BE USED FOR THE VSIE CHARGE ITERATION       FORT0159
C     PARAMETERS.                                                       FORT0160
C                                                                       FORT0161
           DIMENSION AS1(MXUSER),BS1(MXUSER),CS1(MXUSER),AP1(MXUSER),
     .     BP1(MXUSER),CP1(MXUSER),AD1(MXUSER),BD1(MXUSER),CD1(MXUSER)         
           EQUIVALENCE (AS1(1),EXPS(MXUSR2)),(BS1(1),EXPP(MXUSR2)),  
     .     (CS1(1),EXPD(MXUSR2)),(AP1(1),EXPD2(MXUSR2)),
     .     (BP1(1),C1(MXUSR2)),(CP1(1),C2(MXUSR2)),
     .     (AD1(1),COULS(MXUSR2)),(BD1(1),COULP(MXUSR2)), 
     .     (CD1(1),COULD(MXUSR2))                                   
           COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5), FORT0168
     .     BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),     FORT0169
     .     AD3(5),BD3(5),CD3(5)                                         FORT0170
           INTEGER*2 CHANGE(MXUSER)                           
           EQUIVALENCE (CHANGE(1),AS2(1))                               FORT0172
           REAL*4 MADS(MXUSER),MADP(MXUSER),MADD(MXUSER)      
           EQUIVALENCE (MADS(1),NS(MXUSR2)),(MADP(1),NP(MXUSR2)),
     .     (MADD(1),ND(MXUSR2))                              
           COMMON/STARS/STAR,STAR2                                      FORT0176
           INTEGER*2 STAR,STAR2                                         FORT0177
           COMMON/START/NUSER                                           FORT0178
           DIMENSION EXTRA(9)                                           FORT0179
           EQUIVALENCE (X(1),EXTRA(1))                                  FORT0180
           INTEGER*2 CONTIN                                             FORT0181
           EQUIVALENCE (AB(10),CONTIN)                                  FORT0182
           INTEGER*2 HYDROG                                             FORT0183
           DATA HYDROG/' H'/                                            FORT0184
           DATA ATS/' H','HE','LI','BE',' B',' C',' N',' O',' F',
     .      'NE','NA','MG','AL','SI',' P',' S','CL','AR',' K','CA',
     .      'SC','TI',' V','CR','MN','FE','CO','NI','CU','ZN','GA',
     .      'GE','AS','SE','BR','KR','RB','SR',' Y','ZR','NB','MO',
     .      'TC','RU','RH','PD','AG','CD','IN','SN','SB','TE',' I',
     .      'XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD',
     .      'TB','DY','HO','ER','TM','YB','LU','HF','TA',' W','RE',
     .      'OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN',
     .      'FR','RA','AC'/
C                                                                       FORT0185
C     READ AND WRITE TITLE.                                             FORT0186
C     IF CONTIN IS EQUAL TO STAR THEN ANOTHER TITLE CARD WILL           FORT0187
C     BE READ AND PRINTED. HOWEVER ONLY THE FIRST IS STORED             FORT0188
C     FOR PRINTING LATER ON. GET YOUR GOODIES ON THE FIRST.             FORT0189
C                                                                       FORT0190
1000       READ(5,1,END=115) AB                                         FORT0191
1          FORMAT(8A8,A6,A2)                                            FORT0192
C     WRITE TITLE TO DISK FILE 13.   JJN  9-3-90
C
           WRITE(13,2730) AB
2730       FORMAT(8A8,A6,A2)
C
           WRITE(6,2) AB                                                FORT0193
2          FORMAT('1',T10,8A8,A6,A2)                                    FORT0194
11         IF(CONTIN.NE.STAR) GO TO 9                                   FORT0195
           READ(5,1) EXTRA,CONTIN                                       FORT0196
           WRITE(6,12) EXTRA,CONTIN                                     FORT0197
12         FORMAT(T10,8A8,A6,A2)                                        FORT0198
           GO TO 11                                                     FORT0199
9          CONTINUE                                                     FORT0200
C                                                                       FORT0201
C     READ PARAMETER CARD.                                              FORT0202
C                                                                       FORT0203
           READ(5,3) NH,NA,KA,METH,IPRINT,IPUNCH,L1,L2,L3,L4,L5,CON,    FORT0204
     .     PEEP,COULH,(PRT(I),I=1,20),(PUN(J),J=1,20)                   FORT0205
3          FORMAT(6I3,5L1,F5.2,2F6.3,40L1)                              FORT0206
C                                                                       FORT0207
C     INSERT DEFAULT PARAMETERS.                                        FORT0208
C                                                                       FORT0209
           IF(CON.LT.1.E-05) CON=1.75D0                                 FORT0210
           IF(PEEP.LT.1.E-05) PEEP=1.3D0                                FORT0211
           IF(COULH.GT.-1.E-05) COULH=-13.6D0                           FORT0212
           ITERAT=METH.NE.0                                             FORT0213
           NATOM=NH+NA                                                  FORT0214
           NATM=NATOM                                                   FORT0215
C                                                                       FORT0216
C     SET IPRINT OPTION.                                                FORT0217
C                                                                       FORT0218
           IF(IPRINT.GT.1) GO TO 250                                    FORT0219
           PRT(6)=.TRUE.                                                FORT0220
           PRT(7)=.TRUE.                                                FORT0221
           PRT(11)=.TRUE.                                               FORT0222
           PRT(17)=.TRUE.                                               FORT0223
           PRT(19)=.TRUE.                                               FORT0224
           IF(IPRINT.GT.0) GO TO 250                                    FORT0225
           PRT(12)=.TRUE.                                               FORT0226
           PRT(14)=.TRUE.                                               FORT0227
           PRT(20)=.TRUE.                                               FORT0228
           IF(IPRINT.GT.-1) GO TO 250                                   FORT0229
           PRT(13)=.TRUE.                                               FORT0230
           PRT(15)=.TRUE.                                               FORT0231
           PRT(16)=.TRUE.                                               FORT0232
           PRT(18)=.TRUE.                                               FORT0233
           IF(IPRINT.GT.-2) GO TO 250                                   FORT0234
           PRT(10)=.TRUE.                                               FORT0235
C                                                                       FORT0236
C     READ COORDINATES AND HEAVY ATOM CARD.                             FORT0237
C                                                                       FORT0238
250        READ(5,5) (X(I),Y(I),Z(I),I=1,NATOM)                         FORT0239
5          FORMAT(3F15.6)                                               FORT0240
           READ(5,8) (AC(I),I=1,NA)                                     FORT0241
8          FORMAT(40A2)                                                 FORT0242
C                                                                       FORT0243
C     READ AND DECODE ATOM DEFINITION CARDS.                            FORT0244
C                                                                       FORT0245
           JOHN=0
           NDIM=NH                                                      FORT0246
           NTYPE=NH                                                     FORT0247
           NELEC=NH-KA                                                  FORT0248
           K=NUSER                                                      FORT0249
           NUSER2=BB                                              
           IF(METH.GE.2) NUSER2=MXUSER                                   
           DO 100 I=1,NA                                                FORT0252
           IF(NUSER.GT.NUSER2) GO TO 103                                FORT0253
           DO 102 J=NUSER,NUSER2                                        FORT0254
           JSAVE=J                                                      FORT0255
           IF(AC(I) .EQ. SYMBOL(J)) GO TO 101                           FORT0256
102        CONTINUE                                                     FORT0257
C                                                                       FORT0258
C     PROVISION FOR USER SPECIFIED DATA.                                FORT0259
C                                                                       FORT0260
           IF(AC(I).EQ.STAR) GO TO 103                                  FORT0261
           IF(AC(I).EQ.STAR2) GO TO 105                                 FORT0262
           WRITE(6,6) I,AC(I)                                           FORT0263
6          FORMAT(//,T10,'HEAVY ATOM',I3,' NOT RECOGNIZED. SYMBOL0',A2) FORT0264
           IF(METH.GE.2) WRITE(6,13)                                    FORT0265
13         FORMAT(/,T10,'REMEMBER IF USING METH > 1 ALL ATOMIC',        FORT0266
     .     ' PARAMETERS MUST BE DEFINED BY THE USER.')                  FORT0267
115        REWIND 7                                                     FORT0268
           STOP                                                         FORT0269
103        NUSER=NUSER-1                                                FORT0270
105        READ(5,7) SYMBOL(NUSER),VELEC(NUSER),NS(NUSER),EXPS(NUSER),  FORT0271
     .     COULS(NUSER),NP(NUSER),EXPP(NUSER),COULP(NUSER),ND(NUSER),   FORT0272
     .     EXPD(NUSER),COULD(NUSER),C1(NUSER),EXPD2(NUSER),C2(NUSER)    FORT0273
7          FORMAT(A2,I3,3(I3,2F6.3),F6.4,F6.3,F6.4)
C
C    THIS SECTION WRITE OUT THE USER-DEFINED PARAMETERS TO DISK
C    FILE 18. THESE WILL BE USED FOR CONSTRUCTION OF THE PSI1
C    INPUT FILE AS THESE PARAMETERS ARE REQUIRED FOR ELEMENTS
C    GREATER THAN ATOMIC NUMBER 18.    JJN  9-8-90
C
           JOHN=JOHN+1
           WRITE(18,767) SYMBOL(NUSER),VELEC(NUSER),NS(NUSER),
     .      EXPS(NUSER),NP(NUSER),EXPP(NUSER),ND(NUSER),EXPD(NUSER),
     .      C1(NUSER),EXPD2(NUSER),C2(NUSER)
767        FORMAT(A2,I3,I3,F6.3,I3,F6.3,I3,F6.3,F6.4,F6.3,F6.4)
C
C
           JSAVE=NUSER                                                  FORT0275
C                                                                       FORT0276
C     NORMALIZE USER SPECIFIED CONTRACTED D ORBITAL.                    FORT0277
C                                                                       FORT0278
           IF(C2(NUSER).EQ.0.) GO TO 101                                FORT0279
           S=(4.D0*EXPD(NUSER)*EXPD2(NUSER)/(EXPD(NUSER)+EXPD2(NUSER)   FORT0280
     .     )**2)**(ND(NUSER)+.5D0)                                      FORT0281
           S=1.D0/DSQRT(C1(NUSER)**2+C2(NUSER)**2+(S+S)*C1(NUSER)       FORT0282
     .     *C2(NUSER))                                                  FORT0283
           C1(NUSER)=S*C1(NUSER)                                        FORT0284
           C2(NUSER)=S*C2(NUSER)                                        FORT0285
101        NELEC=NELEC+VELEC(JSAVE)                                     FORT0286
C                                                                       FORT0287
C     AC, LATER REFERENCED AS KEY, IS A POINTER TO THE PARAMETER TABLES.FORT0288
C                                                                       FORT0289
111        AC(I)=JSAVE                                                  FORT0290
           NDIM=NDIM+4                                                  FORT0291
           IF(NP(JSAVE).EQ.0) NDIM=NDIM-3                               FORT0292
           IF(ND(JSAVE).NE.0) NDIM=NDIM+5                               FORT0293
           NTYPE=NTYPE+2                                                FORT0294
           IF(NP(JSAVE).EQ.0) NTYPE=NTYPE-1                             FORT0295
           IF(ND(JSAVE).NE.0) NTYPE=NTYPE+1                             FORT0296
100        CONTINUE                                                     FORT0297
C                                                                       FORT0298
C     READ IN CHARGE ITERATION PARAMETERS. SET DEFAULT VALUES.          FORT0299
C                                                                       FORT0300
           IF(.NOT.ITERAT) GO TO 60                                     FORT0301
           IF(K.NE.MXUSR2) GO TO 60                                  
           READ(5,61) DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,     FORT0303
     .     PRTCYC,NCON,PARTIT                                           FORT0304
61         FORMAT(6F10.5,3I5,4X,L1)                                     FORT0305
           IF(.NOT.PARTIT.OR.METH.EQ.2) GO TO 65                        FORT0306
           WRITE(6,66)                                                  FORT0307
66         FORMAT(///,T10,'PARTIAL ITERATION ( PARTIT = TRUE ) MAY',    FORT0308
     .     ' ONLY BE USED IF METH = 2.')                                FORT0309
           STOP                                                         FORT0310
65         IF(DELTAC.EQ.0.0D0) DELTAC=0.0001D0                          FORT0311
           IF(SENSE.EQ.0.0D0) SENSE=2.0D0                               FORT0312
           IF(MAXCYC.EQ.0) MAXCYC=25                                    FORT0313
           IF(PRTCYC.EQ.0) PRTCYC=MAXCYC                                FORT0314
           IF(NCON.EQ.0) NCON=3                                         FORT0315
           IF(DAMP1.EQ.0.0D0) DAMP1=0.1D0                               FORT0316
           IF(METH.GE.3) GO TO 62                                       FORT0317
           IF(DAMP2.EQ.0.0D0) DAMP2=0.25D0                              FORT0318
           IF(LAMPRI.EQ.0.0D0) LAMPRI=0.25D0                            FORT0319
           GO TO 63                                                     FORT0320
62         IF(DAMP2.EQ.0.0D0) DAMP2=0.75D0                              FORT0321
           IF(LAMPRI.EQ.0.0D0) LAMPRI=0.75D0                            FORT0322
63         IF(METH.LT.2) GO TO 60                                       FORT0323
           DO 32 I=1,20                                                 FORT0324
32         ITABLE(I)=.FALSE.                                            FORT0325
           NUSER2=MXUSR2-NUSER                                      
C                                                                       FORT0327
C     READ IN SYMBOLS OF ATOMS ON WHICH CHARGE ITERATION                FORT0328
C     IS TO BE PERFORMED.                                               FORT0329
C                                                                       FORT0330
           IF(.NOT.PARTIT) GO TO 30                                     FORT0331
           READ(5,31) CHANGE                                            FORT0332
31         FORMAT(20A2)                                                 FORT0333
           DO 33 I=1,NUSER2                                             FORT0334
           J=MXUSR2-I                                                   
           DO 33 K=1,NUSER2                                             FORT0336
           IF(SYMBOL(J).EQ.CHANGE(K)) ITABLE(J)=.TRUE.                  FORT0337
33         CONTINUE                                                     FORT0338
           GO TO 34                                                     FORT0339
30         DO 35 I=1,NUSER2                                             FORT0340
           J=MXUSR2-I                                                 
35         ITABLE(J)=.TRUE.                                             FORT0342
C                                                                       FORT0343
C     READ IN VSIE AND MADELUNG PARAMETERS.                             FORT0344
C                                                                       FORT0345
34         DO 36 I=1,NUSER2                                             FORT0346
           J=MXUSR2-I                                                
           IF(.NOT.ITABLE(J)) GO TO 36                                  FORT0348
           READ(5,37) AS1(I),BS1(I),CS1(I),MADS(I)                      FORT0349
37         FORMAT(4F10.8)                                               FORT0350
           IF(NP(J).EQ.0) GO TO 36                                      FORT0351
           IF(ND(J).NE.0) GO TO 38                                      FORT0352
           READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I)                      FORT0353
           GO TO 36                                                     FORT0354
38         IF(NCON.EQ.3) GO TO 39                                       FORT0355
           READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I),AD1(I),BD1(I),       FORT0356
     .     CD1(I),MADD(I)                                               FORT0357
           GO TO 36                                                     FORT0358
39         READ(5,40) AS2(I),BS2(I),CS2(I),AS3(I),BS3(I),CS3(I)         FORT0359
40         FORMAT(3F10.8)                                               FORT0360
           READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I)                      FORT0361
           READ(5,40) AP2(I),BP2(I),CP2(I),AP3(I),BP3(I),CP3(I)         FORT0362
           READ(5,37) AD1(I),BD1(I),CD1(I),MADD(I)                      FORT0363
           READ(5,40) AD2(I),BD2(I),CD2(I),AD3(I),BD3(I),CD3(I)         FORT0364
36         CONTINUE                                                     FORT0365
C                                                                       FORT0366
C     READ IN IOVPOP(I) AND IENRGY(I). INDIVIDULAL OVERLAP POPULATION   FORT0367
C     ANALYSES ARE PERFORMED FROM ORBITAL IOVPOP(N) TO ORBITAL          FORT0368
C     IOVPOP(N+1). INDIVIDUAL ENERGY MATRIX ANALYSES ARE PERFORMED      FORT0369
C     FROM ORBITAL IENRGY(N) TO ORBITAL IENRGY(N+1).                    FORT0370
C                                                                       FORT0371
60         IF(L3) READ(5,67) IOVPOP                                     FORT0372
67         FORMAT(24I3)                                                 FORT0373
           IF(L4) READ(5,67) IENRGY                                     FORT0374
C                                                                       FORT0375
C     PRINT OUT TYPE OF CALCULATION.                                    FORT0376
C                                                                       FORT0377
           IF(METH.EQ.0) WRITE(6,90)                                    FORT0378
           IF(METH.EQ.1) WRITE(6,91)                                    FORT0379
           IF(METH.GE.2) WRITE(6,92)                                    FORT0380
           IF(METH.GT.2) WRITE(6,93)                                    FORT0381
           IF(L5) WRITE(6,94)                                           FORT0382
90         FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION.')               FORT0383
91         FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION WITH CHARGE',    FORT0384
     .     ' ITERATION.',/,T10,'LINEAR CHARGE DEPENDENCE OF SENSE*CHARG'FORT0385
     .     ,'E FOR H(I,I)''S.')                                         FORT0386
92         FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION WITH CHARGE',    FORT0387
     .     ' ITERATION.')                                               FORT0388
93         FORMAT(T10,'MADELUNG CORRECTION INCLUDED.')                  FORT0389
94         FORMAT(T10,'WEIGHTED HIJ FORMULA USED.')                     FORT0390
C                                                                       FORT0391
C     PRINT OUT ATOMIC COORDINATES AND PARAMETERS.                      FORT0392
C                                                                       FORT0393
           IF(PRT(1)) GO TO 80                                          FORT0394
           WRITE(6,74)                                                  FORT0395
74         FORMAT(///,T5,'ATOM',T17,'X',T29,'Y',T41,'Z',T56,'S',T76,'P',FORT0396
     .     T96,'D',T113,'CONTRACTED D'/T47,'N',T50,'EXP',T59,'COUL',    FORT0397
     .     T67,'N',T70,'EXP',T79,'COUL',T87,'N',T90,'EXPD1',T99,'COUL', FORT0398
     .     T109,'C1',T118,'C2',T125,'EXPD2')                            FORT0399
           IF(NH.EQ.0) GO TO 72                                         FORT0400
           J=1                                                          FORT0401
           DO 76 I=1,NH                                                 FORT0402
76         WRITE(6,53) HYDROG,I,X(I),Y(I),Z(I),J,PEEP,COULH             FORT0403
53         FORMAT(T4,A2,I3,3F12.5,3(I3,F8.4,F9.4),2F9.5,F8.4)           FORT0404
C
C  ** THIS CHUNK WRITES ANY HYDROGEN ATOM COORDINATES TO DISK FILE 13
C     THAT ARE DEFINED SPECIFICALLY AS HYDROGENS (I.E., NOT DEFINED AS
C     HEAVY ATOMS). THE HYDROGEN AS A HEAVY ATOM CASE COMES LATER.
C     IT ALSO WRITES THE NUMBER OF ATOMS AND THE NUMBER OF VALENCE
C     ORBITALS (ALSO MOLECULAR ORBITALS) FOR USE IN SETTING UP THE
C     PLOTTING FILES.       JJN   8-28-90
C
           IF (NH.EQ.0) GO TO 72
           WRITE(13,9998) NATOM,NDIM,NELEC,KA,JOHN
9998	   FORMAT(I3,2I4,4I3)
           DO 9991 I=1,NH
9991       WRITE(13,9992) X(I),Y(I),Z(I)
9992       FORMAT(' 1',3F12.6)
C
C  **
C
72         CONTINUE                                                     FORT0405
           DO 151 I=1,NA                                                FORT0406
           KEYI=AC(I)                                                   FORT0407
           INH=I+NH                                                     FORT0408
           IF(NP(KEYI).NE.0) GO TO 152                                  FORT0409
           WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),  FORT0410
     .     EXPS(KEYI),COULS(KEYI)                                       FORT0411
           GO TO 151                                                    FORT0412
152        IF(ND(KEYI).NE.0) GO TO 153                                  FORT0413
           WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),  FORT0414
     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI)       FORT0415
           GO TO 151                                                    FORT0416
153        IF(C2(KEYI).NE.0.0D0) GO TO 154                              FORT0417
           WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),  FORT0418
     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),      FORT0419
     .     ND(KEYI),EXPD(KEYI),COULD(KEYI)                              FORT0420
           GO TO 151                                                    FORT0421
154        WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),  FORT0422
     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),      FORT0423
     .     ND(KEYI),EXPD(KEYI),COULD(KEYI),C1(KEYI),C2(KEYI),EXPD2(KEYI)FORT0424
151        CONTINUE                                                     FORT0425
C
C     THIS CHUNK WRITES THE COORDINATES AND ATOMIC NUMBER 
C     (FOR ELEMENTS LESS THAN OR EQUAL TO 89) TO DISK FILE
C     13 FOR NON-HYDROGEN ATOMS.     JJN  9-8-90
C
C     THE FIRST ITEMS WRITTEN OUT ARE THE NUMBER OF ATOMS AND THE
C     NUMBER OF VALENCE ORBITALS (I.E., MOLECULAR ORBITALS) FOR USE
C     IN CONSTRUCTING THE PLOTTING ROUTINE INPUT FILES.
C
           IF (NH.NE.0) GO TO 10051
10049      WRITE(13,10050) NATOM,NDIM,NELEC,KA,JOHN
10050      FORMAT(I3,2I4,4I3)
10051      DO 9996 I=1,NA
           KEYI=AC(I)
           INH=I+NH
           DO 9995 J=1,89
           IF (SYMBOL(KEYI).EQ.ATS(J)) SYMBL(KEYI)=J
9995       CONTINUE
           WRITE(13,9994) SYMBL(KEYI),X(INH),Y(INH),Z(INH)
9994       FORMAT(I2,3F12.6)
9996       CONTINUE
C
C  
C
           WRITE(6,160) KA,IPRINT,IPUNCH,CON                            FORT0426
160        FORMAT(///,T10,'CHARGE =',I3,8X,'IPRINT =',I3,8X,'IPUNCH =', FORT0427
     .     I3,8X,'HUCKEL CONSTANT =',F7.3)                              FORT0428
C                                                                       FORT0429
C     PRINT OUT ITERATION PARAMETERS.                                   FORT0430
C                                                                       FORT0431
           IF(.NOT.ITERAT) GO TO 80                                     FORT0432
           WRITE(6,81) DAMP1,DAMP2,DAMP3,LAMPRI,MAXCYC,PRTCYC,          FORT0433
     .     SENSE,DELTAC                                                 FORT0434
81         FORMAT(/,T10,'DAMP1 =',F6.3,6X,'DAMP2 =',F6.3,6X,'DAMP3 =',  FORT0435
     .     F6.3,6X,'LAMPRI =',F6.3,//,T10,'MAXCYC =',I3,8X,'PRTCYC =',  FORT0436
     .     I3,8X,'SENSE =',F6.3,6X,'DELTAC =',F10.7)                    FORT0437
C                                                                       FORT0438
C     PRINT OUT VSIE PARAMETERS.                                        FORT0439
C                                                                       FORT0440
           IF(METH.LT.2) GO TO 80                                       FORT0441
           WRITE(6,82)                                                  FORT0442
82         FORMAT(///,' VSIE PARAMETERS',//,T10,'ATOM',T26,'A',T39,'B', FORT0443
     .     T52,'C')                                                     FORT0444
           NUSER2=MXUSR2-NUSER                                         
           DO 83 I=1,NUSER2                                             FORT0446
           J=MXUSR2-I                                             
           IF(.NOT.ITABLE(J)) GO TO 83                                  FORT0448
           WRITE(6,84) SYMBOL(J),AS1(I),BS1(I),CS1(I)                   FORT0449
84         FORMAT(/,T11,A2,4X,3F13.5)                                   FORT0450
           IF(NP(J).EQ.0) GO TO 83                                      FORT0451
           IF(ND(J).NE.0) GO TO 85                                      FORT0452
           WRITE(6,86) AP1(I),BP1(I),CP1(I)                             FORT0453
86         FORMAT(T17,3F13.5)                                           FORT0454
           GO TO 83                                                     FORT0455
85         IF(NCON.EQ.3) GO TO 87                                       FORT0456
           WRITE(6,86) AP1(I),BP1(I),CP1(I),AD1(I),BD1(I),CD1(I)        FORT0457
           GO TO 83                                                     FORT0458
87         WRITE(6,86) AS2(I),BS2(I),CS2(I),AS3(I),BS3(I),CS3(I),AP1(I),FORT0459
     .     BP1(I),CP1(I),AP2(I),BP2(I),CP2(I),AP3(I),BP3(I),CP3(I),     FORT0460
     .     AD1(I),BD1(I),CD1(I),AD2(I),BD2(I),CD2(I),AD3(I),BD3(I),     FORT0461
     .     CD3(I)                                                       FORT0462
83         CONTINUE                                                     FORT0463
80         WRITE(6,99)                                                  FORT0464
99         FORMAT(///)                                                  FORT0465
           RETURN                                                       FORT0466
           END                                                          FORT0467
           SUBROUTINE MOVLAP(H,S,MAD,C,SP,PD,MAXS,MAXP,MAXD,COUL0,SORB, FORT0468
     .     IOCC,HDG,     NDIM, ND1 ,NC,NATOM,NTYPE,NHDG)                FORT0469
C                                                                       FORT0470
C     SUBROUTINE TO CALCULATE INTERATOMIC DISTANCES, OVERLAP            FORT0471
C     INTEGRALS, AND MADELUNG PARAMETERS.                               FORT0472
C                                                                       FORT0473
           IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                  
           PARAMETER (MXUSER=230)
           PARAMETER (MXUSR2=231)
           DIMENSION H(NDIM,NDIM),S(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),
     .     SP(NDIM),PD(NDIM),MAXS(NDIM),MAXP(NDIM),MAXD(NDIM),
     .     COUL0(MXUSER),SORB(NATOM),IOCC(NDIM),HDG(NHDG)              
           REAL*8 MAD                                                   FORT0478
           LOGICAL*4 SP,PD                                              FORT0479
           INTEGER*4 COUL0,SORB           
C                                                                       FORT0482
C     COUL0 DIMENSIONED AT 20 FOR EASY READING DURING PROCESSING        FORT0483
C     OF DELETION INPUT. DELETIONS DONE IN SUBROUTINE DELETS.           FORT0484
C                                                                       FORT0485
	   PARAMETER (MAXATM=500)
	   PARAMETER (BB=250)
           COMMON/TITLE/AB(10)                                          FORT0486
           COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,        FORT0487
     .     IPRINT,IPUNCH,LA,LB,L3,L4,L5,ONEMAT,ITERAT                   FORT0488
           LOGICAL*1 LA,LB,L3,L4,L5,ONEMAT,ITERAT                       FORT0489
           COMMON/OUT/PRT(20),PUN(20)                                   FORT0490
           LOGICAL*1 PRT,PUN                                            FORT0491
           COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),
     .     ND(BB),EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),
     .     C2(BB),COULS(BB),COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),
     .     Z(MAXATM)   
           INTEGER*2 SYMBOL,KEY,VELEC                                   FORT0495
           COMMON/STARS/STAR,STAR2                                      FORT0496
           INTEGER*2 STAR,STAR2                                         FORT0497
           COMMON /LOCLAP/ SK1,SK2,R,L1,L2,M,N1,N2,MAX                  FORT0498
           REAL*4 MADS(MXUSER),MADP(MXUSER),MADD(MXUSER)                 
           EQUIVALENCE (MADS(1),NS(MXUSR2)),(MADP(1),NP(MXUSR2)),       
     .     (MADD(1),ND(MXUSR2))                                         
           DIMENSION PTR(9),DTR(25)                                     FORT0502
           DIMENSION A(MXUSER),B(MXUSER),A1(MXUSER),B1(MXUSER)      
           LOGICAL*1 JGO                                                FORT0504
           EQUIVALENCE (PTR(3),CA),(PTR(8),CB)                          FORT0505
           DATA SQRT3/1.7320508075688770/                               FORT0506
           DATA AUI/1.889644746D0/                                      FORT0507
           DATA PTR(9)/0.D0/,DTR(12)/0.D0/,DTR(22)/0.D0/                FORT0508
           NH1=NH+1                                                     FORT0509
C                                                                       FORT0510
C     HYDROGEN-HYDROGEN OVERLAPS.                                       FORT0511
C                                                                       FORT0512
           IF(NH.LE.1) GO TO 106                                        FORT0513
           DO 107 I=2,NH                                                FORT0514
           IM1=I-1                                                      FORT0515
           DO 107 J=1,IM1                                               FORT0516
           R=DSQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2)        FORT0517
C                                                                       FORT0518
C     STORE OVERLAPS IN UPPER RIGHT TRIANGLE OF S(I,J). PUT             FORT0519
C     DISTANCES IN LOWER RIGHT TRIANGLE.                                FORT0520
C                                                                       FORT0521
           S(I,J)=R                                                     FORT0522
           R=R*PEEP*AUI                                                 FORT0523
           IF(R.GT.50) GO TO 105                                        FORT0524
           SIGMA=(1.D0+R*(1.D0+R/3.D0))/DEXP(R)                         FORT0525
           GO TO 107                                                    FORT0526
105        SIGMA=0.D0                                                   FORT0527
107        S(J,I)=SIGMA                                                 FORT0528
C                                                                       FORT0529
C     HEAVY ATOM-HEAVY ATOM OVERLAPS. LOCAL COORDINATE SYSTEM           FORT0530
C     CENTERED ON ATOM J. FILL IN UPPER RIGHT TRIANGLE OF S(I,J).       FORT0531
C                                                                       FORT0532
106        IORB=NH1                                                     FORT0533
           DO 130 I=1,NA                                                FORT0534
           INH=I+NH                                                     FORT0535
           IM1=I-1                                                      FORT0536
           KEYI=KEY(I)                                                  FORT0537
           MAXD(I)=ND(KEYI)                                             FORT0538
           MAXP(I)=NP(KEYI)                                             FORT0539
           MAXS(I)=NS(KEYI)                                             FORT0540
           SP(I)=EXPS(KEYI) .EQ. EXPP(KEYI)                             FORT0541
           PD(I)=EXPP(KEYI) .EQ. EXPD(KEYI)                             FORT0542
           IF(PD(I)) MAXP(I)=MAX0(MAXP(I),MAXD(I))                      FORT0543
           IF(SP(I)) MAXS(I)=MAX0(NS(KEYI),MAXP(I))                     FORT0544
           SORB(I)=IORB                                                 FORT0545
C                                                                       FORT0546
C     SORB(I) CONTAINS A POINTER TO THE S ORBITAL ON ATOM I.            FORT0547
C                                                                       FORT0548
           IORBS=IORB                                                   FORT0549
           IORB=IORB+4                                                  FORT0550
           IF(NP(KEYI).EQ.0) IORB=IORB-3                                FORT0551
           IF(ND(KEYI).NE.0) IORB=IORB+5                                FORT0552
           IF(NP(KEYI).EQ.0) GO TO 298                                  FORT0553
           JD=IORB-1                                                    FORT0554
           JD1=JD-1                                                     FORT0555
           DO 280 JC = IORBS,JD1                                        FORT0556
           ID=JC+1                                                      FORT0557
           DO 280 IC=ID,JD                                              FORT0558
280        S(JC,IC)=0.D0                                                FORT0559
298        CONTINUE                                                     FORT0560
           IF(I.EQ.1) GO TO 300                                         FORT0561
           DO 131 J=1,IM1                                               FORT0562
           KEYJ=KEY(J)                                                  FORT0563
           JNH=J+NH                                                     FORT0564
           JORBS=SORB(J)                                                FORT0565
           DELX=X(INH)-X(JNH)                                           FORT0566
           DELY=Y(INH)-Y(JNH)                                           FORT0567
           DELZ=Z(INH)-Z(JNH)                                           FORT0568
           RT2=DELX**2+DELY**2                                          FORT0569
           R=DSQRT(RT2+DELZ**2)                                         FORT0570
           S(INH,JNH)=R                                                 FORT0571
C                                                                       FORT0572
C     STORE DISTANCES IN LOWER LEFT TRIANGLE OF S(I,J).                 FORT0573
C                                                                       FORT0574
           IF(R.GT.0.0D0) GO TO 102                                     FORT0575
           ID=IORB-1                                                    FORT0576
           JD=SORB(J+1)-1                                               FORT0577
           DO 103 IC=IORBS,ID                                           FORT0578
           DO 103 JC=JORBS,JD                                           FORT0579
103        S(JC,IC)=0.0D0                                               FORT0580
           GO TO 131                                                    FORT0581
102        IF(RT2.GT.1.E-10) GO TO 135                                  FORT0582
           CB=1.D0                                                      FORT0583
           SB=0.D0                                                      FORT0584
           SA=0.D0                                                      FORT0585
           GOTO 136                                                     FORT0586
135        T=DSQRT(RT2)                                                 FORT0587
           CB=DELX/T                                                    FORT0588
           SB=DELY/T                                                    FORT0589
           SA=T/R                                                       FORT0590
136        CA=DELZ/R                                                    FORT0591
C                                                                       FORT0592
C     THE TRANSFORMATION MATRICES ARE CALCULATED EXPLICITLY.            FORT0593
C     PTR IS THE MATRIX FOR PROJECTING THE X,Y,Z ORBITALS               FORT0594
C     ONTO THE LOCAL SYSTEM. THE ELEMENTS ARE ORDERED SO THAT FIRST     FORT0595
C     X THEN Y THEN Z IS PROJECTED ONTO THE Z' AXIS (SIGMA).            FORT0596
C     THEN THE 3 ARE PROJECTED ONTO THE X' AXIS AND THEN THE Y' (PI).   FORT0597
C     THE D ORBITALS ARE HANDLED SIMILARLY. THE ORDER OF PROJECTION     FORT0598
C     IS X2-Y2,Z2,XY,XZ,YZ FIRST ONTO Z2'(SIGMA)AND THEN ONTO XZ' AND   FORT0599
C     YZ'(PI). FINALLY THE 5 ORBITALS ARE PROJECTED ONTO X'2-Y'2 AND    FORT0600
C     THEN XY' (DELTA).                                                 FORT0601
C                                                                       FORT0602
C     THOSE PTR AND DTR WHICH ARE ZERO ARE INITIALIZED IN A DATA STATE- FORT0603
C     MENT.  CA AND CB HAVE BEEN EQUIVALENCED TO PTR(3) AND PTR(8)      FORT0604
C     RESPECTIVELY TO SAVE TIME.                                        FORT0605
C                                                                       FORT0606
           PTR(1)= SA*CB                                                FORT0607
           PTR(2)= SA*SB                                                FORT0608
C ...      PTR(3)= CA                                                   FORT0609
           PTR(4)= CA*CB                                                FORT0610
           PTR(5)= CA*SB                                                FORT0611
           PTR(6)= -SA                                                  FORT0612
           PTR(7)= -SB                                                  FORT0613
C ...      PTR(8)= CB                                                   FORT0614
C ...      PTR(9)= 0.D0                                                 FORT0615
           IF(ND(KEYI)+ND(KEYJ).EQ.0) GO TO 180                         FORT0616
           CA2=CA**2                                                    FORT0617
           SA2=SA*SA                                                    FORT0618
           CB2=CB*CB                                                    FORT0619
           SB2=SB*SB                                                    FORT0620
           CBSB= CB*SB                                                  FORT0621
           CASA= CA*SA                                                  FORT0622
           CB2SB2= CB2-SB2                                              FORT0623
           DTR(1)= SQRT3*.5D0*SA2*CB2SB2                                FORT0624
           DTR(2)= 1.D0-1.5D0*SA2                                       FORT0625
           DTR(3)= SQRT3*CBSB*SA2                                       FORT0626
           DTR(4)= SQRT3*CASA*CB                                        FORT0627
           DTR(5)= SQRT3*CASA*SB                                        FORT0628
           DTR(6)= CASA*CB2SB2                                          FORT0629
           DTR(7)= -SQRT3*CASA                                          FORT0630
           DTR(8)= 2.D0*CASA*CBSB                                       FORT0631
           DTR(9)= CB*(CA2-SA2)                                         FORT0632
           DTR(10)= SB*(CA2-SA2)                                        FORT0633
           DTR(11)= -2.D0*SA*CBSB                                       FORT0634
C ...      DTR(12)= 0.D0                                                FORT0635
           DTR(13)= SA* CB2SB2                                          FORT0636
           DTR(14)= -PTR(5)                                             FORT0637
           DTR(15)= PTR(4)                                              FORT0638
           IF(ND(KEYI)*ND(KEYJ).EQ.0) GO TO 180                         FORT0639
           DTR(16)=.5D0*(1.D0+CA2)*CB2SB2                               FORT0640
           DTR(17)= .5D0*SQRT3*SA2                                      FORT0641
           DTR(18)= CBSB*(1.D0+CA2)                                     FORT0642
           DTR(19)= -CASA*CB                                            FORT0643
           DTR(20)= -CASA*SB                                            FORT0644
           DTR(21)= -2.D0*CA*CBSB                                       FORT0645
C ...      DTR(22)= 0.D0                                                FORT0646
           DTR(23)= CA*CB2SB2                                           FORT0647
           DTR(24)= PTR(2)                                              FORT0648
           DTR(25)= -PTR(1)                                             FORT0649
180        R=R*AUI                                                      FORT0650
C                                                                       FORT0651
C     (S(I)!S(J)).                                                      FORT0652
C                                                                       FORT0653
           N2=NS(KEYJ)                                                  FORT0654
           N1=NS(KEYI)                                                  FORT0655
           L2=0                                                         FORT0656
           L1=0                                                         FORT0657
           M=0                                                          FORT0658
           MAX=MAXS(I)+MAXS(J)                                          FORT0659
           SK1=EXPS(KEYI)                                               FORT0660
           SK2=EXPS(KEYJ)                                               FORT0661
           CALL ABFNS(A,B)                                              FORT0662
           CALL LOVLAP(SIGMA,A,B)                                       FORT0663
           S(JORBS,IORBS)=SIGMA                                         FORT0664
C                                                                       FORT0665
C     IF THE S EXPONENT OF ATOM I EQUALS THE P EXPONENT WE NEED         FORT0666
C     NOT CALCULATE THE A AND B FUNCTIONS AGAIN.                        FORT0667
C                                                                       FORT0668
C     (P(I)!S(J)).                                                      FORT0669
C                                                                       FORT0670
           JGO=.FALSE.                                                  FORT0671
           IF(KEYI.EQ.KEYJ) GO TO 126                                   FORT0672
           IF((.NOT.SP(I)).OR.(NP(KEYI).EQ.0)) GO TO 126                FORT0673
220        N1=NP(KEYI)                                                  FORT0674
           L1=1                                                         FORT0675
           CALL LOVLAP(SIGMA,A,B)                                       FORT0676
           SIGMA=-SIGMA                                                 FORT0677
           DO 200 IC=1,3                                                FORT0678
200        S(JORBS,IORBS+IC)=PTR(IC)*SIGMA                              FORT0679
           IF(PD(I)) GO TO 221                                          FORT0680
           IF(JGO) GO TO 217                                            FORT0681
           GO TO 137                                                    FORT0682
C                                                                       FORT0683
C     (D(I)!S(J)) CONDITIONALLY AT FIRST CHANCE.                        FORT0684
C                                                                       FORT0685
221        N1=ND(KEYI)                                                  FORT0686
           L1=2                                                         FORT0687
168        CALL LOVLAP(SIGMA,A,B)                                       FORT0688
           IF(C2(KEYI).EQ.0.D0) GO TO 167                               FORT0689
           SK1=EXPD2(KEYI)                                              FORT0690
           CALL ABFNS(A1,B1)                                            FORT0691
           CALL LOVLAP(PART2,A1,B1)                                     FORT0692
           SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2                          FORT0693
           SK1=EXPD(KEYI)                                               FORT0694
167        ID=IORBS+3                                                   FORT0695
           DO 201 IC=1,5                                                FORT0696
201        S(JORBS,ID+IC)=DTR(IC)*SIGMA                                 FORT0697
C                                                                       FORT0698
C     CALCULATE (D(I)!P(J)) IF CAN USE SAME A'S AND B'S.                FORT0699
C                                                                       FORT0700
           IF(SP(J)) GO TO 222                                          FORT0701
           IF(JGO) GO TO 228                                            FORT0702
           GO TO 137                                                    FORT0703
222        N2=NP(KEYJ)                                                  FORT0704
           L2=1                                                         FORT0705
           M=0                                                          FORT0706
           CALL LOVLAP(SIGMA,A,B)                                       FORT0707
           M=1                                                          FORT0708
           CALL LOVLAP(PI,A,B)                                          FORT0709
           IF(C2(KEYI).EQ.0.D0) GO TO 1169                              FORT0710
           SK1=EXPD2(KEYI)                                              FORT0711
           CALL LOVLAP(PART2,A1,B1)                                     FORT0712
           PI=C1(KEYI)*PI+C2(KEYI)*PART2                                FORT0713
           M=0                                                          FORT0714
           CALL LOVLAP(PART2,A1,B1)                                     FORT0715
           SK1=EXPD(KEYI)                                               FORT0716
           SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2                          FORT0717
1169       PI=-PI                                                       FORT0718
           ID=IORBS+3                                                   FORT0719
           DO 195 JC=1,3                                                FORT0720
           DO 195 IC=1,5                                                FORT0721
195        S(JORBS+JC,ID+IC)=PTR(JC)*DTR(IC)*SIGMA+(PTR(JC+3)*DTR(IC+5) FORT0722
     .     +PTR(JC+6)*DTR(IC+10))*PI                                    FORT0723
           IF(JGO) GO TO 131                                            FORT0724
C                                                                       FORT0725
C     NOW TEST FOR DUPLICATE EXPONENTS ON ATOM J.                       FORT0726
C     HOWEVER DO CALCULATIONS ANYHOW.                                   FORT0727
C                                                                       FORT0728
137        N1=NS(KEYI)                                                  FORT0729
           L1=0                                                         FORT0730
C                                                                       FORT0731
C     (S(I)!P(J)).                                                      FORT0732
C                                                                       FORT0733
126        IF(SP(J)) GO TO 138                                          FORT0734
           IF(NP(KEYJ).EQ.0) GO TO 210                                  FORT0735
           MAX=MAXS(I)+MAXP(J)                                          FORT0736
           SK2=EXPP(KEYJ)                                               FORT0737
           CALL ABFNS(A,B)                                              FORT0738
138        N2=NP(KEYJ)                                                  FORT0739
           L2=1                                                         FORT0740
           M=0                                                          FORT0741
           CALL LOVLAP(SIGMA,A,B)                                       FORT0742
           DO 202 IC=1,3                                                FORT0743
202        S(JORBS+IC,IORBS)=PTR(IC)*SIGMA                              FORT0744
           IF(SP(I)) GO TO 156                                          FORT0745
           JGO=.TRUE.                                                   FORT0746
           IF(ND(KEYJ).NE.0) GO TO 149                                  FORT0747
C                                                                       FORT0748
C     BRANCH TO TEST FOR EXPP(J) .EQ. EXPD(J). CALCULATE (S!D) ANYHOW.  FORT0749
C     RETURN WILL BE MADE TO THE NEXT STATEMENT.                        FORT0750
C                                                                       FORT0751
C     (P(I)!P(J))   EXPP(I) EQ,NE EXPS(I).                              FORT0752
C                                                                       FORT0753
           GO TO 646                                                    FORT0754
146        N2=NP(KEYJ)                                                  FORT0755
           L2=1                                                         FORT0756
           SK2=EXPP(KEYJ)                                               FORT0757
646        IF(NP(KEYI).EQ.0) GO TO 210                                  FORT0758
           SK1=EXPP(KEYI)                                               FORT0759
C                                                                       FORT0760
C     THESE STATEMENTS USED ONLY IF HAVE ALREADY CALCULATED (S(I)!D(J)) FORT0761
C     WHICH MEANS THAT SP(I) IS FALSE.                                  FORT0762
C                                                                       FORT0763
           MAX=MAXP(I)+MAXP(J)                                          FORT0764
           CALL ABFNS(A,B)                                              FORT0765
156        N1=NP(KEYI)                                                  FORT0766
           L1=1                                                         FORT0767
148        M=0                                                          FORT0768
           CALL LOVLAP(SIGMA,A,B)                                       FORT0769
           SIGMA=-SIGMA                                                 FORT0770
           M=1                                                          FORT0771
           CALL LOVLAP(PI,A,B)                                          FORT0772
           DO 204 JC=1,3                                                FORT0773
           DO 204 IC=JC,3                                               FORT0774
           S(JORBS+JC,IORBS+IC)=PTR(JC)*PTR(IC)*SIGMA + (PTR(JC+3)*     FORT0775
     .     PTR(IC+3)+PTR(JC+6)*PTR(IC+6))*PI                            FORT0776
204        S(JORBS+IC,IORBS+JC)=S(JORBS+JC,IORBS+IC)                    FORT0777
147        IF(ND(KEYJ).EQ.0) GO TO 210                                  FORT0778
C                                                                       FORT0779
C     BRANCH AROUND (S(I)!D(J)) IF ALREADY DONE.                        FORT0780
C                                                                       FORT0781
           IF(JGO) GO TO 160                                            FORT0782
C                                                                       FORT0783
C     (S(I)!D(J)).                                                      FORT0784
C                                                                       FORT0785
           N1=NS(KEYI)                                                  FORT0786
           L1=0                                                         FORT0787
149        N2=ND(KEYJ)                                                  FORT0788
           L2=2                                                         FORT0789
           IF(PD(J)) GO TO 142                                          FORT0790
           SK2=EXPD(KEYJ)                                               FORT0791
           MAX=MAXS(I)+MAXD(J)                                          FORT0792
           CALL ABFNS(A,B)                                              FORT0793
142        M=0                                                          FORT0794
           CALL LOVLAP(SIGMA,A,B)                                       FORT0795
           IF(C2(KEYJ).EQ.0.D0) GO TO 151                               FORT0796
           SK2=EXPD2(KEYJ)                                              FORT0797
           CALL ABFNS(A1,B1)                                            FORT0798
           CALL LOVLAP(PART2,A1,B1)                                     FORT0799
           SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART2                          FORT0800
           SK2=EXPD(KEYJ)                                               FORT0801
151        JD=JORBS+3                                                   FORT0802
           DO 205 IC=1,5                                                FORT0803
205        S(JD+IC,IORBS)=DTR(IC)*SIGMA                                 FORT0804
150        IF(JGO) GO TO 146                                            FORT0805
C                                                                       FORT0806
C     SP(I) IS TRUE IF HERE SO BRANCH AS WE ALSO HAVE D ON ATOM J.      FORT0807
C                                                                       FORT0808
           GO TO 170                                                    FORT0809
160        JGO=.FALSE.                                                  FORT0810
C                                                                       FORT0811
C          (P(I)!D(J)).                                                 FORT0812
C                                                                       FORT0813
           N2=ND(KEYJ)                                                  FORT0814
           L2=2                                                         FORT0815
           IF(PD(J)) GO TO 178                                          FORT0816
           SK2=EXPD(KEYJ)                                               FORT0817
           MAX=MAXP(I)+MAXD(J)                                          FORT0818
           CALL ABFNS(A,B)                                              FORT0819
178        IF(C2(KEYJ).EQ.0.D0) GO TO 170                               FORT0820
           SK2=EXPD2(KEYJ)                                              FORT0821
           CALL ABFNS(A1,B1)                                            FORT0822
           SK2=EXPD(KEYJ)                                               FORT0823
170        N1=NP(KEYI)                                                  FORT0824
           L1=1                                                         FORT0825
           M=0                                                          FORT0826
           CALL LOVLAP(SIGMA,A,B)                                       FORT0827
           M=1                                                          FORT0828
           CALL LOVLAP(PI,A,B)                                          FORT0829
           IF(C2(KEYJ).EQ.0.D0) GO TO 171                               FORT0830
           SK2=EXPD2(KEYJ)                                              FORT0831
           CALL LOVLAP(PART2,A1,B1)                                     FORT0832
           PI=C1(KEYJ)*PI+C2(KEYJ)*PART2                                FORT0833
           M=0                                                          FORT0834
           CALL LOVLAP(PART2,A1,B1)                                     FORT0835
           SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART2                          FORT0836
           SK2=EXPD(KEYJ)                                               FORT0837
171        SIGMA=-SIGMA                                                 FORT0838
           DO 206 IC=1,3                                                FORT0839
           DO 206 JC=1,5                                                FORT0840
206        S(JD+JC,IORBS+IC)=DTR(JC)*PTR(IC)*SIGMA+(DTR(JC+5)*PTR(IC+3) FORT0841
     .     +DTR(JC+10)*PTR(IC+6))*PI                                    FORT0842
C                                                                       FORT0843
C     (D(I)!D(J)).                                                      FORT0844
C                                                                       FORT0845
           IF(ND(KEYI).EQ.0) GO TO 210                                  FORT0846
           MAX=MAXD(I)+MAXD(J)                                          FORT0847
           IF(PD(I)) GO TO 208                                          FORT0848
           SK1=EXPD(KEYI)                                               FORT0849
           CALL ABFNS(A,B)                                              FORT0850
           IF(C2(KEYJ).EQ.0.D0) GO TO 208                               FORT0851
           SK2=EXPD2(KEYJ)                                              FORT0852
           CALL ABFNS(A1,B1)                                            FORT0853
           SK2=EXPD(KEYJ)                                               FORT0854
208        N1=ND(KEYI)                                                  FORT0855
           L1=2                                                         FORT0856
           M=0                                                          FORT0857
           CALL LOVLAP(SIGMA,A,B)                                       FORT0858
           M=1                                                          FORT0859
           CALL LOVLAP(PI,A,B)                                          FORT0860
           M=2                                                          FORT0861
           CALL LOVLAP(DELTA,A,B)                                       FORT0862
           CC=C2(KEYI)                                                  FORT0863
           IF(C2(KEYJ).EQ.0.D0) GO TO 173                               FORT0864
           CC=C1(KEYJ)*CC                                               FORT0865
           SK2=EXPD2(KEYJ)                                              FORT0866
           CALL LOVLAP(PART2,A1,B1)                                     FORT0867
           DELTA=C1(KEYJ)*DELTA+C2(KEYJ)*PART2                          FORT0868
           M=1                                                          FORT0869
           CALL LOVLAP(PART3,A1,B1)                                     FORT0870
           PI=C1(KEYJ)*PI+C2(KEYJ)*PART3                                FORT0871
           M=0                                                          FORT0872
           CALL LOVLAP(PART4,A1,B1)                                     FORT0873
           SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART4                          FORT0874
           SK2=EXPD(KEYJ)                                               FORT0875
           M=2                                                          FORT0876
173        IF(C2(KEYI).EQ.0.D0) GO TO 172                               FORT0877
           IF(KEYI.EQ.KEYJ) GO TO 176                                   FORT0878
           SK1=EXPD2(KEYI)                                              FORT0879
           CALL ABFNS(A1,B1)                                            FORT0880
           CALL LOVLAP(PART2,A1,B1)                                     FORT0881
           M=1                                                          FORT0882
           CALL LOVLAP(PART3,A1,B1)                                     FORT0883
           M=0                                                                  
           CALL LOVLAP(PART4,A1,B1)                                     FORT0884
176        SIGMA=C1(KEYI)*SIGMA+CC*PART4                                FORT0885
           PI =C1(KEYI)*PI+CC*PART3                                     FORT0886
           DELTA=C1(KEYI)*DELTA+CC*PART2                                FORT0887
           IF(C2(KEYJ).EQ.0.D0) GO TO 172                               FORT0888
           SK1=EXPD2(KEYI)                                              FORT0889
           SK2=EXPD2(KEYJ)                                              FORT0890
           CALL ABFNS(A1,B1)                                            FORT0891
           M=0                                                          FORT0892
           CALL LOVLAP(PART2,A1,B1)                                     FORT0893
           CC=C2(KEYI)*C2(KEYJ)                                         FORT0894
           SIGMA=SIGMA+CC*PART2                                         FORT0895
           M=1                                                          FORT0896
           CALL LOVLAP(PART2,A1,B1)                                     FORT0897
           PI=PI+CC*PART2                                               FORT0898
           M=2                                                          FORT0899
           CALL LOVLAP(PART2,A1,B1)                                     FORT0900
           DELTA=DELTA+CC*PART2                                         FORT0901
172        PI=-PI                                                       FORT0902
           JD=JORBS+3                                                   FORT0903
           DO 211 IC=1,5                                                FORT0904
           ID=IORBS+3                                                   FORT0905
           DO 211 JC=1,5                                                FORT0906
           S(JD+JC,ID+IC) = DTR(IC)*DTR(JC)*SIGMA+(DTR(IC+5)*DTR(JC+5)  FORT0907
     .     +DTR(IC+10)*DTR(JC+10))*PI+(DTR(IC+15)*DTR(JC+15)+DTR(IC+20) FORT0908
     .     *DTR(JC+20))*DELTA                                           FORT0909
211        S(JD+IC,ID+JC)=S(JD+JC,ID+IC)                                FORT0910
C                                                                       FORT0911
C     FILLING IN OTHER HALF OF OVERLAPS FOR (J!I) AS NEEDED.            FORT0912
C                                                                       FORT0913
210        IF(KEYI.EQ.KEYJ) GO TO 213                                   FORT0914
           N2=NS(KEYJ)                                                  FORT0915
           L2=0                                                         FORT0916
           SK2=EXPS(KEYJ)                                               FORT0917
           M=0                                                          FORT0918
           JGO=.TRUE.                                                   FORT0919
           IF(NP(KEYI).EQ.0) GO TO 131                                  FORT0920
           IF(SP(I)) GO TO 215                                          FORT0921
           MAX=MAXP(I)+MAXS(J)                                          FORT0922
           SK1=EXPP(KEYI)                                               FORT0923
           CALL ABFNS(A,B)                                              FORT0924
           GO TO 220                                                    FORT0925
215        IF(PD(I)) GO TO 227                                          FORT0926
217        IF(ND(KEYI).EQ.0) GO TO 131                                  FORT0927
           MAX=MAXD(I)+MAXS(J)                                          FORT0928
           SK1=EXPD(KEYI)                                               FORT0929
           CALL ABFNS(A,B)                                              FORT0930
           GO TO 221                                                    FORT0931
227        IF(SP(J)) GO TO 131                                          FORT0932
           N1=ND(KEYI)                                                  FORT0933
           L1=2                                                         FORT0934
           SK1=EXPD(KEYI)                                               FORT0935
228        IF(NP(KEYJ).EQ.0) GO TO 131                                  FORT0936
           SK2=EXPP(KEYJ)                                               FORT0937
           MAX=MAXD(I)+MAXP(J)                                          FORT0938
           CALL ABFNS(A,B)                                              FORT0939
           IF(C2(KEYI).EQ.0.D0) GO TO 222                               FORT0940
           SK1=EXPD2(KEYI)                                              FORT0941
           CALL ABFNS(A1,B1)                                            FORT0942
           SK1=EXPD(KEYI)                                               FORT0943
           GO TO 222                                                    FORT0944
213        IF(NP(KEYI).EQ.0) GO TO 131                                  FORT0945
           DO 237 IC=1,3                                                FORT0946
237        S(JORBS,IORBS+IC)=-S(JORBS+IC,IORBS)                         FORT0947
           IF(ND(KEYI).EQ.0) GO TO 131                                  FORT0948
           DO 238 IC=4,8                                                FORT0949
           S(JORBS,IORBS+IC)=S(JORBS+IC,IORBS)                          FORT0950
           DO 238 JC=1,3                                                FORT0951
238        S(JORBS+JC,IORBS+IC)=-S(JORBS+IC,IORBS+JC)                   FORT0952
131        CONTINUE                                                     FORT0953
300        IF(NH.EQ.0) GO TO 130                                        FORT0954
           N2=1                                                         FORT0955
           L2=0                                                         FORT0956
           M=0                                                          FORT0957
           SK2=PEEP                                                     FORT0958
           DO 301 J=1,NH                                                FORT0959
           DELX=X(J)-X(INH)                                             FORT0960
           DELY=Y(J)-Y(INH)                                             FORT0961
           DELZ=Z(J)-Z(INH)                                             FORT0962
           RT2=DELX**2+DELY**2                                          FORT0963
           R=DSQRT(RT2+DELZ**2)                                         FORT0964
C                                                                       FORT0965
C     STORE DISTANCES IN LOWER LEFT TRIANGLE OF S(I,J).                 FORT0966
C                                                                       FORT0967
           S(INH,J)=R                                                   FORT0968
           IF(RT2.GT.1.D-10) GO TO 303                                  FORT0969
           CB=1.D0                                                      FORT0970
           SB=0.D0                                                      FORT0971
           SA=0.D0                                                      FORT0972
           GO TO 302                                                    FORT0973
303        T=DSQRT(RT2)                                                 FORT0974
           CB=DELX/T                                                    FORT0975
           SB=DELY/T                                                    FORT0976
           SA=T/R                                                       FORT0977
302        CA=DELZ/R                                                    FORT0978
           R=R*AUI                                                      FORT0979
C                                                                       FORT0980
C     H(J)!S(I)).                                                       FORT0981
C                                                                       FORT0982
           N1=NS(KEYI)                                                  FORT0983
           L1=0                                                         FORT0984
           MAX=1+MAXS(I)                                                FORT0985
           SK1=EXPS(KEYI)                                               FORT0986
           CALL ABFNS(A,B)                                              FORT0987
           CALL LOVLAP(SIGMA,A,B)                                       FORT0988
           S(J,IORBS)=SIGMA                                             FORT0989
           IF(NP(KEYI).EQ.0) GO TO 301                                  FORT0990
           IF(SP(I)) GO TO 304                                          FORT0991
           SK1=EXPP(KEYI)                                               FORT0992
           MAX=1+MAXP(I)                                                FORT0993
           CALL ABFNS(A,B)                                              FORT0994
304        N1=NP(KEYI)                                                  FORT0995
           L1=1                                                         FORT0996
           CALL LOVLAP(SIGMA,A,B)                                       FORT0997
           S(J,IORBS+3)=CA*SIGMA                                        FORT0998
           SIGMA=SIGMA*SA                                               FORT0999
           S(J,IORBS+2)=SB*SIGMA                                        FORT1000
           S(J,IORBS+1)=CB*SIGMA                                        FORT1001
           IF(ND(KEYI).EQ.0) GO TO 301                                  FORT1002
           IF(PD(I)) GO TO 305                                          FORT1003
           SK1=EXPD(KEYI)                                               FORT1004
           MAX=1+ND(KEYI)                                               FORT1005
           CALL ABFNS(A,B)                                              FORT1006
305        N1=ND(KEYI)                                                  FORT1007
           L1=2                                                         FORT1008
           CALL LOVLAP(SIGMA,A,B)                                       FORT1009
           IF(C2(KEYI).EQ.0.D0) GO TO 181                               FORT1010
           SK1=EXPD2(KEYI)                                              FORT1011
           CALL ABFNS(A1,B1)                                            FORT1012
           CALL LOVLAP(PART2,A1,B1)                                     FORT1013
           SK1=EXPD(KEYI)                                               FORT1014
           SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2                          FORT1015
181        CONTINUE                                                     FORT1016
           S(J,IORBS+5)=(1.D0-1.5D0*SA*SA)*SIGMA                        FORT1017
           SIGMA=SIGMA*SQRT3*SA                                         FORT1018
           S(J,IORBS+4)=.5D0*SA*(CB*CB-SB*SB)*SIGMA                     FORT1019
           S(J,IORBS+6)=CB*SB*SA*SIGMA                                  FORT1020
           SIGMA=SIGMA*CA                                               FORT1021
           S(J,IORBS+7)=CB*SIGMA                                        FORT1022
           S(J,IORBS+8)=SB*SIGMA                                        FORT1023
301        CONTINUE                                                     FORT1024
130        CONTINUE                                                     FORT1025
C                                                                       FORT1026
C     CALL DELETS TO SET CERTAIN OVERLAP INTEGRALS = 0.                 FORT1027
C                                                                       FORT1028
           IF(.NOT.LB) GO TO 835                                        FORT1029
           CALL DELETS(S,COUL0,SORB,NDIM)                               FORT1030
           WRITE(6,2010)                                                FORT1031
2010       FORMAT(///)                                                  FORT1032
C                                                                       FORT1033
C     CALCULATE INTERATOMIC MADELUNG PARAMETERS.                        FORT1034
C                                                                       FORT1035
835        IF(METH.LT.3) GO TO 450                                      FORT1036
           IC=1                                                         FORT1037
           DO 401 I=1,NA                                                FORT1038
           KEYI=KEY(I)                                                  FORT1039
           RS=0.0D0                                                     FORT1040
           ID=IC                                                        FORT1041
           MAD(ID,ID)=DBLE(MADS(MXUSR2-KEYI))                         
           IF(NP(KEYI).EQ.0) GO TO 402                                  FORT1043
           ID=IC+1                                                      FORT1044
           MAD(ID,ID)=DBLE(MADP(MXUSR2-KEYI))                           
           IF(ND(KEYI).EQ.0) GO TO 403                                  FORT1046
           ID=IC+2                                                      FORT1047
           MAD(ID,ID)=DBLE(MADD(MXUSR2-KEYI))                           
403        M=IC+1                                                       FORT1049
           DO 404 K=M,ID                                                FORT1050
           K1=K-1                                                       FORT1051
           DO 404 L=IC,K1                                               FORT1052
           CA=MAD(K,K)                                                  FORT1053
           CB=MAD(L,L)                                                  FORT1054
           SA=VALMAD(CA,CB,RS)                                          FORT1055
           MAD(K,L)=SA                                                  FORT1056
404        MAD(L,K)=SA                                                  FORT1057
402        IF(I.EQ.1) GO TO 401                                         FORT1058
           IM1=I-1                                                      FORT1059
           JC=1                                                         FORT1060
           DO 406 J=1,IM1                                               FORT1061
           KEYJ=KEY(J)                                                  FORT1062
           RS=S(I,J)*AUI/27.21D0                                        FORT1063
           JD=JC                                                        FORT1064
           IF(NP(KEYJ).NE.0) JD=JC+1                                    FORT1065
           IF(ND(KEYJ).NE.0) JD=JC+2                                    FORT1066
           DO 407 K=IC,ID                                               FORT1067
           CA=MAD(K,K)                                                  FORT1068
           DO 407 L=JC,JD                                               FORT1069
           CB=MAD(L,L)                                                  FORT1070
           SA=VALMAD(CA,CB,RS)                                          FORT1071
           MAD(K,L)=SA                                                  FORT1072
407        MAD(L,K)=SA                                                  FORT1073
406        JC=JD+1                                                      FORT1074
401        IC=ID+1                                                      FORT1075
C                                                                       FORT1076
C     SET UP DISTANCE MATRIX FOR PRINTING.                              FORT1077
C     STUFF ELEMENTS OF S INTO C TO GET THEM OUT OF THE WAY.            FORT1078
C                                                                       FORT1079
450        ISUB=1                                                       FORT1080
C                                                                       FORT1081
C     ZERO DISTANCE ALONG DIAGONAL.                                     FORT1082
C                                                                       FORT1083
           S(1,1)=0.D0                                                  FORT1084
           DO 1010 I=2,NATOM                                            FORT1085
           S(I,I)=0.D0                                                  FORT1086
           IM1=I-1                                                      FORT1087
           DO 1005 J=1,IM1                                              FORT1088
           C(ISUB)=S(J,I)                                               FORT1089
           ISUB=ISUB+1                                                  FORT1090
1005       S(J,I)=S(I,J)                                                FORT1091
1010       CONTINUE                                                     FORT1092
           IF(PRT(3)) GO TO 2004                                        FORT1093
           WRITE(6,2000)                                                FORT1094
2000       FORMAT('DISTANCE MATRIX')                                    FORT1095
           CALL PEGLEG(S,NATOM,NDIM)                                    FORT1096
2004       IF(PUN(3)) WRITE(7,2050) ((S(I,J),I=1,NATOM),J=1,NATOM)      FORT1097
2050       FORMAT(8F9.6)                                                FORT1098
C                                                                       FORT1099
C     SET UP OVERLAP MATRIX FOR PRINTING.                               FORT1100
C     REPLACE ELEMENTS IN OVERLAP MATRIX FROM C.                        FORT1101
C                                                                       FORT1102
1015       S(1,1)=1.D0                                                  FORT1103
           ISUB=1                                                       FORT1104
           DO 1025 I=2,NDIM                                             FORT1105
           S(I,I)=1.D0                                                  FORT1106
           IM1=I-1                                                      FORT1107
           DO 1020 J=1,IM1                                              FORT1108
           IF(I.GT.NATOM) GO TO 1020                                    FORT1109
           S(J,I)=C(ISUB)                                               FORT1110
           ISUB=ISUB+1                                                  FORT1111
1020       S(I,J)=S(J,I)                                                FORT1112
1025       CONTINUE                                                     FORT1113
           IF(PRT(4)) GO TO 2005                                        FORT1114
           WRITE(6,2001)                                                FORT1115
2001       FORMAT('OVERLAP MATRIX')                                     FORT1116
           CALL PEGLEG(S,NDIM,NDIM)                                     FORT1117
2005       IF(PUN(4)) WRITE(7,2050) S                                   FORT1118
C                                                                       FORT1119
C     PRINT OUT MADELUNG PARAMETERS.                                    FORT1120
C                                                                       FORT1121
           IF(METH.LT.3) GO TO 460                                      FORT1122
           IF(PRT(5)) GO TO 2006                                        FORT1123
           WRITE(6,2002)                                                FORT1124
2002       FORMAT('MADELUNG PARAMETERS')                                FORT1125
           CALL PEGLEG(MAD,NTYPE,NTYPE)                                 FORT1126
2006       IF(PUN(5)) WRITE(7,2050) MAD                                 FORT1127
460        RETURN                                                       FORT1128
           END                                                          FORT1129
      SUBROUTINE DELETS(S,COUL0,SORB,NDIM)                              FORT1130
C                                                                       FORT1131
C  SUBROUTINE FOR SETTING CERTAIN OVERLAP INTEGRALS EQUAL TO ZERO.      FORT1132
C                                                                       FORT1133
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT1134
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      DIMENSION S(NDIM,NDIM),COUL0(MXUSER),SORB(NDIM)                  
      INTEGER*4 COUL0                                                   FORT1136
      INTEGER*2 SORB                                                    FORT1137
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,      FORT1138
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                              FORT1139
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                            FORT1140
      COMMON/OUT/PRT(20),PUN(20)                                        FORT1141
      LOGICAL*1 PRT,PUN                                                 FORT1142
      LOGICAL*1 IERR                                                    FORT1143
      DATA ORBTL/' ORBITAL'/,ATMPR/'  ATOM  '/                          FORT1144
      SORB(NA+1)=NDIM+1                                                 FORT1145
      IERR=.FALSE.                                                      FORT1146
C                                                                       FORT1147
C  READ IN NUMBERS INDICATING WHICH OVERLAP INTEGRALS ARE TO BE SET     FORT1148
C  TO ZERO. A POSITIVE NUMBER REFERS TO AN ATOM, A NEGATIVE ONE TO      FORT1149
C  AN ORBITAL.                                                          FORT1150
C                                                                       FORT1151
10    READ(5,1000,END=300) COUL0                                        FORT1152
1000  FORMAT(20I3)                                                      FORT1153
      DO 200 IDEL=1,19,2                                                FORT1154
      I=COUL0(IDEL)                                                     FORT1155
      J=COUL0(IDEL+1)                                                   FORT1156
C                                                                       FORT1157
C  TERMINATE ON ENCOUNTERING A ZERO.                                    FORT1158
C                                                                       FORT1159
      IF(I.EQ.0.OR.J.EQ.0) GO TO 400                                    FORT1160
      IF(IERR) GO TO 200                                                FORT1161
      IABSV=IABS(I)                                                     FORT1162
      JABSV=IABS(J)                                                     FORT1163
      IF(I.GT.NH) GO TO 20                                              FORT1164
C                                                                       FORT1165
C  I REFERS TO ORBITAL (NEGATIVE) OR H ATOM (POSITIVE,LE NH).           FORT1166
C                                                                       FORT1167
      ILOW=IABSV                                                        FORT1168
      IHIGH=IABSV                                                       FORT1169
C                                                                       FORT1170
C  ERROR IF I OUT OF RANGE OF ORBITAL NUMBERS.                          FORT1171
C                                                                       FORT1172
      IF(IABSV.GT.NDIM) GO TO 160                                       FORT1173
      GO TO 30                                                          FORT1174
C                                                                       FORT1175
C  ERROR IF I OUT OF RANGE OF ATOMS.                                    FORT1176
C                                                                       FORT1177
20    IF(I.GT.NATM) GO TO 160                                           FORT1178
      ILOW=SORB(IABSV-NH)                                               FORT1179
      IHIGH=SORB(IABSV-NH+1)-1                                          FORT1180
30    IF(J.GT.NH) GO TO 40                                              FORT1181
      JLOW=JABSV                                                        FORT1182
      JHIGH=JABSV                                                       FORT1183
C                                                                       FORT1184
C  CHECK TO SEE IF J IS IN RANGE.                                       FORT1185
C                                                                       FORT1186
      IF(JABSV.GT.NDIM) GO TO 160                                       FORT1187
      GO TO 50                                                          FORT1188
40    IF(J.GT.NATM) GO TO 160                                           FORT1189
      JLOW=SORB(JABSV-NH)                                               FORT1190
      JHIGH=SORB(JABSV-NH+1)-1                                          FORT1191
50    X1=ATMPR                                                          FORT1192
      IF(I.LT.0) X1=ORBTL                                               FORT1193
      X2=ATMPR                                                          FORT1194
      IF(J.LT.0) X2=ORBTL                                               FORT1195
      IF(.NOT.PRT(2)) WRITE(6,1002) X1,IABSV,X2,JABSV                   FORT1196
1002  FORMAT('0ALL S(I,J) SET TO ZERO BETWEEN',A8,I4,' AND',A8,I4,'.')  FORT1197
C                                                                       FORT1198
C  J MUST BE LESS THAN OR EQUAL TO I SINCE WE ONLY HAVE A HALF          FORT1199
C  MATRIX AT THIS POINT.                                                FORT1200
C                                                                       FORT1201
      IF(JHIGH.LE.IHIGH) GO TO 60                                       FORT1202
      I=JHIGH                                                           FORT1203
      JHIGH=IHIGH                                                       FORT1204
      IHIGH=I                                                           FORT1205
      I=JLOW                                                            FORT1206
      JLOW=ILOW                                                         FORT1207
      ILOW=I                                                            FORT1208
60    DO 100 I=ILOW,IHIGH                                               FORT1209
      DO 100 J=JLOW,JHIGH                                               FORT1210
100   S(J,I)=0.D0                                                       FORT1211
      GO TO 200                                                         FORT1212
160   IERR=.TRUE.                                                       FORT1213
      WRITE(6,1003) COUL0                                               FORT1214
1003  FORMAT('0NUMBER OUT OF RANGE IN FOLLOWING DELETION CARD'/'0',20I5/
     *'0NO FURTHER DELETIONS, BUT SCANNING FOR ZERO TO TERMINATE CARD RE
     *ADING.')                                                          FORT1217
200   CONTINUE                                                          FORT1218
      GO TO 10                                                          FORT1219
300   WRITE(6,1004)                                                     FORT1220
1004  FORMAT('0END OF FILE IN DELETION CARDS, NO FURTHER DELETIONS.')   FORT1221
400   RETURN                                                            FORT1222
      END                                                               FORT1223
           DOUBLE PRECISION FUNCTION VALMAD(A,B,R)                      FORT1224
C                                                                       FORT1225
C     FUNCTION ROUTINE FOR CALCULATING MADELUNG PARAMETERS.             FORT1226
C                                                                       FORT1227
           IMPLICIT REAL*8(A-H,O-Z)                                     FORT1228
           IF(A.LT.0.01D0.OR.B.LT.0.01D0) GO TO 1                       FORT1229
           AB=(A+B)/(2.0D0*A*B)                                         FORT1230
           VALMAD=1.0D0/DSQRT(R*R+AB*AB)                                FORT1231
           GO TO 2                                                      FORT1232
1          VALMAD=0.0D0                                                 FORT1233
           IF(R.GT.0.001D0) VALMAD=1.0D0/R                              FORT1234
2          RETURN                                                       FORT1235
           END                                                          FORT1236
      SUBROUTINE PEGLEG(A,N,NL)                                         FORT1237
C                                                                       FORT1238
C  SUBROUTINE TO PRINT OUT MATRICES IN READABLE FORMAT.                 FORT1239
C                                                                       FORT1240
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT1241
      DIMENSION A(NL,NL)                                                FORT1242
      NROW=N                                                            FORT1243
      NCOL=N                                                            FORT1244
      GO TO 10                                                          FORT1245
      ENTRY OUTMAT(A,NL,NR,NC)                                          FORT1246
      NROW=NR                                                           FORT1247
      NCOL=NC                                                           FORT1248
10    KITE=0                                                            FORT1249
20    LOW=KITE+1                                                        FORT1250
      KITE=KITE+14                                                      FORT1251
      IF(KITE.GT.NCOL) KITE=NCOL                                        FORT1252
      WRITE(6,1000) (I,I=LOW,KITE)                                      FORT1253
1000  FORMAT(/5X,14I8,//)                                               FORT1254
      DO 30 I=1,NROW                                                    FORT1255
30    WRITE(6,1001) I,(A(I,J),J=LOW,KITE)                               FORT1256
1001  FORMAT(I5,2X,14F8.4)                                              FORT1257
      IF(KITE.LT.NCOL) GO TO 20                                         FORT1258
      RETURN                                                            FORT1259
      END                                                               FORT1260
      SUBROUTINE PEGLEG2(AA,NN,NLL)
C
C     ROUTINE FOR WRITING OUT THE EIGENVECTORS TO DISK FILE 13 FOR
C     PLOTTING ROUTINES.   JJN  8-8-90
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AA(NLL,NLL)
      NROWW=NN
      NCOLL=NN
      WRITE(13,9998) ((AA(J,I),J=1,NROWW), I=1,NCOLL)
9998  FORMAT(8F10.6)
      RETURN
      END
C
C  THIS CONCLUDES THE EIGENVECTOR WRITE ROUTINE TO DISK FILE 13.
C
      SUBROUTINE HUCKEL(H,S,MAD,C,SP,PD,MAXS,MAXP,MAXD,COUL0,SORB,IOCC, FORT1261
     1 HDG,     NDIM, ND1 ,NC,NATOM,NTYPE,NHDG)                         FORT1262
C                                                                       FORT1263
C  SUBROUTINE TO 1) DETERMINE ORBITAL OCCUPATION NUMBERS, AND 2) SETUP  FORT1264
C  HUCKEL MATRIX.                                                       FORT1265
C                                                                       FORT1266
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               FORT1267
      DIMENSION H(NDIM,NDIM),S(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),       FORT1268
     1 SP(NDIM),PD(NDIM),MAXS(NDIM),MAXP(NDIM),MAXD(NDIM),COUL0(NATOM), FORT1269
     2 SORB(NDIM),IOCC(NDIM),HDG(NHDG)                                  FORT1270
      REAL*8 MAD                                                        FORT1271
      REAL*4 IOCC                                                       FORT1272
      PARAMETER (MAXATM=500)
      PARAMETER (BB=250)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      COMMON/TITLE/AB(10)                                               FORT1273
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,      FORT1274
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                              FORT1275
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                            FORT1276
      COMMON/OUT/PRT(20),PUN(20)                                        FORT1277
      LOGICAL*1 PRT,PUN                                                 FORT1278
      COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),ND(BB),
     1 EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),COULS(BB), 
     2 COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)                 
      INTEGER*2 SYMBOL,KEY,VELEC                                        FORT1282
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)                             FORT1284
      REAL*8 LAMPRI                                                     FORT1285
      INTEGER*4 PRTCYC                                                  FORT1286
      LOGICAL*1 PARTIT,PRINTX,ITABLE                                    FORT1287
      COMMON/STARS/STAR,STAR2                                           FORT1288
      INTEGER*2 STAR,STAR2                                              FORT1289
      INTEGER*4 ONE,TWO,STAR1,TRUE                                      FORT1290
      DATA ONE,TWO,STAR1/'1','2','*'/                                   FORT1291
      IF(ONEMAT.AND.IPRINT.LT.-2) RETURN                                FORT1292
      ICYCLE=1                                                          FORT1293
C                                                                       FORT1294
C  SETUP DEFAULT ORBITAL OCCUPATIONS.                                   FORT1295
C                                                                       FORT1296
      IF(L1) GO TO 3                                                    FORT1297
      IH=NELEC/2                                                        FORT1298
      JH=NDIM+1-IH                                                      FORT1299
      DO 1 I=1,JH                                                       FORT1300
1     IOCC(I)=0.0                                                       FORT1301
      DO 2 I=JH,NDIM                                                    FORT1302
2     IOCC(I)=2.0                                                       FORT1303
      IF(IH+IH.NE.NELEC) IOCC(JH-1)=1.0                                 FORT1304
      GO TO 500                                                         FORT1305
C                                                                       FORT1306
C  PROVISION FOR READING IN USER SPECIFIED OCCUPATIONS.                 FORT1307
C  ALSO PROVISION FOR NON-INTEGER OCCUPATIONS.                          FORT1308
C                                                                       FORT1309
3     READ(5,2000) (MAXD(I),I=1,NDIM)                                   FORT1310
2000  FORMAT(80A1)                                                      FORT1311
      TEL=0.0D0                                                         FORT1312
      DO 4 I=1,NDIM                                                     FORT1313
      J=NDIM+1-I                                                        FORT1314
      IOCC(J)=0.0                                                       FORT1315
      IF(MAXD(I).EQ.ONE) IOCC(J)=1.0                                    FORT1316
      IF(MAXD(I).EQ.TWO) IOCC(J)=2.0                                    FORT1317
      IF(MAXD(I).EQ.STAR1) READ(5,900) IOCC(J)                          FORT1318
900   FORMAT(F15.8)                                                     FORT1319
4     TEL=TEL+DBLE(IOCC(J))                                             FORT1320
      TEL=DABS(TEL-DFLOAT(NELEC))                                       FORT1321
      IF(TEL.LT.0.001D0) GO TO 500                                      FORT1322
      WRITE(6,2001)                                                     FORT1323
2001  FORMAT('0*** WARNING ****    ORBITAL POPULATIONS INCONSISTENT WITHFORT1324
     1 ASSUMED CHARGE ON MOLECULE',////,T10,'I',T25,'IOCC(I)',/)        FORT1325
      DO 99 I=1,NDIM                                                    FORT1326
99    WRITE(6,2002) I,IOCC(I)                                           FORT1327
2002  FORMAT(T8,I3,T22,F12.8)                                           FORT1328
      STOP                                                              FORT1329
C                                                                       FORT1330
C  CALL GRMSCH TO ORTHOGONALISE BASIS SET.                              FORT1331
C                                                                       FORT1332
500   CALL GRMSCH(S,C,NDIM)                                             FORT1333
      CON=.5D0*CON                                                      FORT1334
      IF(.NOT.ITERAT) GO TO 15                                          FORT1335
      DO 5 I=1,NDIM                                                     FORT1336
5     SORB(I)=0.D0                                                      FORT1337
      DO 10 I=1,NATOM                                                   FORT1338
      X(I)=0.D0                                                         FORT1339
      Z(I)=0.D0                                                         FORT1340
10    COUL0(I)=0.D0                                                     FORT1341
15    IF(.NOT.ONEMAT) GO TO 25                                          FORT1342
C                                                                       FORT1343
C  IN ONE MATRIX CASE, STUFF DIAGONAL ELEMENTS OF S INTO SP.            FORT1344
C                                                                       FORT1345
      DO 20 I=1,NDIM                                                    FORT1346
20    SP(I)=S(I,I)                                                      FORT1347
C                                                                       FORT1348
C  SETUP DIAGONAL ELEMENTS OF HUCKEL MATRIX IN H(I,J).                  FORT1349
C                                                                       FORT1350
25    IF(NH.EQ.0) GO TO 35                                              FORT1351
      ET=COULH                                                          FORT1352
      DO 30 I=1,NH                                                      FORT1353
      IF(ITERAT) ET=COULH-COUL0(I)                                      FORT1354
30    C(I)=ET                                                           FORT1355
35    IH=NH+1                                                           FORT1356
      ET=0.D0                                                           FORT1357
      DO 50 I=1,NA                                                      FORT1358
      KEYI=KEY(I)                                                       FORT1359
      IF(ITERAT) ET=COUL0(I+NH)                                         FORT1360
      C(IH)=COULS(KEYI)-ET                                              FORT1361
      IH=IH+1                                                           FORT1362
      IF(NP(KEYI).EQ.0) GO TO 50                                        FORT1363
      HH=COULP(KEYI)-ET                                                 FORT1364
      JH=IH+2                                                           FORT1365
      ASSIGN 40 TO IL                                                   FORT1366
      GO TO 42                                                          FORT1367
40    IF(ND(KEYI).EQ.0) GO TO 50                                        FORT1368
      HH=COULD(KEYI)-ET                                                 FORT1369
      JH=IH+4                                                           FORT1370
      ASSIGN 50 TO IL                                                   FORT1371
42    DO 45 J=IH,JH                                                     FORT1372
45    C(J)=HH                                                           FORT1373
      IH=JH+1                                                           FORT1374
      GO TO IL,(40,50)                                                  FORT1375
50    CONTINUE                                                          FORT1376
      DO 55 I=1,NDIM                                                    FORT1377
55    H(I,I)=C(I)                                                       FORT1378
      IF(NHDG.EQ.1) GO TO 59                                            FORT1379
      DO 56 I=1,NDIM                                                    FORT1380
56    HDG(I)=C(I)                                                       FORT1381
C                                                                       FORT1382
C  SETUP OFF-DIAGONAL ELEMENTS OF HUCKEL MATRIX.                        FORT1383
C                                                                       FORT1384
59    CNST=CON                                                          FORT1385
      DO 58 I=2,NDIM                                                    FORT1386
      IL=I-1                                                            FORT1387
      DO 58 J=1,IL                                                      FORT1388
      HH=C(I)+C(J)                                                      FORT1389
      IF(.NOT.L5) GO TO 58                                              FORT1390
      ET=(C(I)-C(J))/HH                                                 FORT1391
      ET=ET*ET                                                          FORT1392
      CNST=CON+ET/2.0D0+ET*ET*(0.5D0-CON)                               FORT1393
58    H(I,J)=CNST*HH*S(I,J)                                             FORT1394
      IF(ONEMAT) GO TO 100                                              FORT1395
      DO 60 I=2,NDIM                                                    FORT1396
      IL=I-1                                                            FORT1397
      DO 60 J=1,IL                                                      FORT1398
60    H(J,I)=H(I,J)                                                     FORT1399
C                                                                       FORT1400
C  PRINT OUT HUCKEL MATRIX. PRINT OUT TITLE IF METH IS NOT              FORT1401
C  EQUAL TO ZERO.                                                       FORT1402
C                                                                       FORT1403
806   IF(ICYCLE.GT.1) GO TO 800                                         FORT1404
      IF(PRT(6)) GO TO 805                                              FORT1405
      IF(METH.EQ.0) GO TO 801                                           FORT1406
      GO TO 802                                                         FORT1407
800   IF(ICYCLE.GE.10000) WRITE(6,701) AB                               FORT1408
701   FORMAT('RESULTS OF CALCULATION  ',8A8,A6,A2,//)                   FORT1409
      IF(PRT(7).OR..NOT.PRINTX) GO TO 805                               FORT1410
801   WRITE(6,803)                                                      FORT1411
803   FORMAT('HUCKEL MATRIX')                                           FORT1412
      CALL PEGLEG(H,NDIM,NDIM)                                          FORT1413
      GO TO 805                                                         FORT1414
802   WRITE(6,804)                                                      FORT1415
804   FORMAT('INPUT HUCKEL MATRIX')                                     FORT1416
      CALL PEGLEG(H,NDIM,NDIM)                                          FORT1417
805   IF(ICYCLE.EQ.1.AND.PUN(6)) WRITE(7,825) H                         FORT1418
      IF(ICYCLE.GT.1.AND.PUN(7).AND.PRINTX) WRITE(7,825) H              FORT1419
825   FORMAT(8F9.5)                                                     FORT1420
C                                                                       FORT1421
C  IF CALCULATING ENERGY MATRIX, STORE H(I,I) IN X(I),Y(I),Z(I).        FORT1422
C                                                                       FORT1423
      IF(ICYCLE.LE.MAXCYC.AND.METH.NE.0) GO TO 100                      FORT1424
      IF(NH.EQ.0) GO TO 369                                             FORT1425
      DO 370 I=1,NH                                                     FORT1426
370   X(I)=H(I,I)                                                       FORT1427
369   IH=NH+1                                                           FORT1428
      JH=NH+1                                                           FORT1429
      DO 371 I=1,NA                                                     FORT1430
      KEYI=KEY(I)                                                       FORT1431
      X(IH)=H(JH,JH)                                                    FORT1432
      JH=JH+1                                                           FORT1433
      IF(NP(KEYI).EQ.0) GO TO 371                                       FORT1434
      Y(IH)=H(JH,JH)                                                    FORT1435
      JH=JH+3                                                           FORT1436
      IF(ND(KEYI).EQ.0) GO TO 371                                       FORT1437
      Z(IH)=H(JH,JH)                                                    FORT1438
      JH=JH+5                                                           FORT1439
371   IH=IH+1                                                           FORT1440
C                                                                       FORT1441
C  CALL TRNFRM TO TRANSFORM HUCKEL MATRIX TO ORTHOGONAL BASIS SET.      FORT1442
C  THEN CALL GIVENS TO PERFORM DIAGONALIZATION.                         FORT1443
C                                                                       FORT1444
100   IH=1                                                              FORT1445
      IF(ONEMAT) IH=2                                                   FORT1446
      CALL TRNFRM(S,H,C,COUL0,NDIM,SP,IH)                               FORT1447
      IF(ONEMAT) GO TO 110                                              FORT1448
      IH=NDIM                                                           FORT1449
      GO TO 120                                                         FORT1450
110   IH=-NDIM                                                          FORT1451
120   CALL GIVENS(NDIM,IH,NDIM,C,SP,COUL0,H)                            FORT1452
130   IF(ICYCLE.GE.10000) ITERAT=.FALSE.                                FORT1453
C                                                                       FORT1454
C  PRINT OUT TITLE, ENERGY LEVELS, AND OCCUPATION NUMBERS.              FORT1455
C                                                                       FORT1456
      IF(ITERAT) GO TO 700                                              FORT1457
      IF(METH.EQ.0) WRITE(6,701) AB                                     FORT1458
      IF(PRT(8)) GO TO 710                                              FORT1459
      WRITE(6,702) (I,COUL0(I),IOCC(I),I=1,NDIM)                        FORT1460
702   FORMAT(////,57X,'ENERGY LEVELS (EV)',/,(/52X,'E(',I3,') =',F12.5, FORT1461
     18X,F6.4))                                                         FORT1462
710   IF(PUN(8)) WRITE(7,825) (COUL0(I),I=1,NDIM)                       FORT1463
700   IF(ONEMAT) GO TO 200                                              FORT1464
C                                                                       FORT1465
C  DIDDLE WITH C,H.                                                     FORT1466
C                                                                       FORT1467
      DO 160 J=1,NDIM                                                   FORT1468
      DO 140 K=1,NDIM                                                   FORT1469
140   C(K)=H(K,IH)                                                      FORT1470
      DO 155 I=1,NDIM                                                   FORT1471
      ET=0.D0                                                           FORT1472
      DO 150 K=I,NDIM                                                   FORT1473
150   ET=ET+S(I,K)*C(K)                                                 FORT1474
155   H(I,IH)=ET                                                        FORT1475
160   IH=IH-1                                                           FORT1476
      K=1                                                               FORT1477
      DO 180 I=2,NDIM                                                   FORT1478
      IL=I-1                                                            FORT1479
      DO 180 J=1,IL                                                     FORT1480
      C(K)=S(I,J)                                                       FORT1481
180   K=K+1                                                             FORT1482
200   IF(METH.GT.1.AND.ITERAT) GO TO 210                                FORT1483
C                                                                       FORT1484
C  CALL OUTPUT FOR FINAL PRINT OUT OF RESULTS.                          FORT1485
C                                                                       FORT1486
205   CALL OUTPUT(H,S,MAD,C,COUL0,SORB,IOCC,HDG,NDIM,NTYPE,NC,NHDG)     FORT1487
      IF(.NOT.ITERAT) GO TO 999                                         FORT1488
      GO TO 220                                                         FORT1489
C                                                                       FORT1490
C  IF DOING CHARGE ITERATIVE CALCULATION ( METH >1 ), CALL ITRATE       FORT1491
C  TO SETUP HUCKEL MATRIX.                                              FORT1492
C                                                                       FORT1493
210   CALL ITRATE(H,S,MAD,C,COUL0,SORB,IOCC,HDG,NDIM,NTYPE,NC,NHDG)     FORT1494
      IF(.NOT.ITERAT) GO TO 205                                         FORT1495
220   IF(ICYCLE.GT.MAXCYC) ICYCLE=10000                                 FORT1496
      IF(METH.GT.1) GO TO 806                                           FORT1497
      GO TO 15                                                          FORT1498
  999 RETURN                                                            FORT1499
      END                                                               FORT1500
      SUBROUTINE GIVENS (NX,NROOTX,NJX,A,B,ROOT,VECT)                   FORT1501
C                                                                       FORT1502
C      SUBROUTINE TO CALCULATE THE EIGENVALUES AND EIGENVECTORS         FORT1503
C      OF A REAL SYMMETRIC MATRIX.                                      FORT1504
C                                                                       FORT1505
C                                                                       FORT1506
C      THE PARAMETERS FOR THE ROUTINE ARE0                              FORT1507
C                                                                       FORT1508
C          NX     ORDER OF MATRIX.                                      FORT1509
C                                                                       FORT1510
C          NROOTX NUMBER OF ROOTS FOR WHICH EIGENVECTORS ARE WANTED.    FORT1511
C                 IF NO VECTORS ARE WANTED, MAKE NROOTX NEGATIVE.       FORT1512
C                                                                       FORT1513
C          NJX    ROW DIMENSION OF VECT ARRAY.  SEE 'VECT' BELOW.       FORT1514
C                 NJX MUST BE NOT LESS THAN NX.                         FORT1515
C                                                                       FORT1516
C          A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR   FORT1517
C                 FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE          FORT1518
C                 LOCATIONS.                                            FORT1519
C                                                                       FORT1520
C          B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST       FORT1521
C                 NX*6 CELLS.                                           FORT1522
C                                                                       FORT1523
C          ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST      FORT1524
C                 NX CELLS LONG.  THE ROOTS ARE ORDERED LARGEST FIRST   FORT1525
C                 IN THIS ARRAY.                                        FORT1526
C                                                                       FORT1527
C          VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN          FORT1528
C                 EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE      FORT1529
C                 DIMENSIONED WITH 'NJX' ROWS AND AT LEAST 'NJX'        FORT1530
C                 COLUMNS, UNLESS NO VECTORS ARE REQUESTED (NEGATIVE    FORT1531
C                 NROOTX).  IN THIS LATTER CASE, THE ARGUMENT VECT      FORT1532
C                 IS JUST A DUMMY, AND THE STORAGE IS NOT USED.         FORT1533
C                 THE EIGENVECTORS ARE NORMALIZED TO UNIT LENGTH.       FORT1534
C                                                                       FORT1535
C      THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE        FORT1536
C      RESULTS APPEAR IN ROOT AND VECT.                                 FORT1537
C                                                                       FORT1538
C      FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING FORT1539
C      POINT UNDERFLOW SHOULD BE A ZERO.                                FORT1540
C                                                                       FORT1541
C      THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE   FORT1542
C      REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.            FORT1543
C                                                                       FORT1544
C      THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS0  FORT1545
C                                                                       FORT1546
C      FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE    FORT1547
C      HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).  FORT1548
C      THE EIGENVALUES OF THE TRIDIAGONAL MATRIX ARE THEN FOUND USING   FORT1549
C      THE QR TRANSFORM METHOD.  SEE J. H. WILKINSON, THE ALGEBRAIC     FORT1550
C      EIGENVALUE PROBLEM(1965) FOR A DESCRIPTION OF THIS ALGORITHM.    FORT1551
C      THE EIGENVECTORS OF THE TRIDIAGONAL FORM ARE THEN EVALUATED      FORT1552
C      (J. H. WILKINSON, COMP. J. 1, 90 (1958)), BY THE METHOD OF       FORT1553
C      INVERSE ITERATION, FOR NONDEGENERATE MATRICES.                   FORT1554
C      FOR MATRICES WITH DEGENERATE OR NEAR-DEGENERATE EIGENVALUES,     FORT1555
C      THE EIGENVECTORS ARE EVALUATED INSTEAD BY FURTHER QR TRANSFORMS. FORT1556
C      THIS METHOD GIVES ORTHOGONAL VECTORS EVEN FOR DEGENERATE ROOTS.  FORT1557
C      FINALLY THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE    FORT1558
C      ORIGINAL ARRAY (FIRST REFERENCE).                                FORT1559
C                                                                       FORT1560
C      THE INVERSE ITERATION PORTION OF THIS PROGRAM WAS ADAPTED        FORT1561
C      FROM THE QUANTUM CHEMISTRY PROGRAM EXCHANGE NUMBER 62.1, BY      FORT1562
C      FRANKLIN PROSSER.  THE EIGENVALUE SUBROUTINE (EVQR) WAS WRITTEN  FORT1563
C      BY WALTER NIELSEN.                                               FORT1564
C                                      ROY GORDON, SEPT. 1969           FORT1565
C                                                                       FORT1566
C      AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN    FORT1567
C      J. M. ORTEGA'S ARTICLE IN 'MATHEMATICS FOR DIGITAL COMPUTERS,'   FORT1568
C      VOLUMD 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.        FORT1569
C                                                                       FORT1570
      IMPLICIT DOUBLE PRECISION(A-H,R-Z)                                FORT1571
      COMMON /VECTOR/ FACT,IDIF                                         FORT1572
      DIMENSION B(NX,6),A(1),ROOT(NX),VECT(NJX,NROOTX)                  FORT1573
C                                                                       FORT1574
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  FORT1575
C  *                                                                 *  FORT1576
C  *   USERS PLEASE NOTE0 TWO PARAMETERS, ETA AND THETA, SHOULD BE   *  FORT1577
C  *   ADJUSTED BY THE USER FOR HIS PARTICULAR MACHINE.              *  FORT1578
C  *                                                                 *  FORT1579
C  *   ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT   *  FORT1580
C  *   REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),  *  FORT1581
C  *   WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).             *  FORT1582
C  *                                                                 *  FORT1583
C  *   THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE    *  FORT1584
C  *   EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE   *  FORT1585
C  *   LARGEST NUMBER).                                              *  FORT1586
C  *                                                                 *  FORT1587
C  *   SOME RECOMMENDED VALUES FOLLOW.                               *  FORT1588
C  *                                                                 *  FORT1589
C  *   FOR IBM 7094, UNIVAC 1108, ETC. (27-BIT BINARY FRACTION,      *  FORT1590
C  *   8-BIT BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.               *  FORT1591
C  *   FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY  *  FORT1592
C  *   EXPONENT), ETA=1.E-11, THETA=1.E307.                          *  FORT1593
C  *   FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY  *  FORT1594
C  *   EXPONENT), ETA=1.E-14, THETA=1.E307.                          *  FORT1595
C  *   FOR IBM 360/50 AND 360/65 DOUBLE PRECISION (56-BIT HEXA-      *  FORT1596
C  *   DECIMAL FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16,    *  FORT1597
C  *   THETA=1.E75.                                                  *  FORT1598
C  *                                                                 *  FORT1599
C  *   OTHER PARAMETERS WHICH MUST BE ADJUSTED ARE0                  *  FORT1600
C  *                                                                 *  FORT1601
C  *   DEL1 = ETA/1.D2, DELTA = ETA**2*1.D2, SMALL = ETA**2/1.D2,    *  FORT1602
C  *   DELBIG = THETA*DELTA/1.D3, THETA1 = 1.D3/THETA, EMAG = ETA,   *  FORT1603
C  *   TOLER = 1.D2*DSQRT(ETA)                                       *  FORT1604
C  *                                                                 *  FORT1605
C  *   TOLER IS A FACTOR USED TO DETERMINE IF ANY ROOTS ARE CLOSE    *  FORT1606
C  *   ENOUGH TOGETHER TO BE CONSIDERED DEGENERATE FOR PURPOSES OF   *  FORT1607
C  *   CALCULATING EIGENVECTORS.  FOR THE MATRIX NORMED TO UNITY, IF *  FORT1608
C  *   THE DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN THE *  FORT1609
C  *   QR TRANSFORMATION IS USED TO FORM THE EIGENVECTORS.           *  FORT1610
C  *                                                                 *  FORT1611
C  *   EMAG IS A TOLERANCE FOR NEGLIGIBLE ELEMENTS IN THE QR         *  FORT1612
C  *   ITERATION FOR EIGENVECTORS FOR DEGENERATE EIGENVALUES.        *  FORT1613
C  *                                                                 *  FORT1614
C  *   IN THE FOLLOWING ROUTINE, ETA = 1.D-16 AND THETA = 1.D75.     *  FORT1615
C  *                                                                 *  FORT1616
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  FORT1617
C                                                                       FORT1618
      DATA ETA/1.D-16/,THETA/1.D38/,DEL1/1.D-18/,DELTA/1.D-30/,         FORT1619
     *SMALL/1.D-34/,DELBIG/1.D05/,THETA1/1.D-38/,TOLER/1.D-6/,          FORT1620
     *EMAG/1.D-16/                                                      FORT1621
      N = NX                                                            FORT1622
      FLOATN = DFLOAT(N)                                                FORT1623
      NROOT = IABS(NROOTX)                                              FORT1624
      IF (NROOT.EQ.0) GO TO 1001                                        FORT1625
      IF (N-1) 1001,1003,105                                            FORT1626
 1003 ROOT(1) = A(1)                                                    FORT1627
      IF (NROOTX.GT.0) VECT(1,1) = 1.D0                                 FORT1628
      GO TO 1001                                                        FORT1629
  105 CONTINUE                                                          FORT1630
C                                                                       FORT1631
C  NSIZE IS THE NUMBER OF ELEMENTS IN THE PACKED ARRAY.                 FORT1632
C                                                                       FORT1633
      NSIZE = (N*(N+1))/2                                               FORT1634
      NM1 = N-1                                                         FORT1635
      NM2 = N-2                                                         FORT1636
      NP1 = N+1                                                         FORT1637
C                                                                       FORT1638
C  COMPUTE TRACE.                                                       FORT1639
C                                                                       FORT1640
      TRACE = 0.D0                                                      FORT1641
      JUMP = 1                                                          FORT1642
      DO 1 J=2,NP1                                                      FORT1643
      TRACE = TRACE + A(JUMP)                                           FORT1644
    1 JUMP = JUMP + J                                                   FORT1645
      TRACE = TRACE/FLOATN                                              FORT1646
C                                                                       FORT1647
C  SUBTRACT TRACE FROM DIAGONAL ELEMENTS TO GIVE A MORE RELIABLE NORM   FORT1648
C  WHEN THERE ARE LARGE DIAGONAL ELEMENTS.                              FORT1649
C                                                                       FORT1650
      JUMP = 1                                                          FORT1651
      DO 2 J=2,NP1                                                      FORT1652
      A(JUMP) = A(JUMP) - TRACE                                         FORT1653
    2 JUMP = JUMP + J                                                   FORT1654
C                                                                       FORT1655
C  SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.         FORT1656
C                                                                       FORT1657
      FACTOR = 0.D0                                                     FORT1658
      DO 70 I=1,NSIZE                                                   FORT1659
   70 FACTOR = DMAX1(FACTOR,DABS(A(I)))                                 FORT1660
      IF (FACTOR.NE.0.D0) GO TO 72                                      FORT1661
C                                                                       FORT1662
C  NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.                   FORT1663
C                                                                       FORT1664
      DO 78 I=1,NROOT                                                   FORT1665
      IF (NROOTX.LT.0) GO TO 78                                         FORT1666
      DO 77 J=1,N                                                       FORT1667
   77 VECT(J,I) = 0.D0                                                  FORT1668
      VECT(I,I) = 1.D0                                                  FORT1669
   78 ROOT(I) = 0.D0                                                    FORT1670
      GO TO 1001                                                        FORT1671
   72 ANORM = 0.D0                                                      FORT1672
   86 SCALE = 1.D0/FACTOR                                               FORT1673
      DO 80 I=1,NSIZE                                                   FORT1674
   80 ANORM = ANORM + (A(I)*SCALE)**2                                   FORT1675
      ANORM = ANORM+ANORM                                               FORT1676
C                                                                       FORT1677
C  SUBTRACT DIAGONAL CONTRIBUTIONS WHICH WERE COUNTED TWICE.            FORT1678
C                                                                       FORT1679
      JUMP = 1                                                          FORT1680
      DO 81 J=2,NP1                                                     FORT1681
      ANORM = ANORM -(A(JUMP)*SCALE)**2                                 FORT1682
   81 JUMP = JUMP + J                                                   FORT1683
   83 ANORM = FACTOR*DSQRT(ANORM)                                       FORT1684
      SCALE = 1.D0/ANORM                                                FORT1685
      DO 91 I=1,NSIZE                                                   FORT1686
   91 A(I) = A(I)*SCALE                                                 FORT1687
      ALIMIT = 1.D0                                                     FORT1688
C                                                                       FORT1689
C  TRIDIAGONALIZATION OF SYMMETRIC MATRIX.                              FORT1690
C                                                                       FORT1691
      ID = 0                                                            FORT1692
      IA = 1                                                            FORT1693
      IF (NM2.EQ.0) GO TO 201                                           FORT1694
      DO 200 J=1,NM2                                                    FORT1695
C                                                                       FORT1696
C  J COUNTS ROW OF A MATRIX TO BE DIAGONALIZED. IA INDICATES START OF   FORT1697
C  NON-CODIAGONAL ELEMENTS IN THE ROW. ID IS THE INDEX OF CODIAGONAL    FORT1698
C  ELEMENT ON THE ROW BEING CODIAGONALIZED.                             FORT1699
C                                                                       FORT1700
      IA = IA+J+2                                                       FORT1701
      ID = ID+J+1                                                       FORT1702
      JP2 = J + 2                                                       FORT1703
      J1 = J + 1                                                        FORT1704
C                                                                       FORT1705
C  FIND LIMITS FOR BAND OF SIGNIFICANT MATRIX ELEMENTS.                 FORT1706
C                                                                       FORT1707
      LIMIT = J1                                                        FORT1708
      II = IA                                                           FORT1709
      DO 99 I=JP2,N                                                     FORT1710
      B(I,5) = A(II)                                                    FORT1711
      IF (DABS(B(I,5)).GT.DEL1) LIMIT = I                               FORT1712
   99 II = II + I                                                       FORT1713
      DTEMP = A(ID)                                                     FORT1714
      IF (LIMIT.GT.J1) GO TO 110                                        FORT1715
C                                                                       FORT1716
C  NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL                FORT1717
C  ELEMENTS ARE TINY.                                                   FORT1718
C                                                                       FORT1719
  120 B(J,1) = DTEMP                                                    FORT1720
      A(ID) = 0.D0                                                      FORT1721
      GO TO 200                                                         FORT1722
C                                                                       FORT1723
C  SUM SQUARES OF SIGNIFICANT NON-CODIAGONAL ELEMENTS OF ROW J.         FORT1724
C                                                                       FORT1725
  110 IDIF = LIMIT -JP2                                                 FORT1726
      SUM = DOT(B(JP2,5),B(JP2,5))                                      FORT1727
C                                                                       FORT1728
C  NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES.                        FORT1729
C                                                                       FORT1730
      SUM = DSQRT(SUM + DTEMP**2)                                       FORT1731
C                                                                       FORT1732
C  NEW CODIAGONAL ELEMENT.                                              FORT1733
C                                                                       FORT1734
      B(J,1) = -DSIGN(SUM,DTEMP)                                        FORT1735
C                                                                       FORT1736
C  FIRST NON-ZERO ELEMENT OF THIS W-VECTOR.                             FORT1737
C                                                                       FORT1738
      B(J+1,2) = DSQRT((1.D0 + DABS(DTEMP)/SUM)*5.D-1)                  FORT1739
C                                                                       FORT1740
C  FORM REST OF THE W-VECTOR ELEMENTS.                                  FORT1741
C                                                                       FORT1742
      TEMP = DSIGN(5.D-1/(B(J+1,2)*SUM),DTEMP)                          FORT1743
      II = IA                                                           FORT1744
      DO 130 I=JP2,LIMIT                                                FORT1745
      B(I,2) = A(II)*TEMP                                               FORT1746
  130 II = II + I                                                       FORT1747
C                                                                       FORT1748
C  FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.             FORT1749
C  SCALAR = W-VECTOR*P-VECTOR.                                          FORT1750
C                                                                       FORT1751
      DAK = 0.D0                                                        FORT1752
C                                                                       FORT1753
C  IC IS THE LOCATION OF THE NEXT DIAGONAL ELEMENT. I RUNS OVER THE     FORT1754
C  NON-ZERO P-ELEMENTS. CASES FOR I LESS THAN LIMIT.                    FORT1755
C                                                                       FORT1756
      IC = ID + 1                                                       FORT1757
      LIMLES = LIMIT - 1                                                FORT1758
      DO 188 I=J1,LIMLES                                                FORT1759
C                                                                       FORT1760
C  FORM FIRST PART OF P ELEMENT THEN MOVE IC TO TOP OF NEXT             FORT1761
C  A-MATRIX 'ROW'.                                                      FORT1762
C                                                                       FORT1763
      IDIF = I - J1                                                     FORT1764
      DTEMP = DOT(B(J1,2),A(IC))                                        FORT1765
      IC = IC + I                                                       FORT1766
C                                                                       FORT1767
C  COMPLETE P ELEMENT. CHANGE INCREMENTING MODE AT DIAGONAL ELEMENT.    FORT1768
C                                                                       FORT1769
  178 IP1 = I + 1                                                       FORT1770
      JJ = IC + IDIF                                                    FORT1771
      DTEMP = DTEMP + DSUM(B(N,1),A(JJ),IP1,LIMIT)                      FORT1772
C                                                                       FORT1773
C  BUILD UP THE K-SCALAR (AK).                                          FORT1774
C                                                                       FORT1775
      DAK = DAK + DTEMP*B(I,2)                                          FORT1776
  188 B(I,1) = DTEMP                                                    FORT1777
C                                                                       FORT1778
C  CASE FOR I = LIMIT.                                                  FORT1779
C                                                                       FORT1780
      IDIF = LIMIT - J1                                                 FORT1781
      DTEMP = DOT(B(J1,2),A(IC))                                        FORT1782
      DAK = DAK + DTEMP*B(LIMIT,2)                                      FORT1783
      B(LIMIT,1) = DTEMP                                                FORT1784
      IDIF = LIMIT - J1                                                 FORT1785
C                                                                       FORT1786
C  TEST TO SEE IF ANY I VALUES REMAIN. DO REMAINING VALUES.             FORT1787
C                                                                       FORT1788
      IF (LIMIT.EQ.N) GO TO 190                                         FORT1789
      IC = IC + LIMIT                                                   FORT1790
      LIMLO = LIMIT + 1                                                 FORT1791
      DO 189 I=LIMLO,N                                                  FORT1792
      B(I,1) = DOT(B(J1,2),A(IC))                                       FORT1793
      B(I,2) = 0.D0                                                     FORT1794
  189 IC = IC + I                                                       FORT1795
C                                                                       FORT1796
C  FORM THE Q-VECTOR.                                                   FORT1797
C                                                                       FORT1798
  190 FACT = -DAK                                                       FORT1799
      CALL VECSUM(B(J1,1),B(J1,2))                                      FORT1800
C                                                                       FORT1801
C  TRANSFORM THE REST OF THE A-MATRIX. JJ INDICATES START-1 OF THE      FORT1802
C  REST OF THE A-MATRIX. MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS  FORT1803
C  TO SAVE SPACE. I RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR. FORT1804
C                                                                       FORT1805
      JJ = ID                                                           FORT1806
      DO 160 I=J1,N                                                     FORT1807
      A(JJ) = B(I,2)                                                    FORT1808
      IF (I.GT.LIMIT) GO TO 161                                         FORT1809
      B2 = B(I,2)                                                       FORT1810
      FACT = -B2 - B2                                                   FORT1811
      IDIF = I - J1                                                     FORT1812
      CALL VECSUM(A(JJ+1),B(J1,1))                                      FORT1813
  161 B1 = B(I,1)                                                       FORT1814
      FACT = -B1 - B1                                                   FORT1815
      IDIF = MIN0(I,LIMIT) - J1                                         FORT1816
      CALL VECSUM(A(JJ+1),B(J1,2))                                      FORT1817
  160 JJ = JJ + I                                                       FORT1818
C                                                                       FORT1819
C  STORE AWAY LIMIT FOR LATER USE IN BACK TRANSFORMATION. MOVE LAST     FORT1820
C  CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE.                        FORT1821
C                                                                       FORT1822
  200 B(J,6) = LIMIT                                                    FORT1823
  201 CONTINUE                                                          FORT1824
      B(NM1,1) = A(NSIZE-1)                                             FORT1825
      A(NSIZE-1) = 0.D0                                                 FORT1826
C                                                                       FORT1827
C  USE QR TRANSFORM METHOD TO FIND EIGENVALUES OF THE TRIDIAGONAL       FORT1828
C  MATRIX. MOVE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX INTO        FORT1829
C  ROOT ARRAY. THIS IS A MORE CONVENIENT INDEXING POSITION. ALSO,       FORT1830
C  PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.               FORT1831
C                                                                       FORT1832
      JUMP = 0                                                          FORT1833
      DO 320 J=1,NM1                                                    FORT1834
      JUMP = JUMP + J                                                   FORT1835
      ROOT(J) = A(JUMP)                                                 FORT1836
  320 B(J,3) = B(J,1)**2                                                FORT1837
      ROOT(N) = A(NSIZE)                                                FORT1838
      CALL EVQR(ROOT,B(1,3),N,30,SMALL)                                 FORT1839
C                                                                       FORT1840
C  ROOT NOW CONTAINS THE SHIFTED AND SCALED EIGENVALUES. STORE          FORT1841
C  EIGENVALUES FOR POSSIBLE LATER USE AS SHIFTS IN EVALUATING           FORT1842
C  EIGENVECTORS FOR DEGENERATE MATRICES.                                FORT1843
C                                                                       FORT1844
      DO 325 J=1,N                                                      FORT1845
  325 B(J,2) = ROOT(J)                                                  FORT1846
C                                                                       FORT1847
C  SORT THE EIGENVALUES INTO DESCENDING ALGEBRAIC ORDER.                FORT1848
C                                                                       FORT1849
      DO 330 I=1,NM1                                                    FORT1850
      IP1 = I + 1                                                       FORT1851
      DO 330 J=IP1,N                                                    FORT1852
      IF (ROOT(I).GE.ROOT(J)) GO TO 330                                 FORT1853
      TEMP = ROOT(I)                                                    FORT1854
      ROOT(I) = ROOT(J)                                                 FORT1855
      ROOT(J) = TEMP                                                    FORT1856
  330 CONTINUE                                                          FORT1857
C                                                                       FORT1858
C  QUIT NOW IF NO VECTORS WERE REQUESTED. OTHERWISE, TEST FOR           FORT1859
C  DEGENERACY OR NEAR DEGENERACY OF EIGENVALUES FOR WHICH EIGENVECTORS  FORT1860
C  WERE REQUESTED. IF ONLY ONE VECTOR REQUESTED, DEGENERACY DOESN'T     FORT1861
C  MATTER.                                                              FORT1862
C                                                                       FORT1863
      IF (NROOTX.LT.0) GO TO 1002                                       FORT1864
      IF (NROOTX.EQ.1) GO TO 807                                        FORT1865
      NTOP = NROOT - 1                                                  FORT1866
      DO 400 I=1,NTOP                                                   FORT1867
      IF (DABS(ROOT(I+1)-ROOT(I)).LE.TOLER) GO TO 410                   FORT1868
  400 CONTINUE                                                          FORT1869
C                                                                       FORT1870
C  NEXT STATEMENT IS REACHED IF ALL EIGENVALUES FOR WHICH EIGENVECTORS  FORT1871
C  WERE REQUESTED ARE WELL SEPARATED.                                   FORT1872
C                                                                       FORT1873
      GO TO 807                                                         FORT1874
C                                                                       FORT1875
C  THE FOLLOWING IS REACHED IF THERE ARE ANY DEGENERATE CLUSTERS        FORT1876
C  OF EIGENVALUES.  USE FURTHER QR TRANSFORMS TO EVALUATE               FORT1877
C  THE EIGENVECTORS OF THE TRIDIAGONAL MATRIX.  THIS METHOD             FORT1878
C  GIVES ORTHOGONAL EIGENVECTORS EVEN WHEN THE EIGENVALUES ARE          FORT1879
C  DEGENERATE.  HOWEVER, IT TAKES MORE ARITHMETIC THAN THE METHOD       FORT1880
C  OF INVERSE ITERATION (AT LEAST AT LARGE N).                          FORT1881
C                                                                       FORT1882
C  PUT DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX INTO ROOT.  PUT OFF-     FORT1883
C  DIAGONAL ELEMENTS INTO B(I,3).                                       FORT1884
C                                                                       FORT1885
  410 JUMP = 0                                                          FORT1886
      DO 440 J=1,NM1                                                    FORT1887
      JUMP = JUMP + J                                                   FORT1888
      ROOT(J) = A(JUMP)                                                 FORT1889
  440 B(J,3) = B(J,1)                                                   FORT1890
C                                                                       FORT1891
C  LAST DIAGONAL ELEMENT.                                               FORT1892
C                                                                       FORT1893
      ROOT(N) = A(NSIZE)                                                FORT1894
C                                                                       FORT1895
C  INITIALIZE VECTORS TO A UNIT MATRIX.                                 FORT1896
C                                                                       FORT1897
      DO 450 I=1,N                                                      FORT1898
      DO 445 J=1,N                                                      FORT1899
  445 VECT(J,I) = 0.D0                                                  FORT1900
  450 VECT(I,I) = 1.D0                                                  FORT1901
C                                                                       FORT1902
C  FORM EIGENVECTORS OF TRIDIAGONAL MATRIX FOR DEGENERATE               FORT1903
C  MATRICES AND TRANSPOSE THE VECTORS.                                  FORT1904
C                                                                       FORT1905
      CALL QRTN(ROOT,B(1,3),VECT,B(1,2),N,25,EMAG,NJX)                  FORT1906
      DO 456 I=1,NM1                                                    FORT1907
      IP1 = I + 1                                                       FORT1908
      DO 455 J=IP1,N                                                    FORT1909
      FLIP = VECT(I,J)                                                  FORT1910
      VECT(I,J) = VECT(J,I)                                             FORT1911
  455 VECT(J,I) = FLIP                                                  FORT1912
  456 CONTINUE                                                          FORT1913
C                                                                       FORT1914
C  IF ROOTS WERE NOT LOCATED IN DESCENDING ORDER, INTERCHANGE ROOTS     FORT1915
C  AND VECTORS.                                                         FORT1916
C                                                                       FORT1917
      ITOP = NROOT-1                                                    FORT1918
      DO 480 I=1,ITOP                                                   FORT1919
      IP1 = I+1                                                         FORT1920
      DO 480 J=IP1,N                                                    FORT1921
      IF (ROOT(I).GE.ROOT(J)) GO TO 480                                 FORT1922
      TEMP = ROOT(I)                                                    FORT1923
      ROOT(I) = ROOT(J)                                                 FORT1924
      ROOT(J) = TEMP                                                    FORT1925
      DO 470 K=1,N                                                      FORT1926
      TEMP = VECT(K,I)                                                  FORT1927
      VECT(K,I) = VECT(K,J)                                             FORT1928
  470 VECT(K,J) = TEMP                                                  FORT1929
  480 CONTINUE                                                          FORT1930
C                                                                       FORT1931
C  DEGENERATE VECTORS ARE NOW COMPLETE.                                 FORT1932
C                                                                       FORT1933
      GO TO 940                                                         FORT1934
C                                                                       FORT1935
C  EIGENVECTORS OF TRIDIAGONAL MATRIX FOR NONDEGENERATE MATRICES.       FORT1936
C  INITIALIZE VECTOR ARRAY.                                             FORT1937
C                                                                       FORT1938
  807 CONTINUE                                                          FORT1939
      DO 705 I=1,NROOT                                                  FORT1940
      DO 15 J=1,N                                                       FORT1941
   15 VECT(J,I) = 1.D0                                                  FORT1942
  705 CONTINUE                                                          FORT1943
      DO 700 I=1,NROOT                                                  FORT1944
C                                                                       FORT1945
C  USE INVERSE ITERATION TO FIND VECTORS.                               FORT1946
C                                                                       FORT1947
  701 AROOT = ROOT(I)                                                   FORT1948
      ELIM1 = A(1) - AROOT                                              FORT1949
      ELIM2 = B(1,1)                                                    FORT1950
      JUMP = 1                                                          FORT1951
      DO 750 J=1,NM1                                                    FORT1952
      JUMP = JUMP + J + 1                                               FORT1953
C                                                                       FORT1954
C  GET THE CORRECT PIVOT EQUATION FOR THIS STEP.                        FORT1955
C                                                                       FORT1956
      IF (DABS(ELIM1).LE.DABS(B(J,1))) GO TO 760                        FORT1957
C                                                                       FORT1958
C  FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.              FORT1959
C                                                                       FORT1960
      B(J,2) = ELIM1                                                    FORT1961
      B(J,3) = ELIM2                                                    FORT1962
      B(J,4) = 0.D0                                                     FORT1963
      TEMP = B(J,1)/ELIM1                                               FORT1964
      ELIM1 = A(JUMP) - AROOT - TEMP*ELIM2                              FORT1965
      ELIM2 = B(J+1,1)                                                  FORT1966
      GO TO 755                                                         FORT1967
C                                                                       FORT1968
C  SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.                     FORT1969
C                                                                       FORT1970
  760 B(J,2) = B(J,1)                                                   FORT1971
      B(J,3) = A(JUMP) - AROOT                                          FORT1972
      B(J,4) = B(J+1,1)                                                 FORT1973
      TEMP = 1.D0                                                       FORT1974
      IF (DABS(B(J,1)).GT.THETA1) TEMP = ELIM1/B(J,1)                   FORT1975
      ELIM1 = ELIM2 - TEMP*B(J,3)                                       FORT1976
      ELIM2 = -TEMP*B(J+1,1)                                            FORT1977
C                                                                       FORT1978
C  SAVE FACTOR FOR THE SECOND ITERATION.                                FORT1979
C                                                                       FORT1980
  755 B(J,5) = TEMP                                                     FORT1981
  750 CONTINUE                                                          FORT1982
      B(N,2) = ELIM1                                                    FORT1983
      B(N,3) = 0.D0                                                     FORT1984
      B(N,4) = 0.D0                                                     FORT1985
      B(NM1,4) = 0.D0                                                   FORT1986
      ITER = 1                                                          FORT1987
C                                                                       FORT1988
C  BACK SUBSTITUTE TO GET THIS VECTOR.                                  FORT1989
C                                                                       FORT1990
  790 L = N + 1                                                         FORT1991
      DO 780 J=1,N                                                      FORT1992
      L = L - 1                                                         FORT1993
  786 CONTINUE                                                          FORT1994
      ELIM1 = VECT(L,I)-VECT(L+1,I)*B(L,3)-VECT(L+2,I)*B(L,4)           FORT1995
C                                                                       FORT1996
C  IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN. THIS APPROACH     FORT1997
C  IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-DEPENDENT CALLS TO     FORT1998
C  OVERFLOW ROUTINES.                                                   FORT1999
C                                                                       FORT2000
      IF (DABS(ELIM1).GT.DELBIG) GO TO 782                              FORT2001
      TEMP = B(L,2)                                                     FORT2002
      IF (DABS(B(L,2)).LT.DELTA) TEMP = DELTA                           FORT2003
      VECT(L,I) = ELIM1/TEMP                                            FORT2004
      GO TO 780                                                         FORT2005
  782 DO 784 K=1,N                                                      FORT2006
  784 VECT(K,I) = VECT(K,I)/DELBIG                                      FORT2007
      GO TO 786                                                         FORT2008
  780 CONTINUE                                                          FORT2009
      GO TO (820,900), ITER                                             FORT2010
C                                                                       FORT2011
C  SECOND ITERATION.                                                    FORT2012
C                                                                       FORT2013
  820 ITER = ITER + 1                                                   FORT2014
  890 ELIM1 = VECT(1,I)                                                 FORT2015
      DO 830 J=1,NM1                                                    FORT2016
      IF (B(J,2).EQ.B(J,1)) GO TO 840                                   FORT2017
C                                                                       FORT2018
C  CASE ONE.                                                            FORT2019
C                                                                       FORT2020
  850 VECT(J,I) = ELIM1                                                 FORT2021
      ELIM1 = VECT(J+1,I) - ELIM1*B(J,5)                                FORT2022
      GO TO 830                                                         FORT2023
C                                                                       FORT2024
C  CASE TWO.                                                            FORT2025
C                                                                       FORT2026
  840 VECT(J,I) = VECT(J+1,I)                                           FORT2027
      ELIM1 = ELIM1 - VECT(J+1,I)*TEMP                                  FORT2028
  830 CONTINUE                                                          FORT2029
      VECT(N,I) = ELIM1                                                 FORT2030
      GO TO 790                                                         FORT2031
C                                                                       FORT2032
C  NORMALIZE THE VECTOR.                                                FORT2033
C                                                                       FORT2034
  900 ELIM1 = 0.D0                                                      FORT2035
      DO 904 J=1,N                                                      FORT2036
  904 ELIM1 = DMAX1(DABS(VECT(J,I)),ELIM1)                              FORT2037
      TEMP = 0.D0                                                       FORT2038
      DO 910 J=1,N                                                      FORT2039
      ELIM2 = VECT(J,I)/ELIM1                                           FORT2040
  910 TEMP = TEMP + ELIM2**2                                            FORT2041
      TEMP = 1.D0/(DSQRT(TEMP)*ELIM1)                                   FORT2042
      DO 920 J=1,N                                                      FORT2043
      VECT(J,I) = VECT(J,I)*TEMP                                        FORT2044
      IF (DABS(VECT(J,I)).LT.DEL1) VECT(J,I) = 0.D0                     FORT2045
  920 CONTINUE                                                          FORT2046
  700 CONTINUE                                                          FORT2047
C                                                                       FORT2048
C  ROTATE THE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY.        FORT2049
C  LOOP OVER ALL THE TRANSFORMATION VECTORS.                            FORT2050
C                                                                       FORT2051
  940 IF (NM2.EQ.0) GO TO 1002                                          FORT2052
      JUMP = NSIZE - NP1                                                FORT2053
      IM = NM1                                                          FORT2054
      DO 950 I=1,NM2                                                    FORT2055
      LIMIT = IDINT(B(IM-1,6))                                          FORT2056
      J1 = JUMP                                                         FORT2057
C                                                                       FORT2058
C  MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.      FORT2059
C                                                                       FORT2060
      DO 955 J=IM,LIMIT                                                 FORT2061
      B(J,2) = A(J1)                                                    FORT2062
  955 J1 = J1 + J                                                       FORT2063
      IDIF = LIMIT - IM                                                 FORT2064
C                                                                       FORT2065
C  MODIFY ALL REQUESTED VECTORS.                                        FORT2066
C                                                                       FORT2067
      DO 960 K=1,NROOT                                                  FORT2068
C                                                                       FORT2069
C  FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR.       FORT2070
C                                                                       FORT2071
      TMP = DOT(B(IM,2),VECT(IM,K))                                     FORT2072
      FACT = -TMP - TMP                                                 FORT2073
      CALL VECSUM(VECT(IM,K),B(IM,2))                                   FORT2074
  960 CONTINUE                                                          FORT2075
      JUMP = JUMP - IM                                                  FORT2076
  950 IM = IM - 1                                                       FORT2077
 1002 CONTINUE                                                          FORT2078
C                                                                       FORT2079
C  RESTORE ROOTS TO THEIR PROPER SIZE AND ADD BACK TRACE.               FORT2080
C                                                                       FORT2081
      DO 95 I=1,N                                                       FORT2082
   95 ROOT(I) = ROOT(I)*ANORM + TRACE                                   FORT2083
 1001 RETURN                                                            FORT2084
      END                                                               FORT2085
      SUBROUTINE EVQR(A,B,N,M,TOL)                                      FORT2086
C                                                                       FORT2087
C  SUBROUTINE FOR PERFORMING A QR TRANSFORM ON A REAL SYMMETRIC         FORT2088
C  TRIDIAGONAL MATRIX.                                                  FORT2089
C                                                                       FORT2090
C  ON ENTERING, A CONTAINS THE N DIAGONAL ELEMENTS OF THE TRIDIAGONAL   FORT2091
C  MATRIX. B CONTAINS THE SQUARES OF THE N-1 OFF-DIAGONAL ELEMENTS.     FORT2092
C  ITERATION IS CONTINUED UNTIL THE SQUARES OF THE OFF-DIAGONAL         FORT2093
C  ELEMENTS ARE LESS THAN TOL. TYPICALLY LESS THAN TWO ITERATIONS PER   FORT2094
C  EIGENVALUE ARE REQUIRED. THUS THE UPPER LIMIT M TO THE NUMBER OF     FORT2095
C  ITERATIONS PER EIGENVALUE MAY BE SAFELY SET AT 20 OR SO.             FORT2096
C                                                                       FORT2097
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT2098
      DIMENSION A(1),B(1)                                               FORT2099
      NX = N                                                            FORT2100
      NI = 1                                                            FORT2101
      SH = 0.D0                                                         FORT2102
C                                                                       FORT2103
C  K COUNTS THE NUMBER OF ITERATIONS PER EIGENVALUE.                    FORT2104
C                                                                       FORT2105
      K = 0                                                             FORT2106
      IF (NX-2) 50,60,85                                                FORT2107
C                                                                       FORT2108
C  EACH NEW ITERATION BEGINS HERE.                                      FORT2109
C                                                                       FORT2110
  100 K = K + 1                                                         FORT2111
      IF (K.LT.M) GO TO 101                                             FORT2112
      WRITE (6,1000) K                                                  FORT2113
 1000 FORMAT(35H NO CONVERGENCE OF QR ALGORITHM IN   ,I4,11H ITERATIONS)FORT2114
      CALL EXIT                                                         FORT2115
C                                                                       FORT2116
C  SOLVE THE TWO BY TWO IN THE LOWER RIGHT CORNER AND USE SMALLER ROOT  FORT2117
C  AS THE NEXT SHIFT.                                                   FORT2118
C                                                                       FORT2119
  101 AT = A(NX) + A(NX-1)                                              FORT2120
      ST = AT*5.D-1                                                     FORT2121
      DISC = AT**2 - 4.D0*(A(NX)*A(NX-1)-B(NX-1))                       FORT2122
      IF (DISC.LE.0.D0) GO TO 15                                        FORT2123
      ST = ST - DSIGN(DSQRT(DISC),ST)*5.D-1                             FORT2124
C                                                                       FORT2125
C  INCREASE THE TOTAL SHIFT BY THE TEMPORARY SHIFT.                     FORT2126
C                                                                       FORT2127
   15 SH = SH + ST                                                      FORT2128
C                                                                       FORT2129
C  THIS LOOP SUBTRACTS THE TEMPORARY SHIFT FROM THE DIAGONAL ELEMENTS.  FORT2130
C                                                                       FORT2131
      DO 20 I=1,NX                                                      FORT2132
   20 A(I) = A(I) - ST                                                  FORT2133
C                                                                       FORT2134
C  INITIALIZE.                                                          FORT2135
C                                                                       FORT2136
      G = A(NI)                                                         FORT2137
      PS = G**2                                                         FORT2138
      RS = PS + B(NI)                                                   FORT2139
      SX = B(NI)/RS                                                     FORT2140
      CXS = 1.D0                                                        FORT2141
      CX = PS/RS                                                        FORT2142
      U = SX*(G+A(NI+1))                                                FORT2143
      A(NI) = G + U                                                     FORT2144
      NTOP = NX - 2                                                     FORT2145
C                                                                       FORT2146
C  THIS LOOP COMPLETES ONE ITERATION, THAT IS ONE QR TRANSFORM.         FORT2147
C                                                                       FORT2148
      DO 10 I=NI,NTOP                                                   FORT2149
C                                                                       FORT2150
C  G IS THE GAMMA IN THE NOTATION OF WILKINSON.                         FORT2151
C                                                                       FORT2152
      G = A(I+1) - U                                                    FORT2153
      IF (CX.GT.TOL) GO TO 12                                           FORT2154
      PS = B(I)*CXS                                                     FORT2155
      GO TO 16                                                          FORT2156
   12 PS = G**2/CX                                                      FORT2157
   16 RS = PS + B(I+1)                                                  FORT2158
C                                                                       FORT2159
C  ROTATE AN OFF-DIAGONAL ELEMENT.                                      FORT2160
C                                                                       FORT2161
      B(I) = SX*RS                                                      FORT2162
      SX = B(I+1)/RS                                                    FORT2163
      CXS = CX                                                          FORT2164
      CX = PS/RS                                                        FORT2165
      U = SX*(G+A(I+2))                                                 FORT2166
C                                                                       FORT2167
C  ROTATE A DIAGONAL ELEMENT.                                           FORT2168
C                                                                       FORT2169
      A(I+1) = G + U                                                    FORT2170
   10 CONTINUE                                                          FORT2171
C                                                                       FORT2172
C  COMPUTE THE LAST DIAGONAL ELEMENT.                                   FORT2173
C                                                                       FORT2174
      A(NX) = A(NX) - U                                                 FORT2175
C                                                                       FORT2176
C  COMPUTE THE LAST OFF-DIAGONAL ELEMENT.                               FORT2177
C                                                                       FORT2178
      IF (CX.GT.TOL) GO TO 112                                          FORT2179
      PS = B(NTOP+1)*CXS                                                FORT2180
      GO TO 116                                                         FORT2181
  112 PS = ((A(NX))**2)/CX                                              FORT2182
  116 B(NTOP+1) = SX*PS                                                 FORT2183
C                                                                       FORT2184
C  END OF ONE ITERATION.                                                FORT2185
C                                                                       FORT2186
   85 IT = NX                                                           FORT2187
C                                                                       FORT2188
C  CHECK UPWARD THROUGH THE OFF-DIAGONAL ELEMENTS TO FIND THOSE LESS    FORT2189
C  THAN TOL. IF NO OFF-DIAGONAL ELEMENTS LESS THAN TOL ARE FOUND,       FORT2190
C  PERFORM ANOTHER ITERATION.                                           FORT2191
C                                                                       FORT2192
   30 IT = IT-1                                                         FORT2193
      IF (DABS(B(IT)).LE.TOL) GO TO 40                                  FORT2194
      IF (IT-NI) 100,100,30                                             FORT2195
C                                                                       FORT2196
C  BRANCH ACCORDING TO WHETHER THE MATRIX ISOLATED BY THE SMALL         FORT2197
C  OFF-DIAGONAL ELEMENT IS OF DIMENSION ONE, TWO, OR MORE.              FORT2198
C                                                                       FORT2199
   40 IF (NX-IT-2) 50,60,70                                             FORT2200
C                                                                       FORT2201
C  EXTRACT THE EIGENVALUE OF A ONE BY ONE MATRIX BY ADDING BACK         FORT2202
C  THE SHIFT.                                                           FORT2203
C                                                                       FORT2204
   50 A(NX) = A(NX) + SH                                                FORT2205
C                                                                       FORT2206
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY           FORT2207
C  LATER ITERATIONS.                                                    FORT2208
C                                                                       FORT2209
      NX = NX-1                                                         FORT2210
C                                                                       FORT2211
C  RESET THE ITERATION COUNTER.                                         FORT2212
C                                                                       FORT2213
      K = 1                                                             FORT2214
      GO TO 80                                                          FORT2215
C                                                                       FORT2216
C  EXTRACT THE EIGENVALUES FROM A TWO BY TWO MATRIX.                    FORT2217
C                                                                       FORT2218
   60 AL = B(NX-1)                                                      FORT2219
      AM = 5.D-1*(A(NX-1)-A(NX))                                        FORT2220
      AMS = AM**2                                                       FORT2221
      SAM = DSIGN(1.D0,AM)                                              FORT2222
      AN = DSQRT(AL+AM**2)                                              FORT2223
      CX = (AN+DABS(AM))/(2.D0*AN)                                      FORT2224
      SX = B(NX-1)/(4.D0*AN**2*CX)                                      FORT2225
      TA = A(NX-1)                                                      FORT2226
      TB = A(NX)                                                        FORT2227
      TC = B(NX-1)                                                      FORT2228
      I = NX                                                            FORT2229
C                                                                       FORT2230
C  ROTATE THE DIAGONAL ELEMENTS AND THE OFF-DIAGONAL ELEMENTS.          FORT2231
C                                                                       FORT2232
      A(NX-1) = TA*CX+TB*SX+TC*SAM/AN+SH                                FORT2233
      A(NX) = TA*SX+TB*CX-TC*SAM/AN+SH                                  FORT2234
      B(NX-1) = 4.D0*AMS*CX*SX-DABS(AM)*TC/AN+TC*(CX-SX)**2             FORT2235
C                                                                       FORT2236
C  RESET THE ITERATION COUNTER.                                         FORT2237
C                                                                       FORT2238
      K = 1                                                             FORT2239
C                                                                       FORT2240
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY THE LATER FORT2241
C  ITERATIONS.                                                          FORT2242
C                                                                       FORT2243
      NX = NX-2                                                         FORT2244
      GO TO 80                                                          FORT2245
C                                                                       FORT2246
C  THE NEXT STATEMENT IS REACHED WHEN THE PORTION OF THE MATRIX         FORT2247
C  ISOLATED IS GREATER THAN TWO BY TWO.  IT CHANGES THE LOWER LIMIT     FORT2248
C  OF THE ITERATION SO THAT ONLY THIS PORTION WILL BE AFFECTED BY       FORT2249
C  SUBSEQUENT ROTATIONS UNTIL ALL ITS EIGENVALUES ARE FOUND.            FORT2250
C                                                                       FORT2251
   70 NI = IT + 1                                                       FORT2252
C                                                                       FORT2253
C  TRANSFER TO BEGINNING OF ANOTHER ITERATION.                          FORT2254
C                                                                       FORT2255
      GO TO 85                                                          FORT2256
C                                                                       FORT2257
C  NEXT STATEMENT IS REACHED AFTER EITHER ONE OR TWO EIGENVALUES        FORT2258
C  HAVE JUST BEEN FOUND.  IT TRANSFERS IF ALL THE EIGENVALUES IN        FORT2259
C  THIS PORTION OF THE MATRIX HAVE BEEN FOUND.                          FORT2260
C                                                                       FORT2261
   80 IF (NX.LT.NI) GO TO 90                                            FORT2262
C                                                                       FORT2263
C  BRANCH ACCORDING TO WHETHER ONE, TWO, OR MORE EIGENVALUES REMAIN     FORT2264
C  TO BE FOUND.                                                         FORT2265
C                                                                       FORT2266
   95 IF (NX-NI-1) 50,60,85                                             FORT2267
C                                                                       FORT2268
C  THE NEXT STATEMENT IS REACHED WHEN ALL EIGENVALUES IN THIS PART      FORT2269
C  OF THE MATRIX HAVE BEEN FOUND. IT RETURNS IF THIS IS THE LAST        FORT2270
C  PART OF THE MATRIX.                                                  FORT2271
C                                                                       FORT2272
   90 IF (NI.EQ.1) RETURN                                               FORT2273
C                                                                       FORT2274
C  ENLARGE THE PORTION OF THE MATRIX BEING TREATED TO INCLUDE THE       FORT2275
C  BEGINNING OF THE MATRIX.                                             FORT2276
C                                                                       FORT2277
      NI = 1                                                            FORT2278
      GO TO 95                                                          FORT2279
      END                                                               FORT2280
      SUBROUTINE QRTN (A,B,V,EIG,N,M,TOL,NJX)                           FORT2281
C                                                                       FORT2282
C  SUBROUTINE FOR PERFORMING A QR TRANSFORM ON A REAL SYMMETRIC         FORT2283
C  TRIDIAGONAL MATRIX.                                                  FORT2284
C                                                                       FORT2285
C  N IS THE DIMENSION OF THE MATRIX. A CONTAINS THE N DIAGONAL          FORT2286
C  ELEMENTS OF THE TRIDIAGONAL MATRIX. B CONTAINS THE N-1 OFF-          FORT2287
C  DIAGONAL ELEMENTS. M IS THE MAXIMUM NUMBER OF ITERATIONS ( SAY       FORT2288
C  20 ). NJX IS THE PHYSICAL ROW DIMENSION OF THE EIGENVECTOR           FORT2289
C  MATRIX V.                                                            FORT2290
C                                                                       FORT2291
C  THE EIGENVALUES ARE ASSUMED KNOWN AND PLACED IN THE FIRST N          FORT2292
C  ELEMENTS OF EIG.                                                     FORT2293
C                                                                       FORT2294
C  N VECTORS V ( EACH OF LENGTH N ) ARE TRANSFORMED INTO THE            FORT2295
C  BASIS IN WHICH THE MATRIX IS DIAGONAL.                               FORT2296
C                                                                       FORT2297
      IMPLICIT REAL*8 (A-H,O-Z)                                         FORT2298
      DIMENSION A(1),B(1),V(1),EIG(1)                                   FORT2299
      NX = N                                                            FORT2300
      NNM1 = NJX*(NX-1)                                                 FORT2301
      NI = 1                                                            FORT2302
C                                                                       FORT2303
C  SET INITIAL TOTAL SHIFT.                                             FORT2304
C                                                                       FORT2305
      SH = 0                                                            FORT2306
      IF (NX-2) 50,60,1                                                 FORT2307
C                                                                       FORT2308
C  K COUNTS THE NUMBER OF ITERATIONS PER EIGENVALUE.                    FORT2309
C                                                                       FORT2310
    1 K=0                                                               FORT2311
C                                                                       FORT2312
C  SET INITIAL TEMPORARY SHIFT.                                         FORT2313
C                                                                       FORT2314
   98 ST = EIG(NX)-SH                                                   FORT2315
C                                                                       FORT2316
C  CHECK FOR SMALL OFF-DIAGONAL ELEMENTS.                               FORT2317
C                                                                       FORT2318
      IT = NX                                                           FORT2319
   99 IT = IT-1                                                         FORT2320
      IF (DABS(B(IT)).LE.TOL) GO TO 40                                  FORT2321
      IF (IT.GT.NI) GO TO 99                                            FORT2322
C                                                                       FORT2323
C  NO SMALL OFF-DIAGONAL ELEMENTS FOUND.  ITERATE.  EACH NEW            FORT2324
C  ITERATION BEGINS HERE.                                               FORT2325
C                                                                       FORT2326
  100 K = K + 1                                                         FORT2327
      IF (K.LT.M) GO TO 11                                              FORT2328
      WRITE (6,1000) K                                                  FORT2329
 1000 FORMAT(35H NO CONVERGENCE OF QR ALGORITHM IN   ,I4,11H ITERATIONS)
      CALL EXIT                                                         FORT2331
   11 IF (K.EQ.1) GO TO 15                                              FORT2332
C                                                                       FORT2333
C  SOLVE THE TWO BY TWO IN THE LOWER RIGHT CORNER AND USE SMALLER ROOT  FORT2334
C  AS THE NEXT SHIFT.                                                   FORT2335
C                                                                       FORT2336
   12 AT = A(NX) + A(NX-1)                                              FORT2337
      ST = AT*5.D-1                                                     FORT2338
      DISC = AT**2-4.D0*(A(NX)*A(NX-1)-B(NX-1)**2)                      FORT2339
      IF (DISC.LE.0.D0) GO TO 15                                        FORT2340
      ST = ST-DSIGN(DSQRT(DISC),ST)*5.D-1                               FORT2341
C                                                                       FORT2342
C  INCREASE THE TOTAL SHIFT BY THE TEMPORARY SHIFT.                     FORT2343
C                                                                       FORT2344
   15 SH = SH + ST                                                      FORT2345
C                                                                       FORT2346
C  THIS LOOP SUBTRACTS THE TEMPORARY SHIFT FROM THE DIAGONAL ELEMENTS.  FORT2347
C                                                                       FORT2348
      DO 20 I=1,NX                                                      FORT2349
   20 A(I) = A(I) - ST                                                  FORT2350
      R = DSQRT(A(NI)**2+B(NI)**2)                                      FORT2351
      S = B(NI)/R                                                       FORT2352
      CS = S                                                            FORT2353
      C = A(NI)/R                                                       FORT2354
      U = (S**2)*(A(NI)+A(NI+1))                                        FORT2355
      A(NI) = A(NI) + U                                                 FORT2356
      CALL ROTATE(V(NI),C,S,NJX,NNM1)                                   FORT2357
      NTOP = NX - 2                                                     FORT2358
C                                                                       FORT2359
C  THIS LOOP COMPLETES ONE ITERATION, THAT IS ONE QR TRANSFORM.         FORT2360
C                                                                       FORT2361
      DO 10 I=NI,NTOP                                                   FORT2362
C                                                                       FORT2363
C  G IS THE GAMMA AND Q IS THE P IN THE NOTATION OF WILKINSON.          FORT2364
C                                                                       FORT2365
      G = A(I+1) - U                                                    FORT2366
      Q = C*A(I+1) - CS*B(I)                                            FORT2367
      R = DSQRT(Q**2+B(I+1)**2)                                         FORT2368
C                                                                       FORT2369
C  ROTATE AN OFF-DIAGONAL ELEMENT.                                      FORT2370
C                                                                       FORT2371
      B(I) = S*R                                                        FORT2372
C                                                                       FORT2373
C  FIND THE NEW SINE AND COSINE FOR THE JACOBI ROTATION. THEN           FORT2374
C  COMPUTE A NEW U.                                                     FORT2375
C                                                                       FORT2376
      S = B(I+1)/R                                                      FORT2377
      CS = C*S                                                          FORT2378
      C = Q/R                                                           FORT2379
      U = (S**2)*(G+A(I+2))                                             FORT2380
C                                                                       FORT2381
C  ROTATE A DIAGONAL ELEMENT.                                           FORT2382
C                                                                       FORT2383
      A(I+1) = G + U                                                    FORT2384
C                                                                       FORT2385
C  ROTATE THE VECTORS.                                                  FORT2386
C                                                                       FORT2387
      CALL ROTATE (V(I+1),C,S,NJX,NNM1)                                 FORT2388
   10 CONTINUE                                                          FORT2389
C                                                                       FORT2390
C  COMPUTE THE LAST OFF DIAGONAL ELEMENT.                               FORT2391
C                                                                       FORT2392
      B(NTOP+1) = S*(C*A(NX)-CS*B(NTOP+1))                              FORT2393
C                                                                       FORT2394
C  COMPUTE THE LAST DIAGONAL ELEMENT.                                   FORT2395
C                                                                       FORT2396
      A(NX) = A(NX) - U                                                 FORT2397
C                                                                       FORT2398
C  END OF ONE ITERATION.                                                FORT2399
C                                                                       FORT2400
   85 IT = NX                                                           FORT2401
C                                                                       FORT2402
C  CHECK UPWARD THROUGH THE OFF DIAGONAL ELEMENTS TO FIND THOSE LESS    FORT2403
C  THAN TOL. IF NO OFF-DIAGONAL ELEMENTS LESS THAN TOL ARE FOUND,       FORT2404
C  PERFORM ANOTHER ITERATION.                                           FORT2405
C                                                                       FORT2406
   30 IT = IT - 1                                                       FORT2407
      IF (DABS(B(IT)).LE.TOL) GO TO 40                                  FORT2408
      IF (IT-NI) 100,100,30                                             FORT2409
C                                                                       FORT2410
C  BRANCH ACCORDING TO WHETHER THE MATRIX ISOLATED BY THE SMALL         FORT2411
C  OFF-DIAGONAL ELEMENT IS OF DIMENSION ONE, TWO, OR MORE.              FORT2412
C                                                                       FORT2413
   40 IF (NX-IT-2) 50,60,70                                             FORT2414
C                                                                       FORT2415
C  EXTRACT THE EIGENVALUE OF A ONE BY ONE MATRIX BY ADDING BACK         FORT2416
C  THE SHIFT.                                                           FORT2417
C                                                                       FORT2418
   50 A(NX) = A(NX) + SH                                                FORT2419
C                                                                       FORT2420
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY           FORT2421
C  LATER ITERATIONS.                                                    FORT2422
C                                                                       FORT2423
      NX = NX - 1                                                       FORT2424
C                                                                       FORT2425
C  RESET THE ITERATION COUNTER.                                         FORT2426
C                                                                       FORT2427
      K = 0                                                             FORT2428
      GO TO 80                                                          FORT2429
C                                                                       FORT2430
C  EXTRACT THE EIGENVALUES FROM A TWO BY TWO MATRIX AND PERFORM THE     FORT2431
C  CORRESPONDING ROTATIONS ON THE VECTORS.                              FORT2432
C                                                                       FORT2433
   60 AL = -B(NX-1)                                                     FORT2434
      AM = 5.D-1*(A(NX-1)-A(NX))                                        FORT2435
      AN = DSQRT(AL**2+AM**2)                                           FORT2436
      C = DSQRT((AN+DABS(AM))/(2.D0*AN))                                FORT2437
      S = DSIGN(5.D-1,AM)*AL/(AN*C)                                     FORT2438
      TA = A(NX-1)                                                      FORT2439
      TB = A(NX)                                                        FORT2440
      TC = B(NX-1)                                                      FORT2441
      CX = C**2                                                         FORT2442
      SX = S**2                                                         FORT2443
      CS = C*S                                                          FORT2444
C                                                                       FORT2445
C  ROTATE THE DIAGONAL ELEMENTS, THE OFF-DIAGONAL ELEMENTS, AND         FORT2446
C  THE VECTORS.                                                         FORT2447
C                                                                       FORT2448
      A(NX-1) = TA*CX+TB*SX-2.D0*TC*CS+SH                               FORT2449
      A(NX) = TA*SX+TB*CX+2.D0*TC*CS+SH                                 FORT2450
      B(NX-1) = 2.D0*AM*CS+TC*(CX-SX)                                   FORT2451
      I = NX-1                                                          FORT2452
      S = -S                                                            FORT2453
      CALL ROTATE (V(I),C,S,NJX,NNM1)                                   FORT2454
C                                                                       FORT2455
C  RESET THE ITERATION COUNTER.                                         FORT2456
C                                                                       FORT2457
      K = 0                                                             FORT2458
C                                                                       FORT2459
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY THE       FORT2460
C  LATER ITERATIONS.                                                    FORT2461
C                                                                       FORT2462
      NX = NX-2                                                         FORT2463
      GO TO 80                                                          FORT2464
C                                                                       FORT2465
C  THE NEXT STATEMENT IS REACHED WHEN THE PORTION OF THE MATRIX         FORT2466
C  ISOLATED IS GREATER THAN TWO BY TWO.  IT CHANGES THE LOWER LIMIT     FORT2467
C  OF THE ITERATION SO THAT ONLY THIS PORTION WILL BE AFFECTED BY       FORT2468
C  SUBSEQUENT ROTATIONS UNTIL ALL ITS EIGENVALUES ARE FOUND.            FORT2469
C                                                                       FORT2470
   70 NI = IT + 1                                                       FORT2471
C                                                                       FORT2472
C  TRANSFER TO BEGINNING OF ANOTHER ITERATION.                          FORT2473
C                                                                       FORT2474
      GO TO 100                                                         FORT2475
C                                                                       FORT2476
C  NEXT STATEMENT IS REACHED AFTER EITHER ONE OR TWO EIGENVALUES        FORT2477
C  HAVE JUST BEEN FOUND.  IT TRANSFERS IF ALL THE EIGENVALUES IN        FORT2478
C  THIS PORTION OF THE MATRIX HAVE BEEN FOUND.                          FORT2479
C                                                                       FORT2480
   80 IF (NX.LT.NI) GO TO 90                                            FORT2481
C                                                                       FORT2482
C  BRANCH ACCORDING TO WHETHER ONE, TWO, OR MORE EIGENVALUES REMAIN     FORT2483
C  TO BE FOUND.                                                         FORT2484
C                                                                       FORT2485
   95 IF (NX-NI-1) 50,60,98                                             FORT2486
C                                                                       FORT2487
C  THE NEXT STATEMENT IS REACHED WHEN ALL EIGENVALUES IN THIS PART      FORT2488
C  OF THE MATRIX HAVE BEEN FOUND. IT RETURNS IF THIS IS THE LAST        FORT2489
C  PART OF THE MATRIX.                                                  FORT2490
C                                                                       FORT2491
   90 IF (NI.EQ.1) RETURN                                               FORT2492
C                                                                       FORT2493
C  ENLARGE THE PORTION OF THE MATRIX BEING TREATED TO INCLUDE THE       FORT2494
C  BEGINNING OF THE MATRIX.                                             FORT2495
C                                                                       FORT2496
      NI = 1                                                            FORT2497
      GO TO 95                                                          FORT2498
      END                                                               FORT2499
      SUBROUTINE ITRATE(H,U,MAD,C,E,W,IOCC,HDG,NDIM,NTYPE,NC,NHDG)      FORT2500
C                                                                       FORT2501
C  SUBROUTINE TO SETUP HUCKEL MATRIX WHEN USING CHARGE ITERATION        FORT2502
C  OPTION ( METH = 2 OR 3 ).                                            FORT2503
C                                                                       FORT2504
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT2505
      DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),       FORT2506
     1 E(NDIM),W(NDIM),IOCC(NDIM),HDG(NHDG)                             FORT2507
      REAL*8 MAD                                                        FORT2508
      REAL*4 IOCC                                                       FORT2509
      PARAMETER (MAXATM=500)
      PARAMETER (BB=250)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      COMMON/TITLE/AB(10)                                               FORT2510
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,      FORT2511
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                              FORT2512
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                            FORT2513
      COMMON/OUT/PRT(20),PUN(20)                                        FORT2514
      LOGICAL*1 PRT,PUN                                                 FORT2515
      COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),ND(BB), 
     1 EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),COULS(BB),
     2 COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)                      
      INTEGER*2 SYMBOL,KEY,VELEC                                        FORT2519
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)                             FORT2521
      REAL*8 LAMPRI                                                     FORT2522
      INTEGER*4 PRTCYC                                                  FORT2523
      LOGICAL*1 PARTIT,PRINTX,ITABLE                                    FORT2524
C                                                                       FORT2525
C  SINCE THE INTERNAL ATOMIC PARAMETERS ( EXPS,EXPP, ETC. ) ARE         FORT2526
C  NOT USED WHEN DOING CHARGE ITERATION, THE SPACE ALLOCATED TO         FORT2527
C  THEM HAS BEEN USED FOR THE VSIE PARAMETERS.                          FORT2528
C                                                                       FORT2529
      DIMENSION AS1(MXUSER),BS1(MXUSER),CS1(MXUSER),AP1(MXUSER),
     . BP1(MXUSER),CP1(MXUSER),AD1(MXUSER),BD1(MXUSER),CD1(MXUSER)  
      EQUIVALENCE (AS1(1),EXPS(MXUSR2)),(BS1(1),EXPP(MXUSR2)),
     . (CS1(1),EXPD(MXUSR2)),(AP1(1),EXPD2(MXUSR2)),
     . (BP1(1),C1(MXUSR2)),(CP1(1),C2(MXUSR2)),
     . (AD1(1),COULS(MXUSR2)),(BD1(1),COULP(MXUSR2)),
     . (CD1(1),COULD(MXUSR2))
      COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5),      FORT2535
     1 BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),AD3(5),  FORT2536
     2 BD3(5),CD3(5)                                                    FORT2537
      REAL*4 ZS(MXUSER),ZP(MXUSER),ZD(MXUSER)                       
      EQUIVALENCE (ZS(1),NS(MXUSR2)),(ZP(1),NP(MXUSR2)),
     . (ZD(1),ND(MXUSR2)) 
      DIMENSION IL(3),JL(3)                                             FORT2540
      EQUIVALENCE (PEEP,ADJUST),(COULH,SIGNN),(SENSE,DENSE),            FORT2541
     1 (VELEC(81),DCHG),(PRT(1),NMIN),(PUN(1),NLN)                      FORT2542
      DIMENSION F(6)                                                    FORT2543
      EQUIVALENCE (A1,F(1)),(A2,F(3)),(A3,F(5))                         FORT2544
      REAL*8 LAMBDA                                                     FORT2545
      QON=DFLOAT(KA)/DFLOAT(NELEC)                                      FORT2546
      IF(ICYCLE.GT.1) GO TO 696                                         FORT2547
      WRITE(6,750)                                                      FORT2548
750   FORMAT('1')                                                       FORT2549
      DO 5 I=1,NDIM                                                     FORT2550
      NMIN=I                                                            FORT2551
      IF(IOCC(I).GT.0.0001) GO TO 697                                   FORT2552
5     CONTINUE                                                          FORT2553
C                                                                       FORT2554
C  IF CHARGE ITERATION WITH MADELUNG CORRECTION ( METH>2 ) IS BEING     FORT2555
C  USED, COMPUTE THE TOTAL ORBITAL OCCUPATIONS FOR THE FREE ATOMS.      FORT2556
C  FILL ATOMIC ORBITALS IN ORDER0 1S 2S 2P 3S 3P 3D 4S 4P 4D 5S 5P      FORT2557
C  5D 6S 6P.                                                            FORT2558
C                                                                       FORT2559
697   IF(METH.LT.3) GO TO 696                                           FORT2560
      PRTCYC=MAXCYC
      DO 700 I=1,NA                                                     FORT2561
      KEYI=KEY(I)                                                       FORT2562
      J=NP(KEYI)                                                        FORT2563
      IF(J.EQ.0) GO TO 701                                              FORT2564
      IL(2)=J*10+1                                                      FORT2565
      JL(2)=2                                                           FORT2566
      J=NS(KEYI)                                                        FORT2567
      IL(1)=J*10                                                        FORT2568
      JL(1)=1                                                           FORT2569
      J=ND(KEYI)                                                        FORT2570
      IF(J.EQ.0) GO TO 702                                              FORT2571
      IL(3)=J*10+2                                                      FORT2572
      JL(3)=3                                                           FORT2573
      DO 710 J=2,3                                                      FORT2574
      IJ=4-J                                                            FORT2575
      K=IJ+1                                                            FORT2576
      DO 710 L=1,IJ                                                     FORT2577
      N=K-L                                                             FORT2578
      IF(IL(K).GT.IL(N)) GO TO 710                                      FORT2579
      I1=IL(K)                                                          FORT2580
      J1=JL(K)                                                          FORT2581
      IL(K)=IL(N)                                                       FORT2582
      JL(K)=JL(N)                                                       FORT2583
      IL(N)=I1                                                          FORT2584
      JL(N)=J1                                                          FORT2585
710   CONTINUE                                                          FORT2586
      GO TO 720                                                         FORT2587
701   JL(1)=1                                                           FORT2588
      JL(2)=2                                                           FORT2589
      JL(3)=3                                                           FORT2590
      GO TO 720                                                         FORT2591
702   JL(3)=3                                                           FORT2592
      IF(IL(1).LT.IL(2)) GO TO 720                                      FORT2593
      JL(2)=1                                                           FORT2594
      JL(1)=2                                                           FORT2595
720   L=VELEC(KEYI)                                                     FORT2596
      DO 700 J=1,3                                                      FORT2597
      K=JL(J)                                                           FORT2598
      J1=4*(K-1)+2                                                      FORT2599
      IF(L.LT.J1) J1=L                                                  FORT2600
      L=L-J1                                                            FORT2601
      IF(K-2) 725,726,727                                               FORT2602
725   ZS(KEYI)=FLOAT(J1)                                                FORT2603
      GO TO 700                                                         FORT2604
726   ZP(KEYI)=FLOAT(J1)                                                FORT2605
      GO TO 700                                                         FORT2606
727   ZD(KEYI)=FLOAT(J1)                                                FORT2607
700   CONTINUE                                                          FORT2608
696   IF(PRINTX) WRITE(6,751)                                           FORT2609
751   FORMAT(////)                                                      FORT2610
      PRINTX=((ICYCLE/PRTCYC)*PRTCYC.EQ.ICYCLE)                         FORT2611
      IF(ICYCLE.EQ.MAXCYC) PRINTX=.TRUE.                                FORT2612
C                                                                       FORT2613
C  CALCULATE SUM OF ONE-ELECTRON ENERGIES.                              FORT2614
C                                                                       FORT2615
753   SUM=0.0D0                                                         FORT2616
      DO 13 I=1,NDIM                                                    FORT2617
      W(I)=DBLE(IOCC(I))                                                FORT2618
13    SUM=SUM+E(I)*W(I)                                                 FORT2619
C                                                                       FORT2620
C  COMPUTE ATOMIC ORBITAL OCCUPATIONS AND STORE IN E(I).                FORT2621
C                                                                       FORT2622
3004  IJ=1                                                              FORT2623
      DO  319  I=1,NDIM                                                 FORT2624
      E(I)=0.0D0                                                        FORT2625
      DO  319  J=1,I                                                    FORT2626
      UB=0.0D0                                                          FORT2627
      DO 320 K=NMIN,NDIM                                                FORT2628
320   UB=UB+H(I,K)*H(J,K)*W(K)                                          FORT2629
      UB=UB*0.50D0                                                      FORT2630
      IF(I.EQ.J)  GO  TO  28                                            FORT2631
      UB=(UB+UB)*C(IJ)                                                  FORT2632
      IJ=IJ+1                                                           FORT2633
28    E(I)=E(I)+UB                                                      FORT2634
      E(J)=E(J)+UB                                                      FORT2635
319   CONTINUE                                                          FORT2636
C                                                                       FORT2637
C  COMPUTE THE ORBITAL OCCUPATION OF A GIVEN TYPE (S,P,D) WHICH         FORT2638
C  VARIES MOST FROM THE LAST CYCLE.                                     FORT2639
C                                                                       FORT2640
3005  DENOM=0.0D0                                                       FORT2641
      DENSE2=DENSE                                                      FORT2642
      DCHG2=DCHG                                                        FORT2643
      J=1                                                               FORT2644
      DO  5000 I=1,NA                                                   FORT2645
      KEYI=KEY(I)                                                       FORT2646
      SDENSE=E(J)                                                       FORT2647
      SIGNS=SDENSE-X(I)                                                 FORT2648
      C(J)=SIGNS                                                        FORT2649
      DIFFS=DABS(SIGNS)                                                 FORT2650
      N=J                                                               FORT2651
      IF(DENOM.GT.DIFFS) GO TO 5001                                     FORT2652
      DENSE=X(I)                                                        FORT2653
      DCHG=SIGNS                                                        FORT2654
      DENOM=DIFFS                                                       FORT2655
      NLF=10*I                                                          FORT2656
5001  IF(NP(KEYI).EQ.0) GO TO 5000                                      FORT2657
      PDENSE=E(J+1)+E(J+2)+E(J+3)                                       FORT2658
      SIGNP=PDENSE-Y(I)                                                 FORT2659
      C(J+3)=SIGNP                                                      FORT2660
      DIFFP=DABS(SIGNP)                                                 FORT2661
      N=J+3                                                             FORT2662
      IF(DENOM.GT.DIFFP) GO TO 5002                                     FORT2663
      DENSE=Y(I)                                                        FORT2664
      DCHG=SIGNP                                                        FORT2665
      DENOM=DIFFP                                                       FORT2666
      NLF=10*I+1                                                        FORT2667
5002  IF(ND(KEYI).EQ.0) GO TO 5000                                      FORT2668
      DDENSE=E(J+4)+E(J+5)+E(J+6)+E(J+7)+E(J+8)                         FORT2669
      SIGND=DDENSE-Z(I)                                                 FORT2670
      C(J+8)=SIGND                                                      FORT2671
      DIFFD=DABS(SIGND)                                                 FORT2672
      N=J+8                                                             FORT2673
      IF(DENOM.GT.DIFFD) GO TO 5000                                     FORT2674
      DENSE=Z(I)                                                        FORT2675
      DCHG=SIGND                                                        FORT2676
      DENOM=DIFFD                                                       FORT2677
      NLF=10*I+2                                                        FORT2678
5000  J=N+1                                                             FORT2679
      SIGNF=DCHG/DENOM                                                  FORT2680
C                                                                       FORT2681
C  DETERMINE LAMBDA, THE DAMPING FACTOR.                                FORT2682
C                                                                       FORT2683
      IF(ICYCLE.LE.2) SIGNN=SIGNF                                       FORT2684
      IF(ICYCLE.LE.2) NLN=NLF                                           FORT2685
      ISIGNN=SIGNN                                                      FORT2686
      ISIGNF=SIGNF                                                      FORT2687
      IF(ICYCLE.EQ.2) ADJUST=DAMP1*DENOM                                FORT2688
      IF(ISIGNN.EQ.ISIGNF.OR.NLN.NE.NLF) GO TO 8001                     FORT2689
      IF(DAMP3.NE.0.0D0) GO TO 401                                      FORT2690
      LAMBDA=(DENSE2-DENSE)/(DCHG-DCHG2)                                FORT2691
      GO TO 402                                                         FORT2692
401   ADJUST=DAMP3*ADJUST                                               FORT2693
8001  IF((ADJUST/DENOM).GE.LAMPRI) ADJUST=DAMP2*DENOM                   FORT2694
      IF(ICYCLE.EQ.1) ADJUST=DENOM                                      FORT2695
      LAMBDA=ADJUST/DENOM                                               FORT2696
402   SIGNN=SIGNF                                                       FORT2697
      NLN=NLF                                                           FORT2698
C                                                                       FORT2699
C                                                                       FORT2700
C  IN THIS SECTION OF THE PROGRAM THREE CALCULATIONS ARE PERFORMED0     FORT2701
C                                                                       FORT2702
C    1. THE TOTAL ORBITAL OCCUPATIONS OF A GIVEN TYPE (S,P,D) ARE       FORT2703
C       DAMPED AND STORED (IN X(I),Y(I),Z(I) RESPECTIVELY).             FORT2704
C    2. THE NET CHARGES ARE CALCULATED AND STORED IN C(I).              FORT2705
C    3. -VSIE'S ARE CALCULATED AND STORED IN W(I).                      FORT2706
C                                                                       FORT2707
C                                                                       FORT2708
      J=1                                                               FORT2709
      DO 803 I=1,NA                                                     FORT2710
      KEYI=KEY(I)                                                       FORT2711
      SDENSE=X(I)+LAMBDA*C(J)                                           FORT2712
      X(I)=SDENSE                                                       FORT2713
      UB=SDENSE                                                         FORT2714
      N=J                                                               FORT2715
      IF(NP(KEYI).EQ.0) GO TO 804                                       FORT2716
      PDENSE=Y(I)+LAMBDA*C(J+3)                                         FORT2717
      Y(I)=PDENSE                                                       FORT2718
      UB=UB+PDENSE                                                      FORT2719
      N=J+3                                                             FORT2720
      IF(ND(KEYI).EQ.0) GO TO 805                                       FORT2721
      N=J+8                                                             FORT2722
      DDENSE=Z(I)+LAMBDA*C(J+8)                                         FORT2723
      Z(I)=DDENSE                                                       FORT2724
      UB=UB+DDENSE                                                      FORT2725
      GO TO 806                                                         FORT2726
804   Q=VELEC(KEYI)-UB                                                  FORT2727
      IF(.NOT.PARTIT) GO TO 1111                                        FORT2728
      IF(ITABLE(KEYI)) GO TO 1111                                       FORT2729
      W(J)=COULS(KEYI)                                                  FORT2730
      GO TO 807                                                         FORT2731
1111  KEYI=MXUSR2-KEY(I)                                            
      W(J)=-((AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI))                       FORT2733
      GO TO 807                                                         FORT2734
805   Q=VELEC(KEYI)-UB                                                  FORT2735
      IF(.NOT.PARTIT) GO TO 1113                                        FORT2736
      IF(ITABLE(KEYI)) GO TO 1113                                       FORT2737
      W(J)=COULS(KEYI)                                                  FORT2738
      W(J+1)=COULP(KEYI)                                                FORT2739
      W(J+2)=COULP(KEYI)                                                FORT2740
      W(J+3)=COULP(KEYI)                                                FORT2741
      GO TO 807                                                         FORT2742
1113  KEYI=MXUSR2-KEY(I)                                               
      VSIES1=(AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI)                        FORT2744
      W(J)=-VSIES1                                                      FORT2745
      VSIEP1=(AP1(KEYI)*Q+BP1(KEYI))*Q+CP1(KEYI)                        FORT2746
      W(J+1)=-VSIEP1                                                    FORT2747
      W(J+2)=-VSIEP1                                                    FORT2748
      W(J+3)=-VSIEP1                                                    FORT2749
      GO TO 807                                                         FORT2750
806   Q=VELEC(KEYI)-UB                                                  FORT2751
      IF(.NOT.PARTIT) GO TO 1115                                        FORT2752
      IF(ITABLE(KEYI)) GO TO 1115                                       FORT2753
      W(J)=COULS(KEYI)                                                          
      W(J+1)=COULP(KEYI)                                                FORT2754
      W(J+2)=COULP(KEYI)                                                FORT2755
      W(J+3)=COULP(KEYI)                                                FORT2756
      W(J+4)=COULD(KEYI)                                                FORT2757
      W(J+5)=COULD(KEYI)                                                FORT2758
      W(J+6)=COULD(KEYI)                                                FORT2759
      W(J+7)=COULD(KEYI)                                                FORT2760
      W(J+8)=COULD(KEYI)                                                FORT2761
      GO TO 807                                                         FORT2762
1115  IF(ND(KEYI).EQ.NP(KEYI)) GO TO 1116                               FORT2763
      IF(NCON.NE.3) GO TO 1116                                          FORT2764
      KEYI=MXUSR2-KEY(I)                                      
      VSIES1=(AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI)                        FORT2766
      VSIES2=(AS2(KEYI)*Q+BS2(KEYI))*Q+CS2(KEYI)                        FORT2767
      VSIES3=(AS3(KEYI)*Q+BS3(KEYI))*Q+CS3(KEYI)                        FORT2768
      W(J)=(SDENSE+PDENSE-2.0D0)*VSIES1+(1.0D0-SDENSE)*VSIES2-          FORT2769
     *PDENSE*VSIES3                                                     FORT2770
      VSIEP1=(AP1(KEYI)*Q+BP1(KEYI))*Q+CP1(KEYI)                        FORT2771
      VSIEP2=(AP2(KEYI)*Q+BP2(KEYI))*Q+CP2(KEYI)                        FORT2772
      VSIEP3=(AP3(KEYI)*Q+BP3(KEYI))*Q+CP3(KEYI)                        FORT2773
      W(J+1)=(SDENSE+PDENSE-2.0D0)*VSIEP1+(1.0D0-PDENSE)*VSIEP2-        FORT2774
     *SDENSE*VSIEP3                                                     FORT2775
      W(J+2)=W(J+1)                                                     FORT2776
      W(J+3)=W(J+1)                                                     FORT2777
      VSIED1=(AD1(KEYI)*Q+BD1(KEYI))*Q+CD1(KEYI)                        FORT2778
      VSIED2=(AD2(KEYI)*Q+BD2(KEYI))*Q+CD2(KEYI)                        FORT2779
      VSIED3=(AD3(KEYI)*Q+BD3(KEYI))*Q+CD3(KEYI)                        FORT2780
      W(J+4)=(SDENSE+PDENSE-1.0D0)*VSIED1-SDENSE*VSIED2-PDENSE*VSIED3   FORT2781
      W(J+5)=W(J+4)                                                     FORT2782
      W(J+6)=W(J+4)                                                     FORT2783
      W(J+7)=W(J+4)                                                     FORT2784
      W(J+8)=W(J+4)                                                     FORT2785
      GO TO 807                                                         FORT2786
1116  KEYI=MXUSR2-KEY(I)                                                
      VSIES1=(AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI)                        FORT2788
      W(J)=-VSIES1                                                      FORT2789
      VSIEP1=(AP1(KEYI)*Q+BP1(KEYI))*Q+CP1(KEYI)                        FORT2790
      W(J+1)=-VSIEP1                                                    FORT2791
      W(J+2)=W(J+1)                                                     FORT2792
      W(J+3)=W(J+1)                                                     FORT2793
      VSIED1=(AD1(KEYI)*Q+BD1(KEYI))*Q+CD1(KEYI)                        FORT2794
      W(J+4)=-VSIED1                                                    FORT2795
      W(J+5)=W(J+4)                                                     FORT2796
      W(J+6)=W(J+4)                                                     FORT2797
      W(J+7)=W(J+4)                                                     FORT2798
      W(J+8)=W(J+4)                                                     FORT2799
807   J=N+1                                                             FORT2800
      C(I)=Q                                                            FORT2801
803   CONTINUE                                                          FORT2802
C                                                                       FORT2803
C  IF CHARGE ITERATION WITHOUT MADELUNG CORRECTION ( METH=2 ) IS        FORT2804
C  BEING USED, SETUP HUCKEL MATRIX. OTHERWISE SKIP THIS SECTION.        FORT2805
C                                                                       FORT2806
      IF(METH.GT.2) GO TO 999                                           FORT2807
      H(1,1)=W(1)                                                       FORT2808
      CNST=CON                                                          FORT2809
      DO 760 I=2,NDIM                                                   FORT2810
      H(I,I)=W(I)                                                       FORT2811
      J1=I-1                                                            FORT2812
      DO 760 J=1,J1                                                     FORT2813
      UB=W(I)+W(J)                                                      FORT2814
      IF(.NOT.L5) GO TO 761                                             FORT2815
      UC=(W(I)-W(J))/UB                                                 FORT2816
      UC=UC*UC                                                          FORT2817
      CNST=CON+UC/2.0D0+UC*UC*(0.5D0-CON)                               FORT2818
761   UB=CNST*U(I,J)*UB                                                 FORT2819
      H(J,I)=UB                                                         FORT2820
760   H(I,J)=UB                                                         FORT2821
      GO TO 850                                                         FORT2822
C                                                                       FORT2823
C  IF CHARGE ITERATION WITH MADELUNG CORRECTION ( METH>2 ) IS BEING     FORT2824
C  USED, SETUP HUCKEL MATRIX. OTHERWISE SKIP THIS SECTION.              FORT2825
C                                                                       FORT2826
999   DGSUM=0.0D0                                                       FORT2827
      N=1                                                               FORT2828
      M=1                                                               FORT2829
      DO 880 I=1,NA                                                     FORT2830
      KEYI=KEY(I)                                                       FORT2831
      N1=N                                                              FORT2832
      IF(NP(KEYI).NE.0) N1=N+1                                          FORT2833
      IF(ND(KEYI).NE.0) N1=N+2                                          FORT2834
      DO 881 J=N,N1                                                     FORT2835
      DG1=0.0D0                                                         FORT2836
      DG2=0.0D0                                                         FORT2837
      L=1                                                               FORT2838
      DO 882 K=1,NA                                                     FORT2839
      KEYK=KEY(K)                                                       FORT2840
      UB=(X(K)-DBLE(ZS(KEYK)))*MAD(J,L)                                 FORT2841
      UC=X(K)*MAD(J,L)                                                  FORT2842
      L=L+1                                                             FORT2843
      IF(NP(KEYK).EQ.0) GO TO 883                                       FORT2844
      UB=UB+(Y(K)-DBLE(ZP(KEYK)))*MAD(J,L)                              FORT2845
      UC=UC+Y(K)*MAD(J,L)                                               FORT2846
      L=L+1                                                             FORT2847
      IF(ND(KEYK).EQ.0) GO TO 883                                       FORT2848
      UB=UB+(Z(K)-DBLE(ZD(KEYK)))*MAD(J,L)                              FORT2849
      UC=UC+Z(K)*MAD(J,L)                                               FORT2850
      L=L+1                                                             FORT2851
883   IF(K.NE.I) DG1=DG1+UB                                             FORT2852
882   DG2=DG2+UC                                                        FORT2853
      J1=M+2*(J-N)                                                      FORT2854
      DO 884 L=M,J1                                                     FORT2855
      H(L,L)=W(L)+DG1                                                   FORT2856
      IF(L5) HDG(L)=H(L,L)+QON*DG2                                      FORT2857
      W(L)=DG2                                                          FORT2858
884   DGSUM=DGSUM+DG2                                                   FORT2859
881   M=J1+1                                                            FORT2860
880   N=N1+1                                                            FORT2861
      CNST=CON                                                          FORT2862
      DO 885 I=2,NDIM                                                   FORT2863
      J1=I-1                                                            FORT2864
      DO 885 J=1,J1                                                     FORT2865
      IF(.NOT.L5) GO TO 886                                             FORT2866
      UB=(HDG(I)-HDG(J))/(HDG(I)+HDG(J))                                FORT2867
      UB=UB*UB                                                          FORT2868
      CNST=CON+UB/2.0D0+UB*UB*(0.5D0-CON)                               FORT2869
886   UB=U(I,J)*(CNST*(H(I,I)+H(J,J))-QON*(0.5D0-CNST)*(W(I)+W(J)))     FORT2870
      H(I,J)=UB                                                         FORT2871
885   H(J,I)=UB                                                         FORT2872
      DGSUM=-(QON*DGSUM)/DFLOAT(NDIM)                                   FORT2873
C                                                                       FORT2874
C  IF DOING LAST CYCLE CALCULATE ENERGY CORRECTIONS.                    FORT2875
C                                                                       FORT2876
      IF(ICYCLE.NE.MAXCYC) GO TO 850                                    FORT2877
      N=1                                                               FORT2878
      UB=0.0D0                                                          FORT2879
      UC=0.0D0                                                          FORT2880
      DO 887 I=1,NA                                                     FORT2881
      KEYI=KEY(I)                                                       FORT2882
      K=N                                                               FORT2883
      A1=X(I)                                                           FORT2884
      E(1)=DBLE(ZS(KEYI))                                               FORT2885
      UB=UB-0.5D0*(A1*A1-A1)*MAD(N,N)                                   FORT2886
      IF(NP(KEYI).EQ.0) GO TO 888                                       FORT2887
      N=N+1                                                             FORT2888
      A2=Y(I)                                                           FORT2889
      E(3)=DBLE(ZP(KEYI))                                               FORT2890
      UB=UB-0.5D0*(A2*A2-A2)*MAD(N,N)-A1*A2*MAD(N-1,N)                  FORT2891
      IF(ND(KEYI).EQ.0) GO TO 888                                       FORT2892
      N=N+1                                                             FORT2893
      A3=Z(I)                                                           FORT2894
      E(5)=DBLE(ZD(KEYI))                                               FORT2895
      UB=UB-0.5D0*(A3*A3-A3)*MAD(N,N)-A1*A3*MAD(N-2,N)-A2*A3*MAD(N-1,N) FORT2896
888   M=1                                                               FORT2897
      I1=I-1                                                            FORT2898
      IF(I1.EQ.0) GO TO 887                                             FORT2899
      DO 889 J=1,I1                                                     FORT2900
      KEYJ=KEY(J)                                                       FORT2901
      L=M                                                               FORT2902
      F(2)=X(J)                                                         FORT2903
      E(2)=DBLE(ZS(KEYJ))                                               FORT2904
      IF(NP(KEYJ).EQ.0) GO TO 890                                       FORT2905
      M=M+1                                                             FORT2906
      F(4)=Y(J)                                                         FORT2907
      E(4)=DBLE(ZP(KEYJ))                                               FORT2908
      IF(ND(KEYJ).EQ.0) GO TO 890                                       FORT2909
      M=M+1                                                             FORT2910
      F(6)=Z(J)                                                         FORT2911
      E(6)=DBLE(ZD(KEYJ))                                               FORT2912
890   DO 891 IJ=K,N                                                     FORT2913
      DO 891 JK=L,M                                                     FORT2914
      N1=2*(IJ-K)+1                                                     FORT2915
      M1=2*(JK-L)+2                                                     FORT2916
891   UC=UC-(F(N1)*F(M1)-E(N1)*E(M1))*MAD(IJ,JK)                        FORT2917
889   M=M+1                                                             FORT2918
887   N=N+1                                                             FORT2919
      A1=UB                                                             FORT2920
      A2=UC                                                             FORT2921
      A3=UB+UC                                                          FORT2922
C                                                                       FORT2923
C  SAVE MADELUNG TERMS FOR USE IN SUBROUTINE OUTPUT IF DOING            FORT2924
C  THE LAST CYCLE.                                                      FORT2925
C                                                                       FORT2926
      K=1                                                               FORT2927
      L=1                                                               FORT2928
      DO 792 I=1,NA                                                     FORT2929
      KEYI=KEY(I)                                                       FORT2930
      MAD(K,K)=W(L)                                                     FORT2931
      K=K+1                                                             FORT2932
      L=L+1                                                             FORT2933
      IF(NP(KEYI).EQ.0) GO TO 792                                       FORT2934
      MAD(K,K)=W(L)                                                     FORT2935
      K=K+1                                                             FORT2936
      L=L+3                                                             FORT2937
      IF(ND(KEYI).EQ.0) GO TO 792                                       FORT2938
      MAD(K,K)=W(L)                                                     FORT2939
      K=K+1                                                             FORT2940
      L=L+5                                                             FORT2941
792   CONTINUE                                                          FORT2942
C                                                                       FORT2943
C  PRINT OUT RESULTS TO SHOW PROGRESS OF ITERATION PROCEDURE.           FORT2944
C                                                                       FORT2945
850   IF(PRTCYC.GT.0) WRITE(6,793) ICYCLE                               FORT2946
793   FORMAT('0CYCLE NO.',I3,'0')                                       FORT2947
      IF(PRTCYC.LT.0) WRITE(6,794)                                      FORT2948
794   FORMAT('0CONVERGENCE REACHED - FINAL CYCLE FOLLOWS0',///)         FORT2949
      PRTCYC=IABS(PRTCYC)                                               FORT2950
      J=NLF/10                                                          FORT2951
      K=NLF-10*J                                                        FORT2952
      WRITE(6,795) SUM,LAMBDA,J,ISIGNF,DENOM,ADJUST,K                   FORT2953
795   FORMAT('+',T25,'ENERGY =',F15.8,T52,'LAMBDA =',F8.5,T78,'ATOM =', FORT2954
     1 I3,T92,'SIGN =',I3,/,T25,'DENOM =',D16.8,T52,'ADJUST =',D14.7,   FORT2955
     2 T78,'NL   =',I3)                                                 FORT2956
C                                                                       FORT2957
C  PRINT OUT ATOMIC CHARGES, ORBITAL OCCUPATIONS, AND CORRECTED         FORT2958
C  H(I,I)'S IF PRINTX IS TRUE.                                          FORT2959
C                                                                       FORT2960
      IF(.NOT.PRINTX) GO TO 500                                         FORT2961
      WRITE(6,600)                                                      FORT2962
600   FORMAT(////,T16,'ATOM',T30,'NET CHG.-DAMPED',T60,'SUMMED ORBITAL O
     1CCUPATIONS-DAMPED'/T65,'S',T75,'P',T85,'D'/)                      FORT2964
      DO 650 I=1,NA                                                     FORT2965
      KEYI=KEY(I)                                                       FORT2966
      WRITE(6,601) X(I)                                                 FORT2967
601   FORMAT(T60,F10.5)                                                 FORT2968
      IF(NP(KEYI).EQ.0) GO TO 625                                       FORT2969
      WRITE(6,602) Y(I)                                                 FORT2970
602   FORMAT('+',T70,F10.5)                                             FORT2971
      IF(ND(KEYI).EQ.0) GO TO 625                                       FORT2972
      WRITE(6,603) Z(I)                                                 FORT2973
603   FORMAT('+',T80,F10.5)                                             FORT2974
625   UB=C(I)                                                           FORT2975
650   WRITE(6,604) SYMBOL(KEYI),I,UB                                    FORT2976
604   FORMAT('+',T15,A2,I3,T30,F10.5)                                   FORT2977
      WRITE(6,809)                                                      FORT2978
809   FORMAT(///,T16,'ATOM',T70,'CORRECTED H(I,I)''S',/T40,'S',T50,     FORT2979
     1 'X',T60,'Y',T70,'Z',T80,'X2-Y2',T90,'Z2',T100,'XY',T110,'XZ',    FORT2980
     2 T120,'YZ'/)                                                      FORT2981
      J=1                                                               FORT2982
      DO 810 I=1,NATM                                                   FORT2983
      KEYI=KEY(I)                                                       FORT2984
      N=J                                                               FORT2985
      IF(NP(KEYI).NE.0) N=J+3                                           FORT2986
      IF(ND(KEYI).NE.0) N=J+8                                           FORT2987
815   WRITE(6,816) SYMBOL(KEYI),I,(H(K,K),K=J,N)                        FORT2988
816   FORMAT(T15,A2,I3,T35,9F10.5)                                      FORT2989
810   J=N+1                                                             FORT2990
      IF(METH.GE.3.AND.DABS(QON).GT.0.0001D0) WRITE(6,870) DGSUM        FORT2991
870   FORMAT(///,T15,'AVERAGE SHIFT OF MO''S DUE TO NON-ZERO TOTAL CHARGFORT2992
     1E =',F12.8,' EV.')                                                FORT2993
      IF(METH.GE.3.AND.ICYCLE.EQ.MAXCYC) WRITE(6,808) A1,A2,A3          FORT2994
808   FORMAT(///,T15,'ENERGY CORRECTIONS0',//,T20,'ONE-CENTER',T40,     FORT2995
     1F16.8,' EV.',//,T20,'TWO-CENTER',T40,F16.8,' EV.',//,T20,'TOTAL', FORT2996
     2T40,F16.8,' EV.')                                                 FORT2997
      WRITE(6,752)                                                      FORT2998
752   FORMAT(///)                                                       FORT2999
500   ICYCLE=ICYCLE+1                                                   FORT3000
C                                                                       FORT3001
C  CHECK FOR CONVERGENCE ( IE. DENOM LESS THAN DELTAC ).                FORT3002
C                                                                       FORT3003
      IF(ICYCLE.GE.MAXCYC) GO TO 433                                    FORT3004
      IF(ICYCLE.LE.2) GO TO 433                                         FORT3005
      IF(DENOM.GE.DELTAC) GO TO 433                                     FORT3006
      ICYCLE=MAXCYC                                                     FORT3007
      PRTCYC=-PRTCYC                                                    FORT3008
433   RETURN                                                            FORT3009
      END                                                               FORT3010
      SUBROUTINE OUTPUT(H,U,MAD,C,E,W,IOCC,HDG,NDIM,NTYPE,NC,NHDG)      FORT3011
C                                                                       FORT3012
C  SUBROUTINE TO ANALYSE AND PRINT OUT RESULTS.                         FORT3013
C                                                                       FORT3014
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3015
      DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),       FORT3016
     1 E(NDIM),W(NDIM),IOCC(NDIM),HDG(NHDG)                             FORT3017
      REAL*8 MAD                                                        FORT3018
      REAL*4 IOCC                                                       FORT3019
      PARAMETER (MAXATM=500)
      PARAMETER (BB=250)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      COMMON/TITLE/AB(10)                                               FORT3020
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,      FORT3021
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                              FORT3022
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                            FORT3023
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)                  FORT3024
      LOGICAL*1 PRT,PUN                                                 FORT3025
      INTEGER*2 IOVPOP,IENRGY                                           FORT3026
      COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),ND(BB),  
     1 EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),COULS(BB), 
     2 COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)                
      INTEGER*2 SYMBOL,KEY,VELEC                                        FORT3030
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)                             FORT3032
      REAL*8 LAMPRI                                                     FORT3033
      INTEGER*4 PRTCYC                                                  FORT3034
      LOGICAL*1 PARTIT,PRINTX,ITABLE                                    FORT3035
      INTEGER*2 SYMB,HYDROG                                             FORT3036
      DATA HYDROG/' H'/                                                 FORT3037
      DO 5 I=1,NDIM                                                     FORT3038
      NMIN=I                                                            FORT3039
      IF(IOCC(I).GT.0.0001) GO TO 7                                     FORT3040
5     CONTINUE                                                          FORT3041
7     IF(ITERAT) GO TO 38                                               FORT3042
C                                                                       FORT3043
C  CALCULATE AND PRINT OUT SUM OF ONE-ELECTRON ENERGIES.                FORT3044
C                                                                       FORT3045
10    SUM=0.0D0                                                         FORT3046
      DO 13 I=1,NDIM                                                    FORT3047
      W(I)=DBLE(IOCC(I))                                                FORT3048
13    SUM=SUM+E(I)*W(I)                                                 FORT3049
      IF(.NOT.PRT(9)) WRITE(6,2001) SUM                                 FORT3050
2001  FORMAT(T10,'SUM OF ONE-ELECTRON ENERGIES =',F16.8,' EV.',///)     FORT3051
      IF(PUN(9)) WRITE(7,2002) SUM                                      FORT3052
2002  FORMAT(F20.8)                                                     FORT3053
      IF(ONEMAT) GO TO 9999                                             FORT3054
C                                                                       FORT3055
C  PRINT OUT WAVE FUNCTIONS.                                            FORT3056
C                                                                       FORT3057
      IF(PRT(10)) GO TO 1003                                            FORT3058
      WRITE(6,1002)                                                     FORT3059
1002  FORMAT('WAVE FUNCTIONS'/'MO''S IN COLUMNS, AO''S IN ROWS')        FORT3060
      CALL PEGLEG(H,NDIM,NDIM)                                          FORT3061
C
C  ** THIS CALL IS TO A SIMILAR ROUTINE TO WRITE THE MO COEFFICIENTS
C     TO DISK FILE 13. THE FORMAT IS DIFFERENT.
C
      CALL PEGLEG2(H,NDIM,NDIM)
C
C  **
C
1003  IF(PUN(10)) WRITE(7,2003) H                                       FORT3062
2003  FORMAT(8F9.6)                                                     FORT3063
C                                                                       FORT3064
C  CALCULATE AND PRINT OUT DENSITY MATRIX.                              FORT3065
C                                                                       FORT3066
500   IF(PRT(11).AND..NOT.PUN(11)) GO TO 38                             FORT3067
      DO 300 I=1,NDIM                                                   FORT3068
      DO 300 J=1,I                                                      FORT3069
      U(I,J)=0.0D0                                                      FORT3070
      DO 310 K=1,NDIM                                                   FORT3071
310   U(I,J)=U(I,J)+H(I,K)*H(J,K)*W(K)                                  FORT3072
300   U(J,I)=U(I,J)                                                     FORT3073
      IF(PRT(11)) GO TO 360                                             FORT3074
      WRITE(6,350)                                                      FORT3075
350   FORMAT('DENSITY MATRIX')                                          FORT3076
      CALL PEGLEG(U,NDIM,NDIM)                                          FORT3077
360   IF(PUN(11)) WRITE(7,2003) U                                       FORT3078
C                                                                       FORT3079
C  CALCULATE ATOMIC ORBITAL OCCUPATIONS AND STORE IN E(I).              FORT3080
C  CALCULATE OVERLAP POPULATION MATRIX.                                 FORT3081
C                                                                       FORT3082
38    IJ=1                                                              FORT3083
      DO 60 I=1,NDIM                                                    FORT3084
      E(I)=0.0D0                                                        FORT3085
      DO 60 J=1,I                                                       FORT3086
      UB=0.0D0                                                          FORT3087
      DO 41 K=NMIN,NDIM                                                 FORT3088
41    UB=UB+H(I,K)*H(J,K)*DBLE(IOCC(K))                                 FORT3089
      UB=UB*0.5D0                                                       FORT3090
      IF(I.EQ.J) GO TO 50                                               FORT3091
      UB=(UB+UB)*C(IJ)                                                  FORT3092
      IJ=IJ+1                                                           FORT3093
50    E(I)=E(I)+UB                                                      FORT3094
      E(J)=E(J)+UB                                                      FORT3095
      IF(ITERAT) GO TO 60                                               FORT3096
      UB=UB+UB                                                          FORT3097
      U(I,J)=UB                                                         FORT3098
      U(J,I)=UB                                                         FORT3099
60    CONTINUE                                                          FORT3100
C                                                                       FORT3101
C  IF DOING CHARGE ITERATION ( METH=1 ) CALL LITER.                     FORT3102
C                                                                       FORT3103
      IF(.NOT.ITERAT) GO TO 80                                          FORT3104
      K=ICYCLE                                                          FORT3105
      PRINTX=((ICYCLE/PRTCYC)*PRTCYC.EQ.ICYCLE)                         FORT3106
      IF(ICYCLE.EQ.MAXCYC) PRINTX=.TRUE.                                FORT3107
      CALL LITER(NH,NA,E,W,J)                                           FORT3108
      IF(ICYCLE.EQ.15000) PRINTX=.TRUE.                                 FORT3109
      IF(K.EQ.1) WRITE(6,355)                                           FORT3110
355   FORMAT('      ATOMIC CHARGES',//)                                 FORT3111
      GO TO (81,82,83,84),J                                             FORT3112
81    WRITE(6,375) K,(X(I),I=1,NATM)                                    FORT3113
375   FORMAT(/,T3,'CYCLE NO.',I3,(T20,10F10.5))                         FORT3114
      GO TO 9000                                                        FORT3115
82    WRITE(6,375) K,(Y(I),I=1,NATM)                                    FORT3116
      GO TO 9000                                                        FORT3117
83    WRITE(6,375) K,(Z(I),I=1,NATM)                                    FORT3118
      GO TO 9000                                                        FORT3119
84    WRITE(6,375) K,(W(I),I=1,NATM)                                    FORT3120
9000  RETURN                                                            FORT3121
C                                                                       FORT3122
C  PRINT OUT OVERLAP POPULATION MATRIX.                                 FORT3123
C                                                                       FORT3124
80    IF(PRT(12)) GO TO 1009                                            FORT3125
      WRITE(6,1006) NELEC                                               FORT3126
1006  FORMAT('OVERLAP POPULATION MATRIX FOR',I4,' ELECTRONS')           FORT3127
      CALL PEGLEG(U,NDIM,NDIM)                                          FORT3128
1009  IF(PUN(12)) WRITE(7,2003) U                                       FORT3129
C                                                                       FORT3130
C  CALL REDUCE TO CALCULATE REDUCED OVERLAP MATRIX. PRINT OUT           FORT3131
C  REDUCED OVERLAP MATRIX.                                              FORT3132
C                                                                       FORT3133
      IF(PRT(13).AND..NOT.PUN(13)) GO TO 2005                           FORT3134
      CALL REDUCE(U,NDIM,NA,NH)                                         FORT3135
      DO 100 I=2,NATM                                                   FORT3136
      K=I-1                                                             FORT3137
      DO 100 J=1,K                                                      FORT3138
100   U(J,I)=U(I,J)                                                     FORT3139
      IF(PRT(13)) GO TO 1015                                            FORT3140
      WRITE(6,1007)                                                     FORT3141
1007  FORMAT('0REDUCED OVERLAP POPULATION MATRIX, ATOM BY ATOM')        FORT3142
      CALL PEGLEG(U,NATM,NDIM)                                          FORT3143
1015  IF(PUN(13)) WRITE(7,2003) ((U(I,J),I=1,NATM),J=1,NATM)            FORT3144
C                                                                       FORT3145
C  IF L3 IS TRUE, CALCULATE AND PRINT OUT OVERLAP POPULATION ANALYSIS,  FORT3146
C  ORBITAL BY ORBITAL, FOR EACH MOLECULAR ORBITAL SPECIFIED.            FORT3147
C                                                                       FORT3148
2005  IF(.NOT.L3) GO TO 21                                              FORT3149
      DO 600 N=1,23,2                                                   FORT3150
      IF(IOVPOP(N).EQ.0) GO TO 25                                       FORT3151
      KMIN=IOVPOP(N)                                                    FORT3152
      KMAX=IOVPOP(N+1)                                                  FORT3153
      DO 600 K=KMIN,KMAX                                                FORT3154
      WRITE(6,2004) K,W(K)                                              FORT3155
2004  FORMAT(///'0OVERLAP POPULATION MATRIX, ORBITAL BY ORBITAL, FOR MOLFORT3156
     1ECULAR ORBITAL',I4,5X,'OCCUPATION IS',F7.4,' ELECTRONS')          FORT3157
      IF(IOCC(K).LT.0.0001) GO TO 600                                   FORT3158
      SUM=0.5D0*W(K)                                                    FORT3159
      IJ=1                                                              FORT3160
      DO 460 I=1,NDIM                                                   FORT3161
      DO 460 J=1,I                                                      FORT3162
      UB=H(I,K)*H(J,K)                                                  FORT3163
      IF(I.EQ.J) GO TO 450                                              FORT3164
      UB=(UB+UB)*C(IJ)                                                  FORT3165
      IJ=IJ+1                                                           FORT3166
450   UB=(UB+UB)*SUM                                                    FORT3167
      U(J,I)=UB                                                         FORT3168
      U(I,J)=UB                                                         FORT3169
460   CONTINUE                                                          FORT3170
      CALL PEGLEG(U,NDIM,NDIM)                                          FORT3171
600   CONTINUE                                                          FORT3172
25    WRITE(6,7003)                                                     FORT3173
7003  FORMAT(///)                                                       FORT3174
C                                                                       FORT3175
C  CALL FULCHM TO CALCULATE COMPLETE CHARGE MATRIX. PRINT OUT           FORT3176
C  COMPLETE CHARGE MATRIX.                                              FORT3177
C                                                                       FORT3178
21    L1=PRT(14).AND..NOT.PUN(14)                                       FORT3179
      L2=PRT(15).AND..NOT.PUN(15)                                       FORT3180
      IF(L1.AND.L2) GO TO 1020                                          FORT3181
      CALL FULCHM(W,U,C,H,NDIM)                                         FORT3182
      IF(PRT(14)) GO TO 2021                                            FORT3183
      WRITE(6,1008)                                                     FORT3184
1008  FORMAT('0COMPLETE CHARGE MATRIX FOR EACH MO, NORMALIZED TO TWO ELEFORT3185
     1CTRONS REGARDLESS OF OCCUPATION')                                 FORT3186
      CALL PEGLEG(U,NDIM,NDIM)                                          FORT3187
2021  IF(PUN(14)) WRITE(7,2003) U                                       FORT3188
C                                                                       FORT3189
C  CALL REDCHM TO CALCULATE REDUCED CHARGE MATRIX. PRINT OUT            FORT3190
C  REDUCED CHARGE MATRIX.                                               FORT3191
C                                                                       FORT3192
      IF(L2) GO TO 1020                                                 FORT3193
      CALL REDCHM(U,NDIM)                                               FORT3194
      IF(PRT(15)) GO TO 1022                                            FORT3195
      WRITE(6,1019)                                                     FORT3196
1019  FORMAT('0REDUCED CHARGE MATRIX, MO''S IN COLUMNS, ATOMS IN ROWS') FORT3197
      CALL OUTMAT(U,NDIM,NATM,NDIM)                                     FORT3198
1022  IF(PUN(15)) WRITE(7,2003) ((U(I,J),I=1,NATM),J=1,NDIM)            FORT3199
C                                                                       FORT3200
C  PRINT OUT ATOMIC CHARGES AND ORBITAL OCCUPATIONS.                    FORT3201
C                                                                       FORT3202
1020  IF(PRT(16)) GO TO 40                                              FORT3203
      WRITE(6,1010)                                                     FORT3204
1010  FORMAT(/'0ATOM',T12,'NET CHG.',T35,'ATOMIC ORBITAL OCCUPATION FOR FORT3205
     *GIVEN MO OCCUPATION'/T35,'S',T45,'X',T55,'Y',T65,'Z',T75,'X2-Y2', FORT3206
     *T85,'Z2',T95,'XY',T105,'XZ',T115,'YZ'/)                           FORT3207
      SYMB=HYDROG                                                       FORT3208
      J=1                                                               FORT3209
      DO 140 I=1,NATM                                                   FORT3210
      IF(I.GT.NH) GO TO 120                                             FORT3211
      UB=1.0-E(I)                                                       FORT3212
      N=I                                                               FORT3213
      GO TO 130                                                         FORT3214
120   KITE=I-NH                                                         FORT3215
      KEYI=KEY(KITE)                                                    FORT3216
      SYMB=SYMBOL(KEYI)                                                 FORT3217
      UB=VELEC(KEYI)-E(J)                                               FORT3218
      N=J                                                               FORT3219
      IF(NP(KEYI).EQ.0) GO TO 130                                       FORT3220
      UB=UB-E(J+1)-E(J+2)-E(J+3)                                        FORT3221
      N=J+3                                                             FORT3222
      IF(ND(KEYI).EQ.0) GO TO 130                                       FORT3223
      UB=UB-E(J+4)-E(J+5)-E(J+6)-E(J+7)-E(J+8)                          FORT3224
      N=J+8                                                             FORT3225
130   WRITE(6,1011) SYMB,I,UB,(E(K),K=J,N)                              FORT3226
1011  FORMAT(1X,A2,I3,T10,F10.5,T30,9F10.5)                             FORT3227
140   J=N+1                                                             FORT3228
C                                                                       FORT3229
C  IF CALCULATING ENERGY MATRIX, PUT DIAGONAL ELEMENTS OF HUCKEL        FORT3230
C  MATRIX IN W(I).                                                      FORT3231
C                                                                       FORT3232
40    PRT(1)=PRT(17).AND..NOT.PUN(17)                                   FORT3233
      PRT(2)=PRT(18).AND..NOT.PUN(18)                                   FORT3234
      PRT(3)=PRT(19).AND..NOT.PUN(19)                                   FORT3235
      PRT(4)=PRT(20).AND..NOT.PUN(20)                                   FORT3236
      L1=PRT(1).AND.PRT(2)                                              FORT3237
      L2=PRT(3).AND.PRT(4)                                              FORT3238
      IF(L1.AND.L2.AND..NOT.L4) GO TO 9999                              FORT3239
      IF(NH.EQ.0) GO TO 23                                              FORT3240
      DO 22 I=1,NH                                                      FORT3241
22    W(I)=X(I)                                                         FORT3242
23    J=NH+1                                                            FORT3243
      K=NH+1                                                            FORT3244
      DO 24 I=1,NA                                                      FORT3245
      KEYI=KEY(I)                                                       FORT3246
      W(J)=X(K)                                                         FORT3247
      J=J+1                                                             FORT3248
      IF(NP(KEYI).EQ.0) GO TO 24                                        FORT3249
      UB=Y(K)                                                           FORT3250
      W(J)=UB                                                           FORT3251
      W(J+1)=UB                                                         FORT3252
      W(J+2)=UB                                                         FORT3253
      J=J+3                                                             FORT3254
      IF(ND(KEYI).EQ.0) GO TO 24                                        FORT3255
      UB=Z(K)                                                           FORT3256
      W(J)=UB                                                           FORT3257
      W(J+1)=UB                                                         FORT3258
      W(J+2)=UB                                                         FORT3259
      W(J+3)=UB                                                         FORT3260
      W(J+4)=UB                                                         FORT3261
      J=J+5                                                             FORT3262
24    K=K+1                                                             FORT3263
C                                                                       FORT3264
C  IF DOING CHARGE ITERATION WITH MADELUNG CORRECTION ON MOLECULE       FORT3265
C  WITH NON-ZERO CHARGE, PUT MADELUNG TERMS IN E(I).                    FORT3266
C                                                                       FORT3267
      QON= DFLOAT(KA)/DFLOAT(NELEC)                                     FORT3268
      IF(METH.LT.3.OR.DABS(QON).LT.0.0001D0) GO TO 180                  FORT3269
      IF(NH.EQ.0) GO TO 181                                             FORT3270
      DO 182 I=1,NH                                                     FORT3271
182   E(I)=MAD(I,I)                                                     FORT3272
181   J=NH+1                                                            FORT3273
      K=NH+1                                                            FORT3274
      DO 183 I=1,NA                                                     FORT3275
      KEYI=KEY(I)                                                       FORT3276
      E(J)=MAD(K,K)                                                     FORT3277
      J=J+1                                                             FORT3278
      K=K+1                                                             FORT3279
      IF(NP(KEYI).EQ.0) GO TO 183                                       FORT3280
      UB=MAD(K,K)                                                       FORT3281
      E(J)=UB                                                           FORT3282
      E(J+1)=UB                                                         FORT3283
      E(J+2)=UB                                                         FORT3284
      J=J+3                                                             FORT3285
      K=K+1                                                             FORT3286
      IF(ND(KEYI).EQ.0) GO TO 183                                       FORT3287
      UB=MAD(K,K)                                                       FORT3288
      E(J)=UB                                                           FORT3289
      E(J+1)=UB                                                         FORT3290
      E(J+2)=UB                                                         FORT3291
      E(J+3)=UB                                                         FORT3292
      E(J+4)=UB                                                         FORT3293
      J=J+5                                                             FORT3294
      K=K+1                                                             FORT3295
183   CONTINUE                                                          FORT3296
C                                                                       FORT3297
C  CALCULATE AND PRINT OUT ENERGY MATRIX.                               FORT3298
C                                                                       FORT3299
180   ONEMAT=.TRUE.                                                     FORT3300
      IF(L1) GO TO 7000                                                 FORT3301
170   SUM=1.0D0                                                         FORT3302
      IF(ONEMAT) SUM=0.0D0                                              FORT3303
      IJ=1                                                              FORT3304
      CNST=2.0D0*CON                                                    FORT3305
      CN2=CNST                                                          FORT3306
      DO 28 I=1,NDIM                                                    FORT3307
      U(I,I)=0.0D0                                                      FORT3308
      DO 28 J=1,I                                                       FORT3309
      UB=0.0D0                                                          FORT3310
      DO 26 K=NMIN,NDIM                                                 FORT3311
26    UB=UB+H(I,K)*H(J,K)*DBLE(IOCC(K))                                 FORT3312
      IF(I.EQ.J) GO TO 28                                               FORT3313
      IF(.NOT.L5) GO TO 35                                              FORT3314
      UC=W(I)                                                           FORT3315
      ET=W(J)                                                           FORT3316
      IF(NHDG.EQ.1) GO TO 36                                            FORT3317
      UC=HDG(I)                                                         FORT3318
      ET=HDG(J)                                                         FORT3319
36    UC=(UC-ET)/(UC+ET)                                                FORT3320
      UC=UC*UC                                                          FORT3321
      CNST=CN2+UC+UC*UC*(1.0D0-CN2)                                     FORT3322
35    IF(METH.LT.3.OR.DABS(QON).LT.0.0001D0) GO TO 31                   FORT3323
      UC=UB*C(IJ)*((CNST-SUM)*(W(I)+W(J))-QON*(1.0D0-CNST)*(E(I)+E(J))) FORT3324
      GO TO 32                                                          FORT3325
31    UC=UB*(CNST-SUM)*C(IJ)*(W(I)+W(J))                                FORT3326
32    UB=SUM*UB*C(IJ)                                                   FORT3327
      IJ=IJ+1                                                           FORT3328
      U(J,I)=UC                                                         FORT3329
      U(I,J)=UC                                                         FORT3330
      U(J,J)=U(J,J)+UB*W(J)                                             FORT3331
28    U(I,I)=U(I,I)+UB*W(I)                                             FORT3332
      IF(PRT(17)) GO TO 3020                                            FORT3333
      IF(ONEMAT) WRITE(6,3005)                                          FORT3334
3005  FORMAT(///,'0ENERGY MATRIX')                                      FORT3335
      IF(.NOT.ONEMAT) WRITE(6,3006)                                     FORT3336
3006  FORMAT('0ENERGY PARTITIONING')                                    FORT3337
      CALL PEGLEG(U,NDIM,NDIM)                                          FORT3338
3020  IF(PUN(17)) WRITE(7,2800) U                                       FORT3339
2800  FORMAT(8F9.5)                                                     FORT3340
C                                                                       FORT3341
C  CALL REDUCE TO CALCULATE REDUCED ENERGY MATRIX. PRINT OUT            FORT3342
C  REDUCED ENERGY MATRIX.                                               FORT3343
C                                                                       FORT3344
      IF(PRT(2)) GO TO 7000                                             FORT3345
      CALL REDUCE(U,NDIM,NA,NH)                                         FORT3346
      DO 700 I=2,NATM                                                   FORT3347
      K=I-1                                                             FORT3348
      DO 700 J=1,K                                                      FORT3349
700   U(J,I)=U(I,J)                                                     FORT3350
      IF(PRT(18)) GO TO 3021                                            FORT3351
      IF(ONEMAT) WRITE(6,3007)                                          FORT3352
3007  FORMAT('0REDUCED ENERGY MATRIX, ATOM BY ATOM')                    FORT3353
      IF(.NOT.ONEMAT) WRITE(6,3008)                                     FORT3354
3008  FORMAT('0REDUCED ENERGY PARTITIONING, ATOM BY ATOM')              FORT3355
      KITE=0                                                            FORT3356
701   LOW=KITE+1                                                        FORT3357
      KITE=KITE+13                                                      FORT3358
      IF(KITE.GT.NATM) KITE=NATM                                        FORT3359
      WRITE(6,702) (I,I=LOW,KITE)                                       FORT3360
702   FORMAT(/5X,13I9,//)                                               FORT3361
      DO 703 I=1,NATM                                                   FORT3362
703   WRITE(6,704) I,(U(I,J),J=LOW,KITE)                                FORT3363
704   FORMAT(I5,2X,13F9.4)                                              FORT3364
      IF(KITE.LT.NATM) GO TO 701                                        FORT3365
3021  IF(PUN(18)) WRITE(7,2111) ((U(I,J),I=1,NATM),J=1,NATM)            FORT3366
2111  FORMAT(7F10.5)                                                    FORT3367
C                                                                       FORT3368
C  IF L4 IS TRUE, CALCULATE AND PRINT OUT ENERGY MATRIX ANALYSIS,       FORT3369
C  ORBITAL BY ORBITAL, FOR EACH MOLECULAR ORBITAL SPECIFIED.            FORT3370
C                                                                       FORT3371
7000  IF(.NOT.ONEMAT) GO TO 9999                                        FORT3372
      IF(.NOT.L4) GO TO 71                                              FORT3373
      DO 246 N=1,23,2                                                   FORT3374
      IF(IENRGY(N).EQ.0) GO TO 72                                       FORT3375
      KMIN=IENRGY(N)                                                    FORT3376
      KMAX=IENRGY(N+1)                                                  FORT3377
      DO 246 K=KMIN,KMAX                                                FORT3378
      WRITE(6,3004) K,IOCC(K)                                           FORT3379
3004  FORMAT(///'0ENERGY MATRIX, ORBITAL BY ORBITAL, FOR MOLECULAR ORBITFORT3380
     1AL',I4,5X,'OCCUPATION IS',F7.4,' ELECTRONS')                      FORT3381
      IF(IOCC(K).LT.0.0001) GO TO 246                                   FORT3382
      IJ=1                                                              FORT3383
      CNST=CON                                                          FORT3384
      EX=DBLE(IOCC(K))                                                  FORT3385
      DO 244 I=1,NDIM                                                   FORT3386
      DO 244 J=1,I                                                      FORT3387
      UB=H(I,K)*H(J,K)                                                  FORT3388
      IF(I.EQ.J) GO TO 242                                              FORT3389
      IF(.NOT.L5) GO TO 236                                             FORT3390
      UC=W(I)                                                           FORT3391
      ET=W(J)                                                           FORT3392
      IF(NHDG.EQ.1) GO TO 237                                           FORT3393
      UC=HDG(I)                                                         FORT3394
      ET=HDG(J)                                                         FORT3395
237   UC=(UC-ET)/(UC+ET)                                                FORT3396
      UC=UC*UC                                                          FORT3397
      CNST=CON+UC/2.0D0+UC*UC*(0.5D0-CON)                               FORT3398
236   IF(METH.LT.3.OR.DABS(QON).LT.0.0001D0) GO TO 247                  FORT3399
      UB=UB*C(IJ)*(CNST*(W(I)+W(J))-QON*(0.5D0-CNST)*(E(I)+E(J)))       FORT3400
      GO TO 248                                                         FORT3401
247   UB=UB*CNST*C(IJ)*(W(I)+W(J))                                      FORT3402
248   IJ=IJ+1                                                           FORT3403
      UB=2.0D0*UB*EX                                                    FORT3404
      U(I,J)=UB                                                         FORT3405
      U(J,I)=UB                                                         FORT3406
      GO TO 244                                                         FORT3407
242   U(I,I)=UB*W(I)*EX                                                 FORT3408
244   CONTINUE                                                          FORT3409
      CALL PEGLEG(U,NDIM,NDIM)                                          FORT3410
246   CONTINUE                                                          FORT3411
72    WRITE(6,7003)                                                     FORT3412
C                                                                       FORT3413
C  CALCULATE AND PRINT OUT ENERGY PARTITIONING AND REDUCED              FORT3414
C  ENERGY PARTITIONING.                                                 FORT3415
C                                                                       FORT3416
71    IF(L2) GO TO 9999                                                 FORT3417
      ONEMAT=.FALSE.                                                    FORT3418
      PRT(17)=PRT(19)                                                   FORT3419
      PUN(17)=PUN(19)                                                   FORT3420
      PRT(18)=PRT(20)                                                   FORT3421
      PUN(18)=PUN(20)                                                   FORT3422
      PRT(2)=PRT(4)                                                     FORT3423
      GO TO 170                                                         FORT3424
9999  RETURN                                                            FORT3425
      END                                                               FORT3426
      SUBROUTINE LITER(NH,NA,E,W,NCYC)                                  FORT3427
C                                                                       FORT3428
C  SUBROUTINE FOR CALCULATING Q*SENSE WHEN USING CHARGE ITERATION       FORT3429
C  OPTION ( METH = 1 ).                                                 FORT3430
C                                                                       FORT3431
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3432
      DIMENSION E(1),W(1)                                               FORT3433
      PARAMETER (MAXATM=500)
      PARAMETER (BB=250)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),ND(BB), 
     1 EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),COULS(BB),
     2 COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)                   
      INTEGER*2 SYMBOL,KEY,VELEC                                        FORT3437
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)                             FORT3439
      REAL*8 LAMPRI                                                     FORT3440
      INTEGER*4 PRTCYC                                                  FORT3441
      LOGICAL*1 PARTIT,PRINTX,ITABLE                                    FORT3442
      NCYC=ICYCLE+1-(ICYCLE/4)*4                                        FORT3443
      ICYCLE=ICYCLE+1                                                   FORT3444
      DELTA=0.D0                                                        FORT3445
      NATOM=NH+NA                                                       FORT3446
      INDEX=1                                                           FORT3447
      DO 60 I=1,NATOM                                                   FORT3448
      INCR=1                                                            FORT3449
      IF(I.GT.NH) GO TO 100                                             FORT3450
      CHG=1.0D0-E(INDEX)                                                FORT3451
      GO TO 150                                                         FORT3452
100   KEYI=KEY(I-NH)                                                    FORT3453
      CHG=E(INDEX)                                                      FORT3454
      IF(ND(KEYI).EQ.0) GO TO 120                                       FORT3455
      INCR=9                                                            FORT3456
      CHG=CHG+E(INDEX+4)+E(INDEX+5)+E(INDEX+6)+E(INDEX+7)+E(INDEX+8)    FORT3457
      GO TO 130                                                         FORT3458
120   IF(NP(KEYI).EQ.0) GO TO 140                                       FORT3459
      INCR=4                                                            FORT3460
130   CHG=CHG+E(INDEX+1)+E(INDEX+2)+E(INDEX+3)                          FORT3461
140   CHG=VELEC(KEYI)-CHG                                               FORT3462
150   INDEX=INDEX+INCR                                                  FORT3463
      GO TO (10,20,30,40),NCYC                                          FORT3464
10    CHG=0.25D0*(CHG+Y(I)+Z(I)+W(I))                                   FORT3465
      X(I)=CHG                                                          FORT3466
      Y(I)=CHG                                                          FORT3467
      Z(I)=CHG                                                          FORT3468
      W(I)=CHG                                                          FORT3469
      GO TO 60                                                          FORT3470
20    DELTA=DELTA+DABS(CHG-X(I))                                        FORT3471
      Y(I)=CHG                                                          FORT3472
      GO TO 60                                                          FORT3473
30    DELTA=DELTA+DABS(CHG-Y(I))                                        FORT3474
      Z(I)=CHG                                                          FORT3475
      GO TO 60                                                          FORT3476
40    DELTA=DELTA+DABS(CHG-Z(I))                                        FORT3477
      W(I)=CHG                                                          FORT3478
60    E(I)=CHG*SENSE                                                    FORT3479
      IF(DELTAC.LT.DELTA.OR.NCYC.EQ.1) RETURN                           FORT3480
      ICYCLE=15000                                                      FORT3481
      RETURN                                                            FORT3482
      END                                                               FORT3483
      SUBROUTINE MATRIX(N,NS1,NS2,NS3,NS4,NS5,NS6,NS7,NS8,NS9,NS10,     FORT3484
     X  NS11,NS12,NS13,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9,ND10,        FORT3485
     X  ND11,ND12,ND13)                                                 FORT3486
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3487
      DIMENSION A(2000000)           
      NSIZE=2000000                                                
C                                                                       FORT3490
C    SUBROUTINE TO ALLOCATE STORAGE FOR MATRICES.                       FORT3491
C                                                                       FORT3492
      CALL ERRSET(74,0,0,0,1,0)                                         FORT3493
      I1=1                                                              FORT3494
      I2=I1+NS1                                                         FORT3495
      I3=I2+NS2                                                         FORT3496
      I4=I3+NS3                                                         FORT3497
      I5=I4+NS4                                                         FORT3498
      I6=I5+NS5                                                         FORT3499
      I7=I6+NS6                                                         FORT3500
      I8=I7+NS7                                                         FORT3501
      I9=I8+NS8                                                         FORT3502
      I10=I9+NS9                                                        FORT3503
      I11=I10+NS10                                                      FORT3504
      I12=I11+NS11                                                      FORT3505
      I13=I12+NS12                                                      FORT3506
      IF=I13+NS13                                                       FORT3507
      IF(IF.GT.NSIZE) GO TO 200                                         FORT3508
      CALL MOVLAP (A(I1),A(I2),A(I3),A(I4),A(I5),A(I6),A(I7),           FORT3509
     X  A(I8),A(I9),A(I10),A(I11),A(I12),A(I13),ND1,ND2,ND3,ND4,        FORT3510
     X  ND5,ND6)                                                        FORT3511
      CALL HUCKEL (A(I1),A(I2),A(I3),A(I4),A(I5),A(I6),A(I7),           FORT3512
     X  A(I8),A(I9),A(I10),A(I11),A(I12),A(I13),ND1,ND2,ND3,ND4,        FORT3513
     X  ND5,ND6)                                                        FORT3514
      RETURN                                                            FORT3515
 200  WRITE(6,1001) IF,NSIZE                                            FORT3516
 1001 FORMAT('0*** INSUFFICIENT SPACE FOR MATRICES'/'0PROGRAM REQUESTS' FORT3517
     1,I7,' DOUBLE WORDS BUT ONLY',I7,' ARE AVAILABLE'/'0RECOMPILE WITH FORT3518
     2LARGER DIMENSION FOR A AND INCREASED VALUE OF NSIZE')             FORT3519
      RETURN                                                            FORT3520
      END                                                               FORT3521
      SUBROUTINE ABFNS(A,B)                                             FORT3522
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3523
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      DIMENSION A(MXUSER),B(MXUSER)                               
      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL                     FORT3525
C                                                                       FORT3526
C    SUBROUTINE FOR CALCULATING A AND B FUNCTIONS FOR USE IN LOVLAP.    FORT3527
C                                                                       FORT3528
      J=MAXCAL+1                                                        FORT3529
      RHO1=0.5D0*(SK1+SK2)*RR                                           FORT3530
      RHO2=0.5D0*(SK1-SK2)*RR                                           FORT3531
      IF(DABS(RHO1).GT.165.D0) GO TO 100                                FORT3532
      IF(DABS(RHO2).GT.165.D0) GO TO 100                                FORT3533
      C=DEXP(-RHO1)                                                     FORT3534
      A(1)=C/RHO1                                                       FORT3535
      DO 15 I=2,J                                                       FORT3536
 15   A(I)=(DFLOAT(I-1)*A(I-1)+C)/RHO1                                  FORT3537
      IX=J                                                              FORT3538
      IR=DABS(2.*RHO2)                                                  FORT3539
      IS=MIN0(IR+1,19)                                                  FORT3540
      IF(RHO2) 25,35,25                                                 FORT3541
 25   D=DEXP(RHO2)                                                      FORT3542
      H=1.D0/D                                                          FORT3543
C                                                                       FORT3544
C    USE THE DSINH ROUTINE INSTEAD OF SUMMING THE INFINITE SERIES.      FORT3545
C                                                                       FORT3546
      R=2.D0*DSINH(RHO2)                                                FORT3547
C                                                                       FORT3548
C    AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE       FORT3549
C    RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.                       FORT3550
C                                                                       FORT3551
      B(1)=R/RHO2                                                       FORT3552
      DO 51 I=2,IX,IS                                                   FORT3553
      IF(IR.EQ.0) GO TO 40                                              FORT3554
      IL=IS-1                                                           FORT3555
C                                                                       FORT3556
C    MODIFICATION TO AVOID EXCEEDING STORAGE LIMITS.                    FORT3557
C    D. WALLACE 04/14/71                                                FORT3558
C                                                                       FORT3559
      DO 31 K=I,IX                                                      FORT3560
      IF((-1)**K) 29,29,30                                              FORT3561
 29   B(K)=(R+DFLOAT(K-1)*B(K-1))/RHO2                                  FORT3562
      GO TO 31                                                          FORT3563
 30   B(K)=-(D+H-DFLOAT(K-1)*B(K-1))/RHO2                               FORT3564
 31   CONTINUE                                                          FORT3565
 40   IN=I+IS-1                                                         FORT3566
C                                                                       FORT3567
C    AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE     FORT3568
C    NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION       FORT3569
C    OF THE INFINITE SERIES.                                            FORT3570
C                                                                       FORT3571
      IF(IN-IX) 39,39,38                                                FORT3572
 39   IF((-1)**IN) 44,44,42                                             FORT3573
 42   TR=RHO2                                                           FORT3574
 105  B(IN)=-2.*TR/DFLOAT(IN+1)                                         FORT3575
      DO 43 J=1,500                                                     FORT3576
      TR=TR*RHO2**2/DFLOAT((2*J)*(2*J+1))                               FORT3577
      IF(DABS(TR/B(IN))-1.0D-7 ) 51,51,43                               FORT3578
 43   B(IN)=B(IN)-2.*TR/DFLOAT(IN+1+2*J)                                FORT3579
      GO TO 51                                                                  
 44   TR=1.                                                             FORT3580
 107  B(IN)=2.*TR/DFLOAT(IN)                                            FORT3581
      DO 46 J=1,500                                                     FORT3582
      TR=TR*RHO2**2/DFLOAT((2*J)*(2*J-1))                               FORT3583
      IF(DABS(TR/B(IN))-1.0D-7 ) 51,51,46                               FORT3584
 46   B(IN)=B(IN)+2.*TR/DFLOAT(IN+2*J)                                  FORT3585
 51   CONTINUE                                                          FORT3586
C                                                                       FORT3587
C    IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.           FORT3588
C                                                                       FORT3589
      GO TO 38                                                          FORT3590
 35   DO 36 I=1,IX,2                                                    FORT3591
      B(I)=2.D0/DFLOAT(I)                                               FORT3592
 36   B(I+1)=0.D0                                                       FORT3593
 38   RETURN                                                            FORT3594
 100  DO 101 I=1,MXUSER                                            
      A(I)=0.D0                                                         FORT3596
 101  B(I)=0.D0                                                         FORT3597
      GO TO 38                                                          FORT3598
      END                                                               FORT3599
      SUBROUTINE LOVLAP (STRAD,A,B)                                     FORT3600
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3601
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      DIMENSION A(MXUSER),B(MXUSER)                                   
      DIMENSION FACT(25)                                                FORT3603
      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX                        FORT3604
      DIMENSION BINCOE(7,7)                                             FORT3605
      DATA BINCOE/7*1.D0,   0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,         FORT3606
     1 2*0.D0,1.D0,3.D0,6.D0,10.D0,15.D0,   3*0.D0,1.D0,4.D0,10.D0,     FORT3607
     2 20.D0,   4*0.D0,1.D0,5.D0,15.D0,   5*0.D0,1.D0,6.D0,   6*0.D0,   FORT3608
     3 1.D0/                                                            FORT3609
      LOGICAL*1 JGO                                                     FORT3610
      DATA JGO/.FALSE./                                                 FORT3611
C                                                                       FORT3612
C    SUBROUTINE TO CALCULATE OVERLAP INTEGRALS IN A LOCAL               FORT3613
C    COORDINATE SYSTEM.                                                 FORT3614
C                                                                       FORT3615
C    INTEGRALS ARE CALCULATED BY TRANSFORMATION TO ELLIPSOIDAL          FORT3616
C    COORDINATES AND THEREBY EXPRESSED IN TERMS OF C-FUNCTIONS.         FORT3617
C    SEE J.C.P.,24,201. ORIGINALLY WRITTEN BY R.M.STEVENS.              FORT3618
C                                                                       FORT3619
C                                                                       FORT3620
C    GENERATE FACTORIALS ONLY ONCE.                                     FORT3621
C                                                                       FORT3622
      IF(JGO) GO TO 10                                                  FORT3623
      JGO=.TRUE.                                                        FORT3624
      FACT(1)=1.D0                                                      FORT3625
      DO 5 I=2,25                                                       FORT3626
 5    FACT(I)=FACT(I-1)*DFLOAT(I-1)                                     FORT3627
 10   CONTINUE                                                          FORT3628
      M2=M1                                                             FORT3629
      STRAD=0.D0                                                        FORT3630
      RHOA=R*SK1                                                        FORT3631
      RHOB=R*SK2                                                        FORT3632
      TERMA=0.5D0**(L1+L2+1) * DSQRT(DFLOAT((L1+L1+1)*(L2+L2+1))*       FORT3633
     1 FACT(L1-M1+1)*FACT(L2-M1+1)/(FACT(N1+N1+1)*FACT(N2+N2+1)*        FORT3634
     2 FACT(L1+M1+1)*FACT(L2+M1+1))*RHOA**(N1+N1+1)*RHOB**(N2+N2+1))    FORT3635
      JEND=1+((L1-M1)/2)                                                FORT3636
      KEND=1+((L2-M2)/2)                                                FORT3637
      IEB=M1+1                                                          FORT3638
      DO 50 J=1,JEND                                                    FORT3639
      JU=J-1                                                            FORT3640
      IAB=N1-L1+JU+JU+1                                                 FORT3641
      ICB=L1-M1-JU-JU+1                                                 FORT3642
      CON1=FACT(L1+L1-JU-JU+1)/(FACT(L1-M1-JU-JU+1)*FACT(JU+1)*         FORT3643
     1 FACT(L1-JU+1))                                                   FORT3644
      DO 50 K=1,KEND                                                    FORT3645
      KU=K-1                                                            FORT3646
      CON12=CON1*FACT(L2+L2-KU-KU+1)/(FACT(L2-M2-KU-KU+1)*FACT(KU+1)*   FORT3647
     1 FACT(L2-KU+1))                                                   FORT3648
      IEV=JU+KU+L2                                                      FORT3649
      IF(2*(IEV/2).NE.IEV) CON12=-CON12                                 FORT3650
      IBB=N2-L2+KU+KU+1                                                 FORT3651
      IDB=L2-M2-KU-KU+1                                                 FORT3652
      VALUE=0.D0                                                        FORT3653
      DO 90 I6=1,IEB                                                    FORT3654
      DO 90 I5=1,IEB                                                    FORT3655
      VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)                              FORT3656
      IEV=I5+I6                                                         FORT3657
      IF(2*(IEV/2).NE.IEV) VALUE1=-VALUE1                               FORT3658
      DO 90 I4=1,IDB                                                    FORT3659
      VALUE1=-VALUE1                                                    FORT3660
      VALUE2=BINCOE(IDB,I4)*VALUE1                                      FORT3661
      DO 90 I3=1,ICB                                                    FORT3662
      VALUE3=BINCOE(ICB,I3)*VALUE2                                      FORT3663
      DO 90 I2=1,IBB                                                    FORT3664
      VALUE3=-VALUE3                                                    FORT3665
      VALUE4=BINCOE(IBB,I2)*VALUE3                                      FORT3666
      DO 90 I1=1,IAB                                                    FORT3667
      TERM=VALUE4*BINCOE(IAB,I1)                                        FORT3668
      IR=I1+I2+IEB+IEB-I6-I6-I3+IDB-I4+ICB-1                            FORT3669
      IP=IAB-I1+IBB-I2+IEB+IEB-I5-I5+ICB-I3+IDB-I4+1                    FORT3670
 90   VALUE=VALUE+A(IP)*B(IR)*TERM                                      FORT3671
 50   STRAD=STRAD+VALUE*CON12                                           FORT3672
      STRAD=STRAD*TERMA                                                 FORT3673
      RETURN                                                            FORT3674
      END                                                               FORT3675
      SUBROUTINE GRMSCH(U,C,NDIM)                                       FORT3676
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3677
      DIMENSION U(NDIM,NDIM),C(NDIM)                                    FORT3678
C                                                                       FORT3679
C    SUBROUTINE TO CARRY OUT A GRAM SCHMIDT ORTHOGONALIZATION ON A      FORT3680
C    SET OF VECTORS X(I). THEY ARE CONVERTED INTO A SET OF ORTHONORMAL  FORT3681
C    VECTORS Y(J).  THE UNITARY TRANSFORMATION IS DEVELOPED.            FORT3682
C                                                                       FORT3683
C    U   CONTAINS THE OVERLAP MATRIX IN THE UPPER RIGHT TRIANGULAR      FORT3684
C        PART. THE LOWER TRIANGULAR PART IS NOT USED. THE UNITARY       FORT3685
C        TRANSFORMATION IS RETURNED IN THE UPPER RIGHT TRIANGLE,        FORT3686
C        INCLUDING THE DIAGONAL.                                        FORT3687
C    C   IS A WORK VECTOR OF LENGTH NDIM.                               FORT3688
C                                                                       FORT3689
C    NDIM IS THE DIMENSION OF C AND U (IE. C(NDIM), U(NDIM,NDIM)).      FORT3690
C                                                                       FORT3691
C    ACTUALLY THE FIRST COLUMN OF U MAY BE USED AS THE WORK AREA        FORT3692
C    PROVIDED THAT THE ELEMENT U(1,1) IS SET EQUAL TO 1.D0 AFTER        FORT3693
C    RETURN IS MADE TO THE CALLING PROGRAM.                             FORT3694
C                                                                       FORT3695
      U(1,1)=1.D0                                                       FORT3696
      DO 100 I=2,NDIM                                                   FORT3697
      I1=I-1                                                            FORT3698
      XNORM=1.D0                                                        FORT3699
      DO 40 J=1,I1                                                      FORT3700
      SUM=0.D0                                                          FORT3701
      DO 30 K=1,J                                                       FORT3702
      SUM=SUM+U(K,J)*U(K,I)                                             FORT3703
 30   CONTINUE                                                          FORT3704
      C(J)=SUM                                                          FORT3705
      XNORM=XNORM-SUM**2                                                FORT3706
 40   CONTINUE                                                          FORT3707
      XNORM=DSQRT(1.D0/XNORM)                                           FORT3708
      DO 70 J=1,I1                                                      FORT3709
      SUM=0.D0                                                          FORT3710
      DO 60 K=J,I1                                                      FORT3711
      SUM=SUM+C(K)*U(J,K)                                               FORT3712
 60   CONTINUE                                                          FORT3713
      U(J,I)=-SUM*XNORM                                                 FORT3714
 70   CONTINUE                                                          FORT3715
      U(I,I)=XNORM                                                      FORT3716
 100  CONTINUE                                                          FORT3717
      RETURN                                                            FORT3718
      END                                                               FORT3719
      SUBROUTINE TRNFRM(S,H,C,COUL0,NDIM,SP,IEXIT)                      FORT3720
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3721
      DIMENSION S(NDIM,NDIM),H(NDIM,NDIM),C(NDIM),COUL0(NDIM),SP(NDIM)  FORT3722
      LOGICAL*1 ONEMAT                                                  FORT3723
C                                                                       FORT3724
C    SUBROUTINE FOR CARRYING OUT A CHANGE IN BASIS SET ON A MATRIX H    FORT3725
C    BY MEANS OF A MATRIX U.                                            FORT3726
C                                                                       FORT3727
C                C=U'*H*U                                               FORT3728
C                                                                       FORT3729
C    H IS ASSUMED REAL SYMMETRIC AND ONLY THE LOWER LEFT TRIANGLE       FORT3730
C    IS REFERENCED. U IS ASSUMED TO BE ENTIRELY CONTAINED IN ITS        FORT3731
C    UPPER RIGHT TRIANGLE AND ONLY THIS PART IS REFERENCED. THE         FORT3732
C    RESULTING REAL SYMMETRIC MATRIX, C, IS STORED IN PACKED FORM.      FORT3733
C                                                                       FORT3734
C    INDEX1 OF LOOP1 POINTS TO U(1,J), J=1,NDIM                         FORT3735
C    INDEX2 OF LOOP2 POINTS TO H(I,1), I=1,J                            FORT3736
C    INDEX3 OF LOOP3 POINTS TO U(K,J), K=1,I                            FORT3737
C    INDEX4 OF LOOP4 POINTS TO U(K,J), K=I+1,J                          FORT3738
C    INDEX5 OF LOOP5 POINTS TO U(1,I), I=1,J                            FORT3739
C    INDEX6 OF LOOP6 POINTS TO U(K,I), K=1,I                            FORT3740
C                                                                       FORT3741
                                                                        FORT3742
      ONEMAT=IEXIT.EQ.2.OR.IEXIT.EQ.3                                   FORT3743
      ISUB=1                                                            FORT3744
      DO 100 J=1,NDIM                                                   FORT3745
      DO 40 I=1,J                                                       FORT3746
      SUM=0.D0                                                          FORT3747
      ILIM=I                                                            FORT3748
      IF(.NOT.ONEMAT) GO TO 10                                          FORT3749
      IF(I.LT.J) GO TO 10                                               FORT3750
      ILIM=I-1                                                          FORT3751
C                                                                       FORT3752
C    IN THE ONE MATRIX CASE, THE S AND H MATRICES OVERLAP ALONG THE     FORT3753
C    DIAGONAL.  THE DIAGONAL OF S IS THEN PUT INTO SP.                  FORT3754
C                                                                       FORT3755
      SUM=SP(I)*H(I,I)                                                  FORT3756
      IF(I.EQ.1) GO TO 25                                               FORT3757
 10   DO 20 K=1,ILIM                                                    FORT3758
 20   SUM=SUM+S(K,J)*H(I,K)                                             FORT3759
 25   IF(I.EQ.J) GO TO 40                                               FORT3760
      I1=I+1                                                            FORT3761
      ILIM=J                                                            FORT3762
      IF(.NOT.ONEMAT) GO TO 27                                          FORT3763
      ILIM=J-1                                                          FORT3764
      SUM=SUM+SP(J)*H(J,I)                                              FORT3765
      IF(ILIM.LT.I1) GO TO 40                                           FORT3766
 27   DO 30 K=I1,ILIM                                                   FORT3767
 30   SUM=SUM+S(K,J)*H(K,I)                                             FORT3768
 40   COUL0(I)=SUM                                                      FORT3769
      DO 60 I=1,J                                                       FORT3770
      SUM=0.D0                                                          FORT3771
      ILIM=I                                                            FORT3772
      IF(.NOT.ONEMAT) GO TO 45                                          FORT3773
      ILIM=I-1                                                          FORT3774
      SUM=COUL0(I)*SP(I)                                                FORT3775
      IF(ILIM.EQ.0) GO TO 55                                            FORT3776
 45   DO 50 K=1,ILIM                                                    FORT3777
 50   SUM=SUM+COUL0(K)*S(K,I)                                           FORT3778
 55   C(ISUB)=SUM                                                       FORT3779
 60   ISUB=ISUB+1                                                       FORT3780
 100  CONTINUE                                                          FORT3781
      RETURN                                                            FORT3782
      END                                                               FORT3783
      DOUBLE PRECISION FUNCTION DSUM(B,A,IP1,LIMIT)                     FORT3784
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3785
      DIMENSION B(1),A(1)                                               FORT3786
C                                                                       FORT3787
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.                 FORT3788
C                                                                       FORT3789
      JJ=1                                                              FORT3790
      DSUM=0.D0                                                         FORT3791
      DO 180 II=IP1,LIMIT                                               FORT3792
      DSUM=DSUM+B(II+1)*A(JJ)                                           FORT3793
 180  JJ=JJ+II                                                          FORT3794
      RETURN                                                            FORT3795
      END                                                               FORT3796
      SUBROUTINE ROTATE(V,C,S,NJX,JTOP)                                 FORT3797
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3798
      DIMENSION V(1)                                                    FORT3799
C                                                                       FORT3800
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.                 FORT3801
C                                                                       FORT3802
      JLIM=JTOP+1                                                       FORT3803
      DO 10 J=1,JLIM,NJX                                                FORT3804
      TA=V(J)                                                           FORT3805
      TB=V(J+1)                                                         FORT3806
      V(J)=TA*C+TB*S                                                    FORT3807
 10   V(J+1)=TB*C-TA*S                                                  FORT3808
      RETURN                                                            FORT3809
      END                                                               FORT3810
      DOUBLE PRECISION FUNCTION DOT(A,B)                                FORT3811
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3812
      DIMENSION A(1),B(1)                                               FORT3813
      COMMON /VECTOR/ FACTOR,LIMIT                                      FORT3814
C                                                                       FORT3815
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.                 FORT3816
C                                                                       FORT3817
      ITOP=LIMIT+1                                                      FORT3818
      DOT=0.D0                                                          FORT3819
      DO 10 I=1,ITOP                                                    FORT3820
 10   DOT=DOT+A(I)*B(I)                                                 FORT3821
      RETURN                                                            FORT3822
      END                                                               FORT3823
      SUBROUTINE VECSUM(A,B)                                            FORT3824
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3825
      DIMENSION A(1),B(1)                                               FORT3826
      COMMON /VECTOR/ FACTOR,LIMIT                                      FORT3827
C                                                                       FORT3828
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.                 FORT3829
C                                                                       FORT3830
      ITOP=LIMIT+1                                                      FORT3831
      DO 10 I=1,ITOP                                                    FORT3832
 10   A(I)=A(I)+FACTOR*B(I)                                             FORT3833
      RETURN                                                            FORT3834
      END                                                               FORT3835
      SUBROUTINE REDUCE(U,NDIM)                                         FORT3836
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3837
      DIMENSION U(NDIM,NDIM)                                            FORT3838
      PARAMETER (MAXATM=500)
      PARAMETER (BB=250)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,      FORT3839
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                              FORT3840
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                            FORT3841
      COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),ND(BB),
     1 EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),       
     2 COULS(BB),COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)           
      INTEGER*2 SYMBOL,KEY,VELEC                                        FORT3845
C                                                                       FORT3846
C    SUBROUTINE TO CALCULATE REDUCED MATRIX.                            FORT3847
C                                                                       FORT3848
      CALL REDCHM(U,NDIM)                                               FORT3849
      ISUB=NH+1                                                         FORT3850
      IAO=ISUB                                                          FORT3851
      DO 100 I=1,NA                                                     FORT3852
      KEYI=KEY(I)                                                       FORT3853
      IF(ND(KEYI)) 10,20,10                                             FORT3854
 10   IATOP=IAO+8                                                       FORT3855
      GO TO 50                                                          FORT3856
 20   IF(NP(KEYI)) 30,40,30                                             FORT3857
 30   IATOP=IAO+3                                                       FORT3858
      GO TO 50                                                          FORT3859
 40   IATOP=IAO                                                         FORT3860
 50   DO 70 J=1,NATM                                                    FORT3861
      SUM=0.D0                                                          FORT3862
      DO 60 K=IAO,IATOP                                                 FORT3863
      SUM=SUM+U(J,K)                                                    FORT3864
 60   CONTINUE                                                          FORT3865
      U(J,ISUB)=SUM                                                     FORT3866
 70   CONTINUE                                                          FORT3867
      ISUB=ISUB+1                                                       FORT3868
 100  IAO=IATOP+1                                                       FORT3869
      RETURN                                                            FORT3870
      END                                                               FORT3871
      SUBROUTINE FULCHM(E,U,C,H,N5)                                     FORT3872
      REAL*8 E(N5),C(N5),U(N5,N5),H(N5,N5)                              FORT3873
      LOGICAL*1 JGO                                                     FORT3874
C                                                                       FORT3875
C    SUBROUTINE TO CALCULATE COMPLETE CHARGE MATRIX.                    FORT3876
C                                                                       FORT3877
      KJ=1                                                              FORT3878
      DO 31 I=1,N5                                                      FORT3879
      JGO=.FALSE.                                                       FORT3880
      IJ=KJ                                                             FORT3881
      DO 32 K=1,N5                                                      FORT3882
      IF(JGO) GO TO 34                                                  FORT3883
 33   IF(I.NE.K) GO TO 35                                               FORT3884
      JGO=.TRUE.                                                        FORT3885
      E(K)=1.0D0                                                        FORT3886
      GO TO 32                                                          FORT3887
 35   E(K)=C(IJ)                                                        FORT3888
      IJ=IJ+1                                                           FORT3889
      GO TO 32                                                          FORT3890
 34   IJ=IJ+K-2                                                         FORT3891
      E(K)=C(IJ)                                                        FORT3892
 32   CONTINUE                                                          FORT3893
      KJ=KJ+I-1                                                         FORT3894
      DO 31 J=1,N5                                                      FORT3895
      UB=0.0D0                                                          FORT3896
      DO 36 K=1,N5                                                      FORT3897
 36   UB=UB+H(K,J)*E(K)                                                 FORT3898
 31   U(I,J)=2.0D0*H(I,J)*UB                                            FORT3899
      RETURN                                                            FORT3900
      END                                                               FORT3901
      SUBROUTINE REDCHM(U,N5)                                           FORT3902
      IMPLICIT REAL*8(A-H,O-Z)                                          FORT3903
      DIMENSION U(N5,N5)                                                FORT3904
      PARAMETER (MAXATM=500)
      PARAMETER (BB=250)
      PARAMETER (MXUSER=230)
      PARAMETER (MXUSR2=231)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,      FORT3905
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT                              FORT3906
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT                            FORT3907
      COMMON/ATOM/KEY(BB),SYMBOL(BB),VELEC(BB),NS(BB),NP(BB),ND(BB),
     . EXPS(BB),EXPP(BB),EXPD(BB),EXPD2(BB),C1(BB),C2(BB),COULS(BB),
     . COULP(BB),COULD(BB),X(MAXATM),Y(MAXATM),Z(MAXATM)
      INTEGER*2 SYMBOL,KEY,VELEC                                        FORT3911
C                                                                       FORT3912
C    SUBROUTINE TO CALCULATE REDUCED CHARGE MATRIX.                     FORT3913
C                                                                       FORT3914
      IF(NA.EQ.0) RETURN                                                FORT3915
      NH1=NH+1                                                          FORT3916
      IAO=NH1                                                           FORT3917
      ISUB=NH1                                                          FORT3918
      DO 128 I=1,NA                                                     FORT3919
      KEYI=KEY(I)                                                       FORT3920
      IF(NP(KEYI)) 122,121,122                                          FORT3921
 121  IATOP=IAO                                                         FORT3922
      GO TO 125                                                         FORT3923
 122  IF(ND(KEYI)) 124,123,124                                          FORT3924
 123  IATOP=IAO+3                                                       FORT3925
      GO TO 125                                                         FORT3926
 124  IATOP=IAO+8                                                       FORT3927
 125  DO 127 J=1,N5                                                     FORT3928
      SUM=0.D0                                                          FORT3929
      DO 126 K=IAO,IATOP                                                FORT3930
 126  SUM=SUM+U(K,J)                                                    FORT3931
 127  U(ISUB,J)=SUM                                                     FORT3932
      ISUB=ISUB+1                                                       FORT3933
 128  IAO=IATOP+1                                                       FORT3934
      RETURN                                                            FORT3935
      END                                                               FORT3936
      SUBROUTINE CORECT(A)                                              FORT3937
      RETURN                                                            FORT3938
      END                                                                       
