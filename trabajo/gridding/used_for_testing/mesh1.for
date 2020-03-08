C*********************************************************************
C
      SUBROUTINE SRV(NCS,IMNCS,NCN,IMNCN,NCH,IMNCH,NN,IMNN,NUM1,
     *               IMNUM1, NUM2  , IMNUM2, NUM3  , IMNUM3,
     *               NP123 , MNP123, PXY   , IMPXY ,
     *               POLAR , MPOLAR, INDEX3, MNDEX3, JWKA  ,
     *               IMJWKA, NG    , IMNG  ,
     *               CORD  , IMCORD, IREGO)
C
C*********************************************************************
C
C ------------------------------------------------
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 NP(3000),N(3000)
C
      INTEGER *4 TEXT1(2),TEXT2(2),TEXT3(2),TEXT4(2),TEXT5(2)
      DIMENSION NCS(2,IMNCS),NCN(2,IMNCN),NCH(2,IMNCH),NN(4,IMNN)
      DIMENSION NUM1(IMNUM1),NUM2(IMNUM2),NUM3(IMNUM3),NP123(3,MNP123)
      DIMENSION PXY(2,IMPXY),POLAR(2,MPOLAR)
      DIMENSION INDEX3(MNDEX3),JWKA(IMJWKA),NG(IMNG)
      DIMENSION CORD(2,IMCORD),NAME(20)
      EQUIVALENCE (N(1),NP(1))
C
C ------------------------------------------------
      COMMON /DAT/ EPS
      COMMON/LPUNIT/LP0,IOUT,IO,IN
C ------------------------------------------------
      NAMELIST / SET /IPRINT,NOPTIM,JK,NE,NUSER,NREG,NDL,MATUSE
      NAMELIST/REGION/I,MATIAL
      NAMELIST/   SYS/NP
      NAMELIST/ POINT/N
      NAMELIST/   SEG/N
C ------------------------------------------------
      DATA           NBN/3/
      DATA           MPOL/50/,MPT/450/,MSEG/300/,MRFL/20/,MCH/100/
      DATA           INP/153/,IND/1000/,INNM/100/,NHM/4000/
      DATA           TEXT1/4HINPU,2HT /,TEXT2/4HMIDD,4HLE P/
      DATA           TEXT3/4HREFL,4HACT /,TEXT4/4HSUPP,4HER E/
      DATA           TEXT5/4HSUPP,4HER C/
      DATA           PI/3.141592653589793D0/
C ------------------------------------------------
      EPS    = 1.E-6
      NE     = 400
      MPOL3  = MPOL*3+3
      NDL=8
      NREG=1
      NUSER  = 0
      MATUSE = 1
      MPT5   = MPT*4
      MSEG6  = MSEG*6
      MRFL2  = MRFL*2+2
      MCH2   = MCH*2+2
C
      REWIND IN
      READ(IN,600 ) NAME
C
C**********************************************************C
C         GENERAL INFORMATION FOR THIS GRID                C
C**********************************************************C
C
      READ(IN,SET)
 3000 FORMAT(4I5,F10.5)
      WRITE(LP0,201) NAME
C
      REWIND IO
      WRITE(IO) NAME
      WRITE(IO) IPRINT,NOPTIM,JK,NE,NUSER,MATUSE,NREG,NDL
C
C**********************************************************C
C         GLOBAL CYCLE FOR REGIONS                         C
C**********************************************************C
C
      IREGO = 1
      MATIAL = 0
C
      DO 2 KREG = 1 , NREG
C
C**********************************************************C
C         GENERAL INFORMATION FOR THIS REGION              C
C**********************************************************C
C
      READ(IN,REGION)
C
C**********************************************************C
C
C         POLAR SYSTEMS
C
C**********************************************************C
C
      IREGF=I
      IREGO = MAX0(IREGO,IREGF)
C
      IF (I.EQ.0)                                GOTO 1000
C
      IF(IPRINT.GE.2)       WRITE(LP0,223) MATIAL,IREGF
      IF(IPRINT.GE.1)       WRITE(LP0,200)TEXT1
C..........................................................C
C     READ POLAR  COORDINATE   SYSTEMS                     C
C..........................................................C
      DO 5 I = 1,MPOL3
   5  NP(I)  = 0.D0
      DO 6 I = 1,MPOL
   6  NUM1(I) = 0
      READ(IN,SYS)
      IPOL    = 0
      IPOLX   = 0
      DO 10 I = 1,INP,3
      J       =IFIX(SNGL(NP(I)))
      IF (J.GT.MPOL)GOTO 700
      IF (J.LT.0   )GOTO 701
      IF (J.EQ.0   )GOTO  20
      IPOL    = IPOL+1
      IF (IPOL.GT.MPOL)GOTO 700
      IPOLX   = MAX0(J,IPOLX)
      NUM1(J) = IPOL
      POLAR(1,IPOL) = NP(I+1)
  10  POLAR(2,IPOL) = NP(I+2)
  20  IF ( IPRINT.GE.1)
     *  CALL LP PS(IPOL,IPOLX,NUM1,POLAR,LP0)
C..........................................................C
C   READ IFORMATION ABOUT POINTS                           C
C..........................................................C
C
      DO 21 I       = 1,MPT5
  21  N(I) = 0.
      DO 22 I = 1,MPT
      NCS(1,I)=0
      NCS(2,I)=0
  22  NUM2(I) = 0
      READ(IN,POINT)
      IPT     = 0
      IPTX    = 0
      DO 30 I = 1,IND,4
      J       =IFIX(SNGL(N(I)))
      IF (J.GT.MPT)GOTO 702
      IF (J.LT.0  )GOTO 703
      IF (J.EQ. 0 )GOTO 40
      IPT     = IPT+1
      IF (IPT.GT.MPT)GOTO 702
      IPTX    = MAX0(J,IPTX)
      NUM2(J) = IPT
      PXY(1,IPT) =N(I+1)
      PXY(2,IPT) =N(I+2)
      NCS(1,IPT)  = 0
      NCS(2,IPT)  =IFIX(SNGL(N(I+3)))
      IF(NCS(2,IPT).NE.0)PXY(2,IPT)=PXY(2,IPT)/180.
  30  CONTINUE
  40  IF (IPT.LT.4)GOTO 704
      IF (IPRINT.GE.1)
     *   CALL LP PT(IPT,IPTX,NUM2,PXY,NCS,LP0)
C..........................................................C
C    READ INFORMATION ABOUT SEGMENTS                       C
C..........................................................C
C
      DO 41 I    = 1,MSEG6
  41  N(I)       = 0.
      DO 42 I    = 1,MSEG
      NCN(1,I)=0
      NCN(2,I)=0
  42  NUM3(I)    = 0
      READ(IN,SEG)
      ISEG       = 0
      ISEGX      = 0
      DO 50 I    = 1,IND,6
      J =IFIX(SNGL(N(I)))
      IF (J.GT.MSEG   )GOTO 705
      IF (J.LT.0      )GOTO 706
      IF (J.EQ.0      )GOTO  60
      ISEG = ISEG+1
      IF (ISEG.GT.MSEG)GOTO 705
      ISEGX = MAX0(J,ISEGX)
      NUM3(J) = ISEG
      NP123(1,ISEG) =IFIX(SNGL(N(I+1)))
      NP123(2,ISEG) =IFIX(SNGL(N(I+2)))
      NP123(3,ISEG) =IFIX(SNGL(N(I+3)))
      NCN(1,ISEG)   =IFIX(SNGL(N(I+4)))
  50  NCN(2,ISEG)   =IFIX(SNGL(N(I+5)))
  60  IF (ISEG.LT.4    )GOTO 707
      IF( IPRINT.GE.1)
     *     CALL LP SG(ISEG,ISEGX,NUM3,NP123,NCN,LP0)
C
C..........................................................C
C  SOME CHECKS OF INPUT DATA                               C
C..........................................................C
C
      CALL CHECK1(ISEG,ISEGX,MPT,MPOL,NUM1,NUM2,NUM3,
     *                                   NP123,NCN,NCS)
C
C..........................................................C
C  DO MIDDLE POINT FOR SOME  SEGMENTS                      C
C..........................................................C
C
      CALL MIDLEP(ISEG,MSEG,MPT,IPT,IPTX,MPOL,NUM1,NUM2,PXY,
     *            NCS,NCN,NP123,POLAR)
      IF(IPRINT.GE.2)WRITE(LP0,200)TEXT2
      IF(IPRINT.GE.2)
     * CALL LP PT(IPT,IPTX,NUM2,PXY,NCS,LP0)
      IF(IPRINT.GE.2)
     * CALL LP SG(ISEG,ISEGX,NUM3,NP123,NCN,LP0)
C
C
C..........................................................C
C     DO SUPERELEMENTS FROM SEGMENTS                       C
C..........................................................C
C
      CALL CSELEM(INN,INNM,ISEG,NP123,NN,INDEX3)
      IF(IPRINT.GE.4) WRITE(LP0,200) TEXT4
      IF(IPRINT.GE.4) CALL LP SUP(INN,NN,LP0)
      CALL CHSIDE(INN,ISEG,NN,NP123)
      IF(IPRINT.GE.1) WRITE(LP0,200) TEXT5
      IF(IPRINT.GE.1) CALL LP SUP(INN,NN,LP0)
C.................................................C
C   FREE NUM1,NUM2,NUM3 AND CHANGE CS,P123,CN     C
C   DO COMPACT INFORMATION                        C
C.................................................C
      DO 100 I1=1,IPT
      J1=NCS(2,I1)
      IF(J1.EQ.0)GOTO 100
      NCS(2,I1)=NUM1(J1)
 100  CONTINUE
      DO 105 I1=1,ISEG
C
      JN=NCN(2,I1)
      IF(NUSER .EQ. 0 )   JN=NCN(2,I1)/2*2+1
C
      NCN(2,I1)=JN
      DO 105 K1=1,3
      J1=NP123(K1,I1)
 105  NP123(K1,I1)=NUM2(J1)
      DO 120 I1=1,INN
      JWKA(I1)=1
      DO 115 J1=1,4
      JKK=NN(J1,I1)
      IF(JK.LT.0)   JKK=-JKK
      JP1=NP123(1,JKK)
      JP2=NP123(2,JKK)
      JP3=NP123(3,JKK)
      JSP1=NCS(2,JP1)
      JSP2=NCS(2,JP2)
      JSP3=NCS(2,JP3)
      IF(JSP1.NE.0.AND.JSP3.NE.0.AND.JSP1.EQ.JSP3) GOTO 115
      JWKA(I1)=0
 115  CONTINUE
 120  CONTINUE
C.................................................C
C DO ALL POINTS ON SEGMENTS                       C
C.................................................C
      NH = 0
C                                                 C
      CALL MPNTS(IPOL,IPT,ISEG,POLAR,PXY,NCS,NP123,NCN,
     *    NHM,NH,CORD,NG,NBN)
      IF(IPRINT.GE.5)CALL LP CORD(NH,ISEG,NP123,CORD,NG,LP0)
      CALL WSCLOK(NH,INN,ISEG,CORD,NN,NP123)
C************************************************************C
C    OUTPUT ON MT0                                           C
C************************************************************C
      IF(IO.LE.0)RETURN
      NSTEP=0
  121 CONTINUE
C
      WRITE(IO)NH,ISEG,INN,NBN
      WRITE(IO)((CORD(I,J),I=1,2),J=1,NH)
      WRITE(IO)(NG(I),I=1,NH)
      WRITE(IO)((NP123(I,J),I=1,3),J=1,ISEG)
      WRITE(IO)((NN(I,J),I=1,4),J=1,INN)
      WRITE(IO)(JWKA(I),I=1,INN)
      WRITE(IO) MATIAL
      IF(IPRINT.GE.6)       WRITE(LP0,335) MATIAL
      WRITE(LP0,202)IREGF
C
    2 CONTINUE
C
 1000 CONTINUE
C
      RETURN
 700   CALL LP PS(IPOL,IPOLX,NUM1,POLAR,LP0)
       WRITE(LP0,500)
       STOP
 701   CALL LP PS(IPOL,IPOLX,NUM1,POLAR,LP0)
       WRITE(LP0,501)
       STOP
 702   CALL LP PT(IPT,IPTX,NUM2,PXY,NCS,LP0)
       WRITE(LP0,502)
       STOP
 703   CALL LP PT(IPT,IPTX,NUM2,PXY,NCS,LP0)
       WRITE(LP0,503)
       STOP
 704   CALL LP PT(IPT,IPTX,NUM2,PXY,NCS,LP0)
       WRITE(LP0,504)
       STOP
 705   CALL LP SG(ISEG,ISEGX,NUM3,NP123,NCN,LP0)
       WRITE(LP0,505)
       STOP
 706   CALL LP SG(ISEG,ISEGX,NUM3,NP123,NCN,LP0)
       WRITE(LP0,506)
       STOP
 707   CALL LP SG(ISEG,ISEGX,NUM3,NP123,NCN,LP0)
       WRITE(LP0,507)
       STOP
 709   WRITE(LP0,509)
       STOP
 710   WRITE(LP0,510)
       STOP
 712   WRITE(LP0,512)
       STOP
 713   WRITE(LP0,513)
       STOP
 714   WRITE(LP0,514)
       STOP
 202  FORMAT(/3X,'*** END OF OUTPUT DATA OF REGION',I4,' OK  ***'/)
 201  FORMAT(/10X,'START OF SERV-EC EXECUTION '
     *        /10X,'PROBLEM : ',20A4/)
 200  FORMAT(/11X,'==== ',2A4,' ===='/)
 223  FORMAT(/5X,'*** MATERIAL      ',I5, ' IN REGION',I5,' ***')
 224  FORMAT(/5X,'*** POLAR SYSTEMS ',I5, ' IN REGION',I5,' ***')
 335  FORMAT ('**MATERIAL IN SUPERELEMENTS OF REGION IS: ',I3)
 500  FORMAT('  ERROR(I1):',
     *       '  SORRY] TOO MANY POLAR SYSTEMS OR TOO BIG NUMER       ')
 501  FORMAT('  ERROR(I2): POLAR SYSTEM HAVE NEGATIVE NUMER          ')
 502  FORMAT('  ERROR(I3):',
     *       '  SORRY] TOO MANY POINTS OR TOO BIG POINT NUMER        ')
 503  FORMAT('  ERROR(I4): POINT HAVE  NEGATIVE NUMER                ')
 504  FORMAT('  ERROR(I5): TOTAL NUMBER OF POINTS LESS THEN 4        ')
 505  FORMAT('  ERROR(I6):',
     *       '  SORRY] TOO MANY SEGMENTS OR TOO BIG SEGMENT NUMER    ')
 506  FORMAT('  ERROR(I7): SEGMENT HAVE NEGATIVE NUMER               ')
 507  FORMAT('  ERROR(I8): TOTAL NUMBER OF SEGMENTS LESS THEN 4      ')
 508  FORMAT('  ERROR(I9): INCORRECT NUMBER OF REFLECTION(<0 OR >100)')
 509  FORMAT('  ERROR(I10): ILLEGAL SEGMENT NUMER DURING REFLACTION  ')
 510  FORMAT('  ERROR(I11): ILLEGAL NODE NUMER IN REFLACTION LINE    ')
 511  FORMAT('  ERROR(I11): INCORRECT NUMBER OF SHIFT (<0 OR >100)   ')
 512  FORMAT('  ERROR(I12): ILLEGAL SEGMENT NUMER DURING SHIFT       ')
 513  FORMAT('  ERROR(I13): ILLEGAL NODE NUMER IN SHIFT      LINE    ')
 514  FORMAT('  ERROR(I14): IRREGULAR SHIFT: ZERO SHIFT              ')
 600  FORMAT(20A4)
      END
C
C****************************************************************
C
      SUBROUTINE CHECK1(ISEG,ISEGX,MPT,MPOL,NUM1,NUM2,NUM3,
     *                                  NP123,NCN,NCS)
C
C****************************************************************
C
C ------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION        NUM1(MPOL),NUM2(MPT),NP123(3,ISEG),NCN(2,ISEG)
      DIMENSION        NCS(2,MPT),NUM3(ISEGX)
C ------------------------------------------------
      COMMON /DAT/ EPS
      COMMON/LPUNIT/LP0,IOUT,MT,IN
C ------------------------------------------------
      JERR = 0
      DO 15 II1 = 1,ISEGX
      I1   = NUM3(II1)
      IF(I1.EQ.0)  GOTO 15
      K1   = NP123(1,I1)
      K2   = NP123(2,I1)
      K3   = NP123(3,I1)
      KK   = (K1-K2)*(K1-K3)*(K2-K3)
      IF(KK.NE.0) GOTO 10
      JERR = -1
      WRITE(LP0,505)II1
  10  CONTINUE
C
C       FIRST POINT
C
      JFL=1
      J1       = NP123(1,I1)
      IF (J1.LE.0.OR.J1.GT.MPT )GOTO 100
      IF (NUM2(J1).EQ.0        )GOTO 101
      IF (NCN(2,I1).LE.0        )GOTO 102
      JJ1      = NUM2(J1)
      IS       = NCS(2,JJ1)
      IF (IS.LT.0.OR.IS.GT.MPOL)GOTO 103
      IF (IS.EQ.0              )GOTO  11
      IF (NUM1(IS).EQ.0        )GOTO 104
      GOTO 11
 100  JERR=-1
      WRITE(LP0,500)II1,JFL
      GOTO 11
 101  JERR=-1
      WRITE(LP0,501)II1,JFL
      GOTO 11
 102  JERR=-1
      WRITE(LP0,502)II1,NCN(2,I1)
      GOTO 11
 103  JERR=-1
      WRITE(LP0,503)II1,JFL
      GOTO 11
 104  JERR=-1
      WRITE(LP0,504)II1,JFL
  11  CONTINUE
C
C       THIRD POINT
C
      JFL=3
      J1       = NP123(3,I1)
      IF (J1.LE.0.OR.J1.GT.MPT )GOTO 300
      IF (NUM2(J1).EQ.0        )GOTO 301
      JJ1      = NUM2(J1)
      IS       = NCS(2,JJ1)
      IF (IS.LT.0.OR.IS.GT.MPOL)GOTO 303
      IF (IS.EQ.0              )GOTO  12
      IF (NUM1(IS).EQ.0        )GOTO 304
      GOTO 12
 300  JERR=-1
      WRITE(LP0,500)II1,JFL
      GOTO 12
 301  JERR=-1
      WRITE(LP0,501)II1,JFL
      GOTO 12
 303  JERR=-1
      WRITE(LP0,503)II1,JFL
      GOTO 12
 304  JERR=-1
      WRITE(LP0,504)II1,JFL
  12  CONTINUE
C
C        SECOND POINT
C
      JFL=2
      J1       = NP123(2,I1)
      IF (J1.EQ.0              )GOTO 15
      IF (J1.LT.0.OR.J1.GT.MPT )GOTO 200
      IF (NUM2(J1).EQ.0        )GOTO 201
      JJ1      = NUM2(J1)
      IS       = NCS(2,JJ1)
      IF (IS.LT.0.OR.IS.GT.MPOL)GOTO 203
      IF (IS.EQ.0              )GOTO 15
      IF (NUM1(IS).EQ.0        )GOTO 204
      GOTO 15
 200  JERR=-1
      WRITE(LP0,500)II1,JFL
      GOTO 15
 201  JERR=-1
      WRITE(LP0,501)II1,JFL
      GOTO 15
 203  JERR=-1
      WRITE(LP0,503)II1,JFL
      GOTO 15
 204  JERR=-1
      WRITE(LP0,504)II1,JFL
 15   CONTINUE
      IF(JERR.LT.0)STOP
      RETURN
 500  FORMAT(//'  ERROR(C1): ILLEGAL NODE NUMER IN SEGMENTS       ',
     *      '  ====== SEGMENT ',I4,' == NODE  ',I4,' ======'//      )
 501  FORMAT(//'  ERROR(C2): NOT EXISTING NODE NUMER IN SEGMENTS  ',
     *      '  ====== SEGMENT ',I4,' == NODE  ',I4,' ======'//        )
 502  FORMAT(//'  ERROR(C3): ILLEGAL SUBDIVISION NUMBER OF SEGMENT',
     *      '  ====== SEGMENT ',I4,' == NDEV  ',I4,' ======'//        )
 503  FORMAT(//'  ERROR(C4): ILLEGAL NUMER COORDINAT SYSTEM OF NODE',
     *      '  ====== SEGMENT ',I4,' == NODE  ',I4,' ======'//        )
 504  FORMAT(//'  ERROR(C5): ILLEGAL NUMER COORDINAT SISTEM OF NODE',
     *      '  ====== SEGMENT ',I4,' == NODE  ',I4,' ======'//        )
 505  FORMAT(//'  ERROR(C6):',
     *                    ' ILLEGAL NODES COMBINATION IN SEGMENT',I4//)
      END
C
C****************************************************************
C
      SUBROUTINE MIDLEP(ISEG,MSEG,MPT,IPT,IPTX,MPOL,NUM1,NUM2,
     *                  PXY,NCS,NCN,NP123,POLAR)
C
C****************************************************************
C
C ------------------------------------------------
      IMPLICIT REAL *8(A-H,O-Z)
C ------------------------------------------------
      DIMENSION NUM1(MPOL),NUM2(MPT),NCS(2,MPT),NCN(2,MSEG),
     * NP123(3,MSEG)
C ------------------------------------------------
      COMMON /DAT/ EPS
      COMMON/LPUNIT/LP0,IOUT,MT,IN
      DIMENSION PXY(2,MPT),POLAR(2,MPOL)
C ------------------------------------------------
      DATA           PI/3.14159 26535/
C ------------------------------------------------
         PI2      = PI*2
      DO 15 I1 = 1,ISEG
         J2       = NP123(2,I1)
      IF (J2.GT.0)GOTO 15
         J1       = NP123(1,I1)
         J3       = NP123(3,I1)
         JJ1      = NUM2(J1)
         JJ3      = NUM2(J3)
         IC1      = NCS(2,JJ1)
         IC3      = NCS(2,JJ3)
         II       = IC1+IC3
      IF (IC1.EQ.IC3.AND.II.NE.0)GOTO 12
  10     X1       = PXY(1,JJ1)
         Y1       = PXY(2,JJ1)
         V1       = X1
         V2       = Y1
      IF (IC1.NE.0)V1 = X1 * DCOS(Y1*PI)+POLAR(1,NUM1(IC1))
      IF (IC1.NE.0)V2 = X1 * DSIN(Y1*PI)+POLAR(2,NUM1(IC1))
         X1 = V1
         Y1 = V2
         X3 = PXY(1,JJ3)
         Y3 = PXY(2,JJ3)
         V1 = X3
         V2 = Y3
      IF (IC3.NE.0)V1 = X3*DCOS(Y3*PI)+POLAR(1,NUM1(IC3))
      IF (IC3.NE.0)V2 = X3*DSIN(Y3*PI)+POLAR(2,NUM1(IC3))
      X3 = V1
      Y3 = V2
      X2 = (X1+X3)/2
      Y2 = (Y1+Y3)/2
      IPT = IPT+1
      IPTX = IPTX+1
      IF (IPT.GT.MPT.OR.IPTX.GT.MPT)GOTO 700
      NUM2(IPTX) = IPT
      NP123(2,I1) = IPTX
      PXY(1,IPT) = X2
      PXY(2,IPT) = Y2
      NCS(1,IPT)  = NCN(1,I1)
      NCS(2,IPT)  = 0
      GOTO 15
  12  J1 = NP123(1,I1)
      J3 = NP123(3,I1)
      JJ1 = NUM2(J1)
      JJ3 = NUM2(J3)
      IC1 = NCS(2,JJ1)
      IC3 = NCS(2,JJ3)
      X1  = PXY(1,JJ1)
      X3  = PXY(1,JJ3)
      IF (DABS(X1-X3).GT.EPS)GOTO 10
      Y1  = PXY(2,JJ1)
      Y3  = PXY(2,JJ3)
      IF (DABS(Y1).GT.2.OR.DABS(Y3).GT.2)GOTO 701
      Y1  = Y1+2+EPS/1000.
      Y3  = Y3+2+EPS/1000.
      K1  = Y1/2
      K3  = Y3/2
      Y1  = Y1-K1*2
      Y3  = Y3-K3*2
      PXY(2,JJ1)=Y1
      PXY(2,JJ3)=Y3
      FMAX =DMAX1(Y1,Y3)
      FMIN =DMIN1(Y1,Y3)
      F    = FMAX-FMIN
      IF (F.GT.1)Y2 = F/2+1.+FMIN
      IF (F.LE.1)Y2 = F/2+FMIN
      IPT = IPT+1
      IPTX = IPTX+1
      IF (IPT.GT.MPT.OR.IPTX.GT.MPT)GOTO 700
      NUM2(IPTX) = IPT
      NP123(2,I1) = IPTX
      PXY(1,IPT) = (X1+X3)/2
      PXY(2,IPT) = Y2
         NCS(1,IPT)  = NCN(1,I1)
         NCS(2,IPT)  = IC1
  15  CONTINUE
      RETURN
 700  WRITE(LP0,500)
      STOP
 701  WRITE(LP0,501)
      STOP
 500  FORMAT('  ERROR(M1): TOO MANY POINTS DURING GENERATION      ')
 501  FORMAT('  ERROR(M2): ILLEGAL  POLAR COORDINATE, MOD(FI)>2*PI')
      END
C

