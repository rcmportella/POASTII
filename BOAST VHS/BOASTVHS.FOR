C      BOASTVHS.FOR
C
C    Main section of BOAST-VHS that sets the size of each array and reads
C   in the II, JJ, and KK values in order to set up the NX, NY, and NZ
C   parameters.
C   The main section call only one subroutine MAINSUB with many arguments.
C
C    For dynamic re-dirensioning the former main section of BOASTVHS is now
C    a subroutine with an argument list.
C
C        (HORIZONTAL/SLANT WELL & CONVENTIONAL VERTICAL WELL MODELS)
C
C        NNX = NO. OF GRID BLOCKS IN X-DIRECTION
C        NNY = NO. OF GRID BLOCKS IN Y-DIRECTION
C        NNZ = NO. OF GRID BLOCKS IN Z-DIRECTION
C	MAX = MAXIMUM BUFFER SIZE FOR THE 3 DIMENSION X,Y, & Z
C
       PARAMETER(NNX=10,NNY=10,NNZ=7)

       PARAMETER(NXP=NNX+1,NYP=NNY+1,NZP=NNZ+1,
     & MAX=NNX*NNY*NNZ,NXPYZ=NXP*NNY*NNZ,NXYPZ=NNX*NYP*NNZ,
     & NXYZP=NNX*NNY*NZP)
C
      PARAMETER (NUM1=20,NUM2=20,NUM4=16,NXY=400)
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      REAL*4 KROT,KRWT,KRGT,MUOT,MUWT,MUGT,MCFG,MCFG1,MCFGT,
     & MCFGI,MBEO,MBEW,MBEG
      CHARACTER*2 TITLE(80)
      CHARACTER*5 WELLID
      CHARACTER*20 INFILE,XAMPLFIL
      COMMON /IDEN/ WELLID(NW)
      COMMON /IO/ NI,NO
      COMMON /IO2/ NO2
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT

      DIMENSION PBOT(MAX),AW(MAX),AE(MAX),AN(MAX),
     & AS(MAX),AB(MAX),AT(MAX),E(MAX),B(MAX)
      DIMENSION PN(MAX),SON(MAX),SWN(MAX),SGN(MAX),
     & A1(MAX),A2(MAX),A3(MAX),
     & SUM(MAX),GAM(MAX)
      REAL KX(MAX),KY(MAX),KZ(MAX)
      DIMENSION EL(MAX)
      DIMENSION P(MAX)
      DIMENSION SO(MAX),SW(MAX),SG(MAX)
      DIMENSION BO(MAX),BW(MAX),BG(MAX)
      DIMENSION QO(MAX),QW(MAX),QG(MAX)
      DIMENSION GOWT(MAX),GWWT(MAX),GGWT(MAX)
      DIMENSION QOWG(MAX)
      DIMENSION CT(MAX)
      DIMENSION VP(MAX)
      DIMENSION DX(MAX),DY(MAX),DZ(MAX)
      DIMENSION TX(NXPYZ),TY(NXYPZ),TZ(NXYZP)
      DIMENSION OW(NXPYZ),OE(NXPYZ),WW(NXPYZ),WE(NXPYZ),
     & OS(NXYPZ),ON(NXYPZ),WS(NXYPZ),WN(NXYPZ),
     & OT(NXYZP),OB(NXYZP),WT(NXYZP),WB(NXYZP)
      DIMENSION OOIP(NZ1),OWIP(NZ1),ODGIP(NZ1),OFGIP(NZ1)
      DIMENSION VAREL(NXY),RDXL(NX1),RDYL(NY1),RDZL(NZ1)
      DIMENSION UM(NX1),AZL(NX1),BZL(NX1),CZL(NX1),DZL(NX1),UZL(NX1)
      REAL*8 DPP(MAX),DPSO(MAX),DPSW(MAX)

      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /SRATE/ PID(NW,NZ1),PWF(NW,NZ1),PWFC(NW,NZ1),KIP(NW),
     & GMO(NW,NZ1),GMW(NW,NZ1),GMG(NW,NZ1),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NZ1),CUMW(NW,NZ1),CUMG(NW,NZ1)

C      DIMENSION BETA(NX1),GAMMA(NX1),W(NX1)

      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NZ1),IQH2(NW,NZ1),IQH3(NW,NZ1),COND(NW)

      COMMON /CODE/ KSM1,KSN1,KCO1,NN,FACT1,FACT2,TMAX,KSOL,MITER,
     & OMEGA,TOL,TOL1,KSN,KSM,KCO,KTR,KCOFF,DSMAX,DPMAX,WORMAX,
     & GORMAX,PAMIN,PAMAX
      COMMON /ADD1/ IM,JM,KM,ETI,FT,FTMAX,DELT
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      COMMON /PSCNTL/ KPI,KSI
      COMMON /MBE/ MBEO,MBEW,MBEG
      COMMON /NUMBER/ II,JJ,KK
      COMMON /VOL/ SCFO,SCFW,SCFG,SCFG1,STBO,STBW,MCFG,MCFG1,MCFGT,
     & MCFGI,STBOI,STBWI,RESVOL
      COMMON /BAL/ OP,WP,GP,WI,GI,PAVG0,PAVG,OPR,WPR,GPR,WIR,GIR,CWOR,
     & WOR,CGOR,GOR
      COMMON /PRTCNT/ IWLCNG,ICHANG,IWLREP,ISUMRY,
     & IPMAP,ISOMAP,ISWMAP,ISGMAP,IPBMAP
C
C      DIMENSION SUM1(MAX),VAREL(MAX),RDXL(NX1),RDYL(NY1),RDZL(NZ1)
C
      CHARACTER*2 FNAME,MNAME,TAG
      CHARACTER*4 LNAME
      COMMON /CHRCTR/ FNAME,MNAME(NUM),LNAME,TAG
      COMMON /TTEST/ NUMPRD,SONTVL
      COMMON /SWTCH/ NGRSW,NTRSW,NRESTART,NTS
      COMMON /R1XYZ/ R11X,R21X,R31X,R41X,R11Y,R21Y,R31Y,R41Y,
     & R11Z,R21Z,R31Z,R41Z
      COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
      COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
      COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K      
      COMMON /COUNT/ N1READ,N2READ
C
C      DATA TX,TY,TZ/NXPYZ*0.0,NXYPZ*0.0,NXYZP*0.0/
C      DATA AW,AS,AT/MAX*0.0,MAX*0.0,MAX*0.0/
C      DATA AE,AN,AB/MAX*0.0,MAX*0.0,MAX*0.0/
C      DATA OW,OE,WW,WE/NXPYZ*0.0,NXPYZ*0.0,NXPYZ*0.0,NXPYZ*0.0/
C      DATA OS,ON,WS,WN/NXYPZ*0.0,NXYPZ*0.0,NXYPZ*0.0,NXYPZ*0.0/
C      DATA OT,OB,WT,WB/NXYZP*0.0,NXYZP*0.0,NXYZP*0.0,NXYZP*0.0/
C      DATA UM,AZL,BZL,CZL,DZL,UZL/NX1*0,NX1*0,NX1*0,NX1*0,NX1*0,NX1*0/
C
C*******************************************************************************
      LOGICAL THERE
      CHARACTER*63 MESSAGE1,MESSAGE2
      CHARACTER*79 DISCLAIM
      DIMENSION DISCLAIM(NUM4)
      DIMENSION MESSAGE1(NUM1),MESSAGE2(NUM2)
      DATA (MESSAGE1(I),I=1,NUM1) /
     &'ษอออออออออออออออออออออ   DOE presents  อออออออออออออออออออออออป',
     &'บ                                                             บ',
     &'บ          BLACK OIL ENGINEERING AND SIMULATION TOOL          บ',
     &'บ                                                             บ',
     &'ฬอออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออน',
     &'บ                                                             บ',
     &'บ                                                             บ',
     &'บ                                                             บ',
     &'บ  ÛÛÛÛÛ  ÛÛÛÛ ÛÛÛÛÛ  ÛÛÛÛ ÛÛÛÛÛ    ÛÛ     ÛÛ ÛÛ   Û   ÛÛÛÛ   บ',
     &'บ  ÛÛ  ÛÛ Û  Û ÛÛ  Û ÛÛ      Û      ÛÛ     ÛÛ ÛÛ   Û  ÛÛ      บ',
     &'บ  ÛÛÛÛÛ  Û  Û ÛÛÛÛÛ  ÛÛÛ    Û  ÛÛÛ  ÛÛ   ÛÛ  ÛÛÛÛÛÛ   ÛÛÛ    บ',
     &'บ  ÛÛ  ÛÛ Û  Û ÛÛ  Û    ÛÛ   Û        ÛÛ ÛÛ   ÛÛ   Û     ÛÛ   บ',
     &'บ  ÛÛÛÛÛ  ÛÛÛÛ ÛÛ  Û ÛÛÛÛ    Û         ÛÛÛ    ÛÛ   Û  ÛÛÛÛ    บ',
     &'บ                                                             บ',
     &'บ                                                             บ',
     &'บ               Vertical/Horizontal/Slant Well                บ',
     &'บ                    Reservoir Simulator                      บ',
     &'บ                        Version 1.3                          บ',
     &'บ                                                             บ',
     &'ศอออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออผ'/
      DATA (MESSAGE2(I),I=1,NUM2) /
     &'ษออออออออออออออออออออออออออ DOE ออออออออออออออออออออออออออออออป',
     &'บ                                                             บ',
     &'บ                 Work Performed for the                      บ',
     &'บ               U.S. Department of Energy                     บ',
     &'บ       Under Cooperative Agreement DE-FC22-83FE60149         บ',
     &'บ               Tom Reid, Project Manager                     บ',
     &'บ              Fred Burtch, Program Manager                   บ',
     &'บ               U.S. Department of Energy                     บ',
     &'บ                      June 1991                              บ',
     &'บ                                                             บ',
     &'ฬออออออออออออออออออออออออ   NIPER   ออออออออออออออออออออออออออน',
     &'บ                                                             บ',
     &'บ    Ming-Ming Chang, Partha Sarathi, Raymond J. Heemstra,    บ',
     &'บ                Aaron Cheng and James F. Pautz               บ',
     &'บ    National Institute for Petroleum and Energy Research     บ',
     &'บ                      P.O. Box 2128                          บ',
     &'บ              Bartlesville, Oklahoma 74005                   บ',
     &'บ                      (918)336-2400                          บ',
     &'บ                                                             บ',
     &'ศอออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออผ'/
C      DATA (MESSAGE3(I),I=1,NUM3) /
C     &'ษออออออออออออออออออออออออ    DOE    ออออออออออออออออออออออออออป',
C     &'บ                                                             บ',
C     &'บ          BLACK OIL ENGINEERING AND SIMULATION TOOL          บ',
C     &'บ                                                             บ',
C     &'ฬอออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออน',
C     &'บ                                                             บ',
C     &'บ                Work Performed for the                       บ',
C     &'บ               U.S. Department of Energy                     บ',
C     &'บ       Under Cooperative Agreement DE-FC22-83FE60149         บ',
C     &'บ               Tom Reid, Project Manager                     บ',
C     &'บ              Fred Burtch, Program Manager                   บ',
C     &'บ               U.S. Department of Energy                     บ',
C     &'บ                      June 1991                              บ',
C     &'บ                                                             บ',
C     &'ศอออออออออออออออออออออออออออออออออออออออออออออออออออออออออออออผ'/
      DATA (DISCLAIM(I),I=1,NUM4)/
     &'ฺฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤ
     &ฤฤฤฤฤฤฤฤฤฤฤฤฤฟ',
     &'ณ                             DISCLAIMER                         
     &             ณ',
     &'ณ                                                                
     &             ณ',
     &'ณ  This software was prepared as an account of work sponsored by 
     &an agency of ณ',
     &'ณthe United States Government.  Neither the United States Governm
     &ent nor any  ณ',
     &'ณof their employees, make any warranty expressed or implied, or a
     &ssumes any   ณ',
     &'ณlegal liability or responsibility for the accuracy, completeness
     &, of useful- ณ',
     &'ณness of any information, apparatus, product, or process disclose
     &d, or repre- ณ',
     &'ณsents that its use would not infringe privately owned rights.  R
     &eference     ณ',
     &'ณherein to any specific commercial product, process, or service b
     &y trade name,ณ',
     &'ณtrademark, manufacturer, or otherwise, does not necessarily cons
     &titute or    ณ',
     &'ณimply its endorsement, recommendation, or favoring by the United
     & States      ณ',
     &'ณGovernment or any agency thereof.  The views and opinions of aut
     &hors express-ณ',
     &'ณed herein do not necessarily state or reflect those of the Unite
     &d States     ณ',
     &'ณGovernment or any agency thereof.                               
     &             ณ',
     &'ภฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤฤ
     &ฤฤฤฤฤฤฤฤฤฤฤฤฤู'/
C
      XAMPLFIL='C:\VHS\BOASTVHS.DAT'
C      VERSION=1.3.1
      MESSAGE1(NUM1-2)(33:37)='1.3.1'
      NI=5
      NO=6
      N1READ=0
      N2READ=0
	CALL CLRSCR
      DO 3 I=1,NUM1
      WRITE(*,501) MESSAGE1(I)
501   FORMAT(6X,A63)
3     CONTINUE
      CALL WAIT(200)
      CALL CLRSCR
      DO 4 I=1,NUM2
      WRITE(*,501) MESSAGE2(I)
4     CONTINUE
      CALL WAIT(400)
      CALL CLRSCR
	WRITE(*,'(//)' )
      DO 6 I=1,NUM4
      WRITE(*,502) DISCLAIM(I)
502   FORMAT(1X,A79)
6     CONTINUE
      PAUSE
      CALL CLRSCR
      WRITE (*,569) MAX,NX1,NY1,NZ1
569   FORMAT(///,
     &T9,'ษออออออออออออออออออออออออออออออออออออออออออออออออออออออออป',/,
     &T9,'บ                                                        บ',/,
     &T9,'บ  THE MAXIMUM NUMBER OF GRID BLOCKS PERMITTED IS ',I4,
     &  '   บ',/,
     &T9,'บ  THE LIMITS IN THE X, Y, & Z AXES ARE : ',3I4,',  บ',/,
     &T9,'บ  RESPECTIVELY.                                         บ',/,
     &T9,'บ                                                        บ',/,
     &T9,'ศออออออออออออออออออออออออออออออออออออออออออออออออออออออออผ')
504   WRITE (*,671) XAMPLFIL
 671  FORMAT(1X,14(/),
     &   1X,'PLEASE ENTER YOUR DISK DRIVE AND INPUT FILE NAME,',/,
     &   '     SUCH AS "',A20,'" : ')
      READ (*,672) INFILE
672   FORMAT(A20)
	IF(INFILE(1:4) .EQ. '    ') INFILE=XAMPLFIL
      INQUIRE(FILE=INFILE,EXIST=THERE)
      IF(.NOT. THERE) WRITE(*,185) INFILE
      IF(.NOT. THERE) GO TO 504
185   FORMAT(1X,'NON-EXISTANT DATA INPUT FILE -- ',A20)
      OPEN(NI,FILE=INFILE,STATUS='OLD')
	READ(NI,69)(TITLE(I),I=1,40)
69	FORMAT(40A2)
	READ(NI,*) II,JJ,KK
	WRITE(*,573) II,JJ,KK,(TITLE(I),I=1,37)
573	FORMAT(1X,'Input coordinates are ',
     &      I3,' by',I3,' by',I2,/,' for ',37A2)
	REWIND (NI)
C
C    SET UP DYNAMIC DIMENSIONS
C
	NX=II
	NY=JJ
	NZ=KK
C
C    GRID BLOCK LIMIT CHECKS
	IJK=II*JJ*KK
	IF(IJK .GT. MAX) GO TO 800
	IF (II .GT. NX1) GO TO 810
	IF (JJ .GT. NY1) GO TO 820
	IF (KK .GT. NZ1) GO TO 830


C*******************************************************************************
C
C      CALL MAINSUB(INFILE,NX,NY,NZ,PBOT,UM,AZL,BZL,CZL,DZL,UZL,
C     & SUM1,VAREL,RDXL,RDYL,RDZL,BETA,GAMMA,W)


	CALL MAINVHS(NX,NY,NZ,PBOT,AW,AE,AN,AS,AB,AT,
     &  E,B,PN,SON,SWN,SGN,A1,A2,A3,
     &  SUM,GAM,AZL,BZL,CZL,DZL,UZL,
     &  KX,KY,KZ,TX,TY,TZ,EL,P,SO,SW,SG,BO,BW,BG,QO,QW,QG,
     &  GOWT,GWWT,GGWT,QOWG,OW,OE,WW,WE,OS,ON,WS,WN,OT,OB,WT,WB,
     &  CT,VP,DX,DY,DZ,OOIP,OWIP,ODGIP,OFGIP,VAREL,UM,
     &  RDXL,RDYL,RDZL,INFILE,TITLE,DPP,DPSO,DPSW)
C
C*******************************************************************************
C      IF(NRESTART.EQ.1)CALL RESTART(INFILE,RSTFILE)     
C        GO TO 840
800	WRITE(NO,900) IJK,MAX
	GO TO 840
810	WRITE(NO,910) II,NX1
	GO TO 840
820	WRITE(NO,920) JJ,NY1
	GO TO 840
830	WRITE(NO,930) KK,NZ1
840	CONTINUE
	CLOSE(NI)
	CLOSE(NO)
	STOP
 900  FORMAT (//1X,'FATAL ERROR;TOTAL NUMBER OF GRID BLOCKS ',
     &'SPECIFIED =',I4,'; EXCEEDS',I4/3X,'****** RUN TERMINATED',
     &' ******'/)
 910  FORMAT (//1X,'FATAL ERROR;TOTAL NUMBER OF X-DIRECTION ',
     &'GRID BLOCKS SPECIFIED =',I4,'; EXCEEDS',I4/3X,'****** RUN ',
     &'TERMINATED *****'/)
 920  FORMAT (//1X,'FATAL ERROR;TOTAL NUMBER OF Y-DIRECTION ',
     &'GRID BLOCKS SPECIFIED =',I4,'; EXCEEDS',I4/3X,'****** RUN ',
     &'TERMINATED ******'/)
 930  FORMAT (//1X,'FATAL ERROR;TOTAL NUMBER OF Z-DIRECTION ',
     &'GRID BLOCKS SPECIFIED =',I4,'; EXCEEDS',I4/3X,'****** RUN ',
     &'TERMINATED ******'/)
C
C	END OF BESTVHS MAIN PROGRAM
C
      END
C       BOAST-VHS.FOR
C        (BOASTVHS.F)
C
C      (HORIZONTAL/SLANT WELL & CONVENTIONAL VERTICAL WELL MODELS)
C       --- AN UPDATED SLANT WELL MODEL  INSTALLED
C       --- OPTION OF RATE AND PRESS. CONSTRAINTS AT THE SAME TIME
C       --- 'RESTART' SUBROUTINE  INSTALLED
C       MAY 1988
C          DIMENSIONED AT 'VHSPARAMA.FOR' & 'VHSPARMB.FOR',
C          REQUIRES "BOASTVHS.DAT" AS INPUT,
C           AND GENERATES "BOASTVHS.OUT" & "BOASTVHS.SUM"
C
C*******************************************************************************
C
	SUBROUTINE MAINVHS(NX,NY,NZ,PBOT,AW,AE,AN,AS,AB,AT,
     &  E,B,PN,SON,SWN,SGN,A1,A2,A3,
     &  SUM,GAM,AZL,BZL,CZL,DZL,UZL,
     &  KX,KY,KZ,TX,TY,TZ,EL,P,SO,SW,SG,BO,BW,BG,QO,QW,QG,
     &  GOWT,GWWT,GGWT,QOWG,OW,OE,WW,WE,OS,ON,WS,WN,OT,OB,WT,WB,
     &  CT,VP,DX,DY,DZ,OOIP,OWIP,ODGIP,OFGIP,
     &  VAREL,UM,RDXL,RDYL,RDZL,INFILE,TITLE,DPP,DPSO,DPSW)
C
      PARAMETER (NUM1=20,NUM2=20,NUM4=16)
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
C	NEW DP LOCAL VARIABLES **RJH**
C       MODIFICATION FOR DOUBLE PRECISION READ OF SO & SW IN RESTART  **RJH**

      CHARACTER*1 IANS
      CHARACTER*20 INFILE,OUTFILE,NEWFILE,RATFILE,RSTFILE
      CHARACTER*20 XAMPLFIL
      REAL KROT,KRWT,KRGT,MUOT,MUWT,MUGT,MCFG,MCFG1,MCFGT,
     & MCFGI,MBEO,MBEW,MBEG
      CHARACTER*2 TITLE(80),FNAME,MNAME,TAG
      CHARACTER*5 WELLID
      CHARACTER*4 LNAME
      COMMON /IDEN/ WELLID(NW)
      COMMON /IO/ NI,NO
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT

      DIMENSION PBOT(NX,NY,NZ)
      DIMENSION AW(NX,NY,NZ),AE(NX,NY,NZ),AN(NX,NY,NZ),
     & AS(NX,NY,NZ),AB(NX,NY,NZ),AT(NX,NY,NZ),E(NX,NY,NZ),B(NX,NY,NZ)
      DIMENSION PN(NX,NY,NZ),
     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      REAL KX(NX,NY,NZ),KY(NX,NY,NZ),KZ(NX,NY,NZ)
      DIMENSION TX(NX+1,NY,NZ),TY(NX,NY+1,NZ),TZ(NX,NY,NZ+1)
      DIMENSION EL(NX,NY,NZ)
      DIMENSION P(NX,NY,NZ)
      DIMENSION SO(NX,NY,NZ),SW(NX,NY,NZ),SG(NX,NY,NZ)
      DIMENSION VAREL(NX,NY),RDXL(NX),RDYL(NY),RDZL(NZ)
      DIMENSION UM(NX+1),AZL(NX+1),BZL(NX+1),CZL(NX+1),
     &          DZL(NX+1),UZL(NX+1)
      DIMENSION BO(NX,NY,NZ),BW(NX,NY,NZ),BG(NX,NY,NZ)
      DIMENSION QO(NX,NY,NZ),QW(NX,NY,NZ),QG(NX,NY,NZ)
      DIMENSION GOWT(NX,NY,NZ),GWWT(NX,NY,NZ),GGWT(NX,NY,NZ),
     & OW(NX+1,NY,NZ),OE(NX+1,NY,NZ),WW(NX+1,NY,NZ),WE(NX+1,NY,NZ),
     & OS(NX,NY+1,NZ),ON(NX,NY+1,NZ),WS(NX,NY+1,NZ),WN(NX,NY+1,NZ),
     & OT(NX,NY,NZ+1),OB(NX,NY,NZ+1),WT(NX,NY,NZ+1),WB(NX,NY,NZ+1)
      DIMENSION QOWG(NX,NY,NZ)
      DIMENSION CT(NX,NY,NZ)
      DIMENSION VP(NX,NY,NZ)
      DIMENSION DX(NX,NY,NZ),DY(NX,NY,NZ),DZ(NX,NY,NZ)
      DIMENSION OOIP(NZ),OWIP(NZ),ODGIP(NZ),OFGIP(NZ)
      REAL*8 DPP(NX,NY,NZ),DPSO(NX,NY,NZ),DPSW(NX,NY,NZ)

      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)

      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
      COMMON /CODE/ KSM1,KSN1,KCO1,NN,FACT1,FACT2,TMAX,KSOL,MITER,
     & OMEGA,TOL,TOL1,KSN,KSM,KCO,KTR,KCOFF,DSMAX,DPMAX,WORMAX,
     & GORMAX,PAMIN,PAMAX
      COMMON /ADD1/ IM,JM,KM,ETI,FT,FTMAX,DELT
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      COMMON /PSCNTL/ KPI,KSI
      COMMON /MBE/ MBEO,MBEW,MBEG
      COMMON /NUMBER/ II,JJ,KK
      COMMON /VOL/ SCFO,SCFW,SCFG,SCFG1,STBO,STBW,MCFG,MCFG1,MCFGT,
     & MCFGI,STBOI,STBWI,RESVOL
      COMMON /BAL/ OP,WP,GP,WI,GI,PAVG0,PAVG,OPR,WPR,GPR,WIR,GIR,CWOR,
     & WOR,CGOR,GOR
      COMMON /PRTCNT/ IWLCNG,ICHANG,IWLREP,ISUMRY,
     & IPMAP,ISOMAP,ISWMAP,ISGMAP,IPBMAP
      COMMON /IO2/ NO2
      COMMON /CHRCTR/ FNAME,MNAME(NUM),LNAME,TAG
      COMMON /TTEST/ NUMPRD,SONTVL
      COMMON /SWTCH/ NGRSW,NTRSW,NRESTART,NTS
      COMMON /R1XYZ/ R11X,R21X,R31X,R41X,R11Y,R21Y,R31Y,R41Y,
     & R11Z,R21Z,R31Z,R41Z
      COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
      COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
      COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K
      COMMON /COUNT/ N1READ,N2READ
      LOGICAL THERE
C      DIMENSION BETA(NX1),GAMMA(NX1),W(NX1)


	CALL VARINIT(NX,NY,NZ,TX,TY,TZ,AW,AS,AT,AE,AN,AB,E,B,
     &   GOWT,GWWT,GGWT,OW,OE,WW,WE,OS,ON,WS,WN,OT,OB,WT,WB)
C
      VERSION=1.3
	ILAST=0
	NGAUGE=1
	IM=NX
	JM=NY
	KM=NZ
505	CALL EXTEN(INFILE,OUTFILE,'OUT')
      WRITE(*,673) OUTFILE
  673 FORMAT(1X,'PLEASE ENTER YOUR DISK DRIVE AND OUTPUT FILE NAME,',/,
     &       1X,'( Another CR will default to "',A20,'" ) : ' )
      READ (*,672) NEWFILE
672   FORMAT(A20)
      IF(NEWFILE(1:2) .NE. '  ') OUTFILE = NEWFILE
      INQUIRE(FILE=OUTFILE,EXIST=THERE)
      IF(THERE) THEN
       WRITE(*,186) OUTFILE
 186  FORMAT(1X,'AN OUTPUT FILE ALREADY EXISTS NAMED : ',A20 )
      WRITE(*,187)
 187  FORMAT(1X,'DO YOU WISH TO DELETE THE OLD FILE? (Y/N) : ' )
       READ (*,'(A1)' ) IANS
       IF( (IANS .NE. 'Y') .AND. (IANS .NE. 'y') ) GO TO 505
      ENDIF
      CALL EXTEN(OUTFILE,RATFILE,'SUM')
      CALL EXTEN(OUTFILE,RSTFILE,'RST')
      IF( INFILE .EQ. RSTFILE )  CALL EXTEN(OUTFILE,RSTFILE,'NEW')
      OPEN(NO,FILE=OUTFILE,STATUS='UNKNOWN')
      OPEN(14,FILE=RATFILE,STATUS='UNKNOWN')
      WRITE(*,3004) RATFILE
3004  FORMAT(10X,'Writing a summary file to : ',A20)
      WRITE(*,3003) OUTFILE
3003  FORMAT(10X,'Writing a complete output file to : ',A20)
      WRITE(14,3007) VERSION
      WRITE(NO,3005)
      WRITE(NO,3007) VERSION
      WRITE(NO,3006)
 3007 FORMAT(T30,'*',T45,'U.S. Department of Energy: ',T98,'*',/T30,
     &'*',T43,'VERTICAL/HORIZONTAL/SLANT WELL RESERVOIR SIMULATOR',T98,
     & '*',/T30,'*',T59,'(VERSION ',F3.1,')',T98,'*')
 3005 FORMAT(//T30,69('*'),/T30,'*',T98,'*'/T30,'*',T98,'*')
 3006 FORMAT(T30,'*',T98,'*'/T30,'*',T98,'*'/T30,69('*'),///)
      WRITE (*,519)
519   FORMAT(/10X,' *** PROGRAM IN PROGRESS, PLEASE DO NOT DISTURB ***')
C
C     ********************************************************
C
C
C**** ESTABLISH RESERVOIR AND BLOCK DIMENSIONS
C
      CALL GRID1(NX,NY,NZ,EL,DX,DY,DZ,SUM,VAREL,RDXL,RDYL,RDZL)
C
C**** ESTABLISH POROSITY AND PERMEABILITY AT EACH ZONE
C
      CALL PARM1(NX,NY,NZ,VP,PV,KX,KY,KZ,DX,DY,DZ)
C
C**** CALCULATE INTERBLOCK TRANSMISSIBILITIES
C
      CALL TRAN1(NX,NY,NZ,A1,A2,A3,KX,KY,KZ,TX,TY,TZ,DX,DY,DZ)
C
C**** EMPIRICAL DATA
C
      CALL TABLE
C
C**** ESTABLISH INITIAL CONDITIONS
C
      CALL UINIT1(NX,NY,NZ,PBOT,PN,SON,SWN,SGN,EL,P,SO,SW,SG,
     &   DPP,DPSO,DPSW)
C
C**** SOLUTION METHOD, DEBUG PRINT, AND TIME-STEP CONTROL PARAMETERS.
C
      CALL CODES
      READ(NI,69)(TITLE(IH),IH=1,40)
69    FORMAT(40A2)
      D5615 = 1./5.615
      D288  = 1./288.
C      D114  = 1./144.
      NMAX  = NN+1
      DO 1000 N=1,NMAX
C
C*** NLOOP DEFINED TO AVOID INDEX MODIFICATION WARNINGS IN
C    CALLS TO MATBAL, SOLMAT, LSOR.
C      NLOOP=N
C
C**** RECURRENT DATA.
C
      IF(FT.LT.FTMAX) GO TO 46
      READ(NI,*) IWLCNG,ICHANG,IWLREP,ISUMRY,
     & IPMAP,ISOMAP,ISWMAP,ISGMAP,IPBMAP
      READ(NI,*) DAY,DTMIN,DTMAX,HOUR,MIN,SEC
      IF(IWLCNG.EQ.0) GO TO 45
      CALL NODES(NVQN,NVQNH,NVQNS,NX,NY,NZ,PN,EL,
     &    KX,KY,KZ,DX,DY,DZ)
   45 DELT=DAY+HOUR/24.+MIN/1440.+SEC/86400.
      FTMAX=ETI+ICHANG*DELT
   46 IF(N.EQ.1)DELT0=DELT
      IF(N.EQ.1)GO TO 1050
      IF(DSMC.LT.DSMAX.AND.DPMC.LT.DPMAX)DELT=DELT*FACT1
      IF(DSMC.GT.DSMAX.OR.DPMC.GT.DPMAX)DELT=DELT*FACT2
      IF(DELT.LT.DTMIN)DELT=DTMIN
      IF(DELT.GT.DTMAX)DELT=DTMAX
      IF(ETI+DELT.GT.FTMAX)DELT=FTMAX-ETI
1050  FT=ETI+DELT
C      WRITE(*,*) FT
C     TIME STEP CHANGE FOR RESTART ONLY **RJH**
      IF(ETI+DELT*0.5.GE.TMAX.AND.NRESTART.EQ.0) GO TO 1010
      IF(ETI+DELT.GT.TMAX.AND.NRESTART.EQ.1) GO TO 1010

C****ITFLAG COUNTS THE NUMBER OF TIME STEP REPITITIONS.
      ITFLAG=0
C****REENTRY POINT FOR REPEATED TIME STEP.
1060  CONTINUE
      DIV1 = 1./DELT
      IF(N.GT.1.OR.ITFLAG.GT.0) GO TO 105
      RESVOL=0.0
      SCFO=0.0
      SCFG=0.0
      SCFG1=0.0
      DO 102 K=1,KK
      OOIP(K)=0.
      OWIP(K)=0.
      ODGIP(K)=0.
      OFGIP(K)=0.
      DO 100 J=1,JJ
      DO 100 I=1,II
      PP=P(I,J,K)
      BPT=PBOT(I,J,K)
      VP(I,J,K)=VP(I,J,K)*DX(I,J,K)*DY(I,J,K)*DZ(I,J,K)
      RESVOL=RESVOL+VP(I,J,K)
C**** NOTE: WE ARE ASSUMING INITIAL PHI IS AT INITIAL RESERVOIR PRESSURE
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PP,RSO)
      CALL INTERP(PWT,RSWT,MPWT,PP,RSW)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PP,BO(I,J,K))
      CALL INTERP(PWT,BWT,MPWT,PP,BW(I,J,K))
      CALL INTERP(PGT,BGT,MPGT,PP,BG(I,J,K))
C      WRITE(6,1) RSO,RSW,BO(I,J,K),BW(I,J,K),BG(I,J,K)
C1     FORMAT(5X,'RSO = ',F15.5,'RSW = ',F15.5,'BO = ',F15.5,'BW = ',
C     & F15.5,'BG = ',F15.5)
      FF1=SO(I,J,K)/BO(I,J,K)
      FF2=SW(I,J,K)/BW(I,J,K)
      SCFO=SCFO+VP(I,J,K)*FF1
      SCFW=SCFW+VP(I,J,K)*FF2
      SCFG=SCFG+VP(I,J,K)*SG(I,J,K)/BG(I,J,K)
      SCFG1=SCFG1+VP(I,J,K)*(RSO*FF1+RSW*FF2)
      CALL INTERP(POT,BOPT,MPOT,PP,BODER)
      CALL INTERP(POT,RSOPT,MPOT,PP,RSODER)
      CALL INTERP(PWT,BWPT,MPWT,PP,BWDER)
      CALL INTERP(PWT,RSWPT,MPWT,PP,RSWDER)
      CALL INTERP(PGT,BGPT,MPGT,PP,BGDER)
      IF(PP.GT.PBOT(I,J,K))BODER=BSLOPE
      IF(PP.GT.PBOT(I,J,K))RSODER=RSLOPE
      CO=-(BODER-BG(I,J,K)*RSODER)/BO(I,J,K)
      CW=-(BWDER-BG(I,J,K)*RSWDER)/BW(I,J,K)
      CG=-BGDER/BG(I,J,K)
      CALL INTERP(PGT,CRT,MPGT,PP,CR)
      CT(I,J,K)=CR + CO*SO(I,J,K) +CW*SW(I,J,K) + CG*SG(I,J,K)
      OOIP(K)=OOIP(K)+D5615*.000001*SON(I,J,K)*VP(I,J,K)/BO(I,J,K)
      OWIP(K)=OWIP(K)+D5615*.000001*SWN(I,J,K)*VP(I,J,K)/BW(I,J,K)
      ODGIP(K)=ODGIP(K)+.001*.000001*(RSO*SON(I,J,K)*VP(I,J,K)/BO(I,J,K)
     & +RSW*SWN(I,J,K)*VP(I,J,K)/BW(I,J,K))
      OFGIP(K)=OFGIP(K)+.001*.000001*SGN(I,J,K)*VP(I,J,K)/BG(I,J,K)
      IF(KCO1.EQ.0)GO TO 100
      WRITE(NO,21)I,J,K,VP(I,J,K),CT(I,J,K),BO(I,J,K),SO(I,J,K),
     &BW(I,J,K),SW(I,J,K),BG(I,J,K),SG(I,J,K)
  21  FORMAT(1X,3I3,8E15.6)
  100 CONTINUE
      WRITE(NO,110)K,OOIP(K),OWIP(K),ODGIP(K),OFGIP(K)
  110 FORMAT(/,1X,'LAYER',I3,' INITIAL FLUID VOLUMES:',
     & /,10X,'OIL IN PLACE (MILLION STB)',T60,F10.4,
     & /,10X,'WATER IN PLACE (MILLION STB)',T60,F10.4,
     & /,10X,'SOLUTION GAS IN PLACE (BILLION SCF)',T60,F10.4,
     & /,10X,'FREE GAS IN PLACE (BILLION SCF)',T60,F10.4/)
  102 CONTINUE
      TOOIP=0.
      TOWIP=0.
      TODGIP=0.
      TOFGIP=0.
      DO 103 K=1,KK
      TOOIP=TOOIP+OOIP(K)
      TOWIP=TOWIP+OWIP(K)
      TODGIP=TODGIP+ODGIP(K)
  103 TOFGIP=TOFGIP+OFGIP(K)
      WRITE(NO,115) TOOIP,TOWIP,TODGIP,TOFGIP
  115 FORMAT(/,1X,'TOTAL INITIAL FLUID VOLUMES IN RESERVOIR:',
     & /,10X,'OIL IN PLACE (MILLION STB)',T60,F10.4,
     & /,10X,'WATER IN PLACE (MILLION STB)',T60,F10.4,
     & /,10X,'SOLUTION GAS IN PLACE (BILLION SCF)',T60,F10.4,
     & /,10X,'FREE GAS IN PLACE (BILLION SCF)',T60,F10.4/)
      STBO=SCFO*D5615
      STBW=SCFW*D5615
      MCFG=SCFG*0.001
      MCFG1=SCFG1*0.001
      STBOI=STBO
      STBWI=STBW
      MCFGI=MCFG+MCFG1
      IF(MCFGI.LE.1.D-7.AND.MCFGT.LE.1.D-7)MBEG=0.0
105   IF(N.EQ.1.AND.ITFLAG.LE.0) CALL PRTPS(N,DELT0,NX,NY,NZ,
     &   PBOT,PN,P,SON,SWN,SGN,SO,SW,SG)
      IF(N.EQ.NMAX) GO TO 1010
C
C**** ESTABLISH RATES & CALCULATE BHFP(IF PID IS NONZERO)
C
      IF(NVQN.NE.0) CALL QRATE(1,NVQN,NVQNH,NX,NY,NZ,PBOT,
     &                         PN,P,SO,SW,SG,QO,QW,QG)
      IF(NVQNH.NE.0)CALL QRATE(2,NVQN,NVQNH,NX,NY,NZ,PBOT,
     &                         PN,P,SO,SW,SG,QO,QW,QG)
      IF(NVQNS.NE.0)CALL QRATE(2,NVQN+NVQNH,NVQNS,NX,NY,NZ,PBOT,
     &                         PN,P,SO,SW,SG,QO,QW,QG)
C      CALL QRATE(NVQN)
C
C**** CALCULATE SEVEN DIAGONAL MATRIX FOR PRESSURE SOLUTION
C
1160  CALL SOLMAT(DIV1,D288,N,NX,NY,NZ,PBOT,
     &  AW,AE,AN,AS,AB,AT,E,B,SUM,GAM,EL,TX,TY,TZ,P,SO,
     &  SW,SG,BO,BW,BG,QO,QW,QG,GOWT,GWWT,GGWT,QOWG,OW,OE,
     &  WW,WE,OS,ON,WS,WN,OT,OB,WT,WB,CT,VP)
C
C**** MODIFY MATRIX ELEMENTS FOR WELLS UNDER IMPLICIT CONTROL.
C
      IF(NVQN.NE.0) CALL PRATEI(1,NVQN,NVQNH,NX,NY,NZ,PBOT,E,B,PN)
      IF(NVQNH.NE.0)CALL PRATEI(2,NVQN,NVQNH,NX,NY,NZ,PBOT,E,B,PN)
      IF(NVQNS.NE.0)CALL PRATEI(2,NVQN+NVQNH,NVQNS,NX,NY,NZ,PBOT,
     &    E,B,PN)
C      CALL PRATEI(NVQN)
C
C**** CALCULATE NEW PRESSURE DISTRIBUTION
C
1170  CALL LSOR(N,NX,NY,NZ,AW,AE,AN,AS,AB,AT,E,B,P,
     &          AZL,BZL,CZL,DZL,UZL,UM)
C
C**** CALCULATE IMPLICIT RATES.
C
      IF(NVQN.NE.0) CALL PRATEO(1,NVQN,NVQNH,NX,NY,NZ,PBOT,P,PN,
     &  BO,BW,BG,QO,QW,QG)
      IF(NVQNH.NE.0)CALL PRATEO(2,NVQN,NVQNH,NX,NY,NZ,PBOT,P,PN,
     &  BO,BW,BG,QO,QW,QG)
      IF(NVQNS.NE.0)CALL PRATEO(2,NVQN+NVQNH,NVQNS,NX,NY,NZ,PBOT,P,PN,
     &  BO,BW,BG,QO,QW,QG)
C      CALL PRATEO(NVQN)
 2051 CONTINUE
C
C**** CALCULATE NEW FLUID SATURATIONS
C
      SCFO=0.0
      SCFW=0.0
      SCFG=0.0
      SCFG1=0.0
      RESVOL=0.0
      DO 400 K=1,KK
      DO 400 J=1,JJ
      DO 400 I=1,II
      PPN=PN(I,J,K)
      PP=P(I,J,K)
      BPT=PBOT(I,J,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PP,RSO)
      CALL INTERP(PWT,RSWT,MPWT,PP,RSW)
      CALL INTERP(PGT,CRT,MPGT,PPN,CR)
      BPT=PBOT(I,J,K)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PP,BBO)
      CALL INTERP(PWT,BWT,MPWT,PP,BBW)
      CALL INTERP(PGT,BGT,MPGT,PP,BBG)
      VPP=VP(I,J,K) * (1.0+CR*(P(I,J,K)-PPN))
      RESVOL=RESVOL+VPP
C
      DP1=0.0
      DP2=0.0
      DP3=0.0
      DP4=0.0
      DP5=0.0
      DP6=0.0
      IF((I-1).GT.0)  DP1=P(I-1,J,K)-PP
      IF((I+1).LE.II) DP2=P(I+1,J,K)-PP
      IF((J-1).GT.0)  DP3=P(I,J-1,K)-PP
      IF((J+1).LE.JJ) DP4=P(I,J+1,K)-PP
      IF((K-1).GT.0)  DP5=P(I,J,K-1)-PP
      IF((K+1).LE.KK) DP6=P(I,J,K+1)-PP
C
      DAODP=OW(I,J,K)*DP1+OE(I,J,K)*DP2+OS(I,J,K)*DP3
     &     +ON(I,J,K)*DP4+OT(I,J,K)*DP5+OB(I,J,K)*DP6
      DAWDP=WW(I,J,K)*DP1+WE(I,J,K)*DP2+WS(I,J,K)*DP3
     &     +WN(I,J,K)*DP4+WT(I,J,K)*DP5+WB(I,J,K)*DP6
      SW(I,J,K)=((DAWDP+GWWT(I,J,K)-QW(I,J,K))*DELT +VP(I,J,K)*
     &SWN(I,J,K)/BW(I,J,K)) * BBW/VPP
      SO(I,J,K)=((DAODP+GOWT(I,J,K)-QO(I,J,K))*DELT +VP(I,J,K)*
     &SON(I,J,K)/BO(I,J,K)) * BBO/VPP
      SG(I,J,K)=1.0-SO(I,J,K)-SW(I,J,K)
      IF(SG(I,J,K).GT.0.0) GO TO 405
      SG(I,J,K)=0.0
       SO(I,J,K)=1.0-SW(I,J,K)
 975   FORMAT(F8.2,2I4,2F9.6,F9.1,F9.6)    
  405 CONTINUE
C
      IF(KCOFF.EQ.0) GO TO 397
      RHO1=VPP*SO(I,J,K)/BBO
      RHO2=VP(I,J,K)*SON(I,J,K)/BO(I,J,K)
      DIFFO=RHO1-RHO2
      RHW1=VPP*SW(I,J,K)/BBW
      RHW2=VP(I,J,K)*SWN(I,J,K)/BW(I,J,K)
      DIFFW=RHW1-RHW2
      RHG1=VPP*SG(I,J,K)/BBG
      RHG2=VP(I,J,K)*SGN(I,J,K)/BG(I,J,K)
      DIFFG=RHG1-RHG2
      WRITE(NO,33)
      WRITE(NO,21)I,J,K,P(I,J,K),SO(I,J,K),SON(I,J,K),SW(I,J,K),
     &SWN(I,J,K),SG(I,J,K),SGN(I,J,K),VPP
      WRITE(NO,21)I,J,K,GOWT(I,J,K),QO(I,J,K),GWWT(I,J,K),
     &QW(I,J,K)
      WRITE(NO,21)I,J,K,DAODP,DAWDP,DELT
      WRITE(NO,21)I,J,K,RHO1,RHO2,DIFFO
      WRITE(NO,21)I,J,K,RHW1,RHW2,DIFFW
      WRITE(NO,21)I,J,K,RHG1,RHG2,DIFFG,GGWT(I,J,K),QG(I,J,K)
C
 397  VP(I,J,K)=VPP
      BO(I,J,K)=BBO
      BW(I,J,K)=BBW
      BG(I,J,K)=BBG
      FF1=SO(I,J,K)/BO(I,J,K)
      FF2=SW(I,J,K)/BW(I,J,K)
      SCFO=SCFO+VP(I,J,K)*FF1
      SCFW=SCFW+VP(I,J,K)*FF2
      SCFG=SCFG+VP(I,J,K)*SG(I,J,K)/BG(I,J,K)
      SCFG1=SCFG1+VP(I,J,K)*(RSO*FF1+RSW*FF2)
      CALL INTERP(POT,BOPT,MPOT,PP,BODER)
      CALL INTERP(POT,RSOPT,MPOT,PP,RSODER)
      CALL INTERP(PWT,BWPT,MPWT,PP,BWDER)
      CALL INTERP(PWT,RSWPT,MPWT,PP,RSWDER)
      CALL INTERP(PGT,BGPT,MPGT,PP,BGDER)
      IF(PP.GT.PBOT(I,J,K))BODER=BSLOPE
      IF(PP.GT.PBOT(I,J,K))RSODER=RSLOPE
      CO=-(BODER-BG(I,J,K)*RSODER)/BO(I,J,K)
      CW=-(BWDER-BG(I,J,K)*RSWDER)/BW(I,J,K)
      CG=-BGDER/BG(I,J,K)
      CALL INTERP(PGT,CRT,MPGT,PP,CR)
      CT(I,J,K)=CR + CO*SO(I,J,K) +CW*SW(I,J,K) + CG*SG(I,J,K)
      IF(N.EQ.KCO)WRITE(NO,21)I,J,K,CR,CO,RSO,CW,RSW,CG
C
  400 CONTINUE
C      AUTO. TIME STEP CONTROL CALC. OF PRESSURE AND SAT. MAXIMA.
      PPM=0.
      SOM=0.
      SWM=0.
      SGM=0.
      DO 240 K=1,KK
      DO 240 J=1,JJ
      DO 240 I=1,II
      DPO=P(I,J,K)-PN(I,J,K)
      DSO=SO(I,J,K)-SON(I,J,K)
      DSW=SW(I,J,K)-SWN(I,J,K)
      DSG=SG(I,J,K)-SGN(I,J,K)
      IF(ABS(DPO).GT.ABS(PPM)) PPM=DPO
      IF(ABS(DSO).GT.ABS(SOM)) SOM=DSO
      IF(ABS(DSW).GT.ABS(SWM)) SWM=DSW
      IF(ABS(DSG).GT.ABS(SGM)) SGM=DSG
 240  CONTINUE
      DPMC=ABS(PPM)
      DSMC=ABS(SOM)
      IF(DSMC.LT.ABS(SGM))DSMC=ABS(SGM)
      IF(DSMC.LT.ABS(SWM))DSMC=ABS(SWM)
C****REPEAT TIME STEP?
      IF(DSMC.LT.DSMAX.AND.DPMC.LT.DPMAX) GO TO 402
      IF(DELT.LE.DTMIN.OR.FACT2.GE.1.0) GO TO 402
      ITFLAG=ITFLAG+1
      DELT=DELT*FACT2
      IF(DELT.LT.DTMIN) DELT=DTMIN
      FT=ETI+DELT
C      WRITE(*,*) FT
      IF(FT.GT.FTMAX) DELT=FTMAX-ETI
C****RESET VARIABLES.
      DO 250 I=1,II
      DO 250 J=1,JJ
      DO 250 K=1,KK
      P(I,J,K)=PN(I,J,K)
      SO(I,J,K)=SON(I,J,K)
      SW(I,J,K)=SWN(I,J,K)
      SG(I,J,K)=SGN(I,J,K)
  250 CONTINUE
C
  33  FORMAT(//)
C
      GO TO 1060
 402  CONTINUE
C
C****UNDERSATURATED GRID BLOCK SATURATION CALCULATION.
C
      DO 410 I=1,II
      DO 410 J=1,JJ
      DO 410 K=1,KK
      IF(P(I,J,K).GT.PN(I,J,K)) GO TO 410
      IF(P(I,J,K).LT.PBOT(I,J,K)) GO TO 410
      IP=I+1
      IM=I-1
      JP=J+1
      JM=J-1
      KP=K+1
      KM=K-1
      IF(IP.GT.II) GO TO 412
      IF(SGN(IP,J,K).GT.0.0001) GO TO 410
  412 IF(IM.LT.1) GO TO 414
      IF(SGN(IM,J,K).GT.0.0001) GO TO 410
  414 IF(JP.GT.JJ) GO TO 416
      IF(SGN(I,JP,K).GT.0.0001) GO TO 410
  416 IF(JM.LT.1) GO TO 418
      IF(SGN(I,JM,K).GT.0.0001) GO TO 410
  418 IF(KP.GT.KK) GO TO 420
      IF(SGN(I,J,KP).GT.0.0001) GO TO 410
  420 IF(KM.LT.1) GO TO 422
      IF(SGN(I,J,KM).GT.0.0001) GO TO 410
  422 SG(I,J,K)=0.0
       SO(I,J,K)=1.0-SW(I,J,K)
  410 CONTINUE
C**** REPRESSURIZATION ALGORITHM.
      IF(IREPRS.EQ.1) GO TO 51
      DO 50 I=1,II
      DO 50 J=1,JJ
      DO 50 K=1,KK
      IF(SG(I,J,K).LE.0.0001) GO TO 50
      PP=P(I,J,K)
      IF(P(I,J,K).GT.PBOT(I,J,K)) PP=PBOT(I,J,K)
      CALL INTERP(POT,BOT,MPOT,PP,BBO)
      CALL INTERP(POT,RSOT,MPOT,PP,RSO)
      CALL INTERP(PGT,BGT,MPGT,PP,BBG)
      IF(SO(I,J,K).EQ.0.0) GO TO 50
      RSONEW=RSO + SG(I,J,K)*BBO/(SO(I,J,K)*BBG)
      CALL INTERP(RSOT,POT,MPOT,RSONEW,PBONEW)
      PBOT(I,J,K)=PBONEW
   50 CONTINUE
   51 CONTINUE
C**** UPDATE OLD FLUID VOLUMES FOR MATERIAL BALANCE.
      STBOI=STBO
      STBWI=STBW
      MCFGI=MCFGT
C**** UPDATE  NEW FLUID VOLUMES.
      STBO=SCFO*D5615
      STBW=SCFW*D5615
      MCFG=SCFG*0.001
      MCFG1=SCFG1*0.001
      MCFGT=MCFG+MCFG1
C
C***** DEBUG PRINT OF PRESENT AND FUTURE P,SO,SW,SG VALUES.
      IF(KCOFF.EQ.0)GO TO 291
      DO 290 K=1,KK
      DO 290 J=1,JJ
      DO 290 I=1,II
      WRITE(NO,21)I,J,K,PN(I,J,K),SON(I,J,K),SWN(I,J,K),SGN(I,J,K)
      WRITE(NO,21)I,J,K,P(I,J,K),SO(I,J,K),SW(I,J,K),SG(I,J,K)
290   CONTINUE
291   CONTINUE
C
C**** WELL REPORT (ALL RATE & PRESSURE DATA APPLICABLE THIS STEP)
C
      IJ=0
      TOR=0.
      TGR=0.
      TWR=0.
      TOC=0.
      TGC=0.
      TWC=0.
      N2=NVQN+NVQNH+NVQNS
      DO 2050 J=1,N2
      GOR=0.
      WOR=0.
      IJ=IJ+1
      IF (J.LE.NVQN)THEN 
       IQ1=IQN1(J)
       IQ2=IQN2(J)
       IQ3=IQN3(J)
       LAY=IQ3+(LAYER(J)-1)
      ELSE
       IQ3=1
       LAY=LAYER(J)
      ENDIF
      DO 2050 KKK=IQ3,LAY
      K=KKK
      IF (J.GT.NVQN)THEN 
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      QOO=QO(IQ1,IQ2,K)*D5615
      QWW=QW(IQ1,IQ2,K)*D5615
      QGG=QG(IQ1,IQ2,K)*.001
      CUMO(J,KKK)=CUMO(J,KKK)+QOO*DELT*.001
      CUMW(J,KKK)=CUMW(J,KKK)+QWW*DELT*.001
      CUMG(J,KKK)=CUMG(J,KKK)+QGG*DELT*.001
      IF(IWLREP.EQ.0) GO TO 2050
C      IF(FT.GT.0.999.AND.FT.LT.1.001) GO TO 891
C      IF(ABS(FT-NUMPRD*SONTVL).GT.0.001) GO TO 2050
  891 IF(IJ.EQ.1.AND.KKK.EQ.IQ3)WRITE(NO,591)FT
      IF(IJ.EQ.1.AND.KKK.EQ.IQ3)WRITE(NO,5911)
 5911 FORMAT(/,T56,'------ RATE ------',22X,'--- CUMULATIVE ---',
     & /,13X,'WELL LOCATION',4X,'CALC    SPEC    SPEC',4X,
     & 'OIL     GAS    WATER    GOR     WOR',5X,
     & 'OIL     GAS    WATER',/,
     & 14X,'ID',3X,'I  J  K    BHFP    BHFP     PI',
     & 3X,'STB/D    MCF/D   STB/D',20X,'MSTB    MMCF    MSTB',/)
      IF(QOO.EQ.0.)GO TO 998
      GOR=QGG*1000./QOO
      WOR=QWW/QOO
  998 WRITE(NO,592)WELLID(J),IQ1,IQ2,K,PWFC(J,KKK),PWF(J,KKK),
     & PID(J,KKK),QOO,QGG,QWW,GOR,WOR,CUMO(J,KKK),CUMG(J,KKK),
     & CUMW(J,KKK)
  592 FORMAT(11X,A5,1X,3I3,2F8.0,F7.3,3F9.0,F7.0,F7.3,3F8.0)
      TOR=TOR+QOO
      TGR=TGR+QGG
      TWR=TWR+QWW
      TOC=TOC+CUMO(J,KKK)
      TGC=TGC+CUMG(J,KKK)
      TWC=TWC+CUMW(J,KKK)
 2050 CONTINUE
      IF(IWLREP.EQ.0)GO TO 2052
C      IF(FT.GT.0.999.AND.FT.LT.1.001) GO TO 791
C      IF(ABS(FT-NUMPRD*SONTVL).GT.0.001) GO TO 2052
  791 WRITE(NO,5912)TOR,TGR,TWR,TOC,TGC,TWC
 5912 FORMAT(12X,102('-'),/,
     &12X,'TOTALS',31X,3F9.0,14X,3F8.0,/)
 2052 CONTINUE
C
C*****CALCULATE MATERIAL BALANCE ERRORS & AVERAGE RESERVOIR PRESSURE
C
      DELT0=DELT
      ETI=ETI+DELT
      CALL MATBAL(DELT0,D5615,NX,NY,NZ,PN,P,QO,QW,QG,VP,PV)
      IF(WOR.GT.WORMAX) GO TO 1002
      IF(GOR.GT.GORMAX) GO TO 1003
      IF(PAVG.LT.PAMIN) GO TO 1004
      IF(PAVG.GT.PAMAX) GO TO 1005
C
C***** SUMMARY REPORT.
      IF(MOD(N,NTS).EQ.0) THEN
       WRITE(14,987)N,ETI,PAVG,OPR,COP,GPR,CGP,GOR,WPR,CWP,WOR
 987  FORMAT(I6,1X,F8.2,1X,F8.2,1X,F8.2,1X,F10.0,1X,F10.2,1X,F11.0,
     &    1X,F8.2,1X,F8.2,1X,F8.0,1X,F6.2)
C
C     BUMP THE PROGRESS GAUGE
C
       CALL GAUGE(NGAUGE,ILAST,ETI,TMAX)
	ENDIF
      IF(ISUMRY.EQ.0) GO TO 2057  
      NLP=N+1
  691 CALL PRTPS(NLP,DELT0,NX,NY,NZ,PBOT,PN,P,SON,SWN,SGN,
     &    SO,SW,SG)
 2057 IF(N.NE.KCO.OR.KCO1.EQ.0)GO TO 500
      DO 300 K=1,KK
      DO 300 J=1,JJ
      DO 300 I=1,II
      WRITE(NO,21)I,J,K,VP(I,J,K),CT(I,J,K),BO(I,J,K),SO(I,J,K),
     &BW(I,J,K),SW(I,J,K),BG(I,J,K),SG(I,J,K)
300   CONTINUE
C
  500 CONTINUE
      IF(N.EQ.KSN)KSN=KSN+KSN1
      IF(N.EQ.KSM)KSM=KSM+KSM1
      IF(N.EQ.KCO)KCO=KCO+KCO1
C**** UPDATE ARRAYS.
      DO 1150 K=1,KK
      DO 1150 J=1,JJ
      DO 1150 I=1,II
      QO(I,J,K)=0.0
      QW(I,J,K)=0.0
      QG(I,J,K)=0.0
      PN(I,J,K)=P(I,J,K)
      SON(I,J,K)=SO(I,J,K)
      SWN(I,J,K)=SW(I,J,K)
      SGN(I,J,K)=SG(I,J,K)
1150  CONTINUE
C
591   FORMAT(///T5,10('*'),' WELL REPORT FOR ALL ACTIVE WELLS ',
     &4X,'ELAPSED TIME =',F11.6,' DAYS FROM BEGINNING OF SIMULATION ',
     &10('*'),//)
C
      IF(WOR.GT.WORMAX) GO TO 1002
      IF(GOR.GT.GORMAX) GO TO 1003
      IF(PAVG.LT.PAMIN) GO TO 1004
      IF(PAVG.GT.PAMAX) GO TO 1005
1000  CONTINUE
C
1002  WRITE(NO,2002)
2002  FORMAT(/T15,'MAXIMUM WOR HAS BEEN EXCEEDED --- SIMULATION',
     &' IS BEING TERMINATED'//)
      GO TO 1010
1003  WRITE(NO,2003)
2003  FORMAT(/T15,'MAXIMUM GOR HAS BEEN EXCEEDED --- SIMULATION',
     &' IS BEING TERMINATED'//)
      GO TO 1010
1004  WRITE(NO,2004)
2004  FORMAT(/T15,'MINIMUM AVERAGE RESERVOIR PRESSURE WAS NOT',
     &' ACHIEVED --- SIMULATION IS BEING TERMINATED'//)
      GO TO 1010
1005  WRITE(NO,2005)
2005  FORMAT(/T15,'MAXIMUM AVERAGE RESERVOIR PRESSURE
     &HAS BEEN EXCEEDED --- SIMULATION IS BEING TERMINATED'//)
1010  CLOSE(UNIT=NI)
      CLOSE(UNIT=NO)
      IF(NRESTART.EQ.1) CALL RESTART(INFILE,RSTFILE,NX,NY,NZ,
     &  P,SO,SW,DPP,DPSO,DPSW)
      CLOSE(NI)
      CLOSE(NO,STATUS='KEEP')
      CLOSE(7,STATUS='KEEP')
      STOP
      END

C     WAIT.FOR
      SUBROUTINE WAIT(DELAY)
C       DELAY IN 1/100 th's of a second
      INTEGER*4 DELAY,NOW,NEXT_TIME
      CALL TIMER ( NOW )
      NEXT_TIME = NOW + DELAY
        DO WHILE ( NOW .LE. NEXT_TIME )
           CALL TIMER ( NOW )
	END DO
      RETURN
	END
C	VARINIT.FOR
C
	SUBROUTINE VARINIT(NX,NY,NZ,TX,TY,TZ,AW,AS,AT,AE,AN,AB,E,B,
     &   GOWT,GWWT,GGWT,OW,OE,WW,WE,OS,ON,WS,WN,OT,OB,WT,WB)

      DIMENSION TX(NX+1,NY,NZ),TY(NX,NY+1,NZ),TZ(NX,NY,NZ+1)
      DIMENSION AW(NX,NY,NZ),AE(NX,NY,NZ),AN(NX,NY,NZ),
     & AS(NX,NY,NZ),AB(NX,NY,NZ),AT(NX,NY,NZ),E(NX,NY,NZ),B(NX,NY,NZ) 
      DIMENSION GOWT(NX,NY,NZ),GWWT(NX,NY,NZ),GGWT(NX,NY,NZ),
     & OW(NX+1,NY,NZ),OE(NX+1,NY,NZ),WW(NX+1,NY,NZ),WE(NX+1,NY,NZ),
     & OS(NX,NY+1,NZ),ON(NX,NY+1,NZ),WS(NX,NY+1,NZ),WN(NX,NY+1,NZ),
     & OT(NX,NY,NZ+1),OB(NX,NY,NZ+1),WT(NX,NY,NZ+1),WB(NX,NY,NZ+1)
C      DIMENSION VP(NX,NY,NZ)
C      DIMENSION UM(NX+1),AZL(NX+1),BZL(NX+1),CZL(NX+1),
C     &    DZL(NX+1),UZL(NX+1)
C      DIMENSION BETA(NX+1),GAMMA(NX+1),W(NX+1)
	DO 10 I=1,NX
	DO 10 J=1,NY
	DO 10 K=1,NZ
	AW(I,J,K)=0.0
	AE(I,J,K)=0.0
	AN(I,J,K)=0.0
	AS(I,J,K)=0.0
	AB(I,J,K)=0.0
	AT(I,J,K)=0.0
	E(I,J,K) =0.0
	B(I,J,K) =0.0
	GOWT(I,J,K)=0.0
	GWWT(I,J,K)=0.0
C        VP(I,J,I)=0.0
10	GGWT(I,J,K)=0.0
	DO 20 I=1,NXP
C         UM(I)=0.0
C         AZL(I)=0.0
C         BZL(I)=0.0
C         CZL(I)=0.0
C         DZL(I)=0.0
C         UZL(I)=0.0
C         BETA(I)=0.0
C         GAMMA(I)=0.0
C         W(I)=0.0
	DO 20 J=1,NY
	DO 20 K=1,NZ
	TX(I,J,K)=0.0
	OW(I,J,K)=0.0
	OE(I,J,K)=0.0
	WW(I,J,K)=0.0
20	WE(I,J,K)=0.0
	DO 30 I=1,NX
	DO 30 J=1,NYP
	DO 30 K=1,NZ
	TY(I,J,K)=0.0
	OS(I,J,K)=0.0
	ON(I,J,K)=0.0
	WS(I,J,K)=0.0
30	WN(I,J,K)=0.0
	DO 40 I=1,NX
	DO 40 J=1,NY
	DO 40 K=1,NZP
	TZ(I,J,K)=0.0
	OT(I,J,K)=0.0
	OB(I,J,K)=0.0
	WT(I,J,K)=0.0
40	WB(I,J,K)=0.0
C	IM=NX	! IS THIS IM EQUAL TO I MAX
C	JM=NY
C	KM=NZ
	RETURN
	END
C	EXTEN.FOR
C
C	******************************RJH*********************************
C
C	AN EXTENSION MANAGER  VERSION 2
C
	SUBROUTINE EXTEN ( OLDFIL, NEWFIL, EXT )
        PARAMETER (NL=20)
	CHARACTER PERIOD*1, SPACE*1, EXT*3, OLDFIL*NL, NEWFIL*NL
	DATA PERIOD /'.'/, SPACE /' '/
	DO 10 I = 1,NL
	 IF ( OLDFIL(I:I) .NE. PERIOD ) THEN
	   NEWFIL(I:I)   = OLDFIL(I:I)
	 ELSE
	   NEWFIL(I:I+3) = OLDFIL(I:I)//EXT
	   RETURN
	 ENDIF
10	CONTINUE
        DO 20 I=1,NL
         IF ( OLDFIL(I:I) .NE. SPACE ) THEN
           NEWFIL(I:I) = OLDFIL(I:I)
         ELSE IF ( I+4 .LE. NL ) THEN
           NEWFIL(I:I+4) = PERIOD//EXT
           RETURN
         ELSE
           WRITE(*,'(26H FILE NAME LENGTH EXCEEDS ,I2 )') NL
           STOP 'FERR' 
         ENDIF
20      CONTINUE
	RETURN
	END

C     GAUGE.FOR
C
C     A THERMOMETER OF PROGRAM PROGESS DESIGNED TO GIVE
C     THE USER AN IDEA OF HOW FAR HIS OR HER PROGRAM IS
C     ADVANCED IN ANY PRE-SELECTED TIME PERIOD
C
      SUBROUTINE GAUGE( NGAUGE, ICELL, PERIOD, TMAX )
      CHARACTER*1 STEP
      DATA STEP /'Û'/
C     WHERE:
C         NGAUGE = THE ACTUAL NUMBER STEPS IN THE LOOP
C         ILAST  = THE NUMBER OF THE LAST VISIBLE THERMOMETER CELL
C         PERIOD = IS THE EQUIVALENT TIME PERIOD (I.E. DAYS)
C         TMAX   = IS THE MAXIMUM NUMBER DAYS PERMITTED IN THE PROGRAM
C
      IF(NGAUGE .LE. 1) THEN
      WRITE(*,100) TMAX
100   FORMAT(1X,T22,'SCALE OF PROGRESS',T49,'Maximum =',F7.2,' Days',
     &    /,7X,'0',T22,'25   Percent',T37,'50',T52,'75',T66,'100',
     &    /,7X,1Hร,3( 4(3Hฤฤย),3Hฤฤล ),4(3Hฤฤย),3Hฤฤด,/,7X,1HÛ)
      ENDIF
      NGAUGE=NGAUGE+1
      IF(PERIOD .LE. TMAX) THEN
       INEXT = 60*(PERIOD/TMAX)
         IF (INEXT .GT. ICELL) THEN 
20        WRITE(*,200) STEP
200       FORMAT(1H&,A1)
          ICELL = ICELL + 1
	  IF(INEXT .GT. ICELL) GO TO 20
         ENDIF
      ENDIF
      RETURN
      END
C *********************************************************
	subroutine ClrScr
C
C DESCRIPTION:
C    This procedure uses the ANSI driver to clear the screen.
C
C SAMPLE CALL:
C    CALL ClrScr
C
C SOURCE FILE:
C    CLRSCR.FOR
C
C NOTES:
C    The cursor is placed in the home the position following
C    the procedure call.
C
C ----------------------------------------------------------

      write (*, 1) 27, 27
 1    format (1X, A1, '[H', A1, '[J')
      return
      end

      BLOCK DATA
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      COMMON /IO/ NI,NO
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      COMMON /ADD1/ IM,JM,KM,ETI,FT,FTMAX,DELT
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      DATA NI,NO /5,6/
      DATA CUMO,CUMW,CUMG/NWZ*0.0,NWZ*0.0,NWZ*0.0/
      DATA ETI,FT,FTMAX/0.0,0.0,0.0/
      DATA COP,CWP,CGP,CWI,CGI/0.D0,0.D0,0.D0,0.D0,0.D0/

      END
C
C
C
      
      BLOCK DATA FILNAM
C
C============================================BLOCK DATA FILNAM
C
C    CHARACTER VARIABLE INITIALIZATION FOR GRAPHICS.
C
C=============================================================
C
       INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      CHARACTER*2 FNAME,MNAME,TAG
      CHARACTER*4 LNAME
      COMMON /CHRCTR/ FNAME,MNAME(NUM),LNAME,TAG
      COMMON /TTEST/ NUMPRD,SONTVL
      COMMON /IO2/ NO2
      DATA NO2/7/
      DATA NUMPRD/1/
      DATA FNAME,LNAME/'SO','.DAT'/
      DATA MNAME/'01','02','03','04','05','06','07','08','09',
     & '10','11','12','13','14','15','16','17','18','19','20',
     & '21','22','23','24','25','26','27','28','29','30','31',
     & '32','33','34','35','36','37','38','39','40'/
      END
C
C
C
      SUBROUTINE INTPVT(BPT,RM,X,Y,N,XO,YO)
C
C=================================================================INTPVT
C
C   PERFORM LINEAR INTERPOLATION WITH FLOATING BUBBLE POINT CALCULATION
C   OF UNDERSATURATED OIL PVT PROPERTIES.
C
C=======================================================================
C
	INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      DIMENSION X(NTE),Y(NTE)
C      WRITE(6,1)(BPT,RM,X(I),Y(I),N,XO,YO,I=1,N)
C1    FORMAT(5X,'B = ',4F15.5,I15,2F15.5)
      IF(XO.GT.BPT)GO TO 100
      IF(XO.GE.X(N))YO=Y(N)
      IF(XO.GE.X(N))RETURN
      DO 10 I=2,N
      IF(XO.GE.X(I))GO TO 10
      YO=Y(I-1)+(XO-X(I-1))*(Y(I)-Y(I-1))/(X(I)-X(I-1))
      RETURN
   10 CONTINUE
  100 CONTINUE
      DO 200 I=2,N
      IF(BPT.GE.X(I))GO TO 200
      YOBP=Y(I-1)+(BPT-X(I-1))*(Y(I)-Y(I-1))/(X(I)-X(I-1))
      YO=YOBP+RM*(XO-BPT)
      RETURN
  200 CONTINUE
      END
C
C
C
      SUBROUTINE INTERP(X,Y,N,XO,YO)
C
C=================================================================INTERP
C
C   PERFORM LINEAR INTERPOLATION.
C
C=======================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      DIMENSION X(NTE),Y(NTE)
      IF(XO.GE.X(N)) YO=Y(N)
      IF(XO.GE.X(N)) RETURN
      DO 10 I=2,N
      IF(XO.GE.X(I)) GO TO 10
      YO=Y(I-1) +(XO-X(I-1))*(Y(I)-Y(I-1))/(X(I)-X(I-1))
      RETURN
  10  CONTINUE
      END
C
C      SUBROUTINE GRID1(NX,NY,NZ)
      SUBROUTINE GRID1(NX,NY,NZ,EL,DX,DY,DZ,SUM,VAREL,RDXL,RDYL,RDZL)
C
C==================================================================GRID1
C
C   READ GRID BLOCK SIZES AND COMPUTE NODE ELEVATIONS
C
C=======================================================================
C

      CHARACTER*2 TITLE(80)
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      DIMENSION EL(NX,NY,NZ)
      DIMENSION DX(NX,NY,NZ),DY(NX,NY,NZ),DZ(NX,NY,NZ)
      COMMON /SWTCH/ NGRSW,NTRSW,NRESTART,NTS
      COMMON /COUNT/ N1READ,N2READ
      DIMENSION SUM(NX,NY),VAREL(NX,NY),RDXL(NX),RDYL(NY),RDZL(NZ)
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,70) (TITLE(IH),IH=1,40)
      WRITE(14,70) (TITLE(IH),IH=1,40)
      WRITE(14,3003)
      WRITE(14,3004)
      WRITE(14,3001)
      WRITE(14,3002)
 3003 FORMAT(//,'                                                    ',
     & '                    GAS                     WATER')
 3004 FORMAT('         PROD      AVG      OIL       CUM        GAS   ',
     & '       CUM       OIL      WATER   CUM     OIL')
 3001 FORMAT('    N    TIME     PRESS     RATE      OIL        RATE  ',
     &  '       GAS      RATIO     RATE   WATER   RATIO')
 3002 FORMAT(' (T.S.)  (DAY)    (PSI)    (B/D)     (BBL)      (MCF/D)',
     &  '      (MCF)    (MCF/B)   (B/D)   (BBL)   (B/B)',/)
      READ(NI,*) II,JJ,KK
      N1READ=N1READ+1
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      READ(NI,*) NGRSW,NTRSW,NRESTART,NTS
      N1READ=N1READ+1
      IF(NGRSW.EQ.1) CALL GRINT
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
C
C*****READ INPUT CODES FOR DX,DY,DZ
C
      READ(NI,*)KDX,KDY,KDZ
      N1READ=N1READ+1
C
C*****ESTABLISH GRID BLOCK LENGTH (DX) DISTRIBUTION
C
      IF(KDX.GE.0)GO TO 180
      READ(NI,*)DXC
      N1READ=N1READ+1
      DO 175 K=1,KK
      DO 175 J=1,JJ
      DO 175 I=1,II
175   DX(I,J,K)=DXC
      WRITE(NO,56)
      WRITE(NO,29)DXC
      GO TO 195
  180 IF(KDX.GT.0)GO TO 185
      READ(NI,*)(RDXL(I),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
C 998  format(10f8.1)
      DO 187 K=1,KK
      DO 187 J=1,JJ
      DO 187 I=1,II
  187 DX(I,J,K)=RDXL(I)
      DO 182 I=1,II
  182 WRITE(NO,511)I,RDXL(I)
      GO TO 195
185   WRITE(NO,43)
      K=1
      WRITE(NO,38)K
      DO 190 J=1,JJ
      READ(NI,*)  (DX(I,J,K),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
  190 WRITE(NO,72)(DX(I,J,K),I=1,II)
      DO 194 K=2,KK
      WRITE(NO,38) K
      DO 194 J=1,JJ
      DO 193 I=1,II
  193 DX(I,J,K)=DX(I,J,1)
  194 WRITE(NO,72) (DX(I,J,K),I=1,II)
195   CONTINUE
      WRITE(NO,56)
C
C*****ESTABLISH GRID BLOCK LENGTH (DY) DISTRIBUTION
C
      IF(KDY.GE.0)GO TO 200
      READ(NI,*)DYC
      N1READ=N1READ+1
      DO 202 K=1,KK
      DO 202 J=1,JJ
      DO 202 I=1,II
202   DY(I,J,K)=DYC
      WRITE(NO,56)
      WRITE(NO,33)DYC
      GO TO 220
  200 IF(KDY.GT.0)GO TO 207
      READ(NI,*)(RDYL(J),J=1,JJ)
      N1READ=N1READ+1+INT(JJ/10)
      IF(MOD(JJ,10).EQ.0)N1READ=N1READ-1
      DO 205 K=1,KK
      DO 205 J=1,JJ
      DO 205 I=1,II
  205 DY(I,J,K)=RDYL(J)
      DO 210 J=1,JJ
  210 WRITE(NO,512)J,RDYL(J)
      GO TO 220
207   WRITE(NO,47)
      K=1
      WRITE(NO,38)K
      DO 215 I=1,II
      READ(NI,*)(DY(I,J,K),J=1,JJ)
      N1READ=N1READ+1+INT(JJ/10)
      IF(MOD(JJ,10).EQ.0)N1READ=N1READ-1
  215 WRITE(NO,72)(DY(I,J,K),J=1,JJ)
      DO 214 K=2,KK
      WRITE(NO,38) K
      DO 214 J=1,JJ
      DO 213 I=1,II
  213 DY(I,J,K)=DY(I,J,1)
  214 WRITE(NO,72) (DY(I,J,K),I=1,II)
220   CONTINUE
      WRITE(NO,56)
C
C*****ESTABLISH GRID BLOCK LENGTH (DZ) DISTRIBUTION
C
      IF(KDZ.GE.0)GO TO 225
      READ(NI,*)DZC
      N1READ=N1READ+1
      DO 230 K=1,KK
      DO 230 J=1,JJ
      DO 230 I=1,II
230   DZ(I,J,K)=DZC
      WRITE(NO,56)
      WRITE(NO,36)DZC
      GO TO 245
  225 IF(KDZ.GT.0)GO TO 232
      READ(NI,*)(RDZL(K),K=1,KK)
      N1READ=N1READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N1READ=N1READ-1
      DO 235 K=1,KK
      DO 235 J=1,JJ
      DO 235 I=1,II
  235 DZ(I,J,K)=RDZL(K)
      DO 237 K=1,KK
  237 WRITE(NO,513)K,RDZL(K)
      GO TO 245
232   WRITE(NO,48)
      DO 240 K=1,KK
      WRITE(NO,38)K
      DO 242 J=1,JJ
      READ(NI,*)(DZ(I,J,K),I=1,II)
      N1READ=N1READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N1READ=N1READ-1
  242 WRITE(NO,72)(DZ(I,J,K),I=1,II)
240   CONTINUE
245   CONTINUE
      WRITE(NO,56)
C
C********GRID BLOCK LENGTH MODIFICATIONS
C
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      READ(NI,*) NUMDX,NUMDY,NUMDZ,IDCODE
      N1READ=N1READ+1
      IF(NUMDX.EQ.0) GO TO 8530
      WRITE(NO,31)
      DO 275 L=1,NUMDX
      READ(NI,*)I,J,K,DX(I,J,K)
      N1READ=N1READ+1
  275 WRITE(NO,32)I,J,K,DX(I,J,K)
      IF(IDCODE.NE.1)GO TO 8530
      WRITE(NO,43)
      DO 853 K=1,KK
      WRITE(NO,38)K
      DO 854 J=1,JJ
  854 WRITE(NO,72)(DX(I,J,K),I=1,II)
853   CONTINUE
8530  CONTINUE
      IF(NUMDY.EQ.0) GO TO 8550
      WRITE(NO,34)
      DO 276 L=1,NUMDY
      READ(NI,*)I,J,K,DY(I,J,K)
      N1READ=N1READ+1
  276 WRITE(NO,32)I,J,K,DY(I,J,K)
      IF(IDCODE.NE.1) GO TO 8550
      WRITE(NO,47)
      DO 855 K=1,KK
      WRITE(NO,38)K
      DO 856 J=1,JJ
  856 WRITE(NO,72)(DY(I,J,K),I=1,II)
855   CONTINUE
8550  CONTINUE
      IF(NUMDZ.EQ.0) GO TO 8570
      WRITE(NO,37)
      DO 277 L=1,NUMDZ
      READ(NI,*)I,J,K,DZ(I,J,K)
      N1READ=N1READ+1
  277 WRITE(NO,32)I,J,K,DZ(I,J,K)
      IF(IDCODE.NE.1) GO TO 8570
      WRITE(NO,48)
      DO 857 K=1,KK
      WRITE(NO,38)K
      DO 858 J=1,JJ
  858 WRITE(NO,72)(DZ(I,J,K),I=1,II)
857   CONTINUE
8570  CONTINUE
C
C*****ESTABLISH NODE MID-POINT ELEVATION.
C
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      READ(NI,*)KEL
      N1READ=N1READ+1
      IF(KEL.EQ.1)GO TO 920
      READ(NI,*)ELEV
      N1READ=N1READ+1
      DO 910 J=1,JJ
      DO 910 I=1,II
  910 VAREL(I,J)=ELEV
  920 IF(KEL.NE.1)GO TO 930
      DO 922 J=1,JJ
      READ(NI,*)(VAREL(I,J),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
  922 CONTINUE
  930 CONTINUE
      DO 923 I=1,II
      DO 923 J=1,JJ
  923 SUM(I,J)=0.
      DO 926 K=1,KK
      DO 926 J=1,JJ
      DO 926 I=1,II
      DEL=SUM(I,J)+DZ(I,J,K)*0.5
      EL(I,J,K)=VAREL(I,J)+DEL
  926 SUM(I,J)=DZ(I,J,K)+SUM(I,J)
      WRITE(NO,390)
      DO 600 K=1,KK
      WRITE(NO,38)K
      DO 600 J=1,JJ
      WRITE(NO,72)(EL(I,J,K),I=1,II)
600   CONTINUE
69    FORMAT(40A2)
56    FORMAT(//)
70    FORMAT(T24,82('*'),/,T24,('*'),T105,('*'),/,
     & T24,('*'),40A2,T105,('*'),/,T24,('*'),T105,
     & ('*'),/,T24,82('*'),//)
72    FORMAT(1X,20F6.0)
29    FORMAT(T15,'GRID BLOCK LENGTH (DX) IS INITIALLY',
     &' SET AT',F10.4,' FOR ALL NODES'//)
31    FORMAT(//T15,'********GRID BLOCK LENGTH (DX) NODE MODIFICATIONS',
     &'**********',//T15,'     I    J    K         NEW DX VALUE')
32    FORMAT(15X,3I5,5X,F14.4)
33    FORMAT(T15,'GRID BLOCK WIDTH  (DY) IS INITIALLY',
     &' SET AT',F10.4,' FOR ALL NODES'//)
34    FORMAT(//T15,'********GRID BLOCK WIDTH  (DY) NODE MODIFICATIONS',
     &'**********',//T15,'     I    J    K         NEW DY VALUE')
36    FORMAT(T15,'GRID BLOCK DEPTH  (DZ) IS INITIALLY',
     &' SET AT',F10.4,' FOR ALL NODES'//)
37    FORMAT(//T15,'********GRID BLOCK DEPTH  (DZ) NODE MODIFICATIONS',
     &'**********',//T15,'     I    J    K         NEW DZ VALUE')
38    FORMAT(/1X,'K =',I2/)
390   FORMAT(//T15,'********** NODE MIDPOINT ELEVATIONS **********'/)
43    FORMAT(//T15,'********GRID BLOCK LENGTH (DX) DISTRIBUTION********'
     &/)
47    FORMAT(//T15,'********GRID BLOCK WIDTH  (DY) DISTRIBUTIO*****'
     &/)
48    FORMAT(//T15,'********GRID BLOCK DEPTH  (DZ) DISTRIBUTION********'
     &/)
  511 FORMAT(T15,'GRID SIZE (DX) IN COLUMN',I5,' IS INITIALLY SET AT'
     &,F8.2,' FOR ALL NODES',/)
  512 FORMAT(T15,'GRID SIZE (DY) IN ROW   ',I5,' IS INITIALLY SET AT'
     &,F8.2,' FOR ALL NODES',/)
  513 FORMAT(T15,'GRID SIZE (DZ) IN LAYER ',I5,' IS INITIALLY SET AT'
     &,F8.2,' FOR ALL NODES',/)
      RETURN
      END
C      SUBROUTINE PARM1(NX,NY,NZ)
      SUBROUTINE PARM1(NX,NY,NZ,VP,PV,KX,KY,KZ,DX,DY,DZ)
C
C==================================================================PARM1
C
C   READ POROSITY AND PERMEABILITY ARRAYS.
C
C       When KKX = -1, a single Kx is read for all grid blocks
C            KKX =  0, a constant Kx is read for each of KK layers
C            KKX =  1, a seperate Kx is read for each II x JJ x KK grid block
C            KKX =  2, a constant Kx is read for each II x KK grid block
C       The same rule is applied to KKY and KKZ
C
C=======================================================================
C
	INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      DIMENSION DX(NX,NY,NZ),DY(NX,NY,NZ),DZ(NX,NY,NZ)
      REAL KXC,KYC,KZC
      CHARACTER*2 TITLE(80)
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      DIMENSION VP(NX,NY,NZ)
      REAL KX(NX,NY,NZ),KY(NX,NY,NZ),KZ(NX,NY,NZ)
      COMMON /COUNT/ N1READ,N2READ
      DIMENSION RPHL(NZ1),RKXL(NZ1),RKYL(NZ1),RKZL(NZ1)
      DIMENSION XKX(NX1),YKX(NX1),ZKX(NX1)
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
C
C*****READ INPUT CODES FOR PHI,KX,KY,KZ
C
      READ(NI,*)KPH,KKX,KKY,KKZ
      N1READ=N1READ+1
C
C*****ESTABLISH POROSITY (PHI) DISTRIBUTION
C
      IF(KPH.GE.0)GO TO 135
      READ(NI,*)PHIC
      N1READ=N1READ+1
      DO 140 K=1,KK
      DO 140 J=1,JJ
      DO 140 I=1,II
140   VP(I,J,K)=PHIC
      WRITE(NO,56)
      WRITE(NO,26)PHIC
      GO TO 165
  135 IF(KPH.GT.0)GO TO 145
      READ(NI,*)(RPHL(K),K=1,KK)
      N1READ=N1READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N1READ=N1READ-1
      DO 550 K=1,KK
      DO 550 J=1,JJ
      DO 550 I=1,II
  550 VP(I,J,K)=RPHL(K)
      DO 560 K=1,KK
  560 WRITE(NO,510)K,RPHL(K)
      GO TO 165
145   WRITE(NO,39)
      DO 160 K=1,KK
      WRITE(NO,38)K
      DO 155 J=1,JJ
      READ(NI,*)(VP(I,J,K),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
  155 WRITE(NO,73)(VP(I,J,K),I=1,II)
160   CONTINUE
165   CONTINUE
      WRITE(NO,56)
C
C*****ESTABLISH PERMEABILITY (KX) DISTRIBUTION
C
      IF(KKX.GE.0)GO TO 180
C *** A single constant Kx is read for all grid blocks  
      READ(NI,*)KXC
      N1READ=N1READ+1
      DO 175 K=1,KK
      DO 175 J=1,JJ
      DO 175 I=1,II
175   KX(I,J,K)=KXC
      WRITE(NO,56)
      WRITE(NO,29)KXC
      GO TO 195
C *** a constant value Kx is read for each of the KK layers
  180 IF(KKX.GT.0)GO TO 185
      READ(NI,*)(RKXL(K),K=1,KK)
      N1READ=N1READ+1+INT(KK/10)
      IF(MOD(KK,1).EQ.0)N1READ=N1READ-1
      DO 187 K=1,KK
      DO 187 J=1,JJ
      DO 187 I=1,II
  187 KX(I,J,K)=RKXL(K)
      DO 182 K=1,KK
  182 WRITE(NO,511)K,RKXL(K)
      GO TO 195
C *** a sperate Kx value is read for each (II x JJ x KK) grid block
185   IF(KKX.GT.1)GO TO 193
      WRITE(NO,43)
      DO 192 K=1,KK
      WRITE(NO,38)K
      DO 190 J=1,JJ
      READ(NI,*)(KX(I,J,K),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
  190 WRITE(NO,72)(KX(I,J,K),I=1,II)
192   CONTINUE
      GO TO 195
C *** a Kx value is read for each (II) column and each (KK) layer : (II x KK)
193   WRITE(NO,43)
      DO 194 K=1,KK
      WRITE(NO,38)K
      READ(NI,*)(XKX(I),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
      DO 196 J=1,JJ
      DO 197 I=1,II
  197  KX(I,J,K)=XKX(I)
  196 WRITE(NO,72)(KX(I,J,K),I=1,II)
194   CONTINUE
195   CONTINUE
      WRITE(NO,56)
C
C*****ESTABLISH PERMEABILITY (KY) DISTRIBUTION
C
      IF(KKY.GE.0)GO TO 200
      READ(NI,*)KYC
      N1READ=N1READ+1
      DO 202 K=1,KK
      DO 202 J=1,JJ
      DO 202 I=1,II
202   KY(I,J,K)=KYC
      WRITE(NO,56)
      WRITE(NO,33)KYC
      GO TO 220
  200 IF(KKY.GT.0)GO TO 207
      READ(NI,*)(RKYL(K),K=1,KK)
      N1READ=N1READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N1READ=N1READ-1
      DO 205 K=1,KK
      DO 205 J=1,JJ
      DO 205 I=1,II
  205 KY(I,J,K)=RKYL(K)
      DO 210 K=1,KK
  210 WRITE(NO,512)K,RKYL(K)
      GO TO 220
 207  IF(KKY.GT.1)GO TO 216
      WRITE(NO,47)
      DO 212 K=1,KK
      WRITE(NO,38)K
      DO 215 J=1,JJ
      READ(NI,*)(KY(I,J,K),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
  215 WRITE(NO,72)(KY(I,J,K),I=1,II)
212   CONTINUE
      GO TO 220
C *** a Ky value is read for each (II) column and each (KK) layer : (II x KK)
216   WRITE(NO,47)
      DO 217 K=1,KK
      WRITE(NO,38)K
      READ(NI,*)(YKX(I),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
      DO 218 J=1,JJ
      DO 219 I=1,II
  219  KY(I,J,K)=YKX(I)
  218 WRITE(NO,72)(KY(I,J,K),I=1,II)
217   CONTINUE
220   CONTINUE
      WRITE(NO,56)
C
C*****ESTABLISH PERMEABILITY (KZ) DISTRIBUTION
C
      IF(KKZ.GE.0)GO TO 225
      READ(NI,*)KZC
      N1READ=N1READ+1
      DO 230 K=1,KK
      DO 230 J=1,JJ
      DO 230 I=1,II
230   KZ(I,J,K)=KZC
      WRITE(NO,56)
      WRITE(NO,36)KZC
      GO TO 245
  225 IF(KKZ.GT.0)GO TO 232
      READ(NI,*)(RKZL(K),K=1,KK)
      N1READ=N1READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N1READ=N1READ-1
      DO 235 K=1,KK
      DO 235 J=1,JJ
      DO 235 I=1,II
  235 KZ(I,J,K)=RKZL(K)
      DO 237 K=1,KK
  237 WRITE(NO,513)K,RKZL(K)
      GO TO 245
 232  IF(KKZ.GT.1)GO TO 243
      WRITE(NO,48)
      DO 240 K=1,KK
      WRITE(NO,38)K
      DO 242 J=1,JJ
      READ(NI,*)(KZ(I,J,K),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
  242 WRITE(NO,72)(KZ(I,J,K),I=1,II)
240   CONTINUE
      GO TO 245
C *** a Kz value is read for each (II) column and each (KK) layer : (II x KK)
243   WRITE(NO,48)
      DO 244 K=1,KK
      WRITE(NO,38)K
      READ(NI,*)(ZKX(I),I=1,II)
      N1READ=N1READ+1+INT(II/10)
      IF(MOD(II,10).EQ.0)N1READ=N1READ-1
      DO 246 J=1,JJ
      DO 247 I=1,II
  247  KZ(I,J,K)=ZKX(I)
  246 WRITE(NO,72)(KZ(I,J,K),I=1,II)
244   CONTINUE
245   CONTINUE
      WRITE(NO,56)
C
C**********POROSITY AND PERMEABILITY MODIFICATIONS
C
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      READ(NI,*) NUMP,NUMKX,NUMKY,NUMKZ,IPCODE
      N1READ=N1READ+1
      IF(NUMP.EQ.0) GO TO 8510
      WRITE(NO,27)
      DO 274 L=1,NUMP
      READ(NI,*)I,J,K,VP(I,J,K)
      N1READ=N1READ+1
  274 WRITE(NO,32)I,J,K,VP(I,J,K)
      IF(IPCODE.NE.1) GO TO 8510
      WRITE(NO,39)
      DO 851 K=1,KK
      WRITE(NO,38)K
      DO 852 J=1,JJ
  852 WRITE(NO,73)(VP(I,J,K),I=1,II)
  851 CONTINUE
8510  CONTINUE
      IF(NUMKX.EQ.0) GO TO 8530
      WRITE(NO,31)
      DO 275 L=1,NUMKX
      READ(NI,*)I,J,K,KX(I,J,K)
      N1READ=N1READ+1
  275 WRITE(NO,32)I,J,K,KX(I,J,K)
      IF(IPCODE.NE.1)GO TO 8530
      WRITE(NO,43)
      DO 853 K=1,KK
      WRITE(NO,38)K
      DO 854 J=1,JJ
  854 WRITE(NO,72)(KX(I,J,K),I=1,II)
853   CONTINUE
8530  CONTINUE
      IF(NUMKY.EQ.0) GO TO 8550
      WRITE(NO,34)
      DO 276 L=1,NUMKY
      READ(NI,*)I,J,K,KY(I,J,K)
      N1READ=N1READ+1
  276 WRITE(NO,32)I,J,K,KY(I,J,K)
      IF(IPCODE.NE.1) GO TO 8550
      WRITE(NO,47)
      DO 855 K=1,KK
      WRITE(NO,38)K
      DO 856 J=1,JJ
  856 WRITE(NO,72)(KY(I,J,K),I=1,II)
855   CONTINUE
8550  CONTINUE
      IF(NUMKZ.EQ.0) GO TO 8570
      WRITE(NO,37)
      DO 277 L=1,NUMKZ
      READ(NI,*)I,J,K,KZ(I,J,K)
      N1READ=N1READ+1
  277 WRITE(NO,32)I,J,K,KZ(I,J,K)
      IF(IPCODE.NE.1) GO TO 8570
      WRITE(NO,48)
      DO 857 K=1,KK
      WRITE(NO,38)K
      DO 858 J=1,JJ
  858 WRITE(NO,72)(KZ(I,J,K),I=1,II)
857   CONTINUE
8570  CONTINUE
C
69    FORMAT(40A2)
56    FORMAT(//)
72    FORMAT(1X,20F6.2)
73    FORMAT(1X,20F6.4)
26    FORMAT(T15,'POROSITY (PHI) IS INITIALLY SET AT',F8.4,' FOR ALL',
     &' NODES'//)
27    FORMAT(//T15,'**********POROSITY NODE MODIFICATIONS**********',
     &//T15,'     I    J    K         NEW PHI VALUE')
29    FORMAT(T15,'PERMEABILITY (KX) IS INITIALLY',
     &' SET AT',F10.4,' FOR ALL NODES'//)
31    FORMAT(//T15,'**********PERMEABILITY (KX) NODE MODIFICATIONS',
     &'**********',//T15,'     I    J    K         NEW KX VALUE')
32    FORMAT(15X,3I5,5X,F14.4)
33    FORMAT(T15,'PERMEABILITY (KY) IS INITIALLY',
     &' SET AT',F10.4,' FOR ALL NODES'//)
34    FORMAT(//T15,'**********PERMEABILITY (KY) NODE MODIFICATIONS',
     &'**********',//T15,'     I    J    K         NEW KY VALUE')
36    FORMAT(T15,'PERMEABILITY (KZ) IS INITIALLY',
     &' SET AT',F10.4,' FOR ALL NODES'//)
37    FORMAT(//T15,'**********PERMEABILITY (KZ) NODE MODIFICATIONS',
     &'**********',//T15,'     I    J    K         NEW KZ VALUE')
38    FORMAT(/1X,'K =',I2/)
39    FORMAT(//T15,'**********POROSITY DISTRIBUTION FOLLOWS**********'/)
390   FORMAT(//T15,'********** NODE MIDPOINT ELEVATIONS **********'/)
43    FORMAT(//T15,'**********PERMEABILITY (KX) DISTRIBUTION**********'
     &/)
445   FORMAT(13F10.2)
47    FORMAT(//T15,'**********PERMEABILITY (KY) DISTRIBUTION**********'
     &/)
48    FORMAT(//T15,'**********PERMEABILITY (KZ) DISTRIBUTION**********'
     &/)
510   FORMAT(//T15,'POROSITY IN LAYER',I5,' IS INITIALLY SET AT',F8.5,
     &' FOR ALL NODES',/)
511   FORMAT(T15,'PERMEABILITY (KX) IN LAYER',I5,' IS INITIALLY SET AT'
     &,F8.2,' FOR ALL NODES',/)
512   FORMAT(T15,'PERMEABILITY (KY) IN LAYER',I5,' IS INITIALLY SET AT'
     &,F8.2,' FOR ALL NODES',/)
513   FORMAT(T15,'PERMEABILITY (KZ) IN LAYER',I5,' IS INITIALLY SET AT'
     &,F8.2,' FOR ALL NODES',/)
      RETURN
      END
C      SUBROUTINE TRAN1(NX,NY,NZ)
      SUBROUTINE TRAN1(NX,NY,NZ,A1,A2,A3,KX,KY,KZ,TX,TY,TZ,DX,DY,DZ)
C
C==================================================================TRAN1
C
C   COMPUTE GEOMETRIC PART OF TRANSMISSIBILITIES.
C
C======================================================================
C
      
      CHARACTER*2 TITLE(80)
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      DIMENSION A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ)
      REAL KX(NX,NY,NZ),KY(NX,NY,NZ),KZ(NX,NY,NZ)
      DIMENSION TX(NX+1,NY,NZ),TY(NX,NY+1,NZ),TZ(NX,NY,NZ+1)
      DIMENSION DX(NX,NY,NZ),DY(NX,NY,NZ),DZ(NX,NY,NZ)
      COMMON /CODE/ KSM1,KSN1,KCO1,NN,FACT1,FACT2,TMAX,KSOL,MITER,
     & OMEGA,TOL,TOL1,KSN,KSM,KCO,KTR,KCOFF,DSMAX,DPMAX,WORMAX,
     & GORMAX,PAMIN,PAMAX
      COMMON /SWTCH/ NGRSW,NTRSW,NRESTART,NTS
      COMMON /COUNT/ N1READ,N2READ
C
      DO 30 K=1,KK
      DO 30 J=1,JJ
      DO 30 I=1,II
      FACX=1.0
      FACY=1.0
      FACZ=1.0
      IF(I.GT.1.AND.I.LT.II)
     & FACX=4.*DX(I,J,K)/(2.*DX(I,J,K)+DX(I+1,J,K)+DX(I-1,J,K))
      IF(J.GT.1.AND.J.LT.JJ)
     & FACY=4.*DY(I,J,K)/(2.*DY(I,J,K)+DY(I,J+1,K)+DY(I,J-1,K))
      IF(K.GT.1.AND.K.LT.KK)
     & FACZ=4.*DZ(I,J,K)/(2.*DZ(I,J,K)+DZ(I,J,K+1)+DZ(I,J,K-1))
C
      A1(I,J,K)=FACX*KX(I,J,K)*DY(I,J,K)*DZ(I,J,K)
      A2(I,J,K)=FACY*KY(I,J,K)*DX(I,J,K)*DZ(I,J,K)
      A3(I,J,K)=FACZ*KZ(I,J,K)*DX(I,J,K)*DY(I,J,K)
   30 CONTINUE
      IF(II.EQ.1)GO TO 501
      DO 50 K=1,KK
      DO 50 J=1,JJ
      DO 50 I=2,II
      TX(I,J,K)=.012656*A1(I-1,J,K)*A1(I,J,K)/(DX(I-1,J,K)*A1(I,J,K)
     &+DX(I,J,K)*A1(I-1,J,K))
   50 CONTINUE
  501 CONTINUE
      IF(JJ.EQ.1)GO TO 551
      DO 55 K=1,KK
      DO 55 J=2,JJ
      DO 55 I=1,II
      TY(I,J,K)=.012656*A2(I,J-1,K)*A2(I,J,K)/(DY(I,J-1,K)*A2(I,J,K)
     &+DY(I,J,K)*A2(I,J-1,K))
   55 CONTINUE
  551 CONTINUE
      IF(KK.EQ.1)GO TO 600
      DO 60 K=2,KK
      DO 60 J=1,JJ
      DO 60 I=1,II
      TZ(I,J,K)=.012656*A3(I,J,K-1)*A3(I,J,K)/(DZ(I,J,K-1)*A3(I,J,K)
     &+DZ(I,J,K)*A3(I,J,K-1))
   60 CONTINUE
600   CONTINUE
C
C****** TRANSMISSIBILITY MODIFICATIONS
C
      IF(NTRSW.EQ.1) CALL SHALE(NX,NY,NZ)
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      READ(NI,*) NUMTX,NUMTY,NUMTZ,ITCODE
      N1READ=N1READ+1
      IF(NUMTX.EQ.0) GO TO 100
      WRITE(NO,31)
      DO 275 L=1,NUMTX
      READ(NI,*)I,J,K,TX(I,J,K)
      N1READ=N1READ+1
  275 WRITE(NO,32)I,J,K,TX(I,J,K)
  100 IF(ITCODE.NE.1)GO TO 8530
      WRITE(NO,44)
      DO 853 K=1,KK
      WRITE(NO,38)K
      DO 854 J=1,JJ
  854 WRITE(NO,72)(TX(I,J,K),I=1,II)
  853 CONTINUE
8530  CONTINUE
      IF(NUMTY.EQ.0) GO TO 110
      WRITE(NO,34)
      DO 276 L=1,NUMTY
      READ(NI,*)I,J,K,TY(I,J,K)
      N1READ=N1READ+1
  276 WRITE(NO,32)I,J,K,TY(I,J,K)
  110 IF(ITCODE.NE.1) GO TO 8550
      WRITE(NO,47)
      DO 855 K=1,KK
      WRITE(NO,38)K
      DO 856 J=1,JJ
  856 WRITE(NO,72)(TY(I,J,K),I=1,II)
  855 CONTINUE
8550  CONTINUE
      IF(NUMTZ.EQ.0) GO TO 120
      WRITE(NO,37)
      DO 277 L=1,NUMTZ
      READ(NI,*)I,J,K,TZ(I,J,K)
      N1READ=N1READ+1
  277 WRITE(NO,32)I,J,K,TZ(I,J,K)
  120 IF(ITCODE.NE.1) GO TO 8570
      WRITE(NO,48)
      DO 857 K=1,KK
      WRITE(NO,38)K
      DO 858 J=1,JJ
  858 WRITE(NO,72)(TZ(I,J,K),I=1,II)
  857 CONTINUE
8570  CONTINUE
C
      IF(KTR.EQ.0)RETURN
C
      WRITE(NO,43)
      DO 42 K=1,KK
      DO 42 J=1,JJ
      DO 42 I=1,II
      WRITE(NO,41)I,J,K,A1(I,J,K),A2(I,J,K),A3(I,J,K),
     &TX(I,J,K),TY(I,J,K),TZ(I,J,K)
42    CONTINUE
      RETURN
C
   69 FORMAT(40A2)
   31 FORMAT(//T15,'**********TRANSMISSIBILITY (TX) NODE MODIFICATIONS',
     & '**********',//T15,'     I    J    K         NEW TX VALUE')
   32 FORMAT(15X,3I5,5X,F14.4)
   34 FORMAT(//T15,'**********TRANSMISSIBILITY (TY) NODE MODIFICATIONS',
     & '**********',//T15,'     I    J    K         NEW TY VALUE')
   37 FORMAT(//T15,'**********TRANSMISSIBILITY (TZ) NODE MODIFICATIONS',
     & '**********',//T15,'     I    J    K         NEW TZ VALUE')
   38 FORMAT(/1X,'K =',I2/)
   41 FORMAT(1X,3I5,6F15.3)
   43 FORMAT(//1X,'    I    J    K        A1            A2      ',
     & '      A3            TX            TY            TZ'/)
   44 FORMAT(//T15,'**********TRANSMISSIBILITY (TX) DISTRIBUTION'
     & ,'**********'/)
   47 FORMAT(//T15,'**********TRANSMISSIBILITY (TY) DISTRIBUTION'
     & ,'**********'/)
   48 FORMAT(//T15,'**********TRANSMISSIBILITY (TZ) DISTRIBUTION'
     & ,'**********'/)
   72 FORMAT(1X,20F6.1)
C
      END
      SUBROUTINE TABLE
C
C==================================================================TABLE
C
C   READ PVT AND RELATIVE PERMEABLITY DATA.
C
C=======================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      CHARACTER*2 TITLE(80)
      REAL KROT,KRWT,KRGT,MUOT,MUWT,MUGT
      COMMON /IO/ NI,NO
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /COUNT/ N1READ,N2READ
C
      WRITE(NO,111)
C
C****  RELATIVE PERMEABILITY & CAPILLARY PRESSURE TABLE
C
      READ(NI,1) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,11) (TITLE(IH),IH=1,40)
      DO 5 I=1, NTE
      READ(NI,*) SAT(I),KROT(I),KRWT(I),KRGT(I),PCOWT(I),PCGOT(I)
      N1READ=N1READ+1
      WRITE(NO,21)SAT(I),KROT(I),KRWT(I),KRGT(I),PCOWT(I),PCGOT(I)
      IF(SAT(I).GE.1.0999) GO TO 10
5     CONTINUE
C
C****  BUBBLE POINT & MAXIMUM PRESSURES
C
10    MSAT=I
      READ(NI,1) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,11)(TITLE(IH),IH=1,40)
      READ(NI,*) PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS
      N1READ=N1READ+1
      WRITE(NO,31)PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS
C
C****  OIL DATA TABLE
C
      READ(NI,1) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,11) (TITLE(IH),IH=1,40)
      DO 15 I=1, NTE
      READ(NI,*) POT(I),MUOT(I),BOT(I),RSOT(I)
      N1READ=N1READ+1
      WRITE(NO,23)POT(I),MUOT(I),BOT(I),RSOT(I)
      RSOT(I)=0.17809*RSOT(I)
      IF(POT(I).GE.PMAXT)GO TO 20
15    CONTINUE
C
C****  WATER DATA TABLE
C
20    MPOT=I
      READ(NI,1) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,11) (TITLE(IH),IH=1,40)
      DO 25 I=1,NTE
      READ(NI,*) PWT(I),MUWT(I),BWT(I),RSWT(I)
      N1READ=N1READ+1
      WRITE(NO,23) PWT(I),MUWT(I),BWT(I),RSWT(I)
      RSWT(I)=0.17809*RSWT(I)
      IF(PWT(I).GE.PMAXT)GO TO 30
25    CONTINUE
C
C****  GAS & ROCK DATA TABLE
C
30    MPWT=I
      READ(NI,1) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,11) (TITLE(IH),IH=1,40)
      DO 35 I=1,999
      READ(NI,*) PGT(I),MUGT(I),BGT(I),CRT(I)
      N1READ=N1READ+1
      WRITE(NO,24) PGT(I),MUGT(I),BGT(I),CRT(I)
      IF(PGT(I).GE.PMAXT) GO TO 40
35    CONTINUE
C
40    MPGT=I
      CONTINUE
C
C****  OIL, WATER, & GAS DENSITIES AT STOCK TANK CONDITIONS
C
      READ(NI,1) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
      WRITE(NO,11) (TITLE(IH),IH=1,40)
      READ(NI,*)RHOSCO,RHOSCW,RHOSCG
      N1READ=N1READ+1
      WRITE(NO,4)RHOSCO,RHOSCW,RHOSCG
C
C**** CALCULATE SLOPES DBO/DP, DRSO/DP, DBW/DP, DRSW/DP, DBG/DP
C
      WRITE(NO,222)
      WRITE(NO,223)
      DO 100 I=2,MPOT
      DIV=1.0/(POT(I)-POT(I-1))
      BOPT(I)=(BOT(I)-BOT(I-1))*DIV
      RSOPT(I)=(RSOT(I)-RSOT(I-1))*DIV
100   WRITE(NO,224)POT(I),BOT(I),BOPT(I),RSOT (I),RSOPT(I)
      WRITE(NO,333)
      DO 200 I=2,MPWT
      DIV=1.0/(PWT(I)-PWT(I-1))
      BWPT(I)=(BWT(I)-BWT(I-1))*DIV
      RSWPT(I)=(RSWT(I)-RSWT(I-1))*DIV
200   WRITE(NO,224)PWT(I),BWT(I),BWPT(I),RSWT(I),RSWPT(I)
C
      WRITE(NO,444)
      DO 300 I=2,MPGT
      BGPT(I)=(BGT(I)-BGT(I-1))/(PGT(I)-PGT(I-1))
300   WRITE(NO,445)PGT(I),BGT(I),BGPT(I)
C
C**** NOTE: SLOPE AT POINT "I" IS BASED ON POINTS "I" AND "I-1"
C
1     FORMAT(40A2)
4     FORMAT(2X,8F10.4)
11    FORMAT(//,6X,40A2/)
21    FORMAT(1X,4F10.4,2F10.2)
23    FORMAT(3X,F10.1,1X,F8.4,2X,F8.4,2X,F8.2)
24    FORMAT(3X,F10.1,1X,F8.4,1X,E10.4,1X,E10.3)
31    FORMAT(F10.2,E10.3,1X,E10.3,2F10.2,I5)
111   FORMAT(//T15,'***** EMPIRICAL DATA TABLE *****'//)
222   FORMAT(//T15,'***** SLOPES FOR COMPRESSIBILITY CALCULATIONS',
     &' ****')
223   FORMAT(//T16,'P',T25,'BO',T33,'DBO/DP',T47,'RSO',
     &T55,'DRSO/DP'/)
224   FORMAT( T13,F7.1,F8.4,3X,E11.4,F8.1,2X,E11.4)
333   FORMAT(//T16,'P',T25,'BW',T33,'DBW/DP',T47,'RSW',
     &T55,'DRSW/DP'/)
444   FORMAT(//15X,'P',12X,'BG',12X,'DBGDP')
445   FORMAT(10X,F9.1,1X,2E15.4)
C
      RETURN
      END
      SUBROUTINE UINIT1(NX,NY,NZ,PBOT,PN,SON,SWN,SGN,EL,P,SO,SW,SG,
     &    DPP,DPSO,DPSW)
C
C=================================================================UINIT1
C
C   READ INITIAL PRESSURE AND SATURATION ARRAY.
C
C=======================================================================
C
       INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
C	NEW DP LOCAL VARIABLES **RJH**
C
      REAL MBEO,MBEW,MBEG,KROT,KRWT,KRGT,MUOT,MUWT,MUGT
      CHARACTER*2 TITLE(80)
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      COMMON /MBE/ MBEO,MBEW,MBEG
      COMMON /PSCNTL/ KPI,KSI
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
     
	DIMENSION PBOT(NX,NY,NZ),PN(NX,NY,NZ),
     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ)
      DIMENSION EL(NX,NY,NZ),P(NX,NY,NZ)
      DIMENSION SO(NX,NY,NZ),SW(NX,NY,NZ),SG(NX,NY,NZ)
      REAL*8 DPP(NX,NY,NZ),DPSO(NX,NY,NZ),DPSW(NX,NY,NZ)

      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),
     & RSOPT(NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),
     & RSWPT(NTE),PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /ADD1/ IM,JM,KM,ETI,FT,FTMAX,DELT
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      COMMON /COUNT/ N1READ,N2READ
      DIMENSION PPK(NZ1),SOK(NZ1),SWK(NZ1)
C
      MBEO=0.0
      MBEW=0.0
      MBEG=0.0
C
C***   INITIALIZE BUBBLE POINT PRESSURE ARRAY
C
      DO 8 I=1,II
      DO 8 J=1,JJ
      DO 8 K=1,KK
 8    PBOT(I,J,K) = PBO
C
      READ(NI,69) (TITLE(IH),IH=1,40)
      N2READ=N2READ+1
      READ(NI,*) KPI,KSI
      N2READ=N2READ+1
      IF(KPI.NE.0) GO TO 5601
      READ(NI,*) PI,PGOC,WOC,GOC
      N2READ=N2READ+1
      DO 200 K=1,KK
      DO 200 J=1,JJ
      DO 200 I=1,II
      IF(EL(I,J,K).LT.GOC) GO TO 175
      IF(EL(I,J,K).GT.WOC)GO TO 150
      BPT=PBOT(I,J,K)
C      WRITE(NO,1)(POT(L),MUOT(L),BOT(L),RSOT(L),L=1,NTE)
C 1     FORMAT(5X,'UINIT1 = ',4F15.5)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PI,BBO)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PI,RSO)
      RHOO=(RHOSCO+RSO*RHOSCG)/BBO
      PN(I,J,K)=PI+RHOO*(EL(I,J,K)-WOC)/144.
      GO TO 200
  150 CALL INTERP(PWT,BWT,MPWT,PI,BBW)
      CALL INTERP(PWT,RSWT,MPWT,PI,RSW)
      RHOW=(RHOSCW+RSW*RHOSCG)/BBW
      PN(I,J,K)=PI+RHOW*(EL(I,J,K)-WOC)/144.
      GO TO 200
  175 CALL INTERP(PGT,BGT,MPGT,PGOC,BBG)
      RHOG=RHOSCG/BBG
      PN(I,J,K)=PGOC+RHOG*(EL(I,J,K)-GOC)/144.
  200 CONTINUE
      GO TO 3015
 5601 IF(KPI.NE.1) GO TO 5610
      DO 3010 K=1,KK
      DO 3005 J=1,JJ
      READ(NI,*)(DPP(I,J,K),I=1,II)
      DO 3034 I=1,II
      PN(I,J,K)=DPP(I,J,K)
3034  CONTINUE
      N2READ=N2READ+1+INT(II/10)
 3005 IF(MOD(KK,10).EQ.0)N2READ=N2READ-1
3010  CONTINUE
      GO TO 3015
C***  PRESSURE AT EACH LAYER TO BE READ ....
5610  READ(NI,*)(PPK(K),K=1,KK)
      N2READ=N2READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N2READ=N2READ-1
      DO 3080 K=1,KK
      DO 3080 J=1,JJ
      DO 3080 I=1,II
3080  PN(I,J,K)=PPK(K)      
C
3015  CONTINUE
C**** INITIALIZE N+1 PRESSURE ARRAY.
      DO 3012 I=1,II
      DO 3012 J=1,JJ
      DO 3012 K=1,KK
 3012 P(I,J,K)=PN(I,J,K)
      IF(KSI.NE.0)GO TO 5600
      READ(NI,*) SOI,SWI,SGI
      N2READ=N2READ+1
      DO 30 K=1,KK
      DO 30 J=1,JJ
      DO 30 I=1,II
      SON(I,J,K)=SOI
      SWN(I,J,K)=SWI
      SGI=1.0-SOI-SWI
      SGN(I,J,K)=SGI
      SO(I,J,K)=SOI
      SW(I,J,K)=SWI
      SG(I,J,K)=SGI
      IF(SG(I,J,K).LT.0.0) SG(I,J,K)=0.0
   30 CONTINUE
      RETURN
5600  IF(KSI.NE.1) GO TO 5620
      DO 3011 K=1,KK
      DO 3006 J=1,JJ
      READ(NI,*)(DPSO(I,J,K),I=1,II)
      N2READ=N2READ+1+INT(II/10)
3006  IF(MOD(KK,10).EQ.0)N2READ=N2READ-1
3011  CONTINUE
      DO 3020 K=1,KK
      DO 3007 J=1,JJ
      READ(NI,*)(DPSW(I,J,K),I=1,II)
      N2READ=N2READ+1+INT(II/10)
3007  IF(MOD(KK,10).EQ.0)N2READ=N2READ-1
3020  CONTINUE
C	MODIFICATION FOR DOUBLE PRECISION READ OF SO & SW IN RESTART  **RJH**
      DO 3031 K=1,KK
      DO 3032 J=1,JJ
      DO 3033 I=1,II
      SO(I,J,K)=DPSO(I,J,K)
      SW(I,J,K)=DPSW(I,J,K)
3033  CONTINUE
3032  CONTINUE
3031  CONTINUE
C	END OF MODIFICATION **RJH**
      DO 3030 K=1,KK
      DO 3030 J=1,JJ
      DO 3030 I=1,II
      SG(I,J,K)=1.0-SO(I,J,K)-SW(I,J,K)
      IF(SG(I,J,K).LT.0.0) SG(I,J,K)=0.0
      SON(I,J,K)=SO(I,J,K)
      SWN(I,J,K)=SW(I,J,K)
      SGN(I,J,K)=SG(I,J,K)
 3030 CONTINUE
69    FORMAT(40A2)
C        IF(COP .EQ. 0.D0) RETURN
      IF(KPI .NE. 1) RETURN
      READ(NI,*) ETI,DELT
      READ(NI,*)  COP,CWP,CGP,CWI,CGI
	RETURN
C       END OF MODIFICATION **RJH**
C***  OIL/WATER SATURATIONS AT EACH LAYER TO BE READ ....
5620  READ(NI,*)(SOK(K),K=1,KK)
      N2READ=N2READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N2READ=N2READ-1
      READ(NI,*)(SWK(K),K=1,KK)
      N2READ=N2READ+1+INT(KK/10)
      IF(MOD(KK,10).EQ.0)N2READ=N2READ-1
      DO 3090 K=1,KK
      DO 3090 J=1,JJ
      DO 3090 I=1,II
       SO(I,J,K)=SOK(K)
       SW(I,J,K)=SWK(K)
       SG(I,J,K)=1.-SOK(K)-SWK(K)
       SON(I,J,K)=SOK(K)
       SWN(I,J,K)=SWK(K)
       SGN(I,J,K)=1.-SOK(K)-SWK(K)
3090  CONTINUE
C      READ(NI,*) ETI,DELT
C      READ(NI,*)  COP,CWP,CGP,CWI,CGI
      RETURN
      END
C
      SUBROUTINE CODES
C
C==================================================================CODES
C
C   READ SOLUTION METHOD, SOLUTION PARAMETERS, AND DIAGNOSTIC CONTROLS
C
C=======================================================================
C
      CHARACTER*2 TITLE(80)
      COMMON /IO/ NI,NO
      COMMON /CODE/ KSM1,KSN1,KCO1,NN,
     &FACT1,FACT2,TMAX,KSOL,MITER,OMEGA,TOL,TOL1,
     & KSN,KSM,KCO,KTR,KCOFF,DSMAX,DPMAX,
     &WORMAX,GORMAX,PAMIN,PAMAX
      READ(NI,69) (TITLE(IH),IH=1,40)
      READ(NI,*) KSN1,KSM1,KCO1,KTR,KCOFF
C EVERY KSM1 TH STEP SOLUTION MATRIX WILL BE WRITTEN
C EVERY KSN1 TH STEP LSOR DATA WILL BE WRITTEN
C EVERY KCO1 TH STEP COMPRESSIBILITIES & VOLUME FACTORS WILL BE WRITTEN
C IF KTR IS NONZERO TRANSMISSIBILITIES WILL BE WRITTEN
      READ(NI,69) (TITLE(IH),I=1,40)
      READ(NI,*)NN,FACT1,FACT2,TMAX,WORMAX,GORMAX,PAMIN,PAMAX
      WRITE(NO,59) NN,FACT1,FACT2,TMAX,WORMAX,GORMAX,PAMIN,PAMAX
C NN--MAXIMUM NUMBER OF TIME-STEPS
C FACT1--FACTOR FOR INCREASING TIME-STEP SIZE
C FACT2--FACTOR FOR DECREASING TIME-STEP SIZE
C TMAX--MAXIMUM SIMULATION TIME
      READ(NI,69) (TITLE(IH),I=1,40)
C      READ(NI,*)KSOL,MITER,OMEGA,TOL,TOL1,DSMAX,DPMAX
      READ(NI,*) MITER,OMEGA,TOL,TOL1,DSMAX,DPMAX
C      IF(KSOL.EQ.1) WRITE(NO,73)
C      IF(KSOL.EQ.2) WRITE(NO,75) MITER,OMEGA,TOL,TOL1
C      IF(KSOL.EQ.3) WRITE(NO,77)
        WRITE(NO,75) MITER,OMEGA,TOL,TOL1
      WRITE(NO,79) DSMAX,DPMAX
      KSN=KSN1
      KSM=KSM1
      KCO=KCO1
      RETURN
C
   73 FORMAT(//T15,'SOLUTION METHOD IS BAND')
   75 FORMAT(//T15,'SOLUTION METHOD IS LSOR:',
     &/T20,'MAXIMUM NUMBER OF ITERATIONS      (MITR) = ',5X,I5,
     &/T20,'INITIAL ACCELERATION PARAMETER    (OMEGA) = ',F10.4,
     &/T20,'MAXIMUM PRESSURE RESIDUAL          (TOL) = ',F10.4,
     &/T20,'PARAMETER FOR CHANGING OMEGA      (TOL1) = ',F10.4)
   77 FORMAT(//T15,'SOLUTION METHOD IS D4.')
   79 FORMAT(/T15,'AUTOMATIC TIME STEP CRITERIA:',
     &/T20,'MAXIMUM ALLOWED SATURATION CHANGE (DSMAX) = ',F10.4,
     &/T20,'MAXIMUM ALLOWED PRESSURE CHANGE   (DPMAX) = ',F10.4/)
59    FORMAT(//T15,'MAXIMUM NUMBER OF TIME-STEPS =',I5/
     &T15,'FACTOR FOR INCREASING DELT =',F10.6,3X,
     &'WHEN DSMAX AND DPMAX NOT EXCEEDED.',/,
     &T15,'FACTOR FOR DECREASING DELT =',F10.6,3X,
     &'WHEN DSMAX OR DPMAX IS EXCEEDED.',/,
     &T15,'MAXIMUM SIMULATION TIME =',F11.6/,
     &T15,'MAXIMUM RESERVOIR WOR/TIME-STEP =',F8.0,' STB/STB'/
     &T15,'MAXIMUM RESERVOIR GOR/TIME-STEP =',F8.0,' SCF/STB'/
     &T15,'MINIMUM AVERAGE RESERVOIR PRESSURE/TIME-STEP =',F8.0/
     &T15,'MAXIMUM AVERAGE RESERVOIR PRESSURE/TIME-STEP =',F8.0//)
69    FORMAT(40A2)
311   FORMAT (I5,F10.0,I5)
      END
C
C
C
C      SUBROUTINE NODES(NVQN,NVQNH,NVQNS,NX,NY,NZ)
      SUBROUTINE NODES(NVQN,NVQNH,NVQNS,NX,NY,NZ,PN,EL,
     &    KX,KY,KZ,DX,DY,DZ)
C
C==================================================================NODES
C
C   READ TIME STEP AND VERTICAL WELL INFORMATION.
C
C         NVQN = # OF VERTICAL WELLS
C         NVQNH= # OF HORIZONTAL WELLS
C         NVQNS= # OF SLANT WELLS
C=======================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      CHARACTER*2 TITLE(80)
      CHARACTER*5 WELLID
      COMMON /IDEN/ WELLID(NW)
      COMMON /IO/ NI,NO
      REAL KX(NX,NY,NZ),KY(NX,NY,NZ),KZ(NX,NY,NZ)
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
      COMMON /R1XYZ/ R11X,R21X,R31X,R41X,R11Y,R21Y,R31Y,R41Y,
     & R11Z,R21Z,R31Z,R41Z
      COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
      COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
      COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K
      DIMENSION EL(NX,NY,NZ)
      DIMENSION DX(NX,NY,NZ),DY(NX,NY,NZ),DZ(NX,NY,NZ)
C
C**** ESTABLISH RATE-SPECIFIED & PRESSURE-SPECIFIED WELLS
C
      READ(NI,69) (TITLE(IH),IH=1,40)
      READ(NI,*)NVQN,NVQNH,NVQNS
      IF(NVQN.EQ.0)GO TO 2011
      WRITE(NO,67)
      WRITE(NO,68)
      DO 2000 J=1,NVQN
C**   INPUT RATE SIGNS BY CONVENTION:
C**   INJECTORS ARE NEGATIVE; PRODUCERS ARE POSITIVE.
      READ(NI,1) WELLID(J)
      READ(NI,*) IQN1(J),IQN2(J),IQN3(J),LAYER(J),KIP(J),
     & QVO(J),QVW(J),QVG(J),QVT(J)
      IQ3=IQN3(J)
      LAY=IQ3+(LAYER(J)-1)
      DO 1900 K=IQ3,LAY
      READ(NI,*) PID(J,K),PWF(J,K)
      WRITE(NO,70)IQN1(J),IQN2(J),K,QVO(J),QVW(J),QVG(J),
     & QVT(J),PWF(J,K),PID(J,K)
 1900 CONTINUE
2000  CONTINUE
C
      WRITE(NO,33)
      DO 2999 J=1,NVQN
      IQ3=IQN3(J)
      LAY=IQ3+(LAYER(J)-1)
      DO 2999 K=IQ3,LAY
      IF(KIP(J).EQ.1) WRITE(NO,2995) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.2) WRITE(NO,2996) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.3) WRITE(NO,2997) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.-1) WRITE(NO,3005) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.-2) WRITE(NO,3006) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.-3) WRITE(NO,3007) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.-11) WRITE(NO,3015) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.-12) WRITE(NO,3016) IQN1(J),IQN2(J),K
      IF(KIP(J).EQ.-13) WRITE(NO,3017) IQN1(J),IQN2(J),K
2999  CONTINUE
C
2011  IF(NVQNH.NE.0) CALL NODESH(NVQN,NVQNH,NVQNS)
      IF(NVQNS.NE.0) CALL SLANT(NVQN,NVQNH,NVQNS,NX,NY,NZ,
     &  EL,PN,KX,KY,KZ,DX,DY,DZ)
1     FORMAT(A5)
69    FORMAT(40A2)
721   FORMAT(4F10.0)
67    FORMAT(//T15,'RESERVOIR CONTAINS THE FOLLOWING RATE NODES
     & (FOR VERTICAL WELLS)'/)
68    FORMAT(T2,'   NODE   ',3X,'  OIL(STBD)',
     &3X,'WATER(STBD)',3X,'GAS(MCFD)',3X,'TOTAL(RBD)',3X,'BHFP(PSIA)'
     &,3X,'PID')
70    FORMAT(1X,3I3,3X,F11.2,4F13.2,F10.6)
71    FORMAT(9X,3I4,F10.2,3F5.1,10X,F10.2,3F5.1,F12.3)
33    FORMAT(/)
 2995 FORMAT(T15,'BLOCK ',3I3,
     & ' CONTAINS A RATE SPECIFIED PRODUCING WELL')
 2996 FORMAT(T15,'BLOCK ',3I3,
     & ' CONTAINS A RATE SPECIFIED WATER INJECTION WELL')
 2997 FORMAT(T15,'BLOCK ',3I3,
     & ' CONTAINS A RATE SPECIFIED GAS INJECTION WELL')
 3005 FORMAT(T15,'BLOCK ',3I3,' CONTAINS AN ',
     & 'EXPLICIT PRESSURE SPECIFIED PRODUCING WELL')
 3006 FORMAT(T15,'BLOCK ',3I3, ' CONTAINS AN ',
     & 'EXPLICIT PRESSURE SPECIFIED WATER INJECTION WELL')
 3007 FORMAT(T15,'BLOCK ',3I3,' CONTAINS AN ',
     & 'EXPLICIT PRESSURE SPECIFIED GAS INJECTION WELL')
 3015 FORMAT(T15,'BLOCK ',3I3,' CONTAINS AN ',
     & 'IMPLICIT PRESSURE SPECIFIED PRODUCING WELL')
 3016 FORMAT(T15,'BLOCK ',3I3,' CONTAINS AN ',
     & 'IMPLICIT PRESSURE SPECIFIED WATER INJECTION WELL')
 3017 FORMAT(T15,'BLOCK ',3I3,' CONTAINS AN ',
     & 'IMPLICIT PRESSURE SPECIFIED GAS INJECTION WELL')
      RETURN
      END
C
      SUBROUTINE NODESH(NVQN,NVQNH,NVQNS)
C
C=================================================================NODESH
C
C   READ HORIZONTAL WELL INFORMATION.
C
C      (NVQN,NVQNH,NVQNS ARE # OF VERTICAL, HORIZONTAL, AND SLANT WELLS,
C         RESPECTIVELY)
C=======================================================================
C
	INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      CHARACTER*2 TITLE(80)
      CHARACTER*5 WELLID
      COMMON /IDEN/ WELLID(NW)
      COMMON /IO/ NI,NO
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
C
C**** ESTABLISH RATE-SPECIFIED & PRESSURE-SPECIFIED WELLS
C
      WRITE(NO,607)
      WRITE(NO,608)
      N1=NVQN+1
      N2=NVQN+NVQNH
      DO 2010 J=N1,N2
C**   INPUT RATE SIGNS BY CONVENTION:
C**   INJECTORS ARE NEGATIVE; PRODUCERS ARE POSITIVE.
      READ(NI,1212) WELLID(J)
      READ(NI,*) LAYER(J),KIP(J),QVO(J),QVW(J),QVG(J),QVT(J),COND(J)
      DO 2001 K=1,LAYER(J)
      READ(NI,*) IQH1(J,K),IQH2(J,K),IQH3(J,K),PID(J,K),PWF(J,K)
      WRITE(NO,703)IQH1(J,K),IQH2(J,K),IQH3(J,K),QVO(J),QVW(J),
     & QVG(J),QVT(J),PWF(J,K),PID(J,K)
2001  CONTINUE
2010  CONTINUE
C
      WRITE(NO,313)
      DO 2998 J=N1,N2
      IF(KIP(J).EQ.1) WRITE(NO,2905) WELLID(J)
      IF(KIP(J).EQ.2) WRITE(NO,2906) WELLID(J)
      IF(KIP(J).EQ.3) WRITE(NO,2907) WELLID(J)
      IF(KIP(J).EQ.-1) WRITE(NO,2915) WELLID(J)
      IF(KIP(J).EQ.-2) WRITE(NO,2916) WELLID(J)
      IF(KIP(J).EQ.-3) WRITE(NO,2917) WELLID(J)
      IF(KIP(J).EQ.-11) WRITE(NO,2925) WELLID(J)
      IF(KIP(J).EQ.-12) WRITE(NO,2926) WELLID(J)
      IF(KIP(J).EQ.-13) WRITE(NO,2927) WELLID(J)
      DO 2997 K=1,LAYER(J)
      WRITE(NO,2976) IQH1(J,K),IQH2(J,K),IQH3(J,K)
2997  CONTINUE
2998  CONTINUE
C
1212  FORMAT(A5)
607   FORMAT(//T15,'RESERVOIR CONTAINS FOLLOWING RATE NODES
     & (FOR HORIZONTAL OR SLANT WELLS) '/)
608    FORMAT(T2,'   NODE   ',3X,'  OIL(STBD)',
     &3X,'WATER(STBD)',3X,'GAS(MCFD)',3X,'TOTAL(RBD)',3X,'BHFP(PSIA)'
     &,3X,'PID')
703   FORMAT(1X,3I3,3X,F11.2,4F13.2,F10.6)
313   FORMAT(/)
 2976 FORMAT(T15,'BLOCK ',3I3)
 2905 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,
     & ' IS A RATE SPECIFIED PRODUCING WELL WHICH CONTAINS:')
 2906 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,
     & ' IS A RATE SPECIFIED WATER INJECTION WELL WHICH CONTAINS:')
 2907 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,
     & ' IS A RATE SPECIFIED GAS INJECTION WELL WHICH CONTAINS:')
 2915 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,' IS AN ',
     & 'EXPLICIT PRESSURE SPECIFIED PRODUCING WELL WHICH CONTAINS:')
 2916 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,' IS AN EXPLICIT ',
     & 'PRESSURE SPECIFIED WATER INJECTION WELL WHICH CONTAINS:')
 2917 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,' IS AN ',
     & 'EXPLICIT PRESSURE SPECIFIED GAS INJECTION WELL WHICH CONTAINS:')
 2925 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,' IS AN ',
     & 'IMPLICIT PRESSURE SPECIFIED PRODUCING WELL WHICH CONTAINS:')
 2926 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,' IS AN IMPLICIT ',
     & 'PRESSURE SPECIFIED WATER INJECTION WELL WHICH CONTAINS:')
 2927 FORMAT(T10,'HORIZONTAL OR SLANT WELL ',A5,' IS AN ',
     & 'IMPLICIT PRESSURE SPECIFIED GAS INJECTION WELL WHICH CONTAINS:')
      RETURN
      END

C
C      SUBROUTINE QRATE(IWELTYPE,NVQN,NVQNH,NX,NY,NZ)
      SUBROUTINE QRATE(IWELTYPE,NVQN,NVQNH,NX,NY,NZ,PBOT,
     &                         PN,P,SO,SW,SG,QO,QW,QG)
C
C==================================================================QRATE
C
C   COMPUTE WELL RATES WHEN RATE OR EXPLICIT PRESSURE CONSTRAIN
C   ACTIVATED;
C        PRATEI: REVISES PRESSURE EQUATION COEFFICENTS WHEN IMPLICIT
C                PRESSURE CONSTRAIN ACTIVATED;
C        PRATEO: COMPUTE RATES WHEN IMPLICIT PRESSURE CONSTRAIN
C                ACTIVATE.
C
C        (IWELTYPE = 1 FOR VERTICAL WELL OPTION        )
C        (         = 2 FOR HORIZONTAL/SLANT WELL OPTION) 
C=====================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
C	NEW DP LOCAL VARIABLES **RJH**
C       MODIFICATION FOR DOUBLE PRECISION READ OF SO & SW IN RESTART  **RJH**

      REAL KROT,KRWT,KRGT,MUWT,MUOT,MUGT
     & ,MUO,MUW,MUG,KRO,KRW,KRG
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
      DIMENSION PBOT(NX,NY,NZ)
      DIMENSION PN(NX,NY,NZ)
C     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
C     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      DIMENSION P(NX,NY,NZ)
      DIMENSION SO(NX,NY,NZ),SW(NX,NY,NZ),SG(NX,NY,NZ)
C      REAL*8 DPP(NX,NY,NZ),DPSO(NX,NY,NZ),DPSW(NX,NY,NZ)
      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      DIMENSION QO(NX,NY,NZ),QW(NX,NY,NZ),QG(NX,NY,NZ)
      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
C
      IF(IWELTYPE.EQ.1)THEN
       N1=1
       N2=NVQN
      ELSE 
       N1=NVQN+1
       N2=NVQN+NVQNH 
      ENDIF
C
      DO 105 J=N1,N2
      IF(IWELTYPE.EQ.1)THEN
       IQ1=IQN1(J)
       IQ2=IQN2(J)
       IQ3=IQN3(J)
       LAY=IQ3+(LAYER(J)-1)
      ELSE
       IQ3=1
       LAY=LAYER(J)
      ENDIF
      DO 1170 KKK=IQ3,LAY
      K=KKK
      IF(IWELTYPE.EQ.2)THEN
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      PWFC(J,KKK)=-1.0
      PP=P(IQ1,IQ2,K)
      BPT=PBOT(IQ1,IQ2,K)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,PP,MUO)
      CALL INTERP(PWT,MUWT,MPWT,PP,MUW)
      CALL INTERP(PGT,MUGT,MPGT,PP,MUG)
      SSO=SO(IQ1,IQ2,K)
      SSW=SW(IQ1,IQ2,K)
      SSG=SG(IQ1,IQ2,K)
      CALL INTERP(SAT,KRWT,MSAT,SSW,KRW)
      CALL INTERP(SAT,KROT,MSAT,SSO,KRO)
      CALL INTERP(SAT,KRGT,MSAT,SSG,KRG)
      GMW(J,KKK)=KRW/MUW
      GMO(J,KKK)=KRO/MUO
      GMG(J,KKK)=KRG/MUG
C      IF(TMP = 1/(1+GMO(J,K)/GMW(J,K))
C      WRITE(*,*) J,K,KRW,KRO,GMW(J,K),GMO(J,K)
 1170 CONTINUE
      IF(KIP(J).LT.0) GO TO 105
      COND1=0.
      COND2=0.
      PCOND1=0.
      PCOND2=0.
      IF(KIP(J).NE.1) GO TO 1190
      ITERQ=0
      QDENOM=0.0
      ALPHAO=0.0
      ALPHAW=0.0
      ALPHAG=0.0
      IF (IWELTYPE.EQ.2)THEN
       IQ3=1
       LAY=LAYER(J)
      ELSE
       LAY=IQ3+(LAYER(J)-1)
      ENDIF
 1172 ITERQ=ITERQ+1
C
      DO 1180 KKK=IQ3,LAY
      K=KKK
      IF(IWELTYPE.EQ.2)THEN
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      PP=P(IQ1,IQ2,K)
      BPT=PBOT(IQ1,IQ2,K)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PP,BBO)
      CALL INTERP(PWT,BWT,MPWT,PP,BBW)
      CALL INTERP(PGT,BGT,MPGT,PP,BBG)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PP,RSO)
      CALL INTERP(PWT,RSWT,MPWT,PP,RSW)
      IF(IWELTYPE.EQ.2.AND.COND(J).EQ.1) GO TO 1179
      IF(ITERQ.NE.1) GO TO 1174
      QDENOM=QDENOM+PID(J,KKK)*GMO(J,KKK)/BBO
      GMT=GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK)
      ALPHAO=GMO(J,KKK)/GMT+ALPHAO
      ALPHAW=GMW(J,KKK)/GMT+ALPHAW
      ALPHAG=GMG(J,KKK)/GMT+ALPHAG
      GO TO 1180
 1174 IF(QVT(J).EQ.0.0) GO TO 1176
      TOTOR=QVT(J)*ALPHAO/(ALPHAO+ALPHAW+ALPHAG)
      GO TO 1178
 1176 TOTOR=QVO(J)
 1178 QO(IQ1,IQ2,K)=TOTOR*5.615*PID(J,KKK)*GMO(J,KKK)
     & /(BBO*QDENOM)
      QW(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*GMW(J,KKK)*BBO
     & /(BBW*GMO(J,KKK))
      QG(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*(GMG(J,KKK)*BBO
     & /(BBG*GMO(J,KKK))+RSO)+RSW*QW(IQ1,IQ2,K)
      GO TO 1180
C***** INFINITE CONDUCTIVITY FOR HORIZONTAL/SLANT WELL
C      AT FIXED OIL RATE
 1179 IF(QVT(J).EQ.0.0) GO TO 11792
      IF(ITERQ.NE.1) GO TO 11791
      FACTOR1=PID(J,KKK)*(GMO(J,KKK)/BBO+GMW(J,KKK)/BBW
     & +GMG(J,KKK)/BBG)
      COND1=COND1+FACTOR1
      PCOND1=PCOND1+FACTOR1*PP
      GO TO 1180
11791 PWFC(J,KKK)=(PCOND1-QVT(J))/COND1
      GO TO 11794
11792 IF(ITERQ.NE.1) GO TO 11793
      FACTOR2=PID(J,KKK)*GMO(J,KKK)/BBO
      COND2=COND2+FACTOR2
      PCOND2=PCOND2+FACTOR2*PP
      GO TO 1180
11793 PWFC(J,KKK)=(PCOND2-QVO(J))/COND2
11794 QO(IQ1,IQ2,K)=5.615*PID(J,KKK)*GMO(J,KKK)*(PP-PWFC(J,KKK))/BBO
      QW(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*GMW(J,KKK)*BBO/(BBW*GMO(J,KKK))
      QG(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*(GMG(J,KKK)*BBO/(BBG*GMO(J,KKK))
     & +RSO)+RSW*QW(IQ1,IQ2,K)
 1180 CONTINUE
      IF(ITERQ.EQ.1) GO TO 1172
      GO TO 105
 1190 CONTINUE
      LAY=IQ3+(LAYER(J)-1)
      ITERQ=0
      QDENOM=0.
      IF(IWELTYPE.EQ.2)THEN
       IQ3=1
       LAY=LAYER(J)
      ENDIF
 1192 ITERQ=ITERQ+1
C
      DO 1200 KKK=IQ3,LAY
      K=KKK
      IF (IWELTYPE.NE.2) GO TO 1193
      IQ1=IQH1(J,KKK)
      IQ2=IQH2(J,KKK)
      K=IQH3(J,KKK)
      IF(COND(J).EQ.1) GO TO 1197
 1193 IF(ITERQ.NE.1)GO TO 1194
      QDENOM=QDENOM+PID(J,KKK)*(GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK))
      GO TO 1200
 1194 IF(KIP(J).NE.2) GO TO 1196
      QW(IQ1,IQ2,K)=QVW(J)*5.615*PID(J,KKK)
     & *(GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK))/QDENOM
      GO TO 1200
 1196 QG(IQ1,IQ2,K)=QVG(J)*1000.*PID(J,KKK)
     & *(GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK))/QDENOM
      GO TO 1200
C***** INFINITE CONDUCTIVITY FOR HORIZONTAL/SLANT WELL
C      AT FIXED WATER/GAS RATE
 1197 IF(ITERQ.NE.1) GO TO 11971
      FACTOR1=PID(J,KKK)*(GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK))
      COND1=COND1+FACTOR1
      PCOND1=PCOND1+FACTOR1*PP
      GO TO 1200
11971 IF(KIP(J).NE.2) GO TO 11972
      PWFC(J,KKK)=(PCOND1-QVW(J))/COND1
      QW(IQ1,IQ2,K)=5.615*PID(J,KKK)*GMW(J,KKK)*(PP-PWFC(J,KKK))
      GO TO 1200
11972 PWFC(J,KKK)=(PCOND1-QVG(J)*1000/5.615)/COND1
      QG(IQ1,IQ2,K)=5.615*PID(J,KKK)*GMG(J,KKK)*(PP-PWFC(J,KKK))
 1200 CONTINUE
C
      IF(ITERQ.EQ.1) GO TO 1192
  105 CONTINUE
      IF(IWELTYPE.EQ.1)THEN
       N1=1
       N2=NVQN
      ELSE
       N1=NVQN+1
       N2=NVQN+NVQNH
      ENDIF
      DO 1210 J=N1,N2
      IF(IWELTYPE.EQ.1)THEN
       IQ1=IQN1(J)
       IQ2=IQN2(J)
       IQ3=IQN3(J)
      ENDIF
      IF(KIP(J).LT.0) GO TO 1210
      IF(IWELTYPE.EQ.2)THEN
       IQ3=1
       LAY=LAYER(J)
      ELSE
       LAY=IQ3+(LAYER(J)-1)
      ENDIF
      IF(COND(J).EQ.1) GO TO 1210
C
      DO 1205 KKK=IQ3,LAY
      K=KKK
      IF(IWELTYPE.EQ.2)THEN
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      PWFC(J,KKK)=0.0
      IF(PID(J,KKK).LE.0.0001) GO TO 1205
      PP=P(IQ1,IQ2,K)
      BPT=PBOT(IQ1,IQ2,K)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PP,BBO)
      CALL INTERP(PWT,BWT,MPWT,PP,BBW)
      CALL INTERP(PGT,BGT,MPGT,PP,BBG)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PP,RSO)
      CALL INTERP(PWT,RSWT,MPWT,PP,RSW)
      FAC=PID(J,KKK)*5.615
      GMTB=GMO(J,KKK)/BBO+GMW(J,KKK)/BBW+GMG(J,KKK)/BBG
      SOLN=RSO*QO(IQ1,IQ2,K)+RSW*QW(IQ1,IQ2,K)
      QT=QO(IQ1,IQ2,K)+QW(IQ1,IQ2,K)+QG(IQ1,IQ2,K)
      PWFC(J,KKK)=PP-(QT-SOLN)/(FAC*GMTB)
 1205 CONTINUE
 1210 CONTINUE
C**** PRESSURE CONSTRAINT
      IF(IWELTYPE.EQ.1)THEN
       N1=1
       N2=NVQN
      ELSE
       N1=NVQN+1
       N2=NVQN+NVQNH
      ENDIF
C
      DO 1350 J=N1,N2
      IF(KIP(J).GE.0) GO TO 1345
      IF(IWELTYPE.EQ.1)THEN
       IQ1=IQN1(J)
       IQ2=IQN2(J)
       IQ3=IQN3(J)
       LAY=IQ3+(LAYER(J)-1)
      ELSE
       IQ3=1
       LAY=LAYER(J)
      ENDIF
      QOSUM=0.
      QTSUM=0.
      QGSUM=0.
      DO 1341 KKK=IQ3,LAY
       K=KKK
      IF(IWELTYPE.EQ.2)THEN
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      PPN=PN(IQ1,IQ2,K)
      BPT=PBOT(IQ1,IQ2,K)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PPN,BBO)
      CALL INTERP(PWT,BWT,MPWT,PPN,BBW)
      CALL INTERP(PGT,BGT,MPGT,PPN,BBG)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PPN,RSO)
      CALL INTERP(PWT,RSWT,MPWT,PPN,RSW)
C**** OIL PRODUCER
      IF(KIP(J).NE.-1) GO TO 1310
      QO(IQ1,IQ2,K)=PID(J,KKK)*5.615*GMO(J,KKK)
     & *(PPN-PWF(J,KKK))/BBO
      IF(PPN.LE.PWF(J,KKK)) QO(IQ1,IQ2,K)=0.0
      QW(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*GMW(J,KKK)*BBO
     & /(BBW*GMO(J,KKK))
      QG(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*(GMG(J,KKK)*BBO
     & /(BBG*GMO(J,KKK))+RSO)+RSW*QW(IQ1,IQ2,K)
      QOSUM=QOSUM+QO(IQ1,IQ2,K)
      QTSUM=QTSUM+QO(IQ1,IQ2,K)+QW(IQ1,IQ2,K)+QG(IQ1,IQ2,K)
      GO TO 1340
C**** WATER INJECTOR
 1310 IF(KIP(J).NE.-2) GO TO 1320
      QW(IQ1,IQ2,K)=PID(J,KKK)*5.615*(GMO(J,KKK)
     & +GMW(J,KKK)+GMG(J,KKK))*(PPN-PWF(J,KKK))/BBW
      IF(PPN.GE.PWF(J,KKK)) QW(IQ1,IQ2,K)=0.0
      GO TO 1340
C**** GAS INJECTOR
 1320 IF(KIP(J).NE.-3) GO TO 1340
      QG(IQ1,IQ2,K)=PID(J,KKK)*5.615*(GMO(J,KKK)
     & +GMW(J,KKK)+GMG(J,KKK))*(PPN-PWF(J,KKK))/BBG
      IF(PPN.GE.PWF(J,KKK)) QG(IQ1,IQ2,K)=0.0
      QGSUM=QGSUM+QG(IQ1,IQ2,K)
 1340 CONTINUE
 1341 CONTINUE
C**** CHECK RATE CONSTRAINTS ALSO ****************
      IF(QVO(J).GT.0.AND.QOSUM.GT.QVO(J)*5.615)THEN
        QORATIO=QVO(J)*5.615/QOSUM
        DO 1351 KKK=IQ3,LAY
          K=KKK
         IF(IWELTYPE.EQ.2)THEN
          IQ1=IQH1(J,KKK)
          IQ2=IQH2(J,KKK)
          K=IQH3(J,KKK)
         ENDIF
         QO(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*QORATIO
         PWFC(J,KKK)=PP-(PP-PWF(J,KKK))*QORATIO
 1351   CONTINUE
       ENDIF
      IF(QVT(J).GT.0.AND.QTSUM.GT.QVT(J)*5.615)THEN
        QTRATIO=QVT(J)*5.615/QTSUM
        DO 1352 KKK=IQ3,LAY
          K=KKK
         IF(IWELTYPE.EQ.2)THEN
          IQ1=IQH1(J,KKK)
          IQ2=IQH2(J,KKK)
          K=IQH3(J,KKK)
         ENDIF
         QO(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*QTRATIO
         QW(IQ1,IQ2,K)=QW(IQ1,IQ2,K)*QTRATIO
         QG(IQ1,IQ2,K)=QG(IQ1,IQ2,K)*QTRATIO
         PWFC(J,KKK)=PP-(PP-PWF(J,KKK))*QTRATIO
 1352   CONTINUE
       ENDIF
      IF(QVG(J).EQ.0.OR.QGSUM.EQ.0)GO TO 1345
        QGRATIO=QVG(J)*1000/QGSUM
       IF(QGRATIO.LT.1.)THEN
        DO 1353 KKK=IQ3,LAY
          K=KKK
         IF(IWELTYPE.EQ.2)THEN
          IQ1=IQH1(J,KKK)
          IQ2=IQH2(J,KKK)
          K=IQH3(J,KKK)
         ENDIF
         QG(IQ1,IQ2,K)=QG(IQ1,IQ2,K)*QGRATIO
         PWFC(J,KKK)=PP-(PP-PWF(J,KKK))*QGRATIO
 1353   CONTINUE
       ENDIF
 1345 CONTINUE
 1350 CONTINUE
      RETURN
      END
C
C      SUBROUTINE SOLMAT(DIV1,D288,N,NX,NY,NZ)
      SUBROUTINE SOLMAT(DIV1,D288,N,NX,NY,NZ,PBOT,
     &  AW,AE,AN,AS,AB,AT,E,B,SUM,GAM,EL,TX,TY,TZ,P,SO,
     &  SW,SG,BO,BW,BG,QO,QW,QG,GOWT,GWWT,GGWT,QOWG,OW,OE,
     &  WW,WE,OS,ON,WS,WN,OT,OB,WT,WB,CT,VP)
C
C=================================================================SOLMAT
C
C   COMPUTE PRESSURE EQUATION COEFFICENTS.
C
C=======================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
C	NEW DP LOCAL VARIABLES **RJH**
C       MODIFICATION FOR DOUBLE PRECISION READ OF SO & SW IN RESTART  **RJH**

      REAL KROT,KRWT,KRGT,MUOT,MUWT,MUGT    
      REAL MUO4,MUW4,MUG4,MUO5,MUW5,MUG5,MUO6,MUW6,MUG6
     &,KRO1,KRW1,KRG1,KRO2,KRW2,KRG2,KRO3,KRW3,KRG3
     &,MUO1,MUW1,MUG1,MUO2,MUW2,MUG2,MUO3,MUW3,MUG3
     &,KRO4,KRW4,KRG4,KRO5,KRW5,KRG5,KRO6,KRW6,KRG6
     &,MO1,MW1,MG1,MO2,MW2,MG2,MO3,MW3,MG3
     &,MO4,MW4,MG4,MO5,MW5,MG5,MO6,MW6,MG6
     &,MUO,MUW,MUG
C
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      COMMON /CODE/ KSM1,KSN1,KCO1,NN,FACT1,FACT2,TMAX,KSOL,MITER,
     & OMEGA,TOL,TOL1,KSN,KSM,KCO,KTR,KCOFF,DSMAX,DPMAX,WORMAX,
     & GORMAX,PAMIN,PAMAX
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
      DIMENSION PBOT(NX,NY,NZ)
      DIMENSION AW(NX,NY,NZ),AE(NX,NY,NZ),AN(NX,NY,NZ),
     & AS(NX,NY,NZ),AB(NX,NY,NZ),AT(NX,NY,NZ),E(NX,NY,NZ),B(NX,NY,NZ)
C      DIMENSION PN(NX,NY,NZ),
C     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
      DIMENSION SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      DIMENSION EL(NX,NY,NZ)
      DIMENSION TX(NX+1,NY,NZ),TY(NX,NY+1,NZ),TZ(NX,NY,NZ+1)
      DIMENSION P(NX,NY,NZ)
      DIMENSION SO(NX,NY,NZ),SW(NX,NY,NZ),SG(NX,NY,NZ)
C      REAL*8 DPP(NX,NY,NZ),DPSO(NX,NY,NZ),DPSW(NX,NY,NZ)
      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      DIMENSION BO(NX,NY,NZ),BW(NX,NY,NZ),BG(NX,NY,NZ)
      DIMENSION QO(NX,NY,NZ),QW(NX,NY,NZ),QG(NX,NY,NZ)
      DIMENSION GOWT(NX,NY,NZ),GWWT(NX,NY,NZ),GGWT(NX,NY,NZ),
     & OW(NX+1,NY,NZ),OE(NX+1,NY,NZ),WW(NX+1,NY,NZ),WE(NX+1,NY,NZ),
     & OS(NX,NY+1,NZ),ON(NX,NY+1,NZ),WS(NX,NY+1,NZ),WN(NX,NY+1,NZ),
     & OT(NX,NY,NZ+1),OB(NX,NY,NZ+1),WT(NX,NY,NZ+1),WB(NX,NY,NZ+1)
      DIMENSION QOWG(NX,NY,NZ)
      DIMENSION CT(NX,NY,NZ)
      DIMENSION VP(NX,NY,NZ)
      DATA RSO1,RSO2,RSO3,RSO4,RSO5,RSO6/6*0.0/
      DATA RSW1,RSW2,RSW3,RSW4,RSW5,RSW6/6*0.0/
C
      DO 200 K=1,KK
      DO 200 J=1,JJ
      DO 200 I=1,II
C
      PP=P(I,J,K)
      BPT=PBOT(I,J,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PP,RSO)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,PP,MUO)
      CALL INTERP(PWT,RSWT,MPWT,PP,RSW)
      CALL INTERP(PWT,MUWT,MPWT,PP,MUW)
      CALL INTERP(PGT,MUGT,MPGT,PP,MUG)
      SSO=SO(I,J,K)
      SSW=SW(I,J,K)
      SSG=SG(I,J,K)
      CALL INTERP(SAT,PCOWT,MSAT,SSW,PCOW)
      CALL INTERP(SAT,PCGOT,MSAT,SSG,PCGO)
C
      RO=(RHOSCO + RSO*RHOSCG)/BO(I,J,K)
      RW=(RHOSCW + RSW*RHOSCG)/BW(I,J,K)
      RG=RHOSCG/BG(I,J,K)
C
      IF(I.EQ.1)GO TO 115
      P1=P(I-1,J,K)
      BPT=PBOT(I-1,J,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,P1,RSO1)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,P1,MUO1)
      CALL INTERP(PWT,RSWT,MPWT,P1,RSW1)
      CALL INTERP(PWT,MUWT,MPWT,P1,MUW1)
      CALL INTERP(PGT,MUGT,MPGT,P1,MUG1)
      SO1S=SO(I-1,J,K)
      SW1S=SW(I-1,J,K)
      SG1S=SG(I-1,J,K)
      CALL INTERP(SAT,PCOWT,MSAT,SW1S,PCOW1)
      CALL INTERP(SAT,PCGOT,MSAT,SG1S,PCGO1)
C
C
      RO1=(RHOSCO + RSO1*RHOSCG)/BO(I-1,J,K)
      RW1=(RHOSCW + RSW1*RHOSCG)/BW(I-1,J,K)
      RG1=RHOSCG/BG(I-1,J,K)
C
      FACT=-D288*(EL(I-1,J,K)-EL(I,J,K))
      GOW1=(RO1+RO)*FACT
      GWW1=(RW1+RW)*FACT + PCOW-PCOW
      GGW1=(RG1+RG)*FACT + PCGO1-PCGO
C
      P11=P1-PP
      HO1=P11+GOW1
      HW1=P11+GWW1
      HG1=P11+GGW1
C
      IF(HO1.GE.0.)CALL INTERP(SAT,KROT,MSAT,SO1S,KRO1)
      IF(HO1.LT.0.)CALL INTERP(SAT,KROT,MSAT,SSO,KRO1)
      IF(HW1.GE.0.)CALL INTERP(SAT,KRWT,MSAT,SW1S,KRW1)
      IF(HW1.LT.0.)CALL INTERP(SAT,KRWT,MSAT,SSW,KRW1)
      IF(HG1.GE.0.)CALL INTERP(SAT,KRGT,MSAT,SG1S,KRG1)
      IF(HG1.LT.0.)CALL INTERP(SAT,KRGT,MSAT,SSG,KRG1)
      MO1=4.0*KRO1/((BO(I-1,J,K)+BO(I,J,K)) * (MUO1+MUO))
      MW1=4.0*KRW1/((BW(I-1,J,K)+BW(I,J,K)) * (MUW1+MUW))
      MG1=4.0*KRG1/((BG(I-1,J,K)+BG(I,J,K)) * (MUG1+MUG))
C
115   AOW=TX(I,J,K)*MO1
      AWW=TX(I,J,K)*MW1
      AGW=TX(I,J,K)*MG1
C
      IF(I.EQ.II)GO TO 125
      P2=P(I+1,J,K)
      BPT=PBOT(I+1,J,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,P2,RSO2)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,P2,MUO2)
      CALL INTERP(PWT,RSWT,MPWT,P2,RSW2)
      CALL INTERP(PWT,MUWT,MPWT,P2,MUW2)
      CALL INTERP(PGT,MUGT,MPGT,P2,MUG2)
      SO2=SO(I+1,J,K)
      SW2=SW(I+1,J,K)
      SG2=SG(I+1,J,K)
      CALL INTERP(SAT,PCOWT,MSAT,SW2,PCOW2)
      CALL INTERP(SAT,PCGOT,MSAT,SG2,PCGO2)
C
      RO2=(RHOSCO + RSO2*RHOSCG)/BO(I+1,J,K)
      RW2=(RHOSCW + RSW2*RHOSCG)/BW(I+1,J,K)
      RG2=RHOSCG/BG(I+1,J,K)
C
      FACT=-D288*(EL(I+1,J,K)-EL(I,J,K))
      GOW2=(RO2+RO)*FACT
      GWW2=(RW2+RW)*FACT + PCOW-PCOW2
      GGW2=(RG2+RG)*FACT + PCGO2-PCGO
C
      P22=P2-PP
      HO2=P22+GOW2
      HW2=P22+GWW2
      HG2=P22+GGW2
      IF(HO2.GE.0.)CALL INTERP(SAT,KROT,MSAT,SO2,KRO2)
      IF(HO2.LT.0.)CALL INTERP(SAT,KROT,MSAT,SSO,KRO2)
      IF(HW2.GE.0.)CALL INTERP(SAT,KRWT,MSAT,SW2,KRW2)
      IF(HW2.LT.0.)CALL INTERP(SAT,KRWT,MSAT,SSW,KRW2)
      IF(HG2.GE.0.)CALL INTERP(SAT,KRGT,MSAT,SG2,KRG2)
      IF(HG2.LT.0.)CALL INTERP(SAT,KRGT,MSAT,SSG,KRG2)
      MO2=4.0*KRO2/((BO(I+1,J,K)+BO(I,J,K)) * (MUO2+MUO))
      MW2=4.0*KRW2/((BW(I+1,J,K)+BW(I,J,K)) * (MUW2+MUW))
      MG2=4.0*KRG2/((BG(I+1,J,K)+BG(I,J,K)) * (MUG2+MUG))
C
125   AOE=TX(I+1,J,K)*MO2
      AWE=TX(I+1,J,K)*MW2
      AGE=TX(I+1,J,K)*MG2
C
C
      IF(J.EQ.1)GO TO 135
      P3=P(I,J-1,K)
      BPT=PBOT(I,J-1,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,P3,RSO3)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,P3,MUO3)
      CALL INTERP(PWT,RSWT,MPWT,P3,RSW3)
      CALL INTERP(PWT,MUWT,MPWT,P3,MUW3)
      CALL INTERP(PGT,MUGT,MPGT,P3,MUG3)
      SO3=SO(I,J-1,K)
      SW3=SW(I,J-1,K)
      SG3=SG(I,J-1,K)
      CALL INTERP(SAT,PCOWT,MSAT,SW3,PCOW3)
      CALL INTERP(SAT,PCGOT,MSAT,SG3,PCGO3)
C
      RO3=(RHOSCO + RSO3*RHOSCG)/BO(I,J-1,K)
      RW3=(RHOSCW + RSW3*RHOSCG)/BW(I,J-1,K)
      RG3=RHOSCG/BG(I,J-1,K)
C
      FACT=-D288*(EL(I,J-1,K)-EL(I,J,K))
      GOW3=(RO3+RO)*FACT
      GWW3=(RW3+RW)*FACT + PCOW-PCOW3
      GGW3=(RG3+RG)*FACT + PCGO3-PCGO
C
      P33=P3-PP
      HO3=P33+GOW3
      HW3=P33+GWW3
      HG3=P33+GGW3
      IF(HO3.GE.0.)CALL INTERP(SAT,KROT,MSAT,SO3,KRO3)
      IF(HO3.LT.0.)CALL INTERP(SAT,KROT,MSAT,SSO,KRO3)
      IF(HW3.GE.0.)CALL INTERP(SAT,KRWT,MSAT,SW3,KRW3)
      IF(HW3.LT.0.)CALL INTERP(SAT,KRWT,MSAT,SSW,KRW3)
      IF(HG3.GE.0.)CALL INTERP(SAT,KRGT,MSAT,SG3,KRG3)
      IF(HG3.LT.0.)CALL INTERP(SAT,KRGT,MSAT,SSG,KRG3)
      MO3=4.0*KRO3/((BO(I,J-1,K)+BO(I,J,K)) * (MUO3+MUO))
      MW3=4.0*KRW3/((BW(I,J-1,K)+BW(I,J,K)) * (MUW3+MUW))
      MG3=4.0*KRG3/((BG(I,J-1,K)+BG(I,J,K)) * (MUG3+MUG))
C
135   AOS=TY(I,J,K)*MO3
      AWS=TY(I,J,K)*MW3
      AGS=TY(I,J,K)*MG3
C
      IF(J.EQ.JJ)GO TO 140
      P4=P(I,J+1,K)
      BPT=PBOT(I,J+1,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,P4,RSO4)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,P4,MUO4)
      CALL INTERP(PWT,RSWT,MPWT,P4,RSW4)
      CALL INTERP(PWT,MUWT,MPWT,P4,MUW4)
      CALL INTERP(PGT,MUGT,MPGT,P4,MUG4)
      SO4=SO(I,J+1,K)
      SW4=SW(I,J+1,K)
      SG4=SG(I,J+1,K)
      CALL INTERP(SAT,PCOWT,MSAT,SW4,PCOW4)
      CALL INTERP(SAT,PCGOT,MSAT,SG4,PCGO4)
C
      RO4=(RHOSCO + RSO4*RHOSCG)/BO(I,J+1,K)
      RW4=(RHOSCW + RSW4*RHOSCG)/BW(I,J+1,K)
      RG4=RHOSCG/BG(I,J+1,K)
C
      FACT=-D288*(EL(I,J+1,K)-EL(I,J,K))
      GOW4=(RO4+RO)*FACT
      GWW4=(RW4+RW)*FACT + PCOW-PCOW4
      GGW4=(RG4+RG)*FACT + PCGO4-PCGO
C
      P44=P4-PP
      HO4=P44+GOW4
      HW4=P44+GWW4
      HG4=P44+GGW4
      IF(HO4.GE.0.)CALL INTERP(SAT,KROT,MSAT,SO4,KRO4)
      IF(HO4.LT.0.)CALL INTERP(SAT,KROT,MSAT,SSO,KRO4)
      IF(HW4.GE.0.)CALL INTERP(SAT,KRWT,MSAT,SW4,KRW4)
      IF(HW4.LT.0.)CALL INTERP(SAT,KRWT,MSAT,SSW,KRW4)
      IF(HG4.GE.0.)CALL INTERP(SAT,KRGT,MSAT,SG4,KRG4)
      IF(HG4.LT.0.)CALL INTERP(SAT,KRGT,MSAT,SSG,KRG4)
      MO4=4.0*KRO4/((BO(I,J+1,K)+BO(I,J,K)) * (MUO4+MUO))
      MW4=4.0*KRW4/((BW(I,J+1,K)+BW(I,J,K)) * (MUW4+MUW))
      MG4=4.0*KRG4/((BG(I,J+1,K)+BG(I,J,K)) * (MUG4+MUG))
C
140   AON=TY(I,J+1,K)*MO4
      AWN=TY(I,J+1,K)*MW4
      AGN=TY(I,J+1,K)*MG4
C
C
      IF(K.EQ.1)GO TO 145
      P5=P(I,J,K-1)
      BPT=PBOT(I,J,K-1)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,P5,RSO5)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,P5,MUO5)
      CALL INTERP(PWT,RSWT,MPWT,P5,RSW5)
      CALL INTERP(PWT,MUWT,MPWT,P5,MUW5)
      CALL INTERP(PGT,MUGT,MPGT,P5,MUG5)
      SO5=SO(I,J,K-1)
      SW5=SW(I,J,K-1)
      SG5=SG(I,J,K-1)
      CALL INTERP(SAT,PCOWT,MSAT,SW5,PCOW5)
      CALL INTERP(SAT,PCGOT,MSAT,SG5,PCGO5)
      RO5=(RHOSCO + RSO5*RHOSCG)/BO(I,J,K-1)
      RW5=(RHOSCW + RSW5*RHOSCG)/BW(I,J,K-1)
      RG5=RHOSCG/BG(I,J,K-1)
C
      FACT=-D288*(EL(I,J,K-1)-EL(I,J,K))
      GOW5=(RO5+RO)*FACT
      GWW5=(RW5+RW)*FACT + PCOW-PCOW5
      GGW5=(RG5+RG)*FACT + PCGO5-PCGO
      P55=P5-PP
      HO5=P55+GOW5
      HW5=P55+GWW5
      HG5=P55+GGW5
      IF(HO5.GE.0.)CALL INTERP(SAT,KROT,MSAT,SO5,KRO5)
      IF(HO5.LT.0.)CALL INTERP(SAT,KROT,MSAT,SSO,KRO5)
      IF(HW5.GE.0.)CALL INTERP(SAT,KRWT,MSAT,SW5,KRW5)
      IF(HW5.LT.0.)CALL INTERP(SAT,KRWT,MSAT,SSW,KRW5)
      IF(HG5.GE.0.)CALL INTERP(SAT,KRGT,MSAT,SG5,KRG5)
      IF(HG5.LT.0.)CALL INTERP(SAT,KRGT,MSAT,SSG,KRG5)
      MO5=4.0*KRO5/((BO(I,J,K-1)+BO(I,J,K)) * (MUO5+MUO))
      MW5=4.0*KRW5/((BW(I,J,K-1)+BW(I,J,K)) * (MUW5+MUW))
      MG5=4.0*KRG5/((BG(I,J,K-1)+BG(I,J,K)) * (MUG5+MUG))
C
C
145   AOT=TZ(I,J,K)*MO5
      AWT=TZ(I,J,K)*MW5
      AGT=TZ(I,J,K)*MG5
C
      IF(K.EQ.KK)GO TO 150
      P6=P(I,J,K+1)
      BPT=PBOT(I,J,K+1)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,P6,RSO6)
      CALL INTPVT(BPT,VSLOPE,POT,MUOT,MPOT,P6,MUO6)
      CALL INTERP(PWT,RSWT,MPWT,P6,RSW6)
      CALL INTERP(PWT,MUWT,MPWT,P6,MUW6)
      CALL INTERP(PGT,MUGT,MPGT,P6,MUG6)
      SO6=SO(I,J,K+1)
      SW6=SW(I,J,K+1)
      SG6=SG(I,J,K+1)
C
      CALL INTERP(SAT,PCOWT,MSAT,SW6,PCOW6)
      CALL INTERP(SAT,PCGOT,MSAT,SG6,PCGO6)
      RO6=(RHOSCO + RSO6*RHOSCG)/BO(I,J,K+1)
      RW6=(RHOSCW + RSW6*RHOSCG)/BW(I,J,K+1)
      RG6=RHOSCG/BG(I,J,K+1)
C
      FACT=-D288*(EL(I,J,K+1)-EL(I,J,K))
      GOW6=(RO6+RO)*FACT
      GWW6=(RW6+RW)*FACT + PCOW-PCOW6
      GGW6=(RG6+RG)*FACT + PCGO6-PCGO
C
      P66=P6-PP
      HO6=P66+GOW6
      HW6=P66+GWW6
      HG6=P66+GGW6
      IF(HO6.GE.0.)CALL INTERP(SAT,KROT,MSAT,SO6,KRO6)
      IF(HO6.LT.0.)CALL INTERP(SAT,KROT,MSAT,SSO,KRO6)
      IF(HW6.GE.0.)CALL INTERP(SAT,KRWT,MSAT,SW6,KRW6)
      IF(HW6.LT.0.)CALL INTERP(SAT,KRWT,MSAT,SSW,KRW6)
      IF(HG6.GE.0.)CALL INTERP(SAT,KRGT,MSAT,SG6,KRG6)
      IF(HG6.LT.0.)CALL INTERP(SAT,KRGT,MSAT,SSG,KRG6)
      MO6=4.0*KRO6/((BO(I,J,K+1)+BO(I,J,K)) * (MUO6+MUO))
      MW6=4.0*KRW6/((BW(I,J,K+1)+BW(I,J,K)) * (MUW6+MUW))
      MG6=4.0*KRG6/((BG(I,J,K+1)+BG(I,J,K)) * (MUG6+MUG))
C
150   AOB=TZ(I,J,K+1)*MO6
      AWB=TZ(I,J,K+1)*MW6
      AGB=TZ(I,J,K+1)*MG6
C
C
      RSO1A=0.5*(RSO1+RSO)
      RSO2A=0.5*(RSO2+RSO)
      RSO3A=0.5*(RSO3+RSO)
      RSO4A=0.5*(RSO4+RSO)
      RSO5A=0.5*(RSO5+RSO)
      RSO6A=0.5*(RSO6+RSO)
      RSW1A=0.5*(RSW1+RSW)
      RSW2A=0.5*(RSW2+RSW)
      RSW3A=0.5*(RSW3+RSW)
      RSW4A=0.5*(RSW4+RSW)
      RSW5A=0.5*(RSW5+RSW)
      RSW6A=0.5*(RSW6+RSW)
      AO1=AOW*GOW1
      AO2=AOE*GOW2
      AO3=AOS*GOW3
      AO4=AON*GOW4
      AO5=AOT*GOW5
      AO6=AOB*GOW6
      AW1=AWW*GWW1
      AW2=AWE*GWW2
      AW3=AWS*GWW3
      AW4=AWN*GWW4
      AW5=AWT*GWW5
      AW6=AWB*GWW6
C
      GOWT(I,J,K)= AO1 + AO2 + AO3 + AO4 + AO5 + AO6
      GWWT(I,J,K)= AW1 + AW2 + AW3 + AW4 + AW5 + AW6
      GGWT(I,J,K)=AGW*GGW1+AGE*GGW2+AGS*GGW3+AGN*GGW4+AGT*GGW5+AGB*GGW6
     &+RSO1A*AO1+RSO2A*AO2+RSO3A*AO3+RSO4A*AO4+RSO5A*AO5+RSO6A*AO6
     &+RSW1A*AW1+RSW2A*AW2+RSW3A*AW3+RSW4A*AW4+RSW5A*AW5+RSW6A*AW6
C
      QOWG(I,J,K)=(BO(I,J,K)-BG(I,J,K)*RSO)*(-GOWT(I,J,K)+QO(I,J,K)) +
     &           (BW(I,J,K)-BG(I,J,K)*RSW)*(-GWWT(I,J,K)+QW(I,J,K)) +
     &            BG(I,J,K)*(-GGWT(I,J,K)+QG(I,J,K))
C
      AW(I,J,K)=(BO(I,J,K) + 0.5*BG(I,J,K)*(RSO1-RSO)) * AOW +
     &         (BW(I,J,K) + 0.5*BG(I,J,K)*(RSW1-RSW)) * AWW +
     &          BG(I,J,K)*AGW
      AE(I,J,K)=(BO(I,J,K) + 0.5*BG(I,J,K)*(RSO2-RSO)) * AOE +
     &         (BW(I,J,K) + 0.5*BG(I,J,K)*(RSW2-RSW)) * AWE +
     &          BG(I,J,K)*AGE
      AS(I,J,K)=(BO(I,J,K) + 0.5*BG(I,J,K)*(RSO3-RSO)) * AOS +
     &         (BW(I,J,K) + 0.5*BG(I,J,K)*(RSW3-RSW)) * AWS +
     &          BG(I,J,K)*AGS
      AN(I,J,K)=(BO(I,J,K) + 0.5*BG(I,J,K)*(RSO4-RSO)) * AON +
     &         (BW(I,J,K) + 0.5*BG(I,J,K)*(RSW4-RSW)) * AWN +
     &          BG(I,J,K)*AGN
      AT(I,J,K)=(BO(I,J,K) + 0.5*BG(I,J,K)*(RSO5-RSO)) * AOT +
     &         (BW(I,J,K) + 0.5*BG(I,J,K)*(RSW5-RSW)) * AWT +
     &          BG(I,J,K)*AGT
      AB(I,J,K)=(BO(I,J,K) + 0.5*BG(I,J,K)*(RSO6-RSO)) * AOB +
     &         (BW(I,J,K) + 0.5*BG(I,J,K)*(RSW6-RSW)) * AWB +
     &          BG(I,J,K)*AGB
      OW(I,J,K)=AOW
      OE(I,J,K)=AOE
      OS(I,J,K)=AOS
      ON(I,J,K)=AON
      OT(I,J,K)=AOT
      OB(I,J,K)=AOB
      WW(I,J,K)=AWW
      WE(I,J,K)=AWE
      WS(I,J,K)=AWS
      WN(I,J,K)=AWN
      WT(I,J,K)=AWT
      WB(I,J,K)=AWB
      IF(KCOFF.EQ.0)GO TO 200
      WRITE(NO,33)
      WRITE(NO,2)I,J,K,MO1,MO2,MO3,MO4,MO5,MO6
      WRITE(NO,2)I,J,K,MW1,MW2,MW3,MW4,MW5,MW6
      WRITE(NO,2)I,J,K,MG1,MG2,MG3,MG4,MG5,MG6
      WRITE(NO,2)I,J,K,AOW,AOE,AOS,AON,AOT,AOB,BO(I,J,K),RSO
      WRITE(NO,2)I,J,K,AWW,AWE,AWS,AWN,AWT,AWB,BW(I,J,K),RSW
      WRITE(NO,2)I,J,K,AGW,AGE,AGS,AGN,AGT,AGB,BG(I,J,K)
      WRITE(NO,2)I,J,K,GOWT(I,J,K),QO(I,J,K),GWWT(I,J,K),QW(I,J,K),
     &            GGWT(I,J,K),QG(I,J,K),QOWG(I,J,K)
200   CONTINUE
C
C**** CALCULATE MAIN DIAGONAL AND RHS VECTOR
C
      DO 300 K=1,KK
      DO 300 J=1,JJ
      DO 300 I=1,II
C
      SUM(I,J,K)=AW(I,J,K)+AE(I,J,K)+AS(I,J,K)+AN(I,J,K)+
     &    AT(I,J,K)+AB(I,J,K)
      GAM(I,J,K)=VP(I,J,K)*CT(I,J,K)*DIV1
      E(I,J,K)=-SUM(I,J,K) - GAM(I,J,K)
      B(I,J,K)= QOWG(I,J,K) - GAM(I,J,K)*P(I,J,K)
300   CONTINUE
      IF(KSM1.EQ.0)RETURN
      IF(N.NE.1.AND.N.NE.NN.AND.N.NE.KSM)RETURN
      WRITE(NO,4)
      DO 404 K=1,KK
      DO 404 J=1,JJ
      DO 404 I=1,II
      WRITE(NO,2)I,J,K,AT(I,J,K),AS(I,J,K),AW(I,J,K),E(I,J,K),
     & AE(I,J,K)   ,AN(I,J,K),AB(I,J,K),B(I,J,K)
404   CONTINUE
4     FORMAT(//T3,'NODE      AT(I,J,K)     AS(I,J,K)      AW(I,J,K)',
     &'      E(I,J,K)     AE(I,J,K)      AN(I,J,K)      AB(I,J,K)',
     &'      B(I,J,K)'/)
2     FORMAT(1X,3I3,8E15.6)
33    FORMAT(//)
C
      RETURN
      END

C      SUBROUTINE PRATEI(IWELTYPE,NVQN,NVQNH,NX,NY,NZ)
      SUBROUTINE PRATEI(IWELTYPE,NVQN,NVQNH,NX,NY,NZ,PBOT,E,B,PN)
C
C==================================================================QRATE
C
C   COMPUTE WELL RATES WHEN RATE OR EXPLICIT PRESSURE CONSTRAIN
C   ACTIVATED;
C        PRATEI: REVISES PRESSURE EQUATION COEFFICENTS WHEN IMPLICIT
C                PRESSURE CONSTRAIN ACTIVATED;
C        PRATEO: COMPUTE RATES WHEN IMPLICIT PRESSURE CONSTRAIN
C                ACTIVATE.
C
C       (IWELTYPE = 1 FOR VERTICAL WELL, 2 FOR HORIZON./SLANT WELL)
C=====================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      REAL KROT,KRWT,KRGT,MUWT,MUOT,MUGT
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
      DIMENSION PBOT(NX,NY,NZ)
C      DIMENSION AW(NX,NY,NZ),AE(NX,NY,NZ),AN(NX,NY,NZ),
C     & AS(NX,NY,NZ),AB(NX,NY,NZ),AT(NX,NY,NZ)
      DIMENSION E(NX,NY,NZ),B(NX,NY,NZ)
      DIMENSION PN(NX,NY,NZ)
C     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
C     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
C
      IF(IWELTYPE.EQ.1)THEN
       N1=1
       N2=NVQN
      ELSE
       N1=NVQN+1
       N2=NVQN+NVQNH
      ENDIF
C
      DO 205 J=N1,N2
      IF(KIP(J).GE.-10) GO TO 205
      IF(IWELTYPE.EQ.1)THEN
       IQ1=IQN1(J)
       IQ2=IQN2(J)
       IQ3=IQN3(J)
       LAY=IQ3+(LAYER(J)-1)
      ELSE
       IQ3=1
       LAY=LAYER(J)
      ENDIF
      DO 203 KKK=IQ3,LAY
      K=KKK
      IF(IWELTYPE.EQ.2)THEN
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      P56=PID(J,KKK)*5.615
      PPN=PN(IQ1,IQ2,K)
      BPT=PBOT(IQ1,IQ2,K)
      CALL INTPVT(BPT,BSLOPE,POT,BOT,MPOT,PPN,BBO)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PPN,RSO)
      CALL INTERP(PWT,BWT,MPWT,PPN,BBW)
      CALL INTERP(PGT,BGT,MPGT,PPN,BBG)
      CALL INTERP(PWT,RSWT,MPWT,PPN,RSW)
      CPIO=GMO(J,KKK)*P56*(BBO-BBG*RSO)/BBO
      CPIW=GMW(J,KKK)*P56*(BBW-BBG*RSW)/BBW
      CPIG=GMG(J,KKK)*P56
      CPI=CPIO+CPIW+CPIG
      B(IQ1,IQ2,K)=B(IQ1,IQ2,K)-CPI*PWF(J,KKK)
      E(IQ1,IQ2,K)=E(IQ1,IQ2,K)-CPI
  203 CONTINUE
  205 CONTINUE
      RETURN
      END
C
C      SUBROUTINE PRATEO(IWELTYPE,NVQN,NVQNH,NX,NY,NZ)
      SUBROUTINE PRATEO(IWELTYPE,NVQN,NVQNH,NX,NY,NZ,PBOT,P,PN,
     &  BO,BW,BG,QO,QW,QG)
C
C==================================================================QRATE
C
C   COMPUTE WELL RATES WHEN RATE OR EXPLICIT PRESSURE CONSTRAIN
C   ACTIVATED;
C        PRATEI: REVISES PRESSURE EQUATION COEFFICENTS WHEN IMPLICIT
C                PRESSURE CONSTRAIN ACTIVATED;
C        PRATEO: COMPUTE RATES WHEN IMPLICIT PRESSURE CONSTRAIN
C                ACTIVATE.
C
C        (IWELTYPE = 1 FOR VERTICAL WELL, 2 FOR HORIZON/SLANT WELL)
C=====================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      REAL KROT,KRWT,KRGT,MUWT,MUOT,MUGT
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
      DIMENSION PBOT(NX,NY,NZ)
C      DIMENSION AW(NX,NY,NZ),AE(NX,NY,NZ),AN(NX,NY,NZ),
C     & AS(NX,NY,NZ),AB(NX,NY,NZ),AT(NX,NY,NZ),E(NX,NY,NZ),B(NX,NY,NZ)
      DIMENSION PN(NX,NY,NZ)
C     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
C     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      DIMENSION P(NX,NY,NZ)
      COMMON /SPVT/ SAT(NTE),KROT(NTE),KRWT(NTE),KRGT(NTE),PCOWT(NTE),
     & PCGOT(NTE),POT(NTE),MUOT(NTE),BOT(NTE),BOPT(NTE),RSOT(NTE),RSOPT
     & (NTE),PWT(NTE),MUWT(NTE),BWT(NTE),BWPT(NTE),RSWT(NTE),RSWPT(NTE),
     & PGT(NTE),MUGT(NTE),BGT(NTE),BGPT(NTE),CRT(NTE)
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      DIMENSION BO(NX,NY,NZ),BW(NX,NY,NZ),BG(NX,NY,NZ)
      DIMENSION  QO(NX,NY,NZ),QW(NX,NY,NZ),QG(NX,NY,NZ)
      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
C
      IF(IWELTYPE.EQ.1)THEN 
       N1=1
       N2=NVQN
      ELSE
       N1=NVQN+1
       N2=NVQN+NVQNH
      ENDIF
C
      DO 2059 J=N1,N2
      IF(KIP(J).GE.-10) GO TO 2059
      IF(IWELTYPE.EQ.1)THEN
       IQ1=IQN1(J)
       IQ2=IQN2(J)
       IQ3=IQN3(J)
       LAY=IQ3+(LAYER(J)-1)
      ELSE
       IQ3=1
       LAY=LAYER(J)
      ENDIF
       QOSUM=0.
       QTSUM=0.
       QGSUM=0.
      DO 2057 KKK=IQ3,LAY
      K=KKK
      IF(IWELTYPE.EQ.2)THEN
       IQ1=IQH1(J,KKK)
       IQ2=IQH2(J,KKK)
       K=IQH3(J,KKK)
      ENDIF
      PP=P(IQ1,IQ2,K)
      PPN=PN(IQ1,IQ2,K)
      BPT=PBOT(IQ1,IQ2,K)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PPN,RSON)
      CALL INTPVT(BPT,RSLOPE,POT,RSOT,MPOT,PP,RSO)
      CALL INTERP(PWT,RSWT,MPWT,PPN,RSWN)
      CALL INTERP(PWT,RSWT,MPWT,PP,RSW)
      RSOAV=0.5*(RSO+RSON)
      RSWAV=0.5*(RSW+RSWN)
      FACTOR=PID(J,KKK)*5.615*(PP-PWF(J,KKK))
      IF(KIP(J).EQ.-13) GO TO 2053
      QW(IQ1,IQ2,K)=GMW(J,KKK)/BW(IQ1,IQ2,K)*FACTOR
      IF(KIP(J).EQ.-12)QW(IQ1,IQ2,K)=
     & (GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK))/BW(IQ1,IQ2,K)*FACTOR
      IF(KIP(J).EQ.-12) GO TO 2057
      QO(IQ1,IQ2,K)=GMO(J,KKK)/BO(IQ1,IQ2,K)*FACTOR
      IF(FACTOR.LE.0.)QO(IQ1,IQ2,K)=0.
      QG(IQ1,IQ2,K)=GMG(J,KKK)/BG(IQ1,IQ2,K)*FACTOR
     & +RSOAV*QO(IQ1,IQ2,K)+RSWAV*QW(IQ1,IQ2,K)
      QOSUM=QOSUM+QO(IQ1,IQ2,K)
      QTSUM=QTSUM+QO(IQ1,IQ2,K)+QW(IQ1,IQ2,K)+QG(IQ1,IQ2,K)
      GO TO 2057
 2053 QG(IQ1,IQ2,K)=(GMO(J,KKK)+GMW(J,KKK)+GMG(J,KKK))
     & /BG(IQ1,IQ2,K)*FACTOR
      QGSUM=QGSUM+QG(IQ1,IQ2,K)
 2057 CONTINUE
C**** CHECK RATE CONSTRAINTS ALSO ****************
      IF(QVO(J).GT.0.AND.QOSUM.GT.QVO(J)*5.615)THEN
        QORATIO=QVO(J)*5.615/QOSUM
        DO 1351 KKK=IQ3,LAY
          K=KKK
         IF(IWELTYPE.EQ.2)THEN
          IQ1=IQH1(J,KKK)
          IQ2=IQH2(J,KKK)
          K=IQH3(J,KKK)
         ENDIF
         QO(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*QORATIO
         PWFC(J,KKK)=PP-(PP-PWF(J,KKK))*QORATIO
 1351   CONTINUE
       ENDIF
      IF(QVT(J).GT.0.AND.QTSUM.GT.QVT(J)*5.615)THEN
        QTRATIO=QVT(J)*5.615/QTSUM
        DO 1352 KKK=IQ3,LAY
          K=KKK
         IF(IWELTYPE.EQ.2)THEN
          IQ1=IQH1(J,KKK)
          IQ2=IQH2(J,KKK)
          K=IQH3(J,KKK)
         ENDIF
         QO(IQ1,IQ2,K)=QO(IQ1,IQ2,K)*QTRATIO
         QW(IQ1,IQ2,K)=QW(IQ1,IQ2,K)*QTRATIO
         QG(IQ1,IQ2,K)=QG(IQ1,IQ2,K)*QTRATIO
         PWFC(J,KKK)=PP-(PP-PWF(J,KKK))*QTRATIO
 1352   CONTINUE
       ENDIF
      IF(QVG(J).EQ.0.OR.QGSUM.EQ.0)GO TO 2059
        QGRATIO=QVG(J)*1000/QGSUM
       IF(QGRATIO.LT.1.)THEN
        DO 1353 KKK=IQ3,LAY
          K=KKK
         IF(IWELTYPE.EQ.2)THEN
          IQ1=IQH1(J,KKK)
          IQ2=IQH2(J,KKK)
          K=IQH3(J,KKK)
         ENDIF
         QG(IQ1,IQ2,K)=QG(IQ1,IQ2,K)*QGRATIO
         PWFC(J,KKK)=PP-(PP-PWF(J,KKK))*QGRATIO
 1353   CONTINUE
       ENDIF
 2059 CONTINUE
      RETURN
      END
C
C      SUBROUTINE MATBAL(DELT0,D5615,NX,NY,NZ)
      SUBROUTINE MATBAL(DELT0,D5615,NX,NY,NZ,PN,P,QO,QW,QG,VP,PV)
C
C=================================================================MATBAL
C
C   COMPUTE SIMPLE MATERIAL BALANCE.
C
C=======================================================================
C
      REAL MBEO,MBEW,MBEG,MCFG,MCFG1,MCFGT,MCFGI
      COMMON /NUMBER/ II,JJ,KK
      COMMON /VOL/ SCFO,SCFW,SCFG,SCFG1,STBO,STBW,MCFG,MCFG1,MCFGT,
     & MCFGI,STBOI,STBWI,RESVOL
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      COMMON /MBE/ MBEO,MBEW,MBEG
      DIMENSION PN(NX,NY,NZ)
C     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
C     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      DIMENSION P(NX,NY,NZ)
      DIMENSION QO(NX,NY,NZ),QW(NX,NY,NZ),QG(NX,NY,NZ)
      DIMENSION VP(NX,NY,NZ)
      COMMON /BAL/ OP,WP,GP,WI,GI,PAVG0,PAVG,OPR,WPR,GPR,WIR,GIR,CWOR,
     & WOR,CGOR,GOR
C
      FACT=D5615*DELT0
      PAVG0=0.0
      PAVG=0.0
      OP=0.0
      WP=0.0
      GP=0.0
      WI=0.0
      GI=0.0
      DO 100 K=1,KK
      DO 100 J=1,JJ
      DO 100 I=1,II
      PAVG0=PAVG0+PN(I,J,K)*VP(I,J,K)
      PAVG=PAVG+P(I,J,K)*VP(I,J,K)
      OP=OP+QO(I,J,K)*FACT
      IF(QW(I,J,K).GT.0.0)WP=WP+QW(I,J,K)*FACT
      IF(QW(I,J,K).LT.0.0)WI=WI+QW(I,J,K)*FACT
      IF(QG(I,J,K).GT.0.0)GP=GP+QG(I,J,K)*DELT0
      IF(QG(I,J,K).LT.0.0)GI=GI+QG(I,J,K)*DELT0
100   CONTINUE
      COP=COP+OP
      CWP=CWP+WP
      CGP=CGP+GP*.001
      CWI=CWI+WI
      CGI=CGI+GI*.001
C****CONVERT SCF TO MCF.
      GP=GP*.001
      GI=GI*.001
      DIV=1.0/DELT0
      OPR=OP*DIV
      WPR=WP*DIV
      GPR=GP*DIV
      WIR=WI*DIV
      GIR=GI*DIV
      PAVG=PAVG/RESVOL
      PAVG0=PAVG0/RESVOL
C
      DENOM1=STBOI - OP
      IF(DENOM1.LT.1.E-7) GO TO 200
      MBEO=(STBO/(STBOI-OP) - 1.0)*100.0
 200  CONTINUE
      DENOM2=STBWI-WP-WI
      IF(DENOM2.LT.1.E-7) GO TO 201
      MBEW=(STBW/(STBWI-WP-WI) - 1.0)*100.0
 201  CONTINUE
      DENOM3=MCFGI-GP-GI
      IF(DENOM3.LT.1.E-7) GO TO 203
      MBEG=(MCFGT/(MCFGI-GP-GI) - 1.0)*100.0
 203  CONTINUE
C
      IF(OP.EQ.0.0)GOR=0.0
      IF(OP.EQ.0.0)WOR=0.0
      IF(OP.EQ.0.0)GO TO 333
      GOR=GP/OP
      WOR=WP/OP
333   IF(COP.EQ.0.D0)CGOR=0.0
      IF(COP.EQ.0.D0)CWOR=0.0
      IF(COP.EQ.0.D0)GO TO 666
      CGOR=CGP/COP*1000.0D0
      CWOR=CWP/COP
666   CONTINUE
C
      GP=GP*0.001
      GI=GI*0.001
      RETURN
      END
C      SUBROUTINE PRTPS(N,DELT0,NX,NY,NZ)
      SUBROUTINE PRTPS(N,DELT0,NX,NY,NZ,PBOT,PN,P,SON,SWN,SGN,
     &    SO,SW,SG)
C
C==================================================================PRTPS
C
C   WRITE SUMMARY REPORT.
C
C======================================================================
C
C      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
C	NEW DP LOCAL VARIABLES **RJH**
C       MODIFICATION FOR DOUBLE PRECISION READ OF SO & SW IN RESTART  **RJH**

      REAL MBEO,MBEW,MBEG,MCFG,MCFG1,MCFGT,MCFGI
      COMMON /IO/ NI,NO
      COMMON /ADD1/ IM,JM,KM,ETI,FT,FTMAX,DELT
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      COMMON /PSCNTL/ KPI,KSI
      COMMON /MBE/ MBEO,MBEW,MBEG
      COMMON /NUMBER/ II,JJ,KK
      COMMON /VOL/ SCFO,SCFW,SCFG,SCFG1,STBO,STBW,MCFG,MCFG1,MCFGT,
     & MCFGI,STBOI,STBWI,RESVOL
      COMMON /BAL/ OP,WP,GP,WI,GI,PAVG0,PAVG,OPR,WPR,GPR,WIR,GIR,CWOR,
     & WOR,CGOR,GOR
      COMMON /PRTCNT/ IWLCNG,ICHANG,IWLREP,ISUMRY,
     & IPMAP,ISOMAP,ISWMAP,ISGMAP,IPBMAP
      COMMON /BUBBLE/ PBO,VSLOPE,BSLOPE,RSLOPE,PMAXT,IREPRS,
     & RHOSCO,RHOSCG,RHOSCW,MSAT,MPOT,MPWT,MPGT
      DIMENSION PBOT(NX,NY,NZ)
      DIMENSION PN(NX,NY,NZ),
     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ)
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
C     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
      DIMENSION P(NX,NY,NZ)
      DIMENSION SO(NX,NY,NZ),SW(NX,NY,NZ),SG(NX,NY,NZ)
C      REAL*8 DPP(NX,NY,NZ),DPSO(NX,NY,NZ),DPSW(NX,NY,NZ)

      PPM=0.
      SOM=0.
      SWM=0.
      SGM=0.
      IF(N.EQ.1) GO TO 300
      DO 240 K=1,KK
      DO 240 J=1,JJ
      DO 240 I=1,II
      DPO=P(I,J,K)-PN(I,J,K)
      DSO=SO(I,J,K)-SON(I,J,K)
      DSW=SW(I,J,K)-SWN(I,J,K)
      DSG=SG(I,J,K)-SGN(I,J,K)
      IF(ABS(DPO).LE.ABS(PPM))GO TO 210
      PPM=DPO
      IPM=I
      JPM=J
      KPM=K
  210 IF(ABS(DSO).LE.ABS(SOM))GO TO 220
      SOM=DSO
      IOM=I
      JOM=J
      KOM=K
  220 IF(ABS(DSW).LE.ABS(SWM))GO TO 230
      SWM=DSW
      IWM=I
      JWM=J
      KWM=K
  230 IF(ABS(DSG).LE.ABS(SGM))GO TO 240
      SGM=DSG
      IGM=I
      JGM=J
      KGM=K
  240 CONTINUE
      WRITE(NO,5)
      WRITE(NO,105)
      WRITE(NO,6)
      NLM=N-1
      GORM=1000.*GOR
      WRITE(NO,110)ETI,NLM,DELT0,PAVG,PAVG0,IPM,JPM,KPM,PPM,
     & IOM,JOM,KOM,SOM,IGM,JGM,KGM,SGM,IWM,JWM,KWM,SWM,
     & MBEO,MBEG,MBEW,OPR,COP,GPR,CGP,WPR,CWP
      WRITE(NO,111) GIR,CGI,WIR,CWI,WOR,CWOR,GORM,CGOR
  110 FORMAT(/,1X,'ELAPSED TIME (DAYS)        =',F9.2,3X,
     & 'TIME STEP NUMBER            =',I5,7X,
     & 'TIME STEP SIZE (DAYS)       =',F9.2,//,1X,
     & 'CURRENT AVG RES PRESSURE    =',F9.1,3X,
     & 'PREVIOUS AVG RES PRESSURE    =',F9.1,3X,
     & 'PRESSURE DPMAX(',I3,',',I3,',',I3,')  =',F9.1,/,1X,
     & 'OIL DSMAX(',I3,',',I3,',',I3,')      =',F9.5,3X,
     & 'GAS DSMAX(',I3,',',I3,',',I3,')      =',F9.5,3X,
     & 'WATER DSMAX(',I3,',',I3,',',I3,')     =',F9.5,/,1X,
     & 'OIL MATERIAL BALANCE (%)    =',F9.6,3X,
     & 'GAS MATERIAL BALANCE (%)     =',F9.6,3X,
     & 'WATER MATERIAL BALANCE (%)   =',F9.6,//,1X,
     & 'OIL PRODUCTION RATE (STB/D) =',F9.1,3X,
     & 'CUM. OIL PRODUCTION (STB)    =',E10.4,/,1X,
     & 'GAS PRODUCTION RATE (MSCF/D)=',F9.1,3X,
     & 'CUM. GAS PRODUCTION (MSCF)   =',E10.4,/,1X,
     & 'WATER PRODUCTION RATE(STB/D)=',F9.1,3X,
     & 'CUM. WATER PRODUCTION (STB)  =',E10.4)
  111  FORMAT(/1X,'GAS INJECTION RATE (MSCF/D) =',F9.1,3X,
     & 'CUM. GAS INJECTION (MSCF)    =',E10.4,/,1X,
     & 'WATER INJECTION RATE (STB/D)=',F9.1,3X,
     & 'CUM. WATER INJECTION (STB)   =',E10.4,//,1X,
     & 'PRODUCING WOR (STB/STB)     =',F9.3,3X,
     & 'CUM. WOR (STB/STB)          =',F9.3,/,1X,
     & 'PRODUCING GOR (SCF/STB)     =',F9.0,3X,
     & 'CUM. GOR (SCF/STB)          =',F9.0,//)
  300 IF(N.EQ.1.AND.KPI.EQ.0) WRITE(NO,304)
  304 FORMAT(//T15,7('*'),' INITIAL ARRAYS ',7('*')//)
      IF(IPMAP.EQ.0.AND.N.NE.1)GO TO 315
      WRITE(NO,61)
      DO 310 K=1,KK
      WRITE(NO,51)K
      DO 310 J=1,JJ
      WRITE(NO,41)(P(I,J,K),I=1,II)
310   CONTINUE
C
  315 IF(ISOMAP.EQ.0.AND.N.NE.1)GO TO 422
      WRITE(NO,71)
      DO 420 K=1,KK
      WRITE(NO,51)K
      DO 420 J=1,JJ
      WRITE(NO,101)(SO(I,J,K),I=1,II)
420   CONTINUE
C
  422 IF(ISWMAP.EQ.0.AND.N.NE.1)GO TO 432
      WRITE(NO,81)
      DO 430 K=1,KK
      WRITE(NO,51)K
      DO 430 J=1,JJ
      WRITE(NO,101)(SW(I,J,K),I=1,II)
430   CONTINUE
C
  432 IF(ISGMAP.EQ.0.AND.N.NE.1)GO TO 4400
      WRITE(NO,91)
      DO 440 K=1,KK
      WRITE(NO,51)K
      DO 440 J=1,JJ
      WRITE(NO,101)(SG(I,J,K),I=1,II)
  440 CONTINUE
4400  CONTINUE
      IF(IPBMAP.EQ.0.AND.N.NE.1)GO TO 4500
      WRITE(NO,102)
  102 FORMAT(//T15,'***** BUBBLE POINT PRESSURE DISTRIBUTION *****'/)
      DO 450 K=1,KK
      WRITE(NO,51)K
      DO 450 J=1,JJ
      WRITE(NO,41)(PBOT(I,J,K),I=1,II)
  450 CONTINUE
4500  CONTINUE
      IF(N.NE.1)WRITE(NO,7)
C
5     FORMAT(//T30,69('*'),/T30,'*',T98,'*'/T30,'*',T98,'*')
6     FORMAT(T30,'*',T98,'*'/T30,'*',T98,'*'/T30,69('*'),///)
7     FORMAT(///T5,54('*'),'  END OF REPORT  ',54('*'),6(/))
105   FORMAT(T30,'*',T45,'SUMMARY REPORT: BOAST - VHS     '
     &,T98,'*')
41    FORMAT(5X,20F6.0)
51    FORMAT(/1X,'K =',I2/)
61    FORMAT(//T15,'***** RESERVOIR PRESSURE DISTRIBUTION *****'/)
71    FORMAT(//T15,'*********  OIL SATURATION  *********'/)
81    FORMAT(//T15,'******** WATER SATURATION *********'/)
91    FORMAT(//T15,'*********  GAS SATURATION  *********'/)
101   FORMAT(5X,20F6.3)
      RETURN
      END
      
C      SUBROUTINE LSOR(N,NX,NY,NZ,AZL,BZL,CZL,DZL,UZL,UM)
      SUBROUTINE LSOR(N,NX,NY,NZ,AW,AE,AN,AS,AB,AT,E,B,P,
     &  AZL,BZL,CZL,DZL,UZL,UM)
C
C===================================================================LSOR
C
C     THIS PROGRAM SOLVES A LINEAR SYSTEM OF THREE-DIMENSIONAL FINITE-
C     DIFFERENCE EQUATIONS.
C     ( LINEAR SUCCESSIVE OVERRELAXATION METHOD )
C
C======================================================================
C
C      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      COMMON /CODE/ KSM1,KSN1,KCO1,NN,FACT1,FACT2,TMAX,KSOL,MITER,
     & OMEGA,TOL,TOL1,KSN,KSM,KCO,KTR,KCOFF,DSMAX,DPMAX,WORMAX,
     & GORMAX,PAMIN,PAMAX
      DIMENSION AW(NX,NY,NZ),AE(NX,NY,NZ),AN(NX,NY,NZ),
     & AS(NX,NY,NZ),AB(NX,NY,NZ),AT(NX,NY,NZ),E(NX,NY,NZ),B(NX,NY,NZ)
      DIMENSION P(NX,NY,NZ)
      
      DIMENSION UM(NX),AZL(NX),BZL(NX),CZL(NX),DZL(NX),UZL(NX)
C
C      DIV=DELT/DELT0
      NITER=0
      DMAX=1.0
      RHO1=0.0
      THETA=0.0
C
   11 CONTINUE
      TW = 1.0 - OMEGA
      DMAX0=DMAX
      THETA0=THETA
      IF(NITER.GE.MITER)WRITE(NO,30)NITER,TOL,DMAX
      IF(NITER.GE.MITER)RETURN
      NITER = NITER +1
      DMAX=0.0
C
   2  DO 20 K=1,KK
      DO 20 J=1,JJ
      DO 15 I=1,II
      UM(I) = P(I,J,K)
C
27    AZL(I) = AW(I,J,K)
      BZL(I)=E(I,J,K)
      CZL(I)=AE(I,J,K)
      DZL(I)=B(I,J,K)
      IF(JJ.EQ.1) GO TO 14
      JM=J-1
      JP=J+1
      IF(J.EQ.1) JM=1
      IF(J.EQ.JJ) JP=JJ
      DZL(I)=DZL(I) - AS(I,J,K)*P(I,JM,K) - AN(I,J,K)*P(I,JP,K)
   14 IF(KK.EQ.1) GO TO 15
      KM=K-1
      KP=K+1
      IF(K.EQ.1)KM=1
      IF(K.EQ.KK) KP=KK
      DZL(I)=DZL(I) - AT(I,J,K)*P(I,J,KM) - AB(I,J,K)*P(I,J,KP)
15    CONTINUE
C
      CALL LTRI(II,AZL,BZL,CZL,DZL,UZL)
C
      DO 16 I=1,II
      GSLSOR=UZL(I)
      P(I,J,K) = TW*UM(I) + OMEGA*GSLSOR
      ARG=P(I,J,K) - UM(I)
      DM=ABS(ARG)
      IF(DM.GT.DMAX) DMAX = DM
16    CONTINUE
C
20    CONTINUE
   3  IF(TOL1.EQ.0.0)GO TO 25
      THETA=DMAX/DMAX0
      DELTA=THETA-THETA0
      ARG=DELTA
      ARG=ABS(ARG)
      IF(ARG.GT.TOL1)GO TO 25
      OM=OMEGA-1.0
      RHO1=(THETA+OM)*(THETA+OM)/(THETA*OMEGA*OMEGA)
      IF(RHO1.GE.1.0)GO TO 25
      ARG=1.0-RHO1
      OMEGA=2.0/(1.0+SQRT(ARG))
C
25    IF(DMAX.GT.TOL)GO TO 11
      IF(N.NE.KSN) GO TO 300
      WRITE(NO,40)NITER,OMEGA,DMAX,THETA,RHO1
 300  CONTINUE
40    FORMAT(T5,'CONVERGENCE(LSOR)  HAS BEEN REACHED AFTER ',I3,
     &' ITERATIONS',5X,'OMEGA = ',F6.3/
     &T5,'DMAX = ',F10.6,5X,'THETA = ',F10.6,5X,'RHO1 = ',F10.6/)
30    FORMAT(T15,'CONVERGENCE(LSOR)  WAS NOT REACHED IN ',I5,
     &' ITERATIONS'/T15,'TOL = ',F10.7,10X,'DMAX = ',F15.7)
C
      RETURN
      END
      SUBROUTINE LTRI(N,A,BI,C,D,X)
C
C===================================================================LTRI
C
C     THIS PR0GRAM SOLVES THE TRIDIAGONAL SYSTEM
C     GENERATED BY THE SYSTEM OF N EQUATIONS
C
C     A(I)*U(I-1) + BI(I)*U(I) + C(I)*U(I+1) = D(I)
C
C=======================================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'

      DIMENSION A(N),BI(N),C(N),D(N),
     & X(N),BETA(NX1),GAMMA(NX1),W(NX1)

      BETA(1)=BI(1)
      GAMMA(1)=D(1)/BI(1)
      NM=N-1
C      COMPUTE FORWARD SOLUTION.
      DO 10 I=1,NM
      W(I)=C(I)/BETA(I)
      IP=I+1
   10 BETA(IP)=BI(IP)-A(IP)*W(I)
      DO 20 I=2,N
      IM=I-1
   20 GAMMA(I)=(D(I)-A(I)*GAMMA(IM))/BETA(I)
C      COMPUTE BACK SOLUTION
      X(N)=GAMMA(N)
      DO 30 J=1,NM
      I=N-J
      IP=I+1
   30 X(I)=GAMMA(I)-W(I)*X(IP)
      RETURN
      END
C
C
C
      SUBROUTINE GRINT
C
C======================================================GRINT
C
C   READS DATA FOR GRAPHICS INITIALIZATION PURPOSES.
C
C===========================================================
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      CHARACTER*2 TITLE(80),FNAME,MNAME,TAG
      CHARACTER*4 LNAME
      COMMON /CHRCTR/ FNAME,MNAME(NUM),LNAME,TAG
      COMMON /TTEST/ NUMPRD,SONTVL
      COMMON /IO/ NI,NO
      COMMON /COUNT/ N1READ,N2READ
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
  69  FORMAT(40A2)
      READ(NI,70) TAG
      N1READ=N1READ+1
  70  FORMAT(A2)
      READ(NI,*) SONTVL
      N1READ=N1READ+1
      RETURN
      END
C
C
C      SUBROUTINE SHALE(NX,NY,NZ)
      SUBROUTINE SHALE(NX,NY,NZ,TZ)
C
C========================================================SHALE
C
C    ALLOWS BROAD ZEROING OF TRANSMISSIBILITIES FOR THE
C    PURPOSE OF SIMULATING SHALES.
C
C=============================================================
C
C      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      CHARACTER*2 TITLE(80)
      COMMON /IO/ NI,NO
      COMMON /NUMBER/ II,JJ,KK
      COMMON /COUNT/ N1READ,N2READ
C      DIMENSION TX(NX+1,NY,NZ),TY(NX,NY+1,NZ)
      DIMENSION TZ(NX,NY,NZ+1)
      READ(NI,69) (TITLE(IH),IH=1,40)
      N1READ=N1READ+1
  69  FORMAT(40A2)
      READ(NI,*) NSHOPS
      N1READ=N1READ+1
      READ(NI,*) NTLAYR
      N1READ=N1READ+1
      IF(NSHOPS.EQ.1) THEN
      DO 1 LZ=1,NTLAYR
      READ(NI,*) KT
      N1READ=N1READ+1
      DO 1 JLZ=1,JJ
      DO 1 ILZ=1,II
      TZ(ILZ,JLZ,KT)=0.0
  1   CONTINUE
      ELSE
      DO 2 LZ=1,NTLAYR
      READ(NI,*) KT,NUTM
      N1READ=N1READ+1
      DO 2 LZF=1,NUTM
      READ(NI,*) IGBL
      N1READ=N1READ+1
      DO 2 JLZ=1,JJ
      TZ(IGBL,JLZ,KT)=0.0
  2   CONTINUE
      ENDIF
      RETURN
      END
C
C      SUBROUTINE SLANT(NVQN,NVQNH,NVQNS,NX,NY,NZ)
      SUBROUTINE SLANT(NVQN,NVQNH,NVQNS,NX,NY,NZ,EL,PN,KX,KY,KZ,
     &   DX,DY,DZ )
C
C=================================================================SLANT
C
C   READ WELL INFORMATION AND CALCULATE PRODUCTIVITY INDEX FOR EACH
C     WELL BLOCK
C
C   (NVQN, NVQNH, NVQNS ARE # OF VERTICAL, HORIZON., & SLANT WELL)
C=======================================================================
C
C   TWO SLANT WELL OPTIONS ARE AVAILABLE :
C        OPTION 1 -- X,Y,Z COORDINATES OF STARTING WELLBLOCK, 
C                    WELL LENGTH, VERTICAL & HORIZONTAL ANGLES
C        OPTION 2 -- X,Y,Z COORDINATES OF STARTING WELLBLOCK
C                    AND X,Y,Z COORDINATES OF ENDING WELLBLOCK
C
      INCLUDE 'D:\F77L\BOASTVHS\VHSPARM.FOR'
      REAL KXN,KYN,KZN
      CHARACTER*5 WELLID
      COMMON /IO/ NI,NO
      DIMENSION ELS(NZ1),PNS(NZ1)
      COMMON /IDEN/ WELLID(NW)
      COMMON /SRATE/ PID(NW,NL),PWF(NW,NL),PWFC(NW,NL),KIP(NW),
     & GMO(NW,NL),GMW(NW,NL),GMG(NW,NL),LAYER(NW),QVO(NW),
     & QVW(NW),QVG(NW),QVT(NW),CUMO(NW,NL),CUMW(NW,NL),CUMG(NW,NL)
      COMMON /IQN/ IQN1(NW),IQN2(NW),IQN3(NW)
      COMMON /IQH/ IQH1(NW,NL),IQH2(NW,NL),IQH3(NW,NL),COND(NW)
      REAL KX(NX,NY,NZ),KY(NX,NY,NZ),KZ(NX,NY,NZ)
      DIMENSION DX(NX,NY,NZ),DY(NX,NY,NZ),DZ(NX,NY,NZ)
      COMMON /NUMBER/ II,JJ,KK
      COMMON /R1XYZ/ R11X,R21X,R31X,R41X,R11Y,R21Y,R31Y,R41Y,
     & R11Z,R21Z,R31Z,R41Z
      COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
      COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
      COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K
      DIMENSION EL(NX,NY,NZ)
      DIMENSION PN(NX,NY,NZ)
C     & SON(NX,NY,NZ),SWN(NX,NY,NZ),SGN(NX,NY,NZ),
C     & A1(NX,NY,NZ),A2(NX,NY,NZ),A3(NX,NY,NZ),
C     & SUM(NX,NY,NZ),GAM(NX,NY,NZ)
C
       WRITE(6,607)
       WRITE(6,608)
       N1=NVQN+NVQNH+1
       N2=NVQN+NVQNH+NVQNS
C
       DO 1000 N=N1,N2
       M=0
       READ(NI,10) WELLID(N)
 10    FORMAT(A5)
       READ(NI,*)IFLAG,KIP(N),QVO(N),QVW(N),QVG(N),QVT(N),COND(N)
       IF(IFLAG.EQ.1)THEN
         READ(NI,*)IS,JS,KS,WELENGTH,THETA,ALPHA,IS1,JS1,KS1
        ELSE
         READ(NI,*)IS,JS,KS,IE,JE,KE,IS1,JS1,KS1
        ENDIF
       READ(NI,*)RW,S,PWF(N,1),JFLAG
       IF(IFLAG.EQ.1)GO TO 108
C ***** CALCULATE WELENGTH, THETA, ALPHA FOR SLANT WELL OPTION 2
        ISTEP=1
        JSTEP=1
        KSTEP=1
        IF(IS.GT.IE)ISTEP=-1
        IF(JS.GT.JE)JSTEP=-1
        IF(KS.GT.KE)KSTEP=-1
        DXWEL=0.
        DYWEL=0.
        DZWEL=0.
        DO 102 I=IS,IE,ISTEP
 102    DXWEL=DXWEL+DX(I,JS,KS)
        DO 104 J=JS,JE,JSTEP
 104    DYWEL=DYWEL+DY(IS,J,KS)
        DO 106 K=KS,KE,KSTEP
 106    DZWEL=DZWEL+DZ(IS,JS,K)
        WELSUM=0.
        IF(IS.NE.IE)WELSUM=WELSUM+DXWEL**2
        IF(JS.NE.JE)WELSUM=WELSUM+DYWEL**2
        IF(KS.NE.KE)WELSUM=WELSUM+DZWEL**2
        WELENGTH=WELSUM**.5
        IF(WELENGTH.EQ.0.)THEN
         WELENGTH=(DXWEL**2+DYWEL**2+DZWEL**2)**.5
        ENDIF
        THETA=ACOS(KSTEP*DZWEL/WELENGTH)*180./3.1416
        ALPHA=ATAN(-ISTEP*JSTEP*DYWEL/DXWEL)*180./3.1416 
        IF(ALPHA.LT.0.)ALPHA=ALPHA+180.
        IF(JSTEP.EQ.1)ALPHA=ALPHA+180.
C
       IF(IS.EQ.IE.AND.JS.EQ.JE.AND.KS.LT.KE)THETA=0.
       IF(IS.EQ.IE.AND.JS.EQ.JE.AND.KS.GT.KE)THETA=180.
       IF(KS.EQ.KE)THETA=90.
       IF(JS.EQ.JE.AND.IS.LE.IE)ALPHA=0.
       IF(JS.EQ.JE.AND.IS.GT.IE)ALPHA=180.
       IF(IS.EQ.IE.AND.JS.LT.JE)ALPHA=270.
       IF(IS.EQ.IE.AND.JS.GE.JE)ALPHA=90.
C ***** DETERMINE THE WELL DIRECTION IN (I,J,K) COORDINATE
 108   IF(THETA.LT.90.)KD=1
       IF(THETA.EQ.90.)KD=0
       IF(THETA.GT.90.)KD=-1
       IF(ALPHA.NE.0.) GO TO 110
       ID=1
       JD=0
 110   IF(ALPHA.NE.90)GO TO 120
       ID=0
       JD=-1
 120   IF(ALPHA.NE.180)GO TO 130
       ID=-1
       JD=0
 130   IF(ALPHA.NE.270)GO TO 140
       ID=0
       JD=1
 140   IF(ALPHA.LE.0.OR.ALPHA.GE.90)GO TO 150
       ID=1
       JD=-1
 150   IF(ALPHA.LE.90.OR.ALPHA.GE.180)GO TO 160
       ID=-1
       JD=-1
 160   IF(ALPHA.LE.180.OR.ALPHA.GE.270)GO TO 170
       ID=-1
       JD=1
 170   IF(ALPHA.LE.270.OR.ALPHA.GE.360)GO TO 180
       ID=1
       JD=1
C  **** DEFINE THE STARTING WELL COORDINATE
 180   RL=WELENGTH
       DXI=DX(IS,JS,KS)
       DYI=DY(IS,JS,KS)
       DZI=DZ(IS,JS,KS)
       TG=THETA*3.1416/180.
       AG=ALPHA*3.1416/180.
       I=IS
       J=JS
       K=KS
C ***** 
       DP=PN(I,J,K)-PWF(N,1)
       DO 190 NK=1,KK
        ELS(NK)=EL(I,J,NK)
        PNS(NK)=PN(I,J,NK)
 190   CONTINUE  
C ***** LOOP TO CALCULATE PID VALUE IN Mth WELL BLOCK ********
 210   M=M+1
       LAYER(N)=M
       KXN=1
       KYN=1
       KZN=1
C       KYN=KY(I,J,K)/KX(I,J,K)
C       KZN=KZ(I,J,K)/KX(I,J,K)
       DX0=DX(I,J,K)
       DY0=DY(I,J,K)
       DZ0=DZ(I,J,K)
       RLX=0.
       RLY=0.
       RLZ=0.
       IF(ID.NE.0.AND.TG.NE.0.) RLX=ABS(DXI/(SIN(TG)*COS(AG)))
       IF(JD.NE.0.AND.TG.NE.0..AND.AG.NE.0.)
     &  RLY=ABS(DYI/(SIN(TG)*SIN(AG)))
       IF(KD.NE.0) RLZ=ABS(DZI/COS(TG))
C  **** CALCULATE THE X,Y,Z COMPONENT OF WELL LENGTH WITHIN WELL BLOCK
       RMIN=99999.
       IF(RLX.NE.0.AND.ABS(RLX).LT.RMIN) RMIN=ABS(RLX)
       IF(RLY.NE.0.AND.ABS(RLY).LT.RMIN) RMIN=ABS(RLY)
       IF(RLZ.NE.0.AND.ABS(RLZ).LT.RMIN) RMIN=ABS(RLZ)
       DXX=ABS(RMIN*SIN(TG)*COS(AG))
       DYY=ABS(RMIN*SIN(TG)*SIN(AG))
       DZZ=ABS(RMIN*COS(TG))
       DL=(DXX*DXX+DYY*DYY+DZZ*DZZ)**.5
       IF(RL.GT.DL)GO TO 250
       DXX=DXX*RL/DL
       DYY=DYY*RL/DL
       DZZ=DZZ*RL/DL
       DL=RL
 250   IS2=0
       JS2=0
       KS2=0
       IF(ID.EQ.0)IS2=1
       IF(ID.EQ.1.AND.I.GE.IS1)IS2=1
       IF(ID.EQ.-1.AND.I.LE.IS1)IS2=1
       IF(JD.EQ.0)JS2=1
       IF(JD.EQ.1.AND.J.GE.JS1)JS2=1
       IF(JD.EQ.-1.AND.J.LE.JS1)JS2=1
       IF(KD.EQ.0)KS2=1
       IF(KD.EQ.1.AND.K.GE.KS1)KS2=1
       IF(KD.EQ.-1.AND.K.LE.KS1)KS2=1
       IF((IS2*JS2*KS2).EQ.1)GO TO 300
C 250   IS2=(I-IE)*(I-IS1)
C       JS2=(J-JE)*(J-JS1)
C       KS2=(K-KE)*(K-KS1)
C       IF(IS2.LE.0.AND.JS2.LE.0.AND.KS2.LE.0)GO TO 300
       M=M-1
       LAYER(N)=M
       GO TO 890
C ***** CONSTANT WELLBORE PRESSURE IF JFLAG = 1, OR 
C       CLACULATE WELLBORE PRESSURE FROM PWF(N,1) AND HYDROLIC GRADIENT 
 300    IF(KIP(N).GT.0.OR.JFLAG.EQ.1)THEN
         PWF(N,M)=PWF(N,1)
        ELSE
         ELSW=EL(I,J,K)+KD*(DZ(I,J,K)/2+DZZ/2-DZI)
         CALL INTERP(ELS,PNS,KK,ELSW,PWFS)
         PWF(N,M)=PWFS-DP
        ENDIF
C
C ***** PID CALCULATION STARTS --- BASED ON ABOU-KASSEM'S PAPER 
C       WHICH CALCULATES PID IN ONE DIRECTION (SPEJ 8/85) ......
       XW=ID*(DX(I,J,K)/2-(DXI-DXX/2))
       YW=JD*(DY(I,J,K)/2-(DYI-DYY/2))
       ZW=KD*(DZ(I,J,K)/2-(DZI-DZZ/2))
C ****  COUNT THE NUMBER OF BOUNDARY THIS WELL BLOCK HAS
       NBOUND=0
       IF(I.EQ.1)THEN
        T1Y=0.
        T4Z=0.
        NBOUND=NBOUND+1
       ELSE
        DXM=(DX(I,J,K)+DX(I-1,J,K))/2
        T1Y=DZ0/DXM*KXN
        T4Z=DY0/DXM*KXN
C        T1Y=TX(I,J,K)
C        T4Z=TX(I,J,K)
       ENDIF
       IF(J.EQ.1)THEN
        NBOUND=NBOUND+1
        T4X=0.
        T1Z=0.
       ELSE
        DYM=(DY(I,J,K)+DY(I,J-1,K))/2
        T4X=DZ0/DYM*KYN
        T1Z=DX0/DYM*KYN
C        T4X=TY(I,J,K)
C        T1Z=TY(I,J,K)
       ENDIF
       IF(K.EQ.1)THEN
        NBOUND=NBOUND+1
        T1X=0.
        T4Y=0.
       ELSE
        DZM=(DZ(I,J,K)+DZ(I,J,K-1))/2
        T1X=DY0/DZM*KZN
        T4Y=DX0/DZM*KZN
C        T1X=TZ(I,J,K)
C        T4Y=TZ(I,J,K)
       ENDIF
       IF(I.EQ.II)THEN
        NBOUND=NBOUND+1
        T3Y=0.
        T2Z=0.
       ELSE
        DXP=(DX(I,J,K)+DX(I+1,J,K))/2
        T3Y=DZ0/DXP*KXN
        T2Z=DY0/DXP*KXN
C        T3Y=TX(I+1,J,K)
C        T2Z=TX(I+1,J,K)
       ENDIF
       IF(J.EQ.JJ)THEN
        NBOUND=NBOUND+1
        T2X=0.
        T3Z=0.
       ELSE
        DYP=(DY(I,J,K)+DY(I,J+1,K))/2
        T2X=DZ0/DYP*KYN
        T3Z=DX0/DYP*KYN
C        T2X=TY(I,J+1,K)
C        T3Z=TY(I,J+1,K)
       ENDIF
       IF(K.EQ.KK)THEN
        NBOUND=NBOUND+1
        T3X=0.
        T2Y=0.
       ELSE
        DZP=(DZ(I,J,K)+DZ(I,J,K+1))/2
        T3X=DY0/DZP*KZN
        T2Y=DX0/DZP*KZN
C        T3X=TZ(I,J,K+1)
C        T2Y=TZ(I,J,K+1)
       ENDIF
C
       BBX=1/(T1X+T2X+T3X+T4X)
       BBY=1/(T1Y+T2Y+T3Y+T4Y)
       BBZ=1/(T1Z+T2Z+T3Z+T4Z)
C
       R11X=YW**2 + (DZM+ZW)**2
       R21X=(DYP-YW)**2 + ZW**2
       R31X=YW**2 + (DZP-ZW)**2
       R41X=(DYM+YW)**2 + ZW**2
       R11Y=ZW**2 + (DXM+XW)**2
       R21Y=(DZP-ZW)**2 + XW**2
       R31Y=ZW**2 + (DXP-XW)**2
       R41Y=(DZM+ZW)**2 + XW**2
       R11Z=XW**2 + (DYM+YW)**2
       R21Z=(DXP-XW)**2 + YW**2
       R31Z=XW**2 + (DYP-YW)**2
       R41Z=(DXM+XW)**2 + YW**2
C  **** RX, RY, RZ FOR NON-BOUNDARY WELL BLOCK
       IF(R11X.EQ.0..AND.T1X.EQ.0.)R11X=1.
       IF(R21X.EQ.0..AND.T2X.EQ.0.)R21X=1.
       IF(R31X.EQ.0..AND.T3X.EQ.0.)R31X=1.
       IF(R41X.EQ.0..AND.T4X.EQ.0.)R41X=1.
       IF(R11Y.EQ.0..AND.T1Y.EQ.0.)R11Y=1.
       IF(R21Y.EQ.0..AND.T2Y.EQ.0.)R21Y=1.
       IF(R31Y.EQ.0..AND.T3Y.EQ.0.)R31Y=1.
       IF(R41Y.EQ.0..AND.T4Y.EQ.0.)R41Y=1.
       IF(R11Z.EQ.0..AND.T1Z.EQ.0.)R11Z=1.
       IF(R21Z.EQ.0..AND.T2Z.EQ.0.)R21Z=1.
       IF(R31Z.EQ.0..AND.T3Z.EQ.0.)R31Z=1.
       IF(R41Z.EQ.0..AND.T4Z.EQ.0.)R41Z=1.
       RX=T1X*ALOG(R11X)+T2X*ALOG(R21X)+T3X*ALOG(R31X)+T4X*ALOG(R41X)
       RY=T1Y*ALOG(R11Y)+T2Y*ALOG(R21Y)+T3Y*ALOG(R31Y)+T4Y*ALOG(R41Y)
       RZ=T1Z*ALOG(R11Z)+T2Z*ALOG(R21Z)+T3Z*ALOG(R31Z)+T4Z*ALOG(R41Z)
       AX=1.
       AY=1.
       AZ=1.
C ***** IF THIS WELL BLOCK HAS NO BOUNDARY ....
       IF(NBOUND.EQ.0)GO TO 800
C ***** IF THIS WELL BLOCK HAS ONE BOUNDARY ....
 400   IF(NBOUND.GT.1)GO TO 500
       IF(I.NE.1.AND.I.NE.II)THEN
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        ENDIF
      IF(J.NE.1.AND.J.NE.JJ)THEN
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        ENDIF
       IF(K.NE.1.AND.K.NE.KK)THEN
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
       GO TO 800
 500   IF(NBOUND.EQ.2)CALL TWOBOUND
       IF(NBOUND.EQ.3)CALL THREEBOUND
       IF(NBOUND.EQ.4)CALL FOURBOUND
C
C ***** PID CALCULATION IN X, Y, AND Z COMPONENT
C     IF "NOT KX=KY=KZ", WELLBORE EQUIVALENT RADIUS IS MODIFIED BASED ON --
C       --- PEACEMAN'S FORMULA (SPEJ, JUNE 1983)
 800   RKXY=KX(I,J,K)/KY(I,J,K)
       RKYZ=KY(I,J,K)/KZ(I,J,K)
       RKXZ=KX(I,J,K)/KZ(I,J,K)
       ROX=EXP(-ALOG(AX)-2*3.1416*BBX+BBX*RX/2)
       IF(RKYZ.NE.1)THEN
        RX1=2*(DY0**2/RKYZ**.5+DZ0**2*RKYZ**.5)**.5
        RX2=((DY0**2+DZ0**2)**.5)*(1/(RKYZ**.25)+RKYZ**.25)
        ROX=ROX*RX1/RX2
       ENDIF
       PIDX=.00708*(KY(I,J,K)*KZ(I,J,K))**.5*DXX/(ALOG(ROX/RW)+S)
C
       ROY=EXP(-ALOG(AY)-2*3.1416*BBY+BBY*RY/2)
       IF(RKXZ.NE.1)THEN
        RY1=2*(DX0**2/RKXZ**.5+DZ0**2*RKXZ**.5)**.5
        RY2=((DX0**2+DZ0**2)**.5)*(1/(RKXZ**.25)+RKXZ**.25)
        ROY=ROY*RY1/RY2
       ENDIF
       PIDY=.00708*(KX(I,J,K)*KZ(I,J,K))**.5*DYY/(ALOG(ROY/RW)+S)
C
       ROZ=EXP(-ALOG(AZ)-2*3.1416*BBZ+BBZ*RZ/2)
       IF(RKXY.NE.1)THEN
        RZ1=2*(DX0**2/RKXY**.5+DY0**2*RKXY**.5)**.5
        RZ2=((DX0**2+DY0**2)**.5)*(1/(RKXY**.25)+RKXY**.25)
        ROZ=ROZ*RZ1/RZ2
       ENDIF
       PIDZ=.00708*(KX(I,J,K)*KY(I,J,K))**.5*DZZ/(ALOG(ROZ/RW)+S)
C
       PID(N,M)=(PIDX**2+PIDY**2+PIDZ**2)**.5
       IQH1(N,M)=I
       IQH2(N,M)=J
       IQH3(N,M)=K
      WRITE(6,703)IQH1(N,M),IQH2(N,M),IQH3(N,M),WELLID(N),DL,QVO(N),
     & QVW(N),QVG(N),QVT(N),PWF(N,M),PID(N,M),PIDX,PIDY,PIDZ
C
 890   RL=RL-DL
       IF(RL.LE.0.001)GO TO 1000
       DXI=DXI-DXX
       DYI=DYI-DYY
       DZI=DZI-DZZ
       IF(DXI.GT.0.1)GO TO 900
       I=I+ID
       IF(I.GT.II.OR.I.LT.1)GO TO 1000
       DXI=DX(I,J,K)
 900   IF(DYI.GT.0.1)GO TO 910
       J=J+JD
       IF(J.GT.JJ.OR.J.LT.1)GO TO 1000
       DYI=DY(I,J,K)
 910   IF(DZI.GT.0.1)GO TO 920
       K=K+KD
       IF(K.GT.KK.OR.K.LT.1)GO TO 1000
       DZI=DZ(I,J,K)
 920   GO TO 210
1000   CONTINUE
C
607   FORMAT(//T15,'RESERVOIR CONTAINS FOLLOWING RATE NODES
     & (FOR SLANT WELLS) '/)
608   FORMAT(T2,'   NODE    WELL   LONG(FT)',2X,' OIL(STBD)',
     & 2X,'WATER(STBD)',2X,'GAS(MCFD)',2X,'TOTAL(RBD)',2X,'BHFP(PSIA)'
     & ,3X,'PID',5X,'PID-X',3X,'PID-Y',3X,'PID-Z')
703   FORMAT(1X,3I3,2X,A5,F8.2,2X,F10.2,4F12.2,F10.4,3F8.3)
       RETURN
       END 
C ======================================================================
       SUBROUTINE RONE(M,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RR,A2)
C ======================================================================
      COMMON /R1XYZ/ R11X,R21X,R31X,R41X,R11Y,R21Y,R31Y,R41Y,
     & R11Z,R21Z,R31Z,R41Z
       COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z

       IF(M.EQ.1)THEN
         R11=R11X
         R21=R21X
         R31=R31X
         R41=R41X
         T1=T1X
         T2=T2X
         T3=T3X
         T4=T4X
        ENDIF
       IF(M.EQ.2)THEN
         R11=R11Y
         R21=R21Y
         R31=R31Y
         R41=R41Y
         T1=T1Y
         T2=T2Y
         T3=T3Y
         T4=T4Y
        ENDIF
       IF(M.EQ.3)THEN
         R11=R11Z
         R21=R21Z
         R31=R31Z
         R41=R41Z
         T1=T1Z
         T2=T2Z
         T3=T3Z
         T4=T4Z
        ENDIF
C
       IF(I.EQ.1)THEN
        A2Z=DX0+2*XW
        R12Z=(DX0+XW)**2+(DYM+YW)**2
        R22Z=(DXP+DX0+XW)**2+YW**2
        R32Z=(DX0+XW)**2+(DYP-YW)**2
        R42Z=1
       ENDIF
       IF(J.EQ.1)THEN
        A2Z=DY0+2*YW
        R12Z=1
        R22Z=(DXP-XW)**2+(DY0+YW)**2
        R32Z=XW**2+(DYP+DY0+YW)**2
        R42Z=(DXM+XW)**2+(DY0+YW)**2  
       ENDIF
       IF(I.EQ.II)THEN
        A2Z=DX0-2*XW
        R12Z=(DX0-XW)**2+(DYM+YW)**2
        R22Z=1
        R32Z=(DX0-XW)**2+(DYP-YW)**2
        R42Z=(DXM+DX0-XW)**2+YW**2
       ENDIF
       IF(J.EQ.JJ)THEN
        A2Z=DY0-2*YW
        R12Z=XW**2+(DY0+DYM-YW)**2
        R22Z=(DXP-XW)**2+(DY0-YW)**2
        R32Z=1
        R42Z=(DXM+XW)**2+(DY0-YW)**2
       ENDIF
 450    R1Z=R11*R12Z
        R2Z=R21*R22Z
        R3Z=R31*R32Z
        R4Z=R41*R42Z
        RZ=T1*ALOG(R1Z)+T2*ALOG(R2Z)+T3*ALOG(R3Z)+T4*ALOG(R4Z)
        RR=RZ
        A2=A2Z
C        WRITE(6,1)R11,R21,R31,R41,R12Z,R22Z,R32Z,R42Z
  1     FORMAT(8E10.2)

        RETURN
        END
C ==============================================================
        SUBROUTINE RTWO1(DX0,DXP,XW,DY0,DYP,YW,T2,T3,RR,A)
C ==============================================================
        A2=DY0+2*YW
        A4=DX0+2*XW
        A3=(A2**2+A4**2)**.5
        A=A2*A3*A4
C        WRITE(6,1)A2,A3,A4,A
    1   FORMAT(4F10.1)
C
        R21=(DXP-XW)**2+YW**2
        R22=(DXP-XW)**2+(A2-YW)**2
        R23=(A4+DXP-XW)**2+(A2-YW)**2
        R24=(A4+DXP-XW)**2+YW**2
        R2=R21*R22*R23*R24
C
        R31=XW**2+(DYP-YW)**2
        R32=XW**2+(DYP+A2-YW)**2
        R33=(DX0+XW)**2+(DY0+DYP+YW)**2
        R34=(DX0+XW)**2+(DYP-YW)**2
        R3=R31*R32*R33*R34
        RR=T2*ALOG(R2)+T3*ALOG(R3)
C
        RETURN
        END
C
C =============================================================
        SUBROUTINE RTWO2(DX0,XW,DYM,DYP,YW,T1,T3,RR,A)
C =============================================================
C      
        A2=DX0-2*XW
        A3=DX0+2*XW
        A=A2*A3
C
        R11=XW**2+(DYM+YW)**2
        R12=(DX0-XW)**2+(DYM+YW)**2
        R13=(DX0+XW)**2+(DYM+YW)**2
        R1=R11*R12*R13
C
        R31=XW**2+(DYP-YW)**2
        R32=(DX0-XW)**2+(DYP-YW)**2
        R33=(DX0+XW)**2+(DYP-YW)**2
        R3=R31*R32*R33
C
        RR=T1*ALOG(R1)+T3*ALOG(R3)
        RETURN
        END
C
C =========================================================
        SUBROUTINE RTHREE(DX0,DXP,XW,DY0,DYP,YW,T3,RR,A)
C ==========================================================
C
        A2=DY0+2*YW
        A3=DX0-2*XW
        A4=DX0+2*XW
        A5=(A2*A2+A3*A3)**.5
        A6=(A2*A2+A4*A4)**.5
        A=A2*A3*A4*A5*A6
C
        R31=XW**2+(DYP-YW)**2
        R32=XW**2+(DY0+DYP+YW)**2
        R33=(DX0-XW)**2+(DYP-YW)**2
        R34=(DX0+XW)**2+(DYP-YW)**2
        R35=(DX0-XW)**2+(DY0+DYP+YW)**2
        R36=(DX0+XW)**2+(DY0+DYP+YW)**2
C 
        R3=R31*R32*R33*R34*R35*R36
        RR=T3*ALOG(R3)
C
        RETURN
        END
C
C ===================================================================
       SUBROUTINE TWOBOUND
C ===================================================================
      COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
      COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
      COMMON /NUMBER/ II,JJ,KK
      COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K
       IF(I.EQ.1.AND.J.EQ.1)THEN
        CALL RTWO1(DX0,DXP,XW,DY0,DYP,YW,T2Z,T3Z,RZ,AZ)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        ENDIF
C
       IF(I.EQ.1.AND.K.EQ.1)THEN
        CALL RTWO1(DZ0,DZP,ZW,DX0,DXP,XW,T2Y,T3Y,RY,AY)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II)THEN
        CALL RTWO2(DX0,XW,DYM,DYP,YW,T1Z,T3Z,RZ,AZ)
        CALL RTWO2(DX0,XW,DZP,DZM,-ZW,T2Y,T4Y,RY,AY)
        ENDIF
C
       IF(I.EQ.1.AND.J.EQ.JJ)THEN
        CALL RTWO1(DY0,DYM,-YW,DX0,DXP,XW,T1Z,T2Z,RZ,AZ)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        ENDIF
C
       IF(I.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO1(DX0,DXP,XW,DZ0,DZM,-ZW,T3Y,T4Y,RY,AY)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.K.EQ.1)THEN
        CALL RTWO1(DY0,DYP,YW,DZ0,DZP,ZW,T2X,T3X,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.I.EQ.II)THEN
        CALL RTWO1(DY0,DYP,YW,DX0,DXM,-XW,T3Z,T4Z,RZ,AZ)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ)THEN
        CALL RTWO2(DY0,YW,DXP,DXM,-XW,T2Z,T4Z,RZ,AZ)
        CALL RTWO2(DY0,YW,DZM,DZP,ZW,T1X,T3X,RX,AX)
        ENDIF
C
       IF(J.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO1(DZ0,DZM,-ZW,DY0,DYP,YW,T1X,T2X,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF

C
       IF(K.EQ.1.AND.I.EQ.II)THEN
        CALL RTWO1(DX0,DXM,-XW,DZ0,DZP,ZW,T1Y,T2Y,RY,AY)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       IF(K.EQ.1.AND.J.EQ.JJ)THEN
        CALL RTWO1(DZ0,DZP,ZW,DY0,DYM,-YW,T3X,T4X,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO2(DZ0,ZW,DYP,DYM,-YW,T2X,T4X,RX,AX)
        CALL RTWO2(DZ0,ZW,DXM,DXP,XW,T1Y,T3Y,RY,AY)
      ENDIF
C
       IF(I.EQ.II.AND.J.EQ.JJ)THEN
        CALL RTWO1(DX0,DXM,-XW,DY0,DYM,-YW,T4Z,T1Z,RZ,AZ)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        ENDIF
C
       IF(I.EQ.II.AND.K.EQ.KK)THEN
        CALL RTWO1(DZ0,DZM,-ZW,DX0,DXM,-XW,T4Y,T1Y,RY,AY)
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       IF(J.EQ.JJ.AND.K.EQ.KK)THEN
        CALL RTWO1(DY0,DYM,-YW,DZ0,DZM,-ZW,T4X,T1X,RX,AX)
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        ENDIF
C
       RETURN
       END
C
C ===================================================================
       SUBROUTINE THREEBOUND
C ===================================================================
      COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
      COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
      COMMON /NUMBER/ II,JJ,KK
      COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K
       IF(I.EQ.1.AND.J.EQ.1.AND.K.EQ.1)THEN
        CALL RTWO1(DY0,DYP,YW,DZ0,DZP,ZW,T2X,T3X,RX,AX)
        CALL RTWO1(DZ0,DZP,ZW,DX0,DXP,XW,T2Y,T3Y,RY,AY)
        CALL RTWO1(DX0,DXP,XW,DY0,DYP,YW,T2Z,T3Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.J.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO1(DZ0,DZM,-ZW,DY0,DYP,YW,T1X,T2X,RX,AX)
        CALL RTWO1(DX0,DXP,XW,DZ0,DZM,-ZW,T3Y,T4Y,RY,AY)
        CALL RTWO1(DY0,DYP,YW,DX0,DXP,XW,T3Z,T2Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.J.EQ.JJ.AND.K.EQ.1)THEN
        CALL RTWO1(DZ0,DZP,ZW,DY0,DYM,-YW,T3X,T4X,RX,AX)
        CALL RTWO1(DX0,DXP,XW,DZ0,DZP,ZW,T3Y,T2Y,RY,AY)
        CALL RTWO1(DY0,DYM,-YW,DX0,DXP,XW,T1Z,T2Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.J.EQ.JJ.AND.K.EQ.KK)THEN
        CALL RTWO1(DY0,DYM,-YW,DZ0,DZM,-ZW,T4X,T1X,RX,AX)
        CALL RTWO1(DZ0,DZM,-ZW,DX0,DXP,XW,T4Y,T3Y,RY,AY)
        CALL RTWO1(DX0,DXP,XW,DY0,DYM,-YW,T2Z,T1Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.II.AND.J.EQ.1.AND.K.EQ.1)THEN
        CALL RTWO1(DZ0,DZP,ZW,DY0,DYP,YW,T3X,T2X,RX,AX)
        CALL RTWO1(DX0,DXM,-XW,DZ0,DZP,ZW,T1Y,T2Y,RY,AY)
        CALL RTWO1(DY0,DYP,YW,DX0,DXM,-XW,T3Z,T4Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.II.AND.J.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO1(DY0,DYP,YW,DZ0,DZM,-ZW,T2X,T1X,RX,AX)
        CALL RTWO1(DZ0,DZM,-ZW,DX0,DXM,-XW,T4Y,T1Y,RY,AY)
        CALL RTWO1(DX0,DXM,-XW,DY0,DYP,YW,T4Z,T3Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.II.AND.J.EQ.JJ.AND.K.EQ.1)THEN
        CALL RTWO1(DY0,DYM,-YW,DZ0,DZP,ZW,T4X,T3X,RX,AX)
        CALL RTWO1(DZ0,DZP,ZW,DX0,DXM,-XW,T2Y,T1Y,RY,AY)
        CALL RTWO1(DX0,DXM,-XW,DY0,DYM,-YW,T4Z,T1Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.II.AND.J.EQ.JJ.AND.K.EQ.KK)THEN
        CALL RTWO1(DZ0,DZM,-ZW,DY0,DYM,-YW,T1X,T4X,RX,AX)
        CALL RTWO1(DX0,DXM,-XW,DZ0,DZM,-ZW,T1Y,T4Y,RY,AY)
        CALL RTWO1(DY0,DYM,-YW,DX0,DXM,-XW,T1Z,T4Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.J.EQ.1)THEN
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RTWO2(DX0,XW,DZP,DZM,-ZW,T2Y,T4Y,RY,AY)
        CALL RTHREE(DX0,DXP,XW,DY0,DYP,YW,T3Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.J.EQ.JJ)THEN
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RTWO2(DX0,XW,DZP,DZM,-ZW,T2Y,T4Y,RY,AY)
        CALL RTHREE(DX0,DXP,XW,DY0,DYM,-YW,T1Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.K.EQ.1)THEN
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RTWO2(DX0,XW,DYM,DYP,YW,T1Z,T3Z,RZ,AZ)
        CALL RTHREE(DX0,DXP,XW,DZ0,DZP,ZW,T2Y,RY,AY)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.K.EQ.KK)THEN
        CALL RONE(1,J,K,JJ,KK,YW,DY0,DYM,DYP,ZW,DZ0,DZM,DZP,RX,AX)
        CALL RTWO2(DX0,XW,DYM,DYP,YW,T1Z,T3Z,RZ,AZ)
        CALL RTHREE(DX0,DXP,XW,DZ0,DZM,-ZW,T4Y,RY,AY)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.I.EQ.1)THEN
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RTWO2(DY0,YW,DZM,DZP,ZW,T1X,T3X,RX,AX)
        CALL RTHREE(DY0,DYP,YW,DX0,DXP,XW,T2Z,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.I.EQ.II)THEN
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RTWO2(DY0,YW,DZM,DZP,ZW,T1X,T3X,RX,AX)
        CALL RTHREE(DY0,DYP,YW,DX0,DXM,-XW,T4Z,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.K.EQ.1)THEN
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RTWO2(DY0,YW,DXP,DXM,-XW,T2Z,T4Z,RZ,AZ)
        CALL RTHREE(DY0,DYP,YW,DZ0,DZP,ZW,T3X,RX,AX)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.K.EQ.KK)THEN
        CALL RONE(2,K,I,KK,II,ZW,DZ0,DZM,DZP,XW,DX0,DXM,DXP,RY,AY)
        CALL RTWO2(DY0,YW,DXP,DXM,-XW,T2Z,T4Z,RZ,AZ)
        CALL RTHREE(DY0,DYP,YW,DZ0,DZM,-ZW,T1X,RX,AX)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.I.EQ.1)THEN
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        CALL RTWO2(DZ0,ZW,DYP,DYM,-YW,T2X,T4X,RX,AX)
        CALL RTHREE(DZ0,DZP,ZW,DX0,DXP,XW,T3Y,RY,AY)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.I.EQ.II)THEN
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        CALL RTWO2(DZ0,ZW,DYP,DYM,-YW,T2X,T4X,RX,AX)
        CALL RTHREE(DZ0,DZP,ZW,DX0,DXM,-XW,T1Y,RY,AY)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.J.EQ.1)THEN
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        CALL RTWO2(DZ0,ZW,DXM,DXP,XW,T1Y,T3Y,RY,AY)
        CALL RTHREE(DZ0,DZP,ZW,DY0,DYP,YW,T2X,RX,AX)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.J.EQ.JJ)THEN
        CALL RONE(3,I,J,II,JJ,XW,DX0,DXM,DXP,YW,DY0,DYM,DYP,RZ,AZ)
        CALL RTWO2(DZ0,ZW,DXM,DXP,XW,T1Y,T3Y,RY,AY)
        CALL RTHREE(DZ0,DZP,ZW,DY0,DYM,-YW,T4X,RX,AX)
        ENDIF
C 
       RETURN
       END
C
C ===================================================================
       SUBROUTINE FOURBOUND
C ===================================================================
       COMMON /T1XYZ/ T1X,T2X,T3X,T4X,T1Y,T2Y,T3Y,T4Y,T1Z,T2Z,T3Z,T4Z
       COMMON /DXYZ/ DX0,DXP,DXM,XW,DY0,DYP,DYM,YW,DZ0,DZP,DZM,ZW
       COMMON /RXYZ/ RX,RY,RZ,AX,AY,AZ,I,J,K
       COMMON /NUMBER/ II,JJ,KK
C
       IF(I.EQ.1.AND.I.EQ.II.AND.J.EQ.1.AND.K.EQ.1)THEN
        CALL RTWO1(DY0,DYP,YW,DZ0,DZP,ZW,T2X,T3X,RX,AX)
        CALL RTHREE(DX0,DXP,XW,DZ0,DZP,ZW,T2Y,RY,AY)
        CALL RTHREE(DX0,DXP,XW,DY0,DYP,YW,T3Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.J.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO1(DZ0,DZM,-ZW,DY0,DYP,YW,T1X,T2X,RX,AX)
        CALL RTHREE(DX0,DXP,XW,DZ0,DZM,-ZW,T4Y,RY,AY)
        CALL RTHREE(DX0,DXP,XW,DY0,DYP,YW,T3Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.J.EQ.JJ.AND.K.EQ.1)THEN
        CALL RTWO1(DZ0,DZP,ZW,DY0,DYM,-YW,T3X,T4X,RX,AX)
        CALL RTHREE(DX0,DXP,XW,DZ0,DZP,ZW,T2Y,RY,AY)
        CALL RTHREE(DX0,DXP,XW,DY0,DYM,-YW,T1Z,RZ,AZ)
        ENDIF
C
       IF(I.EQ.1.AND.I.EQ.II.AND.J.EQ.JJ.AND.K.EQ.KK)THEN
        CALL RTWO1(DY0,DYM,-YW,DZ0,DZM,-ZW,T4X,T1X,RX,AX)
        CALL RTHREE(DX0,DXP,XW,DZ0,DZM,-ZW,T4Y,RY,AY)
        CALL RTHREE(DX0,DXP,XW,DY0,DYM,-YW,T1Z,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.I.EQ.1.AND.K.EQ.1)THEN
        CALL RTWO1(DZ0,DZP,ZW,DX0,DXP,XW,T2Y,T3Y,RY,AY)
        CALL RTHREE(DY0,DYP,YW,DZ0,DZP,ZW,T3X,RX,AX)
        CALL RTHREE(DY0,DYP,YW,DX0,DXP,XW,T2Z,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.I.EQ.1.AND.K.EQ.KK)THEN
        CALL RTWO1(DX0,DXP,XW,DZ0,DZM,-ZW,T3Y,T4Y,RY,AY)
        CALL RTHREE(DY0,DYP,YW,DZ0,DZM,-ZW,T1X,RX,AX)
        CALL RTHREE(DY0,DYP,YW,DX0,DXP,XW,T2Z,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.I.EQ.II.AND.K.EQ.1)THEN
        CALL RTWO1(DX0,DXM,-XW,DZ0,DZP,ZW,T1Y,T2Y,RY,AY)
        CALL RTHREE(DY0,DYP,YW,DZ0,DZP,ZW,T3X,RX,AX)
        CALL RTHREE(DY0,DYP,YW,DX0,DXM,-XW,T4Z,RZ,AZ)
        ENDIF
C
       IF(J.EQ.1.AND.J.EQ.JJ.AND.I.EQ.II.AND.K.EQ.KK)THEN
        CALL RTWO1(DZ0,DZM,-ZW,DX0,DXM,-XW,T4Y,T1Y,RY,AY)
        CALL RTHREE(DY0,DYP,YW,DZ0,DZM,-ZW,T1X,RX,AX)
        CALL RTHREE(DY0,DYP,YW,DX0,DXM,-XW,T4Z,RZ,AZ)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.I.EQ.1.AND.J.EQ.1)THEN
        CALL RTWO1(DX0,DXP,XW,DY0,DYP,YW,T2Z,T3Z,RZ,AZ)
        CALL RTHREE(DZ0,DZP,ZW,DY0,DYP,YW,T2X,RX,AX)
        CALL RTHREE(DZ0,DZP,ZW,DX0,DXP,XW,T3Y,RY,AY)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.I.EQ.1.AND.J.EQ.JJ)THEN
        CALL RTWO1(DY0,DYM,-YW,DX0,DXP,XW,T1Z,T2Z,RZ,AZ)
        CALL RTHREE(DZ0,DZP,ZW,DY0,DYM,-YW,T4X,RX,AX)
        CALL RTHREE(DZ0,DZP,ZW,DX0,DXP,XW,T3Y,RY,AY)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.I.EQ.II.AND.J.EQ.1)THEN
        CALL RTWO1(DY0,DYP,YW,DX0,DXM,-XW,T3Z,T4Z,RZ,AZ)
        CALL RTHREE(DZ0,DZP,ZW,DY0,DYP,YW,T2X,RX,AX)
        CALL RTHREE(DZ0,DZP,ZW,DX0,DXM,-XW,T1Y,RY,AY)
        ENDIF
C
       IF(K.EQ.1.AND.K.EQ.KK.AND.I.EQ.II.AND.J.EQ.JJ)THEN
        CALL RTWO1(DX0,DXM,-XW,DY0,DYM,-YW,T4Z,T1Z,RZ,AZ)
        CALL RTHREE(DZ0,DZP,ZW,DY0,DYM,-YW,T4X,RX,AX)
        CALL RTHREE(DZ0,DZP,ZW,DX0,DXM,-XW,T1Y,RY,AY)
        ENDIF

       RETURN
       END
C
C ====================================================================
C
C      SUBROUTINE RESTART(INFILE,RSTFILE)
C
C      SUBROUTINE RESTART(INFILE,RSTFILE,NX,NY,NZ)
      SUBROUTINE RESTART(INFILE,RSTFILE,NX,NY,NZ,P,SO,SW,DPP,
     &  DPSO,DPSW)
C
C ====================================================================
C     CREATING A "RESTART.DAT" FILE FOR A FOLLOWED RUN AFTER A SPECIFIED
C      SIMULATION TIME. THE DISTRIBUTIONS OF RESERVOIR PRESSURE AND
C      OIL/WATER SATURATION ARE SAVED. 
C      INCLUDE 'D:\F77L\BOASTVHS\VHSPARMA.FOR'
C      INCLUDE 'D:\F77L\BOASTVHS\VHSPARMB.FOR'
C	NEW DP LOCAL VARIABLES **RJH**
C       MODIFICATION FOR DOUBLE PRECISION READ OF SO & SW IN RESTART  **RJH**

      CHARACTER DATA(80)
      DIMENSION P(NX,NY,NZ)
      DIMENSION SO(NX,NY,NZ),SW(NX,NY,NZ)
      REAL*8 DPP(NX,NY,NZ),DPSO(NX,NY,NZ),DPSW(NX,NY,NZ)
      COMMON /NUMBER/ II,JJ,KK
      COMMON /IO/ NI,NO
      COMMON /ADD1/ IM,JM,KM,ETI,FT,FTMAX,DELT
      REAL*8 COP,CWP,CGP,CWI,CGI
      COMMON /ADD2/ COP,CWP,CGP,CWI,CGI
      COMMON /COUNT/ N1READ,N2READ
      CHARACTER*20 INFILE,RSTFILE

      DO 3031 K=1,KK
      DO 3032 J=1,JJ
      DO 3033 I=1,II
      DPP(I,J,K)=P(I,J,K)
      DPSO(I,J,K)=SO(I,J,K)
      DPSW(I,J,K)=SW(I,J,K)
3033  CONTINUE
3032  CONTINUE
3031  CONTINUE
C	END OF MODIFICATION **RJH**
      OPEN(10,FILE=RSTFILE,STATUS='UNKNOWN')
      WRITE(*,3003) RSTFILE
3003  FORMAT(10X,'Writing a restart file to : ',A20)
      OPEN(NI,FILE=INFILE,STATUS='OLD')
      DO 100 N=1, N1READ
       READ(NI,110) (DATA(I),I=1,80)
       WRITE(10,110) (DATA(I),I=1,80)
 100  CONTINUE
      WRITE(10,70)
      WRITE(10,80)
      DO 5 K=1,KK
      DO 5 J=1,JJ
C       WRITE(10,60) (P(I,J,K),I=1,II)
       WRITE(10,*) (DPP(I,J,K),I=1,II)
   5  CONTINUE
      DO 10 K=1,KK
      DO 10 J=1,JJ
C       WRITE(10,50) (SO(I,J,K),I=1,II)
       WRITE(10,*) (DPSO(I,J,K),I=1,II)
  10  CONTINUE
      DO 20 K=1,KK
      DO 20 J=1,JJ
C       WRITE(10,50) (SW(I,J,K),I=1,II)
       WRITE(10,*) (DPSW(I,J,K),I=1,II)
  20  CONTINUE
      WRITE(10,*) ETI,DELT
      WRITE(10,*)  COP,CWP,CGP,CWI,CGI
      DO 115 N=1,N2READ
       READ(NI,110) (DATA(I),I=1,80)
 115  CONTINUE
      DO 120 N=1,9999
       READ(NI,110,END=200,ERR=200) (DATA(I),I=1,80)
       WRITE(10,110) (DATA(I),I=1,80)
 120  CONTINUE
C  50  FORMAT(10F8.5)
C  60  FORMAT(10F8.2)
C  50  FORMAT(7(F10.7,1X))
C  60  FORMAT(7(F10.4,1X))
  70  FORMAT('RESTART RUN -- PRESSURE & OIL/WATER SATURATION
     & OBTAINED FROM PREVIOUS RUN')
  80  FORMAT('    1    1')
 110  FORMAT(80A1)
 200  CLOSE(NI)
      CLOSE(10,STATUS='KEEP')
      STOP
      END 
