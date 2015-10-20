c -- $Id: init.F,v 1.2 2001/07/03 21:29:19 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/init.F,v $
c -- $Name:  $

C ***************************** BASIC *********************************
C
      SUBROUTINE BASIC
C----------------------------------------------------------------------
C
CL            INITIALIZE BASIC CONTROL DATA
C
C---------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMDIF.h>
C---------------------------------------------------------------------
#include <COMFLG.h>
C---------------------------------------------------------------------
#include <COMHED.h>
C---------------------------------------------------------------------
#include <COMHIS.h>
C---------------------------------------------------------------------
#include <COMIMP.h>
C---------------------------------------------------------------------
#include <COMINI.h>
C---------------------------------------------------------------------
#include <COML10.h>
C---------------------------------------------------------------------
#include <COML30.h>
C---------------------------------------------------------------------
#include <COML50.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C---------------------------------------------------------------------
#include <COMMAP.h>
C---------------------------------------------------------------------
#include <COMMUN.h>
C---------------------------------------------------------------------
#include <COMSDT.h>
C---------------------------------------------------------------------
#include <COMTIM.h>
C---------------------------------------------------------------------
#include <CRDCTL.h>
C----------------------------------------------------------------------
#include <PARA_I.h>
C----------------------------------------------------------------------
C
CL            0.1   GENESIS MODEL (VERSGEN), HISTORY-FILE (VERSHIS)
CL                  AND RESTART-FILE (VERSRES) VERSION NUMBERS.
CL                  ALSO MACHINE-REPRESENTATION TYPES (TYPEHIS,TYPERES).
C
      VERSGEN = '3.0'
      VERSHIS = 1.02
#if defined (LKAEROSOL)
      VERSRES = 3.11
#else
      VERSRES = 3.10
#endif

#ifdef isotrac
      VERSRES_ISOT = 1.00
#else
      VERSRES_ISOT = 0.00
#endif

#ifdef cray
      TYPEHIS = 'CRAY'
      TYPERES = 'CRAY'
#else
      TYPEHIS = 'IEEE'
      TYPERES = 'IEEE'
#endif
C
CL            1.0   SET LOGICAL UNIT NUMBERS IN COMLUN. LOCIU IS
CL                  EQUIVALENCED TO UNIT NUMBERS FOR DATA INPUT FILES
CL                  AND RESTART/HISTORY FILES (SEE COMMONS), AND IS
CL                  USED FOR LOOKUP OF DATA INPUT FILES IN FINDIU.
C                     5,6   = NIN, NOUT
C                     11-.. = LOCIU...
C                     50-55 = LSX/LAND, 54 = DYNICE (scratch files)
C                     60-65 = LSX/WATER (scratch files)
C                     66-68 = NUSER* (68 for point output in surfctl)
C                     70    = PLUME (scratch file)
C                     71-76 = TRABUDG (scratch files)
C                     77-79 = used by debplot*, conv, access
C                     80-99 = reserved for OGCM
C
      NIN = 5
#if ! defined (coupogcm) 
      NOUT = 6
#endif
      NUSER0 = 66
      NUSER1 = 67
      NUSER2 = 68
      DO 10 I=1,PLUN
        LOCIU(I) = 10+I
   10 CONTINUE
C
CL            2.0   WRITE MODEL BANNER
C
      CALL BANNER
C
CL            3.8   INITIALIZE VERTICAL GRID (SO READERV CAN USE SIG)
C
      CALL SETSIG
C
CL            4.0   READ NAMELIST INPUT AND INITIAL/RESTART FILE
C
      CALL DATGET
C
CL            4.2   INITIALIZE COMCON, REST OF COMMAP
C
      CALL SETCOM
C
CL            4.4   INITIALIZE CONSTANTS
C
      CALL SETINI
C
CL            5.0   READ TOPOGRAPHY FIELD AND S.D. FIELDS
C
      CALL TOPINI
C
CL            6.0   SET AGCM FIELDS FOR COLD START
C
      IF (NSREST.EQ.-1) CALL COLDSTRT
C
CL            7.0   SET YTSIN, TO BE USED FOR INITIALIZATION OF LSX FOR
CL                  COLD STARTS (NSREST=-1). (DOESN'T *HAVE* TO EQUAL
CL                  LOWEST-LEVEL YT, BUT USED HERE FOR CONVENIENCE.)
C
      CALL SCOPY (NLON*NOREC, YT(1,1,NLEV), 1, YTSIN, 1)
C
CL            7.5   FOR COLD STARTS OR OLD RESTART-FILE VERSIONS,
CL                  CRUDELY SET NEEDED FIELDS (OTHERWISE THEY HAVE
CL                  BEEN READ OR SET IN RSDS)
C
      IF (NSREST.EQ.-1 .OR. VERSRESIN.LE.2.07) THEN
         CALL ZERO (TMPBUD,  NBUDG)
         NBUDGPRE = 0
         CALL ZERO (TRACER,  NLON*NOREC*NLEV*NTRACE)
         CALL ZERO (FRAINC,  NLON*NOREC)
         CALL ZERO (FSNOWC,  NLON*NOREC)
         CALL ZERO (FRAINS,  NLON*NOREC)
         CALL ZERO (FSNOWS,  NLON*NOREC)
         IF (VERSRESIN.LT.1.05) THEN
            CALL ZERO (FPLUMT,  NLON*NOREC*NLEV)
            CALL ZERO (FPLUMQ,  NLON*NOREC*NLEV)
            CALL ZERO (FPLUMU,  NLON*NOREC*NLEV)
            CALL ZERO (FPLUMV,  NLON*NOREC*NLEV)
         ENDIF
C
         DO 710 JJ=1,NOREC
            DO 712 JI=1,NLON
               FIREMS(JI,JJ) = STEBOL*YTSIN(JI,JJ)**4
  712       CONTINUE
  710    CONTINUE
C
         DO 720 JK=1,NLEV
            DO 722 JJ=1,NOREC
               DO 724 JI=1,NLON
                  ALONSAV(JI,JJ,JK) = ALONG(JI)
                  ALATSAV(JI,JJ,JK) = ALATG(JJ)
                  SIGSAV (JI,JJ,JK) = SIG(JK)
  724          CONTINUE
  722       CONTINUE
  720    CONTINUE
         CALL ZERO (WSLTSAV, NLON*NOREC*NLEV)
      ENDIF

c        Initialization of isotope's and reg. tracers for cold start

c>>>>>>>>>>>>>
#ifdef isotrac
c>>>>>>>>>>>>>

c Either:
c     call findiu (iusourc_i, loc)
c#ifdef gcm
c     write(ioterm2,800) locname(loc)(1:lenchr(locname(loc)))
c 800 format(/' Reading isotope source region data file'/4x,a)
c#endif
c     call openf (iusourc_i, locname(loc), ' ', 0, 0)
c     call readdata (iusourc_i,sourc_i,nlon,norec,'SOURC_I ',1,1)
c     call closef(iusourc_i)
c Or:
      call zero (sourc_i, nlon*norec)

      if (nsrest_isot.eq.-1) then
        call zero (rtrq,   nlon*norec*nlev*nrtr)
        call zero (rtrqm1, nlon*norec*nlev*nrtr)
        call zero (rtrfq,  nlon*norec*nlev*nrtr)
        call zero (rtrprc, nlon*norec*nrtr)
        call zero (rtrprs, nlon*norec*nrtr)
        call zero (rtrpsc, nlon*norec*nrtr)
        call zero (rtrpss, nlon*norec*nrtr)
        v_o = ( -80./1000.+1.)*smoj_o
        v_h = (-600./1000.+1.)*smoj_h
        do 100 ji=1,nlon
          do 110 jj=1,norec
            frainc_o(ji,jj)    = v_o*frainc(ji,jj)
            frainc_h(ji,jj)    = v_h*frainc(ji,jj)
            frains_o(ji,jj)    = v_o*frains(ji,jj)
            frains_h(ji,jj)    = v_h*frains(ji,jj)
            fsnowc_o(ji,jj)    = v_o*fsnowc(ji,jj)
            fsnowc_h(ji,jj)    = v_h*fsnowc(ji,jj)
            fsnows_o(ji,jj)    = v_o*fsnows(ji,jj)
            fsnows_h(ji,jj)    = v_h*fsnows(ji,jj)
            do 120 jk=1,nlev
              if (yt(ji,jj,jk).ge.t_r_i) then
                w_o = v_o*avl_o(yt(ji,jj,jk))
                w_h = v_h*avl_h(yt(ji,jj,jk))
              else
                w_o = v_o*avi_o(yt(ji,jj,jk))
                w_h = v_h*avi_h(yt(ji,jj,jk))
              endif
              yq_o(ji,jj,jk)     = v_o*yq(ji,jj,jk)
              yq_h(ji,jj,jk)     = v_h*yq(ji,jj,jk)
              yqm1_o(ji,jj,jk)   = v_o*yqm1(ji,jj,jk)
              yqm1_h(ji,jj,jk)   = v_h*yqm1(ji,jj,jk)
              ylwc_c_o(ji,jj,jk) = w_o*ylwc_c(ji,jj,jk)
              ylwc_c_h(ji,jj,jk) = w_h*ylwc_c(ji,jj,jk)
              ylwc_a_o(ji,jj,jk) = w_o*ylwc_a(ji,jj,jk)
              ylwc_a_h(ji,jj,jk) = w_h*ylwc_a(ji,jj,jk)
              ylwc_s_o(ji,jj,jk) = w_o*ylwc_s(ji,jj,jk)
              ylwc_s_h(ji,jj,jk) = w_h*ylwc_s(ji,jj,jk)
              fplumq_o(ji,jj,jk) = v_o*fplumq(ji,jj,jk)
              fplumq_h(ji,jj,jk) = v_h*fplumq(ji,jj,jk)
 120        continue
 110      continue
 100    continue
      endif
c>>>>>
#endif
c>>>>>

C
CL            8.0   HARD-CODED INITIALIZATION OF TRACERS AND
CL                  SURFACE FLUXES #NTRACA+1 TO #NTRACE
CL                  (#1-3 ARE LIQUID WATER CONTENTS, #4-7 ARE FLUX
CL                  RESERVOIRS, NTRACA = 7)
 
      CALL ZERO (TRAFLUX, NLON*NOREC*(NTRACE-NTRACA))
 
c--------------------------
      IF (NSREST.LE.0) THEN
c--------------------------
 
       DO N=NTRACA+1,NTRACE
         CALL ZERO(TRACER(1,1,1,N),NLON*NOREC*NLEV)
ctest    DO JK=1,6
ctest      DO JJ=1,NOREC
ctest        DO JI=1,NLON
ctest          TRACER(JI,JJ,JK,N) = 1.
ctest        ENDDO
ctest      ENDDO
ctest    ENDDO
       ENDDO
 
c>>>>>>>>>>>>>
#ifndef cosmog
c>>>>>>>>>>>>>
c=====
cc     El Chicon location (R15 lon,lat indices!):
c      DO 800 JK=1,NLEV
c       IF(SIG(JK).GT..015.AND.SIG(JK).LE..055) TRACER(24,36,JK, 8)=999.
c       IF(SIG(JK).GT..055.AND.SIG(JK).LE..095) TRACER(24,36,JK, 9)=999.
c       IF(SIG(JK).GT..095.AND.SIG(JK).LE..135) TRACER(24,36,JK,10)=999.
c 800  CONTINUE
c=====
cc     Horizontally uniform:
c      DO 810 N=8,10
c        IF (N.EQ. 8) JKZ = 2
c        IF (N.EQ. 9) JKZ = 3
c        IF (N.EQ.10) JKZ = 4
c        DO 812 JK = MAX(JKZ-2,1), MIN(JKZ+2,NLEV)
c          ZVAL = 1. * EXP( -(((SIG(JK)-SIG(JKZ))/.020)**2) )
c          DO 814 JJ=1,NOREC
c            DO 816 JI=1,NLON
c              TRACER(JI,JJ,JK,N) = ZVAL
c 816        CONTINUE
c 814      CONTINUE
c 812    CONTINUE
c 810  CONTINUE
 
c      IF (NTRACA.LT.NTRACE) THEN
cc       Upper 4 levels:
c        DO JK = 1,4
c          DO JJ=1,NOREC
c            DO JI=1,NLON
c              TRACER(JI,JJ,JK,NTRACA+1) = 1.
c            ENDDO
c          ENDDO
c        ENDDO
 
cc       Above "stratopause":
c        DO JJ=1,NOREC
c          ZCOS = COS ( 3.1415927*(-0.5 + (JJ-.5)/NOREC) )
c          ZSIGPAUS = .100*(ZCOS**2) + .280*(1.-ZCOS**2)
c          DO JK = 1,NLEV
c            IF (SIG(JK).LE.ZSIGPAUS) THEN
c              DO JI=1,NLON
c                TRACER(JI,JJ,JK,NTRACA+1) = 1.
c              ENDDO
c            ENDIF
c          ENDDO
c        ENDDO
c      ENDIF
c>>>>>
#endif
c>>>>>
 
c----------
      ENDIF
c----------
 
      RETURN
      END
C
C ***************************** DATGET ********************************
C
      SUBROUTINE DATGET
C----------------------------------------------------------------------
C
CL            READS AND PROCESSES NAMELIST INPUT, INITIAL/RESTART FILE
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMDIF.h>
C ---------------------------------------------------------------------
#include <COMDYR.h>
C ---------------------------------------------------------------------
#include <COMFLG.h>
C ---------------------------------------------------------------------
#include <COMHED.h>
C ---------------------------------------------------------------------
#include <COMHIS.h>
C ---------------------------------------------------------------------
#include <COMIMP.h>
C ---------------------------------------------------------------------
#include <COMLUN.h>
C ---------------------------------------------------------------------
#include <COML30.h>
C ---------------------------------------------------------------------
#include <COMMUN.h>
C ---------------------------------------------------------------------
#include <COMSDT.h>
C ---------------------------------------------------------------------
#include <COMTIM.h>
C ---------------------------------------------------------------------
#include <CRDCTL.h>
C ---------------------------------------------------------------------
#include <LGASES.h>
C ---------------------------------------------------------------------
      CHARACTER*80 RUNTITLEX
      CHARACTER*8 CTMP
C ---------------------------------------------------------------------

c        locname holds the names used in the program (see subr findiu,
c        corresponding to list of unit numbers, both in common comlun).
c        Set default filenames here, which can be overrridden by 
c        Namelist input:

      data (locname(i),i=1,presoil) /
     *     'TOPOG',    'OZONE',    'GWD',      'TOPOG2',
     *     'DSBIOME',  'WISCVEG',  'SURFTYP',  'FWATER',
     *     'OCEAN',    'SSTICE', 
     *     'OCEAN_I',  'SOURC_I',  'interm' /
      data (locname(i),i=presoil+1, presoil+psoilun) /
     *     'SOIT01', 'SOIT02', 'SOIT03', 'SOIT04', 'SOIT05', 'SOIT06',
     *     'SOIT07', 'SOIT08', 'SOIT09', 'SOIT10', 'SOIT11', 'SOIT12' /
c     do i=1,psoilun
c       write (locname(presoil+i),"('SOIT',i2.2)") i
c     enddo

c        Flag variables for inputting data input filenames in Namelist,
c        corresponding to locname(1:presoil) in comlun  

      character*96
     *     TOPOG_FILE,     OZONE_FILE,     GWD_FILE,      TOPOG2_FILE,
     *     VEG_FILE,       WISCVEG_FILE,   SURFTYPE_FILE, FWATER_FILE,
     *     OCEANUV_FILE,   SST_FILE, 
     *     ISOT_O_FILE,    ISOT_S_FILE,    DUMMY_INTERM,
     *     SOILTEX01_FILE, SOILTEX02_FILE, SOILTEX03_FILE,
     *     SOILTEX04_FILE, SOILTEX05_FILE, SOILTEX06_FILE,
     *     SOILTEX07_FILE, SOILTEX08_FILE, SOILTEX09_FILE,
     *     SOILTEX10_FILE, SOILTEX11_FILE, SOILTEX12_FILE
      parameter (nblank=presoil+psoilun)
      data 
     *     TOPOG_FILE,     OZONE_FILE,     GWD_FILE,      TOPOG2_FILE,
     *     VEG_FILE,       WISCVEG_FILE,   SURFTYPE_FILE, FWATER_FILE,
     *     OCEANUV_FILE,   SST_FILE, 
     *     ISOT_O_FILE,    ISOT_S_FILE,    DUMMY_INTERM,
     *     SOILTEX01_FILE, SOILTEX02_FILE, SOILTEX03_FILE,
     *     SOILTEX04_FILE, SOILTEX05_FILE, SOILTEX06_FILE,
     *     SOILTEX07_FILE, SOILTEX08_FILE, SOILTEX09_FILE,
     *     SOILTEX10_FILE, SOILTEX11_FILE, SOILTEX12_FILE
     *    / nblank * ' ' /

C ---------------------------------------------------------------------

      NAMELIST /INPUT/
C
     *  RUNTITLE,    BASEDATE,
C
     *  NSREST,      ENDDATE,
     *  NSREST_ISOT,
C
     *  MSPATHIN,    MSNAMIN,
     *  MSPATH,      MSNAMHIS,   MSNAMRES,   MSNAMRE2,
     *  PASSRD,      PASSWD,      RETPD,      FATALWMS,   HISNET,
C
     *  CUSHNAM,     NDENS,      NDENS_STORM,
     *  NHISI,       NDHIS,      HISTLAST,
C
     *  REGPATH,     REGHIS,     REGHIS_A,   REGFREQ,
C
     *  STORMPATH,   STORMHIS,               STORMFREQ,
C
     *  IBISPATHIN,  IBISRESIN,  IBISRES,
C
     *  NSAVE,       NNEXF,      DTIME,      EPS,
     *  TRADLW,      TRADLWA,    TRADSW,
C
     *  ECCU,        OBLU,       PRECU,
     *  TEQCOLD,     TPOCOLD,    PERTURB,
     *  PRIESTQ,     STRATTIME,  STRATFLAG,  SMOLAR,    TRAVDIF,TRACONV,
     *  VEGTYPE,     EVEDELTA,   DSBIOME,    WISCVEG,   
     *  PHOTOTYPE,   TAUBIOME,
     *  NSREST_IBIS, DYNAIBIS,   ASYNIBIS,
     *  SOILTEX,     BUDFREQ, 
     *  STOCHPR,     SHOWMAPS,   LONOLSX,    LATOLSX,   TIMOLSX,
     *  DEPTHML,     DEPTHLAKE,  DEPTHICE,   DEPTHSNO,
     *  OCEANTYPE,   SEQSST,     DYNAMICE,   DYNAMHEM,PRESCOUV,
     *  QFACTOR,     QNORWEG,    QICES,      QICEN,   QLIMIT,
     *  TRIASSIC,    HYBTYPE,    GRAVFAC,    GWROUGH,
     *  ENVFAC,      DELEVTYPE,
     *  LEVRES,      HORRES,     HORRESLSX,  SCOURANT,     HOTLIMIT,
     *  DIF2,        DIF4,       TIMRESERV,
     *  PLUMERAD1, PLUMEAREA1_L, PLUMEAREA1_O, PLUMERAD2,  PLUMEAREA2,
C
     *  ALW,
     *  BLW_C,     BLW_A,     BLW_S,
     *  CLWFAC_C,  CLWFAC_A,  CLWFAC_S,  
     *  CLW_P,
     *  DLW_LL,    DLW_OL,    DLW_LI,   DLW_OI, 
     *  DLWFAC,    
     *  REFFAC_L,  REFFAC_O,  REFFAC_I, RADLIQ,
     *  ELW,       
     *  FLW_C,     FLW_A,     FLW_S,
     *  GLW,
     *  TLIQA,     TLIQB,     SIGANVIL,
     *  FACCLWP_SOL, FACCLWP_IR,
     *  SSACLD,
     *  DAYLENCOR,  DAYLENRAD, SSLON,

     *  TOPOG_FILE,     OZONE_FILE,     GWD_FILE,      TOPOG2_FILE,
     *  VEG_FILE,       WISCVEG_FILE,   SURFTYPE_FILE, FWATER_FILE,
     *  OCEANUV_FILE,   SST_FILE, 
     *  ISOT_O_FILE,    ISOT_S_FILE, 
     *  SOILTEX01_FILE, SOILTEX02_FILE, SOILTEX03_FILE,
     *  SOILTEX04_FILE, SOILTEX05_FILE, SOILTEX06_FILE,
     *  SOILTEX07_FILE, SOILTEX08_FILE, SOILTEX09_FILE,
     *  SOILTEX10_FILE, SOILTEX11_FILE, SOILTEX12_FILE

#ifdef isotrac
     * ,t_s_i,     s_s_i,     hi_s_i,    hs_s_i,    r_s_i,
     *  t_c_i,     s_c_i,     hi_c_i,    hs_c_i,    r_c_i,
     *  t_r_i,     equ_i,
     *  pii_i,     pjj_i,     pdn_i,     hrst_i
#endif
 
C  NAMELIST FOR OTG CODE
 
      NAMELIST /RADOTG/
     *  SOLFAC,  SOLPAR,
     *  CO2PPV, CH4PPV, N2OPPV, F11PPV, F12PPV, O3FAC, O2PPM,
     *  AERTAUVIS, MXAERL, GGUNIF, NOCLDTOP, NOCLDBOT
C
C----------------------------------------------------------------------
C
CL        1.0        SET MUST-BE-ENTERED VARIABLES TO NULL
C
      ENDDATE  = ' '
      NSREST = -9999
      NSREST_ISOT = -9999
      MSNAMIN = ' '
      MSPATH = 'null'
C
CL        1.1        SET OPTIONALLY-ENTERED VARIABLES THAT MAY BE
CL                   OBTAINED FROM THE RESTART FILE TO NULL
C
      RUNTITLE = ' '
      BASEDATE = ' '
C
CL        1.2        SET OPTIONALLY-ENTERED ARRAY VARIABLES TO NULL
CL                   FOR NAMELIST INPUT. (SINCE TRAILING WORDS ARE LEFT
CL                   ALONE BY NAMELIST, CANNOT SET DEFAULTS TILL AFTER
CL                   NAMELIST IS READ)
C
      DO 120 I=1,PFILE
        MSNAMHIS(I) = ' '
        MSNAMRES(I) = ' '
  120 CONTINUE
      MSNAMRE2 = ' '
      DO 122 I=1,NMARK
        NHISI(I) = 0
        NSAVE(I) = 0
        NNEXF(I) = 0
  122 CONTINUE
      DO 124 I=1,NUSHMX
        CUSHNAM(I) = ' '
  124 CONTINUE
C
CL        1.3        SET OPTIONALLY-ENTERED SINGLE VARIABLES TO DEFAULTS
CL                   FOR NAMELIST INPUT
C
      MSPATHIN = 'null'
      PASSRD = ' '
      PASSWD = ' '
      RETPD = 30.
      FATALWMS = .FALSE.
      HISNET = .TRUE.
      NDENS = 1
      NDENS_STORM = 1
      NDHIS = 0
      HISTLAST = .FALSE.
      REGPATH = 'null'
      REGHIS = ' '
      REGHIS_A = ' '
      REGFREQ = 12.
      STORMPATH = 'null'
      STORMHIS = ' '
      STORMFREQ = 12.
      IBISPATHIN = 'null'
      IBISRESIN = ' '
      IBISRES = ' '
      EPS = .06
      TRADLW  = 1.5
      TRADLWA = 24.
      TRADSW  = 1.5
      ECCU = -9999.
      OBLU = -9999.
      PRECU =-9999.
      TEQCOLD = 303.
      TPOCOLD = 260.
      PERTURB = 0.
      PRIESTQ   = .FALSE. ! Priestley fixer for q, else +/- method for q
      STRATTIME = 100.    ! wfilter time scale in days
      STRATFLAG = 0       ! = 0 is off, > 0 is on, < 0 is use zero w's
      SMOLAR    = .FALSE. ! Smolarkiewicz "anti-diffusion",vertical slt
      TRAVDIF   = .TRUE.  ! vertically diffuse tracers
      TRACONV   = .TRUE.  ! plume-convect tracers
      VEGTYPE = 0
      EVEDELTA = 0
      DSBIOME = 0
      WISCVEG = 0
      PHOTOTYPE = 0
      TAUBIOME = 4.
      NSREST_IBIS = 1
      DYNAIBIS = 0
      ASYNIBIS = 1
      SOILTEX = 0
      BUDFREQ = 365.
      STOCHPR = .FALSE.
      SHOWMAPS = .FALSE.
      LONOLSX = -999 ! 295
      LATOLSX = -999 ! -3
      TIMOLSX = 3.
      DEPTHML = 0.
      DEPTHLAKE = -1.
      DEPTHICE = 0.
      DEPTHSNO = 0.
      OCEANTYPE = 1
      SEQSST = .FALSE.
      DYNAMICE = .TRUE.
      DYNAMHEM = .TRUE.
      PRESCOUV = .TRUE.
      QFACTOR = 1.0
      QNORWEG = 500.
      QICES = 15.
      QICEN = 0.
      QLIMIT = 0.
      TRIASSIC = .FALSE.
      HYBTYPE = 2
      SCOURANT = -999.
      HOTLIMIT = -999.
      IF (NLON.EQ.32 .AND. NOREC.EQ.16) THEN            ! T10
        DTIME = 3600.
        GRAVFAC = 2.
        DIF2 = 0.5E5
        DIF4 = 0.4E16
      ELSE IF (NLON.EQ.48 .AND. NOREC.EQ.40) THEN       ! R15
        DTIME = 1800.   ! v1.02
        GRAVFAC = 0.    ! v1.02
        DIF2 = 2.5E5    ! v1.02
        DIF4 = 2.0E16   ! v1.02
      ELSE IF (NLON.EQ.96 .AND.NOREC.EQ.48) THEN        ! T31
        DTIME = 1800.   ! v2.0
        GRAVFAC = 2.    ! v2.0
        DIF2 = 0.5E5    ! v2.0
        DIF4 = 0.4E16   ! v2.0
      ELSE IF (NLON.EQ.128 .AND.NOREC.EQ.64) THEN       ! T42
        DTIME = 1200.
        GRAVFAC = 2.
        DIF2 = 0.5E5
        DIF4 = 0.4E16
      ELSE IF (NLON.EQ.192 .AND.NOREC.EQ.96) THEN       ! T63
        DTIME = 900.
        GRAVFAC = 2.
        DIF2 = 0.35E5
        DIF4 = 0.2E16
      ELSE IF (NLON.EQ.320 .AND.NOREC.EQ.160) THEN      ! T106
        DTIME = 450. ! 600.
        GRAVFAC = 2.
        DIF2 = 0.25E5
        DIF4 = 0.1E16
      ELSE
        WRITE (NOUT,1300) NLON, NOREC
 1300   FORMAT(/'*** Program Error: unknown AGCM horizontal resolution'
     *         /'    NLON=',I6,'   NOREC=',I6
     *         /'    Stop in DATGET setting DIF2,DIF4.')
        CALL ENDRUN (-1)
      ENDIF
c77777GWROUGH = -999.      !-ve = data file + abs for land ice
      GWROUGH = -200.      !-ve = data file + abs for land ice
      ENVFAC  = -999.      !-ve = Himal only, +ve = all, if data file 
      DELEVTYPE = 0
      LEVRES = 0
      HORRES(1) = 0
      HORRES(2) = 0
      HORRESLSX(1) = 0
      HORRESLSX(2) = 0
 
      TIMRESERV = 3.       !time scale(hrs) of T,u,v reservoirs (reserv)
c v2,3:
      PLUMERAD1    = 200.  !pbl  plume radius
      PLUMEAREA1_L = .005  !pbl  plume fractional area over land
      PLUMEAREA1_O = 0.2   !pbl  plume fractional area over ocean
      PLUMERAD2    = 2000. !free plume radius
      PLUMEAREA2   = .01   !free plume fractional area
cexp
c     PLUMERAD1    = 150.  !pbl  plume radius
c     PLUMEAREA1_L = .0075 !pbl  plume fractional area over land
c     PLUMEAREA1_O = 0.3   !pbl  plume fractional area over ocean
c     PLUMERAD2    = 1500. !free plume radius
c     PLUMEAREA2   = .015  !free plume fractional area
 
      ALW     = .95   !strat. condense if q > alw*qsat(above sfc,stabco)
      BLW_C   = 1.    !time scale (hrs) for conv cloud evap. (in reevap)
      BLW_A   = 1.    !same for anvil
      BLW_S   = 1.    !same for stratus
      CLWFAC_C= .00006 !factor for cloud water content in clwcalc, conv
      CLWFAC_A= .00002 !factor for cloud water content in clwcalc, anvil
      CLWFAC_S= .00002 !factor for cloud water content in clwcalc, strat
      CLW_P   =.0010  !within-plume maximum lwd (kg/kg) (in plume)
      DLW_LL   = 6.0  !time scale(hrs)for cloud->orig.precip,land-liq
      DLW_OL   = 6.0  !same for ocn-liq
      DLW_LI   = 0.4  !same for land-ice
      DLW_OI   = 0.4  !same for ocn-ice
      DLWFAC = 1.
      REFFAC_L= 1.
      REFFAC_O= 1.
      REFFAC_I= 1.
      RADLIQ = 0.
      ELW     = 0.3   !time scale(hrs) of sfc-depos(lowest layer,reevap)
      FLW_C   = 0.3   !"full" conv cloud fraction (reevap)
      FLW_A   = 0.5   !same for anvil 
      FLW_S   = 1.0   !same for stratus
      GLW     = 0.3   !time scale (hrs) of accretion by 100 mm/dy precip
      TLIQA   = -15.  !temp (C) for cloud liq-ice (reevap,cldcmp,optdep)
      TLIQB   = 10.   !same as tliqa (transition: tliqa to tliqa+tliqb)
      SIGANVIL= .60   !min sigma level for anvil lwc source (in plume)
      FACCLWP_SOL = 1.6  !factor multiplying sol rad's clwp (radctl)
      FACCLWP_IR  = 1.6  !factor multiplying ir  rad's clwp (radctl)
      SSACLD(1) = .99120 ! cloud single scattering albedo, vis(radsolar)
      SSACLD(2) = .98240 ! cloud single scattering albedo, nir(radsolar)
      DAYLENCOR = 1.E20  !used only for Coriolis (default set below)
      DAYLENRAD = 1.E20  !used only with zencal (default set below)
      SSLON     = 0.     !used only with zencal if daylenrad=0

#ifdef isotrac
c     stratiform clouds:
      t_s_i  = -10.+273.16     ! Liq/sol transition temperature
      s_s_i  = 0.0035 ! 0.004  ! Supersat. function slope
      hi_s_i = 0.5    ! 0.     ! Rel. humidity linear function intercept
      hs_s_i = 1.              ! Rel. humidity linear function slope 
      r_s_i  = 1.              ! Mean drop radius factor 
c     convective clouds:
      t_c_i  = -20.+273.16     ! Liq/sol transition temperature 
      s_c_i  = 0.              ! Supersat. function slope
      hi_c_i = 0.7    ! 0.     ! Rel. humidity linear function intercept
      hs_c_i = 1.              ! Rel. humidity linear function slope 
      r_c_i  = 100.   ! 1.     ! Mean drop radius factor 
c     postcondensation processes:
      t_r_i  = 0.+273.16       ! Liq/sol transition temperature
      equ_i  = 1      ! 0      ! Liq-sol/vap re-equil. flag, 0 no re-eq
      pii_i  = 0      ! 80     ! (i,j) point output: Manaus, Brazil, T31
      pjj_i  = 0      ! 23     !   "      "      "      "
      pdn_i  = 24              ! Point output frequency (time steps)
      hrst_i = 0.              ! Reset reservoirs to value (0=no reset)
#endif
c
      WRITE(NOUT,2000)
 2000 FORMAT(/' Reading namelist input (group name INPUT)...'/)

c        Open namelist file. leave open for lsx namelist read
c        in datgelsx, where it will be closed.

      call findiu (iunamel, loc)
      call openf (iunamel, locname(loc), ' ', 0, 0)
C
CL        1.4        READ NAMELIST INPUT
C
      READ (IUNAMEL,INPUT)
      WRITE(NOUT,*)

c        Set data input filenames if entered 
c        (defaults for locname(...) have been set in data stmt above)

      IF (   TOPOG_FILE.NE.' ') LOCNAME(1)  =    TOPOG_FILE
      IF (   OZONE_FILE.NE.' ') LOCNAME(2)  =    OZONE_FILE
      IF (     GWD_FILE.NE.' ') LOCNAME(3)  =      GWD_FILE
      IF (  TOPOG2_FILE.NE.' ') LOCNAME(4)  =   TOPOG2_FILE
      IF (     VEG_FILE.NE.' ') LOCNAME(5)  =      VEG_FILE
      IF ( WISCVEG_FILE.NE.' ') LOCNAME(6)  =  WISCVEG_FILE
      IF (SURFTYPE_FILE.NE.' ') LOCNAME(7)  = SURFTYPE_FILE
      IF (  FWATER_FILE.NE.' ') LOCNAME(8)  =   FWATER_FILE
      IF ( OCEANUV_FILE.NE.' ') LOCNAME(9)  =  OCEANUV_FILE
      IF (     SST_FILE.NE.' ') LOCNAME(10) =      SST_FILE
      IF (  ISOT_O_FILE.NE.' ') LOCNAME(11) =   ISOT_O_FILE
      IF (  ISOT_S_FILE.NE.' ') LOCNAME(12) =   ISOT_S_FILE

      IF (SOILTEX01_FILE.NE.' ') LOCNAME(PRESOIL+1) = SOILTEX01_FILE
      IF (SOILTEX02_FILE.NE.' ') LOCNAME(PRESOIL+2) = SOILTEX02_FILE
      IF (SOILTEX03_FILE.NE.' ') LOCNAME(PRESOIL+3) = SOILTEX03_FILE
      IF (SOILTEX04_FILE.NE.' ') LOCNAME(PRESOIL+4) = SOILTEX04_FILE
      IF (SOILTEX05_FILE.NE.' ') LOCNAME(PRESOIL+5) = SOILTEX05_FILE
      IF (SOILTEX06_FILE.NE.' ') LOCNAME(PRESOIL+6) = SOILTEX06_FILE
      IF (SOILTEX07_FILE.NE.' ') LOCNAME(PRESOIL+7) = SOILTEX07_FILE
      IF (SOILTEX08_FILE.NE.' ') LOCNAME(PRESOIL+8) = SOILTEX08_FILE
      IF (SOILTEX09_FILE.NE.' ') LOCNAME(PRESOIL+9) = SOILTEX09_FILE
      IF (SOILTEX10_FILE.NE.' ') LOCNAME(PRESOIL+10)= SOILTEX10_FILE
      IF (SOILTEX11_FILE.NE.' ') LOCNAME(PRESOIL+11)= SOILTEX11_FILE
      IF (SOILTEX12_FILE.NE.' ') LOCNAME(PRESOIL+12)= SOILTEX12_FILE
C
CL        2.00       SET DEFAULTS FOR NAMELIST RADOTG
C
      SOLFAC = 1.
      SOLPAR = 1367.     ! W m-2 (was 1370 for Genesis v1.02, 1365 v2)

c     ipcc, 1991-1995, ppv:
      CO2PPV = 355.0e-6
      CH4PPV = 1.714e-6
      N2OPPV = 0.311e-6
      F11PPV = 0.280e-9
      F12PPV = 0.503e-9

      O3FAC    = 1.     ! factor multiplying o3 from data file (radoz1)
                        ! If < 0, average vs lat and month, * by |o3fac|
      O2PPM    = .23143 ! O2 mass mixing ratio
      AERTAUVIS= 0.14   ! backgrd aerosol col vis optical depth (aermix)
      MXAERL   = 3      ! # bottom layers aertauvis distrib over(aermix)
      GGUNIF   = .FALSE. !true for unif distrib ch4,n2o,cfc's in trcmix
      NOCLDTOP = 0      ! number of layers at top with no clds for rad
      NOCLDBOT = 0      ! number of layers at bot with no clds for rad
C
CL        2.04       READ NAMELIST RADOTG
C
      WRITE(NOUT,2002)
 2002 FORMAT(/' Reading namelist input (group name RADOTG)...'/)
      READ (IUNAMEL,RADOTG)
      WRITE(NOUT,*)

c     radiation can't have exactly zero trace gas amounts:
      if (co2ppv.eq.0.) co2ppv = 1.e-15
      if (ch4ppv.eq.0.) ch4ppv = 1.e-15
      if (n2oppv.eq.0.) n2oppv = 1.e-15
      if (f11ppv.eq.0.) f11ppv = 1.e-15
      if (f12ppv.eq.0.) f12ppv = 1.e-15
      if (o3fac.eq.0.)  o3fac  = 1.e-15
      if (o2ppm.eq.0.)  o2ppm  = 1.e-15
C
CL        FOR REST OF THIS SUBR, ACCUMULATE # OF NAMELIST FATAL ERRORS,
CL        AND DON'T CALL ENDRUN UNTIL THE END (SO USERS CAN SEE ALL
CL        ERROR MESSAGES AT ONCE)
C
      IER = 0
C
CL        2.1        CHECK THAT MUST-BE-ENTERED VARS HAVE BEEN ENTERED
C
#ifdef mss
      IF (MSPATH.EQ.'null') THEN
         WRITE(NOUT,2100)
 2100    FORMAT(' *** Must specify a Mass Store directory (MSPATH)',
     *          ' for saving history and restart files.')
         IER = IER + 1
      ELSE IF (MSPATH.EQ.' ') THEN
         WRITE(NOUT,2102)
 2102    FORMAT(/' *** Warning: Since entered MSPATH is blank, any',
     *                 ' history or restart files'
     *          /'     generated by this run will *NOT* be saved to',
     *                 ' the Mass Store !' /)
      ENDIF
#else
      IF (MSPATH.EQ.'null') MSPATH = ' '
#endif
C
      IF (NSREST.EQ.-9999) THEN
         WRITE(NOUT,*) '*** Must specify whether this run is a cold',
     *                 ' start, initial start or restart (NSREST).'
         NSREST = -1
         IER = IER + 1
      ELSE IF (NSREST.LT.-1 .OR. NSREST.GT.1) THEN
         WRITE(NOUT,*) '*** NSREST must be -1, 0 or 1.',
     *                 ' Entered value = ',NSREST,'.'
         NSREST = -1
         IER = IER + 1
      ENDIF
C
#ifdef isotrac
      IF (NSREST_ISOT.EQ.-9999) THEN
        WRITE(NOUT,*)'*** Must specify if isotopes are cold start',
     *               ' or restart file (NSREST_ISOT=-1,1).'
        NSREST_ISOT = -1
        IER = IER + 1
      ELSE
         IF (NSREST_ISOT.NE.-1 .AND. NSREST_ISOT.NE.1) THEN
           WRITE(NOUT,*) '*** NSREST_ISOT must be -1 or 1.',
     *                   ' Entered value = ',NSREST_ISOT,'.'
           NSREST_ISOT = -1
           IER = IER + 1
        else if (nsrest.eq.-1 .and. nsrest_isot.ne.-1) then
          write(nout,*) '*** Incompatible NSREST, NSREST_ISOT.',
     *                  ' Entered values = ',NSREST, NSREST_ISOT,'.'
          ier = ier + 1
        endif
      endif

      if (priestq) then
        write(nout,*) '*** Cannot have PRIESTQ(=.TRUE.) with an',
     *                ' isotopic program.'
        ier = ier + 1
      endif
#else
      if (nsrest_isot.eq.-9999) then
        nsrest_isot = -1
      else
        write(nout,*) '*** Cannot enter NSREST_ISOT with a',
     *                ' non-isotopic program.'
        ier = ier + 1
      endif
#endif
C
      IF (BASEDATE.EQ.' ') THEN
         IF (NSREST.LE.0) THEN
           WRITE(NOUT,*) '*** For NSREST = -1 or 0, a base date',
     *                   ' (BASEDATE) must be specified.'
           IER = IER + 1
         ENDIF
      ELSE
         CALL DATECONV2 (MBDATE, BASEDATE, 'BASEDATE')
      ENDIF
C
      IF (ENDDATE.EQ.' ') THEN
         WRITE(NOUT,*) '*** Must specify an end date for this run',
     *                 ' (ENDDATE).'
         IER = IER + 1
      ELSE
         CALL DATECONV2 (MEDATE, ENDDATE, 'ENDDATE')
      ENDIF
C
      IF (MSPATHIN.EQ.'null') MSPATHIN = MSPATH
C
      IF (NSREST.GE.0 .AND. MSNAMIN.EQ.' ') THEN
         WRITE(NOUT,*) '*** For NSREST = 0 or 1, an initial/restart',
     *                 ' file (MSNAMIN) must be specified.'
         IER = IER + 1
      ELSE IF (NSREST.EQ.-1 .AND. MSNAMIN.NE.' ') THEN
         WRITE(NOUT,*) '*** Warning: Initial/restart file MSNAMIM',
     *                 ' is ignored for cold start NSREST=-1.'
      ENDIF
C
CL        2.2        SET DEFAULTS FOR REGIONAL-MODEL DRIVER FILES
C
      IF (REGPATH.EQ.'null') REGPATH = MSPATH
C
      IF (REGHIS_A.EQ.' ' .AND. REGHIS.NE.' ')
     *   REGHIS_A = REGHIS(1:LENCHR(REGHIS)) // 'LSX'
C
CL        2.25       SET DEFAULT FOR STORM-TRACK HISTORY FILE PATH
C
      IF (STORMPATH.EQ.'null') STORMPATH = MSPATH
C
CL        2.25       SET DEFAULTS FOR IBIS RESTART FILE(S)
C
      IF (DYNAIBIS.EQ.1) THEN
        IF (IBISPATHIN.EQ.'null') IBISPATHIN = MSPATHIN
        IF (NSREST_IBIS.GE.0 .AND. IBISRESIN.EQ.' ') THEN
           WRITE(NOUT,*)
     *        '*** For NSREST_IBIS = 0 or 1, an initial/restart IBIS',
     *        ' file (IBISRESIN) must be specified.'
           IER = IER + 1
        ELSE IF (NSREST_IBIS.EQ.-1 .AND. IBISRESIN.NE.' ') THEN
           WRITE(NOUT,*)
     *        '*** Warning: initial/restart IBIS file IBISRESIN',
     *        ' is ignored for IBIS cold start (NSREST_IBIS=-1).'
        ENDIF
      ENDIF
C
CL        2.3        PROCESS FILENAME LISTS
C
      CALL DATNAMES
C
CL        2.4        PROCESS HISTORY-WRITE,SAVE,NEXT-FILE TIME LISTS
C
      CALL DATTIMES
C
CL        2.5        PROCESS HISTORY FIELD LIST
C
      CALL DATFIELD
C
CL        3.0        READ RESTART FILE (EXCEPT FOR COLD START)
C
      IF (NSREST.GE.0) THEN
         CALL RSDS (MSNAMIN, MSPATHIN, RUNTITLEX,MBDATEX,NSTEPX,DTIMEX)
      ELSE
         VERSRESIN = 0.
         VERSRESIN_ISOT = 0.
      ENDIF
C
CL        3.1        IF A TRUE RESTART, WARN IF INITIAL/RESTART FILE
CL                   DOESN'T MATCH CURRENT RESTART FILE LIST
C
      IF (NSREST.EQ.1) THEN
         IF (MSNAMIN.NE.MSNAMRES(IFILEPRE) .AND.
     *       MSNAMIN.NE.MSNAMRE2) WRITE(NOUT,3150)
     *      MSNAMIN, IFILEPRE, IFILEPRE, MSNAMRES(IFILEPRE)
 3150    FORMAT(
     *     /' *** Warning ***'
     *     /'     The entered restart filename does not agree with the'
     *     /'     expected position in the restart filename list:'
     *     /'     Entered restart filename MSNAMIN           = ',A
     *     /'     Expected list position                     = ',I8
     *     /'     Current filename there MSNAMRES (',I4,')     = ',A /)
C
         IF (IFILE.GT.NFILE) THEN
            WRITE(NOUT,*)'*** Not enough filenames entered',
     *                   ' to complete this run. Stop in DATGET.'
            IER = IER + 1
         ENDIF
      ENDIF
C
CL        3.2        SET VARIABLES THAT FOR TRUE RESTARTS CAN DEPEND
C                    ON RESTART FILE
C
      IF (NSREST.LE.0) THEN
         NSTEP = 1
      ELSE
         IF (RUNTITLE.EQ.' ') RUNTITLE = RUNTITLEX
         IF (BASEDATE.NE.' ') THEN
           IF (MBDATE.NE.MBDATEX) WRITE(NOUT,3200) MBDATE, MBDATEX
 3200      FORMAT(/' *** Warning: Ignoring entered base date',I12
     *            /'     since this is a restart run (NSREST=1).'
     *            /'     Using base date from restart file: ',I12/)
         ENDIF
         MBDATE = MBDATEX
         NSTEP = NINT (NSTEPX*DTIMEX/DTIME) + 1
      ENDIF
C
CL        3.5        SET MAIN RUN-TIME VARIABLES
C
      IF (IABS(MEDATE).GE.100000) THEN
        CALL DATEDIF (NDAYS, MBDATE, MEDATE)
        NSTOP = NINT (NDAYS*86400./DTIME)
      ELSE
        IF (MEDATE.GT.0) THEN
          NSTOP = NINT (MEDATE*86400./DTIME)
        ELSE
          NSTOP = -MEDATE
        ENDIF
        CALL NUTIME (MEDATE, MESEC, MBDATE, NSTOP, DTIME, 0)
      ENDIF
C
      TWODT = 2.0*DTIME
      NRSTRT = NSTEP
C
C        ESTIMATE SIZE OF EACH FILE IN WORDS NUSIZ (FOR OPENF).
C
C     iwpf = estimated # of history writes on each history file
C     (reduced if the length of run is less than 1 year)
      IWPF = (NHISTN/MAX(NNEXFN,1))
     *       * MIN (1., (NSTOP-NRSTRT+1)*DTIME/(86400.*365.))
      IF (NULDAT.GT.0) THEN
C       Add 1000 for an (over)estimate of one set of header records
        NUSIZ = (1000 + NURECL*NULDAT) * MAX(IWPF,1)
      ELSE
        NUSIZ = 0
      ENDIF
C
CL        4.0        OPEN HISTORY FILE
C
      IF (NSREST.LE.0) THEN
C
C        Open new history file for cold start or initial start
C
         IFILE = 1
         IFILEPRE = 1
         CALL HISSET(1)
         CALL OPENF (NDATA, MSNAMHIS(IFILE), MSPATH, 1, NUSIZ)
C
      ELSE
C
         IF (NUOFF.EQ.0) THEN
C
C           Open new history file for true restart.
C           (NUOFF is zero either if last timestep of previous run was
C           a next-file, or if no data records were written since this
C           file was started.)
C
            CALL HISSET(1)
            CALL OPENF (NDATA, MSNAMHIS(IFILE), MSPATH, 1, NUSIZ)
C
         ELSE
C
C           Open existing history file to append to for true restart.
C           (In this case, IFILEPRE must equal IFILE - see LINEMS).
C
            IF (MSNAMHISX.NE.MSNAMHIS(IFILE)) WRITE(NOUT,4000)
     *         IFILE, IFILE, MSNAMHIS(IFILE), MSNAMHISX
 4000       FORMAT(
     *        /' *** Warning ***'
     *        /'     The current history file to be appended to by',
     *               ' this run'
     *        /'     does not agree with its previously recorded name:'
     *        /'     Current position in history filename list = ',I8
     *        /'     Current history file MSNAMHIS (',I4,')      = ',A
     *        /'      Previous filename in that position        = ',A/)
C
            if (.not. hisnet) then
              CALL ACQUIREF (MSNAMHIS(IFILE), MSPATH)
              CALL OPENF (NDATA, MSNAMHIS(IFILE), MSPATH, 0, 0)
              ctmp = 'xxxxxxxx'
              read (ndata,err=4100,iostat=ios) ctmp
 4100         if (ctmp(5:8).eq.' ') ctmp(5:8) = 'CRAY'
              if (ctmp(5:8).ne.typehis) then
                write(nout,4102) ifile,msnamhis(ifile),ctmp(5:8),typehis
 4102           format(
     *            /' *** The current history file to be',
     *                   ' appended to by this run'
     *            /'     has a different machine numerical format',
     *                   ' than the current program:'
     *            /'     MSNAMHIS (',I4,')            = ',a
     *            /'     Machine format of file       = ',a
     *            /'     Machine format of program    = ',a/)
                ier = ier + 1
              endif
 
              IF (CTMP(1:4).NE.'CCM') THEN
                WRITE(NOUT,4104) IFILE, MSNAMHIS(IFILE), CTMP
 4104           FORMAT(
     *            /' *** The current history file to be',
     *                   ' appended to by this run'
     *            /'     has an incorrect first-4-byte label:'
     *            /'     MSNAMHIS (',I4,')            = ',A
     *            /'     First word (should be CCM) = ',A /)
                IER = IER + 1
              ENDIF
              CALL CLOSEF (NDATA)
            endif
C
            CALL OPENF (NDATA, MSNAMHIS(IFILE), MSPATH, -1, 0)
            CALL HISSET(-1)
C
         ENDIF
C
      ENDIF
C
CL        7.0        REMAINING CONSTRAINTS/UNITS CHANGES FOR USER INPUT
C
      REGFREQ   = REGFREQ  *3600.                         !hours to secs
      STORMFREQ = STORMFREQ*3600.                         !hours to secs
C
c     Ensure sw/lw time intervals are multiples of agcm timestep
      IRADLW  = MAX (1, NINT(TRADLW *3600./DTIME))        !hours to secs
      IRADLWA = MAX (1, NINT(TRADLWA*3600./DTIME))
      IRADSW  = MAX (1, NINT(TRADSW *3600./DTIME))
      TRADLW  = IRADLW *DTIME
      TRADLWA = IRADLWA*DTIME
      TRADSW  = IRADSW *DTIME
C
      BUDFREQ = BUDFREQ*86400.                            !days  to secs
      TIMOLSX = TIMOLSX*3600.                             !hours to secs
C
      STRATTIME = STRATTIME*86400.                        !days  to secs
C
      TIMRESERV = TIMRESERV*3600.                         !hours to secs
C
C        nb: daylencor only affects the Coriolis force via omega and ez
C        (setcom/init.f, spect.f), *not* the diurnal solar period.
C
      IF (DAYLENCOR.EQ.1.E20) THEN
        DAYLENCOR = 24.
      ENDIF
      DAYLENCOR  = DAYLENCOR *3600.                       !hours to secs
C
C        nb: daylenrad only affects solar radiation diurnal cycle
C        in radctl, call to zencal
C
      IF (DAYLENRAD.EQ.1.E20) THEN
        DAYLENRAD = 24.
      ENDIF
      DAYLENRAD = DAYLENRAD *3600.                        !hours to secs
C
      RETPD = MIN (32767., MAX (1., RETPD))
C
      IF (NSTOP.LT.NSTEP) THEN
         CALL NUTIME (MCDATE, MCSEC, MBDATE, NSTEP, DTIME, 1)
         WRITE(NOUT,7000) MEDATE, MCDATE
 7000    FORMAT( ' *** Specified end date for this run    :',I12
     *          /'     is earlier than its (re)start date :',I12 )
         IER = IER + 1
      ENDIF
C
      IF (ECCU.NE.-9999. .AND. (ECCU.LT.0..OR.ECCU.GE.1.)) THEN
         WRITE(NOUT,7200) ECCU
 7200    FORMAT(' *** Specified eccentricity ECCU is',f8.5,'.',
     *          ' Must be ge 0 and lt 1.')
         IER = IER + 1
      ENDIF
C
      IF (OBLU.NE.-9999. .AND. (OBLU.LT.0..OR.OBLU.GE.90.)) THEN
         WRITE(NOUT,7300) OBLU
 7300    FORMAT (' *** Specified obliquity OBLU is',f8.3,'.',
     *           ' Must be ge 0 and lt 90.')
         IER = IER + 1
      ENDIF
C
      IF (PRECU.NE.-9999. .AND. (PRECU.LT.0..OR.PRECU.GT.360.)) THEN
         WRITE(NOUT,7400) PRECU
 7400    FORMAT (' *** Specified precession PRECU is',f8.3,'.',
     *           ' Must be ge 0 and le 360.')
         IER = IER + 1
      ENDIF
C
      IF (HYBTYPE.LT.0 .OR. HYBTYPE.GT.2) THEN
         WRITE(NOUT,7500) HYBTYPE
 7500    FORMAT (' *** Specified hybrid-coord-type HYBTYPE is',i6,'.',
     *           ' Must be 0, 1 or 2.')
         IER = IER + 1
      ENDIF
C
      IF ( LEVRES.NE. 0 .AND. 
     *     LEVRES.NE. 6  .AND.LEVRES.NE.12 .AND. LEVRES.NE.16 .AND.
     *     LEVRES.NE.18 .AND. LEVRES.NE.20 .AND. LEVRES.NE.36 ) THEN
         WRITE(NOUT,7600) LEVRES
 7600    FORMAT(' *** Specified restart-file AGCM vertical',
     *                ' resolution LEVRES is',I6,'.'
     *         /'     Must be one of the preset values in subr',
     *                ' READERV: 6, 12, 16, 18, 20 or 36.')
         IER = IER + 1
      ENDIF
C
      IF ( HORRES(1).LT.0 .OR. HORRES(1).GT.360 .OR.
     *     HORRES(2).LT.0 .OR. HORRES(2).GT.180 ) THEN
         WRITE(NOUT,7700) HORRES(1), HORRES(2)
 7700    FORMAT(' *** Specified restart-file AGCM horizontal',
     *                ' resolutions HORRES are (',I8,',',I8,').'
     *         /'     Maximum allowable resolutions are (360,180).')
         IER = IER + 1
      ENDIF
C
      IF ( HORRESLSX(1).LT.0 .OR. HORRESLSX(1).GT.360 .OR.
     *     HORRESLSX(2).LT.0 .OR. HORRESLSX(2).GT.180 ) THEN
         WRITE(NOUT,7800) HORRESLSX(1), HORRESLSX(2)
 7800    FORMAT(' *** Specified restart-file LSX horizontal',
     *                ' resolutions HORRESLSX are (',I8,',',I8,').'
     *         /'     Maximum allowable resolutions are (360,180).')
         IER = IER + 1
      ENDIF
C
      IF (OCEANTYPE.LT.0 .OR. OCEANTYPE.GT.2) THEN
         WRITE(NOUT,7900) OCEANTYPE
 7900    FORMAT (' *** Specified sst/slab/ogcm flag OCEANTYPE',
     *           ' is',i6,'. Must be 0, 1 or 2.')
         IER = IER + 1
      ENDIF
C
      IF (OCEANTYPE.EQ.0) THEN           ! prescribed SSTs
         DYNAMICE = .FALSE.
      ENDIF
C
      IF (VEGTYPE.LT.0 .OR. VEGTYPE.GT.7) THEN
         WRITE(NOUT,7910) VEGTYPE
 7910    FORMAT (' *** Specified vegetation flag VEGTYPE',
     *           ' is',i6,'. Must be 0, 1, 2, 3, 4, 5, 6 or 7.')
         IER = IER + 1
      ENDIF
C
      IF (EVEDELTA.NE.0 .AND. VEGTYPE.NE.3) THEN
         WRITE(NOUT,7912) 
 7912    FORMAT (' *** EVE deltas (EVEDELTA=1) are meaningful only',
     *           ' with interactive EVE (VEGTYPE=3).')
         IER = IER + 1
      ENDIF
C
      IF (DSBIOME.LT.0 .OR. DSBIOME.GT.12) THEN
         WRITE(NOUT,7915) DSBIOME
 7915    FORMAT (' *** Specified Dorman-Sellers biome DSBIOME',
     *           ' is',i6,'. Must be in range 1-12.')
         IER = IER + 1
      ENDIF
C
      IF (WISCVEG.NE.0 .AND. (WISCVEG.LT.1100.OR.WISCVEG.GT.7399)) THEN
         WRITE(NOUT,7917) WISCVEG
 7917    FORMAT (' *** Specified Wisconsin vegetation category WISCVEG',
     *           ' is',i8,'. Must be in range 1100-7399.')
         IER = IER + 1
      ENDIF
C
      IF (PHOTOTYPE.LT.0 .OR. PHOTOTYPE.GT.1) THEN
         WRITE(NOUT,7920) PHOTOTYPE
 7920    FORMAT (' *** Specified photosynthesis flag PHOTOTYPE',
     *           ' is',i6,'. Must be 0 or 1.')
         IER = IER + 1
      ENDIF
C
      IF (DYNAIBIS.LT.0 .OR. DYNAIBIS.GT.1) THEN
         WRITE(NOUT,7922) DYNAIBIS
 7922    FORMAT (' *** Specified dynamic IBIS vegetation flag DYNAIBIS',
     *           ' is',i6,'. Must be 0 or 1.')
         IER = IER + 1
      ENDIF
C
#ifndef ibis
      IF (PHOTOTYPE.EQ.1 .OR. VEGTYPE.EQ.1) THEN
         WRITE(NOUT,7923)
 7923    FORMAT (' *** This program has been compiled without',
     *           ' IBIS photosynthesis (-D ibis),'
     *          /'     which is needed for PHOTOTYPE=1 or VEGTYPE=1.')
         IER = IER + 1
      ENDIF
#endif
C
#ifndef ibisint
      IF (DYNAIBIS.EQ.1) THEN
         WRITE(NOUT,7924)
 7924    FORMAT (' *** This program has been compiled without',
     *           ' dynamic IBIS vegetation (-D ibisint),'
     *          /'     which is needed for DYNAIBIS=1.')
         IER = IER + 1
      ENDIF
#endif
C
      IF (DYNAIBIS.EQ.1 .AND..NOT. (VEGTYPE.EQ.1 .AND. PHOTOTYPE.EQ.1))
     *   THEN
         WRITE(NOUT,7925) VEGTYPE, PHOTOTYPE
 7925    FORMAT (' *** For dynamic IBIS vegetation (DYNAIBIS=1),',
     *           ' VEGTYPE and PHOTOTYPE must be 1.'
     *           /'    VEGTYPE =',i3,'    PHOTOTYPE=',i3)
         IER = IER + 1
      ENDIF
C
      IF (SOILTEX.NE.0) THEN
        IZSAND = SOILTEX/10000
        IZSILT = (SOILTEX - 10000*IZSAND) / 100
        IZCLAY =  SOILTEX - 10000*IZSAND - 100*IZSILT
        IF (IZSAND+IZSILT+IZCLAY.NE.100) THEN
          WRITE(NOUT,7926) SOILTEX
 7926     FORMAT (' *** Specified soil sand-silt-clay texture SOILTEX',
     *            ' is',i8
     *           /'     Must unpack to sand+silt+clay = 100.')
          IER = IER + 1
        ENDIF
      ENDIF
C
#ifndef NETCDF
      IF (HISNET) THEN
         WRITE(NOUT,7960)
 7960    FORMAT (' *** Program has been compiled without -D netcdf,',
     *           ' so HISNET cannot be .TRUE.')
         IER = IER + 1
      ENDIF
#endif
C
      IF (IER.GT.0) THEN
        CALL ENDRUN (-1)
      ENDIF
C
      RETURN
      END
C
C ***************************** DATNAMES ******************************
C
      SUBROUTINE DATNAMES
C----------------------------------------------------------------------
C
CL            PROCESS USER-ENTERED FILENAME LISTS
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      CHARACTER CTMP*4
C----------------------------------------------------------------------
C
C
CL        1.0        POSSIBLE "*" EXPANSION FOR AGCM HISTORY FILE LIST
C
      CALL STARLIST (MSNAMHIS, PFILE, 'AGCM history')
C
CL        1.1        COUNT NUMBER OF AGCM HISTORY FILES NFILE
C
      NFILE = 0
      DO 110 I=PFILE,1,-1
         IF(MSNAMHIS(I).NE.' ') THEN
            NFILE = I
            GOTO 112
         ENDIF
  110 CONTINUE
  112 CONTINUE
C
CL        1.2        SET DEFAULT CCM HISTORY FILENAMES IF NONE ENTERED,
CL                   OR CHECK FOR ANY ENTERED EMBEDDED BLANKS
C
      IF (NFILE.EQ.0) THEN
        DO 120 I=1,PFILE
          WRITE(CTMP,'(I4)') 1000+I
          MSNAMHIS(I) = 'HIST' // CTMP(2:4)
  120   CONTINUE
        NFILE = PFILE
      ELSE
        DO 125 I=1,NFILE
          IF (MSNAMHIS(I).EQ.' ') THEN
            WRITE(NOUT,*) '*** Embedded blank(s) in agcm history file',
     *                    ' list MSNAMHIS. First occurrence at # ',I,'.'
            CALL ENDRUN (-1)
          ENDIF
  125   CONTINUE
      ENDIF
C
CL        1.3        POSSIBLE "*" EXPANSION FOR RESTART FILE LIST
C
      CALL STARLIST (MSNAMRES, PFILE, 'Restart')
C
CL        1.4        SET DEFAULT RESTART FILENAMES, DERIVED FROM
CL                   CCM HISTORY FILENAMES.
C
      DO 140 I=1,NFILE
        IF (MSNAMRES(I).EQ.' ') THEN
          DO 142 L=LEN(MSNAMHIS(I))-3,1,-1
            LENH = L
            IF (MSNAMHIS(I)(L:L).NE.' ') GOTO 144
  142     CONTINUE
  144     MSNAMRES(I) = MSNAMHIS(I)(1:LENH)//'RES'
        ENDIF
  140 CONTINUE
C
      RETURN
      END
C
C ***************************** DATTIMES ******************************
C
      SUBROUTINE DATTIMES
C----------------------------------------------------------------------
C
CL            PROCESS USER-ENTERED HISTORY-WRITE, SAVE, NEXT-FILE TIMES
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      DIMENSION NDAYDEF(12)
      DATA NDAYDEF /  31,  59,  90, 120, 151, 181,
     *               212, 243, 273, 304, 334, 365  /
C----------------------------------------------------------------------
C
      NTPY = NINT(365*86400./DTIME)
C
C--------------------
C HISTORY WRITE TIMES
C--------------------
C
CL        1.0       SET DEFAULT HISTORY-WRITE END TIMES
C
      IF (NHISI(1).EQ.0) THEN
         DO 100 I=1,12
            NHISI(I) = NDAYDEF(I)
  100    CONTINUE
      ENDIF
C
CL        1.1       COUNT NUMBER OF HISTORY-WRITE PERIODS NHISTN, AND
CL                  CHANGE HIST-WRITE END TIMES TO UNITS OF TIMESTEPS.
C
      NHISTN = 0
      KOUNT = 0
      DO 110 I=1,NMARK
         IF (NHISI(I).EQ.0) GOTO 112
C        Convert any "mm/dd" input to day #
         CALL DATECONV (NHISI(I), I, 'NHISI', KOUNT)
         NHISTN = I
         IF (NHISI(I).GT.0) THEN
            NHISI(I) = MAX ( NINT(NHISI(I)*86400./DTIME), 1 )
         ELSE IF (NHISI(I).LT.0) THEN
            NHISI(I) = -NHISI(I)
         ENDIF
  110 CONTINUE
  112 CONTINUE
C
CL       1.14       IF JUST ONE (NON-STRING) WRITE TIME HAS BEEN
CL                  ENTERED, CONVERT TO MULTIPLES OF THAT TIME FROM
CL                  BEGINNING OF YEAR
C
      IF (NHISTN.EQ.1 .AND. KOUNT.EQ.0) THEN
         NHISTN = NTPY / NHISI(1)
         NHISTN = MIN (NMARK, MAX (1, NHISTN))
         IF (NHISTN.GT.1) THEN
            DO 114 I=2,NHISTN
               NHISI(I) = NHISI(1)*I
  114       CONTINUE
         ENDIF
      ENDIF
C
CL       1.2        SET END TIMES NHIST(*,2) DEPENDING ON NHISI,
CL                  CONSTRAINING THEM TO BE WITHIN 1-YEAR SPAN
C
      DO 120 I=1,NHISTN
         NHIST(I,2)= MAX (1, MIN (NTPY, NHISI(I)))
  120 CONTINUE
C
CL       1.3        SET START TIMES NHIST(*,1) DEPENDING ON NDHIS.
CL                  IF NDHIS = 0 (DEFAULT), START TIME = PREVIOUS END
CL                  TIME + 1 TIMESTEP, IE, ACCUMULATED FIELDS WILL BE
CL                  ACCUMULATING 100% OF THE TIME. FIRST CONVERT
CL                  USER-ENTERED NDHIS TO UNITS OF TIMESTEPS.
C
      IF (NDHIS.GT.0) THEN
         NDHIS = NINT (NDHIS*86400./DTIME)
      ELSE
         NDHIS = -NDHIS
      ENDIF
C
      IF (NDHIS.EQ.0) THEN
         NHIST(1,1) = NHIST(NHISTN,2) + 1
         DO 130 I=2,NHISTN
            NHIST(I,1) = NHIST(I-1,2) + 1
  130    CONTINUE
      ELSE
         DO 132 I=1,NHISTN
            NHIST(I,1) = NHIST(I,2) - NDHIS + 1
  132    CONTINUE
      ENDIF
 
c        Constrain start times to not overlap into previous periods
 
      IF (NHIST(1,1).LE.0) NHIST(1,1) = NHIST(1,1) + NTPY
      IF (NHIST(1,1).GT.NHIST(1,2))
     *   NHIST(1,1) = MAX ( NHIST(NHISTN,2) + 1, NHIST(1,1) )
      DO 134 I=2,NHISTN
         NHIST(I,1) = MAX ( NHIST(I-1,2) + 1, NHIST(I,1) )
  134 CONTINUE
C
CL       1.4        CHECK REASONABLENESS OF HISTORY PERIODS
C
      DO 140 I=2,NHISTN
        IF (NHIST(I,1).GT.NHIST(I,2)) THEN
          WRITE(NOUT,1400) I-1, NHISI(I-1)*DTIME/86400.,
     *                     I,   NHISI(I)  *DTIME/86400.
 1400     FORMAT(/' *** AGCM history-write times NHISI must increase',
     *                  ' in time,'
     *           /'     and fall between 00:00 Jan 1 and 24:00 Dec 31:'
     *           /'     NHISI(',I3,') =',F8.3,'   NHISI(',I3,') =',F8.3)
          CALL ENDRUN (-1)
        ENDIF
  140 CONTINUE
C
C-----------
C SAVE TIMES
C-----------
C
      IF (NSAVE(1).EQ.0) THEN
C
CL        2.0        SET DEFAULT SAVE TIMES
C
 
         DO 200 I=1,NHISTN
            NSAVE(I) = NHIST(I,2)
  200    CONTINUE
         NSAVEN = NHISTN
C
      ELSE
C
CL        2.1        COUNT NUMBER OF SAVE TIMES NSAVEN, AND CHANGE
CL                   TO UNITS OF TIMESTEPS.
C
        NSAVEN = 0
        KOUNT = 0
        DO 210 I=1,NMARK
           IF (NSAVE(I).EQ.0) GOTO 212
C          Convert any "mm/dd" input to day #
           CALL DATECONV (NSAVE(I), I, 'NSAVE', KOUNT)
           NSAVEN = I
           IF (NSAVE(I).GT.0) THEN
              NSAVE(I) = MAX ( NINT(NSAVE(I)*86400./DTIME), 1 )
           ELSE IF (NSAVE(I).LT.0) THEN
              NSAVE(I) = -NSAVE(I)
           ENDIF
  210   CONTINUE
  212   CONTINUE
C
CL       2.2        IF JUST ONE (NON-STRING) SAVE TIME HAS BEEN
CL                  ENTERED, CONVERT TO MULTIPLES OF THAT TIME FROM
CL                  BEGINNING OF YEAR
C
        IF (NSAVEN.EQ.1 .AND. KOUNT.EQ.0) THEN
           NSAVEN = NTPY / NSAVE(1)
           NSAVEN = MIN (NMARK, MAX (1, NSAVEN))
           IF (NSAVEN.GT.1) THEN
              DO 224 I=2,NSAVEN
                 NSAVE(I) = NSAVE(1)*I
  224         CONTINUE
           ENDIF
        ENDIF
C
      ENDIF
C
CL       2.3        CONSTRAIN SAVE TIMES TO FALL WITHIN 1-YEAR SPAN
C
      DO 230 I=1,NSAVEN
        NSAVE(I) = MIN (NTPY, NSAVE(I))
  230 CONTINUE
C
C----------------
C NEXT-FILE TIMES
C----------------
C
CL        3.0        SET DEFAULT NEXT-FILE TIME(S)
C
      IF (NNEXF(1).EQ.0) THEN
         NNEXF(1) = 365
         NNEXF(2) = 0
      ENDIF
C
CL        3.1        COUNT NUMBER OF NEXT-FILE TIMES NNEXFN, AND CHANGE
CL                   TO UNITS OF TIMESTEPS.
C
      NNEXFN = 0
      KOUNT = 0
      DO 310 I=1,NMARK
         IF (NNEXF(I).EQ.0) GOTO 312
C        Convert any "mm/dd" input to day #
         CALL DATECONV (NNEXF(I), I, 'NNEXF', KOUNT)
         NNEXFN = I
         IF (NNEXF(I).GT.0) THEN
            NNEXF(I) = MAX ( NINT(NNEXF(I)*86400./DTIME), 1 )
         ELSE IF (NNEXF(I).LT.0) THEN
            NNEXF(I) = -NNEXF(I)
         ENDIF
  310 CONTINUE
  312 CONTINUE
C
CL       3.2        IF JUST ONE (NON-STRING) BEXT-FILE TIME HAS BEEN
CL                  ENTERED, CONVERT TO MULTIPLES OF THAT TIME FROM
CL                  BEGINNING OF YEAR
C
      IF (NNEXFN.EQ.1 .AND. KOUNT.EQ.0) THEN
         NNEXFN = NTPY / NNEXF(1)
         NNEXFN = MIN (NMARK, MAX (1, NNEXFN))
         IF (NNEXFN.GT.1) THEN
            DO 324 I=2,NNEXFN
               NNEXF(I) = NNEXF(1)*I
  324       CONTINUE
         ENDIF
      ENDIF
C
C
CL       3.3        CONSTRAIN NEXT-FILE TIMES TO FALL WITHIN 1-YEAR SPAN
C
      DO 330 I=1,NNEXFN
        NNEXF(I) = MIN (NTPY, NNEXF(I))
  330 CONTINUE
C
      RETURN
      END
C
C ***************************** DATFIELD ******************************
C
      SUBROUTINE DATFIELD
C----------------------------------------------------------------------
C
CL            PROCESS USER-ENTERED HISTORY FIELD LIST
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
C
CL        1.0        INITIALIZE HARD-CODED PROGRAM'S HISTORY FIELD LIST
C
      CALL FLDLST
C
CL        2.0        COUNT NUMBER OF HIST FIELDS ENTERED, AND TRANSFER
CL                   FROM INPUT ARRAY CUSHNAM TO CUSHIS, MUSHFL.
CL                   IF NO FIELDS ENTERED, RELY ON "NO-TRIP" DO-LOOPS
CL                   BELOW AND IN WRITH FOR "FIELDS" RECORD, AND JUST
CL                   SET NULDAT=0 BELOW FOR THAT RECORD. THEN PROGRAM
CL                   WILL RUN, WRITING HISTORY HEADERS BUT NO DATA RECS.
C
      NUSHIS = 0
      DO 200 I=1,NUSHMX
         IF (CUSHNAM(I).EQ.' ') GOTO 202
         CALL UPPER (CUSHNAM(I))
         NUSHIS = I
         N = LENCHR(CUSHNAM(I))
         IF (CUSHNAM(I)(N:N).EQ.'*') THEN
           MUSHFL(I) = 0                              ! instantaneous
           CUSHIS(I) = CUSHNAM(I)(1:N-1)
         ELSE
           MUSHFL(I) = 1                              ! time-accumulated
           CUSHIS(I) = CUSHNAM(I)
         ENDIF
  200 CONTINUE
  202 CONTINUE
C
CL        3.0        FORM INDEX FROM USER-ENTERED LIST TO PROGRAM LIST,
CL                   AND CHECK FOR DUPLICATE ENTRIES
C
      IER = 0
C
      DO 300 I=1,NUSHIS
         DO 302 IP = 1,NINDEX
            IF(CUSHIS(I).EQ.CINDEX(1,IP)) THEN
               MUSHIS(I) = IP
               GOTO 300
            ENDIF
  302    CONTINUE
         IER = IER + 1
         WRITE(NOUT,*) '*** Requested history field ',CUSHIS(I),
     *                 ' is not a valid field label.'
  300 CONTINUE
C
      DO 310 I=1,NUSHIS-1
         DO 312 J=I+1,NUSHIS
            IF (CUSHIS(I).EQ.CUSHIS(J)) THEN
               IER = IER + 1
               WRITE(NOUT,*) '*** Duplicate history fields: ',
     *                       CUSHIS(I),' at entry numbers ',I,' and ',J
            ENDIF
  312    CONTINUE
  310 CONTINUE
C
      IF (IER.GT.0) CALL ENDRUN (-1)
C
CL        4.0        SET "RECORD" POSITIONS IN AGCM HISTORY FILE FOR
CL                   EACH FIELD (RELATIVE TO END OF CURRENT HEADER RECS)
C
      DO 400 I = 1,NUSHIS
         IP = MUSHIS(I)
         IF(I.EQ.1) THEN
            MINDEX(3,IP) = 1
         ELSE
            MINDEX(3,IP) = MINDEX(4,MUSHIS(I-1)) + 1
         ENDIF
         MINDEX(4,IP) = MINDEX(3,IP) + MINDEX(2,IP) - 1
         MINDFL(IP) = MUSHFL(I)
  400 CONTINUE
C
CL        4.5        SET NUMBER OF HEADER AND DATA RECORDS FOR EACH
CL                   HISTORY WRITE
C
      IF (NUSHIS.GT.0) THEN
         NULHDR = 4
         NULDAT = MINDEX(4,MUSHIS(NUSHIS))
      ELSE
         NULHDR = 0
         NULDAT = 0
      ENDIF
C
CL        5.0        CONSTRAIN DATA PACKING DENSITY NDENS, SET
CL                   HIST-FILE RECORD LENGTH NURECL (DETERMINED BY
CL                   LENGTH OF DATA RECORDS, = 1 LABEL WORD +
CL                   UNPACKED/PACKED HORIZONTAL GRIDSIZE).
C
      IZ = NDENS
      NDENS = MAX (1, MIN(4,NDENS))
#ifdef nopack
      NDENS = 1
#endif
      IF (NDENS.NE.IZ) WRITE(NOUT,500) IZ, NDENS
  500 FORMAT(/' *** Warning: Entered history-file packing density',
     *       /'     NDENS =',I2,' has been reset to NDENS =',I2/)
 
      IF (NDENS.EQ.1) THEN
        NURECL = 1 + NLON*NOREC
      ELSE
        NURECL = 1 + (NLON*NOREC + NDENS - 1)/NDENS + 2
      ENDIF
C
CL        5.1        CONSTRAIN DATA PACKING DENSITY NDENS_STORM
C
      IZ = NDENS_STORM
      NDENS_STORM = MAX (1, MIN(4,NDENS_STORM))
#ifdef nopack
      NDENS_STORM = 1
#endif
      IF (NDENS_STORM.NE.IZ) WRITE(NOUT,510) IZ, NDENS_STORM
  510 FORMAT(
     *  /' *** Warning: Entered storm-track history-file packing',
     *   ' density'
     *  /'     NDENS_STORM =',I2,' has been reset to NDENS_STORM =',I2/)
C
      RETURN
      END
C
C ***************************** SETSIG ********************************
C
      SUBROUTINE SETSIG
C----------------------------------------------------------------------
C
CL            INITIALIZE VERTICAL SIGMA GRID (SIG, SIGKMH)
CL            (CALLED BEFORE DATGET SO READERV CAN USE MODEL'S SIG)
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
      DIMENSION SIG6(6),
     *          SIG12(12),
     *          SIG16(16), SIGL16(16+1),
     *          SIG18(18), SIGL18(18+1),
     *          SIG20(20),
     *          SIG36(36), SIGL36(36+1)
C----------------------------------------------------------------------
C     If change vert resol (NLEV in commons,paramnew(slt)), don't need
C     to change anything else...other vertical resolution issues are
C     taken care of automatically: klo,khi in cloud.f.
C     *** Except should add any additional vertical grids to readerv
c     *** in util.f.
C----------------------------------------------------------------------
C      6 layers:
      DATA SIG6  /
     *  .020, .100, .250, .500, .850, .991 /
C----------------------------------------------------------------------
C     12 layers:
      DATA SIG12 /
     *  .009, .025, .060, .110, .165, .245,
     *  .355, .500, .664, .811, .926, .991 /
C----------------------------------------------------------------------
C     16 layers (Martin Wild):
C     DATA SIG16 /
C    *  .010, .030, .060, .110, .175, .255, .350, .450,
C    *  .550, .650, .745, .830, .900, .950, .980, .995 /
C     sigkmh for Martin Wild 16-layer (different from equispaced)
C     DATA SIGL16 /
C    *  .003, .020, .040, .080, .140, .210, .300, .400,
C    *  .500, .600, .700, .790, .870, .930, .970, .990, 1.000 /
C----------------------------------------------------------------------
C     16 layers (med-res stratosphere):
      DATA SIG16 /
     *  .010, .030, .060, .090, .120,
     *  .150, .180, .210, .250, .290,
     *  .355, .500, .664, .811, .926, .991 /
C----------------------------------------------------------------------
C     18 layers (ccm2):
      DATA SIG18 /
     *  .005, .013, .033, .064, .099, .139,
     *  .189, .251, .325, .409, .501, .598,
     *  .695, .787, .866, .929, .970, .993 /
      DATA SIGL18 /
     *  .003, .008, .022, .049, .083, .118, ! top not at 0 for radiation
     *  .163, .219, .288, .367, .456, .551,
     *  .649, .744, .831, .903, .956, .985, 1.000 /
C----------------------------------------------------------------------
C     20 layers (hi-res stratosphere):
      DATA SIG20 /
     *  .005, .010, .020, .040, .060, .080, .100,
     *  .120, .140, .160, .180, .200, .230, .280,
     *  .355, .500, .664, .811, .926, .991 /
C----------------------------------------------------------------------
C     30 layers (very-hi-res stratosphere):
C     DATA SIG30 /
C    *  .003, .006, .010, .015, .020, .030,
C    *  .040, .050, .060, .070, .080, .090,
C    *  .100, .110, .120, .130, .140, .150,
C    *  .160, .170, .180, .200, .230, .280,
C    *  .355, .500, .664, .811, .926, .991 /
C----------------------------------------------------------------------
C     36 layers (doubled sig18):
      data sig36 /
     *    .0050, .0090, .0130, .0230, .0330, .0485,
     *    .0640, .0815, .0990, .1190, .1390, .1640,
     *    .1890, .2200, .2510, .2880, .3250, .3670,
     *    .4090, .4550, .5010, .5495, .5980, .6465,
     *    .6950, .7410, .7870, .8265, .8660, .8975,
     *    .9290, .9495, .9700, .9815, .9930, .9965 /

      data sigl36 /
     *    .0030, .0070, .0110, .0180, .0280, .0408,
     *    .0563, .0728, .0903, .1090, .1290, .1515,
     *    .1765, .2045, .2355, .2695, .3065, .3460,
     *    .3880, .4320, .4780, .5252, .5738, .6222,
     *    .6708, .7180, .7640, .8067, .8462, .8817,
     *    .9133, .9392, .9598, .9758, .9872, .9948, 1.0000 /
C----------------------------------------------------------------------
C
CL            1.0       SET VERTICAL-GRID LAYER MID-POINTS
C
      IF (NLEV.EQ.6) THEN
        CALL SCOPY (NLEV, SIG6,  1, SIG, 1)
      ELSE IF (NLEV.EQ.12) THEN
        CALL SCOPY (NLEV, SIG12, 1, SIG, 1)
      ELSE IF (NLEV.EQ.16) THEN
        CALL SCOPY (NLEV, SIG16, 1, SIG, 1)
      ELSE IF (NLEV.EQ.18) THEN
        CALL SCOPY (NLEV, SIG18, 1, SIG, 1)
      ELSE IF (NLEV.EQ.20) THEN
        CALL SCOPY (NLEV, SIG20, 1, SIG, 1)
      ELSE IF (NLEV.EQ.36) THEN
        CALL SCOPY (NLEV, SIG36, 1, SIG, 1)
      ELSE
        WRITE(NOUT,1000) NLEV
 1000   FORMAT(/' *** Error: AGCM vertical resolution NLEV = ',I6,
     *         /'     Must be one of the pre-set values in subr SETSIG.'
     *         /)
        CALL ENDRUN (-1)
      ENDIF
C
CL            1.1       SET VERTICAL-GRID LAYER BOUNDARIES (EQUISPACED
CL                      BETWEEN MID-POINTS, EXCEPT FOR NLEV=18)
C
      IF (NLEV.EQ.18) THEN
        CALL SCOPY (NLEVP, SIGL18, 1, SIGKMH, 1)
      ELSE IF (NLEV.EQ.36) THEN
        CALL SCOPY (NLEVP, SIGL36, 1, SIGKMH, 1)
      ELSE
        DO 1010 JK=2,NLEV
          SIGKMH(JK) = 0.5*(SIG(JK-1)+SIG(JK))
 1010   CONTINUE
        SIGKMH(1) = min (.001, 0.5*sig(1))
        SIGKMH(NLEVP) = 1.
      ENDIF
 
      RETURN
      END
C
C ***************************** SETCOM ********************************
C
      SUBROUTINE SETCOM
C----------------------------------------------------------------------
C
CL            INITIALIZE COMCON AND COMMAP
CL            (EXCEPT FOR SIG AND SIGKMH, ALREADY DONE IN SETSIG)
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMUL.h>
C----------------------------------------------------------------------
      DIMENSION ZALS(NLEV)
      PARAMETER (PI = 3.14159265358979)
C----------------------------------------------------------------------
C
CL            2.1       INITIALIZE COMCON
C
      CALUNI = 4.1868
      CARMAN = 0.35
      CPAIR  = 1004.64
C
      GRAVIT = 9.80616
      LATVAP = 2.5104E06
      LATICE = 3.336E5
      LATSUB = LATVAP + LATICE
C
C     OMEGA  = 7.292E-5                                     ! ccm1 value
C     OMEGA  =    3.14159265358979 /  43082.                ! ogcm value
      IF (DAYLENCOR.EQ.0.) THEN
        OMEGA = 0.
      ELSE IF (NINT(DAYLENCOR).EQ.86400) THEN
        OMEGA = 2.*3.14159265358979 / (DAYLENCOR*365.25/366.25)!inertial
      ELSE 
        OMEGA = 2.*3.14159265358979 / DAYLENCOR
      ENDIF
      RAIR   = 287.04
      REARTH = 6.37122E+6
      RH2O   = 4.61E2
      RHOH2O = 1.E3
      CH2O = 4.218E3
      CICE = 2.106E3
C
      STEBOL = 5.66961E-8
      TMELT  = 273.16
C
      CAPPA  = RAIR/CPAIR
      EZ     = OMEGA/SQRT(0.375)
      RA     = 1./REARTH
      RGA    = 1./GRAVIT
C
      ZVIR   = RH2O/RAIR -1.
      CPWV   = 1.81E3
      CPVIR  = CPWV/CPAIR - 1.

#ifdef isotrac
      dr_o   = 1./0.9723
      dr_h   = 1./0.9755
      smoj_o = 2.2149e-3
      smoj_h = 3.2883e-4
c     smoj_o = 1.9934e-3 * (20./18.)
c     smoj_h = 0.1557e-3 * (19./18.) * 2.
#endif
C
CL            2.3       INITIALISE COMMAP.
C
CL            2.31      VERTICAL LEVEL DEPENDENT ARRAYS
C
      DO 2310 JK=1,NLEV
        DSIGMA(JK) = SIGKMH(JK+1)- SIGKMH(JK)
        SIGKPK(JK) = SIG(JK)**CAPPA
        SIGKPL(JK) = SIGKMH(JK)**CAPPA
        R2DSIG(JK) = 0.5/DSIGMA(JK)
 2310 CONTINUE
        SIGKPL(NLEVP) = 1.
C
CL            2.32      SPECTRAL ARRAYS.
C
      ZRA2=RA*RA
      DO 2320 J=2,NPMAX
         SQ(J)=J*(J-1)*ZRA2
         RSQ(J)=1./SQ(J)
 2320 CONTINUE
      SQ(1)=0.
      RSQ(1)=0.
C
      DO 2330 J=1,MMAX
         XM(J)=J-1
 2330 CONTINUE
C
CL            2.34      GAUSSIAN LATITUDE DEPENDENT ARRAYS,
CL                      AND LONGITUDE GRID ARRAY.
C
      CALL GAUAW (SIT, W, NOREC)
C
      DO 2340 J=1,NOREC
         IF (SIT(J).GE.0.) THEN
            ALATG(J) = ASIN(SIT(J))
         ELSE
            ALATG(J) = -ASIN(-SIT(J))
         ENDIF
         CS(J)=1.-SIT(J)*SIT(J)
         RCS(J)=1./CS(J)
 2340 CONTINUE
C
      DO 2345 I=1,NLON
c       agcm longitude #1 has edge at Greenwich, as in lsx
c       (old ccm1 longitude #1 had center at Greenwich)
        ALONG(I) = (I-.5) * 2.*PI / NLON
 2345 CONTINUE
C
CL            2.35      INTEGRATION MATRICES OF HYDROSTATIC EQUATION(G)
CL                      AND CONVERSION TERM(A).
CL                      G COMPUTED AS IN CCM0 BUT ISOTHERMAL BOTTOM
CL                      AG CALCULATED TO CONSERVE ENERGY
C
      ILN2=NLEV*NLEV
      DO 2350 JK=1,ILN2
         G(JK)=0.
 2350 CONTINUE
C
      DO 2351 JK=1,NLEV-1
         ZALS(JK)=-.5*LOG(SIG(JK)/SIG(JK+1))
 2351 CONTINUE
C
      IG=0
C
      DO 2354 JK=1,NLEV-2
         G(IG+JK)=ZALS(JK)
         IF(JK.EQ.NLEV-2) GOTO 2353
         DO 2352 JL=JK+1,NLEV-2
            G(IG+JL)=ZALS(JL)+ZALS(JL-1)
 2352    CONTINUE
 2353    CONTINUE
         JL=NLEV
         G(IG+JL-1)=ZALS(JL-2)+ZALS(JL-1)
         G(IG+JL)=-.5*(LOG(SIG(JL))+LOG(SIG(JL-1)))
         IG = IG + NLEV
 2354 CONTINUE
C
      JL=NLEV
      G(IG+JL-1)=ZALS(JL-1)
      G(IG+JL)=-.5*(LOG(SIG(JL))+LOG(SIG(JL-1)))
      IG = IG + NLEV
      JK=NLEV
      G(IG+JK-1)=0.
      G(IG+JK)=-1.*LOG(SIG(JK))
C
      IA=1-NLEV
      IG=0
      DO 2355 JK=1,NLEV
         IA=IA+NLEV
         IG=IG+1
      DO 2355 JL=1,NLEV
         IAA=IA+JL-1
         IGG=IG+(JL-1)*NLEV
         AG(IAA)=G(IGG)*DSIGMA(JL)/DSIGMA(JK)
 2355 CONTINUE
C
C      WRITE(NOUT,7000)
C      WRITE(NOUT,7010) G
C      WRITE(NOUT,7020)
C      WRITE(NOUT,7010) AG
C7000  FORMAT(/' HYDROSTATIC MATRIX G IN SETCOM')
C7020  FORMAT(/' THERMODYNAMIC MATRIX AG IN SETCOM')
C7010  FORMAT(1X,9F13.9)
C
CL            2.5       MULTIPLY G BY R
C
      DO 255 IL=1,ILN2
         G(IL) = G(IL)*RAIR
  255 CONTINUE
c
c        Display resolutions, latitudes, gaussian weights, sigma levels
c
      zsum = 0.
      do 8000 j=1,norec
         zsum = zsum + w(j)
 8000 continue
c
      write (nout,8010)
 8010 format(/' -----------------------------------'
     *       /' *** Model Resolutions and Grids ***'
     *       /' -----------------------------------')
      write (nout,8012) nlon, norec, nlev
 8012 format(/' Agcm    longitude resolution:',i4
     *       /' Agcm    latitude  resolution:',i4
     *       /' Agcm    sigma levels        :',i4)
      if (nsrest.ge.0) write(nout,8013) nlonlsx, nlatlsx
 8013 format( ' Surface longitude resolution:',i4
     *       /' Surface latitude  resolution:',i4)
      write (nout,8030) (alatg(j)*180./3.1415927, j=1,norec)
 8030 format(/' Agcm latitude grid:' / (10(f8.3)) )
      write (nout,8040) zsum, w
 8040 format(/' Agcm gaussian weights (sum =',f8.4,'):' / (10(f8.4)))
      write (nout,8050) sig
 8050 format(/' Agcm sigma mid-levels:' / (10(f8.4)) )
      write (nout,8060) sigkmh
 8060 format(/' Agcm sigma level boundaries:' / (10(f8.4)) )
      write (nout,*)
 
c     Write gaussian latitudes and weights to scratch file for make*
c     write (99,8031) 'alatg'
c     write (99,8032) alatg
c     write (99,8031) 'wgaus'
c     write (99,8032) w
c8031 format(/'      data ',a,' /')
c8032 format(('     * ',4(f14.10,',')))
c     call endrun (0)
 
      RETURN
      END
C
C ***************************** SETINI ********************************
C
      SUBROUTINE SETINI
C----------------------------------------------------------------------
C
CL            INITIALIZE COMMON VARIABLES
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C---------------------------------------------------------------------
#include <COMFFT.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COMLEG.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C---------------------------------------------------------------------
#include <COMMP2.h>
C---------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMSDT.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
C
      DIMENSION ZALP(NMAX,MMAX), ZDALP(NMAX,MMAX)
C
C----------------------------------------------------------------------
C
CL            1.35 COMPUTE TRUNCATION PARAMETERS
C
C     WRITE(NOUT,912)
C     WRITE(NOUT,920)NTRM,NTRN,NTRK
C
      CALL TRUNCAT
C
CL            2.2  SET SVP TABLE AND CONSTANTS
C
      CALL CONSTS
C
CL            2.25 FILL REST OF COMMON BLOCK COMBUD
C
      DO 228 JJ=1,NOREC
         COSBUD(JJ) = W(JJ)/2.
  228 CONTINUE
C     WRITE(NOUT,950)
C     WRITE(NOUT,951) (COSBUD(I),I=1,NOREC)
C
C----------------------------------------------------------------------
C
CL            3.        COMPUTE SEMI-IMPL.TIMESTEP CONSTANTS (COMIMP)
C
      CALL SETTAU
C
C----------------------------------------------------------------------
C
CL            5.        COMPUTE LEGENDRE POLYNOMIALS AND DERIVATIVES
C                       ------- -------- ----------- --- -----------
C
C     NB: LAT INDEX FOR YALP AND YDALP IS FOR NORTHERN HEMISPHERE ONLY
C
      DO 550 JJ = 1,NOREC/2
         JJN = NOREC-JJ+1
         ZSIT = SIT(JJN)
C
         CALL PHCS (ZALP, ZDALP, MMAX, NMAX, ZSIT)
C
CL            5.3       MULTIPLY BY SQRT(2.) IN ORDER TO GET PROPER
CL                      NORMALIZATION. YDALP IS MULTIPLIED BY -1 TO
CL                      CORRECT FOR -SIGN IN COPENHAGEN DEFINITION.
CL                      NB: LAT INDEX FOR YALP AND YDALP IS FOR
CL                      NORTHERN HEMISPHERE ONLY.
C
         ZSQ2 = SQRT(2.)
         DO 531 JN=1,NMAX
            DO 532 JM=1,MMAX
               YALP(JM,JN,JJ) = ZALP(JN,JM)*ZSQ2
               YDALP(JM,JN,JJ) = -ZDALP(JN,JM)*ZSQ2
  532       CONTINUE
  531    CONTINUE
C
  550 CONTINUE
C
CL      5.6     CALCULATE ARRAYS NROW2,3,4(NOREC) FOR MAPPING
CL              BETWEEN NEW, OLD AND LEGENDRE LATITUDE INDICES
C
      DO 560 JJ = 1,NOREC/2
         NROW2(2*JJ) = JJ
         NROW2(2*JJ - 1) = NOREC + 1 - JJ
         NROW3(JJ) = 2*JJ
         NROW3(NOREC + 1 - JJ) = 2*JJ - 1
         NROW4(JJ) = JJ
         NROW4(NOREC + 1 - JJ) = JJ
  560 CONTINUE
C
C----------------------------------------------------------------------
C
CL            6.        COMPUTE CONSTANTS FOR FFT
CL                      ------- --------- --- ---
C
CL            6.1       SET UP TRIGONOMETRIC TABLES
C
C        Old Initialization (for old library ecmfft)
C
C     CALL FAX(IFAX,NLON,3)
C     CALL FFTRIG(TRIG,NLON,3)
C
C        New initialization (for library ecmfftn, SCD News 3/91, p.22)
C
      CALL SET99 (TRIG,IFAX,NLON)
C
C
      RETURN
C
C----------------------------------------------------------------------
C
CL            9.        FORMATS
C
C 912 FORMAT(/' MODEL HEMISPHERIC DOMAIN IS GLOBAL')
C 920 FORMAT(/'     TRUNCATION PARAMETERS'/'     NTRM = 'I4/
C    *        '     NTRN = 'I4/'     NTRK = 'I4/)
C 950   FORMAT(/' LATITUDINAL WEIGHTS FOR BUDGET INTEGRALS',
C    *          ' - COSBUD IN SETINI')
C 951 FORMAT(1X,10E12.3)
C
C----------------------------------------------------------------------
C
      END
C
C ***************************** CONSTS ********************************
C
      SUBROUTINE CONSTS
C----------------------------------------------------------------------
C
CL            INITIALIZE SVP TABLE AND CONSTANTS
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMBUD.h>
C---------------------------------------------------------------------
#include <COMCON.h>
C---------------------------------------------------------------------
#include <COMEST.h>
C---------------------------------------------------------------------
#include <COMMAP.h>
C---------------------------------------------------------------------
#include <COMTIM.h>
C---------------------------------------------------------------------
C
C
      DIMENSION TABL1(81),TABL2(123)
      EQUIVALENCE (TABL1(1),TABLE(1)),(TABL2(1),TABLE(82))
C
C        TABLE OF SATURATION VAPOR PRESSURE ES FROM 173.16 K  TO
C        375.16 K  IN ONE DEGREE INCREMENTS.
C
      DATA TMIN,TMAX/173.16,375.16/
C
      DATA TABL1     /.001403,.001719,.002101,.002561,.003117,.003784,
     1.004584,.005542,.006685,.008049,.009672, .01160, .01388, .01658,
     2 .01977, .02353, .02796, .03316, .03925, .04638, .05472, .06444,
     3 .07577, .08894,  .1042,  .1220,  .1425,  .1662,  .1936,  .2252,
     4  .2615,  .3032,  .3511,  .4060,  .4688,  .5406,  .6225,  .7159,
     5  .8223,  .9432,  1.080,  1.236,  1.413,  1.612,  1.838,  2.092,
     6  2.380,  2.703,  3.067,  3.476,  3.935,  4.449,  5.026,  5.671,
     7  6.393,  7.198,  8.097,  9.098,  10.21,  11.45,  12.83,  14.36,
     8  16.06,  17.94,  20.02,  22.33,  24.88,  27.69,  30.79,  34.21,
     9  37.98,  42.13,  46.69,  51.70,  57.20,  63.23,  69.85,  77.09,
     A  85.02,  93.70, 103.20/
      DATA TABL2             / 114.66, 127.20, 140.81, 155.67, 171.69,
     1 189.03, 207.76, 227.96, 249.67, 272.98, 298.00, 324.78, 353.41,
     2 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78, 656.62,
     3 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2, 1147.4,
     4 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3, 1936.7,
     5 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1, 3167.1,
     6 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1, 5030.7,
     7 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7, 7780.2,
     8 8201.5, 8642.3, 9103.4, 9585.5, 10089., 10616., 11166., 11740.,
     9 12340., 12965., 13617., 14298., 15007., 15746., 16516., 17318.,
     A 18153., 19022., 19926.,
     X                         20867., 21845., 22861., 23918., 25016.,
     1 26156., 27340., 28570., 29845., 31169., 32542., 33965., 35441.,
     2 36971., 38556., 40198., 41898., 43659., 45481., 47367., 49317.,
     3 51335., 53422., 55580., 57809., 60113., 62494., 64953., 67492.,
     4 70113., 72819., 75611., 78492., 81463., 84528., 87688., 90945.,
     5 94302., 97761., 101325., 104994., 108774., 108774./
C
      RETURN
      END
C
C ***************************** TOPINI ********************************
C
      SUBROUTINE TOPINI
C---------------------------------------------------------------------
C
CL            READS IN TOPOGRAPHY FIELD AND S.D. FIELD (FOR ENVELOPE
CL            OROGRAPHY AND GRAV WAVE DRAG), ADDS S.D. TO TOPOG FOR
CL            ENVELOPE OROGRAPHY, THEN RE-TRUNCATES TOPOG TO SPECTRAL
CL            RESOLUTION.
C
C---------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMBUD.h>
C---------------------------------------------------------------------
#include <COMCON.h>
C---------------------------------------------------------------------
#include <COMDRY.h>
C---------------------------------------------------------------------
#include <COMFLG.h>
C---------------------------------------------------------------------
#include <COMGWD.h>
C---------------------------------------------------------------------
#include <COMHED.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C---------------------------------------------------------------------
#include <COML10.h>
C---------------------------------------------------------------------
#include <COMMAP.h>
C---------------------------------------------------------------------
#include <COMMUN.h>
C---------------------------------------------------------------------
      dimension imapz(nlon,norec),  worka(nlon,norec), workb(nlon,norec)
      equivalence (imapz,workb)
C---------------------------------------------------------------------
 
c        Read in topography
 
      call findiu (ntopog, loc)
      write(nout,10) locname(loc)(1:lenchr(locname(loc)))
   10 format(/' Reading agcm topography data file'/4x,a)
      call openf (ntopog, locname(loc), ' ', 0, 0)
      call readdata (ntopog, yfi, nlon, norec, 'TOPO', 2, 1)
      call closef(ntopog)
 
c        Read in s.d. topog for envelope orography and gravity wave drag
 
      if (gwrough.lt.0.) then
        call findiu (ngwdrag, loc)
        write(nout,20) locname(loc)(1:lenchr(locname(loc)))
   20   format(/' Reading agcm standard-deviation topography file'/4x,a)
        call openf (ngwdrag, locname(loc), ' ', 0, 0)
        call readdata (ngwdrag, sdgravxy, nlon, norec, 'GWDH', 2, 1)
        call closef(ngwdrag)
        if (gwrough.ne.-999.) write(nout,22) -gwrough
   22   format(/' *** Single sub-grid gravity-wave roughness',
     *             f8.2,' m set for ice-sheet points')
      else
c       for global uniform single type over all non-ocean:
        write(nout,30) gwrough
   30   format(/' *** Using single sub-grid gravity-wave roughness',
     *          f8.2,' m for all land+ice points')
      endif
c
c        Add s.d (* envfac) to topog field for envelope orography
c        (envfac is a namelist param)
 
      if (envfac.ne.-999. .and. gwrough.lt.0.) then
        if (envfac.ge.0.) then
          do jj=1,norec
            do ji=1,nlon
              yfi(ji,jj) = yfi(ji,jj) + envfac*sdgravxy(ji,jj)
            enddo
          enddo
          write(nout,110) envfac
  110     format(/' Adding envfac*s.d. to topography (envfac=',f5.2,')')
        else
          do jj=1,norec
            do ji=1,nlon
              zlat = alatg(jj)*180./3.1415927
              zlon = along(ji)*180./3.1415927
              if (zlat.gt.  22. and. zlat.lt. 40. .and.      ! Himalayas
     *            zlon.gt.  65. and. zlon.lt. 95.) then      ! Himalayas
                yfi(ji,jj) = yfi(ji,jj) + abs(envfac)*sdgravxy(ji,jj)
              endif
            enddo
          enddo
          write(nout,112) envfac
  112     format(/' Adding |envfac|*s.d. to Himalayan topography only',
     *            ' (envfac=',f5.2,')')
        endif
      endif
 
c        Set topography at specific lon,lat points
 
c     do 160 jj=1,norec
c       zlat = alatg(jj)*180./3.1415927
c       do 162 ji=1,nlon
c         zlon = along(ji)*180./3.1415927
c         if (zlon.gt.180.) zlon = zlon-360.
cc        Central America:
c         if (zlat.gt.   4. and. zlat.lt. 26. .and.
c    *        zlon.gt.-140. and. zlon.lt.-40.) then
c           if (yfi(ji,jj).gt.700.) yfi(ji,jj) = 0.5*(yfi(ji,jj)+700.)
c         endif
cc        Greenland:
c         if (zlat.gt.  60. and. zlat.lt. 86. .and.
c    *        zlon.gt. -90. and. zlon.lt.-10.) then
c           yfi(ji,jj) = min (yfi(ji,jj), 100.)
c         endif
cc        Himalayas:
c         if (zlat.gt.  22. and. zlat.lt. 40. .and.
c    *        zlon.gt.  65. and. zlon.lt. 95.) then
c           yfi(ji,jj) = 1.3*yfi(ji,jj)
c         endif
c         Zero all:
c         yfi(ji,jj) = 0.
c 162   continue
c 160 continue
 
c       Save untruncated topog in topogun (commun) for history,
c       print out untruncated topog if requested
 
      do 170 jj=1,norec
        do 172 ji=1,nlon
          topogun(ji,jj) = yfi(ji,jj)
          imapz(ji,jj) = nint(topogun(ji,jj))
  172   continue
  170 continue
      if (showmaps) then
        call mapout (nout, imapz, nlon, norec, along, alatg, 999999,
     *               'Untruncated topography (m)')
      endif
 
c        Spectrally filter and truncate the topography field
 
c-----
c     call filter (yfi, 'isotropic')
c     call filter (yfi, 'none')
c-----
c        "Physical-space" filter of Navarra et al,1994,J.Climate,7,
c        1169-1183,section 6b.
 
      call scopy (nlon*norec, yfi, 1, worka, 1)
      call filter (yfi, 'none')
      call scopy (nlon*norec, yfi, 1, workb, 1)
      call filter (workb, 'isotropic')
      do 210 jj=1,norec
        do 212 ji=1,nlon
         yfi(ji,jj) = cvmgt (workb(ji,jj),yfi(ji,jj),worka(ji,jj).le.0.)
  212   continue
  210 continue
      call filter (yfi, 'none')
c-----
 
c        Calculate average elevation of "ocean" areas of
c        filtered/truncated field, apply global shift
 
      ztopoc = 0.
      zaroc  = 0.
      do 300 jj=1,norec
        do 302 ji=1,nlon
          if (worka(ji,jj).le.0.) then
            ztopoc = ztopoc + yfi(ji,jj)*cosbud(jj)
            zaroc = zaroc + cosbud(jj)
          endif
  302   continue
  300 continue
      if (ztopoc.ne.0. and. zaroc.gt.0.) then
        ztopoc = ztopoc/zaroc
        do 310 jj=1,norec
          do 312 ji=1,nlon
            yfi(ji,jj) = yfi(ji,jj) - ztopoc
  312     continue
  310   continue
      endif
 
c        Print out filtered and truncated topography
 
      if (showmaps) then
        do 400 jj=1,norec
          do 402 ji=1,nlon
            imapz(ji,jj) = nint(yfi(ji,jj))
c           Offset to show special symbols for "land" in mapout
            if (worka(ji,jj).gt.0.) imapz(ji,jj)= imapz(ji,jj)+1000000
  402     continue
  400   continue
        call mapout (nout, imapz, nlon, norec, along, alatg, -999999,
     *               'Truncated topography (m)')
      endif
 
      if (ztopoc.ge.0.) write(nout,456) 'increased', ztopoc
      if (ztopoc.lt.0.) write(nout,456) 'decreased', ztopoc
  456 format (/' Truncated topography ',a,    ' globally by',f9.2,'  m',
     *         ' (to make ocean mean = zero)')
 
c        Apply correction to dry mass of atmosphere (drymass, in
c        combud) so that global mean sea-level pressure (diagnostic)
c        will be approximately as observed (1013 mb) despite global
c        mean topography different from observed (237.3 m). The present
c        observed dry mass is 98222 N/m2, and the mean surface density
c        used to compute the correction below is 1.29 Kg/m3. The
c        correction to the model's surface pressure to make the
c        implicit dry mass = drymass is done in linems at each timestep.
 
c        For cold starts or older restart files, arbitrarily
c        set previous run's drymassin (otherwise has been read
c        from restart file). Drymass[in] is on the restart file just to
c        detect a change in a new run and so to re-initialiaze budget
c        (see linems).
 
      if (nsrest.eq.-1 .or. versresin.lt.2.01) drymassin = 0.
 
      ztop = 0.
      do 500 jj=1,norec
        do 502 ji=1,nlon
          ztop = ztop + yfi(ji,jj)*cosbud(jj)
  502   continue
  500 continue
      ztop = ztop/nlon
c     drymass = 98222. - gravit*1.29*(ztop-237.3)
c     Use dry mass 98400 to make model's slp close to observed (1013 mb)
c     drymass = 98400. - gravit*1.29*(ztop-237.3)
c     Or use lower value to agree with Shea's slp dataset (1011 mb)
      drymass = 98200. - gravit*1.29*(ztop-237.3)
 
      write(nout,504) ztop, drymass
  504 format (/' Global mean truncated surface topography =',f9.2,'  m'
     *        ,' (cf: present observed 237.3  m)'
     *        /' Global mean dry atmospheric mass         =',f9.2,' Pa'
     *        ,' (cf: present observed 98400 Pa)' )
 
      if (drymassin.ne.0. .and. abs(drymass-drymassin).gt.0.1) then
        write(nout,506) drymassin
  506   format("                     previous run's value =",f9.2," Pa")
c     else if (drymassin.eq.0. .and. nsrest.eq.1) then
c       write(nout,508)
c 508   format(' *** changed from previous run!')
      endif
 
c        Finally, change yfi from elevation to surface geopotential
 
      do 600 jj=1,norec
        do 602 ji=1,nlon
          yfi(ji,jj) = gravit*yfi(ji,jj)
  602   continue
  600 continue
 
      return
      end
C
C *************************** COLDSTRT ********************************
C
      SUBROUTINE COLDSTRT
C----------------------------------------------------------------------
C
CL            SETS AGCM FIELDS FOR COLD START. CURRENTLY TEMPERATURES
CL            DEPEND ON LATITUDE BUT NOT HEIGHT. EASTWARD VELOCITIES
CL            ARE GEOSTROPHIC. [SURFACE FIELDS WILL BE INITIALIZED BY
CL            SURFACE-MODEL ROUTINES BASED ON LOWEST-LEVEL AGCM TEMPS,
CL            PASSED IN FIELD YTSIN = YT(*,*,NLEV). THIS INCLUDES
CL            SEAICE, WHICH WILL EXIST INITIALLY AT OCEAN POINTS WHERE
CL            YTSIN .LE. 271.2 DEG K (SEE SUBR INILSX).]
CL
CL
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMEST.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C---------------------------------------------------------------------
#include <COMIMP.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C---------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COML30.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
      dimension zt(norec,nlev), zq(norec,nlev),
     *          zphi(norec,nlev), zu(norec,nlev), zf(norec)
      parameter (pi = 3.14159265358979)
c
c     ztprof is a mid-latitude "std-atmos" profile, not currently used
c     dimension ztprof(nlev)
c     data ztprof /233., 223., 216., 216., 216., 224.,
c    *             241., 257., 269., 278., 284., 288./
c
c     zrelhum is (global) relative humidity, zpsl is sea-level pressure
c     (zrelhum = 0.5 causes blowup - 12/15/91)
      data zrelhum, zpsl /0.0, 98222./
c----------------------------------------------------------------------
#include <SFEST.h>
c----------------------------------------------------------------------
c
c
      write(nout,50)
   50 format (/' Initializing for cold start (NSREST=-1)')
c
      do 100 jk=1,nlev
        do 102 jj=1,norec
 
c         zt(jj,jk) = teqcold - (teqcold-tpocold)*sin(alatg(jj))**2
          zx = abs(alatg(jj))/(50.*3.1415927/180.)             !for prot
          zt(jj,jk) = tpocold + (teqcold-tpocold)*exp(-zx**4)  !for prot
 
          zvp = zrelhum*estblf(zt(jj,jk))
          zq(jj,jk) = .622*zvp / (zpsl*sig(jk)-(1.-.622)*zvp)
  102   continue
  100 continue
c
c        Calculate geopotential heights (see CCM1 Description, pg. 40),
c        excluding surface geopotential since it will just be used below
c        for geostrophic velocities
c
      do 500 jj=1,norec
        do 502 jk=1,nlev
          zphi(jj,jk) = 0.
          do 504 jm=1,nlev
            ztv = zt(jj,jm) * (1.+zvir*zq(jj,jm))
            zphi(jj,jk) = zphi(jj,jk) + g((jk-1)*nlev+jm)*ztv
  504     continue
  502   continue
  500 continue
c
c        Calculate geostrophic eastward velocity
c
      do 510 jj=1,norec
        jjp = min (norec, jj+1)
        jjm = max (    1, jj-1)
        zf(jj) = 2. * omega * sin(alatg(jj))
        do 512 jk=1,nlev
          zu(jj,jk) = -(zphi(jjp,jk)-zphi(jjm,jk))
     *               / (zf(jj)*rearth*(alatg(jjp)-alatg(jjm)))
  512   continue
  510 continue
c
c        Set log_e(surface pressure), assuming isothermal in height
c
      zlpsl = log(zpsl)
      do 520 jj=1,norec
        do 522 ji=1,nlon
          ylps(ji,jj) = zlpsl - yfi(ji,jj) / (rair*zt(jj,nlev))
  522   continue
  520 continue
c
c        Set other full agcm fields
c
      do 530 jj=1,norec
        jjp = min (norec, jj+1)
        jjm = max (    1, jj-1)
        do 532 jk=1,nlev
          zz = zu(jj,jk)*tan(alatg(jj))/rearth
     *       - (zu(jjp,jk)-zu(jjm,jk))/(rearth*(alatg(jjp)-alatg(jjm)))
          do 534 ji=1,nlon
            yt(ji,jj,jk) = zt(jj,jk)
            yq(ji,jj,jk) = zq(jj,jk)
            yu(ji,jj,jk) = zu(jj,jk)
            yz(ji,jj,jk) = zf(jj) + zz
  534     continue
  532   continue
  530 continue
      call zero (yv, nlon*norec*nlev)
      call zero (yd, nlon*norec*nlev)
      call zero (ydpsl, nlon*norec)
      call zero (ydpsm, nlon*norec)
c
c        Copy from coml30 to coml10
c
      call scopy (nlon*norec,      ylps, 1, ylpsm1, 1)
      call scopy (nlon*norec*nlev, yu,   1, yum1,   1)
      call scopy (nlon*norec*nlev, yv,   1, yvm1,   1)
      call scopy (nlon*norec*nlev, yt,   1, ytm1,   1)
      call scopy (nlon*norec*nlev, yq,   1, yqm1,   1)
      call scopy (nlon*norec*nlev, yz,   1, yzm1,   1)
      call scopy (nlon*norec*nlev, yd,   1, ydm1,   1)
c
c        Display various cold-start fields versus level and latitude
c
      write (nout,900) (sig(jk),jk=nlev,1,-1)
  900 format( /' Cold-start temperatures (deg C):'
     *       //1x,5x,3x,(18f8.3)/)
      do 902 jj=norec,1,-1
        write(nout,904) alatg(jj)*180./pi,
     *                  (yt(1,jj,jk)-tmelt, jk=nlev,1,-1)
  904   format(1x,f5.1,3x,(18f8.2))
  902 continue
c
      write (nout,910) (sig(jk),jk=nlev,1,-1)
  910 format( /' Cold-start eastward velocities (m/s):'
     *       //1x,5x,3x,(18f8.3)/)
      do 912 jj=norec,1,-1
        write(nout,914) alatg(jj)*180./pi,
     *                  (yu(1,jj,jk), jk=nlev,1,-1)
  914   format(1x,f5.1,3x,(18f8.2))
  912 continue
c
      write (nout,920) (sig(jk),jk=nlev,1,-1)
  920 format( /' Cold-start (abs.vort.)/(plan.vort):'
     *       //1x,5x,3x,(18f8.3)/)
      do 922 jj=norec,1,-1
        write(nout,924) alatg(jj)*180./pi,
     *                  (yz(1,jj,jk)/zf(jj), jk=nlev,1,-1)
  924   format(1x,f5.1,3x,(18f8.3))
  922 continue
c
      write (nout,930) (sig(jk),jk=nlev,1,-1)
  930 format( /' Cold-start specific humidity:'
     *       //1x,5x,3x,(18f8.3)/)
      do 932 jj=norec,1,-1
        write(nout,934) alatg(jj)*180./pi,
     *                  (yq(1,jj,jk), jk=nlev,1,-1)
  934   format(1x,f5.1,3x,(18f8.5))
  932 continue
c
      return
      end
C
C ***************************** FILTER  ********************************
C
      subroutine filter (field, filtertype)
c
c       Spectrally truncates field to ccm spectral resolution,
c       exactly parallelling spectral code in old inidat.
c
c       filterype = 'isotropic' or 'exponential' to apply a spectral
c                   filter to reduce gibbs oscillations
c                   (nb: isotropic seems to give better results-9/94)
c
C---------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMLEG.h>
C---------------------------------------------------------------------
#include <COMMAP.h>
C---------------------------------------------------------------------
#include <COMTRU.h>
C---------------------------------------------------------------------
      dimension field(nlon,norec), xfield(mmax2,nmax)
      character*(*) filtertype
C---------------------------------------------------------------------
 
      call zero (xfield, 2*mmax*nmax)
 
c-----
 
      call fftccm( field(1,1), norec, -1)
 
c-----
 
      do 100 jjs=1,norec/2
        jjn = norec-jjs+1
        do 102 ji = 1,nlon
          field(ji,jjn) = (field(ji,jjn) - field(ji,jjs)) * 0.5
          field(ji,jjs) =  field(ji,jjn) + field(ji,jjs)
  102   continue
  100 continue
 
c-----
 
      do 120 jj=1,norec/2
        jjs = jj
        jjn = norec-jj+1
        zw = 2.*w(jj)
 
        do 122 jn=1,nmax
          if(mod(jn,2).eq.1) then
            jj1 = jjs
            jj2 = jjn
          else
            jj1 = jjn
            jj2 = jjs
          endif
          do 124 jma=1,2*mtrunc(jn)
            jm=(jma+1)/2
            xfield(jma,jn) = xfield(jma,jn)
     *                     + zw*yalp(jm,jn,jj)*field(jma,jj1)
  124     continue
  122   continue
 
  120 continue
 
c-----
 
c        If requested, apply a spectral filter to reduce Gibbs
c        oscillations (Navarra et al,1994,J.Climate,7,1169-1183)
 
      if (filtertype.eq.'isotropic') then
 
        do 300 jn=1,nmax
          do 302 jma=1,2*mtrunc(jn)
            jm=(jma+1)/2
            ntrue = jn+jm-2
            zfilt = 1. - float(ntrue)/float(ntrk+1)
            xfield(jma,jn) = zfilt * xfield(jma,jn)
  302     continue
  300   continue
 
      else if (filtertype.eq.'exponential') then
 
        do 310 jn=1,nmax
          do 312 jma=1,2*mtrunc(jn)
            jm=(jma+1)/2
            ntrue = jn+jm-2
            zfilt = exp(-32.*(float(ntrue)/float(ntrk))**(2*2))
            xfield(jma,jn) = zfilt * xfield(jma,jn)
  312     continue
  310   continue
 
      endif
 
c-----
 
      call zero( field(1,1), nlon*norec)
 
c-----
 
      do 500 jj=1,norec/2
        jjs = jj
        jjn = norec-jj+1
 
        do 502 jn=1,nmax
          if(mod(jn,2).eq.1) then
            jj1 = jjs
            jj2 = jjn
          else
            jj1 = jjn
            jj2 = jjs
          endif
          do 504 jma=1,2*mtrunc(jn)
            jm=(jma+1)/2
            field(jma,jj1) = field(jma,jj1)
     *                     + xfield(jma,jn)*yalp(jm,jn,jj)
  504     continue
  502   continue
 
  500 continue
 
c-----
 
      do 520 jjs = 1,norec/2
        jjn = norec-jjs + 1
        do 522 ji = 1,nlon
          field(ji,jjn) = field(ji,jjs) + field(ji,jjn)
          field(ji,jjs) = 2.*field(ji,jjs) - field(ji,jjn)
  522   continue
  520 continue
 
c-----
 
      call fftccm( field(1,1), norec, 1)
 
c-----
 
      return
      end
C
C ***************************** RSDS **********************************
C
      subroutine rsds (cname, cpath, runtitlex, mbdatex, nstepx, dtimex)
 
c             Acquire and read a Genesis initial/restart file
 
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMDRY.h>
C----------------------------------------------------------------------
#include <COMDYR.h>
C---------------------------------------------------------------------
#include <COMFLG.h>
C---------------------------------------------------------------------
#include <COMHED.h>
C---------------------------------------------------------------------
#include <COMHIS.h>
C---------------------------------------------------------------------
#include <COMHIX.h>
C---------------------------------------------------------------------
#include <COMHIY.h>
C---------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COML30.h>
C----------------------------------------------------------------------
#include <COML50.h>
C----------------------------------------------------------------------
#include <COMMUL.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
#include <COMSDT.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C---------------------------------------------------------------------
 
      character*(*) cname, cpath
      character*80 runtitlex
      parameter (pi = 3.14159265358979)
      character cdum16*16, cdum*8, cduma*8, cdumb*8
      equivalence (zdum,cdum), (zduma,cduma), (zdumb,cdumb)
 
c overall restart versions:
      parameter (nverspre=29)
      dimension verspre(nverspre)
      save verspre
      data verspre /
     *   1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08,     ! 9
     *   2.00, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07,           ! 8
     *   3.00, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08,     ! 9
     *   3.09, 3.10, 3.11 /                                        ! 3
 
c isotope restart versions:
      parameter (nverspre_isot=2)
      dimension verspre_isot(nverspre_isot)
      save verspre_isot
      data verspre_isot /
     *   0.00, 1.00 /                                              ! 1

      integer*8 i8a, i8b, i8c, i8d, i8e, i8f

      dimension workskip(nlon,norec)
 
c----------------------------------------------------------------------
 
       write(nout,10) mspathin(1:lenchr(mspathin)), msnamin
   10  format (/' ----------------------------'
     *         /' Opening initial/restart file ',a,a
     *         /' ----------------------------')
 
c         Acquire from permanent location
 
       call acquiref (cname, cpath)
       call openf (nsre, cname, cpath, 0, 0)
 
c=====
 
c        Nb: Must read exact # of words in each record
c        (see wsds) since file is a "pure data" stream.
c
      write(nout,1000)
 1000 format (/4x,'Reading header records...')

c        Determine if Cray restart file (if so, integer*8) 

      krec = 1
      read (nsre, err=8000, iostat=ios) cduma
      call closef (nsre)
      call openf (nsre, cname, cpath, 0, 0)

      if (cduma(1:4).eq.'CRAY') then
        read (nsre, err=8000, iostat=ios)
     *     cduma, runtitlex, i8a, i8b, dtimex, cdumb,
     *     i8c, i8d, i8e
        mbdatex = i8a
        nstepx  = i8b
        nlonx   = i8c
        norecx  = i8d
        nlevx   = i8e
      else
        read (nsre, err=8000, iostat=ios)
     *     cduma, runtitlex, mbdatex, nstepx, dtimex, cdumb,
     *     nlonx, norecx, nlevx
      endif
 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c overall restart version number (in cduma)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (cduma(1:4).eq.'CRAY' .or. cduma(1:4).eq.'IEEE') then
        typeresin = cduma(1:4)
        versresin = 0.
        read (cduma(5:8),'(f4.2)',err=1002) versresin
 1002   continue
      else
        typeresin = 'CRAY'
        call convertype (zduma, 1, 'R', 'CRAY')
        versresin = zduma
      endif

cisoe if (abs(versresin-3.11).lt.1.e-8) versresin = 3.05 ! isoe
 
      do 1010 i=1,nverspre
        if (abs(versresin-verspre(i)).lt.1.e-8) then
          versresin = verspre(i)
          goto 1016
        endif
 1010 continue
      write(nout,1012) msnamin(1:lenchr(msnamin)), versresin,
     *                 (verspre(i),i=1,nverspre)
 1012 format (
     * //'*** Incompatible version number of initial/restart file:'
     *  /4x,'Version number of initial/restart file ',a,':'
     *  /4x, f11.7
     *  /4x,'Version numbers accepted by program:'
     *  /(4x,9(f6.2)) )
      write(nout,1014)
 1014 format(
     *  /4x,'Stop in rsds.'/)
      call endrun (-1)
 1016 continue

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c isotope restart version number (in cdumb)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (cdumb(1:4).eq.'ISOT') then
        versresin_isot = 0.
        read (cdumb(5:8),'(f4.2)',err=1502) versresin_isot
 1502   continue
      else
        versresin_isot = 0.
      endif
cisoe versresin_isot = 1.0 ! isoe
 
      do 1510 i=1,nverspre_isot
        if (abs(versresin_isot-verspre_isot(i)).lt.1.e-8) then
          versresin_isot = verspre_isot(i)
          goto 1516
        endif
 1510 continue
      write(nout,1512) msnamin(1:lenchr(msnamin)), versresin_isot,
     *                 (verspre_isot(i),i=1,nverspre_isot)
 1512 format (
     * //'*** Incompatible isotope version # of initial/restart file:'
     *  /4x,'Isotope version number of initial/restart file ',a,':'
     *  /4x, f11.7
     *  /4x,'Isotope version numbers accepted by program:'
     *  /(4x,9(f6.2)) )
      write(nout,1014)
      call endrun (-1)
 1516 continue

      if (nsrest_isot.ne.-1 .and. versresin_isot.eq.0.) then
        write(nout,1518) msnamin(1:lenchr(msnamin)), versresin_isot
 1518   format (
     *   //' ***Isotope restart requested (NSREST_ISOT=1),'
     *    /4x,'but no isotope fields in initial/restart file ',a
     *    /4x,'(file isotope version number=',f4.2,')' )
        write(nout,1014)
        call endrun (-1)
      endif
c>>>>>>>>>>
 
      call convertype (mbdatex, 1, 'I', typeresin)
      call convertype (nstepx,  1, 'I', typeresin)
      call convertype (dtimex,  1, 'R', typeresin)
      call convertype (nlonx,   1, 'I', typeresin)
      call convertype (norecx,  1, 'I', typeresin)
      call convertype (nlevx,   1, 'I', typeresin)
 
      iyr = mbdatex/10000
      if (versresin.le.3.10) then
c       convert mbdatex years 0-200 to 1900-2100 A.D (<= v2.3, not v3) 
        if (iyr.ge.0.and.iyr.le.200) mbdatex=mbdatex + 1900*10000 !v2.3
      endif
 
      write(nout,2020)
     *  runtitlex(1:lenchr(runtitlex)), versresin, typeresin,
     *  mbdatex, nstepx*dtimex/86400., nstepx*dtimex/(86400.*365.),
     *  nlonx, norecx, nlevx
 2020 format ( 7x,'run title           : ',a
     *        /7x,'file version number :',f14.2
     *        /7x,'file machine type   : ',a
     *        /7x,'base date           :',i14
     *        /7x,'elapsed days        :',f14.3
     *        /7x,'elapsed years       :',f14.3
     *        /7x,'agcm resolutions    :',3i6)
 
 
c       Set expected resols (namelist params horres,levres) to model
c       resols if not entered. Must be done before first call to reader.
 
      if (horres(1).eq.0) horres(1) = nlon
      if (horres(2).eq.0) horres(2) = norec
      if (levres   .eq.0) levres    = nlev
 
      if ( nlonx.ne.horres(1) .or. norecx.ne.horres(2) .or.
     *     nlevx.ne.levres ) then
         write(nout,2030) nlonx,horres(1),norecx,horres(2),nlevx,levres
 2030    format (
     *    //' Agcm resolutions of initial/restart file do not match',
     *      ' the expected values.'
     *    / ' (Expected values are the current model''s, or',
     *      ' as specified by namelist parameters horres[],levres)'
     *    /4x,'Initial/restart file NLON:',i8,'   expected NLON:',i8
     *    /4x,'Initial/restart file NLAT:',i8,'   expected NLAT:',i8
     *    /4x,'Initial/restart file NLEV:',i8,'   expected NLEV:',i8
     *    /4x,'Stop in rsds.'/ )
         call endrun (-1)
      endif
 
c=====
 
      krec = krec + 1
      if (versresin.le.1.02) then
        read (nsre, err=8000, iostat=ios) ifile, ifilepre,
     *                                    msnamhisx, nuoff,
     *                                    msnam_ax,  nuoff_a,
     *                                    cdum16,    cdum
      else
        if (cduma(1:4).eq.'CRAY') then
          read (nsre, err=8000, iostat=ios) i8a, i8b,
     *                                      msnamhisx, i8c,   i8d,
     *                                      msnam_ax,  i8e,   i8f,
     *                                      cdum16,    cdum,    cdum
          ifile    = i8a
          ifilepre = i8b
          nuoff    = i8c
          nwrit    = i8d
          nuoff_a  = i8e
          nwrit_a  = i8f
        else
          read (nsre, err=8000, iostat=ios) ifile, ifilepre,
     *                                      msnamhisx, nuoff,   nwrit,
     *                                      msnam_ax,  nuoff_a, nwrit_a,
     *                                      cdum16,    cdum,    cdum
        endif
      endif
 
      call convertype (ifile,     1, 'I', typeresin)
      call convertype (ifilepre,  1, 'I', typeresin)
      call convertype (nuoff,     1, 'I', typeresin)
      call convertype (nuoff_a,   1, 'I', typeresin)
      if (versresin.gt.1.02) then
        call convertype (nwrit,   1, 'I', typeresin)
        call convertype (nwrit_a, 1, 'I', typeresin)
      endif
 
c=====
 
      krec = krec + 1
 
      if (versresin.le.3.03) then
 
        if (versresin.le.1.01) then
          nz = 40
        else if (versresin.le.1.03) then
          nz = 41
        else if (versresin.le.1.06) then
          nz = 43
        else if (versresin.le.1.08) then
          nz = 45
        else if (versresin.le.2.01) then
          nz = 49
        else if (versresin.le.3.03) then
          nz = 63
        endif
        read (nsre, err=8000, iostat=ios) (cdum, i=1,nz)
        call zero (tmpbud, nbudg)
        nbudgpre = 0
 
      else
 
        if (cduma(1:4).eq.'CRAY') then
          read (nsre, err=8000, iostat=ios) i8a
          nbudgpre = i8a
        else
          read (nsre, err=8000, iostat=ios) nbudgpre
        endif
        call convertype (nbudgpre, 1, 'I', typeresin)
        if (nsrest.eq.1 .and. nbudg.eq.nbudgpre) then
          read (nsre, err=8000, iostat=ios) (tmpbud(i),i=1,nbudg)
          call convertype (tmpbud, nbudg, 'R', typeresin)
        else
          read (nsre, err=8000, iostat=ios) (cdum, i=1,nbudgpre)
          call zero (tmpbud, nbudg)
          nbudgpre = 0
        endif
 
      endif
 
c=====
 
      krec = krec + 1
      if (nsrest.le.0) then
        if (cduma(1:4).eq.'CRAY') then
          read (nsre, err=8000, iostat=ios) xdum, i8a, (xdum, i=1,6)
        else
          read (nsre, err=8000, iostat=ios) xdum, idum, (xdum, i=1,6)
        endif
      else
        if (cduma(1:4).eq.'CRAY') then
          read (nsre, err=8000, iostat=ios)
     *      calday, i8a, ecc, obliq, prec, vern, dist, eccf
            jmnth = i8a
        else
          read (nsre, err=8000, iostat=ios)
     *      calday, jmnth, ecc, obliq, prec, vern, dist, eccf
        endif
        call convertype (calday, 1, 'R', typeresin)
        call convertype (jmnth,  1, 'I', typeresin)
        call convertype (ecc,    1, 'R', typeresin)
        call convertype (obliq,  1, 'R', typeresin)
        call convertype (prec,   1, 'R', typeresin)
        call convertype (vern,   1, 'R', typeresin)
        call convertype (dist,   1, 'R', typeresin)
        call convertype (eccf,   1, 'R', typeresin)
      endif
 
c=====
 
      if (versresin.le.1.08) then                 ! tmass, tmass0, qmass
        krec = krec + 1
        read (nsre, err=8000, iostat=ios) (cdum, i=1,3)
      endif
 
c=====
 
      if (versresin.ge.2.01) then
        krec = krec + 1
        read (nsre, err=8000, iostat=ios) drymassin
        call convertype (drymassin, 1, 'R', typeresin)
      endif
 
c=====
 
c        Read agcm field arrays
 
      write(nout,6000)
 6000 format (/4x,'Reading agcm    fields...')
 
c        common coml30:
 
      call reader (nsre, yf,      nlon,norec,1,    1)
      call reader (nsre, yu,      nlon,norec,nlev, 1)
      call reader (nsre, yv,      nlon,norec,nlev, 1)
      call reader (nsre, ylps,    nlon,norec,1,    1)
      call reader (nsre, yt,      nlon,norec,nlev, 1)
      call reader (nsre, yq,      nlon,norec,nlev, 1)
      call reader (nsre, yz,      nlon,norec,nlev, 1)
      call reader (nsre, yd,      nlon,norec,nlev, 1)
      if (versresin.le.1.08) then
        call reader (nsre, histmp, nlon,norec,nlev, 1)    ! yduh
        call reader (nsre, histmp, nlon,norec,nlev, 1)    ! ydvh
        call reader (nsre, histmp, nlon,norec,nlev, 1)    ! ydth
        call reader (nsre, histmp, nlon,norec,nlev, 1)    ! ydqh
      endif
      call reader (nsre, ydpsl,   nlon,norec,1, 1)
      call reader (nsre, ydpsm,   nlon,norec,1, 1)
 
c        common coml10:
 
      call reader (nsre, yfi,     nlon,norec,1,    1)
      call reader (nsre, yum1,    nlon,norec,nlev, 1)
      call reader (nsre, yvm1,    nlon,norec,nlev, 1)
      call reader (nsre, ylpsm1,  nlon,norec,1,    1)
      call reader (nsre, ytm1,    nlon,norec,nlev, 1)
      call reader (nsre, yqm1,    nlon,norec,nlev, 1)
      call reader (nsre, yzm1,    nlon,norec,nlev, 1)
      call reader (nsre, ydm1,    nlon,norec,nlev, 1)
 
c        Add "mean" temperature for earlier file versions, and divide
c        u,v by cos(lat) to get actual velocities. use approximate
c        latitudes since gaussian lats not calculated yet.
 
      if (versresin.le.1.08) then
        do 6005 jk=1,nlev
          do 6006 jj=1,norec
            zclat = cos ( -pi/2. + pi*(jj-.5)/norec )
            do 6007 ji=1,nlon
              ytm1(ji,jj,jk) = ytm1(ji,jj,jk) + t0(jk)
              yum1(ji,jj,jk) = yum1(ji,jj,jk)/zclat
              yvm1(ji,jj,jk) = yvm1(ji,jj,jk)/zclat
              yt(ji,jj,jk) = yt(ji,jj,jk) + t0(jk)
              yu(ji,jj,jk) = yu(ji,jj,jk)/zclat
              yv(ji,jj,jk) = yv(ji,jj,jk)/zclat
 6007       continue
 6006     continue
 6005   continue
      endif
 
c        Skip records of earlier restart file versions (use levres,
c        not nlev, to calculate number of horiz slices to skip)
 
      if (versresin.le.1.08) then
        if (versresin.le.1.00) then
          nz = 28 + 15*levres +  4*(levres+1) + 1*(levres+1)*(levres+1)
        else if (versresin.le.1.02) then
          nz = 16 + 5*levres
        else if (versresin.le.1.03) then
          nz = 16 + 6*levres
        else if (versresin.le.1.06) then
          nz = 17 + 5*levres
        else if (versresin.le.1.08) then
          nz = 16 + 6*levres
        endif
        do 6010 iz=1,nz
          call reader (nsre, histmp, nlon,norec,1, 1)
 6010   continue
      endif
 
c        Precip, surface upward ir flux, slt 1st-guess pts,
c        vert-mix reservoirs (for cold starts or old versresin,
c        these are set in basic):
 
      if (versresin.ge.2.00) then
        if (versresin.eq.2.00) then
          call reader (nsre, histmp,nlon,norec,nlev, 1)    ! ycld
          call reader (nsre, histmp,nlon,norec,nlev, 1)    ! ycldc
          call reader (nsre, histmp,nlon,norec,nlev, 1)    ! yclda
          call reader (nsre, histmp,nlon,norec,nlev, 1)    ! ycldl
        else if (versresin.eq.2.01) then
          call reader (nsre, histmp,nlon,norec,nlev, 1)    ! ycld
          call reader (nsre, histmp,nlon,norec,nlev, 1)    ! ycldl
        endif
        call reader (nsre, frainc,  nlon,norec,1,    1)
        call reader (nsre, fsnowc,  nlon,norec,1,    1)
        call reader (nsre, frains,  nlon,norec,1,    1)
        call reader (nsre, fsnows,  nlon,norec,1,    1)
        call reader (nsre, firems,  nlon,norec,1,    1)
        call reader (nsre, alonsav, nlon,norec,nlev, 1)
        call reader (nsre, alatsav, nlon,norec,nlev, 1)
        call reader (nsre, sigsav,  nlon,norec,nlev, 1)
        call reader (nsre, wsltsav, nlon,norec,nlev, 1)
      endif
 
c        common coml50 (version 2.05 has liquid water contents in
c        tracers #1-3, versions 2.06+ have flux reservoirs in tracer
c        #4-7). Read one tracer at a time, so readerv can work ok.
 
      call zero (tracer, nlon*norec*nlev*ntrace)
 
      if (versresin.ge.1.05 .and. versresin.le.2.05) then
        call reader (nsre, fplumt,  nlon,norec,nlev, 1)
        call reader (nsre, fplumq,  nlon,norec,nlev, 1)
        call reader (nsre, fplumu,  nlon,norec,nlev, 1)
        call reader (nsre, fplumv,  nlon,norec,nlev, 1)
      endif
 
      if (versresin.ge.1.06 .and. versresin.le.2.01) then
 
        read (nsre) ntracex                               ! # of tracers
        call convertype (ntracex, 1, 'I', typeresin)
        do 6030 n=1,2*ntracex
          call reader (nsre,histmp,nlon,norec,nlev,1) ! tracers,fplumtra
 6030   continue
 
      else if (versresin.ge.2.02 .and. versresin.le.2.04) then
 
        read (nsre) ntracex
        call convertype (ntracex, 1, 'I', typeresin)
        do 6040 n=1,ntracex
          call reader (nsre,histmp,nlon,norec,nlev,1)      ! old tracers
 6040   continue
 
      else if (versresin.eq.2.05) then
 
c       get liquid water amounts (equived to tracers #1-3):
        read (nsre) ntracex
        call convertype (ntracex, 1, 'I', typeresin)
        do 6050 n=1,min(ntracex,3)
          call reader (nsre, tracer(1,1,1,n), nlon,norec,nlev, 1)
 6050   continue
c       skip over any excess tracers on restart file
        if (ntracex.gt.3) then
          do 6052 n=1,ntracex-3
            call reader (nsre,histmp,nlon,norec,nlev,1)    ! old tracers
 6052     continue
        endif
 
      else if (versresin.ge.2.06) then
 
c       get liquid water amounts (#1-3), flux reservoirs(#4-7),
c       and actual tracers (#8-ntracex):
        if (cduma(1:4).eq.'CRAY') then
          read (nsre) i8a
          ntracex = i8a
        else
          read (nsre) ntracex
        endif
        call convertype (ntracex, 1, 'I', typeresin)
        do 6060 n=1, min(ntracex,ntrace)
          call reader (nsre, tracer(1,1,1,n), nlon,norec,nlev, 1)
 6060   continue
c       Skip over any excess tracers on restart file
        if (ntracex.gt.ntrace) then
          do 6062 n=1,ntracex-ntrace
            call reader (nsre,histmp,nlon,norec,nlev,1)    ! old tracers
 6062     continue
        endif
 
c       warn if inconsistent number of tracers:
c       if (ntracex.ne.ntrace) then
c         write(nout,6064) ntracex,ntrace
c6064     format (
c    *      /' *** Warning ***'
c    *      /'     Number of tracer fields on initial/restart file',
c    *             ' differs from program number.'
c    *      /'     Initial/restart file NTRACE=',i6,
c    *             '   Program NTRACE=',i6)
c         if (ntracex.lt.ntrace) then
c           write(nout,6066)
c6066       format(' Excess program tracer fields set to zero.'/)
c         else
c           write(nout,6068)
c6068       format(' Excess initial/restart tracer fields ignored.'/)
c         endif
c       endif
 
      endif
 
c=====
#if defined (LKAEROSOL)
      if (versresin.eq.3.11) call reader (nsre,sfcirrad,nlon,norec,1,1)
#else
      if (versresin.eq.3.11) call reader (nsre,workskip,nlon,norec,1,1)
#endif
c         Read isotopic and regional tracers fields.
c         If nsrest_isot = -1, these will be reset in subr basic

#ifdef isotrac
      if (versresin_isot.ge.1.00) then
        write(nout,6002)
 6002   format (/4x,'Reading agcm isotopic  fields...')
        call reader (nsre, yq_o,     nlon, norec, nlev, 1)
        call reader (nsre, yq_h,     nlon, norec, nlev, 1)
        call reader (nsre, yqm1_o,   nlon, norec, nlev, 1)
        call reader (nsre, yqm1_h,   nlon, norec, nlev, 1)
        call reader (nsre, frainc_o, nlon, norec, 1,    1)
        call reader (nsre, frainc_h, nlon, norec, 1,    1)
        call reader (nsre, fsnowc_o, nlon, norec, 1,    1)
        call reader (nsre, fsnowc_h, nlon, norec, 1,    1)
        call reader (nsre, frains_o, nlon, norec, 1,    1)
        call reader (nsre, frains_h, nlon, norec, 1,    1)
        call reader (nsre, fsnows_o, nlon, norec, 1,    1)
        call reader (nsre, fsnows_h, nlon, norec, 1,    1)
        call reader (nsre, fplumq_o, nlon, norec, nlev, 1)
        call reader (nsre, fplumq_h, nlon, norec, nlev, 1)
        call reader (nsre, ylwc_c_o, nlon, norec, nlev, 1)
        call reader (nsre, ylwc_c_h, nlon, norec, nlev, 1)
        call reader (nsre, ylwc_a_o, nlon, norec, nlev, 1)
        call reader (nsre, ylwc_a_h, nlon, norec, nlev, 1)
        call reader (nsre, ylwc_s_o, nlon, norec, nlev, 1)
        call reader (nsre, ylwc_s_h, nlon, norec, nlev, 1)
        do it=1,nrtr
          call reader (nsre, rtrq  (1,1,1,it), nlon, norec, nlev, 1)
          call reader (nsre, rtrqm1(1,1,1,it), nlon, norec, nlev, 1)
          call reader (nsre, rtrprc(1,1,it)  , nlon, norec, 1,    1)
          call reader (nsre, rtrprs(1,1,it)  , nlon, norec, 1,    1)
          call reader (nsre, rtrpsc(1,1,it)  , nlon, norec, 1,    1)
          call reader (nsre, rtrpss(1,1,it)  , nlon, norec, 1,    1)
          call reader (nsre, rtrfq (1,1,1,it), nlon, norec, nlev, 1)
        enddo
      endif
#else
      if (versresin_isot.ge.1.00) then
        write(nout,6004)
 6004   format (/4x,'Skipping agcm isotopic  fields...')
        do iskip=1, 12*nlev + 8 + nrtr*(3*nlev+4)
          call reader (nsre, workskip, nlon, norec, 1, 1)
        enddo
      endif
#endif

c=====
 
c        Read lsx field arrays
c
      write(nout,6500)
 6500 format (/4x,'Reading surface fields...')
c
      call readlsx (nsre, nout, versresin, typeresin,
     *              horreslsx(1), horreslsx(2), nlonlsx, nlatlsx,
     *              versresin_isot, nsrest_isot, cduma)
 
c=====
 
c        Last-word check
 
      zdum = -999.
      if (versresin.ge.2.07) then
        read (nsre, err=7002, iostat=ios) zdum
        call convertype (zdum, 1, 'R', typeresin)
cisoe   if (abs(zdum-3.11).lt.1.e-8) zdum = 3.05 ! isoe
        do 7000 i=1,nverspre
          if (abs(zdum-verspre(i)).lt.1.e-8) zdum = verspre(i)
 7000   continue
 7002   if (zdum.ne.versresin) then
          write(nout,7004) versresin, zdum, cdum, ios
 7004     format(/' *** Error: incorrect last word on',
     *                  ' initial/restart file'
     *           /'     expected     = ',f8.3
     *           /'     found (real) = ',f25.15
     *           /'     found (char) = ',a8
     *           /'     iostat       = ',i8
     *           /'     Stop in rsds.'/)
          call endrun (-1)
        endif
      endif
 
      call closef (nsre)
 
      write (nout,"(' ----------------------------')")
 
      return
 
c=====
 
c        Fatal errors for reading Genesis restart-file header records
 
 8000 write(nout,8002) krec, ios
 8002 format(/' *** Error or EOF reading header record number',i3,
     *              ' of initial/restart file.'
     *       /'     iostat = ',i8
     *       /'     Stop in rsds.'/)
      call endrun (-1)
      end
C
C ***************************** WSDS **********************************
C
      subroutine wsds (cname, cpath)
 
c             Write and dispose a Genesis restart file.
c             cname = file name (supplied)
c             cpath = path name (supplied)
 
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMDRY.h>
C----------------------------------------------------------------------
#include <COMDYR.h>
C---------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMHED.h>
C----------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMHIY.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COML30.h>
C----------------------------------------------------------------------
#include <COML50.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
#include <COMSDT.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
 
      character*(*) cname, cpath
 
      character*16 cdum16*16,  cdum*8, cduma*8, cdumb*8, runtitlx*80
      save cdum16, cdum
      data cdum16, cdum /' ', ' '/
c----------------------------------------------------------------------

c        If requested, put current file counter at end of runtitle

      runtitlx = runtitle
      i = lenchr(runtitle)
      if (runtitle(i-2:i).eq.'...') 
     *  write (runtitlx(i-2:i),'(i3.3)') ifilepre
c
      call openf (nsre, cname, cpath, 1, 0)
c
c        Restart file version #, run title, time and grid info
c
      write (cduma, '(a4,f4.2)') typeres(1:4), versres
      write (cdumb, '(a4,f4.2)') 'ISOT', versres_isot 
      write (nsre) cduma, runtitlx, mbdate, nstep, dtime, cdumb,
     *             nlon, norec, nlev
 
c        History file information
 
      write (nsre) ifile, ifilepre,
     *             msnamhis(ifilepre), nuoff,    nwrit,
     *             msnam_a(ifilepre),  nuoff_a,  nwrit_a,
     *             cdum16,             cdum,     cdum
 
c        Budget quantitites in combud
 
      write (nsre) nbudg
      write (nsre) tmpbud
 
c        Current calendar and orbit
 
      write (nsre) calday, jmnth, ecc, obliq, prec, vern, dist, eccf
 
c        Dry atmospheric mass
 
      write (nsre) drymass
 
c        Agcm field arrays
 
c     common coml30:
      call writer (nsre, yf,      nlon,norec,1,    1)
      call writer (nsre, yu,      nlon,norec,nlev, 1)
      call writer (nsre, yv,      nlon,norec,nlev, 1)
      call writer (nsre, ylps,    nlon,norec,1,    1)
      call writer (nsre, yt,      nlon,norec,nlev, 1)
      call writer (nsre, yq,      nlon,norec,nlev, 1)
      call writer (nsre, yz,      nlon,norec,nlev, 1)
      call writer (nsre, yd,      nlon,norec,nlev, 1)
      call writer (nsre, ydpsl,   nlon,norec,1,    1)
      call writer (nsre, ydpsm,   nlon,norec,1,    1)
 
c     common coml10:
      call writer (nsre, yfi,     nlon,norec,1,    1)
      call writer (nsre, yum1,    nlon,norec,nlev, 1)
      call writer (nsre, yvm1,    nlon,norec,nlev, 1)
      call writer (nsre, ylpsm1,  nlon,norec,1,    1)
      call writer (nsre, ytm1,    nlon,norec,nlev, 1)
      call writer (nsre, yqm1,    nlon,norec,nlev, 1)
      call writer (nsre, yzm1,    nlon,norec,nlev, 1)
      call writer (nsre, ydm1,    nlon,norec,nlev, 1)
 
c     precip, surface upward ir flux, slt 1st-guess pts:
      call writer (nsre, frainc,  nlon,norec,1,    1)
      call writer (nsre, fsnowc,  nlon,norec,1,    1)
      call writer (nsre, frains,  nlon,norec,1,    1)
      call writer (nsre, fsnows,  nlon,norec,1,    1)
      call writer (nsre, firems,  nlon,norec,1,    1)
      call writer (nsre, alonsav, nlon,norec,nlev, 1)
      call writer (nsre, alatsav, nlon,norec,nlev, 1)
      call writer (nsre, sigsav,  nlon,norec,nlev, 1)
      call writer (nsre, wsltsav, nlon,norec,nlev, 1)
c
c     common coml50 (liquid water contents(#1-3), flux reservoirs(#4-7),
c     and other tracers):
c
      write (nsre) ntrace
      call writer (nsre, tracer,  nlon,norec,nlev*ntrace, 1)

#if defined (LKAEROSOL)
      call writer (nsre, sfcirrad, nlon,norec,1,   1)
#endif

c        isotopic and regional tracer fields:

#ifdef isotrac
        call writer (nsre, yq_o,     nlon, norec, nlev, 1)
        call writer (nsre, yq_h,     nlon, norec, nlev, 1)
        call writer (nsre, yqm1_o,   nlon, norec, nlev, 1)
        call writer (nsre, yqm1_h,   nlon, norec, nlev, 1)
        call writer (nsre, frainc_o, nlon, norec, 1,    1)
        call writer (nsre, frainc_h, nlon, norec, 1,    1)
        call writer (nsre, fsnowc_o, nlon, norec, 1,    1)
        call writer (nsre, fsnowc_h, nlon, norec, 1,    1)
        call writer (nsre, frains_o, nlon, norec, 1,    1)
        call writer (nsre, frains_h, nlon, norec, 1,    1)
        call writer (nsre, fsnows_o, nlon, norec, 1,    1)
        call writer (nsre, fsnows_h, nlon, norec, 1,    1)
        call writer (nsre, fplumq_o, nlon, norec, nlev, 1)
        call writer (nsre, fplumq_h, nlon, norec, nlev, 1)
        call writer (nsre, ylwc_c_o, nlon, norec, nlev, 1)
        call writer (nsre, ylwc_c_h, nlon, norec, nlev, 1)
        call writer (nsre, ylwc_a_o, nlon, norec, nlev, 1)
        call writer (nsre, ylwc_a_h, nlon, norec, nlev, 1)
        call writer (nsre, ylwc_s_o, nlon, norec, nlev, 1)
        call writer (nsre, ylwc_s_h, nlon, norec, nlev, 1)
        do it=1,nrtr
          call writer (nsre, rtrq  (1,1,1,it), nlon, norec, nlev, 1)
          call writer (nsre, rtrqm1(1,1,1,it), nlon, norec, nlev, 1)
          call writer (nsre, rtrprc(1,1,it)  , nlon, norec, 1,    1)
          call writer (nsre, rtrprs(1,1,it)  , nlon, norec, 1,    1)
          call writer (nsre, rtrpsc(1,1,it)  , nlon, norec, 1,    1)
          call writer (nsre, rtrpss(1,1,it)  , nlon, norec, 1,    1)
          call writer (nsre, rtrfq (1,1,1,it), nlon, norec, nlev, 1)
        enddo
#endif
 
c        Lsx field arrays
 
      call writlsx (nsre)
 
c        Last-word check
 
      write (nsre) versres
 
c        Pad non-Cray restart files so can be read on Crays with
c        "assign -s bin" (otherwise get EOF due to Cray buffering?)
 
#ifndef cray
      do 800 i=1,4096
        write(nsre) char(0)
  800 continue
#endif
 
c        Dispose to permanent location
 
      call closef (nsre)
      call disposef (cname, cpath)
 
      return
      end
