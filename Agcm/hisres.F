c -- $Id: hisres.F,v 1.2 2001/07/03 21:29:18 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/hisres.F,v $
c -- $Name:  $

C
C ***************************** BUDCLC ********************************
C
      SUBROUTINE BUDCLC (ZTIME, YPTHIC, YPTHICP)
C----------------------------------------------------------------------
C
CL            CALCULATE BUDGET INTEGRALS, AND SET BUDGET END TIME
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COML30.h>
C----------------------------------------------------------------------
#include <COML50.h>
C----------------------------------------------------------------------
#include <COMMP2.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
      DIMENSION YPTHIC(NLON,NOREC,NLEV), YPTHICP(NLON,NOREC,NLEV)
C----------------------------------------------------------------------
C
CL            1.        SET END TIME FOR THIS BUDGET PERIOD
CL                      (SECONDS SINCE BASE DATE)
C
      BTIME    = ZTIME
      PE       = 0.
      SE       = 0.
      KE       = 0.
      WAVA     = 0.
      WLWC     = 0.
      FPLUMT_G = 0.
      FPLUMQ_G = 0.
C
CL    1.0       COMPUTE GLOBAL MEANS OF PE,KE,SE,WAVA,WLWC.
CL              UNITS ARE J M-2 (PE,KE,SE) OR KG M-2 (WAVA,WLWC)
CL
CL              ERRORS ARE CAUSED BY USING YQ AND YQM1 (ETA GRID)
CL              FOR PE (SIGMA GRID), USING FPLUM[U,V] (ETA GRID) for KE
CL              (SIGMA GRID).
C
      DO 100 JK = 1,NLEV
        DO 102 JJ = 1,NOREC
          ZWEI = COSBUD(JJ)/NLON
          DO 104 JI = 1,NLON
            ZM1 = YPSM1(JI,JJ)*DSIGMA(JK)/GRAVIT
            ZM  = YPS  (JI,JJ)*DSIGMA(JK)/GRAVIT
            ZME1= YPTHIC (JI,JJ,JK)/GRAVIT
            ZME = YPTHICP(JI,JJ,JK)/GRAVIT
 
            ZLWC= YLWC_C(JI,JJ,JK) + YLWC_A(JI,JJ,JK) + YLWC_S(JI,JJ,JK)
            ZQ1 = YQM1(JI,JJ,JK) + FPLUMQ(JI,JJ,JK)/ZME1
            ZQ  = YQ  (JI,JJ,JK) + FPLUMQ(JI,JJ,JK)/ZME
            ZT1 = YTM1(JI,JJ,JK)
            ZT  = YT  (JI,JJ,JK)
            ZU1 = YUM1(JI,JJ,JK) + FPLUMU(JI,JJ,JK)/ZM1
            ZV1 = YVM1(JI,JJ,JK) + FPLUMV(JI,JJ,JK)/ZM1
            ZU  = YU  (JI,JJ,JK) + FPLUMU(JI,JJ,JK)/ZM
            ZV  = YV  (JI,JJ,JK) + FPLUMV(JI,JJ,JK)/ZM
 
            PE = PE + 0.5*ZWEI
     *      * (   (  (1.-ZQ1)*CPAIR*ZT1
     *              +    ZQ1 *(CPWV*(ZT1-TMELT)+LATSUB)
     *              +    ZLWC*LATICE
     *            ) * ZM1
     *          + (  (1.-ZQ )*CPAIR*ZT
     *              +    ZQ  *(CPWV*(ZT -TMELT)+LATSUB)
     *              +    ZLWC*LATICE
     *            ) * ZM
     *        )
     *      + ZWEI*FPLUMT(JI,JJ,JK)
 
            WAVA = WAVA + 0.5*ZWEI*(ZQ1 *ZME1 + ZQ  *ZME)
            WLWC = WLWC + 0.5*ZWEI*(ZLWC*ZME1 + ZLWC*ZME)
 
            KE = KE
     *         + 0.5*ZWEI*0.5*((ZU1**2+ZV1**2)*ZM1 + (ZU**2+ZV**2)*ZM)
 
  104     CONTINUE
  102   CONTINUE
  100 CONTINUE
C
      DO 120 JJ = 1,NOREC
        ZWEI = COSBUD(JJ)/NLON
        DO 122 JI=1,NLON
          SE = SE + ZWEI*0.5*(YPSM1(JI,JJ)+YPS(JI,JJ))*YFI(JI,JJ)/GRAVIT
  122   CONTINUE
  120 CONTINUE
C
CL    2.0       COMPUTE GLOBAL MEANS OF RESERVOIRS FPLUMT AND FPLUMQ
CL              (WHOSE UNITS ARE FLUX*TIME PER LAYER, SO UNITS OF THE
CL              GLOBAL MEANS ARE J M-2 and KG M-2).
C
      DO 200 JK = 1,NLEV
        DO 202 JJ = 1,NOREC
          ZWEI = COSBUD(JJ)/NLON
          DO 204 JI = 1,NLON
            FPLUMT_G = FPLUMT_G + FPLUMT(JI,JJ,JK)*ZWEI
            FPLUMQ_G = FPLUMQ_G + FPLUMQ(JI,JJ,JK)*ZWEI
  204     CONTINUE
  202   CONTINUE
  200 CONTINUE
 
      RETURN
      END
C
C ***************************** BUDINA ********************************
C
      SUBROUTINE BUDINA
C----------------------------------------------------------------------
C
CL            INITIALIZE BUDGET ACCUMULATORS AND TRACKERS, AND SAVE
CL            OLD BUDGET INTEGRALS AND TIME
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
c     DTTOL is tolerance for initial budget test in linems. BTIME0 is
c     set here just to avoid indef error in that test for cold starts.
      DATA BTIME0, DTTOL /0., 30./
C----------------------------------------------------------------------
C
CL            1.        INITIALIZE FLUX ACCUMULATORS, SLT FIX (DWFIX*)
CL                      AND NEGATIVE WATER VAPOR FIX (DWNEG IN VDIF)
C
      DO 10 JJ=1,NOREC
        DHCON(JJ)  = 0.
        DHRADS(JJ) = 0.
        DHRADL(JJ) = 0.
        DRADTP(JJ) = 0.
        DSOLIN(JJ) = 0.
        DSOLRF(JJ) = 0.
        DFIRTP(JJ) = 0.

        DTAIR(JJ)   = 0.
        DTOTCLD(JJ) = 0.
C
        DMCLOUF(JJ) = 0.
        DMCLOUC(JJ) = 0.
        DMCLOUE(JJ) = 0.
        DMCLOUS(JJ) = 0.
        DHCLOUD(JJ) = 0.
        DMREEV(JJ) = 0.
        DHREEV(JJ) = 0.
        DMPRECC(JJ) = 0.
        DMPRECA(JJ) = 0.
        DMPRECS(JJ) = 0.
        DMCLOUP(JJ) = 0.
        DMCLOUA(JJ) = 0.
        DHVDIF(JJ) = 0.
        DHHDIF(JJ) = 0.
C
        DKEVDF(JJ) = 0.
        DKEHDF(JJ) = 0.
C
        DRFLUX(JJ) = 0.
        DSFLUX(JJ) = 0.
        DLFLUX(JJ) = 0.
        DEFLUX(JJ) = 0.
        DPFLUX(JJ) = 0.
        DWNEG(JJ)  = 0.
   10 CONTINUE
C
      DWFIX  = 0.
      DWFIXA = 0.
C
CL            2.        INITIALIZE MIN AND MAX TEMPERATURE TRACKERS
C
      TMPMIN = 999.
      MINTLO = 0
      MINTLA = 0
      MINTLV = 0
      TMPMAX = -999.
      MAXTLO = 0
      MAXTLA = 0
      MAXTLV = 0
      ITMMIN = -999
      ITMMAX = -999
C
CL            3.        INITIALIZE ATMOSPHERIC DRY-MASS ADJUSTMENT
CL                      VARIABLES
C
      DMASS  = 0.
      DMASSA = 0.
C
CL            4.        SAVE OLD VALUES OF THE BUDGET INTEGRALS AND
CL                      START TIME FOR NEXT BUDGET PERIOD
C
      BTIME0    = BTIME
      PE0       = PE
      SE0       = SE
      KE0       = KE
      WAVA0     = WAVA
      WLWC0     = WLWC
      FPLUMT_G0 = FPLUMT_G
      FPLUMQ_G0 = FPLUMQ_G
C
      RETURN
      END
C
C ***************************** BUDPRN ********************************
C
      SUBROUTINE BUDPRN (NBUD)
C----------------------------------------------------------------------
C
CL            PRINT OUT BUDGET SUMMARY TO UNIT NUMBER NBUD.
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C ----------------------------------------------------------------------
C
c        Calculate global integrals from latitude strips in common
C
      ZHCON  = 0.
      ZHRADS = 0.
      ZHRADS = 0.
      ZHRADL = 0.
      ZRADTP = 0.
      ZSOLIN = 0.
      ZSOLRF = 0.
      ZFIRTP = 0.
C
      ZTAIR   = 0.
      ZTOTCLD = 0.
C
      ZMCLOUF = 0.
      ZMCLOUC = 0.
      ZMCLOUE = 0.
      ZMCLOUS = 0.
      ZHCLOUD = 0.
      ZMREEV = 0.
      ZHREEV = 0.
      ZMPRECC = 0.
      ZMPRECA = 0.
      ZMPRECS = 0.
      ZMCLOUP = 0.
      ZMCLOUA = 0.
      ZHVDIF = 0.
      ZHHDIF = 0.
C
      ZKEVDF = 0.
      ZKEHDF = 0.
C
      ZRFLUX = 0.
      ZSFLUX = 0.
      ZLFLUX = 0.
      ZEFLUX = 0.
      ZPFLUX = 0.
      ZWNEG  = 0.
C
      DO 10 JJ=1,NOREC
        ZHCON  = ZHCON  + DHCON(JJ)*COSBUD(JJ)
        ZHRADS = ZHRADS + DHRADS(JJ)*COSBUD(JJ)
        ZHRADL = ZHRADL + DHRADL(JJ)*COSBUD(JJ)
        ZRADTP = ZRADTP + DRADTP(JJ)*COSBUD(JJ)
        ZSOLIN = ZSOLIN + DSOLIN(JJ)*COSBUD(JJ)
        ZSOLRF = ZSOLRF + DSOLRF(JJ)*COSBUD(JJ)
        ZFIRTP = ZFIRTP + DFIRTP(JJ)*COSBUD(JJ)
C
        ZTAIR   = ZTAIR   + DTAIR(JJ)*COSBUD(JJ)
        ZTOTCLD = ZTOTCLD + DTOTCLD(JJ)*COSBUD(JJ)
C
        ZMCLOUF = ZMCLOUF + DMCLOUF(JJ)*COSBUD(JJ)
        ZMCLOUC = ZMCLOUC + DMCLOUC(JJ)*COSBUD(JJ)
        ZMCLOUE = ZMCLOUE + DMCLOUE(JJ)*COSBUD(JJ)
        ZMCLOUS = ZMCLOUS + DMCLOUS(JJ)*COSBUD(JJ)
        ZHCLOUD = ZHCLOUD + DHCLOUD(JJ)*COSBUD(JJ)
        ZMREEV  = ZMREEV  + DMREEV(JJ)*COSBUD(JJ)
        ZHREEV  = ZHREEV  + DHREEV(JJ)*COSBUD(JJ)
        ZMPRECC = ZMPRECC + DMPRECC(JJ)*COSBUD(JJ)
        ZMPRECA = ZMPRECA + DMPRECA(JJ)*COSBUD(JJ)
        ZMPRECS = ZMPRECS + DMPRECS(JJ)*COSBUD(JJ)
        ZMCLOUP = ZMCLOUP + DMCLOUP(JJ)*COSBUD(JJ)
        ZMCLOUA = ZMCLOUA + DMCLOUA(JJ)*COSBUD(JJ)
        ZHVDIF = ZHVDIF + DHVDIF(JJ)*COSBUD(JJ)
        ZHHDIF = ZHHDIF + DHHDIF(JJ)*COSBUD(JJ)
C
        ZKEVDF = ZKEVDF + DKEVDF(JJ)*COSBUD(JJ)
        ZKEHDF = ZKEHDF + DKEHDF(JJ)*COSBUD(JJ)
C
        ZRFLUX = ZRFLUX + DRFLUX(JJ)*COSBUD(JJ)
        ZSFLUX = ZSFLUX + DSFLUX(JJ)*COSBUD(JJ)
        ZLFLUX = ZLFLUX + DLFLUX(JJ)*COSBUD(JJ)
        ZEFLUX = ZEFLUX + DEFLUX(JJ)*COSBUD(JJ)
        ZPFLUX = ZPFLUX + DPFLUX(JJ)*COSBUD(JJ)
        ZWNEG  = ZWNEG  + DWNEG(JJ) *COSBUD(JJ)
   10 CONTINUE
C
C     Set time normalization factors
C
C     ZD = # of days, ZS = # of secs, ZSTEP = # of timesteps
      ZD = (BTIME-BTIME0)/86400.
      ZS = ZD*86400.
      ZSTEP = ZS/DTIME
C
      WRITE (NBUD,*) CHAR(12)
      WRITE(NBUD,50) BTIME0/86400., BTIME/86400.
   50 FORMAT(
     *  /' ------------------------------------------------------------'
     *  /' *** AGCM BUDGET BETWEEN ELAPSED DAYS',  F10.3,   F10.3,' ***'
     *  /' ------------------------------------------------------------'
     *  )
 
c----------------------------------------------------------------------
 
c     Subtract dhcloud,dhreev from dhvdif since these are included in
c     the reservoir fplumt in stabco,plume,reevap but are not diffusive
c     (conservative) processes. Also subtract DHRADS and DHRADL
c     since radiative heating is included in t-forcing in vdif.
 
      ZERDHV = ZHVDIF - ZHRADS - ZHRADL + (FPLUMT_G-FPLUMT_G0)
     *         - (ZSFLUX + ZHCLOUD + ZHREEV)
      ZERPD= PE+SE-PE0-SE0 - (  ZHRADS+ZHRADL+ZSFLUX
     *                        + ZHCLOUD+ZHREEV
     *                        + ZHCON+ZKEVDF+ZKEHDF)
c                             *approx* correction for d(vapor,lwc mass)
     *                     - (CPWV*(280.-TMELT)+LATSUB)*(WAVA-WAVA0)
     *                     - (                  LATICE)*(WLWC-WLWC0)
      ZERP = PE+SE-PE0-SE0 - (  ZHRADS+ZHRADL+ZSFLUX+ZLFLUX
     *                        + ZHCON+ZKEVDF+ZKEHDF)
C
      WRITE(NBUD,110)
     *  (PE+SE-PE0-SE0)/ZS,
     *  ZHRADS/ZS,
     *  ZHRADL/ZS,
     *  ZSFLUX/ZS,
     *  ZLFLUX/ZS,
     *  ZHCON/ZS,
     *  ZKEVDF/ZS,
     *  ZKEHDF/ZS,
     *  ZHVDIF/ZS,
     *  ZHHDIF/ZS,
     *  ZERDHV/ZS,
     *  ZERPD/ZS,
     *  ZERP/ZS
C
  110 FORMAT(/' AGCM THERMAL+POTENTIAL ENERGY:'
     *       /' ------------------------------'/
     */'  D(PE+SE)   change in thermal+potential energy  (W/m2) ',F24.15
     */'  DHRADS     solar radiative heating rate        (W/m2) ',F24.15
     */'  DHRADL     IR    radiative heating rate        (W/m2) ',F24.15
     */'  DSFLUX     upward surface sensible heat flux   (W/m2) ',F24.15
     */'  DLFLUX     upward surface latent   heat flux   (W/m2) ',F24.15
     */'  DHCON      conversion from KE to PE            (W/m2) ',F24.15
     */'  DKEVDF     vert-diffusion  frictional heating  (W/m2) ',F24.15
     */'  DKEHDF     horiz-diffusion frictional heating  (W/m2) ',F24.15
     */'  DHVDIF     vert-mixing/condens net heating     (W/m2) ',F24.15
     */'  DHHDIF     error in horiz-difn net heating     (W/m2) ',F24.15
     */'  ERRDHV     error in vert-mixing heat budget    (W/m2) ',F24.15
     */'  ERRPDRY    error in "dry" PE budget            (W/m2) ',F24.15
     */'  ERRPTOT    error in total PE budget            (W/m2) ',F24.15
     *)
 
c-----
 
      ZERK = KE-KE0 + (ZHCON+ZKEVDF+ZKEHDF)
 
      WRITE(NBUD,120)
     *  (KE-KE0)/ZS,
     *  -ZHCON/ZS,
     *  -ZKEVDF/ZS,
     *  -ZKEHDF/ZS,
     *  ZERK/ZS
 
  120 FORMAT(/' AGCM KINETIC ENERGY:'
     *       /' --------------------'/
     */'  D(KE)      change in kinetic energy            (W/m2) ',F24.15
     */' -DHCON      conversion from PE to KE            (W/m2) ',F24.15
     */' -DKEVDF     vert-diffusion  KE dissipation      (W/m2) ',F24.15
     */' -DKEHDF     horiz-diffusion KE dissipation      (W/m2) ',F24.15
     */'  ERRKE      error in KE budget                  (W/m2) ',F24.15
     *)
 
c-----
 
      WRITE(NBUD,130)
     *  (PE+SE+KE-PE0-SE0-KE0)/ZS, (ZERPD+ZERK)/ZS, (ZERP+ZERK)/ZS
 
  130 FORMAT(/' AGCM TOTAL ENERGY:'
     *       /' ------------------'/
     */'  D(TOTE)    change in PE plus KE                (W/m2) ',F24.15
     */'  ERRTDRY    error in "dry" PE plus KE budget    (W/m2) ',F24.15
     */'  ERRTOT     error in total PE plus KE budget    (W/m2) ',F24.15
     *)
 
c-----
 
      WRITE(NBUD,210)
     *  (ZHRADS+ZHRADL)/ZS,
     *  ZRADTP/ZS,
     *  -ZRFLUX/ZS,
     *  (ZHRADS+ZHRADL-ZRADTP-ZRFLUX)/ZS
 
  210 FORMAT(/' AGCM RADIATIVE FLUXES, CLOUDINESS, TEMPERATURE'
     *       /' ----------------------------------------------'/
     */'  DHRAD      solar+IR heating rate of AGCM       (W/m2) ',F24.15
     */'  DRADTP     net downward solar+IR at AGCM top   (W/m2) ',F24.15
     */' -DRFLUX     net downward solar+IR at surface    (W/m2) ',F24.15
     */'  ERRRAD     error (DHRAD-DRADTP-DRFLUX)         (W/m2) ',F24.15
     *)
 
      WRITE(NBUD,220)
     *  ZSOLIN/ZS,
     *  ZSOLRF/ZS,
     *  (ZSOLIN-ZSOLRF)/ZS,
     *  ZSOLRF/ZSOLIN,
     *  ((ZSOLIN-ZSOLRF)-(ZRADTP+ZFIRTP))/ZS,
     *  ZFIRTP/ZS,

     *  ZTOTCLD/ZS,
     *  ZTAIR/ZS - TMELT
 
  220 FORMAT(
     */'  DSOLIN     solar insolation      at atmos top  (W/m2) ',F24.15
     */'  DSOLRF     solar reflected       at atmos top  (W/m2) ',F24.15
     */'  DSOLABS    solar insol.-refl.    at atmos top  (W/m2) ',F24.15
     */'  DALBTP     planetary albedo      at atmos top  (W/m2) ',F24.15
     */'  DSOLOZ     solar absorbed     above agcm  top  (W/m2) ',F24.15
     */'  DFIRTP     upward IR flux        at agcm  top  (W/m2) ',F24.15
     */
     */'  DTOTCLD    total cloudiness                    (frac) ',F24.15
     */'  DTAIR      air temperature       lowest level  (C)    ',F24.15
     *)
 
c-----
 
      WRITE(NBUD,310)
     *  WAVA,
     *  WLWC,
     *  (WAVA-WAVA0)/ZD,
     *  (WLWC-WLWC0)/ZD,
     *  ZPFLUX/ZD,
     *  ZEFLUX/ZD,
     *  (WAVA-WAVA0+WLWC-WLWC0)/ZD + (ZPFLUX-ZEFLUX)/ZD
 
  310 FORMAT(/' AGCM HYDROLOGY:'
     *       /' ---------------'/
     */'  WAVA       water vapor mass                   (mm)    ',F24.15
     */'  WLWC       cloud liquid water mass            (mm)    ',F24.15
     */'  D(WAVA)    change in AGCM water vapor mass    (mm/day)',F24.15
     */'  D(WLWC)    change in AGCM cloud liquid mass   (mm/day)',F24.15
     */'  DPFLUX     surface precipitation rate         (mm/day)',F24.15
     */'  DEFLUX     surface evaporation rate           (mm/day)',F24.15
     */'  ERRVAP     error (DWAVA+DWLWC+DPFLUX-DEFLUX)  (mm/day)',F24.15
     *)
 
      WRITE(NBUD,320)
     *  (ZPFLUX-ZEFLUX)/ZD,
     *  ZMPRECC/ZD,
     *  ZMPRECA/ZD,
     *  ZMPRECS/ZD,
     *  (ZMPRECC+ZMPRECA+ZMPRECS)/ZD
 
  320 FORMAT(
     */'  P-E        surface precip minus evaporation   (mm/day)',F24.15
     */'  DMPRECC    surface convective   precipitation (mm/day)',F24.15
     */'  DMPRECA    surface anvil-cirrus precipitation (mm/day)',F24.15
     */'  DMPRECS    surface stratiform   precipitation (mm/day)',F24.15
     */'  DMPREC*    surface total        precipitation (mm/day)',F24.15
     *)
 
      WRITE(NBUD,330)
     *  ZMCLOUF/ZD,
     *  ZMCLOUE/ZD,
     *  ZMCLOUC/ZD,
     *  ZMCLOUP/ZD,
     *  ZMCLOUA/ZD,
     *  ZMREEV/ZD,
     *  ZMCLOUS/ZD,
     *  (WLWC-WLWC0)/ZD - (ZMCLOUF-ZMCLOUE-ZMCLOUC-ZMCLOUS)/ZD,
     *  (ZMCLOUC+ZMCLOUS-ZMREEV-ZPFLUX)/ZD
 
  330 FORMAT(
     */'  DMCLOUF    vapor-to-cloud condensation rate   (mm/day)',F24.15
     */'  DMCLOUE    cloud-to-vapor evaporation rate    (mm/day)',F24.15
     */'  DMCLOUC    cloud-to-precip conversion rate    (mm/day)',F24.15
     */'  DMCLOUP    part of DMCLOUC within plumes      (mm/day)',F24.15
     */'  DMCLOUA    part of DMCLOUC due to accretion   (mm/day)',F24.15
     */'  DMREEV     evaporation of falling precip rate (mm/day)',F24.15
     */'  DMCLOUS    cloud-to-surface direct deposition (mm/day)',F24.15
     */'  ERRLWC     error (DWLWC-DMCLOU[F-E-C-S])      (mm/day)',F24.15
     */'  ERRPREC    error (DMCLOU[C+S]-DMREEV-DPFLUX)  (mm/day)',F24.15
     *)
 
      WRITE(NBUD,340)
     *  ZHCLOUD/ZS,
     *  ZHREEV/ZS
 
  340 FORMAT(
     */'  DHCLOUD    all cloud processes                (W/m2)  ',F24.15
     */'  DHREEV     evaporation of falling precip rate (W/m2)  ',F24.15
     *)
 
c-----
 
      WRITE(NBUD,400)
     *  -(ZRFLUX + ZSFLUX + ZLFLUX)/ZS,
     *  FPLUMT_G/PE,
     *  FPLUMQ_G/WAVA
 
  400 FORMAT(/' AGCM MISC INFORMATION:'
     *       /' ----------------------'/
     */'  SURFLX     net downward surface heat flux     (W/m2)  ',F24.15
     */'  FPLUMT/PE  heat reservoir / total heat        (ratio) ',F24.15
     */'  FPLUMQ/WA  vapor reservoir / total vapor      (ratio) ',F24.15
     *)
 
      WRITE(NBUD,410)
     *  DWFIX/ZD,
     *  DWFIXA/ZD,
     *  ZWNEG/ZD,
     *  DMASS/ZSTEP,
     *  DMASSA/ZSTEP
  410 FORMAT(
     */'  DWFIX      slt h2o conservation adjust        (mm/day)',F24.15
     */'  DWFIXA     abs [slt h2o conservation adjust]  (mm/day)',F24.15
     */'  DWNEG      negative h2o fix (vdif)            (mm/day)',F24.15
     */'  DMASS      dry mass adjustment                (Pa/stp)',F24.15
     */'  DMASSA     abs [dry mass adjustment]          (Pa/stp)',F24.15
     *)
 
C----------------------------------------------------------------------
 
CL            5.0       PRINT THE MINIMUM AND MAXIMUM TEMPERATURES,
CL                      AND THE TIMES OF OCCURENCE
C
      WRITE(NBUD,500) TMPMIN, NINT(MINTLO),NINT(MINTLA),NINT(MINTLV),
     *                ITMMIN*DTIME/86400.,
     *                TMPMAX, NINT(MAXTLO),NINT(MAXTLA),NINT(MAXTLV),
     *                ITMMAX*DTIME/86400.
C
  500 FORMAT(//' Minimum temperature for this budget period = ',f24.15,
     *        /'    at (i,j,k) = (',i3,',',i3,',',i3,')',
     *         '    on day =',f10.3
     *        /' Maximum temperature for this budget period = ',f24.15,
     *        /'    at (i,j,k) = (',i3,',',i3,',',i3,')',
     *         '    on day =',f10.3 )
C----------------------------------------------------------------------
C
      WRITE(NBUD,600)
  600 FORMAT(/' ----------------------------------------------'/)
C
      RETURN
      END
C
C ***************************** FLDLST ********************************
C
      SUBROUTINE FLDLST
C----------------------------------------------------------------------
C
CL            INITIALIZE HARD-CODED PROGRAM LIST OF ALL HISTORY FIELDS
C
C---------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMHIS.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C---------------------------------------------------------------------
      CHARACTER*8 COSMONAME(12)
      DATA COSMONAME
     *  /'M0S3    ', 'M0S6    ','M1S3    ', 'M1S6    ', 'M5S3    ',
     *   'M5S6    ',
     *   '7M0S3   ', '7M0S6   ','7M1S3   ', '7M1S6   ', '7M5S3   ',
     *   '7M5S6   '/
C---------------------------------------------------------------------
 
C
      DO 10 I=1,NUSHMX
        MINDEX(1,I) = 0
        MINDEX(2,I) = 0
        MINDEX(3,I) = 0
        MINDEX(4,I) = 0
   10 CONTINUE
C
      IP = 0
C
C01            TOPOG           SURFACE TOPOGRAPHY (FILTERED, TRUNCATED)
      IP = IP + 1
      CINDEX(1,IP) = 'TOPOG'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C02            PS              SURFACE PRESSURE
      IP = IP + 1
      CINDEX(1,IP) = 'PS'
      CINDEX(2,IP) = 'Pa'
      MINDEX(2,IP) = 1
C03            T               TEMPERATURE
      IP = IP + 1
      CINDEX(1,IP) = 'T'
      CINDEX(2,IP) = 'deg K'
      MINDEX(2,IP) = NLEV
C04            U               EAST-WEST WIND
      IP = IP + 1
      CINDEX(1,IP) = 'U'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = NLEV
C05            V               NORTH-SOUTH WIND
      IP = IP + 1
      CINDEX(1,IP) = 'V'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = NLEV
C06            Q               SPECIFIC HUMIDITY
      IP = IP + 1
      CINDEX(1,IP) = 'Q'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C07            DUH             U HORIZONTAL DIFFUSIVE HEATING
      IP = IP + 1
      CINDEX(1,IP) = 'DUH'
      CINDEX(2,IP) = 'deg K/s'
      MINDEX(2,IP) = NLEV
C08            DVH             V HORIZONTAL DIFFUSIVE HEATING
      IP = IP + 1
      CINDEX(1,IP) = 'DVH'
      CINDEX(2,IP) = 'deg K/s'
      MINDEX(2,IP) = NLEV
C09            DTH             T HORIZONTAL DIFFUSION
      IP = IP + 1
      CINDEX(1,IP) = 'DTH'
      CINDEX(2,IP) = 'deg K/s'
      MINDEX(2,IP) = NLEV
C10            TOPOGUN          SURFACE TOPOG, UNFILTERED, UNTRUNCATED
      IP = IP + 1
      CINDEX(1,IP) = 'TOPOGUN'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C11            ORO             CCM SURFACE TYPE (DERIVED FROM LSX)
      IP = IP + 1
      CINDEX(1,IP) = 'ORO'
      CINDEX(2,IP) = ' '
      MINDEX(2,IP) = 1
C12            PRECL           STABLE PRECIPITATION RATE
      IP = IP + 1
      CINDEX(1,IP) = 'PRECL'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C13            PRECC           CONVECTIVE PRECIPITATION RATE
      IP = IP + 1
      CINDEX(1,IP) = 'PRECC'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C14            RAINF           RAINFALL RATE (LIQ EQUIV DEPTH)
      IP = IP + 1
      CINDEX(1,IP) = 'RAINF'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C15            SNOWF           SNOWFALL RATE (LIQ EQUIV DEPTH)
      IP = IP + 1
      CINDEX(1,IP) = 'SNOWF'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C16            HFLR            SENSIBLE HEAT FLUX RATE
      IP = IP + 1
      CINDEX(1,IP) = 'HFLR'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C17            QFLR            LATENT HEAT FLUX RATE
      IP = IP + 1
      CINDEX(1,IP) = 'QFLR'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C18            DUS             HEATING FROM E-W SURFACE STRESS
      IP = IP + 1
      CINDEX(1,IP) = 'DUS'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C19            DVS             HEATING FROM N-S SURFACE STRESS
      IP = IP + 1
      CINDEX(1,IP) = 'DVS'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C20            Z200            200 mb GEOPOTENTIAL HEIGHT
      IP = IP + 1
      CINDEX(1,IP) = 'Z200'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C21            Z500            500 mb GEOPOTENTIAL HEIGHT
      IP = IP + 1
      CINDEX(1,IP) = 'Z500'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C22            Z850            850 mb GEOPOTENTIAL HEIGHT
      IP = IP + 1
      CINDEX(1,IP) = 'Z850'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C23            FRLD            INCIDENT LONGWAVE FLUX ABSORBED AT SFC
      IP = IP + 1
      CINDEX(1,IP) = 'FRLD'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C24            FRSA            NET DOWNWARD SOLAR FLUX AT SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'FRSA'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C25            FRLA            NET UPWARD LONGWAVE FLUX AT SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'FRLA'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C26            FIRTP           NET UPWARD LONGWAVE FLUX AT MODEL TOP
      IP = IP + 1
      CINDEX(1,IP) = 'FIRTP'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C27            SABTP           ABSORBED SOLAR FLUX BY MODEL ATMOS + SFC
      IP = IP + 1
      CINDEX(1,IP) = 'SABTP'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C28            SABAT           ABSORBED SOLAR FLUX ABOVE TOP OF MODEL
      IP = IP + 1
      CINDEX(1,IP) = 'SABAT'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C29            CLOUD           3-D FRACTIONAL TOTAL CLOUD AMOUNT
      IP = IP + 1
      CINDEX(1,IP) = 'CLOUD'
      CINDEX(2,IP) = 'fraction'
      MINDEX(2,IP) = NLEV
C30            PRECIP          PRECIPITATION RATE (LIQ EQUIV DEPTH)
      IP = IP + 1
      CINDEX(1,IP) = 'PRECIP'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C31            CLRLT           CLEAR-SKY LONGWAVE FLUX AT TOP OF ATMOS
      IP = IP + 1
      CINDEX(1,IP) = 'CLRLT'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C32            CLRST           CLEAR-SKY SOLAR FLUX AT TOP OF ATMOS
      IP = IP + 1
      CINDEX(1,IP) = 'CLRST'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C33            CLRLS           CLEAR-SKY LONGWAVE FLUX AT SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'CLRLS'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C34            CLRSS  (=65)    CLEAR-SKY SOLAR FLUX AT SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'CLRSS'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C35            OMEGA           VERTICAL PRESSURE VELOCITY (DP/DT)
      IP = IP + 1
      CINDEX(1,IP) = 'OMEGA'
      CINDEX(2,IP) = 'Pa/s'
      MINDEX(2,IP) = NLEV
C36            SIGDOT          VERTICAL COORDINATE VELOCITY (DSIG/DT)
      IP = IP + 1
      CINDEX(1,IP) = 'SIGDOT'
      CINDEX(2,IP) = 's-1'
      MINDEX(2,IP) = NLEV
C37            TAUX            EAST-WEST WIND STRESS ON SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'TAUX'
      CINDEX(2,IP) = 'N/m2'
      MINDEX(2,IP) = 1
C38            TAUY            NORTH-SOUTH WIND STRESS ON SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'TAUY'
      CINDEX(2,IP) = 'N/m2'
      MINDEX(2,IP) = 1
C39            SRAD            NET DOWNWARD RADIATION FLUX AT SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'SRAD'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C40            QRS             SHORTWAVE ATMOS DOWN-UP FLUX
      IP = IP + 1
      CINDEX(1,IP) = 'QRS'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = NLEV+1
C41            QRL             LONGWAVE ATMOS DOWN-UP FLUX 
      IP = IP + 1
      CINDEX(1,IP) = 'QRL'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = NLEV+1
C42            SOLIN           SOLAR INSOLATION
      IP = IP + 1
      CINDEX(1,IP) = 'SOLIN'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C43            TBOT            AGCM LOWEST-LEVEL TEMPERATURE
      IP = IP + 1
      CINDEX(1,IP) = 'TBOT'
      CINDEX(2,IP) = 'deg K'
      MINDEX(2,IP) = 1
C44            TOTCLD          TOTAL CLOUDINESS (RANDOM OVERLAP)
      IP = IP + 1
      CINDEX(1,IP) = 'TOTCLD'
      CINDEX(2,IP) = 'fraction'
      MINDEX(2,IP) = 1
C45            TOTCLDC         TOTAL CONVECTIVE CLOUDINESS
      IP = IP + 1
      CINDEX(1,IP) = 'TOTCLDC'
      CINDEX(2,IP) = 'fraction'
      MINDEX(2,IP) = 1
C46            TOTCLDA         TOTAL ANVIL-CIRRUS CLOUDINESS
      IP = IP + 1
      CINDEX(1,IP) = 'TOTCLDA'
      CINDEX(2,IP) = 'fraction'
      MINDEX(2,IP) = 1
C47            TOTCLDS         TOTAL STRATUS CLOUDINESS
      IP = IP + 1
      CINDEX(1,IP) = 'TOTCLDS'
      CINDEX(2,IP) = 'fraction'
      MINDEX(2,IP) = 1
C48            WPLUME0         "PBL" PLUME VERTICAL VELOCITY
      IP = IP + 1
      CINDEX(1,IP) = 'WPLUME0'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = NLEV
C49            WPLUME1         "FREE" PLUME VERTICAL VELOCITY
      IP = IP + 1
      CINDEX(1,IP) = 'WPLUME1'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = NLEV
C50            FRSAV           SOLAR FLUX ABSORBED AT SFC, VIS WAVEBAND
      IP = IP + 1
      CINDEX(1,IP) = 'FRSAV'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C51            FRSIV           SOLAR FLUX INCIDENT AT SFC, VIS WAVEBAND
      IP = IP + 1
      CINDEX(1,IP) = 'FRSIV'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C52            FRSAN           SOLAR FLUX ABSORBED AT SFC, NIR WAVEBAND
      IP = IP + 1
      CINDEX(1,IP) = 'FRSAN'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C53            DQASLT          CHANGE IN Q FROM SLT ADVECTION
      IP = IP + 1
      CINDEX(1,IP) = 'DQASLT'
      CINDEX(2,IP) = 'kg/kg/s'
      MINDEX(2,IP) = NLEV
C54            CONVLWC         3-D CONVECTIVE CLOUD LIQUID WATER CONTENT
      IP = IP + 1
      CINDEX(1,IP) = 'CONVLWC'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C55            ANVILWC         3-D ANVIL-CIRRUS CLOUD LIQ. WATER CONTENT
      IP = IP + 1
      CINDEX(1,IP) = 'ANVILWC'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C56            FRSIN           SOLAR FLUX INCIDENT AT SFC, NIR WAVEBAND
      IP = IP + 1
      CINDEX(1,IP) = 'FRSIN'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C57            STRALWC         3-D STRATUS CLOUD LIQUID WATER CONTENT
      IP = IP + 1
      CINDEX(1,IP) = 'STRALWC'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C58            PSLEV           SEA-LEVEL PRESSURE
      IP = IP + 1
      CINDEX(1,IP) = 'PSLEV'
      CINDEX(2,IP) = 'Pa'
      MINDEX(2,IP) = 1
C59            LWC             3-D TOTAL CLOUD LIQUID WATER CONTENT
      IP = IP + 1
      CINDEX(1,IP) = 'LWC'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C60            SURWIND         WIND SPEED AT LOWEST AGCM LEVEL
      IP = IP + 1
      CINDEX(1,IP) = 'SURWIND'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C61            RHUMID          RELATIVE HUMIDITY
      IP = IP + 1
      CINDEX(1,IP) = 'RHUMID'
      CINDEX(2,IP) = 'fraction'
      MINDEX(2,IP) = NLEV
C62            TOTLWC          TOTAL CLOUD LIQUID WATER CONTENT
      IP = IP + 1
      CINDEX(1,IP) = 'TOTLWC'
      CINDEX(2,IP) = 'kg/m2'
      MINDEX(2,IP) = 1
C63            PBL03           PBL HEIGHT AT 03:00 LOCAL TIME
      IP = IP + 1
      CINDEX(1,IP) = 'PBL03'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C64            PBL15           PBL HEIGHT AT 15:00 LOCAL TIME
      IP = IP + 1
      CINDEX(1,IP) = 'PBL15'
      CINDEX(2,IP) = 'm'
      MINDEX(2,IP) = 1
C65            FRSAC  (=34)    NET DWNWD CLEAR-SKY SOLAR FLX AT SURFACE
      IP = IP + 1
      CINDEX(1,IP) = 'FRSAC'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C66            TAIR            TEMPERATURE (SAME AS T, FOR FERRET)
      IP = IP + 1
      CINDEX(1,IP) = 'TAIR'
      CINDEX(2,IP) = 'deg K'
      MINDEX(2,IP) = NLEV
#if defined (LKAEROSOL)
C67            SFCIRRAD
      IP = IP + 1
      CINDEX(1,IP) = 'SFCIRRAD'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
#else
C67            spare
      IP = IP + 1
      CINDEX(1,IP) = 'spare'
      CINDEX(2,IP) = 'spare'
      MINDEX(2,IP) = 1
#endif
C68            spare
      IP = IP + 1
      CINDEX(1,IP) = 'spare'
      CINDEX(2,IP) = 'spare'
      MINDEX(2,IP) = 1

c*************
#ifdef isotrac
c*************
C69            Q_O             SPECIFIC HUMIDITY - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'Q_O'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C70            Q_H             SPECIFIC HUMIDITY - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'Q_H'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C71            CONVLW_O        3-D CONV. CLOUD LIQ. WATER CONTENT - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'CONVLW_O'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C72            CONVLW_H        3-D CONV. CLOUD LIQ. WATER CONTENT - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'CONVLW_H'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C73            ANVILW_O        3-D ANVIL-CIRRUS LIQ. WATER CONT. - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'ANVILW_O'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C74            ANVILW_H        3-D ANVIL-CIRRUS LIQ. WATER CONT. - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'ANVILW_H'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C75            STRALW_O        3-D STRATUS CLOUD LIQ. WATER CONT. - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'STRALW_O'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C76            STRALW_H        3-D STRATUS CLOUD LIQ. WATER CONT. - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'STRALW_H'
      CINDEX(2,IP) = 'kg/kg'
      MINDEX(2,IP) = NLEV
C77            QFLR_O          LATENT HEAT FLUX RATE - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'QFLR_O'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C78            QFLR_H          LATENT HEAT FLUX RATE - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'QFLR_H'
      CINDEX(2,IP) = 'W/m2'
      MINDEX(2,IP) = 1
C79            PRECL_O         STABLE PRECIPITATION RATE - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'PRECL_O'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C80            PRECL_H         STABLE PRECIPITATION RATE - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'PRECL_H'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C81            PRECC_O         CONVECTIVE PRECIPITATION RATE - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'PRECC_O'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C82            PRECC_H         CONVECTIVE PRECIPITATION RATE - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'PRECC_H'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C83            RAINF_O         RAINFALL RATE - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'RAINF_O'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C84            RAINF_H         RAINFALL RATE - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'RAINF_H'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C85            SNOWF_O         SNOWFALL RATE - 18O
      IP = IP + 1
      CINDEX(1,IP) = 'SNOWF_O'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C86            SNOWF_H         SNOWFALL RATE - 2H
      IP = IP + 1
      CINDEX(1,IP) = 'SNOWF_H'
      CINDEX(2,IP) = 'm/s'
      MINDEX(2,IP) = 1
C87-89         spare
      DO 870 I=87,89
        IP = IP + 1
        CINDEX(1,IP) = 'spare'
        CINDEX(2,IP) = 'spare'
        MINDEX(2,IP) = 1
  870 CONTINUE
c****
#else
c****
C69-89         spare
      DO 690 I=69,89
        IP = IP + 1
        CINDEX(1,IP) = 'spare'
        CINDEX(2,IP) = 'spare'
        MINDEX(2,IP) = 1
  690 CONTINUE
c*****
#endif
c*****

c90            UVFLUX          3-D U*V FLUX
      ip = ip + 1
      cindex(1,ip) = 'UVFLUX'
      cindex(2,ip) = 'm2/s2'
      mindex(2,ip) = nlev
c91            WUFLUX        3-D OMEGA*U FLUX
      ip = ip + 1
      cindex(1,ip) = 'WUFLUX'
      cindex(2,ip) = 'Pa m/s2'
      mindex(2,ip) = nlev
c92            WTFLUX        3-D OMEGA*T FLUX
      ip = ip + 1
      cindex(1,ip) = 'WTFLUX'
      cindex(2,ip) = 'Pa K/s'
      mindex(2,ip) = nlev
c93            WVOFLUX        3-D OMEGA*VORTICITY FLUX
      ip = ip + 1
      cindex(1,ip) = 'WVOFLUX'
      cindex(2,ip) = 'Pa/s2'
      mindex(2,ip) = nlev
c94            WPVOFLUX        3-D OMEGA* POTENTIAL VORTICITY FLUX
      ip = ip + 1
      cindex(1,ip) = 'WPVOFLUX'
      cindex(2,ip) = 'Pa/s Ert'
      mindex(2,ip) = nlev
C95-100        spare
      DO 950 I=95,100
        IP = IP + 1
        CINDEX(1,IP) = 'spare'
        CINDEX(2,IP) = 'spare'
        MINDEX(2,IP) = 1
  950 CONTINUE

c>>>>>>>>>>>>>>>>
#ifdef isotracreg
c>>>>>>>>>>>>>>>>
C101-1??       RTR_Qn           3-D REG. TRAC. SPEC. HUMID. #1 TO #NRTR
      do it=1,nrtr
        ip=ip+1
        if (it.le.9) write(cindex(1,ip),'(a,i1)') 'RTR_Q',it
        if (it.gt.9) write(cindex(1,ip),'(a,i2)') 'RTR_Q',it
        cindex(2,ip) = 'kg/kg'
        mindex(2,ip) = nlev
      enddo
C1??-1??       RTR_Pn           REGIO. TRACERS PRECIP. RATES #1 TO #NRTR
      do it=1,nrtr
        ip=ip+1
        if (it.le.9) write(cindex(1,ip),'(a,i1)') 'RTR_P',it
        if (it.gt.9) write(cindex(1,ip),'(a,i2)') 'RTR_P',it
        cindex(2,ip) = 'm/s'
        mindex(2,ip) = 1
      enddo
c>>>>
#else
c>>>>
C101-140       TRACERn          3-D TRUE TRACERS #1 TO #NTRACE
      IHTRACE = IP + 1          ! used in linems
      DO 1010 N=1,40
        IP = IP + 1
#if defined ( cosmog )
        IF (N.LE.NTRACB) THEN
          CINDEX(1,IP) = COSMONAME(N)
          CINDEX(2,IP) = 'ats/kg'
#else
        IF (N.LE.NTRACB) THEN
          IF (N.LE.9) WRITE(CINDEX(1,IP),'(A,I1)') 'TRACER',N
          IF (N.GT.9) WRITE(CINDEX(1,IP),'(A,I2)') 'TRACER',N
          CINDEX(2,IP) = 'kg/kg'
#endif
        ELSE
          CINDEX(1,IP) = 'spare'
          CINDEX(2,IP) = 'spare'
        ENDIF
        MINDEX(2,IP) = NLEV
 1010 CONTINUE
 
C141-180       TRATOTn         VERT-INTEG. TRUE  TRACERS #1 TO #NTRACE
C                              OR SURFACE FLUX FOR COSMOGENIC TRACERS
      DO 1210 N=1,40
        IP = IP + 1
#if defined ( cosmog )
        IF (N.LE.NTRACB) THEN
          LENC = MIN (LENCHR(COSMONAME(N)), 7)
          CINDEX(1,IP) = COSMONAME(N)(1:LENC) // 'FLX'
          CINDEX(2,IP) = 'ats/m2/s'
#else
        IF (N.LE.NTRACB) THEN
          IF (N.LE.9) WRITE(CINDEX(1,IP),'(A,I1)') 'TRATOT',N
          IF (N.GT.9) WRITE(CINDEX(1,IP),'(A,I2)') 'TRATOT',N
          CINDEX(2,IP) = 'kg/m2'
#endif
        ELSE
          CINDEX(1,IP) = 'spare'
          CINDEX(2,IP) = 'spare'
        ENDIF
        MINDEX(2,IP) = 1
 1210 CONTINUE
c>>>>>
#endif
c>>>>>
 
      NINDEX = IP
      IF (NINDEX.GT.NUSHMX) THEN
         WRITE (NOUT,9005) NINDEX, NUSHMX
 9005    FORMAT(' NOT ENOUGH ROOM IN INDICES CINDEX, MINDEX FOR',
     *          ' PROGRAM HISTORY FIELD LIST.'
     *         /' PROGRAM ERROR IN FLDLST. NINDEX,NUSHMX =',2I10)
         CALL ENDRUN (-1)
      ENDIF
C
      RETURN
      END
C
C ***************************** HISSET ********************************
C
      SUBROUTINE HISSET(IFLAG)
C----------------------------------------------------------------------
C
CL            SET RECORD-COUNTER NUOFF AND WRITE-COUNTER NWRIT,
CL            AND ZERO ACCUMULATION COUNTERS MINDCO, FOR AGCM HIST FILE.
CL
CL            IF IFLAG = 1,  STARTING NEW FILE.
CL            IF IFLAG = 0,  ONGOING FILE, INCREMENT FOR NEXT WRITE.
CL            IF IFLAG =-1,  APPEND TO OLD FILE (NUOFF,NWRIT HAVE
CL                           BEEN READ FROM RESTART FILE).
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMHIS.h>
C ---------------------------------------------------------------------
#include <COMLUN.h>
C ---------------------------------------------------------------------
C
      IF (IFLAG.EQ.1) THEN
         NUOFF = 0
         NWRIT = 0
      ELSE IF (IFLAG.EQ.0) THEN
         NUOFF = NUOFF + NULHDR + NULDAT
         NWRIT = NWRIT + 1
      ENDIF
C
      DO 10 I=1,NUSHMX
        MINDCO(I) = 0
   10 CONTINUE
C
      RETURN
      END
C
C ***************************** PUTHIS ********************************
C
      SUBROUTINE PUTHIS (ARR, IV)
C
C----------------------------------------------------------------------
C
CL            WRITES A HISTORY FIELD IN ARRAY ARR TO THE AGCM HISTORY
CL            STORAGE AREA. IV IS THE PROGRAM'S FIELD INDEX. EACH
CL            "RECORD" IN THE STORAGE AREA HOLDS NLON*NOREC VALUES (FOR
CL            ONE  LEVEL), AND THE STORAGE AREA (ACCESSED BY GETHR AND
CL            PUTHR) HOLDS DATA FOR ONE HISTORY "WRITE". ARRAY ARR IS
CL            NOT CHANGED.
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMP2.h>
C----------------------------------------------------------------------
#include <COMSDT.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      DIMENSION ARR(NLON,NOREC,NLEVP)
      DIMENSION TMPX(NLON,NOREC), TMPAV(NLON,NOREC)
C----------------------------------------------------------------------
C
C     DO NOTHING UNLESS USER HAS REQUESTED THIS FIELD
      IF(MINDEX(3,IV).EQ.0) RETURN
C
C     DO NOTHING IF NOT IN ACCUMULATION WINDOW
      IF (.NOT.FAHIST) RETURN
C
C     DO NOTHING FOR INSTANTANEOUS FIELD AND NOT AT WRITE TIMESTEP)
      IF(MINDFL(IV).EQ.0 .AND. .NOT.FWHIST) RETURN
C
C     INCREMENT COUNTER FOR ACCUMULATED FIELDS
      IF(MINDFL(IV).EQ.1) MINDCO(IV) = MINDCO(IV) + 1
C
      DO 100 JK = 1, MINDEX(2,IV)
C
         IREC = MINDEX(3,IV) + JK -1
C
C        FOR NON-FIRST-TIME ACCUMULATED FIELDS, READ EXISTING VALUES
         IF(MINDFL(IV).EQ.1 .AND. MINDCO(IV).GT.1) THEN
            CALL GETHR (IREC, TMPAV, IVDUM)
         ENDIF
C
C        COPY CURRENT FIELD TO TEMPORARY ARRAY
         CALL SCOPY(NLON*NOREC, ARR(1,1,JK), 1, TMPX, 1)
C
C        FOR NON-FIRST-TIME ACCUMULATED FIELDS, ADD EXISTING VALUES
         IF(MINDFL(IV).EQ.1 .AND. MINDCO(IV).GT.1) THEN
            DO 60 JJ=1,NOREC
               DO 62 JI=1,NLON
                  TMPX(JI,JJ) = TMPAV(JI,JJ) + TMPX(JI,JJ)
   62          CONTINUE
   60       CONTINUE
         ENDIF
C
         CALL PUTHR (IREC, TMPX, IV)
C
  100 CONTINUE
C
      RETURN
      END
C
C ***************************** TESTHIST ******************************
C
      LOGICAL FUNCTION TESTHIST (IV, ICODE)
C
C----------------------------------------------------------------------
C        Test if need to calculate a history field (program index IV)
C        If so (ie, need to call PUTHIS), returns true.
C          If ICODE=0, only test if user-requested (indep of timestep).
C          If ICODE=1, also test if need to calculate for this timestep.
C        The distinction is useful since some fields (eg, radiation)
C        are only calculated if FRADSW/LW are true, which may not
C        coincide with user's instantaneous history intervals.
C        This is a "shell" of puthis, with "testhist = .true." replacing
C        any actual puthis operations.
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
 
      TESTHIST = .FALSE.
 
c     Not a requested field
      IF(MINDEX(3,IV).EQ.0) RETURN
 
c        Only do remaining timestep-specific tests if ICODE ne 0
 
      IF (ICODE.NE.0) THEN
c        No accum/write this timestep
         IF (.NOT.FAHIST) RETURN
 
c        Instantaneous field and no write this timestep
         IF(MINDFL(IV).EQ.0 .AND. .NOT.FWHIST) RETURN
      ENDIF
 
      TESTHIST = .TRUE.
 
      RETURN
      END
C
C ***************************** GETHR *********************************
C
      SUBROUTINE GETHR (IREC, ARR, IV)
C----------------------------------------------------------------------
C
CL            GET A HISTORY FIELD FROM THE HISTORY STORAGE AREA.
CL            (HISSAV ARRAY IN COMHIZ).
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      DIMENSION ARR(NLON*NOREC)
C----------------------------------------------------------------------
#if defined ( isotrac )
c     PARAMETER (NHIREC = 300+20*(18+1))! =678 20(3-D tracers + precip.)
      PARAMETER (NHIREC = 60 + 10*NLEV) 
#else
c     PARAMETER (NHIREC = 200)
      PARAMETER (NHIREC = 50 + 8*NLEV )
#endif
      COMMON /COMHIZ/ HISSAV(NLON*NOREC,NHIREC), INDIVSAV(NHIREC)

C----------------------------------------------------------------------
 
      DO 200 I=1,NLON*NOREC
        ARR(I) = HISSAV(I,IREC)
  200 CONTINUE

      iv = indivsav(irec)
 
      RETURN
      END
C
C ***************************** PUTHR *********************************
C
      SUBROUTINE PUTHR (IREC, ARR, iv)
C----------------------------------------------------------------------
C
CL            PUT A HISTORY FIELD INTO THE HISTORY STORAGE AREA.
CL            (HISSAV ARRAY IN COMHIZ).
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      DIMENSION ARR(NLON*NOREC)
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST
C----------------------------------------------------------------------

#if defined ( isotrac )
c     PARAMETER (NHIREC = 300+20*(18+1))! =678 20(3-D tracers + precip.)
      PARAMETER (NHIREC = 60 + 10*NLEV) 
#else
c     PARAMETER (NHIREC = 200)
      PARAMETER (NHIREC = 50 + 8*NLEV )
#endif
      COMMON /COMHIZ/ HISSAV(NLON*NOREC,NHIREC), INDIVSAV(NHIREC)
C----------------------------------------------------------------------
 
      IF (FIRST) THEN
        IF (NULDAT.GT.NHIREC) THEN
          WRITE(NOUT,2000) NULDAT, NHIREC
 2000     FORMAT(/' *** The amount of data for one AGCM history write',
     *            ' is too large for the program buffer.'
     *           /'     Number of data records =',i5
     *           /'     Program maximum        =',i5 /)
          CALL ENDRUN (-1)
        ENDIF
        do i=1,nhirec
          indivsav(i) = 0
        enddo
        FIRST = .FALSE.
      ENDIF
 
      DO 200 I=1,NLON*NOREC
        HISSAV(I,IREC) = ARR(I)
  200 CONTINUE

      indivsav(irec) = iv
 
      RETURN
      END
C
C ***************************** THILO *********************************
C
      SUBROUTINE THILO
C----------------------------------------------------------------------
C
CL            FINDS MIN AND MAX TEMPERATURES FOR BUDGET REPORTING.
CL            CALLED FROM NEAR END OF LINEMS FOR FILTERED TIME N,
CL            BEFORE BUDGET DONE. ALSO WRITES WARING IF TEMPERATURE
CL            EXTREMES ARE OUTSIDE [100 to 400] DEG K.
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      DIMENSION ZT(NLON*NOREC*NLEV)
      EQUIVALENCE (ZT,YTM1)
C----------------------------------------------------------------------
C
      II = ISMIN (NLON*NOREC*NLEV, ZT, 1)
      IF (ZT(II) .LT. MAX(TMPMIN,100.)) THEN
         ILEV = (II-1)/(NLON*NOREC) + 1
         IIR = II - (NLON*NOREC)*(ILEV-1)
         ILAT = (IIR-1)/NLON + 1
         ILON = IIR - NLON*(ILAT-1)
         IF (ZT(II) .LT. TMPMIN) THEN
           TMPMIN = ZT(II)
           ITMMIN = NSTEP
           MINTLV = ILEV
           MINTLA = ILAT
           MINTLO = ILON
         ENDIF
         IF (ZT(II) .LT. 100.) THEN
           WRITE (NOUT,500) 'minimum',
     *       ZT(II), ILON, ILAT, ILEV, NSTEP*DTIME/86400.
  500      FORMAT (/' *** Warning: ',a,' air temperature = ',f9.2,
     *              ' deg K  at (i,j,k,day) = ',3i4,1x,f9.3)
         ENDIF
      ENDIF
C
      II = ISMAX (NLON*NOREC*NLEV, ZT, 1)
      IF (ZT(II) .GT. MIN(TMPMAX,400.)) THEN
         ILEV = (II-1)/(NLON*NOREC) + 1
         IIR = II - (NLON*NOREC)*(ILEV-1)
         ILAT = (IIR-1)/NLON + 1
         ILON = IIR - NLON*(ILAT-1)
         IF (ZT(II) .GT. TMPMAX) THEN
           TMPMAX = ZT(II)
           ITMMAX = NSTEP
           MAXTLV = ILEV
           MAXTLA = ILAT
           MAXTLO = ILON
         ENDIF
         IF (ZT(II) .GT. 400.) THEN
           WRITE (NOUT,500) 'maximum',
     *       ZT(II), ILON, ILAT, ILEV, NSTEP*DTIME/86400.
         ENDIF
      ENDIF
C
      RETURN
      END
C
C ***************************** WRITH **********************************
C
      subroutine writh
 
c-----------------------------------------------------------------------
c
c If HISNET.EQ..F., (Namelist, in comlun), writes a history "write", 
c ie, one set of header recs and data recs, to the current agcm 
c history (binary stream) file.

c If HISNET.EQ..T., does same but to the current agcm *netCDF*
c history file. Uses variables (cusfielnet, numfielnet, nwrit) 
c in common (COMLUN.h). Adapted from convnet, August '00.
c These files are opened in openf, closed in closef (as for binary).

c-----------------------------------------------------------------------
#include <PARAMS.h>
c-----------------------------------------------------------------------
#include <COMFLG.h>
c-----------------------------------------------------------------------
#include <COMHED.h>
c-----------------------------------------------------------------------
#include <COMHIS.h>
c-----------------------------------------------------------------------
#include <COMLUN.h>
c-----------------------------------------------------------------------
#include <COMMAP.h>
c-----------------------------------------------------------------------
#include <COMSDT.h>
c-----------------------------------------------------------------------
#include <COMTIM.h>
c-----------------------------------------------------------------------
#include <COMTRU.h>
c-----------------------------------------------------------------------
#include <CRDCTL.h>
c-----------------------------------------------------------------------
      dimension arr(nlon*norec), brr(nlon*norec)
      character crec*8, runtitlx*80
      character cdum*8, cdum16*16
      save cdum,cdum16
      data cdum,cdum16 /' ',' '/
c-----------------------------------------------------------------------
c for netCDF:
#ifdef NETCDF
#include <netcdf.inc>
#ifdef cray
      integer
#else
      integer*4 
#endif
     * lenattr, nlon4, norec4, nlev4, ival4,     ! ncid in comlun
     * varid, rcode,
     * lonid, latid, levid, timid, minmaxid,
     * minmax(2), ndim, idim(4), start(4), count(4),
     * n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls. 
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)

      save varid, rcode,
     *     lonid, latid, levid, timid, minmaxid,
     *     minmax, idim, start, count,
     *     n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
      data minmax /1,2/
      data n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
     *     /1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13/
      character*8 choriz
      dimension ztmp(1000)
      parameter (pi = 3.14159265358979)
#endif
c-----------------------------------------------------------------------
c
c        If requested, put current file counter at end of runtitle
c
      runtitlx = runtitle
      i = lenchr(runtitle)
c     if (runtitle(i-2:i).eq.'...') write(runtitlx(i-2:i),'(i3.3)')ifile
      write(runtitlx(i+2:i+4),'(i3.3)') ifile

c==========================
      if (.not.hisnet) then
c==========================

c        Regular binary files
 
c---------------------------------
c        Write "CCM" header record
c---------------------------------
 
      crec = 'CCM ' // typehis(1:4)
      nw = 16 + nlon + norec + nlev + nlevp
      write (ndata, err=900, iostat=ios)
     *  crec, nw, nulhdr, nuldat, ndens, nurecl,
     *  vershis, versgen, mbdate, dtime,
     *  nlon, norec, nlev, ntrk, ntrm, ntrn,
     *  along, alatg, sig, sigkmh
 
c----------------------------------
c        Write "FILE" header record
c----------------------------------
 
      crec = 'FILE'
      nw = 6
     *   + (  len(runtitlx) + len(mspathin) + len(msnamin) + len(mspath)
     *      + len(msnamhis(ifile)) + len(msnam_a(ifile))
     *      + len(cdum16) + len(msnamres(ifile))
     *     ) / 8
      write (ndata, err=900, iostat=ios)
     *  crec, nw, nfile, ifile, nsrest, cdum,
     *  runtitlx, mspathin, msnamin, mspath,
     *  msnamhis(ifile), msnam_a(ifile),
     *  cdum16, msnamres(ifile)
 
c-----------------------------------
c        Write "WRITE" header record
c-----------------------------------
 
      crec = 'WRITE'
      nw = 10
      write (ndata, err=900, iostat=ios)
     *  crec, nw, nwrit+1, ndhis, nsteph, nstep,
     *  nint(dtime)*(nsteph-1), nint(dtime)*nstep,
     *  mcdate, mcsec
 
c------------------------------------
c        Write "FIELDS" header record
c------------------------------------
 
      crec = 'FIELDS'
      nw = 3 + nushis*8
      write (ndata, err=900, iostat=ios)
     *  crec, nw, nushis,
     *  ( (cindex(m,mushis(n)),m=1,2),
     *    (mindex(m,mushis(n)),m=1,4),
     *     mindfl(mushis(n)),
     *     mindco(mushis(n)),          n=1,nushis )
 
c---------------------------------
c        Write "DATA" data records
c---------------------------------
 
c        If ndens gt 1, record length nurecl is less than nlon*norec,
c        and the data in arr is first packed into brr.
 
      crec = 'DATA'
 
      do 100 irec = 1,nuldat
        call gethr (irec, arr, iv)

c       for last-time accumulated fields, divide by counter
c       (all counters zero-ed by hisset after return to linems)
        if (mindfl(iv).eq.1) then
           do i=1,nlon*norec
             arr(i) = arr(i) / mindco(iv)
           enddo
        endif
 
        if(ndens.le.1) then
          write (ndata, err=910, iostat=ios) crec, arr
        else
          call packa (arr, nlon*norec, brr, ndens)
          write (ndata, err=910, iostat=ios) crec, (brr(j),j=1,nurecl-1)
        endif
  100 continue
 
      return
 
  900 write(nout,9000) crec, ios
 9000 format(/' *** Error in writing agcm history-file header record.'
     *       /'     crec   = ',a
     *       /'     iostat = ',i8
     *       /'     Stop in WRITH.'/)
      call endrun (-1)
c
 910  write(nout,9100) irec, ios, nurecl, ndens
 9100 format(/' *** Error in writing agcm history-file data record.'
     *       /'     irec   = ',i8
     *       /'     iostat = ',i8
     *       /'     nurecl = ',i8
     *       /'     ndens  = ',i8
     *       /'     Stop in WRITH.'/)
      call endrun (-1)

c=========
      else 
c=========
#ifdef NETCDF

c        netCDF files

      if (nwrit.eq.0) then       ! # of previous writes on file (comlun)

c          Global attributes

        lenattr = lenchr(runtitlx)

        rcode = nf_put_att_text (ncid, ncglobal4, 'title',
     *                           lenattr, runtitlx)

        if (nlon.eq.32 .and. norec.eq.16) then
          choriz = 'T10'
        else if (nlon.eq.48 .and. norec.eq.40) then
          choriz = 'R15'
        else if (nlon.eq.96 .and. norec.eq.48) then
          choriz = 'T31'
        else if (nlon.eq.128 .and. norec.eq.64) then
          choriz = 'T42'
        else
          choriz = 'unknown'
        endif
        lenattr = lenchr(choriz)
        rcode = nf_put_att_text (ncid, ncglobal4, 'horiz_grid',
     *                           lenattr, choriz)

        rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                           n6, 'latlon')

        nlon4  = nlon
        norec4 = norec
        nlev4  = nlev

c          Longitude dimension

        rcode = nf_def_dim (ncid, 'lon', nlon4, lonid)
        rcode = nf_def_var (ncid, 'lon', nf_float4, n1, lonid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n9, 'longitude')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n12, 'degrees east')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3') 
        rcode = nf_enddef (ncid)
        do i=1,nlon
          ztmp(i) = along(i)*180./pi
        enddo
#ifdef cray
        rcode = nf_put_vara_real (ncid, varid, n1, nlon4, ztmp)
#else
        rcode = nf_put_vara_double (ncid, varid, n1, nlon4, ztmp)
#endif
        rcode = nf_redef (ncid)

c          Latitude dimension

        rcode = nf_def_dim (ncid, 'lat', norec4, latid)
        rcode = nf_def_var (ncid, 'lat', nf_float4, n1, latid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n8, 'latitude')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n13, 'degrees north')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3') 
        rcode = nf_enddef (ncid)
        do j=1,norec
          ztmp(j) = alatg(j)*180./pi
        enddo
#ifdef cray
        rcode = nf_put_vara_real (ncid, varid, n1, norec4, ztmp)
#else
        rcode = nf_put_vara_double (ncid, varid, n1, norec4, ztmp)
#endif
        rcode = nf_redef (ncid)

c          Height dimension

        rcode = nf_def_dim (ncid, 'lev', nlev4, levid)
        rcode = nf_def_var (ncid, 'lev', nf_float4, n1, levid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n6, 'height')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n5, 'sigma')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f6.3') 
        rcode = nf_enddef (ncid)
#ifdef cray
        rcode = nf_put_vara_real (ncid, varid, n1, nlev4, sig)
#else
        rcode = nf_put_vara_double (ncid, varid, n1, nlev4, sig)
#endif
        rcode = nf_redef (ncid)

c          Min/Max dimension

        rcode = nf_def_dim (ncid, 'minmax', n2, minmaxid)
        rcode = nf_def_var (ncid, 'minmax', nf_short4, n1, minmaxid,
     *                      varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n7, 'min/max')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n6, '1 or 2')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n2, 'i2') 
        rcode = nf_enddef (ncid)
        rcode = nf_put_vara_int (ncid, varid, n1, n2, minmax)
        rcode = nf_redef (ncid)

c          Time dimension (days since base date)

        rcode = nf_def_dim (ncid, 'time', ncunlim4, timid)
        rcode = nf_def_var (ncid, 'time', nf_float4, n1, timid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n4, 'time')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n4, 'days')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f9.3') 

      endif

c        Add current time (days since base date) to time dimension

      day = dtime * 0.5*((nsteph-1)+nstep) / 86400.
      rcode = nf_enddef (ncid)
      rcode = nf_inq_varid (ncid, 'time', varid)
      ival4 = nwrit + 1
#ifdef cray
      rcode = nf_put_var1_real (ncid, varid, ival4, day)
#else
      rcode = nf_put_var1_double (ncid, varid, ival4, day)
#endif
      rcode = nf_redef (ncid)
 
c        Define variable dimensions and attributes,
c        but only for first encounter of each variable name 
c        for this file (if appending, overwrite previous)
c        (numfielnet zeroed in openf)

      do n=1,nushis

        ifdone = 0
        if (numfielnet.gt.0) then
          do m=1,numfielnet
            if (cusfielnet(m).eq.cindex(1,mushis(n))) ifdone = 1
          enddo
        endif

        if (ifdone.eq.0) then

c         increment field number and name list for this file (in comlun)
          numfielnet = numfielnet + 1
          cusfielnet(numfielnet) = cindex(1,mushis(n))

          idim(1) = lonid
          idim(2) = latid
          if (mindex(2,mushis(n)).eq.1) then
            ndim    = 3
            idim(3) = timid
          else if (mindex(2,mushis(n)).eq.2) then
            ndim    = 4
            idim(3) = minmaxid
            idim(4) = timid
          else
            ndim    = 4
            idim(3) = levid
            idim(4) = timid
          endif
          rcode = nf_def_var (ncid, cindex(1,mushis(n)), nf_float4, 
     *                        ndim, idim, varid)

          ival4 = lenchr(cindex(1,mushis(n)))
          rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                             ival4, cindex(1,mushis(n)))

          ival4 = lenchr(cindex(2,mushis(n)))
          rcode = nf_put_att_text (ncid, varid, 'units',
     *                             ival4, cindex(2,mushis(n)))

          rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                             n5, 'e13.5')

        endif

      enddo

c        Leave define mode
  
      rcode = nf_enddef (ncid)

c        Assign values to variables

      do n=1,nushis

        rcode = nf_inq_varid (ncid, cindex(1,mushis(n)), varid)
        start(1) = 1
        count(1) = nlon
        start(2) = 1
        count(2) = norec


        do j=1,mindex(2,mushis(n))
          if (mindex(2,mushis(n)).eq.1) then
            start(3) = nwrit + 1
            count(3) = 1
          else
            start(3) = j
            count(3) = 1
            start(4) = nwrit + 1
            count(4) = 1
          endif
          irec = mindex(3,mushis(n)) + j - 1
          call gethr (irec, arr, iv)

c         for last-time accumulated fields, divide by counter
c         (all counters zero-ed by hisset after return to linems)
          if(mindfl(iv).eq.1) then
             do i=1,nlon*norec
               arr(i) = arr(i) / mindco(iv)
             enddo
          endif
#ifdef cray
          rcode = nf_put_vara_real (ncid, varid, start, count, arr)
#else
          rcode = nf_put_vara_double (ncid, varid, start, count, arr)
#endif
        enddo

      enddo

c        Return to define mode

      rcode = nf_redef (ncid)
      
      return

#endif
c==========
      endif
c==========
      end
C
C ***************************** PACKA **********************************
C
c------------
#ifdef nopack
c------------
 
      subroutine packa (arr, n, brr, ndens)
 
c        Replaces Cray's PACKA library routine for machines with
c        no such library (just copies betwen arrays).
 
      dimension arr(*), brr(*)
 
      do 10 i=1,n
        brr(i) = arr(i)
   10 continue
 
      return
      end
 
c-----
#endif
c-----
C
C ***************************** REGWRITE *******************************
C
      subroutine regwrite
 
c       Write instantaneous agcm fields at specified interval (namelist
c       param REGFREQ) to a stream of files, one file per month, for
c       driving regional model RegCM2 (one-way nesting).
 
#include <PARAMS.h>
#include <COMCON.h>
#include <COMDYR.h>
#include <COMFLG.h>
#include <COMHED.h>
#include <COMHIX.h>
#include <COMINI.h>
#include <COML10.h>
#include <COMLUN.h>
#include <COMMUN.h>
#include <COMTIM.h>
 
      logical if_bod, if_eod, if_bom, if_eom
      dimension ndaypm(12)
      character cmon(12)*3, cyr*20, cname*200
      save ndaypm, cmon, cname
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     *           'JUL','AUG','SEP','OCT','NOV','DEC' /
 
      dimension regprecip (nlon,norec)
      logical first
      save regprecip, first, kount
      data first /.true./
 
      if (first) then
        kount = 0
        call zero (regprecip, nlon*norec)
        first = .false.
      endif
 
c        If user has not specified permanent filename, do nothing
 
      if (reghis.eq.' ') return
 
 
c       Set beginning/end_of day/month flags
 
      if_bom = .false.
      if_eom = .false.
      if_bod = .false.
      if_eod = .false.
 
      n = nint(86400./dtime)
      if (mod(nstep,n).eq.1) then
        if_bod = .true.
      else if ( mod(nstep,n).eq.0 ) then
        if_eod = .true.
      endif
 
      if (if_bod) then
        if (jday.eq.1) if_bom = .true.
      else if (if_eod) then
        if (jday.eq.ndaypm(jmnth)) if_eom = .true.
      endif
 
c     write (*,888) nstep, if_bod, if_eod, if_bom, if_eom
c888  format ('regwrite: nstep=',i6,'  if[bod,eod,bom,eom]=',4l3)
 
 
c        If at beginning of month, open local file and zero precip
c        accumulators. First assemble permanent name by appending
c        "_year_month" to requested prefix (namelist parameter REGHIS).
 
      if (if_bom .or. nstep.eq.nrstrt) then
         write (cyr,'(i20)') jyear
         call lstrip (cyr)
         cname = reghis(1:lenchr(reghis)) // '_' //
     *          cyr(1:lenchr(cyr))       // '_' //
     *          cmon(jmnth)
         call openf (iureg, cname, regpath, 1, 0)
      endif
 
c        Accumulate agcm precip (not needed by RegCM2, but time-averaged
c        precip rate is included in these files)
 
      kount = kount + 1
      do 20 jj=1,norec
        do 22 ji=1,nlon
          regprecip(ji,jj) = regprecip(ji,jj)
     *                     + (   frainc(ji,jj) + frains(ji,jj)
     *                         + fsnowc(ji,jj) + fsnows(ji,jj) )
   22   continue
   20 continue
 
 
c        If at requested output time, write header and one set of
c        fields; normalize, write and reset precip accumulators.
 
      n = max (nint(regfreq/dtime), 1)
      if (mod(nstep,n).eq.0) then
 
c       Header info:
        write (iureg) typehis, runtitle, nlon, norec, nlev,
     *                jyear, jmnth, jday, mcsec
 
c       land-ocean oro-like mask averaged to agcm grid
c       (0=ocean/seaice, 1=land, 2=icesheet), rounded to ints.
        do 100 jj=1,norec
          do 102 ji=1,nlon
            if (yoroocn(ji,jj).ge.0.5) then
              histmpa(ji,jj) = 0.
            else 
              if (yoroland(ji,jj).ge.yoroice(ji,jj)) then
                histmpa(ji,jj) = 1.
              else 
                histmpa(ji,jj) = 2.
              endif
            endif
  102     continue
  100   continue
        write (iureg) 'ORO     ', histmpa
 
c       Topography (m):
        do 110 jj=1,norec
          do 112 ji=1,nlon
            histmpa(ji,jj) = yfi(ji,jj)/gravit
  112     continue
  110   continue
        write (iureg) 'TOPOG   ', histmpa
 
c       Other fields for RegCM2:
        write (iureg) 'PS      ', ypsm1
        write (iureg) 'T       ', ytm1
        write (iureg) 'Q       ', yqm1
        write (iureg) 'U       ', yum1
        write (iureg) 'V       ', yvm1
 
c       Normalize, write, reset precip accumulators:
        if (kount.ge.1) then
          do 120 jj=1,norec
            do 122 ji=1,nlon
              regprecip(ji,jj) = (regprecip(ji,jj)/kount)/rhoh2o   ! m/s
  122       continue
  120     continue
        endif
        write (iureg) 'PRECIP  ', regprecip
        kount = 0
        call zero (regprecip, nlon*norec)
 
      endif
 
 
c        If end of month, dispose local file to permanent location.
 
      if (if_eom .or. nstep.eq.nstop) then
        call closef (iureg)
        call disposef (cname, regpath)
 
        write (nout,200)
     *    nstep, istepy*dtime/86400., nstep*dtime/86400.,
     *    regpath(1:lenchr(regpath)), cname(1:lenchr(cname))
  200   format(
     *         /' ----------------------------------------------'
     *         /' ******* AGCM RegCM2-DRIVER FILE OUTPUT *******'
     *         /' ----------------------------------------------'
     *         /' timestep                          = ',i10
     *         /' calendar days                     = ',f10.3
     *         /' elapsed days                      = ',f10.3
     *         /' agcm regcm2-driver file name      = ',a,a
     *         /' ----------------------------------------------')
#ifdef mss
        if (regpath.eq.' ') write (nout,202)
  202   format( ' *  File NOT saved to MSS (REGPATH is blank!) *'
     *         /' ----------------------------------------------')
#endif
        write (nout,*)
 
      endif
 
      return
      end
c
c ***************************** STORMWRITE *****************************
c
      subroutine stormwrite
 
c       Write instantaneous agcm fields at specified interval (namelist
c       param STORMFREQ) to a stream of Genesis history files (namelist
c       param STORMHIS), one file per month, for storm-track diagnosis.
 
#include <PARAMS.h>
#include <COMDYR.h>
#include <COMLUN.h>
#include <COMTIM.h>
 
      logical if_bod, if_eod, if_bom, if_eom
      dimension ndaypm(12)
      character cmon(12)*3, cyr*20, cname*200
      save ndaypm, cmon, cname
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     *           'JUL','AUG','SEP','OCT','NOV','DEC' /
      save nwrit_storm
 
 
c        If user has not specified permanent filename, do nothing
 
      if (stormhis.eq.' ') return
 
c       Set beginning/end_of day/month flags
 
      if_bom = .false.
      if_eom = .false.
      if_bod = .false.
      if_eod = .false.
 
      n = nint(86400./dtime)
      if (mod(nstep,n).eq.1) then
        if_bod = .true.
      else if ( mod(nstep,n).eq.0 ) then
        if_eod = .true.
      endif
 
      if (if_bod) then
        if (jday.eq.1) if_bom = .true.
      else if (if_eod) then
        if (jday.eq.ndaypm(jmnth)) if_eom = .true.
      endif
 
c     write (*,888) nstep, if_bod, if_eod, if_bom, if_eom
c888  format ('stormwrite: nstep=',i6,'  if[bod,eod,bom,eom]=',4l3)
 
 
c        If at beginning of month, open local file. First assemble file
c        name by appending "_year_month" to requested prefix (namelist
c        param STORMHIS).
 
      if (if_bom .or. nstep.eq.nrstrt) then
         write (cyr,'(i20)') iabs(jyear)
         call lstrip (cyr)
         cname = stormhis(1:lenchr(stormhis)) // '_' //
     *           cyr(1:lenchr(cyr))       // '_' //
     *           cmon(jmnth)
         nusiz_storm = (1000 + nlon*norec*(3+3*nlev))
     *               * max (nint(31.*86400./stormfreq), 1)
         call openf (iustorm, cname, stormpath, 1, nusiz_storm)
         nwrit_storm = 0
      endif
 
c        If at requested output time, increment write-counter
c        and write one history "write"
 
      n = max (nint(stormfreq/dtime), 1)
 
      if (mod(nstep,n).eq.0) then
        nwrit_storm = nwrit_storm + 1
        call writhstorm (nwrit_storm)
      endif
 
c        If end of month, dispose local file to permanent location.
 
      if (if_eom .or. nstep.eq.nstop) then
        call closef (iustorm)
        call disposef (cname, stormpath)
 
        write (nout,200)
     *    nstep, istepy*dtime/86400., nstep*dtime/86400.,
     *    stormpath(1:lenchr(stormpath)), cname(1:lenchr(cname)),
     *    nwrit_storm
  200   format(
     *         /' ------------------------------------------------'
     *         /' ********* AGCM STORM-TRACK FILE OUTPUT *********'
     *         /' ------------------------------------------------'
     *         /' timestep                            = ',i10
     *         /' calendar days                       = ',f10.3
     *         /' elapsed days                        = ',f10.3
     *         /' storm-track file name               = ',a,a
     *         /' storm-track history writes          = ',i10
     *         /' ------------------------------------------------')
#ifdef mss
        if (stormpath.eq.' ') write (nout,202)
  202   format( ' *  File NOT saved to MSS (STORMPATH is blank!) *'
     *         /' ------------------------------------------------')
#endif
        write (nout,*)
 
      endif
 
      return
      end
c
c ***************************** WRITHSTORM *****************************
c
      subroutine writhstorm (nwrit_storm)
 
c       Writes one set of Genesis-agcm-history-file header records and
c       one set of data records, ie, one  history "write", to the
c       current storm-track history file. Called from stormwrite.
 
c       nwrit_storm = number of the current "write" on this file
 
#include <PARAMS.h>
#include <COMCON.h>
#include <COMFLG.h>
#include <COMHED.h>
#include <COMHIS.h>
#include <COMHIX.h>
#include <COMINI.h>
#include <COMLUN.h>
#include <COML10.h>
#include <COMMAP.h>
#include <COMMUN.h>
#include <COMSDT.h>
#include <COMTIM.h>
#include <COMTRU.h>
#include <CRDCTL.h>
 
      dimension brr(nlon*norec), workh(nlon,norec)
      equivalence (brr, workh)
      character crec*8, runtitlx*80
      character cdum*8, cdum16*16
      save cdum,cdum16
      data cdum,cdum16 /' ',' '/
 
      logical first
      data first /.true./
      save first
 
      parameter (nushis_storm = 5)        ! number of storm-track fields
      character*8 cindex_storm(2,nushis_storm)
      dimension mindex_storm(4,nushis_storm)
      save cindex_storm, mindex_storm,  nuldat_storm, nurecl_storm
 
      data cindex_storm /
     *  'ORO     ',  '        ',
     *  'Z500    ',  'm       ',
     *  'PS      ',  'Pa      ',
     *  'T       ',  'deg K   ',
     *  'OMEGA   ',  'Pa/s    '/
 
      data mindex_storm /
     *  0,   1, 0, 0,
     *  0,   1, 0, 0,
     *  0,   1, 0, 0,
     *  0,nlev, 0, 0,
     *  0,nlev, 0, 0/
 
 
c        Set record-position info (mindex_storm(3:4,*)...(1,*)is spare),
c        number of data records (nuldat_storm), and record length
c        (nurecl_storm). These are local saved variables.
 
      if (first) then
 
        do i=1,nushis_storm
          if (i.eq.1) then
            mindex_storm(3,i) = 1
          else
            mindex_storm(3,i) = mindex_storm(4,i-1) + 1
          endif
          mindex_storm(4,i) = mindex_storm(3,i) + mindex_storm(2,i) - 1
        enddo
 
        nuldat_storm = mindex_storm(4,nushis_storm)
 
        if (ndens_storm.eq.1) then
          nurecl_storm = 1 + nlon*norec
        else
          nurecl_storm =
     *      1 + (nlon*norec + ndens_storm - 1)/ndens_storm + 2
        endif
 
        first = .false.
      endif
 
c        If requested, put current file counter at end of runtitle
 
      runtitlx = runtitle
      i = lenchr(runtitle)
      if (runtitle(i-2:i).eq.'...') write(runtitlx(i-2:i),'(i3.3)')ifile
 
c---------------------------------
c        Write "CCM" header record
c---------------------------------
 
      crec = 'CCM ' // typehis(1:4)
      nw = 16 + nlon + norec + nlev + nlevp
      write (iustorm, err=900, iostat=ios)
     *  crec, nw, nulhdr, nuldat_storm, ndens_storm, nurecl_storm,
     *  vershis, versgen, mbdate, dtime,
     *  nlon, norec, nlev, ntrk, ntrm, ntrn,
     *  along, alatg, sig, sigkmh
 
c----------------------------------
c        Write "FILE" header record
c----------------------------------
 
      crec = 'FILE'
      nw = 6
     *   + (  len(runtitlx) + len(mspathin) + len(msnamin) + len(mspath)
     *      + len(msnamhis(ifile)) + len(msnam_a(ifile))
     *      + len(cdum16) + len(msnamres(ifile))
     *     ) / 8
      write (iustorm, err=900, iostat=ios)
     *  crec, nw, nfile, ifile, nsrest, cdum,
     *  runtitlx, mspathin, msnamin, mspath,
     *  msnamhis(ifile), msnam_a(ifile),
     *  cdum16, msnamres(ifile)
 
c-----------------------------------
c        Write "WRITE" header record
c-----------------------------------
 
      crec = 'WRITE'
      nw = 10
      write (iustorm, err=900, iostat=ios)
     *  crec, nw, nwrit_storm, ndhis, nstep, nstep,
     *  nint(dtime)*(nstep-1), nint(dtime)*nstep,
     *  mcdate, mcsec
 
c------------------------------------
c        Write "FIELDS" header record
c------------------------------------
 
      crec = 'FIELDS'
      nw = 3 + nushis_storm*8
      izkav = 0
      izkac = 1
      write (iustorm, err=900, iostat=ios)
     *  crec, nw, nushis_storm,
     *  ( (cindex_storm(m,n),m=1,2),
     *    (mindex_storm(m,n),m=1,4),
     *    izkav,
     *    izkac,
     *    n=1,nushis_storm
     *  )
 
c---------------------------------
c        Write "DATA" data records
c---------------------------------
 
      crec = 'DATA'
 
c        Loop over storm-track fields
 
c===============================
      do ifield = 1,nushis_storm
c===============================
 
c------
c       land-ocean mask: (yoro was derived from lmask in inisurf by
c       tiling from lsx grid to agcm grid, to 0=ocean/seaice, 1=land,
c       2=icesheet, rounded here to ints.)
 
        if (cindex_storm(1,ifield).eq.'ORO') then
 
          do jj=1,norec
            do ji=1,nlon
              if (yoroocn(ji,jj).ge.0.5) then
                histmp(ji,jj,1) = 0.
              else 
                if (yoroland(ji,jj).ge.yoroice(ji,jj)) then
                  histmp(ji,jj,1) = 1.
                else 
                  histmp(ji,jj,1) = 2.
                endif
              endif
            enddo
          enddo
 
c-----
c       500 mb geopotential height (as in linems):
c       Compute 3-d geopotential on sigma surfaces (see Mod. Proc.
c       User's Guide pg.I-72, and CCM1 Descrip.eq.(3.b.11). Matrix
c       G in commap was set in setcom, already multiplied by rair.
 
        else if (cindex_storm(1,ifield).eq.'Z500') then
 
          do jk=1,nlev
 
            do jj=1,norec
              do ji=1,nlon
                histmp(ji,jj,jk) = yfi(ji,jj)/gravit
              enddo
            enddo
 
            koff = nlev*(jk-1)
            do jl=1,nlev
              zhyd = g(jl+koff) / gravit
              do jj=1,norec
                do ji=1,nlon
c                 use t on sigma, q on eta coords (spatial error!)
                  ztv = ytm1(ji,jj,jl)*(1+zvir*yqm1(ji,jj,jl))
c                 or use unfiltered yq_sig (time error!)
c                 ztv = ytm1(ji,jj,jl)*(1+zvir*yq_sig(ji,jj,jl))
                  histmp(ji,jj,jk) = histmp(ji,jj,jk) + zhyd*ztv
                enddo
              enddo
            enddo
 
          enddo
 
c            Vertically interpolate to 500 mb pressure surface
 
          zgeo = 500*1.e2
 
          zfac = log(sig(nlev)) / log(sig(nlev-1)/sig(nlev))
          zlapse = 6.5e-3
c         zlapse = 10.e-3
          zpow = gravit / (rair*zlapse)
 
          do jj=1,norec
            do ji=1,nlon
 
              zsig = zgeo/ypsm1(ji,jj)
              do jk=2,nlev
                if (sig(jk).ge.zsig) then
                  jkm = jk-1
                  jkp = jk
                  zwei = (sig(jkp)-zsig) / (sig(jkp)-sig(jkm))
                  workh(ji,jj) =     zwei *histmp(ji,jj,jkm)
     *                         + (1.-zwei)*histmp(ji,jj,jkp)
                  goto 1014
                endif
              enddo
 
c             Requested pressure level is below lowest model level,
c             so interpolate between there and sea-level-pressure
c             sigma value (where geopotential=0). For slp, use same
c             extrapolation formula as used for slp history field.
              zts = ytm1(ji,jj,nlev)
     *            + zfac*(ytm1(ji,jj,nlev)-ytm1(ji,jj,nlev-1))
              zts = max (ytm1(ji,jj,nlev-1), zts)
              zsls = (1. + zlapse*yfi(ji,jj)/(zts*gravit))**zpow
 
              zwei = (zsls-zsig) / (zsls-sig(nlev))
              workh(ji,jj) = zwei * histmp(ji,jj,nlev)
 
 1014         continue
            enddo
          enddo
 
          call scopy (nlon*norec, workh, 1, histmp, 1)
 
c-----
c       Surface pressure:
        else if (cindex_storm(1,ifield).eq.'PS') then
          call scopy (nlon*norec, ypsm1, 1, histmp, 1)
 
c-----
c       Air temperature:
        else if (cindex_storm(1,ifield).eq.'T') then
          call scopy (nlon*norec*nlev, ytm1, 1, histmp, 1)
 
c-----
c       Specific humidity:
        else if (cindex_storm(1,ifield).eq.'Q') then
          call scopy (nlon*norec*nlev, yqm1, 1, histmp, 1)
 
c-----
c       Vertical pressure velocity:
        else if (cindex_storm(1,ifield).eq.'OMEGA') then
          call scopy (nlon*norec*nlev, vomega, 1, histmp, 1)
 
        endif
 
 
c          Write data records, one per level.
c          If ndens_storm gt 1, record length nurecl_storm is less than
c          nlon*norec, and the data is first packed into brr.
 
        do jk = 1,mindex_storm(2,ifield)
          if(ndens_storm.le.1) then
            write (iustorm, err=910, iostat=ios)
     *        crec, ((histmp(ji,jj,jk),ji=1,nlon), jj=1,norec)
          else
            call packa (histmp(1,1,jk), nlon*norec, brr, ndens_storm)
            write (iustorm, err=910, iostat=ios)
     *        crec, (brr(j),j=1,nurecl_storm-1)
          endif
        enddo
 
c==========
      enddo
c==========
 
      return
 
c-------------
c Fatal errors
c-------------
 
  900 write(nout,9000) crec, ios
 9000 format(
     *  /' *** Error in writing storm-track history-file header record.'
     *  /'     crec   = ',a
     *  /'     iostat = ',i8
     *  /'     Stop in WRITHSTORM.'/)
      call endrun (-1)
 
 910  write(nout,9100) ifield, jk, ios, nurecl_storm, ndens_storm
 9100 format(
     *  /' *** Error in writing storm-track history-file data record.'
     *  /'     ifield       = ',i8
     *  /'     jk           = ',i8
     *  /'     iostat       = ',i8
     *  /'     nurecl_storm = ',i8
     *  /'     ndens_storm  = ',i8
     *  /'     Stop in WRITHSTORM.'/)
      call endrun (-1)
 
      end
