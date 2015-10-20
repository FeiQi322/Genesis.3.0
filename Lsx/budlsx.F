c -- $Id: budlsx.F,v 1.2 2001/07/03 21:29:22 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/budlsx.F,v $
c -- $Name:  $

c-------
c budlsx
c-------
c
c----------------------------------------------------------------------
c
      SUBROUTINE BUDALL (DT, TIME, NSTEP, NSTOP, LBUDINI)
c
c        Calculates overall surface budget for globe, at end of
c        lsx timestep.
c
c        All energies associated with mass fluxes of h2o are relative
c        to ice at TMELT. This is consistent with the convention in lsx
c        and ccm's CONVAD of having all phase changes occur at TMELT.
c
c        In the flux expressions below, pre-timestep temperature values
c        should be used for exact conservation (especially previous
c        soil-surface temps for runoff) but not for now to save array
c        space.
c
c        DT       = time step between calls to BUDALL
c        TIME     = absolute time (at end of time step just done)
c        NSTEP    = agcm timestep counter
c        NSTOP    = final agcm timestep # for this run
c        LBUDINI  = .T. to force resetting of budget
c
#include <compar.h>
#include <combudl.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
#include <commapl.h>
#include <commud.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <commus.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#include <com1d.h>
c
      DIMENSION
     * DAREA(NLAT), DAREAL(NLON),   DAREAO(NLON),
     * XU(NLON),    XS(NLON),       XL(NLON),    XLB(NLON),
     * XSNO(NLON),  XSNOICE(NLON),  XICE(NLON),  ZICE(NLON),
     * TX(NLON),    TY(NLON)
      LOGICAL LBUDINI, LBUDNET, LFIRST
      SAVE SENETG, SMNETG, DAREA, LFIRST
      DATA TIMENEO, LFIRST /0., .TRUE./        ! (timeneo is in combudl)
      PARAMETER (TAREA = 4.*PI)
 
c
c       Don't do anything for multi-point 2-D icesheet driver
c
#ifdef icedrive
      if (npoi.gt.1) return
#endif
c
c       Set area of grid boxes
c
      IF (LFIRST) THEN
        DO 10 JJ=1,NLAT
#if defined (gcm) && ! defined (rcm) && ! defined (hydmodel)
          DAREA(JJ) = COS(ALAT(JJ)) * (PI/NLAT) * (2.*PI/NLON)
#else
          DAREA(JJ) = TAREA/(NLON*NLAT)
#endif
   10   CONTINUE
        LFIRST = .FALSE.
      ENDIF
c
c===========================
c     Flux section
      IF (.NOT.LBUDINI) THEN
c===========================
c
c          Agcm solar fluxes from atmosphere to surface
c
        DO 110 IB=1,NBAND
          DO 112 JJ=1,NLAT
            DO 114 JI=1,NLON
              FRNETG = FRNETG
     *               + (   (1.-ASURDXY(JI,JJ,IB))*SOLADXY(JI,JJ,IB)
     *                   + (1.-ASURIXY(JI,JJ,IB))*SOLAIXY(JI,JJ,IB) )
     *                 * DT*DAREA(JJ)
  114       CONTINUE
  112     CONTINUE
  110   CONTINUE
c
        DO 120 JJ=1,NLAT
          DO 122 JI=1,NLON
            DAREAL(JI) = DAREA(JJ)*(1.-FWATERXY(JI,JJ))
            DAREAO(JI) = DAREA(JJ)*    FWATERXY(JI,JJ)
 
c              Other agcm fluxes from atmosphere to surface
 
            FRNETG = FRNETG
     *             + (FIRAXY(JI,JJ)-AFIREXY(JI,JJ)) *DT*DAREA(JJ)
            FSNETG = FSNETG
     *             + ASENAXY(JI,JJ) *DT*DAREA(JJ)
            FLNETG = FLNETG
     *             + (   AVAPAXY(JI,JJ)
     *                   *(HVAP+HFUS+(TAXY(JI,JJ)-TMELT)*CVAP)
     *                 + RAINAXY(JI,JJ)
     *                   *((MAX(TAXY(JI,JJ),TMELT)-TMELT)*CH2O + HFUS)
     *                 + SNOWAXY(JI,JJ)
     *                   *((MIN(TAXY(JI,JJ),TMELT)-TMELT)*CICE) )
     *               * DT*DAREA(JJ)
c
            FMNETG = FMNETG
     *             + (AVAPAXY(JI,JJ) + RAINAXY(JI,JJ) + SNOWAXY(JI,JJ))
     *               * DT*DAREA(JJ)

#ifdef hydmodel
            fmhydp = fmhydp + (rainaxy(ji,jj)+snowaxy(ji,jj))
     *                        * dt*darea(jj)
            fmhyde = fmhyde - avapaxy(ji,jj)   * dt*darea(jj)
            fmhydr = fmhydr + grunofxy(ji,jj)  * dt*darea(jj)
            fmhydi = fmhydi + gdrainhxy(ji,jj) * dt*darea(jj)
#endif
c
c              Fluxes implied by icesheet constant thickness despite
c              ablation (GABLATXY has been averaged for land/ocean in
c              putind)
c
            GENETG = GENETG
     *             + ( GABLATXY(JI,JJ)*(TSOIXY(JI,JJ,1)-TMELT)*CICE )
     *               * DT*DAREA(JJ)
c
            GMNETG = GMNETG + GABLATXY(JI,JJ) * DT*DAREA(JJ)
c
c              Flux due to prescribed ocean dynamics
c
            OENETG = OENETG + QFLUXXY(JI,JJ) * DT*DAREAO(JI)
c
c              Implicit fluxes due to land/icesheet runoff/drainage and
c              ocean P-E. For instance, accounts for P-E from ocean yet
c              hocn not changed.
c
c              Implicit fluxes due to land/icesheet runoff/drainage
c              (GRUNOFXY,GDRAINXY,GADJUSTXY have been averaged for
c              land/ocean in putind)
c
            GENETG = GENETG
     *             + ( -GRUNOFXY(JI,JJ)
     *                  *((MAX(TSOIXY(JI,JJ,1),TMELT)-TMELT)*CH2O+HFUS)
     *                 -GDRAINXY(JI,JJ)
     *                  *((TSOIXY(JI,JJ,NSOILAY)-TMELT)*CH2O + HFUS) )
     *               * DT*DAREA(JJ)
c
            GMNETG = GMNETG
     *             + ( - GRUNOFXY(JI,JJ) - GDRAINXY(JI,JJ) )
     *               * DT*DAREA(JJ)
 
            GMNETGA= GMNETGA + GADJUSTXY(JI,JJ) * DT*DAREA(JJ)
c
c              Implicit fluxes due to oceanic P-E
c
            OENETG = OENETG
     *             + ( -FWOCNXY(JI,JJ)
     *                  *((TOCNXY(JI,JJ)-TMELT)*CH2O + HFUS) )
     *               * DT*DAREAO(JI)
c
            OMNETG = OMNETG - FWOCNXY(JI,JJ) * DT*DAREAO(JI)
 
  122     CONTINUE
  120   CONTINUE
c
c========================
c     End of Flux section
      ENDIF
c========================
 
c        Calculate storage if first call this run (LBUDINI=.T.) or
c        if budget-time or if last call this run (LBUDNET=.T.)
 
      LBUDNET = .FALSE.
      N =  MAX ( NINT(xdtbud/DT), 1 )
      IF (MOD(NSTEP,N).EQ.0 .OR. NSTEP.EQ.NSTOP) LBUDNET = .TRUE.
#if ! defined (gcm) && ! defined (hydmodel)
      IF (ABS(TIME-xtrun) .LT. DTTOL) LBUDNET = .TRUE.
#endif
c
c================================================
c     Storage and possible printout/reset section
      IF (LBUDINI .OR. LBUDNET) THEN
c================================================
c
        SENETG = 0.
        SMNETG = 0.
c
        DO 200 JJ=1,NLAT
          DO 202 JI=1,NLON
            DAREAL(JI) = DAREA(JJ)*(1.-FWATERXY(JI,JJ))
            DAREAO(JI) = DAREA(JJ)*    FWATERXY(JI,JJ)
c
c            Set icesheet flag (zice), snow fractional area (xi),
c            seaice fractional area (xice), and land/ocean fractional
c            area (xall).
c
            ZICE(JI) = 0.
            IF (LMASKXY(JI,JJ).EQ.2) ZICE(JI) = 1.
 
            XSNO(JI)    = FIXY(JI,JJ)               *DAREAL(JI)
            XSNOICE(JI) = FI_ICEXY(JI,JJ)*FICEXY(JI,JJ)*DAREAO(JI)
            XICE(JI)    =                 FICEXY(JI,JJ)*DAREAO(JI)
c
c              Vegetation storage
c
            XU(JI)  =  LAIXY(JI,JJ,2) * FUXY(JI,JJ) * DAREAL(JI)
            XS(JI)  =  SAIXY(JI,JJ,2) * FUXY(JI,JJ) * DAREAL(JI)
            XL(JI)  = (LAIXY(JI,JJ,1)+SAIXY(JI,JJ,1))
     *                *FLXY(JI,JJ)*(1.-FIXY(JI,JJ)) * DAREAL(JI)
            XLB(JI) = (LAIXY(JI,JJ,1)+SAIXY(JI,JJ,1))
     *                *FLXY(JI,JJ)*FIXY(JI,JJ)      * DAREAL(JI)
 
            SENETG = SENETG
     *      + XU(JI)*(   TUXY(JI,JJ)*CHU
     *                + WLIQUXY(JI,JJ)*((TUXY(JI,JJ)-TMELT)*CH2O + HFUS)
     *                + WSNOUXY(JI,JJ)* (TUXY(JI,JJ)-TMELT)*CICE )
c
     *      + XS(JI)*(   TSXY(JI,JJ)*CHS
     *                + WLIQSXY(JI,JJ)*((TSXY(JI,JJ)-TMELT)*CH2O + HFUS)
     *                + WSNOSXY(JI,JJ)* (TSXY(JI,JJ)-TMELT)*CICE )
c
     *      + XL(JI)*(   TLXY(JI,JJ)*CHL
     *                + WLIQLXY(JI,JJ)*((TLXY(JI,JJ)-TMELT)*CH2O + HFUS)
     *                + WSNOLXY(JI,JJ)* (TLXY(JI,JJ)-TMELT)*CICE )
c
     *      + XLB(JI) * TLSUBXY(JI,JJ)*CHL
c
            SMNETG = SMNETG
     *               + XU(JI) * (WLIQUXY(JI,JJ) + WSNOUXY(JI,JJ))
     *               + XS(JI) * (WLIQSXY(JI,JJ) + WSNOSXY(JI,JJ))
     *               + XL(JI) * (WLIQLXY(JI,JJ) + WSNOLXY(JI,JJ))
c
c              Snow storage
c
            DO 220 K=1,NSNOLAY
              SENETG = SENETG
     *               + XSNO(JI)   * HSNOXY(JI,JJ,K)*RHOS
     *                            * (TSNOXY(JI,JJ,K)-TMELT)*CICE
     *               + XSNOICE(JI)* HSNO_ICEXY(JI,JJ,K)*RHOS
     *                            * (TSNO_ICEXY(JI,JJ,K)-TMELT)*CICE
              SMNETG = SMNETG
     *               + XSNO(JI)   * HSNOXY(JI,JJ,K)*RHOS
     *               + XSNOICE(JI)* HSNO_ICEXY(JI,JJ,K)*RHOS
  220       CONTINUE
c
c              Soil/icesheet storage
c
            DO 230 K=1,NSOILAY
              TX(JI) = TSOIXY(JI,JJ,K) - ZICE(JI)*TMELT
              TY(JI) = TSOIXY(JI,JJ,K) - TMELT
 
              IF (ITEXXY(JI,JJ,K).GE.1) THEN
c               Soil, not peat:
                MSAND = ITEXXY(JI,JJ,K)/100
                MCLAY = ITEXXY(JI,JJ,K) - 100*MSAND
                ZSAND = .01*MSAND
                ZCLAY = .01*MCLAY
                ZSILT = .01*(100-MSAND-MCLAY)
                ISAND = MAX (1, MIN (NTEX, INT(ZSAND*(NTEX-1) + 1.5)))
                ICLAY = MAX (1, MIN (NTEX, INT(ZCLAY*(NTEX-1) + 1.5)))
                ZPOROS = HYDVAL(ISAND,ICLAY,3)
                ZCSOI  = ZSAND*800. + ZSILT*850. + ZCLAY*890.
                ZRHOSOI= 2660.
              ELSE IF (ITEXXY(JI,JJ,K).EQ.-2) THEN
c               Peat:
                ZPOROS = 0.77
                ZCSOI  = 1920.
                ZRHOSOI = 1300.
              ELSE
c               Icesheet:
                ZPOROS = 0.
                ZCSOI = 2.106E3
                ZRHOSOI = 0.9E3
              ENDIF
 
              SENETG = SENETG
     *               + HSOD(K)
     *                 * ( (1.-ZPOROS)*ZRHOSOI*ZCSOI*TX(JI)
     *                     + ZPOROS*(1.-WISOIXY(JI,JJ,K))
     *                       *WSOIXY(JI,JJ,K)*RHOW*(TY(JI)*CH2O+HFUS)
     *                     + ZPOROS*WISOIXY(JI,JJ,K)
     *                       *RHOW*TY(JI)*CICE ) * DAREAL(JI)
c
              SMNETG = SMNETG
     *               + HSOD(K)
     *                 * (   ZICE(JI)*RHOI
     *                     + ZPOROS*(1.-WISOIXY(JI,JJ,K))
     *                       *WSOIXY(JI,JJ,K)*RHOW
     *                     + ZPOROS*WISOIXY(JI,JJ,K)
     *                       *RHOW ) * DAREAL(JI)
 
  230       CONTINUE
c
c              Seaice storage
c
            DO 240 K=1,NICELAY
              SENETG = SENETG
     *               + XICE(JI)*HICEXY(JI,JJ,K)
     *                 *RHOICE*CSICE*(TICEXY(JI,JJ,K)-TMELT)
              SMNETG = SMNETG
     *               + XICE(JI)*HICEXY(JI,JJ,K)*RHOICE
  240       CONTINUE
c
c              Brine resevoir and ocean storage
c
            SENETG = SENETG
     *             + XICE(JI)*BRINEXY(JI,JJ)
     *             + HOCNXY(JI,JJ)*RHOW
     *               *(CH2O*(TOCNXY(JI,JJ)-TMELT)+HFUS) * DAREAO(JI)
            SMNETG = SMNETG
     *             + HOCNXY(JI,JJ)*RHOW * DAREAO(JI)
 
  202     CONTINUE
  200   CONTINUE
c
c       write(*,888) nstep,icode,istrip,senetg,smnetg
c 888   format(' budall: nstep,icode,istrip,senetg,smnetg=',3i5,2f20.4)
c
c
c       (i)  Printout budget unless budget-reset is forced (LBUDINI)
c       (ii) Zero flux accumulators and reset last budget time
c
c-----------------------------
c       Printout/reset
        IF (.NOT.LBUDINI) THEN
c-----------------------------
 
          Z = (TIME-TIMENEO) * TAREA
          ZD = Z/86400.
          FENETG = FRNETG + FSNETG + FLNETG
c
          DO 9000 ILOOP=1,2
#if defined (gcm)
            IF (ILOOP.EQ.1) IU = IOTERM2
            IF (ILOOP.EQ.2) GOTO 9000
            WRITE(IU,*) CHAR(12)
            WRITE(IU,910)
#elif defined (hydmodel)
            IF (ILOOP.EQ.1) IU = IOTERM
            IF (ILOOP.EQ.2) GOTO 9000
            WRITE(IU,*) CHAR(12)
            WRITE(IU,910)
#else
            IF (ILOOP.EQ.1) IU = IOTERM
            IF (ILOOP.EQ.2) IU = IUTABA
            WRITE(IU, 900)
#endif
     *        TIMENEO/86400., TIME/86400.,
 
     *        FENETG/Z, GENETG/Z, OENETG/Z,
     *        (FENETG+GENETG+OENETG)/Z, (SENETG-SENETGO)/Z,
     *        (FENETG+GENETG+OENETG-(SENETG-SENETGO))/Z,
 
     *        FMNETG/ZD, 
#ifdef hydmodel
     *        fmhydp/zd, fmhyde/zd, fmhydr/zd, fmhydi/zd, 
#endif
     *        GMNETG/ZD, OMNETG/ZD,
     *        (FMNETG+GMNETG+OMNETG)/ZD, (SMNETG-SMNETGO)/ZD,
     *        (FMNETG+GMNETG+OMNETG-(SMNETG-SMNETGO))/ZD,
     *        GMNETGA/ZD
c
  900       FORMAT(
     *     /' -----------------------------------'
     *     /' *** SURFACE-MODEL BUDGET VALUES ***'
     *     /' -----------------------------------'
     *     /' Between  days:',2F9.3
 
     *     /' Energy (W/m2):',
     *      '  agcm :',F10.4,'  runof:',F10.4,'  ocn:',F10.4
     *     /'               ',
     *      '  totfl:',F10.4,'  dstor:',F10.4,'  err:',F10.4
 
     *    //' H2O (mm/day) :',
     *      '  agcm :',F10.6,'  runof:',F10.6,'  ocn:',F10.6
     *     /'               ',
     *      '  totfl:',F10.6,'  dstor:',F10.6,'  err:',F10.6,
     *      '   wadj:',F10.6
     *     /' -----------------------------------'/)
c
  910       FORMAT(
     *     /' ----------------------------------------------', 24('-')
     *     /' *** SURFACE-MODEL BUDGET BETWEEN ELAPSED DAYS ', 2F10.3,
     *      ' ***'
     *     /' ----------------------------------------------', 24('-')
     *    //'  (1) net downward heat flux from agcm     (W/m2)  ',F24.15
     *     /'  (2) icesheet ablation sensible heat      (W/m2)  ',F24.15
     *     /'  (3) ocean q-flux conv. + implicit P-E    (W/m2)  ',F24.15
     *     /'  (4) total surface heat flux (1+2+3)      (W/m2)  ',F24.15
     *     /'  (5) change in surface energy storage     (W/m2)  ',F24.15
     *     /'  (6) error (4-5)                          (W/m2)  ',F24.15
     *     /
     *     /'  (1) precipitation minus evaporation      (mm/day)',F24.15
#ifdef hydmodel
     *     /'  ( ) precipitation                        (mm/day)',F24.15
     *     /'  ( ) evaporation                          (mm/day)',F24.15
     *     /'  ( ) surface runoff                       (mm/day)',F24.15
     *     /'  ( ) infil (wflo5)                        (mm/day)',F24.15
#endif
     *     /'  (2) icesheet ablation                    (mm/day)',F24.15
     *     /'  (3) ocean E-P                            (mm/day)',F24.15
     *     /'  (4) total surface h2o mass flux (1+2+3)  (mm/day)',F24.15
     *     /'  (5) change in surface h2o mass storage   (mm/day)',F24.15
     *     /'  (6) error (4-5)                          (mm/day)',F24.15
     *     /'  (7) soil moisture adjustments (wadjust)  (mm/day)',F24.15
     *    //' --------------------------------------', 24('-') /)
c
 9000     CONTINUE
c
c----------------------------
c       End of printout/reset
        ENDIF
c----------------------------
c
c          Reset budget (zero flux accumulators, set previous storages
c          and time to current values)
c
        FRNETG = 0.
        FSNETG = 0.
        FLNETG = 0.
        FMNETG = 0.
#ifdef hydmodel
        fmhydp = 0.
        fmhyde = 0.
        fmhydr = 0.
        fmhydi = 0.
#endif
        GENETG = 0.
        GMNETG = 0.
        GMNETGA= 0.
        OENETG = 0.
        OMNETG = 0.
        SENETGO = SENETG
        SMNETGO = SMNETG
        TIMENEO = TIME
c
c=======================================================
c     End of storage and possible printout/reset section
      ENDIF
c=======================================================
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE BUDDO (DT, TIME, NSTEP, NSTOP)
c
c        Calculates above-ground storages and writes out budget at the
c        requested interval. Called every time step to check if time.
c        The point to be output is set in FINDPOI.
c
c        DT    = time step in between calls to BUDDO
c        TIME  = absolute time (at end of time step just done)
c        NSTEP = agcm timestep counter
c        NSTOP = final agcm timestep # for this run
c
#include <compar.h>
#include <combudl.h>
#include <comfud.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
      LOGICAL LBUD, LFIRC
      DATA LFIRC / .TRUE. /
c
c        Initialize budget quantities in combudl.
c
      DATA
     *  FSOLIN, FIRIN, FSENIN, FLATIN,  FWATIN,
     *  FSOLAB, FIRAB, FSENAB, FOG12BU, FOG34BU  / 10*0. /
c
c
      CALL FINDPOI(I)
      IF (I.EQ.0) RETURN
c
c        Do budget if time (LBUD) or if first call this run (LFIRC - to
c        set initial storages), or if last call this run (LBUD)
c
      LBUD = .FALSE.
      IF (xdtbud.eq.0.) THEN
        LBUD = .TRUE.
      ELSE
        D = MAX (xdtbud, DT)
        N = MAX (1, INT(TIME/D+.001))
        IF (ABS(TIME-N*D) .LT. DTTOL) LBUD = .TRUE.
#if ! defined (gcm) && ! defined (hydmodel)
        IF (ABS(TIME-xtrun) .LT. DTTOL) LBUD = .TRUE.
#endif
        IF (NSTEP.EQ.NSTOP) LBUD = .TRUE.
      ENDIF
c
c        If lsx is called by a driver for multiple runs over the
c        same years, protect against previous TIMEBUO being used
c
      IF (.NOT.LFIRC) THEN
        IF (TIME.LE.TIMEBUO) LFIRC = .TRUE.
      ENDIF
c
      IF (LBUD .OR. LFIRC) THEN
c
c          Set fractional areas
c
        XU = LAI(I,2) * FU(I)
        XS = SAI(I,2) * FU(I)
        XL = (LAI(I,1)+SAI(I,1)) * FL(I)*(1.-FI(I))
        XLT = LAI(I,1) * FL(I)*(1.-FI(I))
c
c          Calculate current energy and H2O mass storages within lsx
c
        HEATBU = XU*( (CHU+WLIQU(I)*CH2O+WSNOU(I)*CICE) * (TU(I)-TMELT)
     *                + WLIQU(I) * HFUS )
     *         + XS*( (CHS+WLIQS(I)*CH2O+WSNOS(I)*CICE) * (TS(I)-TMELT)
     *                + WLIQS(I) * HFUS )
     *         + XL*( (CHL+WLIQL(I)*CH2O+WSNOL(I)*CICE) * (TL(I)-TMELT)
     *                + WLIQL(I) * HFUS )
     *         + (LAI(I,1)+SAI(I,1)) * FL(I)*FI(I)
     *           * CHL*(TLSUB(I)-TMELT)
c
        WATBU = XU * (WLIQU(I)+WSNOU(I))
     *        + XS * (WLIQS(I)+WSNOS(I))
     *        + XL * (WLIQL(I)+WSNOL(I))
c
c          Don't write if first call this run, since didn't set storages
c          at start of run
c
        IF (.NOT. LFIRC) THEN
          T = TIME - TIMEBUO
          D = T/.864E5
          ERR =
     *          ABS (FSOLIN-FSOLAB) + ABS (FIRIN-FIRAB)
     *        + ABS (FSENIN-FSENAB)
     *        + ABS (FSOLIN+FIRIN+FSENIN+FLATIN - (HEATBU-HEATBUO))
     *        + ABS (FWATIN - (WATBU-WATBUO))
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
          LONO = ALONBUD
          LATO = ALATBUD
#else
          LONO = NINT ( (LOCPOI(I,1)-0.5)*(360./NLON) )
          IF (LONO.GT.180) LONO = LONO - 360
          LATO = NINT ( (LOCPOI(I,2)-0.5)*(180./NLAT) - 90. )
#endif
 
          DO 9000 ILOOP=1,2
            IF (ILOOP.EQ.1) IU = IOTERM
#if defined (gcm) || defined (hydmodel)
            IF (ILOOP.EQ.2) GOTO 9000
#else
            IF (ILOOP.EQ.2) IU = IUTABA
#endif
            WRITE(IU,900)
     *        TIMEBUO/.864E5, TIME/.864E5, LONO, LATO,
     *        FSOLIN/T,  FSOLAB/T,
     *        FIRIN/T,   FIRAB/T,
     *        FSENIN/T,  FSENAB/T,
     *        (FSOLIN+FIRIN+FSENIN+FLATIN)/T,  (HEATBU-HEATBUO)/T,
     *        FWATIN/D,  (WATBU-WATBUO)/D,
     *        FOG12BU/D, FOG34BU/D,
     *        FLATIN/T,  ERR/T
9000      CONTINUE
        ENDIF
c
        FSOLIN = 0.
        FIRIN = 0.
        FSENIN = 0.
        FLATIN = 0.
        FWATIN = 0.
        FSOLAB = 0.
        FIRAB = 0.
        FSENAB = 0.
        FOG12BU = 0.
        FOG34BU = 0.
        HEATBUO = HEATBU
        WATBUO = WATBU
        TIMEBUO = TIME
c
      ENDIF
c
      LFIRC = .FALSE.
      RETURN
c
  900 FORMAT(/' Above-ground budget from',F9.3,'   to',F9.3,
     *        '    at (lon,lat) = ',2i4,':'
     *       /' solar  flux in:',F20.10,5X,'solar absorbed:',F20.10
     *       /' ir     flux in:',F20.10,5X,'ir    absorbed:',F20.10
     *       /' sensib flux in:',F20.10,5X,'sens  absorbed:',F20.10
     *       /' total  heat in:',F20.10,5X,'total d(heat) :',F20.10
     *       /' total  H2O  in:',F20.10,5X,'total d(H2O)  :',F20.10
     *       /' upper  fog flx:',F20.10,5X,'lower fog flux:',F20.10
     *       /' latent flux in:',F20.10,5X,'net   error   :',F20.10)
c
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE BUDAD1 (DT)
c
c        Accumulates above-ground budget fluxes due to lsx model
c        processes (called from SURFCTL). The quantities accumulated
c        here are reported on in BUDDO. The point to be output is set
c        in FINDPOI.
c
c        DT     = time span over which fluxes apply
c
#include <compar.h>
#include <combudl.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
c
c
      CALL FINDPOI(I)
      IF (I.EQ.0) RETURN
c
c        Set fractional areas
c
      XU = LAI(I,2) * FU(I)
      XS = SAI(I,2) * FU(I)
      XL = (LAI(I,1)+SAI(I,1)) * FL(I)*(1.-FI(I))
      XLT = LAI(I,1) * FL(I)*(1.-FI(I))
c
c        Calculate net fluxes across (into) the lsx model boundaries
c
      DO 100 IB=1,NBAND
        FSOLIN = FSOLIN
     *         + ( (1.-ASURD(I,IB))*SOLAD(I,IB)
     *           + (1.-ASURI(I,IB))*SOLAI(I,IB) ) * DT
  100 CONTINUE
c
      FSOLIN = FSOLIN
     *       + ( - (1.-FI(I))*SOLG(I) - FI(I)*SOLI(I) ) * DT
c
      FIRIN = FIRIN
     *      + ( FIRA(I) - FIRB(I) - (1.-FI(I))*FIRG(I)
     *        - FI(I)*FIRI(I) ) * DT
c
      FSENIN = FSENIN
     *       + ( FSENA(I) + (1.-FI(I))*FSENG(I) + FI(I)*FSENI(I) ) * DT
c
c        FLATIN is energy due to all forms of water. Sensible heats are
c        computed relative to ice at TMELT. See subr lsx for fudgbudr/s.
c
      FLATIN = FLATIN
     *       + ( RAINA(I) * (MAX(TA(I)-TMELT,0.)*CH2O + HFUS)
     *         + SNOWA(I) *  MIN(TA(I)-TMELT,0.)*CICE
c
     *         - (1.-FI(I))*( 
     *            (RAING(I)-fudgbudr(i))*((TRAING(I)-TMELT)*CH2O + HFUS)
     *          + (SNOWG(I)-fudgbuds(i))* (TSNOWG(I)-TMELT)*CICE )
c
     *         - FI(I)     *( RAINL(I) * ((TRAINL(I)-TMELT)*CH2O + HFUS)
     *                      + SNOWL(I) *  (TSNOWL(I)-TMELT)*CICE ) ) *DT
c
      FLATIN = FLATIN
     *     + ( (FVAPA(I) + (1.-FI(I))*(FVAPG(I)-fudgbudr(i)-fudgbuds(i))
     *                         + FI(I) *FVAPI(I))
     *          *((TA(I)-TMELT)*CVAP + HFUS + HVAP) ) * DT
c
      FLATIN = FLATIN
     *       + ( XU *FVAPUT(I)*((TU(I)-TMELT)*CH2O + HFUS)
     *         + XLT*FVAPLT(I)*((TU(I)-TMELT)*CH2O + HFUS) )  * DT
c
c        H2O mass budget
c
      FWATIN = FWATIN
     *       + ( FVAPA(I) + RAINA(I) + SNOWA(I)
     *         + (1.-FI(I)) * (FVAPG(I) - RAING(I) - SNOWG(I))
     *         +      FI(I) * (FVAPI(I) - RAINL(I) - SNOWL(I))
     *         + XU*FVAPUT(I)
     *         + XLT*FVAPLT(I) ) * DT
c
c        Calculate solar, ir, sens "convergences" within the lsx model
c
      FSOLAB = FSOLAB
     *       + ( FU(I)*(SOLU(I)+SOLS(I))
     *         + (1.-FI(I))*FL(I)*SOLL(I) ) * DT
c
      FIRAB  = FIRAB
     *       + ( FU(I)*(FIRU(I)+FIRS(I))
     *         + (1.-FI(I))*FL(I)*FIRL(I) ) * DT
c
      FSENAB = FSENAB
     *       + ( - XU*FSENU(I) -XS*FSENS(I) - XL*FSENL(I) ) * DT
c
      FSENAB = FSENAB - FOG12(I) * HFOG12(I) * DT
      FSENAB = FSENAB - FOG34(I) * HFOG34(I) * DT
c
c        Fog fluxes (mass fluxes) for diagnostics only
c
      FOG12BU = FOG12BU + FOG12(I) * DT
      FOG34BU = FOG34BU + FOG34(I) * DT
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE BUDAD2 (DT)
c
c        Accumulates above-ground budget fluxes due to soil and snow
c        quantities (called from SURFCTL every surface time step).
c        The quantities accumulated here are reported on in BUDDO.
c        The point to be output is set in FINDPOI.
c
c        DT     = time span over which fluxes apply
c
#include <compar.h>
#include <combudl.h>
#include <comgrd.h>
#include <comsoi.h>
#include <com1d.h>
c
c
      CALL FINDPOI(I)
      IF (I.EQ.0) RETURN
c
c     HBURY is really sensible ht, but would mess up the FSEN* balance
      FLATIN = FLATIN - HBURY(I) * DT
c
      FLATIN = FLATIN
     *       - ( WBURY(I) * ((TBURY(I)-TMELT)*CH2O + HFUS)
     *         + SBURY(I) *  (TBURY(I)-TMELT)*CICE ) * DT
c
      FWATIN = FWATIN - (WBURY(I)+SBURY(I)) * DT
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE PRINDIAGL (DT, TIME, ISTRIP)
c
c        Writes diagnostic printout for land/icesheet - one line per
c        call per point. The point to be output is set in FINDPOI.
c
c        ISTRIP = current strip number.
c
#include <compar.h>
#include <combudl.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <commapl.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
c
      logical ifirst
      save ifirst,alono,alato,zsol,zsolin,zira,zirb,zsen,zlat,zheat,zvap
      data ifirst /.true./
      data zsol, zsolin, zira, zirb, zsen, zlat, zheat, zvap /8*0./
 
      dimension idaz(12)
      save idaz
      data idaz /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334,365/
c
c
c        If time not a multiple of requested print interval, return
c
      dtp = max (xdtpri, dt)
      if (xdtpri.ne.0.) then
        n = time/dtp + .001
        if (abs(time - n*dtp) .ge. dttol) return
      endif
c
c        Find output point in current strip
c
      call findpoi (i)
      if (i.eq.0) return
c
c        Write long,lat coords
c
      if (ifirst) then
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
        alono = alonbud
        alato = alatbud
#else
        alono = alon(locpoi(i,1))*180./pi
        if (alono.gt.180.) alono = alono - 360.
        alato = alat(locpoi(i,2))*180./pi
#endif
        write(ioterm,700) alono, alato
#if ! defined (gcm) && ! defined (hydmodel)
        write(iutaba,700) alono, alato
#endif
        write(iutabb,700) alono, alato
        write(iutabc,700) alono, alato
        write(iutabd,700) alono, alato
        write(iutabf,700) alono, alato
      endif
c
c        Write header to terminal at start of day, to files once per run
c
      n = max (1, int((time-dtp)/.864e5 + .001))
      if (ifirst .or. abs(time-dtp - n*.864e5).lt.dttol) then
        write(ioterm,800)
#if defined (gcm) || defined (hydmodel)
        ier = flush (ioterm)
#else
        ier = flush (iutaba)
#endif
        ier = flush (iutabb)
        ier = flush (iutabc)
        ier = flush (iutabd)
        ier = flush (iutabf)
      endif
c
      if (ifirst) then
#if ! defined (gcm) && ! defined (hydmodel)
        write(iutaba,800)
#endif
        write(iutabb,801)
        write(iutabc,802)
        write(iutabd,803)
        write(iutabf,805)
      endif
c
      ifirst = .false.
c
      hsx = 0.
      do 50 k=1,nsnolay
        hsx = hsx + hsno(i,k)
   50 continue
c
      hsoix = 0.
      do 55 k=1,nsoilay
        hsoix = hsoix + hsoi(i,k)
   55 continue
c
      sin = 0.
      s = 0.
      do 60 ib=1,nband
        sin = sin + solad(i,ib) + solai(i,ib)
        s   = s   + solad(i,ib)*(1.-asurd(i,ib))
     *            + solai(i,ib)*(1.-asuri(i,ib))
   60 continue
c
      write(ioterm,900)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  ta(i)-tmelt, qa(i), ua(i),
     *  tu(i)-tmelt, ts(i)-tmelt, tl(i)-tmelt,
     *  tsno(i,1)-tmelt, ti(i)-tmelt, hsx, fi(i),
     *  tsoi(i,1)-tmelt, tg(i)-tmelt, ts2(i)-tmelt, richu(i)
#if ! defined (gcm) && ! defined (hydmodel)
      write(iutaba,900)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  ta(i)-tmelt, qa(i), ua(i),
     *  tu(i)-tmelt, ts(i)-tmelt, tl(i)-tmelt,
     *  tsno(i,1)-tmelt, ti(i)-tmelt, hsx, fi(i),
     *  tsoi(i,1)-tmelt, tg(i)-tmelt, ts2(i)-tmelt, richu(i)
#endif
c
      write(iutabb,901)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  ta(i)-tmelt,
     *  tu(i)-tmelt, wliqu(i), wsnou(i),
     *  ts(i)-tmelt, wliqs(i), wsnos(i),
     *  tl(i)-tmelt, wliql(i), wsnol(i),
     *  ts2(i)-tmelt, t12(i)-tmelt, t34(i)-tmelt,
     *  z12(i), z34(i)
c
      write(iutabc,902)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  (tsoi(i,k)-tmelt,k=1,6),
     *  (wsoi(i,k),k=1,6),
     *  (wisoi(i,k),k=1,6)
c
      zsol = s
      zsolin = sin
      zira = fira(i)
      zirb = firb(i)
      zsen = fsena(i)
      zlat = hvap*fvapa(i)
      zheat = heatg(i)
      zvap = fvapa(i)
      zalb = 0.
      if (zsolin.gt.1.e-3) zalb = 1. - zsol/zsolin
c
      write(iutabd,903)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  zalb, zsol, zira, zirb, zsen, zlat, zheat,
     *  raina(i)*.864e5, snowa(i)*.864e5,
     *  zvap*.864e5, (gtransu(i)+gtransl(i))*.864e5, grunof(i)*.864e5,
     *  gdrainh(i)*.864e5, girrig(i),
     *  fu(i),lai(i,2),sai(i,2), fl(i),lai(i,1),sai(i,1), iveg(i)
c
      iyear = time/(86400.*365) + .999999
      id = time/86400. - .000001
      idyr = mod(id,365) + 1
      do 80 im=1,12
        if (idyr.le.idaz(im)) then
          imon = im
          if (imon.eq.1) then
            iday = idyr
          else
            iday = idyr - idaz(imon-1)
          endif
          goto 82
        endif
   80 continue
   82 continue
      hour = (time-id*86400.) / 3600.
      idstep = nint (hour*3600./dt)
c
      zintw = (wliqu(i)+wsnou(i))*lai(i,2)*fu(i)
     *      + (wliqs(i)+wsnos(i))*sai(i,2)*fu(i)
     *      + (wliql(i)+wsnol(i))*(lai(i,1)+sai(i,1))
     *        *fl(i)*(1.-fi(i))
      zsoiw = 0.
      zsoim = 0.
      do 84 k=1,nsoilay
        zsoiw = zsoiw + rhow*poros(i,k)*hsoi(i,k)
     *        * ( (1.-wisoi(i,k))*wsoi(i,k) + wisoi(i,k) )
        zsoim = zsoim + rhow*poros(i,k)*hsoi(i,k)
   84 continue
      if (zsoim.gt.0.) zsoiw = zsoiw/zsoim
c
      write(iutabf,905)
c    *  mod(time/.864e5,365.),
     *  idstep, iday, imon, iyear,
     *  ta(i)-tmelt, ts2(i)-tmelt, tsoi(i,1)-tmelt,
     *  fsena(i), hvap*fvapa(i),
     *  (raina(i)+snowa(i))*.864e5, fvapa(i)*.864e5,
     *  (gtransu(i)+gtransl(i))*.864e5, ginfil(i)*.864e5,
     *  grunof(i)*.864e5, gdrain(i)*.864e5,
     *  zintw, zsoiw,
     *  (1.-wisoi(i,1))*wsoi(i,1)+wisoi(i,1),
     *  (1.-wisoi(i,6))*wsoi(i,6)+wisoi(i,6)
c
#if defined (gcm) || defined (hydmodel)
      ier = flush (ioterm)
#else
      ier = flush (iutaba)
#endif
      ier = flush (iutabb)
      ier = flush (iutabc)
      ier = flush (iutabd)
      ier = flush (iutabf)
c
      return
c
  700 format(/' lon =',f7.1,'   lat =',f6.1)
c
  800 format(/'    day  yr',
     *        '     ta     qa    ua',
     *        '     tu     ts     tl',
     *        '  tsno1     ti  hsno    fi',
     *        '  tsoi1     tg    ts2 richu')
  900 format( f7.3, i4,
     *        f7.2, f7.4, f6.1,
     *        3f7.2,
     *        2f7.2,f6.2,f6.3,
     *        3f7.2,f6.3)
c
  801 format(/'    day  yr',
     *        '     ta',
     *        '     tu wliqu wsnou',
     *        '     ts wliqs wsnos',
     *        '     tl wliql wsnol',
     *        '    ts2   t12   t34',
     *        '    z12   z34')
  901 format( f7.3, i4,
     *        f7.2,
     *        3(f7.2,2f6.2),
     *        f7.2,2f6.2,
     *        f7.2,f6.2)
c
  802 format(/'    day  yr',
     *        '     tsoi1  tsoi2  tsoi3  tsoi4  tsoi5  tsoi6',
     *        '     wsoi1  wsoi2  wsoi3  wsoi4  wsoi5  wsoi6',
     *        '     wiso1  wiso2  wiso3  wiso4  wiso5  wiso6')
  902 format( f7.3, i4,
     *        3x,6f7.2,
     *        3x,6f7.4,
     *        3x,6f7.4)
c
  803 format(/'    day  yr',
     *        '  alb    sol    ird    iru   sens    lat   heatg',
     *        '  raina  snowa  fvapa transp  runof drainh  irrig',
     *        '    fu  laiu  saiu    fl  lail  sail veg')
  903 format( f7.3, i4,
     *        2x,f3.2, 5f7.1, f8.1,
     *        6f7.2,f7.1,
     *        6f6.2, i4)
c
  805 format(/'step day mon  yr',
     *        '     ta   ts2 tsoi1',
     *        '   sens   flat',
     *        ' precip  fvapa transp  infil  runof  drain',
     *        ' inth2o   h2oa   h2o1   h2o6')
  905 format( 4i4,
     *        f7.2,2f6.2,
     *        2f7.1,
     *        6f7.2,
     *        f7.2,
     *        3f7.3)
c
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE PRINDIAGO (DT, TIME, ISTRIP)
c
c        Writes diagnostic printout for ocean/seaice - one line per
c        call per point. The point to be output is set in FINDPOI.
c
c        ISTRIP = current strip number.
c
#include <compar.h>
#include <combudl.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <commapl.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
c
      logical ifirst
      save ifirst, alono,alato,zsol,zsolin,zir,zsen,zlat,zheat,zvap
      data ifirst /.true./
      data zsol, zsolin, zir, zsen, zlat, zheat, zvap /7*0./
 
      dimension idaz(12)
      save idaz
      data idaz /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334,365/
c
c
c        If time not a multiple of requested print interval, return
c
      dtp = max (xdtpri, dt)
      if (xdtpri.ne.0.) then
        n = time/dtp + .001
        if (abs(time - n*dtp) .ge. dttol) return
      endif
c
c        Find output point in current strip
c
      call findpoi (i)
      if (i.eq.0) return
c
c        If grid point is not all water, don't write to terminal on
c        local machine (ioterm, written by prindiagl).
c        Also, note iutabb not used by prindiago.
c
#if defined (gcm) || defined (hydmodel)
      iotermo = ioterm
#else
      if (fwater(i).eq.1.) then
        iotermo = ioterm
      else
        iotermo = -1
      endif
#endif
      iutabao = iutaba
      iutabco = iutabc
      iutabdo = iutabd
      iutabfo = iutabf
c
c        Write long,lat coords
c
      if (ifirst) then
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
        alono = alonbud
        alato = alatbud
#else
        alono = alon(locpoi(i,1))*180./pi
        if (alono.gt.180.) alono = alono - 360.
        alato = alat(locpoi(i,2))*180./pi
#endif
        if (iotermo.ge.0) write(iotermo,700) alono, alato
#if ! defined (gcm) && ! defined (hydmodel)
        write(iutabao,700) alono, alato
#endif
        write(iutabco,700) alono, alato
        write(iutabdo,700) alono, alato
        write(iutabfo,700) alono, alato
      endif
c
c        Write header to terminal at start of day, to files once per run
c
      n = max (1, int((time-dtp)/.864e5 + .001))
      if (ifirst .or. abs(time-dtp - n*.864e5).lt.dttol) then
        if (iotermo.ge.0) then
           write(iotermo,800)
           ier = flush (iotermo)
        endif
#if ! defined (gcm) && ! defined (hydmodel)
        ier = flush (iutabao)
#endif
        ier = flush (iutabco)
        ier = flush (iutabdo)
        ier = flush (iutabfo)
      endif
c
      if (ifirst) then
#if ! defined (gcm) && ! defined (hydmodel)
        write(iutabao,800)
#endif
        write(iutabco,802)
        write(iutabdo,803)
        write(iutabfo,805)
      endif
c
      ifirst = .false.
c
      hsx = 0.
      do 50 k=1,nsnolay
        hsx = hsx + hsno(i,k)
   50 continue
c
      hix = 0.
      brinex = 0.
      xh = 0.
      do 55 k=1,nicelay
        hix = hix + hice(i,k)
        xh = xh + rhoi*hice(i,k) * (cice*(tmelt-tice(i,k)) + hfus)
   55 continue
      if (xh.gt.0.) brinex = brine(i) / xh
c
      sin = 0.
      s = 0.
      do 60 ib=1,nband
        sin = sin + solad(i,ib) + solai(i,ib)
        s   = s   + solad(i,ib)*(1.-asurd(i,ib))
     *            + solai(i,ib)*(1.-asuri(i,ib))
   60 continue
c
      if (iotermo.ge.0) write(iotermo,900)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  ta(i)-tmelt, qa(i), ua(i),
     *  tsno(i,1)-tmelt, hsx, fi(i),
     *  tice(i,1)-tmelt, hix, fice(i),
     *  tocn(i)-tmelt, hocn(i)
#if ! defined (gcm) && ! defined (hydmodel)
      write(iutabao,900)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  ta(i)-tmelt, qa(i), ua(i),
     *  tsno(i,1)-tmelt, hsx, fi(i),
     *  tice(i,1)-tmelt, hix, fice(i),
     *  tocn(i)-tmelt, hocn(i)
#endif
c
      write(iutabco,902)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  ta(i)-tmelt,
     *  (tice(i,k)-tmelt,k=1,3),
     *  (hice(i,k),k=1,3),
     *  hix, fice(i), brinex,
     *  tocn(i)-tmelt, qflux(i), fwater(i)
c
      zsol = s
      zsolin = sin
      zir = fira(i) - firb(i)
      zsen = fsena(i)
      zlat = hvap*fvapa(i)
      zheat = heatg(i)
      zvap = fvapa(i)
      zalb = 0.
      if (zsolin.gt.1.e-3) zalb = 1. - zsol/zsolin
c
      write(iutabdo,903)
     *  mod(time/.864e5,365.), int(1+(time-1.)/(.864e5*365)),
     *  zalb, zsol, zir, zsen, zlat, zheat,
     *  raina(i)*.864e5, snowa(i)*.864e5,
     *  zvap*.864e5, (gtransu(i)+gtransl(i))*.864e5, grunof(i)*.864e5
c
      iyear = time/(86400.*365) + .999999
      id = time/86400. - .000001
      idyr = mod(id,365) + 1
      do 80 im=1,12
        if (idyr.le.idaz(im)) then
          imon = im
          if (imon.eq.1) then
            iday = idyr
          else
            iday = idyr - idaz(imon-1)
          endif
          goto 82
        endif
   80 continue
   82 continue
      hour = (time-id*86400.) / 3600.
      idstep = nint (hour*3600./dt)
c
      write(iutabfo,905)
c    *  mod(time/.864e5,365.),
     *  idstep, iday, imon, iyear,
     *  ta(i)-tmelt, ts2(i)-tmelt, tice(i,1)-tmelt,
     *  fsena(i), hvap*fvapa(i),
     *  (raina(i)+snowa(i))*.864e5, fvapa(i)*.864e5,
     *  (gtransu(i)+gtransl(i))*.864e5, ginfil(i)*.864e5,
     *  grunof(i)*.864e5, gdrain(i)*.864e5
c
#if ! defined (gcm) && ! defined (hydmodel)
      ier = flush (iutabao)
#endif
      ier = flush (iutabco)
      ier = flush (iutabdo)
      ier = flush (iutabfo)
c
      return
c
  700 format(/' lon =',f7.1,'   lat =',f6.1)
c
  800 format(/'    day  yr',
     *        '     ta     qa    ua',
     *        '  tsno1  hsno    fi',
     *        '  tice1   hix  fice',
     *        '   tocn    hocn')
  900 format( f7.3, i4,
     *        f7.2, f7.4, f6.1,
     *        f7.2,f6.2,f6.3,
     *        f7.2,f6.2,f6.3,
     *        f7.2,f8.3)
c
  802 format(/'    day  yr',
     *        '     ta',
     *        '   tice1  tice2  tice3',
     *        '  hice1 hice2 hice3',
     *        '    hix  fice brine',
     *        '   tocn  qflux water')
  902 format( f7.3, i4,
     *        f7.2,
     *        f8.2,2f7.2,
     *        f7.3,2f6.3,
     *        f7.3, f6.3,f6.3,
     *        f7.2,f7.2,f6.2)
c
  803 format(/'    day  yr',
     *        '  alb    sol     ir   sens    lat   heatg',
     *        '  raina  snowa  fvapa transp  runof')
  903 format( f7.3, i4,
     *        2x,f3.2, 4f7.1, f8.1,
     *        5f7.2)
c
  805 format(/'step day mon  yr',
     *        '     ta   ts2 tice1',
     *        '   sens   flat',
     *        ' precip  fvapa transp  infil  runof  drain')
  905 format( 4i4,
     *        f7.2,2f6.2,
     *        2f7.1,
     *        6f7.2)
c
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE PRINDUMP (DT, TIME)
c
c        Writes diagnostic "dump" printout for one land/ice point,
c        every daypr days or every timestep if daypr=0,
c        to unit iunit or no printout if iunit=0.
c        The point to be output is set in FINDPOI.
c
#include <compar.h>
#include <combudl.h>
#include <comfud.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
      save daypr, iunit
      data daypr, iunit /  0.,  88 /
c
c
c        If time not a multiple of the prindump interval, return
c
      dtp = max (daypr*.864e5, dt)
      if (daypr.ne.0.) then
        n = time/dtp + .001
        if (abs(time - n*dtp) .ge. dttol) return
      endif
c
      call findpoi(i)
      if (i.eq.0) return
c
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
      lono = alonbud
      lato = alatbud
#else
      lono = nint ( (locpoi(i,1)-0.5)*(360./nlon) )
      if (lono.gt.180) lono = lono - 360
      lato = nint ( (locpoi(i,2)-0.5)*(180./nlat) - 90. )
#endif
c
      write(iunit,910) time/.864e5, lono, lato
  910 format(/1x,30('-')/
     *       /' Lsx dump for time=',f8.3,'   at (long,lat)= ',2i4 )
c
      write(iunit,920)
  920 format(1x,7x,'solu',7x,'sols',7x,'soll',7x,'solg',7x,'soli')
      write(iunit,924)
     *  solu(i), sols(i), soll(i), solg(i), soli(i)
  924 format(1x,5f11.1)
c
      write(iunit,930)
  930 format( 1x,7x,'firu',7x,'firs',7x,'firl',7x,'firg',7x,'firi',
     *           7x,'firb')
      write(iunit,934)
     *  firu(i), firs(i), firl(i), firg(i), firi(i), firb(i)
  934 format(1x,6f11.1)
c
      write(iunit,940)
  940 format( 1x,6x,'fsenu',6x,'fsens',6x,'fsenl',6x,'fseng',
     *           6x,'fseni',6x,'fsena')
      write(iunit,944)
     *  fsenu(i),fsens(i),fsenl(i),fseng(i),fseni(i),fsena(i)
  944 format(1x,6f11.1)
c
      write(iunit,950)
  950 format( 1x,6x,'fvapu',6x,'fvaps',6x,'fvapl',6x,'fvapg',
     *           6x,'fvapi',6x,'fvapa')
      xla = (lai(i,1)+sai(i,1)) * fl(i) * (1.-fi(i))
      xlt =  lai(i,1)           * fl(i) * (1.-fi(i))
      fvapu = fvapuw(i) + fvaput(i)
      if (xla.ne.0.) then
        fvapl = (xla*fvaplw(i) + xlt*fvaplt(i)) / xla
      else
        fvapl = 0.5 * (fvaplw(i) + fvaplt(i))
      endif
      x = dt
      write(iunit,954)
     *  x*fvapu,x*fvaps(i),x*fvapl,x*fvapg(i),x*fvapi(i),x*fvapa(i)
  954 format(1x,6f11.4)
c
      write(iunit,960)
  960 format( 1x,6x,'wliqu',6x,'wsnou',6x,'wliqs',6x,'wsnos',
     *           6x,'wliql',6x,'wsnol')
      write(iunit,964)
     *  wliqu(i),wsnou(i),wliqs(i),wsnos(i),wliql(i),wsnol(i)
  964 format(1x,6f11.5)
c
      write(iunit,970)
  970 format( 1x,6x,'   ta',6x,'   tu',6x,'   ts',6x,'   tl',
     *           6x,' tsoi',6x,' tsno')
      write(iunit,974)
     *  ta(i), tu(i), ts(i), tl(i), tsoi(i,1), tsno(i,1)
  974 format(1x,6f11.3)
c
      zav = asurd(i,1)*solad(i,1) + asuri(i,1)*solai(i,1)
      ziv = solad(i,1) + solai(i,1)
      if (ziv.gt.0.) zav = zav/ziv
      zan = asurd(i,2)*solad(i,2) + asuri(i,2)*solai(i,2)
      zin = solad(i,2) + solai(i,2)
      if (zin.gt.0.) zan = zan/zin
      write(iunit,982) zav, zan
  982 format(' vis/nir average albedo               :', 2f10.5)
      write(iunit,983) asurd(i,1), asuri(i,1), asurd(i,2), asuri(i,2)
  983 format(' visd/visi/nird/niri albedo           :', 4f10.5)
      write(iunit,984) exp(alogu(i)),exp(alogl(i))
  984 format(' upper/lower roughness length         :', 2f10.5)
      write(iunit,986) ua(i),u1(i),u12(i),u2(i),u3(i),u34(i),u4(i)
  986 format(' ua, u1, u12, u2, u3, u34, u4         :', 7f5.1)
      write(iunit,988) taux(i), tauy(i)
  988 format(' eastward/northward wind stress       :', 2f10.5)
c
      write(iunit,990)
  990 format(/1x,30('-')/)
      return
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE PRINHOCN
c
c        Reports current ocean mixed-layer depth, by finding the first
c        ocean point and reporting its value. (For now, slab thickness
c        is globally uniform.)
c
#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
#include <comocnxy.h>
 
#ifdef gcm
      DO 10 JJ=1,NLAT
        DO 12 JI=1,NLON
          IF (LMASKXY(JI,JJ).GE.3) THEN
            WRITE(IOTERM,900) HOCNXY(JI,JJ)
  900       FORMAT(/' Ocean mixed-layer depth (m) = ',F9.4)
            RETURN
          ENDIF
   12   CONTINUE
   10 CONTINUE
#endif
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE FINDPOI (IP)
c
c        Searches current 1-D strip LOCPOI(NPOI,2) for output point
c        specified by LONBUD,LATBUD. (These are in combudl, set in
c        inilsx). Returns IP=1-D index if found, or IP=0 if not.
c
c
#include <compar.h>
#include <combudl.h>
#include <comgrd.h>
#include <com1d.h>
c
c     135, -25   ! Central Australia
c      53,  19   ! Oman
c     149,  75   ! Offshore Siberian Arctic Ocean
c      13,  47   ! Austria
c     295,  -3   ! Amazon
c     101,  71   ! Northern Siberia
c     263,  43   ! Sioux City
c     113,   1   ! Borneo
c      75,  25   ! N.W. India
c      81,  19   ! E. India
c     127,  61   ! E. Siberia
c     105,  53   ! N. Mongolian border
c      89,  89   ! Central Arctic Ocean
c     265,  73   ! Resolute, Canada
c      75,  55   ! Omsk, USSR
c      91,  33   ! Himalayan peak
c     101, -81   ! Antarctic nlateau
c     275,  61   ! Hudson Bay
c
      IP = 0
      DO 120 I=1,NPOI
        IF ( LOCPOI(I,1).EQ.LONBUD .AND.
     *       LOCPOI(I,2).EQ.LATBUD ) THEN
          IP = I
          RETURN
        ENDIF
  120 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine tscreen (tscr, zscr, za,z1,z12,z2,z3,z34,z4,
     *                    dispu,displ, ta,t12,t34, 
     *                    fu,fl,fi,lai,sai,tg,ti,z0soi,z0sno,richu, np)
 
c        Interpolates diagnostic screen-height temperature tscr at
c        height zscr, for np points.
 
#include <compar.h>
#include <comgrd.h>
      dimension
     *  tscr(np), za(np), z1(np), z12(np), z2(np), z3(np), z34(np),
     *  z4(np), dispu(np), displ(np), ta(np), t12(np), t34(np),
     *  fu(np), fl(np), fi(np), lai(nvec,2), sai(nvec,2), 
     *  tg(np), ti(np), z0soi(np), richu(np)
#include <comrich.h>

      xvmin = 0.2 
 
      do 100 i=1,np

        xvu = fu(i)*           (lai(i,2)+sai(i,2))
        xvl = fl(i)*(1.-fi(i))*(lai(i,1)+sai(i,1))

c------------------------------------------------
        if (xvu.lt.xvmin .and. xvl.lt.xvmin) then
c------------------------------------------------

c         bare ground or snow:
          z0g   = fi(i)*z0sno + (1.-fi(i))*z0soi(i)
          t0g   = fi(i)*ti(i) + (1.-fi(i))*tg(i)
          wlog = log(zscr/z0g) / log(za(i)/z0g)
          wlin = zscr/za(i)
          weip = weiprand (richu(i))
          w    = weip*wlog + (1.-weip)*wlin
          tscr(i) = w*ta(i) + (1.-w)*t0g
 
c----------------------------------
        else if (xvu.ge.xvmin) then
c----------------------------------

          if (zscr.gt.z1(i)) then
c           above upper canopy:
            wlog = log((zscr  -dispu(i)) / (z1(i)-dispu(i)))
     *           / log((za(i) -dispu(i)) / (z1(i)-dispu(i)))
            wlin = (zscr-z12(i))/(za(i)-z12(i))
            weip = weiprand (richu(i))
            w    = weip*wlog + (1.-weip)*wlin
            tscr(i) = w*ta(i) + (1.-w)*t12(i)
          else if (zscr.gt.z12(i)) then
c           within top half of upper canopy:
            tscr(i) = t12(i)
          else if (zscr.gt.z34(i)) then
c           between mid-points of canopies:
            tscr(i) = ( (zscr-z34(i))*t12(i) + (z12(i)-zscr)*t34(i) )
     *                / (z12(i)-z34(i))
          else
c           within lower canopy:
            tscr(i) = t34(i)
          endif
 
c----------------------------------
        else if (xvl.ge.xvmin) then
c----------------------------------

          if (zscr.gt.z3(i)) then
c           above lower canopy:
            wlog = log((zscr  -displ(i)) / (z3(i)-displ(i)))
     *           / log((za(i) -displ(i)) / (z3(i)-displ(i)))
            wlin = (zscr-z34(i))/(za(i)-z34(i))
            weip = weiprand (richu(i))
            w    = weip*wlog + (1.-weip)*wlin
            tscr(i) = w*ta(i) + (1.-w)*t34(i)
          else 
c           within lower canopy:
            tscr(i) = t34(i)
          endif

c------------
        endif
c------------
 
  100 continue
 
      return
      end
