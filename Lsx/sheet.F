c -- $Id: sheet.F,v 1.2 2001/07/03 21:29:30 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/sheet.F,v $
c -- $Name:  $

c------
c sheet
c------
c
c----------------------------------------------------------------------
c
      SUBROUTINE SETSHEET
c
c        Sets diagnostic or prescribed icesheet quantities needed
c        by lsx for current 1-D icesheet strip, returned in commus.
c
#include <compar.h>
#include <comgrd.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsoi.h>
#include <comsat.h>
 
 
c        Set icesheet layer quantities for this strip.
c        Nb: Use strip array hsoi(i,k) in case it varies with position
c        in the future, but don't use "hsoixy" to save memory. So in
c        places (vegdat,budlsx,histlsx) use hsod instead of "hsoixy".
 
      DO 100 K=1,NSOILAY
        DO 102 I=1,NPOI
          HSOI(I,K) = HSOD(K)
          FSAND(I,K) = 0.
          FCLAY(I,K) = 0.
          FSILT(I,K) = 0.
          BEX(I,K) = 1.
          IBEX(I,K) = 1
          SUCTION(I,K) = 0.
          HYDRAUL(I,K) = 0.
          POROS(I,K) = 0.
          SWILT(I,K) = 0.1
          UPCON(I) = 0.
          CSOI(I,K) = 2.106E3
          RHOSOI(I,K) = 0.9E3
          CONSOI(I,K) = 2.2
  102   CONTINUE
  100 CONTINUE
 
c        Set latent heat for icesheet surface, used by lsx.
c        Use prognostic top-layer temperature tsoi, not the diagnostic
c        skin temperature tg, since budgets assume vapor is exchanged
c        with top-layer store at tsoi.
 
      DO 400 I=1,NPOI
        HVASUG(I) = HSUBF(TSOI(I,1),TA(I))
  400 CONTINUE
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SHEETFLX (DTIM)
c
c        Calculates surface fluxes for icesheet, and snow on icesheet.
c        The solar, ir, sensible and latent calculations below parallel
c        those for land surfaces in radlsx and turlsx, but are simpler.
c        Very similar to OCNFLX (ocean.f) for ocean/seaice.
c
#include <compar.h>
#include <comgrd.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comsoi.h>
#include <com1d.h>
      DIMENSION CDLARGE(NVEC), ALOGGH(NVEC),    QOCN(NVEC),
     *          TAURAN(NVEC),  URAN(NVEC),
     *          CDH(NVEC),     CDE(NVEC),       CDM(NVEC),
     *          TAV(NVEC),     QAV(NVEC),
     *          TGOLD(NVEC),   TIOLD(NVEC),
     *          TGPRE(NVEC),   TIPRE(NVEC),
     *          QG(NVEC),      QI(NVEC),
     *          DQG(NVEC),     DQI(NVEC),
     *          COG(NVEC),     COI(NVEC),
     *          WG(NVEC),      WI(NVEC),
     *          EMG(NVEC),     EMI(NVEC),       EMO(NVEC),
     *          ZIRG(NVEC),    ZIRI(NVEC)
      PARAMETER (NITER=3)
#include <comsat.h>
#include <comrich.h>
c
c--------------------
c        Solar fluxes
c--------------------
c
      CALL ZERO (SOLG, NVEC)
      CALL ZERO (SOLI, NVEC)
      CALL SOLSET
      CALL SOLSUR (2)
      IF (NSOL.GT.0) THEN
        DO 100 IB=1,NBAND
!DIR$ IVDEP
          DO 102 J=1,NSOL
            I = INDSOL(J)
            SOLG(I) = SOLG(I) + SOLAD(I,IB)*(1.-ALBSOD(I,IB))
     *                        + SOLAI(I,IB)*(1.-ALBSOI(I,IB))
            SOLI(I) = SOLI(I) + SOLAD(I,IB)*(1.-ALBSND(I,IB))
     *                        + SOLAI(I,IB)*(1.-ALBSNI(I,IB))
  102     CONTINUE
  100   CONTINUE
      ENDIF
c
c-----------------
c        IR fluxes
c-----------------
c
      DO 200 I=1,NPOI
        EMG(I) = EMILICE
        EMI(I) = EMISNOW
 
        ZG = EMG(I) * STEF * (TG(I)**4) + (1.-EMG(I))*FIRA(I)
        ZI = EMI(I) * STEF * (TI(I)**4) + (1.-EMI(I))*FIRA(I)
 
        FIRB(I) = (1.-FI(I))*ZG + FI(I)*ZI
        FIRG(I) = FIRA(I) - ZG
        FIRI(I) = FIRA(I) - ZI
  200 CONTINUE
 
c---------------------------------------------
c        Aerodynamic fluxes over icesheet/snow
c---------------------------------------------
 
c        Set atmos conditions at ZA and other temporary aero vars
c        (like TURINI)
 
      TFAC = 1./(SIGA**CAPPA)
      DO 250 I=1,NPOI
        PA = PSURF(I) * SIGA
        THICA(I) = PSURF(I)*(1.-SIGB) / GRAV
        RHOA(I) = PA / ( RAIR * TA(I) * (1. + (RVAP/RAIR-1.)*QA(I)) )
        CP(I) = CAIR * (1. + (CVAP/CAIR-1.)*QA(I))
        ZA(I) = (PSURF(I) - PA) / (RHOA(I) * GRAV)
#ifdef icedrive
        za(i) = 2. ! impose 2-m height for icedrive
#endif
        ALOGA(I) = LOG (ZA(I))
        UA(I) = MAX ( SQRT(UX(I)**2 + UY(I)**2), 1. )
  250 CONTINUE
 
c--------------------------
      DO 400 ITER = 1,NITER
c--------------------------
 
c        Iterate NITER times for stability of upper icesheet and snow
c        layers (as in turlsx). If first iteration, initialize surface
c        temperatures to icesheet/snow top-layer temps (also set average
c        sfc roughness, first-layer heat conductivities, and slope of
c        ir flux for ice/snow skin).
 
      IF (ITER.EQ.1) THEN
        ZLOGI = LOG (Z0SNO)
        DO 410 I=1,NPOI
          ZLOGG = LOG (Z0SOI(I))
          ALOGAV(I) = (1.-FI(I))*ZLOGG + FI(I)*ZLOGI
          COG(I) = CONSOI(I,1) / (0.5*HSOI(I,1))
          COI(I) = CONSNO(I,1) / (0.5*MAX(HSNO(I,1),HSNOTOP))
          ZIRG(I) = 4.*EMG(I)*STEF*(TG(I)**3)
          ZIRI(I) = 4.*EMI(I)*STEF*(TI(I)**3)
          TGOLD(I) = TG(I)
          TIOLD(I) = TI(I)
  410   CONTINUE
      ENDIF
c
c        Set implicit/explicit factors W* (0 to 1) for this iteration.
c        W* is 1 for fully implicit, 0 for fully explicit. (For first
c        iteration, IMPEXP2 sets W* to 1.)
c
      CALL IMPEXP2 (WG, TG, TGOLD, ITER, 1)
      CALL IMPEXP2 (WI, TI, TIOLD, ITER, 3)
 
c        Adjust T* for this iteration (see comments in turvap),
c        and save current "central" values for final flux calculations.
c        Also set specific humidities, average sfc temp and spec.hum.
c
      DO 420 I=1,NPOI
        TG(I) = WG(I)*TG(I) + (1.-WG(I))*TGOLD(I)
        TI(I) = WI(I)*TI(I) + (1.-WI(I))*TIOLD(I)
 
        TG(I) = MIN (TG(I),TMELT)
        TI(I) = MIN (TI(I),TMELT)
 
        TGPRE(I) = TG(I)
        TIPRE(I) = TI(I)
 
        ES      = ESAT (TG(I))
        QG(I) = QSAT (ES, PSURF(I))
        DQG(I)= DSAT (TG(I), PSURF(I))
        DQG(I)= MIN  (DQG(I), QG(I)/10.)
 
        ES    = ESAT (TI(I))
        QI(I) = QSAT (ES, PSURF(I))
        DQI(I)= DSAT (TI(I), PSURF(I))
        DQI(I)= MIN  (DQI(I), QI(I)/10.)
 
        TAV(I)  = (1.-FI(I))*TG(I) + FI(I)*TI(I)
        QAV(I)  = (1.-FI(I))*QG(I) + FI(I)*QI(I)
  420 CONTINUE
 
c        Calculate stratification factors for icesheet/snow sfc to za
 
      CALL FSTRAT (TAV, TA, TFAC, QAV, QA, VZERO, ZA,
     *             ALOGAV, ALOGAV, ALOGA, UA, RICHU, STRAMU, STRAHU, 1)
 
c        Implicit solution for skin temperatures of icesheet and snow
 
      DO 425 I=1,NPOI
        CDH(I) = UA(I) * (VONK/(ALOGA(I)-ALOGAV(I)))**2 * STRAHU(I)
        CDM(I) = UA(I) * (VONK/(ALOGA(I)-ALOGAV(I)))**2 * STRAMU(I)
        CDH(I) = MIN (CDMAXA, CDH(I) / (1. + CDH(I)/CDMAXB))
        CDM(I) = MIN (CDMAXA, CDM(I) / (1. + CDM(I)/CDMAXB))
 
        SG(I) = RHOA(I)*CDH(I)
 
        TG(I) = (   SOLG(I) + FIRG(I) + ZIRG(I)*TGOLD(I)
     *                + SG(I)*CP(I)*TA(I)*TFAC
     *                - SG(I)*HVASUG(I)*(QG(I)-DQG(I)*TG(I)-QA(I))
     *                + COG(I)*TSOI(I,1) )
     *            / (   SG(I)*(CP(I) + HVASUG(I)*DQG(I))
     *                + COG(I) + ZIRG(I) )
        TI(I)     = (   SOLI(I) + FIRI(I) + ZIRI(I)*TIOLD(I)
     *                + SG(I)*CP(I)*TA(I)*TFAC
     *                - SG(I)*HVASUI(I)*(QI(I)-DQI(I)*TI(I)-QA(I))
     *                + COI(I)*TSNO(I,1) )
     *            / (   SG(I)*(CP(I) + HVASUI(I)*DQI(I))
     *                + COI(I) + ZIRI(I) )
  425 CONTINUE
 
c        Sensible and latent fluxes, and wind stresses over
c        icesheet/snow, and net surface heat fluxes for snow and
c        icesheet models.Only need to compute these for final iteration.
 
      IF (ITER.EQ.NITER) THEN
 
        DO 430 I=1,NPOI
          TGAV = WG(I)*TG(I) + (1.-WG(I))*TGPRE(I)
          FSENG(I) = SG(I) * CP(I) * (TGAV - TA(I)*TFAC)
 
          TIAV = WI(I)*TI(I) + (1.-WI(I))*TIPRE(I)
          FSENI(I) = SG(I) * CP(I) * (TIAV - TA(I)*TFAC)
 
          QGAV = QG(I) + WG(I)*DQG(I)*(TG(I)-TGPRE(I))
          FVAPG(I) = SG(I) * (QGAV - QA(I))
 
          QIAV = QI(I) + WI(I)*DQI(I)*(TI(I)-TIPRE(I))
          FVAPI(I) = SG(I) * (QIAV - QA(I))
 
          TAUX(I) = RHOA(I) * CDM(I) * UX(I)
          TAUY(I) = RHOA(I) * CDM(I) * UY(I)
 
ccc       Partially time-implicit in ta,qa,ux,uy:
ccc       ZMT = THICA(I)/DTIM
ccc       ZFACH = ZMT * (1.-EXP(-RHOA(I)*CDH(I)*TFAC/ZMT))
ccc       ZFACE = ZMT * (1.-EXP(-RHOA(I)*CDH(I)/ZMT))
ccc       ZFACM = ZMT * (1.-EXP(-RHOA(I)*CDM(I)/ZMT))
ccc       FSENG(I) = ZFACH * CP(I) * (TGAV/TFAC - TA(I))
ccc       FSENI(I) = ZFACH * CP(I) * (TIAV/TFAC - TA(I))
ccc       FVAPG(I) = ZFACE * (QGAV - QA(I))
ccc       FVAPI(I) = ZFACE * (QIAV - QA(I))
ccc       TAUX(I) = ZFACM * UX(I)
ccc       TAUY(I) = ZFACM * UY(I)
 
          FIRG(I) = FIRG(I) - WG(I)*ZIRG(I)*(TG(I)-TGOLD(I))
          FIRI(I) = FIRI(I) - WI(I)*ZIRI(I)*(TI(I)-TIOLD(I))
          FIRB(I) = FIRB(I)
     *              + (1.-FI(I))*WG(I)*ZIRG(I)*(TG(I)-TGOLD(I))
     *              +     FI(I) *WI(I)*ZIRI(I)*(TI(I)-TIOLD(I))
 
          HEATG(I) = SOLG(I) + FIRG(I) - FSENG(I) - HVASUG(I)*FVAPG(I)
          HEATI(I) = SOLI(I) + FIRI(I) - FSENI(I) - HVASUI(I)*FVAPI(I)
 
c         Impose final ice/snow constraints on (diagnostic) skin temps
c         TG,TI (eg, for next iterations's ir flux). Don't do earlier
c         since must use matrix result for flux calcs above.
 
          TG(I) = MIN (TG(I),TMELT)
          TI(I) = MIN (TI(I),TMELT)
  430   CONTINUE
 
      ENDIF
 
c----------------------------------------
c     End of aero icesheet/snow iteration
c----------------------------------------
 
 400  CONTINUE
 
c        Set surface rain/snow rates,temps, for snow and icesheet models
 
      DO 500 I=1,NPOI
        RAINL(I)  = RAINA(I)
        TRAINL(I) = MAX (TA(I),TMELT)
        SNOWL(I)  = SNOWA(I)
        TSNOWL(I) = MIN (TA(I),TMELT)
        RAING(I)  = RAINL(I)
        TRAING(I) = TRAINL(I)
        SNOWG(I)  = SNOWL(I)
        TSNOWG(I) = TSNOWL(I)
  500 CONTINUE
 
c        Save overall surface aero fluxes (fsena,fvapa), and evap,
c        ocean roughness length, and 2-m sfc-air temperature, in commud
c        arrays (gsuvap,z0len,ts2) for history fields.
c        (nb: fsena is downward, fseng/i are upward, same for fvap* !)
 
      ZLOG2  = LOG(2.)
      ZLOG10 = LOG(10.)
      DO 600 I=1,NPOI
        FSENA(I) = - ((1.-FI(I))*FSENG(I) + FI(I)*FSENI(I))
        FVAPA(I) = - ((1.-FI(I))*FVAPG(I) + FI(I)*FVAPI(I))
 
        GSUVAP(I)= -FVAPA(I)   ! gsuvap upwards, fvapa downwards
        Z0LEN(I) = EXP (ALOGAV(I))
 
c       Calculation of ts2 is same as in subr tscreen in budlsx
        zweip  = weiprand (richu(i))
        zwlog  = (zlog2-alogav(i)) / (aloga(i)-alogav(i))
        zwlin  = 2./za(i)
        zwg    = zweip*zwlog + (1.-zweip)*zwlin 
        ts2(i) = zwg*ta(i) + (1.-zwg)*tav(i)
#if defined (rcm)
c       calculation of [t,q,u]s10 is same as in subr tscreen in budlsx
        zwlog10 = (zlog10-alogav(i)) / (aloga(i)-alogav(i))
        zwlin10 = 10./za(i)
        zwg10   = zweip*zwlog10 + (1.-zweip)*zwlin10 
        ts10(i) = zwg10*ta(i) + (1.-zwg10)*tav(i)
        qs10(i) = zwg10*qa(i) + (1.-zwg10)*qav(i)
        us10(i) = zwg10*ua(i)
#endif
  600 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SHEETCTL (DTSHEET, TIME, ISTRIP)
 
c        Steps icesheet model through one timestep
 
c        DTSHEET = time step for this call (secs)
c        TIME    = current elapsed time in seconds
c        ISTRIP  = current 1-D strip index (for diagnostics only)
 
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsno.h>
#include <comsoi.h>
#include <com1d.h>
 
      DIMENSION OWSOI(NVEC,NSOILAY),   OWISOI(NVEC,NSOILAY),
     *          OTSOI(NVEC,NSOILAY),
     *          FWTOP(NVEC),           FHTOP(NVEC),
     *          FSQUEEZ(NVEC),
     *          DH(NVEC),
     *          ZPEN(NVEC)
 
 
c        For icesheet, all incoming liquid flux runs off so
c        FWTOP is zero, but is first cooled to TMELT (so all runoff is
c        ge TMELT)
 
!DIR$ IVDEP
      DO 54 I=1,NPOI
        FWTOP(I) = 0.
        FSQUEEZ(I) = 0.
        GINFIL(I) = 0.
        GDRAIN(I) = 0.
        GDRAINH(I) = 0.
        GEVAP(I) = 0.
        GADJUST(I) = 0.
c       rain cooled to zero deg C, immediately runs off:      ! new2
        GABLAT(I) = FVAPG(I)                                  ! new2
        GRUNOF(I) = RAING(I)                                  ! new2
        FHTOP(I) = HEATG(I) + RAING(I)*CH2O*(TRAING(I)-TMELT) ! new2
c or:   rain frozen to ice sheet, doesn't immediately run off):
cold    GABLAT(I) = FVAPG(I) - RAING(I)
cold    FHTOP(I)  = HEATG(I) + RAING(I)*(  CH2O*(TRAING(I)-TMELT) + HFUS
cold *                                   + CICE*(TMELT-TSOI(I,1)) )
   54 CONTINUE
      DO K=1,NSOILAY+1
        DO I=1,NPOI
          WFLO(I,K) = 0.
        ENDDO
      ENDDO

c        Save old soil moisture OWSOI and temperatures OTSOI so
c        implicit SOILH2O and SOILHEAT can aposteriori deduce fluxes.
 
      DO 70 K=1,NSOILAY
!DIR$ IVDEP
        DO 72 I=1,NPOI
          OWSOI(I,K) = WSOI(I,K)
          OWISOI(I,K) = WISOI(I,K)
          OTSOI(I,K) = TSOI(I,K)
   72   CONTINUE
   70 CONTINUE
 
c        Set penetrating solar fluxes across ice boundaries.
c        Don't let any solar flux escape out of bottom,
c        ie, fpen(i,nsoilay+1) = 0.
 
      CALL ZERO (FPEN, NVEC*(NSOILAY+1))
      CALL ZERO (ZPEN, NPOI)
 
      DO 80 K=2,NSOILAY
!DIR$ IVDEP
        DO 82 I=1,NPOI
          ZPEN(I) = ZPEN(I) + HSOI(I,K-1)
c         FPEN(I,K) = PENFR*SOLG(I)*EXP(-ZPEN(I)/PENDE) !not for landice
   82   CONTINUE
   80 CONTINUE
 
c        Zero (for icesheet) soil moisture liquid flow
 
      CALL ZERO (WFLO, NVEC*(NSOILAY+1))
 
c        Step temperatures due to conductive heat transport
 
      CALL SOILHEAT (OTSOI, OWSOI, FHTOP, FPEN,
     *               HFLO, WFLO, CONSOI, DTSHEET)
 
c        Increment icesheet ablation rate GABLAT and runoff GRUNOF due
c        to melting, and set temps to le TMELT.(Don't actually change
c        icesheet layer thicknesses for now.)
 
!DIR$ IVDEP
      DO 300 K=1,NSOILAY
        DO 302 I=1,NPOI
          Z = CVMGT (1., 0., TSOI(I,K).GT.TMELT)
          DH(I) = RHOI*CICE*HSOI(I,K)*(TSOI(I,K)-TMELT)
          GABLAT(I) = GABLAT(I) + Z*(DH(I)/HFUS) / DTSHEET
          GRUNOF(I) = GRUNOF(I) + Z*(DH(I)/HFUS) / DTSHEET
          TSOI(I,K) = Z*TMELT + (1.-Z)*TSOI(I,K)
  302   CONTINUE
  300 CONTINUE
 
      RETURN
      END
