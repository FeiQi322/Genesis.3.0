c -- $Id: lsx.F,v 1.2 2001/07/03 21:29:29 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/lsx.F,v $
c -- $Name:  $

c-----
c lsx
c-----
c
c----------------------------------------------------------------------
c
      SUBROUTINE ALBLSX (ZENGCM,  AGCMD,  AGCMI)
c
c        Calculates overall surface solar albedos for agcm solar calc,
c        one 1D strip of sfc pts at a time. (See comments in SURFCTL.)
c
c        ZENGCM(JI,JJ)   = cosine of solar zenith angle (supplied)
c        AGCMD(JI,JJ,IB) = effective direct  surface albedo (returned)
c        AGCMI(JI,JJ,IB) = effective diffuse surface albedo (returned)
c
#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
      DIMENSION
     *  ZENGCM(NLONG,NLATG),
     *  AGCMD(NLONG,NLATG,NBAND),  AGCMI(NLONG,NLATG,NBAND)

c        Transfer cos(zen angle) from agcm grid to surface grid
 
      CALL GTOL (ZENGCM, COSZENXY, 0)
 
c        Calculate overall surface albedos ASURDXY,ASURIXY,
c        for land pts. Loop over 1D land strips. Subr GETIND
c        copies from global arrs to 1D land strip arrays, and PUTIND
c        copies them back. (Indices have been set up in SETIND, called
c        from inisurf).
 
#if defined ( omp )
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_l, indall_l, npall_l)
!$OMP* private (ib, istrip)
!$OMP* schedule(dynamic,2)
#endif

      DO 100 ISTRIP=1,NSTRIP_L
        CALL GETIND (ISTRIP, INDALL_L, NPALL_L, 1)
        CALL FWETCAL
        CALL SOLSET
        CALL SOLSUR (1)
        DO 102 IB=1,NBAND
          CALL SOLALB (IB)
  102   CONTINUE
        CALL PUTIND (1)
  100 CONTINUE
 
c        Calculate surface albedos for icesheet
 
#if defined ( omp )
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared (nstrip_i, indall_i, npall_i)
!$OMP* private (istrip)
!$OMP* schedule(dynamic,2)
#endif
c
      DO 200 ISTRIP=1,NSTRIP_I
        CALL GETIND (ISTRIP, INDALL_I, NPALL_I, 2)
        CALL SOLSET
        CALL SOLSUR (2)
        CALL PUTIND (2)
  200 CONTINUE
 
c        Calculate surface albedos for ocean/seaice. These are then
c        weighted with ASUR[D,I]XY for land or icesheet pts calculated
c        above (weighting done in PUTIND).
 
#if defined ( omp )
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_o, indall_o, npall_o)
!$OMP* private (istrip)
!$OMP* schedule(dynamic,2)
#endif
c
      DO 300 ISTRIP=1,NSTRIP_O
        CALL GETIND (ISTRIP, INDALL_O, NPALL_O, 3)
        CALL SOLSET
        CALL SOLSUR (3)
        CALL PUTIND (3)
  300 CONTINUE
c
c        Transfer albedos from surface grid to agcm grid
c
      DO 500 IB=1,NBAND
        CALL LTOG (ASURDXY(1,1,IB), AGCMD(1,1,IB), 0)
        CALL LTOG (ASURIXY(1,1,IB), AGCMI(1,1,IB), 0)
  500 CONTINUE

      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE LSX (DTLSX, CO2PPM)
c
c       Steps one 1D strip of lsx model through one timestep.
c
c       DTLSX  = timestep for this call
c       CO2PPM = co2 atmospheric amount, ppmv
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>

      dimension x1(nvec), x2(nvec)

c
c        Calculate areal fractions wetted by intercepted H2O
c        (needed by TWOSTR via SOLALB for effect of intercepted snow
c        on solar radiative properties)
c
      CALL FWETCAL
c
c        Do lsx solar calculations
c
c        Set up for solar calcs (zero SOL[U,S,L,G,I], and set
c        solar index pts with COSZEN gt 0)
c
      CALL SOLSET
c
c        SOLSUR sets surface albedos for soil and snow.
c        SOLALB repeats albedo calculations already done in ALBLSX,
c        and SOLARF uses the unit-incident-flux results from SOLALB
c        to obtain absorbed fluxes SOL[U,S,L,G,I]. SOLARF also
c        estimates incident PARs SUNP[U,L], SHAP[U,L] for TURVAP.
c
      CALL SOLSUR (1)
      DO 100 IB=1,NBAND
        CALL SOLALB (IB)
        CALL SOLARF (IB)
  100 CONTINUE
c
c        Calculate IR fluxes
c
      CALL IRRAD
c
c        Step intercepted h2o due to drip, precip and min/max limits,
c        adjusting intercepted precip to current veg temps, and
c        setting PFLUX implied by that adjustment for TURVAP
c
      CALL CASCADE (DTLSX)
c
c        Re-calculate wetted fractions, changed by CASCADE, needed by
c        TURVAP
c
      CALL FWETCAL
c
c        Step veg temperatures implicitly, and calculate sensible
c        heat and moisture fluxes
c
      CALL TURCTL (DTLSX, CO2PPM)
c
c        Step intercepted h2o due to evaporation, and do final
c        heat-conserving adjustment for snow/liquid above/below
c        melt point
c
      CALL CASCAD2 (DTLSX)
c
c        Set net surface heat fluxes for soil and snow models
c
      DO 500 I=1,NPOI
        HEATG(I) = SOLG(I) + FIRG(I) - FSENG(I) - HVASUG(I)*FVAPG(I)
        HEATI(I) = SOLI(I) + FIRI(I) - FSENI(I) - HVASUI(I)*FVAPI(I)
  500 CONTINUE

c       Convert any condensation (fvapg<0) on bare ground to 
c       rainfall/snowfall (according to qglif(setsoi), at temperature
c       of upper soil layer (according to hvasug(setsoi)). Otherwise,
c       have problem of where to put condensation onto frozen ground
c       with wisoi = 1 in soilctl. To keep above-ground budget working
c       in budad1, which otherwise doesn't know about this change in
c       the fluxes (with heat supplied from the soil as heatg),
c       pass adjusted flux in fudgbudr/s (commus).

      do i=1,npoi
        if (fvapg(i).lt.0.) then
          fudgbudr(i) = -    qglif(i) *fvapg(i)
          fudgbuds(i) = -(1.-qglif(i))*fvapg(i)
          fvapg(i) = 0.
        else
          fudgbudr(i) = 0.
          fudgbuds(i) = 0.
        endif
      enddo
      call mix (raing,traing, raing,traing, fudgbudr,tsoi(1,1),
     *          vzero,vzero)
      call mix (snowg,tsnowg, snowg,tsnowg, fudgbuds,tsoi(1,1), 
     *          vzero,vzero)

      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE FWETCAL
c
c        Calculates FWET[U,S,L], the fractional areas wetted by
c        intercepted H2O (liquid and snow combined). The maximum value
c        FMAX (<1) allows some transpiration even in soaked conditions.
c        Use a linear relation between FWET* and WLIQ*,WSNO* (at least
c        for small values), so that the implied "thickness" is constant
c        (equal to WLIQ*MAX, WSNO*MAX as below) and the typical amount
c        evaporated in one timestep in STEPH2O will not make WLIQ*,WSNO*
c        negative and thus cause a spurious unrecoverable H2O loss.
c        (The MAX(W*MAX,.01) below numericaly allows W*MAX = 0 without
c        blowup.) In fact evaporation in one timestep *does* sometimes
c        exceed WLIQ*MAX (currently 1 Kg/m2), so there is an additional
c        safeguard in TURVAP that limits the wetted-area aerodynamic
c        coefficients SUW,SSW,SLW. If that too fails, there is an
c        ad-hoc adjustment in STEPH2O2 to reset negative WLIQ*,WSNO*
c        amounts to zero by taking some water vapor from the atmosphere.
c
c        Also sets RLIQ[U,S,L], the proportion of FWET[U,S,L] due to
c        liquid alone. FWET,RLIQ are used in TURVAP, RLIQ in STEPH2O.
c        (So RLIQ*FWET, (1-RLIQ)*FWET are the fractional areas wetted
c        by liquid and snow individually.) If FWET is 0, choose RLIQ
c        = 1 if T[U,S,L] ge TMELT or 0 otherwize, for use by TURVAP and
c        STEPH2O in case of initial dew formation on dry surface.
c
#include <compar.h>
#include <comgrd.h>
#include <comveg.h>
#include <com1d.h>
      PARAMETER (FMAX=0.7)    ! was 0.8 in v1.02
c
c        Upper leaves
c
      DO 100 I=1,NPOI
        XLIQ = WLIQU(I) / MAX(WLIQUMAX,.01)
        XTOT = XLIQ + WSNOU(I) / MAX(WSNOUMAX,.01)
        FWETU(I) = MIN (FMAX, XTOT)
        RLIQU(I) = XLIQ / MAX (XTOT, EPSILON)
        RLIQU(I)=CVMGT(1., RLIQU(I), FWETU(I).EQ.0..AND.TU(I).GE.TMELT)
        RLIQU(I)=CVMGT(0., RLIQU(I), FWETU(I).EQ.0..AND.TU(I).LT.TMELT)
  100 CONTINUE
c
c        Upper stems
c
      DO 200 I=1,NPOI
        XLIQ = WLIQS(I) / MAX(WLIQSMAX,.01)
        XTOT = XLIQ + WSNOS(I) / MAX(WSNOSMAX,.01)
        FWETS(I) = MIN (FMAX, XTOT)
        RLIQS(I) = XLIQ / MAX (XTOT, EPSILON)
        RLIQS(I)=CVMGT(1., RLIQS(I), FWETS(I).EQ.0..AND.TS(I).GE.TMELT)
        RLIQS(I)=CVMGT(0., RLIQS(I), FWETS(I).EQ.0..AND.TS(I).LT.TMELT)
  200 CONTINUE
c
c        Lower veg
c
      DO 300 I=1,NPOI
        XLIQ = WLIQL(I) / MAX(WLIQLMAX,.01)
        XTOT = XLIQ + WSNOL(I) / MAX(WSNOLMAX,.01)
        FWETL(I) = MIN (FMAX, XTOT)
        RLIQL(I) = XLIQ / MAX (XTOT, EPSILON)
        RLIQL(I)=CVMGT(1., RLIQL(I), FWETL(I).EQ.0..AND.TL(I).GE.TMELT)
        RLIQL(I)=CVMGT(0., RLIQL(I), FWETL(I).EQ.0..AND.TL(I).LT.TMELT)
  300 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE CASCADE (DTLSX)
c
c        Steps intercepted h2o due to drip, precip, and min/max limits.
c        Calls STEPH2O for upper leaves, upper stems and lower veg in
c        turn, adjusting precips at each level.
c
c        DTLSX = time step for this call
c
c        Local variables:
c        XAI(I)   = LAI and/or SAI for veg component
c                   (allows STEPH2O to work on any veg component)
c        RAIN(I)  = rainfall at appropriate level (modified by STEPH2O)
c        TRAIN(I) = temperature of RAIN (modified by STEPH2O)
c        SNOW(I)  = snowfall at appropriate level (modified by STEPH2O)
c        TSNOW(I) = temperature of SNOW (modified by STEPH2O)
c
#include <compar.h>
#include <combudl.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION
     *  XAI(NVEC),
     *  RAIN(NVEC),  TRAIN(NVEC),   SNOW(NVEC),   TSNOW(NVEC),
     *  X1(NVEC),    X2(NVEC),      X3(NVEC),     X4(NVEC)
 
#ifndef omp
      SAVE WLIQUMIN, WLIQSMIN, WLIQLMIN, WSNOUMIN, WSNOSMIN, WSNOLMIN
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
#endif
 
c        Set WLIQ*MIN and WSNO*MIN (local variables, passed to STEPH2O).
c        Unlike WLIQ*MAX and WSNO*MAX, these are part of the lsx
c        numerical method and not from the vegetation database. The
c        value .0005 should be small compared to typical precip rates
c        times DTLSX to allow any intercepted H2O to be initiated, but
c        not too small to allow evap rates to reduce WLIQ*, WSNO* to
c        that value in a reasonable number of time steps, since
c        fractional area is ~(WLIQ*/WLIQ*MAX). See comments in FWETCAL
c        and STEPH2O.
c
c        Compute every time and don't save if multiprocessing
c        over 1-D strips
c
#ifndef omp
      IF (FIRST) THEN
#endif
        WLIQUMIN = .0005 * (DTLSX/1800.) * (WLIQUMAX/0.1)
        WLIQSMIN = .0005 * (DTLSX/1800.) * (WLIQSMAX/0.1)
        WLIQLMIN = .0005 * (DTLSX/1800.) * (WLIQLMAX/0.1)
        WSNOUMIN = .0005 * (DTLSX/1800.) * (WSNOUMAX/2.0)
        WSNOSMIN = .0005 * (DTLSX/1800.) * (WSNOSMAX/2.0)
        WSNOLMIN = .0005 * (DTLSX/1800.) * (WSNOLMAX/2.0)
#ifndef omp
        FIRST = .FALSE.
      ENDIF
#endif
c
c        Adjust rainfall and snowfall rates at above-tree level,
c        allowing for gcm air temp and upper-level fog formation.
c        Test on HFOG12 = HVAP,HSUB, last set in TURVAP,to distinguish
c        liquid or ice fog (in case T12 crossed TMELT in last
c        iteration in TURVAP). Actually test on le/gt 0.5*(HVAP+HSUB)
c        to allow for small temp-corrections in HVAPF,HSUBF stmt fns.
c
c        Really we should set fog temperature to T12 of last-but-one
c        iter in previous TURVAP, but for simplicity set to MAX or MIN
c        (T12, TMELT) to ensure no rain below or snow above TMELT.
c        (Upper fog adds to precip above tree level.) Note that FOG12,
c        HFOG12, T12 are from previous lsx timestep, so will cause
c        slight energy and h2o mass errors in buddo, since buddo uses
c        fog values at end of this timestep.
c
      HTEST = 0.5*(HVAP+HSUB)
c
      DO 50 I=1,NPOI
        RAINU(I) = RAINA(I)
        TRAINU(I) = MAX (TA(I), TMELT)
        X1(I) = CVMGT (FOG12(I), 0., HFOG12(I).LE.HTEST )
        X2(I) = MAX (T12(I), TMELT)
   50 CONTINUE
      CALL MIX (RAINU,TRAINU, RAINU,TRAINU, X1,X2, VZERO,VZERO)
c
      DO 52 I=1,NPOI
        SNOWU(I) = SNOWA(I)
        TSNOWU(I) = MIN (TA(I), TMELT)
        X1(I) = CVMGT ( FOG12(I), 0., HFOG12(I).GT.HTEST )
        X2(I) = MIN (T12(I), TMELT)
   52 CONTINUE
      CALL MIX (SNOWU,TSNOWU, SNOWU,TSNOWU, X1,X2, VZERO,VZERO)
c
c        Set up for upper leaves.
c
      DO 100 I=1,NPOI
        XAI(I) = LAI(I,2)
        RAIN(I) = RAINU(I)
        TRAIN(I) = TRAINU(I)
        SNOW(I) = SNOWU(I)
        TSNOW(I) = TSNOWU(I)
  100 CONTINUE
c
c        Step upper leaves.
c
      CALL STEPH2O
     *  (TU,  WLIQU,  WSNOU,  XAI,  PFLUXU,  RAIN, TRAIN, SNOW, TSNOW,
     *   TDRIPU, TBLOWU, WLIQUMAX,WSNOUMAX, WLIQUMIN,WSNOUMIN, DTLSX, 3)
c
c        Set up for upper stems. (The upper stems get precip as
c        modified by upper leaves.)
c
      DO 200 I=1,NPOI
        XAI(I) = SAI(I,2)
  200 CONTINUE
c
c        Step upper stems.
c
      CALL STEPH2O
     *  (TS,  WLIQS,  WSNOS,  XAI,  PFLUXS,  RAIN, TRAIN, SNOW, TSNOW,
     *   TDRIPS, TBLOWS, WLIQSMAX,WSNOSMAX, WLIQSMIN,WSNOSMIN, DTLSX, 2)
c
c        Adjust rainfall and snowfall rates at below-tree level,
c        allowing for upper-veg interception/drip/blowoff, and
c        lower-level fog formation. Same comments about HFOG34 apply
c        as for upper-level fog above. (Lower fog adds to precip above
c        the lower veg/snow level.)
c
      DO 300 I=1,NPOI
        X1(I) = FU(I)*RAIN(I)
        X2(I) = (1.-FU(I))*RAINU(I)
        X3(I) = CVMGT ( FOG34(I), 0., HFOG34(I).LE.HTEST )
        X4(I) = MAX (T34(I), TMELT)
  300 CONTINUE
      CALL MIX (RAINL,TRAINL, X1,TRAIN, X2,TRAINU, X3,X4)
c
      DO 310 I=1,NPOI
        X1(I) = FU(I)*SNOW(I)
        X2(I) = (1.-FU(I))*SNOWU(I)
        X3(I) = CVMGT ( FOG34(I), 0., HFOG34(I).GT.HTEST )
        X4(I) = MIN (T34(I), TMELT)
  310 CONTINUE
      CALL MIX (SNOWL,TSNOWL, X1,TSNOW, X2,TSNOWU, X3,X4)
c
c        Set up for lower veg.
c
      DO 400 I=1,NPOI
        XAI(I) = LAI(I,1) + SAI(I,1)
        RAIN(I) = RAINL(I)
        TRAIN(I) = TRAINL(I)
        SNOW(I) = SNOWL(I)
        TSNOW(I) = TSNOWL(I)
  400 CONTINUE
c
c        Step lower veg.
c
      CALL STEPH2O
     *  (TL,  WLIQL,  WSNOL,  XAI,  PFLUXL,  RAIN, TRAIN, SNOW, TSNOW,
     *   TDRIPL, TBLOWL, WLIQLMAX,WSNOLMAX, WLIQLMIN,WSNOLMIN, DTLSX, 1)
c
c        Adjust rainfall and  snowfall rates at soil level,
c        allowing for lower-veg interception/drip/blowoff.
c
      DO 500 I=1,NPOI
        X1(I) = FL(I)*RAIN(I)
        X2(I) = (1.-FL(I))*RAINL(I)
  500 CONTINUE
      CALL MIX (RAING,TRAING, X1,TRAIN, X2,TRAINL, VZERO,VZERO)
c
      DO 510 I=1,NPOI
        X1(I) = FL(I)*SNOW(I)
        X2(I) = (1.-FL(I))*SNOWL(I)
  510 CONTINUE
      CALL MIX (SNOWG,TSNOWG, X1,TSNOW, X2,TSNOWL, VZERO,VZERO)
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE STEPH2O
     *  (TVEG,  WLIQ,  WSNO,  XAI,  PFLUX,  RAIN, TRAIN, SNOW, TSNOW,
     *   TDRIP, TBLOW, WLIQMAX, WSNOMAX, WLIQMIN, WSNOMIN, DTLSX,IWHICH)
c
c        Steps intercepted h2o for one canopy component (upper leaves,
c        upper stems, or lower veg) through one lsx time step, adjusting
c        for h2o sensible heat and phase changes. Also modifies precip
c        due to interception and drip,blowoff.
c
c        All arguments are supplied (unchanged) unless otherwise noted.
c        TVEG(I) = temperature of veg component T[U,S,L]
c        WLIQ(I) = intercepted liquid amount WLIQ[U,S,L] (returned)
c        WSNO(I) = intercepted snow amount WSNO[U,S,L] (returned)
c        XAI(I) = LAI, SAI, LAI+SAI for upper leaves/stems,lower veg
c        PFLUX(I)= ht flux due to adjust of intercep precip (returned)
c        (RAIN, TRAIN, SNOW, TSNOW are modified to below-veg values.)
c        RAIN(I) = rainfall rate above veg component (returned)
c        TRAIN(I)= temperature of RAIN (returned)
c        SNOW(I) = snowfall rate above veg component (returned)
c        TSNOW(I)= temperature of SNOW (returned)
c        TDRIP   = e-folding time of liquid drip  TDRIP[U,S,L]
c        TBLOW   = e-folding time of snow blowoff TBLOW[U,S,L]
c        WLIQMAX = max amount of intercepted liquid WLIQ[U,S,L]MAX
c        WSNOMAX = max amount of intercepted snow   WSNO[U,S,L]MAX
c        WLIQMIN = min amount of intercepted liquid WLIQ*MIN
c        WSNOMIN = min amount of intercepted snow   WSNO*MIN
c        DTLSX   = time step for this call
c        IWHICH  = 3, 2, 1 for upper leaves, upper stems and lower veg
c                  (for diagnostics only)
c
c        Local arrays:
c        FINT(I)= precip fraction intercepted by unit leaf/stem area
c        DRIP(I) = rate of liquid drip
c        BLOW(I) = rate of snow blowoff
c
#include <compar.h>
#include <comgrd.h>
      DIMENSION
     *  TVEG(NVEC), WLIQ(NVEC),  WSNO(NVEC), XAI(NVEC),  PFLUX(NVEC),
     *  RAIN(NVEC), TRAIN(NVEC), SNOW(NVEC), TSNOW(NVEC)
c
      DIMENSION  FINT(NVEC), DRIP(NVEC), BLOW(NVEC)
c
c
c        Calculate FINT, the intercepted precip fraction per unit
c        leaf/stem area. Note 0.5 * LAI or SAI (similar to IRRAD).
c
      DO 50 I=1,NPOI
        FINT(I) =  1.-EXP(-0.5*XAI(I))
        FINT(I) = CVMGT ( FINT(I) / MAX (XAI(I), EPSILON), 0.5,
     *                    XAI(I).GE.EPSILON)
   50 CONTINUE
c
c        Step intercepted liquid and snow amounts due to drip/blow,
c        intercepted rainfall/snowfall, and min/max limits. Also
c        adjust temperature of intercepted precip to current veg temp,
c        storing the heat needed to do this in PFLUX for use in TURVAP.
c
c        Without these PFLUXes, the implicit TURVAP calcs could not
c        account for the heat flux associated with precip adjustments,
c        especially changes of phase (see below), and so could not
c        handle equilibrium situations such as intercepted snowfall
c        being continuously melted by warm atmos fluxes, with the veg
c        temp somewhat lower than the equil atmos temp to supply heat
c        that melts the incoming snow; (TURVAP would just change veg
c        temp to atmos equil, with little sensible heat storage...then
c        final phase adjustment would return veg temp to melt point).
c        The use of the current (ie, previous timestep's) veg temp
c        gives the best estimate of what this timestep's final temp
c        will be, at least for steady conditions.
c
      DO 100 I=1,NPOI
c
c          Liquid
c
        DRIP(I) = XAI(I)*WLIQ(I)/TDRIP
        WLIQ(I) = WLIQ(I) * (1.-DTLSX/TDRIP)
c
        WLIQ(I) = WLIQ(I) + DTLSX*RAIN(I)*FINT(I)
        PFLUX(I) = RAIN(I)*FINT(I) * (TVEG(I)-TRAIN(I))*CH2O
        RAIN(I) = RAIN(I)*(1.-XAI(I)*FINT(I))
c
        X = WLIQ(I)
        WLIQ(I) = MIN (WLIQ(I), WLIQMAX)
        WLIQ(I)  = CVMGT ( 0.,  WLIQ(I), WLIQ(I).LT.WLIQMIN )
        DRIP(I) = DRIP(I) + XAI(I)*(X-WLIQ(I))/DTLSX
c
c          Snow
c
        BLOW(I) = XAI(I)*WSNO(I)/TBLOW
        WSNO(I) = WSNO(I) * (1.-DTLSX/TBLOW)
c
        WSNO(I) = WSNO(I) + DTLSX*SNOW(I)*FINT(I)
        PFLUX(I) = PFLUX(I) + SNOW(I)*FINT(I) * (TVEG(I)-TSNOW(I))*CICE
        SNOW(I) = SNOW(I)*(1.-XAI(I)*FINT(I))
c
        X = WSNO(I)
        WSNO(I) = MIN (WSNO(I), WSNOMAX)
        WSNO(I)  = CVMGT (0.,  WSNO(I), WSNO(I).LT.WSNOMIN)
        BLOW(I) = BLOW(I) + XAI(I)*(X-WSNO(I))/DTLSX
c
  100 CONTINUE
c
c        Change phase of liquid/snow below/above melt point, and add
c        required heat to PFLUX (see comments above). This will only
c        affect the precip intercepted in this timestep, since original
c        WLIQ, WSNO must have been ge/le melt point (ensured in later
c        call to CASCAD2/STEPH2O2).
c
      DO 300 I=1,NPOI
c
c          Liquid below freezing
c
        DW = CVMGT (WLIQ(I), 0., TVEG(I).LT.TMELT)
        PFLUX(I) = PFLUX(I)
     *           + DW * ((CH2O-CICE)*(TMELT-TVEG(I)) - HFUS) / DTLSX
        WLIQ(I) = WLIQ(I) - DW
        WSNO(I) = WSNO(I) + DW
c
c          Snow above freezing
c
        DW = CVMGT (WSNO(I), 0., TVEG(I).GT.TMELT)
        PFLUX(I) = PFLUX(I)
     *           + DW * ((CH2O-CICE)*(TVEG(I)-TMELT) + HFUS) / DTLSX
        WSNO(I) = WSNO(I) - DW
        WLIQ(I) = WLIQ(I) + DW
c
  300 CONTINUE
c
c          Adjust rainfall, snowfall below veg for interception
c          and drip, blowoff.
c
      CALL MIX (RAIN,TRAIN, RAIN,TRAIN, DRIP,TVEG, VZERO,VZERO)
      CALL MIX (SNOW,TSNOW, SNOW,TSNOW, BLOW,TVEG, VZERO,VZERO)
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE CASCAD2 (DTLSX)
c
c        At end of timestep, removes evaporation from intercepted h2o,
c        and does final heat-conserving adjustment for any liquid/snow
c        below/above melt point. Calls STEPH2O2 for upper leaves,
c        upper stems and lower veg in turn.
c
c        DTLSX = time step for this call
c
c        Local variables:
c        (FVEG, XAI allow STEPH2O2 to work on any veg component)
c        FVEG(I)  = fractional areal coverage of veg component
c        XAI(I)   = LAI and/or SAI for veg component
c
#include <compar.h>
#include <combudl.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION FVEG(NVEC), XAI(NVEC)
c
c
c        Set up for upper leaves.
c
      DO 100 I=1,NPOI
        FVEG(I) = FU(I)
        XAI(I) = LAI(I,2)
  100 CONTINUE
c
c        Step upper leaves.
c
      CALL STEPH2O2 (TU,  WLIQU,  WSNOU,  FVEG,  XAI,  RLIQU,  FVAPUW,
     *               CHU,  DTLSX, 3)
c
c        Set up for upper stems.
c
      DO 200 I=1,NPOI
        FVEG(I) = FU(I)
        XAI(I) = SAI(I,2)
  200 CONTINUE
c
c        Step upper stems.
c
      CALL STEPH2O2 (TS,  WLIQS,  WSNOS,  FVEG,  XAI,  RLIQS,  FVAPS,
     *               CHS,  DTLSX, 2)
c
c        Set up for lower veg.
c
      DO 400 I=1,NPOI
        FVEG(I) = (1.-FI(I))*FL(I)
        XAI(I) = LAI(I,1) + SAI(I,1)
  400 CONTINUE
c
c        Step lower veg.
c
      CALL STEPH2O2 (TL,  WLIQL,  WSNOL,  FVEG,  XAI,  RLIQL,  FVAPLW,
     *               CHL,  DTLSX, 1)
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE STEPH2O2
     *  (TVEG,  WLIQ,   WSNO,  FVEG,  XAI,  RLIQ,  FVAPW,
     *   CVEG,  DTLSX,  IWHICH)
c
c        Removes evaporation from intercepted h2o, and does final
c        heat-conserving adjustment for any liquid/snow below/above
c        melt point, for one veg component.
c
c        All arguments are supplied (unchanged) unless otherwise noted.
c        TVEG(I) = temperature of veg component T[U,S,L] (returned)
c        WLIQ(I) = intercepted liquid amount WLIQ[U,S,L] (returned)
c        WSNO(I) = intercepted snow amount WSNO[U,S,L] (returned)
c        FVEG(I) = fractional areal coverage, FU or (1-FI)*FL
c        XAI(I) = LAI, SAI, LAI+SAI for upper leaves/stems,lower veg
c        RLIQ(I) = ratio of area wetted by liquid to total wetted area
c        FVAPW(I)= wetted evap H2O flx per leaf/stem area FVAP[UW,S,LW]
c        CVEG    = specific heat of veg component CH[U,S,L]
c        DTLSX   = time step for this call
c        IWHICH  = 3, 2, 1 for upper leaves, upper stems and lower veg
c                  (for diagnostics only)
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <com1d.h>
      DIMENSION
     *  TVEG(NVEC),  WLIQ(NVEC),  WSNO(NVEC),
     *  FVEG(NVEC),  XAI(NVEC),   RLIQ(NVEC),  FVAPW(NVEC)
      DIMENSION DH(NVEC), DW(NVEC)
#include <comsat.h>
c
c
c        Step intercepted h2o due to evaporation/sublimation.
c        (FVAPW already has been multiplied by FWET factor in TURVAP,
c        so it is per unit leaf/stem area.)
c
c        Due to linear FWET factors (see comments in FWETCAL) and
c        the cap on SUW,SSW,SLW in TURVAP, evaporation in one timestep
c        should hardly ever make WLIQ or WSNO negative. But if this
c        happens, compensate by increasing vapor flux from atmosphere,
c        and decreasing sensib heat flux from atmos. (The former is
c        dangerous since it could suck moisture out of a dry atmos,
c        and both are unphysical but do fix the budget.) TVEG in HVAPF
c        and HSUBF should be pre-TURVAP-timestep values, but are not.
c
      DO 100 I=1,NPOI
        WLIQ(I) = WLIQ(I) - DTLSX*RLIQ(I)*FVAPW(I)
        WSNO(I) = WSNO(I) - DTLSX*(1.-RLIQ(I))*FVAPW(I)
 
#ifndef gcm
#ifndef icedrive
#ifndef hydmodel
        if ( (wliq(i).lt.0. or. wsno(i).lt.0.)
     *       .and. fveg(i)*xai(i).gt.0. )  then
          write (ioterm,102) iwhich, i, wliq(i), wsno(i)
          write (iutaba,102) iwhich, i, wliq(i), wsno(i)
  102   format(' ***Warning: wliq<0 or wsno<0. STEPH2O2 102.',
     *         ' iwhich, i, wliq, wsno:',2i4, 2f12.6)
        endif
#endif
#endif
#endif
 
        ZM = MAX (-WLIQ(I), 0.) *FVEG(I)*XAI(I)/DTLSX
        FVAPA(I) = FVAPA(I) + ZM
        FSENA(I) = FSENA(I) - ZM*HVAPF(TVEG(I),TA(I))
        WLIQ(I) = MAX (WLIQ(I), 0.)
c
        ZM = MAX (-WSNO(I), 0.) *FVEG(I)*XAI(I)/DTLSX
        FVAPA(I) = FVAPA(I) + ZM
        FSENA(I) = FSENA(I) - ZM*HSUBF(TVEG(I),TA(I))
        WSNO(I) = MAX (WSNO(I), 0.)
  100 CONTINUE
c
c        Final heat-conserving correction for liquid/snow below/above
c        melt point
c
      DO 200 I=1,NPOI
c
        CHAV = CVEG + CH2O*WLIQ(I) + CICE*WSNO(I)
c
c          Correct for liquid below melt point.
c          (Nb: if TVEG > TMELT or WLIQ = 0, nothing changes.)
c
        DH(I) = CHAV*(TMELT - TVEG(I))
        DW(I) = MIN (WLIQ(I), MAX (0., DH(I)/HFUS))
        WLIQ(I) = WLIQ(I) - DW(I)
        WSNO(I) = WSNO(I) + DW(I)
        CHAV = CVEG + CH2O*WLIQ(I) + CICE*WSNO(I)
        TVEG(I) = TMELT - (DH(I)-HFUS*DW(I)) / CHAV
c
c          Correct for snow above melt point.
c          (Nb: if TVEG < TMELT or WSNO = 0, nothing changes.)
c
        DH(I) = CHAV*(TVEG(I) - TMELT)
        DW(I) = MIN (WSNO(I), MAX (0., DH(I)/HFUS))
        WSNO(I) = WSNO(I) - DW(I)
        WLIQ(I) = WLIQ(I) + DW(I)
        CHAV = CVEG + CH2O*WLIQ(I) + CICE*WSNO(I)
        TVEG(I) = TMELT + (DH(I)-HFUS*DW(I)) / CHAV
C
  200 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE MIX (XM,TM, X1,T1, X2,T2, X3,T3)
c
c        Calorimetrically mixes masses X1,X2,X3 with temperatures
c        T1,T2,T3 into combined mass XM with temperature TM.
c        XM,TM may be returned into same location as one of X1,T1,..,
c        so hold result temporarily in XTMP,TTMP below. Will work
c        if some of X1,X2,X3 have opposite signs, but may give
c        unphysical TM's.
c
#include <compar.h>
#include <comgrd.h>
      DIMENSION  XM(NVEC), TM(NVEC), X1(NVEC), T1(NVEC),
     *           X2(NVEC), T2(NVEC), X3(NVEC), T3(NVEC)
c
c
      DO 100 I=1,NPOI
        XTMP = X1(I) + X2(I) + X3(I)
 
c       IF (ABS(XTMP).GE.EPSILON) THEN
c         TTMP = (T1(I)*X1(I) + T2(I)*X2(I) + T3(I)*X3(I))
c    *           / XTMP
c       ELSE
c         TTMP = 0.
c         XTMP = 0.
c       ENDIF
        YTMP = SIGN ( MAX(ABS(XTMP),EPSILON), XTMP )
        TTMP = CVMGT( (T1(I)*X1(I) + T2(I)*X2(I) + T3(I)*X3(I))
     *                / YTMP, 0., ABS(XTMP).GE.EPSILON )
        XTMP = CVMGT( XTMP, 0., ABS(XTMP).GE.EPSILON)
 
        XM(I) = XTMP
        TM(I) = TTMP
  100 CONTINUE
c
      RETURN
      END
