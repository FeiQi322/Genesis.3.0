c -- $Id: turlsx.F,v 1.2 2001/07/03 21:29:31 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/turlsx.F,v $
c -- $Name:  $

c--------
c turlsx
c--------
c
c----------------------------------------------------------------------
c
      SUBROUTINE TURCTL (DTLSX, CO2PPM)
c
c        Calculates sensible heat and moisture flux coefficients,
c        and steps canopy temperatures through one timestep.
c        Atmospheric conditions at ZA are supplied in commug
c        arrays TA, QA, UX, UY, PSURF and scalars SIGA,SIGB (p/ps).
c        Downward sensible heat, moisture and wind stress fluxes at ZA
c        are returned in com1d arrays FSENA, FVAPA, TAUX, TAUY.
c        Sensible heat and moisture fluxes from solid objects to air
c        are stored (for other models and budget) in com1d arrays
c        FSEN[U,S,L,G,I], FVAP[U,S,L,G,I].
c
c        The procedure is first to compute wind speeds and aerodynamic
c        transfer coefficients in TURCOF, then call TURVAP to solve an
c        implicit linear system for temperatures and specific
c        humidities and the corresponding fluxes. This is iterated
c        NITER times for non-linearities due to stratification,
c        implicit/explicit (h2o phase), dew, fog, vpd and max soil
c        moisture uptake. T12 and Q12 are changed each iteration,
c        and TU,TS,TL,TG,TI,TAVAP,QAVAP can be adjusted too (see TURVAP)
c
c        DTLSX  = time step for this call
c        CO2PPM = co2 atmospheric amount, ppmv
c
#include <compar.h>
#include <comgrd.h>
#include <commud.h>
#include <commug.h>
#include <comsoi.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
      PARAMETER (NITER=3)
c
c        Initialize aero quantities
c
      CALL TURINI
c
c        Iterate over the whole solution NITER times
c
      DO 100 ITER = 1,NITER
c
c          Calculate wind speeds and aero transfer coeffs
c
        CALL TURCOF (ITER)
c
c         Solve implicit system
c
        CALL TURVAP (DTLSX, CO2PPM, ITER, NITER)
c
  100 CONTINUE
c
c        Calculate the wind stresses at ZA (to surface from agcm)
c
      DO 200 I=1,NPOI
c       CTAU = UA(I) * (VONK/(ALOGA(I)-ALOGU(I)))**2 * STRAMU(I)
c       CTAU = MIN (CDMAXA, CTAU / (1. + CTAU/CDMAXB))
c       TAUX(I) =  RHOA(I) * CTAU * UX(I)
c       TAUY(I) =  RHOA(I) * CTAU * UY(I)
 
c       Partially time-implicit in ux,uy:
        ZCD = RHOA(I)* UA(I) * (VONK/(ALOGA(I)-ALOGU(I)))**2 * STRAMU(I)
        ZMT = THICA(I)/DTLSX
        ZFAC = ZMT * (1.-EXP(-ZCD/ZMT))
        TAUX(I) = ZFAC * UX(I)
        TAUY(I) = ZFAC * UY(I)
  200 CONTINUE
c
c        Save roughness length in commud array for history fields.
c
      DO 300 I=1,NPOI
        Z0LEN(I) = EXP (ALOGU(I))
  300 CONTINUE
 
c        Calculate 2-m surface air temperature (diagnostic, for history)
 
      CALL TSCREEN
     *  (TS2, 2., ZA,Z1,Z12,Z2,Z3,Z34,Z4, DISPU,DISPL, TA,T12,T34,
     *   FU,FL,FI,LAI,SAI,TG,TI,Z0SOI,Z0SNO,RICHU, NPOI)
#if defined (rcm)
      call tscreen
     *  (ts10, 10., za,z1,z12,z2,z3,z34,z4, dispu,displ, ta,t12,t34,
     *   fu,fl,fi,lai,sai,tg,ti,z0soi,z0sno,richu, npoi)
      call tscreen
     *  (qs10, 10., za,z1,z12,z2,z3,z34,z4, dispu,displ, qa,q12,q34,
     *   fu,fl,fi,lai,sai,tg,ti,z0soi,z0sno,richu, npoi)
      call tscreen
     *  (us10, 10., za,z1,z12,z2,z3,z34,z4, dispu,displ, ua,u12,u34,
     *   fu,fl,fi,lai,sai,tg,ti,z0soi,z0sno,richu, npoi)
#endif
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TURINI
c
c        Initializes aerodynamic quantities that remain constant
c        through one gcm timestep. (Some quantities actually are
c        constant as long as the vegetation amounts and fractional
c        coverage remain unchanged, so could re-arrange code for
c        efficiency - currently all arrays initialized here are in
c        com1d which can be overwritten elsewhere in the gcm.)
c
c        Actually some quantities here depend on fractional snow
c        cover FI which may change within one agcm timestep (if the lsx
c        step is shorter), but this effect can be neglected for that
c        short a time.
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
#include <comsat.h>
c
c
      TFAC = 1./(SIGA**CAPPA)
      DO 100 I=1,NPOI
c
c          Atmospheric conditions at ZA. (ZA is variable, although
c          SIGA = p/ps is constant.)
c
        PA = PSURF(I) * SIGA
        RHOA(I) = PA / ( RAIR * TA(I) * (1. + (RVAP/RAIR-1.)*QA(I)) )
        THICA(I) = PSURF(I)*(1.-SIGB) / GRAV
        CP(I) = CAIR * (1. + (CVAP/CAIR-1.)*QA(I))
        ZA(I) = (PSURF(I) - PA) / (RHOA(I) * GRAV)
c       Ensure agcm level higher than canopy top:
        ZA(I) = MAX (ZA(I), ZTOP(I,2)+.50)
c
c          Aerodynamic coefficients for the lower story.
c          CVEGL (drag coeff for momentum) is proportional, and DVEGL
c          (diffusion coeff for momentum) inversely proportional,
c          to X = density of vegetation (without distinction between
c          LAI,SAI and FL*(1-FI)). X is not allowed to be exactly
c          zero to avoid divide-by-zeros, and for X>1 DVEGL is
c          proportional to 1/X**2 so that roughness length tends to
c          zero as X tends to infinity.
c
c          Also the top, bottom and displacement heights Z3(I),Z4(I),
c          DISPL(I) tend to particular values as the density tends to
c          zero, to give same results as equations for no veg at all.
c
        X = FL(I) * (1.-FI(I)) * (LAI(I,1)+SAI(I,1)) / ALAIML
        X  = MIN (X, 3.)               ! for stability with eve/leaf veg
        X1 = MIN (X, 1.)
        CVEGL = (0.4 / MAX(ZTOP(I,1)-ZBOT(I,1),.01)) * MAX(1.E-5,X)
        DVEGL = (0.1 * MAX(ZTOP(I,1)-ZBOT(I,1),.01)) / MAX(1.E-5,X,X**2)
 
        BVEGL = SQRT (2.*CVEGL/DVEGL)  ! 1 / e-folding depth in canopy
        BDL(I) = 0.5 * BVEGL * DVEGL   ! [(tau/rho)/u**2] for inf canopy
        DIL(I) = 1./DVEGL              ! 1 / diffusion coeff
 
        Z3(I) = X1 * ZTOP(I,1) + (1.-X1) * (MAX(Z0SOI(I),Z0SNO) + .01)
        Z4(I) = X1 * ZBOT(I,1) + (1.-X1) * (MAX(Z0SOI(I),Z0SNO) + .01)
        Z34(I) = 0.5*(Z3(I)+Z4(I))
        EXPHL(I) = EXP (0.5 * BVEGL * (Z3(I)-Z4(I)))
        EXPL(I) = EXPHL(I)**2
        DISPL(I) = X1 * 0.7*Z3(I)
c
c          Aerodynamic coefficients for the upper story. Same comments
c          apply as for the lower story.
c
        X = FU(I) * (LAI(I,2)+SAI(I,2)) / ALAIMU
        X  = MIN (X, 3.)               ! for stability with eve/leaf veg
        X1 = MIN (X, 1.)
        CVEGU = (0.4 / MAX(ZTOP(I,2)-ZBOT(I,2),.01)) * MAX(1.E-5,X)
        DVEGU = (0.1 * MAX(ZTOP(I,2)-ZBOT(I,2),.01)) / MAX(1.E-5,X,X**2)
 
        BVEGU = SQRT (2.*CVEGU/DVEGU)
        BDU(I) = 0.5 * BVEGU * DVEGU
        DIU(I) = 1./DVEGU
        Z1(I) = X1 * ZTOP(I,2) + (1.-X1) * (Z3(I)+.01)
        Z2(I) = X1 * ZBOT(I,2) + (1.-X1) * (Z3(I)+.01)
        Z12(I) = 0.5*(Z1(I)+Z2(I))
        EXPHU(I) = EXP (0.5 * BVEGU * (Z1(I)-Z2(I)))
        EXPU(I) = EXPHU(I)**2
        DISPU(I) = X1 * 0.7*Z1(I) + (1.-X1) * DISPL(I)
c
c          Mixing-length logarithms
c
        ALOGG(I) = LOG (Z0SOI(I))
        ALOGI(I) = LOG (Z0SNO)
        ALOGAV(I) = (1.-FI(I))*ALOGG(I) + FI(I)*ALOGI(I)
c       ALOG4 must be > Z0SOI,Z0SNO to avoid possible problems later
        ALOG4(I) = LOG ( MAX (Z4(I), 1.1*Z0SOI(I), 1.1*Z0SNO) )
        ALOG3(I) = LOG (Z3(I)-DISPL(I))
        ALOG2(I) = LOG (Z2(I)-DISPL(I))
        ALOG1(I) = LOG (Z1(I)-DISPU(I))
        ALOGA(I) = LOG (ZA(I)-DISPU(I))
c
c          Set a minimum gcm wind speed to ensure non-zero coefficients
c          and no divide-by-zero problems in TURCOF, TURSEN, TURVAP.
c
#if defined (rcm)
        UA(I) = MAX ( SQRT(UX(I)**2 + UY(I)**2), 0.1 )
#else
        UA(I) = MAX ( SQRT(UX(I)**2 + UY(I)**2), 1. )
#endif
c
c          Initialize U2, ALOGU, ALOGL for first iteration's FSTRAT
c
        U2(I) = UA(I)/EXPHU(I)
        ALOGU(I) = LOG (MAX(.01, .1*(Z1(I)-Z2(I))))
        ALOGL(I) = LOG (MAX(.01, .1*(Z3(I)-Z4(I))))
c
  100 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TURCOF (ITER)
c
c        Solves for wind speeds at various levels, following the
c        method in the LSX Description section 4.e. Also computes
c        upper and lower-region air-air transfer coefficients and
c        saves them in com1d arrays CU and CL for use by TURVAP,
c        and similarly for the solid-air transfer coefficients
c        SU, SS, SL, SG and SI.
c
c        ITER = current iteration number
c
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION YU(NVEC), YL(NVEC)
#include <comsat.h>
c
c        Set stratification factors for lower and upper regions
c        using values from the previous iteration
c
      XFAC = 1.
      CALL FSTRAT (T34, T12, XFAC, Q34, Q12, Z3, Z2,
     *             ALOGL, ALOGL, ALOG2, U2, RICHL, STRAML, STRAHL, ITER)
c
      CALL FSTRAT (T12, TA,  TFAC, Q12, QA,  Z1, ZA,
     *             ALOGU, ALOGU, ALOGA, UA, RICHU, STRAMU, STRAHU, ITER)
c
      DO 100 I=1,NPOI
c
c          Eliminate C/D from eq (28), tau_l/rho from (26),(27), to get
c          lower-story roughness ALOGL. YL/BDL is (tau_l/rho)/(C+D).
c          (Eqn. numbers correspond to LSX Description section 4.e.)
c
        X = ((ALOG4(I)-ALOGAV(I))/VONK)**2 * BDL(I)
        YL(I) = ((X+1)*EXPL(I) + (X-1)/EXPL(I))
     *        / ((X+1)*EXPL(I) - (X-1)/EXPL(I))
        ALOGL(I) = ALOG3(I) - VONK * SQRT(YL(I)/BDL(I))
c
c          Eliminate tau_l/rho from (24),(25), tau_u/rho and A/B from
c          (22),(23), to get upper-story roughness ALOGU. YU/BDU is
c          (tau_u/rho)/(A+B).
c
        X = ((ALOG2(I)-ALOGL(I))/VONK)**2 * BDU(I) / STRAML(I)
        YU(I) = ((X+1)*EXPU(I) + (X-1)/EXPU(I))
     *        / ((X+1)*EXPU(I) - (X-1)/EXPU(I))
        ALOGU(I) = ALOG1(I) - VONK * SQRT(YU(I)/BDU(I))
  100 CONTINUE
c
      DO 200 I=1,NPOI
c
c          Get TAUU (=tau_u/rho) from (21), A and B from (22),(23),
c          TAUL (=tau_u/rho) from (25), C and D from (26),(27).
c
        TAUU = (UA(I) * VONK/(ALOGA(I)-ALOGU(I)))**2 * STRAMU(I)
        A = 0.5 * TAUU*(YU(I)+1)/BDU(I)
        B = 0.5 * TAUU*(YU(I)-1)/BDU(I)
        TAUL = BDU(I) * (A/EXPU(I) - B*EXPU(I))
        C = 0.5 * TAUL*(YL(I)+1)/BDL(I)
        D = 0.5 * TAUL*(YL(I)-1)/BDL(I)
c
c          Evaluate wind speeds at various levels.
c          max (.001,...) is needed in some cases for IBIS dynamic veg.
c
        U1(I)  = SQRT ( MAX (.001, A + B) )
        U12(I) = SQRT ( MAX (.001, A/EXPHU(I) + B*EXPHU(I)) )
        U2(I)  = SQRT ( MAX (.001, A/EXPU(I)  + B*EXPU(I)) )
        U3(I)  = SQRT ( MAX (.001, C + D) )
        U34(I) = SQRT ( MAX (.001, C/EXPHL(I) + D*EXPHL(I)) )
        U4(I)  = SQRT ( MAX (.001, C/EXPL(I)  + D*EXPL(I)) )
c
c          Compute inverse air-air transfer coeffs. Use of inverse
c          individual coeffs CAI, CBI, CCI, CDI, CEI, CFI avoids
c          divide-by-zero as vegetation vanishes. Combine into
c          upper-region coeff CU from ZA to Z12, and lower-region coeff
c          CL from Z34 to Z12.
c
        CA = UA(I)*STRAHU(I)*VONK**2
     *     / ( (ALOGA(I)-ALOGU(I)) * (ALOGA(I)-ALOG1(I)) )
        CA = MIN (CDMAXA, CA / (1. + CA/CDMAXB))
        CAI = 1. / (RHOA(I)*CA)
 
        CBI = DIU(I) * (Z1(I)-Z12(I)) / (RHOA(I) * 0.5*(U1(I)+U12(I)))
        CCI = DIU(I) * (Z12(I)-Z2(I)) / (RHOA(I) * 0.5*(U12(I)+U2(I)))
        CDI = (ALOG2(I)-ALOGL(I)) * (ALOG2(I)-ALOG3(I))
     *      / ( RHOA(I)*U2(I)*STRAHL(I)*VONK**2 )
        CEI = DIL(I) * (Z3(I)-Z34(I)) / (RHOA(I) * 0.5*(U3(I)+U34(I)))
        CFI = DIL(I) * (Z34(I)-Z4(I)) / (RHOA(I) * 0.5*(U34(I)+U4(I)))
        CU(I) = 1. / (CAI+CBI)
        CL(I) = 1. / (CCI+CDI+CEI)
c
c          Compute air-solid transfer coeffs for upper leaves, upper
c          stems, lower story, soil and snow (SU,SS,SL,SG,SI).
c
        SU(I) = RHOA(I) * CLEAF * SQRT (U12(I)/DLEAF(I,2))
        SS(I) = RHOA(I) * CSTEM * SQRT (U12(I)/DSTEM(I,2))
        SL(I) = RHOA(I) * CGRASS* SQRT (U34(I)/DLEAF(I,1))

#ifndef hydmodel
c================
c       v1, 2.0:
        SG0 = RHOA(I) * U4(I) * (VONK/(ALOG4(I)-ALOGG(I)))**2
        SI0 = RHOA(I) * U4(I) * (VONK/(ALOG4(I)-ALOGI(I)))**2
#else
c===============
c       Ibis2.5:
c===============
c replace above formulations which depend on the log-wind profile
c (which may not work well below a canopy), with empirical formulation
c of Norman's. In the original LSX, turcof.f solves for the winds at the
c various levels from the momentum equations. This gives the transfer
c coefficients for heat and moisture. Heat and moisture eqns are then 
c solved in subroutine turvap. Using the empirical formulation of John
c Norman is not consistent with the earlier solution for u4 (based on a
c logarithmic profile just above the ground. However, this is used here 
c because it improved a lot simulations of the sensible heat flux over 
c the HAPEX-MOBILHY and FIFE sites.

        sg0 = rhoa(i) * (0.004 + 0.012 * u4(i))
        si0 = rhoa(i) * (0.003 + 0.010 * u4(i))
c================
#endif

c         Modify the coefficients to deal with cfi

        SG(I) = 1. / (CFI + 1./SG0)
c**********
c       SG(I) = 1. / (CFI + 1./SG0 + 8400.)
c       SG(I) = 1. / (CFI + 1./SG0 + 16000.)
c**********
        SI(I) = 1. / (CFI + 1./SI0)
c
c          Option to crudely increase transfer coeffs if -ve stratif
c          below trees, to avoid unrealistically(?) large -ve stratif.
c
c       cl(i) = cl(i) * max ( 1., 10.*(t34(i)-t12(i)) )
c       sg(i) = sg(i) * max ( 1., 10.*(tg(i,1)-t34(i)) )
c       si(i) = si(i) * max ( 1., 10.*(ti(i,1)-t34(i)) )
c
  200 CONTINUE
c
c888     If no veg, recalculate coefficients appropriately for a
c        single logarithmic profile, and 2 fictitious levels just
c        above soil/snow surface. These levels are arbitrary but are
c        taken as Z2 and Z4, pre-set in vegdat to a few cm height
c        for bare ground. Use STRAHU from above, which used
c        T12 and ALOGU (ok after first iteration).
c
c     DO 888 I=1,NPOI
c       IF (FU(I).EQ.0. .AND. FL(I).EQ.0.) THEN
c         Z = RHOA(I)*UA(I)*STRAHU(I)*VONK**2 / (ALOGA(I)-ALOGAV(I))
c
c         CA    = Z / (ALOGA(I)-ALOG2(I))
c         CU(I) = RHOA(I)*MIN (CDMAXA, CA / (1. + CA/CDMAXB))
c
c         CL(I) = Z / (ALOG2(I)-ALOG4(I))
c
c         SG(I) = Z / (ALOG4(I)-ALOGG(I))
c         SI(I) = Z / (ALOG4(I)-ALOGI(I))
c
c         ALOGU(I) = ALOGAV(I)
c       ENDIF
c 888 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TURVAP (DTLSX, CO2PPM, ITER, NITER)
c
c        Solves canopy system with linearized implicit sensible heat and
c        moisture fluxes. First, assembles matrix ARR of coeffs in
c        linearized eqns for TU,TS,TL,T12,T34,Q12,Q34,TG,TI,TAVAP,QAVAP
c        and assembles the right-hand sides in RHS. Then calls LINSOLVE
c        to solve this system, passing template MPLATE of zero-locns
c        of ARR. Finally calculates the implied fluxes and stores them
c        for the agcm, soil, snow models and budget calcs.
c
c        DTLSX  = time step for this call
c        CO2PPM = co2 atmospheric amount, ppmv
c        ITER   = current iteration number
c        NITER  = total number of iterations
c
c
#include <compar.h>
#include <comgrd.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION
     *  FRADU(NVEC), FRADS(NVEC), FRADL(NVEC),
     *  FWETUX(NVEC),FWETSX(NVEC),FWETLX(NVEC),
     *  QU(NVEC),    QS(NVEC),    QL(NVEC),    QG(NVEC),    QI(NVEC),
     *  DQU(NVEC),   DQS(NVEC),   DQL(NVEC),   DQG(NVEC),   DQI(NVEC),
     *  TUPRE(NVEC), TSPRE(NVEC), TLPRE(NVEC), TGPRE(NVEC), TIPRE(NVEC),
     *  SUW(NVEC),   SSW(NVEC),   SLW(NVEC),
     *  SUT(NVEC),                SLT(NVEC),   SLT0(NVEC),
     *  SUH(NVEC),   SSH(NVEC),   SLH(NVEC),
     *  EMG(NVEC),   EMI(NVEC),
     *  QGFAC(NVEC)
 
c          For single processing, save the following arrays over turctl
c          iterations. For multiprocessing, use task common to do that.
c
      COMMON  /CTURVMP/
     *  XU(NVEC),    XS(NVEC),    XL(NVEC),
     *  CHUX(NVEC),  CHSX(NVEC),  CHLX(NVEC),  CHGX(NVEC),
     *  WLGX(NVEC),  WIGX(NVEC),
     *  COG(NVEC),   COI(NVEC),
     *  ZIRG(NVEC),  ZIRI(NVEC),
     *  TUOLD(NVEC), TSOLD(NVEC), TLOLD(NVEC), TGOLD(NVEC), TIOLD(NVEC),
     *  WU(NVEC),    WS(NVEC),    WL(NVEC),    WG(NVEC),    WI(NVEC),
     *  WGI(NVEC),
     *  QGFACO(NVEC),
     *  TAVAP(NVEC), QAVAP(NVEC)
!$OMP THREADPRIVATE ( /CTURVMP/ )
c--------
#ifndef omp
c--------
      save  /CTURVMP/
c-----
#endif
c-----
c
      PARAMETER (NQN=11)
      DIMENSION ARR(NVEC,NQN,NQN), RHS(NVEC,NQN), VEC(NVEC,NQN),
     *          MPLATE(NQN,NQN)
c                  tu  ts  tl t12 t34 q12 q34  tg  ti  ta  qa
      DATA MPLATE / 1,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,   !tu
     *              0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,   !ts
     *              0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,   !tl
     *              1,  1,  0,  1,  1,  0,  0,  0,  0,  1,  0,   !t12
     *              0,  0,  1,  1,  1,  0,  0,  1,  1,  0,  0,   !t34
     *              1,  1,  0,  0,  0,  1,  1,  0,  0,  0,  1,   !q12
     *              0,  0,  1,  0,  0,  1,  1,  1,  1,  0,  0,   !q34
     *              0,  0,  0,  0,  1,  0,  1,  1,  0,  0,  0,   !tg
     *              0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  0,   !ti
     *              0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,   !ta
     *              0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1    !qa
     *            /
c
#include <comsat.h>
c
c
c        If first iteration, compute various quantities that will not
c        change in subsequent iterations (arrays are "saved" above)
c
c        (T12,T34,Q12,Q34,FOG12,FOG34,HFOG12,HFOG34,TG,TI for the first
c        iteration are saved via global arrays from the previous gcm
c        timestep. That is worth doing only if the agcm forcing is
c        smoothly varying from timestep to timestep.)
c
c------------------------
      IF (ITER.EQ.1) THEN
c------------------------
 
        DO 20 I=1,NPOI
 
c            Weights for canopy coverages
 
          XU(I) = LAI(I,2) * FU(I)
          XS(I) = SAI(I,2) * FU(I)
          XL(I) = (LAI(I,1)+SAI(I,1)) * FL(I) * (1.-FI(I))
c
c            Specific heats per leaf/stem area
c
          CHUX(I) = CHU + CH2O*WLIQU(I) + CICE*WSNOU(I)
          CHSX(I) = CHS + CH2O*WLIQS(I) + CICE*WSNOS(I)
          CHLX(I) = CHL + CH2O*WLIQL(I) + CICE*WSNOL(I)
c
c           Specific heat of top soil layer, and amounts of
c           liquid and ice in it
c
          CHGX(I) = (   (1.-POROS(I,1))*CSOI(I,1)*RHOSOI(I,1)
     *                + POROS(I,1)*(1.-WISOI(I,1))*WSOI(I,1)*CH2O*RHOW
     *                + POROS(I,1)*WISOI(I,1)*CICE*RHOW
     *              ) * HSOI(I,1)
c
          WLGX(I) = POROS(I,1)*(1.-WISOI(I,1))*WSOI(I,1)*RHOW*HSOI(I,1)
c
          WIGX(I) = POROS(I,1)*WISOI(I,1)*RHOW*HSOI(I,1)
c
c            Conductivity coeffs between ground skin and first layer
c
          COG(I) = CONSOI(I,1) / (0.5*    HSOI(I,1))
          COI(I) = CONSNO(I,1) / (0.5*MAX(HSNO(I,1),HSNOTOP))
c
c            d(IR emitted)/dT for soil and snow surfaces
 
          EMG(I) = EMISOIL
          EMI(I) = EMISNOW
 
          ZIRG(I) = 4.*EMG(I)*STEF*(TG(I)**3)
          ZIRI(I) = 4.*EMI(I)*STEF*(TI(I)**3)
 
c            Save original temps in T*OLD for setting of T* before each
c            iteration, for rhs canopy  heat capacity terms in matrix
c            soln, and for adjustment of canopy temps after each iter.
 
          TUOLD(I) = TU(I)
          TSOLD(I) = TS(I)
          TLOLD(I) = TL(I)
          TGOLD(I) = TG(I)
          TIOLD(I) = TI(I)
 
c            Initialize local agcm air temp and specific humidity
c            that will be changed by implicit calcs in turvap
 
          TAVAP(I) = TA(I)
          QAVAP(I) = QA(I)
 
   20   CONTINUE

c----------
      ENDIF
c----------
 
c        Set implicit/explicit factors W* (0 to 1) for this iteration.
c        W* is 1 for fully implicit, 0 for fully explicit. (For first
c        iteration, IMPEXP and IMPEXP2 set W* to 1.)
c
      CALL IMPEXP (WU, TU, CHUX, WLIQU, WSNOU, ITER, 3)
      CALL IMPEXP (WS, TS, CHSX, WLIQS, WSNOS, ITER, 2)
      CALL IMPEXP (WL, TL, CHLX, WLIQL, WSNOL, ITER, 1)
      CALL IMPEXP (WG, TG, CHGX, WLGX,  WIGX,  ITER, 0)
c
c        Do same for points with snow
c
      CALL IMPEXP2 (WI, TI, TIOLD, ITER, 3)
c
c        Adjust T* for this iteration. We are free to choose them,
c        since they are just the central values about which the
c        equations are linearized. (Heat is conserved in the matrix
c        solution because T*OLD are used for the rhs heat capacities.)
c        Here, let T* represent the previous soln if it was fully
c        implicit, but weight towards T*OLD depending on the amount
c        (1-W*) the previous soln was explicit. (This weighting is
c        necessary for melting/freezing surfaces, for which T* is kept
c        at T*OLD, presumably at or near TMELT.)
c
      DO 60 I=1,NPOI
        TU(I) = WU(I)*TU(I) + (1.-WU(I))*TUOLD(I)
        TS(I) = WS(I)*TS(I) + (1.-WS(I))*TSOLD(I)
        TL(I) = WL(I)*TL(I) + (1.-WL(I))*TLOLD(I)
        TG(I) = WG(I)*TG(I) + (1.-WG(I))*TGOLD(I)
        TI(I) = WI(I)*TI(I) + (1.-WI(I))*TIOLD(I)
        TI(I) = MIN (TI(I), TMELT)
   60 CONTINUE
c
c        Save current "central" values for final flux calculations.
c
      DO 70 I=1,NPOI
        TUPRE(I) = TU(I)
        TSPRE(I) = TS(I)
        TLPRE(I) = TL(I)
        TGPRE(I) = TG(I)
        TIPRE(I) = TI(I)
   70 CONTINUE
c
c        Calculate various terms occurring in the linearized eqns,
c        using values of T12, T34, Q12, Q34, FOG12, FOG34 from
c        the previous iteration.
c
      DO 100 I=1,NPOI
c
c          Specific humidities for canopy and ground, and derivs wrt T
c          for canopy. Limit derivs to avoid -ve implicit q's below,
c          as long as d(temp)s in one iteration are le 10 deg K.
c
        E = ESAT(TU(I))
        QU(I) = QSAT (E, PSURF(I))
        DQU(I) = DSAT (TU(I), PSURF(I))
        DQU(I) = MIN (DQU(I), QU(I)/10.)
c
        E = ESAT(TS(I))
        QS(I) = QSAT (E, PSURF(I))
        DQS(I) = DSAT (TS(I), PSURF(I))
        DQS(I) = MIN (DQS(I), QS(I)/10.)
c
        E = ESAT(TL(I))
        QL(I) = QSAT (E, PSURF(I))
        DQL(I) = DSAT (TL(I), PSURF(I))
        DQL(I) = MIN (DQL(I), QL(I)/10.)
c
        E = ESAT(TG(I))
        QG(I) = QSAT (E, PSURF(I))
        DQG(I) = DSAT (TG(I), PSURF(I))
        DQG(I) = MIN (DQG(I), QG(I)/10.)
c
        E = ESAT(TI(I))
        QI(I) = QSAT (E, PSURF(I))
        DQI(I) = DSAT (TI(I), PSURF(I))
        DQI(I) = MIN (DQI(I), QI(I)/10.)
c
  100 CONTINUE

c          Set QGFACO, factor by which soil surface specific humidity is
c          less than saturation, used in TURVAP. The expression should
c          satisfy 30-min CFL criterion for upper-layer soil moisture
c          for small WSOI(I,1). Philip's formula (commented out below)
c          can violate CFL criterion in desert conditions and cause gcm
c          blowup due to negative soil moistures.
c
c          For each iteration, QGFAC is set to QGFACO, or to 1 if
c          condensation onto soil is anticipated (loop 110)
 
#ifndef hydmodel
        if (iter.eq.1) then
c~~~~~~~~~~~~~~~~~~
        do i=1,npoi
c~~~~~~~~~~~~~~~~~~

cc        Soil:
          ZWTOT = (1.-WISOI(I,1))*WSOI(I,1) + WISOI(I,1)
 
cc        Philips:
cc        ZB = MAX ( MAX(0.,MIN(1.,ZWTOT))**IBEX(I,1), 1.E-8 )
cc        ZE = MAX(-50.,  -GRAV*SUCTION(I,1)/(ZB*RVAP*TSOI(I,1)))
cc        QGFACO(I) = EXP(ZE)
 
cc        v1.02:
cc        QGFACO(I) = MIN (1., (ZWTOT/0.4)**2)
 
cc        Noilhan+Planton,1989,MWR,117,536,eq.(28),
cc        also Mahfouf+Noilhan,1991,JAM,30,1354,Table 2,test #1:
cc        ZZ = MAX(0.,MIN(1.,ZWTOT/.75))
cc        QGFACO(I) = 0.5 * ( 1. - COS(PI*SQRT(ZZ)) )
 
cc        v2.0:
cc        "1-exp(-w**2)":
          ZZA = (ZWTOT/WNORM) **2      !wnorm=.3 for v2.0, .5 for v2.1
          ZZB = (   .8/WNORM) **2
          QGFACO(I) = MIN (1., (1.-EXP(-ZZA)) / (1.-EXP(-ZZB)) )
c~~~~~~~~~~~~
        enddo
c~~~~~~~~~~~~
        endif
#else
c=========
c Ibis2.5:
c=========
c Evaporation from bare soil is calculated using the "beta method"
c (e.g., eqns 5 & 7 of Mahfouf and Noilhan 1991, JAM 30 1354-1365),
c but converted to the "alpha method" (eqns 2 & 3 of M&N), to match
c the structure in IBIS. The conversion from the beta to alpha
c method is through the relationship:
c   alpha * qgs - q34 = beta * (hfac * qgs - q34),
c from which one solves for alpha (which is equal to qgfaco):
c   qgfaco = alpha = (beta * hfac) + (1 - beta)*(q34/qgs)
c
c~~~~~~~~~~~~~~~~~~
        do i=1,npoi
c~~~~~~~~~~~~~~~~~~
c
c first calculate the total saturated fraction at the soil surface
c
          zwtot = max ( 0.01, (1.-wisoi(i,1))*wsoi(i,1) + wisoi(i,1) )
c
c next calculate the matric potential (from eqn 9.3 of Campbell and
c Norman), multiply by gravitational acceleration to get in units
c of J/kg, and calculate the relative humidity at the soil water
c surface (i.e., within the soil matrix), based on thermodynamic
c theory (eqn 4.13 of C&N)
c
          psig = -grav * suction(i,1) * (zwtot ** (-ibex(i,1)))
          hfac = exp (max (-100., psig/(rvap*tg(i))))
c
c then calculate the relative humidity of the air (relative to
c saturation at the soil temperature). Note that if hfac2 > 1
c (which would imply condensation), then qgfac is set to 1
c later in the code (to allow condensation to proceed at the
c "potential rate")
c
          hfac2 = q34(i)/qg(i)
c
c set the "beta" factor and then calculate "alpha" (i.e., qgfaco)
c as the beta-weighted average of the soil water RH and the "air RH"
c First calculate beta_w:
c
          zwopt = 1.0
          zwdry = swilt(i,1)
          betaw = max(0.0, min(1., (zwtot - zwdry)/(zwopt - zwdry)) )
c
c Next convert beta_w to beta_s (see Milly 1992, JClim 5 209-226):
c
          emisoil = 0.95
cc        e      = esat(t34(i))
cc        qs1    = qsat (e, psurf(i))
cc        dqs1   = dqsat (t34(i), qs1) ! ibis2.5
          dqs1   = dsat  (t34(i), psurf(i))
          xnumer = hvap * dqs1
          xdenom = cp(i) + (4.0 * emisoil * stef * (t34(i))**3) / sg(i)
          betafac = xnumer / xdenom
          betas = betaw / (1.0 + betafac * (1.0 - betaw))
c
c Combine hfac and hfac2 into qgfaco ("alpha") using beta_s
c
          qgfaco(i) = betas * hfac + (1. - betas) * hfac2
 
          qgfaco(i) = 0.7*qgfaco(i) + 0.3*1.
c~~~~~~~~~~~~
        enddo
c~~~~~~~~~~~~
#endif

c        Set fractions covered by intercepted h2o to 1 if dew forms.
c        (These FWET*X are used only in TURVAP, and are distinct from
c        the real fractions FWET* that are set in FWETCAL.) They must
c        be exactly 1 if Q12 > QU or Q34 > QL, to zero transpiration
c        by the factor 1-FWET[U,L]X below, so preventing "-ve" transp.
c
c        Similarly, set qgfac, allowing for anticipated dew formation
c        to avoid excessive dew formation (which then infiltrates) onto
c        dry soils.
c
      DO 110 I=1,NPOI
        FWETUX(I) = CVMGT (1., FWETU(I),  Q12(I).GT.QU(I))
        FWETSX(I) = CVMGT (1., FWETS(I),  Q12(I).GT.QS(I))
        FWETLX(I) = CVMGT (1., FWETL(I),  Q34(I).GT.QL(I))
        QGFAC(I)  = CVMGT (1., QGFACO(I), Q34(I).GT.QG(I))
  110 CONTINUE
c
c        Set net absorbed radiative fluxes for canopy components
c
      DO 120 I=1,NPOI
c
        FRADU(I) = CVMGT ( (SOLU(I)+FIRU(I)) / MAX(LAI(I,2),EPSILON),
     *                     0., LAI(I,2).GT.0.)
c
        FRADS(I) = CVMGT ( (SOLS(I)+FIRS(I)) / MAX(SAI(I,2),EPSILON),
     *                     0., SAI(I,2).GT.0.)
c
        Z = LAI(I,1)+SAI(I,1)
        FRADL(I) = CVMGT ( (SOLL(I)+FIRL(I)) / MAX(Z,EPSILON),
     *                     0., Z.GT.0.)
c
  120 CONTINUE
 
c--------------------------------------------------------
c     BATS-SiB type stomatal resistance parameterization:
      IF (PHOTOTYPE.EQ.0) THEN
c--------------------------------------------------------
 
c          Compute upper and lower stomatal resistances RTRAU,RTRAL
c          (in comveg)
c
        CALL TURTRA (QU, QL)
c
c          If necessary, increase stomatal resistances to limit
c          transpiration to the maximum soil uptake rates.
c          Increase diagnostic RSUN[U,L],RSHA[U,L] in same proportion.
c          Also calculate canopy-air moisture transfer coeffs for dry
c          (transpiring) leaf areas. Note that SLT0 is transpiration
c          coeff per unit lower leaf area, and SLT is its average over
c          lower leaf+stem area. (uptot[u,l], upsiz[u,l] were set in
c          setsoi).
c
        DO 200 I=1,NPOI
          XOLDU=RTRAU(I)
          RTRAU(I) = MAX ( RTRAU(I), RHOA(I)
     *      * ( LAI(I,2)*(1.-FWETUX(I))*FSTOMA(I,2)*(QU(I)-Q12(I))
     *          / MAX(UPSIZU(I)*UPTOTU(I),1.E-20) - 1./SU(I) ) )
c
          XOLDL=RTRAL(I)
          RTRAL(I) = MAX ( RTRAL(I), RHOA(I)
     *      * ( LAI(I,1)*(1.-FWETLX(I))*FSTOMA(I,1)*(QL(I)-Q34(I))
     *          / MAX(UPSIZL(I)*UPTOTL(I),1.E-20) - 1./SL(I) ) )
c
          RSUNU(I)=RSUNU(I)*(RTRAU(I)/XOLDU)
          RSHAU(I)=RSHAU(I)*(RTRAU(I)/XOLDU)
          RSUNL(I)=RSUNL(I)*(RTRAL(I)/XOLDL)
          RSHAL(I)=RSHAL(I)*(RTRAL(I)/XOLDL)
c
          SUT(I) = (1.-FWETUX(I))*FSTOMA(I,2)
     *             / (1./SU(I) + RTRAU(I)/RHOA(I))
c
          SLT0(I) = (1.-FWETLX(I))*FSTOMA(I,1)
     *              / (1./SL(I) + RTRAL(I)/RHOA(I))
          SLT(I) = SLT0(I) * LAI(I,1) / MAX(LAI(I,1)+SAI(I,1), EPSILON)
  200   CONTINUE
 
c-------------------------------------------------------
c     IBIS-type photosynthesis and stomatal conductance:
      ELSE
c-------------------------------------------------------
 
c>>>>>>>>>>
#ifdef ibis
c>>>>>>>>>>
 
c          Calculate rates of photosynthesis and stomatal conductances
c          gs[sun,sha][ud,uc,l3,l4] (in comveg)
 
        call stomata (dtlsx, co2ppm, iter, niter)
 
c          Calculate canopy-air moisture transfer coeffs for dry
c          (transpiring) leaf areas. Note that SLT0 is transpiration
c          coeff per unit lower leaf area, and SLT is its average over
c          lower leaf+stem area.
 
        do 210 i=1,npoi
 
c            transpiration from sunlit leaves in the upper canopy:
c            conversion factor .029 = (rh2o/rair)*18*1.e-3 from
c            mol_h2o/m2/s per delta(svp/P) to Kg_h2o/m2/s per delta(q),
c            where 18 = molecular weight of h2o.
 
          zssunu =  (1.-confrac(i))/(1./su(i) + 1./(0.029*gssunud(i)))
     *            +     confrac(i) /(1./su(i) + 1./(0.029*gssunuc(i)))
 
c            transpiration from shaded leaves in the upper canopy:
 
          zsshau =  (1.-confrac(i))/(1./su(i) + 1./(0.029*gsshaud(i)))
     *            +     confrac(i) /(1./su(i) + 1./(0.029*gsshauc(i)))
 
c            transpiration from average upper canopy leaves:
 
          sut(i) = max (   (1.-fwetux(i)) * fstoma(i,2)
     *                   * (sunfu(i)*zssunu+(1.-sunfu(i))*zsshau),
     *                   0.
     *                 )
 
c            transpiration from sunlit leaves in the lower canopy:
 
          zssunl =  (1.-c4frac(i))/(1./sl(i) + 1./(0.029*gssunl3(i)))
     *            +     c4frac(i) /(1./su(i) + 1./(0.029*gssunl4(i)))
 
c            transpiration from shaded leaves in the lower canopy:
 
          zsshal =  (1.-c4frac(i))/(1./sl(i) + 1./(0.029*gsshal3(i)))
     *            +     c4frac(i) /(1./su(i) + 1./(0.029*gsshal4(i)))
 
c            transpiration from average lower canopy leaves:
 
          slt0(i) = max (   (1.-fwetlx(i)) * fstoma(i,1)
     *                    * (sunfl(i)*zssunl+(1.-sunfl(i))*zsshal),
     *                    0.
     *                  )
 
c            averaged over stems and lower canopy leaves:
 
          slt(i) = slt0(i) * lai(i,1) / max(lai(i,1)+sai(i,1), epsilon)
 
 210    continue
 
c>>>>>
#endif
c>>>>>
 
c----------
      ENDIF
c----------
c
c        Calculate canopy-air moisture transfer coeffs for wetted
c        leaf/stem areas. The wetted-area coeffs SUW,SSW,SLW are
c        constrained to be less than what would evaporate 0.8 * the
c        intercepted H2O mass in this timestep (using previous
c        iteration's Q* values). This should virtually eliminate
c        evaporation-overshoots and the need for the "negative
c        intercepted h2o"  correction in steph2o2.
c
      DO 300 I=1,NPOI
        SUW(I) = MIN ( FWETUX(I)*SU(I),
     *                 0.8*(WLIQU(I)+WSNOU(I))
     *                 / MAX (DTLSX*(QU(I)-Q12(I)), EPSILON) )
        SUH(I) = SUW(I) * (      RLIQU(I) *HVAPF(TU(I),TA(I))
     *                     + (1.-RLIQU(I))*HSUBF(TU(I),TA(I)) )
     *         + SUT(I) * HVAPF(TU(I),TA(I))
c
        SSW(I) = MIN ( FWETSX(I)*SS(I),
     *                 0.8*(WLIQS(I)+WSNOS(I))
     *                 / MAX (DTLSX*(QS(I)-Q12(I)), EPSILON) )
        SSH(I) = SSW(I) * (      RLIQS(I) *HVAPF(TS(I),TA(I))
     *                     + (1.-RLIQS(I))*HSUBF(TS(I),TA(I)) )
c
        SLW(I) = MIN ( FWETLX(I)*SL(I),
     *                 0.8*(WLIQL(I)+WSNOL(I))
     *                 / MAX (DTLSX*(QL(I)-Q34(I)), EPSILON) )
        SLH(I) = SLW(I) * (      RLIQL(I) *HVAPF(TL(I),TA(I))
     *                     + (1.-RLIQL(I))*HSUBF(TL(I),TA(I)) )
     *         + SLT(I) * HVAPF(TL(I),TA(I))
  300 CONTINUE
c
c        Set latent heats for fog formation
c
      DO 320 I=1,NPOI
        HFOG12(I) = CVMGT ( HVAPF(T12(I),TA(I)),
     *                      HSUBF(T12(I),TA(I)), T12(I).GE.TMELT )
        HFOG34(I) = CVMGT ( HVAPF(T34(I),TA(I)),
     *                      HSUBF(T34(I),TA(I)), T34(I).GE.TMELT )
  320 CONTINUE
c
c        Estimate rates of fog formation, based on a sub-model for
c        each layer giving a linearized iterative correction towards
c        the equilibrium rate at saturation.
c
      DO 340 I=1,NPOI
        E = ESAT (T12(I))
        QSAT12 = QSAT (E, PSURF(I))
        DSAT12 = DSAT (T12(I), PSURF(I))
        E = ESAT (T34(I))
        QSAT34 = QSAT (E, PSURF(I))
        DSAT34 = DSAT (T34(I), PSURF(I))
 
c-------------------------------------------------------------
c Old (pre version 2.0) method, summing all local coefficients
c-------------------------------------------------------------
c       CT12 = ( CU(I) + CL(I) + XU(I)*SU(I) + XS(I)*SS(I) )
c    *         * CP(I) / HFOG12(I)
c       CQ12 = ( CU(I) + CL(I) + XU(I)*(SUW(I)+SUT(I)) + XS(I)*SSW(I) )
c    *         * DSAT12
c       FOG12(I) = MAX ( 0., FOG12(I) + (Q12(I)-QSAT12)*CT12*CQ12
c    *                                  / ((CT12+CQ12)*DSAT12) )
c
c       CT34 = ( CL(I) + XL(I)*SL(I) + (1.-FI(I))*SG(I)+ FI(I)*SI(I) )
c    *         * CP(I) / HFOG34(I)
c       CQ34 = ( CL(I) + XL(I)*(SLW(I)+SLT(I))
c    *         + (1.-FI(I))*SG(I)+ FI(I)*SI(I) ) * DSAT34
c       FOG34(I) = MAX ( 0., FOG34(I) + (Q34(I)-QSAT34)*CT34*CQ34
c    *                                  / ((CT34+CQ34)*DSAT34) )
c
c--------------------------------------------------------------
c New (version 2.0) method, assuming exchanges with vegetation
c and between canopies follow changes in canopy air and so do
c not enter, and with a cap of .01 to avoid iterative overshoot
c--------------------------------------------------------------
        CXX = MIN (CU(I) + (1.-FI(I))*SG(I) + FL(I)*SI(I), .01)
        C12 = CXX / (1. + DSAT12*HFOG12(I)/CP(I))
        FOG12(I) = MAX (0., FOG12(I) + C12*(Q12(I)-QSAT12))
 
        C34 = CXX / (1. + DSAT34*HFOG34(I)/CP(I))
        FOG34(I) = MAX (0., FOG34(I) + C34*(Q34(I)-QSAT34))
 
c       Turn off fog:
c       fog12(i) = 0.
c       fog34(i) = 0.
 
  340 CONTINUE
c
c        Set the matrix of coefficients and the right-hand sides
c        of the linearized eqns. Split into several do loops to avoid
c        Cray "optimization block size exceeded".
c
      CALL ZERO (ARR, NVEC*NQN*NQN)
      CALL ZERO (RHS, NVEC*NQN)
c
      DO 400 I=1,NPOI
c
c          Upper leaf temperature TU
c
        ARR(I,1,1) = CHUX(I)/DTLSX
     *             + WU(I)*SU(I)*CP(I)
     *             + WU(I)*SUH(I)*DQU(I)
        ARR(I,1,4) = -SU(I)*CP(I)
        ARR(I,1,6) = -SUH(I)
        RHS(I,1) = TUOLD(I)*CHUX(I)/DTLSX
     *           - (1.-WU(I))*SU(I)*CP(I)*TU(I)
     *           - SUH(I) * (QU(I)-WU(I)*DQU(I)*TU(I))
     *           + FRADU(I) - PFLUXU(I)
c
c          Upper stem temperature TS
c
        ARR(I,2,2) = CHSX(I)/DTLSX
     *             + WS(I)*SS(I)*CP(I)
     *             + WS(I)*SSH(I)*DQS(I)
        ARR(I,2,4) = -SS(I)*CP(I)
        ARR(I,2,6) = -SSH(I)
        RHS(I,2) = TSOLD(I)*CHSX(I)/DTLSX
     *           - (1.-WS(I))*SS(I)*CP(I)*TS(I)
     *           - SSH(I) * (QS(I)-WS(I)*DQS(I)*TS(I))
     *           + FRADS(I) - PFLUXS(I)
c
c          Lower veg temperature TL
c
        ARR(I,3,3) = CHLX(I)/DTLSX
     *             + WL(I)*SL(I)*CP(I)
     *             + WL(I)*SLH(I)*DQL(I)
        ARR(I,3,5) = -SL(I)*CP(I)
        ARR(I,3,7) = -SLH(I)
        RHS(I,3) = TLOLD(I)*CHLX(I)/DTLSX
     *           - (1.-WL(I))*SL(I)*CP(I)*TL(I)
     *           - SLH(I) * (QL(I)-WL(I)*DQL(I)*TL(I))
     *           + FRADL(I) - PFLUXL(I)
c
  400 CONTINUE
c
      DO 410 I=1,NPOI
c
c          Upper air temperature T12
c
        ARR(I,4,1) = -WU(I)*XU(I)*SU(I)
        ARR(I,4,2) = -WS(I)*XS(I)*SS(I)
        ARR(I,4,4) = CU(I) + CL(I) + XU(I)*SU(I) + XS(I)*SS(I)
        ARR(I,4,5) = -CL(I)
        ARR(I,4,10)= -CU(I)*TFAC
        RHS(I,4) = (1.-WU(I))*XU(I)*SU(I)*TU(I)
     *           + (1.-WS(I))*XS(I)*SS(I)*TS(I)
     *           + HFOG12(I)*FOG12(I)/CP(I)
c
c          Lower air temperature T34
c
        ARR(I,5,3) = -WL(I)*XL(I)*SL(I)
        ARR(I,5,4) = -CL(I)
        ARR(I,5,5) = CL(I) + XL(I)*SL(I)
     *             + (1.-FI(I))*SG(I) + FI(I)*SI(I)
        ARR(I,5,8) = -WG(I)*(1.-FI(I))*SG(I)
        ARR(I,5,9) = -WI(I)*FI(I)*SI(I)
        RHS(I,5) = (1.-WL(I))*XL(I)     *SL(I)*TL(I)
     *           + (1.-WG(I))*(1.-FI(I))*SG(I)*TG(I)
     *           + (1.-WI(I))*FI(I)     *SI(I)*TI(I)
     *           + HFOG34(I)*FOG34(I)/CP(I)
c
c          Upper air specific humidity Q12
c
        ARR(I,6,1) = -WU(I)*XU(I)*(SUW(I)+SUT(I))*DQU(I)
        ARR(I,6,2) = -WS(I)*XS(I)*SSW(I)*DQS(I)
        ARR(I,6,6) = CU(I) + CL(I)
     *             + XU(I)*(SUW(I)+SUT(I)) + XS(I)*SSW(I)
        ARR(I,6,7) = -CL(I)
        ARR(I,6,11)= -CU(I)
        RHS(I,6) = XU(I) * (SUW(I)+SUT(I)) * (QU(I)-WU(I)*DQU(I)*TU(I))
     *           + XS(I) *  SSW(I)         * (QS(I)-WS(I)*DQS(I)*TS(I))
     *           - FOG12(I)
c
c          Lower air specific humidity Q34
c
        ARR(I,7,3) = -WL(I)*XL(I)*(SLW(I)+SLT(I))*DQL(I)
        ARR(I,7,6) = -CL(I)
        ARR(I,7,7) = CL(I) + XL(I)*(SLW(I)+SLT(I))
     *             + (1.-FI(I))*SG(I) +FI(I)*SI(I)
        ARR(I,7,8) = -WG(I)*(1.-FI(I))*SG(I)*QGFAC(I)*DQG(I)
        ARR(I,7,9) = -WI(I)*FI(I)*SI(I)*DQI(I)
        RHS(I,7)= XL(I)     *(SLW(I)+SLT(I))*(QL(I)-WL(I)*DQL(I)*TL(I))
     *          + (1.-FI(I))*SG(I)*QGFAC(I) *(QG(I)-WG(I)*DQG(I)*TG(I))
     *          + FI(I)     *SI(I)          *(QI(I)-WI(I)*DQI(I)*TI(I))
     *          - FOG34(I)
c
  410 CONTINUE
c
      DO 420 I=1,NPOI
c
c          Soil skin temperature
 
c          (There is no WG in this eqn since it solves for a fully
c          implicit TG. WG can be thought of as the fractional soil
c          area using a fully implicit soln, and 1-WG as that using a
c          fully explicit soln. The combined soil temperature is felt
c          by the lower air, so WG occurs in the T34,Q34 eqns above.)
c
        ARR(I,8,5) = -SG(I)*CP(I)
        ARR(I,8,7) = -SG(I)*HVASUG(I)
        ARR(I,8,8) = SG(I)*CP(I) + SG(I)*HVASUG(I)*QGFAC(I)*DQG(I)
     *             + COG(I) + ZIRG(I)
        RHS(I,8) = -SG(I)*HVASUG(I)*QGFAC(I)*(QG(I)-DQG(I)*TG(I))
     *           + COG(I)*TSOI(I,1)
     *           + SOLG(I) + FIRG(I) + ZIRG(I)*TGOLD(I)
c
c          Snow skin temperature
c
c          (There is no WI here, for the same reason as for WG above.)
c
        ARR(I,9,5) = -SI(I)*CP(I)
        ARR(I,9,7) = -SI(I)*HVASUI(I)
        ARR(I,9,9) = SI(I)*CP(I) + SI(I)*HVASUI(I)*DQI(I)
     *             + COI(I) + ZIRI(I)
        RHS(I,9) = -SI(I)*HVASUI(I)*(QI(I)-DQI(I)*TI(I))
     *           + COI(I)*TSNO(I,1)
     *           + SOLI(I) + FIRI(I) + ZIRI(I)*TIOLD(I)
c
c          AGCM air temperature TAVAP
c
        ARR(I,10,4)  = -CU(I)
        ARR(I,10,10) = THICA(I)/DTLSX
     *               + CU(I)*TFAC
        RHS(I,10)    = TA(I)*THICA(I)/DTLSX
c
c          AGCM specific humidity QAVAP
c
        ARR(I,11,6)  = -CU(I)
        ARR(I,11,11) = THICA(I)/DTLSX
     *               + CU(I)
        RHS(I,11)    = QA(I)*THICA(I)/DTLSX
c
  420 CONTINUE
c
c        Solve the systems of equations
c
      CALL LINSOLVE (ARR, RHS, VEC, MPLATE, NPOI, NQN, IOTERM2)
c
c        Optional printout to monitor parts of the iterative scheme
c
#ifdef gcm
      zmonit = 0.
#else
c     zmonit = 1.
      zmonit = 0.
#endif
      iumonit = 56
      if (zmonit.ne.0.) then
        call findpoi(ip)
        if (ip.gt.0) then
c
          if (iter.eq.1) write(iumonit,500)
  500     format(/' i',
     *            '  fwu fwux   wu',
     *            '  tavap     tu     tg    t12    t34',
     *            '  strahu  strahl',
     *            '     qavap        qu        ql       q12       q34',
c    *            '   fog12   fog34 hf12 hf34',
     *            '     ua     u1    u12     u2     u3    u34     u4'
     *            )
c
          write(iumonit,502)
     *     iter,
     *     fwetu(ip), fwetux(ip), wu(ip),
     *     tavap(ip)-tmelt, tu(ip)-tmelt, tg(ip)-tmelt,
     *     t12(ip)-tmelt, t34(ip)-tmelt,
     *     strahu(ip), strahl(ip),
     *     qavap(ip), qu(ip), ql(ip), q12(ip), q34(ip),
c    *     fog12(ip)*.864e5, fog34(ip)*.864e5,
c    *     hfog12(ip)/1.e6, hfog34(ip)/1.e6,
     *     ua(ip), u1(ip), u12(ip), u2(ip), u3(ip), u34(ip),  u4(ip)
  502     format(1x, i1,
     *           3f5.2,
     *           5f7.2,
     *           2f8.3,
     *           5f10.6,
c    *           2f8.3, 2f5.2,
     *           7f7.3
     *           )
 
c
          if (iter.eq.niter) write(iumonit,504)
     *     vec(ip,10)-tmelt, vec(ip,1)-tmelt, vec(ip,3)-tmelt,
     *     vec(ip,4)-tmelt, vec(ip,5)-tmelt,
     *     vec(ip,11), vec(ip,6), vec(ip,7)
  504     format(1x, 1x, 3(5x), 5f7.2, 2(8x), f10.6, 2(10x),2f10.6 / )
c
        endif
      endif
c
c        Copy this iteration's solution to T*, Q12, Q34, TAVAP, QAVAP
c
      DO 600 I=1,NPOI
        TU(I)    = VEC(I,1)
        TS(I)    = VEC(I,2)
        TL(I)    = VEC(I,3)
        T12(I)   = VEC(I,4)
        T34(I)   = VEC(I,5)
        Q12(I)   = VEC(I,6)
        Q34(I)   = VEC(I,7)
        TG(I)    = VEC(I,8)
        TI(I)    = VEC(I,9)
        TAVAP(I) = VEC(I,10)
        QAVAP(I) = VEC(I,11)
  600 CONTINUE
c
c        All done except for final flux calculations,
c        so loop back for the next iteration (except the last)
c
      IF (ITER.LT.NITER) RETURN
c
c          Evaluate sensible heat and moisture fluxes (per unit
c          leaf/stem/snow-free/snow-covered area as appropriate).
c
      DO 700 I=1,NPOI
c
c          Sensible heat fluxes
c
        FSENA(I) = CP(I) * CU(I) * (TAVAP(I)*TFAC-T12(I))
c
        TUAV = WU(I)*TU(I) + (1.-WU(I))*TUPRE(I)
        FSENU(I) = CP(I) * SU(I) * (TUAV-T12(I))
c
        TSAV = WS(I)*TS(I) + (1.-WS(I))*TSPRE(I)
        FSENS(I) = CP(I) * SS(I) * (TSAV-T12(I))
c
        TLAV = WL(I)*TL(I) + (1.-WL(I))*TLPRE(I)
        FSENL(I) = CP(I) * SL(I) * (TLAV-T34(I))
c
        TGAV = WG(I)*TG(I) + (1.-WG(I))*TGPRE(I)
        FSENG(I) = CP(I) * SG(I) * (TGAV-T34(I))
c
        TIAV = WI(I)*TI(I) + (1.-WI(I))*TIPRE(I)
        FSENI(I) = CP(I) * SI(I) * (TIAV-T34(I))
c
c          Moisture fluxes
c
        FVAPA(I)  = CU(I) * (QAVAP(I)-Q12(I))
c
        QUAV = QU(I) + WU(I)*DQU(I)*(TU(I)-TUPRE(I))
        FVAPUW(I) = SUW(I) * (QUAV-Q12(I))
        FVAPUT(I) = SUT(I) * (QUAV-Q12(I))
c
        QSAV = QS(I) + WS(I)*DQS(I)*(TS(I)-TSPRE(I))
        FVAPS(I)  = SSW(I) * (QSAV-Q12(I))
c
        QLAV = QL(I) + WL(I)*DQL(I)*(TL(I)-TLPRE(I))
        FVAPLW(I) = SLW(I) * (QLAV-Q34(I))
        FVAPLT(I) = SLT0(I)* (QLAV-Q34(I))
c
        QGAV = QG(I) + WG(I)*DQG(I)*(TG(I)-TGPRE(I))
        FVAPG(I)  = SG(I) * (QGFAC(I)*QGAV - Q34(I))
c
        QIAV = QI(I) + WI(I)*DQI(I)*(TI(I)-TIPRE(I))
        FVAPI(I)  = SI(I) * (QIAV-Q34(I))
c
c          IR fluxes (correct due to implicit terms above for soil/snow
c          skin temperature...give correction directly to agcm)
c
        FIRG(I) = FIRG(I) - WG(I)*ZIRG(I)*(TG(I)-TGOLD(I))
        FIRI(I) = FIRI(I) - WI(I)*ZIRI(I)*(TI(I)-TIOLD(I))
        FIRB(I) = FIRB(I) + (1.-FI(I))*WG(I)*ZIRG(I)*(TG(I)-TGOLD(I))
     *                    +     FI(I) *WI(I)*ZIRI(I)*(TI(I)-TIOLD(I))
c
c          Impose final snow constraints on (diagnostic) skin
c          temp TI (for instance, for next iterations's IRRAD).
c          Don't do earlier since must use matrix result for flux calcs
c
        TI(I) = MIN (TI(I), TMELT)
c
c       Diagnostic warning (local machine only):
#if ! defined (gcm) && ! defined (hydmodel)
        if (quav.lt.0. .or. qsav.lt.0. .or. qlav.lt.0. .or.
     *      qgav.lt.0. .or. qiav.lt.0. .or.
     *      q12(i).lt.0. .or. q34(i).lt.0.) then
          write(ioterm,702)
     *      quav,qsav,qlav,qgav,qiav,q12(i),q34(i)
          write(iutaba,702)
     *      quav,qsav,qlav,qgav,qiav,q12(i),q34(i)
  702     format(' Warning TURVAP 702:',
     *      ' q[u,s,l,g,i]av, q12, q34: ',7e12.4)
        endif
#endif
c
  700 CONTINUE
c
c        Do small energy-conserving correction since latent heats for
c        evap of intercepted h2o and transpiration were calculated
c        by HVAPF,HSUBF using pre-temps T[U,S,L]PRE, but energy eqns
c        have changed temps to T[U,S,L] including the to-be-removed h2o.
c        (Intercepted evap will be removed in STEPH2O, and transpired
c        h2o will be removed from soil in SOILH2O assuming transpiration
c        occurs at T[U,L].)
c
      DO 750 I=1,NPOI
        X  = CHU + CH2O*WLIQU(I) + CICE*WSNOU(I)
        DX = (   FVAPUW(I)*(RLIQU(I)*CH2O+(1.-RLIQU(I))*CICE)
     *         + FVAPUT(I)*CH2O ) * DTLSX
        DX = MIN (DX, 0.5*X)
        TU(I) = (X*TU(I) - DX*TUPRE(I)) / (X-DX)
c
        X  = CHS + CH2O*WLIQS(I) + CICE*WSNOS(I)
        DX = (   FVAPS(I)*(RLIQS(I)*CH2O+(1.-RLIQS(I))*CICE) ) * DTLSX
        DX = MIN (DX, 0.5*X)
        TS(I) = (X*TS(I) - DX*TSPRE(I)) / (X-DX)
c
        X  = CHL + CH2O*WLIQL(I) + CICE*WSNOL(I)
        DX = (   FVAPLW(I)*(RLIQL(I)*CH2O+(1.-RLIQL(I))*CICE)
     *         + FVAPLT(I)*CH2O ) * DTLSX
        DX = MIN (DX, 0.5*X)
        TL(I) = (X*TL(I) - DX*TLPRE(I)) / (X-DX)
  750 CONTINUE
c
c       Set UPSOI[U,L], the actual soil water uptake rates from each
c       soil layer due to transpiration in the upper and lower stories,
c       for the soil model.
c
      DO 800 K=1,NSOILAY
        DO 802 I=1,NPOI
          UPSOIU(I,K) = FVAPUT(I) * LAI(I,2) * FU(I)
     *                * UPLAYU(I,K)/MAX(UPTOTU(I),EPSILON)
          UPSOIL(I,K) = FVAPLT(I) * LAI(I,1) * FL(I) * (1.-FI(I))
     *                * UPLAYL(I,K)/MAX(UPTOTL(I),EPSILON)
  802   CONTINUE
  800 CONTINUE
c
c       Set net evaporation from intercepted H2O, net evap from surface,
c       and net transpiration rate (GINVAP,GSUVAP,GTRANS*) in commud for
c       history and diagnostics.
c
      DO 810 I=1,NPOI
       GINVAP(I) = FVAPUW(I) * LAI(I,2) * FU(I)
     *           + FVAPS (I) * SAI(I,2) * FU(I)
     *           + FVAPLW(I) * (LAI(I,1)+SAI(I,1)) * FL(I) * (1.-FI(I))
       GSUVAP(I) = FVAPG(I)  * (1.-FI(I))
     *           + FVAPI(I)  * FI(I)
       GTRANSU(I)= FVAPUT(I) * LAI(I,2) * FU(I)
       GTRANSL(I)= FVAPLT(I) * LAI(I,1) * FL(I) * (1.-FI(I))
  810 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TURTRA (QU, QL)
c
c        Computes stomatal resistance coefficients RTRAU, RTRAL in
c        comveg for the upper and lower story, in sec/meter.
c
c        QU    = saturated specific humidity of upper leaves (supplied)
c        QL    = saturated specific humidity of lower leaves (supplied)
c
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION QU(NVEC), QL(NVEC)
c
c
      IB = 1
      DO 100 I=1,NPOI
c
c          PAR dependence (use SUNP[U,L],SHAP[U,L], SUNF[U,L] estimated
c          in SOLARF (IB=1 for visible waveband)
c
c       upper-story sunlit and shaded leaves:
        RPARSUNU = (SUNPU(I)/PARCON(I,2) + 1.)
     *           / (SUNPU(I)/PARCON(I,2) + RSMIN(I,2)/RSMAX(I,2))
        RPARSHAU = (SHAPU(I)/PARCON(I,2) + 1.)
     *           / (SHAPU(I)/PARCON(I,2) + RSMIN(I,2)/RSMAX(I,2))
c       lower-story sunlit and shaded leaves:
        RPARSUNL = (SUNPL(I)/PARCON(I,1) + 1.)
     *           / (SUNPL(I)/PARCON(I,1) + RSMIN(I,1)/RSMAX(I,1))
        RPARSHAL = (SHAPL(I)/PARCON(I,1) + 1.)
     *           / (SHAPL(I)/PARCON(I,1) + RSMIN(I,1)/RSMAX(I,1))
c
c          Leaf temperature dependence (none in v2)
c
c       RTEMU =
c    *    1. / MAX (1. - TEMCON(I,2) * (TEMOPT(I,2)-TU(I))**2, EPSILON)
c       RTEML =
c    *    1. / MAX (1. - TEMCON(I,1) * (TEMOPT(I,1)-TL(I))**2, EPSILON)
        RTEMU = 1.
        RTEML = 1.
c
c          Vapor pressure deficit dependence (none in v2)
c
c       VPDU = MAX ( 0., (PSURF(I)/0.622)*(QU(I)-Q12(I)) )
c       VPDL = MAX ( 0., (PSURF(I)/0.622)*(QL(I)-Q34(I)) )
c       RVPDU = 1. / MAX ( 1. - VPDU/VPDCON(I,2), EPSILON )
c       RVPDL = 1. / MAX ( 1. - VPDL/VPDCON(I,1), EPSILON )
        RVPDU = 1.
        RVPDL = 1.
c
c          Combine the above dependencies (RSUN[U,L],RSHA[U,L] are
c          diagnostic 1D arrays)
c
c       upper story sunlit and shaded leaves:
        RSUNU(I) = MIN ( RSMAX(I,2), RSMIN(I,2)*RPARSUNU*RTEMU*RVPDU )
        RSHAU(I) = MIN ( RSMAX(I,2), RSMIN(I,2)*RPARSHAU*RTEMU*RVPDU )
c       lower story sunlit and shaded leaves:
        RSUNL(I) = MIN ( RSMAX(I,1), RSMIN(I,1)*RPARSUNL*RTEML*RVPDL )
        RSHAL(I) = MIN ( RSMAX(I,1), RSMIN(I,1)*RPARSHAL*RTEML*RVPDL )
c
c          Average over sunlit and shaded leaf fractions (averaging
c          conductances in parallel, not resistances)
c
        RTRAU(I) = 1. / ( SUNFU(I)/RSUNU(I) + (1.-SUNFU(I))/RSHAU(I) )
        RTRAL(I) = 1. / ( SUNFL(I)/RSUNL(I) + (1.-SUNFL(I))/RSHAL(I) )
c
c          Multiply by user-specified (lsx-namelist) factor
c
        RSUNU(I) = STOMFAC*RSUNU(I)
        RSHAU(I) = STOMFAC*RSHAU(I)
        RSUNL(I) = STOMFAC*RSUNL(I)
        RSHAL(I) = STOMFAC*RSHAL(I)
c
        RTRAU(I) = STOMFAC*RTRAU(I)
        RTRAL(I) = STOMFAC*RTRAL(I)
c
  100 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE FSTRAT (TB, TT, TTFAC, QB, QT, ZB, ZT,
     *                   ALBM, ALBH, ALT, U, RICH, STRAM, STRAH, ITER)
c
c        Computes mixing-length stratification correction factors
c        for momentum and heat/vapor, for current 1D strip, using
c        parameterizations in Louis (1979),BLM,17,187. First computes
c        Richardson Numbers. Sets an upper limit to Richardson Numbers
c        so lower-veg winds don't become vanishingly small in very
c        stable conditions (cf, Carson and Richards,1978,BLM,14,68).
c
c        System (I) is as in Louis(1979). System (VI) is improved as
c        described in Louis(1982), ECMWF Workshop on Planetary Boundary
c        Layer Parameterizations,November 1981,59-79 (QC880.4 B65W619).
c
c        TB    = bottom temperature (supplied)
c        TT    = top temperature (supplied)
c        TTFAC = pot. temp factor for TTOP (relative to bottom,supplied)
c        QB    = bottom specific humidity (supplied)
c        QT    = top specific humidity (supplied)
c        ZB    = height of bottom (supplied)
c        ZT    = height of top (supplied)
c        ALBM  = log (bottom roughness length) for momentum (supplied)
c        ALBH  = log (bottom roughness length) for heat/h2o (supplied)
c        ALT   = log (z at top) (supplied)
c        U     = wind speed at top (supplied)
c        RICH  = Richardson Number (returned)
c        STRAM = stratification factor for momentum (returned)
c        STRAH = stratification factor for heat/vap (returned)
c
#include <compar.h>
#include <comgrd.h>
      DIMENSION
     * TB(NVEC),   TT(NVEC),   QB(NVEC),    QT(NVEC),
     * ZB(NVEC),   ZT(NVEC),   ALBM(NVEC),  ALBH(NVEC),
     * ALT(NVEC),  U(NVEC),    RICH(NVEC),  STRAM(NVEC), STRAH(NVEC)
      DIMENSION
     * INDP(NVEC), INDQ(NVEC), STRAMX(NVEC), STRAHX(NVEC)
c
c        Calculate Richardson Numbers, with lower and upper bounds
c
      DO 100 I=1,NPOI
        ZHT = TT(I)*TTFAC*(1.+.622*QT(I))
        ZHB = TB(I)*      (1.+.622*QB(I))
        RICH(I) = GRAV * MAX (ZT(I)-ZB(I), 0.)
     *            * (ZHT-ZHB) / (0.5*(ZHT+ZHB) * (MAX(U(I),.01))**2)
        RICH(I) = MAX (-2., MIN (1., RICH(I)))
  100 CONTINUE
c
c        Set up indices for points with negative or positive Ri
c
      CALL WHENFLE (NPOI, RICH, 1, 0., INDP, NP)
      CALL WHENFGT (NPOI, RICH, 1, 0., INDQ, NQ)
c
c        Calculate momentum and heat/vapor factors for negative Ri.
c
      IF (NP.GT.0) THEN
        DO 300 J=1,NP
          I = INDP(J)
          XM = MAX (ALT(I)-ALBM(I), .5)
          XH = MAX (ALT(I)-ALBH(I), .5)
          YM = (VONK/XM)**2 * EXP (0.5*XM) * SQRT(-RICH(I))
          YH = (VONK/XH)**2 * EXP (0.5*XH) * SQRT(-RICH(I))
 
c         System (I):
c         STRAMX(I) =   1. - 9.4*RICH(I) / (1. + 9.4*7.4*YM)
c         STRAHX(I) = ( 1. - 9.4*RICH(I) / (1. + 9.4*5.3*YH) ) / 0.74
 
c         System (VI):
          STRAMX(I) =   1. - 2*5*RICH(I) / (1. + 75*YM)
          STRAHX(I) =   1. - 3*5*RICH(I) / (1. + 75*YH)
  300   CONTINUE
      ENDIF
c
c        Calculate momentum and heat/vapor factors for positive Ri.
c
      IF (NQ.GT.0) THEN
!DIR$ IVDEP
        DO 400 J=1,NQ
          I = INDQ(J)
 
c         System (I):
c         STRAMX(I) = 1. / (1.+4.7*RICH(I))**2
c         STRAHX(I) = STRAMX(I) / 0.74
 
c         System (VI):
          Z = SQRT(1.+5*RICH(I))
          STRAMX(I) = 1. / (1.+ 2*5*RICH(I)/Z)
          STRAHX(I) = 1. / (1.+ 3*5*RICH(I)*Z)
  400   CONTINUE
      ENDIF
c
c        Except for the first iteration, weight results with the
c        previous iteration's values. This improves convergence by
c        avoiding flip-flop between stable/unstable stratif, eg,
c        with cold upper air and the lower surface being heated by
c        solar radiation.
c
      IF (ITER.EQ.1) THEN
        DO 500 I=1,NPOI
          STRAM(I) = STRAMX(I)
          STRAH(I) = STRAHX(I)
  500   CONTINUE
      ELSE
        W = 0.5
        DO 510 I=1,NPOI
          STRAM(I) = W*STRAMX(I) + (1.-W)*STRAM(I)
          STRAH(I) = W*STRAHX(I) + (1.-W)*STRAH(I)
  510   CONTINUE
      ENDIF
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE IMPEXP (WIMP, T, CH, WLIQ, WSNO, ITER, IWHICH)
c
c        Sets the implicit vs explicit fraction in TURVAP calcs for
c        upper leaves, upper stems or lower veg. This is to account for
c        temperatures of freezing/melting intercepted h2o constrained
c        at the melt point. If a purely implicit calc is used for such
c        a surface, the predicted temperature would be nearly the atmos
c        equil temp with little sensible heat input, so the amount of
c        freezing or melting is underestimated. However, if a purely
c        explicit calc is used with only a small amount of intercepted
c        h2o, the heat exchange can melt/freeze all the h2o and cause
c        an unrealistic huge change in the veg temp. The algorithm
c        below attempts to avoid both pitfalls.
c
c        Also used for soil-skin equation, with CH,WLIQ,WSNO being the
c        heat capacity, liquid amount and ice amount in the top soil
c        layer
c
c        WIMP   = implicit/explicit fraction (0 to 1) (returned)
c        T      = temperature of veg/soil (prev iteration's soln) (supp)
c        CH     = heat capacity of veg/soil (supplied)
c        WLIQ   = veg-intercepted/soil liquid (supplied)
c        WSNO   = veg-intercepted/soil snow (supplied)
c        ITER   = current iteration number (supplied)
c        IWHICH : 3=upper leaves, 2=upper stems, 1=lower veg, 0=soil
c                 (supplied, just for debug/diagnostics)
c
#include <compar.h>
#include <comgrd.h>
      DIMENSION
     *  WIMP(NVEC), T(NVEC), CH(NVEC), WLIQ(NVEC), WSNO(NVEC)
c
c
c        For first iteration, set WIMP to fully implicit, and return
c
      IF (ITER.EQ.1) THEN
        DO 50 I=1,NPOI
          WIMP(I) = 1.
   50   CONTINUE
        RETURN
      ENDIF
c
c        For second and subsequent iterations, estimate WIMP based on
c        the previous iterations's WIMP and its resulting T.
c
c        Calculate H, the "overshoot" heat available to melt any snow
c        or freeze any liquid. Then the explicit fraction is taken to
c        be the ratio of H to the existing h2o's latent heat (ie, 100%
c        explicit calculation if not all of the h2o would be melted or
c        frozen). So WINEW, the implicit amount, is 1 - that ratio.
c        But since we are using the previous iteration's T* results
c        for the next iteration, to ensure convergence we need to damp
c        the returned estimate WIMP by averaging WINEW with the
c        previous estimate. This works reasonably well even with a
c        small number of iterations (3), since for instance with large
c        amounts of h2o so that WIMP should be 0., a good amount of
c        h2o is melted or frozen with WIMP = .25. (WIMP's are "saved"
c        in TURVAP).
c
      DO 100 I=1,NPOI
c
        H = CH(I) * (T(I)-TMELT)
 
c       IF (H.GT. EPSILON) THEN
c         WINEW = 1. - MIN(1., HFUS*WSNO(I)/H)
c       ELSE IF (H.LT.-EPSILON) THEN
c         WINEW = 1. - MIN(1., HFUS*WLIQ(I)/(-H))
c       ELSE
c         WINEW = 1.
c       ENDIF
        Z = MAX (ABS(H), EPSILON)
        WINEW = CVMGT (1. - MIN(1., HFUS*WSNO(I)/Z), 1., H.GT.EPSILON)
        WINEW = CVMGT (1. - MIN(1., HFUS*WLIQ(I)/Z), 1., H.LT.-EPSILON)
c
        WIMP(I) = 0.5*(WIMP(I) + WINEW)
c
  100 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE IMPEXP2 (WIMP, T, TOLD, ITER, IWHICH)
c
c        Sets the implicit vs explicit fraction in TURVAP calcs for
c        icesheet, seaice or snow skin temperatures, to account for
c        temperatures of freezing/melting surfaces being constrained at
c        the melt point. Unlike IMPEXP, don't have to allow for all h2o
c        vanishing within the timestep.
c
c        WIMP = implicit fraction (0 to 1) (returned)
c        T    = temperature of skin layer (supplied)
c        TOLD = temperature of first ice/snow layer (supplied)
c        ITER = current iteration number (supplied)
c        IWHICH : 1=icesheet, 2=seaice, 3=snow
c                 (supplied, just for debug/diagnostics)
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
      DIMENSION
     *  WIMP(NVEC), T(NVEC), TOLD(NVEC)
c
c        For first iteration, set WIMP to fully implicit, and return.
c
      IF (ITER.EQ.1) THEN
        DO 50 I=1,NPOI
          WIMP(I) = 1.
   50   CONTINUE
        RETURN
      ENDIF
c
c        For second and subsequent iterations, estimate WIMP based on
c        the previous iterations's T result...how much T is above
c        TMELT and how much TOLD is below TMELT.
c        (Note that we can assume TOLD le TMELT.)
c
      DO 100 I=1,NPOI
        WIMP(I) = CVMGT ( (TMELT-TOLD(I)) / MAX(T(I)-TOLD(I),EPSILON),
     *                    WIMP(I),
     *                    T(I)-TOLD(I).GT.EPSILON )
 
        WIMP(I) = MAX (0., MIN (1., WIMP(I)))
  100 CONTINUE
c
      RETURN
      END
