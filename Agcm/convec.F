c -- $Id: convec.F,v 1.2 2001/07/03 21:29:17 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/convec.F,v $
c -- $Name:  $

C
C ***************************** CONVAD ********************************
C
      SUBROUTINE CONVAD 
     *  (YU,     YV,      YT,     YQ,
     *   YLWC_C, YLWC_A,  YLWC_S, YCLD,
     *   FPLUMT, FPLUMQ,  FPLUMU, FPLUMV,
     *   TRACER, TRAFLUX, YPM,    YPTHIC, YPMCAP,
     *   PLUMEPREC, HEIMIN, HEIMAX,  JJ
#ifdef isotrac
     *  ,yq_o, yq_h, ylwc_c_o, ylwc_c_h, ylwc_a_o,
     *   ylwc_a_h, ylwc_s_o, ylwc_s_h, fplumq_o, fplumq_h,
     *   plumepr_o, plumepr_h,
     *   rtrq, rtrfq, plumepr_t
#endif
     *  )

C----------------------------------------------------------------------
C
CL            DO PLUME CONVECTION, STRATIFORM ADJUSTMENT,
CL            AND CLOUD CONVERSION, RE-EVAPORATION OF FALLING PRECIP,
CL            FOR ONE LATITUDE STRIP (JJ = LATITUDE INDEX)
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C---------------------------------------------------------------------
#include <COMMUN.h>
C---------------------------------------------------------------------
#include <COMSDT.h>
C---------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      DIMENSION
     *   YU(NLON,NOREC,NLEV),      YV(NLON,NOREC,NLEV),
     *   YT(NLON,NOREC,NLEV),      YQ(NLON,NOREC,NLEV),
     *   YLWC_C(NLON,NOREC,NLEV),  YLWC_A(NLON,NOREC,NLEV),
     *   YLWC_S(NLON,NOREC,NLEV),  YCLD(NLON,NOREC,NLEV),
     *   FPLUMT(NLON,NOREC,NLEV),  FPLUMQ(NLON,NOREC,NLEV),
     *   FPLUMU(NLON,NOREC,NLEV),  FPLUMV(NLON,NOREC,NLEV),
     *   TRACER(NLON,NOREC,NLEV,NTRACE), TRAFLUX(NLON,NOREC,NTRACB),
     *   YPM(NLON,NOREC,NLEV),
     *   YPTHIC(NLON,NOREC,NLEV),  YPMCAP(NLON,NOREC,NLEV),
     *   PLUMEPREC(NLON,NOREC,NLEV),
     *   HEIMIN(NLON,NOREC),       HEIMAX(NLON,NOREC)
#ifdef isotrac
      dimension
     *   yq_o(nlon,norec,nlev),     yq_h(nlon,norec,nlev),
     *   ylwc_c_o(nlon,norec,nlev), ylwc_c_h(nlon,norec,nlev),
     *   ylwc_a_o(nlon,norec,nlev), ylwc_a_h(nlon,norec,nlev),
     *   ylwc_s_o(nlon,norec,nlev), ylwc_s_h(nlon,norec,nlev),
     *   fplumq_o(nlon,norec,nlev), fplumq_h(nlon,norec,nlev),
     *   plumepr_o(nlon,norec,nlev),plumepr_h(nlon,norec,nlev),
     *   za_i(nlon,nlev),           zh_i(nlon,nlev),
     *   rtrq(nlon,norec,nlev,nrtr),rtrfq(nlon,norec,nlev,nrtr),
     *   plumepr_t(nlon,norec,nlev,nrtr)
#endif
C----------------------------------------------------------------------
 
CL              PLUME CONVECTION
 
      CALL PLUME (YTFLX, YQFLX, YUFLX, YVFLX,  YT, YQ, YU, YV,
     *            YLWC_C, YLWC_A, YLWC_S,
     *            FPLUMT, FPLUMQ, FPLUMU, FPLUMV,
     *            TRACER, TRAFLUX, YPM, YPTHIC, YPMCAP,
     *            WPLUME0, WPLUME1, PLUMEPREC,
     *            HEIMIN, HEIMAX, DTIME, JJ
#ifdef isotrac
     *           ,yqflx_o, yqflx_h, yq_o, yq_h,
     *            ylwc_c_o, ylwc_c_h, ylwc_a_o, ylwc_a_h,
     *            ylwc_s_o, ylwc_s_h,
     *            fplumq_o, fplumq_h, plumepr_o, plumepr_h,
     *            za_i, zh_i, rtrfv, rtrq, rtrfq, plumepr_t
#endif
     *           )
c
c-----
c     Turn off convection, just give surface fluxes to lowest layer
c     do 100 ji=1,nlon
c       fplumt(ji,jj,nlev) = fplumt(ji,jj,nlev) + ytflx(ji,jj)*dtime
c       fplumq(ji,jj,nlev) = fplumq(ji,jj,nlev) + yqflx(ji,jj)*dtime
c       fplumu(ji,jj,nlev) = fplumu(ji,jj,nlev) + yuflx(ji,jj)*dtime
c       fplumv(ji,jj,nlev) = fplumv(ji,jj,nlev) + yvflx(ji,jj)*dtime
c 100 continue
c-----

C
CL              LOCAL STRATIFORM CONDENSATION
C
        CALL STABCO (YT, YQ, YPM, YPTHIC, FPLUMT, FPLUMQ, YLWC_S,YCLD,JJ
#ifdef isotrac
     *              ,yq_o, yq_h, fplumq_o, fplumq_h, ylwc_s_o, ylwc_s_h,
     *               rtrq, rtrfq
#endif
     *              )
C
CL              LOCAL EVAPORATION OF CLOUDS, CONVERSION TO PRECIP,
CL              AND REEVAP OF FALLING PRECIP (SEPARATELY FOR
CL              CONV/ANVIL/STRAT DUE TO DIFFERENT PARAM VALUES AND FOR
CL              DIAGNOSTICS). GIVE ANVIL-CIRRUS PRECIP TO FRAINC,FSNOWC.
 
c       convective:
        CALL REEVAP (YT, YQ, YPM, YPTHIC, YLWC_C, YCLD, PLUMEPREC,
     *               FPLUMT, FPLUMQ, FRAINC, FSNOWC, PRECFLX,
     *               DMPRECC, JJ, 1
#ifdef isotrac
     *              ,yq_o, yq_h, ylwc_c_o, ylwc_c_h,
     *               plumepr_o, plumepr_h, fplumq_o, fplumq_h,
     *               frainc_o, frainc_h, fsnowc_o, fsnowc_h, za_i, zh_i,
     *               rtrq, plumepr_t, rtrfq, rtrprc, rtrpsc
#endif
     *              )
 
c       anvil:
        CALL REEVAP (YT, YQ, YPM, YPTHIC, YLWC_A, YCLD, ZDUM,
     *               FPLUMT, FPLUMQ, FRAINC, FSNOWC, PRECFLX,
     *               DMPRECA, JJ, 2
#ifdef isotrac
     *              ,yq_o, yq_h, ylwc_a_o, ylwc_a_h,
     *               zdum_o, zdum_h, fplumq_o, fplumq_h,
     *               frainc_o, frainc_h, fsnowc_o, fsnowc_h, duma, dumh,
     *               rtrq, zdum_t, rtrfq, rtrprc, rtrpsc
#endif
     *              )
 
c       stratus:
        CALL REEVAP (YT, YQ, YPM, YPTHIC, YLWC_S, YCLD, ZDUM,
     *               FPLUMT, FPLUMQ, FRAINS, FSNOWS, PRECFLX,
     *               DMPRECS, JJ, 3
#ifdef isotrac
     *              ,yq_o, yq_h, ylwc_s_o, ylwc_s_h,
     *               zdum_o, zdum_h, fplumq_o, fplumq_h,
     *               frains_o, frains_h, fsnows_o, fsnows_h, duma, dumh,
     *               rtrq, zdum_t, rtrfq, rtrprs, rtrpss
#endif
     *              )
 
      RETURN
      END
C
C ***************************** ESTABV ********************************
C
      SUBROUTINE ESTABV (YES, YQSLEV, YTLEV, YPRES, NP)
C----------------------------------------------------------------------
C
CL            CALCULATES SATURATION VAPOR PRESSURE (YES) AND SATURATION
CL            SPECIFIC HUMIDITY (YQSLEV) FOR GIVEN TEMPERATURE (YTLEV)
CL            AND MID-LAYER PRESSURE (YPRES). DOES CALCULATIONS FOR ONE
CL            STRIP OF NP (<= NLON) POINTS AT A TIME.
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMEST.h>
C----------------------------------------------------------------------
      DIMENSION YES(NP), YQSLEV(NP), YTLEV(NP), YPRES(NP)
C----------------------------------------------------------------------
      DIMENSION IT1(NLON), IT2(NLON), ZTMP1(NLON), ZTMP2(NLON)
C----------------------------------------------------------------------
 
      DO 10 JI=1,NP
        YQSLEV(JI) = MAX (TMIN, MIN (TMAX, YTLEV(JI)))
        ZTMP2(JI) = TMIN + INT(YQSLEV(JI)-TMIN)
        IT1(JI) = YQSLEV(JI) - TMIN + 1
        IT2(JI) = IT1(JI) + 1
   10 CONTINUE
 
c     CALL GATHER (NP, YES,   TABLE, IT1)
c     CALL GATHER (NP, ZTMP1, TABLE, IT2)
      DO 15 JI=1,NP
        YES(JI)   = TABLE(IT1(JI))
        ZTMP1(JI) = TABLE(IT2(JI))
   15 CONTINUE
 
      DO 20 JI=1,NP
        YES(JI) = (ZTMP2(JI) + 1. - YQSLEV(JI)) * YES(JI)
     *          + (YQSLEV(JI) - ZTMP2(JI)) * ZTMP1(JI)
        YQSLEV(JI) = .622*YES(JI)
     *             / MAX ( YPRES(JI)-(1.-.622)*YES(JI), .622*YES(JI) )
   20 CONTINUE
 
      RETURN
      END
C
C ***************************** stabco *********************************
C
      subroutine stabco (
     *             yt, yq, ypm, ypthic, fplumt, fplumq, ylwc_s, ycld, jj
#ifdef isotrac
     *            ,yq_o, yq_h, fplumq_o, fplumq_h, ylwc_s_o, ylwc_s_h,
     *             rtrq, rtrfq
#endif
     *                  )
 
c     Does stable condensation, iterating niter times, for one latitude
c     strip.
 
c     yt       = temperature field (modified)
c     yq       = specific humidity (modified)
c     ypm      = mid-layer pressure (supplied)
c     ypthic   = layer pressure thickness (supplied)
c     fplumt,q = reservoirs of "flux convergence*time" (3-D) (modified)
c     ylwc_s   = stratiform liquid water content (3-D) (modified)
c     ycld     = 3-D cloud fraction (from cldcmp) (supplied)
c     jj       = latitude index (supplied)
 
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMBUD.h>
C ---------------------------------------------------------------------
#include <COMCON.h>
C ---------------------------------------------------------------------
#include <COMFLG.h>
C ---------------------------------------------------------------------
#include <COMLUN.h>
C ---------------------------------------------------------------------
#include <COMMAP.h>
C ---------------------------------------------------------------------
#include <COMTIM.h>
C ---------------------------------------------------------------------
      dimension
     *  yt(nlon,norec,nlev),      yq(nlon,norec,nlev),
     *  ypm(nlon,norec,nlev),     ypthic(nlon,norec,nlev),
     *  fplumt(nlon,norec,nlev),  fplumq(nlon,norec,nlev),
     *  ylwc_s(nlon,norec,nlev),  ycld(nlon,norec,nlev)
#ifdef isotrac
      dimension
     *  yq_o(nlon,norec,nlev),    yq_h(nlon,norec,nlev),
     *  fplumq_o(nlon,norec,nlev),fplumq_h(nlon,norec,nlev),
     *  ylwc_s_o(nlon,norec,nlev),ylwc_s_h(nlon,norec,nlev),
     *  rtrq(nlon,norec,nlev,nrtr),rtrfq(nlon,norec,nlev,nrtr)
#endif
 
      dimension
     *  zt(nlon),      zq(nlon),
     *  delq(nlon),    za(nlon),
     *  zsvp(nlon),    zsvq(nlon)
#ifdef isotrac
      dimension
     *  zq_o(nlon),    zq_h(nlon),
     *  delq_o(nlon),  delq_h(nlon)
#endif
      parameter (niter=1)
C ---------------------------------------------------------------------
 
c        Overall loop over level
 
c=======================
      do 1100 jk=1,nlev
c=======================
 
c          Copy t,q to temporary arrays for this longitude circle
c          (so don't change coml30 arrays for time n)
 
        call scopy (nlon, yt(1,jj,jk), 1, zt, 1)
        call scopy (nlon, yq(1,jj,jk), 1, zq, 1)
#ifdef isotrac
        call scopy (nlon, yq_o(1,jj,jk), 1, zq_o, 1)
        call scopy (nlon, yq_h(1,jj,jk), 1, zq_h, 1)
#endif
 
c          Iterate for this longitude circle
 
c---------------------------
        do 2000 iter=1,niter
c---------------------------
 
          call estabv (zsvp, zsvq, zt, ypm(1,jj,jk), nlon)
 
          do 200 ji=1,nlon
 
c              za is latent heat of condensation, determined by ambient
c              temperature (ie, rain or snow), plus correction
c              for specific heats to be consistent with lsx (imagining
c              all phase changes occur at tmelt)
 
c           za(ji) = cvmgt ( latvap + cpwv*(zt(ji)-tmelt)
c    *                              - ch2o*(zt(ji)-tmelt),
c    *                       latsub + cpwv*(zt(ji)-tmelt)
c    *                              - cice*(zt(ji)-tmelt),
c    *                       zt(ji).ge.tmelt )
c           For now, all cloud water is liquid with no specific heat
            za(ji) = latvap + cpwv*(zt(ji)-tmelt)
 
c              Test if exceed zqmax, which ranges from alw*saturation
c              for zstrat=0 and 100%saturation for zstrat=1, where
c              zstrat = existing (all) cloud fractional area (as
c              in cldcmp). This assumes uniformly distributed specific
c              humidities in non-cloud areas with range of +/- 1.-alw,
c              and saturation in cloud areas. Condense amount that
c              yields zqmax, allowing fopr linearized change in svp with
c              temperature. Store condensed amount in delq.
c              (alw is a namelist param).
c----
c           zstrat = 0.
            zstrat = min (1., max (0., ylwc_s(ji,jj,jk)/clwfac_s))
c           zstrat = ycld(ji,jj,jk)  ! alternate
            zqmax = zsvq(ji) * ( zstrat  + (1.-zstrat)*alw )
 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if (zq(ji).gt.zqmax) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              dz = (za(ji)/rh2o) * zqmax / (zt(ji)**2)
              cpz = cpair*(1.+cpvir*zq(ji))
              delq(ji) = (zq(ji) - zqmax) / (1.+dz*za(ji)/cpz)
#ifdef isotrac
              zqt_i  = zq(ji)
              zqt_o  = zq_o(ji)
              zqt_h  = zq_h(ji)
              ztt_i  = zt(ji)
#endif
              zq(ji) = zq(ji) - delq(ji)
              zt(ji) = zt(ji) + delq(ji)*za(ji)
     *                                   / (cpair*(1.+cpvir*zq(ji)))
#ifdef isotrac

c             If the temperature is above the transition temperature
c             for isotope's in stratus, equilibrate first the liquid
c             with the vapour, then call prm_i to calculate the new
c             isotopic compositions of the liquid and the vapour.
c             If not, call prm_i to calculate the isotopic composition
c             of the newly formed solid independently of the solid
c             already present and deduce the new isotopic composition
c             of the solid and the vapour.

              zlwt_i  = ylwc_s(ji,jj,jk)
              zlwt_o  = ylwc_s_o(ji,jj,jk)
              zlwt_h  = ylwc_s_h(ji,jj,jk)
              if ((ztt_i.gt.t_s_i).and.(zt(ji).ge.t_s_i)) then
c               liquid condensate:
                call equil_i (1., ztt_i, ypm(ji,jj,jk), zqt_i,
     *                      zq_o(ji), zq_h(ji), zlwt_i, zlwt_o, zlwt_h)
                call prm_i (s_s_i, t_s_i, ztt_i, zt(ji),
     *                   zqt_i,  zq(ji),          zq_o(ji), zq_h(ji),
     *                   zlwt_i, zlwt_i+delq(ji), zlwt_o,   zlwt_h)
              else
c               solid condensate:
                call prm_i (s_s_i, t_s_i, ztt_i, zt(ji),
     *                          zqt_i,  zq(ji),   zq_o(ji), zq_h(ji),
     *                          0.,     delq(ji), zlwt_o,   zlwt_h)
              endif
              delq_o(ji) = zqt_o-zq_o(ji)
              delq_h(ji) = zqt_h-zq_h(ji)
#endif
c>>>>>>>>>>>>>>>
            else
c>>>>>>>>>>>>>>>
              delq(ji) = 0.
#ifdef isotrac
              delq_o(ji)= 0.
              delq_h(ji)= 0.
#endif
c>>>>>>>>>>>>>>>>
            endif
c>>>>>>>>>>>>>>>>
c----
c           "Soft" reduction above alw*zsvq:
c           zx = zq(ji)-alw*zsvq(ji)
c           zd = min (0.3*zx + .01*zsvq(ji), zx)
c           delq(ji) = max (zd, zq(ji)-zsvq(ji), 0.)
c           zq(ji) = zq(ji) - delq(ji)
c           zt(ji) = zt(ji) + delq(ji)*za(ji) /(cpair*(1.+cpvir*zq(ji)))
c----
 
c              Accumulate reservoirs fplum[t,q] (flux convergences*time
c              per layer). Use factor of 0.5 since reservoirs are not
c              leapfrog variables (see comments in reserv). Also
c              increment stratiform liquid water content (Kg/Kg) and
c              increment global budget terms (d[m,h]clouf)
 
            zz = delq(ji) * ypthic(ji,jj,jk)/gravit
            fplumt(ji,jj,jk) = fplumt(ji,jj,jk) + zz*za(ji)*0.5
            fplumq(ji,jj,jk) = fplumq(ji,jj,jk) - zz*0.5
            ylwc_s(ji,jj,jk) = ylwc_s(ji,jj,jk) + delq(ji)*0.5
#ifdef isotrac
            fplumq_o(ji,jj,jk) = fplumq_o(ji,jj,jk) -
     *                           0.5*delq_o(ji)*ypthic(ji,jj,jk)/gravit
            fplumq_h(ji,jj,jk) = fplumq_h(ji,jj,jk) -
     *                           0.5*delq_h(ji)*ypthic(ji,jj,jk)/gravit
            ylwc_s_o(ji,jj,jk) = ylwc_s_o(ji,jj,jk) + delq_o(ji)*0.5
            ylwc_s_h(ji,jj,jk) = ylwc_s_h(ji,jj,jk) + delq_h(ji)*0.5
#endif
#ifdef isotracreg
            if (yq(ji,jj,jk).ne.0.) then
              do it=1,nrtr
                rtrfq(ji,jj,jk,it) = rtrfq(ji,jj,jk,it) - 0.5*delq(ji)*
     *            rtrq(ji,jj,jk,it)/yq(ji,jj,jk)*ypthic(ji,jj,jk)/gravit
              enddo
            endif
#endif
            dmclouf(jj) = dmclouf(jj) + zz       *0.5/nlon
            dhcloud(jj) = dhcloud(jj) + zz*za(ji)*0.5/nlon
  200     continue
 
c---------------
 2000   continue
c---------------
 
c=============
 1100 continue
c=============
 
      return
      end
C
C ***************************** plume ********************************
C
      subroutine plume ( ytflx, yqflx, yuflx, yvflx, yt, yq, yu, yv,
     *                   ylwc_c, ylwc_a, ylwc_s,
     *                   fplumt, fplumq, fplumu, fplumv,
     *                   tracer, traflux, ypm, ypthic, ypmcap,
     *                   wplume0, wplume1, plumeprec,
     *                   heimin, heimax, dtplume, jj
#ifdef isotrac
     *                  ,yqflx_o, yqflx_h, yq_o, yq_h,
     *                   ylwc_c_o, ylwc_c_h, ylwc_a_o, ylwc_a_h,
     *                   ylwc_s_o, ylwc_s_h,
     *                   fplumq_o, fplumq_h, plumepr_o, plumepr_h,
     *                   za_i, zh_i, rtrfv, rtrq, rtrfq, plumepr_t
#endif
     *                  )
 
c     Does convection (both pbl and free), by solving a subgrid
c     plume model. For pbl convection, the plume is initiated
c     at the center of the lowest layer (ie, top of constant-flux layer)
c     with w,t,etc values scaled depending on the surface heat/momentum
c     fluxes based on constant-flux scaling, and is integrated
c     vertically until it stops. Doe sone latitude strip at a time.
 
c     For free plumes, there is an option (commented by "csing" or
c     "cmult") to do just one vertical integration from the bottom to
c     the top (csing) or to do multiplte integrations stating from
c     each layer in turn (cmult). In both cases, a free plume is
c     initiated at the center of a model layer if it can rise to
c     (or beyond) the center of the next layer. Plumes are integrated
c     vertically upwards from layer midpoint to midpoint until the
c     plume stops.
 
c     Condensation can occur within the plumes, forming wthin-plume
c     liquid cloud water, which can precipitate out if it exceeds a
c     maximum value.
 
c     The plume equations explicitly solve for plume w,theta,q,u,v where
c     all except w are perturbations from the ambient quantities. The
c     plume vertical finite-difference step is upwards from layer
c     midpoint to layer midpoint. These quantities, along with the plume
c     condensation amounts, then imply the large-scale fluxes at layer
c     midpoints and the large-scale precip and latent heat release.
c     These fluxes are accumulated over pbl and free
c     calculations, and then "shifted" downwards to the next lowest
c     layer interface to calculate large-scale convergences.
c     The large-scale flux convergences for heat, water vapor and
c     momentum are actually put in reservoirs fplum* and released
c     (see reserv, called from vdif) with a time scale of a few hours.
c     Liquid water cloud amounts and tracers are updated here with
c     no reservoir.
 
c     Lines beginning with "cloc" are to convert to standalone mode
c     for use with ./Plume/plume.f
 
c     supplied:
c     ytflx   = upward sensible heat flux from surface (W/m2)
c     yqflx   = upward h2o mass flux from surface (Kg/m2/s)
c     yuflx   = upward u-momentum flux from surface (N/m2)
c     yvflx   = upward v-momentum flux from surface (N/m2)
c     yt      = agcm temperatures (deg K)
c     yq      = agcm specific humidity (Kg/Kg)
c     yu      = agcm eastward  velocity (m/s)
c     yv      = agcm northward velocity  (m/s)
c     ypm     = agcm mid-layer pressures (N/m2)
c     ypthic  = agcm layer pressure-thicknesses (N/m2)
c     ypmcap  = (ypm/surface pressure)**cappa (N/m2)
c     dtplume = plume model timestep (s)
c     jj      = latitude index
 
c     modified:
c     ylwc_c  = convective liquid water content (3-D)
c     ylwc_a  = anvil-cirrus liquid water content (3-D)
c     ylwc_s  = stratiform   liquid water content (3-D)
c     fplum*  = reservoirs of "flux convergence*time" in each agcm layer
c     tracer  = agcm tracer fields. ylwc_[c,a,s] are equivalenced to
c               tracers #1-3, and are changed explicitly below.
c               fplum* are equivalenced to tracers #4-7
c               (ntraca=7, set in commons). Only tracers
c               #ntraca+1 to #ntrace are mixed as "tracers" below.
c     traflux = surface fluxes for tracers ntraca+1 to ntrace
c     wplume0,1 = zonal mean plume vertical velocities (for diagnostics)
c     plumeprec = within-plume precipitation formation rate
c                 (3-D,kg/m2/s, passed from phys, for reevap)
c     heimin    = 03:00 (was daily min) pbl plume heights (for history)
c     heimax    = 15:00 (was daily max) pbl plume heights (for history)
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMBUD.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMEST.h>
c ---------------------------------------------------------------------
#include <COMFLG.h>
c ---------------------------------------------------------------------
#include <COMINI.h>
c ---------------------------------------------------------------------
#include <COMLUN.h>
c ---------------------------------------------------------------------
#include <COMMAP.h>
c ---------------------------------------------------------------------
#include <COMTIM.h>
c-----------------------------------------------------------------------
c#include <complume.h>                                             !cloc
c-----------------------------------------------------------------------
      dimension
     *  ytflx(nlon,norec),        yqflx(nlon,norec),
     *  yuflx(nlon,norec),        yvflx(nlon,norec),
 
     *  yt(nlon,norec,nlev),      yq(nlon,norec,nlev),
     *  yu(nlon,norec,nlev),      yv(nlon,norec,nlev),
     *  ylwc_c(nlon,norec,nlev),  ylwc_a(nlon,norec,nlev),
     *  ylwc_s(nlon,norec,nlev),
     *  fplumt(nlon,norec,nlev),  fplumq(nlon,norec,nlev),
     *  fplumu(nlon,norec,nlev),  fplumv(nlon,norec,nlev),
     *  tracer(nlon,norec,nlev,ntrace), traflux(nlon,norec,ntracb),
 
     *  ypm(nlon,norec,nlev),
     *  ypthic(nlon,norec,nlev),  ypmcap(nlon,norec,nlev),
 
     *  wplume0(norec,nlev),      wplume1(norec,nlev),
     *  plumeprec(nlon,norec,nlev),
     *  heimin(nlon,norec),       heimax(nlon,norec)
#ifdef isotrac
      dimension
     *  yqflx_o(nlon,norec),      yqflx_h(nlon,norec),
     *  yq_o(nlon,norec,nlev),    yq_h(nlon,norec,nlev),
     *  ylwc_c_o(nlon,norec,nlev),ylwc_c_h(nlon,norec,nlev),
     *  ylwc_a_o(nlon,norec,nlev),ylwc_a_h(nlon,norec,nlev),
     *  ylwc_s_o(nlon,norec,nlev),ylwc_s_h(nlon,norec,nlev),
     *  fplumq_o(nlon,norec,nlev),fplumq_h(nlon,norec,nlev),
     *  plumepr_o(nlon,norec,nlev),plumepr_h(nlon,norec,nlev),
     *  rtrfv(nlon,norec,nrtr),   rtrq(nlon,norec,nlev,nrtr),
     *  rtrfq(nlon,norec,nlev,nrtr),plumepr_t(nlon,norec,nlev,nrtr)
#endif
 
      dimension
     *  yh(nlon,nlev),         yrho(nlon,nlev),
     *  ylwc_t(nlon,nlev),
 
     *  dzu(nlon,nlev),        dzl(nlon,nlev),
     *  dztot(nlon,nlev),      dz(nlon,nlev),
     *  dhdz(nlon,nlev),       dqdz(nlon,nlev),
     *  dudz(nlon,nlev),       dvdz(nlon,nlev),
     *  dlwcdz(nlon,nlev),     dadz(nlon,nlev,ntrace),
     *  zelev(nlon,nlev)
#ifdef isotrac
      dimension
     *  ylwc_t_o(nlon,nlev),   ylwc_t_h(nlon,nlev),
     *  dqdz_o(nlon,nlev),     dqdz_h(nlon,nlev),
     *  dlwcdz_o(nlon,nlev),   dlwcdz_h(nlon,nlev),
     *  dqdz_t(nlon,nlev,nrtr)
#endif
 
      dimension
     *  zh(nlon,nlev),        zq(nlon,nlev),
     *  zu(nlon,nlev),        zv(nlon,nlev),
     *  zw(nlon,nlev),        zws(nlon,nlev),
     *  zt(nlon,nlev),        zrho(nlon,nlev),
     *  zlw(nlon,nlev),       za(nlon,nlev,ntrace),
     *  zcond(nlon,nlev),     zlath(nlon,nlev),
     *  zprec(nlon,nlev),
 
     *  zconda(nlon,0:nlev),  zheata(nlon,0:nlev),
     *  zpreca(nlon,0:nlev),  zmfa(nlon,0:nlev),
     *  tplu(nlon,0:nlev),    qplu(nlon,0:nlev),
     *  uplu(nlon,0:nlev),    vplu(nlon,0:nlev),
     *  cpluc(nlon,0:nlev),   cplua(nlon,0:nlev),
     *  cplus(nlon,0:nlev),   aplu(nlon,0:nlev,ntrace)

#ifdef isotrac
      dimension
     *  zq_o(nlon,nlev),      zq_h(nlon,nlev),
     *  zlw_o(nlon,nlev),     zlw_h(nlon,nlev),
     *  zcond_o(nlon,nlev),   zcond_h(nlon,nlev),
     *  zprec_o(nlon,nlev),   zprec_h(nlon,nlev),
     *  zconda_o(nlon,0:nlev),zconda_h(nlon,0:nlev),
     *  zpreca_o(nlon,0:nlev),zpreca_h(nlon,0:nlev),
     *  qplu_o(nlon,0:nlev),  qplu_h(nlon,0:nlev),
     *  cpluc_o(nlon,0:nlev), cpluc_h(nlon,0:nlev),
     *  cplua_o(nlon,0:nlev), cplua_h(nlon,0:nlev),
     *  cplus_o(nlon,0:nlev), cplus_h(nlon,0:nlev),
     *  yl_c_i(nlon,nlev),    yl_a_i(nlon,nlev),    yl_s_i(nlon,nlev),
     *  zapbl(nlon,nlev),     zrpbl(nlon,nlev),
     *  za_i(nlon,nlev),      zh_i(nlon,nlev),

     *  zq_t(nlon,nlev,nrtr),       zcond_t(nlon,nlev,nrtr),
     *  zprec_t(nlon,nlev,nrtr),    zconda_t(nlon,0:nlev,nrtr),
     *  zpreca_t(nlon,0:nlev,nrtr), qplu_t(nlon,0:nlev,nrtr),
     *  delq_t(nrtr),               delw_t(nrtr)
#endif

 
c     tracers 1-3 are used for liquid water contents (mixed by plume),
c     and tracers 4-7 are used for flux reservoirs (not mixed by plume).
      equivalence
     *  (dlwcdz,  dadz(1,1,1)),
     *  (zlw,     za(1,1,1)),
     *  (cpluc,   aplu(1,0,1)),
     *  (cplua,   aplu(1,0,2)),
     *  (cplus,   aplu(1,0,3))
 
      dimension
     *  zarea(nlon,nlev),      zrad(nlon,nlev),
     *  zent(nlon),            zdrag(nlon),
     *  zbuoy(nlon),
 
     *  zcfl(nlon),            zmf(nlon,nlev),
     *  ztotc1(nlon),          ztotc2(nlon),
     *  ztota1(nlon),          ztota2(nlon),
     *  ztots1(nlon),          ztots2(nlon),
     *  zlwcold(nlon,nlev)
 
      dimension
     *  heipbl(nlon),
     *  laypbl(nlon),           numpbl(nlev),
     *  laytop(nlon),           numfret(nlev),
     *  laybot(nlon),           numfreb(nlev),
     *  wfre1(nlon,nlev),       wfre2(nlon,nlev)
 
      logical ifpbl,            ifanv(nlon)
      character*3 chaplu(nlev)
 
c     Diagnostic only:
      dimension qrat(nlev), nqrat(nlev), qneg(nlev), nqneg(nlev)
      parameter (nprinz=1)
cloc  parameter (nprinz=2)                                         !cloc
      dimension iuprinz(nprinz)
 
      parameter (pi = 3.14159265358979)
 
c-----------------------------------------------------------------------
#include <SFEST.h>
c-----------------------------------------------------------------------
cloc  include 'comsatplume'                                        !cloc
c-----------------------------------------------------------------------
c     statement function for plume radius vs. elevation (not used):
c     plumerad (z) = plumerad1 + min (0.5*z, plumerad2)
c ---------------------------------------------------------------------
 
c        Initialize local diagnostics
 
      zptot = 0.
      zetot = 0.
      zqtot = 0.
      do 10 jk=1,nlev
        numpbl(jk) = 0
        numfret(jk) = 0
        numfreb(jk) = 0
   10 continue
 
c        Find longitudes at or spanning 03:00 and 15:00 local time for
c        pbl height history. If timesteps_per_day <= nlon,
c        then ji03a=ji03b and ji15a=ji15b. Else must set more than
c        one longitude per timestep, and ji03a<ji03b, ji15a<ji15b.
c        ji03a, etc, may be < 1 or > nlon...dealt with when used below.
 
      izmin03 = 1000000
      izmin15 = 1000000
      do 20 ji=1,nlon
        izsec = mod ( nint(nstep*dtime + 86400.*along(ji)/(2.*pi)),
     *                86400 )
c       Subtract 1 sec to break "ties" if 03:00,15:00 fall on box edge
        iz03 = iabs(izsec - 1 - 03*3600)
        iz15 = iabs(izsec - 1 - 15*3600)
        if (iz03.lt.izmin03) then
          ji03 = ji
          izmin03 = iz03
        endif
        if (iz15.lt.izmin15) then
          ji15 = ji
          izmin15 = iz15
        endif
  20  continue
 
      izspan = (nlon-1)/nint(86400./dtime) + 1
      ji03a = ji03  - izspan/2
      ji03b = ji03a + izspan -1
      ji15a = ji15  - izspan/2
      ji15b = ji15a + izspan -1
 
 
c         Calculate ambient pot.temp, density and total cloud liquid
c         water at layer midpoints
 
      do 100 jk=1,nlev
        do 102 ji=1,nlon
          yh(ji,jk)     = yt(ji,jj,jk) / ypmcap(ji,jj,jk)
          yrho(ji,jk)   = ypm(ji,jj,jk)
     *                  / (rair*(1.+zvir*yq(ji,jj,jk))*yt(ji,jj,jk))
          ylwc_t(ji,jk) = ylwc_c(ji,jj,jk) + ylwc_a(ji,jj,jk)
     *                  + ylwc_s(ji,jj,jk)
#ifdef isotrac
          ylwc_t_o(ji,jk) = ylwc_c_o(ji,jj,jk) + ylwc_a_o(ji,jj,jk)
     *                  + ylwc_s_o(ji,jj,jk)
          ylwc_t_h(ji,jk) = ylwc_c_h(ji,jj,jk) + ylwc_a_h(ji,jj,jk)
     *                  + ylwc_s_h(ji,jj,jk)
#endif
  102   continue
  100 continue
 
c        Compute ambient dz, d*/dz for full layers. For bottom layer,
c        just use top half of layer.
 
      do 120 jk=1,nlev
        do 122 ji=1,nlon
          dztot(ji,jk) = ypthic(ji,jj,jk) / (gravit*yrho(ji,jk))
          dzu(ji,jk) = dztot(ji,jk) * (sig(jk)-sigkmh(jk))/dsigma(jk)
          dzl(ji,jk) = dztot(ji,jk) - dzu(ji,jk)
  122   continue
  120 continue
 
      do 124 jk=1,nlev-1
        do 126 ji=1,nlon
          dz(ji,jk) = dzl(ji,jk) + dzu(ji,jk+1)
          dhdz(ji,jk)   = (yh(ji,jk)     - yh(ji,jk+1)    ) / dz(ji,jk)
          dqdz(ji,jk)   = (yq(ji,jj,jk)  - yq(ji,jj,jk+1) ) / dz(ji,jk)
          dudz(ji,jk)   = (yu(ji,jj,jk)  - yu(ji,jj,jk+1) ) / dz(ji,jk)
          dvdz(ji,jk)   = (yv(ji,jj,jk)  - yv(ji,jj,jk+1) ) / dz(ji,jk)
          dlwcdz(ji,jk) = (ylwc_t(ji,jk) - ylwc_t(ji,jk+1)) / dz(ji,jk)
#ifdef isotrac
          dqdz_o(ji,jk) = (yq_o(ji,jj,jk)- yq_o(ji,jj,jk+1))/ dz(ji,jk)
          dqdz_h(ji,jk) = (yq_h(ji,jj,jk)- yq_h(ji,jj,jk+1))/ dz(ji,jk)
          dlwcdz_o(ji,jk)=(ylwc_t_o(ji,jk)-ylwc_t_o(ji,jk+1))/dz(ji,jk)
          dlwcdz_h(ji,jk)=(ylwc_t_h(ji,jk)-ylwc_t_h(ji,jk+1))/dz(ji,jk)
#endif
#ifdef isotracreg
          do it=1,nrtr
            dqdz_t(ji,jk,it) =
     *              (rtrq(ji,jj,jk,it)- rtrq(ji,jj,jk+1,it))/ dz(ji,jk)
          enddo
#endif
  126   continue
 
        if (ntraca.lt.ntrace) then
          do 1260 n=ntraca+1,ntrace
            do 1262 ji=1,nlon
             dadz(ji,jk,n) = (tracer(ji,jj,jk,n)-tracer(ji,jj,jk+1,n))
     *                       /dz(ji,jk)
 1262       continue
 1260     continue
        endif
  124 continue
 
c        Compute elevations of layer midpoints (for plume radii calc)
 
      do 130 ji=1,nlon
        zelev(ji,nlev) = dzl(ji,nlev)
  130 continue
      do 132 jk=nlev-1,1,-1
        do 134 ji=1,nlon
          zelev(ji,jk) = zelev(ji,jk+1) + dz(ji,jk)
  134   continue
  132 continue
 
c        Zero plume fluxes accumulated over loop 500
 
      call zero (tplu,  nlon*nlevp)
      call zero (qplu,  nlon*nlevp)
      call zero (uplu,  nlon*nlevp)
      call zero (vplu,  nlon*nlevp)
      call zero (cpluc, nlon*nlevp)
      call zero (cplua, nlon*nlevp)
      call zero (cplus, nlon*nlevp)
      if (ntraca.lt.ntrace)
     *  call zero (aplu(1,0,ntraca+1), nlon*nlevp*(ntrace-ntraca))
      call zero (zconda, nlon*nlevp)
      call zero (zheata, nlon*nlevp)
      call zero (zpreca, nlon*nlevp)
      call zero (zmfa,   nlon*nlevp)
#ifdef isotrac
      call zero (qplu_o,  nlon*nlevp)
      call zero (qplu_h,  nlon*nlevp)
      call zero (cpluc_o, nlon*nlevp)
      call zero (cpluc_h, nlon*nlevp)
      call zero (cplua_o, nlon*nlevp)
      call zero (cplua_h, nlon*nlevp)
      call zero (cplus_o, nlon*nlevp)
      call zero (cplus_h, nlon*nlevp)
      call zero (zconda_o, nlon*nlevp)
      call zero (zconda_h, nlon*nlevp)
      call zero (zpreca_o, nlon*nlevp)
      call zero (zpreca_h, nlon*nlevp)
      call zero (zapbl, nlon*nlev)
      call zero (zrpbl, nlon*nlev)
      call zero (za_i, nlon*nlev)
      call zero (zh_i, nlon*nlev)
#endif
#ifdef isotracreg
      call zero (qplu_t,   nlon*nlevp*nrtr)
      call zero (zconda_t, nlon*nlevp*nrtr)
      call zero (zpreca_t, nlon*nlevp*nrtr)
#endif

c        Loop twice, for "pbl" plumes (ifpbl=.true.) and for "free"
c        plumes (ifpbl=.false., start at second lowest layer)
 
c*********************
      do 500 itype=0,1
c*********************
 
      if (itype.eq.0) then
        ifpbl = .true.
        jkinita = nlev
        jkinitb = nlev
      else
        ifpbl = .false.
        jkinita = nlev
        jkinitb = nlev                                            !csing
cmult   jkinitb = 2                                               !cmult
      endif
 
c         Initialize diagnostics for setting of plume velocities
c         (wfre*,wplume*) and indices of bottom,top layers (lay*, num*)
 
      call zero (wfre1, nlon*nlev)
      call zero (wfre2, nlon*nlev)
 
      do 140 ji=1,nlon
        heipbl(ji) = dzl(ji,nlev)
        laypbl(ji) = nlevp
        laytop(ji) = nlevp
        laybot(ji) = 0
  140 continue
 
c************************************
      do 510 jkinit = jkinita,jkinitb
c************************************
 
c        Assume no new plumes initiated above 150 mb
 
      if (sig(jkinit).lt.0.150) goto 510
 
c        Initialize all plume variables to zero or ambient, and
c        initialize diagnostic plume top/bottom locations
 
      call zero (zw,   nlon*nlev)
      call zero (zws,  nlon*nlev)
      call zero (zt,   nlon*nlev)
      call zero (zh,   nlon*nlev)
      call zero (zq,   nlon*nlev)
      call zero (zu,   nlon*nlev)
      call zero (zv,   nlon*nlev)
      call zero (zlw,  nlon*nlev)
      call scopy (nlon*nlev, yrho, 1, zrho, 1)
      if (ntraca.lt.ntrace)
     *  call zero (za(1,1,ntraca+1), nlon*nlev*(ntrace-ntraca))
      call zero (zarea,nlon*nlev)
      call zero (zcond, nlon*nlev)
      call zero (zlath, nlon*nlev)
      call zero (zprec, nlon*nlev)
#ifdef isotrac
      call zero (zq_o,   nlon*nlev)
      call zero (zq_h,   nlon*nlev)
      call zero (zlw_o,  nlon*nlev)
      call zero (zlw_h,  nlon*nlev)
      call zero (zcond_o, nlon*nlev)
      call zero (zcond_h, nlon*nlev)
      call zero (zprec_o, nlon*nlev)
      call zero (zprec_h, nlon*nlev)
#endif
#ifdef isotracreg
      call zero (zq_t,    nlon*nlev*nrtr)
      call zero (zcond_t, nlon*nlev*nrtr)
      call zero (zprec_t, nlon*nlev*nrtr)
#endif
 
c       For "pbl" plumes, initialize plume conditions
c       at midpoint of bottom layer based on surface fluxes from
c       constant-flux layer. If net buoyancy flux
c       (ytflx/(cp*yt)+.622*yqflx) is < 0, leave initial plume z* = zero
 
c--------------------
      if (ifpbl) then
c--------------------
 
        do 150 ji=1,nlon
 
cloc      zarea(ji,nlev) = xare1                                   !cloc
cloc      zrad (ji,nlev) = xrad1                                   !cloc
 
          zbuoy(ji) = ytflx(ji,jj)/(cpair*yt(ji,jj,nlev))
     *              + .622*yqflx(ji,jj)
 
c         Smaller plumearea1 over land than ocean/ice
          zarea(ji,nlev) = plumearea1_l*    yoroland(ji,jj)     ! .005
     *                   + plumearea1_o*(1.-yoroland(ji,jj))    ! .2
cc        zz = max (0., min (1., (zbuoy(ji)-100.)/200. ))
cc        zarea(ji,nlev) = plumearea1_l*zz + plumearea1_o*(1.-zz)
          zrad(ji,nlev) = plumerad1                            ! 100.
 
          if (zbuoy(ji).gt. 0.) then
 
            ztau = sqrt (yuflx(ji,jj)**2 + yvflx(ji,jj)**2)
            ztauv = max (.01, sqrt(ztau/yrho(ji,nlev)) )
            zconv = ( (max(ytflx(ji,jj),0.)/(yrho(ji,nlev)*cpair))
     *                * rair * (1.-sig(nlev))
     *              ) ** (1./3.)
c-----
c old:
c           zwt = ytflx(ji,jj) / (cpair*yrho(ji,nlev)*zarea(ji,nlev))
c           zw(ji,nlev) = max (ztauv, (abs(zwt))**(1./3.))
c           zt(ji,nlev) = zwt / zw(ji,nlev)
c           if (abs(zt(ji,nlev)).gt.10.) then
c             zt(ji,nlev) = min (10., max(-10., zt(ji,nlev)))
c             zw(ji,nlev) = zwt / zt(ji,nlev)
c           endif
c           zmfbot = zarea(ji,nlev) * yrho(ji,nlev) * zw(ji,nlev)
c           zq(ji,nlev) = yqflx(ji,jj) / zmfbot
c           zu(ji,nlev) = yuflx(ji,jj) / zmfbot
c           zv(ji,nlev) = yvflx(ji,jj) / zmfbot
c-----
c new:
            zfac = 1. / zarea(ji,nlev)
            zwt = ytflx(ji,jj) / (yrho(ji,nlev)*cpair)
            zw(ji,nlev) = sqrt(zfac) * max(ztauv,zconv)
            zt(ji,nlev) = zfac*zwt / zw(ji,nlev)
            if (abs(zt(ji,nlev)).gt.10.) then
              zt(ji,nlev) = min (10., max(-10., zt(ji,nlev)))
              zw(ji,nlev) = zfac*zwt / zt(ji,nlev)
            endif
            zq(ji,nlev)= zfac*yqflx(ji,jj)/(yrho(ji,nlev)*zw(ji,nlev))
            zu(ji,nlev)= zfac*yuflx(ji,jj)/(yrho(ji,nlev)*zw(ji,nlev))
            zv(ji,nlev)= zfac*yvflx(ji,jj)/(yrho(ji,nlev)*zw(ji,nlev))
c-----
            zh(ji,nlev) = zt(ji,nlev) / ypmcap(ji,jj,nlev)
            zws(ji,nlev)= zw(ji,nlev)**2
            zrho(ji,nlev)= ypm(ji,jj,nlev)
     *                  / ( rair*(1.+zvir*(zq(ji,nlev)+yq(ji,jj,nlev)))
     *                      *(zt(ji,nlev)+yt(ji,jj,nlev)) )
#ifdef isotrac
            if (yqflx(ji,jj).ne.0.) then
               zq_o(ji,nlev) = zq(ji,nlev)*yqflx_o(ji,jj)/yqflx(ji,jj)
               zq_h(ji,nlev) = zq(ji,nlev)*yqflx_h(ji,jj)/yqflx(ji,jj)
            else
               zq_o(ji,nlev) = 0.
               zq_h(ji,nlev) = 0.
            endif
#endif
#ifdef isotracreg
            do it=1,nrtr
              if (yqflx(ji,jj).ne.0.) then
                zq_t(ji,nlev,it) = zq(ji,nlev)*
     *                                    rtrfv(ji,jj,it)/yqflx(ji,jj)
              else
                zq_t(ji,nlev,it) = 0.
              endif
            enddo
#endif

          endif

  150   continue

c          In same way, initialize pbl plume tracer concentrations
 
        if (ntraca.lt.ntrace) then
          do 1500 n=ntraca+1,ntrace
            do 1502 ji=1,nlon
              if (zbuoy(ji).gt. 0.) then
                zfac = 1. / zarea(ji,nlev)
                za(ji,nlev,n) = zfac*traflux(ji,jj,n-ntraca)
     *                          / (yrho(ji,nlev)*zw(ji,nlev))
              endif
 1502       continue
 1500     continue
        endif
 
c---------
      else
c---------
 
c          For "free" plumes, set initial plume area and radius.
c          Nb: plume radius affects the physics only via the
c          entrainment and drag coefficients zent,zdrag.
c          Plume area does not affect within-plume physics, only the
c          net flux convergences for changing large-scale quantities.
 
c       jk = jkinit
c       do 160 ji=1,nlon
cloc      zarea(ji,jk) = xare2                                     !cloc
cloc      zrad (ji,jk) = xrad2                                     !cloc
c 160   continue
 
c       Free plume radius dep on column precipitable water vapor:
c       call zero (zpwv, nlon)
c       do 160 jk=1,nlev
c         do 162 ji=1,nlon
c           zpwv(ji) = zpwv(ji) + yq(ji,jj,jk)*ypthic(ji,jj,jk)
c 162     continue
c 160   continue
c       do 164 ji=1,nlon
c         zarea(ji,jkinit) = plumearea2
c         zrad(ji,jkinit) = max (plumerad1, plumerad2*(zpwv(ji)/gravit))
c 164   continue
 
        do 160 ji=1,nlon
          zarea(ji,jkinit) = plumearea2   ! .010
          zrad(ji,jkinit) = plumerad2     ! 2000.
  160   continue
 
c----------
      endif
c----------
 
c        Integrate plumes upwards from initiating layer to top
 
c================================
      do 200 jk = jkinit-1, 1, -1
c================================
 
        jkp = jk+1
 
c          Step plume from midpoint of layer jkp to midpoint of layer jk
 
        do 210 ji=1,nlon
 
c            Force for free plume, or non-zero pbl plume from below
 
          if( (.not.ifpbl)                   .or.zw(ji,jkp).gt.0.)!csing
cmult     if(((.not.ifpbl).and.jkp.eq.jkinit).or.zw(ji,jkp).gt.0.)!cmult
     *      then
 
            zent(ji)  = 0.20 / zrad(ji,jkp)
            zdrag(ji) = 0.5  / (pi*zrad(ji,jkp))
 
            zexp = exp(-zent(ji)*dz(ji,jk))
            zh(ji,jk) = zh(ji,jkp)*zexp
     *                  - (dhdz(ji,jk)/zent(ji))  *(1.-zexp)
            zq(ji,jk) = zq(ji,jkp)*zexp
     *                  - (dqdz(ji,jk)/zent(ji))  *(1.-zexp)
            zu(ji,jk) = zu(ji,jkp)*zexp
     *                  - (dudz(ji,jk)/zent(ji))  *(1.-zexp)
            zv(ji,jk) = zv(ji,jkp)*zexp
     *                  - (dvdz(ji,jk)/zent(ji))  *(1.-zexp)
            zlw(ji,jk)= zlw(ji,jkp)*zexp
     *                  - (dlwcdz(ji,jk)/zent(ji))*(1.-zexp)
            zt(ji,jk) = zh(ji,jk)*ypmcap(ji,jj,jk)
#ifdef isotrac
            zq_o(ji,jk) = zq_o(ji,jkp)*zexp
     *                    - (dqdz_o(ji,jk)/zent(ji))  *(1.-zexp)
            zq_h(ji,jk) = zq_h(ji,jkp)*zexp
     *                    - (dqdz_h(ji,jk)/zent(ji))  *(1.-zexp)
            zlw_o(ji,jk)= zlw_o(ji,jkp)*zexp
     *                    - (dlwcdz_o(ji,jk)/zent(ji))*(1.-zexp)
            zlw_h(ji,jk)= zlw_h(ji,jkp)*zexp
     *                    - (dlwcdz_h(ji,jk)/zent(ji))*(1.-zexp)
#endif
#ifdef isotracreg
            do it=1,nrtr
              zq_t(ji,jk,it) = zq_t(ji,jkp,it)*zexp
     *                    - (dqdz_t(ji,jk,it)/zent(ji))  *(1.-zexp)
            enddo
#endif
 
c              Test for parcel saturation, and condense amount that
c              results in exactly saturated air (zcond, zlath).
c              Add to plume liquid water content, converting excess
c              over maximum (namelist param clw_p) to within-plume
c              precip (zprec).
c              zlath is latent heat of condensation, determined by
c              ambient level temperature (ie, rain or snow), plus
c              correction for specific heats to be consistent with lsx
c              (imagining all phase changes occur at tmelt).
 
            ztabs = zt(ji,jk) + yt(ji,jj,jk)
            zqabs = zq(ji,jk) + yq(ji,jj,jk)
            zqsat = qstblf (ztabs, ypm(ji,jj,jk))
c===================================
            if (zqabs.gt.zqsat) then
c===================================
c             zlath(ji,jk) = cvmgt ( latvap + cpwv*(ztabs-tmelt)
c    *                               - ch2o*(yt(ji,jj,jk)-tmelt),
c    *                               latsub + cpwv*(ztabs-tmelt)
c    *                               - cice*(yt(ji,jj,jk)-tmelt),
c    *                               yt(ji,jj,jk).ge.tmelt )
c             For now, all cloud water is liquid with no specific heat
              zlath(ji,jk) = latvap + cpwv*(ztabs-tmelt)
              dqz = (zlath(ji,jk)/rh2o) * zqsat / (ztabs**2)
              cpz = cpair*(1.+cpvir*zqabs)
              delq = (zqabs-zqsat) / (1.+dqz*zlath(ji,jk)/cpz)
              zqnew = zqabs - delq
              delt = delq*zlath(ji,jk) / (cpair*(1.+cpvir*zqnew))
              ztnew = ztabs + delt
              zq(ji,jk) = zqnew - yq(ji,jj,jk)
              zt(ji,jk) = ztnew - yt(ji,jj,jk)
              zh(ji,jk) = zt(ji,jk)/ypmcap(ji,jj,jk)
#ifdef isotrac
              zqabs_o = zq_o(ji,jk)+yq_o(ji,jj,jk)
              zqabs_h = zq_h(ji,jk)+yq_h(ji,jj,jk)
              zqnew_o = zqabs_o
              zqnew_h = zqabs_h
              zlwt_i  = zlw(ji,jk)+ylwc_t(ji,jk)
              zlwt_o  = zlw_o(ji,jk)+ylwc_t_o(ji,jk)
              zlwt_h  = zlw_h(ji,jk)+ylwc_t_h(ji,jk)
              zlwtt_o = zlwt_o
              zlwtt_h = zlwt_h
#endif
 
              zlwnew = zlw(ji,jk) + ylwc_t(ji,jk) + delq
c             zclw_p = clw_p                               ! (kg/kg dep)
              zclw_p = clw_p / yrho(ji,jk)                 ! (kg/m3 dep)
              delw = max (0., zlwnew-zclw_p)
              zlw(ji,jk) = min (zlwnew, zclw_p) - ylwc_t(ji,jk)

#ifdef isotrac
              zlwtt_i  = zlw(ji,jk)+ylwc_t(ji,jk)

c             If the temperature is above the transition temperature
c             for isotope's in conv. clouds, equilibrate first the liq.
c             with the vapour, then call prm_i to calculate the new
c             isotopic compositions of the liquid and the vapour.
c             If not, call prm_i to calculate the isotopic composition
c             of the newly formed solid independently of the solid
c             already present and deduce the new isotopic composition
c             of the solid and the vapour.

              if ((ztabs.gt.t_c_i).and.(ztnew.ge.t_c_i)) then
c               liquid condensate:
                call equil_i (1., ztabs, ypm(ji,jj,jk), zqabs,
     *                      zqnew_o, zqnew_h, zlwt_i, zlwtt_o, zlwtt_h)
                call prm_i (s_c_i, t_c_i, ztabs, ztnew,
     *                             zqabs,  zqnew,   zqnew_o, zqnew_h,
     *                             zlwt_i, zlwtt_i, zlwtt_o, zlwtt_h)
              else
c               solid condensate:
                if (zlwtt_i.ge.zlwt_i) then
                  call prm_i (s_c_i, t_c_i, ztabs, ztnew,
     *                       zqabs, zqnew,          zqnew_o, zqnew_h,
     *                       0.,    zlwtt_i-zlwt_i, zlwtt_o, zlwtt_h)
                  zlwtt_o = zlwt_o+zlwtt_o
                  zlwtt_h = zlwt_h+zlwtt_h
                else
                  call prm_i (s_c_i, t_c_i, ztabs, ztnew,
     *                             zqabs,  zqnew,   zqnew_o, zqnew_h,
     *                             0.,     0.,      zlwtt_o, zlwtt_h)
                  zlwtt_o = zlwtt_i*zlwt_o/zlwt_i
                  zlwtt_h = zlwtt_i*zlwt_h/zlwt_i
                endif
              endif

              delq_o = zqabs_o-zqnew_o
              delq_h = zqabs_h-zqnew_h
              delw_o = cvmgt(delq_o-(zlwtt_o-zlwt_o), 0., delw.gt.0.)
              delw_h = cvmgt(delq_h-(zlwtt_h-zlwt_h), 0., delw.gt.0.)
              zq_o(ji,jk) = zq_o(ji,jk) -delq_o
              zq_h(ji,jk) = zq_h(ji,jk) -delq_h
              zlw_o(ji,jk)= zlw_o(ji,jk)+delq_o-delw_o
              zlw_h(ji,jk)= zlw_h(ji,jk)+delq_h-delw_h
#endif
#ifdef isotracreg
              do it=1,nrtr
                zrq_t      = (zq_t(ji,jk,it)+rtrq(ji,jj,jk,it))/zqabs
                delq_t(it) = delq*zrq_t
                delw_t(it) = delw*zrq_t
                zq_t(ji,jk,it) = zq_t(ji,jk,it)-delq_t(it)
              enddo
#endif
c===============
            else
c===============
              delq = 0.
              delt = 0.
              delw = 0.
#ifdef isotrac
              delq_o = 0.
              delq_h = 0.
              delw_o = 0.
              delw_h = 0.
#endif
#ifdef isotracreg
              do it=1,nrtr
                delq_t(it) = 0.
                delw_t(it) = 0.
              enddo
#endif
c================
            endif
c================
 
            zrho(ji,jk) = ypm(ji,jj,jk)
     *                  / ( rair*(1.+zvir*(zq(ji,jk)+yq(ji,jj,jk)))
     *                      *(zt(ji,jk)+yt(ji,jj,jk)) )
 
            drhog0= gravit * (yrho(ji,jkp)-zrho(ji,jkp))/zrho(ji,jkp)
            drhog1= gravit * (yrho(ji,jk) -zrho(ji,jk) )/zrho(ji,jk)
            a_rho =  (drhog0 + (drhog1-drhog0)/(1.-zexp)) / zent(ji)
            b_rho = -(         (drhog1-drhog0)/(1.-zexp)) / zent(ji)
 
            zexp2 = zexp*zexp
c           or: include vertical drag
c           zexp2 = exp(-2.*(zent(ji)+zdrag(ji))*dz(ji,jk))
c           a_rho = a_rho* zent(ji) / (zent(ji)+zdrag(ji))
c           b_rho = b_rho* zent(ji) / (2.*(zent(ji)+zdrag(ji))-zent(ji))
 
            zws(ji,jk) = zws(ji,jkp)*zexp2
     *                 + a_rho*(1.-zexp2) + 2.*b_rho*(zexp-zexp2)
            zw(ji,jk) = sqrt ( max(zws(ji,jk),0.) )
 
c             Test if detrainment >=0 is violated. If so, do corrections
 
            if (zw(ji,jk).gt.0.) then
              drho0 = yrho(ji,jkp)-zrho(ji,jkp)
              drho1 = yrho(ji,jk) -zrho(ji,jk)
              zwsmax = (   zrho(ji,jkp)*zws(ji,jkp)
     *                   + gravit*0.5*(drho0+drho1)*dz(ji,jk)
     *                 ) / zrho(ji,jk)
              if (zwsmax.gt.0. .and. zws(ji,jk).gt.zwsmax) then
                zws(ji,jk) = zwsmax
                zw(ji,jk) = sqrt(zws(ji,jk))
                zma = zrho(ji,jkp)*zw(ji,jkp)
                zmb = zrho(ji,jk) *zw(ji,jk)
                zent(ji) = (zmb-zma)/(0.5*(zma+zmb)*dz(ji,jk))
                zexp = exp(-zent(ji)*dz(ji,jk))
                zh(ji,jk) = zh(ji,jkp)*zexp
     *                      - (dhdz(ji,jk)/zent(ji))  *(1.-zexp)
     *                      + delt/ypmcap(ji,jj,jk)
                zq(ji,jk) = zq(ji,jkp)*zexp
     *                      - (dqdz(ji,jk)/zent(ji))  *(1.-zexp)
     *                      - delq
                zu(ji,jk) = zu(ji,jkp)*zexp
     *                      - (dudz(ji,jk)/zent(ji))  *(1.-zexp)
                zv(ji,jk) = zv(ji,jkp)*zexp
     *                      - (dvdz(ji,jk)/zent(ji))  *(1.-zexp)
                zlw(ji,jk)= zlw(ji,jkp)*zexp
     *                      - (dlwcdz(ji,jk)/zent(ji))*(1.-zexp)
     *                      + delq - delw
                zt(ji,jk) = zh(ji,jk)*ypmcap(ji,jj,jk)
#ifdef isotrac
                zq_o(ji,jk) = zq_o(ji,jkp)*zexp
     *                        - (dqdz_o(ji,jk)/zent(ji))  *(1.-zexp)
     *                        - delq_o
                zq_h(ji,jk) = zq_h(ji,jkp)*zexp
     *                        - (dqdz_h(ji,jk)/zent(ji))  *(1.-zexp)
     *                        - delq_h
                zlw_o(ji,jk)= zlw_o(ji,jkp)*zexp
     *                        - (dlwcdz_o(ji,jk)/zent(ji))*(1.-zexp)
     *                        + delq_o - delw_o
                zlw_h(ji,jk)= zlw_h(ji,jkp)*zexp
     *                        - (dlwcdz_h(ji,jk)/zent(ji))*(1.-zexp)
     *                        + delq_h - delw_h
#endif
#ifdef isotracreg
                do it=1,nrtr
                  zq_t(ji,jk,it) = zq_t(ji,jkp,it)*zexp
     *                        - (dqdz_t(ji,jk,it)/zent(ji))*(1.-zexp)
     *                        - delq_t(it)
                enddo
#endif
              endif
            endif
 
c              If vertical velocity = 0, end of plume;
c              else, re-compute zu and zv with (non-linear) form drag
 
            if (zw(ji,jk).eq.0.) then
 
              zh(ji,jk) = 0.
              zq(ji,jk) = 0.
              zu(ji,jk) = 0.
              zv(ji,jk) = 0.
              zlw(ji,jk) = 0.
              zt(ji,jk) = 0.
              zrho(ji,jk) = yrho(ji,jk)
              zcond(ji,jk) = 0.
              zlath(ji,jk) = 0.
              zprec(ji,jk) = 0.
#ifdef isotrac
              zq_o(ji,jk) = 0.
              zq_h(ji,jk) = 0.
              zlw_o(ji,jk) = 0.
              zlw_h(ji,jk) = 0.
              zcond_o(ji,jk) = 0.
              zcond_h(ji,jk) = 0.
              zprec_o(ji,jk) = 0.
              zprec_h(ji,jk) = 0.
#endif
#ifdef isotracreg
              do it=1,nrtr
                zq_t(ji,jk,it) = 0.
                zcond_t(ji,jk,it) = 0.
                zprec_t(ji,jk,it) = 0.
              enddo
#endif
c             Either:
c             heipbl(ji)= heipbl(ji) + dz(ji,jk) * zws(ji,jkp)
c    *                / (max(zws(ji,jkp),1.e-10) - min(zws(ji,jk),0.))
c             Or:
c             Solve quadratic eqn for exp(-zent*z),ie,for z where w=0,
c             for accurate increment to pbl height (diagnostic only).
c             Only valid for no vertical drag(if commented out above).
c             cvmgt(zz1) mod is fudge to avoid quadratic singularity.
              if (ifpbl) then
                zz1 = zws(ji,jkp) - a_rho - 2.*b_rho
                zz1 = cvmgt (zz1, zz1+1.e-9, abs(zz1).gt.1.e-10)
                zz2 = (-b_rho + sqrt(max(0.,b_rho**2-zz1*a_rho))) / zz1
                heipbl(ji) = heipbl(ji) - log(zz2)/zent(ji)
                laypbl(ji) = jkp
              else
                if (zw(ji,jkp).gt.0.) laytop(ji) = min (laytop(ji),jkp)
              endif
 
              zws(ji,jk) = 0.
 
            else
 
              zenu= zent(ji) + zdrag(ji)*abs(0.5*(zu(ji,jkp)+zu(ji,jk)))
     *                         / (0.5*(zw(ji,jkp)+zw(ji,jk)))
              zenv= zent(ji) + zdrag(ji)*abs(0.5*(zv(ji,jkp)+zv(ji,jk)))
     *                         / (0.5*(zw(ji,jkp)+zw(ji,jk)))
              zexu  = exp(-zenu*dz(ji,jk))
              zexv  = exp(-zenv*dz(ji,jk))
              zu(ji,jk) = zu(ji,jkp)*zexu - (dudz(ji,jk)/zenu)*(1.-zexu)
              zv(ji,jk) = zv(ji,jkp)*zexv - (dvdz(ji,jk)/zenv)*(1.-zexv)
 
              zcond(ji,jk) = zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delq
              zprec(ji,jk) = zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delw
#ifdef isotrac
              zcond_o(ji,jk)=zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delq_o
              zcond_h(ji,jk)=zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delq_h
              zprec_o(ji,jk)=zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delw_o
              zprec_h(ji,jk)=zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delw_h
#endif
#ifdef isotracreg
              do it=1,nrtr
                zcond_t(ji,jk,it)=zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)
     *                            *delq_t(it)
                zprec_t(ji,jk,it)=zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)
     *                            *delw_t(it)
              enddo
#endif
              if (ifpbl) then
                heipbl(ji) = heipbl(ji) + dz(ji,jk)
              else
                laybot(ji) = max (laybot(ji),jkp)
              endif
 
            endif
 
          endif
 
  210   continue
 
c          In same way, integrate plume tracers through curr layer
 
        if (ntraca.lt.ntrace) then
          do 2100 n=ntraca+1,ntrace
            do 2102 ji=1,nlon
              if (zw(ji,jk).gt.0.) then
                zexp = exp(-zent(ji)*dz(ji,jk))
                za(ji,jk,n) = za(ji,jkp,n)*zexp
     *                      - (dadz(ji,jk,n)/zent(ji))*(1.-zexp)
              else
                za(ji,jk,n) = 0.
              endif
 2102       continue
 2100     continue
        endif
 
c          Set new plume radius and area.
c          Also count number of plumes, skip out of vertical
c          integration for pbl or multiple-initiation plumes
c          if none left for this strip.
 
        nplu = 0
        do 230 ji=1,nlon
          zrad(ji,jk)  = zrad(ji,jkp)
          zarea(ji,jk) = zarea(ji,jkp)
c         zw(ji,jk)  = zw(ji,jk) * zarea(ji,jkp)/zarea(ji,jk)
c         zws(ji,jk) = zw(ji,jk)**2
          if (zw(ji,jk).gt.0.) nplu = nplu + 1
  230   continue
        if (ifpbl .and.  nplu.eq.0) goto 202                      !csing
cmult   if (             nplu.eq.0) goto 202                      !cmult
 
c=============
  200 continue
  202 continue
c=============
 
c        Accumulate large-scale fluxes at layer midpoints.
c        Also accumulate mass flux zmfa for cfl check below
 
      do 300 jk=1,nlev
 
        do 302 ji=1,nlon
          zmf(ji,jk)  = zarea(ji,jk)*zrho(ji,jk)*zw(ji,jk)
          zmfa(ji,jk) = zmfa(ji,jk) + zmf(ji,jk)
  302   continue
 
        if (ifpbl .and. jk.eq.nlev) then
          do 304 ji=1,nlon
            tplu(ji,jk) = tplu(ji,jk) + ytflx(ji,jj)
            qplu(ji,jk) = qplu(ji,jk) + yqflx(ji,jj)
            uplu(ji,jk) = uplu(ji,jk) + yuflx(ji,jj)
            vplu(ji,jk) = vplu(ji,jk) + yvflx(ji,jj)
#ifdef isotrac
            qplu_o(ji,jk) = qplu_o(ji,jk) + yqflx_o(ji,jj)
            qplu_h(ji,jk) = qplu_h(ji,jk) + yqflx_h(ji,jj)
#endif
#ifdef isotracreg
            do it=1,nrtr
              qplu_t(ji,jk,it)= qplu_t(ji,jk,it) + rtrfv(ji,jj,it)
            enddo
#endif
  304     continue
        else
          do 306 ji=1,nlon
            tplu(ji,jk) = tplu(ji,jk) + zt(ji,jk)*zmf(ji,jk)*cpair
            qplu(ji,jk) = qplu(ji,jk) + zq(ji,jk)*zmf(ji,jk)
c           comment out next two lines to turn off "cumulus friction"
            uplu(ji,jk) = uplu(ji,jk) + zu(ji,jk)*zmf(ji,jk)
            vplu(ji,jk) = vplu(ji,jk) + zv(ji,jk)*zmf(ji,jk)
 
            cpluc(ji,jk) = cpluc(ji,jk)
     *                   + ( zlw(ji,jk) + ylwc_t(ji,jk)
     *                                  - ylwc_c(ji,jj,jk) )*zmf(ji,jk)
            cplua(ji,jk) = cplua(ji,jk) - ylwc_a(ji,jj,jk)  *zmf(ji,jk)
            cplus(ji,jk) = cplus(ji,jk) - ylwc_s(ji,jj,jk)  *zmf(ji,jk)
 
            zconda(ji,jk)= zconda(ji,jk) + zcond(ji,jk)
            zheata(ji,jk)= zheata(ji,jk) + zcond(ji,jk)*zlath(ji,jk)
            zpreca(ji,jk)= zpreca(ji,jk) + zprec(ji,jk)
#ifdef isotrac
            qplu_o(ji,jk) = qplu_o(ji,jk)+zq_o(ji,jk)*zmf(ji,jk)
            qplu_h(ji,jk) = qplu_h(ji,jk)+zq_h(ji,jk)*zmf(ji,jk)
            cpluc_o(ji,jk)= cpluc_o(ji,jk)
     *                   +( zlw_o(ji,jk)+ylwc_t_o(ji,jk)
     *                                  -ylwc_c_o(ji,jj,jk) )*zmf(ji,jk)
            cpluc_h(ji,jk)= cpluc_h(ji,jk)
     *                   +( zlw_h(ji,jk)+ylwc_t_h(ji,jk)
     *                                  -ylwc_c_h(ji,jj,jk) )*zmf(ji,jk)
            cplua_o(ji,jk)= cplua_o(ji,jk)-ylwc_a_o(ji,jj,jk)*zmf(ji,jk)
            cplua_h(ji,jk)= cplua_h(ji,jk)-ylwc_a_h(ji,jj,jk)*zmf(ji,jk)
            cplus_o(ji,jk)= cplus_o(ji,jk)-ylwc_s_o(ji,jj,jk)*zmf(ji,jk)
            cplus_h(ji,jk)= cplus_h(ji,jk)-ylwc_s_h(ji,jj,jk)*zmf(ji,jk)
            zconda_o(ji,jk)= zconda_o(ji,jk)+zcond_o(ji,jk)
            zconda_h(ji,jk)= zconda_h(ji,jk)+zcond_h(ji,jk)
            zpreca_o(ji,jk)= zpreca_o(ji,jk)+zprec_o(ji,jk)
            zpreca_h(ji,jk)= zpreca_h(ji,jk)+zprec_h(ji,jk)
#endif
#ifdef isotracreg
            do it=1,nrtr
              qplu_t(ji,jk,it)=qplu_t(ji,jk,it)+
     *                                         zq_t(ji,jk,it)*zmf(ji,jk)
              zconda_t(ji,jk,it)=zconda_t(ji,jk,it)+zcond_t(ji,jk,it)
              zpreca_t(ji,jk,it)=zpreca_t(ji,jk,it)+zprec_t(ji,jk,it)
            enddo
#endif
  306     continue
        endif
 
        if (ntraca.lt.ntrace) then
          do 3040 n=ntraca+1,ntrace
            if (ifpbl .and. jk.eq.nlev) then
              do 3042 ji=1,nlon
                aplu(ji,jk,n) = aplu(ji,jk,n) + traflux(ji,jj,n-ntraca)
 3042         continue
            else
              do 3062 ji=1,nlon
                aplu(ji,jk,n) = aplu(ji,jk,n) + za(ji,jk,n)*zmf(ji,jk)
 3062         continue
            endif
 3040     continue
        endif
 
  300 continue
 
c        For "pbl" plumes,set history 03:00.15:00 (was daily min/max)
c        pbl heights, accumulate diagnostic number of pbl-tops in each
c        layer(numpbl), and save diagnostic zonal mean plume vertical
c        velocities in wplume0. For "free" plumes, accumulate wfre1,2
c        for later setting (after all initiating layers done) of wplume1
 
c--------------------
      if (ifpbl) then
c--------------------
 
        do 400 ji=1,nlon
          if (laypbl(ji).ne.nlevp)
     *      numpbl(laypbl(ji)) = numpbl(laypbl(ji)) + 1
c         For daily min/max pbl heights:
c         heimin(ji,jj) = min (heipbl(ji), heimin(ji,jj))
c         heimax(ji,jj) = max (heipbl(ji), heimax(ji,jj))
  400   continue
 
c       For 03:00, 15:00 pbl heights (ji03a,etc, may be < 1 or > nlon):
        do 404 ji = ji03a,ji03b
          jix = mod (ji+nlon-1, nlon) + 1
          heimin(jix,jj) = heipbl(jix)
  404   continue
        do 406 ji = ji15a,ji15b
          jix = mod (ji+nlon-1, nlon) + 1
          heimax(jix,jj) = heipbl(jix)
  406   continue
 
        do 410 jk=1,nlev
          wplume0(jj,jk) = 0.
          nzw = 0
          do 412 ji=1,nlon
            if (zw(ji,jk).gt.0.) then
              wplume0(jj,jk) = wplume0(jj,jk) + zw(ji,jk)
              nzw = nzw + 1
            endif
  412     continue
          if (nzw.gt.0) wplume0(jj,jk) = wplume0(jj,jk)/nzw
  410   continue
 
c---------
      else
c---------
 
        do 420 jk=1,nlev
          do 422 ji=1,nlon
c           wfre1(ji,jk) = wfre1(ji,jk) + zarea(ji,jk)*zrho(ji,jk)
c    *                                    *zw(ji,jk)
c           wfre2(ji,jk) = wfre2(ji,jk) + zarea(ji,jk)*zrho(ji,jk)
            wfre1(ji,jk) = max (wfre1(ji,jk),zw(ji,jk))
            wfre2(ji,jk) = 1.
  422     continue
  420   continue
 
c----------
      endif
c----------
 
c        End of loop over initiating layers
 
c*************
  510 continue
c*************
 
c        Set diagnostic zonal "free" plume bottom and top locations
c        and vertical velocity wplume1
 
      if (.not.ifpbl) then
 
        do 520 ji=1,nlon
          if (laytop(ji).ne.nlevp)
     *      numfret(laytop(ji)) = numfret(laytop(ji)) + 1
          if (laybot(ji).ne.0)
     *      numfreb(laybot(ji)) = numfreb(laybot(ji)) + 1
  520   continue
 
        do 530 jk=1,nlev
          wplume1(jj,jk) = 0.
          nzw = 0
          do 532 ji=1,nlon
            if (wfre1(ji,jk).gt.0. .and. wfre2(ji,jk).gt.0.) then
              nzw = nzw + 1
              wplume1(jj,jk)= wplume1(jj,jk) + wfre1(ji,jk)/wfre2(ji,jk)
            endif
  532     continue
          if (nzw.gt.0) wplume1(jj,jk) = wplume1(jj,jk)/nzw
  530   continue
 
      endif

#ifdef isotrac
      if (ifpbl) then
        do 540 jk=1,nlev
          do 542 ji=1,nlon
            zapbl(ji,jk) = zarea(ji,jk)
            zrpbl(ji,jk) = yrho(ji,jk)*zw(ji,jk)
 542      continue
 540    continue
      endif
#endif
 
c        End of "pbl" or "free" loop
 
c*************
  500 continue
c*************
 
c        Check CFL mass-flux  criterion, find worst violation for each
c        column, and reduce fluxes uniformly for a column to make CFL
c        ok. (But not for bottom layer, otherwise wouldn't conserve
c        surface fluxes).
 
      call zero (zcfl, nlon)
 
      do 600 jk=1,nlev-1
        do 602 ji=1,nlon
          zcfl(ji) = max ( zcfl(ji),
     *                     zmfa(ji,jk)*dtplume*gravit
     *                     / min (ypthic(ji,jj,jk),ypthic(ji,jj,jk+1))
     *                   )
  602   continue
  600 continue
 
c     cflmax = 1.e6
      cflmax = 0.75
      cflmax_t = 0.50
 
      do 610 jk=1,nlev-1
        do 612 ji=1,nlon
          zfac = cflmax / max (cflmax, zcfl(ji))
          zmfa(ji,jk)  = zmfa(ji,jk)  * zfac
          tplu(ji,jk)  = tplu(ji,jk)  * zfac
          qplu(ji,jk)  = qplu(ji,jk)  * zfac
          uplu(ji,jk)  = uplu(ji,jk)  * zfac
          vplu(ji,jk)  = vplu(ji,jk)  * zfac
          cpluc(ji,jk) = cpluc(ji,jk) * zfac
          cplua(ji,jk) = cplua(ji,jk) * zfac
          cplus(ji,jk) = cplus(ji,jk) * zfac
          zconda(ji,jk)= zconda(ji,jk)* zfac
          zheata(ji,jk)= zheata(ji,jk)* zfac
          zpreca(ji,jk)= zpreca(ji,jk)* zfac
#ifdef isotrac
          qplu_o(ji,jk)  = qplu_o(ji,jk)  * zfac
          qplu_h(ji,jk)  = qplu_h(ji,jk)  * zfac
          cpluc_o(ji,jk) = cpluc_o(ji,jk) * zfac
          cpluc_h(ji,jk) = cpluc_h(ji,jk) * zfac
          cplua_o(ji,jk) = cplua_o(ji,jk) * zfac
          cplua_h(ji,jk) = cplua_h(ji,jk) * zfac
          cplus_o(ji,jk) = cplus_o(ji,jk) * zfac
          cplus_h(ji,jk) = cplus_h(ji,jk) * zfac
          zconda_o(ji,jk)= zconda_o(ji,jk)* zfac
          zconda_h(ji,jk)= zconda_h(ji,jk)* zfac
          zpreca_o(ji,jk)= zpreca_o(ji,jk)* zfac
          zpreca_h(ji,jk)= zpreca_h(ji,jk)* zfac
#endif
#ifdef isotracreg
          do it=1,nrtr
            qplu_t(ji,jk,it)  = qplu_t(ji,jk,it)  * zfac
            zconda_t(ji,jk,it)= zconda_t(ji,jk,it)* zfac
            zpreca_t(ji,jk,it)= zpreca_t(ji,jk,it)* zfac
          enddo
#endif
  612   continue
 
        if (ntraca.lt.ntrace) then
          do 6120 n=ntraca+1,ntrace
            do 6122 ji=1,nlon
              zfac = cflmax_t / max (cflmax_t, zcfl(ji))
              aplu(ji,jk,n) = aplu(ji,jk,n) * zfac
 6122       continue
 6120     continue
        endif
  610 continue
 
c     zz = 0.
c     do 8880 ji=1,nlon
c       zz = zz +  cflmax / max (cflmax, zcfl(ji))
c8880 continue
c     write(*,8881) nstep, jj, ji, cflmax, zz/nlon
c8881 format('plume: zonal mean cfl: nstep,jj,cflmax,zfac=',2i4,2f6.3)
 
c        Compute flux convergence into each layer and add to reservoirs
c        fplum* ("flux*time per layer", in commun) for use in vdif.
c        Also increment tracer fields themselves (no reservoir).
c        Use dtplume since reservoirs are not leapfrog variables
c        (see comments in reserv).
 
c        For simplicity, shift fluxes at layer midpoints down to next
c        lowest interface, (so give zconda,zheata,zpreca to layer below)
 
      do 700 jk=1,nlev
        do 702 ji=1,nlon
          fplumt(ji,jj,jk) = fplumt(ji,jj,jk)
     *                     + (  tplu(ji,jk)-tplu(ji,jk-1)
     *                        + zheata(ji,jk-1)           ) * dtplume
          fplumq(ji,jj,jk) = fplumq(ji,jj,jk)
     *                     + (  qplu(ji,jk)-qplu(ji,jk-1)
     *                        - zconda(ji,jk-1)           ) * dtplume
          fplumu(ji,jj,jk) = fplumu(ji,jj,jk)
     *                     + (  uplu(ji,jk)-uplu(ji,jk-1) ) * dtplume
          fplumv(ji,jj,jk) = fplumv(ji,jj,jk)
     *                     + (  vplu(ji,jk)-vplu(ji,jk-1) ) * dtplume

#ifdef isotrac
          fplumq_o(ji,jj,jk)=fplumq_o(ji,jj,jk)
     *                     + (  qplu_o(ji,jk)-qplu_o(ji,jk-1)
     *                      - zconda_o(ji,jk-1)           ) * dtplume
          fplumq_h(ji,jj,jk)=fplumq_h(ji,jj,jk)
     *                     + (  qplu_h(ji,jk)-qplu_h(ji,jk-1)
     *                      - zconda_h(ji,jk-1)           ) * dtplume
#ifdef isotracreg
          do it=1,nrtr
            rtrfq(ji,jj,jk,it)=rtrfq(ji,jj,jk,it)
     *                       + (  qplu_t(ji,jk,it)-qplu_t(ji,jk-1,it)
     *                        - zconda_t(ji,jk-1,it)      ) * dtplume
          enddo
#endif
c         Compute and save plume (pbl+free) fractional area and relative
c         humidity for the parameterization of precip,evap in reevap:
          ctw = zrpbl(ji,jk)+yrho(ji,jk)*zw(ji,jk)
          cta = zapbl(ji,jk)*zrpbl(ji,jk)+
     *          zarea(ji,jk)*yrho(ji,jk)*zw(ji,jk)
          if (ctw.gt.0.) za_i(ji,jk) = cta/ctw
          if (cta.ne.0.) then
            tc =           yt(ji,jj,jk) + tplu(ji,jk)/cta/cpair
            zh_i(ji,jk) = (yq(ji,jj,jk) + qplu(ji,jk)/cta)/
     *                                   qstblf(tc, ypm(ji,jj,jk))
          endif
#endif
 
          zqtot = zqtot
     *          + (qplu(ji,jk)-qplu(ji,jk-1))*dtplume*cosbud(jj)/nlon
          zptot = zptot + zconda(ji,jk-1)*dtplume*cosbud(jj)/nlon
          if (ifpbl .and. jk.eq.nlev)
     *      zetot = zetot + yqflx(ji,jj)*dtplume*cosbud(jj)/nlon
  702   continue
 
c          Don't change passive tracers (>= #8) if turned off
c          by namelist traconv (default = .true.)
 
        if (ntraca.lt.ntrace .and. traconv) then
          do 7020 n=ntraca+1,ntrace
            do 7022 ji=1,nlon
              zm = ypthic(ji,jj,jk)/gravit
              tracer(ji,jj,jk,n) = tracer(ji,jj,jk,n)
     *                      + (aplu(ji,jk,n)-aplu(ji,jk-1,n))*dtplume/zm
 7022       continue
 7020     continue
        endif
  700 continue
 
c        Increment large-scale liquid water contents (ylwc_*, kg/kg)
c        set within-plume precipitation (plumeprec, kg/m2/s), increment
c        global budget terms (d[m,h]clouf, kg/m2)
 
      do 720 ji=1,nlon
        ifanv(ji) = .true.
  720 continue
 
      do 730 jk=1,nlev                                   ! top to bottom
        do 732 ji=1,nlon
 
          zm = ypthic(ji,jj,jk)/gravit
          zc = cpluc (ji,jk)   - cpluc(ji,jk-1)
     *       + zconda(ji,jk-1) - zpreca(ji,jk-1)
#ifdef isotrac
          zc_o = cpluc_o (ji,jk)   - cpluc_o(ji,jk-1)
     *         + zconda_o(ji,jk-1) - zpreca_o(ji,jk-1)
          zc_h = cpluc_h (ji,jk)   - cpluc_h(ji,jk-1)
     *         + zconda_h(ji,jk-1) - zpreca_h(ji,jk-1)
#endif
          zlwcold(ji,jk) = ylwc_c(ji,jj,jk)
 
c            Give to anvil if plume top is high enough, else to conv

          if (zc.ne.0.) then
            if (ifanv(ji) .and. sig(jk).lt.siganvil) then
              ylwc_a(ji,jj,jk)= ylwc_a(ji,jj,jk) + dtplume*zc/zm
#ifdef isotrac
              ylwc_a_o(ji,jj,jk)=ylwc_a_o(ji,jj,jk)+dtplume*zc_o/zm
              ylwc_a_h(ji,jj,jk)=ylwc_a_h(ji,jj,jk)+dtplume*zc_h/zm
#endif
            else
              ylwc_c(ji,jj,jk)= ylwc_c(ji,jj,jk) + dtplume*zc/zm
#ifdef isotrac
              ylwc_c_o(ji,jj,jk)=ylwc_c_o(ji,jj,jk)+dtplume*zc_o/zm
              ylwc_c_h(ji,jj,jk)=ylwc_c_h(ji,jj,jk)+dtplume*zc_h/zm
#endif
            endif
            ifanv(ji) = .false.
          endif
 
          ylwc_a(ji,jj,jk) = ylwc_a(ji,jj,jk)
     *                     + (cplua(ji,jk)-cplua(ji,jk-1)) * dtplume/zm
          ylwc_s(ji,jj,jk) = ylwc_s(ji,jj,jk)
     *                     + (cplus(ji,jk)-cplus(ji,jk-1)) * dtplume/zm
#ifdef isotrac
          ylwc_a_o(ji,jj,jk)=ylwc_a_o(ji,jj,jk)
     *                     +(cplua_o(ji,jk)-cplua_o(ji,jk-1))*dtplume/zm
          ylwc_a_h(ji,jj,jk)=ylwc_a_h(ji,jj,jk)
     *                     +(cplua_h(ji,jk)-cplua_h(ji,jk-1))*dtplume/zm
          ylwc_s_o(ji,jj,jk)=ylwc_s_o(ji,jj,jk)
     *                     +(cplus_o(ji,jk)-cplus_o(ji,jk-1))*dtplume/zm
          ylwc_s_h(ji,jj,jk)=ylwc_s_h(ji,jj,jk)
     *                     +(cplus_h(ji,jk)-cplus_h(ji,jk-1))*dtplume/zm
#endif
 
          plumeprec(ji,jj,jk) = zpreca(ji,jk-1)
#ifdef isotrac
          plumepr_o(ji,jj,jk) = zpreca_o(ji,jk-1)
          plumepr_h(ji,jj,jk) = zpreca_h(ji,jk-1)
#endif
#ifdef isotracreg
          do it=1,nrtr
            plumepr_t(ji,jj,jk,it)=zpreca_t(ji,jk-1,it)
          enddo
#endif
          dmclouf(jj) = dmclouf(jj) + zconda(ji,jk-1)*dtplume/nlon
          dhcloud(jj) = dhcloud(jj) + zheata(ji,jk-1)*dtplume/nlon
  732   continue
  730 continue
 
c        Crudely fix very slightly negative cloud liquid water contents
c        (otherwise causes blowups in cldcmp/radctl)
 
      call zero (ztotc1, nlon)
      call zero (ztota1, nlon)
      call zero (ztots1, nlon)
      call zero (ztotc2, nlon)
      call zero (ztota2, nlon)
      call zero (ztots2, nlon)
      do 750 jk=1,nlev
        do 752 ji=1,nlon
          zm = ypthic(ji,jj,jk)/gravit
          ztotc1(ji) = ztotc1(ji) + ylwc_c(ji,jj,jk)*zm
          ztota1(ji) = ztota1(ji) + ylwc_a(ji,jj,jk)*zm
          ztots1(ji) = ztots1(ji) + ylwc_s(ji,jj,jk)*zm
#ifdef isotrac
          yl_c_i(ji,jk) = ylwc_c(ji,jj,jk)
          yl_a_i(ji,jk) = ylwc_a(ji,jj,jk)
          yl_s_i(ji,jk) = ylwc_s(ji,jj,jk)
#endif
          ylwc_c(ji,jj,jk) = max (ylwc_c(ji,jj,jk), 0.)
          ylwc_a(ji,jj,jk) = max (ylwc_a(ji,jj,jk), 0.)
          ylwc_s(ji,jj,jk) = max (ylwc_s(ji,jj,jk), 0.)
          ztotc2(ji) = ztotc2(ji) + ylwc_c(ji,jj,jk)*zm
          ztota2(ji) = ztota2(ji) + ylwc_a(ji,jj,jk)*zm
          ztots2(ji) = ztots2(ji) + ylwc_s(ji,jj,jk)*zm
  752   continue
  750 continue
 
      do 754 ji=1,nlon
c       Since z*2 > z*1, if z*1>0, so is z*2:
        ztotc1(ji) = cvmgt ( ztotc1(ji)/max(ztotc2(ji),1.e-20), 0.,
     *                       ztotc1(ji).gt.0. )
        ztota1(ji) = cvmgt ( ztota1(ji)/max(ztota2(ji),1.e-20), 0.,
     *                       ztota1(ji).gt.0.)
        ztots1(ji) = cvmgt ( ztots1(ji)/max(ztots2(ji),1.e-20), 0.,
     *                       ztots1(ji).gt.0.)
  754 continue
 
      do 756 jk=1,nlev
        do 758 ji=1,nlon
          ylwc_c(ji,jj,jk) = ztotc1(ji) * ylwc_c(ji,jj,jk)
          ylwc_a(ji,jj,jk) = ztota1(ji) * ylwc_a(ji,jj,jk)
          ylwc_s(ji,jj,jk) = ztots1(ji) * ylwc_s(ji,jj,jk)
#ifdef isotrac
          if (yl_c_i(ji,jk).ne.0.) then
            ylwc_c_o(ji,jj,jk)=ylwc_c_o(ji,jj,jk)*ylwc_c(ji,jj,jk)
     *                         /yl_c_i(ji,jk)
            ylwc_c_h(ji,jj,jk)=ylwc_c_h(ji,jj,jk)*ylwc_c(ji,jj,jk)
     *                         /yl_c_i(ji,jk)
          else
            ylwc_c_o(ji,jj,jk)= 0.
            ylwc_c_h(ji,jj,jk)= 0.
          endif
          if (yl_a_i(ji,jk).ne.0.) then
            ylwc_a_o(ji,jj,jk)=ylwc_a_o(ji,jj,jk)*ylwc_a(ji,jj,jk)
     *                         /yl_a_i(ji,jk)
            ylwc_a_h(ji,jj,jk)=ylwc_a_h(ji,jj,jk)*ylwc_a(ji,jj,jk)
     *                         /yl_a_i(ji,jk)
          else
            ylwc_a_o(ji,jj,jk)= 0.
            ylwc_a_h(ji,jj,jk)= 0.
          endif
          if (yl_s_i(ji,jk).ne.0.) then
            ylwc_s_o(ji,jj,jk)=ylwc_s_o(ji,jj,jk)*ylwc_s(ji,jj,jk)
     *                         /yl_s_i(ji,jk)
            ylwc_s_h(ji,jj,jk)=ylwc_s_h(ji,jj,jk)*ylwc_s(ji,jj,jk)
     *                         /yl_s_i(ji,jk)
          else
            ylwc_s_o(ji,jj,jk)= 0.
            ylwc_s_h(ji,jj,jk)= 0.
          endif
#endif
  758   continue
  756 continue
 
c*****
c     Diagnostic printout: cloud liquid water
c     do 8882 ji=1,nlon
c       do 8883 jk=1,nlev
c         if (ylwc_c(ji,jj,jk).lt.0.) then
c           write (*,*) 'plume: ji=',ji,' jj=',jj,' jk=',jk
c           write (*,8887)
c           do 8884 k=1,nlev
c             zm = ypthic(ji,jj,k)/gravit
c             zca = (cpluc(ji,k)-cpluc(ji,k-1))*dtplume/zm
c             zcb = zconda(ji,k-1)*dtplume/zm
c             zcc = -zpreca(ji,k-1)*dtplume/zm
c             write(*,8888) k,
c    *                      zw(ji,k),
c    *                      zmf(ji,k)*dtplume*gravit/ypthic(ji,jj,k),
c    *                      ylwc_t(ji,k)*1.e3,
c    *                      (zlw(ji,k)+ylwc_t(ji,k))*1.e3,
c    *                      zlwcold(ji,k)*1.e3,
c    *                      ylwc_c(ji,jj,k)*1.e3,
c    *                      zca*1.e3,
c    *                      zcb*1.e3,
c    *                      zcc*1.e3
c8887         format(' k', 2x, '         zw', '        zmf',
c    *                        '      ylwc_t', '        zlw',
c    *                        '    ylw_cold', '    ylw_cnew',
c    *                        '     d(cplu)',
c    *                        '   d(zconda)', '   d(zpreca)')
c8888         format(i2, 2x, 9f12.8)
c8884       continue
c           call endrun (-1)
c         endif
c8883   continue
c8882 continue
c*****
 
      iuprinz(1) = nout
cloc  iuprinz(1) = 6
cloc  iuprinz(2) = iuout
 
c*****
c     Diagnostic printout: vertical profile at one point
c     nz = nint(86400./dtime)
c     if (nstep.gt.nstop-nz .and. jj.eq.norec/2)
c       then
c       zzmax = -1.e20
c       do 8000 ji=1,nlon
c         zz = 0.
c         do 8002 jk=1,nlev
c           if (zconda(ji,jk-1).gt.0.) zz = zz + 1.
c8002     continue
c         if (zz.gt.zzmax) then
c           jim =ji
c           zzmax = zz
c         endif
c8000   continue
c
c      if (jj.eq.14) then
c       jim = 18
cloc   if (jj.eq.1) then                                           !cloc
cloc    jim = 1                                                    !cloc
c       do 8010 ipz=1,nprinz
c         write(iuprinz(ipz),8012) nstep*dtime/86400.,jim,jj
c8012     format(/' day=',f8.3,'  ji=',i3,'  jj=',i3)
c         write(iuprinz(ipz),8014)
c    *     (jk, sigkmh(jk),
c    *          zarea(jim,jk), zw(jim,jk),
c    *          zt(jim,jk),    zq(jim,jk), zrho(jim,jk),
c    *          yh(jim,jk),    yt(jim,jj,jk),
c    *          yq(jim,jj,jk), yrho(jim,jk),
c    *          tplu(jim,jk),  qplu(jim,jk)*.864e5,
c    *          zconda(jim,jk)*.864e5,
c    *          fplumt(jim,jj,jk)/(cpair*ypthic(jim,jj,jk)/gravit),
c    *          fplumq(jim,jj,jk)/(ypthic(jim,jj,jk)/gravit),
c    *      jk=nlev,1,-1)
c8014     format(' jk',4x,'sig',
c    *                  5x,'zarea',8x,'zw',
c    *                  8x,'zt',8x,'zq',6x,'zrho',
c    *                  8x,'yh',8x,'yt',
c    *                  8x,'yq',6x,'yrho',
c    *                  6x,'tplu',6x,'qplu',
c    *                  4x,'zconda',
c    *                  4x,'fplumt',4x,'fplumq'
c    *           /(i3,f7.3,14f10.5))
c8010   continue
c      endif
c     endif
c*****
 
c*****
c     Diagnostic printout: one line for one point
c     nz = nint(86400./dtime)
c     if (nstep.gt.nstop-2*nz .and. nstep.le.nstop-nz.and.jj.eq.norec/2)
c       then
c       ji = nlon/2
cloc  if (jj.eq.1)                                                 !cloc
cloc    ji = 1                                                     !cloc
c       if (mod(nstep,nz).eq.1) write(iuout,8020) ji,jj
c8020   format(/' ji,jj=',2i3)
c       do 8022 jk=1,nlev
c         chaplu(jk) = '   '
c8022   continue
c       chaplu(laytop(ji)) = '***'
c
c       do 8030 ipz=1,nprinz
c         write(iuprinz(ipz),8032)
c8032     format(5x,'day',2x,'heipbl',3x,'ytflx',
c    *           5x,5(6x,'yh',1x),5x,5(6x,'yq',1x) )
c         write(iuprinz(ipz),8034) nstep*dtime/86400.,
c    *      heipbl(ji), ytflx(ji,jj),
c    *      (yh(ji,jk),chaplu(jk)(1:1),
c    *      jk=nlev,nlev-4,-1),
c    *      (yq(ji,jj,jk),chaplu(jk)(1:1), jk=nlev,nlev-4,-1)
c8034     format(f8.3,f8.1,f8.1,5x,5(f8.2,a1),5x,5(f8.5,a1))
c8030   continue
c     endif
c*****
 
c     Rest of diagnostics are global (don't work if multiprocessed)
#ifndef omp
c--------------------------
      if (jj.eq.norec) then
c--------------------------
c*****
c     Diagnostic printout: number of negative q's globally vs level
c     zz = dtime / (3600.)
c     qnegmax = 1.e20
c     nqnegtot = 0
c     do 8040 jk=1,nlev
c       qrat(jk) = 0.
c       qneg(jk) = 0.
c       nqrat(jk) = 0
c       nqneg(jk) = 0
c       do 8042 jz=1,norec
c       do 8042 ji=1,nlon
c         dqnew = zz * fplumq(ji,jz,jk) / (ypthic(ji,jz,jk)/gravit)
c         if (abs(yq(ji,jz,jk)).gt.0.00005) then
c           qrat(jk) = qrat(jk) + abs(dqnew) / abs(yq(ji,jz,jk))
c           nqrat(jk) = nqrat(jk) + 1
c         endif
c         qnew = yq(ji,jz,jk) + dqnew
c         if (qnew.lt.0.) then
c           qneg(jk) = qneg(jk) + qnew
c           nqneg(jk) = nqneg(jk) + 1
c           nqnegtot = nqnegtot + 1
c           if (qnew.lt.qnegmax) then
c             qnegmax = qnew
c             jkmax = jk
c             jimax = ji
c             jjmax = jz
c           endif
c         endif
c8042   continue
c       if (nqrat(jk).gt.0) qrat(jk) = qrat(jk)/nqrat(jk)
c       if( nqneg(jk).gt.0) qneg(jk) = qneg(jk)/nqneg(jk)
c8040 continue
c
c     if (nqnegtot.gt.0) then
c       do 8050 ipz=1,nprinz
c         write(iuprinz(ipz),8052)
c    *       nstep*dtime/86400., zptot, zetot, zqtot,
c    *       qnegmax, jimax,jjmax,jkmax,
c    *       nqrat(jkmax),qrat(jkmax),nqneg(jkmax),qneg(jkmax)
c8052     format(/' day=',f8.3,
c    *            '   zptot=',f10.5,'   zetot=',f10.5,'   zqtot=',f10.5
c    *           /' qnegmax=',f10.5,'  at=',3i4
c    *           /' nqrat=',i4,'   qrat=',f9.6,'  nqneg=',i4,
c    *            '  qneg=',f9.6)
c8050   continue
c     endif
c*****
 
c*****
c     Diagnostic printout: max heimax, numfret
c     nz = nint(86400./dtime)
c     if (nstep.gt.nstop-nz) then
c       zzm = -1.e20
c       do 8060 jz=1,norec
c       do 8060 ji=1,nlon
c         if (heimax(ji,jz).gt.zzm) then
c           jim = ji
c           jzm = jz
c           zzm = heimax(jim,jzm)
c         endif
c8060   continue
c       write(*,8062) nstep*dtime/86400., zzm, jim, jzm,
c    *                (numfret(jk),jk=nlev,1,-1)
c8062   format(' day=',f8.3,'  heimax=',f8.1,'  at ji,jj=',2i3,
c    *         '  numfret=',(18i5))
c     endif
c*****
 
c*****
c     Diagnostic printout: number of pbl and plume tops in each layer
c     (written to fort.70)
c     if (mod(nstep,nint(1.*86400./dtime)) .eq. 0) then
c       iu = 70
c       write(iu,8070) (nstep-nrstrt)*dtime/86400.
c       write(iu,8072) 'pbl  ', (numpbl (jk),jk=nlev,1,-1)
c       write(iu,8072) 'fret ', (numfret(jk),jk=nlev,1,-1)
c       write(iu,8072) 'freb ', (numfreb(jk),jk=nlev,1,-1)
c       ier = flush (iu)
c8070   format('day:',f8.3)
c8072   format(4x,'num',a,'=',(18i5))
c     endif
c*****
c----------
      endif
c----------
#endif
 
      return
      end
c
c ***************************** reevap *********************************
c
      subroutine reevap (yt, yq, ypm, ypthic, ylwc, ycld, plumeprec,
     *                   fplumt, fplumq, frain, fsnow, precflx,
     *                   dmprec, jj, icode
#ifdef isotrac
     *                  ,yq_o, yq_h, ylwc_o, ylwc_h,
     *                   plumepr_o, plumepr_h, fplumq_o, fplumq_h,
     *                   frain_o, frain_h, fsnow_o, fsnow_h, za_i, zh_i,
     *                   rtrq, plumepr_t, rtrfq, rtrpr, rtrps
#endif
     *                  )

c     Does (one latitude strip at a time):
c     (i) local evaporation of cloud liquid water to large-scale vapor,
c     (ii) re-evaporation of falling precipitation, and (iii) local
c     conversion of cloud liquid water to originating precip.
c     Called separately for convective/anvil/stratiform precip, to
c     allow (i) different param values, (ii) diagnostics, and (iii)
c     different stochastic params in lsx for conv vs strat (so need
c     frainc/s and fsnowc/s separately).
c
c     As falling precip enters a layer from above, its sensible heat and
c     phase are adjusted to that layer's temperature and the associated
c     heat flux given to fplumt for that layer. Then a fraction of the
c     precip re-evaporates or sublimates within that layer proportional
c     to the ambient dryness (qsat-q), using a conceptual model with
c     fixed fall velocities, transfer coeffs and drop radii.
 
c     Supplied:
c     yt     = agcm temperatures (deg K)
c     yq     = agcm specific humidity (kg/kg)
c     ypm    = agcm mid-layer pressure (N/m2)
c     ypthic = agcm layer pressure thicknesses (N/m2)
c     ycld   = agcm total cloud fraction (3-D)
c     plumeprec = within-plume precip source (kg/m2/s), only for conv.
c     jj     = latitude index
c     icode  = 1 for convective, 2 for anvil, 3 for stratus
 
c     Modified:
c     ylwc   = cloud liquid water content (conv/anvil/strat) (kg/kg)
c     fplum* = reservoirs of "flux convergence*time" in each agcm layer
c              (J/m2 for fplumt, kg/m2 for fplumq)
c     frain  = rainfall rate at ground (kg/m2/s)
c     fsnow  = snowfall rate at ground (kg/m2/s)
c     precflx= precip flux through interfaces (kg/m2/s)
c     dmprec = precip_at_ground budget term (dmprec[c,a,s])
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMBUD.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMEST.h>
c ---------------------------------------------------------------------
#include <COMFLG.h>
c ---------------------------------------------------------------------
#include <COMINI.h>
c ---------------------------------------------------------------------
#include <COMLUN.h>
c ---------------------------------------------------------------------
#include <COMMAP.h>
c ---------------------------------------------------------------------
#include <COMTIM.h>
c ---------------------------------------------------------------------
      dimension
     *  yt(nlon,norec,nlev),      yq(nlon,norec,nlev),
     *  ypm(nlon,norec,nlev),     ypthic(nlon,norec,nlev),
 
     *  ylwc(nlon,norec,nlev),    ycld(nlon,norec,nlev),
     *  plumeprec(nlon,norec,nlev),
     *  fplumt(nlon,norec,nlev),  fplumq(nlon,norec,nlev),
     *  frain(nlon,norec),        fsnow(nlon,norec),
     *  precflx(nlon,norec,nlevp),
     *  dmprec(norec)
#ifdef isotrac
      dimension
     *  yq_o(nlon,norec,nlev),        yq_h(nlon,norec,nlev),
     *  ylwc_o(nlon,norec,nlev),      ylwc_h(nlon,norec,nlev),
     *  plumepr_o(nlon,norec,nlev),   plumepr_h(nlon,norec,nlev),
     *  fplumq_o(nlon,norec,nlev),    fplumq_h(nlon,norec,nlev),
     *  frain_o(nlon,norec),          frain_h(nlon,norec),
     *  fsnow_o(nlon,norec),          fsnow_h(nlon,norec),
     *  za_i(nlon,nlev),              zh_i(nlon,nlev),
     *  rtrq(nlon,norec,nlev,nrtr),   plumepr_t(nlon,norec,nlev,nrtr),
     *  rtrfq(nlon,norec,nlev,nrtr),  rtrpr(nlon,norec,nrtr),
     *  rtrps(nlon,norec,nrtr)
#endif
 
      dimension
     *  zqsat(nlon),   zdqsat(nlon),  zrho(nlon),
     *  zdz(nlon),     zm(nlon),      zch(nlon),
     *  zfacp(nlon),   zfaco(nlon),   zlath(nlon),
     *  zprec(nlon),   zdprec(nlon),  zcoolp(nlon),
     *  zorig(nlon),   za(nlon),      zb(nlon),
     *  zplump(nlon),  clw(nlon)
#ifdef isotrac
      dimension
     *  zprec_o(nlon),    zprec_h(nlon),
     *  zdprec_o(nlon),   zdprec_h(nlon),
     *  zorig_o(nlon),    zorig_h(nlon),
     *  zplump_o(nlon),   zplump_h(nlon),
     *  zprec_t(nlon,nrtr),zdprec_t(nlon,nrtr)
      PARAMETER (PI = 3.14159265358979)
#endif
c ---------------------------------------------------------------------
c     zevp is a factor in the re-evap fraction of falling precip
c     parameter (zevp=0.0) !turn off
      parameter (zevp=.05)
      save zdiagpa, zdiagpg, kount
      data zdiagpa, zdiagpg, kount /0., 0., 0/
c ---------------------------------------------------------------------
#include <SFEST.h>
c ---------------------------------------------------------------------
#include <PARA_I.h>
#ifdef isotrac
c     xd_i: exponent of the kinetic fractionation term as a function of
c     the drop radius for an uniform population
      xd_i(R)  = min(1.,0.5112+1./(20350.*R+1.905))
#endif
c ---------------------------------------------------------------------
 
c       Set namelist params depending on cloud type
 
      if (icode.eq.1) then
        blw    = blw_c
        flw    = flw_c
        flwmin = .05
      else if (icode.eq.2) then
        blw    = blw_a
        flw    = flw_a
        flwmin = .05
      else if (icode.eq.3) then
        blw    = blw_s
        flw    = flw_s
c       flwmin = .05 
        flwmin = .10 ! 777777
      endif
 
c       Algorithm does one pass from top to bottom layer, one longitude
c       strip at a time.
c         zprec = falling precip (kg/m2/s)
c         zorig = precip originating in current layer (kg/m2/s)
 
      call zero (zprec,nlon)
#ifdef isotrac
      call zero (zprec_o,nlon)
      call zero (zprec_h,nlon)
#endif
#ifdef isotracreg
      call zero (zprec_t,nlon*nrtr)
#endif
 
c=====
      do 100 jk=1,nlev
c=====

c          Set ambient quantities

        do ji=1,nlon
          zm(ji) = ypthic(ji,jj,jk)/gravit
          zqsat(ji) = qstblf (yt(ji,jj,jk),ypm(ji,jj,jk))
          zrho(ji) = ypm(ji,jj,jk)
     *             / (rair*(1.+zvir*yq(ji,jj,jk))*yt(ji,jj,jk))
          zdz(ji) = ypthic(ji,jj,jk) / (gravit*zrho(ji))
        enddo

c          Set cloud water content clw, * by "full" fraction flw 

        call clwcalc (clw, yt(1,jj,jk), zrho, icode)
        do ji=1,nlon
          clw(ji) = clw(ji)*flw
        enddo

c-------------------------------------------------------------------
c       Local evaporation of cloud liquid water to large-scale vapor
c-------------------------------------------------------------------
 
        do 110 ji=1,nlon
 
c         finite-difference:
c         zde = min ( ylwc(ji,jj,jk),
c    *                max (0., (zqsat(ji)-yq(ji,jj,jk))/.01)
c    *                * max (ylwc(ji,jj,jk), .02*clw(ji))
c    *                * dtime/(blw*3600.)
c    *              )
c         analytic:
          zz = max (0.,(zqsat(ji)-yq(ji,jj,jk))/.01) * dtime/(blw*3600.)
          zde = ylwc(ji,jj,jk) * (1. - exp(-zz))
          zde = cvmgt (zde,ylwc(ji,jj,jk),ylwc(ji,jj,jk).ge..02*clw(ji))

#ifdef isotrac
c         This process is taken as isotopicaly non fractionating since
c         it can be seen as the total evaporation of droplets at the
c         periphery of clouds.
          if (ylwc(ji,jj,jk).ne.0.) then
            fplumq_o(ji,jj,jk)=fplumq_o(ji,jj,jk)+zde*zm(ji)*
     *                         ylwc_o(ji,jj,jk)/ylwc(ji,jj,jk)
            fplumq_h(ji,jj,jk)=fplumq_h(ji,jj,jk)+zde*zm(ji)*
     *                         ylwc_h(ji,jj,jk)/ylwc(ji,jj,jk)
            ylwc_o(ji,jj,jk) = ylwc_o(ji,jj,jk)*(1-zde/ylwc(ji,jj,jk))
            ylwc_h(ji,jj,jk) = ylwc_h(ji,jj,jk)*(1-zde/ylwc(ji,jj,jk))
          endif
#endif
#ifdef isotracreg
          if (yq(ji,jj,jk).ne.0.) then
            do it=1,nrtr
              rtrfq(ji,jj,jk,it)= rtrfq(ji,jj,jk,it)+zde*zm(ji)*
     *                            rtrq(ji,jj,jk,it)/yq(ji,jj,jk)
            enddo
          endif
#endif
 
c         For now, all cloud water is liquid with no specific heat
          zl = latvap + cpwv*(yt(ji,jj,jk)-tmelt)
          ylwc(ji,jj,jk) = ylwc(ji,jj,jk) - zde
          fplumq(ji,jj,jk) = fplumq(ji,jj,jk) +    zde*zm(ji)
          fplumt(ji,jj,jk) = fplumt(ji,jj,jk) - zl*zde*zm(ji)
#ifdef isotrac
          if (ylwc(ji,jj,jk).eq.0.) then
            ylwc_o(ji,jj,jk) = 0.
            ylwc_h(ji,jj,jk) = 0.
          endif
#endif
 
          dmcloue(jj) = dmcloue(jj) +    zde*zm(ji)/nlon
          dhcloud(jj) = dhcloud(jj) - zl*zde*zm(ji)/nlon
  110   continue
 
c----------------------------------------------------------------------
c       Temperature-adjust and evaporation of precip falling from above
c----------------------------------------------------------------------
 
c       If any precip falling into this layer from above (for any
c       longitude), (i) adjust it to ambient temperature and (ii)
c       evaporate fraction zfacp. Give associated heat and vapor fluxes
c       from both processes to the large-scale reservoirs fplum[t,q],
c       and increment budget diagnostics. First set:
c         zlath = sensible+latent heat for reevap within layer
c                 (as if all phase changes occur at tmelt as in LSX)
c         zfacp = fraction of falling precip reevaped in this layer
c         zch   = sensible+latent heat for adjust to layer temperature
c                 (as if all phase changes occur at tmelt as in LSX)
 
        zprecsum = ssum(nlon,zprec,1)
 
c-----
        if (zprecsum.ne.0.) then
c-----
          do 120 ji=1,nlon
 
            zlath(ji) = cvmgt ( latvap + cpwv*(yt(ji,jj,jk)-tmelt)
     *                                 - ch2o*(yt(ji,jj,jk)-tmelt),
     *                          latsub + cpwv*(yt(ji,jj,jk)-tmelt)
     *                                 - cice*(yt(ji,jj,jk)-tmelt),
     *                          yt(ji,jj,jk).ge.tmelt )
            zdqsat(ji) = zqsat(ji) - yq(ji,jj,jk)
            zz = zevp * zdz(ji) * max(zdqsat(ji),0.)
            zfacp(ji) = 1.-exp(-zz)
c           Fudge to make reevap happen more in lower layers:
c           zfacp(ji) = zfacp(ji) * max (0., min (1., (sig(jk)-.7)/.3))
 
            za(ji) = cvmgt (0.5, -0.5, yt(ji,jj,jk-1).ge.tmelt)
            zb(ji) = cvmgt (0.5, -0.5, yt(ji,jj,jk) .ge.tmelt)
            zch(ji) = - cvmgt ( ch2o, cice, yt(ji,jj,jk-1).ge.tmelt )
     *                        * (yt(ji,jj,jk-1)-tmelt)
     *                + cvmgt ( ch2o, cice, yt(ji,jj,jk) .ge.tmelt )
     *                        * (yt(ji,jj,jk)-tmelt)
     *                + latice* (zb(ji)-za(ji))
 
            zdprec(ji) = zprec(ji)*zfacp(ji)

c+++++++++++++
#ifdef isotrac
c+++++++++++++
c        If the temperature is below the transition temperature defined
c        for reevap then the process is the sublimation of solid, which
c        is isotopicaly non fractionating. Liquid precip will exchange
c        with the vapour in saturated atmospheres (h>0.99) and evaporate
c        in unsaturated ones. Two different formulations based on
c        Stewart (1975) are used for the calculation of the respective
c        isotopic changes. The kinetic fractionation is parametrized by
c        the drop radius of an uniform population, itself deduced from
c        the precipitation rate.
 
c        For convective precip only: A coefficient, fa_i, allows to
c        consider that the precip falls on an area proportional to the
c        plume area, which can vary with height, computed in PLUME, or
c        on all the grid box. Another coeff., fh_i, determines if the
c        relative humidity is taken as for the grid box, as in the
c        plume, saturated, or any intermediate profile.
 
           if (zdprec(ji).eq.0. .or. (yq(ji,jj,jk).eq.0.)) then
             zdprec_o(ji) = 0.
             zdprec_h(ji) = 0.
#ifdef isotracreg
             do it=1,nrtr
               zdprec_t(ji,it) = 0.
             enddo
#endif
           else
c        Solid
             if (yt(ji,jj,jk-1).lt.t_r_i) then
               zdprec_o(ji) = zdprec(ji)*zprec_o(ji)/zprec(ji)
               zdprec_h(ji) = zdprec(ji)*zprec_h(ji)/zprec(ji)
#ifdef isotracreg
               do it=1,nrtr
                 zdprec_t(ji,it) = zdprec(ji)*zprec_t(ji,it)/zprec(ji)
               enddo
#endif
c        Liquid
             else
c          Effective relative humidity and mean drop radius
               hr_i = yq(ji,jj,jk)/zqsat(ji)
c            For convective clouds:
               if (icode.eq.1) then
                 hr_i = max(0.01, min(1., hi_c_i+hs_c_i*
     *                  (cvmgt(zh_i(ji,jk), hr_i, za_i(ji,jk).gt.0.))))
                 if (r_c_i.gt.0.) then
                   z_i  = cvmgt(za_i(ji,jk), 1., za_i(ji,jk).gt.0.)
                   rd_i = min(0.003, (zprec(ji)*sqrt(zrho(ji))/
     *                                    5.85e11/r_c_i/z_i)**.22222)
                 else
                   rd_i = -r_c_i
                 endif
c            For stratiform and anvil clouds:
               else
                 hr_i = max(0.01, min(1., hi_s_i+hs_s_i*hr_i))
                 if (r_s_i.gt.0.) then
                   rd_i = min(0.003, (zprec(ji)*sqrt(zrho(ji))/
     *                                    5.85e11/r_s_i)**.22222)
                 else
                   rd_i = -r_s_i
                 endif
               endif
c          Exchange
               if (hr_i.gt.0.99) then
                 vd_i = 2.13*sqrt(rhoh2o*gravit*rd_i/zrho(ji))
                 ff_i = 1+fr_i(rd_i)*sqrt(vd_i*rd_i/2/pi/
     *                              dv_i(yt(ji,jj,jk),ypm(ji,jj,jk)))
                 cte_i= 3.*zrho(ji)*dv_i(yt(ji,jj,jk),ypm(ji,jj,jk))*
     *                         ff_i*zqsat(ji)/(rd_i*rd_i*rhoh2o*vd_i)
                 cte_o=cte_i/(avl_o(yt(ji,jj,jk))*dr_o**xd_i(rd_i))
                 cte_h=cte_i/(avl_h(yt(ji,jj,jk))*dr_h**xd_i(rd_i))
                 zdprec_o(ji)=zprec_o(ji)-(zprec(ji)-zdprec(ji))*(
     *                  avl_o(yt(ji,jj,jk))*yq_o(ji,jj,jk)/yq(ji,jj,jk)
     *                  +(zprec_o(ji)/zprec(ji)-avl_o(yt(ji,jj,jk))*
     *                yq_o(ji,jj,jk)/yq(ji,jj,jk))*exp(-cte_o*zdz(ji)))
                 zdprec_h(ji)=zprec_h(ji)-(zprec(ji)-zdprec(ji))*(
     *                  avl_h(yt(ji,jj,jk))*yq_h(ji,jj,jk)/yq(ji,jj,jk)
     *                  +(zprec_h(ji)/zprec(ji)-avl_h(yt(ji,jj,jk))*
     *                yq_h(ji,jj,jk)/yq(ji,jj,jk))*exp(-cte_h*zdz(ji)))
#ifdef isotracreg
                 do it=1,nrtr
                   zdprec_t(ji,it)=zprec_t(ji,it)-(zprec(ji)-zdprec(ji))
     *                  *(rtrq(ji,jj,jk,it)/yq(ji,jj,jk)+(zprec_t(ji,it)
     *                  /zprec(ji)-rtrq(ji,jj,jk,it)/yq(ji,jj,jk))
     *                  *exp(-cte_i*zdz(ji)))
                 enddo
#endif
c          Evaporation
               else
                 cte_o= max(0.001,
     *                  1-(1-hr_i)*avl_o(yt(ji,jj,jk))*dr_o**xd_i(rd_i))
                 cte_h= max(0.001,
     *                  1-(1-hr_i)*avl_h(yt(ji,jj,jk))*dr_h**xd_i(rd_i))
                 zdprec_o(ji)=zprec_o(ji)-(zprec(ji)-zdprec(ji))*(avl_o
     *                  (yt(ji,jj,jk))*hr_i*yq_o(ji,jj,jk)/yq(ji,jj,jk)+
     *                  (cte_o*zprec_o(ji)/zprec(ji)-avl_o(yt(ji,jj,jk))
     *                  *hr_i*yq_o(ji,jj,jk)/yq(ji,jj,jk))*(1-zdprec(ji)
     *                  /zprec(ji))**(cte_o/(1-cte_o)))/cte_o
                 zdprec_h(ji)=zprec_h(ji)-(zprec(ji)-zdprec(ji))*(avl_h
     *                  (yt(ji,jj,jk))*hr_i*yq_h(ji,jj,jk)/yq(ji,jj,jk)+
     *                  (cte_h*zprec_h(ji)/zprec(ji)-avl_h(yt(ji,jj,jk))
     *                  *hr_i*yq_h(ji,jj,jk)/yq(ji,jj,jk))*(1-zdprec(ji)
     *                  /zprec(ji))**(cte_h/(1-cte_h)))/cte_h
#ifdef isotracreg
                 do it=1,nrtr
                   zdprec_t(ji,it)=zprec_t(ji,it)-(zprec(ji)-zdprec(ji))
     *                  *(rtrq(ji,jj,jk,it)/yq(ji,jj,jk)+(zprec_t(ji,it)
     *                  /zprec(ji)-rtrq(ji,jj,jk,it)/yq(ji,jj,jk))
     *                  *(1-zdprec(ji)/zprec(ji))**(hr_i/(1-hr_i)))
                 enddo
#endif
               endif
             endif
           endif
 
           fplumq_o(ji,jj,jk)=fplumq_o(ji,jj,jk)+zdprec_o(ji)*dtime
           fplumq_h(ji,jj,jk)=fplumq_h(ji,jj,jk)+zdprec_h(ji)*dtime
           zprec_o(ji)=zprec_o(ji)-zdprec_o(ji)
           zprec_h(ji)=zprec_h(ji)-zdprec_h(ji)
#ifdef isotracreg
           do it=1,nrtr
             rtrfq(ji,jj,jk,it) = rtrfq(ji,jj,jk,it)
     *                                      +zdprec_t(ji,it)*dtime
             zprec_t(ji,it) = zprec_t(ji,it)-zdprec_t(ji,it)
           enddo
#endif
c+++++
#endif
c+++++

            zcoolp(ji) = zprec(ji)*zch(ji) + zdprec(ji)*zlath(ji)
            fplumq(ji,jj,jk) = fplumq(ji,jj,jk) + zdprec(ji)*dtime
            fplumt(ji,jj,jk) = fplumt(ji,jj,jk) - zcoolp(ji)*dtime
            zprec(ji) = zprec(ji) - zdprec(ji)
  120     continue
 
          dmreev(jj)= dmreev(jj) + ssum(nlon,zdprec,1)*dtime/nlon
          dhreev(jj)= dhreev(jj) - ssum(nlon,zcoolp,1)*dtime/nlon
c-----
        endif
c-----
 
c------------------------------------------------------------------
c       Conversion of large-scale cloud water to originating precip
c------------------------------------------------------------------
 
c          If convective, copy within-plume precip from 3-D plumeprec
 
        if (icode.eq.1) then
           call scopy (nlon, plumeprec(1,jj,jk), 1, zplump, 1)
#ifdef isotrac
           call scopy (nlon,plumepr_o(1,jj,jk),1,zplump_o,1)
           call scopy (nlon,plumepr_h(1,jj,jk),1,zplump_h,1)
#endif
        else
           call zero (zplump, nlon)
#ifdef isotrac
           call zero (zplump_o, nlon)
           call zero (zplump_h, nlon)
#endif
        endif
 
        do 130 ji=1,nlon
 
c            Autoconversion
 
          zliq  = max (0., min (1., (yt(ji,jj,jk)-tmelt-tliqa)/tliqb ))
          zland = yoroland(ji,jj)
          zdlw =     zland *    zliq *dlw_ll*dlwfac
     *         + (1.-zland)*    zliq *dlw_ol*dlwfac
     *         +     zland *(1.-zliq)*dlw_li
     *         + (1.-zland)*(1.-zliq)*dlw_oi
          zdlw = zdlw*dlwfac
 
c         Dependence on ambient density (geometric droplet density
c         decreases, but fall speed increases, as zrho decreases):
c         zdlw  = zdlw *      zrho(ji)
          zdlw  = zdlw * sqrt(zrho(ji))
c         zdlw  = zdlw / sqrt(zrho(ji))
c         zdlw  = zdlw /      zrho(ji)
 
c----
c         Sundqvist, finite-difference:
c         zdpa = ylwc(ji,jj,jk)
c    *         * min ( 1.,
c    *                 (1. - exp(-(ylwc(ji,jj,jk)/clw(ji))**2) )
c    *                 * dtime/(zdlw*3600.)
c    *               )
c----
c         linear, quadratic or cubic, finite-difference:
c         zdpa = ylwc(ji,jj,jk)
clin *         * min (1.0,                            dtime/(zdlw*3600))
cquad*         * min (1.0,(ylwc(ji,jj,jk)/clw(ji))   *dtime/(zdlw*3600))
ccubi*         * min (1.0,(ylwc(ji,jj,jk)/clw(ji))**2*dtime/(zdlw*3600))

c         linear, quadratic or cubic, analytic:
clin      zz =      ( 1.+                             dtime/(zdlw*3600))
cquad     zz =      ( 1.+ (ylwc(ji,jj,jk)/clw(ji))   *dtime/(zdlw*3600))
ccubi     zz = sqrt(1.+2.*(ylwc(ji,jj,jk)/clw(ji))**2*dtime/(zdlw*3600))
          zz =      ( 1.+ (ylwc(ji,jj,jk)/clw(ji))   *dtime/(zdlw*3600))
          zdpa = ylwc(ji,jj,jk) * (1. -  1./zz)
c----
 
c         autoconversion stops for ylwc below some small amount:
          zdpa = min (zdpa, max (0., ylwc(ji,jj,jk)-clw(ji)*flwmin))
#ifdef isotrac
c         autoconversion is isotopically non-fractionating:
          zdpa_o= 0.
          zdpa_h= 0.
          zdpb_o= 0.
          zdpb_h= 0.
          if (ylwc(ji,jj,jk).ne.0.) then
            zdpa_o=zdpa*ylwc_o(ji,jj,jk)/ylwc(ji,jj,jk)
            zdpa_h=zdpa*ylwc_h(ji,jj,jk)/ylwc(ji,jj,jk)
            ylwc_o(ji,jj,jk)=ylwc_o(ji,jj,jk)-zdpa_o
            ylwc_h(ji,jj,jk)=ylwc_h(ji,jj,jk)-zdpa_h
          endif
#endif
 
          ylwc(ji,jj,jk) = ylwc(ji,jj,jk) - zdpa
 
c            Accretion by falling precip
 
c         finite-difference:
c         zdpb = ylwc(ji,jj,jk)
c    *           * min (1.0, (zprec(ji)*86400./100.)*dtime/(glw*3600.) )
c         analytic:
          zz = (zprec(ji)*86400./100.) * dtime/(glw*3600.)
          zdpb = ylwc(ji,jj,jk) * (1. - exp(-zz))
#ifdef isotrac
c         accretion is isotopically non-fractionating:
          if (ylwc(ji,jj,jk).ne.0.) then
            zdpb_o=zdpb/ylwc(ji,jj,jk)*ylwc_o(ji,jj,jk)
            zdpb_h=zdpb/ylwc(ji,jj,jk)*ylwc_h(ji,jj,jk)
            ylwc_o(ji,jj,jk)=ylwc_o(ji,jj,jk)-zdpb_o
            ylwc_h(ji,jj,jk)=ylwc_h(ji,jj,jk)-zdpb_h
          endif
#endif
 
          ylwc(ji,jj,jk) = ylwc(ji,jj,jk) - zdpb
#ifdef isotrac
          if (ylwc(ji,jj,jk).eq.0.) then
            ylwc_o(ji,jj,jk) = 0.
            ylwc_h(ji,jj,jk) = 0.
          endif
#endif
 
c            Transfer converted+accreted cloud water to originating
c            precip (plus within-plume convective precip source),
c            add to zprec for next layer below, give associated heat
c            fluxes to the large-scale reservoir fplumt, and increment
c            budget diagnostics.
 
          zorig(ji) = (zdpa+zdpb)*zm(ji)/dtime + zplump(ji)
          zprec(ji) = zprec(ji) + zorig(ji)
#ifdef isotrac
          zorig_o(ji)=(zdpa_o+zdpb_o)*zm(ji)/dtime + zplump_o(ji)
          zorig_h(ji)=(zdpa_h+zdpb_h)*zm(ji)/dtime + zplump_h(ji)
          zprec_o(ji)=zprec_o(ji)+zorig_o(ji)
          zprec_h(ji)=zprec_h(ji)+zorig_h(ji)
#endif
#ifdef isotracreg
          if (yq(ji,jj,jk).ne.0.) then
            do it=1,nrtr
              zorig_t=(zdpa+zdpb)*rtrq(ji,jj,jk,it)/yq(ji,jj,jk)*zm(ji)/
     *               dtime + cvmgt(plumepr_t(ji,jj,jk,it),0.,icode.eq.1)
              zprec_t(ji,it)=zprec_t(ji,it)+zorig_t
            enddo
          endif
#endif
 
c            Accumulate 3-D precip flux through interfaces (kg/m2/s)
c            (passed to cosmogenic isotope subr cosmog in phys)
#ifdef cosmog
          precflx(ji,jj,jk+1) = precflx(ji,jj,jk+1) + zprec(ji)
#endif
 
c         For now, all cloud water is liquid with no specific heat
          zl = cvmgt (        - ch2o*(yt(ji,jj,jk)-tmelt),
     *                 latice - cice*(yt(ji,jj,jk)-tmelt),
     *                 yt(ji,jj,jk).ge.tmelt )
          fplumt(ji,jj,jk) = fplumt(ji,jj,jk) + zl*zorig(ji)*dtime
 
          dmcloup(jj) = dmcloup(jj) + zplump(ji)  *dtime/nlon
          dmclouc(jj) = dmclouc(jj) + zorig(ji)   *dtime/nlon
          dmcloua(jj) = dmcloua(jj) + zdpb*zm(ji)       /nlon
          dhcloud(jj) = dhcloud(jj) + zl*zorig(ji)*dtime/nlon
          zdiagpa     = zdiagpa     +    zorig(ji)*cosbud(jj)/nlon

#ifdef isotrac
c         liquid-vapor isotopic equilibrium (optional) forcing:
          if ((equ_i.eq.1).and.(yt(ji,jj,jk).ge.t_r_i)) then
            q_o = yq_o(ji,jj,jk)
            q_h = yq_h(ji,jj,jk)
            call equil_i (1.,           yt(ji,jj,jk),    ypm(ji,jj,jk),
     *                yq(ji,jj,jk),    q_o,              q_h,
     *              ylwc(ji,jj,jk), ylwc_o(ji,jj,jk), ylwc_h(ji,jj,jk))
            fplumq_o(ji,jj,jk) = fplumq_o(ji,jj,jk)+
     *                           (q_o-yq_o(ji,jj,jk))*zm(ji)
            fplumq_h(ji,jj,jk) = fplumq_h(ji,jj,jk)+
     *                           (q_h-yq_h(ji,jj,jk))*zm(ji)
          endif
#endif
 
  130   continue
 
c=====
  100 continue
c=====
 
c------------------------------------------------------------
c        Turbulent deposition of cloud particles onto surface
c        (only for lowest layer)
c------------------------------------------------------------
 
      do 300 ji=1,nlon
        zm(ji) = ypthic(ji,jj,nlev)/gravit
 
c       Don't need to recompute zrho,zdz, since last calculated
c       in loop 100 above for jk=nlev(!)
c       zrho(ji) = ypm(ji,jj,nlev)
c    *           / (rair*(1.+zvir*yq(ji,jj,nlev))*yt(ji,jj,nlev))
c       zdz(ji) = ypthic(ji,jj,nlev) / (gravit*zrho(ji))

c       finite difference:
c       zds = ylwc(ji,jj,nlev)
c    *      * min ( 1., (100./zdz(ji))*dtime/(elw*3600.) )
c       analytic:
        zz = (100./zdz(ji))*dtime/(elw*3600.)
        zds = ylwc(ji,jj,nlev) * (1. - exp(-zz))

#ifdef isotrac
c       turbulent deposition is isotopically non-fractionating:
        if (ylwc(ji,jj,nlev).ne.0.) then
          zds_o = zds*ylwc_o(ji,jj,nlev)/ylwc(ji,jj,nlev)
          zds_h = zds*ylwc_h(ji,jj,nlev)/ylwc(ji,jj,nlev)
          ylwc_o(ji,jj,nlev) = ylwc_o(ji,jj,nlev)-zds_o
          ylwc_h(ji,jj,nlev) = ylwc_h(ji,jj,nlev)-zds_h
          zprec_o(ji) = zprec_o(ji)+zds_o*zm(ji)/dtime
          zprec_h(ji) = zprec_h(ji)+zds_h*zm(ji)/dtime
        endif
#endif
#ifdef isotracreg
        if (yq(ji,jj,nlev).ne.0.) then
          do it=1,nrtr
            zds_t = zds*rtrq(ji,jj,nlev,it)/yq(ji,jj,nlev)
            zprec_t(ji,it) = zprec_t(ji,it)+zds_t*zm(ji)/dtime
          enddo
        endif
#endif
 
c       For now, all cloud water is liquid with no specific heat
        zl = cvmgt (        - ch2o*(yt(ji,jj,nlev)-tmelt),
     *               latice - cice*(yt(ji,jj,nlev)-tmelt),
     *               yt(ji,jj,nlev).ge.tmelt )
 
        ylwc(ji,jj,nlev) = ylwc(ji,jj,nlev) - zds

#ifdef isotrac
        if (ylwc(ji,jj,nlev).eq.0.) then
          ylwc_o(ji,jj,nlev) = 0.
          ylwc_h(ji,jj,nlev) = 0.
        endif
c       liquid-vapor isotopic equilibrium (optional) forcing:
        if ((equ_i.eq.1).and.(yt(ji,jj,nlev).ge.t_r_i)) call equil_i
     *                     (1.,     yt(ji,jj,nlev),    ypm(ji,jj,nlev),
     *          yq(ji,jj,nlev),   yq_o(ji,jj,nlev),   yq_h(ji,jj,nlev),
     *        ylwc(ji,jj,nlev), ylwc_o(ji,jj,nlev), ylwc_h(ji,jj,nlev))
#endif

        fplumt(ji,jj,nlev) = fplumt(ji,jj,nlev) + zl*zds*zm(ji)
        zprec(ji) = zprec(ji) + zds*zm(ji)/dtime
 
c          Accumulate 3-D precip flux through interfaces (kg/m2/s)
c          (passed to cosmogenic isotope subr cosmog in phys)
#ifdef cosmog
        precflx(ji,jj,nlev+1) = precflx(ji,jj,nlev+1) + zprec(ji)
#endif
 
        dmclous(jj) = dmclous(jj) +    zds*zm(ji)/nlon
        dhcloud(jj) = dhcloud(jj) + zl*zds*zm(ji)/nlon
  300 continue
 
c        At this point zprec = precip falling out of bottom of lowest
c        layer, with temperature of lowest layer and phase depending
c        on that temp >= or < tmelt. Increment model precip quantities
c        (frain,fsnow) and budget diagnostics.
 
      do 400 ji=1,nlon
        frain(ji,jj) = cvmgt ( frain(ji,jj)+zprec(ji), frain(ji,jj),
     *                         yt(ji,jj,nlev).ge.tmelt )
        fsnow(ji,jj) = cvmgt ( fsnow(ji,jj)+zprec(ji), fsnow(ji,jj),
     *                         yt(ji,jj,nlev).lt.tmelt )
#ifdef isotrac
        frain_o(ji,jj)=cvmgt (frain_o(ji,jj)+zprec_o(ji),frain_o(ji,jj),
     *                         yt(ji,jj,nlev).ge.tmelt )
        frain_h(ji,jj)=cvmgt (frain_h(ji,jj)+zprec_h(ji),frain_h(ji,jj),
     *                         yt(ji,jj,nlev).ge.tmelt )
        fsnow_o(ji,jj)=cvmgt (fsnow_o(ji,jj)+zprec_o(ji),fsnow_o(ji,jj),
     *                         yt(ji,jj,nlev).lt.tmelt )
        fsnow_h(ji,jj)=cvmgt (fsnow_h(ji,jj)+zprec_h(ji),fsnow_h(ji,jj),
     *                         yt(ji,jj,nlev).lt.tmelt )
#endif
#ifdef isotracreg
        do it=1,nrtr
          rtrpr(ji,jj,it)=cvmgt (rtrpr(ji,jj,it)+zprec_t(ji,it),
     *                       rtrpr(ji,jj,it), yt(ji,jj,nlev).ge.tmelt )
          rtrps(ji,jj,it)=cvmgt (rtrps(ji,jj,it)+zprec_t(ji,it),
     *                       rtrps(ji,jj,it), yt(ji,jj,nlev).lt.tmelt )
        enddo
#endif
        dmprec(jj)  = dmprec(jj)  + zprec(ji)*dtime/nlon
        zdiagpg     = zdiagpg     + zprec(ji)*(cosbud(jj)/nlon)
  400 continue
 
c     Write diagnostics (doesn't work for multiprocessing)
#ifndef omp
c     kount = kount + 1
c     if (kount.eq.3*norec) then
c       if ( mod(nstep,nint(1.*86400./dtime)) .eq. 0 ) then
c         write(*,*)
c         zz = 8.64e7/rhoh2o
c         write(*,500) nstep*dtime/86400.,
c    *      zz*zdiagpa, zz*zdiagpg, zdiagpg/max(zdiagpa,1.e-10)
c 500     format(' reevap:  day=',f9.3,'  source=',f8.3,
c    *           '  precip=',f8.3,'  p/s=',f6.3)
c         zdiagpa = 0.
c         zdiagpg = 0.
c         kount = 0
c       endif
c     endif
#endif
 
      return
      end
c
c ***************************** clwcalc ********************************
c
      subroutine clwcalc (clw, ta, rho, icode)

c     Sets cloud water content (kg/kg), given air temperature
c     ta (deg K) and density rho (kg/m3), for one longitude circle.
c     icode = 1,2,3 for convective, anvil, cirrus.
c     Needs COMCON.h (tmelt) and  COMFLG (clwfac_[c,a,s]).
c     Multiply by 1.e-3 for g/m3 to kg/m3, and by 1/density for 
c     kg/m3 to kg/kg.
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMFLG.h>
c ---------------------------------------------------------------------

      dimension clw(nlon), ta(nlon), rho(nlon)

      if (icode.eq.1) then
        clwfac = clwfac_c
      else if (icode.eq.2) then
        clwfac = clwfac_a
      else if (icode.eq.3) then
        clwfac = clwfac_s
      endif

c        Lemus et al.,1997,J.Clim.,10, 52-64, Eq. (1) + max:
c
c     do i=1,nlon
c       clw(i) = max (   0.127 
c    *                   + 6.78e-3*(ta(i)-tmelt) 
c    *                   + 1.29e-4*(ta(i)-tmelt)**2 
c    *                   + 8.68e-7*(ta(i)-tmelt)**3,  1.e-3 )
c    *           * 1.e-3 * clwfac / rho(i)
c     enddo

c        Lemus et al.,1997, parabolic fit to Fig. 1, log(CLW) vs T,
c        using parabsolve.f:
c
c     clw7:
c     do i=1,nlon
c       clw(i) = 10. ** ( - 0.910457 
c    *                    + 0.0161058   * (ta(i)-tmelt) 
c    *                    - 0.000360577 * (ta(i)-tmelt)**2 )
c    *           * 1.e-3 * clwfac / rho(i)
c     enddo

c        Lemus et al.,1997, 2-piece linear fit to Fig. 1, log(CLW) vs T,
c        (and constant/rho for convective): 
c
c     if (icode.eq.1) then                                  ! convective
c       do i=1,nlon
c         clw(i) = .10 * 1.e-3 * clwfac / rho(i)
c       enddo
c     else
c       do i=1,nlon
c         zt = ta(i)-tmelt
c         zlog = -1. + 0.5*max (0., (zt +15.)/(20. +15.))
c    *               - 2.0*max (0., (-15.-zt)/(-15.+60.))  ! lemus 
c    *               - 2.0*max (0., (-15.-zt)/(-15.+40.))  ! a17a
c    *               - 3.0*max (0., (-15.-zt)/(-15.+40.))  ! a16,a17
c         clw(i) = (10.**zlog) * 1.e-3 * clwfac / rho(i)
c       enddo
c     endif
        
      if (clwfac.lt.0.) then
        do i=1,nlon
c         Lemus et al.,1997,J.Clim.,10, 52-64, Eq. (1) + max:
          clw(i) = max (   0.127 
     *                     + 6.78e-3*(ta(i)-tmelt) 
     *                     + 1.29e-4*(ta(i)-tmelt)**2 
     *                     + 8.68e-7*(ta(i)-tmelt)**3,  1.e-3 )
     *             * 1.e-3 * abs(clwfac) / rho(i)
        enddo
      else
        do i=1,nlon
          clw(i) = clwfac 
        enddo
      endif

      return
      end
 
c=============
#ifdef isotrac
c=============

c ***************************** fr_i **********************************
 
      real function fr_i(r)

c   Ventilation factor as a function of the drop radius, fitted from
c   Beard & Pruppacher and Kinzer & Gunn experimental results.

      A = exp(-1000.*r)
      if (r.lt.0.000346) then
        fr_i = max (0., 52.313-199.45*A+257.73*A*A-110.61*A*A*A)
      else
        fr_i =          2.0198-6.5565*A+12.537*A*A-7.3916*A*A*A
      endif

      return
      end

c ***************************** equil_i *******************************
 
      subroutine equil_i (h, T, p, qv, qv_o, qv_h, ql, ql_o, ql_h)
 
c   Equil_i does the reequilibration of the cloud condensed phase with
c   the vapour, following Jouzel et al, 1975. The calculations consider
c   that only small droplets are present (r < 0.03 mm), which allows to
c   neglect the kinetic component of the fractionation, and to take into
c   account pseudo-equilibrium in an unsaturated atmosphere.
c       Supplied:
c        - h, T, p   : Relative humidity, temperature and pressure
c        - qv, ql    : vapour and liquid phase specific humidities
c       Modified:
c        - qv_*, ql_*: vapour and liquid phase isotopic specific
c                       humidities (18O, 2H)
 
c   For now, equil_i is called from Stabco, Plume, and optionaly from
c   Reevap, with (h=1.), and so apply simply equilibrium fractionation.
c   Code lines relevant to unsaturated air are commented.
c---------
#include <COMFLG.h>
#include <PARA_I.h>
c---------
c--
c needs comsat.h, comcon.h?
c       cth    = 7.45e-5*(T-273.16)+0.0241885
c       rvs    = 0.018*esat(T)/8.31/T
c       A      = 1/(1+dv_i(T,p)*rvs*2.5e6/cth/T*(2.5e6*0.018/8.31/T-1))
c--
        ak_o = avl_o(T)
c    *           *abs(h)/(abs(h)+A*(1-abs(h))*(1-avl_o(T)/dr_o))
        ak_h = avl_h(T)
c    *           *abs(h)/(abs(h)+A*(1-abs(h))*(1-avl_h(T)/dr_h))
        if (qv.gt.1.e-20) then 
          qv_o   = (qv_o+ql_o)/(1+ak_o*ql/qv)
          qv_h   = (qv_h+ql_h)/(1+ak_h*ql/qv)
          ql_o   = qv_o*ak_o*ql/qv
          ql_h   = qv_h*ak_h*ql/qv
        else
          qv_o   = 0.                             ! dp 6/04
          qv_h   = 0.                             ! dp 6/04
        endif
 
      return
      end

c ****************************** prm_i ********************************
 
      subroutine prm_i (s, tli, t, ta, qv, qva, qv_o, qv_h,
     *                                 ql_in, qla_in, ql_o, ql_h) ! dp
cdp  *                                 ql,    qla,    ql_o, ql_h) ! dp
 
c   Prm_i does the isotopic fractionation for condensation following
c   the model of Merlivat and Jouzel, 1979. The formation of liquid
c   (precip or cloud) is made at the equilibrium while the formation
c   of solid involves a kinetic fractionation parametrized by a
c   temperature dependent supersaturation function. The model iterates
c   several times if the relative variation of the vapour specific
c   humidity exceeds 15%.
c       Supplied:
c        - s      : slope of the supersaturation function
c        - tli    : liquid/solid transition temperature
c        - t, ta  : initial and final temperatures
c        - qv, qva: initial and final vapour phase specific humidities
c        - ql, qla: initial and final liquid/solid phase specific humid.
c       Modified:
c        - qv_*   : vapour phase isotopic specific humidities (18O, 2H)
c        - ql_*   : liquid/solid phase isotopic specif. humid. (18O, 2H)
 
c---------
#include <COMCON.h>
#include <PARA_I.h>
c---------

      ql = ql_in      ! dp
      qla = qla_in    ! dp

      tt  = t
      qvt = qv
      qlt = ql
      n   = 1
      if (qv.gt.0.0) n = 1+int(abs((qva-qv)/qv/0.15))
      dt  = (ta-t)/n
      dqv = (qva-qv)/n
      dql = (qla-ql)/n
                          do 1 i=1,n
      ta  = t+dt
      qva = qv+dqv
      qla = ql+dql
      if ((t.gt.tli).and.(ta.ge.tli)) then
         al_o = avl_o(t)
         al_h = avl_h(t)
         ala_o= avl_o(ta)
         ala_h= avl_h(ta)
      else
         ss   = 1.-s*(t -273.15)
         ssa  = 1.-s*(ta-273.15)
         al_o = avi_o(t) *ss /((ss -1.)*avi_o(t) *dr_o+1.)
         al_h = avi_h(t) *ss /((ss -1.)*avi_h(t) *dr_h+1.)
         ala_o= avi_o(ta)*ssa/((ssa-1.)*avi_o(ta)*dr_o+1.)
         ala_h= avi_h(ta)*ssa/((ssa-1.)*avi_h(ta)*dr_h+1.)
      endif
 
      alm_o=(ala_o+al_o)/2.
      alm_h=(ala_h+al_h)/2.
      qvm=(qva+qv)/2.
      qlm=(qla+ql)/2.
      if (qvm.ne.0.) then
         vk_o=(alm_o*(qva-qv)*(1+qlm/qvm)-qlm*(ala_o-al_o))
     *       /(qvm+alm_o*qlm)
         vk_h=(alm_h*(qva-qv)*(1+qlm/qvm)-qlm*(ala_h-al_h))
     *       /(qvm+alm_h*qlm)
         qv_o=qv_o*(1+vk_o/(1-vk_o/2.))
         qv_h=qv_h*(1+vk_h/(1-vk_h/2.))
         ql_o=ala_o*qv_o*qla/qva
         ql_h=ala_h*qv_h*qla/qva
      endif
      t  = ta
      qv = qva
      ql = qla
 1                        continue
 
      t   = tt
      qv  = qvt
      ql  = qlt
      return
      end
c===== 
#endif
c===== 
