c -- $Id: surfctl.F,v 1.2 2001/07/03 21:29:30 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/surfctl.F,v $
c -- $Name:  $

c--------
c surfctl
c--------
c
c----------------------------------------------------------------------
c
      SUBROUTINE SURFCTL
     *    ( NSTEP,    DTGCM,    FRADSW,
     *      ZENGCM,
     *      SRGCMD,   SRGCMI,  
     *      FIRABS,   UGCM,     VGCM,     TGCM,    QGCM,
     *      RAICGCM,  SNOCGCM,  RAISGCM, SNOSGCM, STOCHPR,PGCM,  CLDGCM,
     *      FIREMS,   FSEGCM,   FVAGCM,  TAXGCM,  TAYGCM, SIGCMA,SIGCMB,
     *      OMEGA,    SEQSST, 
     *      DYNAMICE, DYNAMHEM,
     *      QFACTOR, QNORWEG, QICES,  QICEN, QLIMIT,
     *      CO2PPM,   HOTLIMIT,
     *      NSREST,   NSREST_IBIS,VERSRESIN,FSAVE,   FNEXF,   ISTEPY,
     *      NRSTRT,   NSTOP,
     *      MBDATE,   MCDATE,   MCSEC,   JYEAR,   JMNTH,  JDAY,
     *      RUNTITLE, SHOWMAPS,
     *      BASEDATE, ENDDATE,  NSAVE,   NHISI,   NDHIS
#ifdef rcm
     *     ,SSTRCM,   IDATEFILE
#endif
#ifdef isotrac
     *     ,qgcm_o, qgcm_h, raicgcm_o, raicgcm_h, snocgcm_o, snocgcm_h,
     *      raisgcm_o, raisgcm_h, snosgcm_o, snosgcm_h,
     *      fvagcm_o, fvagcm_h
#endif
     *     )
c
c        Driver for surface models. Takes global arrays from gcm, calls
c        surface model subrs to step them through one agcm timestep,
c        one 1D strip of (land/ocean) points at a time. The time loop
c        uses the shortest model (lsx) timestep. At present this is the
c        same as the agcm step, but if it were shorter, the 1D-strip
c        loop is outside the lsx time loop so we wouldn't have to
c        transfer to global arrays every lsx timestep - also if the
c        other surface model (snow,soil) timesteps were integer divisors
c        of the agcm timestep, we wouldn't need global arrays for
c        communication between them and lsx.
c
c        All arguments are supplied unless otherwise noted:
c        NSTEP    = agcm timestep counter
c        DTGCM    = length of time to run lsx this call
c        FRADSW   = .TRUE. if agcm/rcm just done solar calculations
c        ZENGCM(JI,JJ)   = cosine of solar zenith angle
c        SRGCMD(JI,JJ,IB)= dwnwd direct  solar flux at SIGCMA
c                          This is per unit horizontal area, ie, it has
c                          already been multiplied by cos (zen angle).
c        SRGCMI(JI,JJ,IB)= dwnwd diffuse solar flux at SIGCMA
c        FIRABS(JI,JJ)   = absorbed IR flux at surface
c        UGCM(JI,JJ)     = eastward air velocity at SIGCMA
c        VGCM(JI,JJ)     = northward air velocity at SIGCMA
c        TGCM(JI,JJ)     = air temperature (deg K) at SIGCMA
c        QGCM(JI,JJ)     = specific humidity at SIGCMA
c        RAICGCM(JI,JJ)  = convective rainfall rate (mass flux)
c        SNOCGCM(JI,JJ)  = convective snowfall rate (mass flux)
c        RAISGCM(JI,JJ)  = stratiform rainfall rate (mass flux)
c        SNOSGCM(JI,JJ)  = stratiform snowfall rate (mass flux)
c        STOCHPR         = T/F flag for stochastic precip
c        PGCM(JI,JJ)     = surface pressure (N/m2)
c        CLDGCM(JI,JJ)   = total cloudiness (fraction)
c        FIREMS(JI,JJ)   = emitted IR flux at surface (returned)
c        FSEGCM(JI,JJ)   = dwnwd sensible heat flux at SIGCMA (returned)
c        FVAGCM(JI,JJ)   = downward H2O vapor flux at SIGCMA (returned)
c        TAXGCM(JI,JJ)   = eastward stress on surface (returned)
c        TAYGCM(JI,JJ)   = northward stress on surface (returned)
c        SIGCMA          = sigma (p/ps) at lowest gcm *level*
c        SIGCMB          = sigma (p/ps) at top of lowest gcm *layer*
c        OMEGA           = 2*pi/day length for dynam.ice Coriolis(rad/s)
c        SEQSST          = .TRUE. for sequential (eg, AMIP) prescribed
c                          SST/seaice , else (default) cyclic 12 months
c        DYNAMICE        = .TRUE. for dynamic seaice
c        DYNAMHEM        = .TRUE. for 2 hemispheric passes, .FALSE. for
c                           global (needed if sea ice at equator)
c        QFACTOR         = oceanic q-flux global factor
c        QNORWEG         = max value (W m-2) for Norwegian-Sea q-flux
c        QICES           = heat convergence under 100% sea ice, s.h.
c        QICEN           = heat convergence under 100% sea ice, n.h.
c        QLIMIT          = limit on SST gradients in SETQFLUX (also
c                          imposes uniform diffn coeff)
c        CO2PPM          = atmospheric co2 amount, ppmv
c        HOTLIMIT        = max limit for temps (deg C, -999.=do nothing)
c        NSREST          = -1 for cold start, 0 for initial start from a
c                          restart file, 1 for true restart
c        NSREST_IBIS     = -1 for IBIS dynamic veg cold start, 0 for
c                          initial start, 1 for restart
c        VERSRESIN       = version number of initial/restart file used
c                          to start this run
c        FSAVE           = .TRUE. if history file save this timestep
c        FNEXF           = .TRUE. if history next-file this timestep
c        ISTEPY          = gcm timestep # within current year
c        NRSTRT          = first timestep # for this run
c        NSTOP           = final timestep #
c        MBDATE          = base date for curr. sequence of runs (yymmdd)
c        MCDATE          = current date (at end of time step) (yymmdd)
c        MCSEC           = seconds into current day (at end of timestep)
c        JYEAR           = current year (2001 AD,etc, or -21000 BP,etc)
c        JMNTH           = current month number (1-12)
c        JDAY            = current day number (1-31)
c        RUNTITLE        = run title
c        SHOWMAPS        = .TRUE. to print 2-D maps of prescribed fields
c        BASEDATE        = basedate of run (passed to iniogcm)
c        ENDDATE         = end date of run (passed to iniogcm)
c        NSAVE           = save times for run (passed to iniogcm)
c        NHISI           = agcm history write times (passed to iniogcm)
c        NDHIS           = agcm history accum period (passed to iniogcm)
c        SSTRCM          = prescribed SSTs, deg K (for RCM only)
c
#include <compar.h>
#include <comgrd.h>
#include <combudl.h>
#include <comeve.h>
#include <comhisl.h>
#include <comice.h>
#include <comicexy.h>
#ifdef gcm
#include <COMLUN.h>
#endif
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <commus.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comvegxy.h>
#include <com1d.h>
      LOGICAL FRADSW, STOCHPR, SEQSST, DYNAMICE, DYNAMHEM, 
     *        FSAVE,  FNEXF,    SHOWMAPS
      CHARACTER*80 RUNTITLE
      DIMENSION
     *  ZENGCM(NLONG,NLATG),
     *  SRGCMD(NLONG,NLATG,NBAND), SRGCMI(NLONG,NLATG,NBAND),
     *  FIRABS(NLONG,NLATG),
     *  UGCM(NLONG,NLATG),         VGCM(NLONG,NLATG),
     *  TGCM(NLONG,NLATG),         QGCM(NLONG,NLATG),
     *  RAICGCM(NLONG,NLATG),      SNOCGCM(NLONG,NLATG),
     *  RAISGCM(NLONG,NLATG),      SNOSGCM(NLONG,NLATG),
     *  PGCM(NLONG,NLATG),         CLDGCM(NLONG,NLATG),
     *  FIREMS(NLONG,NLATG),
     *  FSEGCM(NLONG,NLATG),       FVAGCM(NLONG,NLATG),
     *  TAXGCM(NLONG,NLATG),       TAYGCM(NLONG,NLATG)
#ifdef rcm
      DIMENSION SSTRCM(NLONG,NLATG)
#endif
#ifdef isotrac
      dimension
     *      qgcm_o(nlong,nlatg),    qgcm_h(nlong,nlatg),
     *      raicgcm_o(nlong,nlatg), raicgcm_h(nlong,nlatg),
     *      snocgcm_o(nlong,nlatg), snocgcm_h(nlong,nlatg),
     *      raisgcm_o(nlong,nlatg), raisgcm_h(nlong,nlatg),
     *      snosgcm_o(nlong,nlatg), snosgcm_h(nlong,nlatg),
     *      fvagcm_o(nlong,nlatg),  fvagcm_h(nlong,nlatg),
     *      workl_o(nlon,nlat),     workl_h(nlon,nlat),
     *      workm_o(nlon,nlat),     workm_h(nlon,nlat)
#endif

 
      CHARACTER*20 BASEDATE,       ENDDATE
      DIMENSION NSAVE(NMARK),      NHISI(NMARK)
c
c     Local flags for global budget, above-ground budget and lsx dump:
      LOGICAL LBUDINI, ZBUDABV, ZDUMP
      SAVE ZBUDABV, ZDUMP
c     Local arrays:
      DIMENSION CONLOW(NVEC),      DFIRAXY(NLON,NLAT)

c     Local, for dynice (h<->u grids) 
      dimension 
     *  uocnxy_u (0:nlon+1,nlat+1),    vocnxy_u (0:nlon+1,nlat+1),
     *  fuocnet_u(0:nlon+1,nlat+1),    fvocnet_u(0:nlon+1,nlat+1),
     *  lmaskdice(nlon,nlat)
 
#ifdef gcm
      DATA ZBUDABV, ZDUMP /.FALSE.,.FALSE./
#else
      DATA ZBUDABV, ZDUMP /.TRUE., .FALSE./
#endif
 
#include <comsat.h>

#ifdef dotdot
c this allows progress to be monitored for short test runs
         print 453, NSTEP
  453    format (/I4,$)
#endif
c
c        Set run time at start of this time step
c
      TIMGCM = (NSTEP-1) * DTGCM
c
c        Transfer gcm downward solar fluxes to lsx grid, if gcm has
c        just done solar calcs  (FRADSW=.T.). If FRADSW=.F., downward
c        solar fluxes (SRGCM*) are unchanged from previous, so
c        SOLA*XY don't need to be recalculated.
c
c        Subr GTOL transfers gcm field to lsx grid, and GTOLSOL adjusts
c        SOLAD, SOLAI so that lsx global-mean absorbed solar flx equals
c        that implied by gcm rad calcs.
c
c     IF (FRADSW) THEN
       CALL GTOL (ZENGCM, COSZENXY, 0)

       DO 10 IB=1,NBAND
         CALL GTOL (SRGCMD(1,1,IB), SOLADXY(1,1,IB), 0)
#if defined (gcm) && ! defined (rcm)
         CALL GTOLSOL(SRGCMD(1,1,IB), SOLADXY(1,1,IB), ASURDXY(1,1,IB))
#endif
c
         CALL GTOL (SRGCMI(1,1,IB), SOLAIXY(1,1,IB), 0)
#if defined (gcm) && ! defined (rcm)
         CALL GTOLSOL(SRGCMI(1,1,IB), SOLAIXY(1,1,IB), ASURIXY(1,1,IB))
#endif
   10  CONTINUE
c     ENDIF
c
c        Transfer agcm sfc downward (absorbed by sfc) ir flx
c
      CALL GTOL (FIRABS, FIRAXY, 1)
c
c        Transfer agcm conditions at ZA to lsx global grid arrays.
c
      SIGA = SIGCMA
      SIGB = SIGCMB
      CALL GTOL (UGCM, UXXY, 0)
      CALL GTOL (VGCM, UYXY, 0)
      CALL GTOL (TGCM, TAXY, 0)
      CALL GTOL (QGCM, QAXY, 0)
      CALL GTOL (PGCM, PSURFXY, 0)
      CALL GTOL (CLDGCM, TOTCLDXY, 0)
#ifdef isotrac
      call gtol (qgcm_o, qaxy_o, 0)
      call gtol (qgcm_h, qaxy_h, 0)
      call zero (avapaxy_o, nlon*nlat)
      call zero (avapaxy_h, nlon*nlat)
#endif

c        If stochastic precip requested, subr stoch does stochastic
c        perturbs on lsx grid (different stochastic params for conv and
c        strat), then returns conv+strat. If not requested, add GCM
c        conv+strat then interpolate to lsx grid.

c+++++
c rain
c+++++

      CALL GTOL  (RAICGCM, WORKL, -1)
      CALL GTOL  (RAISGCM, WORKM, -1)
#ifdef isotrac
      call gtol  (raicgcm_o, workl_o, -1)
      call gtol  (raicgcm_h, workl_h, -1)
      call gtol  (raisgcm_o, workm_o, -1)
      call gtol  (raisgcm_h, workm_h, -1)
#endif
      if (lonbud.ne.-999 .and. latbud.ne.-999) then
        raic_p = workl(lonbud,latbud)
        rais_p = workm(lonbud,latbud)
      endif
c>>>>>>>>>>>>>>>>>>>>>>
      IF (STOCHPR) THEN
c>>>>>>>>>>>>>>>>>>>>>>
        CALL STOCH (WORKL, WORKM, RAINAXY)
#ifdef isotrac
        do 20 jj=1,nlat
          do 22 ji=1,nlon
            if ((workl(ji,jj)+workm(ji,jj)).ne.0.) then
              rainaxy_o(ji,jj) = rainaxy  (ji,jj)*
     *       (workl_o(ji,jj)+workm_o(ji,jj))/(workl(ji,jj)+workm(ji,jj))
              rainaxy_h(ji,jj) = rainaxy  (ji,jj)*
     *       (workl_h(ji,jj)+workm_h(ji,jj))/(workl(ji,jj)+workm(ji,jj))
            else
              rainaxy_o(ji,jj) = 0.
              rainaxy_h(ji,jj) = 0.
            endif
 22       continue
 20     continue
#endif
c>>>>>>>>>
      ELSE
c>>>>>>>>>
        DO JJ=1,NLAT
          DO JI=1,NLON
            RAINAXY(JI,JJ) = WORKL(JI,JJ) + WORKM(JI,JJ)
#ifdef isotrac
            rainaxy_o(ji,jj) = workl_o(ji,jj)+workm_o(ji,jj)
            rainaxy_h(ji,jj) = workl_h(ji,jj)+workm_h(ji,jj)
#endif
          ENDDO
        ENDDO
c>>>>>>>>>>
      ENDIF
c>>>>>>>>>>

c       Augment rainfall from irrigation reservoir (soilctl,commudxy)

#ifndef gcm
c     zirrig_rate = 100./86400.              ! irrigation rate (kg/m2/s)
c     girrig_max  = 500.                     ! max capacity (kg/m2)
c     do jj=1,nlat
c       do ji=1,nlon
c         if (lmaskxy(ji,jj).eq.1) then
c           if (taxy(ji,jj).gt.tmelt) then
c             dirrig = min (girrigxy(ji,jj), zirrig_rate*dtgcm)
c             dirrig = max (dirrig, girrigxy(ji,jj) - girrig_max)
c             rainaxy(ji,jj) = rainaxy(ji,jj) + dirrig/dtgcm
c             girrigxy(ji,jj) = girrigxy(ji,jj) - dirrig
c           endif
c         endif
c       enddo
c     enddo
#endif

c+++++
c snow
c+++++

      CALL GTOL  (SNOCGCM, WORKL, -1)
      CALL GTOL  (SNOSGCM, WORKM, -1)
#ifdef isotrac
      call gtol  (snocgcm_o, workl_o, -1)
      call gtol  (snocgcm_h, workl_h, -1)
      call gtol  (snosgcm_o, workm_o, -1)
      call gtol  (snosgcm_h, workm_h, -1)
#endif
      if (lonbud.ne.-999 .and. latbud.ne.-999) then
        snoc_p = workl(lonbud,latbud)
        snos_p = workm(lonbud,latbud)
      endif

c~~~~~~~~~~~~~~~~~~~~~~
      IF (STOCHPR) THEN
c~~~~~~~~~~~~~~~~~~~~~~
        CALL STOCH (WORKL, WORKM, SNOWAXY)
#ifdef isotrac
        do 30 jj=1,nlat
          do 32 ji=1,nlon
            if ((workl(ji,jj)+workm(ji,jj)).ne.0.) then
              snowaxy_o(ji,jj) = snowaxy  (ji,jj)*
     *       (workl_o(ji,jj)+workm_o(ji,jj))/(workl(ji,jj)+workm(ji,jj))
              snowaxy_h(ji,jj) = snowaxy  (ji,jj)*
     *       (workl_h(ji,jj)+workm_h(ji,jj))/(workl(ji,jj)+workm(ji,jj))
            else
              snowaxy_o(ji,jj) = 0.
              snowaxy_h(ji,jj) = 0.
            endif
 32      continue
 30   continue
#endif
c~~~~~~~~~
      ELSE
c~~~~~~~~~
        DO JJ=1,NLAT
          DO JI=1,NLON
            SNOWAXY(JI,JJ) = WORKL(JI,JJ) + WORKM(JI,JJ)
#ifdef isotrac
            snowaxy_o(ji,jj) = workl_o(ji,jj)+workm_o(ji,jj)
            snowaxy_h(ji,jj) = workl_h(ji,jj)+workm_h(ji,jj)
#endif
          ENDDO
        ENDDO
c~~~~~~~~~~
      ENDIF
c~~~~~~~~~~

#ifdef ibisint
c     perform deltas for interactive ibis veg (temp,qa/rh, and
c     rain,snow; nb: *after* stochpr).
      call ibis_delta (taxy,    qaxy,    psurfxy, lmaskxy,
     *                 siga, istepy, dtgcm, 1)
      call ibis_delta (rainaxy, snowaxy, xdum,    lmaskxy,
     *                 siga, istepy, dtgcm, 2)
#endif

c------------------------------------------------------
c       Single-point driving meteorol values to fort.68
c------------------------------------------------------

      if ((lonbud.ne.-999 .and. latbud.ne.-999) .and. .false.) then
c     if ((lonbud.ne.-999 .and. latbud.ne.-999)              ) then

c        Save single-point (lsx-grid) driving meteorol values, before
c        delevtype corrections, for tabular output to fort.68 below
c
        coszen_p = coszenxy(lonbud,latbud)
        solad1_p = soladxy(lonbud,latbud,1)
        solad2_p = soladxy(lonbud,latbud,2)
        solai1_p = solaixy(lonbud,latbud,1)
        solai2_p = solaixy(lonbud,latbud,2)
        fira_p = firaxy(lonbud,latbud)
        ux_p = uxxy(lonbud,latbud)
        uy_p = uyxy(lonbud,latbud)
        ta_p = taxy(lonbud,latbud)
        qa_p = qaxy(lonbud,latbud)
        psurf_p = psurfxy(lonbud,latbud)
        cld_p = totcldxy(lonbud,latbud)
 
c        Write single-point (lsx-grid) driving meteorol values
c        every timestep to fort.68 (header at start of each a year)
 
        if (istepy.eq.1) then
          zlon = ((lonbud-.5)/nlon)*360.
          if (zlon.gt.180.) zlon = zlon - 360.
          zlat = -90. + ((latbud-.5)/nlat)*180.
          write (68,6800) lonbud, latbud, zlon, zlat
 6800     format(
     *      'ilon,jlat=',2i5,'     lon,lat=',2f8.2
     *      /'      date   ts',
     *       '   coszen (frac)',
     *       '   solad1 (W/m2)', '   solad2 (W/m2)',
     *       '   solai1 (W/m2)', '   solai2 (W/m2)',
     *       '   firabs (W/m2)', '   totcld (frac)',
     *       '        ux (m/s)', '        uy (m/s)',
     *       '          ta (K)', '      qa (kg/kg)',
     *       '      psurf (mb)',
     *       '   raic (mm/day)', '   rais (mm/day)',
     *       '   snoc (mm/day)', '   snos (mm/day)',
     *       '       delev (m)'
     *      /)
        endif
 
        write (68,6802)
     *    mcdate,        nint(mcsec/dtgcm),
     *    coszen_p,
     *    solad1_p,      solad2_p,
     *    solai1_p,      solai2_p,
     *    fira_p,        cld_p,
     *    ux_p,          uy_p,
     *    ta_p,          qa_p,
     *    psurf_p*.01,
     *    raic_p*.864e5, rais_p*.864e5,
     *    snoc_p*.864e5, snos_p*.864e5,
     *    delevxy(lonbud,latbud)
 6802   format(i10,i4,17f16.9)
        ier = flush(68)

      endif

c---------------------------
c     end of fort.68 section
c---------------------------
c
c        Compute and print histogram of precip rates
c
c     CALL HISTOGRAM (NSTEP, DTGCM, 69, STOCHPR,  RAINAXY, SNOWAXY)
c
c        Report the current mixed-layer depth, first timestep only
c        (no longer, since ml depths are constant in v2)
c
c     IF (NSTEP.EQ.NRSTRT) CALL PRINHOCN
c
c        Set lsx history-write flags FWHIST_A,FAHIST_A,FRANGE,FRANGE24
c
#ifdef gcm
      CALL HISTPRE (DTGCM, FSAVE, FNEXF, ISTEPY, NSTEP, NRSTRT, NSTOP)
#endif
 
c        Apply requested elevation corrections to surface air
c        temperature (taxy) and downward incident IR flux (firaxy).
c        Save dnwd IR correction (dfiraxy) to apply later to returned
c        upward ir (firems), for energy conservation.
c        Also adjust surface pressure (psurfxy), and specific humidity
c        (qaxy) assuming constant relative humidity (7/18/95).
c
      IF (DELEVTYPE.EQ.0) THEN
 
        CALL ZERO (DFIRAXY, NLON*NLAT)
 
      ELSE
 
cv20a   ZLAPSE = .007                                     !7.0 deg K/km
        ZLAPSE = .008
        ZEPS =  RVAP/RAIR
 
        DO 40 JJ=1,NLAT
          DO 42 JI=1,NLON
 
            ZT1 = TAXY(JI,JJ)
            ZT2 = TAXY(JI,JJ) - ZLAPSE*DELEVXY(JI,JJ)
            ZRAT = 1. + (ZEPS-1.)*QAXY(JI,JJ)
 
            ZS1 = ESAT(ZT1)
            ZS2 = ESAT(ZT2)
            ZE1 = SIGA*PSURFXY(JI,JJ)*QAXY(JI,JJ)*ZEPS/ZRAT
c           Assume constant relative humidity:
            ZE2 = ZE1 * (ZS2/ZS1)
 
c           Vertically integrate hydrostatic eqn with const lapse rate:
            ZPOW = GRAV/(RAIR*ZRAT*ZLAPSE)
            ZP2 = SIGA*PSURFXY(JI,JJ) * (ZT2/ZT1)**ZPOW
            ZQ2 = ZE2 / (ZE2*(1.-ZEPS) + ZP2*ZEPS)
 
            ZE1 = .01*ZE1                                 !convert to mb
            ZE2 = .01*ZE2                                 !convert to mb
 
c           Sfc-air temperature and specific humidity, sfc pressure:
#ifdef isotrac
            if (qaxy(ji,jj).ne.0.) then
              qaxy_o(ji,jj) = zq2*qaxy_o(ji,jj)/qaxy(ji,jj)
              qaxy_h(ji,jj) = zq2*qaxy_h(ji,jj)/qaxy(ji,jj)
            else
              qaxy_o(ji,jj) = 0.
              qaxy_h(ji,jj) = 0.
            endif
#endif
            TAXY(JI,JJ) = ZT2
            QAXY(JI,JJ) = ZQ2
            PSURFXY(JI,JJ) = ZP2/SIGA
 
c           Total cloudiness:
            ZCLD1 = TOTCLDXY(JI,JJ)
c           ZCLD2 = ZCLD1                                     !no change
c           ZCLD3 = ZCLD1                                     !no change
c           protect against very negative delevxy:
            ZCLD2 = MIN (1., ZCLD1*EXP(-DELEVXY(JI,JJ)/2000.)) ! for ir
            ZCLD3 = ZCLD2                                      ! for sol
 
c           Downward IR radiation:
c           Oke,1987,p.373;Idso,1981,WRR,17,295-304;Sellers,1965,p.58:
            ZIR1 = (0.7 + 5.95E-5*ZE1*EXP(1500./ZT1)) * STEF*ZT1**4
     *           * (1. + 0.2 * ZCLD1**2)
            ZIR2 = (0.7 + 5.95E-5*ZE2*EXP(1500./ZT2)) * STEF*ZT2**4
     *           * (1. + 0.2 * ZCLD2**2)
            DFIRAXY(JI,JJ) = FIRAXY(JI,JJ) * (ZIR2/ZIR1 - 1.)
            FIRAXY(JI,JJ) = FIRAXY(JI,JJ) + DFIRAXY(JI,JJ)
 
c           Solar radiation: use simple expression for multiple refls,
c           giving surface *incident* flux = C(1-r)/(1-a*r) + 1-C,
c           where r = cloud transmission/reflection = 0.5, a = surface
c           albedo = 0.4, and C = cloudiness fraction. This does *not*
c           conserve global energy.
            ZSOL = (1.-.375*ZCLD3) / (1.-.375*ZCLD1)
            DO IB=1,NBAND
              SOLADXY(JI,JJ,IB) = SOLADXY(JI,JJ,IB)*ZSOL
              SOLAIXY(JI,JJ,IB) = SOLAIXY(JI,JJ,IB)*ZSOL
            ENDDO
 
c           Correct precip for (i) rain vs. snow depending on corrected
c           TA, and (ii) total amount dep. on SVP. The latter does
c           *not* conserve global h2o mass.
c           ZPREC = 1.                            ! no change
c           ZPREC = EXP(.05**(ZT2-ZT1))           ! Huy.,PPP,1991,89,399
            ZPREC = MAX (0.7, MIN (1.3, ZS2/ZS1)) ! dep on svp
            IF (ZT2.GE.TMELT) THEN
              RAINAXY(JI,JJ) = (RAINAXY(JI,JJ)+SNOWAXY(JI,JJ))*ZPREC
              SNOWAXY(JI,JJ) = 0.
#ifdef isotrac
              rainaxy_o(ji,jj)=(rainaxy_o(ji,jj)+snowaxy_o(ji,jj))*zprec
              snowaxy_o(ji,jj) = 0.
              rainaxy_h(ji,jj)=(rainaxy_h(ji,jj)+snowaxy_h(ji,jj))*zprec
              snowaxy_h(ji,jj) = 0.
#endif
            ELSE
              SNOWAXY(JI,JJ) = (RAINAXY(JI,JJ)+SNOWAXY(JI,JJ))*ZPREC
              RAINAXY(JI,JJ) = 0.
#ifdef isotrac
              snowaxy_o(ji,jj)=(rainaxy_o(ji,jj)+snowaxy_o(ji,jj))*zprec
              rainaxy_o(ji,jj) = 0.
              snowaxy_h(ji,jj)=(rainaxy_h(ji,jj)+snowaxy_h(ji,jj))*zprec
              rainaxy_h(ji,jj) = 0.
#endif
            ENDIF
 
   42     CONTINUE
   40   CONTINUE
 
      ENDIF
c
c        If first timestep of initial/cold start, set flag to initialize
c        global budget in budall. Even if restart, do if earlier
c        restart-file version (without full lsx budget info), or if
c        previous run ended with a budget calculation (in case model
c        has changed).
c
      LBUDINI = .FALSE.
      IF ( NSTEP.EQ.NRSTRT .AND.
     *     ( NSREST.LE.0 .OR. 
#if ! defined (rcm) && ! defined (hydmodel)
     *       VERSRESIN.LE.1.08 .OR.
#endif
     *       ABS((NSTEP-1)*DTGCM-TIMENEO).LT.DTTOL ) ) LBUDINI = .TRUE.
c
c        Loop over 1D land strips. Subr GETIND copies from global
c        arrs to 1D land strip arrays, and PUTIND copies them back.
c        (Indices have been set up in SETIND, called from inisurf).
c
#ifdef omp
#ifdef cray
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_l, indall_l, npall_l,
!$OMP*          dtgcm, co2ppm, timgcm, zbudabv, zdump, nstep, nstop)
!$OMP* private (i, istrip, npoi,
!$OMP*          afire, asena, avapa, ataux, atauy,
!$OMP*          firb,  fsena, fvapa, taux,  tauy,
!$OMP*          tsoi,  hsoi, consoi)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_l, indall_l, npall_l,
!$OMP*          dtgcm, co2ppm, timgcm, zbudabv, zdump, nstep, nstop)
!$OMP* private (i, istrip)
!$OMP* schedule(dynamic,2)
#endif
#endif
c
c-----------------------------
      DO 100 ISTRIP=1,NSTRIP_L
c-----------------------------
c
c          Get strip arrays from global arrays
c
        CALL GETIND (ISTRIP, INDALL_L, NPALL_L, 1)
c
c          Set physical soil quantities
c
        CALL SETSOI
c
c          Set physical snow quantities
c
        CALL SETSNO
c
c          Step lsx model
c
        CALL LSX (DTGCM, CO2PPM)
c
c          Arbitrarily set veg temps & intercepted h2o for no-veg locns
c
        CALL NOVEG
c
c          Save commug[xy] agcm fluxes from com1d strip arrays
c
        DO 140 I=1,NPOI
          AFIRE(I) = FIRB(I)
          ASENA(I) = FSENA(I)
          AVAPA(I) = FVAPA(I)
          ATAUX(I) = TAUX(I)
          ATAUY(I) = TAUY(I)
  140   CONTINUE
c
c          Above-ground budget accumulation (call before snow model
c          since snow model can change fractional snow cover FI).
c
        IF (ZBUDABV) CALL BUDAD1 (DTGCM)
c
c          Step snow model
c
        CALL SNOWCTL (DTGCM, TSOI(1,1), HSOI(1,1), CONSOI(1,1), 1)
c
c          Step soil model
c
        CALL SOILCTL (DTGCM, TIMGCM+DTGCM, ISTRIP)
c
c          Diagnostic printout for land/icesheet
c
        CALL PRINDIAGL (DTGCM, TIMGCM+DTGCM, ISTRIP)
c
c          Diagnostic "dump" for above-ground environment
c
        IF (ZDUMP) CALL PRINDUMP (DTGCM, TIMGCM+DTGCM)
c
c          Above-ground budget accum for effects of snow,soil,icesheet,
c          and compute and write out budget and daily means. Doesn't
c          work for ocean pts on Cray since some land/icesheet variables
c          are not set.
c
        IF (ZBUDABV) CALL BUDAD2 (DTGCM)
        IF (ZBUDABV) CALL BUDDO (DTGCM, TIMGCM+DTGCM, NSTEP, NSTOP)
c
c          Save npp-related results
c
c       CALL NPPCTL (DTGCM, TIMGCM+DTGCM, ISTRIP)
c
c          Copy strip arrays back to global arrays
c
        CALL PUTIND (1)
c
c-------------
  100 CONTINUE
c-------------
c
c        Loop over 1D icesheet strips.
c
#ifdef omp
#ifdef cray
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_i, indall_i, npall_i,
!$OMP*          dtgcm, timgcm, zbudabv, zdump, nstep, nstop)
!$OMP* private (i, istrip, npoi,
!$OMP*          afire, asena, avapa, ataux, atauy,
!$OMP*          firb,  fsena, fvapa, taux,  tauy,
!$OMP*          tsoi,  hsoi, consoi)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_i, indall_i, npall_i,
!$OMP*          dtgcm, timgcm, zbudabv, zdump, nstep, nstop)
!$OMP* private (i, istrip)
!$OMP* schedule(dynamic,2)
#endif
#endif
c
c-----------------------------
      DO 200 ISTRIP=1,NSTRIP_I
c-----------------------------
c
c          Get strip arrays from global arrays
c
        CALL GETIND (ISTRIP, INDALL_I, NPALL_I, 2)
c
c          Set physical icesheet quantities
c
        CALL SETSHEET
c
c          Set physical snow quantities
c
        CALL SETSNO
c
c          Calculate radiative and aerodynamic fluxes over icesheet
c
        CALL SHEETFLX (DTGCM)
c
c          Save commug[xy] agcm fluxes from com1d strip arrays
c
        DO 240 I=1,NPOI
          AFIRE(I) = FIRB(I)
          ASENA(I) = FSENA(I)
          AVAPA(I) = FVAPA(I)
          ATAUX(I) = TAUX(I)
          ATAUY(I) = TAUY(I)
  240   CONTINUE
c
c          Above-ground budget accumulation (call before snow model
c          since snow model can change fractional snow cover FI).
c
        IF (ZBUDABV) CALL BUDAD1 (DTGCM)
c
c          Step snow model
c
        CALL SNOWCTL (DTGCM, TSOI(1,1), HSOI(1,1), CONSOI(1,1), 2)
c
c          Step icesheet model
c
        CALL SHEETCTL (DTGCM, TIMGCM+DTGCM, ISTRIP)
c
c          Diagnostic printout for land/icesheet
c
        CALL PRINDIAGL (DTGCM, TIMGCM+DTGCM, ISTRIP)
c
c          Diagnostic "dump" for above-ground environment
c
        IF (ZDUMP) CALL PRINDUMP (DTGCM, TIMGCM+DTGCM)
c
c          Above-ground budget accum for effects of snow,soil,icesheet,
c          and compute and write out budget and daily means. Doesn't
c          work for ocean pts on Cray since some land/icesheet variables
c          are not set.
c
        IF (ZBUDABV) CALL BUDAD2 (DTGCM)
        IF (ZBUDABV) CALL BUDDO (DTGCM, TIMGCM+DTGCM, NSTEP, NSTOP)
c
c          Copy strip arrays back to global arrays
c
        CALL PUTIND (2)
c
c-------------
  200 CONTINUE
c-------------
c
c        If slab ocean requested (oceantype=1), set oceanic
c        heat flux. Do before FICEXY changed in loop 200. Also do
c        if ogcm (oceantype=2), for lsx ocean points that are land in
c        the ogcm.
c
#ifdef rcm
      CALL ZERO (QFLUXXY, NLON*NLAT)
#else
      IF (OCEANTYPE.EQ.0) THEN
        CALL ZERO (QFLUXXY, NLON*NLAT)
      ELSE
        CALL SETQFLUX (QFLUXXY, TOCNXY,  FICEXY, TFREEZXY, LMASKXY,
     *                 QFACTOR, QNORWEG, QICES,  QICEN, QLIMIT,
     *                 MCDATE,  SHOWMAPS)
      ENDIF
#endif
c
c        Loop over ocean/seaice strips.
c
#ifdef omp
#ifdef cray
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_o, indall_o, npall_o,
!$OMP*          dtgcm, timgcm, conice)
!$OMP* private (i, istrip, npoi, npice,
!$OMP*          afire, asena, avapa, ataux, atauy,
!$OMP*          firb,  fsena, fvapa, taux,  tauy,
!$OMP*          tice, hice, conlow)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (nstrip_o, indall_o, npall_o,
!$OMP*          dtgcm, timgcm, conice)
!$OMP* private (i, istrip, conlow)
!$OMP* schedule(dynamic,2)
#endif
#endif
c
c-----------------------------
      DO 300 ISTRIP=1,NSTRIP_O
c-----------------------------
c
c          Get strip arrays from global arrays
c
        CALL GETIND (ISTRIP, INDALL_O, NPALL_O, 3)
c
c          Set physical sea-ice quantities
c
        CALL SETICE
c
c          Set physical snow quantities
c
        CALL SETSNO
c
c          Calculate radiative and aerodynamic fluxes over ocean/seaice.
c
        CALL OCNFLX (DTGCM)
c
c          Save commug[xy] agcm fluxes from com1d strip arrays
c
        DO 340 I=1,NPOI
          AFIRE(I) = FIRB(I)
          ASENA(I) = FSENA(I)
          AVAPA(I) = FVAPA(I)
          ATAUX(I) = TAUX(I)
          ATAUY(I) = TAUY(I)
  340   CONTINUE
c
c          Step snow on seaice (if any ice pts in current strip - npice
c          is in comice and was set in ocnflx)
c
        IF (NPICE.GT.0.) THEN
          DO 350 I=1,NPOI
            CONLOW(I) = CONICE
  350     CONTINUE
          CALL SNOWCTL (DTGCM, TICE(1,1), HICE(1,1), CONLOW, 3)
        ENDIF
c
c          Step sea-ice thermodynamic model (vertical processes,
c          no change to sea-ice fraction FICE)
c
        CALL SEAICE (DTGCM)

#ifdef dotdot
c this allows progress to be monitored for short test runs
            print 454,"."
 454    format (A,$)
#endif

c
c          OCNML steps the ocean mixed-layer temperature, and
c          sea-ice processes that can change the fraction FICE.
c
        CALL OCNML (DTGCM)
c
c          Diagnostic printout for ocean
c
        CALL PRINDIAGO (DTGCM, TIMGCM+DTGCM, ISTRIP)
c
c          Copy strip arrays back to global arrays
c
        CALL PUTIND (3)
c
c-------------
  300 CONTINUE
c-------------

c        Adjust upward IR returned to the gcm by -1 * elev-correction
c        applied earlier to the downward IR at each point, so energy is
c        conserved

      DO 400 JJ=1,NLAT
        DO 402 JI=1,NLON
          AFIREXY(JI,JJ) = AFIREXY(JI,JJ) - DFIRAXY(JI,JJ)
  402   CONTINUE
  400 CONTINUE


c     isotopic surface processes:

#ifdef isotrac
#ifdef omp
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared  (dtgcm)
!$OMP* private (jj)
!$OMP* schedule(dynamic,2)
#endif
      do jj=1,nlat
        call surf_isot (jj, dtgcm)
      enddo
#endif

c        Set surface/canopy temperature (for rcm only)
c
c#ifdef rcm
c     do jj=1,nlat
c       do ji=i,nlon
c         soil and lower veg:
c         ztl = tlxy(ji,jj)*flxy(ji,jj) + tgxy(ji,jj)*(1.-flxy(ji,jj))
c
cc        lower surface including snow:
c         ztg = tsnoxy(ji,jj,1)*fixy(ji,jj) + ztl*(1.-fixy(ji,jj))
c
cc        upper canopy:
c
c         ztu =(tuxy(ji,jj)*laixy(ji,jj,2) + tsxy(ji,jj)*saixy(ji,jj,2))
c    *         / max (1.e-4, laixy(ji,jj,2) + saixy(ji,jj,2))
c
cc        land:
c         x = ztu*fuxy(ji,jj) + ztg*(1.-fuxy(ji,jj))
c
cc        sea ice/lake ice:
c         y = (1. - fi_icexy(ji,jj))*ticexy(ji,jj,1)
c    *            + fi_icexy(ji,jj) *tsno_icexy(ji,jj,1)
c
cc        open water:
c         z = tocnxy(ji,jj)
c
c         tsurfxy(ji,jj) = (1.-fwaterxy(ji,jj))*                      x
c    *                   +     fwaterxy(ji,jj) *(      ficexy(ji,jj) *y
c    *                                           + (1.-ficexy(ji,jj))*z)
c       enddo
c     enddo
c#endif

c        Transfer lsx fluxes to agcm arrays. LTOG transfers an lsx
c        field to the gcm grid.

      CALL LTOG (AFIREXY, FIREMS, 1)
      CALL LTOG (ASENAXY, FSEGCM, 1)
      CALL LTOG (AVAPAXY, FVAGCM, 1)
      CALL LTOG (ATAUXXY, TAXGCM, 1)
      CALL LTOG (ATAUYXY, TAYGCM, 1)
#ifdef isotrac
c     call ltog (avapaxy,   workg,    0)
c     call ltog (avapaxy_o, fvagcm_o, 0)
c     call ltog (avapaxy_h, fvagcm_h, 0)
c     do 650 jj=1,nlat
c       do 652 ji=1,nlon
c         if (workg(ji,jj).ne.0.) then
c            fvagcm_o(ji,jj)=fvagcm_o(ji,jj)*fvagcm(ji,jj)/workg(ji,jj)
c            fvagcm_h(ji,jj)=fvagcm_h(ji,jj)*fvagcm(ji,jj)/workg(ji,jj)
c         endif
c 652   continue
c 650 continue
      call ltog (avapaxy_o, fvagcm_o, 1)
      call ltog (avapaxy_h, fvagcm_h, 1)
#endif
c
c        Accumulate global budget fluxes and possibly compute global
c        storage and printout
c
#ifndef icedrive
      CALL BUDALL (DTGCM, TIMGCM+DTGCM, NSTEP, NSTOP, LBUDINI)
#endif
c
c        If requested, prescribe sst and seaice, overriding
c        ocean temps and seaice fractions and thicknesses set in OCNML.
c
c-----------------------------
      IF (OCEANTYPE.EQ.0) THEN
c-----------------------------
c
#ifdef rcm
        CALL DOPRESST_RCM (SSTRCM)
#else
        CALL DOPRESST (MBDATE, MCDATE, MCSEC, SEQSST, SHOWMAPS)
#endif
c
c----------
      ENDIF
c----------
c
c        If requested (DYNAMICE), step dynamic sea ice model dynice.
c        Call dynice in any case if OCEANTYPE=2 (ogcm) to set and filter
c        net ocean stresses F[U,V]OCNET. 
c
c------------------------------------------------------------------
      IF ((DYNAMICE.AND.(OCEANTYPE.EQ.1)) .OR. OCEANTYPE.EQ.2) THEN
c------------------------------------------------------------------
c
#if defined (gcm) && ! defined (rcm)

c       transfer [u,v]ocnxy (comocnxy,h-grid) to local u-grid arrays
        do jj=1,nlat+1
          jjc = min (jj,nlat)
          jjm = max (jj-1,1)
          do ji=1,nlon
            jim = mod (ji+nlon-2,nlon) + 1
c           if (lmaskxy(ji,jjc).eq.3 .and. lmaskxy(jim,jjc).eq.3 .and.
c    *          lmaskxy(ji,jjm).eq.3 .and. lmaskxy(jim,jjm).eq.3) then
            if (fwaterxy(ji,jjc).eq.1..and.fwaterxy(jim,jjc).eq.1..and.
     *          fwaterxy(ji,jjm).eq.1..and.fwaterxy(jim,jjm).eq.1.) then

              if (oceantype.eq.1) then
c               from lsx (h) grid:
                uocnxy_u(ji,jj)=0.25*(  uocnxy(ji,jjc)+uocnxy(jim,jjc)
     *                                + uocnxy(ji,jjm)+uocnxy(jim,jjm))
                vocnxy_u(ji,jj)=0.25*(  vocnxy(ji,jjc)+vocnxy(jim,jjc)
     *                                + vocnxy(ji,jjm)+vocnxy(jim,jjm))
              else
c               from mom2,4's u-grid(at NE corner of h grid,move to SW):
                uocnxy_u(ji,jj)= uocnxy(jim,jjm)
                vocnxy_u(ji,jj)= vocnxy(jim,jjm)
              endif

            else
              uocnxy_u(ji,jj) = 0.
              vocnxy_u(ji,jj) = 0.
            endif
          enddo
c         wraparound:
          uocnxy_u(0,jj)      = uocnxy_u(nlon,jj)
          vocnxy_u(0,jj)      = vocnxy_u(nlon,jj)
          uocnxy_u(nlon+1,jj) = uocnxy_u(1,jj)
          vocnxy_u(nlon+1,jj) = vocnxy_u(1,jj)
        enddo

c       Set lmaskdice so dynamic ice is done over 100% "lakes" as well 
c       as true lsx ocean (so for coupled ogcm, include true (and lsx)
c       ocean areas that are not ocean in ogcm, eg, Arctic)
        do jj=1,nlat
          do ji=1,nlon+1
            if (fwaterxy(ji,jj).eq.1.0) then
              lmaskdice(ji,jj) = 3
            else
              lmaskdice(ji,jj) = 0
            endif
          enddo
        enddo

        CALL DYNICE (DYNAMICE, DYNAMHEM,
     *               LMASKDICE, FICEXY, HICEXY, TICEXY, BRINEXY,
     *               FIXY, HSNOXY, TSNOXY,
     *               FUICEXY,FVICEXY, FUOCNXY,FVOCNXY,UOCNXY_U,VOCNXY_U,
     *               FUICEXYM,FVICEXYM,FUOCNXYM,FVOCNXYM,
     *               FUOCNET_U, FVOCNET_U,
     *               EPSILON, EPSILON,
     *               OMEGA,   DTGCM,   DTGCM*ISTEPY, DTGCM*NSTEP)

c       transfer f[u,v]ocnet_u (u-grid) to mom2,4 f[u,v]ocnet (comocnxy)
        do jj=1,nlat
          jjp = jj+1
          jjpa = min (jj+1, nlat)
          do ji=1,nlon
            jip = ji+1
            jipa = mod (ji,nlon) + 1

c           to lsx (h) grid:
c           if (lmaskxy(ji,jj).eq.3) then
c             fuocnet(ji,jj)=0.25*
c    *                       ( fuocnet_u(ji,jj)  + fuocnet_u(jip,jj)
c    *                       + fuocnet_u(ji,jjp) + fuocnet_u(jip,jjp))
c             fvocnet(ji,jj)=0.25*
c    *                       ( fvocnet_u(ji,jj)  + fvocnet_u(jip,jj)
c    *                       + fvocnet_u(ji,jjp) + fvocnet_u(jip,jjp))

c           to mom2,4's u grid, at NE corner of h grid
c           (not at SW corner as used by dynice):
c           if (lmaskxy(ji,jj).  eq.3 .or. lmaskxy(jipa,jj)  .eq.3 .or.
c    *          lmaskxy(ji,jjpa).eq.3 .or. lmaskxy(jipa,jjpa).eq.3) then
            if (lmaskdice(ji,jj)  .eq.3 .or. lmaskdice(jipa,jj).eq.3.or.
     *          lmaskdice(ji,jjpa).eq.3 .or. lmaskdice(jipa,jjpa).eq.3)
     *        then
              fuocnet(ji,jj) = fuocnet_u(jip,jjp)
              fvocnet(ji,jj) = fvocnet_u(jip,jjp)
            else
              fuocnet(ji,jj) = 0.
              fvocnet(ji,jj) = 0.
            endif
          enddo
        enddo
        call puthilsx (fuocnet, 56)
        call puthilsx (fvocnet, 57)
#endif
c
c---------
      ELSE
c---------
c
c          If dynice not called, write null lsx history fields
c
#ifdef gcm
        CALL ZERO (WORKL, NLON*NLAT)
        CALL PUTHILSX (WORKL, 38)
        CALL PUTHILSX (WORKL, 39)
        CALL PUTHILSX (WORKL, 44)
        CALL PUTHILSX (WORKL, 45)
        call puthilsx (workl, 56)
        call puthilsx (workl, 57)
#endif
c
c----------
      ENDIF
c----------

c---------
#ifdef gcm
c---------
c
c        Process daily-mean and monthly-mean t,p,rh for interactive
c        Eve/Leaf vegetation
c
      CALL EVEWRITE (DTGCM, NSTEP, NRSTRT, FSAVE, JYEAR, JMNTH, JDAY)
c
c        Process daily-mean and monthly-mean t,p,rh for interactive
c        Biome model
c
      CALL BIOWRITE (DTGCM, NSTEP, NRSTRT, FSAVE,
     *               JYEAR, JMNTH, JDAY)
c
c        Process climate and vegetation statistics for interactive
c        IBIS vegetation
c
#ifdef ibis
      CALL IBIS_STAT (DTGCM, NSREST_IBIS, NSTEP, NRSTRT, NSTOP,
     *                ISTEPY, MBDATE, JYEAR, JMNTH, JDAY)
#endif

#if defined (WEATHERING)
      call weathwrite (dtgcm, nstep, nrstrt, fsave,
     *                 jyear, jmnth, jday)
#endif
c
c        Write lsx fields to regional-model driver files
c
      CALL REGWRILSX (DTGCM, NSTEP, NRSTRT, ISTEPY, MBDATE,
     *                JYEAR, JMNTH, JDAY, MCSEC, RUNTITLE)
c
c        Write lsx history variables (either for ccm history file,
c        first averaged onto agcm grid, or directly to lsx history file)
c
      CALL HISTWRIT
c
c        Do next-write,save and next-file processing for lsx hist files
c
      CALL HISTPOST (DTGCM,  NSREST, FSAVE, FNEXF,
     *               NSTEP,  NRSTRT, NSTOP, ISTEPY,
     *               MBDATE, MCDATE, MCSEC, RUNTITLE
#ifdef rcm
     *              ,IDATEFILE
#endif
     *              )
c
c        Report the current mixed-layer depth, last timestep only
c        (no longer, since ml depths are constant in v2)
c
c     IF (NSTEP.EQ.NSTOP) CALL PRINHOCN

c-----
#endif
c-----

c        Impose upper limit (agcm namelist) on soil and veg temperatures

      if (hotlimit.ne.-999.) then
        do j=1,nlat
          do i=1,nlon
            tocnxy(i,j) = min (tocnxy(i,j),tmelt+hotlimit)
            tuxy(i,j) = min (tuxy(i,j),tmelt+hotlimit)
            tsxy(i,j) = min (tsxy(i,j),tmelt+hotlimit)
            tlxy(i,j) = min (tlxy(i,j),tmelt+hotlimit)
            do k=1,nsoilay
              tsoixy(i,j,k) = min (tsoixy(i,j,k), tmelt+hotlimit)
            enddo
          enddo
        enddo
      endif

      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine dopresst (mbdate, mcdate, mcsec, seqsst, showmaps)
 
c     Reads mid-month SST fields from prescribed sst/ice data file
c     into arrays sstprev and sstnext. These values are in 10 * deg C,
c     with 10% seaice indicated by value -20. Prescribes current ocean
c     mixed-layer temperature and sea-ice fractions by linearly
c     interpolating in time, and reads file to update sstprev,sstnext
c     whenever a mid-month time point is crossed.
 
c     If namelist SEQSST is false (default), cycles through 12
c     mid-months (Jan-Dec) on the data file. If true, reads the file
c     sequentially assuming the first field on the data file
c     represents the mid-month preceding the basedate of this run.
 
c     For sequential: sets ibm1 and icm1, the indices of the mid-month
c     points preceding the basedate and current date respectively.
c     Then the field #s in the file surrounding the current date are:
c       im1 = (jcy-jby)*nmon + icm1 - ibm1 + 1, and im2 = im1 + 1,
c     where jby and jcy are the basedate year # and current year #.
 
c     For cyclic: much the same logic, except all basedate info is
c     ignored, and im1=icm1, im2=im1+1, reset appropriately if 0 or 13.
 
c mbdate   = base date for curr. sequence of runs (yymmdd)
c mcdate   = current date (at end of time step) (yymmdd)
c mcsec    = seconds into current day (at end of timestep)
c seqsst   = .TRUE. for sequential SSTs, .FALSE. (default) for cyclic
c showmaps = .TRUE. to print 2-D maps of prescribed fields
 
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commapl.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
 
      logical seqsst, showmaps
 
      parameter (nmon=12)
      dimension sstprev(nlon,nlat), sstnext(nlon,nlat)
      dimension tz(nlon), fz(nlon)
      character ctitle*80
      logical first
      data first /.true./
 
      dimension daybeg(nmon), daydata(0:nmon+1)
      character cmon(nmon)*3
c     Calendar days at start of each month:
      data daybeg /            0,    31,    59,    90,   120,   151,
     *                       181,   212,   243,   273,   304,   334/
c     Mid-month calendar days which SST data fields represent:
      data daydata / -15.5, 15.5,  45.0,  74.5, 105.0, 135.5, 166.0,
     *                     196.5, 227.5, 258.0, 288.5, 319.0, 349.5,
     *               380.5 /
      data cmon /          'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     *                     'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /
 
      save first, loc, im2sav, daybeg, daydata, cmon, sstprev, sstnext,
     *     ibm1, jby, jby_data, jbm_data
 
 
c        If first call and if sequential SSTs are specified (seqsst),
c        calculate calendar info (jby,jbm,jbd,zbcal) for basedate,
c        and find the index ibm1 of the mid-month preceding basedate.
c        Nb: ibm1 and jby are saved variables.
 
c---------------------------------
      if (first .and. seqsst) then
c---------------------------------
c       if >= 0, relative to 0 A.D; if < 0, relative to 1950 A.D. (B.P.)
        jby = mbdate/10000
        jbm = mod(iabs(mbdate)/100,100)
        jbd = mod(iabs(mbdate),100)
        zbcal = daybeg(jbm) + jbd-1 + 0.
 
        do im=1,nmon+1
          if (zbcal.le.daydata(im)) then
            ibm1 = im-1
            goto 100
          endif
        enddo
  100   continue
 
c          Set year and month indices of the data field preceding
c          basedate, used later **for messages only**. These can
c          differ from jby and jbm since they apply to the mid-month
c          fields on the file, not to the basedate itself.
c          Nb: these are saved variables.
 
        if (ibm1.eq.0) then
          jby_data = jby-1
          jbm_data = 12
        else
          jby_data = jby
          jbm_data = ibm1
        endif
 
c       first is set to false in another if block below
c----------
      endif
c----------
 
c        Calculate calendar info (jcy,jcm,jcd,zccal) from current
c        mcdate and mcsec, and find the index icm1 and weight wm1 of
c        the mid-month preceding the current date and time.
 
c     if >= 0, relative to 0 A.D; if < 0, relative to 1950 A.D. (B.P.)
      jcy = mcdate/10000
      jcm = mod(iabs(mcdate)/100,100)
      jcd = mod(iabs(mcdate),100)
      zccal = daybeg(jcm) + jcd-1 + mcsec/86400.
 
      do im=1,nmon+1
        if (zccal.le.daydata(im)) then
          icm1 = im-1
          wm1 = (daydata(im)-zccal) / (daydata(im)-daydata(icm1))
          goto 150
        endif
      enddo
  150 continue
 
c        Set year and month indices of data field following the
c        current date, used later **for messages only**. These can
c        differ from jcy and jcm since they apply to the mid-month
c        fields on the file, not to the current date itself.
 
      if (icm1.eq.0) then
        jcy_data = jcy
        jcm_data = 1
      else if (icm1.eq.12) then
        jcy_data = jcy+1
        jcm_data = 1
      else
        jcy_data = jcy
        jcm_data = icm1+1
      endif
 
c        Set indices corresponding to the data fields preceding and
c        following the current time. If namelist SEQSST parameter is
c        .true. (ie, sequential SSTs), the first field on the file is
c        assumed to correspond to the mid-month immediately preceding
c        basedate. If not (default), there are assumed to be 12 monthly
c        data fields on the file (Jan-Dec) and cycle through them
c        just by rewinding the file when mid-December is crossed.
 
      if (seqsst) then
c       "+ 1" here since field #1 on file immediately precedes basedate
        im1 = (jcy-jby)*nmon + icm1 - ibm1 + 1
        im2 = im1 + 1
      else
        im1 = icm1
        if (im1.eq.0) im1 = nmon
        im2 = im1 + 1
        if (im2.eq.nmon+1) im2 = 1
      endif
 
 
c        If first call, read through file to "prev" position, then read
c        in "next". Also print first "next" field if requested.
 
c====================
      if (first) then
c====================
 
#ifdef gcm
        call findiu (nsstice, loc)
        if (seqsst) then
          write(ioterm2,200) jby,jbm,jbd, jcy,jcm,jcd,
     *                       locname(loc)(1:lenchr(locname(loc))),
     *                       cmon(jbm_data),jby_data,
     *                       cmon(jcm_data),jcy_data
  200     format(
     *      /' The prescribed SST/sea-ice data file will be read',
     *       ' sequentially (SEQSST = .T.).'
     *      /' The first field on this file should correspond to the',
     *       ' mid-month immediately'
     *      /' preceding the base date of this run.'
     *      /'   Base date of this run (BASEDATE)=',i6,'/',i2.2,'/',i2.2
     *      /'   Current date of this run        =',i6,'/',i2.2,'/',i2.2
     *      /' The file will now be read up to the two mid-months',
     *       ' spanning the current date.'
     *     //' Reading prescribed SST/sea-ice data file'/4x,a
     *      /' from ',a,i6,'  to ',a,i6,' ...')
        else
          write(ioterm2,202) locname(loc)(1:lenchr(locname(loc))),
     *                       cmon(1), cmon(im1), cmon(im2)
  202     format(/' Reading prescribed SST/sea-ice data file'/4x,a
     *           /' for months ',a,'-',a,' and ',a,'...')
        endif
 
        call openf (nsstice, locname(loc), ' ', 0, 0)
 
        do 210 imon=1,im1
          call readdata (nsstice, sstprev, nlon, nlat, 'SSTI', 2, imon)
  210   continue
        if ((.not.seqsst) .and. (im2.eq.1)) rewind (nsstice)
        call readdata (nsstice, sstnext, nlon, nlat, 'SSTI', 2, im2)
#endif
 
        if (showmaps) then
          do j=1,nlat
            do i=1,nlon
              if (abs(hocnxy(i,j)-hdeep).lt.0.01) then
                iworkl(i,j) = nint(sstnext(i,j))
              else
                iworkl(i,j) = 1000000
              endif
            enddo
          enddo
 
          if (seqsst) then
            write (ctitle,'(a,a,i6,a)')
     *        'prescribed SST for ',cmon(jcm_data),jcy_data,
     *        '  (10*deg C, seaice = -20)'
          else
            ctitle = 'prescribed SST for '//cmon(im2)
     *               //' (10*deg C, seaice = -20)'
          endif
          call mapout (ioterm2, iworkl, nlon, nlat,
     *                 alon, alat, 1000000, ctitle)
        endif
 
        im2sav = im2
 
        first=.false.
 
c==========
      endif
c==========
 
c        If next-month time, copy "next" to "prev", read new "next"
 
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      if (im2.ne.im2sav) then
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
        call scopy (nlon*nlat, sstnext, 1, sstprev, 1)
 
#ifdef gcm
        if (seqsst) then
          write(ioterm2,400) locname(loc)(1:lenchr(locname(loc))),
     *                       cmon(jcm_data),jcy_data
  400     format(/' Reading prescribed SST/sea-ice data file'/4x,a
     *           /' for ',a,i6,' ...')
        else
          write(ioterm2,402) locname(loc)(1:lenchr(locname(loc))),
     *                       cmon(im2)
  402     format(/' Reading prescribed SST/sea-ice data file'/4x,a
     *           /' for month ',a,'...')
        endif
        if ((.not.seqsst) .and. (im2.eq.1)) rewind (nsstice)
        call readdata (nsstice, sstnext, nlon, nlat, 'SSTI', 2, im2)
#endif
 
        im2sav = im2
 
c<<<<<<<<<<
      endif
c<<<<<<<<<<
 
c        Set snow, seaice and ocean fields depending on interpolated
c        value between sstprev and sstnext. Convert temperatures to
c        deg K. If seaice is present, set thickness totz and fraction
c        fracz as linear ramp functions of latitude.
c
c        Set only for "open-ocean" points with hocnxy=hdeep, set in
c        iniocn. "Inland lakes" (hocnxy=hshallow) use predicted
c        ml-slab calculations in ocean.f. Points with no open water
c        have hocnxy = 0.
 
      do j=1,nlat
 
        al = abs ((j-.5)*180./nlat - 90.)
 
        if (j.le.nlat/2) then
c         Southern hemisphere:
          zwei = max (0., min (1., (al-60.)/(70.-60.)))
          totz  = (1.-zwei)*1.0  + zwei*1.0
c         fracz = (1.-zwei)*0.98 + zwei*0.98
          fracz = (1.-zwei)*0.90 + zwei*0.90
        else
c         Northern hemisphere:
          zwei = max (0., min (1., (al-70.)/(90.-70.)))
          totz  = (1.-zwei)*1.0  + zwei*4.0
c         fracz = (1.-zwei)*0.95 + zwei*0.998
          fracz = (1.-zwei)*0.92 + zwei*0.998
        endif
        totz  = max (nicelay*hicetop, min (5.0, totz))
        fracz = max (0., min (1., fracz))

        if (ticrit.eq.-999.) then    
c         sea ice fraction based on latitude (set above)
          do i=1,nlon
           tz(i) = wm1*sstprev(i,j) + (1.-wm1)*sstnext(i,j)
           tz(i) = max ( 0.1*tz(i) + tmelt, tmelt )
           fz(i) =     wm1 * cvmgt(1.,0.,sstprev(i,j).lt.-19.999)
     *           + (1.-wm1)* cvmgt(1.,0.,sstnext(i,j).lt.-19.999)
           fz(i) = cvmgt (fz(i), 0., fz(i).gt.epsilon)
           fz(i) = min (fz(i), fracz)
         enddo
        else
c         sea ice fraction based on sst (namelist ticrit)
          do i=1,nlon
            tz(i) = wm1*sstprev(i,j) + (1.-wm1)*sstnext(i,j)
            tz(i) = max ( 0.1*tz(i) + tmelt, 271.2 )
            fz(i) = max (0., min (1., 
     *              (ticrit+tmelt - tz(i))/(ticrit+tmelt - 271.2)))
            fz(i) = cvmgt (fz(i), 0., fz(i).gt.epsilon)
          enddo
        endif

        do i=1,nlon
 
          if (abs(hocnxy(i,j)-hdeep).lt.0.01) then   ! not for "shallow"
 
            if (fz(i).eq.0.) then
 
c                Prescribed open ocean (set mask, sst, snow=0, seaice=0)
c                (sst includes global-shift namelist variable delsst)
 
crd           next line commented out for Rob DeConto fixsst runs:
              tocnxy(i,j) = tz(i) + delsst
 
              do k=1,nsnolay
                hsno_icexy(i,j,k) = 0.
                tsno_icexy(i,j,k) = tmelt
              enddo
              fi_icexy(i,j) = 0.
              ti_icexy(i,j) = tmelt
 
              do k=1,nicelay
                hicexy(i,j,k) = 0.
                ticexy(i,j,k) = tmelt
              enddo
              ficexy(i,j)   = 0.
              tg_icexy(i,j) = tmelt
 
            else
 
c                Prescribed seaice (set mask, sst, ice fraction,
c                ice thickness, ice temp if new ice; don't set snow)
 
crd           next if test added for Rob DeConto fixsst runs:
crd           if (fz(i).gt.0.5) then
              tocnxy(i,j) = tz(i) + delsst
crd           endif
 
              hicexy(i,j,1) = hicetop
              do k=2,nicelay
                hicexy(i,j,k) = (totz-hicetop) / (nicelay-1)
              enddo
              if (ficexy(i,j).eq.0.) then
                ticexy(i,j,1) = min (taxy(i,j), tmelt)
                do k=2,nicelay
                  ticexy(i,j,k) = ( (nicelay-k)*ticexy(i,j,1)
     *                              + (k-1)*tmelt )
     *                            / (nicelay-1)
                enddo
                tg_icexy(i,j) = ticexy(i,j,1)
              endif
              ficexy(i,j) = fz(i)
 
            endif
 
          endif
 
        enddo
 
      enddo
 
      return
      end
c
c----------------------------------------------------------------------
c
#ifdef rcm
      subroutine dopresst_rcm (sstrcm)
 
c     Sets SSTs and sea ice fields for RCM. Stripped-down version 
c     of subr dopresst, using current sst field sstrcm (deg K).
c     Values < tzcrit indicate sea ice. 
 
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
 
      dimension sstrcm(nlon,nlat)
 
      parameter (totz  = 1.0)    ! sea ice thickness
      parameter (fracz = 1.0)    ! sea ice fraction
      parameter (tzcrit = 272.0) ! cutoff value in sstrcm => sea ice
 
c        Set only for "open-ocean" points with hocnxy=hdeep, set in
c        iniocn. "Inland lakes" (hocnxy=hshallow) use predicted
c        ml-slab calculations in ocean.f. Points with no open water
c        have hocnxy = 0.

      do j=1,nlat
        do i=1,nlon

          if (lmaskxy(i,j).ge.3) then
c         if (abs(hocnxy(i,j)-hdeep).lt.0.01) then   ! not for "shallow"
 
            if (sstrcm(i,j).ge.tzcrit) then
 
c                Prescribed open ocean (set mask, sst, snow=0, seaice=0)
 
              tocnxy(i,j) = sstrcm(i,j)
 
              do k=1,nsnolay
                hsno_icexy(i,j,k) = 0.
                tsno_icexy(i,j,k) = tmelt
              enddo
              fi_icexy(i,j) = 0.
              ti_icexy(i,j) = tmelt
 
              do k=1,nicelay
                hicexy(i,j,k) = 0.
                ticexy(i,j,k) = tmelt
              enddo
              ficexy(i,j)   = 0.
              tg_icexy(i,j) = tmelt
 
            else
 
c                Prescribed seaice (set mask, sst, ice fraction,
c                ice thickness, ice temp if new ice; don't set snow)
 
              tocnxy(i,j) = 271.2
 
              hicexy(i,j,1) = hicetop
              do k=2,nicelay
                hicexy(i,j,k) = (totz-hicetop) / (nicelay-1)
              enddo
              if (ficexy(i,j).eq.0.) then
                ticexy(i,j,1) = min (taxy(i,j), tmelt)
                do k=2,nicelay
                  ticexy(i,j,k) = ( (nicelay-k)*ticexy(i,j,1)
     *                              + (k-1)*tmelt )
     *                            / (nicelay-1)
                enddo
                tg_icexy(i,j) = ticexy(i,j,1)
              endif
              ficexy(i,j) = fracz
 
            endif
 
          endif
 
        enddo
      enddo
 
      return
      end
#endif
c
c----------------------------------------------------------------------
c
      subroutine stoch (pconv, pstrat, pstoch)
c
c       Applies stochastic perturbations to convective and stratiform
c       precip rates on lsx grid, then adds the two. Called separately
c       for rain and snow to maintain lsx's distinction between the two.
c
c       The probability density function (prob of precip between x
c       and dx) for locations with x > 0 is:
c         f(x)*dx = (q/P) * exp(-q*x/P) * dx
c       where q is the fraction of the grid square with x > 0,
c       and P is the mean (GCM) precip rate over the whole grid square.
c       (eg, Thomas and Henderson-Sellers,1991,J.Climate,4,898-910).
c       Then the integral of f(x) from 0 to 1 = 1, and the integral
c       of x*f(x) from 0 to 1 = P/q.
c
c       Then given a uniformly distributed random variable r in the
c       range [0,1], the stochastic precip rate x is:
c          such that g(x) = r/q if 0<r<q,
c          or          x  = 0   if q<r<1
c       where g(x) is the integral of f() from 0 to x.
 
c pconv  = mean convective precip (rain or snow) rate (supplied)
c pstrat = mean stratiform precip (rain or snow) rate (supplied)
c pstoch = combined stochastic precip (rain or snow) rate (returned)
 
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
 
      dimension
     *  pconv(nlon,nlat), pstrat(nlon,nlat), pstoch(nlon,nlat)
 
      logical first
      save first
      data first /.true./
      data ztot, ntoz /0., 0/
 
c-----------------------------------------------------------------------
c     qconv,qstrat are fractions of lsx grid square that get convective
c     or stratiform precip > 0.
#ifdef gcm
      parameter (qconv = 0.3, qstrat = 0.7)  !qstrat was 1.0 in early v2
#else
      parameter (qconv = 1.0, qstrat = 1.0)
#endif
c-----------------------------------------------------------------------
 
#ifdef sun
      real*4 rand
#endif
 
c       Set random seed for ranf() if first call (then identical runs
c       with the identical model will reproduce results exactly).
 
      if (first) then
 
#if   defined ( aix )
        call random_seed()
#elif defined ( cray )
        zdum = ranset(12345)
#elif defined ( sgi )
        call srand(12345)
#elif defined ( sun )
        zdum = rand(12345)
#elif defined ( linux )
        zdum = 0
        call random_seed()
#else
        write (*,*) 'Error (stoch): unrecognized machine type' 
        call endrun (-1)
#endif

        write(*,10) qconv, qstrat
   10   format(/'*** Using stochastic precipitation perturbations:',
     *          '  qconv=',f6.3,'  qstrat=',f6.3)
        first = .false.
      endif
 
      do 100 j=1,nlat
        do 102 i=1,nlon
 
c         Get stochastic convective rate
 
#if   defined ( aix )
          call random_number(zran)
#elif defined ( cray )
          zran = ranf()
#elif defined ( sgi )
          zran = rand(0)
#elif defined ( sun )
          zran = rand(0)
#elif defined ( linux )
          call random_number(zran)
#else
        write (*,*) 'Error (stoch): unrecognized machine type' 
        call endrun (-1)
#endif
 
c         qq = qconv
c         Eltahir and Bras, 1993, J.Climate, 6, 641, eq.(4):
          qq = max (.001, min (0.7, pconv(i,j)*86400./100.))
 
          if (zran.lt.qq-.0001) then
            zc = -(pconv(i,j)/qq) * log(1.-(zran/qq))
          else
            zc = 0.
          endif
 
c         Get stochastic stratiform rate
 
#if   defined ( aix )
          call random_number(zran)
#elif defined ( cray )
          zran = ranf()
#elif defined ( sgi )
          zran = rand(0)
#elif defined ( sun )
          zran = rand(0)
#elif defined ( linux )
          call random_number(zran)
#else
        write (*,*) 'Error (stoch): unrecognized machine type' 
        call endrun (-1)
#endif
c         qq = qstrat
c         Eltahir and Bras, 1993, J.Climate, 6, 641, eq.(4):
          qq = max (.001, min (1.0, pstrat(i,j)*86400./25.))
 
          if (zran.lt.qq-.0001) then
            zs = -(pstrat(i,j)/qq) * log(1.-(zran/qq))
          else
            zs = 0.
          endif
 
c         Add convective and stratiform rates
          pstoch(i,j) = zc + zs
 102    continue
 100  continue
 
#ifndef gcm
c     Diagnostics for local machine:
c     ztot = ztot + pstoch(1,1)
c     ntoz = ntoz + 1
c     write(*,200) pconv(1,1)*86400., pstoch(1,1)*86400.
c    *             (ztot*86400.)/ntoz
c 200 format('Stoch:  pconv=',f8.3,'  pstoch=',f8.3,
c    *       '  <pstoch>=',f8.3)
#endif
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine histogram (nstep, dtgcm, iu, stochpr, rainaxy, snowaxy)
c
c       Accumulates and prints global histogram of time-mean precip
c       rates. The bin divisions are set by data stmt below for zbin
c       (mm/day). The averaging period is ztimav (secs).
 
c Parameters (all supplied):
c nstep   = gcm timestep counter
c dtgcm   = gcm timestep, sec
c iu      = logical unit number
c stochpr = stochastic-precip T/F flag (for display only)
c rainaxy = instantaneous rainfall rate (kg/m2/s)
c snowaxy = instantaneous snowfall rate (kg/m2/s)
 
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
 
      dimension rainaxy(nlon,nlat), snowaxy(nlon,nlat)
      logical stochpr
 
      dimension phistog(nlon,nlat)
      parameter (nbin=24)
      dimension zbin(nbin+1), zfreq(nbin)
      character cfreq(nbin)*4
      logical first, firstpr
      save first, firstpr, phistog, nhistog, zbin, ztimav
      data first, firstpr /.true., .true./
      data zbin /
     *    0,  .1,   5,  10,  15,  20,  25,  30,  35,  40,
     *   45,  50,  55,  60,  65,  70,  75,  80,  85,  90,
     *   95, 100, 105, 110, 999/
      data ztimav /86400./          ! for daily means
c     data ztimav /0./              ! for instantaneous
 
 
c        Initialize daily mean accumulators
 
      if (first) then
        call zero (phistog,  nlon*nlat)
        nhistog = 0
        first = .false.
      endif
 
c         Accumulate time-mean precip rates
 
      do 10 j=1,nlat
       do 12 i=1,nlon
         phistog(i,j) = phistog(i,j) + (rainaxy(i,j)+snowaxy(i,j))
   12   continue
   10 continue
      nhistog = nhistog + 1
 
c        If time, compute+print histogram, reset accumulators
 
c------------------------------------------------------------
      if ( mod(nstep,nint(max(ztimav/dtgcm,1.))) .eq.0 ) then
c------------------------------------------------------------
 
c          Compute fractional areas for each precip bin
 
        call zero (zfreq, nbin)
        zptot = 0.
        zweitot = 0.
 
        do 20 j=1,nlat
          zwei = cos(alat(j))
          do 22 i=1,nlon
            zz = phistog(i,j) * (86400.*1.e3/rhow) / nhistog
            do 24 ib=1,nbin
              if (zz.ge.zbin(ib) .and. zz.lt.zbin(ib+1)) then
                zfreq(ib) = zfreq(ib) + zwei
                goto 26
              endif
   24       continue
            if (zz.lt.zbin(1)     ) zfreq(1)    = zfreq(1)    + zwei
            if (zz.ge.zbin(nbin+1)) zfreq(nbin) = zfreq(nbin) + zwei
   26       continue
            zptot   = zptot   + zz*zwei
            zweitot = zweitot +    zwei
   22     continue
   20   continue
        zptot = zptot/zweitot
 
        do 30 ib=1,nbin
          zfreq(ib) = zfreq(ib)/zweitot
   30   continue
 
c          Transfer to output character array (to show zeros as blank)
 
        do 40 ib=1,nbin
          if (zfreq(ib).eq.0.) then
            cfreq(ib) = ' '
          else
            write(cfreq(ib),'(i4)') nint(100.*zfreq(ib))
          endif
   40   continue
 
c          Print out histogram for this timestep
 
        if (firstpr) then
          write (iu,100) stochpr,
     *                   zbin(1), zbin(2), (nint(zbin(ib)),ib=3,nbin+1)
          write (iu,*)
          firstpr = .false.
        endif
        zday = nstep*dtgcm/86400.
        write (iu,102) zday, (cfreq(ib),ib=1,nbin), zptot
 
  100   format( 'Histogram of precipitation rates:',
     *          ' % global area vs mm/day.  (stochpr = ',l1,')'
     *        //'      day', 2x, 2f4.1, 23i4)
  102   format(       f9.3 , 4x, 24a4, '     mean=',f5.2)
 
c          Flush file buffer
 
        ier = flush (iu)
 
c          Reset daily mean accumulators
 
        call zero (phistog, nlon*nlat)
        nhistog = 0
 
c----------
      endif
c----------
 
      return
      end
c
c----------------------------------------------------------------------
c
#ifdef ibisint
      subroutine ibis_delta (ta, qa, ps, mask, siga, istepy,dtgcm,icode)

c     Performs climate deltas on lowest-level air temp, rh, precip,
c     for interactive IBIS vegetation. If first call, reads
c     files for observed and control monthly t, rh, precip, and saves
c     deltas (or ratios for precip) for each month. Interpolates in 
c     time between mid-months using istepy*dtgcm.

c     Only do this over land, since don't want to influence gcm/lsx over
c     ocean or sea ice. (Also, in case data files are nearest-neigbor 
c     interpolated from land, don't want to use garbage values over 
c     mid-oceans!)

c     ta = air temp (icode=1), or 1st rain/snow rate (icode=2)   (modif)
c     qa = air spec hum(icode=1), or 2nd rain/snow rate(icode=2) (modif)
c     ps = sfc pressure (only used for icode=1)
c     mask = land/ocean mask (1=land,2=icesheet,3=ocean)
c     siga  = pressure (agcm lowest level) / surface pressure
c     istepy = time step within current year
c     dtgcm = time step (s)
c     icode = 1 for (ta,qa); = 2 for rain/snow rates (still in ta,qa)

#include <compar.h>
#include <comgrd.h>

      parameter (nmon=12)

      dimension
     *  ta(nlon,nlat), qa(nlon,nlat), ps(nlon,nlat), mask(nlon,nlat)

      dimension 
     *  tmp1(nlon,nlat,nmon), tmp2(nlon,nlat,nmon),
     *  tomc(nlon,nlat,nmon), romc(nlon,nlat,nmon), pomc(nlon,nlat,nmon)

c     mid-month calendar days which obs,ctl data fields represent:
      dimension daydata(0:nmon+1)
      data daydata / -15.5, 15.5,  45.0,  74.5, 105.0, 135.5, 166.0,
     *                     196.5, 227.5, 258.0, 288.5, 319.0, 349.5,
     *               380.5 /

      logical first
      data first /.true./

      save tomc, romc, pomc, daydata, first

#include <comsat.h>

c        If first call this run, read in control and observed fields

c--------------------
      if (first) then
c--------------------
        write (6,*) '  ibisdelta: reading observed,control temp. files'
        call readibisfile('dat/dat_temp.2x2', tmp1, nlon,nlat,nmon,90)
        call readibisfile('ctl/ctl_temp.2x2', tmp2, nlon,nlat,nmon,90)
        do k=1,nmon
          do j=1,nlat
            do i=1,nlon
              tomc(i,j,k) = (tmp1(i,j,k)-tmp2(i,j,k))*.01   ! 100*C to K
            enddo
          enddo
        enddo

        write (6,*) '  ibisdelta: reading observed,control r.h. files'
        call readibisfile('dat/dat_relh.2x2', tmp1, nlon,nlat,nmon,90)
        call readibisfile('ctl/ctl_relh.2x2', tmp2, nlon,nlat,nmon,90)
        do k=1,nmon
          do j=1,nlat
            do i=1,nlon
              romc(i,j,k) = (tmp1(i,j,k)-tmp2(i,j,k))*.0001  !100*% to 1
            enddo
          enddo
        enddo

        write (6,*) '  ibisdelta: reading observed,control prec. files'
        call readibisfile('dat/dat_prec.2x2', tmp1, nlon,nlat,nmon,90)
        call readibisfile('ctl/ctl_prec.2x2', tmp2, nlon,nlat,nmon,90)
        do k=1,nmon
          do j=1,nlat
            do i=1,nlon
              pomc(i,j,k) = max (tmp1(i,j,k),10.)          ! 10*mm/month
     *                    / max (tmp2(i,j,k),10.)
            enddo
          enddo
        enddo

        first = .false.
c----------
      endif
c----------

c        Calculate indices and weights for interpolation from mid-months

      zday = istepy*dtgcm / 86400.
      do k=1,nmon+1
        if (zday.le.daydata(k)) then
          km1 = k-1
          km2 = k
          wm1 = (daydata(km2)-zday) / (daydata(km2)-daydata(km1))
          wm2 = 1.-wm1
          if (km1.eq.0) km1 = nmon
          if (km2.eq.nmon+1) km2 = 1
          goto 50
        endif
      enddo
   50 continue

c        Do deltas, using (tomc,romc) or pomc time-interp from mid mnths

      if (icode.eq.1) then

c          temperature and humidity (additive for t and rh):

        zeps = rvap/rair

        do j=1,nlat
          do i=1,nlon
            if (mask(i,j).eq.1 .or. mask(i,j).eq.2) then

c             convert from spec hum to rel hum:
              zrat = 1. + (zeps-1.)*qa(i,j)
              zp = siga*ps(i,j)
              zs1 = esat(ta(i,j))
              zr1 = qa(i,j)*zp*zeps / (zs1 * (1.+(zeps-1.)*qa(i,j)))

c             do deltas:
              ta(i,j) = ta(i,j) + (wm1*tomc(i,j,km1)+wm2*tomc(i,j,km2))
              zr2     = zr1     + (wm1*romc(i,j,km1)+wm2*romc(i,j,km2))
              zr2 = max (.0001, min (.9999, zr2))

c             convert from rel hum to spec hum:
              zs2 = esat(ta(i,j))
              qa(i,j) = zr2*zs2 / (zp*zeps - zr2*zs2*(zeps-1.))

            endif
          enddo
        enddo

      else  if (icode.eq.2) then

c          rainfalls/snowfalls (ratios, not additive):

        do j=1,nlat
          do i=1,nlon
            if (mask(i,j).eq.1 .or. mask(i,j).eq.2) then

              zratio =  wm1*pomc(i,j,km1)+wm2*pomc(i,j,km2)
              ta(i,j) = ta(i,j)*zratio
              qa(i,j) = qa(i,j)*zratio
              ta(i,j) = min (500./86400., ta(i,j))        ! < ... mm/day
              qa(i,j) = min (500./86400., qa(i,j))        ! < ... mm/day
              ta(i,j) = max (0., ta(i,j))
              qa(i,j) = max (0., qa(i,j))

            endif
          enddo
        enddo

      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine readibisfile (cfile, arr, nlon, nlat, nval, iu)
 
c        Read a temp/prec/relh file for IBIS model veg (for climate
c        deltas) in 12i6 format. These are files "borrowed" from
c        the EVE veg model, and have longitude cuts at dateline,
c        which must be shifted here to Greenwich for LSX.

      character*(*) cfile
      dimension arr(nlon,nlat,nval)
 
      write (6,*) '    reading file ',cfile
 
      open (iu, file=cfile, status='old')

      do j=nlat,1,-1
        read (iu,'(12f6.0)') ((arr(i,j,k),k=1,nval),i=nlon/2 + 1,nlon),
     *                       ((arr(i,j,k),k=1,nval),i=1,nlon/2)
        read (iu,'()')
      enddo

      close (iu)

      return
      end

#endif
c
c----------------------------------------------------------------------
c
      SUBROUTINE NOVEG
c
c       If no veg surfaces exist, set prog vars to nominal values.
c       (Sensible fluxes FSEN[U,S,L], latent fluxes FVAP[U,S,L]*,
c       temperature T[U,S,L], and intercepted liquid, snow amounts
c       WLIQ[U,S,L], WSNO[U,S,L] have been calculated for a unit
c       leaf/stem surface, whether or not one exists.)
c
#include <compar.h>
#include <comgrd.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
c
      DO 100 I=1,NPOI
        TAV = (1.-FI(I))*TG(I) + FI(I)*TI(I)
c
        Z = CVMGT (0., 1., LAI(I,2).EQ.0.OR.FU(I).EQ.0.)
        TU(I) = Z*TU(I) + (1.-Z)*TAV
        WLIQU(I) = Z*WLIQU(I)
        WSNOU(I) = Z*WSNOU(I)
c
        Z = CVMGT (0., 1., SAI(I,2).EQ.0.OR.FU(I).EQ.0.)
        TS(I) = Z*TS(I) + (1.-Z)*TAV
        WLIQS(I) = Z*WLIQS(I)
        WSNOS(I) = Z*WSNOS(I)
c
        X = LAI(I,1) + SAI(I,1)
        Y = FL(I)*(1.-FI(I))
        Z = CVMGT (0., 1., X.EQ.0.OR.Y.EQ.0.)
        TL(I) = Z*TL(I) + (1.-Z)*TAV
        WLIQL(I) = Z*WLIQL(I)
        WSNOL(I) = Z*WSNOL(I)
c
  100 CONTINUE
c
      RETURN
      END

c-----------------------------------------------------------------------
 
c>>>>>>>>>>>>>
#ifdef isotrac
c>>>>>>>>>>>>>

      subroutine surf_isot (jj, dt)
 
c Does isotopic computations for the surface, modifies surface
c reservoirs and returns the isotopic vapor flux (downward)
c * = isotope
 
c  GCM<-->LSX
c p(_*)      = precipitation at the ground
c avapaxy(_*)= downward total vapor flux
c rainaxy(_*)= downward rain flux
c snowaxy(_*)= downward snow flux
c qaxy(_*)   = specific humidity at the lowest GCM level
c  LAND+ICESHEET
c fvc(_*)    = vapor flux over land
c rct_*      = isotopic ratio of the top soil reservoir
c hct_i      = amount of water in the top soil reservoir
c hctm_i     = max amount of water in the top soil reservoir
c rcb_*      = isotopic ratio of the lower soil reservoir
c hcb_i      = amount of water in the bottom soil reservoir
c  OCEAN+LAKES
c fw         = surface fraction of ocean and lakes
c fvo(_*)    = vapor flux over ocean and lakes
c ro_*       = isotopic ratio of ocean and lakes
c qoxy_i     = spec. sat. humid. at the ocean surf. (computed in OCNFLX)
c xoxy_i     = exponent of the kinetic fractionation factor ( idem )
c zfxy_i     = aerodynamic coefficient for the vapor flux ( idem )
c hl_i       = amount of water in the lake reservoir
c hldm_i     = max amount of water for deep lakes
c hlsm_i     = max amount of water for shallow lakes
c  SEA-ICE
c fi         = surface fraction of sea-ice
c fvi(_*)    = vapor flux over sea-ice
c ri_*       = isotopic ratio of the sea-ice reservoir (snow)
c hi_i       = amount of water in the sea-ice reservoir
c him_i      = max amount of water in the sea-ice reservoir
c ----------
#include <compar.h>
#include <comgrd.h>
#include <comicexy.h>
#include <commudxy.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsnoxy.h>
#include <comsoixy.h>
#include <COMLUN.h>
#include <COMTIM.h>
#include <PARA_I.h>


      iua = 82                ! logical unit number
      iub = 83                ! logical unit number
 1    format ('### [',i2,',',i2,'] # ocean_i= ',i4,' # fwater = ',f6.4)
 2    format (i6,2(1x,f6.4),12(1x,g13.7))
 3    format (i6,1x,f5.1,3(1x,g13.7,1x,f5.1,1x,f6.1),
     *                                4(1x,f7.3,1x,f5.1,1x,f6.1))
 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (nstep.eq.nrstrt .and. jj.eq.pjj_i) then 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        open (iua, form='formatted', status='unknown', file='surfa_i')
        write (iua,1)
     *    pii_i,pjj_i, mocean_i(pii_i,pjj_i), fwaterxy(pii_i,pjj_i)
        write (iua,'(a)')
     *   '^^Step ^^fi,n fi,n+1 ^^^^^^^^^^fvo '//
     *   '^^^^^^^^fvo_O ^^^^^^^^fvo_H ^^^^^^^^^^fvc ^^^^^^^^fvc_O '//
     *   '^^^^^^^^fvc_H ^^^^^^^^^^fvi ^^^^^^^^fvi_O ^^^^^^^^fvi_H '//
     *   '^^^^^^^^^^fvt ^^^^^^^^fvt_O ^^^^^^^^fvt_H'
        open (iub, form='formatted', status='unknown', file='surfb_i')
        write (iub,1)
     *    pii_i,pjj_i, mocean_i(pii_i,pjj_i), fwaterxy(pii_i,pjj_i)
        write (iub,'(a)') '^^Step ^^^Ta ^^^^^^^^^rain ^d18O ^^^d2H '//
     *      '^^^^^^^^^snow ^d18O ^^^d2H ^^^^^^^^^^^Qa ^d18O ^^^d2H '//
     *                  '^^^SoiT ^d18O ^^^d2H ^^^SoiB ^d18O ^^^d2H '//
     *                  'Openwat ^d18O ^^^d2H ^SeaIce ^d18O ^^^d2H'
c~~~~~~~~~~
      endif
c~~~~~~~~~~
 
c++++++++++++++++++++
      do 20 ji=1,nlon
c++++++++++++++++++++

          fw    = fwaterxy(ji,jj)
          fi    = ficexy_i(ji,jj)
          if (fi.eq.0.) hi_i(ji,jj) = 0.
          fvc   = fvcxy_i(ji,jj)
          fvc_o = 0.
          fvc_h = 0.
          fvo   = fvoxy_i(ji,jj)
          fvo_o = 0.
          fvo_h = 0.
          fvi   = fvixy_i(ji,jj)
          fvi_o = 0.
          fvi_h = 0.
          p_o   = 0.
          p_h   = 0.

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (fwaterxy(ji,jj).lt.1.) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c           LAND+ICE-SHEET:

c           precipitation:
            p = max ( 0., rainaxy(ji,jj) + snowaxy(ji,jj)
     *                                   - abs(grunofxy(ji,jj)) )
            if (p.ne.0.) then
              p_o = p * (rainaxy_o(ji,jj)+snowaxy_o(ji,jj))
     *                / (rainaxy  (ji,jj)+snowaxy  (ji,jj))
              p_h = p * (rainaxy_h(ji,jj)+snowaxy_h(ji,jj))
     *                / (rainaxy  (ji,jj)+snowaxy  (ji,jj))
            endif

c------------------------------
            if (fvc.le.0.) then
c------------------------------
c             evaporation:
              if (-fvc*dt.le.hct_i(ji,jj)) then
                fvc_o = rct_o(ji,jj)*fvc
                fvc_h = rct_h(ji,jj)*fvc
                hct_i(ji,jj) = hct_i(ji,jj) + fvc*dt
              else
                fvc_o = -(   hct_i(ji,jj)*rct_o(ji,jj) / dt
     *                     + (-fvc-hct_i(ji,jj)/dt)*rcb_o(ji,jj) )
                fvc_h = -(   hct_i(ji,jj)*rct_h(ji,jj) / dt
     *                     + (-fvc-hct_i(ji,jj)/dt)*rcb_h(ji,jj) )
                hct_i(ji,jj) = 0.
                hcb_i(ji,jj) = max(0., hct_i(ji,jj)+hcb_i(ji,jj)+fvc*dt)
              endif
c---------------
            else
c---------------
c             condensation (liquid or solid):
              if ((tgxy(ji,jj).ge.tmelt).or.(qaxy(ji,jj).le.0.)) then
                fvc_o = fvc*rct_o(ji,jj)
                fvc_h = fvc*rct_h(ji,jj)
              else
                fvc_o = ( fvc*avi_o(tgxy(ji,jj))*1.1
     *                    / (1.+ avi_o(tgxy(ji,jj))*0.1*dr_o)
     *                  )
     *                * ( qaxy_o(ji,jj) / qaxy(ji,jj) )
                fvc_h = ( fvc*avi_h(tgxy(ji,jj))*1.1
     *                    / (1.+ avi_h(tgxy(ji,jj))*0.1*dr_h)
     *                  )
     *                * ( qaxy_h(ji,jj)/qaxy(ji,jj) )
              endif
              p   = p   + fvc
              p_o = p_o + fvc_o
              p_h = p_h + fvc_h
c----------------
            endif
c----------------

            hpp = hct_i(ji,jj) + p*dt

            if (hpp.ne.0.) then
              rct_o(ji,jj) = (hct_i(ji,jj)*rct_o(ji,jj) + p_o*dt) / hpp
              rct_h(ji,jj) = (hct_i(ji,jj)*rct_h(ji,jj) + p_h*dt) / hpp
              hct_i(ji,jj) = min (hpp,hctm_i)
            endif

            if (hpp.gt.hctm_i) then
              rcb_o(ji,jj) = (   hcb_i(ji,jj)*rcb_o(ji,jj)
     *                         + (hpp-hctm_i)*rct_o(ji,jj) )
     *                       / (hcb_i(ji,jj)+hpp-hctm_i)
              rcb_h(ji,jj) = (   hcb_i(ji,jj)*rcb_h(ji,jj)
     *                         + (hpp-hctm_i)*rct_h(ji,jj) )
     *                       / (hcb_i(ji,jj)+hpp-hctm_i)
              hcb_i(ji,jj) = max (0., hcb_i(ji,jj)+hpp-hctm_i)
            endif

c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>
 
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          if (fwaterxy(ji,jj).gt.0.) then
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c--------------------------------------------------------------------
            if (mocean_i(ji,jj).ne.8 .and. mocean_i(ji,jj).ne.9) then
c--------------------------------------------------------------------

c             OCEANS AND SEAS:

c             evaporation and condensation:
              fvo_o = -zfxy_i(ji,jj)
     *              * ( qoxy_i(ji,jj)*ro_o(ji,jj) / avl_o(tocnxy(ji,jj))
     *                  - qaxy_o(ji,jj) )
     *              / ( dr_o**xoxy_o(ji,jj))
              fvo_h = -zfxy_i(ji,jj)
     *              * ( qoxy_i(ji,jj)*ro_h(ji,jj) / avl_h(tocnxy(ji,jj))
     *                  - qaxy_h(ji,jj) )
     *              / ( dr_h**xoxy_h(ji,jj))

c---------------
            else
c---------------
 
c             LAKES:

              if (fvo.le.0. .or. qaxy(ji,jj).eq.0.) then
c               evaporation:
                fvo_o = ro_o(ji,jj)*fvo
                fvo_h = ro_h(ji,jj)*fvo
              else
c               condensation, liquid or solid
                if (tgxy(ji,jj).ge.tmelt) then
                  fvo_o = fvo*avl_o(tgxy(ji,jj))
     *                  * qaxy_o(ji,jj)/qaxy(ji,jj)
                  fvo_h = fvo*avl_h(tgxy(ji,jj))
     *                  * qaxy_h(ji,jj)/qaxy(ji,jj)
                else
                  fvo_o = ( fvo*avi_o(tgxy(ji,jj))*1.1
     *                      / (1. + avi_o(tgxy(ji,jj))*0.1*dr_o)
     *                    )
     *                    * ( qaxy_o(ji,jj)/qaxy(ji,jj) )
                  fvo_h = (  fvo*avi_h(tgxy(ji,jj))*1.1
     *                       / (1. +  avi_h(tgxy(ji,jj))*0.1*dr_h)
     *                    )
     *                    * ( qaxy_h(ji,jj)/qaxy(ji,jj) )
                endif
              endif

              hpp = hl_i(ji,jj) + (rainaxy(ji,jj)+snowaxy(ji,jj)+fvo)*dt

              ro_o(ji,jj) = 
     *                 (  hl_i(ji,jj)*ro_o(ji,jj)
     *                    + (rainaxy_o(ji,jj)+snowaxy_o(ji,jj)+fvo_o)*dt
     *                 ) / hpp
              ro_h(ji,jj) =
     *                 (  hl_i(ji,jj)*ro_h(ji,jj)
     *                    + (rainaxy_h(ji,jj)+snowaxy_h(ji,jj)+fvo_h)*dt
     *                 ) / hpp
              zzz = cvmgt ( hldm_i, hlsm_i, mocean_i(ji,jj).eq.8 )
              hl_i(ji,jj) = max (10., min (hpp,zzz))
c----------------
            endif
c----------------

c------------------------------------------
            if (ficexy_i(ji,jj).gt.0.) then
c------------------------------------------
 
c             SEA-ICE:

              p   = snowaxy(ji,jj)
              p_o = snowaxy_o(ji,jj)
              p_h = snowaxy_h(ji,jj)

              if (fvi.lt.0.) then
c               evaporation:
                fvi_o = cvmgt (ri_o(ji,jj)*fvi,
     *                        -(hi_i(ji,jj)/dt)*ri_o(ji,jj)
     *                        -(-fvi-hi_i(ji,jj)/dt)*1.003*ro_o(ji,jj),
     *                        -fvi.le.hi_i(ji,jj)/dt
     *                        )
                fvi_h = cvmgt (ri_h(ji,jj)*fvi,
     *                        -(hi_i(ji,jj)/dt)*ri_h(ji,jj)
     *                        -(-fvi-hi_i(ji,jj)/dt)*1.020*ro_h(ji,jj),
     *                        -fvi.le.hi_i(ji,jj)/dt)
                hi_i(ji,jj) = max (0., hi_i(ji,jj)+fvi*dt)

              else

c               condensation, liquid or solid:
                if ( (tg_icexy(ji,jj).ge.tmelt).or.
     *               (qaxy(ji,jj).le.0.) ) then
                  fvi_o = fvi * cvmgt (ri_o(ji,jj), 1.003*ro_o(ji,jj),
     *                                              hi_i(ji,jj).gt.0.)
                  fvi_h = fvi * cvmgt (ri_h(ji,jj), 1.020*ro_h(ji,jj),
     *                                              hi_i(ji,jj).gt.0.)
                else
                  fvi_o = ( fvi*avi_o(tg_icexy(ji,jj))*1.1
     *                      / (1.+avi_o(tg_icexy(ji,jj))*0.1*dr_o)
     *                    )
     *                    * ( qaxy_o(ji,jj)/qaxy(ji,jj) )
                  fvi_h = ( fvi*avi_h(tg_icexy(ji,jj))*1.1
     *                      / (1.+avi_h(tg_icexy(ji,jj))*0.1*dr_h)
     *                    )
     *                    * ( qaxy_h(ji,jj)/qaxy(ji,jj) )
                endif

                p   = p  +fvi
                p_o = p_o+fvi_o
                p_h = p_h+fvi_h

              endif

              hpp = hi_i(ji,jj)+p*dt
              if (hpp.ne.0.) then
                ri_o(ji,jj) = (hi_i(ji,jj)*ri_o(ji,jj)+p_o*dt)/hpp
                ri_h(ji,jj) = (hi_i(ji,jj)*ri_h(ji,jj)+p_h*dt)/hpp
                hi_i(ji,jj) = min(hpp,him_i)
              endif
c----------------
            endif
c----------------
 
c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>
 
c         MERGE:

          avapaxy_o(ji,jj) = (fvi_o*fi + fvo_o*(1.-fi)) * fw
     *                     +  fvc_o                     * (1.-fw)
          avapaxy_h(ji,jj) = (fvi_h*fi + fvo_h*(1.-fi)) * fw
     *                     +  fvc_h                     * (1.-fw)
 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if ((nstep.eq.nrstrt.or.mod(nstep,pdn_i).eq.0).and.
     *                            (ji.eq.pii_i).and.(jj.eq.pjj_i)) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write (iua,2)
     *        nstep,fi,ficexy(ji,jj),
     *        fvo,fvo_o,fvo_h,fvc,fvc_o,fvc_h,fvi,fvi_o,fvi_h,
     *        avapaxy(ji,jj),avapaxy_o(ji,jj),avapaxy_h(ji,jj)
            if (rainaxy(ji,jj).ne.0.) then
              rain_o = (rainaxy_o(ji,jj)/rainaxy(ji,jj)/smoj_o-1.)*1e3
              rain_h = (rainaxy_h(ji,jj)/rainaxy(ji,jj)/smoj_h-1.)*1e3
            else
              rain_o = -1000.
              rain_h = -1000.
            endif
            if (snowaxy(ji,jj).ne.0.) then
              snow_o = (snowaxy_o(ji,jj)/snowaxy(ji,jj)/smoj_o-1.)*1e3
              snow_h = (snowaxy_h(ji,jj)/snowaxy(ji,jj)/smoj_h-1.)*1e3
            else
              snow_o = -1000.
              snow_h = -1000.
            endif
            qa_o = (qaxy_o(ji,jj)/qaxy(ji,jj)/smoj_o-1.)*1e3
            qa_h = (qaxy_h(ji,jj)/qaxy(ji,jj)/smoj_h-1.)*1e3
            write (iub,3)
     *        nstep,taxy(ji,jj)-273.16,
     *        rainaxy(ji,jj),rain_o,rain_h,
     *        snowaxy(ji,jj),snow_o,snow_h,
     *        qaxy(ji,jj),qa_o,qa_h,
     *        hct_i(ji,jj),(rct_o(ji,jj)/smoj_o-1.)*1e3,
     *        (rct_h(ji,jj)/smoj_h-1.)*1e3,hcb_i(ji,jj),
     *        (rcb_o(ji,jj)/smoj_o-1.)*1e3,(rcb_h(ji,jj)/smoj_h-1.)*1e3,
     *        hl_i(ji,jj),(ro_o(ji,jj)/smoj_o-1.)*1e3,
     *        (ro_h(ji,jj)/smoj_h-1.)*1e3,hi_i(ji,jj),
     *        (ri_o(ji,jj)/smoj_o-1.)*1e3,(ri_h(ji,jj)/smoj_h-1.)*1e3
c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~

c+++++++++++++++
 20     continue
c+++++++++++++++

      return
      end

c>>>>>
#endif
c>>>>>
