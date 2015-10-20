c -- $Id: comice.h,v 1.2 2001/07/03 21:29:24 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comice.h,v $
c -- $Name:  $

c-------
c comice
c-------
c
c     1-D arrays for sea-ice model. Most correspond to global arrays in
c     comicexy (mostly prognostic, except for diagnostic fuice, fvice)
c
c        HICE   = sea-ice layer thicknesses
c        TICE   = temperature for each sea-ice layer
c        FICE   = fractional seaice cover
c        BRINE  = amount of heat in brine resevoir (J m-2)
c        TG_ICE = skin temperature of sea ice
c        HITOTO = total seaice thickness at start of each timestep
c        FICEOLD= FICE before adjustments in OCNML (no corresp 2-D arr)
c        FUICE = eastwd wind stress imparted from atmos to seaice sfc
c        FVICE = northwd wind stress imparted from atmos to seaice sfc
c        INDICE = index into current ocean strip of pts with sea ice
c        NPICE  = number of points in INDICE index
c
      COMMON  /SILSX/
     *  HICE(NVEC,NICELAY),              TICE(NVEC,NICELAY),
     *  FICE(NVEC),                      BRINE (NVEC),
     *  TG_ICE(NVEC),
     *  HITOTO(NVEC),                    FICEOLD(NVEC),
     *  FUICE (NVEC),                    FVICE (NVEC),
#if defined TWOSTRICE
     *  fluxiced(nvec,nicelay+1),       fluxicei(nvec,nicelay+1),
#endif
     *  INDICE(NVEC),                    NPICE
#ifdef isotrac
c       Surface fraction of sea-ice and vapor flux above sea-ice
     * ,fice_i(NVEC),                    fvi_i(NVEC)
#endif
!$OMP THREADPRIVATE ( /SILSX/ )
c
c----------------------------------------------------------------------
c
c     Scalars and non-prognostic 1-D arrays for sea-ice thermodynamics
c     and dynamics, the latter corresponding to global arrays in 
c     comicexy.
c
c        BASEFAC= coefficient for vertical heat flux into seaice base
c        SIDEFAC= coefficient for lateral heat flux into seaice sides
c        SIDEZER= lateral heat flux factor value for FICE=0.
c
c        HICEMIN  = minimum vertically integrated seaice thickness
c        HICETOP  = "adaptive grid" thickness of top seaice layer
c        FICEMAXN = maximum fractional seaice cover, northern hemis
c        FICEMAXS = maximum fractional seaice cover, southern hemis
c                   (ficemax* are namelist, may be set in datgelsx)
c        TICRIT   = temperature for setting prescribed sea ice
c                   fraction in DOPRESST as fn of SST (namelist,
c                   default = not used and ice fraction = f(latitude))
c        HICEEXP  = e-folding thickness for ice strength (m)
c        Z0ICE    = roughness length of sea-ice surface (m)
c        RHOICE   = sea-ice density (Kg m-3)
c        CSICE    = specific heat of sea ice (J Kg-1 deg-1)
c        CONICE   = thermal conductivity of seaice (W m-1 K-1)
c        PENFR    = fraction of absorbed solar flux that penetrates ice
c        PENDE    = e-folding depth of penetrating solar flux into ice
c        HICECRE  = thicknesses for newly formed seaice layers
c        ZMCRE    = column mass of newly formed seaice
c
c        FBRINMAX = maximum brine resevoir amount as a fraction of the
c                   heat needed to completely melt the seaice column
c        FBRINSTO = fraction of seaice melt that is retained as "brine"
c                   (rest "leaks" immediately to ocean as meltwater)
c
      COMMON /SI2LSX/
     *  BASEFAC,             SIDEFAC,         SIDEZER,
c
     *  HICEMIN,             HICETOP,
     *  FICEMAXN,            FICEMAXS,        TICRIT,     HICEEXP,
     *  Z0ICE,               RHOICE,          
     *  CSICE,               CONICE,
c
     *  PENFR,               PENDE,
     *  HICECRE(NICELAY),    ZMCRE,
c
     *  FBRINMAX,            FBRINSTO
c
