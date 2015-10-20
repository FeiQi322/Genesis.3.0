c -- $Id: comvegxy.h,v 1.2 2001/07/03 21:29:27 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comvegxy.h,v $
c -- $Name:  $

c---------
c comvegxy
c---------
c
c     Global arrays for vegetation.
c     See comveg for variable descriptions, and veg scalars.
c
c----------------------------------------------------------------------
c
      COMMON /VEPLSXXY/
     *  TUXY (NLON,NLAT),     TSXY (NLON,NLAT),     TLXY (NLON,NLAT),
     *  TLSUBXY (NLON,NLAT),
     *  WLIQUXY (NLON,NLAT),  WLIQSXY (NLON,NLAT),  WLIQLXY(NLON,NLAT),
     *  WSNOUXY (NLON,NLAT),  WSNOSXY (NLON,NLAT),  WSNOLXY(NLON,NLAT),
c
     *  T12XY (NLON,NLAT),    T34XY (NLON,NLAT),
     *  Q12XY (NLON,NLAT),    Q34XY (NLON,NLAT),
     *  FOG12XY (NLON,NLAT),  FOG34XY (NLON,NLAT),
     *  HFOG12XY (NLON,NLAT), HFOG34XY (NLON,NLAT)
c
c----------------------------------------------------------------------
c
      COMMON /VEGLSXXY/
     *  IVEGXY (NLON,NLAT),
     *  FLXY (NLON,NLAT),               FUXY (NLON,NLAT),
     *  LAIXY (NLON,NLAT,2),            SAIXY (NLON,NLAT,2), 
c
     *  RHOVEGXY (NLON,NLAT,NBAND,2),   TAUVEGXY (NLON,NLAT,NBAND,2),
     *  ORIEVXY (NLON,NLAT,2),          ORIEHXY (NLON,NLAT,2),
c
     *  ZTOPXY (NLON,NLAT,2),           ZBOTXY (NLON,NLAT,2),
     *  DLEAFXY (NLON,NLAT,2),          DSTEMXY (NLON,NLAT,2),
     *  FSTOMAXY (NLON,NLAT,2),
c
     *  RSMINXY (NLON,NLAT,2),          RSMAXXY (NLON,NLAT,2),
     *  PARCONXY (NLON,NLAT,2),         VPDCONXY (NLON,NLAT,2),
     *  TEMCONXY (NLON,NLAT,2),         TEMOPTXY (NLON,NLAT,2),
     *  FROOTXY (NLON,NLAT,NSOILAY,2),
     *  IBIOMEXY(NLON,NLAT)                   !  biome # for BIOME model
c
      REAL LAIXY
c
c     For vegdat_eve in vegdat.f (so FLXY,FUXY must be contiguous):
c
      DIMENSION FALLXY(NLON,NLAT,2)
      EQUIVALENCE (FALLXY, FLXY)
c
c----------------------------------------------------------------------
c
#ifdef ibis
c     For IBIS photosynthesis (see comments in comveg for equiv 1D arrs)

      COMMON /CIBISXY/ 
     *  VMAXSUNUXY(NLON,NLAT,2),        VMAXSUNLXY(NLON,NLAT,2),
     *  VMAXSHAUXY(NLON,NLAT,2),        VMAXSHALXY(NLON,NLAT,2),

     *  CISUNUDXY(NLON,NLAT),           CISUNUCXY(NLON,NLAT),
     *  CISUNL3XY(NLON,NLAT),           CISUNL4XY(NLON,NLAT),
     *  CISHAUDXY(NLON,NLAT),           CISHAUCXY(NLON,NLAT),
     *  CISHAL3XY(NLON,NLAT),           CISHAL4XY(NLON,NLAT),

     *  CSSUNUDXY(NLON,NLAT),           CSSUNUCXY(NLON,NLAT),
     *  CSSUNL3XY(NLON,NLAT),           CSSUNL4XY(NLON,NLAT),
     *  CSSHAUDXY(NLON,NLAT),           CSSHAUCXY(NLON,NLAT),
     *  CSSHAL3XY(NLON,NLAT),           CSSHAL4XY(NLON,NLAT),

     *  C4FRACXY(NLON,NLAT),            CONFRACXY(NLON,NLAT)

c        For IBIS photosynthesis/stomatal conductance history fields:
      COMMON /CIBISHXY/
     *  ANUDXY(NLON,NLAT),              ANUCXY(NLON,NLAT),
     *  ANL3XY(NLON,NLAT),              ANL4XY(NLON,NLAT),
     *  AGUDXY(NLON,NLAT),              AGUCXY(NLON,NLAT),
     *  AGL3XY(NLON,NLAT),              AGL4XY(NLON,NLAT),
     *  GSUDXY(NLON,NLAT),              GSUCXY(NLON,NLAT),
     *  GSL3XY(NLON,NLAT),              GSL4XY(NLON,NLAT)

c     For Wisconsin veg and dynamic IBIS phenology
c     (no corresponding 1-D arrays)
c        T15RUNXY    = 15-day exponential running mean 2-m air temp (K)
c        AGDD0XY     = degree days above 0 C based on t15runxy,
c                      reset to 0 when T15RUNXY drops below threshold
c        AGDD5XY     = degree days above 5 C based on t15runxy,
c                      reset to 0 when T15RUNXY drops below threshold

      COMMON /CIBISPXY/
     *  T15RUNXY(NLON,NLAT),
     *  AGDD0XY(NLON,NLAT),             AGDD5XY(NLON,NLAT)
#endif

#ifdef ibisint
c     For IBIS dynamic vegetation (no corresp 1-D arrs)

c        NPFTU = number of vegetation types in upper canopy
c        NPFTL = number of vegetation types in lower canopy

c        FRAC[U,L}XY = fraction of upper/lower canopy for each veg type
c        PLAI[U,L]XY = "potential" one-sided LAI (before phenology)
c        CL[U,L]XY = carbon stores for leaves, all u/l   types (kg/m2)
c        CWUXY     = carbon stores for wood,   all upper types (kg/m2)
c        CR[U,L]XY = carbon stores for roots,  all u/l   types (kg/m2)
c        EXIST[U,L] = mask for climate-controlled range limits (0 or 1)
c        SAPFRAC = fraction of (upper) wood that is sapwood

c        TDAILXY     = 24-hr mean 2-m temp (K) (diagnostic only)
c        TMONXY      = monthly mean 2-m air temperature (K)
c        TCOLDXY     = minimum monthly mean 2-m temperature (K)
c        TCOLDPREVXY = TCOLDXY from previous year (K)
c        TWARMXY     = maximum monthly mean 2-m temperature (K)
c        TWARMPREVXY = TWARMXY from previous year (K)
c        ANPP[U,L]XY = annual npp, without yet adjusting for stem and
c                      root respiration, or growth respiration in subr
c                      carbon (kg/m2/yr)
c        ANPPXY      = annual npp (all veg combined) (kg/n2/yr)
c        AGPPXY      = annual gpp (all veg combined) (kg/n2/yr)
c        GDD[0,5]    = degree days above 0 or 5 deg C (C*day)
c        GDD[0,5]PREV= GDD[0,5] from previous year (C*day)
c        AVGFRESP[A,B]XY = temperature-dependence term for respiration
c                          for upper stems and top ~1m of soil (non-dim)

      PARAMETER (NPFTU = 7, NPFTL = 2)

      COMMON /CIBISIXY/ 

     *  FRACUXY(NLON,NLAT,NPFTU),
     *  PLAIUXY(NLON,NLAT,NPFTU),
     *  CLUXY(NLON,NLAT,NPFTU),
     *  CWUXY(NLON,NLAT,NPFTU),
     *  CRUXY(NLON,NLAT,NPFTU),

     *  FRACLXY(NLON,NLAT,NPFTL),
     *  PLAILXY(NLON,NLAT,NPFTL),
     *  CLLXY(NLON,NLAT,NPFTL),
     *  CRLXY(NLON,NLAT,NPFTL),

     *  EXISTU(NLON,NLAT,NPFTU),
     *  EXISTL(NLON,NLAT,NPFTL),

     *  SAPFRAC(NLON,NLAT),

     *  TDAILXY(NLON,NLAT),             TMONXY(NLON,NLAT),
     *  TCOLDXY(NLON,NLAT),             TCOLDPREVXY(NLON,NLAT),
     *  TWARMXY(NLON,NLAT),             TWARMPREVXY(NLON,NLAT),
     *  ANPPUXY(NLON,NLAT,NPFTU),       ANPPLXY(NLON,NLAT,NPFTL),
     *  ANPPXY(NLON,NLAT),              AGPPXY(NLON,NLAT),
     *  GDD0(NLON,NLAT),                GDD0PREV(NLON,NLAT),
     *  GDD5(NLON,NLAT),                GDD5PREV(NLON,NLAT),
     *  AVGFRESPAXY(NLON,NLAT),         AVGFRESPBXY(NLON,NLAT)
#endif
c
