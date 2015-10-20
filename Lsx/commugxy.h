c -- $Id: commugxy.h,v 1.2 2001/07/03 21:29:25 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/commugxy.h,v $
c -- $Name:  $

c---------
c commugxy
c---------
c
c     Global arrays for communicating between lsx and agcm.
c     See commug for variable descriptions (except delevxy).
c
c     DELEVXY = elevation correction field (m) for adjusting taxy 
c               and firaxy in surfctl, due to subgrid or agcm-filtered 
c               topography (difference between "true" topog on LSX
c               grid read from data file in inilsx, and passed
c               gcm topog interpolated to lsx grid).       
c
c     TOTCLDXY = total cloudiness (used only for delevtype,and biome3/4)
c
c     TOPNODXY = uncorrected (no DELEVXY) topography passed from gcm 
c                and interpolated to lsx grid, for history only)  
c
c----------------------------------------------------------------------
c
      COMMON /GCMLSXXY/
     *  COSZENXY (NLON,NLAT),
     *  SOLADXY (NLON,NLAT,NBAND),    SOLAIXY (NLON,NLAT,NBAND), 
     *  ASURDXY (NLON,NLAT,NBAND),    ASURIXY (NLON,NLAT,NBAND),
c 
     *  FIRAXY (NLON,NLAT),
     *  UXXY (NLON,NLAT),             UYXY (NLON,NLAT),  
     *  TAXY (NLON,NLAT),             QAXY (NLON,NLAT),
     *  RAINAXY (NLON,NLAT),          SNOWAXY (NLON,NLAT),
     *  PSURFXY (NLON,NLAT), 
c
     *  AFIREXY (NLON,NLAT),
     *  ASENAXY (NLON,NLAT),          AVAPAXY (NLON,NLAT),
     *  ATAUXXY (NLON,NLAT),          ATAUYXY (NLON,NLAT),
c
     *  LMASKXY (NLON,NLAT),          FWATERXY (NLON,NLAT),
c
     *  DELEVXY (NLON,NLAT),          TOTCLDXY (NLON,NLAT),
     *  TOPNODXY(NLON,NLAT)
#ifdef isotrac
c       Isotopic specific humidity at the lowest GCM level; isotopic
c       rainfall and snowfall rates; isotopic evaporation/condensation
c       flux to the surface; ocean/land mask for isotopes
     * ,qaxy_o(NLON,NLAT),            qaxy_h(NLON,NLAT),
     *  rainaxy_o(NLON,NLAT),         rainaxy_h(NLON,NLAT),
     *  snowaxy_o(NLON,NLAT),         snowaxy_h(NLON,NLAT),
     *  avapaxy_o(NLON,NLAT),         avapaxy_h(NLON,NLAT),
     *  mocean_i(NLON,NLAT)
#endif


#ifdef rcm
c       Global RCM arrays for easily transfering veg/land map and
c       topography map to where they are used (in inilsx/vegdat_dorm)
c
      COMMON /RCMLSXXY/
     *  VEGRCM(NLON,NLAT),            TOPRCM(NLON,NLAT),
     *  ISURFTYPE

#elif hydmodel

c       Global hydrol arrays for easily transfering veg/land map and
c       topography map to where they are used (in inilsx/vegdat_dorm/
c       inisoi). Also specified soil moistures for categories
c       (thetahyd), their uppermost soil-layer numbers (nlayhyd),
c       standing water depths if any (hwaterhyd), and number
c       of categories (ncateghyd, must be <= 100).
c
      COMMON /HYDLSXXY/
     *  LMASKHYD(NLON,NLAT),          TOPHYD(NLON,NLAT),
     *  IVEGHYD(NLON,NLAT),           ISOILHYD(NLON,NLAT,NSOILAY),
     *  THETAHYD(100),                NLAYHYD(100),
     *  HWATERHYD(100),               NCATEGHYD
#endif
