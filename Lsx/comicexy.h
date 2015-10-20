c -- $Id: comicexy.h,v 1.2 2001/07/03 21:29:24 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comicexy.h,v $
c -- $Name:  $

c---------
c comicexy
c---------
c
c     Global arrays for sea-ice model. See comice for variable descrips.
c
      COMMON /SILSXXY/
     *  HICEXY (NLON,NLAT,NICELAY),   TICEXY (NLON,NLAT,NICELAY),
     *  FICEXY (NLON,NLAT),           BRINEXY(NLON,NLAT),
     *  TG_ICEXY(NLON,NLAT)
#ifdef isotrac
c       Sea-ice surface reservoir (snow) for isotope's (liquid equiv.
c       depth), max capacity and isotopic ratios
     * ,hi_i(NLON,NLAT), him_i, ri_o(NLON,NLAT), ri_h(NLON,NLAT),
c       see comice
     *  ficexy_i(NLON,NLAT), fvixy_i(NLON,NLAT)
#endif
c
c----------------------------------------------------------------------
c
c     Global arrays for sea-ice (thermodynamic and dynamic) model.
c     See comice for variable descriptions and other scalars.
c
      COMMON /SI2LSXXY/
     *  FUICEXY(NLON,NLAT),   FVICEXY(NLON,NLAT)
c
c     Global arrays for running time means of wind stress on sea ice and
c     on open ocean, for dynamic sea-ice model
c
      COMMON /SIRUNXY/
     *  FUICEXYM(NLON,NLAT),  FVICEXYM(NLON,NLAT),
     *  FUOCNXYM(NLON,NLAT),  FVOCNXYM(NLON,NLAT)
c
