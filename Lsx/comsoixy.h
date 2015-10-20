c -- $Id: comsoixy.h,v 1.2 2001/07/03 21:29:27 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comsoixy.h,v $
c -- $Name:  $

c---------
c comsoixy
c---------
c
c     Global arrays for soil model. See comsoi for variable descriptions
c
      COMMON /SOILSXXY/
     *  ITEXXY (NLON,NLAT,NSOILAY),   TSOIXY (NLON,NLAT,NSOILAY),
     *  WSOIXY (NLON,NLAT,NSOILAY),   WISOIXY (NLON,NLAT,NSOILAY),
     *  WPUDXY (NLON,NLAT),           WIPUDXY (NLON,NLAT),
     *  TGXY   (NLON,NLAT),
c
     *                                Z0SOIXY (NLON,NLAT), 
     *  ALBSAVXY (NLON,NLAT),         ALBSANXY (NLON,NLAT)
#ifdef isotrac
c       Top soil reservoir for isotope's (liquid depth equivalent) and
c       isotopic ratios and maximum capacity; same for bottom reservoir
c       (infinite capacity); vapor flux
     * ,hct_i(NLON,NLAT), rct_o(NLON,NLAT), rct_h(NLON,NLAT), hctm_i,
     *  hcb_i(NLON,NLAT), rcb_o(NLON,NLAT), rcb_h(NLON,NLAT),
     *  fvcxy_i(NLON,NLAT)
#endif
c
c----------------------------------------------------------------------
