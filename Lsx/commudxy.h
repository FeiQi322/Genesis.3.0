c -- $Id: commudxy.h,v 1.2 2001/07/03 21:29:24 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/commudxy.h,v $
c -- $Name:  $

c---------
c commudxy
c---------
c
c----------------------------------------------------------------------
c
c     Global arrays for various history-file and budget-related 
c     diagnostic fields. See commud for variable descriptions.
c     *** Except POROSXY, UPTOTUXY, UPTOTLXY which do 
c     not have corresponding strip arrays in commud (POROS is in 
c     comsoi, UPTOT[U,L] are in commus).
c
      COMMON /DIALSXXY/
     *  GINVAPXY(NLON,NLAT),  GSUVAPXY(NLON,NLAT),
     *  GTRANSUXY(NLON,NLAT), GTRANSLXY(NLON,NLAT),
     *  GRUNOFXY(NLON,NLAT),  GDRAINXY(NLON,NLAT), GDRAINHXY(NLON,NLAT),
     *  WFLOXY(NLON,NLAT,NSOILAY+1),
     *  GABLATXY(NLON,NLAT),
     *  GADJUSTXY(NLON,NLAT), TS2XY(NLON,NLAT),    Z0LENXY(NLON,NLAT),
     *  GINFILXY(NLON,NLAT),  GIRRIGXY(NLON,NLAT),

     *  POROSXY(NLON,NLAT,NSOILAY),
     *  UPTOTUXY(NLON,NLAT),  UPTOTLXY(NLON,NLAT)

#if defined (rcm)
     * ,ts10xy(nlon,nlat),    qs10xy(nlon,nlat),   us10xy(nlon,nlat)
#endif

c
