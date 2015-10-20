c -- $Id: comeve.h,v 1.2 2001/07/03 21:29:23 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comeve.h,v $
c -- $Name:  $

c-------
c comeve
c-------
c
c     Global arrays for daily and monthly mean T, precip, r.h.,
c     to be accumulated for interactive EVE model.
c
      common /evelsxxy/
     *  teve_d (nlon,nlat),  teve_m (nlon,nlat),
     *  peve_d (nlon,nlat),  peve_m (nlon,nlat),
     *  reve_d (nlon,nlat),  reve_m (nlon,nlat),
c
c     Global arrays for monthly mean T, precip, and sunshine hours
c     to be accumulated for interactive BIOME model.
c
     *  tbiome (nlon,nlat), pbiome (nlon,nlat),
     *  cbiome (nlon,nlat)

