c -- $Id: comsnoxy.h,v 1.2 2001/07/03 21:29:26 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comsnoxy.h,v $
c -- $Name:  $

c---------
c comsnoxy
c---------
c
c     Global arrays for snow model.
c     See comsno for variable descriptions, and snow scalars.
c
c----------------------------------------------------------------------
c
      COMMON /SNOLSXXY/
     *  TSNOXY (NLON,NLAT,NSNOLAY),    HSNOXY (NLON,NLAT,NSNOLAY),
     *  FIXY (NLON,NLAT),              TIXY (NLON,NLAT)
c
c----------------------------------------------------------------------
c
c     Global arrays for snow on sea ice. Uses same 1D arrays in comsno
c     and scalars as for snow on land or icesheet.
c
      COMMON /SNOLSXXY/
     *  TSNO_ICEXY(NLON,NLAT,NSNOLAY), HSNO_ICEXY(NLON,NLAT,NSNOLAY),
     *  FI_ICEXY(NLON,NLAT),           TI_ICEXY(NLON,NLAT)
c
