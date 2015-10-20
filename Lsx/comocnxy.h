c -- $Id: comocnxy.h,v 1.2 2001/07/03 21:29:26 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comocnxy.h,v $
c -- $Name:  $

c---------
c comocnxy
c---------
c
c     Global arrays for ocean-surface-flux model. See comocn for
c     variable descriptions and other scalars, except:
c
c     SOCNXY = surface ogcm salinity (fraction,for history,no 1-D array)
c     TFREEZXY = surface freezing pt (K, dep on salinity)
c
c----------------------------------------------------------------------
c
      COMMON /OCNLSXXY/
     *  TOCNXY(NLON,NLAT),    HOCNXY(NLON,NLAT),
     *  FTOCNXY(NLON,NLAT),   FWOCNXY(NLON,NLAT),
     *  FUOCNXY(NLON,NLAT),   FVOCNXY(NLON,NLAT),
     *  QFLUXXY(NLON,NLAT),  
     *  SOCNXY(NLON,NLAT),    TFREEZXY(NLON,NLAT),
c       for mom2,4:
     *  FHOCNXY (NLON,NLAT),   FQOCNXY (NLON,NLAT),
     *  FSOLOCNXY(NLON,NLAT),  FIROCNXY(NLON,NLAT),
     *  FPRLOCNXY (NLON,NLAT), FPRFOCNXY(NLON,NLAT),
     *  FSALOCNXY (NLON,NLAT), PRESOCNXY(NLON,NLAT)

#ifdef isotrac
c       See comocn.h for description. vapor flux from ocean and lakes;
c       isotopic ratios for ocean and lake surface water; capacity of
c       lake reservoir, maximum for deep and shallow lakes
     * ,xoxy_o(NLON,NLAT),  xoxy_h(NLON,NLAT),
     *  qoxy_i(NLON,NLAT),  zfxy_i(NLON,NLAT),
     *  fvoxy_i(NLON,NLAT), ro_o(NLON,NLAT),   ro_h(NLON,NLAT),
     *  hl_i(NLON,NLAT),    hldm_i,            hlsm_i
#endif
c
c----------------------------------------------------------------------
c
c     Global arrays for OGCM coupling.
c     nb: F[U,V]OCNET are time filtered in dynice 
c
c        UOCNXY  = eastward  surface velocity (m/s)
c        VOCNXY  = northward surface velocity (m/s)
c        FUOCNET = net eastward  stress on ocean (N m-2)
c        FVOCNET = net northward stress on ocean (N m-2)
c
      COMMON /OGCMLSX/
     *  UOCNXY (NLON,NLAT),    VOCNXY (NLON,NLAT),
     *  FUOCNET(NLON,NLAT),    FVOCNET(NLON,NLAT)
c
