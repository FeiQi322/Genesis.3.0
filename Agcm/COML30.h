c -- $Id: COML30.h,v 1.2 2001/07/03 21:29:08 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COML30.h,v $
c -- $Name:  $

C
CL    THIS COMMON BLOCK REPLACES OLD NLINE3
C
      COMMON /COML30/
     +       YF(NLON,NOREC),              YU(NLON,NOREC,NLEV),
     +       YV(NLON,NOREC,NLEV),         YLPS(NLON,NOREC),
     +       YT(NLON,NOREC,NLEV),
     +       YZ(NLON,NOREC,NLEV),         YD(NLON,NOREC,NLEV),
     +       YQ(NLON,NOREC,NLEV),
     +       YDUH(NLON,NOREC,NLEV),       YDVH(NLON,NOREC,NLEV),
     +       YDTH(NLON,NOREC,NLEV),
     +       YDPSL(NLON,NOREC),           YDPSM(NLON,NOREC)
#ifdef isotrac
c    Isotopic and rtr specific humidity for the vapor phase
     +      ,yq_o(nlon,norec,nlev),       yq_h(nlon,norec,nlev),
     +       rtrq (nlon,norec,nlev,nrtr)
#endif

