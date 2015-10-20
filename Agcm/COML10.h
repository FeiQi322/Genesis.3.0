c -- $Id: COML10.h,v 1.2 2001/07/03 21:29:08 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COML10.h,v $
c -- $Name:  $

C
CL    THIS COMMON BLOCK REPLACES OLD NLINE1
C
C     DYNAMICAL VARIABLES
C
      COMMON /COML10/
     +       YFI(NLON,NOREC),             YUM1(NLON,NOREC,NLEV),
     +       YVM1(NLON,NOREC,NLEV),       YLPSM1(NLON,NOREC),
     +       YTM1(NLON,NOREC,NLEV),
     +       YZM1(NLON,NOREC,NLEV),       YDM1(NLON,NOREC,NLEV),
     +       YQM1(NLON,NOREC,NLEV)
#ifdef isotrac
c    Isotopic and rtr specific humidity for the vapor phase
     +      ,yqm1_o(nlon,norec,nlev),     yqm1_h(nlon,norec,nlev),
     +       rtrqm1(nlon,norec,nlev,nrtr)
#endif

