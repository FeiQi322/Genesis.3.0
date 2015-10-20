c -- $Id: COMHIY.h,v 1.2 2001/07/03 21:29:07 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMHIY.h,v $
c -- $Name:  $

C
CL            TEMPORARY ARRAYS FOR STORING HISTORY RADIATION QUANTITIES
CL            LATITUDE BY LATITUDE, AND THEN TRANSMITTING THEM TO THE
CL            HISTORY STORAGE AREA.
C
      COMMON /COMHIY/
     +   HFRSA(NLON,NOREC),   HSABTP(NLON,NOREC),  HSABAT(NLON,NOREC),
     +   HCLRSS(NLON,NOREC),  HCLRST(NLON,NOREC),  HSOLIN(NLON,NOREC),
     +   HFRLA(NLON,NOREC),   HFIRTP(NLON,NOREC),  HCLRLS(NLON,NOREC),
     +   HCLRLT(NLON,NOREC),  HFRSAC(NLON,NOREC)

