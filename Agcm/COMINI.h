c -- $Id: COMINI.h,v 1.2 2001/07/03 21:29:07 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMINI.h,v $
c -- $Name:  $

C
CL            AGCM-GRID FIELDS PASSED TO INITIALIZE LSX FROM SCRATCH,
CL            OR PASSED BACK JUST ONCE FROM LSX AT INITIALIZATION
C
      COMMON /COMINI/
     +   YTSIN (NLON,NOREC),
     +   YOROOCN(NLON,NOREC), YOROLAND(NLON,NOREC), YOROICE(NLON,NOREC)

