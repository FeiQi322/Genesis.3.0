c -- $Id: COMSTA.h,v 1.2 2001/07/03 21:29:10 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMSTA.h,v $
c -- $Name:  $

C
CL            STATISTICS INTEGRALS
C
      COMMON /COMSTA/
     +   RMSZ,     RMSD,     RMST,     STPS,     STPE,
     +   STKE,     STTE,
     +   VZ2(NLEV),D2(NLEV), T2(NLEV), PE(NLEV), PT(NLEV),
     +   PS,       PFI

