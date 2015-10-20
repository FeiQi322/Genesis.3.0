c -- $Id: COMMAP.h,v 1.2 2001/07/03 21:29:08 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMMAP.h,v $
c -- $Name:  $

C
CL            MAP FACTORS AND HYDROSTATIC AND ENERGY CONVERSION MATRICES
C
      COMMON /COMMAP/
     +   SIG(NLEV),     SIGKMH(NLEVP), DSIGMA(NLEV),
     +   SIGKPK(NLEV),  SIGKPL(NLEVP), R2DSIG(NLEV),
     +   SQ(NPMAX),     RSQ(NPMAX),    XM(MMAX),
     +   SIT(NOREC),    W(NOREC),      CS(NOREC),     RCS(NOREC),
     +   ALATG(NOREC),  ALONG(NLON),
     +   G(NLEV2),      AG(NLEV2),     BHYD(NLEV2)

