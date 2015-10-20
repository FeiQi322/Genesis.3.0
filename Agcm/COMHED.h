c -- $Id: COMHED.h,v 1.2 2001/07/03 21:29:06 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMHED.h,v $
c -- $Name:  $

C
CL            RUN TITLE AND VERSION NUMBERS
C
      COMMON /COMHED/
     +   VERSGEN, 
     +   VERSHIS,   TYPEHIS,
     +   VERSRES,   TYPERES,
     +   VERSRESIN, TYPERESIN,
     +   VERSRES_ISOT, VERSRESIN_ISOT

      CHARACTER*8 VERSGEN, TYPEHIS, TYPERES, TYPERESIN

      COMMON /COMHEDC/
     +   RUNTITLE
      CHARACTER*80 RUNTITLE

