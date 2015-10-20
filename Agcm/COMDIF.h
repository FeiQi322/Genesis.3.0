c -- $Id: COMDIF.h,v 1.2 2001/07/03 21:29:05 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMDIF.h,v $
c -- $Name:  $

C
CL            HORIZONTAL DIFFUSION 
C
      COMMON /COMDIF/
     +   DIF2,               DIF4,
     +   KDEL4,              KDEL2INC,    
     +   KCOURANT,           SCOURANT,
     +   HDIFT2(NPMAX),      HDIFZ2(NPMAX),
     +   HDIFT4(NPMAX),      HDIFZ4(NPMAX),
     +   HDIFT(NPMAX,NLEV),  HDIFZ(NPMAX,NLEV),
     +   HFACT(NPMAX,NLEV),  HFACZ(NPMAX,NLEV)

