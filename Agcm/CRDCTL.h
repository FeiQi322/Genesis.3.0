c -- $Id: CRDCTL.h,v 1.2 2001/07/03 21:29:12 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/CRDCTL.h,v $
c -- $Name:  $

C
CL            RADIATION CONTROL VARIABLES
C
C             FRADSW   = .T. IFF FULL SHORTWAVE COMPUTATION
C             FRADLW   = .T. IFF FULL LONGWAVE COMPUTATION
C
      COMMON /CRDCTL/
     +   TRADLW,    TRADLWA,    TRADSW,
     +   IRADLW,    IRADLWA,    IRADSW,
     +   FRADLW,    FRADLWA,    FRADSW
      LOGICAL 
     +   FRADLW,    FRADLWA,    FRADSW

