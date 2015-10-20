c -- $Id: COMTIM.h,v 1.2 2001/07/03 21:29:10 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMTIM.h,v $
c -- $Name:  $

C
CL            MODEL TIME VARIABLES
C
      COMMON /COMTIM/
     +   DTIME,    TWODT,
     +   NRSTRT,   NSTEP,    NSTEPH,   NSTOP,  ISTEPY,
     +   MCDATE,   MCSEC,
     +   MBDATE,   MEDATE,
     +   BASEDATE, ENDDATE
      CHARACTER*20 BASEDATE, ENDDATE

