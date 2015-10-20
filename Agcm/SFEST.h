c -- $Id: SFEST.h,v 1.2 2001/07/03 21:29:17 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/SFEST.h,v $
c -- $Name:  $

C     STATEMENT FTN USED IN SATURATION VAPOR PRESSURE TABLE LOOKUP
C
      TLIM(S) = MAX (TMIN, MIN (TMAX, S))
C
      ESTBLF(T) = ((TMIN + INT(TLIM(T)-TMIN)) + 1. - TLIM(T))
     *            * TABLE(INT(TLIM(T)-TMIN)+1)
     *          + (TLIM(T) - (TMIN + INT(TLIM(T)-TMIN)))
     *            * TABLE(INT(TLIM(T)-TMIN)+2)
C
      QSTBLF(T,P) = .622 / MAX ( P/ESTBLF(T) - (1.-.622), .622 )

