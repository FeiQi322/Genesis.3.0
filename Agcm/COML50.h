c -- $Id: COML50.h,v 1.2 2001/07/03 21:29:08 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COML50.h,v $
c -- $Name:  $

C
CL    THIS COMMON BLOCK IS FOR NON-LEAPFROG PROGNOSTIC VARIABLES
CL    (LIQUID WATER CONTENTS EQUIVALENCED TO TRACERS 1-3,
CL    FLUX RESERVOIRS EQUIVALENCED TO TRACERS 4-7; SEE NTRACA ABOVE),
CL    AND GLOBAL DIAGNOSTIC TRACER QUANTITIES TOTT*.
C
      COMMON /COML50/
     +       TRACER (NLON,NOREC,NLEV,NTRACE),
     +       TRAFLUX(NLON,NOREC,NTRACB),
     +       TOTTO(NTRACE), TOTTN(NTRACE), TOTTA(NTRACE)
#ifdef isotrac
c    Isotopic liquid water contents for convective, anvil-cirrus and
c    stratus clouds / isotopic and rtr vapor flux reservoir
     +      ,ylwc_c_o(nlon,norec,nlev),      ylwc_c_h(nlon,norec,nlev),
     +       ylwc_a_o(nlon,norec,nlev),      ylwc_a_h(nlon,norec,nlev),
     +       ylwc_s_o(nlon,norec,nlev),      ylwc_s_h(nlon,norec,nlev),
     +       fplumq_o(nlon,norec,nlev),      fplumq_h(nlon,norec,nlev),
     +       rtrfq (nlon,norec,nlev,nrtr)
#endif
      DIMENSION
     +       YLWC_C(NLON,NOREC,NLEV),
     +       YLWC_A(NLON,NOREC,NLEV),
     +       YLWC_S(NLON,NOREC,NLEV),
     +       FPLUMT(NLON,NOREC,NLEV),        FPLUMQ(NLON,NOREC,NLEV),
     +       FPLUMU(NLON,NOREC,NLEV),        FPLUMV(NLON,NOREC,NLEV)
      EQUIVALENCE 
     +       (TRACER(1,1,1,1), YLWC_C),
     +       (TRACER(1,1,1,2), YLWC_A),
     +       (TRACER(1,1,1,3), YLWC_S),
     +       (TRACER(1,1,1,4), FPLUMT),
     +       (TRACER(1,1,1,5), FPLUMQ),
     +       (TRACER(1,1,1,6), FPLUMU),
     +       (TRACER(1,1,1,7), FPLUMV)

