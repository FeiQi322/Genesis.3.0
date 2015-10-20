c -- $Id: COMMUN.h,v 1.2 2001/07/03 21:29:09 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMMUN.h,v $
c -- $Name:  $

C
CL    THESE COMMON BLOCKS ARE FOR COMMUNICATION BETWEEN SEPARATE PARTS
CL    OF THE CODE, MOSTLY WITHIN ONE TIMESTEP OR BETWEEN FRADSW,FRADLW
CL    CALCULATIONS, BUT SOME NEEDED ACROSS EVERY TIMESTEP BOUNDARY,
CL    EG, PRECIP, SLT FIRST GUESS, SO MUST BE WRITTEN TO RESTART FILE
CL    (SEE RSDS,WSDS).
C
      COMMON /COMMUN1/
C
C       RADIATION HEATING RATES, INSOLATION PARAMS
C       (COMPUTED ONLY IF FRADSW = .TRUE. (EXCEPT YQRL))
     +  YQRS(NLON,NOREC,NLEVP),         YQRL(NLON,NOREC,NLEVP),
     +  COSQ(NLON,NOREC),               FRAQ(NLON,NOREC),
C
C       SURFACE PRESSURE, VERTICAL VELOCITY, SURFACE PRECIP RATES
     +  YPSM1(NLON,NOREC),              YPS(NLON,NOREC),
     +  TOPOGUN(NLON,NOREC),            VOMEGA(NLON,NOREC,NLEV),       
     +  WPLUME0(NOREC,NLEV),            WPLUME1(NOREC,NLEV),
     +  FRAINC(NLON,NOREC),             FRAINS(NLON,NOREC), 
     +  FSNOWC(NLON,NOREC),             FSNOWS(NLON,NOREC), 
C
C       FOR COSMOGENIC ISOTOPES: PRECIPITATION FLUX THROUGH LEVEL
C       INTERFACES (kg/m2/s)
#ifdef cosmog
     +  PRECFLX(NLON,NOREC,NLEVP),
#endif
C
C       HEATING RATE AND MOMENTUM TENDENCIES
     +  YT2(NLON,NOREC,NLEV),           
     +  YFU(NLON,NOREC,NLEV),           YFV(NLON,NOREC,NLEV)
#ifdef isotrac
c    Isotopic precipitation rates for convective and stable rainfall
c    and snowfall. Idem for rtr.
     + ,frainc_o(nlon,norec),           frainc_h(nlon,norec),
     +  frains_o(nlon,norec),           frains_h(nlon,norec),
     +  fsnowc_o(nlon,norec),           fsnowc_h(nlon,norec),
     +  fsnows_o(nlon,norec),           fsnows_h(nlon,norec),
     +  rtrprc(nlon,norec,nrtr),        rtrprs(nlon,norec,nrtr),
     +  rtrpsc(nlon,norec,nrtr),        rtrpss(nlon,norec,nrtr)
#endif
C
      COMMON /COMMUN2/
C
C       PASSED TO AND FROM LSX MODULE
     +  ALBSUD(NLON,NOREC,2),           ALBSUI(NLON,NOREC,2),
     +  SOLSUD(NLON,NOREC,2),           SOLSUI(NLON,NOREC,2),
     +  FIRABS(NLON,NOREC),             FIREMS(NLON,NOREC),
     +  YUFLX(NLON,NOREC),              YVFLX(NLON,NOREC),
     +  YTFLX(NLON,NOREC),              YQFLX(NLON,NOREC),
     +  SLDFLX(NLON,NOREC),

#if defined (LKAEROSOL)
     +  SFCIRRAD(NLON,NOREC),
#endif
C
C       PASSED TO AND FROM SLT (PREVIOUS DEPARTURE POINTS,
C       RUNNING-MEAN VERTIAL VELOCITY,  FOR SLT)
     +  ALONSAV(NLON,NOREC,NLEV),       ALATSAV(NLON,NOREC,NLEV),
     +  SIGSAV(NLON,NOREC,NLEV),        WSLTSAV(NLON,NOREC,NLEV)

#ifdef isotrac
c    Isotopic and rtr evaporation/condensation flux from the surface
     + ,yqflx_o(nlon,norec),            yqflx_h(nlon,norec),
     +  rtrfv (nlon,norec,nrtr),
c    Sources spatial distribution for regional tracers
     +  sourc_i(nlon,norec)
#endif

