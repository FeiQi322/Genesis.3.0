c -- $Id: combudl.h,v 1.2 2001/07/03 21:29:23 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/combudl.h,v $
c -- $Name:  $

c--------
c combudl 
c--------
c
c    For lsx budget and daily means.
c
c----------------------------------------------------------------------
c
c
c        DTTOL   = timing tolerance for budget and print outputs (sec)
c
c        Following are used in BUDALL:
c        FRNETG  = net radiative energy flux from atmosphere to surface
c        FSNETG  = net sensible  energy flux from atmosphere to surface
c        FLNETG  = net latent    energy flux from atmosphere to surface
c        FMNETG  = net h2o mass flux from atmosphere to surface
c        GENETG  = net energy   flux from surface to runoff,drain,ablat
c        GMNETG  = net h2o mass flux from surface to runoff,drain,ablat
c        GMNETGA = net h2o mass flux from surface to subr wadjust
c        OENETG  = net energy   flux from surface to deep ocean
c        OMNETG  = net h2o mass flux from surface to deep ocean
c        SENETGO = net energy   storage in surface, previous BUDALL set
c        SMNETGO = net h2o mass storage in surface, previous BUDALL set
c        TIMENEO = absolute time at previous BUDALL set
c
c        Following are used in FINDPOI,BUDDO,BUDAD1/2,PRINDIAG,PRINDUMP
c        LONBUD  = global lsx long index of lsx output/dump point
c        LATBUD  = global lsx lat  index of lsx output/dump point
c
c        Following are used in BUDDO and/or BUDAD[1,2]:
c        F*IN    = fluxes across lsx external boundaries
c        F*AB    = fluxes absorbed by veg and intercepted H2O 
c        FOG*BU  = net fog fluxes (diagnostic only)
c        HEATBU  = net thermodynamic energy of veg and intercepted H2O
c        HEATBUO = previous budget value of HEATBU
c        WATBU   = net mass of intercepted H2O
c        WATBUO  = previous budget value of WATBU
c        TIMEBUO = absolute time at previous lsx budget set
c
      COMMON /BUDLSX/
c
     *  DTTOL,
c
     *  FRNETG,   FSNETG,   FLNETG,
     *  FMNETG,   GENETG,   GMNETG,   GMNETGA,
#ifdef hydmodel
     *  fmhydp,   fmhyde,   fmhydr,   fmhydi,
#endif
     *  OENETG,   OMNETG,
     *  SENETGO,  SMNETGO,  TIMENEO,
c
     *  LONBUD,   LATBUD,
c
     *  FSOLIN,   FIRIN,    FSENIN,    FLATIN,    FWATIN,
     *  FSOLAB,   FIRAB,    FSENAB,
     *  FOG12BU,  FOG34BU,
     *  HEATBU,   HEATBUO,  WATBU,     WATBUO,    TIMEBUO
c
c     For convenient writing/reading to/from restart file:
#ifdef hydmodel
      PARAMETER (NBUDGL = 16)
#else
      PARAMETER (NBUDGL = 12)
#endif
      DIMENSION TMPBUDL(NBUDGL)
      EQUIVALENCE (TMPBUDL, FRNETG)
c
c----------------------------------------------------------------------
c
c        In separate common to allow multiprocessing:
c        Set in SNOW, used in BUDAD2:
c        WBURY = lower-veg intercepted liquid lost by snow burial
c        SBURY = lower-veg intercepted snow lost by snow burial
c        TBURY = temperature of WBURY and SBURY
c        HBURY = heat lost as buried lower veg is exhumed
c
      COMMON  /BUD2LSX/
     *  WBURY(NVEC), SBURY(NVEC), TBURY(NVEC), HBURY(NVEC)
!$OMP THREADPRIVATE ( /BUD2LSX/ )
c
c----------------------------------------------------------------------
c
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
c
c     long and lat (degrees) of output point (corresp. to lonbud,latbud)
c
      common /cicedrive/
     *  alonbud, alatbud
#endif
c
