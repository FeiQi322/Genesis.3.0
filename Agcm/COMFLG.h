c -- $Id: COMFLG.h,v 1.2 2001/07/03 21:29:06 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMFLG.h,v $
c -- $Name:  $

C
CL            LOGICAL CONTROL FLAGS AND OTHER USER-SETTABLE VALUES
C
      COMMON /COMFLG/
     +   NSREST,    NSREST_ISOT, TEQCOLD,   TPOCOLD,    PERTURB,
     +   PRIESTQ,   STRATTIME, STRATFLAG,  SMOLAR,   TRAVDIF,  TRACONV,
     +   VEGTYPE,   EVEDELTA, DSBIOME,   WISCVEG,    PHOTOTYPE,TAUBIOME,
     +   DYNAIBIS,  ASYNIBIS,  NSREST_IBIS,SOILTEX,   
     +   CO2PPM,    SOLFAC,     SOLPAR,
     +   STOCHPR,   SHOWMAPS,  LONOLSX,    LATOLSX,   TIMOLSX,
     +   DEPTHML,   DEPTHLAKE, DEPTHICE,   DEPTHSNO,
     +   OCEANTYPE, SEQSST,    DYNAMICE,   DYNAMHEM,  PRESCOUV,
     +   QFACTOR,   QNORWEG,   QICES,      QICEN,     QLIMIT,

     +   TRIASSIC,  HYBTYPE,   GRAVFAC,    GWROUGH,   ENVFAC,
     +   DELEVTYPE, LEVRES,    HORRES(2),  HORRESLSX(2), HOTLIMIT,
     +   TIMRESERV,
     +   PLUMERAD1, PLUMEAREA1_L, PLUMEAREA1_O, PLUMERAD2, PLUMEAREA2,

     +   ALW,
     +   BLW_C,     BLW_A,     BLW_S,
     +   CLWFAC_C,  CLWFAC_A,  CLWFAC_S, 
     +   CLW_P,
     +   DLW_LL,    DLW_OL,    DLW_LI,    DLW_OI,
     +   DLWFAC,    REFFAC_L,  REFFAC_O,  REFFAC_I, RADLIQ,
     +   ELW, 
     +   FLW_C,     FLW_A,     FLW_S,
     +   GLW,
     +   TLIQA,     TLIQB,     SIGANVIL,
     +   AERTAUVIS, MXAERL,    GGUNIF,    NOCLDTOP,  NOCLDBOT, 
     +   FACCLWP_SOL,  FACCLWP_IR,
     +   SSACLD(2)

#ifdef isotrac
c    Liquid-solid transition temperature for condensation, slope of the
c    supersaturation function used to determine the kinetic
c    fractionation for solid condensate formation, intercept and slope
c    of the linear function used in the determination of the effective
c    relative humidity for precip. evaporation+exchange, factor used in
c    the computation of the mean drop radius for precip.,
c      For stratiform and anvil clouds:
     +  ,t_s_i,     s_s_i,     hi_s_i,    hs_s_i,    r_s_i,
c      For convective clouds:
     +   t_c_i,     s_c_i,     hi_c_i,    hs_c_i,    r_c_i,
c    Liquid-solid transition temperature for evaporation+exchange and
c    flag for liquid-vapor equilibrium forcing.
     +   t_r_i,     equ_i
#endif

      LOGICAL
     +   PRIESTQ,   SMOLAR,    TRAVDIF,   TRACONV,
     +   STOCHPR,   SHOWMAPS,
     +   SEQSST,    DYNAMICE,  DYNAMHEM,  PRESCOUV,
     +   TRIASSIC,  RADSTRIP,  GGUNIF,    VERS21

      INTEGER 
     +   STRATFLAG,
     +   VEGTYPE,  EVEDELTA,    DSBIOME,   WISCVEG,
     +   PHOTOTYPE, DYNAIBIS,ASYNIBIS,
     +   SOILTEX,
     +   OCEANTYPE, HYBTYPE,   DELEVTYPE,  HORRES,    HORRESLSX
#ifdef isotrac
     +  ,equ_i
#endif

