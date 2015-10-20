c -- $Id: comsno.h,v 1.2 2001/07/03 21:29:26 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comsno.h,v $
c -- $Name:  $

c-------
c comsno
c-------
c
c     Scalars and 1D arrays for snow model. The 1D arrays correspond to
c     global arrays in comsnoxy.
c
c----------------------------------------------------------------------
c
c        Z0SNO   = roughness length of snow surface
c        RHOS    = density of snow
c        HSNOTOP = thickness of top snow layer (meters)
c        HSNOMIN = minimum total thickness of snow (meters)
c        HSNOMAX = maximum total thickness of snow (meters)
c        FIMIN   = minimum fractional snow cover
c        FIMAX   = maximum fractional snow cover
c
      COMMON /SNOLSX/
     *  Z0SNO,            RHOS,
     *  HSNOTOP,          HSNOMIN,          HSNOMAX,
     *  FIMIN,            FIMAX
c
c----------------------------------------------------------------------
c
c        In separate common to allow multiprocessing:
c        TSNO    = temperature of snow layers
c        HSNO    = thickness of snow layers (meters)
c        FI      = fractional snow cover
c        TI      = skin temperature of snow
c
      COMMON  /SNO2LSX/
     *  TSNO (NVEC,NSNOLAY),  HSNO (NVEC,NSNOLAY),
     *  FI (NVEC),            TI(NVEC)
!$OMP THREADPRIVATE ( /SNO2LSX/ )
c
