c -- $Id: comsoi.h,v 1.2 2001/07/03 21:29:27 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comsoi.h,v $
c -- $Name:  $

c-------
c comsoi
c-------

c----------------------------------------------------------------------

c     soil model prescribed constants, set in INISOD or INISOI
c     and Namelist soil params

c        HSOD     = soil or land-ice layer thicknesses vs nsoilay only
c                   (see comments in setsoi)
c        WSOI_INIT = initial soil moisture

      COMMON /SODLSX/
     *  HSOD (NSOILAY), WSOI_INIT

c----------------------------------------------------------------------

c     1D arrays for soil model. These 1D arrays correspond to global
c     arrays in comsoixy, except hsoi (see comments in setsoi)
c
c        ITEX     = 100*sand fraction + clay fraction (prescribed),
c                   or: -2 = peat, -1 = land-ice, 0 = ocean/seaice
c        TSOI     = soil temperature for each soil layer
c        WSOI     = liquid content of soil, relative to pore space 
c                   remaining AFTER ice content
c        WISOI    = ice content of soil, relative to total pore space
c        WPUD     = liquid content of puddles (Kg m-2) per soil area
c        WIPUD    = ice    content of puddles (Kg m-2) per soil area
c        TG       = skin temperature of soil
c
c        HSOI    = soil or land-ice layer thicknesses (m)
c        Z0SOI   = roughness length of soil/landice surface (m)
c        ALBSAV  =saturated soil/landice surface albedo(visible wavebnd)
c        ALBSAN  =saturated soil/landice surface albedo(near-ir wavebnd)

      COMMON  /SOILSX/
     *  ITEX(NVEC,NSOILAY),         TSOI(NVEC,NSOILAY),
     *  WSOI(NVEC,NSOILAY),         WISOI(NVEC,NSOILAY),
     *  WPUD(NVEC),                 WIPUD(NVEC),
     *  TG(NVEC),
c
     *  HSOI (NVEC,NSOILAY),        Z0SOI (NVEC), 
     *  ALBSAV (NVEC),              ALBSAN (NVEC)
!$OMP THREADPRIVATE ( /SOILSX/ )

c----------------------------------------------------------------------

c     Soil 1-D arrays, set in SETSOI from texture database
c     (no correspondig 2-D arrays)

c        HYDRAUL = saturated hydraulic conductivity (m/s)
c        SUCTION = saturated matric potential (m)
c        BEX     = Clapp-Hornberger exponent "B" in soil water potential
c        IBEX    = nint(bex), used for cpu speed
c        SWILT   = wilting soil moisture value (must be > 0, < 1)
c        UPCON   = const in soil water uptake rate param (Kg m-2 s-1)
c        CSOI    = specific heat of soil, no pore spaces (J Kg-1 deg-1)
c        RHOSOI  = soil density (without pores, not bulk) (Kg m-3)
c        POROS   = porosity (mass of water per unit vol at sat / rhow)
c        POROSFLO= porosity after reduction by ice content
c        FSAND   = fraction of sand (unpacked from ITEX)
c        FSILT   = fraction of silt (unpacked from ITEX)
c        FCLAY   = fraction of clay (unpacked from ITEX)
c
      COMMON  /SOTLSX/
     *  HYDRAUL (NVEC,NSOILAY),          SUCTION (NVEC,NSOILAY),
     *  BEX (NVEC,NSOILAY),              IBEX (NVEC,NSOILAY),
     *  SWILT (NVEC,NSOILAY),            UPCON (NVEC),
     *  CSOI (NVEC,NSOILAY),             RHOSOI (NVEC,NSOILAY),
     *  POROS (NVEC,NSOILAY),            POROSFLO (NVEC,NSOILAY),
     *  FSAND (NVEC,NSOILAY),            FSILT (NVEC,NSOILAY),
     *  FCLAY (NVEC,NSOILAY)
!$OMP THREADPRIVATE ( /SOTLSX/ )
c
c----------------------------------------------------------------------
c
c     Soil look-up table, set in claphorn:
      PARAMETER (NHYD=4, NTEX=51)

      COMMON /CHYDVAL/
     *  HYDVAL(NTEX,NTEX,NHYD)
c
