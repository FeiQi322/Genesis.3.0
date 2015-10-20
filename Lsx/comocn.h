c -- $Id: comocn.h,v 1.2 2001/07/03 21:29:25 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comocn.h,v $
c -- $Name:  $

c-------
c comocn
c-------
c
c     Variables for ocean surface flux model.
c     1D arrays correspond to global arrays in comocnxy.
c
c----------------------------------------------------------------------
c
c        TOCN  = mixed-layer temperature
c        HOCN  = mixed-layer depth
c        FTOCN = net heat flux into overall ocn sfc (open + under ice)
c        FWOCN = net H2O  flux into overall ocn sfc (open + under ice)
c        FUOCN = eastwd wind stress imparted from atm to open ocn sfc
c        FVOCN = northwd wind stress imparted from atm to open ocn sfc
c        QFLUX = ocean heat flux convergence (W/m2)
c    *           (prescribed in version 1, diffusive in version 2)
c        SOCN  = mixed-layer salinity (fraction)
c        TFREEZ= mixed-layer freezing point (K, dep on salinity)

c       for mom2,4 (also fwocn,fuocn,fvocn):
c        FHOCN = upward sensible heat flux at ocn sfc (open + under ice)
c        FSOLOCN= net dwnwd solar rad flux at ocn sfc (open + under ice)
c        FIROCN = net dwnwd ir rad flux at ocn sfc (open + under ice)
c        FPRLOCN= liquid precip at ocn sfc (open + under ice)
c        FPRFOCN= frozen precip at ocn sfc (open + under ice)
c        FSALOCN= salt flux (see comments in genmom2,4) (open+underice)
c        PRESOCN = hydrostatic pressure due to atmos + sea ice on ocean
c
      COMMON  /OCNLSX/
     *  TOCN (NVEC),         HOCN (NVEC),
     *  FTOCN (NVEC),        FWOCN (NVEC),
     *  FUOCN (NVEC),        FVOCN (NVEC),
     *  QFLUX (NVEC),
     *  SOCN  (NVEC),        TFREEZ(NVEC),
c       for mom2,4:
     *  FHOCN   (NVEC),       FQOCN(NVEC),
     *  FSOLOCN (NVEC),       FIROCN(NVEC),
     *  FPRLOCN (NVEC),       FPRFOCN(NVEC),
     *  FSALOCN (NVEC),       PRESOCN(NVEC)

#ifdef isotrac
c       Exponent for isotopic kinetic effect, saturation spec. humidity
c       at ocean surface and aerodynamic coefficient for ocean water
c       evaporation/condensation
     * ,xo_o(NVEC), xo_h(NVEC), qo_i(NVEC), zf_i(NVEC)
#endif
!$OMP THREADPRIVATE ( /OCNLSX/ )
c
c----------------------------------------------------------------------
c
c        In separate common to allow multiprocessing:
c        Z0OCNH = roughness length for open ocean (T,q)
c                 (set in datgelsx/namelist or in iniocn)
c        HDEEP = mixed layer depth for open ocean (set in iniocn)
c        HSHALLOW = mixed layer depth for inland lakes (set in iniocn)
c
      COMMON /OCN2LSX/
     *  Z0OCNH, HDEEP, HSHALLOW
