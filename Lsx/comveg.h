c -- $Id: comveg.h,v 1.2 2001/07/03 21:29:27 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comveg.h,v $
c -- $Name:  $

c-------
c comveg
c-------
c
c     1D arrays and scalars for vegetation. The 1D arrays correspond to
c     global arrays in comvegxy.
c
c----------------------------------------------------------------------
c
c     Prognostic upper and lower vegetation arrays, and those saved
c     between timesteps (T12, Q12, FOG12, HFOG12, and ...34).
c
c        TU    = temperature of upper-story leaves
c        TS    = temperature of upper-story stems
c        TL    = temperature of lower-story leaves and stems
c        TLSUB = temperature of lower-story veg buried by snow
c        WLIQU = mass of intercepted liquid H2O on upper-leaf unit area
c        WLIQS = mass of intercepted liquid H2O on upper-stem unit area
c        WLIQL = mass of intercepted liquid H2O on lower-veg unit area
c        WSNOU = mass of intercepted frozen H2O on upper-leaf unit area
c        WSNOS = mass of intercepted frozen H2O on upper-stem unit area
c        WSNOL = mass of intercepted frozen H2O on lower-veg unit area
c
c        T[12,34]    = air temperatures at Z12 and Z34
c        Q[12,34]    = specific humidity of air at Z12 and Z34
c        FOG[12,34]  = rate of fog formation (mass flux) at Z12 and Z34
c        HFOG[12,34] = latent heat of vap or subl for fog formation
c
      COMMON  /VEPLSX/
     *  TU (NVEC),           TS (NVEC),           TL (NVEC),
     *  TLSUB (NVEC),
     *  WLIQU (NVEC),        WLIQS (NVEC),        WLIQL (NVEC),
     *  WSNOU (NVEC),        WSNOS (NVEC),        WSNOL (NVEC),
c
     *  T12 (NVEC),          T34 (NVEC),
     *  Q12 (NVEC),          Q34 (NVEC),
     *  FOG12 (NVEC),        FOG34 (NVEC),
     *  HFOG12 (NVEC),       HFOG34 (NVEC)
!$OMP THREADPRIVATE ( /VEPLSX/ )
c
c----------------------------------------------------------------------
c
c     1D upper and lower vegetation arrays (not saved btwn timesteps),
c     and scalars.
c
c        SUNPU = absorbed PAR on *sunlit* upper-story leaves (W/m2)
c        SUNPL = absorbed PAR on *sunlit* lower-story leaves (W/m2)
c        SHAPU = absorbed PAR on *shaded* upper-story leaves (W/m2)
c        SHAPL = absorbed PAR on *shaded* lower-story leaves (W/m2)
c        SUNFU = fraction of sunlit LAI to total LAI for upper story
c        SUNFL = fraction of sunlit LAI to total LAI for lower story
c        S[U,S,L} = air-veg transfer coeffs for upper leaves,
c                   upper stems, lower veg.
c        For BATS-SiB stomatal resistance parameterization:
c        RTRA[U,L] = stomatal resistances of upper and lower leaves(s/m)
c
      COMMON  /VEP1LSX/
     *  SUNPU (NVEC),        SUNPL (NVEC),
     *  SHAPU (NVEC),        SHAPL (NVEC),
     *  SUNFU (NVEC),        SUNFL (NVEC),
     *  SU (NVEC),           SS (NVEC),           SL (NVEC),
     *  RTRAU(NVEC),         RTRAL(NVEC)
!$OMP THREADPRIVATE ( /VEP1LSX/ )
c
c----------------------------------------------------------------------
c
c     Vegetation quantities obtained from vegetation database,
c     or from stomatal/photosynthesis model that need to be retained
c     between timesteps.
c
c        IVEG   = index into veg type database
c        FU     = fraction of overall area covered by upper story
c        FL     = fraction of snow-free area covered by lower story
c        LAI    = current two-sided leaf area index
c        SAI    = current "total-bark" stem area index
c
c        RHOVEG = reflectance of an average leaf/stem
c        TAUVEG = transmittance of an average leaf/stem
c        ORIEV  = fraction of leaf/stems with vertical   orientation
c        ORIEH  = fraction of leaf/stems with horizontal orientation
c        (1-ORIEV-ORIEH = fraction with spherical (random) orientations)
c       
c        ZTOP   = height of plant top above ground
c        ZBOT   = height of lowest branches above ground
c        DLEAF  = typical leaf dimension (used for grass too)
c        DSTEM  = typical stem dimension
c        FSTOMA = 0.5 if one-sided stomata, 1 if both-sided or needles
c
c        RSMIN  = minimum stomatal resistance
c        RSMAX  = maximum stomatal resistance
c        PARCON = constant flux in PAR param, times (1-rhoveg-tauveg)
c                 (visible) to allow use of absorbed PAR, not incident
c        VPDCON = constant vapor pressure deficit in vpd param
c        TEMCON = constant multiplier in leaf temp param
c        TEMOPT = optimum temperature (constant) in leaf temp param
c        FROOT  = fraction of roots in each soil layer
c
      COMMON  /VEGLSX/
     *  IVEG (NVEC),
     *  FU (NVEC),                      FL (NVEC),
     *  LAI (NVEC,2),                   SAI (NVEC,2), 
c
     *  RHOVEG(NVEC,NBAND,2),           TAUVEG(NVEC,NBAND,2),
     *  ORIEV(NVEC,2),                  ORIEH(NVEC,2),
c
     *  ZTOP (NVEC,2),                  ZBOT (NVEC,2),
     *  DLEAF (NVEC,2),                 DSTEM (NVEC,2),
     *  FSTOMA (NVEC,2),
c
     *  RSMIN (NVEC,2),                 RSMAX (NVEC,2),
     *  PARCON (NVEC,2),                VPDCON (NVEC,2),
     *  TEMCON (NVEC,2),                TEMOPT (NVEC,2),
     *  FROOT (NVEC,NSOILAY,2)
!$OMP THREADPRIVATE ( /VEGLSX/ )
c
      REAL LAI 
c
c----------------------------------------------------------------------
c
c        In separate common to allow multiprocessing:
c
c        CHU    = specific heat of upper leaves per unit leaf area
c        CHS    = specific heat of upper stems per unit stem area
c        CHL    = specific heat of lower story per unit veg area
c
c        ALAIM[U,L] = mature 2-sided leaf+stem area index (normalizing
c                     constants for aero diffusion coeffs in TURINI)
c        BLAIM      = mature 2-sided leaf+stem area index (normalizing
c                     constant for root amount in VEGDAT)
c
c        CLEAF  = constant in upper leaf-air aero transfer param
c        CSTEM  = constant in upper stem-air aero transfer param
c        CGRASS = constant in lower story-air aero transfer param
c
c        WLIQ[U,S,L]MAX  = max mass of intercepted liquid per veg area
c        WSNO[U,S,L]MAX  = max mass of intercepted snow per veg area
c
c        TDRIPU = decay time for upper-leaf intercepted liquid dripoff
c        TDRIPS = decay time for upper-stem intercepted liquid dripoff
c        TDRIPL = decay time for lower-veg  intercepted liquid dripoff
c        TBLOWU = decay time for upper-leaf intercepted snow blowoff
c        TBLOWS = decay time for upper-stem intercepted snow blowoff
c        TBLOWL = decay time for lower-veg  intercepted snow blowoff
c
      COMMON /VEG1LSX/
     *  CHU,      CHS,      CHL,
     *  ALAIMU,   ALAIML,   BLAIM,
     *  CLEAF,    CSTEM,    CGRASS,
     *  WLIQUMAX, WLIQSMAX, WLIQLMAX, WSNOUMAX, WSNOSMAX, WSNOLMAX,
     *  TDRIPU,   TDRIPS,   TDRIPL,   TBLOWU,   TBLOWS,   TBLOWL
c
c-----------------------------------------------------------------------
c
c        For IBIS photosynthesis and stomatal conductance:
c
c        No corresponding global arrays:
c        GSSUNUD = stomatal conductance, sunlit, upper broadleaf leaves
c        GSSUNUC = stomatal conductance, sunlit, upper conifer leaves
c        GSSUNL3 = stomatal conductance, sunlit, lower c3      leaves
c        GSSUNL4 = stomatal conductance, sunlit, lower c4      leaves
c        GSSHA*  = same for shaded ...

      COMMON  /CIBIS1D/ 
     *  GSSUNUD(NVEC),                  GSSUNUC(NVEC),
     *  GSSUNL3(NVEC),                  GSSUNL4(NVEC),
     *  GSSHAUD(NVEC),                  GSSHAUC(NVEC),
     *  GSSHAL3(NVEC),                  GSSHAL4(NVEC)
!$OMP THREADPRIVATE ( /CIBIS1D/ )

c----------
#ifdef ibis
c----------

c        With corresponding global arrays in comvegxy:
c        VMAXSUNU= optimum Rubisco capacity,sunlit,upper (mol/m2/s)
c                  (...,1)=broadleaf, (...,2)=conifer
c        VMAXSUNL= optimum Rubisco capacity,sunlit,lower (mol/m2/s)
c                  (...,1)=c3, (...,2)=c4
c        VMAXSHA*= same except for shaded
c        VMAXSHA = optimum Rubisco capacity, shaded (mol/m2/s)
c        CISUNUD = intracellular co2, sunlit, upper, broadleaf (mol/mol)
c        CISUNUC = intracellular co2, sunlit, upper, conifer   (mol/mol)
c        CISUNL3 = intracellular co2, sunlit, lower, c3        (mol/mol)
c        CISUNL4 = intracellular co2, sunlit, lower, c4        (mol/mol)
c        CISHA*  = same except for shaded, ...,                (mol/mol)
c        CS*     = leaf-surface co2, ...,                      (mol/mol)
c        CONFRAC = fraction of upper story LAI that is conifer 
c        C4FRAC  = fraction of lower story LAI that is c4
c
      COMMON  /CIBIS/ 
     *  VMAXSUNU(NVEC,2),               VMAXSUNL(NVEC,2),
     *  VMAXSHAU(NVEC,2),               VMAXSHAL(NVEC,2),

     *  CISUNUD(NVEC),                  CISUNUC(NVEC),
     *  CISUNL3(NVEC),                  CISUNL4(NVEC),
     *  CISHAUD(NVEC),                  CISHAUC(NVEC),
     *  CISHAL3(NVEC),                  CISHAL4(NVEC),

     *  CSSUNUD(NVEC),                  CSSUNUC(NVEC),
     *  CSSUNL3(NVEC),                  CSSUNL4(NVEC),
     *  CSSHAUD(NVEC),                  CSSHAUC(NVEC),
     *  CSSHAL3(NVEC),                  CSSHAL4(NVEC),

     *  CONFRAC(NVEC),                  C4FRAC(NVEC)
!$OMP THREADPRIVATE ( /CIBIS/ )

c        For IBIS photosynthesis/stomatal conductance history fields:
c        ANUD = npp of C per upper broadleaf    area (mol_co2/m2/s)
c        ANUC = npp of C per upper conifer leaf area (mol_co2/m2/s)
c        ANL3 = npp of C per lower c3      leaf area (mol_co2/m2/s)
c        ANL4 = npp of C per lower c4      leaf area (mol_co2/m2/s)
c        AG*  = same except for gpp
c        GS*  = stomatal comductance of leaves (mol_h2o/m2/s)
c
      COMMON  /CIBISH/
     *  ANUD(NVEC),                     ANUC(NVEC),
     *  ANL3(NVEC),                     ANL4(NVEC),
     *  AGUD(NVEC),                     AGUC(NVEC),
     *  AGL3(NVEC),                     AGL4(NVEC),
     *  GSUD(NVEC),                     GSUC(NVEC),
     *  GSL3(NVEC),                     GSL4(NVEC)
!$OMP THREADPRIVATE ( /CIBISH/ )

c----
#else
c----

c     Dummy commons for -Xlocal declaration in jclink for sgi
      COMMON /CIBIS/ DUMIBIS
      COMMON /CIBISH/ DUMIBISH

c-----
#endif
c-----
c
