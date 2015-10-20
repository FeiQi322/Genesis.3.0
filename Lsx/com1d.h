c -- $Id: com1d.h,v 1.2 2001/07/03 21:29:22 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/com1d.h,v $
c -- $Name:  $

c------
c com1d  
c------
c
c        Contains temporary 1-D lsx "strip" arrays, lost when next
c        1D strip is processed.
c
c----------------------------------------------------------------------
c
c        Internal lsx fluxes, and passed to soil/ice and snow models.
c
c        SOLU   = solar flux absorbed by upper-story leaves
c        SOLS   = solar flux absorbed by upper-story stems
c        SOLL   = solar flux absorbed by lower-story leaves and stems
c        SOLG   = solar flux absorbed by snow-free soil/ice
c        SOLI   = solar flux absorbed by snow
c        SOLO   = solar flux absorbed by open ocean
c
c        PAREX  = extinction coefficient for direct PAR (lower/upper)
c        RSUNU  = stomatal resistance for upper-story sunlit leaves
c        RSHAU  = stomatal resistance for upper-story shaded leaves
c        RSUNL  = stomatal resistance for lower-story sunlit leaves
c        RSHAL  = stomatal resistance for lower-story shaded leaves
c
c        FIRB   = upward IR flux at ZA (net emitted by surface)
c        FIRU   = IR flux absorbed by upper-story leaves
c        FIRS   = IR flux absorbed by upper-story stems
c        FIRL   = IR flux absorbed by lower-story leaves and stems
c        FIRG   = IR flux absorbed by snow-free soil/ice
c        FIRI   = IR flux absorbed by snow
c        FIRO   = IR flux absorbed by open ocean
c
c        FSENA  = downward sensible heat flux at ZA
c        FSENU  = sensible heat flux from upper leaves to air
c        FSENS  = sensible heat flux from upper stems to air
c        FSENL  = sensible heat flux from lower story to air
c        FSENG  = sensible heat flux from soil/ice surface to air
c        FSENI  = sensible heat flux from snow surface to air
c        FSENO  = sensible heat flux from open ocean surface to air
c
c        FVAPA  = downward water vapor mass flux at ZA
c        FVAPUW = non-trans  H2O vapor flux from upper leaves to air
c        FVAPUT = transpired H2O vapor flux from upper leaves to air
c        FVAPS  = (non-trans.) H2O vapor flux from upper stems to air
c        FVAPLW = non-trans  H2O vapor flux from lower leaves to air
c        FVAPLT = transpired H2O vapor flux from lower story to air
c        FVAPG  = H2O vapor flux from soil/ice surface to air
c        FVAPI  = H2O vapor flux from snow surface to air
c        FVAPO  = H2O vapor flux from open ocean surface to air
c
c        TAUX   = eastward wind stress to sfc(below ZA) from atmos
c        TAUY   = poleward wind stress to sfc(below ZA) from atmos
c        TAUXO  = eastward wind stress to open ocn (below ZA) from atmos
c        TAUYO  = poleward wind stress to open ocn (below ZA) from atmos
c        TAUXI  = eastward wind stress to sea ice  (below ZA) from atmos
c        TAUYI  = poleward wind stress to sea ice  (below ZA) from atmos
c
c        RAINU  = adjusted rainfall rate at above-tree level
c        TRAINU = temperature of RAINU
c        SNOWU  = adjusted snowfall rate at above-tree level
c        TSNOWU = temperature of SNOWU
c        RAINL  = adjusted rainfall rate at below-tree level
c        TRAINL = temperature of RAINL
c        SNOWL  = adjusted snowfall rate at below-tree level
c        TSNOWL = temperature of SNOWL
c        RAING  = adjusted rainfall rate at soil/ice level
c        TRAING = temperature of RAING
c        SNOWG  = adjusted snowfall rate at soil/ice level
c        TSNOWG = temperature of SNOWG
c
      COMMON  /FLULSX/
     *  SOLU (NVEC),   SOLS (NVEC),    SOLL (NVEC),
     *  SOLG (NVEC),   SOLI (NVEC),    SOLO (NVEC),
c
     *  PAREX(NVEC,2), 
     *  RSUNU(NVEC),   RSHAU(NVEC),    RSUNL(NVEC),   RSHAL(NVEC),
c
     *  FIRB (NVEC),   FIRU (NVEC),    FIRS (NVEC),   FIRL (NVEC),
     *  FIRG (NVEC),   FIRI (NVEC),    FIRO(NVEC),
c
     *  FSENA (NVEC),  FSENU (NVEC),   FSENS (NVEC),  FSENL (NVEC),
     *  FSENG (NVEC),  FSENI(NVEC),    FSENO(NVEC),
c
     *  FVAPA (NVEC),  FVAPUW (NVEC),  FVAPUT (NVEC), FVAPS (NVEC),
     *  FVAPLW (NVEC), FVAPLT (NVEC),
     *  FVAPG (NVEC),  FVAPI (NVEC),   FVAPO(NVEC),
c
     *  TAUX (NVEC),   TAUY (NVEC),  
     *  TAUXO(NVEC),   TAUYO(NVEC),  
     *  TAUXI(NVEC),   TAUYI(NVEC),  
c
     *  RAINU(NVEC), TRAINU(NVEC),  SNOWU(NVEC),  TSNOWU(NVEC),
     *  RAINL(NVEC), TRAINL(NVEC),  SNOWL(NVEC),  TSNOWL(NVEC),
     *  RAING(NVEC), TRAING(NVEC),  SNOWG(NVEC),  TSNOWG(NVEC)
!$OMP THREADPRIVATE ( /FLULSX/ )
c
c----------------------------------------------------------------------
c        
c       For solar calculations.
c
c       ALBSOD  = direct albedo for soil/land-ice surface
c       ALBSOI  = diffuse albedo for soil/land-ice surface
c       ALBSND  = direct albedo for snow surface
c       ALBSNI  = diffuse albedo for snow surface
c       ALBICD  = direct albedo for sea-ice surface
c       ALBICI  = diffuse albedo for sea-ice surface
c       ALBOCD  = direct albedo for open ocean surface
c       ALBOCI  = diffuse albedo for open ocean surface
c
c       In the following, ** = LO for lower story, UP for upper story.
c       Also, all direct fluxes are per unit horizontal area,
c       ie, already including a factor of cos(zen angle).
c
c       AB**D/I = direct or diffuse flux absorbed  by 2-strm layer 
c       RE**D/I = direct or diffuse flux reflected by 2-strm layer 
c       F**DD   = direct  flux below 2-strm layer, for direct incident
c       F**DI   = diffuse flux below 2-strm layer, for direct incident
c       DUMMY   = (no direct flux produced for diffuse incident)
c       F**II   = diffuse flux below 2-strm layer,for diffuse incident
c       SOL2D/I = direct or diffuse downward flux below trees
c       SOL3D/I = direct or diffuse downward flux above surface
c       INDSOL  = index of current strip for pts with positive COSZEN
c       NSOL    = number of points in INDSOL
c      
      COMMON  /WSOLSX/
     *  ALBSOD(NVEC,NBAND), ALBSOI(NVEC,NBAND),
     *  ALBSND(NVEC,NBAND), ALBSNI(NVEC,NBAND),
     *  ALBICD(NVEC,NBAND), ALBICI(NVEC,NBAND),
     *  ALBOCD(NVEC,NBAND), ALBOCI(NVEC,NBAND),
c
     *  ABLOD(NVEC),  ABLOI(NVEC), 
     *  RELOD(NVEC),  RELOI(NVEC),
     *  FLODD(NVEC),  FLODI(NVEC),  
     *  DUMMY(NVEC),  FLOII(NVEC),
c
     *  ABUPD(NVEC),  ABUPI(NVEC),
     *  REUPD(NVEC),  REUPI(NVEC),
     *  FUPDD(NVEC),  FUPDI(NVEC),  
     *                FUPII(NVEC),
c
     *  SOL2D(NVEC),  SOL2I(NVEC),
     *  SOL3D(NVEC),  SOL3I(NVEC),
c
     *  INDSOL(NVEC), NSOL
!$OMP THREADPRIVATE ( /WSOLSX/ )
c
c----------------------------------------------------------------------
c
c     IR emissivities of soil,snow,landice,seaice,water (set in irrad)
c 
      COMMON /EMILSX/
     *  EMISOIL, EMISNOW, EMILICE, EMISICE, EMIOCN
c
c----------------------------------------------------------------------
c
c      Mostly for aerodynamic calculations.
c 
c      RHOA      = air density at ZA (allowing for h20 vapor)
c      CP        = specific heat of air at ZA (allowing for h20 vap)
c      UA        = wind speed at ZA
c      THICA     = thickness of lowest agcm *layer* (kg/m2)
c      Z1,Z12,Z2 = effective top, mid, and bottom levels of upper story
c      Z3,Z34,Z4 = effective top, mid, and bottom levels of lower story
c      ZA        = height above surface of lowest agcm level
c      DISP*     = zero-plane displacement height for each story
c      DI[U,L]   = 1. / momentum diffusion coefficient within each story
c      BD[U,L]   = aero. coeff ([(tau/rho)/u**2] for an infinite canopy)
c      EXPH[U,L] = exp(...) for each story's half thickness
c      EXP[U,L]  = EXPH[U,L]**2
c      ALOG*     = log(...) for various levels or soil/ice roughness
c      ALOG*     = log (roughness length of upper,lower story)
c      U*        = wind speed at varios levels
c      C*        = air transfer coeffs for upper,lower air regions
c      S[G,I]    = air-soil and air-snow transfer coeffs
c      RICHU     = Richardson Number above trees (Z1 to ZA)
c      STRAMU    = momentum correction factor for stratif above trees
c      STRAHU    = heat/vap correction factor for stratif above trees
c      RICHL     = Richardson Number below trees (Z3 to Z2)
c      STRAML    = momentum correction factor for stratif below trees
c      STRAHL    = heat/vap correction factor for stratif below trees
c      PFLUX*    = heat fluxes on leaves/stems due to intercepted h2o
c      TFAC      = (ps/p) ** (RAIR/CAIR) for lowest agcm level (const)
c
      COMMON  /WTULSX/
     *  RHOA(NVEC),   CP(NVEC),     UA(NVEC),     THICA(NVEC),
     *  Z1(NVEC),     Z12(NVEC),    Z2(NVEC), 
     *  Z3(NVEC),     Z34(NVEC),    Z4(NVEC),     ZA(NVEC),
     *  DISPU(NVEC),  DISPL(NVEC),  DIU(NVEC),    DIL(NVEC), 
     *  BDU(NVEC),    BDL(NVEC),    EXPHU(NVEC),  EXPHL(NVEC),
     *  EXPU(NVEC),   EXPL(NVEC),
     *  ALOGG(NVEC),  ALOGI(NVEC),  ALOGAV(NVEC), ALOG4(NVEC),
     *  ALOG3(NVEC),  ALOG2(NVEC),  ALOG1(NVEC),  ALOGA(NVEC),
     *  ALOGU(NVEC),  ALOGL(NVEC),
     *  U1(NVEC),     U12(NVEC),    U2(NVEC),
     *  U3(NVEC),     U34(NVEC),    U4(NVEC),
     *  CU(NVEC),     CL(NVEC),
     *  SG(NVEC),     SI(NVEC),
     *  RICHU(NVEC),  STRAMU(NVEC), STRAHU(NVEC),
     *  RICHL(NVEC),  STRAML(NVEC), STRAHL(NVEC),
     *  PFLUXU(NVEC), PFLUXS(NVEC), PFLUXL(NVEC),
     *  TFAC
!$OMP THREADPRIVATE ( /WTULSX/ )
c
c----------------------------------------------------------------------
c
c      Miscellaneous strip arrays.
c
c        FWET[U,S,L] = areal fraction of veg components wetted by
c                      intercepted liquid and snow combined.
c        RLIQ[U,S,L] = proportion of FWET[U,S,L] due to liquid.
c
      COMMON  /VARLSX/
     *  FWETU (NVEC),   FWETS (NVEC),   FWETL (NVEC),
     *  RLIQU (NVEC),   RLIQS (NVEC),   RLIQL (NVEC)
!$OMP THREADPRIVATE ( /VARLSX/ )
c
