c -- $Id: comgrd.h,v 1.2 2001/07/03 21:29:23 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comgrd.h,v $
c -- $Name:  $

c-------
c comgrd   
c-------
c
c     Pervasive surface model quantities.
c
c----------------------------------------------------------------------
c
c     Misc, including gcm-lsx and lsx-strip grid mapping
c
c        INTERM       = standard input unit number (5)
c        IOTERM       = standard output unit number (6 if local machine,
c                       or formatted file # if called from Cray)
c        IOTERM2      = standard output unit number (6 in all cases)
c        IUTAB*       = unit numbers for tabular output
c        VZERO        = a real array of zeros, of length NVEC
c        EPSILON      = small quantity to avoid zero-divides and other
c                       truncation or machine-limit troubles with small
c                       values. Should be slightly greater than O(1) 
c                       machine precision.
c        INDALL_L     = array of lon/lat indices for land 1D strips
c        INDALL_I     = array of lon/lat indices for icesheet 1D strips
c        INDALL_O     = array of lon/lat indices for ocean 1D strips
c        NPALL_L      = array of number of pts in each land 1D strip
c        NPALL_I      = array of number of pts in each icesheet 1D strip
c        NPALL_O      = array of number of pts in each ocean 1D strip
c        NSTRIP_L     = number of land 1D strips
c        NSTRIP_I     = number of icesheet 1D strips
c        NSTRIP_O     = number of ocean 1D strips
c
      COMMON /GRDLSX/
     *  INTERM,     IOTERM,     IOTERM2,
     *  IUTABA,     IUTABB,     IUTABC,    IUTABD,   IUTABE,  IUTABF,
     *  VZERO (NVEC),           EPSILON,
c
     *  INDALL_L(NVEC,NALL),   INDALL_I(NVEC,NALL), INDALL_O(NVEC,NALL),
     *  NPALL_L(NALL),         NPALL_I(NALL),       NPALL_O(NALL),  
     *  NSTRIP_L,              NSTRIP_I,            NSTRIP_O
c
c----------------------------------------------------------------------
c
c        In separate common to allow multiprocessing:
c        NSTRIP       = current number of strips (NSTRIP_[L,I,O])
c        LOCPOI       = lon and lat indices of current strip
c        NPOI         = number of points in current strip
c
      COMMON  /GRD2LSX/
     *  NSTRIP,
     *  LOCPOI (NVEC,2),        NPOI
!$OMP THREADPRIVATE ( /GRD2LSX/ )
c
c----------------------------------------------------------------------
c
c     Physical constants, et al.
c
c        STEF  = Stefan-Boltzmann constant
c        VONK  = von Karman constant
c        GRAV  = gravitational acceleration
c        TMELT = freezing point of water (except seaice bottom)
c        HFUS  = latent heat of fusion of water (except seaice bottom)
c        HVAP  = latent heat of vaporization of water 
c        HSUB  = latent heat of sublimation of ice 
c        CH2O  = specific heat of liquid water 
c        CICE  = specific heat of ice
c        CAIR  = specific heat of dry air at constant pressure
c        CVAP  = specific heat of water vapor at constant pressure
c        RAIR  = gas constant for dry air
c        RVAP  = gas constant for water vapor
c        CAPPA = RAIR/CAIR
c        RHOW  = density of liquid water (all types)
c        RHOI  = density of ice (all types)
c        REARTH= radius of earth 
c
c        CDMAXA= maximum  sfc-agcm drag coefficient * wind speed.
c        CDMAXB= critical sfc-agcm drag coefficient * wind speed,
c                based on implicit agcm layer response. These
c                ensure CFL stability of AGCM lowest layer.
c
      COMMON /PHYLSX/
     *  STEF,  VONK,  GRAV,
     *  TMELT, HFUS,  HVAP,  HSUB,  CH2O,  CICE,
     *  CAIR,  CVAP,  RAIR,  RVAP,  CAPPA,
     *  RHOW,  RHOI,  REARTH,
     *  CDMAXA,CDMAXB
#ifdef isotrac
c       ratios of isotopic species and pure water diffusion coefficients
c       in the vapor phase / VSMOW standard ratios (time 2*19/18 for 2H
c       and time 20/18 for 18O)
     * ,dr_o, dr_h, smoj_o, smoj_h,
     *  ietest, ! choice of wind-speed expt in ocnflx
     *  oceano18 !unif ocn sfc o18 (def=use ocean_i, 100=d18ocret(lat))
#endif
c
c----------------------------------------------------------------------
c
c        Temporary work space. iwork.. equivalenced to work.. only
c        to save space...on sgi, work is real*8 and iwork is int*4, 
c        so no correspondence between elements.
c
      COMMON /WRKLSX/
     *  WORKL(NLON,NLAT),  WORKM(NLON,NLAT),
     *  WORKG(NLONG,NLATG)
c
      DIMENSION WORKL1(NLON*NLAT),   IWORKL(NLON,NLAT)
      EQUIVALENCE (WORKL,WORKL1,IWORKL)
c
      DIMENSION WORKM1(NLON*NLAT),   IWORKM(NLON,NLAT)
      EQUIVALENCE (WORKM,WORKM1,IWORKM)

      DIMENSION WORKG1(NLONG*NLATG), IWORKG(NLONG,NLATG)
      EQUIVALENCE (WORKG,WORKG1,IWORKG)
c
c----------------------------------------------------------------------
c
c     Namelist parameters specific to lsx
c     Physical constants.
c
c        STOMFAC   = final multiplier for computed stomatal resistance
c        DELSST    = increment for prescribed SSTs, applied in DOPRESST
c        CD*FAC    = factors in ocean heat/vap/drag flux coeffs
c
c        DELEVTYPE = flag (0 to 4) for type of elevation correction
c        OCEANTYPE = flag (0/1/2) for prescr.SSTs/slab/ogcm
c        VEGTYPE   = 0 to call vegdat_dorm (Dorman+Sellers),
c                    1 to call vegdat_wisc (Wisconsin),
c                    2 to call vegdat_eve  (Eve/Leaf non-interactive),
c                    3 to call vegdat_eve  (Eve/Leaf     interactive)
c                    4 to call vegdat_biome(Biome non-interactive)
c                    5 to call vegdat_biome(Biome     interactive)
c                    6 to call vegdat_sdvm (SDVM non-interactive)
c                    7 to call vegdat_sdvm (SDVM     interactive)
c        EVEDELTA  = 0 for no      delta's when calling interactive EVE
c                  = 1 for climate delta's when calling interactive EVE
c        DSBIOME   = global prescribed Dorman-Sellers single biome type
c                    (1-12, default=0), ignored unless VEGTYPE = 0
c        WISCVEG   = global prescribed Wisconsin single veg category,
c                    (default=0), ignored unless VEGTYPE = 1
c        PHOTOTYPE = 0 for old stomatal resistance param (BATS/SiB)
c                  = 1 for new photosynthesis/stom cond.(Ibis)
c        TAUBIOME  = time scale for time filter of biome forcing (years)
c        DYNAIBIS  = 1 for dynamical IBIS (Wisconsin) vegetation, else 0
c        ASYNIBIS  = number of dynamic IBIS yearly iterations (dynaveg)
c        SOILTEX   = global prescribed soil sand-silt-clay texture
c                    (10000*sand + 100*silt + clay, default=0)
c
c        DSLONS(2) = longitudes (-180 to 180) between which to set veg
c                    to DSVEGS, soil to DSSOIS, in Southern Hemisphere
c        DSLATS(2) = latitudes (-90 to 0) within which to set veg
c                    to DSVEGS, soil to DSSOIS, in Southern Hemisphere
c        DSVEGS    = D+S veg type (integer) for prescribed SH rectangle
c        DSSOIS    = D+S soil sand/silt/clay code for prescr.SH rectang.
c        DSLONN,DSLATN,DSVEGN,DSSOIN = as above for Northern Hemisphere

c
      COMMON /NALLSX/
     *  STOMFAC, DELSST,
     *  CDHFAC,  CDEFAC,    CDMFAC,
c    *  WNORM,
     *  DELEVTYPE, OCEANTYPE, VEGTYPE,  EVEDELTA, DSBIOME,
     *  WISCVEG,   PHOTOTYPE, TAUBIOME, DYNAIBIS, ASYNIBIS,
     *  SOILTEX,
     *  DSLONS(2), DSLATS(2), DSVEGS, DSSOIS,
     *  DSLONN(2), DSLATN(2), DSVEGN, DSSOIN,
     *  imagefile, ifimage(nlon,nlat)

c        WNORM     = normalizing value of surface soil moisture for
c                    evaporation parameterization in turvap
c        Was in common and namelist - caused bomb in turvap with 
c        lion ifort 11.1, Feb 2010 (!)
      parameter (wnorm=0.3)
c
      INTEGER
     *  DELEVTYPE, OCEANTYPE, VEGTYPE, EVEDELTA, 
     *  DSBIOME, WISCVEG, PHOTOTYPE, DYNAIBIS, ASYNIBIS,
     *  SOILTEX,
     *  DSVEGS, DSVEGN, DSSOIS, DSSOIN

      character*200 imagefile
c
      PARAMETER (PI = 3.14159265358979)

