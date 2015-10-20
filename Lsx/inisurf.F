c -- $Id: inisurf.F,v 1.2 2001/07/03 21:29:28 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/inisurf.F,v $
c -- $Name:  $

c--------
c inisurf
c--------
c
c----------------------------------------------------------------------
c
      SUBROUTINE INISURF
     *(
     * MLON, MLAT, ALATG, WGAUSS, MBAND, DTGCM,
     * YTSIN,      YOROOCN,    YOROLAND, YOROICE,
#ifdef rcm
     * vegrcmin,   toprcmin,  isurftypin,
#elif hydmodel
     * lmaskhydin, tophydin,  iveghydin,   isoilhydin, 
     * thetahydin, nlayhydin, hwaterhydin, ncateghydin,
#endif
     * DELEVTYPIN, YFI,       VEGTYPIN, EVEDELTAIN,  DSBIOMIN,
     * WISCVEGIN,  PHOTOTYPIN,TAUBIOMIN,DYNAIBISIN,ASYNIBISIN,
     * SOILTEXIN,  DEPTHML,   DEPTHLAKE,DEPTHICE,DEPTHSNO,
     * OCEANTYPIN, DYNAMICE,  PRESCOUV, SHOWMAPS,
     * LONOLSX,    LATOLSX,   NSREST,   NSREST_ISOT,
     * NNEXFN_DUM, NRSTRT,    NSTOP,
     * NHISI,      NDHIS,     HISTLAST
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
     * ,alonbudin, alatbudin
#endif
#if defined (rcm)
     * ,alon2din, alat2din
#endif
     * )
c
c        Calls the various surface initialization subrs.
c        The order in which they are called matters (see below).
c
c        All arguments are supplied from the ccm except 
c        YOROOCN, YOROLAND, YOROICE which are returned.
c
c        MLON   = gcm's longitude dimension (should = lsx's value)
c        MLAT   = gcm's latitude  dimension (should = lsx's value)
c        ALATG  = gaussian latitudinal grid (radians)
c        WGAUSS = gaussian weights for each ccm latitude interval
c        MBAND  = gcm's number of solar wavebands (should = lsx's value)
c        DTGCM  = gcm's timestep
c
c        YTSIN  = gcm temperature field used to initialize lsx (deg K)
c        YOROOCN =ocean/no ocean mask averaged to ccm grid
c        YOROLAND=land/no land mask averaged to ccm grid
c        YOROICE =ice/no ice mask averaged to ccm grid

c        VEGRCMIN = vegetation map (rcm only, MM4 landuse categories)
c        TOPRCMIN = surface elevations (m) (rcm only)

c        LMASKHYDIN = land-ocean mask (hydrol only)
c        TOPHYDIN   = surface elevations (m) (hydrol only)
c        IVEGHYDIN = Dorman-Sellers veg types (hydrol only)
c        ISOILHYDIN = soil textures (100*sand% + clay%) (hydrol only)
c        THETAHYDIN = prescribed at-depth soil moistures for categs(hyd)
c        NLAYHYDIN  = uppermost soil-layer numbers for thetahyd (hydrol)
c        HWATERHYDIN = depth of standing water, categs, if any (hydrol)
c        NCATEGHYDIN = number of categories (hydrol only)

c        DELEVTYPIN= flag (0 to 4) for type of elevation correction
c        YFI       = gcm truncated sfc geopotential (elev*grav,m2/s2)
c        VEGTYPIN  = 0 to call vegdat_dorm (Dorman+Sellers),
c                    1 to call vegdat_wisc (Wisconsin)
c                    2 to call vegdat_eve  (Eve/Leaf prescribed),
c                    3 to call vegdat_eve  (Eve/Leaf interactive)
c                    4 to call vegdat_biome (Biome3/4 non-interactive),
c                    5 to call vegdat_biome (Biome3/4 interactive)
c                    6 to call vegdat_sdvm  (SDVM non-interactive)
c                    7 to call vegdat_sdvm  (SDVM interactive)
c        EVEDELTAIN = 1 for climate deltas for interactive eve (else 0)
c        DSBIOMIN  = global prescribed Dorman-Sellers single biome type
c                    (1-12; default=0), ignored unless VEGTYPIN = 0
c        WISCVEGIN = global prescribed Wisconsin single veg category,
c                    (default=0), ignored unless VEGTYPE = 1
c        PHOTOTYPIN= 0 for old stomatal resistance param (BATS/SiB)
c                  = 1 for new photosynthesis/stom cond.(Ibis)
c        TAUBIOMIN = time scale for time filter of biome forcing (years)
c        DYNAIBISIN= 1 for dynamic IBIS Wisconsin vegetation, else 0
c        ASYNIBISIN= number of dynamic-IBIS yearly iterations (dynaveg)
c        SOILTEXIN = global prescribed soil sand-silt-clay texture(or 0)
c        DEPTHML  = general ocean mixed layer depth (=0 if not entered),
c                   applied at first timestep
c        DEPTHLAKE= as DEPTHML but for inland open water
c        DEPTHICE = max limit on seaice thickness, applied at first step
c        DEPTHSNO = max limit on snow thickness, applied at first step
c        OCEANTYPIN=flag (0/1/2) for prescr.SSTs/slab/ogcm
c        DYNAMICE = .TRUE. for dynamic sea ice (used by INIOCN)
c        PRESCOUV = .TRUE. for prescribed sfc-ocn currents for dyn.ice.
c        SHOWMAPS = .TRUE. to print 2-D maps of prescribed fields
c        LONOLSX  = lon (deg E) for lsx output pt; index if icedrive/rcm
c        LATOLSX  = lat (deg N) for lsx output pt; index if icedrive/rcm
c        NSREST   = -1 for cold start, 0 for initial start from a
c                   restart file, 1 for true restart
c        NSREST_ISOT = -1 for cold-start isotopic vars,1 from rstrt file
c        NNEXFN_DUM = number of (agcm,lsx) history files per year
c                    (no longer used)
c        NRSTRT   = first timestep number of this run
c        NSTOP    = last timestep number of this run
c        NHISI    = agcm history write times (passed to datgelsx)
c        NDHIS    = agcm history accum period (passed to datgelsx)
c        HISTLAST = agcm flag to force history write last step (""")
c        ALONBUDIN= lon (deg E) of lsx output point if icedrive/rcm
c        ALATBUDIN= lat (deg N) of lsx output point if icedrive/rcm
c        ALON2DIN = lon (deg E) of rcm grid box centers (cross grid)
c        ALAT2DIN = lat (deg N) of rcm grid box centers (cross grid)
c
#include <compar.h>
#include <comeve.h>
#include <comfud.h>
#include <comgrd.h>
#include <comhisl.h>
#include <comice.h>
#ifdef gcm
#include <COMLUN.h>
#endif
#include <commapl.h>
#include <commugxy.h>
#include <comocn.h>
#include <comsoi.h>
      DIMENSION
     *  ALATG(NLATG),           WGAUSS(NLATG),
     *  YTSIN(NLONG,NLATG),
     *  YOROOCN(NLONG,NLATG),   YOROLAND(NLONG,NLATG), 
     *  YOROICE(NLONG,NLATG),
#ifdef rcm
     *  VEGRCMIN(NLONG,NLATG),  TOPRCMIN(NLONG,NLATG),
#elif hydmodel
     *  LMASKHYDIN(NLONG,NLATG),TOPHYDIN(NLONG,NLATG),
     *  IVEGHYDIN(NLONG,NLATG), ISOILHYDIN(NLONG,NLATG,NSOILAY),
     *  THETAHYDIN(*),          NLAYHYDIN(*),
     *  HWATERHYDIN(*),
#endif
     *  YFI(NLONG,NLATG),   NHISI(NMARK)

#if defined (rcm)
     * ,ALON2DIN(NLONG,NLATG),  ALAT2DIN(NLONG,NLATG)
#endif

      INTEGER DELEVTYPIN, VEGTYPIN,   EVEDELTAIN, DSBIOMIN,
     *        WISCVEGIN,  PHOTOTYPIN, DYNAIBISIN, ASYNIBISIN,
     *        SOILTEXIN,  OCEANTYPIN
      LOGICAL DYNAMICE, PRESCOUV, SHOWMAPS, HISTLAST
c
c        Copy gcm flags to lsx flags in comgrd
c
      DYNAIBIS  = DYNAIBISIN
      ASYNIBIS  = ASYNIBISIN
      VEGTYPE   = VEGTYPIN
      EVEDELTA  = EVEDELTAIN
      DSBIOME   = DSBIOMIN
      WISCVEG   = WISCVEGIN
      PHOTOTYPE = PHOTOTYPIN
      TAUBIOME  = TAUBIOMIN
      SOILTEX   = SOILTEXIN
      OCEANTYPE = OCEANTYPIN
      DELEVTYPE = DELEVTYPIN

c       If using table svp as in agcm (and not Lowe polys in comsat.h),
c       initialize table

#ifdef tablesvp
      call svpinit
#endif

c        Set lsx version numbers and machine representation format
c        (needed by datgelsx)

#ifdef gcm
      VERSLSX = '3.0'
      VERSHIS_A = 1.02
#ifdef cray
      TYPEHIS_A = 'CRAY'
#else
      TYPEHIS_A = 'IEEE'
#endif
#endif
c
c        Initialize lsx history file processing variables
c
#ifdef gcm
      CALL DATGELSX  (DTGCM, NSREST, NRSTRT, NSTOP,
     *                NHISI, NDHIS,  HISTLAST)
#else
      STOMFAC = 1.
      DELSST = 0.
      Z0OCNH = 0.
      FICEMAXN = 0.
      FICEMAXS = 0.
      TICRIT = -999.
      HICEEXP = 0.5 ! was 3. in v2.0.a.1
c v2:
      CDHFAC = 1.0
      CDEFAC = 1.35
      CDMFAC = 1.0
c v3:
c     CDHFAC = 0.6
c     CDEFAC = 0.6
c     CDMFAC = 0.6

      dslons(1) = -181.
      dslons(2) = -181.
      dslats(1) = -91.
      dslats(2) = -91.
      dsvegs    = 0
      dssois    = 0
      dslonn(1) = -181.
      dslonn(2) = -181.
      dslatn(1) = -91.
      dslatn(2) = -91.
      dsvegn    = 0
      dssoin    = 0
      imagefile = ' '
      wsoi_init = 0.6
#endif

c        If rcm, copy veg/land and topog arrays to common arrays (in 
c        commugxy), used later instead of data input files
c        SURFTYP,VEGTYPE,TOPOG2. Relies on NLONG=NLON,NLATG=NLAT.
c        Also, set unit numbers for COMLUN (basic_rcm).
c
c        If hydrol, do similarly for hydrol input data.

#ifdef rcm
      call scopy (nlong*nlatg, vegrcmin, 1, vegrcm, 1)
      call scopy (nlong*nlatg, toprcmin, 1, toprcm, 1)
      call scopy (nlong*nlatg, alon2din, 1, alon2d, 1)
      call scopy (nlong*nlatg, alat2din, 1, alat2d, 1)
      isurftype = isurftypin
      call basic_rcm
#elif hydmodel      
      do jj=1,nlatg
        do ii=1,nlong
          lmaskhyd(ii,jj)= lmaskhydin(ii,jj)
          tophyd(ii,jj)  = tophydin(ii,jj)
          iveghyd(ii,jj) = iveghydin(ii,jj)
          do kk=1,nsoilay
            isoilhyd(ii,jj,kk) = isoilhydin(ii,jj,kk)
          enddo
        enddo
      enddo
      ncateghyd = ncateghydin
      do ic =1,ncateghyd
        thetahyd(ic)  = thetahydin(ic)
        nlayhyd(ic)   = nlayhydin(ic)
        hwaterhyd(ic) = hwaterhydin(ic)
      enddo
#endif
c
c        Initialize physical consts, dimensions, lsx model
c
      CALL INILSX (MLON, MLAT, ALATG, WGAUSS, MBAND, DTGCM,
     *             YTSIN, YFI, SHOWMAPS,
     *             LONOLSX, LATOLSX,
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
     *             alonbudin, alatbudin,
#endif
     *             NSREST)
c
c        Initialize ocean mixed layer and seaice models
c
      CALL INIOCN (YTSIN, DEPTHML, DEPTHLAKE, DEPTHICE,
     *             DYNAMICE, PRESCOUV, SHOWMAPS, DTGCM, 
     *             NSREST, NSREST_ISOT)
c
c        Initialize snow model
c
      CALL INISNO (YTSIN, DEPTHSNO, DTGCM, NSREST)
c
c        Initialize soil database
c
      CALL INISOD (SHOWMAPS)
c
c       Initialize soil, seaice models (uses INISOD, INIOCN)
c
      CALL INISOI (YTSIN, NSREST, NSREST_ISOT)
c
c        Set YORO[OCN,LAND,ICE] = ocean/land/ice-only masks,
c        used by agcm for diagnostics, history, cloud, vdif, rad...
c        Note YORO's are continuous(0-1), not discrete(0 or 1)
 
      DO JJ=1,NLAT
        DO JI=1,NLON
          WORKL(JI,JJ) = CVMGT (1., 0., LMASKXY(JI,JJ).EQ.3)
        ENDDO
      ENDDO
      CALL LTOG (WORKL, YOROOCN,    0)

      DO JJ=1,NLAT
        DO JI=1,NLON
          WORKL(JI,JJ) = CVMGT (1., 0., LMASKXY(JI,JJ).EQ.1)
        ENDDO
      ENDDO
      CALL LTOG (WORKL, YOROLAND, 0)

      DO JJ=1,NLAT
        DO JI=1,NLON
          WORKL(JI,JJ) = CVMGT (1., 0., LMASKXY(JI,JJ).EQ.2)
        ENDDO
      ENDDO
      CALL LTOG (WORKL, YOROICE, 0)
c
c        Set strip indices for 1-D strips in alblsx and surfctl
c
      CALL SETIND
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE INILSX (MLON, MLAT, ALATG, WGAUSS, MBAND, DTLSX,
     *                   YTSIN, YFI, SHOWMAPS,
     *                   LONOLSX, LATOLSX,
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
     *                   alonbudin, alatbudin,
#endif
     *                   NSREST)
c
c        Does initialization for lsx model. (See INISURF for descrips
c        of the passed variables.)
c
#include <compar.h>
#include <combudl.h>
#include <comeve.h>
#include <comfud.h>
#include <comgrd.h>
#include <comicexy.h>
#ifdef gcm
#include <COMLUN.h>
#include <comhisl.h>
#endif
#include <commapl.h>
#include <commudxy.h>
#include <commugxy.h>
#include <comocnxy.h>
#ifdef isotrac
#include <comsoixy.h>
#endif
#include <comveg.h>
#include <comvegxy.h>
      DIMENSION
     *  ALATG(NLATG),         WGAUSS(NLATG),
     *  YTSIN(NLONG,NLATG),  YFI(NLONG,NLATG)
      LOGICAL SHOWMAPS
 
c     Greenland:
      PARAMETER (NGREENL=5)
      DIMENSION XGREENL(NGREENL), YGREENL(NGREENL)
      DATA XGREENL / -50., -75., -60,  -40.,  -7./
      DATA YGREENL /  85.,  78.,  70.,  50.,  85./
 
c     Antarctica:
      PARAMETER (NANTARC=4)
      DIMENSION XANTARC(NANTARC), YANTARC(NANTARC)
      DATA XANTARC /-180.,-180., 180., 180. /
      DATA YANTARC / -90., -60., -60., -90. /
 
c     Baffin Island:
      PARAMETER (NBAFFIN=5)
      DIMENSION XBAFFIN(NBAFFIN), YBAFFIN(NBAFFIN)
      DATA XBAFFIN / -84., -82., -62., -66., -70. /
      DATA YBAFFIN /  73.,  74.,  67.,  63.,  68. /
 
c     Labrador-Ungava:
      PARAMETER (NLABRAD=4)
      DIMENSION XLABRAD(NLABRAD), YLABRAD(NLABRAD)
      DATA XLABRAD / -72., -64., -64., -72. /
      DATA YLABRAD /  56.,  56.,  50.,  50. /
 
c     Keewatin #1:
      PARAMETER (NKEEWAT=4)
      DIMENSION XKEEWAT(NKEEWAT), YKEEWAT(NKEEWAT)
      DATA XKEEWAT / -94., -90., -90., -94. /
      DATA YKEEWAT /  72.,  72.,  68.,  68. /
 
c     Keewatin #2:
      PARAMETER (NKEEWA2=4)
      DIMENSION XKEEWA2(NKEEWA2), YKEEWA2(NKEEWA2)
      DATA XKEEWA2 / -92., -88., -88., -92. /
      DATA YKEEWA2 /  66.,  66.,  64.,  64. /
 
c     Melville Hills:
      PARAMETER (NMELVIL=4)
      DIMENSION XMELVIL(NMELVIL), YMELVIL(NMELVIL)
      DATA XMELVIL /-124.,-120.,-120.,-124. /
      DATA YMELVIL /  70.,  70.,  68.,  68. /
 
c     Horn Mtns:
      PARAMETER (NHORNMT=4)
      DIMENSION XHORNMT(NHORNMT), YHORNMT(NHORNMT)
      DATA XHORNMT /-120.,-118.,-118.,-120. /
      DATA YHORNMT /  64.,  64.,  62.,  62. /
c
c        Set lsx physical constants (mks, should be same values as ccm)
c
 
      STEF = 5.66961E-8
      VONK = 0.4
      GRAV = 9.80616
      TMELT = 273.16
      HVAP = 2.5104E6
      HFUS = 0.3336E6
      HSUB = HVAP + HFUS
      CH2O = 4.218E3
      CICE = 2.106E3   !cice,rhoi used for veg,soil ice (+cice for snow)
      CAIR = 1.00464E3 !Other variables are used for icesheet,sea ice
      CVAP = 1.81E3    !(set in inisoi, iniocn)
      RAIR = 287.04
      RVAP = 461.
      CAPPA = RAIR / CAIR
      RHOW = 1.E3
      RHOI = 0.9E3
      REARTH = 6.37122E+6
#ifdef isotrac
      dr_o   = 1./0.9723
      dr_h   = 1./0.9755
      smoj_o = 2.2149e-3
      smoj_h = 3.2883e-4
#endif
c
c        Set max and crit values for "cd" = drag_coeff * wind_speed,
c        for "cfl" stability of agcm lowest layer (on the order of
c        first_layer_thickness / agcm_leapfrog_timestep). Used in
c        TURCTL, TURCOF, OCNFLX to modify "cd" as follows:
c        cd = min (cdmaxa, cd/(1. + cd/cdmaxb))
c
c     CDMAXA = 100. / (2.*DTLSX)
c     CDMAXA = 300. / (2.*DTLSX)
c
c        No longer needed due to implicit treatment of TA,QA, and wind
c        stress in turlsx.f and ocean.f
c
      CDMAXA = 1.E20
      CDMAXB = 1.E20
c
c        Set lsx grid dimensions, iteration numbers, etc
c
      CALL ZERO (VZERO, NVEC)
#ifdef r4
      EPSILON = 1.E-5
#else
      EPSILON = 1.E-10
#endif
      DTTOL = 30.
c
c        Check for gcm-lsx dimension consistency
c
      IF (NLONG.NE.MLON .OR. NLATG.NE.MLAT .OR. NBAND.NE.MBAND) THEN
        WRITE(IOTERM2,900) NLONG, MLON, NLATG, MLAT, NBAND, MBAND
  900   FORMAT(/' *** Fatal Error ***'
     *         /' Dimensions declared in LSX do not agree with GCM.'
     *         /' LSX longitude  =',I8,5X,'GCM longitude  =',I8
     *         /' LSX latitude   =',I8,5X,'GCM latitude   =',I8
     *         /' LSX solar bands=',I8,5X,'GCM solar bands=',I8)
        CALL ENDRUN (-1)
      ENDIF
c
c        Set indices to be used by GTOL,LTOG for mapping between
c        ccm and lsx grids
c
#if defined (rcm)
      do j=1,nlat
c       alat(j) = alatbudin*pi/180.        ! needed in vegdat_dorm
        alat(j) = alat2d(nlon/2,j)*pi/180.
      enddo
      do i=1,nlon
        alon(i) = alon2d(i,nlat/2)*pi/180.
      enddo
#elif defined (hydmodel)
      do j=1,nlat
        alat(j) = alatg(j)                 ! needed in vegdat_dorm
      enddo
#elif defined (gcm)
      CALL SETMAP (ALATG, WGAUSS, .FALSE.)
#endif
c
c        Zero agcm-flux arrays, and global diagnostic arrays, so
c        GETIND and LTOG do not cause Cray "indefinite" errors at first
c
      CALL ZERO (ASURDXY, NLON*NLAT*NBAND)
      CALL ZERO (ASURIXY, NLON*NLAT*NBAND)
      CALL ZERO (AFIREXY, NLON*NLAT)
      CALL ZERO (ASENAXY, NLON*NLAT)
      CALL ZERO (AVAPAXY, NLON*NLAT)
      CALL ZERO (ATAUXXY, NLON*NLAT)
      CALL ZERO (ATAUYXY, NLON*NLAT)
c
      CALL ZERO (GINVAPXY, NLON*NLAT)
      CALL ZERO (GSUVAPXY, NLON*NLAT)
      CALL ZERO (GTRANSUXY, NLON*NLAT)
      CALL ZERO (GTRANSLXY, NLON*NLAT)
      CALL ZERO (GRUNOFXY, NLON*NLAT)
      CALL ZERO (GDRAINXY, NLON*NLAT)
      CALL ZERO (GDRAINHXY, NLON*NLAT)
      CALL ZERO (WFLOXY,   NLON*NLAT*(NSOILAY+1))
      CALL ZERO (GABLATXY, NLON*NLAT)
      CALL ZERO (GADJUSTXY,NLON*NLAT)
      CALL ZERO (TS2XY,    NLON*NLAT)
      CALL ZERO (Z0LENXY,  NLON*NLAT)
      CALL ZERO (GINFILXY, NLON*NLAT)
      CALL ZERO (POROSXY,  NLON*NLAT*NSOILAY)
#if defined (rcm)
      call zero (ts10xy,   nlon*nlat)
      call zero (qs10xy,   nlon*nlat)
      call zero (us10xy,   nlon*nlat)
#endif
c
c        Zero ocean fluxes to avoid Cray "indefinite" errors in first
c        GETIND
c
      CALL ZERO (FTOCNXY, NLON*NLAT)
      CALL ZERO (FWOCNXY, NLON*NLAT)
      CALL ZERO (FUOCNXY, NLON*NLAT)
      CALL ZERO (FVOCNXY, NLON*NLAT)
      CALL ZERO (FHOCNXY,   NLON*NLAT)
      CALL ZERO (FQOCNXY,   NLON*NLAT)
      CALL ZERO (FSOLOCNXY, NLON*NLAT)
      CALL ZERO (FIROCNXY,  NLON*NLAT)
      CALL ZERO (FPRLOCNXY, NLON*NLAT)
      CALL ZERO (FPRFOCNXY, NLON*NLAT)
      CALL ZERO (FSALOCNXY, NLON*NLAT)
      CALL ZERO (PRESOCNXY, NLON*NLAT)
      CALL ZERO (QFLUXXY, NLON*NLAT)
      CALL ZERO (FUICEXY, NLON*NLAT)
      CALL ZERO (FVICEXY, NLON*NLAT)
#ifdef isotrac
      call zero (fvcxy_i, nlon*nlat)
      call zero (xoxy_o, nlon*nlat)
      call zero (xoxy_h, nlon*nlat)
      call zero (qoxy_i, nlon*nlat)
      call zero (zfxy_i, nlon*nlat)
      call zero (ficexy_i, nlon*nlat)
      call zero (fvoxy_i, nlon*nlat)
      call zero (fvixy_i, nlon*nlat)
#endif
c
c       If cold start, zero time-running means of wind stress on ice
c       and open ocean for dynamic sea-ice model. (For initial or
c       restarts from early restart-file versions, these were zeroed in
c       readlsx). Also initialize vegetation prognostic variables and
c       canopy air conditions (used in TURVAP). Also set [t,p,r]eve_d
c       for first call to eve_leaf. 
c
      IF (NSREST.EQ.-1) THEN
c
        CALL ZERO (FUICEXYM, NLON*NLAT)
        CALL ZERO (FVICEXYM, NLON*NLAT)
        CALL ZERO (FUOCNXYM, NLON*NLAT)
        CALL ZERO (FVOCNXYM, NLON*NLAT)
c
        CALL GTOL (YTSIN, TUXY, 0)
        P = 98222.
        DO 50 JJ=1,NLAT
          DO 52 JI=1,NLON
            TSXY(JI,JJ) = TUXY(JI,JJ)
            TLXY(JI,JJ) = TUXY(JI,JJ)
            TLSUBXY(JI,JJ) = TUXY(JI,JJ)
            WLIQUXY(JI,JJ) = 0.
            WLIQSXY(JI,JJ) = 0.
            WLIQLXY(JI,JJ) = 0.
            WSNOUXY(JI,JJ) = 0.
            WSNOSXY(JI,JJ) = 0.
            WSNOLXY(JI,JJ) = 0.
c
            T12XY(JI,JJ) = TUXY(JI,JJ)
            Q12XY(JI,JJ) = 0.
            T34XY(JI,JJ) = TUXY(JI,JJ)
            Q34XY(JI,JJ) = 0.
            FOG12XY(JI,JJ) = 0.
            FOG34XY(JI,JJ) = 0.
            HFOG12XY(JI,JJ) = HVAP
            HFOG34XY(JI,JJ) = HVAP

            TEVE_D(JI,JJ) = 15.       ! deg C
            PEVE_D(JI,JJ) = 0.        ! mm per day
            REVE_D(JI,JJ) = 50.       ! %
c
   52     CONTINUE
   50   CONTINUE

        CALL ZERO (GIRRIGXY, NLON*NLAT)
 
      ENDIF
c
c        Read lsx surface-type data input file. The lsx mask values are
c        1=land, 2=icesheet, 3=ocean. This file has ocean pts = blank(0)
c        so convert to 3. Also older file versions may have ocean
c        pts = 4, so convert to 3.
c
c        If rcm, instead use vegrcm array (in commugxy, set in 
c        inisurf passed from the rcm). Also set fwaterxy here, 
c        not in iniocn, so can set inland lakes (slab).
c
c        If hydrol, instead use lmaskhyd array (in commugxy, set in 
c        inisurf passed from hydrol code). 
c
c---------
#ifdef rcm
c---------
      if (isurftype.eq.0) then       ! BATS landuse categories
        iocean_rcm   = 15
        icesheet_rcm = 12
        ilake_rcm    = 14
      else if (isurftype.eq.1) then  ! MM4 landuse categories
        iocean_rcm   = 7
        icesheet_rcm = 11
        ilake_rcm    = 14
      else                           ! MM5 landuse categories
        iocean_rcm   = 16
        icesheet_rcm = 24
        ilake_rcm    = 25
      endif
   
      do jj=1,nlat
        do ii=1,nlon
          if (nint(vegrcm(ii,jj)).eq.iocean_rcm) then
            lmaskxy(ii,jj) = 3
            fwaterxy(ii,jj) = 1.
          else if (nint(vegrcm(ii,jj)).eq.icesheet_rcm) then
            lmaskxy(ii,jj) = 2
            fwaterxy(ii,jj) = 0.
          else if (nint(vegrcm(ii,jj)).eq.ilake_rcm) then
            lmaskxy(ii,jj) = 1
            fwaterxy(ii,jj) = -1.         !100% inland lake (see iniocn)
          else
            lmaskxy(ii,jj) = 1
            fwaterxy(ii,jj) = 0.
          endif
        enddo
      enddo
c-------------
#elif hydmodel
c-------------
      do jj=1,nlat
        do ii=1,nlon
          lmaskxy(ii,jj) = lmaskhyd(ii,jj)
          if (lmaskxy(ii,jj).eq.0 .or. lmaskxy(ii,jj).gt.3)
     *       lmaskxy(ii,jj) = 3
        enddo
      enddo
c----
#else
c----
      CALL FINDIU (IUSURT, LOC)
#ifdef gcm
      WRITE(IOTERM2,80) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC)))
   80 FORMAT(/' Reading surface-type data file'/4X,A)
#endif
      CALL OPENF (IUSURT, LOCNAME(LOC), ' ', 0, 0)
      CALL READDATA (IUSURT, WORKL, NLON, NLAT, 'SURT', 1, 1)
      CALL CLOSEF(IUSURT)
      DO 82 JJ=1,NLAT
        DO 84 JI=1,NLON
          LMASKXY(JI,JJ) = NINT (WORKL(JI,JJ))
          IF (LMASKXY(JI,JJ).EQ.0 .OR. LMASKXY(JI,JJ).GT.3)
     *       LMASKXY(JI,JJ) = 3
   84   CONTINUE
   82 CONTINUE
c-----
#endif
c-----
c
c       Set elevation correction field delevxy (commugxy) for adjusting
c       taxy,firaxy in surfctl due to agcm truncated vs. sfc topography.
c       "True" surface-grid topography is from topog data file TOPOG2.
c       DELEVTYPE is the lsx elevation-correction-type flag (in comgrd,
c       set by passed param from gcm):
c         0 = no corrections
c         1 = corrections applied everywhere
c         2 = only Greenland + Antarctica polygons (and not ocean)
c         3 = as 2 plus all icesheet points (lmaskxy=2) (eg, 21Kyr)
c         4 = as 2 plus specific NH icesheet initiation sites
 
      CALL ZERO (DELEVXY, NLON*NLAT)
 
c---------
#ifdef gcm
c---------
 
c-----------------------------
      IF (DELEVTYPE.GE.1) THEN
c-----------------------------
 
c         Read in true topography on surface grid.
c         If rcm, instead use toprcm array (in commugxy, set in 
c         inisurf passed from the rcm).
c
c         If hydrol, instead use tophyd array (in commugxy, set in 
c         inisurf passed from hydrol code). 

#ifdef rcm
        do jj=1,nlat
          do ii=1,nlon
            workl(ii,jj) = toprcm(ii,jj)
          enddo
        enddo
#elif hydmodel
        do jj=1,nlat
          do ii=1,nlon
            workl(ii,jj) = tophyd(ii,jj)
          enddo
        enddo

#else
        CALL FINDIU (NTOPOG2, LOC)
        WRITE(IOTERM2,100) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC)))
  100   FORMAT(/' Reading surface topography data file'/4X,A)
        CALL OPENF (NTOPOG2, LOCNAME(LOC), ' ', 0, 0)
        CALL READDATA (NTOPOG2, WORKL, NLON, NLAT, 'TOPO', 2, 1)
        CALL CLOSEF(NTOPOG2)
#endif
 
c         Elevation corrections everywhere:
 
        IF (DELEVTYPE.EQ.1) THEN
 
          CALL SCOPY (NLON*NLAT, WORKL, 1, WORKM, 1)
 
 
c         Elevation corrections for specific sites:
 
        ELSE
 
          CALL RESETR (WORKM, NLON*NLAT, -999999.)
 
c         Greenland:
          CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                  XGREENL, YGREENL, NGREENL, 1.)
 
c         Antarctica:
          CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                  XANTARC, YANTARC, NANTARC, 1.)
 
c         Except ocean points:
          DO 120 JJ=1,NLAT
            DO 122 JI=1,NLON
              IF (WORKM(JI,JJ).EQ.1. .AND. LMASKXY(JI,JJ).NE.3) THEN
                WORKM(JI,JJ) = WORKL(JI,JJ)
              ELSE
                WORKM(JI,JJ) = -999999.
              ENDIF
  122       CONTINUE
  120     CONTINUE
 
c           All icesheet points (lmaskxy=2)
 
          IF (DELEVTYPE.EQ.3) THEN
            DO 130 JJ=1,NLAT
              DO 132 JI=1,NLON
                IF (LMASKXY(JI,JJ).EQ.2) WORKM(JI,JJ) = WORKL(JI,JJ)
  132         CONTINUE
  130       CONTINUE
          ENDIF
 
c           Icesheet initiation sites:
 
          IF (DELEVTYPE.EQ.4) THEN
 
c           Baffin Island:
            CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                    XBAFFIN, YBAFFIN, NBAFFIN, 900.)
 
c           Labrador-Ungava:
            CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                    XLABRAD, YLABRAD, NLABRAD, 800.)
 
c           Keewatin #1:
            CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                    XKEEWAT, YKEEWAT, NKEEWAT, 700.)
 
c           Keewatin #2:
            CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                    XKEEWA2, YKEEWA2, NKEEWA2, 700.)
 
c           Melville Hills:
            CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                    XMELVIL, YMELVIL, NMELVIL, 600.)
 
c           Horn Mts:
            CALL POLYGON (WORKM, ALON, ALAT, NLON, NLAT,
     *                    XHORNMT, YHORNMT, NHORNMT, 700.)
 
          ENDIF
 
        ENDIF
 
c         Interpolate truncated agcm topography to surface grid
 
        DO 140 JJ=1,NLATG
          DO 142 JI=1,NLONG
            WORKG(JI,JJ) = YFI(JI,JJ)/GRAV
  142     CONTINUE
  140   CONTINUE
        CALL GTOL (WORKG, WORKL, 0)
 
c         Set elevation correction (true minus agcm topography)
 
        DO 150 JJ=1,NLAT
          DO 152 JI=1,NLON
            IF (WORKM(JI,JJ).EQ.-999999.) THEN
              DELEVXY(JI,JJ) = 0.
            ELSE
              DELEVXY(JI,JJ) = WORKM(JI,JJ) - WORKL(JI,JJ)
            ENDIF
  152     CONTINUE
  150   CONTINUE

c       for history only (gcm):
        call scopy (nlong*nlatg, workl, 1, topnodxy,1)
 
c---------
      ELSE
c---------

c       for history only (gcm):
        do jj=1,nlatg
          do ji=1,nlong
            workg(ji,jj) = yfi(ji,jj)/grav
           enddo
        enddo
        call gtol (workg, topnodxy, 0)

c----------
      ENDIF
c----------
c
c        Print surface-type and elevation correction maps if requested
c
      IF (SHOWMAPS) THEN
        CALL MAPOUT (IOTERM2, LMASKXY, NLON, NLAT, ALON, ALAT, 4,
     *               'Surface Type')
c
        IF (DELEVTYPE.NE.0) THEN
          DO 510 JJ=1,NLAT
            DO 512 JI=1,NLON
              IWORKL(JI,JJ) = NINT (.01*DELEVXY(JI,JJ))
  512       CONTINUE
  510     CONTINUE
          CALL MAPOUT (IOTERM2, IWORKL, NLON, NLAT, ALON, ALAT, 0,
     *                 'Elevation correction (100s of m)')
        ENDIF
      ENDIF
c----
#else
c----

c       for history only (rcm):
        do jj=1,nlat
          do ji=1,nlon
            topnodxy(ji,jj) = yfi(ji,jj)/grav
           enddo
        enddo

c-----
#endif
c-----

c       Initialize IBIS optimum Rubisco capacity (to 200e-6 mol m-2 s-1)
c       if not obtained from restart file (see readlsx, where vmax*xy
c       are set to 1.e20 if old restart file, and writlsx where
c       1.e20's are written to restart file if phototype=0).
c       Similarly, initialiaze IBIS co2 concs (to 300e-6 mol/mol),
c       not currently on restart file.
c
c       Also initialize IBIS global history fields to avoid indef
c       errors in GETIND.
c
#ifdef ibis
      IF (NSREST.EQ.-1 .OR. VMAXSUNUXY(1,1,1).EQ.1.E20) THEN
        CALL RESETR (VMAXSUNUXY,    NLON*NLAT*2, 100.E-6)
        CALL RESETR (VMAXSUNLXY,    NLON*NLAT*2, 100.E-6)
        CALL RESETR (VMAXSHAUXY,    NLON*NLAT*2,  10.E-6)
        CALL RESETR (VMAXSHALXY,    NLON*NLAT*2,  10.E-6)
      ENDIF
      IF (NSREST.EQ.-1 .OR. UPTOTUXY(1,1).EQ.1.E20) THEN
        CALL RESETR (UPTOTUXY,      NLON*NLAT,        0.)
        CALL RESETR (UPTOTLXY,      NLON*NLAT,        0.)
        CALL GTOL (YTSIN, T15RUNXY, 0)
        CALL RESETR (AGDD0XY,       NLON*NLAT,        0.)
        CALL RESETR (AGDD5XY,       NLON*NLAT,        0.)
      ENDIF
      CALL RESETR (CISUNUDXY,  NLON*NLAT,   350.*0.7*1.E-6)
      CALL RESETR (CISUNUCXY,  NLON*NLAT,   350.*0.7*1.E-6)
      CALL RESETR (CISUNL3XY,  NLON*NLAT,   350.*0.7*1.E-6)
      CALL RESETR (CISUNL4XY,  NLON*NLAT,   350.*0.4*1.E-6)
      CALL RESETR (CISHAUDXY,  NLON*NLAT,   350.*0.7*1.E-6)
      CALL RESETR (CISHAUCXY,  NLON*NLAT,   350.*0.7*1.E-6)
      CALL RESETR (CISHAL3XY,  NLON*NLAT,   350.*0.7*1.E-6)
      CALL RESETR (CISHAL4XY,  NLON*NLAT,   350.*0.4*1.E-6)
      CALL RESETR (CSSUNUDXY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSUNUCXY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSUNL3XY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSUNL4XY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSHAUDXY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSHAUCXY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSHAL3XY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL RESETR (CSSHAL4XY,  NLON*NLAT,   350.*0.85*1.E-6)
      CALL ZERO   (CONFRACXY,  NLON*NLAT)
      CALL ZERO   (C4FRACXY,   NLON*NLAT)
 
      CALL ZERO (ANUDXY, NLON*NLAT)
      CALL ZERO (ANUCXY, NLON*NLAT)
      CALL ZERO (ANL3XY, NLON*NLAT)
      CALL ZERO (ANL4XY, NLON*NLAT)
      CALL ZERO (AGUDXY, NLON*NLAT)
      CALL ZERO (AGUCXY, NLON*NLAT)
      CALL ZERO (AGL3XY, NLON*NLAT)
      CALL ZERO (AGL4XY, NLON*NLAT)
      CALL ZERO (GSUDXY, NLON*NLAT)
      CALL ZERO (GSUCXY, NLON*NLAT)
      CALL ZERO (GSL3XY, NLON*NLAT)
      CALL ZERO (GSL4XY, NLON*NLAT)
#endif
c
c       Zero various land-only quantities over ocean (never used
c       over ocean, but for niceity in case used for history fields)
c
      DO 600 JJ=1,NLAT
        DO 602 JI=1,NLON
          IF (LMASKXY(JI,JJ).EQ.3) THEN
            T12XY(JI,JJ) = TMELT
            Q12XY(JI,JJ) = 0.
            T34XY(JI,JJ) = TMELT
            Q34XY(JI,JJ) = 0.
            FOG12XY(JI,JJ) = 0.
            FOG34XY(JI,JJ) = 0.
            HFOG12XY(JI,JJ) = HVAP
            HFOG34XY(JI,JJ) = HVAP
#ifdef ibis
            VMAXSUNUXY(JI,JJ,1) = 0.
            VMAXSUNUXY(JI,JJ,2) = 0.
            VMAXSUNLXY(JI,JJ,1) = 0.
            VMAXSUNLXY(JI,JJ,2) = 0.
            VMAXSHAUXY(JI,JJ,1) = 0.
            VMAXSHAUXY(JI,JJ,2) = 0.
            VMAXSHALXY(JI,JJ,1) = 0.
            VMAXSHALXY(JI,JJ,2) = 0.
 
            CISUNUDXY(JI,JJ) = 0.
            CISUNUCXY(JI,JJ) = 0.
            CISUNL3XY(JI,JJ) = 0.
            CISUNL4XY(JI,JJ) = 0.
            CISHAUDXY(JI,JJ) = 0.
            CISHAUCXY(JI,JJ) = 0.
            CISHAL3XY(JI,JJ) = 0.
            CISHAL4XY(JI,JJ) = 0.
 
            CSSUNUDXY(JI,JJ) = 0.
            CSSUNUCXY(JI,JJ) = 0.
            CSSUNL3XY(JI,JJ) = 0.
            CSSUNL4XY(JI,JJ) = 0.
            CSSHAUDXY(JI,JJ) = 0.
            CSSHAUCXY(JI,JJ) = 0.
            CSSHAL3XY(JI,JJ) = 0.
            CSSHAL4XY(JI,JJ) = 0.
#endif
          ENDIF
  602   CONTINUE
  600 CONTINUE
c
c       Transfer requested lon,lat of output tables to lsx variables
c       LONBUD,LATBUD, and convert to 2-D lsx indices. This is the
c       output point for BUDDO, PRINDIAG and PRINDUMP (budlsx.f)
c       and VEGDAT_EVE (vegdat.f)
 
      LONBUD = LONOLSX
      LATBUD = LATOLSX
      IF (LONBUD.EQ.-999 .AND. LATBUD.EQ.-999) THEN
c       no output, don't reset

      ELSE IF ((LONBUD.GT.1000 .AND. LATBUD.GT.1000)
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
     *         .OR. .TRUE.
#endif
     *        ) THEN
c       global lsx grid indices
        IF (LONBUD.GT.1000) LONBUD = LONBUD - 1000
        IF (LATBUD.GT.1000) LATBUD = LATBUD - 1000
        LONBUD = MAX (1, MIN (NLON, LONBUD))
        LATBUD = MAX (1, MIN (NLAT, LATBUD))
      ELSE

c       lat,lon (degrees), converted to global lsx grid indices
        IF (LONBUD.LT.0) LONBUD = LONBUD + 360
        LONBUD = LONBUD*(NLON/360.) + 1.001
        LATBUD = (LATBUD+90.)*(NLAT/180.) + 1.001
        LONBUD = MAX (1, MIN (NLON, LONBUD))
        LATBUD = MAX (1, MIN (NLAT, LATBUD))
      ENDIF
c
#if defined (icedrive) || defined (rcm) || defined (hydmodel)
      alonbud = alonbudin
      alatbud = alatbudin
#endif
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine inifud (minterm, mioterm, dtpricray, dtbudcray, icode)
c
c        Gets parameter values for this run, entered from keyboard
c        for local machine, or hard-coded for Cray (which only needs
c        xdtpri and xdtbud).
c
c        Entered values are held as x* in comfud, and assigned to the
c        appropriate variables in various initialization subrs.
c
c        icode = 0 for keyboard entry of x*
c        icode = 1 for hard-coded settings below of x* (for gcm,etc)
c        icode = 2 if x* have already been set
c        dtpricray = requested print  interval (secs) if icode=1
c        dtbudcray = requested budget interval (secs) if icode=1
c        (In all cases, convert final units of x* to mks below.)
c
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
c
c
      interm = minterm
#if defined (gcm) || defined (hydmodel)
      ioterm = 50
#else
      ioterm = mioterm
#endif
      ioterm2 = mioterm

#ifdef rcm
      ioterm = 150
      iutaba = 150
      iutabb = 151
      iutabc = 152
      iutabd = 153
      iutabe = 154
      iutabf = 155
#else
      iutaba = 50
      iutabb = 51
      iutabc = 52
      iutabd = 53
      iutabe = 54
      iutabf = 55
#endif

      if (icode.eq.0) then
 
        write(ioterm,500)
        read(interm,501)
     *    xtrun, xifres, xdtgcm, xdtpri, xdtbud, xsolgcm, xzengcm,
     *    xt1gcm, xt2gcm, xqgcm, xraigcm, xsnogcm, xwind
c
c          If blank (zero) from keyboard, assign non-zero defaults
c
c       Soil experiments 11/94:
c       if(xtrun.eq.0.) xtrun = 1.
c       if(xdtgcm.eq.0.) xdtgcm = 5.
c       if(xsolgcm.eq.0.) xsolgcm = -1.
c       if(xt1gcm.eq.0.) xt1gcm = 25.
c       if(xt2gcm.eq.0.) xt2gcm = 25.
c       if(xraigcm.eq.0.) xraigcm = 250.
 
        if(xtrun.eq.0.) xtrun = 5.
        if(xdtgcm.eq.0.) xdtgcm = 30.
        if(xdtpri.eq.0.) xdtpri = 1440. !. 180. ! 30.
        if(xdtbud.eq.0.) xdtbud = 365.
        if(xsolgcm.eq.0.) xsolgcm = 850.
        if(xzengcm.eq.0.) xzengcm = 70.
        if(xt1gcm.eq.0.) xt1gcm = 20.
        if(xt2gcm.eq.0.) xt2gcm = 30.
        if(xqgcm.eq.0.) xqgcm = 0.8
        if(xwind.eq.0.) xwind = 7.
c
c          If negative from keyboard, assign secondary defaults
c
        if(xdtpri.lt.0.) xdtpri = 0.
        if(xdtbud.lt.0.) xdtbud = 0.
        if(xsolgcm.lt.0.) xsolgcm = 0.
        if(xzengcm.lt.0.) xzengcm = 0.
        if(xqgcm.lt.0) xqgcm = 0.
        if(xraigcm.lt.0) xraigcm = 0.
        if(xsnogcm.lt.0) xsnogcm = 0.
        if(xwind.lt.0.) xwind = 0.
 
      else if (icode.eq.1) then
 
        xtrun = 0.
        xifres = 0.
        xdtgcm = 0.
        xdtpri = dtpricray/60.         ! dtpricray is in seconds
        xdtbud = dtbudcray/.864e5      ! dtbudcray is in seconds
        xsolgcm = 0.
        xzengcm = 0.
        xt1gcm = 0.
        xt2gcm = 0.
        xqgcm = 0.
        xraigcm = 0.
        xsnogcm = 0.
        xwind = 0.
 
      endif
c
#if ! defined (gcm) && ! defined (hydmodel)
      if (icode.eq.0) then
        write(ioterm,502)
     *    xtrun, nint(xifres), xdtgcm, xdtpri, xdtbud, xsolgcm, xzengcm,
     *    xt1gcm, xt2gcm, xqgcm, xraigcm, xsnogcm, xwind
      else
        write(ioterm,503)
     *    xtrun, nint(xifres), xdtgcm, xsolgcm, xzengcm,
     *    xt1gcm, xt2gcm, xqgcm, xraigcm, xsnogcm, xwind
      endif
      do iloop=1,5
        if (iloop.eq.1) iuz = iutaba
        if (iloop.eq.2) iuz = iutabb
        if (iloop.eq.3) iuz = iutabc
        if (iloop.eq.4) iuz = iutabd
        if (iloop.eq.5) iuz = iutabf
        write(iuz,503)
     *    xtrun, nint(xifres), xdtgcm, xsolgcm, xzengcm,
     *    xt1gcm, xt2gcm, xqgcm, xraigcm, xsnogcm, xwind
       enddo
#endif
c
c        Convert units to mks
c
      xtrun = xtrun*365.*.864e5
      xdtgcm = xdtgcm*60.
      xdtpri = xdtpri*60.
      xdtbud = xdtbud*.864e5
      xt1gcm = xt1gcm + 273.16
      xt2gcm = xt2gcm + 273.16
      xraigcm = xraigcm/.864e5
      xsnogcm = xsnogcm/.864e5
c
  500 format(/1x,78('-')/
     * ' trun ifres dtgcm dtpri dtbud solgcm zengcm',
     * ' t1gcm t2gcm qgcm raigcm snogcm wind')
  501 format(
     *   f5.0, f6.0, f6.0, f6.0, f6.0,  f7.0,  f7.0,
     *    f6.0, f6.0,f5.0,  f7.0,  f7.0,f5.0)
  502 format(
     *   f5.0,   i6, f6.0, f6.0, f6.0,  f7.0,  f7.1,
     *    f6.0, f6.0,f5.1,  f7.0,  f7.0,f5.0)
  503 format(
     * ' c1=''',
     * ' yrs=',f7.2,' ifres=',i2,' dtgcm=',f6.2,' solgcm=',f6.1,
     * ' zengcm=',f6.2,''''
     */' c2=''',
     * 't1gcm=',f6.2,' t2gcm=',f6.2,' qgcm=',f6.2,
     * ' raigcm=',f6.2,' snogcm=',f6.2,' wind=',f6.2,'''')
c
c
      if (icode.eq.0) then
c
        write(ioterm,600)
        read(interm,601)
     *    xdate, xmask, xwater, xvegt, xphoto, xclay, xbperm,
     *    xtsoi, xwsoi, xwisoi, xfsno, xhsno, xfice, xhice, xtocn
c
c          If blank (zero) from keyboard, assign non-zero defaults
c
        if(xdate.eq.0.)  xdate = 910101.
        if(xmask.eq.0.)  xmask= 1.
        if(xvegt.eq.0.)  xvegt = 11.  !  1.
        if(xphoto.eq.0.) xphoto = 0.
        if(xclay.eq.0.)  xclay = 3333.
        if(xbperm.eq.0.) xbperm = 0.1
        if(xtsoi.eq.0.) xtsoi = 0.5*(xt1gcm + xt2gcm) - 273.16
        if(xwsoi.eq.0. .and. xwisoi.eq.0.) then
          if (xtsoi.ge.0.) then
            xwsoi  = 0.5
          else
            xwisoi = 0.5
          endif
        endif
c
c          If negative from keyboard, assign secondary defaults
c
        if(xphoto.lt.0.) xphoto = 0.
        if(xclay.lt.0.) xclay = 0.
        if(xbperm.lt.0.) xbperm = 0.
        if(xwsoi.lt.0.) xwsoi = 0.
        if(xwisoi.lt.0.)xwisoi = 0.
c
c          Logical constraints
c
        xwater = .01 * nint(100.*xwater)
        if (xmask.le.2. .and. xwater.lt.1.) then
          if (xwater.eq.0.) then
            xfice = 0.
            xhice = 0.
            xtocn = 271.2 - 273.16
          else
            if (xtocn.eq.0..and.xfice.gt.0.) xtocn = 271.2 - 273.16
            if (xtocn.eq.0..and.xfice.eq.0.) xtocn = 0.5*(xt1gcm+xt2gcm)
     *                                               - 273.16
          endif
        else
          xmask = 3.
          xwater = 1.
          if (xtocn.eq.0..and.xfice.gt.0.) xtocn = 271.2 - 273.16
          if (xtocn.eq.0..and.xfice.eq.0.) xtocn = 0.5*(xt1gcm+xt2gcm)
     *                                               - 273.16
        endif
c
        if (xfice.eq.0.) then
           xhice = 0.
        else
           xhice = max (xhice, 0.2)
        endif
c
        if (xmask.ge.2.) then
          xvegt = 11.
        endif
 
      else if (icode.eq.1) then
 
        xdate = 0.
        xmask = 0.
        xwater = 0.
        xvegt = 0.
        xphoto = 0.
        xclay = 0.
        xbperm = 0.
        xtsoi = 0.
        xwsoi = 0.
        xwisoi= 0.
        xfsno = 0.
        xhsno = 0.
        xfice = 0.
        xhice = 0.
        xtocn = 0.
 
      endif
 
c
#if ! defined (gcm) && ! defined (hydmodel)
      if (icode.eq.0) then
        write(ioterm,602)
     *    nint(xdate),nint(xmask),xwater,nint(xvegt),nint(xphoto),
     *    nint(xclay),xbperm,
     *    xtsoi,xwsoi,xwisoi,xfsno,xhsno,xfice,xhice,xtocn
      else
        write(ioterm,603)
     *    nint(xdate),nint(xmask),xwater,nint(xvegt),nint(xphoto),
     *    nint(xclay),xbperm,
     *    xtsoi,xwsoi,xwisoi,xfsno,xhsno,xfice,xhice,xtocn
      endif
      do iloop=1,5
        if (iloop.eq.1) iuz = iutaba
        if (loopb.eq.1) iuz = iutaba
        if (loopb.eq.2) iuz = iutabb
        if (loopb.eq.3) iuz = iutabc
        if (loopb.eq.4) iuz = iutabd
        if (loopb.eq.5) iuz = iutabf 
        write(iuz,603)
     *    nint(xdate),nint(xmask),xwater,nint(xvegt),nint(xphoto),
     *    nint(xclay),xbperm,
     *    xtsoi,xwsoi,xwisoi,xfsno,xhsno,xfice,xhice,xtocn
      enddo
#endif
c
c        Convert units to mks (and fractional water to percent)
c
      xwater = xwater*100.
      xtsoi = xtsoi + 273.16
      xtocn = xtocn + 273.16
c
  600 format(/
     *'   date mask water vegt photo clay bperm',
     *  ' tsoi wsoi wisoi fsno hsno fice hice  tocn')
  601 format(
     *    f7.0,f5.0, f6.0,f5.0, f6.0,f5.0, f6.0,
     *    f5.0,f5.0, f6.0,f5.0,f5.0,f5.0,f5.0, f6.0)
  602 format(
     *      i7,  i5, f6.1,  i5,  i6,   i5, f6.2,
     *    f5.1,f5.1, f6.1,f5.1,f5.1,f5.1,f5.1, f6.2)
  603 format(
     * ' c3=''',
     * 'date=',i7,' mask=',i2,' water=',f6.2,' vegt=',i3,' photo=',i2,
     * ' clay=',i4,' bperm=',f6.2
     */' c4=''',
     * 'tsoi=',f6.2,' wsoi=',f6.2,' wisoi=',f6.2,' fsno=', f6.2,
     * ' hsno=',f6.2,' fice=',f6.2,' hice=',f6.2,' tocn=',f6.2,'''')
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine polygon (arr, alon, alat, mlon, mlat,
     *                    xvert, yvert, nvert, valin)
 
c Identifies points in a regular 2-D mlon-by-mlat grid that lie
c within a polygon defined by nvert vertices (xvert, yvert).
c The method is to test one regular grid point at a time, by calculating
c the angle swept out by the vector originating at the grid point
c and moving from one vertex to the next (anticlockwise = positive).
c Then the sum of all such arcs will be zero for points ouside the
c polygon, and +/- 360 degrees for points inside.
 
c Currently used only for setting elevation corrections for specific
c areas like Baffin Island, called from inilsx.
 
c arr = regular-grid array whose "inside" pts are set to valin(modified)
c alon = longitudes of regular grid, radians, 0 to 2*pi (supplied)
c alat = latitudes of regular grid,  radians, -pi/2 to pi/2 (supplied)
c mlon = x-dimension of regular grid (supplied)
c mlat = y-dimension of regular grid (supplied)
c xvert = vertex longitudes, degrees, -180 to 180 (supplied)
c yvert = vertex latitudes,  degrees  (supplied)
c nvert = number of vertices (supplied)
c valin = value for "inside" points of arr (supplied)
 
#include <compar.h>
 
      dimension arr (mlon,mlat), alon(mlon), alat(mlat),
     *          xvert(nvert), yvert(nvert)
 
      parameter (pi = 3.14159265358979)
c     degrees per radian
      parameter (dpr = 180./pi)
 
c     eps should be negligible but greater than machine precision
#ifdef r4
      parameter (eps = 1.e-5)
#else
      parameter (eps = 1.e-10)
#endif
 
c-----
 
c        Loop over all points in the regular grid
 
      do 100 j=1,mlat
        ytest = alat(j)
        do 102 i=1,mlon
 
c            Shift grid longitudes to range -pi to pi
c            (so polygons can span greenwich, but can't span dateline)
 
          xtest = alon(i)
          if (xtest.gt.pi) xtest = xtest - 2.*pi
 
c            Loop over vertex pairs (including nvert-to-1),
c            compute arc increment, and accumulate in angle.
 
          angle = 0.
          do 110 k=1,nvert
            kp = k+1
            if (k.eq.nvert) kp = 1
            xv1 = xvert(k) /dpr
            yv1 = yvert(k) /dpr
            xv2 = xvert(kp)/dpr
            yv2 = yvert(kp)/dpr
 
            if ( abs(yv1-ytest) + abs(xv1-xtest) .gt. eps ) then
              ang1 = dpr * atan2 (yv1-ytest, xv1-xtest)
            else
              ang1 = 0.
            endif
 
            if ( abs(yv2-ytest) + abs(xv2-xtest) .gt. eps ) then
              ang2 = dpr * atan2 (yv2-ytest, xv2-xtest)
            else
              ang2 = 0.
            endif
 
            dang = ang2 - ang1
 
c           If grid pt. on a vertex line, set to "in" and skip out
            if (abs(abs(dang)-180.).lt.eps) then
              arr(i,j) = valin
              goto 102
            endif
 
c           Increment total angle (I *think* the following logic is
c           independent of the range returned by atan2)
            if      (dang.le. 360 .and. dang .ge.180.) then
              angle = angle + dang - 360.
            else if (dang.lt. 180 .and. dang .ge. -180.) then
              angle = angle + dang
            else if (dang.lt.-180 .and. dang .ge.-360.) then
              angle = angle + dang + 360.
            endif
  110     continue
 
c         write(*,120) alat(j)*dpr, alon(i)*dpr, angle
c 120     format('alat,alon,angle=',3f10.4)
 
c            If angle is (almost) +/- 360 degrees, set to "in"
 
          if (abs(abs(angle)-360.) .lt. eps) arr(i,j) = valin
  102   continue
  100 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
#ifdef tablesvp

      subroutine svpinit

c  Table of saturation vapor pressure es from 173.16 to 376.16 K
c  in one degree increments. As in convec.F. Alternative to Lowe
c  polynomials in comsat.h. (namelist tablesvp).

      parameter (ntable=203, tmin=173.16)
      common /csatc/ table(0:ntable), dtable(0:ntable)

      data table     /.001403,.001719,.002101,.002561,.003117,.003784,
     1.004584,.005542,.006685,.008049,.009672, .01160, .01388, .01658,
     2 .01977, .02353, .02796, .03316, .03925, .04638, .05472, .06444,
     3 .07577, .08894,  .1042,  .1220,  .1425,  .1662,  .1936,  .2252,
     4  .2615,  .3032,  .3511,  .4060,  .4688,  .5406,  .6225,  .7159,
     5  .8223,  .9432,  1.080,  1.236,  1.413,  1.612,  1.838,  2.092,
     6  2.380,  2.703,  3.067,  3.476,  3.935,  4.449,  5.026,  5.671,
     7  6.393,  7.198,  8.097,  9.098,  10.21,  11.45,  12.83,  14.36,
     8  16.06,  17.94,  20.02,  22.33,  24.88,  27.69,  30.79,  34.21,
     9  37.98,  42.13,  46.69,  51.70,  57.20,  63.23,  69.85,  77.09,
     a  85.02,  93.70, 103.20, 114.66, 127.20, 140.81, 155.67, 171.69,
     1 189.03, 207.76, 227.96, 249.67, 272.98, 298.00, 324.78, 353.41,
     2 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78, 656.62,
     3 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2, 1147.4,
     4 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3, 1936.7,
     5 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1, 3167.1,
     6 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1, 5030.7,
     7 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7, 7780.2,
     8 8201.5, 8642.3, 9103.4, 9585.5, 10089., 10616., 11166., 11740.,
     9 12340., 12965., 13617., 14298., 15007., 15746., 16516., 17318.,
     a 18153., 19022., 19926., 20867., 21845., 22861., 23918., 25016.,
     1 26156., 27340., 28570., 29845., 31169., 32542., 33965., 35441.,
     2 36971., 38556., 40198., 41898., 43659., 45481., 47367., 49317.,
     3 51335., 53422., 55580., 57809., 60113., 62494., 64953., 67492.,
     4 70113., 72819., 75611., 78492., 81463., 84528., 87688., 90945.,
     5 94302., 97761.,101325.,104994.,108774.,108774./

      do i=1,ntable-1
        dtable(i) = (table(i+1)-table(i-1)) / 2.
      enddo
      dtable(0) = dtable(1)
      dtable(ntable) = dtable(ntable-1)

      return
      end
 
c-----------------------------------------------------------------------

      function esat(t)

c     returns svp (N/m2) as a function of temperature t (K)

      parameter (ntable=203, tmin=173.16)
      common /csatc/ table(0:ntable), dtable(0:ntable)
      it1 = int(t-tmin)
      zw  = (t-tmin) - it1
      it1 = max (0, min (ntable-1,it1))
      esat = (1.-zw)* table(it1) + zw* table(it1+1)
      return
      end
 
c-----------------------------------------------------------------------

      function qsat(es,p)

c     returns saturation specific humidity (kg/kg) as a function 
c     of svp es (N/m2) and ambient pressure (N/m2). Impose an upper
c     limit of 1 to avoid spurious values for very high es and/or
c     small p.

      qsat = .622 * es / max (p-(1.-.622)*es, .622*es)
      return
      end
 
c-----------------------------------------------------------------------

      function dsat(t,p)
c
c     returns d(sat. spec. humidity)/dt as a function of t (K)
c     and ambient pressure p (N/m2). Neglects q in denominator of the
c     q(vp) reln.

      parameter (ntable=203, tmin=173.16)
      common /csatc/ table(0:ntable), dtable(0:ntable)
      it1 = int(t-tmin)
      zw  = (t-tmin) - it1
      it1 = max (0, min (ntable-1,it1))
      es = (1.-zw)* table(it1) + zw* table(it1+1)
      ds = (1.-zw)*dtable(it1) + zw*dtable(it1+1)
      dsat = .622 * ds / max (p-(1.-.622)*es, .622*es)
c     write (166,*) 't=',t,'  p=',p,'  it1=',it1,' dtable=',dtable(it1)
      return
      end
 
#endif

c-----------------------------------------------------------------------
