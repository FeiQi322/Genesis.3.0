c -- $Id: compar.h,v 1.2 2001/07/03 21:29:26 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/compar.h,v $
c -- $Name:  $

#ifdef rcm
c for rcm dimensions jrcm,ircm, used below to set nlong,nlatg
#include <reglsx.param>
#endif

c-------
c compar 
c-------
c
c     General save statement for all lsx subrs, to assign all local
c     variables to static storage (not stack) on Cray
#ifdef saveall
      save                                           
#endif
c
c     On sgi, declare nint and int external, and supply own versions
c     in utilsx, to avoid sgi bug with -O3 and -i8.
#ifdef sgi
      external nint, int
#endif

#if defined (NOFLUSH) || defined (OWNFLUSH)
      external flush
      integer flush
#endif
c
c----------------------------------------------------------------------
c
c     All parameters for surface models.
c
c----------------------------------------------------------------------
c        
c        NLONG    = gcm longitude dimension
c        NLATG    = gcm latitude dimension
c        NLON     = lsx longitude dimension
c        NLAT     = lsx latitude dimension
c        NTILE    = max # of lsx boxes touched by one overlayed ccm box
c        NVEC     = max # of lsx pts processed as a 1D strip
c        NALL     = max # of 1D strips
c        NBAND    = number of solar wavebands
c        NSOILAY  = number of soil layers in soil db and model (ge 2)
c        NICELAY  = number of seaice layers in seaice model (ge 2)
c        NCHEM    = number of chemical species in soil chemistry model
c        NSNOLAY  = number of snow layers in snow model (ge 2)
c        NMARK    = maximum number of lsx history write periods
c        NUSHMX   = maximum number of lsx  history fields (program list)
c
      INTEGER  NLONG,  NLATG,   NLON,    NLAT,  NTILE,   NVEC,    NALL,
     *         NBAND,  NSOILAY, NICELAY, NCHEM, NSNOLAY, 
     *         NMARK, NUSHMX
c
      PARAMETER (
c----------
#ifndef gcm
c----------
#ifdef  icedrive
c    *            nlon  =     5,       nlat   =   5,           !surfall
c    *            nlong =  nlon,       nlatg  =  nlat,
c    *            ntile =    10,
c    *            nvec  =   192,       nall   = (nlon*nlat-1)/nvec + 1,

     *            NLON  =  960,        NLAT  =   1,            !driveall
     *            NLONG = NLON,        NLATG =  NLAT,
     *            NTILE = 1,
c    *            NVEC  = 192,
     *            NVEC  =  96,         NALL  = (NLON*NLAT-1)/NVEC + 1,

#elif hydmodel

#if defined (subseta)
c    *            alongd1 =  -95.00,    alongd2= -65.00,  ! Gt Lakes
c    *            alatgd1 =   35.00,    alatgd2=  55.00,
c    *            alongd1 =  -100.00,   alongd2= -80.00,  ! Miss south
c    *            alatgd1 =   25.00,    alatgd2=  40.00,
     *            alongd1 =  -100.00,   alongd2= -80.00,  ! Miss bottom
     *            alatgd1 =   25.00,    alatgd2=  40.00,
#elif defined (subsetb)
     *            alongd1 =  -120.00,   alongd2= -70.00,  ! Mississippi
     *            alatgd1 =   25.00,    alatgd2=  55.00,
#elif defined (subsetc)
     *            alongd1 =  -95.00,    alongd2= -65.00,  ! Gt. Lakes
     *            alatgd1 =   35.00,    alatgd2=  55.00,
#elif defined (subsetd)
     *            alongd1 =  -95.50,    alongd2= -94.50,  ! 45N,95W pt
     *            alatgd1 =   44.50,    alatgd2=  45.50,
#elif defined (subseti)
     *            alongd1 =  -120.00,   alongd2= -100.00, ! Colorado
     *            alatgd1 =   30.00,    alatgd2=  47.00,
#elif defined (subsetj)
     *            alongd1 =  -125.00,   alongd2= -105.00, ! Columbia
     *            alatgd1 =   40.00,    alatgd2=  56.00,
#elif defined (subsetk)
     *            alongd1 =  -80.00,    alongd2= -74.00,  ! Susquehanna
     *            alatgd1 =   37.00,    alatgd2=  44.00,
#else
     *            alongd1 = -168.75,    alongd2= -48.75,  ! full N.Amer.
     *            alatgd1 = -4.76184,   alatgd2=  84.7532,
#endif
c T62:
     *            dlongd  = 360./192,   dlatgd = 180./94,
     *            ncat  = 6,
     *            nlongz= (alongd2-alongd1)/dlongd + 1.501,
     *            nlong = nlongz*ncat, 
     *            nlatg =  (alatgd2-alatgd1)/dlatgd + 1.501,
     *            NLON  =  NLONG,      NLAT  = NLATG,
     *            NVEC  = MIN (192, NLON*NLAT),
     *            NTILE = 1,
     *            NALL  = (NLON*NLAT-1)/NVEC + 1,

#else
c single point:
     *            NLONG =  1,          NLATG =  1,
     *            NLON  =  1,          NLAT  =  1,
     *            NTILE =  1,
     *            NVEC  =  1,          NALL  =  1,
#endif
c-----
#endif
c-----

c                 agcm resolution:

#ifdef A_T10    
     *            NLONG =  32,         NLATG = 16,             !T10   
#endif
#ifdef A_R15_12 
     *            NLONG =  48,         NLATG = 40,             !R15_12
#endif
#ifdef A_R15 
     *            NLONG =  48,         NLATG = 40,             !R15
#endif
#ifdef A_R30
     *            NLONG =  96,         NLATG = 80,             !R30
#endif
#ifdef A_T31_12
     *            NLONG =  96,         NLATG = 48,             !T31_12
#endif
#ifdef A_T31 
     *            NLONG =  96,         NLATG = 48,             !T31
#endif
#ifdef A_T31_36
     *            NLONG =  96,         NLATG = 48,             !T31_36
#endif
#ifdef A_T42_12 
     *            NLONG = 128,         NLATG = 64,             !T42_12
#endif
#ifdef A_T42 
     *            NLONG = 128,         NLATG = 64,             !T42
#endif
#ifdef A_T63_12 
     *            NLONG = 192,         NLATG = 96,             !T63_12
#endif
#ifdef A_T63 
     *            NLONG = 192,         NLATG = 96,             !T63
#endif
#ifdef A_T106_12
     *            NLONG = 320,         NLATG =160,             !T106_12
#endif
#ifdef A_T106
     *            NLONG = 320,         NLATG =160,             !T106
#endif

c                 surface resolution:

#ifdef S_T10 
     *            NLON  =  32,         NLAT  = 16,             !T10
#endif
#ifdef S_R15 
     *            NLON  =  48,         NLAT  = 40,             !R15
#endif
#ifdef S_R30
     *            NLON  =  96,         NLAT  = 80,             !R30
#endif
#ifdef S_T31 
     *            NLON  =  96,         NLAT  = 48,             !T31
#endif
#ifdef S_T42 
     *            NLON  = 128,         NLAT  = 64,             !T42
#endif
#ifdef S_T63 
     *            NLON = 192,          NLAT = 96,              !T63
#endif
#ifdef S_T106
     *            NLON = 320,          NLAT =160,              !T106
#endif
#ifdef S_2X2 
     *            NLON  = 180,         NLAT  = 90,             !2X2
#endif
#ifdef S_1XH
     *            NLON  = 360,         NLAT  =360,             !1XH
#endif

#ifdef rcm  
c    *            NLONG =  30,         NLATG = 30,     !for rcm test
c    *            NLONG =  99,         NLATG = 72,     !for rcm stage3
c    *            NLONG = 135,         NLATG = 92,     !for rcm Nat.As.
     *            NLONG = JRCM,        NLATG = IRCM,   !in reglsx.param
     *            NLON  =  NLONG,      NLAT  = NLATG,
#endif

#ifdef gcm
#ifdef rcm
c    *            NVEC    = NLON,
     *            NVEC    = MIN (192, NLON*NLAT),
     *            NTILE   = 1, 
#else
     *            NVEC    = MIN (192, NLON*NLAT), !3*64,good vector lgth
c    *            NTILE   = (NLON/NLONG+2) * (NLAT/NLATG+2),
     *            NTILE   = (NLON/NLONG+3) * (NLAT/NLATG+3),
#endif
     *            NALL    = (NLON*NLAT-1)/NVEC + 1,
#endif

     *            NBAND   = 2,
#ifdef S_T10 
     *            NSOILAY = 3,
#else
     *            NSOILAY = 6,
#endif
#if defined (TWOSTRICE)
     *            NICELAY =10,         ! 8? 6?
#else
     *            NICELAY = 3,
#endif
     *            NCHEM   = 5,         NSNOLAY = 3,
c    for standalone dynice:
c    *            NSOILAY = 1,         NICELAY = 1,
c    *            NCHEM   = 1,         NSNOLAY = 1,
c
     *            NMARK   = 20000,
     *            NUSHMX  = 200,
c
c     for Dorman-Sellers and SDVM vegetation 
     *            NVED    = 13)

#ifndef gcm
      common /cdumloc/
     *  locname(1) ! to avoid unresolved external without COMLUN
#endif
c
