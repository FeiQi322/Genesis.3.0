c -- $Id: commug.h,v 1.2 2001/07/03 21:29:25 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/commug.h,v $
c -- $Name:  $

c-------
c commug
c
c     1D arrays and scalars for communication between lsx and agcm.
c     The 1D arrays correspond to global arrays in commugxy. Fields in
c     commug[xy] correspond to those in comagcm but on the lsx grid,
c     not on the agcm grid.
c
c----------------------------------------------------------------------
c        
c       COSZEN = cosine of solar zenith angle (supplied)
c       SOLAD  = direct downwd solar flux per unit horiz area(supplied)
c       SOLAI  = diffuse downward solar flux (supplied)
c       ASURD  = direct  albedo of surface system (returned by ALBLSX)
c       ASURI  = diffuse albedo of surface system (returned by ALBLSX)
c      
c       FIRA   = absorbed IR flux by overall surface (supplied)
c       UX     = eastward air velocity (supplied)
c       UY     = northward air velocity (supplied)
c       TA     = air temperature (supplied)
c       QA     = specific humidity (supplied)
c       RAINA  = rainfall rate (mass flux) (supplied)
c       SNOWA  = snowfall rate (mass flux) (supplied)
c       PSURF  = surface pressure (supplied)
c
c       AFIRE  = emitted IR flux by overall surface (returned)
c       ASENA  = downward sensible heat flux (returned)
c       AVAPA  = downward H2O vapor flux (returned)
c       ATAUX  = eastward wind stress to sfc (returned)
c       ATAUY  = northward wind stress to sfc (returned)
c
c       Surface mask (only concerns surface models, not AGCM):
c       LMASK  = land-icesheet/ocean mask (1=land, 2=icesheet, 3=ocean)
c       FWATER = fraction of open water (ocean or lakes) in grid box 
c
      COMMON  /GCMLSX/
     *  COSZEN (NVEC),
     *  SOLAD (NVEC,NBAND),     SOLAI (NVEC,NBAND), 
     *  ASURD (NVEC,NBAND),     ASURI (NVEC,NBAND),
c 
     *  FIRA (NVEC),
     *  UX (NVEC),              UY (NVEC),  
     *  TA (NVEC),              QA (NVEC),
     *  RAINA (NVEC),           SNOWA (NVEC),
     *  PSURF (NVEC),
c
     *  AFIRE(NVEC),
     *  ASENA (NVEC),           AVAPA (NVEC),
     *  ATAUX (NVEC),           ATAUY (NVEC),
c
     *  LMASK(NVEC),            FWATER(NVEC)
!$OMP THREADPRIVATE ( /GCMLSX/ )
c
c----------------------------------------------------------------------
c
c       In separate common to allow multiprocessing:
c       SIGA   = sigma (p/ps) of lowest agcm *level* (supplied, scalar)
c       SIGB   = sigma (p/ps) of top of lowest agcm *layer*(supp,scalar)
c
      COMMON /GCM2LSX/
     *  SIGA,                   SIGB
c
