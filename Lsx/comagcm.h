c -- $Id: comagcm.h,v 1.2 2001/07/03 21:29:23 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comagcm.h,v $
c -- $Name:  $

c--------
c comagcm 
c--------
c
c     "agcm" quantities.
c        
c----------------------------------------------------------------------
c
c        NLONG, NLONG    = agcm longitude dimension
c        NLATG, NLATG    = agcm latitude dimension
c        NBAND, NBAND    = number of solar wavebands
c        INTERM          = unit number for agcm standard input
c        IOTERM          = unit number for agcm standard output
c        IUREST          = unit number for agcm restart file
c        DTGCM           = agcm timestep
c        OMEGDA          = 2 pi / 1 day
c
c        ALATG(JJ)       = agcm latitude grid (radians)
c        WGAUSS(JJ)      = "gaussian" weights for agcm lat intervals
c        ZENGCM(JI,JJ)   = cosine of solar zenith angle
c        SRGCMD(JI,JJ,IB)= dwnwd direct  solar flux at SIGCMA (supplied)
c                          This is per unit horizontal area, ie, it has
c                          already been multiplied by cos (zen angle).
c        SRGCMI(JI,JJ,IB)= dwnwd diffuse solar flux at SIGCMA (supplied)
c        AGCMD(JI,JJ,IB) = direct  overall surface albedo 
c        AGCMI(JI,JJ,IB) = diffuse overall surface albedo 
c        FIRABS(JI,JJ)   = absorbed IR flux at surface (supplied)
c        FIREMS(JI,JJ)   = emitted IR flux at surface (returned)
c        UGCM(JI,JJ)     = eastward air velocity at SIGCMA (supplied)
c        VGCM(JI,JJ)     = northward air velocity at SIGCMA (supplied)
c        TGCM(JI,JJ)     = air temperature (deg K) at SIGCMA (supplied)
c        QGCM(JI,JJ)     = specific humidity at SIGCMA (supplied)
c        RAICGCM(JI,JJ)  = convective rainfall rate(mass flux)(supplied)
c        SNOCGCM(JI,JJ)  = convective snowfall rate(mass flux)(supplied)
c        RAISGCM(JI,JJ)  = stratiform rainfall rate(mass flux)(supplied)
c        SNOSGCM(JI,JJ)  = stratiform snowfall rate(mass flux)(supplied)
c        PGCM(JI,JJ)     = surface pressure (supplied)
c        CLDGCM(JI,JJ)  = total cloudiness (fraction) (supplied)
c        FSEGCM(JI,JJ)   = dwnwd sensible heat flux at SIGCMA (returned)
c        FVAGCM(JI,JJ)   = downward H2O vapor flux at SIGCMA (returned)
c        TAXGCM(JI,JJ)   = eastward stress imparted to gcm (returned)
c        TAYGCM(JI,JJ)   = northward stress imparted to gcm (returned)
c
c        YTSIN           = initial ccm surface temperatures (K)
c        YORO            = land/ocean mask on ccm grid (returned)
c                          (0 = ocean, 1 = land)
c        YOROAER         = like YORO except continuous weight (0-1) and
c                          land-ice lumped with ocean
c        YFI             = gcm truncated surface geopotential
c                          (= elev*grav,m2/s2)
c        SIGCMA          = sigma (p/ps) at lowest gcm level (supplied)
c        SIGCMB          = sigma (p/ps) at top of lowest gcm layer (sup)
c        TMELT           = h2o melt point (deg K) for comsat
c
      INTEGER NLONG, NLATG, NBAND
      PARAMETER (NLONG=1, NLATG=1, NBAND=2)
c
      COMMON /DUMGCM1/
     *  DTGCM,  OMEGDA,
     *  INTERM, IOTERM, IUREST,
c
     *  ALATG(NLATG),              WGAUSS(NLATG),
     *  ZENGCM(NLONG,NLATG),
     *  SRGCMD(NLONG,NLATG,NBAND), SRGCMI(NLONG,NLATG,NBAND), 
     *  AGCMD(NLONG,NLATG,NBAND),  AGCMI(NLONG,NLATG,NBAND),
     *  FIRABS(NLONG,NLATG), 
     *  FIREMS(NLONG,NLATG),       FIREMSA(NLONG,NLATG),
     *  UGCM(NLONG,NLATG),         VGCM(NLONG,NLATG),
     *  TGCM(NLONG,NLATG),         QGCM(NLONG,NLATG),
     *  RAICGCM(NLONG,NLATG),      SNOCGCM(NLONG,NLATG),
     *  RAISGCM(NLONG,NLATG),      SNOSGCM(NLONG,NLATG),
     *  PGCM(NLONG,NLATG),         CLDGCM(NLONG,NLATG),
     *  FSEGCM(NLONG,NLATG),       FVAGCM(NLONG,NLATG),
     *  TAXGCM(NLONG,NLATG),       TAYGCM(NLONG,NLATG)
c
        COMMON /DUMGCM2/
     *  YTSIN(NLONG,NLATG),
     *  YORO(NLONG,NLATG),         YOROAER(NLONG,NLATG),
     *  YFI(NLONG,NLATG),
     *  SIGCMA,                    SIGCMB,
     *  TMELT
c
