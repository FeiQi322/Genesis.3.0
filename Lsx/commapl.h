c -- $Id: commapl.h,v 1.2 2001/07/03 21:29:24 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/commapl.h,v $
c -- $Name:  $

c--------
c commapl   
c--------
c
c----------------------------------------------------------------------
c
c     Indices and weights for mapping between ccm and lsx grids.
c     Nb: 1D indices are used into 2D arrays (ie, lon + nlon*(lat-1)).
c
c        INDG(KL,1-4) = up to 4 ccm points influencing lsx point KL.
c        FACG(KL,1-4) = weights for those ccm points.
c
c        INDL(KG,1-NTILE) = up to NTILE lsx pts influencing ccm pt KG.
c        FACL(KG,1-NTILE) = weights for those lsx points.
c
c        WEIG(KL) = area of each ccm box / global area
c        WEIL(KL) = area of each lsx box / global area
c
c        ALON     = lsx longitudes (radians)
c        ALAT     = lsx latitudes (radians)
c
c        ALON2D   = lon (deg E) of rcm grid box centers (cross grid)
c        ALAT2D   = lat (deg N) of rcm grid box centers (cross grid)
c
      COMMON /MAPLSX/
     *  INDG (NLON*NLAT,4),          FACG (NLON*NLAT,4),
     *  INDL (NLONG*NLATG,NTILE),    FACL (NLONG*NLATG,NTILE),
     *  WEIG (NLONG*NLATG),          WEIL (NLON*NLAT),
c
     *  ALON (NLON),                 ALAT (NLAT)
#ifdef rcm
     * ,ALON2D(NLON,NLAT),           ALAT2D(NLON,NLAT)
#endif
c
