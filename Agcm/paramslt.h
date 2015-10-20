c -- $Id: paramslt.h,v 1.2 2001/07/03 21:29:19 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/paramslt.h,v $
c -- $Name:  $

c paramslt

c-------------------------------
c Parameters and commons for slt
c-------------------------------

c-----------------------------------------------------------------------

c nlon     = gcm longitude dimension
c nlat     = gcm latitude  dimension   
c nlev     = gcm mid-layer sigma dimension   

c nlonset  = dimension of uniform longitude grid (see locset)
c nlatset  = dimension of uniform latitude  grid (see locset)
c nlevset  = dimension of uniform sigma     grid (see locset)

c polelatd = latitude (deg) poleward of which "cartesian" grid is used
c            (see vectorh) and previous dep pts are used (see slt)
c niter    = number of mid-point iterations (see slt)

c=======================
c     AGCM dimensions (sets nlon,norec,nlev...nlon,nlat,nlev used here)
#include <AGCMSIZE.h>
      parameter (nlat=norec)
c=======================

      parameter (nlonp2=nlon+2, nlatp2=nlat+2, 
     *           nlevp2=nlev+2, nlevm=nlev-1)

      parameter (nlonset=50000, nlatset=50000, nlevset=50000)

      parameter (polelatd = 80.,
     *           niter    = 2)

      parameter (pi=3.14159265358979, twopi=2.*pi, pio2=pi/2.,
     *           epsilon = 1.e-10)

c-----------------------------------------------------------------------

c Extended grids:

c lona  = long with two wraparound points on each end
c lata  = latg with two extra points on each end (poles + beyond poles)
c siga  = sigg with two extra points top and bottom (<0, 0, ..., 1, >1)
c slona = sin(lona)
c clona = cos(lona)
c slata = sin(lata)
c clata = cos(lata)

      common /sltexta/
     *  lona (-1:nlonp2),                 lata (-1:nlatp2),
     *  slona(-1:nlonp2),                 clona(-1:nlonp2),
     *  slata(-1:nlatp2),                 clata(-1:nlatp2),
     *  siga (-1:nlevp2)

      real lona, lata

c-----------------------------------------------------------------------

c Uniform (fine-spaced) grid variables:

c iset  = index of closest lona-grid  point to left of nlonset 
c         uniformly spaced values spanning lona(-1:nlonp2)
c jset  = index of closest lata-grid  point to south of nlatset
c         uniformly spaced values spanning lata(-1:nlatp2)
c kset  = index of closest siga-grid  point upwards from nlevset 
c         uniformly spaced values spanning siga(-1:nlevp2)

      common /sltunif/
     *  iset(nlonset), jset(nlatset), kset(nlevset)

c-----------------------------------------------------------------------

c Cubic-interpolation terms that depend only on the grid:

      common /sltfacc/
     * delt_x(-1:nlon+1), 
     * faca_x(0:nlon), facap_x(0:nlon), facb_x(0:nlon), facbp_x(0:nlon),

     * delt_y(-1:nlat+1), 
     * faca_y(0:nlat), facap_y(0:nlat), facb_y(0:nlat), facbp_y(0:nlat),

     * delt_z(-1:nlev+1),
     * faca_z(0:nlev), facap_z(0:nlev), facb_z(0:nlev), facbp_z(0:nlev)

c-----------------------------------------------------------------------
