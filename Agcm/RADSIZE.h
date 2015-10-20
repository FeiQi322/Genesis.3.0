c     Basic radiation-code grid resolution parameters

c====================
#include <AGCMSIZE.h>
c====================

      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     plevp

      parameter (plon   = nlon,
     $           plev   = nlev,
     $           plat   = norec,
     $           plevp  = plev + 1)

