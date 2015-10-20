c     domain-wide buffers for ir absorbtivities and emissivities
c
      real  absnxt,              ! Nearest layer absorptivities
     $     abstot,              ! Non-adjacent layer absorptivites
     $     emstot               ! Total emissivity
      common /radbuf/ absnxt(plon,plev,4,plat)
     &      , abstot(plon,plevp,plevp,plat),emstot(plon,plevp,plat)
c
