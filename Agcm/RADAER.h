c Aerosol mixing ratio 
c
      real  aermmr        ! Radiation level aerosol mass mixing ratio
      common /radaer/ aermmr(plon,plev)

!$OMP THREADPRIVATE ( /RADAER/ )
