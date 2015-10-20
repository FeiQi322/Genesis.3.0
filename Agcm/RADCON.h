C
C Radiation constants
C
      common/radcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,o2mmr   ,
     $              dpfo3   ,dpfco2
C
      real  gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     o2mmr,     ! O2 mass mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2     ! Voigt correction factor for CO2
c
