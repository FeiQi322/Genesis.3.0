c -- $Id: comdice.h,v 1.2 2001/07/03 21:29:17 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/comdice.h,v $
c -- $Name:  $

c--------
c comdice
c--------
c
c     Parameters and commons for the dynamical seaice model.
c
c-----------------------------------------------------------------------
c
c radius of earth (m)
      parameter (earthr=6.37122e6)

c turning angle for water stress (deg)
      parameter (wturn=25.)
c     parameter (wturn=0.)

c drag coefficients for water and air ((N m-2)/(m s-1))
      parameter (csubw=0.6524, csuba=0.01256)

c flag for which cavitating iteration method to use.
c if true,  solve for free u,v just once and increment u,v (cavit)
c if false, solve for "free+pressure" u,v and increment pressure(cavit2)
      logical flato
      parameter (flato=.false.)

c max, min numbers of cavit iterations,and convergence criterion (m s-1)
c     parameter (nloopmax=100, nloopmin=30, convinc=.0030)
      parameter (nloopmax=150, nloopmin=30, convinc=.0010)

c underrelaxation factor for pressure increment in cavit
      parameter (relax=0.5)

c grids 
      common /dicegrid/ 
     *  alon(0:nlon+1),             alonu(0:nlon+1),
     *  slonu(0:nlon+1),            clonu(0:nlon+1),

     *  alat(0:nlat+1),             alatu(nlat+1),
     *  slat(0:nlat+1),             slatu(nlat+1),
     *  deltx(0:nlat+1),            deltxu(nlat+1),
     *  delty(0:nlat+1),            deltyu(nlat+1),
     *  srfarea(nlat),              srfareau(nlat+1)

c     Task common block for Cray-parallelization (in jclink for sgi)

      COMMON  /dicegrid2/
     *  latha, lathb, latua, latub
!$OMP THREADPRIVATE ( /dicegrid2/ )

c water-stess turning angles (for u-grid)
      common /diceturn/
     *  coswturn(nlat+1),           sinwturn(nlat+1)

c u-grid land-ocean mask umask
      common /dicemask/
     *  umask(0:nlon+1,nlat+1)

c diagnostic output 
      common /dicediag/
     *  loop(2), rmsinc(2), rmsvel(2), icolaps(nlon,nlat)

c work space
      common /dicework/
     *  work(nlon,nlat),            work2(nlon,nlat),
     *  workh(0:nlon+1,0:nlat+1),
     *  worku(0:nlon+1,nlat+1),     worku2(0:nlon+1,nlat+1)

