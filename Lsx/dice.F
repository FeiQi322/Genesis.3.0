c -- $Id: dice.F,v 1.2 2001/07/03 21:29:17 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/dice.F,v $
c -- $Name:  $

c
c----------------------------------------------------------------------
c
      subroutine dynice
     *  (dynamice, dynamhem, lmask, fice, hicelay, ticelay, brine,
     *   fsno, hsnolay, tsnolay,
     *   fuwindi, fvwindi, fuwindo, fvwindo, uocn, vocn,
     *   fuwindim,fvwindim,fuwindom,fvwindom,
     *   fuocn,   fvocn,
     *   ficemin, fsnomin,
     *   omega_in,dt, secyear, secrun)
 
c Main driver of the dynamical sea-ice model of Flato and Hibler (1990),
c Annals of Glaciology,14,72-77. Two cavitating-fluid iteration methods
c are provided, selected by logical flag flato in comdice. In one
c method (flato=.true.), the free-drift equations are solved just once,
c then the velocities themselves are iteratively adjusted to remove or
c reduce convergence, incrementing a "pressure" field up to a maximum
c strength. In the other method (flato=.false.), the "free-drift +
c pressure" equations are solved iteratively, again incrementing the
c pressure field where there is convergence up to a maximum strength.
 
c A staggered Arakawa-B grid is used, with thickness and concentration
c on the "h-grid", and velocities on the "u-grid" at the corners of the
c h-grid. The u-grid pt (i,j) is at the SW corner of h-grid pt (i,j).
c Land-ocean grid boxes from lsx coincide with h-grid boxes, so that
c h-grid edges are at the poles. (So u-grid centers are at the poles;
c u and v are set to zero at the poles.) A no-slip boundary condition
c is applied at land boundaries.
 
c All velocities are zeroed over land by a local "umask", so
c field values on land are never changed by this package. This is
c nice since for vectorization this package computes results both 
c over land and ocean.
 
c Generally, arrays on the h-grid run eastward from i = 1 to nlon in
c longitude with no wraparound pts (1E to 359E if 2x2 resol), and from
c j = 1 to nlat northward (89S to 89N if 2x2) in latitude. Occasionally
c these arrays are copied to a wrapped h-grid work-array workh
c dimensioned (0:nlon+1, 0:nlat+1), whose extra latitude points 0 and
c nlat+1 are referenced only in advect. u-grid arrays (including those
c supplied and returned to the ocean model) are wrapped in longitude
c but not in latitude, so are dimensioned (0:nlon+1, nlat+1), with
c lat # 1 at the south pole and lat # nlat+1 at the north pole.
 
c Ice velocities at the poles are fixed at zero, so some u-grid fields
c including wind stresses and ocean currents at j=1 and nlat+1 are
c never used. However for uniformity all u-grid arrays have the same
c dimensions.
 
c The include-file comdice contains variables unique to this model.
c The files compar and comgrd contain basic commons shared with the
c lsx surface package of the agcm. MKS units are used throughout.
 
c All arguments are supplied and on the h-grid unless otherwise noted.
 
c *** The "flato=.false." iterative method is relatively untested ***
c *** Ocean currents uocn,vocn are currently modified by polefilt ***
c *** Code is indep of grid size, except iceprint assumes 2x2 deg ***
 
c dynamice= flag (T/F). If false, no ice is advected (skip advection
c           section)...just set net ocean stresses f[u,v]ocn for ogcm.
c dynamhem= flag (T/F). If T (default), do entire dynamic ice in 2
c           loops, one for each hemisphere (and eliminate non-ice
c           latitudes first). If F, do in one global loop over all
c           lats, necessary to avoid blowup if sea ice reaches equator.
c lmask   = lsx surface-type mask: 1=land, 2=icesheet,
c           3=seaice/ocean (supplied)
c fice    = seaice areal fraction, ie, concentration (modified)
c hicelay = thickness of each seaice layer, meters (modified)
c ticelay = temperature of each seaice layer, deg K (modified)
c brine   = brine resevoir, J m-2 (modified)
c fsno    = snow areal fraction per unit sea-ice area (modified)
c hsnolay = thickness of each snow layer, meters (modified)
c tsnolay = temperature of each snow layer, deg K (modified)
c fuwindi = eastward wind stress on seaice, N/m**2 (modified)
c fvwindi = northward wind stress on seaice, N/m**2 (modified)
c fuwindo = eastward wind stress on open ocean, N/m**2
c fvwindo = northward wind stress on open ocean, N/m**2
c f*wind*m= running time means of f[u,v]wind[i,o], N/m**2
c uocn    = eastward  surface ocean current, m/s (u-grid)
c vocn    = northward surface ocean current, m/s (u-grid)
c fuocn   = net eastward  stress on ocean, N/m**2, (u-grid, returned)
c fvocn   = net northward stress on ocean, N/m**2, (u-grid, returned)
c ficemin = minimum seaice areal fraction below which to kill ice
c fsnomin = minimum snow   areal fraction below which to kill snow
c omega_in= 2*pi/day lenth for coriolis force, rad/s
c dt      = dynamic seaice time step, seconds
c secyear = seconds into current calendar year (for diagnostics only)
c secrun  = seconds into run from base date (for diagnostics only)
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      logical dynamice, dynamhem
 
      dimension
     *  lmask(nlon,nlat),
     *  fice(nlon,nlat),
     *  hicelay(nlon,nlat,nicelay),  ticelay(nlon,nlat,nicelay),
     *  brine(nlon,nlat),            fsno(nlon,nlat),
     *  hsnolay(nlon,nlat,nsnolay),  tsnolay(nlon,nlat,nsnolay),
     *  fuwindi(nlon,nlat),          fvwindi(nlon,nlat),
     *  fuwindo(nlon,nlat),          fvwindo(nlon,nlat),
     *  fuwindim(nlon,nlat),         fvwindim(nlon,nlat),
     *  fuwindom(nlon,nlat),         fvwindom(nlon,nlat),
     *  uocn(0:nlon+1,nlat+1),       vocn(0:nlon+1,nlat+1),
     *  fuocn(0:nlon+1,nlat+1),      fvocn(0:nlon+1,nlat+1)
 
c Local h-grid arrays
 
      dimension
     *  hice(nlon,nlat),
     *  ficeold(nlon,nlat),          hiceold(nlon,nlat),
     *  fsnoold(nlon,nlat),          hsnoold(nlon,nlat),
     *  pres(nlon,nlat),             streng(nlon,nlat),
 
     *  prespd(nlon,nlat),
     *  du1(nlon,nlat),              dv1(nlon,nlat),
     *  du2(nlon,nlat),              dv2(nlon,nlat),
     *  du3(nlon,nlat),              dv3(nlon,nlat),
     *  du4(nlon,nlat),              dv4(nlon,nlat)
 
c Local u-grid arrays
 
      dimension
     *  hiceu(0:nlon+1,nlat+1),      ficeu(0:nlon+1,nlat+1),
     *  u(0:nlon+1,nlat+1),          v(0:nlon+1,nlat+1),
     *  fuwindis(0:nlon+1,nlat+1),   fvwindis(0:nlon+1,nlat+1),
     *  fuwindos(0:nlon+1,nlat+1),   fvwindos(0:nlon+1,nlat+1)
 
      logical first, ifdone
      save first
      data first /.true./
 
c Statement function for ice strength (N m-1)
c     fstreng(h,f) = 0.5e4 * h * exp (-20.*(1.-f)) (Hibler,1979)
c     fstreng(h,f) = 2.7e4 * h * exp (-20.*(1.-f)) (Hibler+Walsh,1982)
c                                                  (and F+H,JPO,1992)
c     fstreng(h,f) = 1.0e4 * h * exp (max(h-5.,0.) - 10.*(1.-f))
c     fstreng(h,f) = 2.0e4 * h * exp (-10.*(1.-f)) (used before 4/92)
      fstreng(h,f) = 2.0e4 * h * exp (-20.*(1.-f))
 
      omega_old =  2.*3.14159265358979 / (86400.*365.25/366.25)
      if ( abs(omega_in - omega_old).lt.1.e-8 ) then
c       for bit agreement with old 2.0.a (was in old comdice):
        omega = (pi/43082.)
      else
        omega = omega_in
      endif

c Initial setup, done once only at the (re)start of each run
 
      if (first) call icesetup (lmask)
 
 
c Initialize some local arrays globally, since ice model will operate on
c reduced latitudinal regions of non-zero ice, but sometimes refer to
c adjacent points. (Also hice, ficeold will be needed globally for
c history field computations, and pres,streng may be needed in iceprint)
 
      call getthick (hicelay, hice)
 
      do 50 j=1,nlat
      do 50 i=1,nlon
        ficeold(i,j) = fice(i,j)
        pres(i,j) = 0.
        streng(i,j) = 0.
        icolaps(i,j) = -1
   50 continue
 
      do 52 j=1,nlat+1
      do 52 i=0,nlon+1
        u(i,j) = 0.
        v(i,j) = 0.
        hiceu(i,j) = 0.
        ficeu(i,j) = 0.
   52 continue
 
 
c Since supplied seaice wind stress is zero if no ice, fill in by
c setting to open-ocean wind stress. This avoids contamination by
c zero values at ice margins when interpolating to u-grid. Multiply
c by a factor to allow approximately for different roughness lengths,
c currently .0024 m for ocean, .01 m for seaice in lsx, with reference
c height at about 70 m (sigma=.991).
 
      do 60 j=1,nlat
      do 60 i=1,nlon
        fuwindi(i,j) = cvmgt ( fuwindo(i,j)*1.348, fuwindi(i,j),
     *                         fice(i,j).eq.0. )
        fvwindi(i,j) = cvmgt ( fvwindo(i,j)*1.348, fvwindi(i,j),
     *                         fice(i,j).eq.0. )
   60 continue
 
 
c Set latha,etc, to global extent for calls to wrapu in polefilt below
 
      latha = 1
      lathb = nlat
      latua = 1
      latub = nlat+1
 
 
c Time-smooth, transfer to u-grid and polar-filter the ice wind stresses
 
      call timefilt (fuwindi, fvwindi, fuwindim, fvwindim, dt)
      call ltou (fuwindim, fuwindis)
      call ltou (fvwindim, fvwindis)
      call polefilt (fuwindis, fvwindis, 86, 1)
 
 
c The net ocean stresses f[u,v]ocn computed for the ocean model will be
c calculated in loop 1000 for the reduced latitudinal regions. But they
c are needed everywhere, so first set them globally to the open-ocean
c wind stresses as if for no ice. (Note that f[u,v]windo from the agcm
c are supplied at all ocean points, even if 100% ice cover.) First
c time-smooth, transfer to u-grid, polar-filter the open-ocean stresses.
 
      call timefilt (fuwindo, fvwindo, fuwindom, fvwindom, dt)
      call ltou (fuwindom, fuwindos)
      call ltou (fvwindom, fvwindos)
      call polefilt (fuwindos, fvwindos, 86, 2)
 
      do 70 j=1,nlat+1
      do 70 i=0,nlon+1
        fuocn(i,j) = fuwindos(i,j)
        fvocn(i,j) = fvwindos(i,j)
   70 continue
 
 
c Polar-filter the prescribed ocean currents (only need do first call
c since these are currenlt presribed and time-invariant)
 
      if (first) call polefilt (uocn, vocn, 86, 3)
 
      first = .false.
 
 
c Overall loop for each hemisphere (1 for SH, 2 for NH)

      if (dynamhem) then
        nhem = 2      ! two hemispheric domains
      else
        nhem = 1      ! global domain, needed for ice at equator
      endif
 
#if defined ( omp )
#ifdef cray
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared   (nhem, coswturn, sinwturn,
!$OMP*           dynamice, dynamhem, work,work2,workh,worku,worku2, dt,
!$OMP*           ficemin, fsnomin, tmelt,
!$OMP*           lmask,
!$OMP*           fice,
!$OMP*           hicelay,          ticelay,
!$OMP*           brine,            fsno,
!$OMP*           hsnolay,          tsnolay,
!$OMP*           fuwindi,          fvwindi,
!$OMP*           fuwindo,          fvwindo,
!$OMP*           fuwindim,         fvwindim,
!$OMP*           fuwindom,         fvwindom,
!$OMP*           uocn,             vocn,
!$OMP*           fuocn,            fvocn,
!$OMP*           hice,
!$OMP*           ficeold,          hiceold,
!$OMP*           fsnoold,          hsnoold,
!$OMP*           pres,             streng,
!$OMP*           prespd,           omega,
!$OMP*           du1,              dv1,
!$OMP*           du2,              dv2,
!$OMP*           du3,              dv3,
!$OMP*           du4,              dv4,
!$OMP*           hiceu,            ficeu,
!$OMP*           u,                v,
!$OMP*           fuwindis,         fvwindis,
!$OMP*           fuwindos,         fvwindos)
!$OMP* private  (i, j, k, iloop, ifdone, ihem,
!$OMP*           latha, lathb, latua, latub)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared   (nhem, coswturn, sinwturn,
!$OMP*           dynamice, dynamhem, work,work2,workh,worku,worku2, dt,
!$OMP*           ficemin, fsnomin, tmelt,
!$OMP*           lmask,
!$OMP*           fice,
!$OMP*           hicelay,          ticelay,
!$OMP*           brine,            fsno,
!$OMP*           hsnolay,          tsnolay,
!$OMP*           fuwindi,          fvwindi,
!$OMP*           fuwindo,          fvwindo,
!$OMP*           fuwindim,         fvwindim,
!$OMP*           fuwindom,         fvwindom,
!$OMP*           uocn,             vocn,
!$OMP*           fuocn,            fvocn,
!$OMP*           hice,
!$OMP*           ficeold,          hiceold,
!$OMP*           fsnoold,          hsnoold,
!$OMP*           pres,             streng,
!$OMP*           prespd,           omega,
!$OMP*           du1,              dv1,
!$OMP*           du2,              dv2,
!$OMP*           du3,              dv3,
!$OMP*           du4,              dv4,
!$OMP*           hiceu,            ficeu,
!$OMP*           u,                v,
!$OMP*           fuwindis,         fvwindis,
!$OMP*           fuwindos,         fvwindos)
!$OMP* private  (i, j, k, iloop, ifdone, ihem)
!$OMP* schedule(dynamic,2)
#endif
#endif
 
c==========================
      do 1000 ihem = 1,nhem
c==========================
 
 
c Find the reduced latitudinal band for this hemisphere, bounded by
c latha,lathb (h-grid) and latua,latub (u-grid). Skip this hemisphere
c if no ice, flagged by latha=0. (lat[ha,hb,ua,ub] are in common
c dicegrid, set as task common for multiprocessing)
 
      call icebound (fice, ihem, dynamhem)
      if (latha.eq.0) goto 1000
 
 
c Interpolate hice and fice to u-grid (hiceu, ficeu)
 
      call ltou (hice, hiceu)
      call ltou (fice, ficeu)
 
 
c Dynamical section for u and v.
c First, initialize pressure to zero and set strength.
 
      do 200 j=latha,lathb
      do 200 i=1,nlon
        pres(i,j) = 0.
        streng(i,j) = fstreng ( hice(i,j), fice(i,j) )
  200 continue
 
 
c If "flato" method, solve for free-drift velocities (and set the ring
c at 90N to pure slab motion). If "not flato" method, initialize u,v to
c zero for the first iteration.
 
      if (flato) then
        call icefree (u,v,hiceu,fice,fuwindis,fvwindis,uocn,vocn,pres,
     *                omega)
        call polefilt (u, v, 90, 4)
      else
        do 210 j=latua,latub
        do 210 i=0,nlon+1
          u(i,j) = 0.
          v(i,j) = 0.
  210   continue
      endif
 
c Set prespd = 1/divergence and du*,dv*= corner velocs, per unit "delta"
c pressure at each h-grid pt with pressure=0 elsewhere, for use in cavit
 
      call cavset (hiceu, prespd, du1,dv1,du2,dv2,du3,dv3,du4,dv4,omega)
 
c Begin cavitating-fluid iteration
 
      do 250 iloop=1,nloopmax
 
c          Save pre-iteration u,v for convergence test in icediag.
c          (So cavit and icefree must not change worku[2]).
 
        do 252 j=latua,latub
        do 252 i=0,nlon+1
          worku(i,j) = u(i,j)
          worku2(i,j) = v(i,j)
  252   continue
 
c          If "flato" method, apply increments to both pressure and u,v
c          in cavit. If "not flato" method, apply increment to pressure
c          then call icefree to solve for the new free-drift u,v.
 
        call cavit (u, v, pres, streng,
     *              prespd,du1,dv1,du2,dv2,du3,dv3,du4,dv4)
        if (.not.flato)
     *    call icefree (u,v,hiceu,fice,fuwindis,fvwindis,uocn,vocn,pres,
     *                  omega)
 
c          Reset the ring of velocities at 90N to pure slab motion
 
        call polefilt (u, v, 90, 4)
 
c          Test if converged - if so, compute diagnostics and skip out
 
        call icediag (u,v, fice,ficeu, pres,streng, iloop, ifdone, ihem)
        if (ifdone) goto 260
 
  250 continue
  260 continue
 
 
c Compute the net stress on the ocean (f[u,v]ocn, on u-grid) for use
c by the dynamic ocean model, by weighting sub-ice stresses computed
c here with open-ocean stresses f[u,v]windos. Note that total stress
c between agcm and ocean is not conserved due to earlier transfer from
c h-grid f[u,v]windo to u-grid f[u,v]windos. (Could fix by weighting
c with 1-fice before transfer, but that would complicate spatial polar
c filtering and time smoothing.)
 
      do 300 j=latua,latub
      do 300 i=1,nlon
        fuocn(i,j) = (1.-ficeu(i,j))*fuwindos(i,j)
     *             + ficeu(i,j)*csubw*( (u(i,j)-uocn(i,j))*coswturn(j)
     *                                 -(v(i,j)-vocn(i,j))*sinwturn(j))
        fvocn(i,j) = (1.-ficeu(i,j))*fvwindos(i,j)
     *             + ficeu(i,j)*csubw*( (v(i,j)-vocn(i,j))*coswturn(j)
     *                                 +(u(i,j)-uocn(i,j))*sinwturn(j))
  300 continue
      call wrapu (fuocn)
      call wrapu (fvocn)
 
 
c Advection section.
c Fields are advected in a particular order so already-advected fields
c can be used to recover thickness from mass and temperature from heat.
c Each layer of the multi-layer thermodynamic ice and snow models is
c advected totally independently of the other layers (!). Any illogical
c results (eg, snow with no ice) should be taken care of by logical
c adjustments in the thermodynamic model(s), but I don't think such
c cases can arise here except possibly in pathological situations(?).
 
c-----------------------
      if (dynamice) then
c-----------------------
 
c        advect ice concentration
 
      do 400 j=latha,lathb
      do 400 i=1,nlon
        ficeold(i,j) = fice(i,j)
        work(i,j) = fice(i,j)
        work2(i,j) = 1.
  400 continue
      call advect (fice, work, work2, u, v, dt, ficemin, 0)
 
      do 420 k=1,nicelay
 
c          advect ice layer mass, recovering thickness
 
        do 422 j=latha,lathb
        do 422 i=1,nlon
          hiceold(i,j) = hicelay(i,j,k)
          work(i,j) = ficeold(i,j) * hicelay(i,j,k)
          work2(i,j) = fice(i,j)
  422   continue
        call advect (hicelay(1,1,k), work, work2, u, v, dt, 0., 1)
 
c          advect ice layer heat, recovering temperature
 
        do 424 j=latha,lathb
        do 424 i=1,nlon
          work(i,j)  = ficeold(i,j) * hiceold(i,j) * ticelay(i,j,k)
          work2(i,j) = fice(i,j) * hicelay(i,j,k)
  424   continue
        call advect (ticelay(1,1,k), work, work2, u, v, dt, 0., 2)
 
  420 continue
 
c        advect ice brine resevoir
 
      do 440 j=latha,lathb
      do 440 i=1,nlon
        work(i,j) = ficeold(i,j)*brine(i,j)
        work2(i,j) = fice(i,j)
  440 continue
      call advect (brine, work, work2, u, v, dt, 0., 3)
 
 
c        advect snow concentration
 
      do 500 j=latha,lathb
      do 500 i=1,nlon
        fsnoold(i,j) = fsno(i,j)
        work(i,j) = ficeold(i,j) * fsno(i,j)
        work2(i,j) = fice(i,j)
  500 continue
      call advect (fsno, work, work2, u, v, dt, fsnomin, 10)
 
      do 520 k=1,nsnolay
 
c          advect snow layer mass, recovering thickness
 
        do 522 j=latha,lathb
        do 522 i=1,nlon
          hsnoold(i,j) = hsnolay(i,j,k)
          work(i,j) = ficeold(i,j) * fsnoold(i,j) * hsnolay(i,j,k)
          work2(i,j) = fice(i,j) * fsno(i,j)
  522   continue
        call advect (hsnolay(1,1,k), work, work2, u, v, dt, 0., 11)
 
c          advect snow layer heat, recovering temperature
 
        do 524 j=latha,lathb
        do 524 i=1,nlon
          work(i,j)  = ficeold(i,j) * fsnoold(i,j) * hsnoold(i,j)
     *                 * tsnolay(i,j,k)
          work2(i,j) = fice(i,j) * fsno(i,j) * hsnolay(i,j,k)
  524   continue
        call advect (tsnolay(1,1,k), work, work2, u, v, dt, 0., 12)
 
  520 continue
 
 
c Kill off small amounts of ice or snow. ficemin and fsnomin (supplied
c by lsx) should be larger than what can typically be advected into a
c previously ice/snow-free grid box in one time step, otherwise the ice
c margin can never advance.
 
      do 600 j = latha,lathb
      do 600 i = 1,nlon
        if (lmask(i,j).ge.3) then
 
          if (fice(i,j).lt.ficemin) then
            fice(i,j) = 0.
            brine(i,j) = 0.
            do 602 k=1,nicelay
              hicelay(i,j,k) = 0.
              ticelay(i,j,k) = tmelt
  602       continue
          endif
 
c         no snow (either because fsno is small, or because no seaice)
 
          if (fsno(i,j).lt.fsnomin .or. fice(i,j).lt.ficemin) then
            fsno(i,j) = 0.
            do 606 k=1,nsnolay
              hsnolay(i,j,k) = 0.
              tsnolay(i,j,k) = tmelt
  606       continue
          endif
 
        endif
  600 continue
 
c----------
      endif
c----------
 
c Loop back for other hemisphere
 
c=============
 1000 continue
c=============
 
c lsx history fields section.
c Since history variables must be on h-grid, transfer u-grid variables
c (u,v) to h-grid first, (so eventual plots will be erroneously shifted
c in longitude and latitude). Shift in latitude to get u,v at north pole
c and lose the south pole.
 
      do 2000 j=1,nlat
        do 2002 i=1,nlon
           if (ficeu(i,j+1)*hiceu(i,j+1).gt.0.) then
             work (i,j) = u(i,j+1)
             work2(i,j) = v(i,j+1)
           else
             work (i,j) = 0.
             work2(i,j) = 0.
           endif
 2002   continue
 2000 continue
      call puthilsx (work,  38)
      call puthilsx (work2, 39)

c fice = new fraction, ficeold = old. hicelay = new layer thicknesses,
c hice = old total thickness, updated here to new. 
 
      do 2010 j=1,nlat
        do 2012 i=1,nlon
           work2(i,j) = hice(i,j)
 2012   continue
 2010 continue
 
      call getthick (hicelay, hice)
 
      do 2014 j=1,nlat
        do 2016 i=1,nlon
           work (i,j) = (fice(i,j)-ficeold(i,j)) / dt
           work2(i,j) = (hice(i,j)-work2(i,j)) / dt
 2016   continue
 2014 continue
      call puthilsx (work,  44)
      call puthilsx (work2, 45)
 
 
c Diagnostic grid printout (assumes 2x2 deg)
 
c     if (nlon.eq.180 .and. nlat.eq.90) 
c    *  call iceprint (hicelay, hice, fice, fuwindis, fvwindis,
c    *                 uocn, vocn, u, v, pres, streng, lmask,
c    *                 dt, secyear, secrun)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine icesetup (lmask)
 
c Initialization for dynamical sea ice, called once at (re)start.
 
c lmask = lsx surface-type mask, on h-grid (supplied)
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension lmask(nlon,nlat)
 
 
c twopi is a local variable (not in common), and is scaled by nlon/180
c so that the longitudinal box width is always 2 degrees independent of
c nlon (so we can run fast test-cases with just a few boxes, not
c extending e-w round the world).
 
      twopi = 2.*pi * (nlon/180.)
 
 
c Set grid metrics. deltx,delty are e-w and n-s widths of h-grid boxes,
c and deltxu,deltyu are e-w and n-s widths of u-grid boxes.
 
      do 10 i=0,nlon+1
        alon(i) =  (i-.5)*twopi/nlon
        alonu(i) = (i-1.)*twopi/nlon
        slonu(i) = sin(alonu(i))
        clonu(i) = cos(alonu(i))
   10 continue
 
      do 12 j=1,nlat
        alat(j)  = ((j-.5)/nlat)*pi - pi/2.
        slat(j) = sin(alat(j))
        deltx(j)  = earthr * cos(alat(j)) * (twopi/nlon)
        delty(j)  = earthr * (pi/nlat)
   12 continue
      alat(0) = 2*alat(1)-alat(2)
      slat(0) = 2*slat(1)-slat(2)
      deltx(0) = deltx(1)
      delty(0) = delty(1)
      alat(nlat+1) = 2*alat(nlat)-alat(nlat-1)
      slat(nlat+1) = 2*slat(nlat)-slat(nlat-1)
      deltx(nlat+1) = deltx(nlat)
      delty(nlat+1) = delty(nlat)
 
      do 14 j=1,nlat+1
        alatu(j) = ((j-1.) /nlat)*pi - pi/2.
        slatu(j) = sin(alatu(j))
        deltxu(j) = earthr * cos(alatu(j)) * (twopi/nlon)
        deltyu(j) = delty(j)
   14 continue
 
 
c Set the surface area of each h-grid box
 
      do 100 j=1,nlat
        srfarea(j) = (slatu(j+1)-slatu(j))*(twopi/nlon)*earthr*earthr
  100 continue
 
 
c Set the surface area of each u-grid box, including triangular pole
c areas
 
      do 120 j=2,nlat
        srfareau(j) = (slat(j)-slat(j-1))*(twopi/nlon) * earthr*earthr
  120 continue
 
      srfareau(1) = (slat(1)-slatu(1))*(twopi/nlon) * earthr*earthr
      srfareau(nlat+1) = (slatu(nlat+1)-slat(nlat))*(twopi/nlon)
     *                   * earthr*earthr
 
 
c Set the water turning parameters on u-grid
 
      c = cos(wturn*pi/180.)
      s = sin(wturn*pi/180.)
      do 200 j=1,nlat+1
        coswturn(j) = c
        sinwturn(j) = s
        if (alatu(j).lt.0.) sinwturn(j) = -s
  200 continue
 
 
c Set u-grid land-ocean mask: 1 for ocean points, 0 for any point
c on the corner of a land h-grid box.
 
      do 300 j=1,nlat+1
      do 300 i=0,nlon+1
        umask(i,j) = 1.
  300 continue
 
      do 310 j=1,nlat
        do 312 i=1,nlon
          if(lmask(i,j).le.2) then
            umask(i,j) = 0.
            umask(i+1,j) = 0.
            umask(i,j+1) = 0.
            umask(i+1,j+1) = 0.
          endif
  312   continue
  310 continue
 
c At this point, umask on the Greenwich meridian (i=1 on u-grid) has
c only been affected by the eastern h-grid box. Remedy this here.
c Also set umask wraparound.
 
      do 320 j=1,nlat+1
        umask(1,j) = min (umask(1,j), umask(nlon+1,j))
  320 continue
 
      do 330 j=1,nlat+1
        umask(0,j) = umask(nlon,j)
        umask(nlon+1,j) = umask(1,j)
  330 continue
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine timefilt (fuh, fvh, fu, fv, dt)
 
c Transfers h-grid fields fuh,fvh to the u-grid, and applies this input
c to fu,fv using a running-mean time smoother.
 
c fuh = eastward  wind stress on h-grid (supplied)
c fvh = northward wind stress on h-grid (supplied)
c fu  = eastward  time-smoothed wind stress on h-grid (returned)
c fv  = northward time-smoothed wind stress on h-grid (returned)
c dt  = dynamic seaice time step, seconds
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  fuh(nlon,nlat),             fvh(nlon,nlat),
     *  fu (nlon,nlat),             fv (nlon,nlat)
 
c Timescale of the running time smoother, in days
 
      parameter (tfilt = 4.)
 
c Transfer to temporary work arrays and apply the time smoother by
c linear weighting with previous values. The weight of the n-th
c timestep in the past turns out to be wt * (1-wt)**n. (The running
c mean arrays fu,fv have been initialized or read from the restart file
c and passed by the calling program.)
 
      wt = dt/(tfilt*86400.)
 
      do 100 j=1,nlat
      do 100 i=1,nlon
        fu(i,j) = (1.-wt)*fu(i,j) + wt*fuh(i,j)
        fv(i,j) = (1.-wt)*fv(i,j) + wt*fvh(i,j)
 100  continue
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine polefilt (u, v, latslab1, icode)
 
c Spatially filters wind stresses or ocean currents near North Pole
c towards pure slab motion.
 
c u        = eastward  stress/current on u-grid (modified)
c v        = northward stress/current on u-grid (modified)
c latslab1 = lat (deg) northward of which to apply filter (inclusive)
c icode    = 1 for seaice wind stresses,
c            2 for open-ocean wind stresses,
c            3 for ocean currents,
c            4 for ice velocities.
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  u(0:nlon+1,nlat+1),         v(0:nlon+1,nlat+1)
 
c     umaxpd = max zonal-mean stress/speed per deg from NP
c     wslab1 = filtered vs unfiltered weight at latslab1, (=1 at NP)
 
      dimension umaxpd(4), utmp(nlon), vtmp(nlon)
      save umaxpd, wslab1
      data umaxpd /.02, .02, .01, .01/
      data wslab1 /0.5/
 
 
c Set latitudinal boundaries for the filter. If outside current
c hemispheric range of ice, do nothing.
 
      jslab1 = max (latua, nint (nlat*(latslab1+90.)/180. + 1.))
      jslab2 = min (latub, nlat+1)
      if (jslab1.gt.jslab2) return
 
 
c Calculate mean slab motion. Weight all pts equally, not by true area,
c since loosely speaking the more polar values are more "significant".
 
      uslab = 0.
      vslab = 0.
 
      do 100 j=jslab1,jslab2
        do 102 i=1,nlon
          utmp(i) = -u(i,j)*slonu(i) - v(i,j)*clonu(i)
          vtmp(i) =  u(i,j)*clonu(i) - v(i,j)*slonu(i)
  102   continue
        uslab = uslab + ssum(nlon, utmp, 1)
        vslab = vslab + ssum(nlon, vtmp, 1)
  100 continue
 
      uslab = uslab / (nlon*(jslab2-jslab1+1))
      vslab = vslab / (nlon*(jslab2-jslab1+1))
 
 
c Calculate zonal mean eastward and northward values uav,vav, and
c constrain them by umax. Set the final field, weighted with the
c unfiltered field by w, which varies linearly with latitude to 1 at NP.
c The filtered field is pure slab motion + uav (rot) + vav(conv).
c Note if latslab1 = 90, umax = uav = vav = 0.
 
      do 200 j=jslab1,jslab2
 
        uav = 0.
        vav = 0.
        do 202 i=1,nlon
          uav = uav + u(i,j)
          vav = vav + v(i,j)
  202   continue
        uav = uav / nlon
        vav = vav / nlon
 
        umax = umaxpd(icode) * (alatu(nlat+1)-alatu(j)) * (180./pi)
        uav = min ( umax, max(-umax,uav) )
        vav = min ( umax, max(-umax,vav) )
 
        if (jslab1.lt.jslab2) then
          w = (1.*(j-jslab1) + wslab1*(jslab2-j)) / (jslab2-jslab1)
        else
          w = 1.
        endif
 
        zslat = abs(slatu(j))
 
        do 204 i=1,nlon
          u(i,j) = w * ((-uslab*slonu(i) + vslab*clonu(i))*zslat + uav)
     *           + (1.-w)*u(i,j)
          v(i,j) = w * ( -uslab*clonu(i) - vslab*slonu(i)        + vav)
     *           + (1.-w)*v(i,j)
  204   continue
 
  200 continue
 
 
      call wrapu(u)
      call wrapu(v)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine icebound (fice, ihem, dynamhem)
 
c For the current hemisphere, finds the latitudinal indices bounding
c the extent of ice. These latitudinal indices are latha,lathb for the
c h-grid, and latua,latub for the u-grid (in common dicegrid). If ice
c doesn't extend to the poles, the southernmost ice point is at latha+1,
c the northernmost ice point is at lathb-1, latua is the southern
c boundary of the ice (=latha+1), and latub is the northernmost boundary
c of the ice (=lathb). If ice does extend to the poles, then latha=1
c (89S) or lathb=nlat(89N), and latua=1 (90S) or latub=nlat+1 (90N).
 
c (Note that extra h-grid points are set in ltoh at latha-1 and lathb+1,
c which for the polar case are at j=0 (91S) or j=nlat+1 (91N). So in all
c cases latua,latub are "inside" the region of defined h-grid points.)
 
c fice = ice concentration, on h-grid (supplied)
c ihem = 1 for southern hemisphere, 2 for northern hemisphere (supp)
c dynamhem = T for hemispheric, F for global
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension fice(nlon,nlat)
      logical dynamhem
      dimension ifany (nlat)
 
      if (dynamhem) then
        if (ihem.eq.1) then   ! SH or NH hemispheric domain
          nla = 1
          nlb = nlat/2
        else
          nla = nlat/2 + 1
          nlb = nlat
        endif
      else
        nla = 1               ! global domain, needed for ice at equator
        nlb = nlat
      endif
 
 
c First set ifany, flag array for whether ice exists at each latitude
 
      do 100 j=nla,nlb
        ifany(j) = 0
        do 102 i=1,nlon
          if (fice(i,j).gt.0.) then
            ifany(j) = 1
            goto 100
          endif
  102   continue
  100 continue
 
 
c Search northward from southernmost latitude, and set latha,latua.
c If no ice in this hemisphere, set latha=0 and return.
 
      do 200 j=nla,nlb
        if (ifany(j).eq.1) then
          latha = max (1, j-1)
          latua = latha + 1
          goto 202
        endif
  200 continue
      latha = 0
      return
  202 if (latha.eq.1) latua = 1
 
 
c Search southward from northernmost latitude, and set lathb,latub
 
      do 300 j=nlb,nla,-1
        if (ifany(j).eq.1) then
          lathb = min (nlat, j+1)
          latub = lathb
          goto 302
        endif
  300 continue
  302 if (lathb.eq.nlat) latub = nlat+1
 
      if ( (ihem.eq.1 .and. lathb.eq.nlat/2+1) .or.
     *     (ihem.eq.2 .and. latha.eq.nlat/2-1) ) then
        write(iutabe,900) ihem
  900   format(/' *** Warning. icebound 900. ice at equator. ihem=',i3/)
      endif
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine icefree (u, v, hiceu, fice, fuwindis, fvwindis,
     *                    uocn, vocn, pres, omega)
 
c Solves for free-drift velocities u and v.
c Note that total stress between agcm and ice is not conserved due
c to earlier transfer from h-grid f[u,v]windi to u-grid f[u,v]windis.
c (Could fix by weighting with fice before transfer and divide by ficeu
c after, but that would complicate polar filtering and time-smoothing.)
 
c u        = eastward  velocity on u-grid (returned)
c v        = northward velocity on u-grid (returned)
c hiceu    = total ice thickness on u-grid (supplied)
c fice     = ice concentration, on h-grid (supplied)
c fuwindis = eastward  wind stress on u-grid (supplied)
c fvwindis = northward wind stress on u-grid (supplied)
c uocn     = eastward  surface ocean velocity on u-grid (supplied)
c vocn     = northward surface ocean velocity on u-grid (supplied)
c pres     = internal ice "pressure" on h-grid (supplied)
c omega    = 2*pi/day lenth for coriolis force, rad/s (supplied)
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  u(0:nlon+1,nlat+1),           v(0:nlon+1,nlat+1),
     *  hiceu(0:nlon+1,nlat+1),       fice(nlon,nlat),
     *  fuwindis(0:nlon+1,nlat+1),    fvwindis(0:nlon+1,nlat+1),
     *  uocn(0:nlon+1,nlat+1),        vocn(0:nlon+1,nlat+1),
     *  pres(nlon,nlat)
 
 
c Transfer internal ice pressure to a wrapped h-grid array, since its
c gradients on the u-grid are computed below
 
      call ltoh (pres, workh)
 
 
c Solve free-drift equations
 
      do 100 j=latua,latub
        a = csubw*coswturn(j)
 
        do 102 i=1,nlon
 
c            Compute gradients of internal ice pressure
 
          prese = 0.5* ( (workh(i,j)  -workh(i-1,j)  )/deltx(j)
     *                  +(workh(i,j-1)-workh(i-1,j-1))/deltx(j-1) )
          presn = 0.5* ( (workh(i,j)  -workh(i,j-1)  )/deltyu(j)
     *                  +(workh(i-1,j)-workh(i-1,j-1))/deltyu(j) )
 
c            Compute gradients of dynamic topography assuming geostrophy
 
          slope =  2.*omega*slatu(j)*vocn(i,j) / grav
          slopn = -2.*omega*slatu(j)*uocn(i,j) / grav
 
c            add up all forcing terms (those indep of u and v)
 
          uforce = fuwindis(i,j)
     *           - prese
     *           - grav*rhoi*hiceu(i,j)*slope
     *           + csubw*coswturn(j)*uocn(i,j)
     *           - csubw*sinwturn(j)*vocn(i,j)
 
          vforce = fvwindis(i,j)
     *           - presn
     *           - grav*rhoi*hiceu(i,j)*slopn
     *           + csubw*coswturn(j)*vocn(i,j)
     *           + csubw*sinwturn(j)*uocn(i,j)
 
c            solve the 2 linearized equations for u and v
 
          b = rhoi*hiceu(i,j)*2.*omega*slatu(j) + csubw*sinwturn(j)
          u(i,j) = umask(i,j) * (vforce*b + uforce*a) / (a*a+b*b)
          v(i,j) = umask(i,j) * (vforce*a - uforce*b) / (a*a+b*b)
 
  102   continue
  100 continue
 
 
c Set wraparound u,v points
 
      call wrapu (u)
      call wrapu (v)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine cavset (hiceu,prespd,du1,dv1,du2,dv2,du3,dv3,du4,dv4,
     *                   omega)
 
c Set prespd = 1/divergence and du*,dv* = velocs at corners, per unit
c pressure at a h-grid box with zero pressure elsewhere, for use in
c cavit. If divergence is zero (on land, or surrounded by land), prespd
c is set to zero. The global array for prespd is only necessary since
c we account for zero veloc at one or more corners (umask). The global
c arrays du*,dv* are only necessary since we include the coriolis terms
c (b[1-4] below) which depend on hiceu (Flato doesn't include them).
 
c hiceu    = total ice thickness on u-grid (supplied)
c prespd   = 1/divergence per unit "delta" pressure on h-grid (returned)
c du1      = eastward  veloc at SW corner per unit pressure (returned)
c dv1      = northward veloc at SW corner per unit pressure (returned)
c du2      = eastward  veloc at SE corner per unit pressure (returned)
c dv2      = northward veloc at SE corner per unit pressure (returned)
c du3      = eastward  veloc at NE corner per unit pressure (returned)
c dv3      = northward veloc at NE corner per unit pressure (returned)
c du4      = eastward  veloc at NW corner per unit pressure (returned)
c dv4      = northward veloc at NW corner per unit pressure (returned)
c omega    = 2*pi/day lenth for coriolis force, rad/s (supplied)
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  hiceu(0:nlon+1,nlat+1),
     *  prespd(nlon,nlat),
     *  du1(nlon,nlat),              dv1(nlon,nlat),
     *  du2(nlon,nlat),              dv2(nlon,nlat),
     *  du3(nlon,nlat),              dv3(nlon,nlat),
     *  du4(nlon,nlat),              dv4(nlon,nlat)
 
      dimension divinc(nlon)
 
 
c Calculate corner-velocities and divergence for a unit "delta" pressure
c at each h-grid point, with zero pressure elsewhere. Algebra is same as
c in icefree. im,ip and jm,jp are indices on u-grid.
 
      do 100 j=latha,lathb
        jm = j
        jp = j+1
        am = csubw*coswturn(jm)
        ap = csubw*coswturn(jp)
        pem =  0.5/deltx(j)
        pep = -0.5/deltx(j)
        pnm =  0.5/deltyu(jm)
        pnp = -0.5/deltyu(jp)
 
        do 102 i=1,nlon
          im = i
          ip = i+1
 
          b1 = rhoi*hiceu(im,jm)*2.*omega*slatu(jm) + csubw*sinwturn(jm)
          b2 = rhoi*hiceu(ip,jm)*2.*omega*slatu(jm) + csubw*sinwturn(jm)
          b3 = rhoi*hiceu(ip,jp)*2.*omega*slatu(jp) + csubw*sinwturn(jp)
          b4 = rhoi*hiceu(im,jp)*2.*omega*slatu(jp) + csubw*sinwturn(jp)
 
          du1(i,j) = umask(i,j)     * (-pnm*b1 - pem*am) / (am*am+b1*b1)
          dv1(i,j) = umask(i,j)     * (-pnm*am + pem*b1) / (am*am+b1*b1)
          du2(i,j) = umask(i+1,j)   * (-pnm*b2 - pep*am) / (am*am+b2*b2)
          dv2(i,j) = umask(i+1,j)   * (-pnm*am + pep*b2) / (am*am+b2*b2)
          du3(i,j) = umask(i+1,j+1) * (-pnp*b3 - pep*ap) / (ap*ap+b3*b3)
          dv3(i,j) = umask(i+1,j+1) * (-pnp*ap + pep*b3) / (ap*ap+b3*b3)
          du4(i,j) = umask(i,j+1)   * (-pnp*b4 - pem*ap) / (ap*ap+b4*b4)
          dv4(i,j) = umask(i,j+1)   * (-pnp*ap + pem*b4) / (ap*ap+b4*b4)
 
          divinc(i) = ( - 0.5 * (du1(i,j) + du4(i,j)) * delty(j)
     *                  + 0.5 * (du2(i,j) + du3(i,j)) * delty(j)
     *                  - 0.5 * (dv1(i,j) + dv2(i,j)) * deltxu(j)
     *                  + 0.5 * (dv3(i,j) + dv4(i,j)) * deltxu(j+1)
     *                ) / srfarea(j)
  102   continue
 
 
c Invert to get 1/divergence, or set to 0 if divinc is 0. (Do this
c here so don't have to do it for every iteration in cavit.)
 
        do 104 i=1,nlon
          prespd(i,j) = cvmgt ( 1./max(divinc(i),1.e-20),
     *                          0., divinc(i).gt.0. )
  104   continue
 
  100 continue
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine cavit (u, v, pres, streng,
     *                  prespd,du1,dv1,du2,dv2,du3,dv3,du4,dv4)
 
c Accounts for ice rheology using cavitating-fluid method. In each
c iteration, nothing is done for divergent h-grid points, and for
c convergent h-grid points a cancelling amount of divergence is added
c by incrementing "pressure" an appropriate amount for each h-grid point
c in isolation. The amount of pressure per divergence and the resulting
c corner-velocs u,v have been pre-computed in cavset. The resulting
c pressure and u,v fields continue to obey the (linear) momentum eqns.
c The pressure is not allowed to exceed a maximum value ("strength")
c that is a function of ice thickness and concentration - if the maximum
c pressure is reached, the flow is allowed to converge ("collapse").
 
c Since the increments in u and v affect convergence at previously
c adjusted h-grid points, the process must be iterated globally.
c Each iteration reduces kinetic energy but conserves momentum (at
c least for uniform hice,fice on a cartesian grid) and so convergence
c is guaranteed.
 
c u      = eastward  velocity on u-grid (modified)
c v      = northward velocity on u-grid (modified)
c pres   = internal ice "pressure" on h-grid (modified)
c streng = internal ice strength on h-grid (supplied)
c prespd = 1/divergence per unit "delta" pressure on h-grid (supplied)
c du1    = eastward  veloc at SW corner per unit pressure (supplied)
c dv1    = northward veloc at SW corner per unit pressure (supplied)
c du2    = eastward  veloc at SE corner per unit pressure (supplied)
c dv2    = northward veloc at SE corner per unit pressure (supplied)
c du3    = eastward  veloc at NE corner per unit pressure (supplied)
c dv3    = northward veloc at NE corner per unit pressure (supplied)
c du4    = eastward  veloc at NW corner per unit pressure (supplied)
c dv4    = northward veloc at NW corner per unit pressure (supplied)
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  u(0:nlon+1,nlat+1),     v(0:nlon+1,nlat+1),
     *  pres(nlon,nlat),        streng(nlon,nlat),
     *  prespd(nlon,nlat),
     *  du1(nlon,nlat),         dv1(nlon,nlat),
     *  du2(nlon,nlat),         dv2(nlon,nlat),
     *  du3(nlon,nlat),         dv3(nlon,nlat),
     *  du4(nlon,nlat),         dv4(nlon,nlat)
 
      dimension
     *  div(nlon),
     *  presinc(nlon,nlat),     presold(nlon),
     *  uwrapsav(nlat+1),       vwrapsav(nlat+1)
 
 
c Calculate the divergence, paralleling calculations in advect,
c i.e, using 2-point differences and using spherical metrics delt*
 
      do 100 j=latha,lathb
        do 102 i=1,nlon
          div(i) = ( - 0.5 * (u(i,j)   + u(i,j+1))   * delty(j)
     *               + 0.5 * (u(i+1,j) + u(i+1,j+1)) * delty(j)
     *               - 0.5 * (v(i,j)   + v(i+1,j))   * deltxu(j)
     *               + 0.5 * (v(i,j+1) + v(i+1,j+1)) * deltxu(j+1)
     *             ) / srfarea(j)
 
c Increment the pressure where the flow is convergent. If the
c accumulated pressure would exceed the maximum, reduce the increment
c and allow some convergence (ie, collapse). relax is an underrelaxation
c factor set in comdice.
 
          presold(i) = pres(i,j)
          pres(i,j)= min( pres(i,j) + relax*max(-div(i),0.)*prespd(i,j),
     *                    streng(i,j) )
          presinc(i,j) = pres(i,j) - presold(i)
  102   continue
  100 continue
 
 
c--------------------
      if (flato) then
c--------------------
 
c If "flato" method, apply corner-velocity increments here using unit
c veloc increments pre-computed in cavset. If "not flato" method,
c this will be done by a call to icefree.
 
c Save u and v values at nlon+1, for wraparound in loop 220 below
 
        do 200 j=latua,latub
          uwrapsav(j) = u(nlon+1,j)
          vwrapsav(j) = v(nlon+1,j)
  200   continue
 
 
c This calculation of the u,v increments is an unfolding of the
c divergence calculation.
 
        do 210 j=latha,lathb
          do 212 i=1,nlon
            u(i,j)     = u(i,j)     + du1(i,j)*presinc(i,j)
            v(i,j)     = v(i,j)     + dv1(i,j)*presinc(i,j)
            u(i+1,j)   = u(i+1,j)   + du2(i,j)*presinc(i,j)
            v(i+1,j)   = v(i+1,j)   + dv2(i,j)*presinc(i,j)
            u(i+1,j+1) = u(i+1,j+1) + du3(i,j)*presinc(i,j)
            v(i+1,j+1) = v(i+1,j+1) + dv3(i,j)*presinc(i,j)
            u(i,j+1)   = u(i,j+1)   + du4(i,j)*presinc(i,j)
            v(i,j+1)   = v(i,j+1)   + dv4(i,j)*presinc(i,j)
  212     continue
  210   continue
 
 
c Allow for the east-west boundary, since the westernmost u-grid points
c (1,j) have not yet been influenced by the easternmost h-grid boxes
c (nlon,j). [u,v]wrapsav hold the previous iteration's u,v at nlon+1.
 
        do 220 j=latua,latub
          u(1,j) = u(1,j) + u(nlon+1,j) - uwrapsav(j)
          v(1,j) = v(1,j) + v(nlon+1,j) - vwrapsav(j)
 
c       also set wraparounds
 
          u(nlon+1,j) = u(1,j)
          v(nlon+1,j) = v(1,j)
 
          u(0,j) = u(nlon,j)
          v(0,j) = v(nlon,j)
  220   continue
 
c----------
      endif
c----------
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine icediag (u,v,fice,ficeu,pres,streng,iloop,ifdone,ihem)
 
c (i)   Calculates rmsinc, the u,v increment for convergence test
c (ii)  Tests if converged or have reached max iteration number
c (iii) If last iteration, calculates diagnostics (in common comdiag)
 
c u      = eastward  velocity on u-grid (supplied)
c v      = northward velocity on u-grid (supplied)
c fice   = ice concentration on h-grid (supplied)
c ficeu  = ice concentration on u-grid (supplied)
c pres   = internal ice "pressure" on h-grid (supplied)
c streng = internal ice strength on h-grid (supplied)
c iloop  = iteration number (supplied)
c ifdone = flag, set to .true. if last iteration (returned)
c ihem   = 1 for SH, 2 for NH
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *   u(0:nlon+1,nlat+1),           v(0:nlon+1,nlat+1),
     *   fice(nlon,nlat),              ficeu(0:nlon+1,nlat+1),
     *   pres(nlon,nlat),              streng(nlon,nlat)
 
      logical ifdone
 
 
      ifdone = .false.
 
 
c Calculate max velocity increment for this iteration (old test on
c rms ice-area-weighted velocity is commented out). Also calculate
c tareau = total ice area.
 
      rmsinc(ihem) = 0.
      tareau = 0.
      do 100 j=latua,latub
        do 102 i=1,nlon
          dvel2 = (u(i,j)-worku(i,j))**2 + (v(i,j)-worku2(i,j))**2
c         local max:
          rmsinc(ihem) = max ( rmsinc(ihem), dvel2 )
c         rms ice-area weighted:
c         rmsinc(ihem) = rmsinc(ihem) + dvel2*srfareau(j)*ficeu(i,j)
 
          tareau = tareau + srfareau(j)*ficeu(i,j)
  102   continue
  100 continue
c     local max:
      rmsinc(ihem) = sqrt(rmsinc(ihem))
c     rms ice-area weighted:
c     if (tareau.gt.0.) rmsinc(ihem) = sqrt(rmsinc(ihem)/tareau)
 
 
c     If last loop, set ifdone flag and calc diagnostics (in comdiag)
 
      if ( (rmsinc(ihem).le.convinc .and. iloop.ge.nloopmin)
     *     .or. iloop.eq.nloopmax ) then
 
        ifdone = .true.
 
        loop(ihem) = iloop
 
        rmsvel(ihem) = 0.
        do 200 j=latua,latub
        do 200 i=1,nlon
          rmsvel(ihem) = rmsvel(ihem)
     *           + (u(i,j)**2 + v(i,j)**2) * srfareau(j)*ficeu(i,j)
  200   continue
        if (tareau.gt.0.) rmsvel(ihem) = sqrt(rmsvel(ihem)/tareau)
 
        do 210 j=latha,lathb
        do 210 i=1,nlon
          icolaps(i,j) = -1
          if (fice(i,j).gt.0.) then
            icolaps(i,j) = 0
            if (pres(i,j).gt.0.) icolaps(i,j) = 1
            if (pres(i,j).ge.streng(i,j)) icolaps(i,j) = 2
          endif
  210   continue
 
 
c       Warn only if failed to get below 2*convinc...below that but
c       above convinc is not too bad
 
        if (rmsinc(ihem) .gt. 2*convinc)
     *  write(iutabe,220) ihem, iloop, nloopmax, rmsinc(ihem), 2*convinc
  220   format(/' *** Warning. icediag 220. Bad cavit convergence.',
     *          ' ihem,iloop,nloopmax,rmsinc,2*convinc=',3i4,2e15.5/)
 
      endif
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine advect (result, field, factor, u, v, dt, vmin, icode)
 
c Advects an h-grid field "field" by u-grid velocities (u,v) through
c time dt. Divides advected field by "factor" to give "result".
 
c The flux-form used here guarantees conservation of the advected
c fields. Note that "upstream" values are advected - which guarantees
c that recovered thicknesses and temperatures will have reasonable
c values (as long as the CFL condition is satisfied so no fields become
c negative). Given the CFL condition, ice concentration (and other
c fields) will approach but never reach zero in a continuously divergent
c grid box - kill off small ice quantities elsewhere in the code.
 
c All arguments are supplied except "result", and on h-grid except u,v.
c result = final advected result (not necessarily the advected field)
c field  = pre-advection field to be advected
c factor = post-advection qty to divide advected field to get result
c u      = eastward  velocity (u-grid)
c v      = northward velocity (u-grid)
c dt     = time step
c vmin   = value below which "result" = 0 (for ice,snow fractions)
c icode  =  0 for ice concentration,
c           1 for ice mass,
c           2 for ice heat,
c           3 for ice brine resevoir,
c          10 for snow concentration,
c          11 for snow mass,
c          12 for snow heat.
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  result(nlon,nlat),     field(nlon,nlat),    factor(nlon,nlat),
     *  u(0:nlon+1,nlat+1),    v(0:nlon+1,nlat+1)
 
 
c Transfer h-grid field to be advected to a wrapped h-grid array
 
      call ltoh (field, workh)
 
 
c Do advection
 
      do 400 j=latha,lathb
      do 400 i=1,nlon
        ul = 0.5 * (u(i,j)   + u(i,j+1))
        ur = 0.5 * (u(i+1,j) + u(i+1,j+1))
        vb = 0.5 * (v(i,j)   + v(i+1,j))
        vt = 0.5 * (v(i,j+1) + v(i+1,j+1))
 
        fl = max(ul,0.)*workh(i-1,j) + min(ul,0.)*workh(i,j)
        fr = max(ur,0.)*workh(i,j)   + min(ur,0.)*workh(i+1,j)
        fb = max(vb,0.)*workh(i,j-1) + min(vb,0.)*workh(i,j)
        ft = max(vt,0.)*workh(i,j)   + min(vt,0.)*workh(i,j+1)
 
        field(i,j) = workh(i,j) + dt
     *             * (  fl*delty(j)  - fr*delty(j)
     *                + fb*deltxu(j) - ft*deltxu(j+1) ) / srfarea(j)
  400 continue
 
 
c Optionally smear latitude bands nearest poles (on h-grid, j=1 or nlat)
 
c     do 450 j=1,nlat,nlat-1
c       if (latha.le.j .and. j.le.lathb) then
c         zf = 0.
c         do 452 i=1,nlon
c           zf = zf + field(i,j)
c 452     continue
c         zf = zf/nlon
c         do 454 i=1,nlon
c           field(i,j) = zf
c 454     continue
c       endif
c 450 continue
 
 
c For the primary (first-advected) field of ice concentration, set
c "result" to zero if below a tiny limit. Small amounts of ice will be
c killed anyway at the end of dynice, but the advection code should
c still guard against machine-limit problems. The limit here (0.1*vmin)
c is smaller than the later "kill-amount" (vmin=ficemin, passed from
c dynice) and so only affects creation in a previously ice-free box;
c since pre-advected fice (for existing ice) is always >= ficemin, the
c advected "field" will be substantially above .1*ficemin as long as
c the "CFL" criterion is well-satisfied.
 
      if (icode.eq.0) then
        do 500 j=latha,lathb
        do 500 i=1,nlon
          result(i,j) = field(i,j)
          if (result(i,j) .lt. .1*vmin) result(i,j) = 0.
  500   continue
 
c For snow concentration, similar action as for fice above, with passed
c vmin=fsnomin. But if "factor" (= advected fice) is zero, don't do
c anything...snow and ice will be killed anyway at end of dynice.
 
      else if (icode.eq.10) then
        do 510 j=latha,lathb
        do 510 i=1,nlon
          if (factor(i,j).gt.0.) then
            result(i,j) = field(i,j) / factor(i,j)
            if (result(i,j) .lt. .1*vmin) result(i,j) = 0.
          endif
  510   continue
 
c For all other fields, just divide by "factor" (=fice*...) if non-zero.
c "factor" is zero for no ice (possibly cut off by small limit in
c loop 500), or for land points, where we mustn't change "result"
c since that would corrupt soil thicknesses and land snow amounts.
 
      else
        do 520 j=latha,lathb
        do 520 i=1,nlon
          if (factor(i,j).gt.0.) result(i,j) = field(i,j) / factor(i,j)
  520   continue
      endif
 
 
c For  ice and snow concentrations, impose maximum of 1.
 
      if (icode.eq.0 .or. icode.eq.10) then
        do 600 j=latha,lathb
        do 600 i=1,nlon
          result(i,j) = min (result(i,j), 1.)
  600   continue
      endif
 
 
c Check - advection should never create negative values
 
#ifndef gcm
      do 900 j=latha,lathb
      do 900 i=1,nlon
        if (field(i,j).lt.0.) write(iutabe,902) i,j,icode,
     *  field(i,j),workh(i,j),factor(i,j),result(i,j)
  902   format(/' *** Warning. advect. i,j,icode,',
     *          'field,workh,factor,result:',3i4,4e13.4/)
  900 continue
#endif
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine iceprint (hicelay, hice, fice,fuwindis,fvwindis,
     *                     uocn, vocn, u, v, pres, streng, lmask,
     *                     dt, secyear, secrun)
 
c Dumps out printed grids of various arrays for a selected region
c at a selected time interval. Local arrays *mean can be used to
c accumulate and output time means instead of instantaneous values.
c Writes to unit iutabe (in comgrd).
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension
     *  hicelay(nlon,nlat,nicelay),
     *  hice(nlon,nlat),             fice(nlon,nlat),
     *  fuwindis(0:nlon+1,nlat+1),   fvwindis(0:nlon+1,nlat+1),
     *  uocn(0:nlon+1,nlat+1),       vocn(0:nlon+1,nlat+1),
     *  u(0:nlon+1,nlat+1),          v(0:nlon+1,nlat+1),
     *  pres(nlon,nlat),             streng(nlon,nlat),
     *  lmask(nlon,nlat)
 
c Local arrays for time-mean accumulations (*mean) and other diagnostic
c calcs. Model global grid (skipping i's, shifted in i and j) will be
c mapped to continuous smaller display grid (ii,jj), so *mean arrays
c don't have to be global.
 
      parameter (nlonmean=40, nlatmean=40)
      dimension
     *  fumean(nlonmean,nlatmean),   fvmean(nlonmean,nlatmean),
     *  umean(nlonmean,nlatmean),    vmean(nlonmean,nlatmean),
     *  dmean(nlonmean,nlatmean)
 
      dimension ndmo(12)
      data ndmo /31,59,90,120,151,181,212,243,273,304,334,365/
 
      logical first, ifdmo
      save fumean, fvmean, umean, vmean, dmean, ndmo, first, kount,
     *     i1, i2, i3, i4, iu1, iu2, iu3, iu4,
     *     j1, j2, j3, j4, ju1, ju2, ju3, ju4
      data first /.true./
 
 
      if (first) then
 
c       Set output-region grid indices from requested long,lat degrees
c       (i* for h-grid, iu* for u-grid, j* for h-grid, ju* for u-grid)
 
c8880   Arctic, Canadian Archipelago
 8881   lonlef = 235
 8882   lonrit = 305
 8884   lattop = 89
 8885   latbot = 71
 
c8880   Arctic, Soviet peninsula
c8881   lonlef = 65
c8882   lonrit = 135
c8884   lattop = 89
c8885   latbot = 71
 
c8880   Cape Adair, Antarctica
c8881   lonlef = 129
c8882   lonrit = 199
c8884   lattop = -71
c8885   latbot = -77
 
        i1 = max (1, min (nlon, nint(( lonlef     /360.)*nlon + 0.5) ))
        i2 = max (1, min (nlon, nint(( lonrit     /360.)*nlon + 0.5) ))
        j1 = max (1, min (nlat, nint(((lattop+90.)/180.)*nlat + 0.5) ))
        j2 = max (1, min (nlat, nint(((latbot+90.)/180.)*nlat + 0.5) ))
 
        i3 = (i2-i1)/36 + 1
        iu1 = i1+1
        iu2 = i2+1
        iu3 = i3
 
        j3 = sign (1,j2-j1)
        ju1 = j1
        ju2 = j2+1
        ju3 = j3
 
        i2  = min (i2,  i1 +i3 *(nlonmean-1))
        iu2 = min (iu2, iu1+iu3*(nlonmean-1))
        j2  = min (j2,  j1 +j3 *(nlatmean-1))
        ju2 = min (ju2, ju1+ju3*(nlatmean-1))
 
        i4  = (i2 -i1 )/i3  + 1
        j4  = (j2 -j1 )/j3  + 1
        iu4 = (iu2-iu1)/iu3 + 1
        ju4 = (ju2-ju1)/ju3 + 1
 
c       Initialize the time-mean accumulation counter
 
        kount = 0
 
        first = .false.
 
      endif
 
 
c Zero accumulators if first call, or if printout on previous call
 
      if (kount.eq.0) then
        do 5 jj=1,ju4
        do 5 ii=1,iu4
          fumean(ii,jj) = 0.
          fvmean(ii,jj) = 0.
          umean(ii,jj) = 0.
          vmean(ii,jj) = 0.
          dmean(ii,jj) = 0.
    5   continue
      endif
 
 
c If not printout time, accumulate time means, increment kount, return.
c If printout time, normalize time means and set kount to zero.
 
      ifdmo = .false.
#ifdef gcm
      do 8 m=1,12
        if ( abs(secyear-86400.*ndmo(m))/dt .lt. 0.5 ) ifdmo = .true.
    8 continue
#else
      if (mod(nint(secyear),86400).eq.0) ifdmo = .true.
#endif
 
      if ((.not.ifdmo) .or. kount.eq.0) then
 
        do 10 j=ju1,ju2,ju3
        jj = (j-ju1)/ju3 + 1
        do 10 i=iu1,iu2,iu3
          ii = (i-iu1)/iu3 + 1
          fumean(ii,jj) = fumean(ii,jj) + fuwindis(i,j)
          fvmean(ii,jj) = fvmean(ii,jj) + fvwindis(i,j)
          umean(ii,jj) = umean(ii,jj) + u(i,j)
          vmean(ii,jj) = vmean(ii,jj) + v(i,j)
   10   continue
        do 12 j=j1,j2,j3
        jj = (j-j1)/j3 + 1
        do 12 i=i1,i2,i3
          ii = (i-i1)/i3 + 1
          dmean(ii,jj) = dmean(ii,jj)
     *                 + ( - 0.5* (u(i,j)   + u(i,j+1))   * delty(j)
     *                     + 0.5* (u(i+1,j) + u(i+1,j+1)) * delty(j)
     *                     - 0.5* (v(i,j)   + v(i+1,j))   * deltxu(j)
     *                     + 0.5* (v(i,j+1) + v(i+1,j+1)) * deltxu(j+1)
     *                   ) / srfarea(j)
   12   continue
        kount = kount + 1
        return
 
      else
 
        do 20 j=ju1,ju2,ju3
        jj = (j-ju1)/ju3 + 1
        do 20 i=iu1,iu2,iu3
          ii = (i-iu1)/iu3 + 1
          fumean(ii,jj) = fumean(ii,jj) / kount
          fvmean(ii,jj) = fvmean(ii,jj) / kount
          umean(ii,jj) = umean(ii,jj) / kount
          vmean(ii,jj) = vmean(ii,jj) / kount
   20   continue
        do 22 j=j1,j2,j3
        jj = (j-j1)/j3 + 1
        do 22 i=i1,i2,i3
          ii = (i-i1)/i3 + 1
          dmean(ii,jj) = dmean(ii,jj) / kount
   22   continue
        kount = 0
 
      endif
 
 
c Calculate hemispheric means of fice and fice*hice, write header
 
      tares = 0.
      thics = 0.
      do 44 j=1,nlat/2
      do 44 i=1,nlon
        tares = tares + fice(i,j)*srfarea(j)
        thics = thics + fice(i,j)*hice(i,j)*srfarea(j)
   44 continue
      tares = tares / (2.*pi*earthr**2)
      thics = thics / (2.*pi*earthr**2)
 
      taren = 0.
      thicn = 0.
      do 46 j=nlat/2+1,nlat
      do 46 i=1,nlon
        taren = taren + fice(i,j)*srfarea(j)
        thicn = thicn + fice(i,j)*hice(i,j)*srfarea(j)
   46 continue
      taren = taren / (2.*pi*earthr**2)
      thicn = thicn / (2.*pi*earthr**2)
 
      write(iutabe,50)
     *  secyear/86400., 1 + int((secrun-1.)/(.864e5*365)),
     *  taren, tares, thicn, thics,
     *  (loop(ihem), rmsinc(ihem), rmsvel(ihem), ihem=1,2)
#ifndef gcm
      write(ioterm,50)
     *  secyear/86400., int(1+(secrun-1)/(.864e5*365)),
     *  taren, tares, thicn, thics,
     *  (loop(ihem), rmsinc(ihem), rmsvel(ihem), ihem=1,2)
#endif
   50 format(/1x,79('-')
     *      //' day =',f8.3,'   year =',i3,
     *                      '   area (nh) =',f8.4,'   area (sh) =',f8.4,
     *                      '   vol (nh) =', f8.4,'   vol (sh) =', f8.4
     *      //' sh:  loop =',i4,'   rmsinc =',e15.5,'   rmsvel =',e15.5
     *       /' nh:  loop =',i4,'   rmsinc =',e15.5,'   rmsvel =',e15.5)
 
      write(iutabe,52) (2*i-1,i=i1,i2,i3)
   52 format(/9x,36i3)
 
 
c Start of grid output for arrays
 
      do 1010 j=1,nlat
      do 1010 i=1,nlon
        if (lmask(i,j).ge.3) then
          work(i,j) = 10.*hice(i,j)
        else
          work(i,j) = 999.
        endif
 1010 continue
      write(iutabe,*) 'hice*10:'
      do 101 j=ju1,ju2,ju3
  101 write(iutabe,800) -91+2*j, (nint(work(i,j)),i=i1,i2,i3)
 
      do 1020 j=1,nlat
      do 1020 i=1,nlon
        if (fice(i,j).gt..005) then
          work(i,j) = 100.*fice(i,j)
        else if (fice(i,j).gt.1.e-8) then
          work(i,j) = -10000.*fice(i,j)
        else
          work(i,j) = 0.
        endif
 1020 continue
      write(iutabe,*) 'fice*100, -fice*10000:'
      do 102 j=j2,j1,j3
 102  write(iutabe,800) -91+2*j, (nint(work(i,j)),i=i1,i2,i3)
 
      write(iutabe,*) 'time-mean fuwindis*100:'
      do 103 jj=ju4,1,-1
      j = ju1 + ju3*(jj-1)
  103 write(iutabe,800) -92+2*j, (nint(fumean(ii,jj)*100.),ii=1,iu4)
 
      write(iutabe,*) 'time-mean fvwindis*100:'
      do 104 j=ju1,ju2,ju3
      jj = (j-ju1)/ju3 + 1
  104 write(iutabe,800) -92+2*j, (nint(fvmean(ii,jj)*100.),ii=1,iu4)
 
      write(iutabe,*) 'uocn*100:'
      do 105 j=ju1,ju2,ju3
  105 write(iutabe,800) -92+2*j, (nint(uocn(i,j)*100.),i=iu1,iu2,iu3)
 
      write(iutabe,*) 'vocn*100:'
      do 106 j=ju1,ju2,ju3
  106 write(iutabe,800) -92+2*j, (nint(vocn(i,j)*100.),i=iu1,iu2,iu3)
 
      write(iutabe,*) 'time-mean u*100:'
      do 107 j=ju1,ju2,ju3
      jj = (j-ju1)/ju3 + 1
  107 write(iutabe,800) -92+2*j, (nint(umean(ii,jj)*100.),ii=1,iu4)
 
      write(iutabe,*) 'time-mean v*100:'
      do 108 j=ju1,ju2,ju3
      jj = (j-ju1)/ju3 + 1
  108 write(iutabe,800) -92+2*j, (nint(vmean(ii,jj)*100.),ii=1,iu4)
 
      write(iutabe,*) 'time-mean div*1.e7:'
      do 109 j=j1,j2,j3
      jj = (j-j1)/j3 + 1
  109 write(iutabe,800) -91+2*j,(nint(dmean(ii,jj)*1.e7),ii=1,i4)
 
      write(iutabe,*) 'streng*1.e-3:'
      do 110 j=j1,j2,j3
  110 write(iutabe,800) -91+2*j,(nint(streng(i,j)*1.e-3),i=i1,i2,i3)
 
      write(iutabe,*) 'pres*1.e-3:'
      do 111 j=j2,j1,j3
  111 write(iutabe,800) -91+2*j,(nint(pres(i,j)*1.e-3),i=i1,i2,i3)
 
      write(iutabe,*) 'icolaps:'
      do 112 j=j2,j1,j3
  112 write(iutabe,800) -91+2*j, (icolaps(i,j),i=i1,i2,i3)
 
      return
 
  800 format(1x,i4,4x,36i3)
 
      end
c
c----------------------------------------------------------------------
c
      subroutine ltoh (al, ah)
 
c Copies lsx-grid (= h-grid with no wraparound) array al to h-grid
c array ah, sets ah's extra latitude points, and sets ah's
c longitudinal wraparound points.
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension al(nlon,nlat), ah(0:nlon+1,0:nlat+1)
 
 
c Copy to ah, including extra points at latha-1 and lathb+1 if not
c next to poles.
 
      do 5 j = max(latha-1,1), min (lathb+1,nlat)
       do 6 i=1,nlon
         ah(i,j) = al(i,j)
    6   continue
    5 continue
 
 
c If latha is next to south pole, set extra points ah(i,0).
 
      if (latha.eq.1) then
        do 10 i=1,nlon
          ah(i,0) = al(i,1)
   10   continue
      endif
 
c If lathb is next to north pole, set extra points ah(i,nlat+1).
 
      if (lathb.eq.nlat) then
        do 15 i=1,nlon
          ah(i,nlat+1) = al(i,nlat)
   15   continue
      endif
 
 
c Set longitudinal wraparound
 
      do 20 j=latha-1,lathb+1
        ah(0,j) = ah(nlon,j)
        ah(nlon+1,j) = ah(1,j)
   20 continue
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine ltou (al, au)
 
c Interpolates lsx-grid array al to u-grid array au, and sets au's
c wraparound points. Note that the interpolation is area weighted,
c by the use of the sine latitude scales slat,slatu.
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension  al(nlon,nlat),  au(0:nlon+1, nlat+1)
      dimension  ah(0:nlon+1,0:nlat+1)
 
      call ltoh (al, ah)
 
      do 10 j=latua,latub
        x1 = slat(j)  - slatu(j)
        x2 = slatu(j) - slat(j-1)
        w1 = 0.5*x1/(x1+x2)
        w2 = 0.5*x2/(x1+x2)
        do 12 i=1,nlon+1
          au(i,j) = w1*ah(i,j)   + w1*ah(i-1,j)
     *            + w2*ah(i,j-1) + w2*ah(i-1,j-1)
   12   continue
   10 continue
 
      call wrapu (au)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine wrapu (au)
 
c Sets longitudinal wraparound points for u-grid array au
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension au(0:nlon+1,nlat+1)
 
      do 10 j=latua,latub
        au(0,j) = au(nlon,j)
        au(nlon+1,j) = au(1,j)
  10  continue
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine getthick (hicelay, hice)
 
c Computes lsx-grid total ice thickness hice from lsx-grid
c layer thicknesses hicelay 
 
#include <compar.h>
#include <comgrd.h>
#include <comdice.h>
 
      dimension hicelay(nlon,nlat,nicelay), hice(nlon,nlat)
 
      do 100 j=1,nlat
      do 100 i=1,nlon
        hice(i,j) = 0.
  100 continue
      do 110 k=1,nicelay
        do 112 j=1,nlat
        do 112 i=1,nlon
          hice(i,j) = hice(i,j) + hicelay(i,j,k)
  112   continue
  110 continue
 
      return
      end
