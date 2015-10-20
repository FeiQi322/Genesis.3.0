c -- $Id: slth.F,v 1.2 2001/07/03 21:29:21 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/slth.F,v $
c -- $Name:  $

c
c-----------------------------------------------------------------------
c
      subroutine slthyb (wgaus, sigg, sigl, dsig,
     *                   gravit, radius, dtime,
     *                   ug, vg, wg, wgsav,
     *                   priestq, strattime, stratflag, smolar,
     *                   pthicm, pthicp,
     *                   qg, qgnew, tracer, ntrace, ntraca,
     *                   lonsav, latsav, sigsav,
     *                   totqo, totqn, totqa 
#ifdef isotrac
     *                  ,qg_o, qgnew_o, qg_h, qgnew_h, tc_o, tc_h,
     *                   ta_o, ta_h, ts_o, ts_h,
     *                   qg_t, qgnew_t, nrtr
#endif
     *                  )
 
c Controller for hybrid-coordinate semi-Lagrangian transport.
 
c Features:
c * Cubic interpolation for all fields and in all 3 dimensions.
c * Operator-splitting (horizontal then vertical).
c * Fields extended by 2 extra points in all dimensions,
c   (sigma <0, =0, =1, >1 for vertical interpolation).
c * Global mass-conservation adjustment (either Priestley or +/- method
c   for q, done in linems for other tracers)
 
c ug, vg, wg, qg and tracers are on hybrid (eta) coords. The actual
c values of the sigma and hybrid vertical grids are assumed to be
c identical (sigg,sigl,dsig); in fact these names involving "sig" are
c often used for hybrid (not sigma) variables.
 
c The slt subrs proper (slthyb and below) are actually independent
c of whether sigma or hybrid coords are used. The required interpolation
c and setup is done by subrs sigeta, etadot and etapres in sltgcm.f
c (called from gcm subrs before slthyb), as determined by the namelist
c flag hybrid (in comflg) passed to those routines. Subr hybini_[a,d]
c (in sltini.f, called from sltini/linems) initially sets up the hybrid
c transformation arrays in paramhyb used in sigeta, etadot and etapres.
 
c Input:
c wgaus = gcm Gaussian weights, south to north (sum=2)
c sigg  = gcm mid-layer sigma/eta grid, high (0) to low (1)
c sigl  = gcm sigma/eta boundaries, high (0) to low (1)
c dsig  = gcm sigma/eta grid thicknesses
c gravit= gravitational acceleration (m s-2)
c radius= earth radius, m
c dtime = gcm time step (2-dt leapfrog), s
c ug    = gcm eastward  wind velocity, m s-1, at time n
c vg    = gcm northward wind velocity, m s-1, at time n
c wg    = gcm d[sigma/eta]/dt vertical velocity, s-1, at time n
c wgsav = time running-mean of wg, for vertical advect. in stratos (s-1)
c priestq  = .true. for Priestley mass fixer for q, else +/- method
c strattime= exponential timescale for w-filter in extfl (s)
c stratflag= 0 for off, > 0 for on, < 0 for zero w's)
c smolar= .true. for Smolarkiewicz anti-diffusion in vertical
c pthicm= gcm pressure-thicknesses of sig/eta layers,at time n-1 (N m-2)
c pthicp= gcm pressure-thicknesses of sig/eta layers,at time n+1 (N m-2)
c qg    = gcm specific humidity before advection, at time n-1
c tracer= gcm tracer fields at time n-1 (*** modified by advection ***)
c ntrace= number of tracer fields
c ntraca= number of non-passive fields (1-3=cloud liq.h2o,
c         4-7 = flux reservoirs for t,q,u,v)
c lonsav= first-guess departure-point longitudes for this call
c latsav= first-guess departure-point latitudes for this call
c sigsav= first-guess departure-point sigma/eta values for this call
 
c Output:
c qgnew = gcm specific humidity after advection, at time n+1
c tracer= gcm tracer fields (modified by advection), at time n
c lonsav= departure-point longitudes to be saved for next call
c latsav= departure-point latitudes to be saved for next call
c sigsav= departure-point eta values to be saved for next call
c totq[o,n,a] = global mean column mass of water vapor before advection,
c               after advection, and after global adjustment (kg m-2)
 
#include <paramslt.h>
 
      dimension
     *  wgaus(nlat),
     *  sigg(nlev),               sigl(nlev+1),          dsig(nlev),
     *  ug(nlon,nlat,nlev),       vg(nlon,nlat,nlev),
     *  wg(nlon,nlat,nlev),       wgsav(nlon,nlat,nlev),
     *  pthicm(nlon,nlat,nlev),   pthicp(nlon,nlat,nlev),
     *  qg(nlon,nlat,nlev),       qgnew(nlon,nlat,nlev),
     *  tracer(nlon,nlat,nlev,ntrace),
     *  lonsav(nlon,nlat,nlev),   latsav(nlon,nlat,nlev),
     *  sigsav(nlon,nlat,nlev)
#ifdef isotrac
      dimension
     *  qg_o(nlon,nlat,nlev),     qgnew_o(nlon,nlat,nlev),
     *  qg_h(nlon,nlat,nlev),     qgnew_h(nlon,nlat,nlev),
     *  tc_o(nlon,nlat,nlev),     tc_h(nlon,nlat,nlev),
     *  ta_o(nlon,nlat,nlev),     ta_h(nlon,nlat,nlev),
     *  ts_o(nlon,nlat,nlev),     ts_h(nlon,nlat,nlev),
     *  qg_t(nlon,nlat,nlev,nrtr),qgnew_t(nlon,nlat,nlev,nrtr)
#endif
      real lonsav, latsav
 
      integer stratflag
      logical priestq, smolar
 
c Extended-grid variables:
 
c ua             = eastward  velocity
c va             = northward velocity
c wa             = eta-dot vertical velocity
c qa[v]          = specific humidity before advection
c uax,vax,qax    = cubic interpolation factors, used in docube
 
      dimension
     *  ua (-1:nlonp2,-1:nlatp2),
     *  uax(-1:nlonp2,-1:nlatp2,2),
 
     *  va (-1:nlonp2,-1:nlatp2),
     *  vax(-1:nlonp2,-1:nlatp2,2),
 
     *  qa (-1:nlonp2,-1:nlatp2),
     *  qax(-1:nlonp2,-1:nlatp2,2),
 
     *  wa  (nlon,-1:nlevp2),
     *  qav (nlon,-1:nlevp2)

#ifdef isotrac
      dimension
     *  qa_o (-1:nlonp2,-1:nlatp2),   qa_h (-1:nlonp2,-1:nlatp2),
     *  qax_o(-1:nlonp2,-1:nlatp2,2), qax_h(-1:nlonp2,-1:nlatp2,2),
     *  qav_o(   nlon,  -1:nlevp2),   qav_h(   nlon,  -1:nlevp2)
#endif
#ifdef isotracreg
      parameter (nrtrlocal = 1)
      dimension
     *  qa_t (-1:nlonp2,-1:nlatp2,  nrtrlocal),
     *  qax_t(-1:nlonp2,-1:nlatp2,2,nrtrlocal),
     *  qav_t(   nlon,  -1:nlevp2,  nrtrlocal)
c    *  qa_t (-1:nlonp2,-1:nlatp2,  nrtr),
c    *  qax_t(-1:nlonp2,-1:nlatp2,2,nrtr),
c    *  qav_t(   nlon,  -1:nlevp2,  nrtr)
#endif
 
c qglin  = gcm array for linearly interp. water vapor field (Priestley)
c qalinv = extended array for linearly interp. fields (for Priestley)
 
      dimension
     *  qglin  (nlon,nlat,nlev),
     *  qalinv (nlon,-1:nlevp2)
 
c wgfilt = gcm time-filtered vertical velocities to use this timestep
 
      dimension
     *  wgfilt(nlon,nlat,nlev)
 
c Mid-point or departure-point variables, for current longitude strip:
 
c lond  = longitude of mid or departure point
c latd  = latitude  of mid or departure point
c sigd  = sigma     of mid or departure point
c id    = index of closest lona-grid  point to left   of lond
c jd    = index of closest lata-grid  point to south  of latd
c kd    = index of closest siga-grid  point upwards from sigd
c fd    = fraction of lona-grid  interval from lona(id)   to lond
c gd    = fraction of lata-grid  interval from lata(jd)   to latd
c hd    = fraction of siga-grid  interval from siga(kd)   to sigd
c ud    = eastward  velocity
c vd    = northward velocity
c wd    = sigma-dot vertical velocity
 
      dimension
     *  lond(nlon), latd(nlon), sigd(nlon),
     *  id(nlon),   jd(nlon),   kd(nlon),
     *  fd(nlon),   gd(nlon),   hd(nlon),
     *  ud(nlon),   vd(nlon),   wd(nlon)
      real lond, latd
 
c Temporay saves of id,jd,kd,fd,gd,hd for last-mid point of departure-pt
c iteration, used for departure pt for tracers
 
      dimension
     *  idsav(nlon,nlat),    fdsav(nlon,nlat),
     *  jdsav(nlon,nlat),    gdsav(nlon,nlat),
     *  kdsav(nlon,nlev),    hdsav(nlon,nlev)
 
c Temporary array used in counteracting numerical dispersion
 
      dimension
     *  zflux(nlon,nlev+1)
 
 
c-------------------------
c Start of executable code
c-------------------------
 
 
c*********************************
c Horizontal slt for q and tracers
c*********************************
 
c---------------------------------------------------------------------
c Do horizontal departure-point algorithm, q and tracer interpolation,
c one horizontal x,y slice at a time
c---------------------------------------------------------------------
 
#if defined ( omp )
#ifndef isotrac
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared   (ug, vg, qg, lona, lonsav, lata, latsav, iset, jset,
!$OMP*           qgnew, qglin, dtime, radius, ntrace, ntraca, tracer)
!$OMP* private  (i, j, lond, latd, iter, id, fd, jd, gd,
!$OMP*           ua, va, uax, vax, qa, qax, ud, vd, dt,
!$OMP*           idsav, fdsav, jdsav, gdsav, n, k)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared   (ug, vg, qg, lona, lonsav, lata, latsav, iset, jset,
!$OMP*           qgnew, qglin, dtime, radius, ntrace, ntraca, tracer,
!$OMP*           qg_o, qg_h, tc_o, tc_h, ta_o, ta_h, ts_o, ts_h,
!$OMP*           qg_t, nrtr,
!$OMP*           qgnew_o, qgnew_h, qgnew_t)
!$OMP* private  (i, j, lond, latd, iter, id, fd, jd, gd,
!$OMP*           ua, va, uax, vax, qa, qax, ud, vd, dt,
!$OMP*           idsav, fdsav, jdsav, gdsav, n, k,
!$OMP*           qa_o, qax_o, qa_h, qax_h, 
!$OMP*           qa_t, qax_t)
!$OMP* schedule(dynamic,2)
#endif
#endif
 
c---------------------
      do 1000 k=1,nlev
c---------------------
 
c          Set extended fields and e-w cubic interp factors for
c          ua,va,qa for this horizontal x,y slice
 
        call extuv   (ug, vg, ua, va, k)
        call extcube (ua, uax, k)
        call extcube (va, vax, k)
 
        call extfl   (qg, qa, k)
        call extcube (qa, qax, k)
#ifdef isotrac
        call extfl   (qg_o, qa_o, k)
        call extcube (qa_o, qax_o, k)
        call extfl   (qg_h, qa_h, k)
        call extcube (qa_h, qax_h, k)
#endif
#ifdef isotracreg
        do it=1,nrtr
          call extfl   (qg_t(1,1,1,it), qa_t(1,1,it), k)
          call extcube (qa_t(1,1,it),   qax_t(1,1,1,it), k)
        enddo
#endif
 
c-----------------------
        do 1002 j=1,nlat
c-----------------------
 
c            Set first-guess mid-pts using previous departure points
c            lonsav, latsav. The simple averaging here would be
c            inaccurate near the poles, so in polar regions just use
c            the arrival points. The cvmgt test for lond takes care of
c            branch cut at Greenwich, relying on lona and lonsav ge 0
c            and le 2*pi (see vector).
 
          if (abs(lata(j)).le.polelatd*pi/180.) then
            do 100 i=1,nlon
              lond(i) = cvmgt (
     *                      0.5*(lona(i)+lonsav(i,j,k)),
     *                      mod(0.5*(lona(i)+lonsav(i,j,k))+pi,twopi),
     *                      abs(lona(i)-lonsav(i,j,k)).le.pi )
              latd(i) = 0.5 * (lata(j) + latsav(i,j,k))
  100       continue
          else
            do 102 i=1,nlon
              lond(i) = lona(i)
              latd(i) = lata(j)
  102       continue
          endif
 
c            Iterate to find mid-point coords
 
          do 150 iter=1,niter
 
c              Compute grid-indices id,fd,jd,gd of current mid-pt coords
 
            call locate (lond,nlon, lona,nlonp2, nlonset,iset, id,fd)
            call locate (latd,nlon, lata,nlatp2, nlatset,jset, jd,gd)
 
c              Estimate velocities ud,vd at current mid-points
 
            call interph (id, jd, fd, gd, ua, uax, ud, 1)
            call interph (id, jd, fd, gd, va, vax, vd, 0)
 
c              Save final-iteration mid-point grid-indices in global
c              arrays for tracer interpolation below. Use mid-points
c              since tracers are advected over 1*dt timestep, not
c              2*dt (leapfrog) as for q.
 
            if (iter.eq.niter) then
              do 160 i=1,nlon
                idsav(i,j) = id(i)
                fdsav(i,j) = fd(i)
                jdsav(i,j) = jd(i)
                gdsav(i,j) = gd(i)
  160         continue
            endif
 
c              Compute new mid-pts = arrival pts - (ud,vd)*0.5*dtime.
c              If final iteration, use dtime instead of 0.5*dtime.
 
            dt = 0.5*dtime
            if (iter.eq.niter) dt = dtime
            call vectorh (j, ud, vd, dt, radius, lond, latd)
 
  150     continue
 
c            Compute grid-indices id,fd,jd,gd of departure-point coords
 
          call locate (lond,nlon, lona,nlonp2, nlonset,iset, id,fd)
          call locate (latd,nlon, lata,nlatp2, nlatset,jset, jd,gd)
 
c            Interpolate q to departure points
 
          call interph (id, jd, fd, gd, qa, qax, qgnew(1,j,k), 1)
#ifdef isotrac
          call interph (id, jd, fd, gd, qa_o, qax_o, qgnew_o(1,j,k), 0)
          call interph (id, jd, fd, gd, qa_h, qax_h, qgnew_h(1,j,k), 0)
#endif
#ifdef isotracreg
          do it=1,nrtr
            call interph (id, jd, fd, gd, qa_t(1,1,it), qax_t(1,1,1,it),
     *                                            qgnew_t(1,j,k,it), 0)
          enddo
#endif
 
c            *Linearly* interpolate q to departure pts (for Priestley)
 
          call inth_lin (id, jd, fd, gd, qa, qglin(1,j,k))
 
c            Save departure-point coords lonsav,latsav for next call
 
          call scopy (nlon, lond, 1, lonsav(1,j,k), 1)
          call scopy (nlon, latd, 1, latsav(1,j,k), 1)
 
c---------------
 1002   continue
c---------------
 
c          Loop over tracers for horizontal tracer slt. Do outside
c          latitude loop 1002 above so can use local arrays qa,qax for
c          each tracer in turn. Use previously saved mid-pts from last
c          iteration in loop 150 above (since tracer timestep is 1*dt
c          not 2*dt).
 
        do 1500 n=1,ntrace
c         Don't do for flux reservoirs of q,u,v (#s 5,6,7) (see linems):
          if (n.le.4 .or. n.ge.ntraca+1) then
            call extfl (tracer(1,1,1,n), qa, k)
            call extcube (qa, qax, k)
#ifdef isotrac
            if (n.eq.1) then
              call extfl (tc_o, qa_o, k)
              call extcube (qa_o, qax_o, k)
              call extfl (tc_h, qa_h, k)
              call extcube (qa_h, qax_h, k)
            else if (n.eq.2) then
              call extfl (ta_o, qa_o, k)
              call extcube (qa_o, qax_o, k)
              call extfl (ta_h, qa_h, k)
              call extcube (qa_h, qax_h, k)
            else if (n.eq.3) then
              call extfl (ts_o, qa_o, k)
              call extcube (qa_o, qax_o, k)
              call extfl (ts_h, qa_h, k)
              call extcube (qa_h, qax_h, k)
            endif
#endif
            do 1502 j=1,nlat
              call interph  (idsav(1,j), jdsav(1,j),
     *                       fdsav(1,j), gdsav(1,j),
     *                       qa, qax, tracer(1,j,k,n), 1)
#ifdef isotrac
              if (n.eq.1) then
                call interph  (idsav(1,j), jdsav(1,j),
     *                         fdsav(1,j), gdsav(1,j),
     *                         qa_o, qax_o, tc_o(1,j,k), 0)
                call interph  (idsav(1,j), jdsav(1,j),
     *                         fdsav(1,j), gdsav(1,j),
     *                         qa_h, qax_h, tc_h(1,j,k), 0)
              else if (n.eq.2) then
                call interph  (idsav(1,j), jdsav(1,j),
     *                         fdsav(1,j), gdsav(1,j),
     *                         qa_o, qax_o, ta_o(1,j,k), 0)
                call interph  (idsav(1,j), jdsav(1,j),
     *                         fdsav(1,j), gdsav(1,j),
     *                         qa_h, qax_h, ta_h(1,j,k), 0)
              else if (n.eq.3) then
                call interph  (idsav(1,j), jdsav(1,j),
     *                         fdsav(1,j), gdsav(1,j),
     *                         qa_o, qax_o, ts_o(1,j,k), 0)
                call interph  (idsav(1,j), jdsav(1,j),
     *                         fdsav(1,j), gdsav(1,j),
     *                         qa_h, qax_h, ts_h(1,j,k), 0)
              endif
#endif
 1502       continue
          endif
 1500   continue
 
c-------------
 1000 continue
c-------------
 
 
c*******************************
c Vertical slt for q and tracers
c*******************************
 
c        Time filter vertical velocity in stratosphere
 
      call wfilter (wg, wgsav, wgfilt,
     *              sigg, dsig, strattime, stratflag, dtime)
 
c-------------------------------------------------------------------
c Do vertical departure-point algorithm, q and tracer interpolation,
c one longitude-level slice at a time
c-------------------------------------------------------------------
 
#if defined ( omp )
#ifndef isotrac
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared   (siga, sigsav, kset,
!$OMP*           wgfilt, qgnew, qglin,
!$OMP*           dtime, radius, dsig, ntrace, ntraca, tracer)
!$OMP* private  (i, k, sigd, iter, kd, hd,
!$OMP*           wa, qav, qalinv, wd, dt, wfud,
!$OMP*           kdsav, hdsav, n, j)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared   (siga, sigsav, kset,
!$OMP*           wgfilt, qgnew, qglin,
!$OMP*           dtime, radius, dsig, ntrace, ntraca, tracer,
!$OMP*           qgnew_o, qgnew_h, 
!$OMP*           tc_o, tc_h, ta_o, ta_h, ts_o, ts_h,
!$OMP*           qgnew_t, nrtr)
!$OMP* private  (i, k, sigd, iter, kd, hd,
!$OMP*           wa, qav, qalinv, wd, dt, wfud,
!$OMP*           kdsav, hdsav, n, j,
!$OMP*           qav_o, qav_h, 
!$OMP*           qav_t)
!$OMP* schedule(dynamic,2)
#endif
#endif
 
c---------------------
      do 2000 j=1,nlat
c---------------------
 
c          Set extended field wa for this lon-lev slice.
c          Reset the extended field qav (from qgnew this time),
c          and set qalinv for linear interp (for Priestley)
 
        call extflv (wgfilt, wa,     j, 1)
        call extflv (qgnew,  qav,    j, 0)
        call extflv (qglin,  qalinv, j, 0)
#ifdef isotrac
        call extflv (qgnew_o, qav_o, j, 0)
        call extflv (qgnew_h, qav_h, j, 0)
#endif
#ifdef isotracreg
        do it=1,nrtr
          call extflv (qgnew_t(1,1,1,it), qav_t(1,1,it), j, 0)
        enddo
#endif
 
c-----------------------
        do 2002 k=1,nlev
c-----------------------
 
c            Set first-guess mid-pts using previous departure pts sigsav
 
          do 200 i=1,nlon
            sigd(i) = 0.5 * (siga(k) + sigsav(i,j,k))
  200     continue
 
c            Iterate to find mid-point coords
 
          do 250 iter=1,niter
 
c              Compute grid-indices kd,hd of current mid-pt coords
 
            call locate (sigd,nlon, siga,nlevp2, nlevset,kset, kd,hd)
 
c              Estimate velocities wd at current mid-points
 
            call interpv (kd, hd, wa, wd)
 
c8880       Fudge upwards   w for lowest layer to zero,
c           or    downwards w for lowest layer by * 2.
c           (Affects tracers as well as q)
c           if (k.eq.nlev) then
c             do 8880 i=1,nlon
c               if (wd(i).lt.0.) wd(i) = 0.
c               if (wd(i).gt.0.) wd(i) = 2.*wd(i)
c8880         continue
c           endif
 
c              Save final-iteration mid-point grid-indices in global
c              arrays for tracer interpolation below. Use mid-points
c              since tracers are advected over 1*dt timestep, not
c              2*dt (leapfrog) as for q.
 
            if (iter.eq.niter) then
              do 260 i=1,nlon
                kdsav(i,k) = kd(i)
                hdsav(i,k) = hd(i)
  260         continue
            endif
 
c              Compute new mid-pts = arrival pts - (wd)*0.5*dtime.
c              If final iteration, use dtime instead of 0.5*dtime.
 
            dt = 0.5*dtime
            if (iter.eq.niter) dt = dtime
            call vectorv (k, wd, dt, sigd)
 
  250     continue
 
c         Fudge departure points in lowest layer upwards.
c         (Doesn't affect tracers since their kd,hd are already saved)
c         if (k.eq.nlev) then
cc          wfud = 0.06
cc          do 8881 i=1,nlon
cc            sigd(i) = wfud*(1.-dsig(nlev)) + (1.-wfud)*sigd(i)
cc8881      continue
c
c           do 8882 i=1,nlon
c             sigd(i) = sigd(i) - .00195   !               jan = -.05310
c             sigd(i) = sigd(i) - .00180   ! ann = -.3383, jan = -.00124
c             sigd(i) = sigd(i) - .00175   !               jan =  .02158
c             sigd(i) = sigd(i) - .00120   ! ann = -.0662
c             sigd(i) = sigd(i) - .00110   ! ann = -.0048
c             sigd(i) = sigd(i) - .00500   ! for new convection
c8882       continue
c         endif
 
c            Compute grid-indices kd,hd of departure-point coords
 
          call locate (sigd,nlon, siga,nlevp2, nlevset,kset, kd,hd)
 
c            Interpolate q to departure points
 
          call interpv (kd, hd, qav, qgnew(1,j,k))
#ifdef isotrac
          call interpv (kd, hd, qav_o, qgnew_o(1,j,k))
          call interpv (kd, hd, qav_h, qgnew_h(1,j,k))
#endif
#ifdef isotracreg
          do it=1,nrtr
            call interpv (kd, hd, qav_t(1,1,it), qgnew_t(1,j,k,it))
          enddo
#endif

c            *Linearly* interpolate q to departure pts (for Priestley)
 
          call intv_lin (kd, hd, qalinv, qglin(1,j,k))
 
c            Save departure-point coords sigsav for next call
 
          call scopy (nlon, sigd, 1, sigsav(1,j,k), 1)
 
c---------------
 2002   continue
c---------------
 
c          Loop over tracers for vertical tracer slt. Do outside
c          level loop 2002 above so can use local array qav for
c          each tracer in turn. Use previously saved mid-pts from last
c          iteration in loop 250 above (since tracer timestep is 1*dt
c          not 2*dt).
 
        do 2500 n=1,ntrace
c         Don't do for flux reservoirs of q,u,v (#s 5,6,7) (see linems):
          if (n.le.4 .or. n.ge.ntraca+1) then
            call extflv (tracer(1,1,1,n), qav,    j, 0)
#ifdef isotrac
            if (n.eq.1) then
              call extflv (tc_o, qav_o, j, 0)
              call extflv (tc_h, qav_h, j, 0)
            else if (n.eq.2) then
              call extflv (ta_o, qav_o, j, 0)
              call extflv (ta_h, qav_h, j, 0)
            else if (n.eq.3) then
              call extflv (ts_o, qav_o, j, 0)
              call extflv (ts_h, qav_h, j, 0)
            endif
#endif
            do 2502 k=1,nlev
              call interpv  (kdsav(1,k), hdsav(1,k),
     *                       qav, tracer(1,j,k,n))
#ifdef isotrac
              if (n.eq.1) then
                call interpv (kdsav(1,k), hdsav(1,k),
     *                        qav_o, tc_o(1,j,k))
                call interpv (kdsav(1,k), hdsav(1,k),
     *                        qav_h, tc_h(1,j,k))
              else if (n.eq.2) then
                call interpv (kdsav(1,k), hdsav(1,k),
     *                        qav_o, ta_o(1,j,k))
                call interpv (kdsav(1,k), hdsav(1,k),
     *                        qav_h, ta_h(1,j,k))
              else if (n.eq.3) then
                call interpv (kdsav(1,k), hdsav(1,k),
     *                        qav_o, ts_o(1,j,k))
                call interpv (kdsav(1,k), hdsav(1,k),
     *                        qav_h, ts_h(1,j,k))
              endif
#endif
 2502       continue
          endif
 2500   continue
 
c-------------
 2000 continue
c-------------
 
 
c******************************************************************
c Global mass-conservation adjustment for q.
c Choice of two methods...namelist param priestq (default = .true.)
c (Done for tracers in linems via subr traglob)
c******************************************************************
 
      if (priestq) then
 
c          Priestley method
 
        call mascalc (qg,    wgaus, pthicm, gravit, totqo)
        call mascalc (qgnew, wgaus, pthicp, gravit, totqn)
 
        call priestley (qgnew, qglin, totqo, totqn, totqa,
     *                  wgaus, pthicp, gravit)
 
      else
 
c          "+/-" method
 
        call traglob (qgnew, qg, pthicp, pthicm,
     *                totqo, totqn, totqa, wgaus, gravit, 0)
#ifdef isotrac
        call traglob (qgnew_o, qg_o, pthicp, pthicm,
     *                dum1,  dum2,  tot,   wgaus, gravit, 0)
        call traglo_i(qgnew, qg, qgnew_o, qg_o,
     *                pthicp, wgaus, gravit, tot)
        call traglob (qgnew_h, qg_h, pthicp, pthicm,
     *                dum1,  dum2,  tot,   wgaus, gravit, 0)
        call traglo_i(qgnew, qg, qgnew_h, qg_h,
     *                pthicp, wgaus, gravit, tot)
#endif
#ifdef isotracreg
        do it=1,nrtr
          call traglob (qgnew_t(1,1,1,it), qg_t(1,1,1,it),
     *                  pthicp, pthicm,
     *                  dum1,  dum2,  tot,   wgaus, gravit, 0)
          call traglo_i(qgnew, qg, qgnew_t(1,1,1,it), qg_t(1,1,1,it),
     *                  pthicp, wgaus, gravit, tot)
        enddo
#endif
 
      endif
 
c****************************************
c Negative vertical diffusion for tracers
c****************************************
 
c=======================
      do 3000 n=1,ntrace
c=======================
 
c---------------------------------------------------------------------
c Apply some negative vertical diffusion (at all levels) to counteract
c slt vertical numerical dispersion in stratosphere. See also
c running-time filter of vertical velocity applied in extfl.
c Only do for passive tracers ( >= ntraca+1),
c and only if requested by namelist smolar (default = .false.)
c---------------------------------------------------------------------
 
c-------------------------------------------
        if (n.ge.ntraca+1 .and. smolar) then
c-------------------------------------------
 
          do 3300 j=1,nlat
 
c             Compute zflux = downward tracer flux across interfaces
 
            call zero (zflux, nlon*(nlev+1))
            do 3306 k=2,nlev
              zd = sigg(k)-sigg(k-1)
c             cfl constraint:
              zk1 = -0.5 * min(dsig(k),dsig(k-1))*zd / (0.5*dtime)
              do 3308 i=1,nlon
 
c               Either uniform:
c               zk = -30.        ! order of (5 km)**2 / 10 days
c               Or Smolarkiewicz-like:
c               zw = max (abs(wgfilt(i,j,k)), abs(wgfilt(i,j,k-1)))
                zw = abs (0.5*(wgfilt(i,j,k)+wgfilt(i,j,k-1)))
c               Use smaller factor than standard Smol's 0.5 to reduce
c               max concentration values greater than previous max:
c               zk = -0.5 * (zw*zd - 0.5*dtime*(zw**2))
                zk = -0.2 * (zw*zd - 0.5*dtime*(zw**2))
 
c               Constrain so -ve diffusion can't generate -ve tracer:
                zk2 = -0.5 * min (tracer(i,j,k,n)*dsig(k),
     *                            tracer(i,j,k-1,n)*dsig(k-1)) * zd
     *                / ( 0.5*dtime
     *                   * max ( abs(tracer(i,j,k,n)-tracer(i,j,k-1,n)),
     *                            epsilon ) )
 
                zflux(i,k) = max(zk,zk1,zk2)
     *                     * (tracer(i,j,k-1,n)-tracer(i,j,k,n))/zd
 3308         continue
 3306       continue
 
c              Change tracer concentrations
 
            do 3310 k=1,nlev
              do 3312 i=1,nlon
                tracer(i,j,k,n) = tracer(i,j,k,n) + 0.5*dtime
     *                          * (zflux(i,k)-zflux(i,k+1))/dsig(k)
 3312         continue
 3310       continue
 
 3300     continue
 
c------------
        endif
c------------
 
c=============
 3000 continue
c=============
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfilter (wg, wgsav, wgfilt,
     *                    sigg, dsig, strattime, stratflag, dtime)
 
 
c  Apply running time filter to vertical velocity wg, and save in wgsav.
c  This helps to reduce vertical numerical dispersion in the
c  stratosphere due to (presumed non-breaking) waves with timescales
c  <= ~strattime. Do this only above a crude latitudinal-dependent
c  tropopause. See also counteracting of numerical dispersion (smolar)
c  in slthyb.
 
c  nb: must only be called once per timestep for correct wgsav.
 
c Input:
c wg    = gcm d[sigma/eta]/dt vertical velocity, s-1, at time n
c sigg  = gcm mid-layer sigma/eta grid, high (0) to low (1)
c dsig  = gcm sigma/eta grid thicknesses
c strattime = exponential timescale for w-filter in extfl (s)
c stratflag = 0 for off, > 0 for on, < 0 for zero w's)
c dtime = gcm time step (2-dt leapfrog) (s)
 
c Output:
c wgsav = time running-mean of gcm's wg (s-1)
c wgfilt= vertical velocities to use for this timestep (s-1)
 
#include <paramslt.h>
 
      dimension
     *  wg(nlon,nlat,nlev),
     *  wgsav(nlon,nlat,nlev),
     *  wgfilt(nlon,nlat,nlev),
     *  sigg(nlev),
     *  dsig(nlev)
 
      integer stratflag
 
c Local array for rough tropopause height
      dimension sigpaus(nlat)
      logical first
      data first /.true./
      save sigpaus, first
 
 
c        Set tropopause height, first call only
c        (neglecting difference btwn gaussian grid and regular spacing)
 
      if (first) then
        do 50 j=1,nlat
          zcos = cos (-.5*pi + pi*(j-.5)/nlat)
          sigpaus(j) = .100*(zcos**2) + .280*(1.-zcos**2)
   50   continue
        first = .false.
      endif
 
c        Update time-filtered vertical velocity field wgsav
c        and set wgfilt field to be used this timestep
 
      zfil = exp(-0.5*dtime/max(strattime,100.))
      do 100 k=1,nlev
        do 102 j=1,nlat
 
          do 104 i=1,nlon
            wgsav(i,j,k) = zfil*wgsav(i,j,k) + (1.-zfil)*wg(i,j,k)
  104     continue
 
c            Use time-filtered wgsav if in our "strastosphere", and
c            stratflag > 0...or if stratflag < 0, set to zero.
c            Else don't modify w's.
 
          if (sigg(k).le.sigpaus(j) .and. stratflag.ne.0) then
            if (stratflag .gt. 0) then
              do 106 i=1,nlon
                wgfilt(i,j,k) = wgsav(i,j,k)
  106         continue
            else
              do 107 i=1,nlon
                wgfilt(i,j,k) = 0.
  107         continue
            endif
          else
            do 108 i=1,nlon
              wgfilt(i,j,k) =wg(i,j,k)
  108       continue
          endif
 
  102   continue
  100 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine extfl (fieldg, fielda, k)
 
c Sets an extended-grid scalar field fielda from gcm field fieldg,
c for current horizontal x,y slice at one level.
 
c Input:
c fieldg= gcm scalar field
c k     = vertical index
 
c Output:
c fielda = extended scalar field
 
#include <paramslt.h>
 
      dimension
     *  fieldg(nlon,nlat,nlev),
     *  fielda(-1:nlonp2,-1:nlatp2)
 
 
c         Copy gcm fieldg to extended fielda
 
      do 112 j=1,nlat
        do 114 i=1,nlon
          fielda(i,j) = fieldg(i,j,k)
  114   continue
  112 continue
 
c        Extend fielda in latitude, by setting pole values to the
c        zonal mean of the last gcm latitude band, and by setting
c        beyond-pole values to last latitude band + 180 deg longitude.
 
      zavs = ssum(nlon, fielda(1,1), 1) / nlon
      zavn = ssum(nlon, fielda(1,nlat), 1) / nlon
      do 202 i=1,nlon
        fielda(i,0)      = zavs
        fielda(i,nlat+1) = zavn
  202 continue
 
!DIR$ ivdep
      do 204 i=1,nlon/2
        fielda(i,       -1    ) = fielda(i+nlon/2,1   )
        fielda(i+nlon/2,-1    ) = fielda(i,       1   )
        fielda(i,       nlat+2) = fielda(i+nlon/2,nlat)
        fielda(i+nlon/2,nlat+2) = fielda(i,       nlat)
  204 continue
 
c        Extend fielda in longitude by two wraparound points on each end
 
      do 402 j=-1,nlatp2
        fielda(-1,j)     = fielda(nlon-1,j)
        fielda(0,j)      = fielda(nlon,j)
        fielda(nlon+1,j) = fielda(1,j)
        fielda(nlon+2,j) = fielda(2,j)
  402 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine extflv (fieldg, fielda, j, ifw)
 
c Sets an extended-grid scalar field fielda from gcm field fieldg,
c for current longitude-level slice at one latitide.
 
c Input:
c fieldg= gcm scalar field
c j     = latitude index
c ifw   = 0 for scalars, 1 for vertical velocity
 
c Output:
c fielda = extended scalar field
 
#include <paramslt.h>
 
      dimension
     *  fieldg(nlon,nlat,nlev),
     *  fielda(nlon,-1:nlevp2)
 
 
c         Copy gcm fieldg to extended fielda
 
      do 112 k=1,nlev
        do 114 i=1,nlon
          fielda(i,k) = fieldg(i,j,k)
  114   continue
  112 continue
 
c-----------------------
      if (ifw.eq.0) then
c-----------------------
 
c        Extend *scalar* fielda in level. Current method simply uses
c        topmost or bottom-most layer value. The second method
c        (commented out) uses linear vertical extrapolation, which seems
c        to avoid distortions in the cubic vertical interpolations in
c        interpv, and reduces global conservation errors (in combination
c        with setting of siga - see extgrd in inislt.f), at least for q.
c        HOWEVER, the second method can cause spurious increase of
c        tracer amounts if there is a maximum in the upper layer and
c        upper-layer sigma-dots are downward.
 
c       second method:
c       xa = (siga(-1)-siga(1)) / (siga(2)-siga(1))
c       xb = (siga(0) -siga(1)) / (siga(2)-siga(1))
c       xc = (siga(nlev+1)-siga(nlevm)) / (siga(nlev)-siga(nlevm))
c       xd = (siga(nlev+2)-siga(nlevm)) / (siga(nlev)-siga(nlevm))
 
        do 302 i=1,nlon
          fielda(i,-1)    = fielda(i,1)
          fielda(i,0)     = fielda(i,1)
          fielda(i,nlev+1) = fielda(i,nlev)
          fielda(i,nlev+2) = fielda(i,nlev)
 
c         second method:
c         fielda(i,-1)    = (1.-xa)*fielda(i,1) + xa*fielda(i,2)
c         fielda(i,0)     = (1.-xb)*fielda(i,1) + xb*fielda(i,2)
c         fielda(i,nlev+1) =
c    *                   (1.-xc)*fielda(i,nlevm) + xc*fielda(i,nlev)
c         fielda(i,nlev+2) =
c    *                   (1.-xd)*fielda(i,nlevm) + xd*fielda(i,nlev)
 
c         Constrain extrapolated values ge 0 (q,tracers):
          fielda(i,-1)     = max (0., fielda(i,-1))
          fielda(i,0)      = max (0., fielda(i,0))
          fielda(i,nlev+1) = max (0., fielda(i,nlev+1))
          fielda(i,nlev+2) = max (0., fielda(i,nlev+2))
  302   continue
 
c---------
      else
c---------
 
c         Extend vertical velocity in level, by setting to zero at
c         boundaries (0,1) and linear vertical extrapolation outside.
 
        xa = (siga(0)-siga(-1)) / (siga(1)-siga(0))
        xd = (siga(nlev+2)-siga(nlev+1)) / (siga(nlev+1)-siga(nlev))
        do 322 i=1,nlon
          fielda(i,-1)     = -xa*fielda(i,1)
          fielda(i,0)      = 0.
          fielda(i,nlev+1) = 0.
          fielda(i,nlev+2) = -xd*fielda(i,nlev)
  322   continue
 
c----------
      endif
c----------
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ext1d (fieldg, fielda)
 
c Sets a 1-D (horiz) extended-grid scalar field fielda from gcm field
c fieldg (used for surface pressure).
 
c Input:
c fieldg = gcm scalar field
 
c Output:
c fielda = extended scalar field
 
#include <paramslt.h>
 
      dimension
     *  fieldg(nlon,nlat),   fielda(-1:nlonp2,-1:nlatp2)
 
 
c        Copy extended fielda from gcm fieldg
 
      do 102 j=1,nlat
        do 104 i=1,nlon
          fielda(i,j) = fieldg(i,j)
  104   continue
  102 continue
 
c        Extend fielda in latitude, by setting pole values to the
c        zonal mean of the last gcm latitude band, and by setting
c        beyond-pole values to last latitude band + 180 deg longitude.
 
      zavs = ssum(nlon, fielda(1,1), 1) / nlon
      zavn = ssum(nlon, fielda(1,nlat), 1) / nlon
      do 202 i=1,nlon
        fielda(i,0)      = zavs
        fielda(i,nlat+1) = zavn
  202 continue
 
!DIR$ ivdep
      do 204 i=1,nlon/2
        fielda(i,       -1)     = fielda(i+nlon/2,1)
        fielda(i+nlon/2,-1)     = fielda(i,       1)
        fielda(i,       nlat+2) = fielda(i+nlon/2,nlat)
        fielda(i+nlon/2,nlat+2) = fielda(i,       nlat)
  204 continue
 
c        Extend fielda in longitude by two wraparound points on each end
 
      do 402 j=-1,nlatp2
        fielda(-1,j)     = fielda(nlon-1,j)
        fielda(0,j)      = fielda(nlon,j)
        fielda(nlon+1,j) = fielda(1,j)
        fielda(nlon+2,j) = fielda(2,j)
  402 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine extuv (ug, vg, ua, va, k)
 
c Sets extended-grid horizontal velocities ua,va from gcm ug,vg,
c for current horizontal x,y slice at one level.
 
c Input:
c ug,vg  = gcm eastward,northward velocities
c k      = vertical index
 
c Output:
c ua,va  = extended eastward,northward velocities
 
#include <paramslt.h>
 
      dimension
     *  ug(nlon,nlat,nlev),       vg(nlon,nlat,nlev)
 
      dimension
     *  ua(-1:nlonp2,-1:nlatp2),
     *  va(-1:nlonp2,-1:nlatp2)
 
      dimension
     *  utmp(nlon), vtmp(nlon)
 
 
c        Copy extended ua,va from ug,vg
 
      do 102 j=1,nlat
        do 104 i=1,nlon
          ua(i,j) = ug(i,j,k)
          va(i,j) = vg(i,j,k)
  104   continue
  102 continue
 
c        Extend ua,va in latitude. For pole, fit polemost gcm latitude
c        circle to pure slab motion, and apply that slab motion at the
c        pole (note sign changes for va between poles). For beyond-pole,
c        use -1 * polemost latitude band + 180 shift in longitude.
 
c        south pole
 
      do 202 i=1,nlon
        utmp(i) = -ua(i,1)*slona(i) + va(i,1)*clona(i)
        vtmp(i) =  ua(i,1)*clona(i) + va(i,1)*slona(i)
  202 continue
      uslab = ssum(nlon, utmp, 1) / nlon
      vslab = ssum(nlon, vtmp, 1) / nlon
      do 204 i=1,nlon
        ua(i,0) = -uslab*slona(i) + vslab*clona(i)
        va(i,0) =  uslab*clona(i) + vslab*slona(i)
  204 continue
 
c        north pole
 
      do 206 i=1,nlon
        utmp(i) = -ua(i,nlat)*slona(i) - va(i,nlat)*clona(i)
        vtmp(i) =  ua(i,nlat)*clona(i) - va(i,nlat)*slona(i)
  206 continue
      uslab = ssum(nlon, utmp, 1) / nlon
      vslab = ssum(nlon, vtmp, 1) / nlon
      do 208 i=1,nlon
        ua(i,nlat+1) = -uslab*slona(i) + vslab*clona(i)
        va(i,nlat+1) = -uslab*clona(i) - vslab*slona(i)
  208 continue
 
c        beyond both poles
 
!DIR$ ivdep
      do 210 i=1,nlon/2
        ua(i,       -1) = -ua(i+nlon/2,1)
        va(i,       -1) = -va(i+nlon/2,1)
        ua(i+nlon/2,-1) = -ua(i,       1)
        va(i+nlon/2,-1) = -va(i,       1)
 
        ua(i,       nlat+2) = -ua(i+nlon/2,nlat)
        va(i,       nlat+2) = -va(i+nlon/2,nlat)
        ua(i+nlon/2,nlat+2) = -ua(i,       nlat)
        va(i+nlon/2,nlat+2) = -va(i,       nlat)
  210 continue
 
c        Extend ua,va in longitude by two wraparound points on each end
 
      do 402 j=-1,nlatp2
        ua(-1,j)     = ua(nlon-1,j)
        ua(0,j)      = ua(nlon,j)
        ua(nlon+1,j) = ua(1,j)
        ua(nlon+2,j) = ua(2,j)
 
        va(-1,j)     = va(nlon-1,j)
        va(0,j)      = va(nlon,j)
        va(nlon+1,j) = va(1,j)
        va(nlon+2,j) = va(2,j)
  402 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine extcube (vala, valax, k)
 
c Sets e-w cubic interpolation factors valax that depend on the field,
c for extended-grid field vala.
 
c Input:
c vala = extended field
c k    = vertical index
 
c Output:
c valax = e-w cubic interpolation factors for each extended grid point,
c         for a horizontal x,y slice at one level,
c         to be used in interph for pts in the x-interval to the right
c         of each x-grid pt. valax(*,*,1) involves the value and
c         derivatives at the right-hand side of the interval (i+1), and
c         valax(*,*,2) involves those values at the left-hand side(i).
c         This pre-computation is feasible since x-computations are
c         done first so always involve grid values, and is worthwhile
c         since 4 times as many e-w interps are done than n-s interps
c         (see interph).
 
#include <paramslt.h>
 
      dimension
     *  vala (-1:nlonp2,-1:nlatp2),
     *  valax(-1:nlonp2,-1:nlatp2,2)
 
      dimension
     *  dm1(nlon), d(nlon), dp1(nlon),
     *  der(nlon), derp(nlon), dermin(nlon), dermax(nlon)
 
 
      do 102 j=-1,nlatp2
        do 104 i=1,nlon
 
c            Calculate the Lagrange derivative estimates der and derp
 
          dm1(i) = (vala(i,  j)-vala(i-1,j)) / delt_x(i-1)
          d(i)   = (vala(i+1,j)-vala(i,  j)) / delt_x(i)
          dp1(i) = (vala(i+2,j)-vala(i+1,j)) / delt_x(i+1)
 
          der(i)  = dm1(i)
     *            + faca_x(i)*(d(i)-dm1(i)) + facap_x(i)*(dp1(i)-d(i))
          derp(i) = dm1(i)
     *            + facb_x(i)*(d(i)-dm1(i)) + facbp_x(i)*(dp1(i)-d(i))
 
c            Apply SCM0 (sufficient condition for monotonicity)
 
          dermin(i) = cvmgt (0., 3*d(i), d(i).gt.0.)
          dermax(i) = cvmgt (3*d(i), 0., d(i).gt.0.)
          der(i)  = max ( dermin(i), min (dermax(i),der(i)) )
          derp(i) = max ( dermin(i), min (dermax(i),derp(i)) )
 
c            Estimate Hermite-cubic interpolation factors
 
          valax(i,j,1) = 3.*vala(i+1,j) - delt_x(i)*derp(i)
          valax(i,j,2) = 3.*vala(i,  j) + delt_x(i)*der(i)
 
  104   continue
  102 continue
 
c       Set wraparound points for valax (but -1 and nlon+2 never used)
 
      do 202 j=-1,nlatp2
        valax(-1,j,1)     = valax(nlon-1,j,1)
        valax(-1,j,2)     = valax(nlon-1,j,2)
        valax( 0,j,1)     = valax(nlon,  j,1)
        valax( 0,j,2)     = valax(nlon,  j,2)
        valax(nlon+1,j,1) = valax(1,j,1)
        valax(nlon+1,j,2) = valax(1,j,2)
        valax(nlon+2,j,1) = valax(2,j,1)
        valax(nlon+2,j,2) = valax(2,j,2)
  202 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine locate (targ,nlon, grid, ngrid, nset, iset, id, fd)
 
c Finds indices id, fd into an extended grid for a longitude strip of
c values targ. This is done by finding the closest grid points in a
c finely-spaced uniform grid of nset values between the extremes of
c the extended grid, and using the uniform grid's indices iset computed
c previously in locset in inislt.f.
 
c Input:
c targ = target coords (x, y or sigma) for one longitude strip
c nlon = dimension of targ (and id,fd)
c grid = extended grid (first and last elements = extremes of unif grid)
c ngrid= dimension of extended grid (-1:ngrid)
c nset = number of uniform grid points
c iset = index of closest extended grid pt to left of each unif grid pt
 
c Output:
c id   = copy of iset for each element of targ
c fd   = fraction of extended grid interval from left edge to target
 
      dimension
     *  targ(nlon), grid(-1:ngrid),
     *  iset(nset), id(nlon), fd(nlon)
 
c Set extremes of uniform grid
 
      seta = grid(-1)
      setb = grid(ngrid)
 
c Locate targ in the uniform grid, copy its iset to id, and compute fd
 
      do 100 i=1,nlon
        m = nset*(targ(i)-seta)/(setb-seta) + 1.001
        m = max (1, min(nset,m))
        zv = seta + (setb-seta)*(m-.5)/nset
        id(i) = iset(m)
        fd(i) = (zv-grid(id(i))) / (grid(id(i)+1)-grid(id(i)))
  100 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine interph (id, jd, fd, gd, fielda, fieldax, fieldd,iflag)
 
c Estimates field values fieldd at a longitude strip of mid or departure
c points, by *horizontal* interpolation from extended-grid field fielda.
c Calls subr docube to actually do 1-D interpolations, first in e-w
c then in n-s directions.
 
c Input:
c id     = index of closest lona-grid point to left  of mid/dep pt
c jd     = index of closest lata-grid point to south of mid/dep pt
c fd     = fraction of lona-grid interval from grid pt to mid/dep pt
c gd     = fraction of lata-grid interval from grid pt to mid/dep pt
c fielda = extended field to be interpolated
c fieldax= cubic interpolation factors for e-w direction
c iflag  = 1 to compute fdpow,gdpow, else use previous values
 
c Output:
c fieldd = interpolated values for current strip of mid/dep points
 
#include <paramslt.h>
 
      dimension
     *  id(nlon),          jd(nlon),
     *  fd(nlon),          gd(nlon),
     *  fielda  (-1:nlonp2,-1:nlatp2),
     *  fieldax (-1:nlonp2,-1:nlatp2, 2),
     *  fieldd(nlon)

c     Task common block for parallelization (local to interph)
      COMMON  /cinterph/
     *  fdpow(nlon,4),     gdpow(nlon,4)
!$OMP THREADPRIVATE ( /cinterph/ )
#ifndef omp
      save  /cinterph/
#endif

      dimension
     *  valy(nlon,-1:2)
 
c        If iflag = 1, pre-compute cubic interpolation terms that depend
c        just on the fractional location within the grid square. If
c        iflag ne 1, use values from previous call to save cpu time
c        (done with task common for multiprocessing)
 
      if (iflag.eq.1) then
        do 100 i=1,nlon
          fdpow(i,1) = fd(i)**3
          fdpow(i,2) = fd(i)**2 * (1.-fd(i))
          fdpow(i,3) = fd(i)    * (1.-fd(i))**2
          fdpow(i,4) =            (1.-fd(i))**3
 
          gdpow(i,1) = gd(i)**3
          gdpow(i,2) = gd(i)**2 * (1.-gd(i))
          gdpow(i,3) = gd(i)    * (1.-gd(i))**2
          gdpow(i,4) =            (1.-gd(i))**3
  100   continue
      endif
 
c        Do cubic interpolation, first in longitude to get valy, then
c        in latitude to get fieldd. The second subscript of valy runs
c        from -1 to 2 for the 4 points on the y grid surrounding the
c        location being interpolated to.
 
      do 200 jsh = -1,2
 
        do 300 i=1,nlon
          valy(i,jsh) = fielda  (id(i)+1, jd(i)+jsh)    * fdpow(i,1)
     *                + fieldax (id(i),   jd(i)+jsh, 1) * fdpow(i,2)
     *                + fieldax (id(i),   jd(i)+jsh, 2) * fdpow(i,3)
     *                + fielda  (id(i),   jd(i)+jsh)    * fdpow(i,4)
  300   continue
 
  200 continue
 
      call docube (jd, gdpow,
     *             delt_y, faca_y, facap_y, facb_y, facbp_y,
     *             valy, fieldd)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine inth_lin (id, jd, fd, gd, fielda, fieldd)
 
c Estimates field values fieldd at a longitude strip of mid or departure
c pts, by **linear** *horizontal* interpolation from extended-grid field
c fielda.
 
c Input:
c id     = index of closest lona-grid point to left  of mid/dep pt
c jd     = index of closest lata-grid point to south of mid/dep pt
c fd     = fraction of lona-grid interval from grid pt to mid/dep pt
c gd     = fraction of lata-grid interval from grid pt to mid/dep pt
c fielda = extended-grid field to be interpolated
 
c Output:
c fieldd = interpolated values for current strip of mid/dep points
 
#include <paramslt.h>
 
      dimension
     *  id(nlon),          jd(nlon),
     *  fd(nlon),          gd(nlon),
     *  fielda (-1:nlonp2,-1:nlatp2),
     *  fieldd(nlon)
 
 
      do 10 i=1,nlon
        fieldd(i) = (1.-fd(i))*(1.-gd(i))*fielda(id(i)  ,jd(i)  )
     *            + (   fd(i))*(1.-gd(i))*fielda(id(i)+1,jd(i)  )
     *            + (1.-fd(i))*(   gd(i))*fielda(id(i)  ,jd(i)+1)
     *            + (   fd(i))*(   gd(i))*fielda(id(i)+1,jd(i)+1)
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine interpv (kd, hd, fielda, fieldd)
 
c Estimates field values fieldd at a longitude strip of mid or departure
c points, by *vertical* interpolation from extended-grid field fielda.
c Calls subr docube to actually do 1-D interpolation in the vertical.
 
c Input:
c kd     = index of closest siga-grid point upwards from mid/dep pt
c hd     = fraction of siga-grid interval from grid pt to mid/dep pt
c fielda = extended field to be interpolated
 
c Output:
c fieldd = interpolated values for current strip of mid/dep points
 
#include <paramslt.h>
 
      dimension
     *  kd(nlon),  hd(nlon),
     *  fielda (nlon,-1:nlevp2),
     *  fieldd (nlon)
 
      dimension
     *  hdpow(nlon,4),
     *  valz(nlon,-1:2)
 
c        Compute cubic interpolation terms that depend just on the
c        fractional location within the grid square.
 
      do 100 i=1,nlon
        hdpow(i,1) = hd(i)**3
        hdpow(i,2) = hd(i)**2 * (1.-hd(i))
        hdpow(i,3) = hd(i)    * (1.-hd(i))**2
        hdpow(i,4) =            (1.-hd(i))**3
  100 continue
 
c        Assemble valz array for current strip as needed by docube.
c        (Second dimension of valz holds the 4 surrounding vertical pts)
 
      do 200 ksh = -1,2
        do 202 i=1,nlon
          valz(i,ksh) = fielda (i, kd(i)+ksh)
  202   continue
  200 continue
 
      call docube (kd, hdpow,
     *             delt_z, faca_z, facap_z, facb_z, facbp_z,
     *             valz, fieldd)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine intv_lin (kd, hd, fielda, fieldd)
 
c Estimates field values fieldd at a longitude strip of mid or departure
c pts, by **linear** *vertical* interpolation from extended-grid field
c fielda.
 
c Input:
c kd     = index of closest siga-grid point upwards from mid/dep pt
c hd     = fraction of siga-grid interval from grid pt to mid/dep pt
c fielda = extended-grid field to be interpolated
 
c Output:
c fieldd = interpolated values for current strip of mid/dep points
 
#include <paramslt.h>
 
      dimension
     *  kd(nlon),  hd(nlon),
     *  fielda (nlon,-1:nlevp2),
     *  fieldd(nlon)
 
 
      do 10 i=1,nlon
        fieldd(i) = (1.-hd(i))* fielda(i,kd(i)  )
     *            +     hd(i) * fielda(i,kd(i)+1)
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine docube (id, fdpow,
     *                   delt, faca, facap, facb, facbp,
     *                   valin, valout)
 
c Does interpolation in one dimension(*) for a longitude strip of mid or
c departure points. 4 surrounding values on the grid in this dimension
c are in valin(*,-1:2), and result is returned in valout(*).
c Interpolation is done by Lagrange/SCM0/Hermite cubics, following
c Williamson and Rasch (1989), MWR, 117, 102-129.
 
c Input:
c id = index of closest ext. grid pt to left/south/above of mid/dep pt.
c fdpow = factors that depend only on fraction of grid interval
c delt = extended grid intervals (dx, dy or dz) (pre-computed for speed)
c fac[a,ap,b,bp] = factors for Lagrangian deriv estimates, that depend
c                  only on the grid (pre-computed for speed)
c valin = field values at 4 surrounding points for each mid/dep point.
 
c Output:
c valout = interpolated values for the current strip of mid/dep points.
 
#include <paramslt.h>
 
c     parameter (ndim=max(nlon,nlat,nlev))
      parameter (ndim=nlon+nlat+nlev)
 
      dimension
     *  id (nlon), fdpow(nlon,4),
     *  delt(-1:ndim+1),
     *  faca(0:ndim), facap(0:ndim), facb(0:ndim), facbp(0:ndim),
     *  valin(nlon,-1:2), valout(nlon)
 
      dimension
     *  dm1(nlon), d(nlon), dp1(nlon),
     *  der(nlon), derp(nlon), dermin(nlon), dermax(nlon)
 
 
      do 100 i=1,nlon
 
c          Calculate the Lagrange derivative estimates der and derp
 
        k = id(i)
 
        dm1(i) = (valin(i,0)-valin(i,-1)) / delt(k-1)
        d(i)   = (valin(i,1)-valin(i,0) ) / delt(k)
        dp1(i) = (valin(i,2)-valin(i,1) ) / delt(k+1)
 
        der(i)  = dm1(i)
     *          + faca(k) * (d(i)-dm1(i)) + facap(k) * (dp1(i)-d(i))
        derp(i) = dm1(i)
     *          + facb(k) * (d(i)-dm1(i)) + facbp(k) * (dp1(i)-d(i))
 
c          Apply SCM0 (sufficient condition for monotonicity)
 
cold    dermin(i) = cvmgt (0., 3*d(i), d(i).gt.0.)
cold    dermax(i) = cvmgt (3*d(i), 0., d(i).gt.0.)
        dermin(i) = min (3.*d(i), 0.)
        dermax(i) = max (3.*d(i), 0.)
        der(i)  = max ( dermin(i), min (dermax(i),der(i)) )
        derp(i) = max ( dermin(i), min (dermax(i),derp(i)) )
 
c          Estimate final values by Hermite cubic interpolation
 
        valout(i)= valin(i,1)                      * fdpow(i,1)
     *           + (3.*valin(i,1)-delt(k)*derp(i)) * fdpow(i,2)
     *           + (3.*valin(i,0)+delt(k)*der(i) ) * fdpow(i,3)
     *           + valin(i,0)                      * fdpow(i,4)
 
  100 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vectorh (j, ud, vd, dt, radius, lond, latd)
 
c Computes a new longitude-strip set of mid or departure points by
c subtracting the horizontal vector (ud,vd)*dt from the arrival points.
 
c Equatorward of +/- polelatd degrees, horizontal motion is done simply
c by a linearized spherical increment. Poleward of polelatd, horizontal
c motion is done by transforming to an approximately cartesian grid
c defined by:
c   x = earth radius * (pi/2-latitude) * cos(longitude)
c   y = earth radius * (pi/2-latitude) * sin(longitude)
c (or ...latitude+pi/2... for southern hemisphere). The resulting
c longitudes can cross Greenwich, but are reset to the range 0 to 2*pi.
 
c Input:
c j      = gcm latitude index of current longitude strip
c ud     = mid or departure-point eastward  velocities
c vd     = mid or departure-point northward velocities
c dt     = time step
c radius = earth radius
 
c nb: this strip's arrival points  = (lona(i),lata(j), i=1,nlon)
 
c Output:
c lond   = longitudes of new mid or departure points
c latd   = latitudes  of new mid or departure points
 
#include <paramslt.h>
 
      dimension lond(nlon), latd(nlon), ud(nlon), vd(nlon)
      real lond, latd
 
      dimension ucart(nlon), vcart(nlon), xcart(nlon), ycart(nlon)
 
 
      polelat = polelatd*pi/180.
 
      If (lata(j).lt.-polelat) then
 
c            Horizontal motion, near south pole
 
        rcart = radius*(lata(j) + pio2)
 
        do 100 i=1,nlon
          ucart(i) = -ud(i)*slona(i) + vd(i)*clona(i)
          vcart(i) =  ud(i)*clona(i) + vd(i)*slona(i)
 
          xcart(i) = rcart*clona(i) - ucart(i)*dt
          ycart(i) = rcart*slona(i) - vcart(i)*dt
 
          latd(i) = -pio2 + sqrt(xcart(i)**2 + ycart(i)**2) / radius
 
          xcart(i) = cvmgt (xcart(i), 1.e-20, abs(xcart(i)).gt.1.e-20)
          lond(i) = atan2 (ycart(i), xcart(i))
          lond(i) = cvmgt (lond(i), lond(i) + twopi, lond(i).ge.0.)
          lond(i) = mod (lond(i),twopi)
  100   continue
 
      else if (lata(j).gt.polelat) then
 
c            Horizontal motion, near north pole.
c            Note sign changes for vd between poles.
 
        rcart = radius*(pio2-lata(j))
 
        do 200 i=1,nlon
          ucart(i) = -ud(i)*slona(i) - vd(i)*clona(i)
          vcart(i) =  ud(i)*clona(i) - vd(i)*slona(i)
 
          xcart(i) = rcart*clona(i) - ucart(i)*dt
          ycart(i) = rcart*slona(i) - vcart(i)*dt
 
          latd(i) = pio2 - sqrt(xcart(i)**2 + ycart(i)**2) / radius
 
          xcart(i) = cvmgt (xcart(i), 1.e-20, abs(xcart(i)).gt.1.e-20)
          lond(i) = atan2 (ycart(i), xcart(i))
          lond(i) = cvmgt (lond(i), lond(i) + twopi, lond(i).ge.0.)
          lond(i) = mod (lond(i),twopi)
  200   continue
 
      else
 
c            Horizontal motion, not near poles
 
        do 300 i=1,nlon
          lond(i) = lona(i) - (ud(i)/(radius*clata(j)))*dt
          lond(i) = cvmgt (lond(i), lond(i) + twopi, lond(i).ge.0.)
          lond(i) = mod (lond(i),twopi)
 
          latd(i) = lata(j) - (vd(i)/radius)*dt
  300   continue
 
      endif
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vectorv (k, wd, dt, sigd)
 
c Computes a new longitude-strip set of mid or departure points by
c subtracting the vector wd*dt from the arrival points.
 
c Vertical motion is done everywhere by a simple increment of sigma,
c with the result constrained to be between sigma=0 and sigma=1.
 
c Input:
c k      = gcm vertical index of current longitude strip
c wd     = mid or departure-point vertical sigma-dot velocities
c dt     = time step
 
c nb: this strip's arrival points  = (siga(k),i=1,nlon)
 
c Output:
c sigd   = sigma coords of new mid or departure points
 
#include <paramslt.h>
 
      dimension wd(nlon), sigd(nlon)
 
 
      do 400 i=1,nlon
        sigd(i) = siga(k) - wd(i)*dt
        sigd(i) = max ( 0., min(1.,sigd(i)) )
  400 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mascalc (fieldg, wgaus, pthic, gravit, tot)
 
c Calculates global mass of fieldg.
 
c Input:
c fieldg= gcm field (kg/kg)
c wgaus = gcm Gaussian weights, south to north (sum=2)
c pthic = gcm pressure-thicknesses of (hybrid) eta-grid boxes
c gravit= gravitational acceleration (m s-2)
 
c Output:
c tot  = global mean column mass of fieldg (kg m-2)
 
#include <paramslt.h>
 
      dimension
     *  fieldg(nlon,nlat,nlev),          wgaus(nlat),
     *  pthic(nlon,nlat,nlev)
 
 
      tot = 0.
 
      do 100 k=1,nlev
        do 102 j=1,nlat
          do 104 i=1,nlon
            tot = tot + fieldg(i,j,k)*wgaus(j)*pthic(i,j,k)
  104     continue
  102   continue
  100 continue
 
c        Normalize (2 for wgaus, nlon for longitudes, gravit for pthic)
 
      tot = tot / (2.*nlon*gravit)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine priestley (fieldg, fieldl, toto, totn, tota,
     *                      wgaus, pthic, gravit)
 
c Does global adjustment of fieldg to conserve total mass (restores
c totn to toto), similar to Priestley (1993, MWR, 121, 621-629).
 
c Input:
c fieldg = gcm field after slt advection (kg/kg)
c fieldl = gcm field after *linear* interpolation (kg/kg)
c toto   = global mean column mass before advection (kg m-2)
c totn   = global mean column mass after advection (kg m-2)
c wgaus  = gcm Gaussian weights, south to north (sum=2)
c pthic  = gcm pressure-thicknesses of (hybrid) eta-grid boxes
c gravit = gravitational acceleration (m s-2)
 
c Output:
c fieldg = gcm *adjusted* field (kg/kg)
c tota   = global mean column mass after global adjustment (kg m-2)
 
#include <paramslt.h>
 
      dimension
     *  fieldg(nlon,nlat,nlev),        fieldl(nlon,nlat,nlev),
     *  wgaus(nlat),                   pthic(nlon,nlat,nlev)
 
 
      surplus = totn - toto
      dir = sign (1., surplus)
 
c        Calculate maximum possible change by shifting all fieldg's
c        to fieldl's wherever this shift has the right sign
 
      avail = 0.
      do 100 k=1,nlev
        do 102 j=1,nlat
          do 104 i=1,nlon
            avail = avail
     *            + max (0., dir*(fieldg(i,j,k)-fieldl(i,j,k)))
     *            * wgaus(j) * pthic(i,j,k)
  104     continue
  102   continue
  100 continue
      avail = avail / (2.*nlon*gravit)
 
c        Set fraction of maximum change needed to restore global mass
c        (if would be > 1, cannot fully restore keeping within fieldl)
 
      alpha = min ( 1., abs(surplus) / max(avail,epsilon) )
 
c        Shift fieldg to restore global mass
 
      do 200 k=1,nlev
        do 202 j=1,nlat
          do 204 i=1,nlon
            fieldg(i,j,k) = fieldg(i,j,k) - alpha*dir
     *                    * max (0., dir*(fieldg(i,j,k)-fieldl(i,j,k)))
  204     continue
  202   continue
  200 continue
 
c        Set adjusted global mass (diagnostic)
 
      beta = min ( 1., avail / max(abs(surplus),epsilon) )
      tota = totn + beta*(toto-totn)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine traglob (trnew, trold, pthic, pthico,
     *                    toto,  totn,  tota, wgaus, gravit, n)
 
c "+/-" method: global mass conservation adjustment for one tracer field
 
c Modified:
c trnew      = gcm tracer field after advection
c tot[o,n,a] = global mean column mass of tracers before advection,
c              after advection, and after global adjustment (kg m-2)
 
c Supplied:
c trold  = gcm tracer field before advection
c pthic  = gcm pressure-thicknesses of sig/eta layers, current (N m-2)
c pthico = gcm pressure-thicknesses before advection (N m-2)
c wgaus  = gcm Gaussian weights, south to north (sum=2)
c gravit = gravitational acceleration (m s-2)
c n      = tracer number (or zero for q), for diagnostics only
 
#include <paramslt.h>
 
      dimension
     *  trnew(nlon,nlat,nlev),             trold(nlon,nlat,nlev),
     *  pthic(nlon,nlat,nlev),             pthico(nlon,nlat,nlev),
     *  wgaus(nlat)
 
 
c        Calculate total mass of tracer before and after advection
 
      call mascalc (trold, wgaus, pthico, gravit, toto)
      call mascalc (trnew, wgaus, pthic,  gravit, totn)
 
c        Use a multiplicative method proportional to the local
c        *change* in the field, one value for +ve changes, another
c        for -ve changes.
 
c        Calculate global totals where change in field is
c        non zero, and weighted change is +ve (zplus) and -ve (zminus)
 
      zsame  = 0.
      zplus  = 0.
      zminus = 0.
      do k=1,nlev
        do j=1,nlat
          zwei = wgaus(j)/(2.*nlon*gravit)
          do i=1,nlon
 
c           First, if change is due to roundoff error in slt,
c           set back to exact pre-advection value
            zz = trnew(i,j,k) - trold(i,j,k)
            if ( abs(zz)/max(abs(trold(i,j,k)),1.e-20) .lt. 1.e-10) then
              trnew(i,j,k)= trold(i,j,k)
              zz = 0.
            endif
 
            zinc = trnew(i,j,k)*pthic (i,j,k)
     *           - trold(i,j,k)*pthico(i,j,k)
 
            if (zz.eq.0.) then
              zsame = zsame + zinc*zwei
            else
              if (zinc.gt.0.) then
                zplus  = zplus  + zinc*zwei
              else
              zminus = zminus - zinc*zwei
              endif
            endif
 
          enddo
        enddo
      enddo
 
c     write (*,500) n, zsame, zplus, zminus
c 500 format('traglob: n=',i2,'  zsame,zplus,zminus=',3f25.12)
 
      if (zplus.gt.1.e-20 .and. zminus.gt.1.e-20) then
        zfacp = 0.5*(zplus+zminus-zsame) / zplus
        zfacm = 0.5*(zplus+zminus+zsame) / zminus
        do k=1,nlev
          do j=1,nlat
            do i=1,nlon
              zz = trnew(i,j,k) - trold(i,j,k)
              if (abs(zz).gt.0.) then
                zinc =   trnew(i,j,k)*pthic (i,j,k)
     *                 - trold(i,j,k)*pthico(i,j,k)
                if (zinc.gt.0.) then
                  zfac = zfacp
                else
                  zfac = zfacm
                endif
                trnew(i,j,k)= (   trold(i,j,k)*pthico(i,j,k)
     *                          + zinc*zfac )
     *                        / pthic(i,j,k)
              endif
 
            enddo
          enddo
        enddo
      endif
 
c        Calculate total mass of tracer after adjustment
 
      call mascalc (trnew, wgaus, pthic, gravit, tota)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sltoff (nstep,  dtime,  mcdate, calday, hybtype,
     *                   along,  alatg,  wgaus,  sigg,   sigl,   dsig,
     *                   gravit, radius, twodt,
     *                   ug,     vg,     wg,     wgsav,
     *                   psurf,  pmid,   pint,   pthicm, pthic,  pthicp,
     *                   albsud, albsui, solsud, solsui, cosq,
     *                   tg,     qg,     qgnew,
     *                   cloudc, clouda, clouds)
 
c Writes variables needed for slt to binary files at each time step,
c to be read by an off-line slt driver. Each file contains one day's
c data, with names sltoff.<mmdd> where mm is the month and dd the day.
 
c nstep = gcm timestep counter
c dtime = gcm timestep (s)
c mcdate= calendar day number (no leap years) (yymmdd)
c calday= day within year (1.0 at 0:0 GMT 1/1, 366. at 24:00 GMT 12/31)
c hybtype= hybrid-coordinate method (0/1/2)
c along = gcm longitude grid, east from Greenwich (radians)
c alatg = gcm latitude grid, south to north (radians)
c wgaus = gcm Gaussian weights, south to north (sum=1)
c sigg  = gcm mid-layer sigma/eta grid, high to low
c sigl  = gcm sigma/eta boundaries, high (0) to low (1)
c dsig  = gcm sigma/eta grid thicknesses
c gravit= gravitational acceleration (m s-2)
c radius= earth radius, m
c twodt = leapfrog 2*dtime timestep (s)
c ug    = gcm eastward  wind velocity, m s-1, at time n
c vg    = gcm northward wind velocity, m s-1, at time n
c wg    = gcm d[sigma/eta]/dt vertical velocity, s-1, at time n
c wgsav = time running-mean of wg, for tracer advection in stratos (s-1)
c         at time n (s-1)
c psurf = gcm surface pressure, at time n+1 (N m-2)
c pmid  = gcm pressure at mid points of sig/eta layers, time n (N m-2)
c pint  = gcm pressure at interfaces of sig/eta layers, time n (N m-2)
c pthicm= gcm pressure-thicknesses of sig/eta layers,at time n-1 (N m-2)
c pthic = gcm pressure-thicknesses of sig/eta layers,at time n (N m-2)
c pthicp= gcm pressure-thicknesses of sig/eta layers,at time n+1 (N m-2)
c albsud= gcm surface albedo for direct solar radiation (0-1)
c albsui= gcm surface albedo for diffuse solar radiation (0-1)
c solsud= gcm downward incident direct solar radiation at sfc (W m-2)
c solsui= gcm downward incident diffuse solar radiation at sfc (W m-2)
c cosq  = cosine of the solar zenith angle (0-1)
c         nb: albsu*, solsu* are for 2 wavebands: #1 = visible
c             (<.7 microns), #2 = near-ir (>.7 microns).
c         nb: albsu*, solsu*, cosq are computed only every ~1.5 hours
c             and are held constant for intervening timesteps.
c tg    = gcm air temperature, at time n (deg K)
c qg    = gcm specific humidity before advection, at time n-1
c qgnew = gcm specific humidity after advection, at time n+1
c cloudc= convective   cloud amount, at time n+1 (kg/kg)
c clouda= anvil-cirrus cloud amount, at time n+1 (kg/kg)
c clouds= stratus      cloud amount, at time n+1 (kg/kg)
 
#include <paramslt.h>
 
      integer hybtype
 
      dimension
     *  along(nlon),              alatg(nlat),
     *  wgaus(nlat),
     *  sigg(nlev),               sigl(nlev+1),
     *  dsig(nlev),
     *  ug(nlon,nlat,nlev),       vg(nlon,nlat,nlev),
     *  wg(nlon,nlat,nlev),       wgsav(nlon,nlat,nlev),
     *  psurf(nlon,nlat),
     *  pmid(nlon,nlat,nlev),     pint(nlon,nlat,nlev+1),
     *  pthicm(nlon,nlat,nlev),   pthic(nlon,nlat,nlev),
     *  pthicp(nlon,nlat,nlev),
     *  albsud(nlon,nlat,2),      albsui(nlon,nlat,2),
     *  solsud(nlon,nlat,2),      solsui(nlon,nlat,2),
     *  cosq(nlon,nlat),          tg(nlon,nlat,nlev),
     *  qg(nlon,nlat,nlev),       qgnew(nlon,nlat,nlev),
     *  cloudc(nlon,nlat,nlev),   clouda(nlon,nlat,nlev),
     *  clouds(nlon,nlat,nlev)
 
      dimension worka(nlon,nlat)
 
      character cname*25
      logical first
      save cname, first, iu, outhours
      data first, iu /.true., 95/
      data outhours /12./        ! outhours = file output interval (hrs)
 
 
      icalday  = int(calday-.0001)
      isecday  = mod (nint(nstep*dtime)-1, 86400) + 1
 
c       Open a new output file
 
      nos = nint(outhours*3600./dtime)
      if (first .or. mod(nstep,nos).eq.1) then
 
        ihourout = nint((isecday + (nos-1)*dtime)/3600.)
        write(cname,"('sltoff.',i3.3,'.',i2.2)") icalday, ihourout
#ifdef cray
        call assign ('assign -a ' //cname// ' -b 256 -s bin ' //cname)
#endif
        open (iu, file=cname, status='unknown',form='unformatted')
        write(iu) nlon, nlat, nlev, hybtype
        write(iu) along, alatg, wgaus, sigg, sigl, dsig, gravit, radius,
     *            twodt
 
        first = .false.
      endif
 
      icalday = int(calday-.0001)
      isecday = mod (nint(nstep*dtime)-1, 86400) + 1
 
c        Set surface albedo, visible waveband (set to zero where night)
 
      iband = 1
      do 10 j=1,nlat
        do 12 i=1,nlon
          zinc = solsud(i,j,iband) + solsui(i,j,iband)
          zref = albsud(i,j,iband)*solsud(i,j,iband)
     *         + albsui(i,j,iband)*solsui(i,j,iband)
          if (zinc.gt.0.1) then
            worka(i,j) = zref/zinc
          else
            worka(i,j) = 0.
          endif
   12   continue
   10 continue
 
      write(iu) icalday, isecday, nstep*dtime
      write(iu) ug
      write(iu) vg
      write(iu) wg
      write(iu) wgsav
      write(iu) psurf
      write(iu) pmid
      write(iu) pint
      write(iu) pthicm
      write(iu) pthic
      write(iu) pthicp
      write(iu) worka
      write(iu) cosq
      write(iu) tg
      write(iu) qg
      write(iu) qgnew
      write(iu) cloudc
      write(iu) clouda
      write(iu) clouds
 
      if (mod(nstep,nos).eq.0) close(iu)
 
      return
      end
c==================================================================
c Routines in this file are called directly from gcm subrs
c to set up arrays for hybrid or sigma coords before call to slthyb
c==================================================================
c
c-----------------------------------------------------------------------
c
      subroutine lpsdot (u, v, div, dpsl, dpsm, dsigma, wslt)
 
c Calculates d[ln(surface pressure)]/dt for use in hybrid-coordinate
c slt, on full (mid-layer) sigma levels.
 
c Input:
c u      = gcm eastward  velocity, m s-1
c v      = gcm northward velocity, m s-1
c div    = gcm divergence, s-1
c dpsl   = d(ln ps)/d(earth radius*longitude), m-1
c dpsm   = (1-mu**2) d(ln ps)/d(earth radius*mu), where mu=sin(lat), m-1
c dsigma = thicknesses of sigma grid layers
 
c Output:
c wslt   = gcm lagrangian "lps_dot" d[ln(surface pressure)]/dt
 
#include <paramslt.h>
 
      dimension
     *   u(nlon,nlat,nlev),  v(nlon,nlat,nlev),  div(nlon,nlat,nlev),
     *   dpsl(nlon,nlat),    dpsm(nlon,nlat),
     *   dsigma(nlev),       wslt(nlon,nlat,nlev)
 
      dimension wrk(nlon)
 
 
c        Loop over latitude
 
      do 500 j=1,nlat
 
         do 5  i=1,nlon
            wrk(i) = 0.0
    5    continue
 
c           Compute total vertical integral
 
         do 100  k=1,nlev
            zdsig = dsigma(k)
            do 102  i=1,nlon
               wrk(i) = wrk(i)
     *                + div(i,j,k)*zdsig
     *                + (u(i,j,k)*dpsl(i,j) + v(i,j,k)*dpsm(i,j))*zdsig
  102       continue
  100    continue
 
c           Add local term
 
         do 200  k=1,nlev
            do 202  i=1,nlon
               wslt(i,j,k) = - wrk(i)
     *              + (u(i,j,k)*dpsl(i,j) + v(i,j,k)*dpsm(i,j))
  202       continue
  200    continue
 
  500 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sigeta (field1, field2, j, ilat, pg, jpg, idir, iflag,
     *                   hybtype)
 
c Interpolates one longitude-level slice of a gcm (non-extended) field
c between sigma and hybrid vertical coordinates, using statement
c functions in paramhyb. Does one lon-lev slice at a time so can be
c called by linems (for 3-D fields u,v,t,w,q) and by vdif (for
c delta-reservoirs dplum*(nlon,nlev)). Calculations that depend
c only on the grids (not the fields) are done only if iflag = 1, else
c they are assumed to be the same as the previous call (saved in
c kd(nlon,nlev),hdpow(nlon,*,nlev) and ksav(nlon,*,nlev).
 
c Input:
c field1 = gcm field before interpolation
c j      = current latitude index for field1, field2
c ilat   = latitude index of field1 and field2
c pg     = gcm 2-D surface pressure (N m-2)
c jpg    = current latitude index for pg
c idir   = 1 for sigma-to-eta, = -1 for eta-to-sigma
c iflag  = 1 to compute grid-dependent terms, else use previous values
c hybtype = 0 for no hybrid, 1 for anal hybrid, 2 for discrete hybrid
 
c Output:
c field2 = gcm field after interpolation (can be same array as field1
c          in calling routine)
 
#include <paramslt.h>
#include <paramhyb.h>
 
      dimension
     *  field1(nlon,ilat,nlev),         field2(nlon,ilat,nlev),
     *  pg(nlon,nlat)
 
      integer hybtype
 
      dimension
     *  ftmp(nlon,nlev),                targ(nlon),
     *  kd(nlon,nlev),                  hd(nlon),
     *  hdpow(nlon,4,nlev),
     *  ksav(nlon,-1:2,nlev),
     *  valz(nlon,-1:2),
     *  indps(nlon)
 
c If not multiprocessed, can save the grid interpolation information
c for a latitude band and use it for subsequent calls for the same band
c (ie, iflag=0). If multiprocessed, cannot do this and must always
c call with iflag=1, and *not* save these arrays (which must be local).
 
#ifndef omp
      save kd, hdpow, ksav
#endif
 
 
c        If non-hybrid, just copy field1 to field2 and return
 
      if (hybtype.eq.0) then
        do 50 k=1,nlev
          do 52 i=1,nlon
            field2(i,j,k) = field1(i,j,k)
   52     continue
   50   continue
        return
      endif
 
c====================
      do 100 k=1,nlev
c====================
 
c       If iflag = 1, set cubic interpolation terms and indices that
c       depend on the grid, not the field. If iflag ne 1, use values
c       from previous call (hdpow, ksav) for this lon-lev slice.
 
c---------------------------
        if (iflag.eq.1) then
c---------------------------
 
c            Calculate fine-spaced surface-pressure indices indps
 
          if (k.eq.1) call getindps (indps, pg(1,jpg))
 
c           Calculate target vertical coords for a longitude strip
 
          if (idir.eq.1) then
c           sigma-to-eta:
            do 120 i=1,nlon
              targ(i) = etaval (k,indps(i))
  120       continue
          else
c           eta-to-sigma:
            do 122 i=1,nlon
              targ(i) = sigval (k,indps(i))
  122       continue
          endif
 
c           Set kd and hd for these targets (via uniform grid)
 
          call locate (targ,nlon, siga,nlevp2, nlevset,kset, kd(1,k),hd)
 
c           Compute cubic interpolation terms, save in hdpow(*,*,nlev).
c           (The following code is similar to interpv)
 
          do 130 i=1,nlon
            hdpow(i,1,k) = hd(i)**3
            hdpow(i,2,k) = hd(i)**2 * (1.-hd(i))
            hdpow(i,3,k) = hd(i)    * (1.-hd(i))**2
            hdpow(i,4,k) =            (1.-hd(i))**3
  130     continue
 
c           Compute interpolation indices, save in ksav(*,*,nlev)
 
          do 140 ksh = -1,2
            do 142 i=1,nlon
              ksav(i,ksh,k) = max (1, min (nlev, kd(i,k)+ksh))
  142       continue
  140     continue
 
c------------
        endif
c------------
 
c         Assemble valz array for current strip as needed by docube.
c         (2nd dimension of valz holds the 4 surrounding vertical pts)
 
        do 160 ksh = -1,2
          do 162 i=1,nlon
            valz(i,ksh) = field1 (i,j,ksav(i,ksh,k))
  162     continue
  160   continue
 
c         Do cubic vertical interpolation
 
        call docube (kd(1,k), hdpow(1,1,k),
     *               delt_z, faca_z, facap_z, facb_z, facbp_z,
     *               valz, ftmp(1,k))
 
c=============
  100 continue
c=============
 
c       After this lon-lev slice is done, copy ftmp to field2 (so
c       field1 and field2 can be the same array in the calling routine)
 
      do 200 k=1,nlev
        call scopy (nlon, ftmp(1,k), 1, field2(1,j,k), 1)
  200 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine etadot (wg, omega, pg, sigg, hybtype)
 
c Computes hybrid-coordinate "eta_dot" vertical velocity d[eta]/dt,
c using statement functions in paramhyb. If slt-hybrid not used,
c computes sigma-coordinate "sig_dot" vertical velocity d[p/ps]/dt.
c Called directly from linems.
 
c Input:
c wg    = gcm d[ln(surface pressure)]/dt (lagrangian), s-1
c omega = gcm d[p]/dt (lagrangian), N m-2 s-1
c pg    = gcm surface pressure, N m-2
c sigg  = gcm mid-layer sigma grid, high to low
c hybtype = 0 for no hybrid, 1 for anal hybrid, 2 for discrete hybrid
 
c Output:
c wg    = gcm d[eta]/dt for hybrid or d[p/ps]/dt for non-hybrid, s-1
 
#include <paramslt.h>
#include <paramhyb.h>
 
      dimension
     *  wg(nlon,nlat,nlev),       omega(nlon,nlat,nlev),
     *  pg(nlon,nlat),            sigg(nlev)
 
      integer hybtype
 
      dimension
     *  indps(nlon)
 
c---------------------------
      if (hybtype.eq.0) then
c---------------------------
 
c       Set to sigma_dot:
        do 100 j=1,nlat
          do 102 k=1,nlev
            do 104 i=1,nlon
              wg(i,j,k) = omega(i,j,k)/pg(i,j) - sigg(k)*wg(i,j,k)
  104       continue
  102     continue
  100   continue
 
c---------
      else
c---------
 
c       Set to eta_dot:
        do 200 j=1,nlat
          call getindps (indps, pg(1,j))
          do 202 k=1,nlev
            do 204 i=1,nlon
              wg(i,j,k)
     *          = detadp (k,indps(i)) *omega(i,j,k)
     *          + detadps(k,indps(i)) *wg(i,j,k)*pg(i,j)
  204       continue
  202     continue
  200   continue
 
c----------
      endif
c----------
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine etapres (pg, etam, etal, pm, pl, pthic, pmcap, plcap,
     *                    cappa, hybtype)
 
c Computes various pressure quantities of hybrid-coordinate grid boxes
c (global gcm arrays), using statement functions in paramhyb.
c Called directly from linems.
 
c Input:
c pg   = gcm surface pressure (N m-2)
c etam = gcm vertical eta-grid mid-points
c etal = gcm vertical eta-grid boundaries
c cappa = R/c_p (exponent for pressure dependence of potential temp)
c hybtype = 0 for no hybrid, 1 for anal hybrid, 2 for discrete hybrid
 
c Output:
c pm    = gcm pressures at eta-grid layer mid-points
c pl    = gcm pressures at eta-grid layer boundaries
c pthic = gcm pressure-thicknesses of eta-grid boxes
c pmcap = (pm/pg)**cappa
c plcap = (pl/pg)**cappa
 
#include <paramslt.h>
#include <paramhyb.h>
 
      dimension
     *  pg(nlon,nlat),           etam(nlev),           etal(nlev+1),
     *  pm(nlon,nlat,nlev),      pl(nlon,nlat,nlev+1),
     *  pthic(nlon,nlat,nlev),
     *  pmcap(nlon,nlat,nlev),   plcap(nlon,nlat,nlev+1)
 
      integer hybtype
 
      dimension
     *  indps(nlon)
 
c---------------------------
      if (hybtype.eq.0) then
c---------------------------
 
c         If non-hybrid, simply set to appropriate sigma-grid values
 
        do 100 k=1,nlev
          zzm = etam(k)**cappa
          zzl = etal(k)**cappa
          do 102 j=1,nlat
            do 104 i=1,nlon
              pm(i,j,k) = pg(i,j)*etam(k)
              pl(i,j,k) = pg(i,j)*etal(k)
              pmcap(i,j,k) = zzm
              plcap(i,j,k) = zzl
  104       continue
  102     continue
  100   continue
 
        do 110 j=1,nlat
          do 112 i=1,nlon
            pl(i,j,nlev+1) = pg(i,j)
            plcap(i,j,nlev+1) = 1.
  112     continue
  110   continue
 
c---------
      else
c---------
 
c         Calculate mid-layer pressures and layer-boundary pressures
 
        do 200 j=1,nlat
          call getindps (indps, pg(1,j))
          do 202 k=1,nlev
            do 204 i=1,nlon
              pm(i,j,k) = pg(i,j) * sigval (k,indps(i))
              pl(i,j,k) = pg(i,j) * sigvall(k,indps(i))
 
              pmcap(i,j,k) = (pm(i,j,k)/pg(i,j))**cappa
              plcap(i,j,k) = (pl(i,j,k)/pg(i,j))**cappa
  204       continue
  202     continue
  200   continue
 
        do 210 j=1,nlat
          do 212 i=1,nlon
            pl(i,j,nlev+1) = pg(i,j)
            plcap(i,j,nlev+1) = 1.
  212     continue
  210   continue
 
c----------
      endif
c----------
 
c       Set layer pressure-thickesses
 
      do 300 k=1,nlev
        do 302 j=1,nlat
          do 304 i=1,nlon
            pthic(i,j,k) = pl(i,j,k+1) - pl(i,j,k)
  304     continue
  302   continue
  300 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine etapres_a (pg, etam, etal, pthic, hybtype)
 
c Stripped-down version of etapres, just to calculate pressure
c thicknesses (global gcm array) of hybrid-coordinate boxes,
c using statement functions in paramhyb. Called directly from linems.
 
c Input:
c pg   = gcm surface pressure (N m-2)
c etam = gcm vertical eta-grid mid-points
c etal = gcm vertical eta-grid boundaries
c hybtype = 0 for no hybrid, 1 for anal hybrid, 2 for discrete hybrid
 
c Output:
c pthic = gcm pressure-thicknesses of eta-grid boxes
 
#include <paramslt.h>
#include <paramhyb.h>
 
      dimension
     *  pg(nlon,nlat),           etam(nlev),           etal(nlev+1),
     *  pthic(nlon,nlat,nlev)
 
      integer hybtype
 
      dimension
     *  indps(nlon),             pl(nlon,nlev+1)
 
 
c        Overall loop over latitude so local array pl need not be global
 
c--------------------
      do 400 j=1,nlat
c--------------------
 
c---------------------------
      if (hybtype.eq.0) then
c---------------------------
 
c         If non-hybrid, simply set to appropriate sigma-grid values
 
        do 100 k=1,nlev
            do 104 i=1,nlon
              pl(i,k) = pg(i,j)*etal(k)
  104       continue
  100   continue
 
          do 112 i=1,nlon
            pl(i,nlev+1) = pg(i,j)
  112     continue
 
c---------
      else
c---------
 
c         Calculate mid-layer pressures and layer-boundary pressures
 
          call getindps (indps, pg(1,j))
          do 202 k=1,nlev
            do 204 i=1,nlon
              pl(i,k) = pg(i,j) * sigvall(k,indps(i))
  204       continue
  202     continue
 
          do 212 i=1,nlon
            pl(i,nlev+1) = pg(i,j)
  212     continue
 
c----------
      endif
c----------
 
c       Set layer pressure-thickesses
 
      do 300 k=1,nlev
          do 304 i=1,nlon
            pthic(i,j,k) = pl(i,k+1) - pl(i,k)
  304     continue
  300 continue
 
c-------------
  400 continue
c-------------
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine getindps (indps, pg)
 
c     Computes fine-spaced surface-pressure indices for sigval[l],
c     etaval, detadp, detadps arrays for hybrid coords, for one
c     longitude strip.
 
c     indps = fine-spaced ps index for this longitude strip (returned)
c     pg    = surface pressure,N m-2, for this longitude strip(supplied)
 
#include <paramslt.h>
#include <paramhyb.h>
 
      dimension indps(nlon), pg(nlon)
 
      do 100 i=1,nlon
        indps(i) = 1 + nint((nps-1)*(pg(i)-psa)/(psb-psa))
  100 continue
 
      return
      end
c================================================================
c sltini is called from gcm to initialize various slt grid arrays
c================================================================
c
c-----------------------------------------------------------------------
c
      subroutine sltini (long, latg, sigg, sigl, hybtype)
 
c Initialization for slt. Called once before any slt or hybrid-sigma
c related calls. Sets extended grids (lona, lata, siga), uniform
c fine-spaced indices (iset,jset,kset), and cubic interpolation factors
c that depend only on the grid and not on the field (delt_*, fac*_*).
c These are all in commons in paramslt.
 
c long  = gcm longitude grid, radians, greenwich+.. to greenwich+2*pi-..
c latg  = gcm latitude  grid, radians, south to north
c sigg  = gcm mid-layer sigma/eta grid, high     to low
c sigl  = gcm interface sigma/eta grid, high (0) to low (1)
c hybtype = 0 for no hybrid, 1 for anal hybrid, 2 for discrete hybrid
 
#include <paramslt.h>
#include <paramhyb.h>
 
      integer hybtype
 
      dimension
     *  long(nlon),  latg(nlat),  sigg(nlev), sigl(nlev+1)
      real long, latg
 
c          Set extended grids *a
 
      call extgrd (long,latg,sigg)
 
c        Set uniform (fine-spaced) grid indices *set
 
      call locset (lona,  nlonp2, nlonset, iset)
      call locset (lata,  nlatp2, nlatset, jset)
      call locset (siga,  nlevp2, nlevset, kset)
 
c        Set verious cubic interpolation terms for
c        extcube that depend only on the grid, not on the field.
 
      call setcube
c
c        Set sigval,etc, arrays in paramhyb for hybrid coords
c        (analytic or discrete type)
 
      iuchek = 0
 
      if (hybtype.eq.1) then
        call hybini_a (sigg, sigl, sigg, sigl)
      else if (hybtype.eq.2) then
        call hybini_d (sigg, sigl, sigg, sigl, iuchek)
      endif
 
c        Write out sigval,etc, arrays as a check (if iuchek ne 0)
 
      if (iuchek.ne.0) then
        npsout = 10
        do 500 ips = 1,npsout+1
          if (ips.le.npsout) then
            ps = psa + (psb-psa)*(ips-1.)/(npsout-1.)
          else
            ps = p0
          endif
          kps = 1 + nint((nps-1)*(ps-psa)/(psb-psa))
          write(iuchek,*)
          write(iuchek,*)
     *      ' k, sig, sigval, etaval, detadp, detadps for ps=',ps,
     *      '  (kps=',kps,')'
          write(iuchek,*)
          do 502 k=1,nlev+1
            write(iuchek,510) k, sigl(k), sigvall(k,kps)
            if (k.le.nlev)
     *        write(iuchek,510)k, sigg(k), sigval(k,kps), etaval(k,kps),
     *                         detadp(k,kps), detadps(k,kps)
  502     continue
  500   continue
  510   format(i5,f8.4,4f15.8)
      endif
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine extgrd (long,latg,sigg)
 
c Sets extended grids *a (in common in paramslt) from gcm grids *g.
 
c Input:
c long  = gcm longitude grid
c latg  = gcm latitude  grid
c sigg  = gcm mid-layer sigma grid
 
#include <paramslt.h>
 
      dimension
     *  long(nlon),               latg(nlat),
     *  sigg(nlev)
      real long, latg
 
c        Copy and extend longitude grid by two wraparound pt on each end
 
      call scopy (nlon, long, 1, lona(1), 1)
      lona(0)      = 2*lona(1)      - lona(2)
      lona(-1)     = 2*lona(0)      - lona(1)
      lona(nlon+1) = 2*lona(nlon)   - lona(nlon-1)
      lona(nlon+2) = 2*lona(nlon+1) - lona(nlon)
 
      do 10 i=-1,nlonp2
        slona(i) = sin(lona(i))
        clona(i) = cos(lona(i))
   10 continue
 
c        Copy and extend latitude grid by two extra points, at pole and
c        one beyond
 
      call scopy (nlat, latg, 1, lata(1), 1)
      lata(0)      = -pi/2.
      lata(-1)     = 2*lata(0)      - lata(1)
      lata(nlat+1) = pi/2.
      lata(nlat+2) = 2*lata(nlat+1) - lata(nlat)
 
      do 20 j=-1,nlatp2
        slata(j) = sin(lata(j))
        clata(j) = cos(lata(j))
   20 continue
 
c        Copy and extend mid-layer sigma grid by two extra points.
 
      call scopy (nlev, sigg, 1, siga(1), 1)
      siga(0)      = 0.
      siga(-1)     = 2*siga(0) - siga(1)
      siga(nlev+1) = 1.
      siga(nlev+2) = 2*siga(nlev+1) - siga(nlev)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine locset (arr, narr, nset, iset)
 
c Sets indices iset into a non-uniform extended grid arr for
c nset uniformly spaced points between arr(-1) and arr(narr). Then in
c subr locate, index values are obtained by first calculating the
c location in the uniform grid and then using iset. This is a lot
c faster than searching in the non-uniform grid directly.
 
c Input:
c arr  = non-uniform extended grid
c narr = max subscript for arr, ie, arr is dimensioned arr(-1:narr)
c nset = number of uniformly spaced points (much larger than narr)
 
c Output:
c iset = index of closest arr grid pt to "left" of each uniform grid pt
 
      dimension arr(-1:narr), iset(nset)
 
 
c        Set uniform grid interval
 
      dtarset = (arr(narr)-arr(-1)) / nset
 
c        For each uniform grid pt (=center of uniform boxes), do simple
c        search in arr. Assumes arr(j) values increase with j.
 
      do 100 i=1,nset
        tarset = (i-0.5)*dtarset + arr(-1)
        do 102 j=0,narr
          if (tarset.le.arr(j)) then
            iset(i) = j-1
            goto 100
          endif
  102   continue
  100 continue
 
c        Impose constraints on final iset values to ensure that the id's
c        from subr locate will be ge 0 and le nlon,nlat or nlev. (This
c        is ok due to the restricted ranges of all departure points
c        calculated by subr vector.)  Then id+ish (and jd+jsh,
c        kd+ksh) in subrs interp[h,v] will always be within -1 to
c        nlonp2,nlatp2 or nlevp2.
 
      do 200 i=1,nset
        if (iset(i).lt.0) then
          iset(i) = 0
        else if (iset(i).gt.narr-2) then
          iset(i) = narr-2
        endif
  200 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setcube
 
c        Set various terms for the cubic interpolation factors in
c        extcube, that depend *only* on the grid, not on the field.
c        (These terms delt*, fac* are in common in paramslt.)
 
#include <paramslt.h>
 
      do 10 i=-1,nlon+1
        delt_x(i) = lona(i+1)-lona(i)
   10 continue
      do 12 j=-1,nlat+1
        delt_y(j) = lata(j+1)-lata(j)
   12 continue
      do 14 k=-1,nlev+1
        delt_z(k)  = siga(k+1) -siga(k)
   14 continue
 
 
      do 20 i=0,nlon
        fa = (lona(i)-lona(i-1)) * (lona(i)-lona(i+1))
     *                           / (lona(i+2)-lona(i-1))
        faca_x(i)  = (lona(i)-lona(i-1)) / (lona(i+1)-lona(i-1))
     *             - fa / (lona(i+1)-lona(i-1))
        facap_x(i) = fa / (lona(i+2)-lona(i))
 
 
        fb = (lona(i+1)-lona(i-1)) * (lona(i+1)-lona(i))
     *                             / (lona(i+2)-lona(i-1))
        facb_x(i)  = (2*lona(i+1)-lona(i)-lona(i-1))
     *               / (lona(i+1)-lona(i-1))
     *             - fb / (lona(i+1)-lona(i-1))
        facbp_x(i) = fb / (lona(i+2)-lona(i))
   20 continue
 
 
      do 22 j=0,nlat
        fa = (lata(j)-lata(j-1)) * (lata(j)-lata(j+1))
     *                           / (lata(j+2)-lata(j-1))
        faca_y(j)  = (lata(j)-lata(j-1)) / (lata(j+1)-lata(j-1))
     *             - fa / (lata(j+1)-lata(j-1))
        facap_y(j) = fa / (lata(j+2)-lata(j))
 
 
        fb = (lata(j+1)-lata(j-1)) * (lata(j+1)-lata(j))
     *                             / (lata(j+2)-lata(j-1))
        facb_y(j)  = (2*lata(j+1)-lata(j)-lata(j-1))
     *               / (lata(j+1)-lata(j-1))
     *             - fb / (lata(j+1)-lata(j-1))
        facbp_y(j) = fb / (lata(j+2)-lata(j))
   22 continue
 
 
      do 24 k=0,nlev
        fa = (siga(k)-siga(k-1)) * (siga(k)-siga(k+1))
     *                           / (siga(k+2)-siga(k-1))
        faca_z(k)  = (siga(k)-siga(k-1)) / (siga(k+1)-siga(k-1))
     *             - fa / (siga(k+1)-siga(k-1))
        facap_z(k) = fa / (siga(k+2)-siga(k))
 
 
        fb = (siga(k+1)-siga(k-1)) * (siga(k+1)-siga(k))
     *                             / (siga(k+2)-siga(k-1))
        facb_z(k)  = (2*siga(k+1)-siga(k)-siga(k-1))
     *             / (siga(k+1)-siga(k-1))
     *             - fb / (siga(k+1)-siga(k-1))
        facbp_z(k) = fb / (siga(k+2)-siga(k))
   24 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine hybini_d (eta, etal, sig, sigl, iuchek)
 
c       Sets hybrid transformation arrays for ***discrete*** hybrid
c       system, defined by sig(k)= p/ps = a(eta(k))*(p0/ps) + b(eta(k)).
c       These arrays are in common in paramhyb, and are used in
c       subrs sigeta, etaval and etapres. They are:
 
c     * sigval, giving sigma (=p/ps) value for given mid-pt eta-grid
c       value and ps index.
 
c     * sigvall, giving sigma (=p/ps) value for given interface eta-grid
c       value and ps index.
 
c     * etaval, giving eta for given mid-pt sigma-grid value and ps
c       index.
 
c     * detadp, giving d[eta]/d[p] (at constant ps) for given mid-point
c       sigma-grid value and ps index.
 
c     * detadps, giving d[eta]/d[ps] (at constant ps) for given
c       mid-point sigma-grid value and ps index.
 
c      etaval, detadp and detadps are dimensioned (nlev,nps), where
c      nlev is for the midpoint grid values and nps is for fine-spaced
c      ps values, equispaced from psa to psb (Pascals).
 
c      The gcm eta grid (=sigma grid) is in arrays eta,etal (=sig,sigl).
 
c      The hybrid transformation is from the 18-level ccm2 grid, read
c      by ~/Hybrid/readccm2.f, which reads a ccm2 initial file and
c      writes the values in the data statements below for araw, braw
c      and etaraw. (These values are an interleaving of the ccm2
c      mid-point and interface-point values.)
c         p/ps = araw(etaraw)*(p0/ps) + braw(etaraw)
c      where ps=surface pressure, p0=gcm reference pressure.
 
c      Then a, b, d[a]/d[eta], d[b]/d[eta] for any given eta are
c      linearly interpolated in eta from the raw values.

c      cf: a18 and b18 settings for ccm3 in subr set_eta in:
c      ~/Csm2.0.2/cam1/models/atm/cam/src/dynamics/fv/dynamics_vars.F90
 
#include <paramslt.h>
#include <paramhyb.h>
 
      parameter (nraw=2*18+1, nloop=25)
 
      dimension
     *  eta (nlev ), etal(nlev+1), sig (nlev), sigl(nlev+1)
 
      dimension
     *  araw(nraw),       braw(nraw),        etaraw(nraw),
     *  aval (nlev ),     bval (nlev),
     *  avall(nlev+1),    bvall(nlev+1)
 
      data araw /
     *   0.002917, 0.004809, 0.007929, 0.013073, 0.021554, 0.032559,
     *   0.049183, 0.063947, 0.083143, 0.081677, 0.079931, 0.078020,
     *   0.075774, 0.073367, 0.070575, 0.067648, 0.064296, 0.060862,
     *   0.056984, 0.053110, 0.048791, 0.044600, 0.039989, 0.035661,
     *   0.030963, 0.026727, 0.022190, 0.018307, 0.014204, 0.010942,
     *   0.007541, 0.005147, 0.002684, 0.001352, 0.000000, 0.000000,
     *   0.000000 /
 
      data braw /
     *   0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     *   0.000000, 0.000000, 0.000000, 0.017366, 0.038054, 0.060693,
     *   0.087309, 0.115824, 0.148904, 0.183592, 0.223300, 0.263985,
     *   0.309941, 0.355846, 0.407010, 0.456676, 0.511298, 0.562588,
     *   0.618247, 0.668443, 0.722193, 0.768203, 0.816817, 0.855465,
     *   0.895759, 0.924129, 0.953314, 0.969094, 0.985112, 0.992528,
     *   1.000000 /
 
      data etaraw /
     *   0.002917, 0.004809, 0.007929, 0.013073, 0.021554, 0.032559,
     *   0.049183, 0.063947, 0.083143, 0.099043, 0.117985, 0.138713,
     *   0.163083, 0.189191, 0.219479, 0.251239, 0.287596, 0.324848,
     *   0.366924, 0.408955, 0.455801, 0.501275, 0.551287, 0.598248,
     *   0.649210, 0.695169, 0.744383, 0.786510, 0.831021, 0.866407,
     *   0.903300, 0.929275, 0.955997, 0.970446, 0.985112, 0.992528,
     *   1.000000 /
 
 
c        Calculate a,b on model eta grid (for "sigval")
 
      do 100 k=1,nlev+1
        if (k.le.nlev) then
          aval(k) = abterp (eta(k), etaraw, araw, nraw)
          bval(k) = abterp (eta(k), etaraw, braw, nraw)
        endif
        avall(k) = abterp (etal(k), etaraw, araw, nraw)
        bvall(k) = abterp (etal(k), etaraw, braw, nraw)

c       5/07: avoid avall(1) < 0 due to roundoff:
        avall(k) = max (avall(k),0.)
        bvall(k) = max (bvall(k),0.)
  100 continue
 
c------------------------------------
c Outer loop on fine-spaced ps values
c------------------------------------
      do 200 kps = 1,nps
        ps = psa + (psb-psa)*(kps-1.)/(nps-1.)
 
c          Find eta that corresponds to given sigma and ps, store those
c          eta's for full-sigma grid and nps fine-spaced ps-values
c          between psa and psb (for "etaval", "detadp", "detadps")
 
        do 210 k=1,nlev
 
c           Do binary search (nloop) for the eta that corresponds to
c           each sigma mid-layer grid value
 
          elo = 0.
          ehi = 1.
          do 212 iloop=1,nloop
            emid=0.5*(elo+ehi)
            sigmid =  abterp(emid,etaraw,araw,nraw)*(p0/ps)
     *              + abterp(emid,etaraw,braw,nraw)
            if (sigmid.le.sig(k)) then
              elo = emid
            else
              ehi = emid
            endif
  212     continue
 
c            Store the eta value for this sigma and ps
 
          etaval(k,kps) = emid
 
c            Compute and store deta/dp, deta/dps for this sigma and ps
 
          einc = 0.01
          ezp = min (emid+einc, 1.)
          ezm = max (emid-einc, 0.)
          dade = (  abterp(ezp,etaraw,araw,nraw)
     *            - abterp(ezm,etaraw,araw,nraw) ) / (ezp-ezm)
          dbde = (  abterp(ezp,etaraw,braw,nraw)
     *            - abterp(ezm,etaraw,braw,nraw) ) / (ezp-ezm)
          bz = abterp(emid,etaraw,braw,nraw)
 
          detadp (k,kps) =   1./ (dade*p0 + dbde*ps)
          detadps(k,kps) = -bz / (dade*p0 + dbde*ps)
 
  210   continue
 
c          Compute and store sigval, sigvall for eta grid and this ps
 
        do 220 k=1,nlev
          sigval (k,kps) = aval (k)*(p0/ps) + bval (k)
          sigvall(k,kps) = avall(k)*(p0/ps) + bvall(k)
  220   continue
        sigvall(nlev+1,kps) = avall(nlev+1)*(p0/ps) + bvall(nlev+1)
 
c-------------
  200 continue
c-------------
 
c        Write out computed arrays as a check (if iuchek ne 0)
 
      if (iuchek.ne.0) then
 
c       Report *raw values used:
        write (iuchek,*)
        write (iuchek,*) 'Echo of ccm2-values for araw, braw, etaraw:'
        write (iuchek,300) (araw(k),braw(k),etaraw(k), k=1,nraw)
  300   format(//(5x,f10.6,2f20.6))
 
        write(iuchek,*)
        write(iuchek,*) 'k, eta, aval, bval (on interfaces and midpts):'
        write(iuchek,*)
        do 400 k=1,nlev+1
          write(iuchek,410)
     *    k, etal(k), avall(k), bvall(k), avall(k)+bvall(k)
          if (k.le.nlev)
     *      write(iuchek,410)
     *      k, eta(k),  aval(k),  bval(k),  aval(k) +bval(k)
  400   continue
  410   format(i5,f10.4,3f20.4)
 
      endif
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine hybini_a (eta, etal, sig, sigl)
 
c       Sets hybrid transformation arrays for *analytic* hybrid
c       system, defined by eta = (p/ps)*p/ps + (1-(p/ps))*p/p0
c       ("linear-transition", Simmons and Burridge,1981,MWR,109,pg.762).
c       These arrays (sigval[l],etaval,detadp[s]) are in common
c       in paramhyb, and are used in subrs sigeta, etaval and etapres.
c       They are described further in comments in hybini_d.
 
c       sigval[l] = sigma value for given eta and ps index
c       etaval = eta value for given sigma and ps index
c       detadp = d(eta)/d(p)  at const ps for given sigma and ps index
c       detadps= d(eta)/d(ps) at const p  for given sigma and ps index
 
 
#include <paramslt.h>
#include <paramhyb.h>
 
      dimension
     *  eta (nlev ), etal(nlev+1), sig (nlev), sigl(nlev+1)
 
c-------------------------------------------------
c     Statement functions (z=sig or eta, a=ps/p0):
 
      fetaval (z,a) = (1.-a)*(z**2) + a*z
 
      fsigval (z,a) =cvmgt( (-.5*a + sqrt((.5*a)**2+(1.-a)*z)) / (1.-a),
     *                      z,
     *                      abs(1.-a).gt.epsilon )
 
      fdetadp (z,a) = ( 2.*z*(1./a - 1.) + 1. ) / p0
 
      fdetadps(z,a) = ( (z**2)*(-2./a + 1.) ) / p0
c-------------------------------------------------
 
 
c        Outer loop on fine-spaced ps values
 
      do 200 kps = 1,nps
        ps = psa + (psb-psa)*(kps-1.)/(nps-1.)
 
        do 210 k=1,nlev
 
c         sig and eta are identical, distinguished here just for clarity
 
          etaval(k,kps) = fetaval (sig(k), ps/p0)
 
          sigval (k,kps) = fsigval (eta (k), ps/p0)
          sigvall(k,kps) = fsigval (etal(k), ps/p0)
 
          detadp (k,kps) = fdetadp (sig(k), ps/p0)
          detadps(k,kps) = fdetadps(sig(k), ps/p0)
 
  210   continue
 
c       Protect against truncation error (-ve fsigval) for etal(1) = 0.
        if (etal(1).lt.epsilon) sigvall(1,kps) = 0.
 
        sigvall(nlev+1,kps) = fsigval (etal(nlev+1), ps/p0)
 
  200 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      function abterp (xval, x, y, n)
 
c     Linearly interpolates y(x), returns interpolated y at x=xval.
c     (x(n) must be monotonically increasing.) (Used by hybini_d).
 
      dimension x(n), y(n)
 
      do 10 i=2,n
        if (xval.le.x(i)) then
          ip = i
          goto 12
        endif
   10 continue
      ip = n
   12 im = ip-1
 
      wei = (xval-x(im))/(x(ip)-x(im))
      abterp = (1.-wei)*y(im) + wei*y(ip)
 
      return
      end

c-----------------------------------------------------------------------

c-------------
#ifdef isotrac
c-------------
      subroutine traglo_i (qnew, qold, qnew_t, qold_t,
     *                     pthic, wgaus, gravit, sqmfx)
 
c TRAGLO_I does tracer ratios extrema fixing and additional mass fixing
c after transport in SLT
 
c qnew   : specific humidity field after transport and mass fixing
c qnew_t : same for tracer (MODIFIED)
c qold   : specific humidity field before transport
c qold_t : same for tracer
c pthic  : gcm pressure-thicknesses of sig/eta layers, current (N m-2)
c wgaus  : gcm Gaussian weights, south to north (sum=2)
c gravit : gravitational acceleration (m s-2)
c sqmfx  : tracer mass after transport and mass fixing

#include <paramslt.h>

      dimension pthic (nlon,nlat,nlev), wgaus (nlat),
     *          qnew  (nlon,nlat,nlev), qold  (nlon,nlat,nlev),
     *          qnew_t(nlon,nlat,nlev), qold_t(nlon,nlat,nlev),
     *          qmin_t(nlon,nlat,nlev), qmax_t(nlon,nlat,nlev),
     *          r     (nlon,nlat,nlev)
 
c Determine tracer ratios before transport
      do k=1,nlev
        do j=1,nlat
          do i=1,nlon
           if (qold(i,j,k).ne.0.) r(i,j,k) = qold_t(i,j,k)/qold(i,j,k)
           if (abs(qold(i,j,k)).gt.1.e6) then
             write (6,*)'Warning:traglo_i: i,j,k=',i,j,k,
     *                  ' qold=',qold(i,j,k)
           endif
          enddo
        enddo
      enddo
 
c Determine min/max ratios and associated specific humidities for each
c point from the point itself and the 26 surrounding points
      sqmmm = 0.
      do k=1,nlev
        do j=1,nlat
          do i=1,nlon
            rmin =  1.e20
            rmax = -1.e20
            do kk=k-1,k+1
              lk = max(min(kk,nlev),1)
              do kj=j-1,j+1
                lj = max(min(kj,nlat),1)
                do ki=i-1,i+1
                  if      (ki.lt.1)    then
                    li = nlon
                  else if (ki.gt.nlon) then
                    li = 1
                  else
                    li = ki
                  endif
                  if (qold(li,lj,lk).ne.0.) then
                    rmin = min(rmin, r(li,lj,lk))
                    rmax = max(rmax, r(li,lj,lk))
                  endif
                enddo
              enddo
            enddo
            if (rmin.le.rmax) then
              qmin_t(i,j,k) = rmin*qnew(i,j,k)
              qmax_t(i,j,k) = rmax*qnew(i,j,k)
            else
              qmin_t(i,j,k) = qnew_t(i,j,k)
              qmax_t(i,j,k) = qnew_t(i,j,k)
            endif
 
c Fix tracer ratios extrema
            qnew_t(i,j,k) = min(max(qnew_t(i,j,k),
     *                              qmin_t(i,j,k)), qmax_t(i,j,k))
 
c Compute tracer global min/max mass difference product
            sqmmm = sqmmm + (qmax_t(i,j,k)-qnew_t(i,j,k))*(qnew_t(i,j,k)
     *                      -qmin_t(i,j,k))*pthic(i,j,k)*wgaus(j)
          enddo
        enddo
      enddo
      sqmmm = sqmmm/(2*nlon*gravit)
 
c Compute tracer global mass after ratios extrema fixing
      call mascalc (qnew_t, wgaus, pthic, gravit, sqrfx)
 
c Fix global mass with a method that conserves tracer ratios extrema
c     if (sqmmm.ne.0. .and. sqmfx.ne.sqrfrx) then
      if (sqmmm.ne.0. .and. sqmfx.ne.sqrfx) then
        sqr = (sqmfx-sqrfx)/sqmmm
        do k=1,nlev
          do j=1,nlat
            do i=1,nlon
              qnew_t(i,j,k) = qnew_t(i,j,k) +
     *                       (qmax_t(i,j,k)-qnew_t(i,j,k))*
     *                       (qnew_t(i,j,k)-qmin_t(i,j,k))*sqr
            enddo
          enddo
        enddo
      endif
 
      return
      end
 
c-----
#endif
c-----

