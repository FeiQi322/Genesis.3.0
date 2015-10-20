c -- $Id: vdif.F,v 1.2 2001/07/03 21:29:22 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/vdif.F,v $
c -- $Name:  $

c
c****************************** VDIF_SIG *******************************
c
      subroutine vdif_sig (yum1, yu, yvm1, yv, ytm1, yt, yq,
     *                     fplumu, fplumv, fplumt, j)
 
c     Does time-implicit vertical diffusion calculations on sigma grid
c     for u,v,t, one latitude band at a time, for leapfrog timestep n-1
c     to n+1. 3-D forcing for each field from surface_fluxes+PBL+plume
c     have been stored in reservoirs FPLUM* on eta coords, and portions
c     of these are applied in this timestep (DPLUM* from subr reserv),
c     interpolated from eta to sigma coords, and passed to subr
c     vertdifn.
c
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C---------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      dimension
     * yum1(nlon,norec,nlev),    yu(nlon,norec,nlev),
     * yvm1(nlon,norec,nlev),    yv(nlon,norec,nlev),
     * ytm1(nlon,norec,nlev),    yt(nlon,norec,nlev),
     * yq(nlon,norec,nlev),
     * fplumu(nlon,norec,nlev),  fplumv(nlon,norec,nlev),
     * fplumt(nlon,norec,nlev)
c----------------------------------------------------------------------
      dimension
     *  dm(nlon,nlev),           vertk(nlon,nlevp),
     *  dz(nlon,nlev),           dz_h(nlon,nlevp),
     *  rho(nlon,nlev),          rho_h(nlon,nlevp),
     *  dmcp(nlon,nlev),         rhocp_h(nlon,nlevp),
     *  zthet(nlon,nlev),
     *  dplumu(nlon,nlev),       dplumv(nlon,nlev),   dplumt(nlon,nlev),
     *  zplumu0(nlon),           zplumv0(nlon),       zplumt0(nlon),
     *  zplumu1(nlon),           zplumv1(nlon),       zplumt1(nlon),
     *  unew(nlon,nlev),         vnew(nlon,nlev),     tnew(nlon,nlev),
     *  conev(nlev),             conel(nlevp),        ztmp(nlon)
c----------------------------------------------------------------------
      save conev, conel
      data conev /nlev*1./
      data conel /nlevp*1./
c     rimp is time-implicit fraction in vertdifn for u,v,t
      save rimp
      data rimp /1./
c----------------------------------------------------------------------
 
c        Compute rho and cp*rho at full levels, and vertical layer
c        thicknesses, using fields at time n
 
      do 100 k=1,nlev
        do 102 i=1,nlon
          rho(i,k) = sig(k)*yps(i,j)
     *             / (rair*(1.+zvir*yq(i,j,k))*yt(i,j,k))
          dm(i,k) = dsigma(k)*yps(i,j) / gravit
          dz(i,k) = dm(i,k) / rho(i,k)
          dmcp(i,k) = dm(i,k)*cpair*(1.+cpvir*yq(i,j,k))
  102   continue
  100 continue
 
c        Compute rho and cp*rho at half levels, and vertical
c        thicknesses between layer centers. Doesn't matter what
c        these values are at top and bottom (k=1 and nlevp), since
c        vertk will be zero there.
 
      do 110 k=1,nlevp
        km1 = max (k-1,1)
        kk  = min (k,nlev)
        do 112 i=1,nlon
          rho_h(i,k) = 0.5*(rho(i,km1)+rho(i,kk))
          dz_h(i,k) = 0.5*(dz(i,km1) + dz(i,kk))
          rhocp_h(i,k) = rho_h(i,k)*cpair
  112   continue
  110 continue
 
c-------------------
 
c        Calculate vertical diffusion coefficients at sigma half levels,
c        using u,v,t,q fields at time n. First set potential temperature
c        at full levels.
 
      do 150 k=1,nlev
        do 152 i=1,nlon
          zthet(i,k) = yt(i,j,k)*(1.+zvir*yq(i,j,k)) / sigkpk(k)
  152   continue
  150 continue
 
      call vertcoef (vertk, yu, yv, zthet, dz_h, j)
 
c-------------------
 
c        Extract some of reservoirs fplum* into flux convergences dplum*
c        (still on eta grid)
 
      call reserv (fplumu, dplumu, dtime, j)
      call reserv (fplumv, dplumv, dtime, j)
      call reserv (fplumt, dplumt, dtime, j)
 
c        Add radiative heating to dplumt (still on eta grid)
 
      do 200 k=1,nlev
        do 202 i=1,nlon
          dplumt(i,k) = dplumt(i,k) + yt2(i,j,k)
  202   continue
  200 continue
 
c        Interpolate delta-reservoirs from eta to sigma grid (sigeta),
c        applying a vertically uniform correction to conserve column
c        totals at each grid point. (If hybrid-slt not used, SIGETA
c        merely copies 1st array to 2nd.)
 
      call zero (zplumu0, nlon)
      call zero (zplumv0, nlon)
      call zero (zplumt0, nlon)
      do 250 k=1,nlev
        do 252 i=1,nlon
          zplumu0(i) = zplumu0(i) + dplumu(i,k)
          zplumv0(i) = zplumv0(i) + dplumv(i,k)
          zplumt0(i) = zplumt0(i) + dplumt(i,k)
  252   continue
  250 continue
 
      call sigeta (dplumu, dplumu, 1, 1, yps, j, -1, 1, hybtype)
#ifdef omp
      call sigeta (dplumv, dplumv, 1, 1, yps, j, -1, 1, hybtype)
      call sigeta (dplumt, dplumt, 1, 1, yps, j, -1, 1, hybtype)
#else
      call sigeta (dplumv, dplumv, 1, 1, yps, j, -1, 0, hybtype)
      call sigeta (dplumt, dplumt, 1, 1, yps, j, -1, 0, hybtype)
#endif
 
      call zero (zplumu1, nlon)
      call zero (zplumv1, nlon)
      call zero (zplumt1, nlon)
      do 260 k=1,nlev
        do 262 i=1,nlon
          zplumu1(i) = zplumu1(i) + dplumu(i,k)
          zplumv1(i) = zplumv1(i) + dplumv(i,k)
          zplumt1(i) = zplumt1(i) + dplumt(i,k)
  262   continue
  260 continue
 
      do 270 k=1,nlev
        do 272 i=1,nlon
          dplumu(i,k) = dplumu(i,k) + (zplumu0(i)-zplumu1(i))*dsigma(k)
          dplumv(i,k) = dplumv(i,k) + (zplumv0(i)-zplumv1(i))*dsigma(k)
          dplumt(i,k) = dplumt(i,k) + (zplumt0(i)-zplumt1(i))*dsigma(k)
  272   continue
  270 continue
 
c-------------------
 
c        Calculate vertical difn for t
 
      call vertdifn (tnew, ytm1, dplumt,
     *               dmcp, dz_h, rhocp_h, vertk, sigkpk, sigkpl,
     *               twodt, rimp, j)
 
c        Store heating rate (W m-2 per layer) in yt2 (don't increment
c        yt2 since radiative heating was included in vertdifn forcing),
c        and compute vertical heat diffusion budget term
 
      do 300 k=1,nlev
        do 302 i=1,nlon
          ztmp(i) = (tnew(i,k)-ytm1(i,j,k)) * dmcp(i,k)
          yt2(i,j,k) = ztmp(i)/twodt
  302   continue
        dhvdif(j) = dhvdif(j) + ssum(nlon,ztmp,1)*0.5/nlon
  300 continue
 
c-------------------
 
c        Calculate vertical diffusion for u and v
 
      call vertdifn (unew, yum1, dplumu,
     *               dm, dz_h, rho_h, vertk, conev, conel,
     *               twodt, rimp, j)
 
      call vertdifn (vnew, yvm1, dplumv,
     *               dm, dz_h, rho_h, vertk, conev, conel,
     *               twodt, rimp, j)
 
c        Store tendencies in yfu and yfv
 
      do 400 k=1,nlev
        do 402 i=1,nlon
          yfu(i,j,k) = (unew(i,k)-yum1(i,j,k))/twodt
          yfv(i,j,k) = (vnew(i,k)-yvm1(i,j,k))/twodt
  402   continue
  400 continue
 
c        Add momentum diffusive heating to net heating rate (W/m2/layer)
c        and compute vertical diffusion ke dissipation budget term
 
      do 410 k=1,nlev
        do 412 i=1,nlon
           ztmp(i) = 0.5*(  yum1(i,j,k)**2 - unew(i,k)**2
     *                    + yvm1(i,j,k)**2 - vnew(i,k)**2 )
     *               *dsigma(k)*yps(i,j)/gravit
           yt2(i,j,k) = yt2(i,j,k) + ztmp(i)/twodt
 412    continue
        dkevdf(j) = dkevdf(j) + ssum(nlon,ztmp,1)*0.5/nlon
 410  continue
 
      return
      end
c
c****************************** VDIF_ETA *******************************
c
      subroutine vdif_eta (yqm1, yq, yu, yv, yt, fplumq, tracer,
     *                     ypm, ypthic, ypmcap, j 
#ifdef isotrac
     *                    ,yqm1_o, yqm1_h, yq_o, yq_h,
     *                     fplumq_o, fplumq_h, ylwc_c_o, ylwc_c_h,
     *                     ylwc_a_o, ylwc_a_h, ylwc_s_o, ylwc_s_h,
     *                     rtrqm1, rtrq, rtrfq
#endif
     *                    )

 
c     Does time-implicit vertical diffusion calculations on eta grid
c     for q, one latitude band at a time, for leapfrog timestep n-1
c     to n+1. Also does tracers, for non-leapfrog timestep.
c     3-D forcing for each field (except tracers) from
c     surface_fluxes+PBL+plume have been stored in reservoirs FPLUM*
c     on eta coords, and portions of these are applied in this timestep
c     (DPLUM* from subr reserv) and passed to subr vertdifn.
c
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C---------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      dimension
     * yqm1(nlon,norec,nlev),    yq(nlon,norec,nlev),
     * yu(nlon,norec,nlev),      yv(nlon,norec,nlev),
     * yt(nlon,norec,nlev),
     * fplumq(nlon,norec,nlev),  tracer(nlon,norec,nlev,ntrace),
     * ypm(nlon,norec,nlev),     ypthic(nlon,norec,nlev),
     * ypmcap(nlon,norec,nlev)
#ifdef isotrac
      dimension
     * yqm1_o(nlon,norec,nlev),  yqm1_h(nlon,norec,nlev),
     * yq_o(nlon,norec,nlev),    yq_h(nlon,norec,nlev),
     * fplumq_o(nlon,norec,nlev),fplumq_h(nlon,norec,nlev),
     * ylwc_c_o(nlon,norec,nlev),ylwc_c_h(nlon,norec,nlev),
     * ylwc_a_o(nlon,norec,nlev),ylwc_a_h(nlon,norec,nlev),
     * ylwc_s_o(nlon,norec,nlev),ylwc_s_h(nlon,norec,nlev),
     * rtrqm1(nlon,norec,nlev,nrtr),rtrq(nlon,norec,nlev,nrtr),
     * rtrfq(nlon,norec,nlev,nrtr)
#endif
c----------------------------------------------------------------------
      dimension
     *  dm(nlon,nlev),           vertk(nlon,nlevp),
     *  dz(nlon,nlev),           dz_h(nlon,nlevp),
     *  rho(nlon,nlev),          rho_h(nlon,nlevp),
     *  zthet(nlon,nlev),
     *  dplumq(nlon,nlev),       dpluma(nlon,nlev),
     *  qnew(nlon,nlev),         anew(nlon,nlev),
     *  conev(nlev),             conel(nlevp),
     *  vzero(nlon,nlev),        ztmp(nlon),
     *  ztot1(nlon),             ztot2(nlon)
#ifdef isotrac
      dimension
     * rq(nlon,nlev),           ktrop(nlon),
     * dplumq_o(nlon,nlev),     dplumq_h(nlon,nlev),
     * qnew_o(nlon,nlev),       qnew_h(nlon,nlev),
     * anew_o(nlon,nlev),       anew_h(nlon,nlev),
     * ztot1_o(nlon),           ztot1_h(nlon),
     * ztot2_o(nlon),           ztot2_h(nlon),
     * rmax_o(nlon),            rmax_h(nlon),
     * rmin_o(nlon),            rmin_h(nlon),
     * dplumq_t(nlon,nlev,nrtr),qnew_t(nlon,nlev,nrtr),
     * ztot1_t(nlon,nrtr),      ztot2_t(nlon,nrtr),
     * ztot2s_t(nrtr),
     * rmax_t(nlon,nrtr),       rmin_t(nlon,nrtr)
#endif
c----------------------------------------------------------------------
      save conev, conel, vzero
      data conev /nlev*1./
      data conel /nlevp*1./
      parameter (izero = nlon*nlev)
      data vzero /izero*0./
c     rimq is time-implicit fraction in vertdifn for q,tracers
      save rimq
      data rimq /.5/
c----------------------------------------------------------------------
 
c        Compute rho at full levels, and vertical layer thicknesses,
c        using fields at time n
 
      do 100 k=1,nlev
        do 102 i=1,nlon
          rho(i,k) = ypm(i,j,k) / (rair*(1.+zvir*yq(i,j,k))*yt(i,j,k))
          dm(i,k) = ypthic(i,j,k) / gravit
          dz(i,k) = dm(i,k) / rho(i,k)
  102   continue
  100 continue
 
c        Compute rho at half levels, and vertical thicknesses between
c        layer centers. Doesn't matter what these values are at top and
c        bottom (k=1 and nlevp), since vertk will be zero there.
 
      do 110 k=1,nlevp
        km1 = max (k-1,1)
        kk  = min (k,nlev)
        do 112 i=1,nlon
          rho_h(i,k) = 0.5*(rho(i,km1)+rho(i,kk))
          dz_h(i,k) = 0.5*(dz(i,km1) + dz(i,kk))
  112   continue
  110 continue
 
c-------------------
 
c        Calculate vertical diffusion coeffs at eta-grid half levels,
c        using u,v,t,q fields at time n. First set potential temperature
c        at full levels.
 
      do 150 k=1,nlev
        do 152 i=1,nlon
          zthet(i,k) = yt(i,j,k)*(1.+zvir*yq(i,j,k)) / ypmcap(i,j,k)
  152   continue
  150 continue
 
      call vertcoef (vertk, yu, yv, zthet, dz_h, j)
 
c        Save Brunt-Vaisala frequency in global array
c
c     do 160 k=2,nlev-1
c       do 162 i=1,nlon
c         vertoz(i,j,k) =
c    *      sqrt ( gravit  * max (zthet(i,k-1)-zthet(i,k+1),0.)
c    *                   / ( 0.5*(zthet(i,k-1)+zthet(i,k+1))
c    *                          *(dz_h (i,k)  +dz_h (i,k+1)) )
c    *           )
c 162   continue
c 160 continue
c     do 164 i=1,nlon
c         vertoz(i,j,1) =
c    *      sqrt ( gravit  * max (zthet(i,1)-zthet(i,2),0.)
c    *                   / ( 0.5*(zthet(i,1)+zthet(i,2))
c    *                          * dz_h(i,2) )
c    *           )
c         vertoz(i,j,nlev) =
c    *      sqrt ( gravit  * max (zthet(i,nlev-1)-zthet(i,nlev),0.)
c    *                   / ( 0.5*(zthet(i,nlev-1)+zthet(i,nlev))
c    *                          * dz_h(i,nlev) )
c    *           )
c 164 continue
 
c-------------------
 
c        Extract some of q-reservoir fplumq into flux convergence dplumq
 
      call reserv  (fplumq, dplumq, dtime, j)
#ifdef isotrac
      call reserv  (fplumq_o, dplumq_o, dtime, j)
      call reserv  (fplumq_h, dplumq_h, dtime, j)
#endif
#ifdef isotracreg
      do it=1,nrtr
        call reserv  (rtrfq(1,1,1,it), dplumq_t(1,1,it), dtime, j)
      enddo
#endif
 
c       Calculate vertical diffusion for q
 
      call vertdifn (qnew, yqm1, dplumq,
     *               dm, dz_h, rho_h, vertk, conev, conel,
     *               twodt, rimq, j)
#ifdef isotrac
      call vertdifn (qnew_o, yqm1_o, dplumq_o,
     *               dm, dz_h, rho_h, vertk, conev, conel,
     *               twodt, rimq, j)
      call vertdifn (qnew_h, yqm1_h, dplumq_h,
     *               dm, dz_h, rho_h, vertk, conev, conel,
     *               twodt, rimq, j)
#endif
#ifdef isotracreg
      do it=1,nrtr
        call vertdifn (qnew_t(1,1,it),rtrqm1(1,1,1,it),dplumq_t(1,1,it),
     *                 dm, dz_h, rho_h, vertk, conev, conel,
     *                 twodt, rimq, j)
      enddo
#endif

c        Fix negative q, since vertfdifn can produce o(10**-4)
c        negative values, and accumulate budget error quantity dwneg.
c        Store resulting q values in global array yqm1.
 
c-----
c     Simple setting back to zero:
c     zneg = 0.
c     do 430 k=1,nlev
c       zwei = 1./(nlon*gravit)
c       do 432 i=1,nlon
c         zneg = zneg + ypthic(i,j,k)*max(-qnew(i,k), 0.)
c         yqm1(i,j,k) = max (qnew(i,k), 0.)
c 432   continue
c 430 continue
c     dwneg(j) = dwneg(j) + 0.5*zneg/(nlon*gravit)
c-----
c     As in plume for cloud liquid water (dwneg stays zero):
      call zero (ztot1, nlon)
      call zero (ztot2, nlon)
c=============
#ifdef isotrac
c=============
      call zero (ztot1_o, nlon)
      call zero (ztot1_h, nlon)
      call zero (ztot2_o, nlon)
      call zero (ztot2_h, nlon)
#ifdef isotracreg
      call zero (ztot1_t, nlon*nrtr)
      call zero (ztot2_t, nlon*nrtr)
#endif

c     Search for tracer ratios min/max for each vertical column, taking
c     into account the change associated with dplumq in vertdifn. Fixing
c     tracer ratios extrema is done below to avoid most of the crazy 
c     ratios that are sometimes generated (often at the tropopause) in
c     vertdifn.
      do i=1,nlon
        rmax_o(i) = -1.e20
        rmax_h(i) = -1.e20
        rmin_o(i) =  1.e20
        rmin_h(i) =  1.e20
          do k=1,nlev
            rq(i,k)= yqm1(i,j,k)  +dplumq(i,k)  *twodt/dm(i,k)
            rq_o   = yqm1_o(i,j,k)+dplumq_o(i,k)*twodt/dm(i,k)
            rq_h   = yqm1_h(i,j,k)+dplumq_h(i,k)*twodt/dm(i,k)
            rmax_o(i) = cvmgt(max(rmax_o(i), rq_o/rq(i,k)),
     *                            rmax_o(i), rq(i,k).ne.0.)
            rmax_h(i) = cvmgt(max(rmax_h(i), rq_h/rq(i,k)),
     *                            rmax_h(i), rq(i,k).ne.0.)
            rmin_o(i) = cvmgt(min(rmin_o(i), rq_o/rq(i,k)),
     *                            rmin_o(i), rq(i,k).ne.0.)
            rmin_h(i) = cvmgt(min(rmin_h(i), rq_h/rq(i,k)),
     *                            rmin_h(i), rq(i,k).ne.0.)
          enddo
      enddo
#ifdef isotracreg
      do it=1,nrtr
        do i=1,nlon
          rmax_t(i,it) = -1.e20
          rmin_t(i,it) =  1.e20
          do k=1,nlev
            rq_t = rtrqm1(i,j,k,it)+dplumq_t(i,k,it)*twodt/dm(i,k)
            rmax_t(i,it) = cvmgt(max(rmax_t(i,it), rq_t/rq(i,k)),
     *                               rmax_t(i,it), rq(i,k).ne.0.)
            rmin_t(i,it) = cvmgt(min(rmin_t(i,it), rq_t/rq(i,k)),
     *                               rmin_t(i,it), rq(i,k).ne.0.)
          enddo
        enddo
      enddo
#endif
c=====
#endif
c=====

      do 430 k=1,nlev
        do 432 i=1,nlon
          zm = ypthic(i,j,k)/gravit
          ztot1(i) = ztot1(i) + qnew(i,k)*zm
          qnew(i,k) = max (qnew(i,k), 0.)
          ztot2(i) = ztot2(i) + qnew(i,k)*zm
c~~~~~~~~~~~~~
#ifdef isotrac
c~~~~~~~~~~~~~
c         Fix negative q for tracers and fix tracer ratio
c         extrema violations:
          ztot1_o(i) = ztot1_o(i) + qnew_o(i,k)*zm
          ztot1_h(i) = ztot1_h(i) + qnew_h(i,k)*zm
          qnew_o(i,k) = max (qnew_o(i,k), 0.)
          qnew_h(i,k) = max (qnew_h(i,k), 0.)
          if (yqm1(i,j,k).ne.0..and.qnew(i,k).ne.0..and.
     *       (qnew_o(i,k).eq.0..or.qnew_h(i,k).eq.0.)) then
            qnew_o(i,k) = qnew(i,k)*yqm1_o(i,j,k)/yqm1(i,j,k)
            qnew_h(i,k) = qnew(i,k)*yqm1_h(i,j,k)/yqm1(i,j,k)
          endif
          if (qnew(i,k).eq.0.) then
            qnew_o(i,k) = 0.
            qnew_h(i,k) = 0.
          else
            qnew_o(i,k) = max(min(qnew_o(i,k),
     *                    rmax_o(i)*qnew(i,k)), rmin_o(i)*qnew(i,k))
            qnew_h(i,k) = max(min(qnew_h(i,k),
     *                    rmax_h(i)*qnew(i,k)), rmin_h(i)*qnew(i,k))
          endif
          ztot2_o(i) = ztot2_o(i) + qnew_o(i,k)*zm
          ztot2_h(i) = ztot2_h(i) + qnew_h(i,k)*zm
#ifdef isotracreg
          do it=1,nrtr
            ztot1_t(i,it) = ztot1_t(i,it) + qnew_t(i,k,it)*zm
            qnew_t(i,k,it) = max (qnew_t(i,k,it), 0.)
            if (qnew(i,k).eq.0.) then
              qnew_t(i,k,it) = 0.
            else
              qnew_t(i,k,it) = max(min(qnew_t(i,k,it),
     *                rmax_t(i,it)*qnew(i,k)), rmin_t(i,it)*qnew(i,k))
            endif
            ztot2_t(i,it) = ztot2_t(i,it) + qnew_t(i,k,it)*zm
          enddo
#endif
c~~~~~
#endif
c~~~~~
  432   continue
  430 continue
 
      do 434 i=1,nlon
c       Since ztot2 > ztot1, if ztot1>0, so is ztot2:
        ztot1(i) = cvmgt ( ztot1(i)/max(ztot2(i),1.e-20), 0.,
     *                     ztot1(i).gt.0. )
c+++++++++++++
#ifdef isotrac
c+++++++++++++
c       Determine stratospheric water vapor tracer mass and tropopause
c       k-level. The "tropopause" is defined by the level where the 
c       mass of water in the atmosphere above it is a 1/1000 of the 
c       total mass of water in the air column.
        ztot     = 0.
        ztot2s_o = 0.
        ztot2s_h = 0.
#ifdef isotracreg
        call zero (ztot2s_t, nrtr)
#endif
        do k=1,nlev
          ztot     = ztot     + qnew(i,k)  *dm(i,k)
          if (ztot/ztot2(i).gt.0.001) goto 435
          ztot2s_o = ztot2s_o + qnew_o(i,k)*dm(i,k)
          ztot2s_h = ztot2s_h + qnew_h(i,k)*dm(i,k)
#ifdef isotracreg
          do it=1,nrtr
            ztot2s_t(it) = ztot2s_t(it) + qnew_t(i,k,it)*dm(i,k)
          enddo
#endif
        enddo
 435    ktrop(i) = k-1
        ztot1_o(i) = cvmgt ( (ztot1_o(i)-ztot1(i)*ztot2s_o)/
     *                   max((ztot2_o(i)-ztot2s_o), 1.e-20), 0.,
     *                       (ztot1_o(i)-ztot1(i)*ztot2s_o).gt.0. )
        ztot1_h(i) = cvmgt ( (ztot1_h(i)-ztot1(i)*ztot2s_h)/
     *                   max((ztot2_h(i)-ztot2s_h), 1.e-20), 0.,
     *                       (ztot1_h(i)-ztot1(i)*ztot2s_h).gt.0. )
#ifdef isotracreg
        do it=1,nrtr
         ztot1_t(i,it)=cvmgt((ztot1_t(i,it)-ztot1(i)*ztot2s_t(it))/
     *                   max((ztot2_t(i,it)-ztot2s_t(it)), 1.e-20), 0.,
     *                     (ztot1_t(i,it)-ztot1(i)*ztot2s_t(it)).gt.0.)
        enddo
#endif
c+++++
#endif
c+++++
  434 continue
 
      do 436 k=1,nlev
        do 438 i=1,nlon
          qnew(i,k) = ztot1(i) * qnew(i,k)
          yqm1(i,j,k) = qnew(i,k)
c.............
#ifdef isotrac
c.............
c         Fix mass, conserving tracer ratios in the tropopause:
          if (k.le.ktrop(i)) then
            qnew_o(i,k) = ztot1(i) * qnew_o(i,k)
            qnew_h(i,k) = ztot1(i) * qnew_h(i,k)
#ifdef isotracreg
            do it=1,nrtr
              qnew_t(i,k,it)   = ztot1(i) * qnew_t(i,k,it)
              rtrqm1(i,j,k,it) = qnew_t(i,k,it)
            enddo
#endif
          else
            qnew_o(i,k) = ztot1_o(i) * qnew_o(i,k)
            qnew_h(i,k) = ztot1_h(i) * qnew_h(i,k)
#ifdef isotracreg
            do it=1,nrtr
              qnew_t(i,k,it)   = ztot1_t(i,it) * qnew_t(i,k,it)
              rtrqm1(i,j,k,it) = qnew_t(i,k,it)
            enddo
#endif
          endif
          yqm1_o(i,j,k) = qnew_o(i,k)
          yqm1_h(i,j,k) = qnew_h(i,k)
c.....
#endif
c.....
  438   continue
  436 continue
 
c-------------------
 
c        Calculate vertical diffusion for tracers (non-leapfrog step).
c        Don't do for flux reservoirs (equivalenced to tracers #4-7).
c        Don't do for passive tracers (>= #8) if turned off
c        by namelist travdif (default = .true.)
 
      do 600 n=1,ntrace
 
        if (n.le.3 .or. (n.ge.ntraca+1 .and. travdif)) then
 
c           Calculate vert diffusion for tracers (1-dt non-leapfrg step)
 
          call vertdifn (anew, tracer(1,1,1,n), vzero,
     *                   dm, dz_h, rho_h, vertk, conev, conel,
     *                   dtime, rimq, j)
#ifdef isotrac
          if (n.eq.1) then
            call vertdifn (anew_o, ylwc_c_o, vzero,
     *                     dm, dz_h, rho_h, vertk, conev, conel,
     *                     dtime, rimq, j)
            call vertdifn (anew_h, ylwc_c_h, vzero,
     *                     dm, dz_h, rho_h, vertk, conev, conel,
     *                     dtime, rimq, j)
          else if (n.eq.2) then
            call vertdifn (anew_o, ylwc_a_o, vzero,
     *                     dm, dz_h, rho_h, vertk, conev, conel,
     *                     dtime, rimq, j)
            call vertdifn (anew_h, ylwc_a_h, vzero,
     *                     dm, dz_h, rho_h, vertk, conev, conel,
     *                     dtime, rimq, j)
          else if (n.eq.3) then
            call vertdifn (anew_o, ylwc_s_o, vzero,
     *                     dm, dz_h, rho_h, vertk, conev, conel,
     *                     dtime, rimq, j)
            call vertdifn (anew_h, ylwc_s_h, vzero,
     *                     dm, dz_h, rho_h, vertk, conev, conel,
     *                     dtime, rimq, j)
          endif
#endif
 
c           Set new global tracer field, fixing negative values
 
          do 620 k=1,nlev
            do 622 i=1,nlon
#ifdef isotrac
              if (n.le.3) then
                anew_o(i,k) = max (anew_o(i,k), 0.)
                anew_h(i,k) = max (anew_h(i,k), 0.)
                if (anew(i,k).le.0.) then
                  anew_o(i,k) = 0.
                  anew_h(i,k) = 0.
                endif
                if ((anew(i,k).gt.0.).and.(tracer(i,j,k,n).gt.0.).and.
     *              ((anew_o(i,k).eq.0.).or.(anew_h(i,k).eq.0.))) then
                  if (n.eq.1) then
                   anew_o(i,k)=anew(i,k)*ylwc_c_o(i,j,k)/tracer(i,j,k,n)
                   anew_h(i,k)=anew(i,k)*ylwc_c_h(i,j,k)/tracer(i,j,k,n)
                  else if (n.eq.2) then
                   anew_o(i,k)=anew(i,k)*ylwc_a_o(i,j,k)/tracer(i,j,k,n)
                   anew_h(i,k)=anew(i,k)*ylwc_a_h(i,j,k)/tracer(i,j,k,n)
                  else if (n.eq.3) then
                   anew_o(i,k)=anew(i,k)*ylwc_s_o(i,j,k)/tracer(i,j,k,n)
                   anew_h(i,k)=anew(i,k)*ylwc_s_h(i,j,k)/tracer(i,j,k,n)
                  endif
                endif
                if (n.eq.1) then
                  ylwc_c_o(i,j,k) = anew_o(i,k)
                  ylwc_c_h(i,j,k) = anew_h(i,k)
                else if (n.eq.2) then
                  ylwc_a_o(i,j,k) = anew_o(i,k)
                  ylwc_a_h(i,j,k) = anew_h(i,k)
                else if (n.eq.3) then
                  ylwc_s_o(i,j,k) = anew_o(i,k)
                  ylwc_s_h(i,j,k) = anew_h(i,k)
                endif
              endif
#endif

c             tracer(i,j,k,n) = anew(i,k)
c             Or crudely fix negative values:
              tracer(i,j,k,n) = max (anew(i,k), 0.)
  622       continue
  620     continue
 
        endif
 
  600 continue
 
      return
      end
c
c****************************** VERTCOEF *******************************
c
      subroutine vertcoef (vertk, u, v, zthet, dz_h, j)
 
c     Calculates vertical-diffusion coefficient on half sigma or eta
c     levels, one latitude band at a time.
 
c     Returned:
c     vertk = diffusion coefficient on half levels (lon-lev slice, m2/s)
 
c     Supplied:
c     u    = eastward velocity at full levels (3-D, m/s)
c     v    = northward velocity at full levels (3-D, m/s)
c     zthet= potential temperature at full levels (lon-lev slice, K)
c     dz_h = vertical distance between half levels (lon-lev slice, m)
c     j    = latitude index.
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMMAP.h>
c ---------------------------------------------------------------------
      dimension
     *  vertk(nlon,nlevp),
     *  u(nlon,norec,nlev),     v(nlon,norec,nlev),
     *  zthet(nlon,nlev),       dz_h(nlon,nlevp)
 
      dimension
     *  zuv(nlon,nlevp),        zrich(nlon,nlevp),
     *  zneut(nlon,nlevp)
c ----------------------------------------------------------------------
      parameter (vertkmin=0.1, vertkmax=1.e10, vertlen=30.)
c ----------------------------------------------------------------------
 
c        Calculate zuv = du**2 + dv**2, zrich = Richardson Number,
c        zneut = neutral difn coeff, and vertk, all at half levels.
c        Impose lower limit of 1 m2/s2 on zuv to allow motionless
c        cold starts.
 
      do 100 k=2,nlev
        do 102 i=1,nlon
          zuv(i,k) =  max (  (u(i,j,k-1)-u(i,j,k))**2
     *                     + (v(i,j,k-1)-v(i,j,k))**2, 1. )
 
          zrich(i,k) = gravit * (zthet(i,k-1)-zthet(i,k)) * dz_h(i,k)
     *                  / ( 0.5*(zthet(i,k-1)+zthet(i,k)) * zuv(i,k) )
 
          zneut(i,k) = (vertlen**2) * sqrt(zuv(i,k)) / dz_h(i,k)
 
          vertk(i,k) = zneut(i,k)
     *               * cvmgt ( max(1.-5.*zrich(i,k), 0.),
     *                         sqrt(max(1.-18.*zrich(i,k), 0.)),
     *                         zrich(i,k).ge.0. )
 
          vertk(i,k) = min (vertkmax, max (vertkmin, vertk(i,k)))
 
  102  continue
  100 continue
 
c        Zero  diffusion coefficient at top and bottom boundaries
 
      call zero (vertk(1,1), nlon)
      call zero (vertk(1,nlevp), nlon)
 
      return
      end
c
c****************************** VERTDIFN *******************************
c
      subroutine vertdifn (fnew, field, forc,
     *                     dm, dz_h, rho_h, vertk, wei, weil,
     *                     dt, rimp, j)
 
c     Steps fields due to vertical diffusion, one latitude band at a
c     time. Grid can be sigma or hybrid. Timestep is leapfrog, from n-1
c     to n+1. Calls TRIDIA for time-implicit linear calculation, with
c     time-explicit diffusion coeffs calculated at time n.
 
c     Returned:
c     fnew  = new value of field after timestep (this latitude slice)
 
c     Supplied:
c     field = field(u,v,t on sigma grid, or q,tracer on hybrid grid)
c     forc  = forcing (flux convergence per layer)
c     dm    = layer masses (kg m-2), or cp*mass for T
c     dz_h  = vertical distance between half levels (m)
c     rho_h = density at half levels (kg m-3), or cp*density for T
c     vertk = vertical diffusion coefficient at half levels (m2/s)
c     wei   = sigma**kappa at full levels for T, otherwise = 1
c     weil  = sigma**kappa at boundaries for T, otherwise = 1
c     dt    = timestep (s)
c     rimp  = implicit fraction of the tridiagonal calculation (0 to 1)
c             (=1 for u,v,t, =~.5 for q,tracers to avoid -ve values)
c     j     = latitude index
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMMAP.h>
c ---------------------------------------------------------------------
      dimension
     *  fnew(nlon,nlev),
     *  field(nlon,norec,nlev),    forc(nlon,nlev),
 
     *  dm(nlon,nlev),             dz_h(nlon,nlevp),
     *  rho_h(nlon,nlevp),         vertk(nlon,nlevp),
     *  wei(nlev),                 weil(nlevp)
 
      dimension
     *  zdt(nlon,nlev),     zdif(nlon,nlevp),
     *  zd1(nlon,nlev),     zd2(nlon,nlev),    zd3(nlon,nlev),
     *  zrhs(nlon,nlev),    zwork1(nlon,nlev), zwork2(nlon,nlev)
c ---------------------------------------------------------------------
 
c          Set full-level dt/(dz*rho)
 
        do 100 k=1,nlev
          do 102 i=1,nlon
            zdt(i,k)= dt / dm(i,k)
  102     continue
  100   continue
 
c          Set half-level K*rho/dz (vertk=zdif=0 at k=1 and nlev+1)
 
        do 200 k=1,nlevp
          do 202 i=1,nlon
            zdif(i,k)= weil(k)*vertk(i,k)*rho_h(i,k)/dz_h(i,k)
  202     continue
  200   continue
 
c        Set matrix diagonals
 
      do 300 k=1,nlev
        weim1 = wei(max(k-1,1))
        weip1 = wei(min(k+1,nlev))
        do 302 i=1,nlon
          zd1(i,k) =    - zdt(i,k)*rimp* zdif(i,k)/weim1
          zd2(i,k) = 1. + zdt(i,k)*rimp*(zdif(i,k)+zdif(i,k+1))/wei(k)
          zd3(i,k) =    - zdt(i,k)*rimp* zdif(i,k+1)/weip1
  302   continue
  300 continue
 
c       Set right-hand sides
 
      if (rimp.eq.1.) then
        do 400 k=1,nlev
          km1 = max (k-1, 1)
          kp1 = min (k+1, nlev)
          do 402 i=1,nlon
            zrhs(i,k) = field(i,j,k) + zdt(i,k)*forc(i,k)
  402     continue
  400   continue
      else
        do 410 k=1,nlev
          km1 = max (k-1, 1)
          kp1 = min (k+1, nlev)
          weim1 = wei(km1)
          weip1 = wei(kp1)
          do 412 i=1,nlon
            zrhs(i,k) = field(i,j,k) + zdt(i,k)*forc(i,k)
     *       + zdt(i,k) * (1.-rimp)
     *       * ( zdif(i,k)  *(field(i,j,km1)/weim1-field(i,j,k)/wei(k))
     *          +zdif(i,k+1)*(field(i,j,kp1)/weip1-field(i,j,k)/wei(k)))
  412     continue
  410   continue
      endif
 
c        Solve the tridiagonal systems
 
      call tridia (nlon,nlon,nlev, zd1,zd2,zd3, zrhs, fnew,
     *             zwork1,zwork2)
 
      return
      end
C
C ***************************** RESERV ********************************
C
      subroutine reserv (freserv, falloc, dtres, j)
 
c       Computes flux falloc from reservoir freserv (one of
c       fplum[t,q,u,v] in commun, which are in units of flux*time per
c       layer), to be applied this timestep in vdif, and
c       decrements freserv by falloc*dtres. The physical timescale for
c       release from all reservoirs is set by namelist parameter
c       timreserv (s). These reservoirs are all on hybrid (eta) coords.
 
c       Called one latitude band at a time from vdif, so falloc does
c       not have to be global.
 
c          freserv = "flux_convergence*time" reservoir (modified)
c          falloc  = "flux_convergence" (returned)
c          dtres   = reservoir timestep for this call (always = dtime)
c          j       = latitude band
 
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
      dimension freserv(nlon,norec,nlev), falloc(nlon,nlev)
C----------------------------------------------------------------------
 
      zfrac = min (1., dtres / max(timreserv,.001))
 
      do 100 k=1,nlev
        do 102 i=1,nlon
          falloc(i,k)  = zfrac * freserv(i,j,k)
          freserv(i,j,k) = freserv(i,j,k) - falloc(i,k)
          falloc(i,k)  = falloc(i,k) / dtres
  102   continue
  100 continue
 
      return
      end
c
c****************************** DIFGRID ********************************
c
c     subroutine difgrid (field, coef, dt)
c
cc       Horizontally diffuses (del**2) each level of 3-D array field
cc       on agcm grid through time dt with difn coeff coef (m s-2).
cc       Used only for fplumt, called at start of vdif_sig.
c
cc--------------------------------------------------------------------
c#include <PARAMS.h>
cc---------------------------------------------------------------------
c#include <COMCON.h>
cc--------------------------------------------------------------------
c#include <COMMAP.h>
cc---------------------------------------------------------------------
c#include <COMTIM.h>
cc--------------------------------------------------------------------
c     dimension field(nlon,norec,nlev), ztmp(nlon,norec)
c     parameter (pi = 3.14159265358979)
cc--------------------------------------------------------------------
c
c     dlon = 2.*pi/nlon
c
c     do 100 k=1,nlev
c
cc      Diagnostic only:
cc      if (nstep.eq.nstop) then
cc        zglob1 = 0.
cc        do 800 j=1,norec
cc          zwei = 0.5*w(j)/nlon
cc          do 802 i=1,nlon
cc            zglob1 = zglob1 + field(i,j,k)*zwei
cc802       continue
cc800     continue
cc      endif
c
c       do 102 j=1,norec
c         jp = min (j+1,norec)
c         jm = max (j-1,1)
c
c         clat  = cos(alatg(j))
c         clatm = cos(0.5*(alatg(j)+alatg(jm)))
c         clatp = cos(0.5*(alatg(j)+alatg(jp)))
c
c         dlm = alatg(j)-alatg(jm)
c         if (j.eq.1) dlm = alatg(j) + 0.5*pi
c         dlp = alatg(jp)-alatg(j)
c         if (j.eq.norec) dlp = 0.5*pi - alatg(j)
c
c         do 104 i=1,nlon
c           ip = mod(i+1,nlon)
c           im = mod(i-2+nlon,nlon) + 1
c           ztmp(i,j) = (field(im,j,k)+field(ip,j,k)-2.*field(i,j,k))
c    *                  / (clat**2 * dlon**2)
c    *                + (   clatp*(field(i,jp,k)-field(i,j,k))/dlp
c    *                    + clatm*(field(i,jm,k)-field(i,j,k))/dlm )
c    *                  / w(j)
c 104     continue
c 102   continue
c
c       do 150 j=1,norec
c         do 152 i=1,nlon
c           field(i,j,k) = field(i,j,k) + ztmp(i,j)*dt*coef/rearth**2
c 152     continue
c 150   continue
c
cc      Diagnostic only:
cc      if (nstep.eq.nstop) then
cc        zglob2 = 0.
cc        do 820 j=1,norec
cc          zwei = 0.5*w(j)/nlon
cc          do 822 i=1,nlon
cc            zglob2 = zglob2 + field(i,j,k)*zwei
cc822       continue
cc820     continue
cc        write(*,830) coef, dt, k, zglob1, zglob2
cc830     format('difgrid: coef, dt, k, zglob1/2 =',2e15.5,i4,2e20.10)
cc      endif
c
c 100 continue
c
c     return
c     end
c
c****************************** GWD ************************************
c
      subroutine gwd (yu, yv, yt, yps, yfu, yfv, yt2, j)
 
c Apply Mcfarlane (1987,JAS,44,1775) stationary orographic gravity wave
c drag scheme. Adapted from ccm2 subr gwd. Does one lon-lev slice at a
c time. The prescribed s.d. of topography is sdgravxy (in comgwd), read
c in subr basic. Tuned by uniform-factor namelist param gravfac (in
c comflg)...do nothing if gravfac=0 or gwrough=0.
 
c yu     = eastward velocity, m/s (supplied)
c yv     = northward velocity, m/s (supplied)
c yt     = temperature, K (supplied)
c yps    = surface pressure, N/m2 (supplied)
c yfu    = layer eastward momentum tendency, m/s2 (modified)
c yfv    = layer northward momentum tendency, m/s2 (modified)
c yt2    = layer heating rate, W/m2 (modified)
c j      = current latitude index (supplied)
 
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMGWD.h>
C----------------------------------------------------------------------
#include <COMINI.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      dimension
     *   yu(nlon,norec,nlev),        yv(nlon,norec,nlev),
     *   yt(nlon,norec,nlev),        yps(nlon,norec),
     *   yfu(nlon,norec,nlev),       yfv(nlon,norec,nlev),
     *   yt2(nlon,norec,nlev)
c----------------------------------------------------------------------
c Local arrays:
      dimension
     *   sdgrav(nlon),
     *   pmid(nlon,nlev),            pint(nlon,nlevp),
     *   ub(nlon),                   vb(nlon),
     *   vmb(nlon),                  vp(nlon,nlev),
     *   bv(nlon,nlevp),             ui(nlon,nlev),
     *   ti(nlon,nlev),              gtau(nlon,nlevp),
     *   indx(nlon),                 ztmp(nlon)
 
c Local arrays vs. latitude, saved for output:
      dimension
     *   zduout(norec,nlev),         zdvout(norec,nlev)
      save zduout, zdvout
 
      logical first(norec)
      save first
      data first /norec*.true./
c-----------------------------------------------------------------------
      parameter ( vmin  = 3.0,     ! minimum reference-level wind
     *            v0    = 1.e-5,   ! approx zero wind
     *            hmin  = 5.0,     ! minimum subgrid s.d. topog
     *            fc2   = 0.5,     ! critical Froude Number squared
     *            fudge = 8.0e-6,  ! tunable param e*mu/2 = e*pi/l
     *            tnlim = 2.9e-3,  ! max |v| tendency, ~250 m/s/day
     *            bv2min= 2.5e-5 ) ! min value of Brunt-Vaisala freq**2
c-----------------------------------------------------------------------
 
c        Do nothing if gravfac=0 or gwrough=0
 
      if (gravfac.le.0. .or. gwrough.eq.0.) return
 
c        Initialize diagnostic accumulators
 
      if (first(j)) then
         do 2 k=1,nlev
           zduout(j,k) = 0.
           zdvout(j,k) = 0.
    2   continue
        first(j) = .false.
      endif
 
      g2or = -gravit*gravit/rair
      g2ocp = gravit*gravit/cpair
      fgor  = gravfac*fudge*gravit/rair
 
c        If requested, set sdgrav to single value over land/ice pts,
c        and to zero over ocean pts (can't set in topini since yoro's
c        not set until inisurf called)
 
      if (gwrough.eq.-999.) then
c       set from data file:
        do i=1,nlon
          sdgrav(i) = sdgravxy(i,j)
        enddo
      else if (gwrough.lt.0.) then
c       set over ice sheet only, else from data file:
        do i=1,nlon
          sdgrav(i) = (1.-yoroice(i,j))*sdgravxy(i,j)
     *              +     yoroice(i,j) *(-gwrough)
        enddo
      else 
c       set uniformly over all land/ice sheet, zero over ocean:
        do i=1,nlon
          sdgrav(i) = (1.-yoroocn(i,j))*gwrough
        enddo
      endif
 
c        Calculate reference-level wind speeds
 
      do 10 i=1,nlon
         vmb(i) = max ( sqrt(yu(i,j,nlev)**2+yv(i,j,nlev)**2), 1. )
   10 continue
 
c        Determine the set of points at which drag will be applied:
c        lowest level wind > vmin, and subgrid height sdgrav > hmin
 
      lengw = 0
      do 20 i=1,nlon
         if (vmb(i).gt.vmin .and. sdgrav(i).gt.hmin) then
            lengw = lengw + 1
            indx(lengw) = i
         end if
   20 continue
 
c        Don't do any calculations if no longitudes in this strip
 
c-------------------------
      if (lengw.ne.0) then
c-------------------------
 
c           Set mid-layer and interface pressures (yu,yv,yt are on sigma
c           coords)
 
         do 100 k=1,nlev
           do 102 ii=1,lengw
              i=indx(ii)
              pmid(i,k) = sig(k)*yps(i,j)
              pint(i,k) = sigkmh(k)*yps(i,j)
  102      continue
  100    continue
 
c           Reset top interface pressure assuming equal increment in
c           ln p between the two interfaces below. Also set bottom
c           interface pressure, and lowest midlevel wind unit vectors
 
         do 110 ii=1,lengw
            i=indx(ii)
            pint(i,1) = exp(2.0*log(pint(i,2)) - log(pint(i,3)))
            pint(i,nlevp) = sigkmh(nlevp)*yps(i,j)
            ub(i) = yu(i,j,nlev)/vmb(i)
            vb(i) = yv(i,j,nlev)/vmb(i)
  110    continue
 
c           Calculate projection of wind in reference level direction.
c           No negative or zero values permitted.
 
         do 120 k=1,nlev
            do 122 ii=1,lengw
               i=indx(ii)
               vp(i,k) = max ( yu(i,j,k)*ub(i) + yv(i,j,k)*vb(i), v0 )
  122       continue
  120    continue
 
         do 200 k=2,nlev
!DIR$ IVDEP
            do 202 ii=1,lengw
               i=indx(ii)
 
c           Compute interior interface temperatures and wind projections
 
               zfac = log(pint(i,k)/pmid(i,k-1))
     *              / log(pmid(i,k)/pmid(i,k-1))
               ti(i,k) = yt(i,j,k-1) + zfac*(yt(i,j,k)-yt(i,j,k-1))
               ui(i,k) = vp(i,  k-1) + zfac*(vp(i,  k)-vp(i,  k-1))
 
c           Compute interior interface Brunt-Vaisala frequencies
 
               bv(i,k) = sqrt ( max ( bv2min,
     *                 g2or*pint(i,k)/(ti(i,k)**2)*
     *                 ( (yt(i,j,k)-yt(i,j,k-1))/(pmid(i,k)-pmid(i,k-1))
     *                   - cappa*ti(i,k)/pint(i,k) ) ) )
  202       continue
  200    continue
 
c           Compute top interface Brunt-Vaisala frequency assuming
c           isothermal profile above top mid-layer. Also calculate
c           lowest level g*stress, using Brunt-Vaisala frequency for
c           first interior interface (top of first layer).
 
         do 240 ii=1,lengw
            i=indx(ii)
            bv(i,1) = sqrt(g2ocp/yt(i,j,1))
            gtau(i,nlevp) = fgor * min ( 4.*sdgrav(i)**2,
     *                                   fc2*(vmb(i)/bv(i,nlev))**2 )
     *                    * pmid(i,nlev)*bv(i,nlev)*vmb(i)/yt(i,j,nlev)
  240    continue
 
c           Calc g*stress at all other interior levels from bottom up
 
         do 290 k=nlev,2,-1
           zw = min (1., max (0., (0.7-sig(k))/0.3 ))
           do 292 ii=1,lengw
             i=indx(ii)
             gmax = fgor*fc2*(ui(i,k)**3)*pint(i,k)/(bv(i,k)*ti(i,k))
             gtau(i,k) = min (gtau(i,k+1), zw*gmax+(1.-zw)*gtau(i,k+1))
  292      continue
  290    continue
 
c           Calculate g*stress at top interface, assuming it is located
c           an equal increment in ln p above the 2 interfaces below,
c           and all variables are constant above top midpoint.
 
         do 300 ii=1,lengw
            i=indx(ii)
            gtau(i,1) = min ( gtau(i,2), fgor*fc2*vp(i,1)**3
     *                                  *pint(i,1)/(bv(i,1)*yt(i,j,1)) )
  300    continue
 
c           Increment 3-D momentum tendencies yfu,yfv by
c           [g*stress divergence / delta_p]. Bound the gravity-wave
c           increment by a prescribed limit. Also add momentum diffusive
c           heating to net heating rate yt2 (W/m2/layer) and increment
c           vertical diffusion ke dissipation budget term dkevdf
c           (nb: u and v are on sigma coords so use dsigma)
 
         do 320 k=1,nlev
            call zero (ztmp,nlon)
!DIR$ IVDEP
            do 322 ii=1,lengw
               i=indx(ii)
               utgw = -ub(i) * (gtau(i,k+1)-gtau(i,k))
     *                       / (pint(i,k+1)-pint(i,k))
               vtgw = -vb(i) * (gtau(i,k+1)-gtau(i,k))
     *                       / (pint(i,k+1)-pint(i,k))
               tausx = max ( sqrt(utgw**2+vtgw**2), 1.e-15 )
               utgw = utgw*min(tnlim,tausx)/tausx
               vtgw = vtgw*min(tnlim,tausx)/tausx
 
               yfu(i,j,k) = yfu(i,j,k) + utgw
               yfv(i,j,k) = yfv(i,j,k) + vtgw
 
               ztmp(i) = -(yu(i,j,k)*utgw + yv(i,j,k)*vtgw)
     *                   *dsigma(k)*yps(i,j)/gravit
               yt2(i,j,k) = yt2(i,j,k) + ztmp(i)
 
c              Diagnostics:
               zduout(j,k) = zduout(j,k) + utgw
               zdvout(j,k) = zdvout(j,k) + vtgw
c              if (i.eq.nint(nlon*(360.-110.)/360.)) then
c                zduout(j,k) = utgw
c                zdvout(j,k) = vtgw
c                zuout(j,k) = yu(i,j,k)
c                zvout(j,k) = yv(i,j,k)
c                zbout(j,k) = bv(i,k)
c              endif
 
  322       continue
            dkevdf(j) = dkevdf(j) + ssum(nlon,ztmp,1)*twodt*0.5/nlon
  320    continue
 
c----------
      endif
c----------
 
      return
 
c----------------------------------------------------------------------
c        Separate entry point for printout section (called after gwd
c        latitude loop in phys, to allow multiprocessing)
 
      entry gwdo (iu)
c----------------------------------------------------------------------
 
c        Do nothing if gravfac=0 or gwrough=0
 
      if (gravfac.le.0. .or. gwrough.eq.0.) return
 
      if (nstep.eq.nstop) then
 
        incz = (norec-1)/40 + 1
        nz = (nstop-nrstrt+1)*nlon
 
        write(iu,7000) char(12),
     *    'gravity wave drag: u tendency (10*m/s/day)',
     *    (nint(alatg(jj)*180./3.14159),jj=norec,1,-incz)
        write(iu,*)
        do 7100 k=1,nlev
          write(iu,7002)
     *      sig(k),(nint(86400.*10.*zduout(jj,k)/nz), jj=norec,1,-incz)
 7100   continue
 
        write(iu,7000) ' ',
     *    'gravity wave drag: v tendency (10*m/s/day)',
     *    (nint(alatg(jj)*180./3.14159),jj=norec,1,-incz)
        write(iu,*)
        do 7200 k=1,nlev
          write(iu,7002)
     *      sig(k), (nint(86400.*10.*zdvout(jj,k)/nz),jj=norec,1,-incz)
 7200   continue
 
 7000   format(/a1,' ',a/(6x,3x,40i3))
 7002   format( f6.3,3x,(40i3))
 
      endif
 
      return
      end
