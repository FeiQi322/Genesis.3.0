c
c***************************** RADCLW *********************************
c
      subroutine radclw(jslc,   ts   ,tnm     ,qnm     ,o3vmr   ,
     $                  pmid    ,pint    ,pmln    ,piln    ,plco2   ,
     $                  plh2o   ,n2o     ,ch4     ,cfc11   ,cfc12   ,
     $                  cld     ,tclrsf  ,qrl     ,flns    ,flnt    ,
     $                  flnsc   ,flntc   ,flwds   ,
     $                  dolwa)
C-----------------------------------------------------------------------
C
C Compute longwave radiation heating rates and boundary fluxes
C
C Uses broad band absorptivity/emissivity method to compute clear sky;
C assumes randomly overlapped clouds with variable cloud emissivity to
C include effects of clouds.
C
C Computes clear sky absorptivity/emissivity at lower frequency (in
C general) than the model radiation frequency; uses previously computed
C and stored values for efficiency
C
C Note: This subr contains vertical indexing which proceeds
C       from bottom to top rather than the top to bottom indexing
C       used in the rest of the model.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radclw.F,v 1.3 1995/03/03 17:44:43 bonan Exp $
c $Author: bonan $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
#include <RADSIZE.h>
C-----------------------------------------------------------------------
      integer plevp2,plevp3,plevp4
      parameter (plevp2=plev+2,plevp3=plev+3,plevp4=plev+4)
C------------------------------Commons----------------------------------
#include <RADCON.h>
#include <RADBUF.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer jslc                 ! Model latitude index
      real  ts(plon)               ! Ground (skin) temperature
C
C Input arguments which are only passed to other routines
C
      real  tnm(plon,plev),        ! Level temperature
     $     qnm(plon,plev),        ! Level moisture field
     $     o3vmr(plon,plev),      ! ozone volume mixing ratio
     $     pmid(plon,plev),       ! Level pressure
     $     pint(plon,plevp),      ! Model interface pressure
     $     pmln(plon,plev),       ! Ln(pmid)
     $     piln(plon,plevp),      ! Ln(pint)
     $     plco2(plon,plevp),     ! Path length co2
     $     plh2o(plon,plevp)      ! Path length h2o
      real  n2o(plon,plev),        ! nitrous oxide mass mixing ratio
     $     ch4(plon,plev),        ! methane mass mixing ratio
     $     cfc11(plon,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plon,plev)       ! cfc12 mass mixing ratio

      logical dolwa
C
C Input/Output arguments
C
      real  cld(plon,plevp),       ! Cloud cover
     $     tclrsf(plon,plevp)     ! Clear sky fraction
C
C Output arguments
C
      real  
     $     qrl(plon,plevp),       ! Longwave dwn-up flx at interfaces
     $     flns(plon),            ! Surface cooling flux
     $     flnt(plon),            ! Net outgoing flux
     $     flnsc(plon),           ! Clear sky surface cooing
     $     flntc(plon),           ! Net clear sky outgoing flux
     $     flwds(plon)            ! Down longwave flux at surface
C
C---------------------------Local variables-----------------------------
C
      integer     i,    ! Longitude index
     $            k,    ! Level index
     $           k1,    ! Level index
     $           k2,    ! Level index
     $           k3,    ! Level index
     $           km,    ! Level index
     $          km1,    ! Level index
     $          km2,    ! Level index
     $          km3,    ! Level index
     $          km4     ! Level index
C
      real  tmp(plon),
     $     tmp1,                ! Temporary 1
     $     absbt(plon)         ! Downward emission at model top
      real  plol(plon,plevp),   ! O3 pressure wghted path length
     $     plos(plon,plevp)    ! O3 path length
C
      real  co2em(plon,plevp),  ! Layer co2 normalized planck funct. derivative
     $     co2eml(plon,plev),  ! Interface co2 normalized planck funct. deriv.
     $     delt(plon),         ! Diff t**4 mid layer to top interface
     $     delt1(plon),        ! Diff t**4 lower intrfc to mid layer
     $     bk1(plon),          ! Absrptvty for vertical quadrature
     $     bk2(plon),          ! Absrptvty for vertical quadrature
     $     ful(plon,plevp),    ! Total upwards longwave flux
     $     fsul(plon,plevp),   ! Clear sky upwards longwave flux
     $     fdl(plon,plevp),    ! Total downwards longwave flux
     $     fsdl(plon,plevp),   ! Clear sky downwards longwv flux
     $     fclb4(plon,plev),   ! Sig t**4 for cld bottom interfc
     $     fclt4(plon,plev),   ! Sig t**4 for cloud top interfc
     $     s(plon,plevp,plevp) ! Flx integral sum
      real  tplnka(plon,plevp), ! Planck fnctn temperature
     $     s2c(plon,plevp),    ! H2o cont amount
     $     s2t(plon,plevp),    ! H2o cont temperature
     $     w(plon,plevp),      ! H2o path
     $     tplnke(plon)        ! Planck fnctn temperature
      real  h2otr(plon,plevp),  ! H2o trnmsn for o3 overlap
     $     co2t(plon,plevp),   ! Prs wghted temperature path
     $     tint(plon,plevp),   ! Interface temperature
     $     tint4(plon,plevp),  ! Interface temperature**4
     $     tlayr(plon,plevp),  ! Level temperature
     $     tlayr4(plon,plevp)  ! Level temperature**4
      real  rtclrsf(plon,plevp),! 1./tclrsf(i,k)
     $     gocp                 ! gravit/cpair
      integer klov(plon),      ! Cloud lowest level index
     $        khiv(plon),      ! Cloud highest level index
     $        khivm(plon)      ! khiv(i) - 1
      integer indx(plon),
     $        nptsc,
     $        ii,
     $        khighest
      logical done(plon),
     $        start(plon)

c
c Trace gas variables
c
      real  ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp),   ! pressure factor for ch4
     $     uptype(plon,plevp)  ! p-type continuum path length
      real  abplnk1(14,plon,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plon,plevp)  ! nearest layer factor
C
C
C------------------------------Externals--------------------------------
C
      integer intmax
      external intmax
c      external radtpl,           ! Compute path lengths
c     $         radems,           ! H2o,co2,o3 emissivity
c     $         radabs            ! H2o,co2,o3 absorptivity
C-----------------------------------------------------------------------
C
C Initialize and recompute the tclrsf array
C
      do i=1,plon
         rtclrsf(i,1) = 1.0/tclrsf(i,1)
      end do
C
      do k=1,plev
         do i=1,plon
            fclb4(i,k) = 0.
            fclt4(i,k) = 0.
            tclrsf(i,k+1) = tclrsf(i,k)*(1. - cld(i,k+1))
            rtclrsf(i,k+1) = 1./tclrsf(i,k+1)
         end do
      end do
C
C Calculate some temperatures needed to derive absorptivity and
C emissivity, as well as some h2o path lengths
C
      call radtpl(tnm     ,ts      ,qnm     ,pint    ,plh2o   ,
     $            tplnka  ,s2c     ,s2t     ,w       ,tplnke  ,
     $            tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    ,
     $            piln    )

c do emissivity and absorptivity calculations
c only if abs/ems computation
c
      if (dolwa) then
C
C Compute ozone path lengths at frequency of a/e calculation.
C
         call radoz2(o3vmr   ,pint    ,plol    ,plos    )
c
c Compute trace gas path lengths
c
         call trcpth(tnm, pint, cfc11, cfc12, n2o, ch4, qnm,
     $               ucfc11, ucfc12, un2o0,  un2o1,  uch4,
     $               uco211, uco212, uco213, uco221, uco222,
     $               uco223, bn2o0,  bn2o1,  bch4,   uptype)
C
C
C Compute total emissivity:
C
         call radems(s2c     ,s2t     ,w       ,tplnke  ,plh2o   ,
     $               pint    ,plco2   ,tint    ,tint4   ,tlayr   ,
     $               tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  ,
     $               un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $               uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $               bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  ,
     $               co2t    ,h2otr   ,abplnk1 ,abplnk2 ,jslc  )
C
C Compute total absorptivity:
C
         call radabs(pmid    ,pint    ,co2em   ,co2eml  ,tplnka  ,
     $               s2c     ,s2t     ,w       ,h2otr   ,plco2   ,
     $               plh2o   ,co2t    ,tint    ,tlayr   ,plol    ,
     $               plos    ,pmln    ,piln    ,ucfc11  ,ucfc12  ,
     $               un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $               uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $               bn2o0   ,bn2o1   ,bch4    ,abplnk1 ,abplnk2 ,
     $               jslc  )

      endif
C
C Find the lowest and highest level cloud for each grid point
C Note: Vertical indexing here proceeds from bottom to top
C
      do i=1,plon
         klov(i) = 0
         done(i) = .false.
      end do
      do k=1,plev
         do i=1,plon
            if (.not.done(i) .and. cld(i,plevp2-k).gt.0.0) then
               done(i) = .true.
               klov(i) = k
            end if
         end do
      end do
      call whenne (plon,klov,1,0,indx,nptsc)
      do i=1,plon
         khiv(i) = klov(i)
         done(i) = .false.
      end do
      do k=plev,1,-1
!DIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (.not.done(i) .and. cld(i,plevp2-k).gt.0.0) then
               done(i) = .true.
               khiv(i) = k
            end if
         end do
      end do
      do i=1,plon
         khivm(i) = khiv(i) - 1
      end do
C
C Note: Vertical indexing here proceeds from bottom to top
C
      do ii=1,nptsc
         i=indx(ii)
         do k=klov(i),khiv(i)
            fclt4(i,plevp-k) = stebol*tint4(i,plevp2-k)
            fclb4(i,plevp-k) = stebol*tint4(i,plevp3-k)
         end do
      end do
C
C Compute sums used in integrals (all longitude points)
C
C Definition of bk1 & bk2 depends on finite differencing.  for
C trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
C layers only.
C
C delt=t**4 in layer above current sigma level km.
C delt1=t**4 in layer below current sigma level km.
C
      do i=1,plon
         delt(i) = tint4(i,plev) - tlayr4(i,plevp)
         delt1(i) = tlayr4(i,plevp) - tint4(i,plevp)
         s(i,plevp,plevp) = stebol*(delt1(i)*absnxt(i,plev,1,jslc) +
     $                              delt (i)*absnxt(i,plev,4,jslc))
         s(i,plev,plevp)  = stebol*(delt (i)*absnxt(i,plev,2,jslc) +
     $                              delt1(i)*absnxt(i,plev,3,jslc))
      end do
      do k=1,plev-1
         do i=1,plon
            bk2(i) = (abstot(i,k,plev,jslc) +
     &                abstot(i,k,plevp,jslc))*0.5
            bk1(i) = bk2(i)
            s(i,k,plevp) = stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
         end do
      end do
C
C All k, km>1
C
      do km=plev,2,-1
         do i=1,plon
            delt(i)  = tint4(i,km-1) - tlayr4(i,km)
            delt1(i) = tlayr4(i,km) - tint4(i,km)
         end do
         do k=plevp,1,-1
            if (k.eq.km) then
               do i=1,plon
                  bk2(i) = absnxt(i,km-1,4,jslc)
                  bk1(i) = absnxt(i,km-1,1,jslc)
               end do
            else if(k.eq.km-1) then
               do i=1,plon
                  bk2(i) = absnxt(i,km-1,2,jslc)
                  bk1(i) = absnxt(i,km-1,3,jslc)
               end do
            else
               do i=1,plon
                  bk2(i) = (abstot(i,k,km-1,jslc) +
     &                      abstot(i,k,km,jslc))*0.5
                  bk1(i) = bk2(i)
               end do
            end if
            do i=1,plon
               s(i,k,km) = s(i,k,km+1) + stebol*
     $                    (bk2(i)*delt(i) + bk1(i)*delt1(i))
            end do
         end do
      end do
C
C Computation of clear sky fluxes always set first level of fsul
C
      do i=1,plon
         fsul(i,plevp) = stebol*(ts(i)**4)
      end do
C
C Downward clear sky fluxes store intermediate quantities in down flux
C Initialize fluxes to clear sky values.
C
      do i=1,plon
         tmp(i) = fsul(i,plevp) - stebol*tint4(i,plevp)
         fsul(i,1) = fsul(i,plevp)-abstot(i,1,plevp,jslc)*tmp(i)
     &             + s(i,1,2)
         fsdl(i,1) = stebol*(tplnke(i)**4)*emstot(i,1,jslc)
         ful(i,1) = fsul(i,1)
         fdl(i,1) = fsdl(i,1)
      end do
C
C fsdl(i,plevp) assumes isothermal layer
C
      do k=2,plev
         do i=1,plon
            fsul(i,k) = fsul(i,plevp)-abstot(i,k,plevp,jslc)*tmp(i)
     $                + s(i,k,k+1)
            ful(i,k) = fsul(i,k)
            fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k,jslc) -
     $                  (s(i,k,2) - s(i,k,k+1))
            fdl(i,k) = fsdl(i,k)
         end do
      end do
C
C Store the downward emission from level 1 = total gas emission * sigma
C t**4.  fsdl does not yet include all terms
C
      do i=1,plon
         ful(i,plevp) = fsul(i,plevp)
         absbt(i) = stebol*(tplnke(i)**4)*emstot(i,plevp,jslc)
         fsdl(i,plevp) = absbt(i) - s(i,plevp,2)
         fdl(i,plevp) = fsdl(i,plevp)
      end do
C
C Modifications for clouds
C
C Further qualify longitude subset for computations.  Select only those
C locations where there are clouds (total cloud fraction <= 1.e-3 treated
C as clear)
C
      call whenflt(plon,tclrsf(1,plevp),1,0.999,indx,nptsc)
C
C Compute downflux at level 1 for cloudy sky
C
      do ii=1,nptsc
         i=indx(ii)
C
C First clear sky flux plus flux from cloud at level 1
C
         fdl(i,plevp) = fsdl(i,plevp)*tclrsf(i,plev)*
     $         rtclrsf(i,plevp-khiv(i)) + fclb4(i,plev-1)*cld(i,plev)
      end do
C
C Flux emitted by other layers
C Note: Vertical indexing here proceeds from bottom to top
C
      khighest = khiv(intmax(plon,khiv,1))
      do km=3,khighest
         km1 = plevp - km
         km2 = plevp2 - km
         km4 = plevp4 - km
!DIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (km.le.khiv(i)) then
               tmp1 = cld(i,km2)*tclrsf(i,plev)*rtclrsf(i,km2)
               fdl(i,plevp) = fdl(i,plevp) +
     $                       (fclb4(i,km1) - s(i,plevp,km4))*tmp1
            end if
         end do
      end do
C
C Note: Vertical indexing here proceeds from bottom to top
C
      do k=1,khighest-1
         k1 = plevp - k
         k2 = plevp2 - k
         k3 = plevp3 - k
!DIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.ge.klov(i) .and. k.le.khivm(i)) then
               ful(i,k2) = fsul(i,k2)*(tclrsf(i,plevp)*rtclrsf(i,k1))
            end if
         end do
         do km=1,k
            km1 = plevp - km
            km2 = plevp2 - km
            km3 = plevp3 - km
!DIR$ IVDEP
            do ii=1,nptsc
               i=indx(ii)
               if (k.le.khivm(i) .and. km.ge.klov(i) .and.
     $             km.le.khivm(i)) then
C
                  ful(i,k2) = ful(i,k2) +
     $                 (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))*
     $                 cld(i,km2)*(tclrsf(i,km1)*rtclrsf(i,k1))
               end if
            end do
         end do            ! km=1,k
      end do               ! k=1,khighest-1
C
      do k=1,plevp
         k2 = plevp2 - k
         k3 = plevp3 - k
         do i=1,plon
            start(i) = .false.
         end do
!DIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.ge.khiv(i)) then
               start(i) = .true.
               ful(i,k2) = fsul(i,k2)*tclrsf(i,plevp)*
     $                     rtclrsf(i,plevp-khiv(i))
            end if
         end do
         do km=1,khighest
            km1 = plevp - km
            km2 = plevp2 - km
            km3 = plevp3 - km
!DIR$ IVDEP
            do ii=1,nptsc
               i=indx(ii)
               if(start(i) .and. km.ge.klov(i) .and. km.le.khiv(i)) then
                  ful(i,k2) = ful(i,k2)  +
     $              (cld(i,km2)*tclrsf(i,km1)*rtclrsf(i,plevp-khiv(i)))*
     $              (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))
               end if
            end do
         end do         ! km=1,khighest
      end do            ! k=1,plevp
C
C Computation of the downward fluxes
C
      do k=2,khighest-1
         k1 = plevp - k
         k2 = plevp2 - k
         k3 = plevp3 - k
!DIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.le.khivm(i)) fdl(i,k2) = 0.
         end do
         do km=k+1,khighest
            km1 = plevp - km
            km2 = plevp2 - km
            km4 = plevp4 - km
!DIR$ IVDEP
            do ii=1,nptsc
               i=indx(ii)
               if (k.le.khiv(i) .and. km.ge.max0(k+1,klov(i)) .and.
     $             km.le.khiv(i)) then
C
                  fdl(i,k2) = fdl(i,k2) +
     $                 (cld(i,km2)*tclrsf(i,k1)*rtclrsf(i,km2))*
     $                 (fclb4(i,km1) - s(i,k2,km4) + s(i,k2,k3))
               end if
            end do
         end do            ! km=k+1,khighest
!DIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.le.khivm(i)) then
               fdl(i,k2) = fdl(i,k2) + fsdl(i,k2)*
     $              (tclrsf(i,k1)*rtclrsf(i,plevp-khiv(i)))
            end if
         end do
      end do               ! k=1,khighest-1
C
C End cloud modification loops
C
C All longitudes: store history tape quantities
C
      do i=1,plon
C
C Downward longwave flux
C
         flwds(i) = fdl(i,plevp)
C
C Net flux
C
         flns(i) = ful(i,plevp) - fdl(i,plevp)
C
C Clear sky flux at top of atmosphere
C
         flntc(i) = fsul(i,1)
         flnsc(i) = fsul(i,plevp) - fsdl(i,plevp)
C
C Outgoing ir
C
         flnt(i) = ful(i,1) - fdl(i,1)
      end do
C
C Computation of longwave heating (flux down-up at interfaces)
C
c     gocp = gravit/cpair
c     do k=1,plev
      do k=1,plevp
         do i=1,plon
c           qrl(i,k) = (ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))*
c    $                 gocp/((pint(i,k) - pint(i,k+1)))
            qrl(i,k) = fdl(i,k) - ful(i,k)
         end do
      end do
C
      return
      end
c
c***************************** BLKDAT *********************************
c
c
      block data blkdat
C-----------------------------------------------------------------------
C
C Contains all data statements for common variables
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      L. Bath, June 1992
C
C-----------------------------------------------------------------------
c
c $Id: blkdat.F,v 1.1.1.1 1995/02/09 23:26:36 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
#include <RADSIZE.h>
C------------------------------Commons----------------------------------
#include <RADCAE.h>
C-----------------------------------------------------------------------
C
C /CRDCAE/
C H2O EMISSIVITY AND ABSORTIVITY COEFFICIENTS
C
      data coefa/1.01400e+00, 6.41695e-03, 2.85787e-05,
     $           1.01320e+00, 6.86400e-03, 2.96961e-05,
     $           1.02920e+00, 1.01680e-02, 5.30226e-05,
     $           1.02743e+00, 9.85113e-03, 5.00233e-05/
C
      data coefb/8.85675e+00,-3.51620e-02, 2.38653e-04,-1.71439e-06,
     $           5.73841e+00,-1.91919e-02, 1.65993e-04,-1.54665e-06,
     $           6.64034e+00, 1.56651e-02,-9.73357e-05, 0.0,
     $           7.09281e+00, 1.40056e-02,-1.15774e-04, 0.0/
C
      data coefc/9.90127e-01, 1.22475e-03, 4.90135e-06,
     $           9.89753e-01, 1.97081e-03, 3.42046e-06,
     $           9.75230e-01, 1.03341e-03, 0.0,
     $           9.77366e-01, 8.60014e-04, 0.0/
C
      data coefd/7.03047e-01,-2.63501e-03,-1.57023e-06, 0.0,
     $           5.29269e-01,-3.14754e-03, 4.39595e-06, 0.0,
     $           7.88193e-02, 1.31290e-03, 4.25827e-06,-1.23982e-08,
     $           1.62744e-01, 2.22847e-03, 2.60102e-06,-4.30133e-08/
C
      data coefe/3.93137e-02,-4.34341e-05, 3.74545e-07,
     $           3.67785e-02,-3.10794e-05, 2.94436e-07,
     $           7.42500e-02, 3.97397e-05, 0.0,
     $           7.52859e-02, 4.18073e-05, 0.0/
C
      data coeff/2.2037 e-01,1.39719e-03,-7.32011e-06,
     $          -1.40262e-08,2.13638e-10,-2.35955e-13,
     $           3.07431e-01,8.27225e-04,-1.30067e-05,
     $           3.49847e-08,2.07835e-10,-1.98937e-12/
C
      data coefg/9.04489e+00,-9.56499e-03,
     $           1.80898e+01,-1.91300e-02,
     $           8.72239e+00,-9.53359e-03,
     $           1.74448e+01,-1.90672e-02/
C
      data coefh/5.46557e+01,-7.30387e-02,
     $           1.09311e+02,-1.46077e-01,
     $           5.11479e+01,-6.82615e-02,
     $           1.02296e+02,-1.36523e-01/
C
      data coefi/3.31654e-01,-2.86103e-04,-7.87860e-06,
     $           5.88187e-08,-1.25340e-10,-1.37731e-12,
     $           3.14365e-01,-1.33872e-03,-2.15585e-06,
     $           6.07798e-08,-3.45612e-10,-9.34139e-15/
C
      data coefj/2.82096e-02,2.47836e-04,1.16904e-06,
     $           9.27379e-02,8.04454e-04,6.88844e-06/
C
      data coefk/2.48852e-01,2.09667e-03,2.60377e-06,
     $           1.03594e+00,6.58620e-03,4.04456e-06/
C
C Narrow band data for H2O
C 200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
C
      data realk/  0.18967069430426e-04, 0.70172244841851e-04   /
      data   st /  0.31930234492350e-03, 0.97907319939060e-03   /
      data   a1 /  0.28775403075736e-01, 0.23236701470511e-01   /
      data   a2 / -0.57966222388131e-04,-0.95105504388411e-04   /
      data   b1 /  0.29927771523756e-01, 0.21737073577293e-01   /
      data   b2 / -0.86322071248593e-04,-0.78543550629536e-04   /
      end
c
c***************************** CLDEMS *********************************
c
      subroutine cldems(clwp, fice, rei, emis)
C-----------------------------------------------------------------------
C
C Compute cloud emissivity using cloud liquid water path (g/m**2)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, J. Kiehl, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: cldems.F,v 1.2 1995/02/17 21:28:09 jhack Exp $
c $Author: jhack $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
      real  kabs                   ! longwave absorption coeff (m**2/g)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real  clwp(plon,plev),      ! cloud liquid water path (g/m**2)
     $     rei(plon,plev),       ! ice effective drop size (microns)
     $     fice(plon,plev)       ! fractional ice content within cloud
C
C Output arguments
C
      real  emis(plon,plev)       ! cloud emissivity (fraction)
C
C---------------------------Local workspace-----------------------------
C
      integer i,k                 ! longitude, level indices
      real  kabsl,                 ! longwave absorption coeff (m**2/g)
     $     kabsi                  ! ice absorption coefficient
      parameter (kabsl = 0.090361)
C
C-----------------------------------------------------------------------
C
      do k=1,plev
         do i=1,plon
            kabsi = 0.005 + 1./rei(i,k)
            kabs = kabsl*(1.-fice(i,k))+kabsi*fice(i,k)
            emis(i,k) = 1. - exp(-1.66*kabs*clwp(i,k))
         end do
      end do
C
      return
      end
c
c***************************** RADABS *********************************
c
      subroutine radabs(pbr    ,pnm     ,co2em    ,co2eml  ,tplnka  ,
     $                  s2c    ,s2t     ,w        ,h2otr   ,plco2   ,
     $                  plh2o  ,co2t    ,tint     ,tlayr   ,plol    ,
     $                  plos   ,pmln    ,piln     ,ucfc11  ,ucfc12  ,
     $                  un2o0  ,un2o1   ,uch4     ,uco211  ,uco212  ,
     $                  uco213 ,uco221  ,uco222   ,uco223  ,uptype  ,
     $                  bn2o0  ,bn2o1   ,bch4    ,abplnk1  ,abplnk2 ,
     $                  jslc  )
C-----------------------------------------------------------------------
C
C Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
C
C h2o  ....  Uses nonisothermal emissivity for water vapor from
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
C
C co2  ....  Uses absorptance parameterization of the 15 micro-meter
C            (500 - 800 cm-1) band system of Carbon Dioxide, from
C            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
C            of the Absorptance Due to the 15 micro-meter Band System
C            of Carbon Dioxide Jouranl of Geophysical Research,
C            vol. 96., D5, pp 9013-9019.
C            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
C            are also included.
C
C o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
C            band system of ozone, from Ramanathan, V. and R.Dickinson,
C            1979: The Role of stratospheric ozone in the zonal and
C            seasonal radiative energy balance of the earth-troposphere
C            system. Journal of the Atmospheric Sciences, Vol. 36,
C            pp 1084-1104
C
C ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
C
C n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
C            bands of nitrous oxide
C
C cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
C            micron bands of CFC11
C
C cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
C            micron bands of CFC12
C
C
C Computes individual absorptivities for non-adjacent layers, accounting
C for band overlap, and sums to obtain the total; then, computes the
C nearest layer contribution.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C
C-----------------------------------------------------------------------
c
c $Id: radabs.F,v 1.1.16.1 1998/08/17 21:41:08 zender Exp $
c
      implicit none
      integer jslc
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCAE.h>
#include <RADCON.h>
C Output arguments in radbuf.cb 
#include <RADBUF.h>
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real  pbr(plon,plev),           ! Prssr at mid-levels (dynes/cm2)
     $     pnm(plon,plevp),          ! Prssr at interfaces (dynes/cm2)
     $     co2em(plon,plevp),        ! Co2 emissivity function
     $     co2eml(plon,plev),        ! Co2 emissivity function
     $     tplnka(plon,plevp),       ! Planck fnctn level temperature
     $     s2c(plon,plevp),          ! H2o continuum path length
     $     s2t(plon,plevp),          ! H2o tmp and prs wghted path
     $     w(plon,plevp),            ! H2o prs wghted path
     $     h2otr(plon,plevp),        ! H2o trnsmssn fnct for o3 overlap
     $     plco2(plon,plevp),        ! Co2 prs wghted path length
     $     plh2o(plon,plevp),        ! H2o prs wfhted path length
     $     co2t(plon,plevp),         ! Tmp and prs wghted path length
     $     tint(plon,plevp),         ! Interface temperatures
     $     tlayr(plon,plevp),        ! K-1 level temperatures
     $     plol(plon,plevp),         ! Ozone prs wghted path length
     $     plos(plon,plevp)          ! Ozone path length
      real  pmln(plon,plev),          ! Ln(pmidm1)
     $     piln(plon,plevp)          ! Ln(pintm1)
c
c   Trace gas variables
c
      real  ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     uptype(plon,plevp), ! continuum path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp)    ! pressure factor for ch4
      real  abplnk1(14,plon,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plon,plevp)  ! nearest layer factor
      real  abstrc(plon)            ! total trace gas absorptivity
      real  bplnk(14,plon,4) ! Planck functions for sub-divided layers
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        k1,           ! Level index
     $        k2,           ! Level index
     $        kn,           ! Nearest level index
     $        iband         ! Band  index
C
      real  pnew(plon),     ! Effective pressure for H2O vapor linewidth
     $     trline(plon,2), ! Transmission due to H2O lines in window
     $     u(plon),        ! Pressure weighted H2O path length
     $     tbar(plon,4),   ! Mean layer temperature
     $     emm(plon,4),    ! Mean co2 emissivity
     $     o3emm(plon,4),  ! Mean o3 emissivity
     $     o3bndi,          ! Ozone band parameter
     $     temh2o(plon,4), ! Mean layer temperature equivalent to tbar
     $     k21,             ! Exponential coefficient used to calculate
C                           !  rotation band transmissvty in the 650-800
C                           !  cm-1 region (tr1)
     $     k22,             ! Exponential coefficient used to calculate
C                           !  rotation band transmissvty in the 500-650
C                           !  cm-1 region (tr2)
     $     uc1(plon)       ! H2o continuum pathlength in 500-800 cm-1
      real  to3h2o(plon),   ! H2o trnsmsn for overlap with o3
     $     pi,              ! For co2 absorptivity computation
     $     sqti(plon),     ! Used to store sqrt of mean temperature
     $     et,              ! Co2 hot band factor
     $     et2,             ! Co2 hot band factor squared
     $     et4,             ! Co2 hot band factor to fourth power
     $     omet,            ! Co2 stimulated emission term
     $     f1co2,           ! Co2 central band factor
     $     f2co2(plon),    ! Co2 weak band factor
     $     f3co2(plon),    ! Co2 weak band factor
     $     t1co2(plon),    ! Overlap factr weak bands on strong band
     $     sqwp,            ! Sqrt of co2 pathlength
     $     f1sqwp(plon)    ! Main co2 band factor
      real  oneme,           ! Co2 stimulated emission term
     $     alphat,          ! Part of the co2 stimulated emission term
     $     wco2,            ! Constants used to define co2 pathlength
     $     posqt,           ! Effective pressure for co2 line width
     $     u7(plon),       ! Co2 hot band path length
     $     u8,              ! Co2 hot band path length
     $     u9,              ! Co2 hot band path length
     $     u13,             ! Co2 hot band path length
     $     rbeta7(plon),   ! Inverse of co2 hot band line width par
     $     rbeta8,          ! Inverse of co2 hot band line width par
     $     rbeta9,          ! Inverse of co2 hot band line width par
     $     rbeta13          ! Inverse of co2 hot band line width par
      real  tpatha(plon),     ! For absorptivity computation
     $     a,                 ! Eq(2) in table A3a of R&D
     $     abso(plon,6),     ! Absorptivity for various gases/bands
     $     dtp(plon),        ! Path temp minus 300 K used in h2o
C                             !  rotation band absorptivity
     $     dtx(plon),        ! Planck temperature minus 250 K
     $     dty(plon),        ! Path temperature minus 250 K
     $     dtz(plon),        ! Planck temperature minus 300 K
     $     term1(plon,4),    ! Equation(5) in table A3a of R&D(1986)
     $     term2(plon,4)     ! Delta a(Te) in table A3a of R&D(1986)
      real  term3(plon,4),    ! DB/dT function for rotation and
C                             !  vibration-rotation band absorptivity
     $     term4(plon,4),    ! Equation(6) in table A3a of R&D(1986)
     $     term5(plon,4),    ! Delta a(Tp) in table A3a of R&D(1986)
     $     term6(plon,plevp),! DB/dT function for window region
     $     term7(plon,2),    ! Kl_inf(i) in eq(8) of table A3a of R&D
     $     term8(plon,2),    ! Delta kl_inf(i) in eq(8)
     $     term9(plon,plevp),! DB/dT function for 500-800 cm-1 region
     $     tr1,               ! Eqn(6) in table A2 of R&D for 650-800
     $     tr10(plon),       ! Eqn (6) times eq(4) in table A2
C                             !  of R&D for 500-650 cm-1 region
     $     tr2                ! Eqn(6) in table A2 of R&D for 500-650
      real  tr5,               ! Eqn(4) in table A2 of R&D for 650-800
     $     tr6,               ! Eqn(4) in table A2 of R&D for 500-650
     $     tr9(plon),        ! Equation (6) times eq(4) in table A2
C                             !  of R&D for 650-800 cm-1 region
     $     uc(plon)          ! Y + 0.002U in eq(8) of table A2 of R&D
      real  sqrtu(plon),      ! Sqrt of pressure weighted h20 pathlength
     $     fwk(plon),        ! Equation(33) in R&D far wing correction
     $     fwku(plon),       ! GU term in eqs(1) and (6) in table A2
     $     r2st(2),           ! 1/(2*beta) in eq(10) in table A2
     $     dtyp15(plon),     ! DeltaTp in eqs(11) & (12) in table A3a
     $     dtyp15sq(plon),   ! (DeltaTp)^2 in eqs(11) & (12) table A3a
     $     to3co2(plon),     ! P weighted temp in ozone band model
     $     dpnm(plon),       ! Pressure difference between two levels
     $     pnmsq(plon,plevp),! Pressure squared
     $     dw(plon),         ! Amount of h2o between two levels
     $     uinpl(plon,4),    ! Nearest layer subdivision factor
     $     winpl(plon,4),    ! Nearest layer subdivision factor
     $     zinpl(plon,4),    ! Nearest layer subdivision factor
     $     pinpl(plon,4),    ! Nearest layer subdivision factor
     $     dplh2o(plon)      ! Difference in press weighted h2o amount
      real  r80257,            ! Conversion factor for h2o pathlength
     $     r293,              ! 1/293
     $     r250,              ! 1/250
     $     r3205,             ! Line width factor for o3 (see R&Di)
     $     r300,              ! 1/300
     $     rsslp,             ! Reciprocal of sea level pressure
     $     r2sslp             ! 1/2 of rsslp
      real   ds2c,     ! Y in eq(7) in table A2 of R&D
     $      a11,      ! A1 in table A3b for rotation band absorptivity
     $      a31,      ! A3 in table A3b for rotation band absorptivity
     $      a21,      ! First part in numerator of A2 in table A3b
     $      a22,      ! Second part in numerator of A2 in table A3b
     $      a23,      ! Denominator of A2 in table A3b (rotation band)
     $      t1t4,     ! Eq(3) in table A3a of R&D
     $      t2t5,     ! Eq(4) in table A3a of R&D
     $      rsum,     ! Eq(1) in table A2 of R&D
     $      a41,      ! Numerator in A2 in Vib-rot abstivity(table A3b)
     $      a51,      ! Denominator in A2 in Vib-rot (table A3b)
     $      a61       ! A3 factor for Vib-rot band in table A3b
      real   phi,      ! Eq(11) in table A3a of R&D
     $      psi,      ! Eq(12) in table A3a of R&D
     $      cf812,    ! Eq(11) in table A2 of R&D
     $      ubar,     ! H2o scaled path see comment for eq(10) table A2
     $      pbar,     ! H2o scaled pres see comment for eq(10) table A2
     $      g4        ! Arguement in exp() in eq(10) table A2
C
      real   dplos,    ! Ozone pathlength eq(A2) in R&Di
     $      dplol,    ! Presure weighted ozone pathlength
     $      tlocal,   ! Local interface temperature
     $      beta,     ! Ozone mean line parameter eq(A3) in R&Di
C                       (includes Voigt line correction factor)
     $      rphat,    ! Effective pressure for ozone beta
     $      tcrfac,   ! Ozone temperature factor table 1 R&Di
     $      tmp1,     ! Ozone band factor see eq(A1) in R&Di
     $      u1,       ! Effective ozone pathlength eq(A2) in R&Di
     $      realnu,   ! 1/beta factor in ozone band model eq(A1)
     $      tmp2,     ! Ozone band factor see eq(A1) in R&Di
     $      u2,       ! Effective ozone pathlength eq(A2) in R&Di
     $      rsqti     ! Reciprocal of sqrt of path temperature
C
      real   tpath,    ! Path temperature used in co2 band model
     $      tmp3,     ! Weak band factor see K&B
     $      rdpnmsq,  ! Reciprocal of difference in press^2
     $      rdpnm,    ! Reciprocal of difference in press
     $      p1,       ! Mean pressure factor
     $      p2,       ! Mean pressure factor
     $      dtym10,   ! T - 260 used in eq(9) and (10) table A3a
     $      dplco2,   ! Co2 pathlength
     $      corfac,   ! Correction factors in table A3b
     $      g2,       ! Part of arguement in eq(10) in table A2
     $      te,       ! A_0 T factor in ozone model table 1 of R&Di
     $      denom     ! Denominator in eq(8) of table A3a of R&D
c
c
      real  th2o(plon),
     $     tco2(plon),
     $     to3(plon)
      integer wvl
C
C
C Transmission terms for various spectral intervals:
C
C     real  trab1(plon)   ! H2o     0 -  800 cm-1
      real  trab2(plon)   ! H2o   500 -  800 cm-1
C     real  trab3(plon)   ! Co2   band system
      real  trab4(plon)   ! H2o   800 - 1000 cm-1
C     real  trab5(plon)   ! 9.6 micrometer band
      real  trab6(plon)   ! H2o  1000 - 1200 cm-1
C     real  trab7(plon)   ! H2o  1200 - 2200 cm-1
C
      real  bndfct, ! Band absorptance parameter for co2
     $     absbnd  ! Proportional to co2 band absorptance
C
      real  dbvtit(plon,plevp),      ! Intrfc drvtv plnck fnctn for o3
     $     dbvtly(plon,plev)        ! Level drvtv plnck fnctn for o3
C
C--------------------------Statement function---------------------------
C
      real  dbvt,t     ! Planck fnctn tmp derivative for o3
C
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
C
C-----------------------------------------------------------------------
C
C Initialize
C
      do k=1,plev
         do i=1,plon
            dbvtly(i,k) = dbvt(tlayr(i,k+1))
            dbvtit(i,k) = dbvt(tint(i,k))
         end do
      end do
      do i=1,plon
         dbvtit(i,plevp) = dbvt(tint(i,plevp))
      end do
C
      r80257  = 1./8.0257e-04
      r293    = 1./293.
      r250    = 1./250.
      r3205   = 1./.3205
      r300    = 1./300.
      rsslp   = 1./sslp
      r2sslp  = 1./(2.*sslp)
      r2st(1) = 1./(2.*st(1))
      r2st(2) = 1./(2.*st(2))
      bndfct  = 2.0*22.18/(sqrt(196.)*300.)
C
C Non-adjacent layer absorptivity:
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C abso(i,6)   co2 15  micrometer band system
C
      do k=1,plevp
         do i=1,plon
            pnmsq(i,k) = pnm(i,k)**2
            dtx(i) = tplnka(i,k) - 250.
            term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                   (1. +  c9*dtx(i)*(1. + c11*dtx(i)*
     $                   (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
            term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                   (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                   (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
         end do
      end do
C
C Non-nearest layer level loops
C
      do 200 k1=plevp,1,-1
         do 100 k2=plevp,1,-1
            if(k1.eq.k2) go to 100
            do i=1,plon
               dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
               u(i)      = abs(dplh2o(i))
               sqrtu(i)  = sqrt(u(i))
               ds2c      = abs(s2c(i,k1) - s2c(i,k2))
               dw(i)     = abs(w(i,k1) - w(i,k2))
               uc1(i)    = (ds2c + 1.7e-3*u(i))*(1. +  2.*ds2c)/
     $                                          (1. + 15.*ds2c)
               uc(i)     = ds2c + 2.e-3*u(i)
            end do
            do i=1,plon
               pnew(i)   = u(i)/dw(i)
               tpatha(i) = (s2t(i,k1) - s2t(i,k2))/dplh2o(i)
               dtx(i)      = tplnka(i,k2) - 250.
               dty(i)      = tpatha(i)    - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
            end do
            do iband=2,4,2
               do i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
               end do
            end do
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C
            do i=1,plon
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac  = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4    = term1(i,2)*term4(i,2)
               t2t5    = term2(i,2)*term5(i,2)
               a       = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)  = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i) = fwk(i)*u(i)
               rsum    = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
C               trab1(i)  = rsum
            end do
C
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C
            do i=1,plon
               a41   = 1.75 - 3.960e-03*dtz(i)
               a51   = 1.00 + 1.3*sqrtu(i)
               a61   = 1.00 + 1.250e-03*dtp(i) + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
C               trab7(i)  = rsum
            end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
            do k=1,2
               do i=1,plon
                  phi   = exp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = exp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(sqrt(g2) - 1.)
                  trline(i,k) = exp(-g4)
               end do
            end do
            do i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
            end do
C
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C
            do i=1,plon
               k21    = term7(i,1) + term8(i,1)/
     $             (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
               k22    = term7(i,2) + term8(i,2)/
     $             (1. + (c28 + c29*(dty(i)-10.))*sqrtu(i))
               tr1    = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2    = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5    = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
               tr6    = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)   = tr1*tr5
               tr10(i)  = tr2*tr6
               th2o(i) = tr10(i)
               trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
               trab4(i) = exp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i) = exp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               abso(i,3) = term6(i,k2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4) = term9(i,k2)*.5*(tr1 - tr9(i) + tr2 - tr10(i))
            end do
            if(k2.lt.k1) then
               do i=1,plon
                  to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
               end do
            else
               do i=1,plon
                  to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
               end do
            end if
C
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C
            do i=1,plon
               dpnm(i)  = pnm(i,k1) - pnm(i,k2)
               to3co2(i)=(pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/
     $                   dpnm(i)
               te       = (to3co2(i)*r293)**.7
               dplos    = plos(i,k1) - plos(i,k2)
               dplol    = plol(i,k1) - plol(i,k2)
               u1       = 18.29*abs(dplos)/te
               u2       = .5649*abs(dplos)/te
               rphat    = dplol/dplos
               tlocal   = tint(i,k2)
               tcrfac   = sqrt(tlocal*r250)*te
               beta     = r3205*(rphat + dpfo3*tcrfac)
               realnu   = te/beta
               tmp1     = u1/sqrt(4. + u1*(1. + realnu))
               tmp2     = u2/sqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*log(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
               to3(i)   = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
C               trab5(i)  = 1.-(o3bndi/(1060-980.))
            end do
C
C abso(i,6)      co2 15  micrometer band system
C
            do i=1,plon
               sqwp      = sqrt(abs(plco2(i,k1) - plco2(i,k2)))
               et        = exp(-480./to3co2(i))
               sqti(i)   = sqrt(to3co2(i))
               rsqti     = 1./sqti(i)
               et2       = et*et
               et4       = et2*et2
               omet      = 1. - 1.5*et2
               f1co2     = 899.70*omet*
     $                   (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i) = f1co2*sqwp
               t1co2(i)  = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme     = 1. - et2
               alphat    = oneme**3*rsqti
               pi        = abs(dpnm(i))
               wco2      =  2.5221*co2vmr*pi*rga
               u7(i)     =  4.9411e4*alphat*et2*wco2
               u8        =  3.9744e4*alphat*et4*wco2
               u9        =  1.0447e5*alphat*et4*et2*wco2
               u13       = 2.8388e3*alphat*et4*wco2
               tpath     = to3co2(i)
               tlocal    = tint(i,k2)
               tcrfac    = sqrt(tlocal*r250*tpath*r300)
               posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp +
     $                     dpfco2*tcrfac)*rsqti
               rbeta7(i) = 1./(5.3228*posqt)
               rbeta8    = 1./(10.6576*posqt)
               rbeta9    = rbeta7(i)
               rbeta13   = rbeta9
               f2co2(i)  = (u7(i)/sqrt(4. + u7(i)*(1. + rbeta7(i)))) +
     $                     (u8/sqrt(4. + u8*(1. + rbeta8))) +
     $                     (u9/sqrt(4. + u9*(1. + rbeta9)))
               f3co2(i)  = u13/sqrt(4. + u13*(1. + rbeta13))
            end do
            if (k2.ge.k1) then
               do i=1,plon
                  sqti(i) = sqrt(tlayr(i,k2))
               end do
            end if
C
            do i=1,plon
               tmp1      = log(1. + f1sqwp(i))
               tmp2      = log(1. + f2co2(i))
               tmp3      = log(1. + f3co2(i))
               absbnd    = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
               abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
               tco2(i)=1./(1.0+10.0*(u7(i)
     &                         /sqrt(4.+u7(i)*(1.+rbeta7(i)))))
C               trab3(i)  = 1. - bndfct*absbnd
            end do
c
c     Calculate absorptivity due to trace gases
c
            call trcab(k1, k2, ucfc11, ucfc12, un2o0,  un2o1,
     $                         uch4,   uco211, uco212, uco213,
     $                         uco221, uco222, uco223, bn2o0,
     $                         bn2o1,  bch4,   to3co2, pnm,
     $                         dw,     pnew,   s2c,    uptype,
     $                         u,      abplnk1,tco2,   th2o,
     $                         to3,    abstrc)
C
C Sum total absorptivity
C
            do i=1,plon
               abstot(i,k1,k2,jslc) =
     &                           abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
     $                           + abstrc(i)
            end do
  100    continue
  200 continue        ! End of non-nearest layer level loops
C
C Non-adjacent layer absorptivity:
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C abso(i,6)   co2 15  micrometer band system
C
C Nearest layer level loop
C
      do 500 k2=plev,1,-1
         do i=1,plon
            tbar(i,1)   = 0.5*(tint(i,k2+1) + tlayr(i,k2+1))
            emm(i,1)    = 0.5*(co2em(i,k2+1) + co2eml(i,k2))
            tbar(i,2)   = 0.5*(tlayr(i,k2+1) + tint(i,k2))
            emm(i,2)    = 0.5*(co2em(i,k2) + co2eml(i,k2))
            tbar(i,3)   = 0.5*(tbar(i,2) + tbar(i,1))
            emm(i,3)    = emm(i,1)
            tbar(i,4)   = tbar(i,3)
            emm(i,4)    = emm(i,2)
            o3emm(i,1)  = 0.5*(dbvtit(i,k2+1) + dbvtly(i,k2))
            o3emm(i,2)  = 0.5*(dbvtit(i,k2) + dbvtly(i,k2))
            o3emm(i,3)  = o3emm(i,1)
            o3emm(i,4)  = o3emm(i,2)
            temh2o(i,1) = tbar(i,1)
            temh2o(i,2) = tbar(i,2)
            temh2o(i,3) = tbar(i,1)
            temh2o(i,4) = tbar(i,2)
            dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
         end do
c----------------------------------------------------------
c  Weighted Planck functions for trace gases
c
      do wvl = 1,14
         do i = 1,plon
            bplnk(wvl,i,1) = 0.5*(abplnk1(wvl,i,k2+1) +
     $                          abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5*(abplnk1(wvl,i,k2) +
     $                          abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
         end do
      end do
c---------------------------------------------------------
         do i=1,plon
            rdpnmsq    = 1./(pnmsq(i,k2+1) - pnmsq(i,k2))
            rdpnm      = 1./dpnm(i)
            p1         = .5*(pbr(i,k2) + pnm(i,k2+1))
            p2         = .5*(pbr(i,k2) + pnm(i,k2  ))
            uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
            uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
            uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
            uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
            winpl(i,1) = (.5*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
            winpl(i,2) = (.5*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
            winpl(i,3) = (.5*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*
     $                   rdpnm
            winpl(i,4) = (.5*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*
     $                   rdpnm
            tmp1       = 1./(piln(i,k2+1) - piln(i,k2))
            tmp2       = piln(i,k2+1) - pmln(i,k2)
            tmp3       = piln(i,k2  ) - pmln(i,k2)
            zinpl(i,1) = (.5*tmp2          )*tmp1
            zinpl(i,2) = (        - .5*tmp3)*tmp1
            zinpl(i,3) = (.5*tmp2 -    tmp3)*tmp1
            zinpl(i,4) = (   tmp2 - .5*tmp3)*tmp1
            pinpl(i,1) = 0.5*(p1 + pnm(i,k2+1))
            pinpl(i,2) = 0.5*(p2 + pnm(i,k2  ))
            pinpl(i,3) = 0.5*(p1 + pnm(i,k2  ))
            pinpl(i,4) = 0.5*(p2 + pnm(i,k2+1))
         end do
         do 400 kn=1,4
            do i=1,plon
               u(i)     = uinpl(i,kn)*abs(plh2o(i,k2) - plh2o(i,k2+1))
               sqrtu(i) = sqrt(u(i))
               dw(i)    = abs(w(i,k2) - w(i,k2+1))
               pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
               ds2c     = abs(s2c(i,k2) - s2c(i,k2+1))
               uc1(i)   = uinpl(i,kn)*ds2c
               uc1(i)   = (uc1(i) + 1.7e-3*u(i))*(1. +  2.*uc1(i))/
     $                                           (1. + 15.*uc1(i))
               uc(i)    = uinpl(i,kn)*ds2c + 2.e-3*u(i)
            end do
            do i=1,plon
               dtx(i)      = temh2o(i,kn) - 250.
               dty(i)      = tbar(i,kn) - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
            end do
            do iband=2,4,2
               do i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
               end do
            end do
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C
            do i=1,plon
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac    = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4      = term1(i,2)*term4(i,2)
               t2t5      = term2(i,2)*term5(i,2)
               a         = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)    = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i)   = fwk(i)*u(i)
               rsum      = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
C               trab1(i) = rsum
            end do
C
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C
            do i=1,plon
               a41   = 1.75 - 3.960e-03*dtz(i)
               a51   = 1.00 + 1.3*sqrtu(i)
               a61   = 1.00 + 1.250e-03*dtp(i) + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
C               trab7(i) = rsum
            end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
            do k=1,2
               do i=1,plon
                  phi   = exp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = exp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*winpl(i,kn)*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(sqrt(g2) - 1.)
                  trline(i,k) = exp(-g4)
               end do
            end do
            do i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
            end do
C
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C
            do i=1,plon
               dtym10     = dty(i) - 10.
               denom      = 1. + (c30 + c31*dtym10*dtym10)*sqrtu(i)
               k21        = term7(i,1) + term8(i,1)/denom
               denom      = 1. + (c28 + c29*dtym10       )*sqrtu(i)
               k22        = term7(i,2) + term8(i,2)/denom
               term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                     (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                      (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
               tr1     = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2     = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5     = exp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
               tr6     = exp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)  = tr1*tr5
               tr10(i) = tr2*tr6
               trab2(i)= 0.65*tr9(i) + 0.35*tr10(i)
               th2o(i) = tr10(i)
               trab4(i)= exp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i)= exp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                     (1. + c9*dtx(i)*(1. + c11*dtx(i)*
     $                     (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
               abso(i,3)  = term6(i,2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4)  = term9(i,2)*.5*(tr1-tr9(i)+tr2-tr10(i))
            end do
C
C abso(i,5)  o3  9.6 micrometer (nu3 and nu1 bands)
C
            do i=1,plon
               te        = (tbar(i,kn)*r293)**.7
               dplos     = abs(plos(i,k2+1) - plos(i,k2))
               u1        = zinpl(i,kn)*18.29*dplos/te
               u2        = zinpl(i,kn)*.5649*dplos/te
               tlocal    = tbar(i,kn)
               tcrfac    = sqrt(tlocal*r250)*te
               beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
               realnu    = te/beta
               tmp1      = u1/sqrt(4. + u1*(1. + realnu))
               tmp2      = u2/sqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*log(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*o3emm(i,kn)*
     $                     (h2otr(i,k2+1)/h2otr(i,k2))
               to3(i)    = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
C               trab5(i) = 1.-(o3bndi/(1060-980.))
            end do
C
C abso(i,6)   co2 15  micrometer band system
C
            do 300 i=1,plon
               dplco2   = plco2(i,k2+1) - plco2(i,k2)
               sqwp     = sqrt(uinpl(i,kn)*dplco2)
               et       = exp(-480./tbar(i,kn))
               sqti(i)  = sqrt(tbar(i,kn))
               rsqti    = 1./sqti(i)
               et2      = et*et
               et4      = et2*et2
               omet     = (1. - 1.5*et2)
               f1co2    = 899.70*omet*
     $                    (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i)= f1co2*sqwp
               t1co2(i) = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme    = 1. - et2
               alphat   = oneme**3*rsqti
               pi       = abs(dpnm(i))*winpl(i,kn)
               wco2     = 2.5221*co2vmr*pi*rga
               u7(i)    = 4.9411e4*alphat*et2*wco2
               u8       = 3.9744e4*alphat*et4*wco2
               u9       = 1.0447e5*alphat*et4*et2*wco2
               u13      = 2.8388e3*alphat*et4*wco2
               tpath    = tbar(i,kn)
               tlocal   = tbar(i,kn)
               tcrfac   = sqrt((tlocal*r250)*(tpath*r300))
               posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
               rbeta7(i)= 1./(5.3228*posqt)
               rbeta8   = 1./(10.6576*posqt)
               rbeta9   = rbeta7(i)
               rbeta13  = rbeta9
               f2co2(i) = u7(i)/sqrt(4. + u7(i)*(1. + rbeta7(i))) +
     $                    u8/sqrt(4. + u8*(1. + rbeta8)) +
     $                    u9/sqrt(4. + u9*(1. + rbeta9))
               f3co2(i) = u13/sqrt(4. + u13*(1. + rbeta13))
               tmp1     = log(1. + f1sqwp(i))
               tmp2     = log(1. + f2co2(i))
               tmp3     = log(1. + f3co2(i))
               absbnd   = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
               abso(i,6)= trab2(i)*emm(i,kn)*absbnd
               tco2(i)=1.0/(1.0+ 10.0*u7(i)
     &                          /sqrt(4.+u7(i)*(1.+rbeta7(i))))
C               trab3(i) = 1. - bndfct*absbnd
  300       continue
c
c   Calculate trace gas absorptivity for nearest layer
c
            call trcabn(k2, kn, ucfc11, ucfc12, un2o0,  un2o1,
     $                          uch4,   uco211, uco212, uco213,
     $                          uco221, uco222, uco223, bn2o0,
     $                          bn2o1,  bch4,   tbar,   bplnk,
     $                          winpl,  pinpl,  tco2,   th2o,
     $                          to3,    uptype, dw,     s2c,
     $                          u,      pnew,   abstrc     ,uinpl)
C
C Total next layer absorptivity:
C
            do i=1,plon
               absnxt(i,k2,kn,jslc) =
     &                           abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
     $                           + abstrc(i)
            end do
  400    continue
  500 continue                  !  end of nearest layer level loop
C
      return
      end
c
c***************************** RADEMS *********************************
c
      subroutine radems(s2c     ,s2t     ,w       ,tplnke  ,plh2o   ,
     $                  pnm     ,plco2   ,tint    ,tint4   ,tlayr   ,
     $                  tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  ,
     $                  un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $                  bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  ,
     $                  co2t    ,h2otr   ,abplnk1 ,abplnk2 ,jslc  )
C-----------------------------------------------------------------------
C
C Compute emissivity for H2O, CO2, O3
C
C H2O  ....  Uses nonisothermal emissivity for water vapor from
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
C
C
C CO2  ....  Uses absorptance parameterization of the 15 micro-meter
C            (500 - 800 cm-1) band system of Carbon Dioxide, from
C            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
C            of the Absorptance Due to the 15 micro-meter Band System
C            of Carbon Dioxide Jouranl of Geophysical Research,
C            vol. 96., D5, pp 9013-9019
C
C O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
C            band system of ozone, from Ramanathan, V. and R. Dickinson,
C            1979: The Role of stratospheric ozone in the zonal and
C            seasonal radiative energy balance of the earth-troposphere
C            system. Journal of the Atmospheric Sciences, Vol. 36,
C            pp 1084-1104
C
C Computes individual emissivities, accounting for band overlap, and
C sums to obtain the total.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radems.F,v 1.2 1995/02/17 21:28:41 jhack Exp $
c $Author: jhack $
c
      implicit none
      integer jslc
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCAE.h>
#include <RADCON.h>
#include <RADBUF.h>
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real  s2c(plon,plevp),    ! H2o continuum path length
     $     s2t(plon,plevp),    ! Tmp and prs wghted h2o path length
     $     w(plon,plevp),      ! H2o path length
     $     tplnke(plon),       ! Layer planck temperature
     $     plh2o(plon,plevp),  ! H2o prs wghted path length
     $     pnm(plon,plevp),    ! Model interface pressure
     $     plco2(plon,plevp),  ! Prs wghted path of co2
     $     tint(plon,plevp),   ! Model interface temperatures
     $     tint4(plon,plevp),  ! Tint to the 4th power
     $     tlayr(plon,plevp),  ! K-1 model layer temperature
     $     tlayr4(plon,plevp), ! Tlayr to the 4th power
     $     plol(plon,plevp),   ! Pressure wghtd ozone path
     $     plos(plon,plevp)    ! Ozone path
c
c   Trace gas variables
c
      real  ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp),   ! pressure factor for ch4
     $     uptype(plon,plevp)  ! p-type continuum path length
C
C Output arguments (also in radbuf.cb)
C
      real  
     $     co2em(plon,plevp),   ! Layer co2 normalzd plnck funct drvtv
     $     co2eml(plon,plev),   ! Intrfc co2 normalzd plnck func drvtv
     $     co2t(plon,plevp),    ! Tmp and prs weighted path length
     $     h2otr(plon,plevp)    ! H2o transmission over o3 band
      real  emplnk(14,plon),        ! emissivity Planck factor
     $     abplnk1(14,plon,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plon,plevp)  ! nearest layer factor
      real  emstrc(plon,plevp)  ! total trace gas emissivity

C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                  ! Longitude index
     $     k,                  ! Level index]
     $     k1,                 ! Level index
     $     iband               ! H2o band index
C
C Local variables for H2O:
C
      real  h2oems(plon,plevp),! H2o emissivity
     $     tpathe(plon),      ! Used to compute h2o emissivity
     $     a(plon),           ! Eq(2) in table A3a of R&D
     $     corfac(plon),      ! Correction factors in table A3b
     $     dtp(plon),         ! Path temperature minus 300 K used in
C                                h2o rotation band absorptivity
     $     dtx(plon),         ! Planck temperature minus 250 K
     $     dty(plon),         ! Path temperature minus 250 K
     $     dtz(plon),         ! Planck temperature minus 300 K
     $     emis(plon,4),      ! Total emissivity (h2o+co2+o3)
     $     rsum(plon),        ! Eq(1) in table A2 of R&D
     $     term1(plon,4),     ! Equation(5) in table A3a of R&D(1986)
     $     term2(plon,4)      ! Delta a(Te) in table A3a of R&D(1986)
      real  term3(plon,4),     ! B(T) function for rotation and
C                                vibration-rotation band emissivity
     $     term4(plon,4),     ! Equation(6) in table A3a of R&D(1986)
     $     term5(plon,4),     ! Delta a(Tp) in table A3a of R&D(1986)
     $     term6(plon,2),     ! B(T) function for window region
     $     term7(plon,2),     ! Kl_inf(i) in eq(8) of table A3a of R&D
     $     term8(plon,2),     ! Delta kl_inf(i) in eq(8)
     $     term9(plon,2),     ! B(T) function for 500-800 cm-1 region
     $     tr1(plon),         ! Equation(6) in table A2 for 650-800
     $     tr2(plon),         ! Equation(6) in table A2 for 500-650
     $     tr3(plon)          ! Equation(4) in table A2 for 650-800
      real  tr4(plon),         ! Equation(4),table A2 of R&D for 500-650
     $     tr7(plon),         ! Equation (6) times eq(4) in table A2
C                              !   of R&D for 650-800 cm-1 region
     $     tr8(plon),         ! Equation (6) times eq(4) in table A2
C                              !   of R&D for 500-650 cm-1 region
     $     uc(plon),          ! Y + 0.002U in eq(8) of table A2 of R&D
     $     pnew(plon),        ! Effective pressure for h2o linewidth
     $     trline(plon,2),    ! Transmission due to H2O lines in window
     $     k21(plon),         ! Exponential coefficient used to calc
C                              !  rot band transmissivity in the 650-800
C                              !  cm-1 region (tr1)
     $     k22(plon),         ! Exponential coefficient used to calc
C                              !  rot band transmissivity in the 500-650
C                              !  cm-1 region (tr2)
     $     u(plon),           ! Pressure weighted H2O path length
     $     uc1(plon),         ! H2o continuum pathlength 500-800 cm-1
     $     r80257              ! Conversion factor for h2o pathlength
      real  a11,                ! A1 in table A3b for rotation band emiss
     $     a31,                ! A3 in table A3b for rotation band emiss
     $     a21,                ! First part in numerator of A2 table A3b
     $     a22,                ! Second part in numertor of A2 table A3b
     $     a23,                ! Denominator of A2 table A3b (rot band)
     $     t1t4,               ! Eq(3) in table A3a of R&D
     $     t2t5,               ! Eq(4) in table A3a of R&D
     $     fwk,                ! Equation(33) in R&D far wing correction
     $     a41,                ! Numerator in A2 in Vib-rot (table A3b)
     $     a51,                ! Denominator in A2 in Vib-rot(table A3b)
     $     a61,                ! A3 factor for Vib-rot band in table A3b
     $     phi,                ! Eq(11) in table A3a of R&D
     $     psi,                ! Eq(12) in table A3a of R&D
     $     ubar,               ! H2o scaled path comment eq(10) table A2
     $     g1,                 ! Part of eq(10) table A2
     $     pbar,               ! H2o scaled pres comment eq(10) table A2
     $     g3,                 ! Part of eq(10) table A2
     $     g2,                 ! Part of arguement in eq(10) in table A2
     $     g4,                 ! Arguement in exp() in eq(10) table A2
     $     cf812               ! Eq(11) in table A2 of R&D
      real  troco2(plon,plevp) ! H2o overlap factor for co2 absorption
C
C Local variables for CO2:
C
      real  co2ems(plon,plevp), ! Co2 emissivity
     $     co2plk(plon),       ! Used to compute co2 emissivity
     $     sum(plon),          ! Used to calculate path temperature
     $     t1i,                 ! Co2 hot band temperature factor
     $     sqti,                ! Sqrt of temperature
     $     pi,                  ! Pressure used in co2 mean line width
     $     et,                  ! Co2 hot band factor
     $     et2,                 ! Co2 hot band factor
     $     et4,                 ! Co2 hot band factor
     $     omet,                ! Co2 stimulated emission term
     $     ex                   ! Part of co2 planck function
      real  f1co2,               ! Co2 weak band factor
     $     f2co2,               ! Co2 weak band factor
     $     f3co2,               ! Co2 weak band factor
     $     t1co2,               ! Overlap factor weak bands strong band
     $     sqwp,                ! Sqrt of co2 pathlength
     $     f1sqwp,              ! Main co2 band factor
     $     oneme,               ! Co2 stimulated emission term
     $     alphat,              ! Part of the co2 stimulated emiss term
     $     wco2,                ! Consts used to define co2 pathlength
     $     posqt,               ! Effective pressure for co2 line width
     $     rbeta7,              ! Inverse of co2 hot band line width par
     $     rbeta8,              ! Inverse of co2 hot band line width par
     $     rbeta9,              ! Inverse of co2 hot band line width par
     $     rbeta13              ! Inverse of co2 hot band line width par
      real  tpath,               ! Path temp used in co2 band model
     $     tmp1,                ! Co2 band factor
     $     tmp2,                ! Co2 band factor
     $     tmp3,                ! Co2 band factor
     $     tlayr5,              ! Temperature factor in co2 Planck func
     $     rsqti,               ! Reciprocal of sqrt of temperature
     $     exm1sq               ! Part of co2 Planck function
      real  u7,    ! Absorber amount for various co2 band systems
     $     u8,    ! Absorber amount for various co2 band systems
     $     u9,    ! Absorber amount for various co2 band systems
     $     u13    ! Absorber amount for various co2 band systems
      real  r250,  ! Inverse 250K
     $     r300,  ! Inverse 300K
     $     rsslp  ! Inverse standard sea-level pressure
C
C Local variables for O3:
C
      real  o3ems(plon,plevp),  ! Ozone emissivity
     $     dbvtt(plon),        ! Tmp drvtv of planck fctn for tplnke
     $     te,                  ! Temperature factor
     $     u1,                  ! Path length factor
     $     u2,                  ! Path length factor
     $     phat,                ! Effecitive path length pressure
     $     tlocal,              ! Local planck function temperature
     $     tcrfac,              ! Scaled temperature factor
     $     beta,                ! Absorption funct factor voigt effect
     $     realnu,              ! Absorption function factor
     $     o3bndi               ! Band absorption factor
C
C Transmission terms for various spectral intervals:
C
C     real     trem1(plon)         ! H2o     0 -  800 cm-1
C     real     trem2(plon)         ! H2o   500 -  800 cm-1
C     real     trem3(plon)         ! Co2   500 -  800 cm-1
      real  trem4(plon)         ! H2o   800 - 1000 cm-1
C     real     trem5(plon)         ! O3     9.6 micro-meter band
      real  trem6(plon)         ! H2o  1000 - 1200 cm-1
C     real     trem7(plon)         ! H2o  1200 - 2200 cm-1
      real  bndfct,              ! Band absorptance parameter for co2
     $     absbnd               ! Proportional to co2 band absorptance
      real  tco2(plon),         ! co2 overlap factor
     $     th2o(plon),         ! h2o overlap factor
     $     to3(plon)           ! o3 overlap factor
C
C---------------------------Statement functions-------------------------
C
C Statement functions
C Derivative of planck function at 9.6 micro-meter wavelength, and
C an absorption function factor:
C
      real  dbvt,fo3,t,ux,vx
C
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.+(-6.1364820707e-3+1.5550319767e-5*t)*t)
C
      fo3(ux,vx)=ux/sqrt(4.+ux*(1.+vx))
C
C-----------------------------------------------------------------------
C
C Initialize
C
      r80257  = 1./8.0257e-04
C
      r250  = 1./250.
      r300  = 1./300.
      rsslp = 1./sslp
C
C Planck function for co2
C
      do i=1,plon
         ex        = exp(960./tplnke(i))
         co2plk(i) = 5.e8/((tplnke(i)**4)*(ex - 1.))
         co2t(i,1) = tplnke(i)
         sum(i)    = co2t(i,1)*pnm(i,1)
      end do
      k = 1
      do k1=plevp,2,-1
         k = k + 1
         do i=1,plon
            sum(i)         = sum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
            ex             = exp(960./tlayr(i,k1))
            tlayr5         = tlayr(i,k1)*tlayr4(i,k1)
            co2eml(i,k1-1) = 1.2e11*ex/(tlayr5*(ex - 1.)**2)
            co2t(i,k)      = sum(i)/pnm(i,k)
         end do
      end do
      bndfct = 2.*22.18/(sqrt(196.)*300.)
C
C Initialize planck function derivative for O3
C
      do i=1,plon
         dbvtt(i) = dbvt(tplnke(i))
      end do
c
c   Calculate trace gas Planck functions
c
      call trcplk(tint, tlayr, tplnke, emplnk, abplnk1, abplnk2)
C
C Interface loop
C
      do 200 k1=1,plevp
C
C H2O emissivity
C
C emis(i,1)     0 -  800 cm-1   rotation band
C emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
C emis(i,3)   800 - 1200 cm-1   window
C emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
C
C For the p type continuum
C
         do i=1,plon
            uc(i)     = s2c(i,k1) + 2.e-3*plh2o(i,k1)
            u(i)      = plh2o(i,k1)
            pnew(i)   = u(i)/w(i,k1)
C
C Apply scaling factor for 500-800 continuum
C
            uc1(i)    = (s2c(i,k1) + 1.7e-3*plh2o(i,k1))*
     $                 (1. + 2.*s2c(i,k1))/(1. + 15.*s2c(i,k1))
            tpathe(i) = s2t(i,k1)/plh2o(i,k1)
         end do
         do i=1,plon
            dtx(i) = tplnke(i) - 250.
            dty(i) = tpathe(i) - 250.
            dtz(i) = dtx(i) - 50.
            dtp(i) = dty(i) - 50.
         end do
         do iband=1,3,2
            do i=1,plon
               term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                          dtx(i)*(1. + c1(iband)*dtx(i))
               term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                          dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                  (1. + c3(iband)*dtx(i)))
               term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                          dtx(i)*(1. +  c4(iband)*dtx(i)*
     $                                  (1. + c5(iband)*dtx(i)))
               term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                          dty(i)*(1. + c6(iband)*dty(i))
               term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                          dty(i)*(1. + c7(iband)*dty(i))
            end do
         end do
C
C emis(i,1)     0 -  800 cm-1   rotation band
C
         do i=1,plon
            a11  = .37 - 3.33e-5*dtz(i) + 3.33e-6*dtz(i)*dtz(i)
            a31  = 1.07 - 1.00e-3*dtp(i) + 1.475e-5*dtp(i)*dtp(i)
            a21  = 1.3870 + 3.80e-3*dtz(i) - 7.8e-6*dtz(i)*dtz(i)
            a22  = 1.0 - 1.21e-3*dtp(i) - 5.33e-6*dtp(i)*dtp(i)
            a23  = 0.9 + 2.62*sqrt(u(i))
            corfac(i) = a31*(a11 + ((a21*a22)/a23))
            t1t4 = term1(i,1)*term4(i,1)
            t2t5 = term2(i,1)*term5(i,1)
            a(i) = t1t4 + t2t5/(1. + t2t5*sqrt(u(i))*corfac(i))
            fwk  = fwcoef + fwc1/(1. + fwc2*u(i))
            rsum(i)   = exp(-a(i)*(sqrt(u(i)) + fwk*u(i)))
            emis(i,1) = (1. - rsum(i))*term3(i,1)
C            trem1(i)  = rsum(i)
C
C emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
C
            a41      = 1.75 - 3.96e-3*dtz(i)
            a51      = 1.00 + 1.3*sqrt(u(i))
            a61      = 1.00 + 1.25e-3*dtp(i) + 6.25e-5*dtp(i)*dtp(i)
            corfac(i)= .3*(1. + (a41)/(a51))*a61
            t1t4     = term1(i,3)*term4(i,3)
            t2t5     = term2(i,3)*term5(i,3)
            a(i)     = t1t4 + t2t5/(1. + t2t5*sqrt(u(i))*corfac(i))
            fwk      = fwcoef + fwc1/(1. + fwc2*u(i))
            rsum(i)  = exp(-a(i)*(sqrt(u(i)) + fwk*u(i)))
            emis(i,2)= (1. - rsum(i))*term3(i,3)
C            trem7(i) = rsum(i)
         end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
         do k=1,2
            do i=1,plon
               phi  = a1(k)*(dty(i) + 15.) + a2(k)*(dty(i) + 15.)**2
               psi  = b1(k)*(dty(i) + 15.) + b2(k)*(dty(i) + 15.)**2
               phi  = exp(phi)
               psi  = exp(psi)
               ubar = w(i,k1)*phi
               ubar = (ubar*1.66)*r80257
               pbar = pnew(i)*(psi/phi)
               cf812 = cfa1 + ((1.-cfa1)/(1. + ubar*pbar*10.))
               g1   = (realk(k)*pbar)/(2.*st(k))
               g2   = 1. + (ubar*4.0*st(k)*cf812)/pbar
               g3   = sqrt(g2) - 1.
               g4   = g1*g3
               trline(i,k) = exp(-g4)
            end do
         end do
         do i=1,plon
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
         end do
C
C emis(i,3)   800 - 1200 cm-1   window
C
         do i=1,plon
            term6(i,1) = coeff(1,1) + coeff(2,1)*dtx(i)*
     $                  (1. +  c8*dtx(i)*(1. + c10*dtx(i)*
     $                  (1. + c12*dtx(i)*(1. + c14*dtx(i)))))
            trem4(i)  = exp(-(coefg(1,1)+coefg(2,1)*dtx(i))*uc(i))
     $                  *trline(i,2)
            trem6(i)  = exp(-(coefg(1,2)+coefg(2,2)*dtx(i))*uc(i))
     $                  *trline(i,1)
            emis(i,3) = term6(i,1)*(1. - .5*trem4(i) -.5*trem6(i))
C
C emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
C
            k21(i) = term7(i,1) + term8(i,1)/
     $           (1.+(c30+c31*(dty(i)-10.)*(dty(i)-10.))*sqrt(u(i)))
            k22(i) = term7(i,2) + term8(i,2)/
     $           (1. + (c28 + c29*(dty(i)-10.))*sqrt(u(i)))
            term9(i,1) = coefi(1,1) + coefi(2,1)*dtx(i)*
     $                  (1. + c18*dtx(i)*(1. + c20*dtx(i)*
     $                   (1. + c22*dtx(i)*(1. + c24*dtx(i)))))
            fwk    = fwcoef + fwc1/(1.+fwc2*u(i))
            tr1(i) = exp(-(k21(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
            tr2(i) = exp(-(k22(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
            tr3(i) = exp(-((coefh(1,1) + coefh(2,1)*dtx(i))*uc1(i)))
            tr4(i) = exp(-((coefh(1,2) + coefh(2,2)*dtx(i))*uc1(i)))
            tr7(i) = tr1(i)*tr3(i)
            tr8(i) = tr2(i)*tr4(i)
            emis(i,4) = term9(i,1)*.5*(tr1(i)-tr7(i) + tr2(i)-tr8(i))
            h2oems(i,k1) = emis(i,1)+emis(i,2)+emis(i,3)+emis(i,4)
            troco2(i,k1) = 0.65*tr7(i) + 0.35*tr8(i)
            th2o(i) = tr8(i)
C            trem2(i)     = troco2(i,k1)
         end do
C
C CO2 emissivity for 15 micron band system
C
         do 100 i=1,plon
            t1i    = exp(-480./co2t(i,k1))
            sqti   = sqrt(co2t(i,k1))
            rsqti  = 1./sqti
            et     = t1i
            et2    = et*et
            et4    = et2*et2
            omet   = 1. - 1.5*et2
            f1co2  = 899.70*omet*(1. + 1.94774*et + 4.73486*et2)*rsqti
            sqwp   = sqrt(plco2(i,k1))
            f1sqwp = f1co2*sqwp
            t1co2  = 1./(1. + 245.18*omet*sqwp*rsqti)
            oneme  = 1. - et2
            alphat = oneme**3*rsqti
            wco2   = 2.5221*co2vmr*pnm(i,k1)*rga
            u7     = 4.9411e4*alphat*et2*wco2
            u8     = 3.9744e4*alphat*et4*wco2
            u9     = 1.0447e5*alphat*et4*et2*wco2
            u13    = 2.8388e3*alphat*et4*wco2
C
            tpath  = co2t(i,k1)
            tlocal = tplnke(i)
            tcrfac = sqrt((tlocal*r250)*(tpath*r300))
            pi     = pnm(i,k1)*rsslp + 2.*dpfco2*tcrfac
            posqt  = pi/(2.*sqti)
            rbeta7 =  1./( 5.3288*posqt)
            rbeta8 = 1./ (10.6576*posqt)
            rbeta9 = rbeta7
            rbeta13= rbeta9
            f2co2  = (u7/sqrt(4. + u7*(1. + rbeta7))) +
     $               (u8/sqrt(4. + u8*(1. + rbeta8))) +
     $               (u9/sqrt(4. + u9*(1. + rbeta9)))
            f3co2  = u13/sqrt(4. + u13*(1. + rbeta13))
            tmp1   = log(1. + f1sqwp)
            tmp2   = log(1. +  f2co2)
            tmp3   = log(1. +  f3co2)
            absbnd = (tmp1 + 2.*t1co2*tmp2 + 2.*tmp3)*sqti
            tco2(i)=1.0/(1.0+10.0*(u7/sqrt(4. + u7*(1. + rbeta7))))
            co2ems(i,k1)  = troco2(i,k1)*absbnd*co2plk(i)
            ex     = exp(960./tint(i,k1))
            exm1sq = (ex - 1.)**2
            co2em(i,k1) = 1.2e11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
C            trem3(i) = 1. - bndfct*absbnd
  100    continue
C
C O3 emissivity
C
         do i=1,plon
            h2otr(i,k1) = exp(-12.*s2c(i,k1))
            te          = (co2t(i,k1)/293.)**.7
            u1          = 18.29*plos(i,k1)/te
            u2          = .5649*plos(i,k1)/te
            phat        = plos(i,k1)/plol(i,k1)
            tlocal      = tplnke(i)
            tcrfac      = sqrt(tlocal*r250)*te
            beta        = (1./.3205)*((1./phat) + (dpfo3*tcrfac))
            realnu      = (1./beta)*te
            o3bndi      = 74.*te*(tplnke(i)/375.)*
     $         log(1. + fo3(u1,realnu) + fo3(u2,realnu))
            o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
            to3(i)=1.0/(1. + 0.1*fo3(u1,realnu) + 0.1*fo3(u2,realnu))
C            trem5(i)    = 1.-(o3bndi/(1060-980.))
         end do
c
c   Calculate trace gas emissivities
c
      call trcems(k1,     co2t,   pnm,    ucfc11, ucfc12, un2o0,
     $            un2o1,  bn2o0,  bn2o1,  uch4,   bch4,   uco211,
     $            uco212, uco213, uco221, uco222, uco223, uptype,
     $            w,      s2c,    u,      emplnk, th2o,   tco2,
     $            to3,    emstrc)
C
C Total emissivity:
C
         do i=1,plon
            emstot(i,k1,jslc) = h2oems(i,k1) + co2ems(i,k1)
     &                         + o3ems(i,k1) + emstrc(i,k1)
         end do
  200 continue          ! End of interface loop
C
      return
      end
c
c***************************** RADOZ2 *********************************
c
      subroutine radoz2(o3vmr   ,pint    ,plol    ,plos    )
C-----------------------------------------------------------------------
C
C Computes the path length integrals to the model interfaces given the
C ozone volume mixing ratio
C
C---------------------------Code history--------------------------------
C
C Original version:     CCM1
C Standardized:         J. Rosinski, June 1992
C Reviewed:             J. Kiehl, B. Briegleb, August 1992
C Mixing ratio version: Bruce Biegleb, September 1992
C
C-----------------------------------------------------------------------
c
c $Id: radoz2.F,v 1.1.1.1 1995/02/09 23:27:03 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C-----------------------------------------------------------------------
c
#include <RADSIZE.h>
c
C------------------------------Commons----------------------------------
#include <RADOZP.h>
C------------------------------Input arguments--------------------------
C
      real  o3vmr(plon,plev)   ! ozone volume mixing ratio
      real  pint(plon,plevp)   ! Model interface pressures
C
C----------------------------Output arguments---------------------------
C
      real  plol(plon,plevp),  ! Ozone prs weighted path length (cm)
     $     plos(plon,plevp)   ! Ozone path length (cm)
C
C---------------------------Local workspace-----------------------------
C
      integer   i,             ! longitude index
     $          k              ! level index
C
C-----------------------------------------------------------------------
C
C Evaluate the ozone path length integrals to interfaces;
C factors of .1 and .01 to convert pressures from cgs to mks:
C
C Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
C
      do i=1,plon
         plos(i,1) = 0.1*cplos*o3vmr(i,1)*pint(i,1)
         plol(i,1) = 0.01*cplol*o3vmr(i,1)*pint(i,1)*pint(i,1)
      end do
      do k=2,plevp
         do i=1,plon
            plos(i,k) = plos(i,k-1) + 0.1*cplos*o3vmr(i,k-1)*
     $                  (pint(i,k) - pint(i,k-1))
            plol(i,k) = plol(i,k-1) + 0.01*cplol*o3vmr(i,k-1)*
     $                  (pint(i,k)*pint(i,k)-pint(i,k-1)*pint(i,k-1))
         end do
      end do
C
      return
      end
c
c***************************** RADTPL *********************************
c
      subroutine radtpl(tnm    ,ts     ,qnm    ,pnm    ,plh2o  ,
     $                  tplnka ,s2c    ,s2t    ,w      ,tplnke ,
     $                  tint   ,tint4  ,tlayr  ,tlayr4 ,pmln   ,
     $                  piln   )
C-----------------------------------------------------------------------
C
C Compute temperatures and path lengths for longwave radiation
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radtpl.F,v 1.1.1.1 1995/02/09 23:27:03 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCON.h>
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real  tnm(plon,plev),     ! Model level temperatures
     $     ts(plon),           ! Surface skin temperature
     $     qnm(plon,plev),     ! Model level specific humidity
     $     pnm(plon,plevp),    ! Pressure at model interfaces (dynes/cm2)
     $     plh2o(plon,plevp)   ! Pressure weighted h2o path
C
C Output arguments
C
      real  tplnka(plon,plevp), ! Level temperature from interface temperatures
     $     s2c(plon,plevp),    ! H2o continuum path length
     $     s2t(plon,plevp),    ! H2o tmp and prs wghtd path length
     $     w(plon,plevp),      ! H2o path length
     $     tplnke(plon),       ! Equal to tplnka
     $     tint(plon,plevp),   ! Layer interface temperature
     $     tint4(plon,plevp),  ! Tint to the 4th power
     $     tlayr(plon,plevp),  ! K-1 level temperature
     $     tlayr4(plon,plevp), ! Tlayr to the 4th power
     $     pmln(plon,plev),    ! Ln(pmidm1)
     $     piln(plon,plevp)    ! Ln(pintm1)
C
C---------------------------Local variables-----------------------------
C
      integer   i,              ! Longitude index
     $          k               ! Level index
      real  r296,                ! Inverse stand temp for h2o continuum
     $     repsil,              ! Inver ratio mol weight h2o to dry air
     $         dy,              ! Thickness of layer for tmp interp
     $       dpnm,              ! Pressure thickness of layer
     $     dpnmsq,              ! Prs squared difference across layer
     $       rtnm               ! Inverse level temperature
C
C-----------------------------------------------------------------------
C
      r296   = 1./296.
      repsil = 1./epsilo
C
C Set the top and bottom intermediate level temperatures,
C top level planck temperature and top layer temp**4.
C
C Tint is lower interface temperature
C (not available for bottom layer, so use ground temperature)
C
      do i=1,plon
         tint(i,plevp)  = ts(i)
         tint4(i,plevp) = tint(i,plevp)**4
         tplnka(i,1)    = tnm(i,1)
         tint(i,1)      = tplnka(i,1)
         tlayr4(i,1)    = tplnka(i,1)**4
         tint4(i,1)     = tlayr4(i,1)
      end do
C
C Intermediate level temperatures are computed using temperature
C at the full level below less dy*delta t,between the full level
C
      do k=2,plev
         do i=1,plon
            dy = (piln(i,k) - pmln(i,k))/(pmln(i,k-1) - pmln(i,k))
            tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
            tint4(i,k) = tint(i,k)**4
         end do
      end do
C
C Now set the layer temp=full level temperatures and establish a
C planck temperature for absorption (tplnka) which is the average
C the intermediate level temperatures.  Note that tplnka is not
C equal to the full level temperatures.
C
      do k=2,plevp
         do i=1,plon
            tlayr(i,k)  = tnm(i,k-1)
            tlayr4(i,k) = tlayr(i,k)**4
            tplnka(i,k) = .5*(tint(i,k) + tint(i,k-1))
         end do
      end do
C
C Calculate tplank for emissivity calculation.
C Assume isothermal tplnke i.e. all levels=ttop.
C
      do i=1,plon
         tplnke(i)  = tplnka(i,1)
         tlayr(i,1) = tint(i,1)
      end do
C
C Now compute h2o path fields:
C
      do i=1,plon
         s2t(i,1) = plh2o(i,1) * tnm(i,1)
C ccm3.2
C        w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
C        s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil

C ccm3.6.6
         w(i,1)   = sslp * (plh2o(i,1)*2.) / pnm(i,1)
         rtnm     = 1./tnm(i,1)
         s2c(i,1) = plh2o(i,1)*
     $            exp(1800.*(rtnm - r296))*qnm(i,1)*repsil
      end do
      do k=1,plev
         do i=1,plon
            dpnm       = pnm(i,k+1) - pnm(i,k)
            dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
            rtnm       = 1./tnm(i,k)
            s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
            w(i,k+1)   = w(i,k)   + rga*qnm(i,k)*dpnm
            s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)*
     $             exp(1800.*(rtnm - r296))*qnm(i,k)*repsil
         end do
      end do
C
      return
      end
c
c***************************** TRCAB *********************************
c
      subroutine trcab(k1, k2, ucfc11, ucfc12, un2o0,  un2o1,
     $                         uch4,   uco211, uco212, uco213,
     $                         uco221, uco222, uco223, bn2o0,
     $                         bn2o1,  bch4,   to3co2, pnm,
     $                         dw,     pnew,   s2c,    uptype,
     $                         dplh2o, abplnk1,tco2,   th2o,
     $                         to3,    abstrc)
c----------------------------------------------------------------------
c Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
c CFC12.
c
c             Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcab.F,v 1.2 1995/02/17 21:28:52 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCON.h>
C
C------------------------------Arguments--------------------------------
      integer k1,k2
      real  to3co2(plon),       ! pressure weighted temperature
     $     pnm(plon,plevp),    ! interface pressures
     $     ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp)    ! pressure factor for ch4
      real  dw(plon),           ! h2o path length
     $     pnew(plon),         ! pressure
     $     s2c(plon,plevp),    ! continuum path length
     $     uptype(plon,plevp), ! p-type h2o path length
     $     dplh2o(plon)        ! p squared h2o path length
      real  abplnk1(14,plon,plevp) ! Planck factor
      real  tco2(plon),         ! co2 transmission factor
     $     th2o(plon),         ! h2o transmission factor
     $     to3(plon)           ! o3 transmission factor
c
c  Output Arguments
c
      real  abstrc(plon)  ! total trace gas absorptivity
c
c  Local Variables
c
      real  sqti(plon),         ! square root of mean temp
     $     du1,                 ! cfc11 path length
     $     du2,                 ! cfc12 path length
     $     acfc1,               ! cfc11 absorptivity 798 cm-1
     $     acfc2,               ! cfc11 absorptivity 846 cm-1
     $     acfc3,               ! cfc11 absorptivity 933 cm-1
     $     acfc4,               ! cfc11 absorptivity 1085 cm-1
     $     acfc5,               ! cfc12 absorptivity 889 cm-1
     $     acfc6,               ! cfc12 absorptivity 923 cm-1
     $     acfc7,               ! cfc12 absorptivity 1102 cm-1
     $     acfc8,               ! cfc12 absorptivity 1161 cm-1
     $     du01,                ! n2o path length
     $     dbeta01,             ! n2o pressure factor
     $     dbeta11,             !         "
     $     an2o1,               ! absorptivity of 1285 cm-1 n2o band
     $     du02,                ! n2o path length
     $     dbeta02,             ! n2o pressure factor
     $     an2o2,               ! absorptivity of 589 cm-1 n2o band
     $     du03                 ! n2o path length
      real  dbeta03,             ! n2o pressure factor
     $     an2o3,               ! absorptivity of 1168 cm-1 n2o band
     $     duch4,               ! ch4 path length
     $     dbetac,              ! ch4 pressure factor
     $     ach4,                ! absorptivity of 1306 cm-1 ch4 band
     $     du11,                ! co2 path length
     $     du12,                !       "
     $     du13,                !       "
     $     dbetc1,              ! co2 pressure factor
     $     dbetc2,              ! co2 pressure factor
     $     aco21,               ! absorptivity of 1064 cm-1 band
     $     du21,                ! co2 path length
     $     du22,                !       "
     $     du23,                !       "
     $     aco22                ! absorptivity of 961 cm-1 band
      real  tt(plon),           ! temp. factor for h2o overlap factor
     $     psi1,                !                 "
     $     phi1,                !                 "
     $     p1,                  ! h2o overlap factor
     $     w1,                  !        "
     $     ds2c(plon),         ! continuum path length
     $     duptyp(plon),       ! p-type path length
     $     tw(plon,6),         ! h2o transmission factor
     $     g1(6),               !         "
     $     g2(6),               !         "
     $     g3(6),               !         "
     $     g4(6),               !         "
     $     ab(6),               ! h2o temp. factor
     $     bb(6),               !         "
     $     abp(6),              !         "
     $     bbp(6)               !         "
      real  tcfc3,               ! transmission for cfc11 band
     $     tcfc4,               ! transmission for cfc11 band
     $     tcfc6,               ! transmission for cfc12 band
     $     tcfc7,               ! transmission for cfc12 band
     $     tcfc8,               ! transmission for cfc12 band
     $     tlw,                 ! h2o transmission
     $     tch4                 ! ch4 transmission
c
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,
     $         0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,
     $         0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
c------------------------------------------------------------------
      real  func, u, b
      func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = sqrt(to3co2(i))
c h2o transmission
         tt(i) = abs(to3co2(i) - 250.0)
         ds2c(i) = abs(s2c(i,k1) - s2c(i,k2))
         duptyp(i) = abs(uptype(i,k1) - uptype(i,k2))
      end do
c
      do l = 1,6
            do i = 1,plon
               psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
               phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
               p1 = pnew(i) * (psi1/phi1) / sslp
               w1 = dw(i) * phi1
               tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                     - g3(l)*ds2c(i)-g4(l)*duptyp(i))
            end do
      end do
c
      do i = 1,plon
            du1 = abs(ucfc11(i,k1) - ucfc11(i,k2))
            du2 = abs(ucfc12(i,k1) - ucfc12(i,k2))
c cfc transmissions
            tcfc3 = exp(-175.005*du1)
            tcfc4 = exp(-1202.18*du1)
            tcfc6 = exp(-5786.73*du2)
            tcfc7 = exp(-2873.51*du2)
            tcfc8 = exp(-2085.59*du2)
c  Absorptivity for CFC11 bands
            acfc1 = 50.0*(1.0-exp(-54.09*du1))*tw(i,1)*abplnk1(7,i,k2)
            acfc2 = 60.0*(1.0-exp(-5130.03*du1))*tw(i,2)
     $                                            *abplnk1(8,i,k2)
            acfc3 = 60.0*(1.0-tcfc3) * tw(i,4)*tcfc6*abplnk1(9,i,k2)
            acfc4 = 100.0*(1.0-tcfc4)* tw(i,5) * abplnk1(10,i,k2)
c  Absorptivity for CFC12 bands
            acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)*
     $                                             abplnk1(11,i,k2)
            acfc6 = 50.0*(1.0 - tcfc6)* tw(i,4)*abplnk1(12,i,k2)
            acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5)*tcfc4*abplnk1(13,i,k2)
            acfc8 = 70.0*(1.0 - tcfc8)* tw(i,6)*abplnk1(14,i,k2)
c  Emissivity for CH4 band 1306 cm-1
            tlw = exp(-1.0*sqrt(dplh2o(i)))
            duch4 = abs(uch4(i,k1) - uch4(i,k2))
            dbetac = abs(bch4(i,k1) - bch4(i,k2))/duch4
            ach4 = 6.00444*sqti(i)*log(1.0 + func(duch4,dbetac)) *
     $           tlw * abplnk1(3,i,k2)
            tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
c  Absorptivity for N2O bands
            du01 = abs(un2o0(i,k1) - un2o0(i,k2))
            du11 = abs(un2o1(i,k1) - un2o1(i,k2))
            dbeta01 = abs(bn2o0(i,k1) - bn2o0(i,k2))/du01
            dbeta11 = abs(bn2o1(i,k1) - bn2o1(i,k2))/du11
c     1285 cm-1 band
            an2o1 = 2.35558*sqti(i)*log(1.0 + func(du01,dbeta01)
     $            +  func(du11,dbeta11)) * tlw* tch4*abplnk1(4,i,k2)
            du02 = 0.100090*du01
            du12 = 0.0992746*du11
            dbeta02 = 0.964282*dbeta01
c     589 cm-1 band
            an2o2 = 2.65581*sqti(i)*log(1.0 + func(du02,dbeta02)
     $            +  func(du12,dbeta02)) * th2o(i) * tco2(i) *
     $               abplnk1(5,i,k2)
            du03 = 0.0333767*du01
            dbeta03 = 0.982143*dbeta01
c     1168 cm-1 band
            an2o3 = 2.54034*sqti(i)*log(1.0 + func(du03,dbeta03)) *
     $           tw(i,6) * tcfc8 * abplnk1(6,i,k2)
c  Emissivity for 1064 cm-1 band of CO2
            du11 = abs(uco211(i,k1) - uco211(i,k2))
            du12 = abs(uco212(i,k1) - uco212(i,k2))
            du13 = abs(uco213(i,k1) - uco213(i,k2))
            dbetc1 = 2.97558*abs(pnm(i,k1) + pnm(i,k2))/
     $                                           (2.0*sslp*sqti(i))
            dbetc2 = 2.0 * dbetc1
            aco21 = 3.7571*sqti(i)*log(1.0 + func(du11,dbetc1)
     $         + func(du12,dbetc2) + func(du13,dbetc2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * abplnk1(2,i,k2)
c  Emissivity for 961 cm-1 band
            du21 = abs(uco221(i,k1) - uco221(i,k2))
            du22 = abs(uco222(i,k1) - uco222(i,k2))
            du23 = abs(uco223(i,k1) - uco223(i,k2))
            aco22 = 3.8443*sqti(i)*log(1.0 + func(du21,dbetc1)
     $         + func(du22,dbetc1) + func(du23,dbetc2))
     $         * tw(i,4) * tcfc3 * tcfc6 * abplnk1(1,i,k2)
c total trace gas absorptivity
            abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6
     $            +  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4
     $            +  aco21 + aco22
      end do
      return
      end
c
c***************************** TRCABN *********************************
c
      subroutine trcabn(k2, kn, ucfc11, ucfc12, un2o0,  un2o1,
     $                          uch4,   uco211, uco212, uco213,
     $                          uco221, uco222, uco223, bn2o0,
     $                          bn2o1,  bch4,   tbar,   bplnk,
     $                          winpl,  pinpl,  tco2,   th2o,
     $                          to3,    uptype, dw,     s2c,
     $                          up2,    pnew,   abstrc     ,uinpl)
c----------------------------------------------------------------------
c Calculate nearest layer absorptivity due to CH4, N2O, CFC11 and CFC12
c
c         Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcabn.F,v 1.2 1995/02/17 21:28:54 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
#include <RADSIZE.h>
C------------------------------Commons----------------------------------
#include <RADCON.h>
C------------------------------Arguments--------------------------------
      integer k2,
     $        kn
      real  tbar(plon,4),       ! pressure weighted temperature
     $     ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp),   ! pressure factor for ch4
     $     bplnk(14,plon,4),   ! weighted Planck function for absorptivity
     $     winpl(plon,4),      ! fractional path length
     $     pinpl(plon,4)       ! pressure factor for subdivided layer
      real  tco2(plon),         ! co2 transmission
     $     th2o(plon),         ! h2o transmission
     $     to3(plon)           ! o3 transmission
      real  dw(plon),           ! h2o path length
     $     pnew(plon),         ! pressure factor
     $     s2c(plon,plevp),    ! h2o continuum factor
     $     uptype(plon,plevp), ! p-type path length
     $     up2(plon)           ! p squared path length
      real  uinpl(plon,4)    ! Nearest layer subdivision factor
c
c  Output Arguments
c
      real  abstrc(plon)        ! total trace gas absorptivity
c
c  Local Variables
c
      real  sqti(plon),         ! square root of mean temp
     $     rsqti(plon),        ! reciprocal of sqti
     $     du1,                 ! cfc11 path length
     $     du2,                 ! cfc12 path length
     $     acfc1,               ! absorptivity of cfc11 798 cm-1 band
     $     acfc2,               ! absorptivity of cfc11 846 cm-1 band
     $     acfc3,               ! absorptivity of cfc11 933 cm-1 band
     $     acfc4,               ! absorptivity of cfc11 1085 cm-1 band
     $     acfc5,               ! absorptivity of cfc11 889 cm-1 band
     $     acfc6,               ! absorptivity of cfc11 923 cm-1 band
     $     acfc7,               ! absorptivity of cfc11 1102 cm-1 band
     $     acfc8,               ! absorptivity of cfc11 1161 cm-1 band
     $     du01,                ! n2o path length
     $     dbeta01,             ! n2o pressure factors
     $     dbeta11              !        "
      real   an2o1,              ! absorptivity of the 1285 cm-1 n2o band
     $     du02,                ! n2o path length
     $     dbeta02,             ! n2o pressure factor
     $     an2o2,               ! absorptivity of the 589 cm-1 n2o band
     $     du03,                ! n2o path length
     $     dbeta03,             ! n2o pressure factor
     $     an2o3,               ! absorptivity of the 1168 cm-1 n2o band
     $     duch4,               ! ch4 path length
     $     dbetac,              ! ch4 pressure factor
     $     ach4,                ! absorptivity of the 1306 cm-1 ch4 band
     $     du11,                ! co2 path length
     $     du12,                !       "
     $     du13,                !       "
     $     dbetc1,              ! co2 pressure factor
     $     dbetc2,              ! co2 pressure factor
     $     aco21,               ! absorptivity of the 1064 cm-1 co2 band
     $     du21,                ! co2 path length
     $     du22,                !       "
     $     du23,                !       "
     $     aco22                ! absorptivity of the 961 cm-1 co2 band
      real  tt(plon),           ! temp. factor for h2o overlap
     $     psi1,                !          "
     $     phi1,                !          "
     $     p1,                  ! factor for h2o overlap
     $     w1,                  !          "
     $     ds2c(plon),         ! continuum path length
     $     duptyp(plon),       ! p-type path length
     $     tw(plon,6),         ! h2o transmission overlap
     $     g1(6),               ! h2o overlap factor
     $     g2(6),               !         "
     $     g3(6),               !         "
     $     g4(6),               !         "
     $     ab(6),               ! h2o temp. factor
     $     bb(6),               !         "
     $     abp(6),              !         "
     $     bbp(6)               !         "
      real  tcfc3,               ! transmission of cfc11 band
     $     tcfc4,               ! transmission of cfc11 band
     $     tcfc6,               ! transmission of cfc12 band
     $     tcfc7,               !         "
     $     tcfc8,               !         "
     $     tlw,                 ! h2o transmission
     $     tch4                 ! ch4 transmission
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,
     $         0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,
     $         0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
c------------------------------------------------------------------
      real  func, u, b
      func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = sqrt(tbar(i,kn))
         rsqti(i) = 1. / sqti(i)
c h2o transmission
         tt(i) = abs(tbar(i,kn) - 250.0)
         ds2c(i) = abs(s2c(i,k2+1) - s2c(i,k2))*uinpl(i,kn)
         duptyp(i) = abs(uptype(i,k2+1)-uptype(i,k2))*uinpl(i,kn)
      end do
c
      do l = 1,6
         do i = 1,plon
            psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
            phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
            p1 = pnew(i) * (psi1/phi1) / sslp
            w1 = dw(i) * winpl(i,kn) * phi1
            tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                  - g3(l)*ds2c(i)-g4(l)*duptyp(i))
         end do
      end do
c
      do i = 1,plon
c
            du1 = abs(ucfc11(i,k2+1) - ucfc11(i,k2)) * winpl(i,kn)
            du2 = abs(ucfc12(i,k2+1) - ucfc12(i,k2)) * winpl(i,kn)
c cfc transmissions
            tcfc3 = exp(-175.005*du1)
            tcfc4 = exp(-1202.18*du1)
            tcfc6 = exp(-5786.73*du2)
            tcfc7 = exp(-2873.51*du2)
            tcfc8 = exp(-2085.59*du2)
c  Absorptivity for CFC11 bands
            acfc1 = 50.0*(1.0-exp(-54.09*du1)) * tw(i,1)*bplnk(7,i,kn)
            acfc2 = 60.0*(1.0-exp(-5130.03*du1))*tw(i,2)*bplnk(8,i,kn)
            acfc3 = 60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6 * bplnk(9,i,kn)
            acfc4 = 100.0*(1.0 - tcfc4)* tw(i,5) * bplnk(10,i,kn)
c  Absorptivity for CFC12 bands
            acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)
     $                                            *bplnk(11,i,kn)
            acfc6 = 50.0*(1.0 - tcfc6)*tw(i,4)*bplnk(12,i,kn)
            acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5)*tcfc4 *bplnk(13,i,kn)
            acfc8 = 70.0*(1.0 - tcfc8)*tw(i,6)*bplnk(14,i,kn)
c  Emissivity for CH4 band 1306 cm-1
            tlw = exp(-1.0*sqrt(up2(i)))
            duch4 = abs(uch4(i,k2+1) - uch4(i,k2)) * winpl(i,kn)
            dbetac = 2.94449 * pinpl(i,kn) * rsqti(i) / sslp
            ach4 = 6.00444*sqti(i)*log(1.0 + func(duch4,dbetac)) *
     $           tlw * bplnk(3,i,kn)
            tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
c  Absorptivity for N2O bands
            du01 = abs(un2o0(i,k2+1) - un2o0(i,k2)) * winpl(i,kn)
            du11 = abs(un2o1(i,k2+1) - un2o1(i,k2)) * winpl(i,kn)
            dbeta01 = 19.399 *  pinpl(i,kn) * rsqti(i) / sslp
            dbeta11 = dbeta01
c     1285 cm-1 band
            an2o1 = 2.35558*sqti(i)*log(1.0 + func(du01,dbeta01)
     $            + func(du11,dbeta11)) * tlw * tch4 * bplnk(4,i,kn)
            du02 = 0.100090*du01
            du12 = 0.0992746*du11
            dbeta02 = 0.964282*dbeta01
c     589 cm-1 band
            an2o2 = 2.65581*sqti(i)*log(1.0 + func(du02,dbeta02)
     $             +  func(du12,dbeta02)) * tco2(i) * th2o(i) *
     $                bplnk(5,i,kn)
            du03 = 0.0333767*du01
            dbeta03 = 0.982143*dbeta01
c     1168 cm-1 band
            an2o3 = 2.54034*sqti(i)*log(1.0 + func(du03,dbeta03)) *
     $           tw(i,6) * tcfc8 * bplnk(6,i,kn)
c  Emissivity for 1064 cm-1 band of CO2
            du11 = abs(uco211(i,k2+1) - uco211(i,k2)) * winpl(i,kn)
            du12 = abs(uco212(i,k2+1) - uco212(i,k2)) * winpl(i,kn)
            du13 = abs(uco213(i,k2+1) - uco213(i,k2)) * winpl(i,kn)
            dbetc1 = 2.97558 * pinpl(i,kn) * rsqti(i) / sslp
            dbetc2 = 2.0 * dbetc1
            aco21 = 3.7571*sqti(i)*log(1.0 + func(du11,dbetc1)
     $         + func(du12,dbetc2) + func(du13,dbetc2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * bplnk(2,i,kn)
c  Emissivity for 961 cm-1 band of co2
            du21 = abs(uco221(i,k2+1) - uco221(i,k2)) * winpl(i,kn)
            du22 = abs(uco222(i,k2+1) - uco222(i,k2)) * winpl(i,kn)
            du23 = abs(uco223(i,k2+1) - uco223(i,k2)) * winpl(i,kn)
            aco22 = 3.8443*sqti(i)*log(1.0 + func(du21,dbetc1)
     $         + func(du22,dbetc1) + func(du23,dbetc2))
     $         * tw(i,4) * tcfc3 * tcfc6 * bplnk(1,i,kn)
c total trace gas absorptivity
            abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6
     $                + acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4
     $                + aco21 + aco22
      end do
      return
      end
c
c***************************** TRCEMS *********************************
c
      subroutine trcems(k,      co2t,   pnm,    ucfc11, ucfc12, un2o0,
     $                  un2o1,  bn2o0,  bn2o1,  uch4,   bch4,   uco211,
     $                  uco212, uco213, uco221, uco222, uco223, uptype,
     $                  w,      s2c,    up2,    emplnk, th2o,   tco2,
     $                  to3,    emstrc )
c----------------------------------------------------------------------
c  Calculate emissivity for CH4, N2O, CFC11 and CFC12 bands.
c            Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcems.F,v 1.2 1995/02/17 21:28:56 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCON.h>
C
C------------------------------Arguments--------------------------------
      real  co2t(plon,plevp),   ! pressure weighted temperature
     $     pnm(plon,plevp),    ! interface pressure
     $     ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     uptype(plon,plevp), ! continuum path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp)    ! pressure factor for ch4
      real  emplnk(14,plon)     ! emissivity Planck factor
      real  th2o(plon),         ! water vapor overlap factor
     $     tco2(plon),         ! co2 overlap factor
     $     to3(plon)           ! o3 overlap factor
      real  s2c(plon,plevp),    ! h2o continuum path length
     $     w(plon,plevp),      ! h2o path length
     $     up2(plon)           ! pressure squared h2o path length
      integer k                 ! level index
c
c  Output Arguments
c
      real  emstrc(plon,plevp)  ! total trace gas emissivity
c
c  Local Variables
c
      real  sqti(plon),         ! square root of mean temp
     $     ecfc1,               ! emissivity of cfc11 798 cm-1 band
     $     ecfc2,               !     "      "    "   846 cm-1 band
     $     ecfc3,               !     "      "    "   933 cm-1 band
     $     ecfc4,               !     "      "    "   1085 cm-1 band
     $     ecfc5,               !     "      "  cfc12 889 cm-1 band
     $     ecfc6,               !     "      "    "   923 cm-1 band
     $     ecfc7,               !     "      "    "   1102 cm-1 band
     $     ecfc8,               !     "      "    "   1161 cm-1 band
     $     u01,                 ! n2o path length
     $     u11,                 ! n2o path length
     $     beta01,              ! n2o pressure factor
     $     beta11,              ! n2o pressure factor
     $     en2o1,               ! emissivity of the 1285 cm-1 N2O band
     $     u02,                 ! n2o path length
     $     u12,                 ! n2o path length
     $     beta02,              ! n2o pressure factor
     $     en2o2,               ! emissivity of the 589 cm-1 N2O band
     $     u03                  ! n2o path length
      real  beta03,              ! n2o pressure factor
     $     en2o3,               ! emissivity of the 1168 cm-1 N2O band
     $     betac,               ! ch4 pressure factor
     $     ech4,                ! emissivity of 1306 cm-1 CH4 band
     $     betac1,              ! co2 pressure factor
     $     betac2,              ! co2 pressure factor
     $     eco21,               ! emissivity of 1064 cm-1 CO2 band
     $     eco22                ! emissivity of 961 cm-1 CO2 band
      real  tt(plon),           ! temp. factor for h2o overlap factor
     $     psi1,                ! narrow band h2o temp. factor
     $     phi1,                !             "
     $     p1,                  ! h2o line overlap factor
     $     w1,                  !          "
     $     tw(plon,6),         ! h2o transmission overlap
     $     g1(6),               ! h2o overlap factor
     $     g2(6),               !          "
     $     g3(6),               !          "
     $     g4(6),               !          "
     $     ab(6),               !          "
     $     bb(6),               !          "
     $     abp(6),              !          "
     $     bbp(6)               !          "
      real  tcfc3,               ! transmission for cfc11 band
     $     tcfc4,               !          "
     $     tcfc6,               ! transmission for cfc12 band
     $     tcfc7,               !          "
     $     tcfc8,               !          "
     $     tlw,                 ! h2o overlap factor
     $     tch4                 ! ch4 overlap factor
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,
     $         0.0540398,0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,
     &         0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,
     &         9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,
     $         0.0182932,0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
      real  func, u, b
      func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = sqrt(co2t(i,k))
c Transmission for h2o
         tt(i) = abs(co2t(i,k) - 250.0)
      end do
c
      do l = 1,6
         do i = 1,plon
            psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
            phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
            p1 = pnm(i,k) * (psi1/phi1) / sslp
            w1 = w(i,k) * phi1
            tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                  - g3(l)*s2c(i,k)-g4(l)*uptype(i,k))
         end do
      end do
c
      do i = 1,plon
c transmission due to cfc bands
            tcfc3 = exp(-175.005*ucfc11(i,k))
            tcfc4 = exp(-1202.18*ucfc11(i,k))
            tcfc6 = exp(-5786.73*ucfc12(i,k))
            tcfc7 = exp(-2873.51*ucfc12(i,k))
            tcfc8 = exp(-2085.59*ucfc12(i,k))
c  Emissivity for CFC11 bands
            ecfc1 = 50.0*(1.0 - exp(-54.09*ucfc11(i,k))) * tw(i,1) *
     $                                                  emplnk(7,i)
            ecfc2 = 60.0*(1.0 - exp(-5130.03*ucfc11(i,k)))* tw(i,2) *
     $                                                  emplnk(8,i)
            ecfc3 = 60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
            ecfc4 = 100.0*(1.0 - tcfc4)*tw(i,5)*emplnk(10,i)
c  Emissivity for CFC12 bands
            ecfc5 = 45.0*(1.0 - exp(-1272.35*ucfc12(i,k)))*tw(i,3)*
     $                                                     emplnk(11,i)
            ecfc6 = 50.0*(1.0 - tcfc6)*tw(i,4)*emplnk(12,i)
            ecfc7 = 80.0*(1.0 - tcfc7)*tw(i,5)* tcfc4 * emplnk(13,i)
            ecfc8 = 70.0*(1.0 - tcfc8)*tw(i,6) * emplnk(14,i)
c  Emissivity for CH4 band 1306 cm-1
            tlw = exp(-1.0*sqrt(up2(i)))
            betac = bch4(i,k)/uch4(i,k)
            ech4 = 6.00444*sqti(i)*log(1.0 + func(uch4(i,k),betac)) *
     $                 tlw * emplnk(3,i)
            tch4 = 1.0/(1.0 + 0.02*func(uch4(i,k),betac))
c  Emissivity for N2O bands
            u01 = un2o0(i,k)
            u11 = un2o1(i,k)
            beta01 = bn2o0(i,k)/un2o0(i,k)
            beta11 = bn2o1(i,k)/un2o1(i,k)
c     1285 cm-1 band
            en2o1 = 2.35558*sqti(i)*log(1.0 + func(u01,beta01) +
     $              func(u11,beta11))*tlw*tch4*emplnk(4,i)
            u02 = 0.100090*u01
            u12 = 0.0992746*u11
            beta02 = 0.964282*beta01
c     589 cm-1 band
            en2o2 = 2.65581*sqti(i)*log(1.0 + func(u02,beta02) +
     $              func(u12,beta02)) * tco2(i) * th2o(i) *
     $              emplnk(5,i)
            u03 = 0.0333767*u01
            beta03 = 0.982143*beta01
c     1168 cm-1 band
            en2o3 = 2.54034*sqti(i)*log(1.0 + func(u03,beta03)) *
     $                 tw(i,6) * tcfc8 * emplnk(6,i)
c  Emissivity for 1064 cm-1 band of CO2
            betac1 = 2.97558*pnm(i,k) / (sslp*sqti(i))
            betac2 = 2.0 * betac1
            eco21 = 3.7571*sqti(i)*log(1.0 + func(uco211(i,k),betac1)
     $         + func(uco212(i,k),betac2) + func(uco213(i,k),betac2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * emplnk(2,i)
c  Emissivity for 961 cm-1 band
            eco22 = 3.8443*sqti(i)*log(1.0 + func(uco221(i,k),betac1)
     $         + func(uco222(i,k),betac1) + func(uco223(i,k),betac2))
     $         * tw(i,4) * tcfc3 * tcfc6 * emplnk(1,i)
c total trace gas emissivity
            emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 +ecfc6
     $                + ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4
     $                + eco21 + eco22
      end do
      return
      end
c
c***************************** TRCPLK *********************************
c
      subroutine trcplk(tint, tlayr, tplnke, emplnk, abplnk1, abplnk2)
c----------------------------------------------------------------------
c   Calculate Planck factors for absorptivity and emissivity of
c   CH4, N2O, CFC11 and CFC12
c
c-----------------------------------------------------------------------
c
c $Id: trcplk.F,v 1.2 1995/02/17 21:29:00 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCON.h>
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real  tint(plon,plevp),  ! interface temperatures
     $     tlayr(plon,plevp), ! k-1 level temperatures
     $     tplnke(plon)       ! Top Layer temperature
c
c output arguments
c
      real  emplnk(14,plon),        ! emissivity Planck factor
     $     abplnk1(14,plon,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plon,plevp)  ! nearest layer factor
c
c local workspace
c
      integer wvl                   ! wavelength index
      integer i,k
      real  f1(14),                  ! Planck function factor
     $     f2(14),                  !        "
     $     f3(14)                   !        "
c
      data f1 /5.85713e8,7.94950e8,1.47009e9,1.40031e9,1.34853e8,
     $         1.05158e9,3.35370e8,3.99601e8,5.35994e8,8.42955e8,
     $         4.63682e8,5.18944e8,8.83202e8,1.03279e9/
      data f2 /2.02493e11,3.04286e11,6.90698e11,6.47333e11,
     $         2.85744e10,4.41862e11,9.62780e10,1.21618e11,
     $         1.79905e11,3.29029e11,1.48294e11,1.72315e11,
     $         3.50140e11,4.31364e11/
      data f3 /1383.,1531.,1879.,1849.,848.,1681.,
     $         1148.,1217.,1343.,1561.,1279.,1328.,
     $         1586.,1671./
c
c Calculate emissivity Planck factor
c
      do wvl = 1,14
         do i = 1,plon
            emplnk(wvl,i) = f1(wvl)/
     $                   (tplnke(i)**4.0*(exp(f3(wvl)/tplnke(i))-1.0))
         end do
      end do
c
c Calculate absorptivity Planck factor for tint and tlayr temperatures
c
      do wvl = 1,14
         do k = 1, plevp
            do i = 1, plon
c non-nearlest layer function
               abplnk1(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tint(i,k)))
     $              /(tint(i,k)**5.0*(exp(f3(wvl)/tint(i,k))-1.0)**2.0)
c nearest layer function
               abplnk2(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tlayr(i,k)))
     $            /(tlayr(i,k)**5.0*(exp(f3(wvl)/tlayr(i,k))-1.0)**2.0)
            end do
         end do
      end do
      return
      end
c
c***************************** TRCPTH *********************************
c
      subroutine trcpth(tnm, pnm, cfc11, cfc12, n2o, ch4, qnm,
     $                  ucfc11, ucfc12, un2o0,  un2o1,  uch4,
     $                  uco211, uco212, uco213, uco221, uco222,
     $                  uco223, bn2o0,  bn2o1,  bch4,   uptype)
c----------------------------------------------------------------------
c Calculate path lengths and pressure factors for CH4, N2O, CFC11
c and CFC12.
c           Coded by J.T. Kiehl, November 21, 1994.
c
c-----------------------------------------------------------------------
c
c $Id: trcpth.F,v 1.2 1995/02/17 21:29:04 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#include <RADSIZE.h>
C
C------------------------------Commons----------------------------------
#include <RADCON.h>
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real  tnm(plon,plev),     ! Model level temperatures
     $     pnm(plon,plevp),    ! Pressure at model interfaces (dynes/cm2)
     $     qnm(plon,plev),     ! h2o specific humidity
     $     cfc11(plon,plev),   ! CFC11 mass mixing ratio
     $     cfc12(plon,plev),   ! CFC12 mass mixing ratio
     $     n2o(plon,plev),     ! N2O mass mixing ratio
     $     ch4(plon,plev)      ! CH4 mass mixing ratio
C
C Output arguments
C
      real  ucfc11(plon,plevp), ! CFC11 path length
     $     ucfc12(plon,plevp), ! CFC12 path length
     $     un2o0(plon,plevp),  ! N2O path length
     $     un2o1(plon,plevp),  ! N2O path length (hot band)
     $     uch4(plon,plevp),   ! CH4 path length
     $     uco211(plon,plevp), ! CO2 9.4 micron band path length
     $     uco212(plon,plevp), ! CO2 9.4 micron band path length
     $     uco213(plon,plevp), ! CO2 9.4 micron band path length
     $     uco221(plon,plevp), ! CO2 10.4 micron band path length
     $     uco222(plon,plevp), ! CO2 10.4 micron band path length
     $     uco223(plon,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plon,plevp),  ! pressure factor for n2o
     $     bn2o1(plon,plevp),  ! pressure factor for n2o
     $     bch4(plon,plevp),   ! pressure factor for ch4
     $     uptype(plon,plevp)  ! p-type continuum path length
C
C---------------------------Local variables-----------------------------
C
      integer   i,              ! Longitude index
     $          k               ! Level index
      real  co2fac(plon,1),     ! co2 factor
     $     alpha1(plon),       ! stimulated emission term
     $     alpha2(plon),       ! stimulated emission term
     $     rt(plon),           ! reciprocal of local temperature
     $     rsqrt(plon),        ! reciprocal of sqrt of temp
     $     pbar(plon),         ! mean pressure
     $     dpnm(plon),         ! difference in pressure
     $     co2mmr
      real  diff                 ! diffusivity factor
      data diff /1.66/
c-----------------------------------------------------------------------
c  Calculate path lengths for the trace gases
c-----------------------------------------------------------------------
      co2mmr = 1.51913 * co2vmr
      do i = 1,plon
         ucfc11(i,1) = 1.8 * cfc11(i,1) * pnm(i,1) * rga
         ucfc12(i,1) = 1.8 * cfc12(i,1) * pnm(i,1) * rga
         un2o0(i,1) = diff * 1.02346e5 * n2o(i,1) * pnm(i,1) * rga
     $                       / sqrt(tnm(i,1))
         un2o1(i,1) = diff * 2.01909 * un2o0(i,1) *
     $                    exp(-847.36/tnm(i,1))
         uch4(i,1) = diff * 8.60957e4 * ch4(i,1) * pnm(i,1) * rga
     $                       / sqrt(tnm(i,1))
         co2fac(i,1) = diff * co2mmr * pnm(i,1) * rga
         alpha1(i) = (1.0 - exp(-1540.0/tnm(i,1)))**3.0/sqrt(tnm(i,1))
         alpha2(i) = (1.0 - exp(-1360.0/tnm(i,1)))**3.0/sqrt(tnm(i,1))
         uco211(i,1) = 3.42217e3 * co2fac(i,1) * alpha1(i) *
     $                             exp(-1849.7/tnm(i,1))
         uco212(i,1) = 6.02454e3 * co2fac(i,1) * alpha1(i) *
     $                             exp(-2782.1/tnm(i,1))
         uco213(i,1) = 5.53143e3 * co2fac(i,1) * alpha1(i) *
     $                             exp(-3723.2/tnm(i,1))
         uco221(i,1) = 3.88984e3 * co2fac(i,1) * alpha2(i) *
     $                             exp(-1997.6/tnm(i,1))
         uco222(i,1) = 3.67108e3 * co2fac(i,1) * alpha2(i) *
     $                             exp(-3843.8/tnm(i,1))
         uco223(i,1) = 6.50642e3 * co2fac(i,1) * alpha2(i) *
     $                             exp(-2989.7/tnm(i,1))
         bn2o0(i,1) = diff * 19.399 * pnm(i,1)**2.0 * n2o(i,1) *
     $                  1.02346e5 * rga / (sslp*tnm(i,1))
         bn2o1(i,1) = bn2o0(i,1) * exp(-847.36/tnm(i,1)) *
     $                  2.06646e5
         bch4(i,1) = diff * 2.94449 * ch4(i,1) * pnm(i,1)**2.0 * rga *
     $                  8.60957e4 / (sslp*tnm(i,1))
         uptype(i,1) = diff * qnm(i,1) * pnm(i,1)**2.0 *
     $                   exp(1800.0*(1.0/tnm(i,1) - 1.0/296.0)) *
     $                   rga / sslp
      end do
      do k = 1,plev
         do i = 1,plon
            rt(i) = 1./tnm(i,k)
            rsqrt(i) = sqrt(rt(i))
            pbar(i) = 0.5 * (pnm(i,k+1) + pnm(i,k)) / sslp
            dpnm(i) = (pnm(i,k+1) - pnm(i,k)) * rga
            alpha1(i) = diff * rsqrt(i) *
     $                         (1.0 - exp(-1540.0/tnm(i,k)))**3.0
            alpha2(i) = diff * rsqrt(i) *
     $                         (1.0 - exp(-1360.0/tnm(i,k)))**3.0
            ucfc11(i,k+1) = ucfc11(i,k) +  1.8 * cfc11(i,k) * dpnm(i)
            ucfc12(i,k+1) = ucfc12(i,k) +  1.8 * cfc12(i,k) * dpnm(i)
            un2o0(i,k+1) = un2o0(i,k) + diff * 1.02346e5 *
     $                                  n2o(i,k) * rsqrt(i) * dpnm(i)
            un2o1(i,k+1) = un2o1(i,k) + diff * 2.06646e5 * n2o(i,k) *
     $           rsqrt(i) * exp(-847.36/tnm(i,k)) * dpnm(i)
            uch4(i,k+1) = uch4(i,k) + diff * 8.60957e4 * ch4(i,k) *
     $           rsqrt(i) * dpnm(i)
            uco211(i,k+1) = uco211(i,k) + 1.15*3.42217e3 * alpha1(i) *
     $            co2mmr * exp(-1849.7/tnm(i,k)) * dpnm(i)
            uco212(i,k+1) = uco212(i,k) + 1.15*6.02454e3 * alpha1(i) *
     $            co2mmr * exp(-2782.1/tnm(i,k)) * dpnm(i)
            uco213(i,k+1) = uco213(i,k) + 1.15*5.53143e3 * alpha1(i) *
     $            co2mmr * exp(-3723.2/tnm(i,k)) * dpnm(i)
            uco221(i,k+1) = uco221(i,k) + 1.15*3.88984e3 * alpha2(i) *
     $            co2mmr * exp(-1997.6/tnm(i,k)) * dpnm(i)
            uco222(i,k+1) = uco222(i,k) + 1.15*3.67108e3 * alpha2(i) *
     $            co2mmr * exp(-3843.8/tnm(i,k)) * dpnm(i)
            uco223(i,k+1) = uco223(i,k) + 1.15*6.50642e3 * alpha2(i) *
     $            co2mmr * exp(-2989.7/tnm(i,k)) * dpnm(i)
            bn2o0(i,k+1) = bn2o0(i,k) + diff * 19.399 * pbar(i) * rt(i)
     $          * 1.02346e5 * n2o(i,k) * dpnm(i)
            bn2o1(i,k+1) = bn2o1(i,k) + diff * 19.399 * pbar(i) * rt(i)
     $          * 2.06646e5 * exp(-847.36/tnm(i,k)) * n2o(i,k)*dpnm(i)
            bch4(i,k+1) = bch4(i,k) + diff * 2.94449 * rt(i) * pbar(i)
     $            * 8.60957e4 * ch4(i,k) * dpnm(i)
            uptype(i,k+1) = uptype(i,k) + diff *qnm(i,k) *
     $                   exp(1800.0*(1.0/tnm(i,k) - 1.0/296.0)) *
     $                   pbar(i) * dpnm(i)
         end do
       end do
       return
       end
