c -- $Id: snow.F,v 1.2 2001/07/03 21:29:30 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/snow.F,v $
c -- $Name:  $

c------
c snow
c------
c
c----------------------------------------------------------------------
c
      SUBROUTINE INISNO (YTSIN, DEPTHSNO, DTSNOW, NSREST)
c
c        Does initialization for snow model
c
c        YTSIN   = ccm surface temperature field for initialization
c        DEPTHSNO= max limit on total snow thickness
c        DTSNOW  = time step for snow model
c        NSREST  = -1 for cold start, 0 for initial start from a
c                  restart file, 1 for true restart
c
#include <compar.h>
#include <comfud.h>
#include <comicexy.h>
#include <commugxy.h>
#include <comgrd.h>
#include <comsno.h>
#include <comsnoxy.h>
      DIMENSION YTSIN(NLONG, NLATG)
c
c
c        RHOS is density of snow
c
      RHOS = 250.
c
c        HSNOTOP is "adaptive-grid" thickness of top snow layer
c
      HSNOTOP = .05
c
c        HSNOMIN is the total snow thickness as snow fraction FI
c        approaches zero (global constant). As FI increases from 0 to 1,
c        snow thickness is assumed to increase linearly from HSNOMIN
c        to HSNORUF. HSNORUF is set in subr snow below and depends on
c        the lower veg cover. (HSNOMIN should be >= NSNOLAY*HSNOTOP for
c        VADAPT to work properly.)
c        Ref: Donald et al,1995,Water Res.Res.,31,4,995-1009.
c
      HSNOMIN = MAX (0.15, NSNOLAY*HSNOTOP)
c
c        HSNOMAX is an arbitrary upper limit to total snow thickness,
c        imposed near the top of subr snow (alternative to DEPTHSNO).
c
      HSNOMAX = 10.
c
c        FIMIN and FIMAX are minimum and maximum snowcover fractions.
c        All snow is removed when FI falls below FIMIN from above.
c        This is necessary to make snow vanish since we have no lateral
c        melt, just vertical melt which would have to melt thickness
c        ~HSNOMIN in one time step to make snow vanish. The value .001
c        should be small compared to typical snowfall rates * DTSNOW
c        / (RHOS*HSNOMIN) to allow creation of snow, but large enough
c        that typical melt rates will reduce FI to FIMIN in a
c        reasonable amount of time.
c
#if defined (SNOWBALL) 
      FIMIN = .000001* (DTSNOW/1800.) * (0.15/HSNOMIN)
#else
c     FIMIN = .00020 * (DTSNOW/1800.) * (0.15/HSNOMIN) ! v2.0.a
      FIMIN = .00002 * (DTSNOW/1800.) * (0.15/HSNOMIN)
#endif
 
#if defined (gcm) && ! defined (SNOWBALL)
      FIMAX = 0.95
#else
      FIMAX = 1.0
#endif
c
c        Z0SNO is roughness length of snow cover
c
#ifdef icedrive
      Z0SNO = 0.0010
#else
      Z0SNO = 0.0005
#endif
c
c        Interpolate initialization temperature field to lsx grid
c
      CALL GTOL (YTSIN, WORKL, 0)
c
c        Initialize snow variables. For Cray, set to no snow for
c        cold start (NSREST=-1), so the following for FIXY > 0 and
c        NSREST=-1 is really only for local machine with xfsno > 0.
c
      DO 100 JJ=1,NLAT
        DO 102 JI=1,NLON
 
c--------------------------------
c         Snow over land/icesheet
c--------------------------------
 
          IF (NSREST.EQ.-1) THEN
#ifdef gcm
            FIXY(JI,JJ) = 0.
#else
            FIXY(JI,JJ) = xfsno
#endif
          ENDIF
 
          IF (FIXY(JI,JJ).EQ.0. .OR. FWATERXY(JI,JJ).EQ.1.) THEN
 
            FIXY(JI,JJ) = 0.
            DO 110 K=1,NSNOLAY
              HSNOXY(JI,JJ,K) = 0.
              TSNOXY(JI,JJ,K) = TMELT
  110       CONTINUE
            TIXY(JI,JJ) = TMELT
c
          ELSE
c
            IF (NSREST.EQ.-1) THEN
#ifdef gcm
              H = MAX (0.3,  HSNOMIN)
#else
              H = MAX (xhsno,HSNOMIN)
#endif
              DO 112 K=1,NSNOLAY
                IF (K.EQ.1) HSNOXY(JI,JJ,K) = HSNOTOP
                IF (K.GT.1) HSNOXY(JI,JJ,K) = (H-HSNOTOP) / (NSNOLAY-1)
                TSNOXY(JI,JJ,K) = MIN (TMELT, WORKL(JI,JJ))
  112         CONTINUE
              TIXY(JI,JJ) = TSNOXY(JI,JJ,1)
            ENDIF
c
c              If namelist param depthsno entered, apply snow
c              thickness limit.
c
            IF (DEPTHSNO.GT.0.) THEN
               H = MAX (DEPTHSNO, HSNOMIN)
               DO 114 K=2,NSNOLAY
                 HSNOXY(JI,JJ,K) = MIN ( HSNOXY(JI,JJ,K),
     *                                   (H-HSNOTOP)/(NSNOLAY-1) )
  114          CONTINUE
            ENDIF
c
          ENDIF
 
c--------------------------------------
c         Snow over sea ice or lake ice
c--------------------------------------
 
          IF (NSREST.EQ.-1) THEN
#ifdef gcm
            FI_ICEXY(JI,JJ) = 0.
#else
            FI_ICEXY(JI,JJ) = xfsno
#endif
          ENDIF
 
          IF ( FI_ICEXY(JI,JJ).EQ.0. .OR. FWATERXY(JI,JJ).EQ.0. .OR.
     *         FICEXY(JI,JJ).EQ.0. ) THEN
 
            FI_ICEXY(JI,JJ) = 0.
            DO 120 K=1,NSNOLAY
              HSNO_ICEXY(JI,JJ,K) = 0.
              TSNO_ICEXY(JI,JJ,K) = TMELT
  120       CONTINUE
            TI_ICEXY(JI,JJ) = TMELT
c
          ELSE
c
            IF (NSREST.EQ.-1) THEN
#ifdef gcm
              H = MAX (0.3,  HSNOMIN)
#else
              H = MAX (xhsno,HSNOMIN)
#endif
              DO 122 K=1,NSNOLAY
                IF (K.EQ.1) HSNO_ICEXY(JI,JJ,K) = HSNOTOP
                IF (K.GT.1) HSNO_ICEXY(JI,JJ,K) =(H-HSNOTOP)/(NSNOLAY-1)
                TSNO_ICEXY(JI,JJ,K) = MIN (TMELT, WORKL(JI,JJ))
  122         CONTINUE
              TI_ICEXY(JI,JJ) = TSNO_ICEXY(JI,JJ,1)
            ENDIF
c
c              If namelist param depthsno entered, apply snow
c              thickness limit.
c
            IF (DEPTHSNO.GT.0.) THEN
               H = MAX (DEPTHSNO, HSNOMIN)
               DO 124 K=2,NSNOLAY
                 HSNO_ICEXY(JI,JJ,K) = MIN ( HSNO_ICEXY(JI,JJ,K),
     *                                       (H-HSNOTOP)/(NSNOLAY-1) )
  124          CONTINUE
            ENDIF
 
          ENDIF
c
  102   CONTINUE
  100 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SETSNO
c
c        Sets diagnostic snow quantities for this strip needed by snow
c        and lsx or ocnflx, returned in commus.
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
      DIMENSION ZDEPTH(NVEC)
#include <comsat.h>
 
 
c        Set thermal conductivity of snow (CONSNO). Depends on depth
c        via assumed increase from 100 Kg/m3 (new snow) towards
c        300 Kg/m3 (settled snow) (Paterson,1969,p.6). However use
c        constant snow density elsewhere for simplicity
c        (rhos, set in inisno).
 
c     CALL ZERO (ZDEPTH, NPOI)
 
      DO 100 K=1,NSNOLAY
        DO 102 I=1,NPOI
 
c         ZDEPTH(I) = ZDEPTH(I) + 0.5*HSNO(I,K)
c         10 cm e-folding increase of density with depth
c         ZRHOS = 100. + 200.*(1.-EXP(-ZDEPTH(I)/.10))
c         Loth,1993,JGR,98,10451,eq.5; cf.Williams+Smith,1989,p.110
c         CONSNO(I,K) = .02 + 2.5E-6*(ZRHOS**2)
c         ZDEPTH(I) = ZDEPTH(I) + 0.5*HSNO(I,K)
 
          CONSNO(I,K) = 0.2   ! v1.02
 
  102   CONTINUE
  100 CONTINUE
 
c        Set latent heat for snow surface, used by lsx and ocnflx.
c        Use prognostic top-layer temperature tsno, not the diagnostic
c        skin temperature ti, since budgets assume vapor is exchanged
c        with top-layer store at tsno.
 
      DO 400 I=1,NPOI
        HVASUI(I) = HSUBF(TSNO(I,1),TA(I))
  400 CONTINUE
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SNOWCTL (DTSNOW, TLOW, HLOW, CONLOW, ICODE)
c
c        Steps snow model through one timestep (over soil, icesheet
c        or seaice depending on ICODE=1,2, or 3)
c
#include <compar.h>
#include <combudl.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
c
c     DTSNOW = time step for this call
c     TLOW   = temperature of first soil/ice layer below snow
c     HLOW   = thickness of first soil/ice layer below snow
c     CONLOW = heat conductivity of first soil/ice layer below snow
c     ICODE  = 1 for snow on land, 2 on icesheet, 3  on seaice
c
      DIMENSION TLOW(NVEC), HLOW(NVEC), CONLOW(NVEC)
c
c     Local arrays:
c     HINIT  = initial layer thicknesses when snow first forms
c     INDSNO = index of points with snow in current 1D strip
c     FIOLD  = old FI at start of this timestep
c     FHTOP  = heat flux into upper snow surface
c     SFLO   = heat flux across snow and buried-lower-veg layer bdries
c     ZMELT  = liquid mass flux increments to soil/ice, at temperature
c              TMELT, due to processes occuring during this step
c     ZHEAT  = heat flux to soil/ice, due to processes occuring this
c              step
c     DFI    = change in FI
c     XL     = lower veg density
c     X*,HT  = temporary arrays
c
      DIMENSION
     *   HINIT(NVEC,NSNOLAY),           HSNORUF(NVEC),
     *   INDSNO(NVEC),
     *   FIOLD(NVEC),    FHTOP(NVEC),   SFLO(NVEC,NSNOLAY+2),
     *   ZMELT(NVEC),    ZHEAT(NVEC),   DFI(NVEC),    XL(NVEC),
     *   XH(NVEC),       XM(NVEC),      HT(NVEC),
     *   X1(NVEC),       X2(NVEC),      X3(NVEC)
c
c        Set local total snow thickness as FI increases to 1 (HSNORUF)
c        based on geometric height and coverage of lower veg,
c        and newly-created snow thickness profile (HINIT).
c        (HSNORUF must be > HSNOMIN; see comments in inisno, and do
c        loop 404 below). Also save pre-timestep snow fraction FIOLD,
c        and set local buried-veg quantities XL,X1 for later.
c
      IF (ICODE.EQ.1) THEN                                    ! veg/soil
        DO 10 I=1,NPOI
          HSNORUF(I) =  MIN (0.70, MAX (HSNOMIN+.05, FL(I)*ZTOP(I,1)))
          XL(I) = FL(I) * (LAI(I,1) + SAI(I,1))
          X1(I) = TLSUB(I)
   10   CONTINUE
      ELSE IF (ICODE.EQ.2 .OR. ICODE.EQ.3) THEN        ! icesheet/seaice
        DO 12 I=1,NPOI
          HSNORUF(I) = HSNOMIN
   12   CONTINUE
        CALL ZERO (XL, NPOI)
        CALL ZERO (X1, NPOI)
      ENDIF
c
      DO 20 I=1,NPOI
        HINIT(I,1) = HSNOTOP
        FIOLD(I) = FI(I)
   20 CONTINUE
      DO 22 K=2,NSNOLAY
        DO 24 I=1,NPOI
          HINIT(I,K) = (HSNOMIN - HSNOTOP) / (NSNOLAY-1)
   24   CONTINUE
   22 CONTINUE
c
c        Initialize local arrays and diagnostics in combudl
c        (in case not set below due to no snow in this strip)
c
      CALL ZERO (SFLO,  NVEC*(NSNOLAY+2))
      CALL ZERO (ZMELT, NVEC)
      CALL ZERO (ZHEAT, NVEC)
      CALL ZERO (WBURY, NVEC)
      CALL ZERO (SBURY, NVEC)
      CALL ZERO (TBURY, NVEC)
      CALL ZERO (HBURY, NVEC)
c
c        Set up index INDSNO, NPN for pts with snow. INDSNO is used
c        only by VADAPT - elsewhere below, just test on NPN > 0. Cray
c        efficiency without index outweighs saving by eliminating
c        non-snow points; also, index can't be used anyway for some
c        loops below.
c        Nb: This index is re-computed below after snowfall on bare grd
c
      CALL WHENFGT (NPOI, FI, 1, 0., INDSNO, NPN)
c
c-----------------------
      IF (NPN.GT.0) THEN
c-----------------------
c
c          Impose arbitrary maximum snow depth of 10 m, mainly for
c          aesthetic purposes for accumulation areas of ice sheets,
c          as an alternative to using namelist DEPTHSNO.
c          *** Nb: This does not conserve h2o mass ***
 
        DO K=2,NSNOLAY
          DO I=1,NPOI
            HSNO(I,K) = MIN (HSNO(I,K), (HSNOMAX-HSNOTOP)/(NSNOLAY-1))
          ENDDO
        ENDDO
c
c          Set surface heat flux FHTOP and increment top layer thickness
c          due to snowfall and sublimation on existing snow. (Account
c          for heat to bring rainfall to melt point...after that,
c          accounted for by percolation/refreezing below)
c
        DO 40 I=1,NPOI
          FHTOP(I) = HEATI(I)
     *             + RAINL(I)*CH2O*(TRAINL(I)-TMELT)
     *             + SNOWL(I)*CICE*(TSNOWL(I)-TSNO(I,1))
          IF (FI(I).GT.0.) HSNO(I,1)= HSNO(I,1)
     *                   + DTSNOW*(SNOWL(I)-FVAPI(I))/RHOS
   40   CONTINUE
c
c          Step temperatures due to heat conduction, including buried
c          lower-veg temperature TLSUB. Record heat lost by buried
c          lower veg in HBURY, for lsx budget.
c
        CALL SNOWHEAT (TLSUB, FHTOP, SFLO, TLOW, HLOW, CONLOW,
     *                 XL, CHL, CONSNO, DTSNOW)
c
        IF (ICODE.EQ.1) THEN
          DO 44 I=1,NPOI
            HBURY(I)= HBURY(I) - FI(I)*XL(I)*CHL*(TLSUB(I)-X1(I))/DTSNOW
   44     CONTINUE
        ENDIF
 
c----------
      ENDIF
c----------
c
c        Put snowfall from 1-FI snow-free area onto side of existing
c        snow, or create new snow if current FI = 0. Then reset index.
c        (Total depth of newly created snow = HSNOMIN.)
c        (FI will not become gt 1 here if one timestep's snowfall
c        <= HSNOMIN, but protect against this anyway.)
c
c        If no adjacent snowfall or FI = 1, DFI = 0, so no effect.

c     adjust snowg, fvapg within 1-fi fraction, to avoid giving all of
c     snowg to side of existing fi snow if snowg<fvapg. Added 9/23/11:
      do i=1,npoi
        if (fi(i).lt.1.) then
          zsme = snowg(i) - fvapg(i)
          if (zsme.ge.0.) then
c           give as snow to fi snow fraction (side):
            snowg(i) = zsme
            fvapg(i) = 0.
          else
c           remove as evap from 1-fi non-snow fraction:
            snowg(i) = 0.
            fvapg(i) = -zsme
          endif
        endif
      enddo
c
      CALL ZERO (HT, NPOI)
      DO 190 K=1,NSNOLAY
        DO 192 I=1,NPOI
          HT(I) = HT(I) + HSNO(I,K)
  192   CONTINUE
  190 CONTINUE
c
      DO 195 I=1,NPOI
        HT(I) = CVMGT (HSNOMIN, HT(I), HT(I).EQ.0.)
  195 CONTINUE
c
      DO 200 I=1,NPOI
        DFI(I) = (1.-FI(I))*DTSNOW*(SNOWG(I)/RHOS) / HT(I)
        DFI(I) = MIN (DFI(I), 1.-FI(I))
  200 CONTINUE
c
      DO 210 K=1,NSNOLAY
        DO 212 I=1,NPOI
          IF (FI(I)+DFI(I).GT.0.)
     *      TSNO(I,K) = (TSNO(I,K)*FI(I) + TSNOWG(I)*DFI(I))
     *                / (FI(I)+DFI(I))
c         Set initial thicknesses for newly created snow
          IF (FI(I).EQ.0. .AND. DFI(I).GT.0.) HSNO(I,K) = HINIT(I,K)
  212   CONTINUE
  210 CONTINUE
c
      DO 215 I=1,NPOI
c       FI(I) = FI(I) + DFI(I)           ! if snowg < 0, fi+dfi can < 0,
        FI(I) = MAX (FI(I) + DFI(I), 0.) ! so avoid regcm3 Nans. 8/06.
 215  CONTINUE
c
c        Recompute snow/no-snow index
c
      CALL WHENFGT (NPOI, FI, 1, 0., INDSNO, NPN)
c
c-----------------------
      IF (NPN.GT.0) THEN
c-----------------------
c
c          Melt from any layer (due to implicit heat conduction, any
c          layer can exceed TMELT, not just the top layer), and reduce
c          thicknesses (even to zero, and give extra heat to soil/ice).
c          Ok to do it for non-snow points, for which XH = XM = 0.

c          4/02: percolation/refreezing added, kept track of by xm()
c          through loop 300 from top to bottom. Initialize xm to
c          rainfall on top of snow (already at melt point, see above)
c
        do i=1,npoi
          xm(i) = cvmgt (rainl(i)*dtsnow, 0., fi(i).gt.0.)
        enddo

        DO 300 K=1,NSNOLAY
          DO 302 I=1,NPOI
c           percolation/refreezing of liquid from above. First, 
c           convert all liquid to snow, add to existing layer at a new
c           temperature that may temporarily be > tmelt:
            zhi = rhos*hsno(i,k)*cice * (tsno(i,k)-tmelt)
            zhw = hfus*xm(i)
            hsno(i,k) = hsno(i,k) + xm(i)/rhos
            tsno(i,k) = tmelt + (zhi+zhw)
     *                          / (rhos*max(hsno(i,k),epsilon)*cice)

c           Second, if > tmelt, convert some ice to liquid (xm, used 
c           for next layer down), at expense of layer thickness:
            XH(I) = RHOS*HSNO(I,K)*CICE * MAX(TSNO(I,K)-TMELT, 0.)
            XM(I) = MIN (RHOS*HSNO(I,K), XH(I)/HFUS)
            HSNO(I,K) = HSNO(I,K) - XM(I)/RHOS
            TSNO(I,K) = MIN (TSNO(I,K),TMELT)
c           extra heat if entire layer is melted: 
            ZHEAT(I) = ZHEAT(I) + FI(I)*(XH(I)-HFUS*XM(I))
  302     CONTINUE
  300   CONTINUE

c       Assign any liquid reaching bottom of snowpack to zmelt:
        do i=1,npoi
          zmelt(i) = zmelt(i) + fi(i)*xm(i)
        enddo
c
c          Adjust FI and thicknesses for coverage-vs-volume relation
c          ie, total thickness = HSNOMIN*(1-FI) + HSNORUF*FI for
c          0<FI<FIMAX. (Ok to do it for no-snow points, for which
c          HT=FI=FINEW=0.)
c
        CALL ZERO (HT, NPOI)
        DO 400 K=1,NSNOLAY
          DO 402 I=1,NPOI
            HT(I) = HT(I) + HSNO(I,K)
  402     CONTINUE
  400   CONTINUE
c
        ZHH = 0.5*HSNOMIN
        DO 404 I=1,NPOI
          if (icode.eq.2 .or. icode.eq.3) then ! changed 9/23/11
            FINEW = FI(I)*HT(I)/HSNORUF(I)     ! constant depth, 0<fi<1
          else
            ZDH = HSNORUF(I)-HSNOMIN           ! linear variation,0<fi<1
            FINEW = ( -ZHH + SQRT(ZHH**2 + ZDH*FI(I)*HT(I)) ) / ZDH
          endif
          FINEW = MAX (0., MIN (FIMAX, FINEW))
          X1(I) =  FI(I) / MAX (FINEW, EPSILON)
          FI(I) =  FINEW
  404   CONTINUE
c
        DO 406 K=1,NSNOLAY
          DO 408 I=1,NPOI
            HSNO(I,K) = HSNO(I,K) * X1(I)
  408     CONTINUE
  406   CONTINUE
c
c----------
      ENDIF
c----------
c
c        Re-adapt snow thickness profile, so top thickness = HSNOTOP
c        and other thicknesses are equal (VADAPT is in ocean.f). Also
c        adjust temperature to conserve sensible heat.
c
      CALL VADAPT (HSNO, TSNO, HSNOTOP, INDSNO, NPN, NSNOLAY)
c
c-----------------------
      IF (NPN.GT.0) THEN
c-----------------------
c
c        If FI is below FIMIN, melt all snow and adjust soil/ice fluxes
c
        CALL SCOPY (NPOI, FI, 1, X1, 1)
        DO 500 K=1,NSNOLAY
          DO 502 I=1,NPOI
            IF (X1(I).LT.FIMIN) THEN
              XM(I) = X1(I) * RHOS * HSNO(I,K)
              ZMELT(I) = ZMELT(I) + XM(I)
              ZHEAT(I) = ZHEAT(I) - XM(I)*(CICE*(TMELT-TSNO(I,K))+HFUS)
              HSNO(I,K) = 0.
              TSNO(I,K) = TMELT
              FI(I) = 0.
            ENDIF
  502     CONTINUE
  500   CONTINUE
c
c          Adjust buried lower veg for FI changes. If FI has increased,
c          incorporate newly buried intercepted h2o into bottom-layer
c          snow, giving associated heat increment to soil/ice, and mix
c          the specific heat of newly buried veg (at TL) into TLSUB.
c          If FI has decreased, change temp of newly exhumed veg to TL,
c          giving assoc heat increment to soil/ice, and smear out
c          intercepted h2o.
c          [W,S,T,H]BURY are for budget accumulation in BUDAD2, called
c          from SURFCTL, since buried lower veg and its intercepted h2o
c          are included in the lsx budget.
c
c          Don't use index for Cray vec. Ok for non-snow points for
c          which DFI = 0.
c
        DO 600 I=1,NPOI
          DFI(I) = FI(I) - FIOLD(I)
c
          IF (DFI(I).GT.0.) THEN
c           Factor of XL*CHL has been divided out of next line
            TLSUB(I)= (TLSUB(I)*FIOLD(I) + TL(I)*DFI(I)) / FI(I)
            ZHEAT(I) = ZHEAT(I) + DFI(I)*XL(I)
     *                 * ( WLIQL(I) * (CH2O*(TL(I)-TMELT) + HFUS
     *                                +CICE*(TMELT-TSNO(I,NSNOLAY)))
     *                   + WSNOL(I) *  CICE*(TL(I)-TSNO(I,NSNOLAY)) )
c
            HSNO(I,NSNOLAY) = HSNO(I,NSNOLAY)
     *                      + DFI(I)*XL(I)*(WLIQL(I)+WSNOL(I))
     *                        / (RHOS*FI(I))
            WBURY(I) = WBURY(I) + DFI(I)*XL(I)*WLIQL(I) / DTSNOW
            SBURY(I) = SBURY(I) + DFI(I)*XL(I)*WSNOL(I) / DTSNOW
            TBURY(I) = TL(I)
          ENDIF
c
          IF (DFI(I).LT.0.) THEN
            ZHEAT(I) = ZHEAT(I) - DFI(I)*XL(I)*CHL*(TLSUB(I)-TL(I))
            WLIQL(I) = WLIQL(I) * (1.-FIOLD(I)) / (1.-FI(I))
            WSNOL(I) = WSNOL(I) * (1.-FIOLD(I)) / (1.-FI(I))
            HBURY(I) = HBURY(I)
     *               - DFI(I)*XL(I)*CHL*(TLSUB(I)-TL(I)) / DTSNOW
          ENDIF
c
  600   CONTINUE
c
c          Areally average fluxes to be used by soil/ice model.
c          (Don't use index due to MIX call, but only need at all
c          if NPN > 0).
c
        DO 700 I=1,NPOI
          HEATG(I) = (1.-FIOLD(I))*HEATG(I)
     *             + FIOLD(I)*SFLO(I,NSNOLAY+2)
     *             + ZHEAT(I)/DTSNOW
          SOLG(I) = (1.-FIOLD(I)) * SOLG(I)
          FVAPG(I) = (1.-FIOLD(I)) * FVAPG(I)
          X1(I) = (1.-FIOLD(I)) * RAING(I)
          X2(I) = ZMELT(I)/DTSNOW
          X3(I) = TMELT
  700   CONTINUE
        CALL MIX (RAING,TRAING, X1,TRAING, X2,X3, VZERO,VZERO)
c
c----------
      ENDIF
c----------
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SNOWHEAT (TLSUB, FHTOP, SFLO, TLOW, HLOW, CONLOW,
     *                     XL, CHL, CONSNO, DTIM)
c
c        Sets up call to TRIDIA to solve implicit snow heat conduction,
c        using snow temperatures in TSNO (in comsno). Adds an extra
c        buried-lower-veg layer to the bottom of the snow with
c        conduction coefficient CONBUR/XL and heat capacity CHL*XL.
c        CONBUR (for XL=1) is chosen to be equiv to 10 cm of snow.
c
c        TLSUB = temperature of buried lower veg (supplied, returned)
c        FHTOP = heat flux into top snow layer from atmos (supplied)
c        SFLO  = downward heat flow across layer boundaries (returned)
c        TLOW  = temperature of first soil/ice layer (supplied)
c        HLOW  = thickness of first soil/ice layer (supplied)
c        CONLOW= heat conductivity of first soil/ice layer (supplied)
c        XL    = (LAI(I,1)+SAI(I,1))*FL(I), lower-veg density(supplied)
c        CHL   = specific heat of lower veg per l/s area (supplied)
c        CONSNO= heat conductivity in each snow layer (supplied)
c        DTIM  = time step (supplied)
c
c        Local arrays and scalars:
c        CON    = conduction coefficents between layers
c        TEMP   = combined snow and buried-veg temperatures
c        D1,2,3 = diagonals of tridiagonal systems of equations
c        RHS    = right-hand sides of systems of equations
c        W1,2   = work arrays needed by TRIDIA
c        RIMP   = implicit fraction of the calculation (0 to 1)
c        CONBUR = conduction coeff of buried lower veg layer
c                 for unit density XL=(LAI+SAI)*FL, in W m-2 K-1
c        HFAKE  = arbitrary small thickness to allow processing
c                 for zero snow. (Doesn't use index since TRIDIA
c                 not set up for index.)
c
#include <compar.h>
#include <comgrd.h>
#include <comsno.h>
      DIMENSION
     * TLSUB(NVEC),   FHTOP(NVEC),    SFLO(NVEC,NSNOLAY+2),
     * TLOW(NVEC),    HLOW(NVEC),     CONLOW(NVEC),           XL(NVEC),
     * CONSNO(NVEC,NSNOLAY)
c
      DIMENSION
     *  CON(NVEC,NSNOLAY+2),  TEMP(NVEC,NSNOLAY+1),
     *  D1(NVEC,NSNOLAY+1),   D2(NVEC,NSNOLAY+1),   D3(NVEC,NSNOLAY+1),
     *  RHS(NVEC,NSNOLAY+1),  W1(NVEC,NSNOLAY+1),   W2(NVEC,NSNOLAY+1)
c
      SAVE RIMP, CONBUR, HFAKE
      DATA RIMP, CONBUR, HFAKE /1.0, 2.0, .01/
c
c
c        Copy snow and buried-lower-veg temperatures into combined
c        array TEMP
c
      CALL SCOPY (NVEC*NSNOLAY, TSNO,  1, TEMP,              1)
      CALL SCOPY (NVEC,         TLSUB, 1, TEMP(1,NSNOLAY+1), 1)
c
c       Set conduction coefficients between layers (CONSNO has been set
c       in SETSNO)
c
      DO 100 K=1,NSNOLAY+2
        IF (K.EQ.1) THEN
          CALL ZERO (CON(1,K), NPOI)
c
        ELSE IF (K.LE.NSNOLAY) THEN
          DO 102 I=1,NPOI
            CON(I,K) = 1. / (   0.5*MAX(HSNO(I,K-1),HFAKE)/CONSNO(I,K-1)
     *                        + 0.5*MAX(HSNO(I,K)  ,HFAKE)/CONSNO(I,K) )
  102     CONTINUE
c
        ELSE IF (K.EQ.NSNOLAY+1) THEN
          DO 104 I=1,NPOI
            CON(I,K) = 1. / (   0.5*MAX(HSNO(I,K-1),HFAKE)/CONSNO(I,K-1)
     *                        + 0.5*XL(I)/CONBUR )
  104     CONTINUE
c
        ELSE IF (K.EQ.NSNOLAY+2) THEN
          DO 106 I=1,NPOI
            CON(I,K) = 1. / (   0.5*XL(I)/CONBUR
     *                        + 0.5*MAX(HLOW(I),HFAKE)  /CONLOW(I) )
  106     CONTINUE
        ENDIF
  100 CONTINUE
c
c        Set matrix diagonals and right-hand side. For layer NSNOLAY+1
c        (buried-lower-veg layer), use explicit contact with soil/ice,
c        and multiply eqn through by XL*CHL/DTIM to allow zero XL.
c
      DO 200 K=1,NSNOLAY+1
        KM1 = MAX (K-1,1)
        KP1 = MIN (K+1,NSNOLAY+1)
c
        IF (K.LE.NSNOLAY) THEN
          DO 202 I=1,NPOI
            DT = DTIM / (RHOS*CICE*MAX(HSNO(I,K),HFAKE))
            D1(I,K) =    - DT*RIMP* CON(I,K)
            D2(I,K) = 1. + DT*RIMP*(CON(I,K)+CON(I,K+1))
            D3(I,K) =    - DT*RIMP* CON(I,K+1)
c
            RHS(I,K) = TEMP(I,K) + DT
     *               * ( (1.-RIMP)*CON(I,K)  *(TEMP(I,KM1)-TEMP(I,K))
     *                 + (1.-RIMP)*CON(I,K+1)*(TEMP(I,KP1)-TEMP(I,K)) )
  202     CONTINUE
c
          IF (K.EQ.1) THEN
            DO 204 I=1,NPOI
              DT = DTIM / (RHOS*CICE*MAX(HSNO(I,K),HFAKE))
              RHS(I,K) = RHS(I,K) + DT*FHTOP(I)
  204       CONTINUE
          ENDIF
c
        ELSE IF (K.EQ.NSNOLAY+1) THEN
c
          DO 206 I=1,NPOI
            DTI = (XL(I)*CHL) / DTIM
            D1(I,K) =     -  RIMP* CON(I,K)
            D2(I,K) = DTI +  RIMP*(CON(I,K)+CON(I,K+1))
            D3(I,K) = 0.
            RHS(I,K) = DTI*TEMP(I,K)
     *               + ( (1.-RIMP)*CON(I,K)*(TEMP(I,KM1)-TEMP(I,K))
     *                 + CON(I,K+1)*(TLOW(I)-(1.-RIMP)*TEMP(I,K)) )
  206     CONTINUE
        ENDIF
  200 CONTINUE
c
c        Solve the tridiagonal systems
c
      CALL TRIDIA (NPOI, NVEC, NSNOLAY+1, D1,D2,D3, RHS, TEMP, W1,W2)
c
c        Deduce downward heat fluxes between layers
c
      CALL SCOPY (NPOI, FHTOP, 1, SFLO(1,1), 1)
c
      DO 400 K=1,NSNOLAY+1
        IF (K.LE.NSNOLAY) THEN
          DO 402 I=1,NPOI
            SFLO(I,K+1) = SFLO(I,K) - RHOS*CICE*HSNO(I,K)
     *                                *(TEMP(I,K)-TSNO(I,K)) / DTIM
  402     CONTINUE
c
        ELSE
          DO 404 I=1,NPOI
            SFLO(I,K+1) = SFLO(I,K)
     *                  - XL(I)*CHL*(TEMP(I,NSNOLAY+1)-TLSUB(I)) / DTIM
  404     CONTINUE
        ENDIF
  400 CONTINUE
c
c        Copy temperature solution to TSNO and TLSUB, but not for
c        points with no snow
c
      DO 500 K=1,NSNOLAY
        DO 502 I=1,NPOI
          TSNO(I,K) = CVMGT (TEMP(I,K), TSNO(I,K), FI(I).GT.0.)
  502   CONTINUE
  500 CONTINUE
c
      DO 510 I=1,NPOI
        TLSUB(I) = CVMGT (TEMP(I,NSNOLAY+1), TLSUB(I), FI(I).GT.0.)
  510 CONTINUE
c
      RETURN
      END
