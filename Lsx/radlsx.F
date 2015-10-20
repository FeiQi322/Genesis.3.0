c -- $Id: radlsx.F,v 1.2 2001/07/03 21:29:29 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/radlsx.F,v $
c -- $Name:  $

c--------
c radlsx
c--------
c
c----------------------------------------------------------------------
c
      SUBROUTINE SOLSET
c
c       Zeros albedos and internal lsx absorbed solar fluxes, and sets
c       index for other solar routines. The index INDSOL, with number
c       of points NSOL, points to current 1D strip arrays whose COSZEN
c       values are gt 0. (INDSOL, NSOL are in com1d.)
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
c
c
c        Zero albedos returned to gcm just as a niceity (will remain
c        zero only for points with COSZEN lt 0).
 
      CALL ZERO (ASURD, NVEC*NBAND)
      CALL ZERO (ASURI, NVEC*NBAND)
c
c        Zeros absorbed solar fluxes SOL[U,S,L,G,I]1 since only points
c        with +ve COSZEN will be set in SOLARF, and since
c        SOL[U,L,S,G,I] are summed over wavebands in SOLARF.
c
      CALL ZERO (SOLU, NVEC)
      CALL ZERO (SOLL, NVEC)
      CALL ZERO (SOLS, NVEC)
      CALL ZERO (SOLG, NVEC)
      CALL ZERO (SOLI, NVEC)
 
c        Similarly zero PAR-related arrays set in SOLARF for TURVAP
 
      CALL ZERO (SUNPU, NVEC)
      CALL ZERO (SUNPL, NVEC)
      CALL ZERO (SHAPU, NVEC)
      CALL ZERO (SHAPL, NVEC)
      CALL ZERO (SUNFU, NVEC)
      CALL ZERO (SUNFL, NVEC)
c
c        Set index of points with positive COSZEN
c
      CALL WHENFGT (NPOI, COSZEN, 1, 0., INDSOL, NSOL)
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SOLSUR (ITYPE)
c
c        Sets surface albedos for soil, icesheet or ocean/seaice, and
c        for snow, prior to other lsx solar calculations.
c        For ocean strips, also sets overall surface albedos
c        ASUR[D,I] (done later in SOLALB for land strips)
c
c        ITYPE = 1 for land, 2 for icesheet, 3 for ocean strips.
c
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <commug.h>
#include <comocn.h>
#include <comsno.h>
#include <comsoi.h>
#include <com1d.h>
 
      DIMENSION X(NVEC), ZFAC(NVEC)
 
#if defined (TWOSTRICE)         
      dimension zfluxd(nicelay+1),  zfluxi(nicelay+1)
#endif

c     Albedos at [LOw,HIgh] threshold temperatures:
c     A7SV[LO,HI] .. Snow, Visible
c     A7SN[LO,HI] .. Snow, Near-ir
c     A7TV[LO,HI] .. Snow over icesheet, Visible
c     A7TN[LO,HI] .. Snow over icesheet, Near-ir
c     A7LV[LO,HI] .. Ice sheet, Visible
c     A7LN[LO,HI] .. Ice sheet, Near-ir
c     A7IV[LO,HI] .. Sea Ice, Visible
c     A7IN[LO,HI] .. Sea Ice, Near-ir
 
      SAVE A7SVLO,A7SVHI,A7SNLO,A7SNHI,
     *     A7TVLO,A7TVHI,A7TNLO,A7TNHI,
     *     A7LVLO,A7LVHI,A7LNLO,A7LNHI,
     *     A7IVLO,A7IVHI,A7INLO,A7INHI
 
c     [LOw,HIgh] threshold temperatures for albedos:
c     T7S[LO,HI] .. Snow
c     T7T[LO,HI] .. Snow over ice sheet
c     T7L[LO,HI] .. Ice sheet
c     T7I[LO,HI] .. Sea ice
 
#ifndef omp
      SAVE T7SLO, T7SHI,
     *     T7TLO, T7THI,
     *     T7LLO, T7LHI,
     *     T7ILO, T7IHI
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
#endif
 
c--------------
#ifdef icedrive
c--------------
      DATA A7SVLO, A7SVHI /0.90, 0.55/        ! snow (non-icesheet)
      DATA A7SNLO, A7SNHI /0.60, 0.35/
c     DATA A7SVLO, A7SVHI /0.01, 0.01/        ! snow (non-icesheet)
c     DATA A7SNLO, A7SNHI /0.01, 0.01/
 
#ifdef testb
      DATA A7TVLO, A7TVHI /0.90, 0.60/        ! snow over icesheet
      DATA A7TNLO, A7TNHI /0.70, 0.45/
      DATA A7LVLO, A7LVHI /0.80, 0.50/        ! ice sheet
      DATA A7LNLO, A7LNHI /0.60, 0.35/
#else
      DATA A7TVLO, A7TVHI /0.90, 0.55/        ! snow over icesheet
      DATA A7TNLO, A7TNHI /0.60, 0.35/
      DATA A7LVLO, A7LVHI /0.90, 0.70/        ! ice sheet
      DATA A7LNLO, A7LNHI /0.70, 0.55/
#endif
 
c----
#else
c----
      DATA A7SVLO, A7SVHI /0.90, 0.55/        ! snow (non-icesheet)
      DATA A7SNLO, A7SNHI /0.60, 0.35/
 
c     Higher cold_snow_visible albedo over Antarctic icesheet
c     (Grenfell,Warren,Mullen,1994,JGR,99,D9,18669-18684):
      DATA A7TVLO, A7TVHI /0.95, 0.70/        ! snow over icesheet
      DATA A7TNLO, A7TNHI /0.70, 0.55/
 
      DATA A7LVLO, A7LVHI /0.90, 0.70/        ! ice sheet
      DATA A7LNLO, A7LNHI /0.70, 0.55/
 
c     DATA A7LVLO, A7LVHI /0.90, 0.85/        ! cona
c     DATA A7LNLO, A7LNHI /0.70, 0.65/        ! cona
c-----
#endif
c-----
 
c     DATA A7IVLO, A7IVHI /0.80, 0.70/        ! sea ice ! v2.3
c     DATA A7INLO, A7INHI /0.50, 0.40/                  ! v2.3
      DATA A7IVLO, A7IVHI /0.85, 0.80/                  ! v3 
      DATA A7INLO, A7INHI /0.55, 0.50/                  ! v3
 
c       Compute T7* every time and don't save if multiprocessing
c       over 1D strips
 
#ifndef omp
      IF (FIRST) THEN
#endif
cv1.02  T7SLO = TMELT - 15.
        T7SLO = TMELT - 5.
        T7SHI = TMELT
        T7TLO = TMELT - 5.
        T7THI = TMELT
        T7LLO = TMELT - 5.
        T7LHI = TMELT
        T7ILO = TMELT - 5.
        T7IHI = TMELT
#ifndef omp
        FIRST = .FALSE.
      ENDIF
#endif
 
c        Do nothing if all points in current strip have COSZEN le 0
 
      IF (NSOL.EQ.0) RETURN
 
c====================
c        Snow albedos
c====================
 
c-------------------------
      IF (ITYPE.NE.2) THEN
c-------------------------
 
c          Land/Ocean strip (snow)
 
!DIR$ IVDEP
        DO 50 J=1,NSOL
          I = INDSOL(J)
 
c         Snow, visible:
          X(I) = (A7SVHI*(TSNO(I,1)-T7SLO) + A7SVLO*(T7SHI-TSNO(I,1)))
     *           / (T7SHI-T7SLO)
          X(I) = MIN (A7SVLO, MAX (A7SVHI, X(I)))
          ZFAC(I) = MAX ( 0., 1.5/(1.+4.*COSZEN(I)) - 0.5 )
          ALBSND(I,1) = MIN (0.99, X(I) + (1.-X(I))*ZFAC(I))
          ALBSNI(I,1) = X(I)
 
c         Snow, near-IR:
          X(I) = (A7SNHI*(TSNO(I,1)-T7SLO) + A7SNLO*(T7SHI-TSNO(I,1)))
     *           / (T7SHI-T7SLO)
          X(I) = MIN (A7SNLO, MAX (A7SNHI, X(I)))
          ALBSND(I,2) = MIN (0.99, X(I) + (1.-X(I))*ZFAC(I))
          ALBSNI(I,2) = X(I)
 
   50   CONTINUE
 
c---------
      ELSE
c---------
 
c          Icesheet strip (snow)
 
!DIR$ IVDEP
        DO 60 J=1,NSOL
          I = INDSOL(J)
 
c         Snow, visible:
          X(I) = (A7TVHI*(TSNO(I,1)-T7TLO) + A7TVLO*(T7THI-TSNO(I,1)))
     *           / (T7THI-T7TLO)
          X(I) = MIN (A7TVLO, MAX (A7TVHI, X(I)))
          ZFAC(I) = MAX ( 0., 1.5/(1.+4.*COSZEN(I)) - 0.5 )
          ALBSND(I,1) = MIN (0.99, X(I) + (1.-X(I))*ZFAC(I))
          ALBSNI(I,1) = X(I)
 
c         Snow, near-IR:
          X(I) = (A7TNHI*(TSNO(I,1)-T7TLO) + A7TNLO*(T7THI-TSNO(I,1)))
     *           / (T7THI-T7TLO)
          X(I) = MIN (A7TNLO, MAX (A7TNHI, X(I)))
          ALBSND(I,2) = MIN (0.99, X(I) + (1.-X(I))*ZFAC(I))
          ALBSNI(I,2) = X(I)
 
   60   CONTINUE
 
c----------
      ENDIF
c----------
 
c=============================================
c        Soil/ice/seaice/ocean surface albedos
c=============================================
 
c-------------------------
      IF (ITYPE.EQ.1) THEN
c-------------------------
 
c          Land strip
 
!DIR$ IVDEP
        DO 100 J=1,NSOL
          I = INDSOL(J)
 
c         ZW = WSOI(I,1)*(1.-WISOI(I,1)) + WISOI(I,1)
          ZW = WSOI(I,1)*(1.-WISOI(I,1))
 
c         Soil, visible:
c         ALBSOD(I,1) = MIN (1., 0.15 * 1.1 / (1.+.2*COSZEN(I)))
c         ALBSOI(I,1) = 0.15
c         BATS formulation (BATS Tech Note, Eqn 2)
c         DINC = MAX (.11-.40*ZW, 0.)
c         ALBSOD(I,1) = MIN (ALBSAV(I)+DINC, 2.*ALBSAV(I))
 
          DINC = 1.0 + 1.0 * MIN (1., MAX (0., 1.-(ZW/.3)    ))
          ALBSOD(I,1) = MIN (ALBSAV(I)*DINC, .80)
          ALBSOI(I,1) = ALBSOD(I,1)
 
c         Soil, near-IR:
c         ALBSOD(I,2) = MIN (1., 0.3 * 1.1 / (1.+.2*COSZEN(I)))
c         ALBSOI(I,2) = 0.3
c         BATS formulation (BATS Tech Note, Eqn 2)
c         ALBSOD(I,2) = MIN (ALBSAN(I)+2.*DINC, 2.*ALBSAN(I))
 
          ALBSOD(I,2) = MIN (ALBSAN(I)*DINC, .80)
          ALBSOI(I,2) = ALBSOD(I,2)
  100   CONTINUE
 
c------------------------------
      ELSE IF (ITYPE.EQ.2) THEN
c------------------------------
 
c          Icesheet strip
 
!DIR$ IVDEP
        DO 200 J=1,NSOL
          I = INDSOL(J)
 
c         Icesheet, visible:
          X(I) = (A7LVHI*(TSOI(I,1)-T7LLO) + A7LVLO*(T7LHI-TSOI(I,1)))
     *           / (T7LHI-T7LLO)
          X(I) = MIN (A7LVLO, MAX (A7LVHI, X(I)))
          ALBSOD(I,1) = X(I)
          ALBSOI(I,1) = X(I)
 
c         Icesheet, near-IR:
          X(I) = (A7LNHI*(TSOI(I,1)-T7LLO) + A7LNLO*(T7LHI-TSOI(I,1)))
     *           / (T7LHI-T7LLO)
          X(I) = MIN (A7LNLO, MAX (A7LNHI, X(I)))
          ALBSOD(I,2) = X(I)
          ALBSOI(I,2) = X(I)
  200   CONTINUE
 
c          Overall icesheet surface albedos
 
        DO 250 IB=1,NBAND
!DIR$ IVDEP
          DO 252 J=1,NSOL
            I = INDSOL(J)
            ASURD(I,IB) = (1.-FI(I))*ALBSOD(I,IB)
     *                       +FI(I) *ALBSND(I,IB)
            ASURI(I,IB) = (1.-FI(I))*ALBSOI(I,IB)
     *                       +FI(I) *ALBSNI(I,IB)
  252     CONTINUE
  250   CONTINUE
 
c------------------------------
      ELSE IF (ITYPE.EQ.3) THEN
c------------------------------

#if defined (TWOSTRICE)         
        call zero (fluxiced, nvec*(nicelay+1))
        call zero (fluxicei, nvec*(nicelay+1))
#endif
 
c        Ocean strip
 
!DIR$ IVDEP
        DO 300 J=1,NSOL
          I = INDSOL(J)
 
c         Sea-ice, visible:
          X(I) = (A7IVHI*(TICE(I,1)-T7ILO) + A7IVLO*(T7IHI-TICE(I,1)))
     *           / (T7IHI-T7ILO)
          X(I) = MIN (A7IVLO, MAX (A7IVHI, X(I)))
          ALBICD(I,1) = X(I)
          ALBICI(I,1) = X(I)

#if defined (TWOSTRICE)         
c         Solve for two-stream visible-band fluxes through (cold) ice, 
c         returning surface albedo (albic[d,i]_lo) and net downward 
c         fluxes across ice layer interfaces (fluxice[d,i],in comice.h).
c         Adapted from Snwoball-Earth ebm. Here, downward incident
c         surface flux is set to 1. Returned fluxes fluxice[d,i] are 
c         scaled by actual dwnd incident fluxes in OCNFLX.
c--------------------------------
          if (fice(i).gt.0.) then
c--------------------------------
            fluxin = 1.
            cos_diffuse = cos (60.*pi/180.)       ! "mean" diffuse angle

c           direct visible flux:
            call twostr_ice_ctl (hice(i,1), nicelay, fluxin, coszen(i), 
     *        albicd_lo, zfluxd)

c           indirect (diffuse) visible flux:
            call twostr_ice_ctl (hice(i,1), nicelay, fluxin,cos_diffuse,
     *        albici_lo, zfluxi)

c           hi/lo temperature weighting as in ebm, hi albedo as above: 
            whi  = max (0., min (1.,  (tice(i,1)-(tmelt-5.))/ 5. ))
            whiv = min (whi, 0.7)
            wlov = 1.-whiv
            albicd(i,1) = a7ivhi*whiv + albicd_lo*wlov
            albici(i,1) = a7ivhi*whiv + albici_lo*wlov

c           reduce net downward penetrating flux by lo-temp factor wlov:
            do k=1,nicelay+1
              fluxiced(i,k) =  wlov*zfluxd(k)
              fluxicei(i,k) =  wlov*zfluxi(k)
            enddo
c--------------
          endif
c--------------
#endif
 
c         Sea-ice, near-IR:
          X(I) = (A7INHI*(TICE(I,1)-T7ILO) + A7INLO*(T7IHI-TICE(I,1)))
     *           / (T7IHI-T7ILO)
          X(I) = MIN (A7INLO, MAX (A7INHI, X(I)))
          ALBICD(I,2) = X(I)
          ALBICI(I,2) = X(I)
 
c         Ocean (no waveband dependence):
          CZ = COSZEN(I)
          ALBOCD(I,1)= .026/(CZ**1.7+.065) + .15*(CZ-.1)*(CZ-.5)*(CZ-1.)
          ALBOCI(I,1) = .060
          ALBOCD(I,2) = ALBOCD(I,1)
          ALBOCI(I,2) = ALBOCI(I,1)
  300   CONTINUE
 
c          Overall ocean/sea-ice surface albedos
 
        DO 350 IB=1,NBAND
!DIR$ IVDEP
          DO 352 J=1,NSOL
            I = INDSOL(J)
            ASURD(I,IB) = (1.-FICE(I))*              ALBOCD(I,IB)
     *                  +     FICE(I) * ( (1.-FI(I))*ALBICD(I,IB)
     *                                       +FI(I) *ALBSND(I,IB) )
            ASURI(I,IB) = (1.-FICE(I))*              ALBOCI(I,IB)
     *                  +     FICE(I) * ( (1.-FI(I))*ALBICI(I,IB)
     *                                       +FI(I) *ALBSNI(I,IB) )
  352     CONTINUE
  350   CONTINUE
 
c----------
      ENDIF
c----------
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SOLALB (IB)
c
c        Calculates effective albedos of the lsx surface system,
c        separately for unit incoming direct and diffuse flux. The
c        incoming direct zenith angles are supplied in commug array
c        COSZEN, and the effective albedos are returned in commug
c        arrays ASURD, ASURI. Also detailed absorbed and reflected flux
c        info is stored in com1d arrays, for later use by SOLARF.
c
c        The procedure is first to calculate the grass+soil albedos,
c        then the tree + (grass+soil+snow) albedos. The labels
c        (a) to (d) correspond to those in the LSX Description doc.
c
c        IB = waveband number.
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
c
c
c        Do nothing if all points in current strip have COSZEN le 0
c
      IF (NSOL.EQ.0) RETURN
c
c        (a) Obtain albedos, etc, for two-stream lower veg + soil
c            system, for direct and diffuse incoming unit flux
c
      DO 100 J=1,NSOL
        I = INDSOL(J)
        ASURD(I,IB) = ALBSOD(I,IB)
        ASURI(I,IB) = ALBSOI(I,IB)
  100 CONTINUE
c
      CALL TWOSTR (ABLOD,  RELOD,  FLODD,  FLODI,
     *             ASURD,  ASURI,     1,   COSZEN, IB, 0)
      CALL TWOSTR (ABLOI,  RELOI,  DUMMY,  FLOII,
     *             ASURD,  ASURI,     1,   COSZEN, IB, 1)
c
c        (b) Areally average surface albedos (lower veg, soil, snow)
c
      DO 200 J=1,NSOL
        I = INDSOL(J)
        ASURD(I,IB) = FL(I)*(1.-FI(I))*RELOD(I)
     *                + (1.-FL(I))*(1.-FI(I))*ALBSOD(I,IB)
     *                + FI(I)*ALBSND(I,IB)
        ASURI(I,IB) = FL(I)*(1.-FI(I))*RELOI(I)
     *                + (1.-FL(I))*(1.-FI(I))*ALBSOI(I,IB)
     *                + FI(I)*ALBSNI(I,IB)
  200 CONTINUE
c
c        (c) Obtain albedos, etc, for two-stream upper veg + surface
c            system, for direct and diffuse incoming unit flux
c
      CALL TWOSTR (ABUPD,  REUPD,  FUPDD,  FUPDI,
     *             ASURD,  ASURI,     2,   COSZEN, IB, 0)
      CALL TWOSTR (ABUPI,  REUPI,  DUMMY,  FUPII,
     *             ASURD,  ASURI,     2,   COSZEN, IB, 1)
c
c        (d) Calculate average overall albedos
c
!DIR$ IVDEP
      DO 400 J=1,NSOL
        I = INDSOL(J)
        ASURD(I,IB) = FU(I)*REUPD(I)
     *              + (1.-FU(I))*ASURD(I,IB)
        ASURI(I,IB) = FU(I)*REUPI(I)
     *              + (1.-FU(I))*ASURI(I,IB)
  400 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SOLARF (IB)
c
c        Calculates solar fluxes absorbed by upper and lower stories,
c        soil and snow. The incoming gcm direct and diffuse fluxes
c        are supplied in commug arrays SOLAD and SOLAI. The direct
c        zenith angles are in commug array COSZEN, and must be the same
c        as supplied earlier to SOLALB. SOLARF uses the detailed
c        unit-incoming-flux results obtained earlier by SOLALB and
c        stored in com1d arrays. The absorbed fluxes are returned in
c        com1d arrays SOL[U,S,L,G,I].
c
c        The procedure is first to calculate the upper-story absorbed
c        fluxes and fluxes below the upper story, then the lower-story
c        absorbed fluxes and fluxes below the lower story, then fluxes
c        absorbed by the soil and snow. The labels (f) to (i)
c        correspond to those in the LSX Description document.
c
c        IB = waveband number.
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
c
c
c        Do nothing if all points in current strip have COSZEN le 0
c
      IF (NSOL.EQ.0) RETURN
c
c        (f) Calculate fluxes absorbed by upper leaves and stems,
c            and downward fluxes below upper veg, using unit-flux
c            results of solalb(c). (Apportion absorbed flux between
c            leaves and stems in proportion to their LAI and SAI.)
c
!DIR$ IVDEP
      DO 600 J=1,NSOL
        I = INDSOL(J)
        X = SOLAD(I,IB)*ABUPD(I) + SOLAI(I,IB)*ABUPI(I)
        Y = LAI(I,2) / MAX (LAI(I,2)+SAI(I,2), EPSILON)
        SOLU(I) = SOLU(I) + X * Y
        SOLS(I) = SOLS(I) + X * (1.-Y)
        SOL2D(I) = SOLAD(I,IB)*FUPDD(I)
        SOL2I(I) = SOLAD(I,IB)*FUPDI(I) + SOLAI(I,IB)*FUPII(I)
  600 CONTINUE
c
c        (g) Areally average fluxes to lower veg, soil, snow
c
      DO 700 J=1,NSOL
        I = INDSOL(J)
        SOL3D(I) = FU(I)*SOL2D(I) + (1.-FU(I))*SOLAD(I,IB)
        SOL3I(I) = FU(I)*SOL2I(I) + (1.-FU(I))*SOLAI(I,IB)
  700 CONTINUE
c
c        (h,i) Calculate fluxes absorbed by lower veg, snow-free soil
c              and snow, using results of (g) and unit-flux results
c              of solalb(a)
c
!DIR$ IVDEP
      DO 800 J=1,NSOL
        I = INDSOL(J)
        SOLL(I) = SOLL(I) + SOL3D(I)*ABLOD(I) + SOL3I(I)*ABLOI(I)
c
        XD = (FL(I)*FLODD(I) + 1.-FL(I)) * SOL3D(I)
        XI = FL(I)*(SOL3D(I)*FLODI(I) + SOL3I(I)*FLOII(I))
     *       + (1.-FL(I)) * SOL3I(I)
        SOLG(I) = SOLG(I)
     *    + (1.-ALBSOD(I,IB))*XD + (1.-ALBSOI(I,IB))*XI
c
        SOLI(I) = SOLI(I)
     *            + (1.-ALBSND(I,IB))*SOL3D(I)
     *            + (1.-ALBSNI(I,IB))*SOL3I(I)
  800 CONTINUE
c
c        Estimate absorbed PARs for sunlit and shaded parts of LAI,
c        SUNP[U,L],SHAP{U,L], and sunlit fraction of LAI, SUNF[U,L]
c        (in comveg), used in TURTRA.
c
c        Direct absorbed PAR is total canopy absorbed direct flux
c        divided by total LAI+SAI (one-sided to correspond to
c        stomatal resistance and photosynthesis param. definitions).
c        Similarly for diffuse absorbed PAR, added to direct for sunlit.
c        This neglects complications due to differing values of dead vs
c        live elements, averaged into RHOVEG,TAUVEG in vegdat, and
c        modifications of omega due to intercepted snow in TWOSET.
c        Fraction of sunlit/total LAI is via SiB's Eq.(40)/Goudriaan's
c        extinction coefficeinet PAREX.
c        Do only for visible band (IB=1).
c
      IF (IB.EQ.1) THEN
!DIR$ IVDEP
        DO 900 J=1,NSOL
          I = INDSOL(J)
c         Upper canopy:
          XAIU  = MAX (0.5*(LAI(I,2)+SAI(I,2)), EPSILON)
          SUNFU(I) = (1.-EXP(-PAREX(I,2)*XAIU)) / (PAREX(I,2)*XAIU)
          SHAPU(I) = SOLAI(I,IB)*ABUPI(I)/XAIU
          SUNPU(I) = SOLAD(I,IB)*ABUPD(I)/(SUNFU(I)*XAIU) + SHAPU(I)
 
c         Lower canopy:
          XAIL  = MAX (0.5*(LAI(I,1)+SAI(I,1)), EPSILON)
          SUNFL(I) = (1.-EXP(-PAREX(I,1)*XAIL)) / (PAREX(I,1)*XAIL)
          SHAPL(I) = SOL3I(I)*ABLOI(I)/XAIL
          SUNPL(I) = SOL3D(I)*ABLOD(I)/(SUNFL(I)*XAIL) + SHAPL(I)
 
c         Account for 0.9 mu (Genesis) vs. 0.7 mu (IBIS) split:
#ifdef ibis
c         SHAPU(I) = SHAPU(I) * 46.0 / 64.7 !slevis
c         SUNPU(I) = SUNPU(I) * 46.0 / 64.7 !slevis
c         SHAPL(I) = SHAPL(I) * 46.0 / 64.7 !slevis
c         SUNPL(I) = SUNPL(I) * 46.0 / 64.7 !slevis
          SHAPU(I) = SHAPU(I) * 0.78        !slevis  ibis2.0
          SUNPU(I) = SUNPU(I) * 0.78        !slevis  ibis2.0
          SHAPL(I) = SHAPL(I) * 0.84        !slevis  ibis2.0
          SUNPL(I) = SUNPL(I) * 0.84        !slevis  ibis2.0
#endif
 
c-----
c888 Old method (no distinction btwn sunlit and shaded, 2-sided LAI+SAI)
c888      SUNPU(I)= (SOLAD(I,IB)*ABUPD(I)+SOLAI(I,IB)*ABUPI(I))
c    *             / (2.*XAIU)
c         SHAPU(I) = SUNPU(I)
c         SUNPL(I) = (SOL3D(I)*ABLOD(I) + SOL3I(I)*ABLOI(I))
c    *             / (2.*XAIL)
c         SHAPL(I) = SUNPL(I)
c-----
  900   CONTINUE
      ENDIF
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TWOSTR (ABVEG,  REFL,   FBELD,  FBELI,
     *                   ASURD,  ASURI,  IV,     COSZEN,  IB,  IFDIF)
c
c        Solves canonical radiative transfer problem of two-stream veg
c        layer + underlying surface of known albedo, for unit incoming
c        direct or diffuse flux. Returns flux absorbed within layer,
c        reflected flux, and downward fluxes below layer. Note that all
c        direct fluxes are per unit horizontal zrea, ie, already
c        including a factor cos (zenith angle).
c
c        ABVEG = flux absorbed by two-stream layer (returned)
c        REFL  = flux reflected above two-stream layer (returned)
c        FBELD = downward direct  flux below two-stream layer(returned)
c        FBELI = downward diffuse flux below two-stream layer(returned)
c        ASURD = direct  albedo of underlying surface (supplied)
c        ASURI = diffuse albedo of underlying surface (supplied)
c        IV    = 1 or 2 for lower or upper story params (supplied)
c        COSZEN= cosine of direct zenith angle (supplied, must be gt 0)
c        IB    = waveband number (supplied)
c        IFDIF = 0 for unit incoming direct flux and no diffuse, or
c                1 for vice-versa (supplied)
c
#include <compar.h>
#include <comgrd.h>
#include <comveg.h>
#include <com1d.h>
c
      DIMENSION ABVEG(NVEC),         REFL(NVEC),
     *          FBELD(NVEC),         FBELI(NVEC),
     *          ASURD(NVEC,NBAND),   ASURI(NVEC,NBAND),
     *          COSZEN(NVEC)
c
      DIMENSION OMEGA(NVEC),         BETAD(NVEC),         BETAI(NVEC),
     *          AVMU(NVEC),          GDIR(NVEC),
     *          TMP0(NVEC),          TMP1(NVEC),          TMP2(NVEC)
c
c
c        Do nothing if all points in current strip have COSZEN le 0
c
      IF (NSOL.EQ.0) RETURN
c
c        Set incoming direct and diffuse intensities DIRIN and DIFIN
c
      IF (IFDIF.EQ.0) THEN
        DIRIN = 1.
        DIFIN = 0.
      ELSE
        DIRIN = 0.
        DIFIN = 1.
      ENDIF
c
c        Calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR
c
      CALL TWOSET (OMEGA, BETAD, BETAI, AVMU, GDIR, COSZEN, IV, IB)
c
!DIR$ IVDEP
      DO 100 J=1,NSOL
        I = INDSOL(J)

c-----------------------------------------------
        if (lai(i,iv)+sai(i,iv).lt.epsilon) then
c-----------------------------------------------
c
c            If no veg, set ABVEG,etc, to exactly 0 or 1. If 
c            this is not done, roundoff error causes small(+/-)
c            SOLS, SOLL values in SOLARF and subsequent problems in 
c            TURVAP via TURTRA (also REFL is significantly non-zero(?),
c            and FBELI has non-zero roundoff error).

          fbeld(i) = 1.
          fbeli(i) = 1.
          abveg(i) = 0.
          refl(i)  = 0.

c-----------
        else
c-----------

c            Set up for general two-strm diffuse up/down fluxes y1, y2:
c                - d(y1)/dx + A y1 - B y2 = C exp ( -P x ),
c                  d(y2)/dx + A y2 - B y1 = D exp ( -P x ).
c            where x is cumulative 0.5*(LAI+SAI) / AVMU,
c            P,Q are e-folding rates of the particular and complementary
c            solutions (guaranteed gt 0 if OMEGA ne 1),
c            Y1, Y2 (*exp(-Px)) is the particular solution.
c
c            The determinant A*A - B*B - P*P for the particular solution
c            can be zero (ie, P = Q). Mathematically, would have other
c            form of particular soln, involving x*exp(-Qx). Handle below
c            by perturbing P, which gives same results to o(delta P)...
c            but the resulting solution involves the difference between
c            two large quantities, so beware of truncation error.
c
          A = 1. - OMEGA(I) * (1.-BETAI(I))
          B = OMEGA(I) * BETAI(I)
          TMP0(I) = A*A-B*B
          Q = SQRT (TMP0(I))
          P = AVMU(I) * GDIR(I) / MAX (COSZEN(I), .01)
c
c         Next line perturbs P if P = Q (as per comment above)
          P = CVMGT ( P, (1.+SIGN(.001,P-Q)) * P, ABS(P-Q).GE..001*P )
c
          C0 = OMEGA(I) * P * DIRIN
          C = C0 * BETAD(I)
          D = C0 * (1.-BETAD(I))
          Y1 = ((A-P)*C + B*D) / (A*A - B*B - P*P)
          Y2 = (B*C + (A+P)*D) / (A*A - B*B - P*P)
c
c            Find CON1 and CON2, the constants multiplying the
c            complementary solns to satisfy the boundary conditions
c                F1*CON1 + F2*CON2 = G1,      F3*CON1 + F4*CON2 = G2.
c            The determinant DET will be non-zero unless Q=0 or B=0
c            (ie, unless OMEGA = 1, OMEGA = 0 or BETAI = 0).
c
          TMP1(I) = P * 0.5*(LAI(I,IV)+SAI(I,IV)) / AVMU(I)
          TMP2(I) = Q * 0.5*(LAI(I,IV)+SAI(I,IV)) / AVMU(I)
          EP = EXP (-TMP1(I))
          EQ = EXP (-TMP2(I))
          F1 = A + Q
          F2 = A - Q
          F3 = ( B - ASURI(I,IB) * (A + Q) ) * EQ
          F4 = ( B - ASURI(I,IB) * (A - Q) ) / EQ
          DET = F1*F4 - F2*F3
          G1 = -Y2 + DIFIN
          G2 =(-Y1 + ASURI(I,IB)*Y2 + ASURD(I,IB)*DIRIN) * EP
          CON1 = ( G1*F4 - G2*F2 ) / DET
          CON2 = (-G1*F3 + G2*F1 ) / DET
c
c            Save downward direct, diffuse fluxes below two-stream layer
c
          FBELD(I) = DIRIN*EP
          FBELI(I) = CON1*(A+Q)*EQ + CON2*(A-Q)/EQ + Y2*EP
c
c            Save reflected flux, and flux absorbed by two-stream layer
c
          REFL(I) = (CON1 + CON2)*B + Y1
          ABSUR = (1.-ASURD(I,IB)) * FBELD(I)
     *          + (1.-ASURI(I,IB)) * FBELI(I)
          ABVEG(I) = DIRIN + DIFIN - REFL(I) - ABSUR

c------------
        endif
c------------

  100 CONTINUE
c
c        Compute PAR extinction coefficient (SiB Eq.40, Goudriaan,1970)
c        for PAR calculations in SOLARF. Do only for visible (IB=1)
c        and direct (IFDIF=0). (GDIR depends on IV and IB, OMEGA
c        depends on IV, neither depend on IFDIF).
c
      IF (IB.EQ.1.AND.IFDIF.EQ.0) THEN
!DIR$ IVDEP
         DO 300 J=1,NSOL
            I = INDSOL(J)
            PAREX(I,IV) = GDIR(I)*SQRT(1.-OMEGA(I))/MAX(COSZEN(I),0.01)
  300    CONTINUE
      ENDIF
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TWOSET (OMEGA, BETAD, BETAI, AVMU, GDIR, COSZEN, IV,IB)
c
c        Sets two-stream parameters, given single-element transmittance
c        and reflectance, leaf orientation weights, and cosine of the
c        zenith angle. Then adjusts for amounts of intercepted snow.
c
c        The two-stream parameters OMEGA,BETAD,BETAI are weighted
c        combinations of the "exact" values for the 3 orientations:
c        all vertical, all horizontal, or all random (ie, spherical).
c        The vertical,horizontal weights are in ORIEV,ORIEH (comveg).
c        The "exact" expressions are as derived in my notes(8/6/91,p.6).
c        Note that values for OMEGA*BETAD and OMEGA*BETAI are calculated
c        and then divided by the new OMEGA, since those products are
c        actually used in TWOSTR. Also those depend *linearly* on the
c        single-element transmittances and reflectances TAUVEG,RHOVEG,
c        which are themselves linear weights of leaf and stem values
c        (set in VEGDAT).
c
c        For random orientation, OMEGA*BETAD depends on COSZEN according
c        to the function in array TABLEMU, generated by a standalone
c        program "leaf" (stored at the end of this file.)
c
c        The procedure is approximate since OMEGA*BETA[D,I] and GDIR
c        should depend non-linearly on the complete leaf-angle
c        distribution. Then we should also treat leaf and stem angle
c        distributions separately, and allow for the cylindrical
c        shape of stems (Norman and Jarvis, App.B; the expressions
c        below are appropriate for flat leaves).
c
c        All quantities are returned unless otherwise noted:
c        OMEGA  = fraction of intercepted radiation that is scattered
c        BETAD  = fraction of scattered *direct* radiation that is
c                 scattered into upwards hemisphere
c        BETAI  = fraction of scattered downward *diffuse* radiation
c                 that is scattered into upwards hemisphere (or fraction
c                 of scattered upward diffuse rad. into downwards hemis)
c        AVMU   = average diffuse optical depth
c        GDIR   = average projected leaf area into solar direction
c        COSZEN = cosine of solar zenith angle (supplied)
c        IV     = 1 or 2 for lower or upper story params (supplied)
c        IB     = waveband number (supplied)
c
#include <compar.h>
#include <comgrd.h>
#include <comveg.h>
#include <com1d.h>
c
      DIMENSION OMEGA(NVEC),         BETAD(NVEC),        BETAI(NVEC),
     *          AVMU(NVEC),          GDIR(NVEC),         COSZEN(NVEC)
c
      DIMENSION OTMP(NVEC)
c
      PARAMETER (NTMU=100)
      DIMENSION TABLEMU(NTMU+1)
      SAVE TABLEMU
c
      DATA TABLEMU /
     *   0.5000, 0.4967, 0.4933, 0.4900, 0.4867, 0.4833, 0.4800, 0.4767,
     *   0.4733, 0.4700, 0.4667, 0.4633, 0.4600, 0.4567, 0.4533, 0.4500,
     *   0.4467, 0.4433, 0.4400, 0.4367, 0.4333, 0.4300, 0.4267, 0.4233,
     *   0.4200, 0.4167, 0.4133, 0.4100, 0.4067, 0.4033, 0.4000, 0.3967,
     *   0.3933, 0.3900, 0.3867, 0.3833, 0.3800, 0.3767, 0.3733, 0.3700,
     *   0.3667, 0.3633, 0.3600, 0.3567, 0.3533, 0.3500, 0.3467, 0.3433,
     *   0.3400, 0.3367, 0.3333, 0.3300, 0.3267, 0.3233, 0.3200, 0.3167,
     *   0.3133, 0.3100, 0.3067, 0.3033, 0.3000, 0.2967, 0.2933, 0.2900,
     *   0.2867, 0.2833, 0.2800, 0.2767, 0.2733, 0.2700, 0.2667, 0.2633,
     *   0.2600, 0.2567, 0.2533, 0.2500, 0.2467, 0.2433, 0.2400, 0.2367,
     *   0.2333, 0.2300, 0.2267, 0.2233, 0.2200, 0.2167, 0.2133, 0.2100,
     *   0.2067, 0.2033, 0.2000, 0.1967, 0.1933, 0.1900, 0.1867, 0.1833,
     *   0.1800, 0.1767, 0.1733, 0.1700, 0.1667 /
c
      SAVE A7SVLO,A7SVHI,A7SNLO,A7SNHI,T7SNO,BETADSNO,BETAISNO
      DATA A7SVLO, A7SVHI /0.90, 0.55/
      DATA A7SNLO, A7SNHI /0.60, 0.35/
      DATA T7SNO /5./
      DATA BETADSNO, BETAISNO /0.5, 0.5/
c
c
c        Set two-stream parameters OMEGA, BETAD, BETAI, GDIR and AVMU
c        as weights of those for 100% vert,horiz,random orientations.
c
!DIR$ IVDEP
      DO 100 J=1,NSOL
        I = INDSOL(J)
c
        ZRHO = RHOVEG(I,IB,IV)
        ZTAU = TAUVEG(I,IB,IV)
c       Weight for random orientation is 1 - those for vert and horiz
        ORAND = 1. - ORIEV(I,IV) - ORIEH(I,IV)
c
        OMEGA(I) = ZRHO + ZTAU
c
c       ZTAB is transmittance coeff. for random-orientation OMEGA*BETAD,
c       given by TABLEMU as a function of COSZEN.
        ITAB = NINT (COSZEN(I)*NTMU + 1)
        ZTAB = TABLEMU(ITAB)
        BETAD(I) = (  ORIEV(I,IV) * 0.5*(ZRHO + ZTAU)
     *              + ORIEH(I,IV) * ZRHO
     *              + ORAND       * ((1.-ZTAB)*ZRHO + ZTAB*ZTAU) )
     *             / OMEGA(I)
c
        BETAI(I) = (  ORIEV(I,IV) * 0.5*(ZRHO + ZTAU)
     *              + ORIEH(I,IV) * ZRHO
     *              + ORAND       * ((2./3.)*ZRHO + (1./3.)*ZTAU) )
     *             / OMEGA(I)
c
        GDIR(I)  = ORIEV(I,IV) * (2./PI) * SQRT (1.-COSZEN(I)*COSZEN(I))
     *           + ORIEH(I,IV) * COSZEN(I)
     *           + ORAND       * 0.5
 
        AVMU(I) = 1.
c
  100 CONTINUE
c
c        Adjust OMEGA, BETAD and BETAI for amounts of intercepted snow
c        (omega of snow ramps between -5 to 0 deg C, as in solsur).
 
      IF (IB.EQ.1) THEN
        ZA7LO = A7SVLO
        ZA7HI = A7SVHI
      ELSE
        ZA7LO = A7SNLO
        ZA7HI = A7SNHI
      ENDIF
 
      IF (IV.EQ.1) THEN
c
c       Lower story
!DIR$ IVDEP
        DO 210 J=1,NSOL
          I = INDSOL(J)
          Y = FWETL(I)*(1.-RLIQL(I))
          Z = MAX (0., MIN (1., (TMELT-TL(I))/T7SNO ))
          O = Z*ZA7LO + (1.-Z)*ZA7HI
          OTMP(I)  = OMEGA(I)
          OMEGA(I) =  (1-Y)*OTMP(I)          + Y*O
          BETAD(I) = ((1-Y)*OTMP(I)*BETAD(I) + Y*O*BETADSNO) / OMEGA(I)
          BETAI(I) = ((1-Y)*OTMP(I)*BETAI(I) + Y*O*BETAISNO) / OMEGA(I)
  210   CONTINUE
c
      ELSE
c
c       Upper story
!DIR$ IVDEP
        DO 220 J=1,NSOL
          I = INDSOL(J)
          X = LAI(I,IV) / MAX (LAI(I,IV)+SAI(I,IV), EPSILON)
          Y = X * FWETU(I)*(1.-RLIQU(I)) + (1-X) *FWETS(I)*(1.-RLIQS(I))
          ZT = X*TU(I) + (1.-X)*TS(I)
          Z = MAX (0., MIN (1., (TMELT-ZT)/T7SNO ))
          O = Z*ZA7LO + (1.-Z)*ZA7HI
          OTMP(I)  = OMEGA(I)
          OMEGA(I) =  (1-Y)*OTMP(I)          + Y*O
          BETAD(I) = ((1-Y)*OTMP(I)*BETAD(I) + Y*O*BETADSNO) / OMEGA(I)
          BETAI(I) = ((1-Y)*OTMP(I)*BETAI(I) + Y*O*BETAISNO) / OMEGA(I)
  220   CONTINUE
c
      ENDIF
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE IRRAD
 
c        Calculates overall emitted IR flux, and net absorbed minus
c        emitted IR fluxes for upper leaves, upper stems, lower story,
c        soil and snow. Assumes upper leaves, upper stems and lower
c        story each form a semi-transparent plane, with the upper-leaf
c        plane just above the upper-stem plane. The solid surfaces
c        (snow, soil) can have emissivities < 1.
c
c        The gcm incoming flux is supplied in commug array FIRA.
c        The emitted IR flux by overall surface system is returned in
c        com1d array FIRB. The IR fluxes absorbed by upper leaves,
c        upper stems, lower veg, soil and snow are returned in com1d
c        arrays FIRU, FIRS, FIRL, FIRG and FIRI.
c
c        Other com1d arrays used are:
c        EMU, EMS, EML  = emissivities of the vegetation planes
c        FUP, FDOWN     = upward and downward fluxes below tree level
 
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <comsno.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
 
c     Set emissivities of soil,snow,icesheet,seaice,water
c     (in com1d for use in ocnflx and sheetflx)
      DATA EMISOIL, EMISNOW, EMILICE, EMISICE, EMIOCN
#ifdef icedrive
     *    /    1.0,     1.0,     1.0,     1.0,    1.0 /
#else
c    *    /    .95,     .95,     1.0,     1.0,    1.0 / ! v2
     *    /    .85,     .98,     .98,     .98,    1.0 / ! v3
#endif
 
c     Local arrays:
c     EMU   = IR emissivity of upper-leaves veg plane
c     EMS   = IR emissivity of upper-stems veg plane
c     EML   = IR emissivity of lower-story veg plane
c     EMG   = IR emissivity (gray) of soil surface
c     EMI   = IR emissivity (gray) of snow surface
c     FDOWN = downward IR flux below tree level per overall area
c     FUP   = upward   IR flux below tree level per overall area
c     FDOWNG= downward IR flux below lower-story veg
c     FUPG  = upward   IR flux below lower-story veg
c     FUPGB = upward   IR flux above bare soil surface
c     FUPI  = upward   IR flux above snow surface
      DIMENSION
     *  EMU(NVEC),   EMS(NVEC),   EML(NVEC),   EMG(NVEC),  EMI(NVEC),
     *  FDOWN(NVEC), FDOWNG(NVEC),
     *  FUP(NVEC),   FUPG(NVEC),  FUPGB(NVEC), FUPI(NVEC)
 
c     Use uniform value 1.0 for average diffuse optical depth
c     (although an array for solar, all values are set to 1 in TWOSET).
      SAVE AVMUIR
      DATA AVMUIR /1./
 
      DO 100 I=1,NPOI
 
        EMU(I) = 1. - EXP ( -0.5*LAI(I,2) / AVMUIR )
        EMS(I) = 1. - EXP ( -0.5*SAI(I,2) / AVMUIR )
        EML(I) = 1. - EXP ( -0.5*(LAI(I,1)+SAI(I,1)) / AVMUIR )
 
        EMG(I) = EMISOIL
        EMI(I) = EMISNOW
 
        FDOWN(I) =  (1.-FU(I)) * FIRA(I)
     *            + FU(I) * ( (1.-EMU(I))*(1.-EMS(I))*FIRA(I)
     *                       +    EMU(I)* (1.-EMS(I))*STEF*(TU(I)**4)
     *                       +    EMS(I)*STEF*(TS(I)**4) )
 
        FDOWNG(I) = (1.-EML(I))*FDOWN(I)  + EML(I)*STEF*(TL(I)**4)
 
        FUPG(I)   = (1.-EMG(I))*FDOWNG(I) + EMG(I)*STEF*(TG(I)**4)
 
        FUPGB(I)  = (1.-EMG(I))*FDOWN(I)  + EMG(I)*STEF*(TG(I)**4)
 
        FUPI(I)   = (1.-EMI(I))*FDOWN(I)  + EMI(I)*STEF*(TI(I)**4)
 
        FUP(I) = (1.-FI(I))*(      FL(I)*(       EML(I) *STEF*(TL(I)**4)
     *                                     + (1.-EML(I))*FUPG(I) )
     *                        +(1.-FL(I))*FUPGB(I)
     *                      )
     *         +     FI(I) * FUPI(I)
 
        FIRB(I) =   (1.-FU(I)) * FUP(I)
     *            + FU(I)  * ( (1.-EMU(I))*(1.-EMS(I))*FUP(I)
     *                        +    EMU(I)*STEF*(TU(I)**4)
     *                        +    EMS(I)*(1.-EMU(I))*STEF*(TS(I)**4) )
 
        FIRU(I) =   EMU(I)*EMS(I)*STEF*(TS(I)**4)
     *            + EMU(I)*(1.-EMS(I))*FUP(I)
     *            + EMU(I)*FIRA(I)
     *            - 2*EMU(I)*STEF*(TU(I)**4)
 
        FIRS(I) =   EMS(I)*EMU(I)*STEF*(TU(I)**4)
     *            + EMS(I)*FUP(I)
     *            + EMS(I)*(1.-EMU(I))*FIRA(I)
     *            - 2*EMS(I)*STEF*(TS(I)**4)
 
        FIRL(I) =   EML(I)*FDOWN(I)
     *            + EML(I)*FUPG(I)
     *            - 2*EML(I)*STEF*(TL(I)**4)
 
        FIRG(I) =       FL(I)  * (FDOWNG(I) - FUPG(I))
     *            + (1.-FL(I)) * (FDOWN(I)  - FUPGB(I))
 
        FIRI(I) =   FDOWN(I) - FUPI(I)
 
  100 CONTINUE
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine leafdum
c     program    leafdum
c
c        This is a separate program, stored here for convenience.
c        It is not to be linked into the lsx program.
 
c        Computes G (projected area), "table" factor, and omega*beta
c        for *direct* incident radiation versus solar zenith angle,
c        for any given leaf angle distribution. (See notes,8/6/91,p.5).
c        Also writes table values in Fortran data statement form to
c        be used in lsx subroutine TWOSET.
 
c        Arrays versus solar zenith angle (dimensioned ntab):
c        cosz  = cosine (solar zenith angle) (1 = overhead)
c        the   = solar zenith angle
c        table = weighting vs zenith angle, such that:
c                omega*beta = trans*table + refl*(1-table)
c        g     = projected area of leaves in solar direction
c        omb   = omega*beta (fraction of intercepted direct radiation
c                            that is scattered into upward hemisphere)
 
c        Arrays versus leaf distribution angle (dimensioned npoi):
c        alf = leaf angle (between leaf normal and vertical)
c        adis = fraction of leaves between alf-.5*dalf and alf+.5*dalf
 
c        Other variables:
c        tadis = integral of leaf angle distribution (should be 1)
c        psi   = pi - azimuth angle beyond which leaves are illuminated
c                on their bottom side
c        trans = single leaf transmittance
c        refl  = single leaf reflectance
 
      parameter (ntab=101, npoi=1000)
      dimension cosz(ntab), the(ntab), table(ntab), g(ntab), omb(ntab),
     *          alf(npoi), adis(npoi)
 
      parameter (pi = 3.14159265358979)
 
      pio2 = pi/2
 
c        Set arrays that depend on solar zenith angle
 
      do 10 i=1,ntab
        cosz(i) = float(i-1) / float(ntab-1)
        the(i) = acos(cosz(i))
        table(i) = 0.
        g(i) = 0.
   10 continue
 
c        Set arrays that depend on leaf distribution angle
 
      tadis = 0.
      dalf = pio2/npoi
      do 20 i=1,npoi
        alf(i) = ((i-0.5)/npoi)*pio2
c       Spherical distribution
        adis(i) = sin(alf(i))
c       Uniform vs angle distribution
c       adis(i) = 1./pio2
        tadis = tadis + adis(i)*dalf
   20 continue
 
c        Loop over solar zenith angles
 
      do 100 it=1,ntab
 
c          Loop (for integrations) over leaf distribution angles
 
        do 200 ia=1,npoi
          if (the(it)+alf(ia).le.pio2+.00001) then
            ga = cos(alf(ia))*cos(the(it))
            gb = 0
          else
            psi = acos (1./(tan(alf(ia))*tan(the(it))))
            gb = (  sin(alf(ia))*sin(the(it))*sin(psi)
     *            - cos(alf(ia))*cos(the(it))*psi      ) / pi
            ga = gb + cos(alf(ia))*cos(the(it))
          endif
 
          g(it) = g(it) + adis(ia) * (ga + gb) * dalf
          table(it) = table(it) + adis(ia) * (  gb*0.5*(1+cos(alf(ia)))
     *                                        + ga*0.5*(1-cos(alf(ia))))
     *                                     * dalf
  200   continue
 
        table(it) = table(it) / g(it)
        g(it) = g(it) / tadis
  100 continue
 
c        Prompt for trans and refl, and printout versus solar angle
 
      do 800 it=1,ntab
        if (it.eq.1) then
          write(6,810)
  810     format(/'Enter trans, refl : ',$)
          read(*,*) trans, refl
          write(6,812) trans, refl, tadis
  812     format(  'trans=',f5.3,'  refl=',f5.3,'  tadis=',f11.8
     *           //'   cosz    the         g     table       omb'/)
        endif
 
        omb(it) = trans*table(it) + refl*(1.-table(it))
        write(6,814) cosz(it), the(it)*180/pi, g(it), table(it), omb(it)
  814   format(2f7.3,3f10.6)
  800 continue
 
c       Compute mean omega*beta, weighted by sin(zenith angle), to
c       get omega*beta for isotropic diffuse downward radiation.
c       Result (ombdif) should be (2/3)*refl + (1/3)*trans.
 
      ombdif = 0.
      dthe = pio2/(ntab-1)
      do 900 it=2,ntab-1
        ombdif = ombdif + sin(the(it)) * omb(it) * dthe
  900 continue
      ombdif = ombdif + sin(the(1))    * omb(1)    * dthe/2.
      ombdif = ombdif + sin(the(ntab)) * omb(ntab) * dthe/2.
      write(6,902) ombdif
  902 format(/'ombdif = ',f10.6)
 
c       Write table values in fortran data statement form
 
      write (6,910) table
  910 format(//'      data table /'
     *       /('     *  ',8(f7.4,',')) )
 
      call endrun (0)
      end

c-----------------------------------------------------------------------

c++++++++++++++++++++++
#if defined (TWOSTRICE)
c++++++++++++++++++++++

      subroutine twostr_ice_ctl (hice, nice, visflux, u0, 
     *                           albavg, fluxice)

c     Calls twostr (Jim Kasting, 01/04) to solve delta-Eddington  
c     visible-waveband solar fluxes in a bubbly ice layer with absorbing
c     substrate (pure ice) with absorption coefficient ka, and bubbles
c     giving a scattering coefficient ks. The "single"-scattering
c     albedo w0s = ks/(ks+ka), and ks is specified via w0s in twostr
c     (and ka, and asymmetry factor of bubbles are hard coded).  
c     See Mullen and Warren,1988,JGR,93,D7,8403, and Warren et al.,
c     2002,JGR,107,C10,3167. Specular surface reflections are ignored.
c
c     Twostr uses NZ layers, and returns
c     flux profiles (DIRECT, FNET(net downwards), FUP, FDN) at 
c     layer boundaries (1=top, nz+1=bottom) for a given direct 
c     incident beam with cos(zenith angle) = U0 . The incident flux at 
c     the surface per unit horizontal area (DIRECT(1)) is set to U0
c     in twostr. All flux quantities are then multiplied by VISFLUX,
c     the supplied downward incident flux at the surface.
c
c     Supplied: 
c     hice    = ice-model layer thicknesses (m)
c     nice    = number of ice-model layers
c     visflux = downward incident flux at surface (W/m2)
c     u0      = cosine (incident angle) of visflux
c     Returned:
c     albavg  = surface albedo
c     fluxice= net downward flux across ice model layer interfaces
c               (interpolated from nz grid to ice-model grid) 

      dimension hice(nice), fluxice(nice+1)

      parameter (nz=50, nz1=nz+1, nz2=2*nz)          ! same as in twostr
      dimension 
     *  direct(nz1), fnet(nz1), fup(nz1), fdn(nz1), dz(nz), z(nz1)

c     single-scattering albedo of ice:
      parameter (w0s = .999)
c     parameter (w0s = .994)

c     absorption coefficient of ice (1/m)
      parameter (viskap = 0.10)
c     parameter (viskap = 0.07)

c        Set total ice thickness (h)

      h = 0.
      do m=1,nice
        h = h + hice(m)
      enddo

c        Set adaptive layer thicknesses (dz(nz)) and boundary levels
c        (z(nz1)). z coordinate increases downwards from surface at
c        z=0 to base at z=h. Thinner layers near top, tending to 
c        h/nz at base. (And index k increases downward too).

      zhtot = 0.
      kc = nint(0.3*nz)
      do k=1,nz
c       zwei = sqrt(exp(-(k-1.)/(0.3*nz)))
c       dz(k) = min(0.5,h/nz)*zwei + (h/nz)*(1.-zwei)
        if (k.le.kc) then
          dz(k) = min(0.5,h/nz)
        else
          zwei = 1. - float(k-kc)/(nz-kc)
          dz(k) = zwei*min(0.5,h/nz) + (1.-zwei)*(h/nz)
        endif
        zhtot = zhtot + dz(k)
      enddo
      z(1) = 0.
      do k=1,nz
        dz(k) = dz(k)*h/zhtot
        z(k+1) = z(k) + dz(k)
      enddo
      z(nz1) = h  ! to avoid possible roundoff errors

      CALL TWOSTR_ICE (U0,W0S,VISKAP,H,TAUB,SRFALB,TRANS,
     *                 direct, fnet, fup, fdn, dz)

c        Incident flux direct(1) has been set to u0 in twostr, 
c        accounting for incidence through a unit horizontal area

      albavg = fup(1)/direct(1)
      traavg = fdn(nz1)/direct(1)
      absavg = (fnet(1)-fnet(nz1))/direct(1)
      erravg = albavg + traavg + absavg - 1.
      solpenb= traavg * visflux

c        Rescale downward flux (fnet) by dwnd incident flux (visflux)

      do k=1,nz1
        fnet(k) = fnet(k)*(visflux/direct(1))
      enddo

c        Interpolate net downward fluxes to ice grid (fluxice)

      fluxice(1) = fnet(1) 

      zice = hice(1)
      do m=2,nice  
        do k=2,nz1
          if (z(k).ge.zice) then
            zwei = (z(k)-zice)/(z(k)-z(k-1))
            fluxice(m) = zwei*fnet(k-1) + (1.-zwei)*fnet(k) 
            goto 100
          endif
        enddo
  100   continue
        zice = zice + hice(m)
      enddo

      fluxice(nice+1) = fnet(nz1) 

c        Print results

c     write (6,'(/a)') 
c    *  '     w0s       h     alb     tra     abs    penb'
c     write (6,'(f8.3, f8.2, 4f8.3)')
c    *  w0s, h, albavg, traavg, absavg, solpenb

c     write (6,'(/a)') 
c    *   '   k  direct     fup     fdn    fnet       z      dz'
c     do k=1,nz1
c       write(6,'(i4, 6f8.3)') 
c    *    k, direct(k), fup(k), fdn(k), fnet(k), z(k), dz(min(k,nz))  
c     enddo

c     write (6,'(/a)') 
c    *   '   m    zice  fluxice'
c     zice = 0.
c     do m=1,nice+1
c       write (6,'(i4,f8.3,f10.3)')  m, zice, fluxice(m)
c       if (m.le.nice) zice = zice + hice(m)
c     enddo

      if (abs(erravg).gt.1.e-6) then
        write (6,'(/a,f15.10)') 
     *    '*** Error (twostr_ice_ctl): erravg=', erravg
        stop
      endif

      return
      end

c-----------------------------------------------------------------------

      SUBROUTINE TWOSTR_ICE (U0,W0S,VISKAP,H,TAUB,SRFALB,TRANS,
     *                       direct, fnet, fup, fdn, dz)

c     Code from Jim Kasting (1/04). 
c     Ref: Toon et al.,1989,JGR,94,16,287-16,301.

c     Significant modifications (dp, 1/04): 
c       Sign of FNET, from +ve upwards to +ve downwards.
c       VISKAP passed as an argument.
c       Non-uniform vertical grid (DZ).

      parameter (nz=50, nz1=nz+1, nz2=2*nz)  ! same as in twostr_ice_ctl
      DIMENSION DIRECT(NZ1),FNET(NZ1),FUP(NZ1),FDN(NZ1), dz(nz)
C
C   Input parameters
C     U0 = cosine of solar zenith angle
C     W0S = single scattering albedo
C     H = ice thickness (in meters)
C
C   Output parameters
C     TAUB = total extinction optical depth
C     SRFALB = surface albedo
C     TRANS = ice transmissivity
C
C   This version of TWOSTR has been modified to include U0 in the
C   definition of DIRECT. 
C
C   This is my version of the Toon et al. 2-stream code.  (Ref.: JGR
C   94, 16287, 1989). It vectorizes over height, rather than wavelength,
C   and is designed to work with PRIMS3 and its companion photochemical 
C   models.
C
C   For now, at least, it is hardwired as the quadrature approximation.
C   The Eddington coefficients are in the code but are commented out.
C
C   Note that I have not included any logic to prevent W0 (the single
C   scattering albedo) from getting close to unity. In Tom's code it has
C   a maximum value of 0.999 .
C
      DIMENSION TAU(NZ), TAUC(NZ), G(NZ), GAM1(NZ), GAM2(NZ), GAM3(NZ),
     1  GAM4(NZ), ALAM(NZ), CGAM(NZ), E1(NZ), E2(NZ), E3(NZ), E4(NZ),
     2  CP0(NZ), CPB(NZ), CM0(NZ), CMB(NZ), Y1(NZ), Y2(NZ), AMEAN(NZ1),
     3  W0(NZ)
      DIMENSION A(NZ2), B(NZ2), D(NZ2), E(NZ2), Y(NZ2)
C
C   Parameters:
C     ALB = albedo at the base of the ice (the ice-water interface)
C     GS = asymmetry factor in the visible region
C     VISKAP = absorption coefficient of ice in the visible (1/m)
      ALB = 0.
      GS = 0.85
c     VISKAP = 0.1     ! dp (passed in as argument)
C
C   Calculate total absorption optical depth (TAUABS) and extinction
C   optical depth (TAUB)
      TAUABS = VISKAP * H
      TAUB = TAUABS/(1.-W0S)
C
C   U0 and TAUB are the zenith angle and total optical depth.  These are
C   just trial values.  ALB is the surface albedo, TAUB is the total 
C   extinction optical depth, W0S and GS are unique values of the single
C   scattering albedo and asymmetry factor that are used to load the 
C   altitude-dependent vectors W0 and G.  This can all be changed.
C
C      U1 = 0.5
C  (Eddington value above)
      SQ3 = SQRT(3.)
      U1 = 1./SQ3
C  (Quadrature value above)
      U0M = 1./U0
      U0M2 = U0M*U0M
      U1M = 1./U1
      NZM1 = NZ - 1
      NZP1 = NZ + 1
      MZ2 = NZ2
      PI = 3.14159265358979
C
C   For now, just make up an optical depth grid.
C   First, apply delta scaling
      F = GS * GS
      GP = GS/(1. + GS)
      TAUP = (1. - W0S*F)*TAUB
      W0P = (1. - F)*W0S/(1. - W0S*F)
c
      DO 1 N=1,NZ
c     TAU(N) = TAUP/NZ                ! dp
      TAU(N) = TAUP * DZ(N)/H         ! dp
      W0(N) = W0P
   1  G(N) = GP
C   Temporary above here
C
C   Calculate the gamma's, lambda's, and e's
      DO 2 N=1,NZ
C      GAM1(N) = (7. - W0(N)*(4.+3.*G(N)))/4.
C      GAM2(N) = - (1. - W0(N)*(4.-3.*G(N)))/4.
C      GAM3(N) = (2. - 3.*G(N)*U0)/4.
C   (Eddington values above -- not currently used)
C
      GAM1(N) = SQ3*(2. - W0(N)*(1.+G(N)))/2.
      GAM2(N) = SQ3*W0(N)*(1.-G(N))/2.
      GAM3(N) = (1. - SQ3*G(N)*U0)/2.
      GAM4(N) = 1. - GAM3(N)
C   (Quadrature values above)
C
      ALAM(N) = SQRT(GAM1(N)*GAM1(N) - GAM2(N)*GAM2(N))
      CGAM(N) = (GAM1(N) - ALAM(N))/GAM2(N)
      EMLT = EXP(-ALAM(N)*TAU(N))
C
      E1(N) = 1. + CGAM(N)*EMLT
      E2(N) = 1. - CGAM(N)*EMLT
      E3(N) = CGAM(N) + EMLT
   2  E4(N) = CGAM(N) - EMLT
C
C   Calculate A, B, and D, i.e. the coefficients of the tridiagonal
C      matrix
C   Top of atmosphere
      A(1) = 0.
      B(1) = E1(1)
      D(1) = -E2(1)
C
C   Odd coefficients
      DO 3 N=1,NZM1
      L = 2*N + 1
      A(L) = E2(N)*E3(N) - E4(N)*E1(N)
      B(L) = E1(N)*E1(N+1) - E3(N)*E3(N+1)
   3  D(L) = E3(N)*E4(N+1) - E1(N)*E2(N+1)
C
C   Even coefficients
      DO 4 N=1,NZM1
      L = 2*N
      A(L) = E2(N+1)*E1(N) - E3(N)*E4(N+1)
      B(L) = E2(N)*E2(N+1) - E4(N)*E4(N+1)
   4  D(L) = E1(N+1)*E4(N+1) - E2(N+1)*E3(N+1)
C
C   Bottom of atmosphere
      A(NZ2) = E1(NZ) - ALB*E3(NZ)
      B(NZ2) = E2(NZ) - ALB*E4(NZ)
      D(NZ2) = 0.
C
C   Now, set up the RHS of the equation:
C   TAUC(N) is the optical depth above layer N
      TAUC(1) = 0.
      DO 5 N=2,NZ
   5  TAUC(N) = TAUC(N-1) + TAU(N-1)
C
      DIRECT(1) = U0

      DO 6 N=1,NZ
      ET0 = EXP(-TAUC(N)/U0)
      ETB = ET0 * EXP(-TAU(N)/U0)
      DIRECT(N+1) = U0*ETB
      DENOM = ALAM(N)*ALAM(N) - U0M2
      FACP = W0(N) * ((GAM1(N)-U0M)*GAM3(N) + GAM4(N)*GAM2(N))
      FACM = W0(N) * ((GAM1(N)+U0M)*GAM4(N) + GAM2(N)*GAM3(N))
C
      CP0(N) = ET0*FACP/DENOM
      CPB(N) = ETB*FACP/DENOM
      CM0(N) = ET0*FACM/DENOM
   6  CMB(N) = ETB*FACM/DENOM
      SSFC = ALB*DIRECT(NZ1)
C   DIRECT(N) is the direct solar flux at the top of layer N
C
C   Odd coefficients
      E(1) = - CM0(1)
      DO 7 N=1,NZM1
      L = 2*N + 1
   7  E(L) = (CP0(N+1)-CPB(N))*E3(N) + (CMB(N)-CM0(N+1))*E1(N)
C
C   Even coefficients
      DO 8 N=1,NZM1
      L = 2*N
   8  E(L) = (CP0(N+1)-CPB(N))*E2(N+1) - (CM0(N+1)-CMB(N))*E4(N+1)
      E(NZ2) = SSFC - CPB(NZ) + ALB*CMB(NZ)
C
C   Call the tridiagonal solver.  Use Numerical Recipes for now.
      CALL TRIDAG(A,B,D,E,Y,MZ2)
C
      DO 9 N=1,NZ
      L = 2*N
      L1 = L-1
      Y1(N) = Y(L1)
   9  Y2(N) = Y(L)
C
C   Calculate the mean intensities, AMEAN(N), at the boundaries between
C   the layers.  AMEAN(N) is the intensity at the top of layer N.
C   (Top of atmosphere is different)
      AMEAN(1) = U1M * (Y1(1)*E3(1) - Y2(1)*E4(1) + CP0(1)) + 1.
      F1P0 = (Y1(1)*E3(1) - Y2(1)*E4(1) + CP0(1)) * PI
      DO 110 N=1,NZ
 110  AMEAN(N+1) = U1M * (Y1(N)*(E1(N)+E3(N)) + Y2(N)*(E2(N)+E4(N)) 
     1  + CPB(N) + CMB(N)) + DIRECT(N+1)/U0
C
C      PRINT 210, AMEAN(1),AMEAN(NZ+1)
 210  FORMAT(/1X,'AMEAN(1) =',1PE10.3,5X,'AMEAN(NZ1) =',E10.3)
C
C  Calculate upward and downward fluxes
      FUP(1) = ((Y1(1)*E3(1) - Y2(1)*E4(1)) + CP0(1))
      FDN(1) = DIRECT(1)
c     FNET(1) = FUP(1) - FDN(1)  ! dp
      fnet(1) = fdn(1) - fup(1)  ! dp
      DO 10 N=1,NZ
         FUP(N+1) = (Y1(N)*E1(N) + Y2(N)*E2(N)
     &     + CPB(N))
         FDN(N+1) = (Y1(N)*E3(N) + Y2(N)*E4(N)
     &     + CMB(N)) + DIRECT(N+1)
c        FNET(N+1) = FUP(N+1) - FDN(N+1)  ! dp
         fnet(n+1) = fdn(n+1) - fup(n+1)  ! dp
 10   CONTINUE
C
C   Calculate surface albedo (SRFALB) and ice transmissivity (TRANS)
      SRFALB = FUP(1)/DIRECT(1)
      TRANS = FDN(NZ1)/DIRECT(1)
C
C   PI times these results. FMT is for diffuse flux only (for comparison
C   with Toon et al., Table 4). FLUXT is total downward flux (for 
C   comparison with Table 5).
      FP0 = PI * FUP(1)
      FMT = PI * (FDN(NZ+1) - DIRECT(NZ+1))
c     FLUXT = - PI * FNET(NZ+1)
      FLUXT =   PI * FNET(NZ+1)
C      PRINT 209, FP0,FMT,FLUXT
 209  FORMAT(1X,'FP0 =',1PE10.3,5X,'FMT =',E10.3,5X,'FLUXT =',E10.3)
C

c     depth distribution of the flux
c
c     WRITE(6,108)
c108  FORMAT(/4X,'DIRECT',5X,'FUP',6X,'FDN',8X,'FNET')
c     WRITE(6,109) (DIRECT(N),FUP(N),FDN(N),FNET(N),N=1,NZ1)
c109  FORMAT(1X,1P3E10.3,1x,E10.3)

      RETURN
      END

c-----------------------------------------------------------------------

      SUBROUTINE TRIDAG(A,B,C,R,U,N)
      PARAMETER (NMAX=200)
      DIMENSION GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
      IF(B(1).EQ.0.)PAUSE
      BET=B(1)
      U(1)=R(1)/BET
      DO 11 J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0.)PAUSE
        U(J)=(R(J)-A(J)*U(J-1))/BET
11    CONTINUE
      DO 12 J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
12    CONTINUE
      RETURN
      END

c+++++
#endif
c+++++
