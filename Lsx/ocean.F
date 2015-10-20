c -- $Id: ocean.F,v 1.2 2001/07/03 21:29:29 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/ocean.F,v $
c -- $Name:  $

c------
c ocean
c------
c
c----------------------------------------------------------------------
c
      SUBROUTINE INIOCN (YTSIN, DEPTHML, DEPTHLAKE, DEPTHICE,
     *                   DYNAMICE, PRESCOUV, SHOWMAPS, DTGCM, 
     *                   NSREST, NSREST_ISOT)
c
c        Does initialization for ocean and seaice.
c
c        YTSIN    = ccm surface temperature field for initialization
c        DEPTHML  = general ocean mixed-layer depth (=0 if not entered)
c        DEPTHLAKE= as DEPTHML but for inland open water(=-1 if not ent)
c        DEPTHICE = max limit on total ice thickness
c        DYNAMICE = .TRUE. for dynamic sea ice, .FALSE. for not.
c        PRESCOUV = .TRUE. for prescribed sfc-ocean currents for dynice
c        SHOWMAPS = .TRUE. to print 2-D maps of prescribed fields
c        DTGCM    = gcm timestep (seconds)
c        NSREST   = -1 for cold start, 0 for initial start from a
c                   restart file, 1 for true restart
 
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commapl.h>
#include <commudxy.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comvegxy.h>
      DIMENSION YTSIN(NLONG,NLATG)
      LOGICAL DYNAMICE, PRESCOUV, SHOWMAPS

#ifdef isotrac
c     delta 18O vs latitude for oceanic surface waters (pacific ocean,
c     indian ocean and atlantic ocean)
      dimension opac_o(19), oind_o(19), oatl_o(19)
c     latitudes:
c       -90,  -80,  -70,  -60,  -50,  -40,  -30,  -20,  -10,    0,
c        10,   20,   30,   40,   50,   60,   70,   80,   90
      data opac_o /
     * 0.00, 0.00, 0.00, 0.00, 0.00, 0.27, 0.56, 0.64, 0.46, 0.41,
     * 0.29, 0.21, 0.19,-0.23,-0.71,-0.64,-0.26, 0.00, 0.00 /
      data oind_o /
     * 0.00, 0.00, 0.00, 0.00, 0.00, 0.27, 0.36, 0.40, 0.19, 0.23,
     * 0.27, 0.68, 0.68, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /
      data oatl_o /
     * 0.00, 0.00, 0.00,-0.28,-0.05, 0.27, 0.82, 0.94, 1.00, 0.90,
     * 0.90, 1.00, 1.10, 0.79,-0.05, 0.19, 0.02, 0.00, 0.00 /

c     for user-provided file with prescribed ocean d18o vs latitude: 
      dimension alat18o(500), d18olat(500), d18olatm(nlat) 
#endif

c
c
c        Physical constants of sea ice
 
      Z0ICE  = .0015  ! .0005
      RHOICE = 0.9E3
      CSICE  = 2.106E3
      CONICE = 2.2
 
c        Coefficient for vertical heat flux from ocean into seaice base
c        (ref: Gordon and Huber,1990,JGR,95,C7,11655-11672).
 
      BASEFAC = 70. ! (100. in early v2)  ! (10. in v1.02)
 
c        Coefficients for lateral heat flux into seaice sides
 
      SIDEFAC = 5.  ! (10. in early v2)
      SIDEZER = 0.1
 
c        HICETOP is "adaptive-grid" thickness of top seaice layer
 
c        HICEMIN is minimum total seaice thickness needed for stability
c        in vertical heat diffusion. (Must be ge NICELAY*HICETOP for
c        VADAPT to work ok). Note that a minimum seaice fraction is
c        not needed due to lateral melt in OCNML, unlike snow and
c        intercepted veg h2o.
 
      HICETOP = .05
      HICEMIN = MAX (.20, NICELAY*HICETOP)
 
c        If not entered in namelist (datgelsx), set open-ocean
c        roughness length for T,s, and max seaice fractional areas
c        in each hemisphere (also depend on ice thickness in OCNML)
 
      IF (Z0OCNH.EQ.0.) Z0OCNH = .0001
 
#ifdef gcm
      IF (FICEMAXS.EQ.0.) FICEMAXS = 0.92  ! 0.92 v2, 0.96  v1
      IF (FICEMAXN.EQ.0.) FICEMAXN = 0.99  ! 0.92 v2, 0.995 v1
#else
      IF (FICEMAXS.EQ.0.) FICEMAXS = 1.0
      IF (FICEMAXN.EQ.0.) FICEMAXN = 1.0
#endif
 
c        Seaice penetrating solar flux parameters (Maykut+Untersteiner)
 
c     penfr now in lsx namelist (6/19/12), default = 0.17.
c     (set to 0 for no penetrating solar radiation): 
c     PENFR = 0.17
      PENDE = 0.666
 
c        Properties of newly formed seaice (frozen in leads).
 
      HICECRE(1) = HICETOP
      DO 60 K=2,NICELAY
        HICECRE(K) = (HICEMIN - HICETOP) / (NICELAY-1)
   60 CONTINUE
 
      ZMCRE = 0.
      DO 62 K=1,NICELAY
        ZMCRE = ZMCRE + RHOICE*HICECRE(K)
   62 CONTINUE
 
c        FBRINMAX is max brine resevoir capacity as a fraction of the 
c        heat needed to melt whole seaice column (Semtner(1976) used 
c        0.3). Even if FBRINMAX is set to 0, whole column must still 
c        reach TMELT before any meltwater discharge occurs (unless 
c        FBRINSTO < 1). Don't set FBRINMAX = 1 or else get 
c        divide-by-zero in SEAICE.
c
c        FBRINSTO is fraction of seaice melt retained as "brine", 
c        as opposed to immediately being discharged (set to 0 for 
c        no brine effects)
c
c        FBRINMAX,STO now in lsx namelist (6/19/12), defaults = 
c        0.25,0.75, or 0,0 if defined RCM (sea ice prescribed, 
c        and brine -> Nan).
c
c     FBRINMAX = 0.25
c     FBRINSTO = 0.75
c
c        Read lsx open-water fraction data input file...except if
c        (entered) depthlake = 0., have no lakes (100% ocean
c        will be set below from lmaskxy=3)
c
c        If rcm, fwaterxy has already been set in inilsx
c
c~~~~~~~~~~
#ifndef rcm
c~~~~~~~~~~
      IF (DEPTHLAKE.EQ.0.) THEN
        CALL ZERO (FWATERXY, NLON*NLAT)
      ELSE
        CALL FINDIU (IUWATER, LOC)
#ifdef gcm
        WRITE(IOTERM2,70) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC)))
   70   FORMAT(/' Reading percent open water data file'/4X,A)
#endif
        CALL OPENF (IUWATER, LOCNAME(LOC), ' ', 0, 0)
        CALL READDATA (IUWATER, FWATERXY, NLON, NLAT, 'FWATER', 2, 1)
        CALL CLOSEF(IUWATER)
      ENDIF
 
c         Change open-water percents to fractions
 
      DO 72 JJ=1,NLAT
        DO 74 JI=1,NLON
#ifdef gcm
c         change from percent (data input file units) to fraction:
          FWATERXY(JI,JJ) = .01*FWATERXY(JI,JJ)
#endif
c         Ensure between -1 and 1, and for safety in later if-tests
c         (...lt.1., etc):
          IF (FWATERXY(JI,JJ).GT. 0.999) FWATERXY(JI,JJ) =  1.
          IF (FWATERXY(JI,JJ).LT.-0.999) FWATERXY(JI,JJ) = -1.
  74    CONTINUE
  72  CONTINUE
c~~~~~
#endif
c~~~~~
c
c        Interpolate initialization temperature field to lsx grid
c        (for setting cold start ocean temperatures and sea ice)
c
      CALL GTOL (YTSIN, WORKL, 0)
 
c        Over land, set ocean and seaice variables to null values
c        (nb: HOCNXY and TOCNXY must be 0 for BUDALL).
c        Over ocean, initialize variables if cold-start (NSREST=-1) or
c        if user has modified the land-ocean mask to add new ocean
c        points (flagged by previous HOCNXY=0).
c
c        First, set HDEEP and HSHALLOW (comocn, used in dopresst).
c        If namelist param depthml entered, set deep depths to
c        max (depthml,0.1 m), else to 50 m. Similarly set shallow
c        depths to 5 m, or max (depthlake, 0.1 m) if depthlake entered
c        (if not entered, default depthlake = -1.).
 
      HDEEP = 50.
      IF (DEPTHML.GT.0.) HDEEP = MAX (DEPTHML, .1)
      HSHALLOW = 5.
      IF (DEPTHLAKE.GT.0.) HSHALLOW = MAX (DEPTHLAKE, .1)
 
c=======================
      DO 100 JJ=1,NLAT
        DO 102 JI=1,NLON
c=======================
 
          HOCNPREV = HOCNXY(JI,JJ)
 
c            Set final open-water fraction fwater and deep/shallow depth
c            hocnxy, distinguished by data-file fwater being +ve (deep)
c            or -ve (shallow).
 
c---------
#if defined (PATRICK)
          IF (LMASKXY(JI,JJ).EQ.2) THEN

c           ice sheet (always 100% of grid cell):
            FWATERXY(JI,JJ) = 0.
            HOCNXY(JI,JJ) = 0.

          ELSE IF (LMASKXY(JI,JJ).EQ.3 .OR. FWATERXY(JI,JJ).GT.0.) THEN
#else
          IF (LMASKXY(JI,JJ).EQ.3 .OR. FWATERXY(JI,JJ).GT.0.) THEN
#endif
 
c           Deep:
            IF (LMASKXY(JI,JJ).EQ.3) FWATERXY(JI,JJ) = 1.
#ifdef hydmodel
c           set standing water depths for hydrol categories
            icat = mod(ji-1,ncateghyd) + 1 
            hocnxy(ji,jj) = hwaterhyd(icat)
#else
            HOCNXY(JI,JJ) = HDEEP
#endif
 
          ELSE IF (FWATERXY(JI,JJ).LT.0.) THEN
 
c           Shallow:
            FWATERXY(JI,JJ) = -FWATERXY(JI,JJ)
            HOCNXY(JI,JJ) = HSHALLOW
 
          ELSE IF (FWATERXY(JI,JJ).EQ.0.) THEN
 
c           None:
            HOCNXY(JI,JJ) = 0.
 
          ENDIF
c---------
 
c            Set ocean temperatures and sea ice if cold start(NSREST=-1)
c            or if user has modified the land-ocean map to add new water
c            points (flagged by previous HOCNXY=0).
 
c----------------------------------------
          IF (FWATERXY(JI,JJ).EQ.0.) THEN
c----------------------------------------
 
            TOCNXY(JI,JJ) = TMELT
            SOCNXY(JI,JJ) = 0.
            TFREEZXY(JI,JJ) = TMELT
            FICEXY(JI,JJ) = 0.
 
c-------------
          ELSE
c-------------
 
            IF (NSREST.EQ.-1. OR. HOCNPREV.EQ.0.) THEN
 
              SOCNXY(JI,JJ) = .035
              TFREEZXY(JI,JJ) = TMELT - .054 * 1000.*SOCNXY(JI,JJ)
              TOCNXY(JI,JJ) = MAX (TFREEZXY(JI,JJ), WORKL(JI,JJ))
#ifdef gcm
              IF (LMASKXY(JI,JJ).EQ.3) THEN
                IF (WORKL(JI,JJ).LT.TFREEZXY(JI,JJ)+.001) THEN
                  TOCNXY(JI,JJ) = TFREEZXY(JI,JJ)
                  FICEXY(JI,JJ) = 1.
                ELSE
                  TOCNXY(JI,JJ) = WORKL(JI,JJ)
                  FICEXY(JI,JJ) = 0.
                ENDIF
              ELSE
                TOCNXY(JI,JJ) = MAX (TMELT, WORKL(JI,JJ))
c~~~
#ifdef rcm
                TOCNXY(JI,JJ) = MAX (TMELT+1., WORKL(JI,JJ)) !1C, no ice
#else
                TOCNXY(JI,JJ) = MAX (TMELT, WORKL(JI,JJ))
                IF (WORKL(JI,JJ).LT.TMELT+.001) FICEXY(JI,JJ) = 1.
#endif
c~~~
              ENDIF
#else
              FICEXY(JI,JJ) = xfice
#endif
 
              IF (FICEXY(JI,JJ).GT.0.) THEN
c               For local machine, allow initial tocnxy>tfreez below ice
#ifdef gcm
                TOCNXY(JI,JJ) = TFREEZXY(JI,JJ)
                H = MAX (2.0,   HICEMIN)
#else
                H = MAX (xhice, HICEMIN)
#endif
                DO 130 K=1,NICELAY
                  IF (K.EQ.1) HICEXY(JI,JJ,K) = HICETOP
                  IF (K.GT.1) HICEXY(JI,JJ,K) = (H-HICETOP)/(NICELAY-1)
                  TICEXY(JI,JJ,K) = MIN (TMELT, WORKL(JI,JJ))
  130           CONTINUE
 
                BRINEXY(JI,JJ) = 0.
                TG_ICEXY(JI,JJ) = TICEXY(JI,JJ,1)
              ELSE
                FICEXY(JI,JJ) = 0.
              ENDIF
 
            ELSE
c              tfreezxy not on restart file
               TFREEZXY(JI,JJ) = TMELT - .054 * 1000.*SOCNXY(JI,JJ)
            ENDIF
 
c              If namelist param depthice entered, apply seaice
c              thickness limit.
 
            IF (DEPTHICE.GT.0. .AND. FICEXY(JI,JJ).GT.0.) THEN
               H = MAX (DEPTHICE, HICEMIN)
               DO 134 K=2,NICELAY
                 HICEXY(JI,JJ,K) = MIN ( HICEXY(JI,JJ,K),
     *                                   (H-HICETOP)/(NICELAY-1) )
  134          CONTINUE
            ENDIF
c--------------
          ENDIF
c--------------
 
c            Set sea ice to null if ice fraction = 0
 
          IF (FICEXY(JI,JJ).EQ.0.) THEN
            DO 140 K=1,NICELAY
              HICEXY(JI,JJ,K) = 0.
              TICEXY(JI,JJ,K) = TMELT
  140       CONTINUE
            FICEXY(JI,JJ) = 0.
            BRINEXY(JI,JJ) = 0.
            TG_ICEXY(JI,JJ) = TMELT
          ENDIF
 
c===============
  102   CONTINUE
  100 CONTINUE
c===============
c
c       Print open-water map if requested (after possible changes above)
c
      IF (SHOWMAPS) THEN
        DO JJ=1,NLAT
          DO JI=1,NLON
            IWORKL(JI,JJ) = NINT(100.*FWATERXY(JI,JJ))
          ENDDO
        ENDDO
        CALL MAPOUT (IOTERM2, IWORKL, NLON, NLAT, ALON, ALAT, 0,
     *               'Open-water percent (-shallow, +deep)')
      ENDIF

c!!!!!!!!!!!!!
#ifdef isotrac
c!!!!!!!!!!!!!

c................................
      if (oceano18.eq.-999.) then
c................................

c          Read indiv-ocean map

        call findiu (iuocean_i, loc)
#ifdef gcm
        write(ioterm2,148) locname(loc)(1:lenchr(locname(loc)))
  148   format(/' Reading isotope surface-type data file'/4x,a)
#endif
        call openf (iuocean_i, locname(loc), ' ', 0, 0)
        call readdata (iuocean_i, workm, nlon, nlat, 'OCEAN_I', 1, 1)
        call closef(iuocean_i)
        do jj=1,nlat
          do ji=1,nlon
            mocean_i(ji,jj) = nint(workm(ji,jj))
          enddo
        enddo

        nerr = 0
        do jj=1,nlat
          do ji=1,nlon

            if (fwaterxy(ji,jj).gt.0.) then
              if ( mocean_i(ji,jj).lt.1.or.mocean_i(ji,jj).gt.10 ) then
                write (ioterm2,*)
     *              'OCEAN_I error: i,j,val,lmask,fwater: ',
     *              ji,jj,mocean_i(ji,jj),lmaskxy(ji,jj),fwaterxy(ji,jj)
                nerr = nerr+1
              endif 
              if ( (lmaskxy(ji,jj).ne.1 .or. fwaterxy(ji,jj).eq.1.)
     *             .and. (mocean_i(ji,jj).eq.8.or.mocean_i(ji,jj).eq.9)
     *           ) then
                write (ioterm2,*)
     *              'OCEAN_I warning: i,j,val,lmask,fwater: ',
     *              ji,jj,mocean_i(ji,jj),lmaskxy(ji,jj),fwaterxy(ji,jj)
c               nerr = nerr+1
              endif 
            else 
              mocean_i(ji,jj) = 0
            endif
          enddo
        enddo
        if (nerr.ge.1) call endrun (-1)

c.....................................
      else if (oceano18.eq.1000.) then
c.....................................

c          Read user file with prescribed (zonal) ocean d18O vs latitude

#ifdef gcm
        write(ioterm2,150) 
  150   format(/' Reading ocean d18O vs latitude file OCEAN_I_LAT')
#endif
        iu = iuocean_i
        open (unit=iu, file='OCEAN_I_LAT',
     *        form='formatted', status='old', err=152)
        goto 154
  152   write(ioterm2,*) ' *** Error opening file OCEAN_I_LAT'
        call endrun (-1)
  154   continue

        np18o = 0
        do jz=1,500 
          read (iu,*,end=156) alat18o(jz), d18olat(jz)
          np18o = np18o + 1
        enddo
  156   continue
        close (iu)

c       linear interp to model lat grid (d18olat to d18olatm):
        do jj=1,nlat
          zlat = alat(jj)*180./pi
          if (zlat.le.alat18o(1)) then
            d18olatm(jj) = d18olat(1)
          else if (zlat.ge.alat18o(np18o)) then
            d18olatm(jj) = d18olat(np18o)
          else
            do jz = 2,np18o
              if (zlat.ge.alat18o(jz-1).and.zlat.le.alat18o(jz)) then
                zwei =(zlat-alat18o(jz-1)) / (alat18o(jz)-alat18o(jz-1))
                d18olatm(jj) =(1.-zwei)*d18olat(jz-1) + zwei*d18olat(jz)
              endif
            enddo
          endif
          write (6,*) 'jj=',jj,'  d18olatm=',d18olatm(jj)
        enddo

c.........
      else
c.........

        do jj=1,nlat
          do ji=1,nlon
            mocean_i(ji,jj) = 0
          enddo
        enddo
c..........
      endif
c..........

c###################
      DO JJ=1,NLAT
        DO JI=1,NLON
c###################

c            Set isotopic compositions of ocean surface water
c            (and ice sheets)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if ( nsrest_isot.eq.-1 .and. fwaterxy(ji,jj).gt.0.) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
            xlj = (jj-.5)*180./nlat-90.
            k   = int((xlj+90.)/10.)+1
            xlk = 10.*(k-10.)

            r_o = -999.
            r_h = -999.

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if (oceano18.eq.-999.) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c              From data input file OCEAN_I

              if (mocean_i(ji,jj).eq.1) then
c               Pacific:
                r_o = opac_o(k)+(xlj-xlk)*(opac_o(k+1)-opac_o(k))/10.
                r_h = 7.0*r_o
              else if (mocean_i(ji,jj).eq.2 .or. mocean_i(ji,jj).eq.10)
     *          then
c               Atlantic or next to ice sheet:
                r_o = oatl_o(k)+(xlj-xlk)*(oatl_o(k+1)-oatl_o(k))/10.
                r_h = 7.3*r_o
              else if (mocean_i(ji,jj).eq.3) then
c               Indian:
                r_o = oind_o(k)+(xlj-xlk)*(oind_o(k+1)-oind_o(k))/10.
                r_h = 7.0*r_o
              else if (mocean_i(ji,jj).eq.4) then
c               Hudson Bay:
                r_o = -4.8
                r_h = -35.
              else if (mocean_i(ji,jj).eq.5) then
c               Mediterranean:
                r_o = 1.6
                r_h = 8.0
              else if (mocean_i(ji,jj).eq.6) then
c               Black Sea:
                r_o = -2.9
                r_h = -24.
              else if (mocean_i(ji,jj).eq.7) then
c               Baltic Sea:
                r_o = -6.1
                r_h = -54.
              endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            else if (oceano18.eq.1000.) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c                User-provided latitudinal profile file (read above)

              r_o = d18olatm(jj)
              r_h = 7.0*r_o

c>>>>>>>>>>>>>>>
            else
c>>>>>>>>>>>>>>>

c                Uniform d18o value

              r_o = oceano18
              r_h = 7.0*r_o

c>>>>>>>>>>>>>>>>
            endif
c>>>>>>>>>>>>>>>>
 
c           convert delta to isotopic ratio
            if (r_o.ne.-999.) then
c             write (6,*)'ocean d180: ji,jj=',ji,jj,' r_o=',r_o
              ro_o(ji,jj) = (r_o/1.e3 + 1.)*smoj_o
              ro_h(ji,jj) = (r_h/1.e3 + 1.)*smoj_h
            endif
c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~
c############
        ENDDO
      ENDDO
c############
c!!!!!
#endif
c!!!!!
 
c        If slab ocean and dynamic ice and PRESCOUV requested, read in
c        prescribed ocean surface currents UOCNXY,VOCNXY from data input
c        file NOCEAN for dynamic seaice model. (Dynamic topography
c        gradients are computed from UOCNXY,VOCNXY by the dynamic seaice
c        model assuming geostrophy.) Do both for initial and restarts.
 
c        If slab ocean, dynamic ice, and PRESCOUV = .F., zero currents.
c        If ogcm and dyanmic ice, currents will be set by ogcm.
c
c        If read data file, allow for:
c        (i)  file latitudes are staggered, omit south pole.
c        (ii) file values are mm/s, convert to m/s.
c
      CALL ZERO (UOCNXY, NLON*NLAT)
      CALL ZERO (VOCNXY, NLON*NLAT)
 
c---------
#ifdef gcm
c---------
 
c-------------------------------------------------------------
      IF (DYNAMICE .AND. (OCEANTYPE.EQ.1) .AND. PRESCOUV) THEN
c-------------------------------------------------------------
 
        CALL FINDIU (NOCEAN, LOC)
        WRITE(IOTERM2,200) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC)))
  200   FORMAT(/' Reading surface-ocean-current data file'/4X,A)
 
        CALL OPENF (NOCEAN, LOCNAME(LOC), ' ', 0, 0)
c
        CALL READDATA (NOCEAN, WORKL, NLON, NLAT, 'OCNU', 2, 1)
        DO 202 JJ=1,NLAT
          JJM = MAX (JJ-1,1)
          DO 204 JI=1,NLON
            UOCNXY(JI,JJ) = 0.5*(WORKL(JI,JJM)+WORKL(JI,JJ))*.001
  204     CONTINUE
  202   CONTINUE
c
        CALL READDATA (NOCEAN, WORKL, NLON, NLAT, 'OCNV', 2, 2)
        DO 212 JJ=1,NLAT
          JJM = MAX (JJ-1,1)
          DO 214 JI=1,NLON
            VOCNXY(JI,JJ) = 0.5*(WORKL(JI,JJM)+WORKL(JI,JJ))*.001
  214     CONTINUE
  212   CONTINUE
c
        CALL CLOSEF(NOCEAN)
 
c          Ensure ocean currents are zero in lsx land interiors
c
        DO 310 JJ=1,NLAT
          DO 312 JI=1,NLON
c           IF (LMASKXY(JI,JJ) .LE.2) THEN
            if (fwaterxy(ji,jj) .lt. 1.) THEN
              UOCNXY(JI,JJ) = 0.
              VOCNXY(JI,JJ) = 0.
            ENDIF
  312     CONTINUE
  310   CONTINUE
c
c        Print ocean current maps if requested
c
        IF (SHOWMAPS) THEN
c
          DO 410 JJ=1,NLAT
            DO 412 JI=1,NLON
c             IF (LMASKXY(JI,JJ).LE.2) THEN
              if (fwaterxy(ji,jj) .lt. 1.) THEN
                IWORKL(JI,JJ) = 999999
              ELSE
                IWORKL(JI,JJ) = NINT(UOCNXY(JI,JJ)*1000.)
              ENDIF
  412       CONTINUE
  410     CONTINUE
          CALL MAPOUT (IOTERM2, IWORKL, NLON, NLAT, ALON, ALAT,999999,
     *                 'Eastward surface ocean current (mm/sec)')
c
          DO 420 JJ=1,NLAT
            DO 422 JI=1,NLON
c             IF (LMASKXY(JI,JJ).LE.2) THEN
              if (fwaterxy(ji,jj) .lt. 1.) THEN
                IWORKL(JI,JJ) = 999999
              ELSE
                IWORKL(JI,JJ) = NINT(VOCNXY(JI,JJ)*1000.)
              ENDIF
  422       CONTINUE
  420     CONTINUE
          CALL MAPOUT (IOTERM2, IWORKL, NLON, NLAT, ALON, ALAT,999999,
     *                 'Northward surface ocean current (mm/sec)')
c
        ENDIF
 
c-----------------------------------------------------------------------
      ELSE IF (DYNAMICE .AND. (OCEANTYPE.EQ.1) .AND. .NOT.PRESCOUV) THEN
c-----------------------------------------------------------------------
 
        WRITE (IOTERM2,500)
  500   FORMAT(/' *** Ocean currents for dynamic seaice set to zero !')
 
c----------
      ENDIF
c----------
 
c-----
#endif
c-----


      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SETICE
c
c        Sets diagnostic seaice quantities (in commus) needed by
c        seaice and lsx.
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsat.h>
 
c        Set latent heat for sea-ice surface, used in ocnflx and by
c        lsx for heatg. Use prognostic top-layer temperature tice, not
c        the diagnostic skin temperature tg_ice, since budgets assume
c        vapor is exchanged with top-layer store at tice.
 
      DO 400 I=1,NPOI
        HVASUG(I) = HSUBF(TICE(I,1),TA(I))
  400 CONTINUE
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE OCNFLX (DTIM)
c
c        Calculates surface fluxes for ocean, seaice and snow on seaice.
c        The solar, ir, sensible and latent calculations below parallel
c        those for land surfaces in radlsx and turlsx, but are simpler.
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION CDLARGE(NVEC), ALOGGH(NVEC),    QOCN(NVEC),
     *          TAURAN(NVEC),  URAN(NVEC),
     *          CDH(NVEC),     CDE(NVEC),       CDM(NVEC),
     *          TAV(NVEC),     QAV(NVEC),
     *          TGOLD(NVEC),   TIOLD(NVEC),
     *          TGPRE(NVEC),   TIPRE(NVEC),
     *          QG(NVEC),      QI(NVEC),
     *          DQG(NVEC),     DQI(NVEC),
     *          COG(NVEC),     COI(NVEC),
     *          WG(NVEC),      WI(NVEC),
     *          EMG(NVEC),     EMI(NVEC),       EMO(NVEC),
     *          ZIRG(NVEC),    ZIRI(NVEC),      richu2(nvec)
      PARAMETER (NITER=3)
#include <comsat.h>
#include <comrich.h>
#ifdef isotrac
c     dv_i: coefficient of diffusion of water vapour in air
c     (T in K, p in Pa):
      dv_i(T,P)= 2.17e-5*(101325./P)*(T/273.15)**1.88
#endif

c        Zero net ocean and sea-ice fluxes since they may be
c        incremented by various processes in OCNFLX, SEAICE, ICEPHASE,
c        and OCNML.
c
      CALL ZERO (FTOCN, NVEC)
      CALL ZERO (FWOCN, NVEC)
      CALL ZERO (FUOCN, NVEC)
      CALL ZERO (FVOCN, NVEC)
c     for mom2,4:
      CALL ZERO (FHOCN,   NVEC)
      CALL ZERO (FQOCN,   NVEC)
      CALL ZERO (FSOLOCN, NVEC)
      CALL ZERO (FIROCN,  NVEC)
      CALL ZERO (FPRLOCN, NVEC)
      CALL ZERO (FPRFOCN, NVEC)
      CALL ZERO (FSALOCN, NVEC)
      CALL ZERO (PRESOCN, NVEC)

      CALL ZERO (FUICE, NVEC)
      CALL ZERO (FVICE, NVEC)
c
c--------------------
c        Solar fluxes
c--------------------
c
      CALL ZERO (SOLG, NVEC)
      CALL ZERO (SOLI, NVEC)
      CALL ZERO (SOLO, NVEC)
      CALL SOLSET
      CALL SOLSUR (3)
      IF (NSOL.GT.0) THEN
        DO 100 IB=1,NBAND
!DIR$ IVDEP
          DO 102 J=1,NSOL
            I = INDSOL(J)
            SOLG(I) = SOLG(I) + SOLAD(I,IB)*(1.-ALBICD(I,IB))
     *                        + SOLAI(I,IB)*(1.-ALBICI(I,IB))
            SOLI(I) = SOLI(I) + SOLAD(I,IB)*(1.-ALBSND(I,IB))
     *                        + SOLAI(I,IB)*(1.-ALBSNI(I,IB))
            SOLO(I) = SOLO(I) + SOLAD(I,IB)*(1.-ALBOCD(I,IB))
     *                        + SOLAI(I,IB)*(1.-ALBOCI(I,IB))
  102     CONTINUE
  100   CONTINUE

#if defined (TWOSTRICE)
c       scale visible two-stream penetrative fluxes by incident fluxes: 
        do j=1,nsol
          i = indsol(j)
          do k=1,nicelay+1
            fluxiced(i,k) = fluxiced(i,k)*solad(i,1)
            fluxicei(i,k) = fluxicei(i,k)*solai(i,1)
          enddo
        enddo
#endif
      ENDIF
c
c-----------------
c        IR fluxes
c-----------------
c
      DO 200 I=1,NPOI
        EMG(I) = EMISICE
        EMI(I) = EMISNOW
        EMO(I) = EMIOCN
 
        ZG = EMG(I) * STEF * (TG_ICE(I)**4) + (1.-EMG(I))*FIRA(I)
        ZI = EMI(I) * STEF * (TI(I)**4)     + (1.-EMI(I))*FIRA(I)
        ZO = EMO(I) * STEF * (TOCN(I)**4)   + (1.-EMO(I))*FIRA(I)
 
        FIRB(I) = (1.-FICE(I))*ZO
     *               +FICE(I) *( (1.-FI(I))*ZG + FI(I)*ZI )
        FIRG(I) = FIRA(I) - ZG
        FIRI(I) = FIRA(I) - ZI
        FIRO(I) = FIRA(I) - ZO
  200 CONTINUE
 
c------------------------------------------
c        Aerodynamic fluxes over open ocean
c------------------------------------------
 
c        Set atmos conditions at ZA and other temporary aero vars
c        (like TURINI)
 
      TFAC = 1./(SIGA**CAPPA)
      DO 250 I=1,NPOI
        PA = PSURF(I) * SIGA
        THICA(I) = PSURF(I)*(1.-SIGB) / GRAV
        RHOA(I) = PA / ( RAIR * TA(I) * (1. + (RVAP/RAIR-1.)*QA(I)) )
        CP(I) = CAIR * (1. + (CVAP/CAIR-1.)*QA(I))
        ZA(I) = (PSURF(I) - PA) / (RHOA(I) * GRAV)
        ALOGA(I) = LOG (ZA(I))
        UA(I) = MAX ( SQRT(UX(I)**2 + UY(I)**2), 1. )
  250 CONTINUE
 
c        Set ocean roughness length for momentum (alogg). Do via
c        10m-drag coeff of Large and Pond (ref: Trenberth,Olson,Large,
c        1982,NCAR TN-338, pg.8).
 
c        Also set ocean roughness length for heat/vap (aloggh). This is
c        different from momentum, indep. of wind speed (ref, Miller
c        et al,J.Climate,5,418-434). (Default for namelist variable
c        z0ocnh is .0001 m, set in iniocn).
 
      ZLOGH  = LOG (Z0OCNH)
      ZLOG10 = LOG (10.)
 
      DO 300 I=1,NPOI
        CDLARGE(I) = CVMGT (1.14, .62+1.56/UA(I), UA(I).GT.3.)
  300 CONTINUE
      DO 310 I=1,NPOI
c       CDLARGE(I) = .62 + 1.56/UA(I)
c       IF (UA(I).GT.3.) CDLARGE(I) = 1.14
c       IF (UA(I).GT.10.) CDLARGE(I) = .49+.065*UA(I)
        CDLARGE(I) = CVMGT (.49+.065*UA(I), CDLARGE(I), UA(I).GT.10.)
c       For u:
        ALOGG(I) = ZLOG10 - VONK / SQRT(.001*CDLARGE(I))
c       For T,q:
        ALOGGH(I) = ZLOGH
 
c       Set ocean surface specific humidity
        ES = ESAT (TOCN(I))
        QOCN(I) = QSAT (ES, PSURF(I))
  310 CONTINUE
 
c        Calculate stratification factors for open ocean surface to za
 
      CALL FSTRAT (TOCN, TA, TFAC, QOCN, QA, VZERO, ZA,
     *             ALOGG, ALOGGH, ALOGA, UA, RICHU, STRAMU, STRAHU, 1)
 
c       Assume triangular sub-grid windspeed distribution between
c       0 (actually 0.1 m/s) and 2*UA, and calculate mean neutral
c       wind stress TAURAN for that wind distribution. (Triangular
c       distribution is similar to Weibull distrib: Pavia+O'Brien,
c       J.Clim.Appl.Meteor.,25,1324-1332).
c
c     CALL ZERO (TAURAN, NPOI)
c     WEIRAN = 0.
c     NRAN = 11
c     DO 320 IRAN = 1,NRAN
c       ZFAC = 2.*(IRAN-1.)/(NRAN-1.)
c       ZWEI = 1. - ABS(1.-ZFAC)
c       DO 322 I=1,NPOI
c         URAN(I) = MAX (0.1, UA(I)*ZFAC)
c         IF (URAN(I).LE.3.) THEN
c           CDLARGE(I) = .62+1.56/URAN(I)
c         ELSE IF (URAN(I).LE.10.) THEN
c           CDLARGE(I) = 1.14
c         ELSE
c           CDLARGE(I) = .49+.065*URAN(I)
c         ENDIF
c         CDRAN = URAN(I) * (VONK/(ALOGA(I)-ALOGRAN))**2
c         TAURAN(I) =  TAURAN(I) + RHOA(I)*CDRAN*URAN(I)*ZWEI
c 322   CONTINUE
c       WEIRAN = WEIRAN + ZWEI
c 320 CONTINUE
c     DO 324 I=1,NPOI
c       TAURAN(I) = TAURAN(I)/WEIRAN
c 324 CONTINUE
c
c        Sensible and latent fluxes, and wind stresses over open ocean
c
      ZLOG2 = LOG(2.)
      DO 330 I=1,NPOI
        CDH(I) = CDHFAC*UA(I)*(VONK/(ALOGA(I)-ALOGGH(I)))**2 * STRAHU(I)
        CDE(I) = CDEFAC* (CDH(I)/CDHFAC)
        CDM(I) = CDMFAC*UA(I)*(VONK/(ALOGA(I)-ALOGG (I)))**2 * STRAMU(I)
 
c       Either: use triangular wind-speed distrib computed above:
c       CDM(I) = CDMFAC* (TAURAN(I)/(RHOA(I)*UA(I))) * STRAMU(I)
c       Or: crudely account for random wind speeds:
        CDM(I) = CDM(I) * (1. + .6 * MIN(UA(I)/10.,1.))
 
c       Crudely account for enhanced evap of ocean spray at high u:
c       (Bortkovskii, in Smith,BLM,1989,47,277-293,Fig.1,Eq.13)
c       cde(i) = cde(i) + 1.5 * sqrt (max(ua(i)-15.,0.)/5. )

c       arb increase, high ua:
c       cde(i) = cde(i) * (1. + max (0., min (2., (ua(i)-8.)/4.)))
c       arb increase, high ts-ta:
c       ztas = tocn(i) - ta(i)
c       cde(i) = cde(i) * (1. + max (0., min (2., (ztas-5.)/2.)))
 
        CDH(I) = MIN (CDMAXA, CDH(I) / (1. + CDH(I)/CDMAXB))
        CDE(I) = MIN (CDMAXA, CDE(I) / (1. + CDE(I)/CDMAXB))
        CDM(I) = MIN (CDMAXA, CDM(I) / (1. + CDM(I)/CDMAXB))
c
c       FSENO(I) = RHOA(I) * CDH(I) * CP(I) * (TOCN(I) - TA(I)*TFAC)
c       FVAPO(I) = RHOA(I) * CDE(I) * (QOCN(I)-QA(I))
c       TAUXO(I) = RHOA(I) * CDM(I) * UX(I)
c       TAUYO(I) = RHOA(I) * CDM(I) * UY(I)
 
c       Partially time-implicit in ta,qa,ux,uy:
        ZMT = THICA(I)/DTIM
        ZFACH = ZMT * (1.-EXP(-RHOA(I)*CDH(I)*TFAC/ZMT))
        ZFACE = ZMT * (1.-EXP(-RHOA(I)*CDE(I)/ZMT))
        ZFACM = ZMT * (1.-EXP(-RHOA(I)*CDM(I)/ZMT))
        FSENO(I) = ZFACH * CP(I) * (TOCN(I)/TFAC - TA(I))
        FVAPO(I) = ZFACE * (QOCN(I)-QA(I))
#ifdef isotrac
c       Save the values of ZFACE and QOCN:
        zf_i(I) = ZFACE
        qo_i(I) = QOCN(I)
#endif
        TAUXO(I) = ZFACM * UX(I)
        TAUYO(I) = ZFACM * UY(I)
c
c          Increment ocean fluxes (open-ocean contribution), combining
c          heat and h2o flxs with under-seaice contribs already accumed.
c          Nb: wind stresses F[U,V]OCN are kept separately from the
c          seaice-sfc wind stresses F[U,V]ICE; (under-seaice stresses
c          will be calced and combined with F[U,V]OCN by the dynamic
c          seaice module).
c
        ZH = SOLO(I) + FIRO(I) - FSENO(I)
     *       - HVAPF(TOCN(I),TA(I))*FVAPO(I)
     *       - RAINA(I)* CH2O*(TOCN(I)-MAX(TA(I),TMELT))
     *       - SNOWA(I)*(HFUS + CSICE*(TMELT-MIN(TA(I),TMELT))
     *                        + CH2O*(TOCN(I)-TMELT))
        FTOCN(I) = FTOCN(I) + ZH * (1.-FICE(I))

        ZW = RAINA(I) + SNOWA(I) - FVAPO(I)
        FWOCN(I) = FWOCN(I) + ZW * (1.-FICE(I))

c       for mom2,4:
        FHOCN(I)   = FHOCN(I)   + FSENO(I)  *(1.-FICE(I))
        FQOCN(I)   = FQOCN(I)   + FVAPO(I)  *(1.-FICE(I))
        FSOLOCN(I) = FSOLOCN(I) + SOLO(I)   *(1.-FICE(I))
        FIROCN(I)  = FIROCN(I)  + FIRO(I)   *(1.-FICE(I))
        FPRLOCN(I) = FPRLOCN(I) + RAINA(I)  *(1.-FICE(I))
        FPRFOCN(I) = FPRFOCN(I) + SNOWA(I)  *(1.-FICE(I))
c       fsalocn only for sea-ice freeze/melt:
c       FSALOCN(I) = FSALOCN(I) + 0.
        PRESOCN(I) = PRESOCN(I) + PSURF(I)  *(1.-FICE(I))

        FUOCN(I) = TAUXO(I)
        FVOCN(I) = TAUYO(I)
c
  330 CONTINUE

c        Isotopic fluxes and fractionation for evap over open ocean
 
c>>>>>>>>>>>>>
#ifdef isotrac
c>>>>>>>>>>>>>
c===============================================
c       ietest is a namelist variable, in comgrd
      zfac = 1.
      if (ietest.le.-2) zfac = -ietest/100.    ! (factor*ua, isot below)
c===============================================
 
c~~~~~~~~~~~~~~~~~~
      do i = 1,npoi
c~~~~~~~~~~~~~~~~~~
 
c         Compute the friction velocity ufr_i;
c                 the kinematic molecular air viscosity rvisc;
c                 the smooth/rough regime transition fric. vel. ufr_t;
c                 the ratio between turb. and diff. resistances rtorm;
c         and finaly the value of the kinetic fractionation exponent xo
c         (Merlivat and Coantic, 75; Merlivat and Jouzel, 79)
 
c=============================
        if (ietest.eq.-1) then
c=============================
 
c         Use gcm's coeffs cdm, cde to account for stability,
c         and to use Brutsaert's (1975) Eq. A6 (non-neutral).
          ufr_i = sqrt(cdm(i)*ua(i))
 
          rvisc = 1.3067e-5+9.8e-8*(ta(i)-273.16)
          ufr_t = (81.1*rvisc*grav)**(1./3.)
 
c smooth regime
          if (ufr_i.le.ufr_t) then
            x     = 2./3.
            rtorm = (sqrt(ua(i)/cde(i)) - 13.5)
     *            / (13.6*(rvisc/dv_i(ta(i),psurf(I)))**x)
c rough regime
          else
            x     = 0.5
            rtorm = (sqrt(ua(i)/cde(i)) - 5.)
     *            / ( 7.3*(ufr_i*ufr_i*ufr_i/81.1/grav/rvisc)**0.25
     *                *(rvisc/dv_i(ta(i),psurf(I)))**x )
          endif
 
          xo_o(i) = log((rtorm+dr_o**x)/(rtorm+1))/log(dr_o)
          xo_h(i) = log((rtorm+dr_h**x)/(rtorm+1))/log(dr_h)
 
c===================================================
        else if (ietest.le.-2 .or. ietest.eq.0) then
c===================================================
 
c         "standard" renaud if ietest = 0, or zfac*ua if ietest <= -2
 
c         friction velocity from ua(i) by 3 iters
          ufr_i = 1.e-4
          do i_i=1,3
c           ufr_i = ufr_i-(0.4     *ua(i)/sqrt(795.6*za(i))
            ufr_i = ufr_i-(0.4*zfac*ua(i)/sqrt(795.6*za(i))
     *            + 2.*ufr_i*log(ufr_i))/2./(log(ufr_i)+1.)
          enddo
          ufr_i = sqrt(795.6*za(i))*ufr_i
 
          rvisc = 1.3067e-5+9.8e-8*(ta(i)-273.16)
          ufr_t = (81.1*rvisc*grav)**(1./3.)
 
c smooth regime
          if (ufr_i.le.ufr_t) then
            x     = 2./3.
            rtorm = log(ufr_i*za(I)/30./rvisc)/
     *                  (0.4*13.6*(rvisc/dv_i(ta(i),psurf(I)))**x)
c rough regime
          else
            x     = 0.5
            rtorm = (log(81.1*grav*za(I)/ufr_i/ufr_i)/0.4-5.)/
     *                  (7.3*(ufr_i*ufr_i*ufr_i/81.1/grav/rvisc)**0.25*
     *                  (rvisc/dv_i(ta(i),psurf(I)))**x)
          endif
 
          xo_o(I) = log((rtorm+dr_o**x)/(rtorm+1))/log(dr_o)
          xo_h(I) = log((rtorm+dr_h**x)/(rtorm+1))/log(dr_h)
 
c=================================
        else if (ietest.eq.1) then
c=================================
 
c etest (Es, 100% smooth):
 
c       xo_o(i) = log(1.0063)   /log(dr_o)   ! 10 m (paper p.34)
c       xo_h(i) = log(1.0050)   /log(dr_h)   ! 10 m (paper p.34)
        xo_o(i) = log(1.0054821)/log(dr_o)   ! 60 m (xls table)
        xo_h(i) = log(1.0048356)/log(dr_h)   ! 60 m (xls table)
c       xo_o(i) = 0.1946                     ! Renaud email
c       xo_h(i) = 0.1945                     ! Renaud email
 
c=================================
        else if (ietest.eq.2) then
c=================================
 
c etest (Er, 100% rough):
 
c       xo_o(i) = log(1.0027)   /log(dr_o)   ! 10 m (paper p.34)
c       xo_h(i) = log(1.0024)   /log(dr_h)   ! 10 m (paper p.34)
        xo_o(i) = log(1.0023763)/log(dr_o)   ! 60 m (xls table)
        xo_h(i) = log(1.0020967)/log(dr_h)   ! 60 m (xls table)
c       xo_o(i) = 0.0845                  ! Renaud email
c       xo_h(i) = 0.0844                  ! Renaud email
 
c=================================
        else if (ietest.eq.3) then
c=================================
 
c etest: Jouzel et al. (1987) GISS, Renaud email 7/18/01 10->60 m:
 
        if (ua(i).lt.8.) then
          zko = 5.25 / 1000.
        else
          zko = (0.210*ua(i) + 0.80) / 1000.
        endif
        zkh = 0.88*zko
 
        zalf_o = 1. / (1.-zko)
        zalf_h = 1. / (1.-zkh)
 
        xo_o(i) = log(zalf_o)/log(dr_o)
        xo_h(i) = log(zalf_h)/log(dr_h)
 
c===========
        else
c===========
 
        write (ioterm2,*) '*** Error: unknown ietest=',ietest
        call endrun (-1)
 
c============
        endif
c============
 
c~~~~~~~~~~
      enddo
c~~~~~~~~~~
c>>>>>
#endif
c>>>>>

c--------------------------------------------
c        Aerodynamic fluxes over sea-ice/snow
c--------------------------------------------
c
c        Set up index INDICE, NPICE (in comice) for pts with seaice.
c        INDICE itself is used only by VADAPT, called from ICEPHASE;
c        but NPICE>0 is tested here, in SURFCTL,SEAICE,ICEPHASE,...
c        Cray efficiency without index outweighs saving by eliminating
c        non-ice points.
c
      CALL WHENFGT (NPOI, FICE, 1, 0., INDICE, NPICE)
 
c-------------------------
      IF (NPICE.EQ.0) THEN
c-------------------------
 
        DO 395 I=1,NPOI
          FSENG(I) = 0.
          FVAPG(I) = 0.
          FSENI(I) = 0.
          FVAPI(I) = 0.
          TAUXI(I) = 0.
          TAUYI(I) = 0.
          HEATI(I) = 0.
          HEATG(I) = 0.
          ALOGAV(I) = 0.
          TAV(I) = TMELT
  395   CONTINUE
 
c---------
      ELSE
c---------
 
c--------------------------
      DO 400 ITER = 1,NITER
c--------------------------
 
c        Iterate NITER times for stability of upper seaice and snow
c        layers (as in turlsx). If first iteration, initialize surface
c        temperatures to sea-ice/snow top-layer temps (also set average
c        sfc roughness, first-layer heat conductivities, and slope of
c        ir flux for ice/snow skin).
 
      IF (ITER.EQ.1) THEN
        ZLOGG = LOG (Z0ICE)
        ZLOGI = LOG (Z0SNO)
        DO 410 I=1,NPOI
          ALOGAV(I) = (1.-FI(I))*ZLOGG + FI(I)*ZLOGI
          COG(I) = CONICE      / (0.5*MAX(HICE(I,1),HICETOP))
          COI(I) = CONSNO(I,1) / (0.5*MAX(HSNO(I,1),HSNOTOP))
          ZIRG(I) = 4.*EMG(I)*STEF*(TG_ICE(I)**3)
          ZIRI(I) = 4.*EMI(I)*STEF*(TI(I)**3)
          TGOLD(I) = TG_ICE(I)
          TIOLD(I) = TI(I)
  410   CONTINUE
      ENDIF
c
c        Set implicit/explicit factors W* (0 to 1) for this iteration.
c        W* is 1 for fully implicit, 0 for fully explicit. (For first
c        iteration, IMPEXP2 sets W* to 1.)
c
      CALL IMPEXP2 (WG, TG_ICE, TGOLD, ITER, 2)
      CALL IMPEXP2 (WI, TI,     TIOLD, ITER, 3)
 
c        Adjust T* for this iteration (see comments in turvap),
c        and save current "central" values for final flux calculations.
c        Also set specific humidities, average sfc temp and spec.hum.
c
      DO 420 I=1,NPOI
        TG_ICE(I) = WG(I)*TG_ICE(I) + (1.-WG(I))*TGOLD(I)
        TI(I)     = WI(I)*TI(I)     + (1.-WI(I))*TIOLD(I)
 
        TG_ICE(I) = MIN (TG_ICE(I),TMELT)
        TI(I)     = MIN (TI(I),    TMELT)
 
        TGPRE(I) = TG_ICE(I)
        TIPRE(I) = TI(I)
 
        ES      = ESAT (TG_ICE(I))
        QG(I) = QSAT (ES, PSURF(I))
        DQG(I)= DSAT (TG_ICE(I), PSURF(I))
        DQG(I)= MIN  (DQG(I), QG(I)/10.)
 
        ES      = ESAT (TI(I))
        QI(I) = QSAT (ES, PSURF(I))
        DQI(I)= DSAT (TI(I), PSURF(I))
        DQI(I)= MIN  (DQI(I), QI(I)/10.)
 
        TAV(I)  = (1.-FI(I))*TG_ICE(I)+ FI(I)*TI(I)
        QAV(I)  = (1.-FI(I))*QG(I)    + FI(I)*QI(I)
  420 CONTINUE
 
c        Calculate stratification factors for seaice/snow surface to za
 
      CALL FSTRAT (TAV, TA, TFAC, QAV, QA, VZERO, ZA,
     *             ALOGAV, ALOGAV, ALOGA, UA, RICHU2, STRAMU, STRAHU, 1)
 
c        Implicit solution for skin temperatures of sea ice and snow
 
      DO 425 I=1,NPOI
        CDH(I) = UA(I) * (VONK/(ALOGA(I)-ALOGAV(I)))**2 * STRAHU(I)
        CDM(I) = UA(I) * (VONK/(ALOGA(I)-ALOGAV(I)))**2 * STRAMU(I)
        CDH(I) = MIN (CDMAXA, CDH(I) / (1. + CDH(I)/CDMAXB))
        CDM(I) = MIN (CDMAXA, CDM(I) / (1. + CDM(I)/CDMAXB))
 
        SG(I) = RHOA(I)*CDH(I)
 
        TG_ICE(I) = (   SOLG(I) + FIRG(I) + ZIRG(I)*TGOLD(I)
     *                + SG(I)*CP(I)*TA(I)*TFAC
     *                - SG(I)*HVASUG(I)*(QG(I)-DQG(I)*TG_ICE(I) - QA(I))
     *                + COG(I)*TICE(I,1) )
     *            / (   SG(I)*(CP(I) + HVASUG(I)*DQG(I))
     *                + COG(I) + ZIRG(I) )
        TI(I)     = (   SOLI(I) + FIRI(I) + ZIRI(I)*TIOLD(I)
     *                + SG(I)*CP(I)*TA(I)*TFAC
     *                - SG(I)*HVASUI(I)*(QI(I)-DQI(I)*TI(I) - QA(I))
     *                + COI(I)*TSNO(I,1) )
     *            / (   SG(I)*(CP(I) + HVASUI(I)*DQI(I))
     *                + COI(I) + ZIRI(I) )
  425 CONTINUE
 
c        Sensible and latent fluxes, and wind stresses over seaice/snow,
c        and net surface heat fluxes for snow and seaice models.
c        Only need to compute these for final iteration.
 
      IF (ITER.EQ.NITER) THEN
 
        DO 430 I=1,NPOI
          TGAV = WG(I)*TG_ICE(I) + (1.-WG(I))*TGPRE(I)
ccc       FSENG(I) = SG(I) * CP(I) * (TGAV - TA(I)*TFAC)
 
          TIAV = WI(I)*TI(I) + (1.-WI(I))*TIPRE(I)
ccc       FSENI(I) = SG(I) * CP(I) * (TIAV - TA(I)*TFAC)
 
          QGAV = QG(I) + WG(I)*DQG(I)*(TG_ICE(I)-TGPRE(I))
ccc       FVAPG(I) = SG(I) * (QGAV - QA(I))
 
          QIAV = QI(I) + WI(I)*DQI(I)*(TI(I)-TIPRE(I))
ccc       FVAPI(I) = SG(I) * (QIAV - QA(I))
 
ccc       TAUXI(I) = RHOA(I) * CDM(I) * UX(I)
ccc       TAUYI(I) = RHOA(I) * CDM(I) * UY(I)
 
ccc       Partially time-implicit in ta,qa,ux,uy:
          ZMT = THICA(I)/DTIM
          ZFACH = ZMT * (1.-EXP(-RHOA(I)*CDH(I)*TFAC/ZMT))
          ZFACE = ZMT * (1.-EXP(-RHOA(I)*CDH(I)/ZMT))
          ZFACM = ZMT * (1.-EXP(-RHOA(I)*CDM(I)/ZMT))
          FSENG(I) = ZFACH * CP(I) * (TGAV/TFAC - TA(I))
          FSENI(I) = ZFACH * CP(I) * (TIAV/TFAC - TA(I))
          FVAPG(I) = ZFACE * (QGAV - QA(I))
          FVAPI(I) = ZFACE * (QIAV - QA(I))
          TAUXI(I) = ZFACM * UX(I)
          TAUYI(I) = ZFACM * UY(I)
 
          FIRG(I) = FIRG(I) - WG(I)*ZIRG(I)*(TG_ICE(I)-TGOLD(I))
          FIRI(I) = FIRI(I) - WI(I)*ZIRI(I)*(TI(I)    -TIOLD(I))
          FIRB(I) = FIRB(I) + FICE(I) *
     *              (  (1.-FI(I))*WG(I)*ZIRG(I)*(TG_ICE(I)-TGOLD(I))
     *               +     FI(I) *WI(I)*ZIRI(I)*(TI(I)    -TIOLD(I)) )
 
          HEATG(I) = SOLG(I) + FIRG(I) - FSENG(I) - HVASUG(I)*FVAPG(I)
          HEATI(I) = SOLI(I) + FIRI(I) - FSENI(I) - HVASUI(I)*FVAPI(I)
 
c         Impose final ice/snow constraints on (diagnostic) skin temps
c         TG,TI (eg, for next iterations's ir flux). Don't do earlier
c         since must use matrix result for flux calcs above.
 
          TG_ICE(I) = MIN (TG_ICE(I),TMELT)
          TI(I)     = MIN (TI(I),    TMELT)
  430   CONTINUE
 
      ENDIF
 
c---------------------------------------
c     End of aero sea-ice/snow iteration
 400  CONTINUE
c---------------------------------------
 
c----------
      ENDIF
c----------
 
c        Set surface rain/snow rates,temps, for snow and seaice models
 
      DO 500 I=1,NPOI
        RAINL(I)  = RAINA(I)
        TRAINL(I) = MAX (TA(I),TMELT)
        SNOWL(I)  = SNOWA(I)
        TSNOWL(I) = MIN (TA(I),TMELT)
        RAING(I)  = RAINL(I)
        TRAING(I) = TRAINL(I)
        SNOWG(I)  = SNOWL(I)
        TSNOWG(I) = TSNOWL(I)
  500 CONTINUE
 
c        Save overall surface aero fluxes (fsena,fvapa), and evap,
c        ocean roughness length, and 2-m sfc-air temperature, in commud
c        arrays (gsuvap,z0len,ts2) for history fields.
c        (nb: fsena is downward, fseno/g/i are upward, same for fvap* !)
 
      ZLOG2  = LOG(2.)
      ZLOG10 = LOG(10.)
      DO 600 I=1,NPOI
        FSENA(I) = - (1.-FICE(I))*            FSENO(I)
     *             -     FICE(I) *((1.-FI(I))*FSENG(I) + FI(I)*FSENI(I))
        FVAPA(I) = - (1.-FICE(I))*            FVAPO(I)
     *             -     FICE(I) *((1.-FI(I))*FVAPG(I) + FI(I)*FVAPI(I))
        TAUX(I)  = (1.-FICE(I))*TAUXO(I)
     *                +FICE(I) *TAUXI(I)
        TAUY(I)  = (1.-FICE(I))*TAUYO(I)
     *                +FICE(I) *TAUYI(I)
 
        GSUVAP(I)= -FVAPA(I)   ! gsuvap upwards, fvapa downwards
        Z0LEN(I) = EXP ( (1.-FICE(I))*ALOGGH(I)
     *                      +FICE(I) *ALOGAV(I) )
 
c       Calculation of ts2 is same as in subr tscreen in budlsx

        zweip  = weiprand (richu(i))
        zweip2 = weiprand (richu2(i))

        zwolog = (zlog2-aloggh(i)) / (aloga(i)-aloggh(i))
        zwolin = 2./za(i)
        zwo    = zweip*zwolog + (1.-zweip)*zwolin

        zwglog = (zlog2-alogav(i)) / (aloga(i)-alogav(i))
        zwglin = 2./za(i)
        zwg    = zweip2*zwglog + (1.-zweip2)*zwglin

        TS2(I) = (1.-FICE(I)) * (ZWO*TA(I) + (1.-ZWO)*TOCN(I))
     *              +FICE(I)  * (ZWG*TA(I) + (1.-ZWG)*TAV(I) )

#if defined (rcm)
c       calculation of [t,q,u]s10 is same as in subr tscreen in budlsx
        zwo10log = (zlog10-aloggh(i)) / (aloga(i)-aloggh(i))
        zwo10lin = 10./za(i)
        zwo10  = zweip*zwo10log + (1.-zweip)*zwo10lin

        zwg10log = (zlog10-alogav(i)) / (aloga(i)-alogav(i))
        zwg10lin = 10./za(i)
        zwg10  = zweip2*zwg10log + (1.-zweip2)*zwg10lin

        ts10(i) = (1.-fice(i)) * (zwo10*ta(i) + (1.-zwo10)*tocn(i))
     *               +fice(i)  * (zwg10*ta(i) + (1.-zwg10)*tav(i) )
        qs10(i) = (1.-fice(i)) * (zwo10*qa(i) + (1.-zwo10)*qocn(i))
     *               +fice(i)  * (zwg10*qa(i) + (1.-zwg10)*qav(i) )
        us10(i) = (1.-fice(i)) *  zwo10*ua(i)
     *               +fice(i)  *  zwg10*ua(i)
#endif

#ifdef isotrac
c       Save the values of FICE and of the vapor flux above sea-ice:
        fice_i(I) = FICE(I)
        fvi_i(I)  = -((1.-FI(I))*FVAPG(I) + FI(I)*FVAPI(I))
#endif

  600 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE SEAICE (DTICE)
 
c        Steps sea-ice vertical-column heat conduction through one step
 
c        DTICE  = thermodynamic seaice timestep in seconds
 
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
#include <comice.h>
#include <commud.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsno.h>
#include <comveg.h>
#include <com1d.h>
      DIMENSION FHTOP(NVEC), ZPEN(NVEC)
 
c        Set various diagnostics (commud) to zero for ocean/seaice,
c        and save seaice fraction (FICEOLD) and total thickness (HITOTO)
c        at beginning of timestep, for later adjustments in OCNML.
 
      DO 5 I=1,NPOI
        GINFIL(I) = 0.
        GINVAP(I) = 0.
        GTRANSU(I) = 0.
        GTRANSL(I) = 0.
        GEVAP(I)  = 0.
        GRUNOF(I) = 0.
        GABLAT(I) = 0.
        GADJUST(I)= 0.
        GDRAIN(I) = 0.
        GDRAINH(I) = 0.
        FICEOLD(I) = FICE(I)
        TFREEZ(I) = TMELT - .054 * 1000.*SOCN(I)  ! as in mom2,4 (mu_ts)
    5 CONTINUE
      DO K=1,NSOILAY+1
        DO I=1,NPOI
          WFLO(I,K) = 0.
        ENDDO
      ENDDO
c
      CALL ZERO (HITOTO, NPOI)
      DO 10 K=1,NICELAY
        DO 12 I=1,NPOI
          HITOTO(I) = HITOTO(I) + HICE(I,K)
   12   CONTINUE
   10 CONTINUE
c
c        Do nothing if no sea-ice points in this strip (initial creation
c        of sea ice in freezing open ocean done later in ocnml)
c
      IF (NPICE.EQ.0) RETURN
c
c        For sea ice, all incoming liquid flux runs off, but is first
c        cooled to TMELT (so all runoff is ge TMELT). Also
c        sublimation/accretion is removed/added from/to top seaice layer
 
!DIR$ IVDEP
      DO 58 I=1,NPOI
        IF (FICE(I).GT.0.) THEN
          FHTOP(I) = HEATG(I) + RAING(I)*CH2O*(TRAING(I)-TMELT)
          HICE(I,1) = HICE(I,1) - DTICE*FVAPG(I)/RHOICE
        ELSE
          FHTOP(I) = 0.
        ENDIF
   58 CONTINUE
 
c        Set penetrating solar fluxes across sea ice boundaries
 
      CALL ZERO (FPENICE, NVEC*(NICELAY+1))
      CALL ZERO (ZPEN, NPOI)
 
      DO 80 K=2,NICELAY+1
!DIR$ IVDEP
        DO 82 I=1,NPOI
          IF (FICE(I).GT.0.) THEN
            ZPEN(I) = ZPEN(I) + HICE(I,K-1)
#if defined (TWOSTRICE)
            fpenice(i,k) = fluxiced(i,k) + fluxicei(i,k)
#else
            FPENICE(I,K) = PENFR*SOLG(I)*EXP(-ZPEN(I)/PENDE)
#endif
          ENDIF
   82   CONTINUE
   80 CONTINUE
 
c        Step seaice temperatures due to conductive heat transport
 
      CALL ICEHEAT (FHTOP, FPENICE, HFLOICE, DTICE)
 
c        Thermodynamically adjust seaice for surface and internal
c        melting, and bottom melt/freeze
 
      CALL ICEPHASE (DTICE)
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE ICEHEAT (FHTOP, FPENICE, HFLOICE, DTIM)
c
c        Sets up call to TRIDIA to solve implicit sea-ice heat
c        conduction, using layer temperatures in TICE (in comice).
c        Lower bc is temperature = TFREEZ(I) at base of seaice.
c
c        All arguments are supplied except HFLOICE (returned).
c        FHTOP = heat flux into top layer from atmos
c        FPENICE= penetrating solar flx across layer bdries (0 at top)
c        HFLOICE=downward heat flow across layer boundaries
c        DTIM  = time step for this call
c
c        Local arrays and scalars:
c        OTICE  = ice temperatures at start of timestep (for HFLOICE)
c        CON    = conduction coefficents between layers
c        D1,2,3 = diagonals of tridiagonal systems of equations
c        RHS    = right-hand sides of systems of equations
c        W1,2   = work arrays needed by TRIDIA
c        RIMP   = implicit fraction of the calculation (0 to 1)
c        HFAKE  = arbitrary small thickness to allow processing
c                 for zero ice.
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <comocn.h>
      DIMENSION
     *  FHTOP(NVEC),           FPENICE(NVEC,NICELAY+1),
     *  HFLOICE(NVEC,NICELAY+1)
 
      DIMENSION
     *  OTICE(NVEC,NICELAY),   CON(NVEC,NICELAY+1),
     *  D1(NVEC,NICELAY),      D2(NVEC,NICELAY),  D3(NVEC,NICELAY),
     *  RHS(NVEC,NICELAY),     W1(NVEC,NICELAY),  W2(NVEC,NICELAY)
c
      SAVE RIMP, HFAKE
      DATA RIMP, HFAKE /1.0, .01/
c
c        Save pre-step temperatures (for calc of HFLOICE below)
c
      CALL SCOPY (NVEC*NICELAY, TICE, 1, OTICE, 1)
c
c        Set conduction coefficient between layers
c
      DO 100 K=1,NICELAY+1
c
        IF (K.EQ.1) THEN
          CALL ZERO (CON(1,K), NPOI)
c
        ELSE IF (K.LE.NICELAY) THEN
          DO 102 I=1,NPOI
            CON(I,K) =  1. / (  0.5*MAX(HICE(I,K-1),HFAKE)/CONICE
     *                        + 0.5*MAX(HICE(I,K),  HFAKE)/CONICE )
  102     CONTINUE
c
        ELSE IF (K.EQ.NICELAY+1) THEN
          DO 104 I=1,NPOI
            CON(I,K) = CONICE/(0.5*MAX(HICE(I,K-1),HFAKE))
  104     CONTINUE
        ENDIF
c
  100 CONTINUE
c
c        Set diagonals of matrix and right-hand side
c
      DO 200 K=1,NICELAY
        KM1 = MAX (K-1,1)
        KP1 = MIN (K+1,NICELAY)
c
        DO 202 I=1,NPOI
          DT = DTIM / (RHOICE*CSICE*MAX(HICE(I,K),HFAKE))
          D1(I,K) =    - DT * RIMP * CON(I,K)
          D2(I,K) = 1. + DT * RIMP * (CON(I,K)+CON(I,K+1))
          D3(I,K) =    - DT * RIMP * CON(I,K+1)
c
          RHS(I,K) = TICE(I,K) + DT
     *               * ( (1.-RIMP)*CON(I,K)  *(TICE(I,KM1)-TICE(I,K))
     *                 + (1.-RIMP)*CON(I,K+1)*(TICE(I,KP1)-TICE(I,K))
     *                 + FPENICE(I,K) - FPENICE(I,K+1) )
  202   CONTINUE
c
        IF (K.EQ.1) THEN
          DO 204 I=1,NPOI
            DT = DTIM / (RHOICE*CSICE*MAX(HICE(I,K),HFAKE))
            RHS(I,K) = RHS(I,K) + DT*FHTOP(I)
  204     CONTINUE
        ENDIF
c
        IF (K.EQ.NICELAY) THEN
          DO 206 I=1,NPOI
            DT = DTIM / (RHOICE*CSICE*MAX(HICE(I,K),HFAKE))
            RHS(I,K) = RHS(I,K) + DT * CON(I,K+1)
     *                 * (TFREEZ(I) - (1.-RIMP)*TICE(I,KP1))
            D3(I,K) = 0.
  206     CONTINUE
        ENDIF
c
  200 CONTINUE
c
c        Solve systems of equations
c
      CALL TRIDIA (NPOI, NVEC, NICELAY, D1,D2,D3, RHS, TICE, W1,W2)
c
c        Deduce downward heat fluxes between layers
c
      CALL SCOPY (NPOI, FHTOP, 1, HFLOICE(1,1), 1)
c
      DO 400 K=1,NICELAY
        DO 402 I=1,NPOI
          HFLOICE(I,K+1) = HFLOICE(I,K) - RHOICE*CSICE*HICE(I,K)
     *                                    *(TICE(I,K)-OTICE(I,K))/DTIM
 
  402   CONTINUE
  400 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE ICEPHASE (DTICE)
c
c        Thermodynamically adjusts seaice model for surface and
c        internal melting, and bottom melt/freeze, through one timestep.
c        Also contributes to the accumulated ocean fluxes F*OCN that
c        will be used by the ocean model.
c
c        DTICE  = time step for this call
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <com1d.h>
      DIMENSION ZBASE(NVEC), XH(NVEC), HITOT(NVEC)
 
c     HFAKE  = arbitrary small thickness to allow processing for no ice
      SAVE HFAKE
      DATA HFAKE /.01/
c
c        If all open ocean points in this strip, do nothing
c        (for safety only - already skipped out at start of SEAICE)
c
      IF (NPICE.EQ.0) RETURN
c
c        If temperature of each layer is gt TMELT, give fraction
c        FBRINSTO of the excess heat to the brine resevoir, and
c        1-FBRINSTO to meltwater production. (Assume each layer cannot
c        be entirely melted in one timestep.)
c
      CALL ZERO (SITOP,NVEC)
      CALL ZERO (SIBASE,NVEC)
c
      DO 100 K=1,NICELAY
!DIR$ IVDEP
        DO 102 I=1,NPOI
          XH(I) = RHOICE*CSICE*HICE(I,K) * MAX(TICE(I,K)-TMELT,0.)
          TICE(I,K) = MIN (TICE(I,K),TMELT)
c
          BRINE(I)  = BRINE(I)  + FBRINSTO*XH(I)
          HICE(I,K) = HICE(I,K) - (1.-FBRINSTO)*XH(I)/(RHOICE*HFUS)
          SITOP(I)  = SITOP(I)  + (1.-FBRINSTO)*XH(I)/(HFUS*DTICE)
  102   CONTINUE
  100 CONTINUE
c
c        As long as there is any heat in the brine resevoir, use it to
c        keep layer temperatures at TMELT. Apply from top to bottom.
c
      DO 120 K=1,NICELAY
!DIR$ IVDEP
        DO 122 I=1,NPOI
          XH(I) = MIN(RHOICE*CSICE*HICE(I,K)*(TMELT-TICE(I,K)),BRINE(I))
          TICE(I,K) = TICE(I,K)
     *                + XH(I) / (RHOICE*CSICE*MAX(HFAKE,HICE(I,K)))
          BRINE(I) = BRINE(I) - XH(I)
  122   CONTINUE
  120 CONTINUE
c
c        If brine resevoir exceeds capacity = FBRINMAX * heat needed
c        to melt entire column, use excess heat to melt top-layer ice.
c        At this point, if BRINE > 0, all layer temperatures must be at
c        TMELT, so it doesn't matter which layer is melted (VADAPT will
c        be called soon below), as long as the layer cannot be entirely
c        melted in one timestep. Actually use more than the excess heat
c        amount by a factor 1/(1-FBRINMAX), which results from a simple
c        calculation allowing for the reduction in the capacity itself
c        due to the ice melted. (This factor corresponds to reducing
c        the latent heat of fusion as in Semtner,1976.)
c
c        First calculate heat needed to melt entire column.
c
      CALL ZERO (XH, NVEC)
      DO 140 K=1,NICELAY
!DIR$ IVDEP
        DO 142 I=1,NPOI
          XH(I)= XH(I) + RHOICE*HICE(I,K)*((TMELT-TICE(I,K))*CSICE+HFUS)
  142   CONTINUE
  140 CONTINUE
c
c        Now any excess brine-heat melts top layer (at TMELT if BRINE>0)
c
!DIR$ IVDEP
      DO 145 I=1,NPOI
        DBRINE = MAX ( 0., (BRINE(I) - FBRINMAX*XH(I)) / (1.-FBRINMAX) )
        DBRINE = MIN ( DBRINE, BRINE(I) )
c
        BRINE(I)  = BRINE(I)  - DBRINE
        HICE(I,1) = HICE(I,1) - DBRINE/(RHOICE*HFUS)
        SITOP(I)  = SITOP(I)  + DBRINE/(HFUS*DTICE)
  145 CONTINUE
c
c        Melt or freeze at seaice base, with adjustment to latent heat
c        (HFREEZ) for difference between temp of bottom ice layer and 
c        TFREEZ, to conserve heat in ice-water-ice cycles.
c        HFLOICE,FPENICE (commus) were set in ICEHEAT and SEAICE resp;
c        must subtract escaping penetrating solar flux FPENICE since
c        HFLOICE includes it but it passes through the seaice base
c        to the ocean (added to FTOCN below) without causing bottom
c        melting.
c
!DIR$ IVDEP
      DO 240 I=1,NPOI
        ZBASE(I) = BASEFAC * (TOCN(I)-TFREEZ(I))
        DH = DTICE * ( HFLOICE(I,NICELAY+1) - FPENICE(I,NICELAY+1)
     *                 + ZBASE(I) )
        HFREEZ = HFUS - (CH2O-CSICE)*(TMELT-TFREEZ(I))
        DW = DH / ( CSICE*(TFREEZ(I)-TICE(I,NICELAY))+HFREEZ )
        HICE(I,NICELAY) = HICE(I,NICELAY) - DW/RHOICE
        SIBASE(I) = DW/DTICE
  240 CONTINUE
c
c        Re-adapt seaice thickness profile, so top thickness = HICETOP
c        and other thicknesses are equal.
c
      CALL VADAPT (HICE, TICE, HICETOP, INDICE, NPICE, NICELAY)
c
c        Accumulate ocean heat and h2o fluxes (under-seaice contrib),
c        and wind stress on ice surface. (Don't call MIX to mix the
c        various liquid fluxes, since SIBASE can be negative due to
c        freezing onto ice bottom.)
c
!DIR$ IVDEP
      DO 400 I=1,NPOI
        FTOCN(I) = FTOCN(I)
     *             + (   FPENICE(I,NICELAY+1)
     *                 - ZBASE(I)
     *                 - RAING(I)*CH2O*(TOCN(I)-TMELT)
     *                 - SITOP(I) *CH2O*(TOCN(I)-TMELT)
     *                 - SIBASE(I)*CH2O*(TOCN(I)-TFREEZ(I)) )
     *               * FICE(I)
        FWOCN(I) = FWOCN(I) +(RAING(I)+SITOP(I)+SIBASE(I))*FICE(I) !csal
        FUICE(I) = TAUXI(I)
        FVICE(I) = TAUYI(I)
  400 CONTINUE

      CALL ZERO (HITOT,NPOI)
      DO K=1,NICELAY
        DO I=1,NPOI
          HITOT(I) = HITOT(I) + HICE(I,K)
        ENDDO
      ENDDO

!DIR$ IVDEP
      DO I=1,NPOI
c       for mom2,4:
        FHOCN(I)   = FHOCN(I)   + ZBASE(I) * FICE(I)
c       FQOCN(I)   = FQOCN(I)   + 0.
        FSOLOCN(I) = FSOLOCN(I) + FPENICE(I,NICELAY+1) * FICE(I)
c       FIROCN(I)  = FIROCN(I)  + 0.
        FPRLOCN(I) = FPRLOCN(I) + RAING(I) * FICE(I)
c       FPRFOCN(I) = FPRFOCN(I) + 0.
        FSALOCN(I) = FSALOCN(I) -(SITOP(I)+SIBASE(I)) * FICE(I) 
        PRESOCN(I) = PRESOCN(I) + (PSURF(I)+RHOI*GRAV*HITOT(I)) *FICE(I)
      ENDDO

      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE OCNML (DTGCM)
c
c        Steps ML temperature through one agcm timestep. Also adjusts
c        seaice in ways that can alter its fractional coverage FICE.
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <commug.h>
#include <commus.h>
#include <comocn.h>
#include <comsno.h>
#include <com1d.h>
c
      DIMENSION
     *  FICETMP(NVEC),
     *  XM(NVEC),       XH(NVEC),        DM(NVEC),       DH(NVEC),
     *  DFI(NVEC),      RATIO(NVEC),     FX(NVEC),
     *  HITOT(NVEC),    Z(NVEC),         ZMSNO(NVEC),    ZTSNO(NVEC),
     *  ZHSNO(NVEC),    ZZICE(NVEC),     ZMICE(NVEC)
c
c        Set flag if any seaice points so can skip some sections below
c        if none in current strip. Don't use an index since it slows
c        Cray speed and makes vectorization harder, and seaice regions
c        are clustered anyway. Still need to protect all code below
c        for non-ice points, since most strips will contain some.
c
      IZ = ISMAX (NPOI, FICE, 1)
      IF (FICE(IZ).GT.0.) IFANY = 1
c
c        Step ocean mixed-layer temp due to net sfc heat flux FTOCN
c        and oceanic heat flux convergence QFLUX (set in SETQFLUX)
c
      DO 100 I=1,NPOI
        TOCN(I) = TOCN(I)
     *          + DTGCM * (FTOCN(I)+QFLUX(I)) / (RHOW*CH2O*HOCN(I))
        tocn(i) = min (tocn(i), tmelt + 100.)    ! for hoigh obliq
  100 CONTINUE
c
c        Decrement ice fraction implied by accumulated vertical melting
c        in SEAICE since last call to OCNML. Previous total thickness
c        saved in HITOTO in SEAICE. Implicitly assumes linear
c        distribution of thicknesses between 0 and 2*HITOTO, as in
c        Hibler(1979) and Harvey(1988), or triangular distribution
c        between 0 and 2*HITOTO, as suggested by Steve Vavrus (U.Wisc)
c
      zdecrem = 1.
      IF (IFANY.EQ.1 .and. zdecrem.ne.0.) THEN
 
        CALL ZERO (HITOT,NPOI)
        DO 150 K=1,NICELAY
          DO 152 I=1,NPOI
            HITOT(I) = HITOT(I) + HICE(I,K)
  152     CONTINUE
  150   CONTINUE
 
        DO 160 I=1,NPOI
          ZZ = MAX(HITOTO(I)-HITOT(I),0.) / MAX(HITOTO(I),EPSILON)
c         RATIO(I) = 1. -  ZZ     / 2.         ! linear distribution
          RATIO(I) = 1. - (ZZ*ZZ) / 2.         ! triangular distribution
          FICE(I) = FICE(I)*RATIO(I)
  160   CONTINUE
 
        DO 165 K=1,NICELAY
          DO 167 I=1,NPOI
            HICE(I,K) = HICE(I,K)/RATIO(I)
  167     CONTINUE
  165   CONTINUE
 
      ENDIF
c
c        Do lateral seaice melt. First calculate XM = seaice column
c        mass, and XH = heat required to melt a seaice column to ocean
c        temperature TOCN. Then increment FICE, FTOCN, FWOCN.
c        (Even if this made all seaice pts vanish, don't reset IFANY
c        since would still need to do final snow adjustments at end.)
c
      IF (IFANY.EQ.1) THEN
c
        CALL ZERO (XH,NVEC)
        CALL ZERO (XM,NVEC)
        DO 200 K=1,NICELAY
          DO 202 I=1,NPOI
            XM(I) = XM(I) + RHOICE*HICE(I,K)
            XH(I) = XH(I) + RHOICE*HICE(I,K)
     *                    * (   CSICE*(TMELT-TICE(I,K))+ HFUS
     *                        + CH2O*(TOCN(I)-TMELT) )
  202     CONTINUE
  200   CONTINUE
c
c          DH(I)=... line says how much ocean heat melts seaice sides.
c          Its dep on FICE is like Park+Wash, but with small value
c          SIDEZER so it doesn't tend to zero at same rate as FICE,
c          and so can make seaice vanish. It depends on ocean temp
c          rather than FTOCN, more physical than Park+Wash.
c
        DO 220 I=1,NPOI
          DH(I) = DTGCM*SIDEFAC*MAX(TOCN(I)-TFREEZ(I), 0.)
     *            *(1.-FICE(I))*(FICE(I)+SIDEZER)
          DH(I) = MIN (DH(I), XH(I)*FICE(I))
          DFI(I) = DH(I) / MAX(XH(I),EPSILON)
          FICE(I) = FICE(I) - DFI(I)
          TOCN(I) = TOCN(I) - DH(I)/(RHOW*CH2O*HOCN(I))
          FTOCN(I) = FTOCN(I) - DH(I)/DTGCM
          FWOCN(I) = FWOCN(I) + DFI(I)*XM(I)/DTGCM        ! csal
c         for mom2,4:
          FHOCN(I)   = FHOCN(I)   + DH(I)/DTGCM
          FSALOCN(I) = FSALOCN(I) - DFI(I)*XM(I)/DTGCM
  220   CONTINUE
c
      ENDIF
c
c        If mixed layer is below freezing, form new seaice in leads
c        with properties *CRE (in comocn)...ie, with a certain fixed
c        thickness and with temperature TFREEZ. Then combine with
c        existing ice by summing existing and new fractional coverages,
c        and finally adjusting layer thicknesses and temperatures.
c
c        Resulting FICE + DFI could exceed 1...taken care of below.
c
      DO 300 I=1,NPOI
        DH(I) = RHOW*CH2O*HOCN(I) * MAX (TFREEZ(I)-TOCN(I), 0.)
        HFREEZ = HFUS - (CH2O-CSICE)*(TMELT-TFREEZ(I))
        DFI(I) = DH(I)/(ZMCRE*HFREEZ)
        TOCN(I) = MAX (TOCN(I), TFREEZ(I))
        FTOCN(I) = FTOCN(I) + DH(I)/DTGCM
        FWOCN(I) = FWOCN(I) - DFI(I)*ZMCRE/DTGCM        ! csal
c       for mom2,4:
        FHOCN(I)   = FHOCN(I)   - DH(I)/DTGCM
        FSALOCN(I) = FSALOCN(I) + DFI(I)*ZMCRE/DTGCM
c
c       Following adjustments needed to protect EPSILON logic below
        FICE(I) = CVMGT (0., FICE(I), FICE(I).LT.EPSILON)
        DFI(I) = CVMGT (0., DFI(I), DFI(I).LT.EPSILON)
c
        FICETMP(I) = FICE(I) + DFI(I)
  300 CONTINUE
c
c        Set IFANY to 1 in case first seaice point(s) for this strip
c        is being formed by mixed layer freezing. If not and IFANY is
c        0, don't have to do any more of this section (loops 320,330).
c
      IZ = ISMAX (NPOI, FICETMP, 1)
      IF (FICETMP(IZ).GT.0.) IFANY = 1
c
      IF (IFANY.EQ.1) THEN
c
        DO 320 K=1,NICELAY
          DO 322 I=1,NPOI
            HICENEW = (FICE(I)*HICE(I,K) + DFI(I)*HICECRE(K))
     *              / MAX (EPSILON, FICETMP(I))
c
            TICE(I,K) = (  FICE(I)*HICE(I,K)*TICE(I,K)
     *                    + DFI(I)*HICECRE(K)*TFREEZ(I) )
     *                / MAX (EPSILON*.001, FICETMP(I)*HICENEW)
c
            HICE(I,K) = HICENEW
  322     CONTINUE
  320   CONTINUE
c
        DO 330 I=1,NPOI
          FICE(I) = FICETMP(I)
  330   CONTINUE
c
      ENDIF
c
c        Ensure total ice thickness is ge minimum value HICEMIN, by
c        shrinking FICE if necessary
c
      IF (IFANY.EQ.1) THEN
c
        CALL ZERO (XH, NVEC)
        DO 400 K=1,NICELAY
          DO 402 I=1,NPOI
             XH(I) = XH(I) + HICE(I,K)
  402     CONTINUE
  400   CONTINUE
c
        DO 410 I=1,NPOI
          RATIO(I) =  MIN (XH(I),HICEMIN) / HICEMIN
          FICE(I) = FICE(I) * RATIO(I)
  410   CONTINUE
c
        DO 420 K=1,NICELAY
          DO 422 I=1,NPOI
            HICE(I,K) = HICE(I,K) / MAX (RATIO(I), EPSILON)
  422     CONTINUE
  420   CONTINUE
c
      ENDIF
c
c        Ensure fractional ice coverage is le maximum value
c        FICEMAX[N,S], by increasing thicknesses if necessary
c
      IF (IFANY.EQ.1) THEN
c
        CALL ZERO (HITOT,NPOI)
        DO 500 K=1,NICELAY
          DO 502 I=1,NPOI
            HITOT(I) = HITOT(I) + HICE(I,K)
  502     CONTINUE
  500   CONTINUE
c
        DO 510 I=1,NPOI
          FX(I)= CVMGT (FICEMAXS, FICEMAXN, LOCPOI(I,2).LE.NLAT/2)
#ifdef gcm
c         Dependence on thickness from Harvey(1988),J.Climate,p.1223.
c         Not for SNOWBALL (allow fice to be 1.0)
#  if ! defined (SNOWBALL)
          FX(I) = 1. - (1.-FX(I))* EXP(-MAX(HITOT(I)-1.,0.)/HICEEXP)
#  endif
#endif
          RATIO(I) = MAX (FICE(I),FX(I)) / FX(I)
          FICE(I) = MIN (FICE(I),FX(I))
  510   CONTINUE
c
        DO 520 K=1,NICELAY
          DO 522 I=1,NPOI
            HICE(I,K) = HICE(I,K)*RATIO(I)
  522     CONTINUE
  520   CONTINUE
c
      ENDIF
c
c        Although we don't have a minimum fractional ice coverage,
c        protect against arbitrarily small machine-values by killing
c        seaice (and snow) if FICE lt EPSILON, similar to adjustments
c        at end of DYNICE. Do before last snow adjustments below, since
c        they don't change FICE and it protects them for no-ice points.
c
      IF (IFANY.EQ.1) THEN
c
        DO 600 I=1,NPOI
         Z(I) = CVMGT (0., 1., FICE(I).LT.EPSILON)
  600   CONTINUE
c
        DO 602 I=1,NPOI
          FICE(I) = FICE(I)*Z(I)
          FI(I) = FI(I)*Z(I)
          BRINE(I) = BRINE(I)*Z(I)
  602   CONTINUE
c
        DO 604 K=1,NICELAY
          DO 606 I=1,NPOI
            HICE(I,K) = HICE(I,K)*Z(I)
            TICE(I,K) = TICE(I,K)*Z(I) + TMELT*(1.-Z(I))
  606     CONTINUE
  604   CONTINUE
c
        DO 608 K=1,NSNOLAY
          DO 610 I=1,NPOI
            HSNO(I,K) = HSNO(I,K)*Z(I)
            TSNO(I,K) = TSNO(I,K)*Z(I) + TMELT*(1.-Z(I))
  610     CONTINUE
  608   CONTINUE
c
      ENDIF
c
c        Compute some column-integrated quantities for snow and seaice,
c        needed for snow adjustments further below
c
      IF (IFANY.EQ.1) THEN
c
        CALL ZERO (ZMSNO, NVEC)
        CALL ZERO (ZTSNO, NVEC)
        CALL ZERO (ZHSNO, NVEC)
        DO 700 K=1,NSNOLAY
          DO 702 I=1,NPOI
            ZMSNO(I) = ZMSNO(I) + RHOS*HSNO(I,K)
            ZTSNO(I) = ZTSNO(I) + RHOS*HSNO(I,K)*TSNO(I,K)
            ZHSNO(I) = ZHSNO(I)
     *                + RHOS*HSNO(I,K)* ( HFUS + CSICE*(TMELT-TSNO(I,K))
     *                                         + CH2O*(TOCN(I)-TMELT) )
  702     CONTINUE
  700   CONTINUE
c
        CALL ZERO (ZZICE, NVEC)
        CALL ZERO (ZMICE, NVEC)
        DO 704 K=1,NICELAY
          DO 706 I=1,NPOI
            ZZICE(I) = ZZICE(I) + HICE(I,K)
            ZMICE(I) = ZMICE(I) + RHOICE*HICE(I,K)
  706     CONTINUE
  704   CONTINUE
c
c       If net increase in ice coverage, stretch snow & brine resevoir.
c       If net decrease in ice coverage, dump some snow+brine into ocn;
c       (could make TOCN(I) < TFREEZ(I) - tolerate till next call).
c       (FICEOLD = FICE before any changes was saved in SEAICE).
c
        DO 720 I=1,NPOI
c
          IF (FICE(I).GT.FICEOLD(I)) THEN
            FI(I) = FI(I) * FICEOLD(I)/FICE(I)
            BRINE(I) = BRINE(I) * FICEOLD(I)/FICE(I)
c
          ELSE IF (FICE(I).LT.FICEOLD(I)) THEN
            DM(I) = (FICEOLD(I)-FICE(I))*   FI(I)*ZMSNO(I)
            DH(I) = (FICEOLD(I)-FICE(I))* (-FI(I)*ZHSNO(I) + BRINE(I))
            TOCN(I) = TOCN(I)   + DH(I)/(RHOW*CH2O*HOCN(I))
            FTOCN(I) = FTOCN(I) + DH(I)/DTGCM
            FWOCN(I) = FWOCN(I) + DM(I)/DTGCM      ! csal
c           for mom2,4:
            FHOCN(I)    = FHOCN(I)   - DH(I)/DTGCM
            FSALOCN(I)  = FSALOCN(I) - DM(I)/DTGCM
c           fsalocn only for sea-ice freeze/melt (not snow)?:
c           FSALOCN(I)  = FSALOCN(I) + 0.
          ENDIF
c
  720   CONTINUE
c
c          If "white ice", ie, snow-ice boundary below waterline,
c          convert some snow to ice to restore bdry to waterline,
c          and give to upper ice layer. (Take from entire snow
c          column and not from bottom snow layer, to avoid making
c          the bottom snow layer thickness negative.)
c
        DO 740 I=1,NPOI
          DM(I) = 0.
          IF ((ZMICE(I)+FI(I)*ZMSNO(I))/RHOW .GT. ZZICE(I)) THEN
            DM(I) = (RHOICE/RHOW)*FI(I)*ZMSNO(I)
     *            - (1.-(RHOICE/RHOW)) *ZMICE(I)
            TICE(I,1)  = ( HICE(I,1)*TICE(I,1)
     *                     + (ZTSNO(I)/ZMSNO(I))*(DM(I)/RHOICE) )
     *                 / ( HICE(I,1) + DM(I)/RHOICE )
            HICE(I,1) = HICE(I,1) + DM(I)/RHOICE
          ENDIF
  740   CONTINUE
c
        DO 742 K=1,NSNOLAY
          DO 744 I=1,NPOI
            HSNO(I,K) = HSNO(I,K)
     *                * ( 1. - DM(I) / MAX(EPSILON,FI(I)*ZMSNO(I)) )
  744     CONTINUE
  742   CONTINUE
c
      ENDIF
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE VADAPT (HCUR, TCUR, HTOP, INDP, NP, NLAY)
c
c        Re-adapt snow or seaice layer thicknesses, so top thickness
c        equals HSNOTOP or HICETOP and other thicknesses are equal.
c        Also adjusts profile of tracer field TCUR so its vertical
c        integral is conserved (eg, temperature).
c
c        HCUR = layer thicknesses (supplied and returned)
c        TCUR = tracer field (supplied and returned)
c        HTOP = prescribed top layer thickness (supplied)
c        INDP = index of snow or seaice pts in current strip (supplied)
c        NP   = number of snow or seaice pts in current strip (supplied)
c        NLAY = number of layers (=NSNOLAY or NICELAY) (supplied)
c
#include <compar.h>
#include <comgrd.h>
c
      DIMENSION HCUR(NVEC,NLAY),       TCUR(NVEC,NLAY),
     *          INDP(NVEC)
c
      INTEGER PLAY
#if defined  ( cray ) || defined ( linux )
      PARAMETER (PLAY = NSNOLAY+NICELAY)
#else
      PARAMETER (PLAY = MAX(NSNOLAY,NICELAY))
#endif
c
      DIMENSION
     *  HNEW(NVEC,PLAY),    TNEW(NVEC,PLAY),
     *  ZOLD(NVEC,PLAY+1),  HT(NVEC),           H1(NVEC),
     *  ZA(NVEC),           ZB(NVEC),           ZHEAT(NVEC)
c
c
c        If no snow or seaice points in current 1D strip, return. Note
c        that the index is not used below (for Cray vec and efficiency)
c        except in the final loop setting the returned values.
c
      IF (NP.EQ.0) RETURN
c
c        Set distances ZOLD from surface to old layer boundaries
c
      CALL ZERO (ZOLD(1,1), NPOI)
c
      DO 300 K=1,NLAY
        DO 302 I=1,NPOI
          ZOLD(I,K+1) = ZOLD(I,K) + HCUR(I,K)
  302   CONTINUE
  300 CONTINUE
c
c        Set new layer thicknesses HNEW (tot thickness is unchanged).
c        If total thickness is less than NLAY*HTOP (which should be
c        le HSNOMIN or HICEMIN), make all new layers equal including
c        top one, so other layers aren't so thin. Use EPSILON to
c        handle no-snow or no-seaice points.
c
      CALL SCOPY (NPOI, ZOLD(1,NLAY+1), 1, HT, 1)
c
      DO 304 I=1,NPOI
        IF (HT(I).GE.NLAY*HTOP) THEN
          H1(I) = (HT(I)-HTOP)/(NLAY-1)
        ELSE
          H1(I) = MAX (HT(I)/NLAY, EPSILON)
        ENDIF
  304 CONTINUE
c
      DO 306 K=1,NLAY
        DO 308 I=1,NPOI
          HNEW(I,K) = H1(I)
  308   CONTINUE
  306 CONTINUE
c
      DO 310 I=1,NPOI
        IF (HT(I).GE.NLAY*HTOP) HNEW(I,1) = HTOP
  310 CONTINUE
c
c        Integrate old temperature profile (loop 410) over each
c        new layer (loop 400), to get new field TNEW
c
      CALL ZERO (ZB, NPOI)
c
      DO 400 K=1,NLAY
c
        DO 402 I=1,NPOI
          ZA(I) = ZB(I)
          ZB(I) = ZA(I) + HNEW(I,K)
          ZHEAT(I) = 0.
  402   CONTINUE
c
        DO 410 KO=1,NLAY
          DO 412 I=1,NPOI
            IF (ZA(I).LT.ZOLD(I,KO+1) .AND. ZB(I).GT.ZOLD(I,KO)) THEN
              DZ = MIN(ZOLD(I,KO+1),ZB(I)) - MAX(ZOLD(I,KO),ZA(I))
              ZHEAT(I) = ZHEAT(I) + TCUR(I,KO)*DZ
            ENDIF
  412     CONTINUE
  410   CONTINUE
c
        DO 420 I=1,NPOI
          TNEW(I,K) = ZHEAT(I) / HNEW(I,K)
  420   CONTINUE
c
  400 CONTINUE
c
c        Use index for final copy to seaice or snow arrays, to avoid
c        changing nominal values for no-seaice or no-snow points
c
      DO 500 K=1,NLAY
        DO 502 J=1,NP
          I = INDP(J)
          HCUR(I,K) = HNEW(I,K)
          TCUR(I,K) = TNEW(I,K)
  502   CONTINUE
  500 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      subroutine setqflux (qfluxxy, tocnxy,  ficexy, tfreezxy, lmaskxy,
     *                     qfactor, qnorweg, qices,  qicen, qlimit,
     *                     mcdate,  showmaps)
c
c        If qfactor >= 0 and <= 100:
c
c        Computes oceanic heat flux convergence (modeled as linear
c        diffusion), and adjusts it in various ways that depend on
c        season (Norwegian Sea) or current conditions (under sea ice).
c        Then makes a uniform correction so that its (instantaneous!)
c        global integral is zero.
c
c        If qfactor < 0:
c
c        Reads OHT.DAT for (lat(deg), oht(Pw)), computes (zonally
c        symmetric) qflux. Does same adjustments as above (multiplying
c        by abs(qfactor).
c
c        If qfactor > 100:
c
c        Reads QFLUX global map for convergence (W/m2). Does same
c        adjustments as above (multiplying by qfactor-100).
c
c        QFLUXXY = oceanic heat flux convergence, W/m2 (returned)
c        TOCNXY  = ocean mixed-layer temperature, deg K (supplied)
c        FICEXY  = sea ice fraction (supplied)
c        TFREEZXY= freezing pt of mixed layer (supplied)
c        LMASKXY = lsx land-icesheet-ocean mask (supplied, for 1st call)
c        QFACTOR = q-flux factor multiplying standard diffusion coeff
c        QNORWEG = max value (W m-2) for Norwegian-Sea q-flux
c        QICES   = heat convergence under 100% sea ice, southern hemis
c        QICEN   = heat convergence under 100% sea ice, northern hemis
c        MCDATE  = current date (at end of time step) (yymmdd)
c        SHOWMAPS= .TRUE. to print 2-D maps of prescribed fields
c----
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <commapl.h>
#include <comocn.h>
c----
      dimension qfluxxy(nlon,nlat),  tocnxy(nlon,nlat),
     *          ficexy(nlon,nlat),   tfreezxy(nlon,nlat),
     *          lmaskxy(nlon,nlat)
      logical showmaps
c-----
      dimension xmask(nlon,nlat),    ymask(nlon,nlat+1),
     *          coeflat(nlat+1),     zonfrac(nlat),
     *          darea(nlat),         zonf(nlat+1),
     *          zonq(nlat),          zont(nlat),
     *          zond(nlat+1)
      save xmask, ymask, coeflat, zonfrac, dlon, dlat, darea
      character*80 ctmp
      logical firsta, firstb
      save firsta, firstb
      data firsta, firstb /.true., .true./
c-----
      dimension alohtin(1000), ohtin(1000), oht(nlat+1), qfluxoht(nlat)
      save qfluxoht
c-----
      dimension qfluxxyin(nlon,nlat)
      save qfluxxyin
      character ckey*8, comment*100
c-----

c      coefstd = 2.0e13                         !  (W K-1)  ! v2
c      coefstd = 4.0e13                         !  (W K-1)  ! v2.3
c      dependent on lsx resolution (runpres vs runpresa, 3/07):
       if (nlon.ge.180 .and. nlat.ge.90) then
         coefstd = 5.2e13                       !  (W K-1)  ! v3
       else
         coefstd = 5.2e13 * 1.2                 !  (W K-1)  ! v3
       endif
 
c++++++++++++++++++++++++++++++++++++++++++++++++++
      if (qfactor.ge.0. .and. qfactor.le.100.) then
c++++++++++++++++++++++++++++++++++++++++++++++++++
 
c        If qfactor >= 0 and <= 100, diffusive qflux
 
c---------------------
      if (firsta) then
c---------------------
 
c          Set masks for ocean-ocean boundaries between h-grid boxes.
c          Only allow diffusion between true ocean points (lmaskxy=3),
c          not lakes (lmaskxy=1)
 
        call zero (xmask, nlon*nlat)
        call zero (ymask, nlon*(nlat+1))
        do 50 j=1,nlat
          jm = max (j-1,1)
          do 52 i=1,nlon
            im = mod(i-2+nlon,nlon) + 1
            if (lmaskxy(im,j).eq.3 .and. lmaskxy(i,j).eq.3)
     *        xmask(i,j) = 1.
            if (lmaskxy(i,jm).eq.3 .and. lmaskxy(i,j).eq.3 .and. j.gt.1)
     *        ymask(i,j) = 1.
   52     continue
   50   continue
 
c          Set latitudinal diffusion coefficient coeflat (on "h-grid"
c          boudaries), which depends on:
c            qfactor (user-settable overall factor)
c            zfaco   (depends on zonfrac = zonal fraction of ocean)
c            zfacl   (depends on latitude)
 
        do 60 j=1,nlat
          zonfrac(j) = 0.
          do 62 i=1,nlon
            if (lmaskxy(i,j).eq.3) zonfrac(j) = zonfrac(j) + 1.
   62     continue
          zonfrac(j) = zonfrac(j)/nlon
   60   continue
 
        coeflat(1) = 0.
        coeflat(nlat+1) = 0.
        do 70 j=2,nlat
 
c         Imagine N-S transport is equally efficient if at least some
c         continental barriers, but for ocean fractions near 1,
c         upper-ocean flow is ~zonal with little N-S currents or eddies:
 
          zf = 0.5*(zonfrac(j-1) + zonfrac(j))
          if (zf.le.0.75) then
            zfaco = 1.
          else
            zfaco = 1. - .9  * (zf-.75)/.25
c           zfaco = 1. - .8  * (zf-.75)/.25
c           zfaco = 1. - .75 * (zf-.75)/.25
          endif
 
c         Latitude dependence to have max transport at ~20 deg N/S
 
          zl = abs ( 0.5*(alat(j-1) + alat(j))*180./pi )
          if (zl.le.20.) then
            zfacl = 1.
          else if (zl.le.30.) then
            zfacl = 1. - .7 * (zl-20.)/10.
          else
            zfacl = .3
          endif

c         Uniform coef if qlimit .ne. 0 (namelist). Also see below.

          if (qlimit.ne.0.) then
            zfacl = 1.
            zfaco = 1.
          endif
 
          coeflat(j) = coefstd * qfactor * zfaco * zfacl

c Drake Passage modification # 2
c         if (alat(j).lt. 0.) then
c           coeflat(j) = coefstd * qfactor * zfaco * zfacl ! SH
c         else
c           coeflat(j) = coefstd * 2.0 * zfaco * zfacl     ! NH
c         endif

   70   continue
 
        dlon = 2.*pi/nlon
        dlat = pi/nlat
 
        firsta = .false.
c----------
      endif
c----------
 
      do 100 j=1,nlat
        jm = max (j-1,1)
        jp = min (j+1,nlat)
 
        alm = 0.5 * (alat(j)+alat(jm))
        alp = 0.5 * (alat(j)+alat(jp))
        if (j.eq.1)    alm = -0.5*pi
        if (j.eq.nlat) alp =  0.5*pi
 
        dx  = dlon * cos(alat(j))
        dxm = dlon * cos(alm)
        dxp = dlon * cos(alp)
 
        dy  = dlat
        dyp = dlat
        dym = dlat
 
        slatm = sin(alm)
        slatp = sin(alp)
 
        darea(j) = dlon * (slatp-slatm)
 
c          Reduce longitudinal diffusion for numerical CFL stability
c          near poles (cos(70) = 0.34202)
 
        coefx  = 0.5*(coeflat(j)+coeflat(j+1))
     *           * min ( 1., (cos(alat(j))/0.34202)**2 )
        coefym = coeflat(j)
        coefyp = coeflat(j+1)
 
        do i=1,nlon
          ip = mod(i,nlon) + 1
          im = mod(i-2+nlon,nlon) + 1

          zt1 = (tocnxy(im,j)-tocnxy(i,j)) / dx
          zt2 = (tocnxy(ip,j)-tocnxy(i,j)) / dx
          zt3 = (tocnxy(i,jm)-tocnxy(i,j)) / dym
          zt4 = (tocnxy(i,jp)-tocnxy(i,j)) / dyp
c         Impose limit on abs(dT/dx or dT/dy), deg K per deg arc,
c         if requested (qlimit > 0, namelist, also see above).
c         Convert to deg K per radian.
          if (qlimit.gt.0.) then
            zt1 =  sign ( min(abs(zt1),qlimit*180./pi), zt1 )
            zt2 =  sign ( min(abs(zt2),qlimit*180./pi), zt2 )
            zt3 =  sign ( min(abs(zt3),qlimit*180./pi), zt3 )
            zt4 =  sign ( min(abs(zt4),qlimit*180./pi), zt4 )
          endif
          qfluxxy(i,j) =
     *      (  coefx * xmask(i, j)*zt1 * dy
     *       + coefx * xmask(ip,j)*zt2 * dy
     *       + coefym* ymask(i, j)*zt3 * dxm
     *       + coefyp* ymask(i,jp)*zt4 * dxp
     *      ) / (darea(j) * rearth**2)
        enddo
 
  100 continue
 
c+++++++++++++++++++++++++++++++++
      else if (qfactor.lt.0.) then
c+++++++++++++++++++++++++++++++++
 
c        If qfactor < 0, get qflux from prescribed OHT vs. latitude file
 
c---------------------
      if (firsta) then
c---------------------
 
c          Read alohtin, ohtin (lat in deg, oht in Pw) from input file
 
        write (ioterm2,140)
  140   format(/' Reading ocean heat transport file OHT.DAT')
        iu = 90
        open (iu,file='OHT.DAT',status='old', form='formatted')
        nohtin = 0
        do i=1,999999
          read (iu,*,end=150) alohtin(nohtin+1), ohtin(nohtin+1)
          nohtin = nohtin + 1
        enddo
  150   close(iu)
 
c         Interpolate to model grid's oht(j)
 
        oht(1) = 0.
        oht(nlat+1) = 0.
        do j=2,nlat
          zl = 0.5*(alat(j-1)+alat(j))*180./pi
          if (zl.le.alohtin(1)) then
            oht(j) = ohtin(1)*(zl+90.)/(alohtin(1)+90.)
          else if (zl.ge.alohtin(nohtin)) then
            oht(j) = ohtin(nohtin)*(90.-zl)/(90.-alohtin(nohtin))
          else
            do jj=nohtin-1,1,-1
              if (zl.ge.alohtin(jj)) then
                zw = (zl-alohtin(jj)) / (alohtin(jj+1)-alohtin(jj))
                oht(j) = (1.-zw)*ohtin(jj) + zw*ohtin(jj+1)
                goto 160
              endif
            enddo
  160       continue
          endif
        enddo
 
c         Convert to convergence (from Pw to W/m2)
 
        dlon = 2.*pi/nlon
        dlat = pi/nlat
 
        do j=1,nlat
          zonfrac(j) = 0.
          do i=1,nlon
            if (lmaskxy(i,j).eq.3) zonfrac(j) = zonfrac(j) + 1.
          enddo
          zonfrac(j) = zonfrac(j)/nlon
 
          if (zonfrac(j).eq.0.) then
            oht(j) = 0.
            oht(j+1) = 0.
          endif
 
          darea(j) = dlon*dlat*cos(alat(j))
        enddo
 
        call zero (qfluxoht, nlat)
        do j=1,nlat
          zareastrip = zonfrac(j)*nlon*darea(j)*(rearth**2)
          if (zareastrip.gt.0.) then
            qfluxoht(j) = -(oht(j+1)-oht(j))*1.e15 / zareastrip
          endif
        enddo
 
c       write (6,775) (j, alohtin(j),ohtin(j), j=1,nohtin)
c 775   format('j,alohtin,ohtin=',i3,2f10.3)
c       write (6,*)
c       write (6,776) (j, alat(j)*180./pi,oht(j),qfluxoht(j), j=1,nlat)
c 776   format('j,al,oht,qflux=',i3,3f10.3)
 
        firsta = .false.
c----------
      endif
c----------
 
      call zero (qfluxxy, nlon*nlat)
      do j=1,nlat
        do i=1,nlon
          if (lmaskxy(i,j).eq.3) qfluxxy(i,j) = abs(qfactor)*qfluxoht(j)
        enddo
      enddo
 
c+++++++++++++++++++++++++++++++++++
      else if (qfactor.gt.100.) then
c+++++++++++++++++++++++++++++++++++
 
c        If qfactor > 100, specified global qflux map
 
c---------------------
      if (firsta) then
c---------------------
 
        dlon = 2.*pi/nlon
        dlat = pi/nlat
        do j=1,nlat
          darea(j) = dlon*dlat*cos(alat(j))
        enddo
 
        write(ioterm2,180) 'QFLUX'
  180   format(/' Reading prescribed oceanic q-flux data file ',a)
 
        iu = 90
        open (iu, file='QFLUX', status='old', form='formatted',err=190)
 
        read(iu,'(a8,2i8,8x,a)', err=190) ckey, mlon, mlat, comment
        if (ckey.ne.'QFLUX' .or. mlon.ne.nlon .or. mlat.ne.nlat) then
           write(ioterm2,182) 'QFLUX',  comment(1:lenchr(comment)),
     *                     ckey,mlon,mlat, 'QFLUX',nlon,nlat
  182      format(/' *** Fatal error.',
     *             ' Incorrect attributes for data file ',a
     *            /'     file comment                        = ',a,
     *            /'     file     label, longitude, latitude = ',a8, 2i6
     *            /'     expected label, longitude, latitude = ',a8, 2i6
     *            /)
           call endrun (-1)
        else
           write(ioterm2,184) comment(1:lenchr(comment)), mlon, mlat
  184      format('    file comment   : ',a
     *           /'    file resolution:',2i6)
        endif
c
c          skip 3 records (blank, longitude grid, blank)
c
        read(iu, '(//)', err=190)
c
c          read "f5.0-format" data records
c
        do j=nlat,1,-1
          read (iu,'(5x,3x,360f5.0)',err=190)
     *      (qfluxxyin(i,j),i=nlon/2+1,nlon),
     *      (qfluxxyin(i,j),i=1,nlon/2)
        enddo
 
        close (iu)
 
        goto 194
  190   write (ioterm2,192) 'QFLUX'
  192   format(/'*** Error opening or reading data file ',a)
        call endrun (-1)
  194   continue
 
        do j=1,nlat
          do i=1,nlon
            if (lmaskxy(i,j).ne.3) qfluxxyin(i,j) = 0.
            if (nint(qfluxxyin(i,j)).eq.9999) qfluxxyin(i,j) = 0. !v1.02
          enddo
        enddo
 
        firsta = .false.
c----------
      endif
c----------
 
      do j=1,nlat
        do i=1,nlon
          qfluxxy(i,j) = (qfactor-100.)*qfluxxyin(i,j)
        enddo
      enddo
 
c++++++++++
      endif
c++++++++++
 
c        Adjust towards (user-settable) qice[s,n] for sea-ice fraction
c        (do only if qice[s,n] >= 0)
 
      do 200 j=1,nlat
        if (j.le.nlat/2) then
          qice = qices
        else
          qice = qicen
        endif
        if (qice.ge.0.) then
          do 202 i=1,nlon
            if (lmaskxy(i,j).eq.3)
     *        qfluxxy(i,j) = (1.-ficexy(i,j)) * qfluxxy(i,j)
     *                          +ficexy(i,j)  * qice
  202     continue
        endif
  200 continue
c
c        Adjust for Norwegian-sea lateral+vertical oceanic heat flux,
c        if requested
c
      if (qnorweg.gt.0.) then
        do 220 j=1,nlat
          zlat = alat(j)*180./pi
          do 222 i=1,nlon
            zlon = alon(i)*180./pi
 
            if (lmaskxy(i,j).eq.3) then
              if ( (zlat.ge. 66 .and. zlat.le.78) .and.
     *             (zlon.ge.350 .or.  zlon.le.56)       )  then
                zq= qnorweg*(1.- min(1.,(tocnxy(i,j)-tfreezxy(i,j))/3.))
                qfluxxy(i,j) = max (zq, qfluxxy(i,j))
              endif
            endif
 
  222     continue
  220   continue
      endif
c
c        Calculate total flux qadj, and adjust qfluxxy by an additive
c        correction so that global integral is zero. Only apply this
c        correction between 56S and 56N, so that fluxes under seaice
c        are unaffected. (Don't do this on local machine, which only
c        does 1 point.)
c
      qadj = 0.
      aadj = 0.
      do 300 j=1,nlat
        zlat = alat(j)*180./pi
        do 302 i=1,nlon
          if (lmaskxy(i,j).eq.3) then
             qadj = qadj + qfluxxy(i,j)*darea(j)
             if (abs(zlat).lt.56.) aadj = aadj + darea(j)
          endif
  302   continue
  300 continue
c
c=========================
#ifdef gcm
      if (aadj.gt.0.) then
#else
      if (.false.) then
#endif
c=========================
c
        qadj = qadj/aadj
        do 310 j=1,nlat
          zlat = alat(j)*180./pi
          do 312 i=1,nlon
          if (lmaskxy(i,j).eq.3 .and. abs(zlat).lt.56.)
     *      qfluxxy(i,j)= qfluxxy(i,j) - qadj
  312     continue
  310   continue
 
c          If first call, check adjusted global average = 0,
c          compute zonally integrated northward flux, print map
 
c======================================
        if (showmaps .and. firstb) then
c======================================
 
          zqtot = 0.
          zatot = 0.
          zonf(1) = 0.
          zond(1) = 0.
 
          do 400 j=1,nlat
            zq = 0.
            zt = 0.
            zd = 0.
            nz = 0
            nzd = 0
            do 402 i=1,nlon
              if (lmaskxy(i,j).eq.3) then
                zq = zq + qfluxxy(i,j)
                zt = zt + tocnxy(i,j)
                nz = nz + 1
                if (lmaskxy(i,min(j+1,nlat)).eq.3) then
                  zd = zd + (tocnxy(i,j) - tocnxy(i,min(j+1,nlat)))
                  nzd = nzd + 1
                endif
              endif
  402       continue
            zqtot = zqtot + zq*darea(j)
            zatot = zatot + nz*darea(j)
            zonf(j+1) = zonf(j) - zq*darea(j)*rearth**2
            zonq(j) = zq / max (nz,1)
            zont(j) = zt / max (nz,1)
            if (nz.eq.0) zont(j) = tmelt
            zond(j+1) = zd / max(nzd,1)
  400     continue
 
 
          write(ioterm2,410) qfactor, -qadj, zqtot/max(zatot,epsilon)
  410     format(/' Adjusted oceanic q-flux (first timestep):',
     *           /'    qfactor           =',f13.3,
     *           /'    adjustment (W/m2) =',f13.6,
     *           /'    final mean (W/m2) =',e13.6)
 
          if (qfactor.ge.0. .and. qfactor.le.100.) then
            write(ioterm2,411)
     *        ( alat(min(j,nlat))*180./pi,
     *          zonfrac(min(j,nlat)),
     *          coeflat(min(j,nlat))*1.e-13,
     *          zont(min(j,nlat))-tmelt,
     *          zond(j),
     *          zonq(min(j,nlat)),
     *          zonf(j)*1.e-15, j=nlat+1,1,-1 )
  411       format(
     *        /'    lat   ocn   dcoef    tocn    delt    conv     flux'
     *        /'         frac   10^13     (C)     (C)  (W/m2)  10^15 W'
     *        /'                W/m/K'
     *        /(   f7.2, f6.3,   f8.3,   f8.3,   f8.3,   f8.3,    f9.3))
 
          else if (qfactor.lt.0.) then
            write(ioterm2,412)
     *        ( alat(min(j,nlat))*180./pi,
     *          zonfrac(min(j,nlat)),
     *          qfluxoht(min(j,nlat)),
     *          zonq(min(j,nlat)),
     *          oht(j),
     *          zonf(j)*1.e-15, j=nlat+1,1,-1 )
  412       format(
     *        /'    lat   ocn  convin    conv    ohtin      oht'
     *        /'         frac    W/m2    W/m2   10^15W   10^15W'
     *        /(   f7.2, f6.3,  f8.3,    f8.3,    f9.3,    f9.3))
          endif
 
          do 420 j=1,nlat
            do 422 i=1,nlon
              iworkl(i,j) = nint (qfluxxy(i,j))
              if (lmaskxy(i,j).ne.3) iworkl(i,j)=1000000
  422       continue
  420     continue
          jcm = mod(iabs(mcdate)/100,100)
          jcd = mod(iabs(mcdate),100)
          write (ctmp,425) jcm, jcd
  425     format('Ocean heat flux (W/m2) for date ',i2.2,'/',i2.2)
          call mapout (ioterm2, iworkl, nlon, nlat, alon, alat,
     *                 1000000, ctmp(1:lenchr(ctmp)))
 
c         do 440 j=1,nlat
c           do 442 i=1,nlon
c             iworkl(i,j) = nint (tocnxy(i,j)-273.16)
c             if (lmaskxy(i,j).ne.3) iworkl(i,j)=1000000
c 442       continue
c 440     continue
c         write (ctmp,445) jcm, jcd
c 445     format('Ocean temperature (deg C) for date ',i2.2,'/',i2,2)
c         call mapout (ioterm2, iworkl, nlon, nlat, alon, alat,
c    *                 1000000, ctmp(1:lenchr(ctmp)))
 
          firstb = .false.
 
c============
        endif
c============
      endif
c============
 
      return
      end
