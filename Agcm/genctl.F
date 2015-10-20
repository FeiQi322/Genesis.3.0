c -- $Id: genctl.F,v 1.2 2001/07/03 21:29:18 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/genctl.F,v $
c -- $Name:  $

C
C ***************************** GENCTL ********************************
C
 
c=========================
#if ! defined ( coupogcm )
c=========================

      PROGRAM GENCTL
C ---------------------------------------------------------------------
C
CL            MAIN PROGRAM FOR THE GENESIS GCM
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      LOGICAL IFEXIST
C----------------------------------------------------------------------

 
C        Initialize model
#if defined ( aixsummary )
      ! - initialize for summary call (echos resource use) at end:
c     call start()
#endif

c        Set floating-point exception handler for sun (solaris)

#if defined (sun)
      external hand1, hand2, hand3
      ier = ieee_handler ('set', 'division', hand1)
      ier = ieee_handler ('set', 'overflow', hand2)
      ier = ieee_handler ('set', 'invalid',  hand3)
#endif

c        Delete stub file that will signify normal end (not abend)

      open (188, file='okend', status='unknown')
      close (188, status='delete')

      CALL BASIC
      IER = FLUSH (NOUT)
 
C        Step through the run
 
      DO 100 IS = NRSTRT, NSTOP
        NSTEP = IS
        CALL LINEMS
 
c          Check for external stop signal
 
        inquire (file='STOP', exist=ifexist)
        if (ifexist) then
          write(nout,*) '------------------------------------'
          write(nout,*) '*** Stopping: file STOP detected ***'
          write(nout,*) '------------------------------------'
          open (99, file='STOP', status='old')
          close (99, status='delete')
          call endrun (0) 
        endif
 
  100 CONTINUE
 
#if defined ( aixsummary )
      ! - initialize for summary call (echos resource use) at end:
      ! - call summary to echo resource use:
      call summary()
#endif
 
      WRITE (NOUT,200)
  200 FORMAT(/' --------------------------------'
     *       /' ******* END OF MODEL RUN *******'
     *       /' --------------------------------'/)

c     create stub file to signify normal end of run:
      open (188, file='okend', status='new')
      close (188)

      CALL ENDRUN (0)

      END
 
c=====
#endif
c=====
C
C ***************************** LINEMS ********************************
C
      SUBROUTINE LINEMS
C----------------------------------------------------------------------
C
CL            CONTROL PHYSICS, DYNAMICS, HISTORY, REGEN/RESTART,
CL            THROUGH ONE TIME STEP
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMDRY.h>
C----------------------------------------------------------------------
#include <COMDYR.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMHED.h>
C----------------------------------------------------------------------
#include <COMHIS.h>
C----------------------------------------------------------------------
#include <COMHIX.h>
C----------------------------------------------------------------------
#include <COMHIY.h>
C----------------------------------------------------------------------
#include <COMIMP.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COML30.h>
C----------------------------------------------------------------------
#include <COML50.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMP2.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
#include <COMSDT.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
#include <CRDCTL.h>
C----------------------------------------------------------------------
c     Local arrays to enable use of sigma-hybrid coords
c     (prognostic u,v,t on sigma, prognostic q and tracers on hybrid):
 
c     yu_eta = eastward  velocity on hybrid coords (time n)
c     yv_eta = northward velocity on hybrid coords (time n)
c     yt_eta = temperature on hybrid coords (time n)
c     yq_sig = specific humidity on sigma coords (time n)
c     ypm    = pressures at mid-layer points on hybrid coords (time n)
c     ypl    = pressures at layer interfaces on hybrid coords (time n)
c     ypthicm= layer pressure-thicknesses on hybrid coords (time n-1)
c     ypthic = layer pressure-thicknesses on hybrid coords (time n)
c     ypthicp= layer pressure-thicknesses on hybrid coords (time n+1)
c     ypmcap = ypm**cappa on hybrid coords (time n)
c     yplcap = ypl**cappa on hybrid coords (time n)
c     wslt   = first d[ln(ps)]/dt, then etadot or sigdot (time n)
 
      DIMENSION
     *   YU_ETA(NLON,NOREC,NLEV),   YV_ETA(NLON,NOREC,NLEV),
     *   YT_ETA(NLON,NOREC,NLEV),   YQ_SIG(NLON,NOREC,NLEV),
     *   YPM(NLON,NOREC,NLEV),      YPL(NLON,NOREC,NLEVP),
     *   YPTHICM(NLON,NOREC,NLEV),  YPTHIC(NLON,NOREC,NLEV),
     *   YPTHICP(NLON,NOREC,NLEV),
     *   YPMCAP(NLON,NOREC,NLEV),   YPLCAP(NLON,NOREC,NLEVP),
     *   WSLT(NLON,NOREC,NLEV)
 
c        For slt global mass-adjustment fixer (traglob)
 
      DIMENSION
     *  TROLD(NLON,NOREC,NLEV,NTRACE), YPTHICOLD(NLON,NOREC,NLEV)
#ifdef isotrac
      dimension trold_o(nlon,norec,nlev,3), trold_h(nlon,norec,nlev,3)
#endif

      DIMENSION WORKH(NLON,NOREC)
      CHARACTER CMD*200, CNAME*96, DNAME*96, CONCAT*96
      LOGICAL TESTHIST
      parameter (pi=3.14159265358979)
C----------------------------------------------------------------------
 
C
CL              CALCULATE CURRENT DATE AND SECONDS, CALENDAR VARIABLES
C
      CALL NUTIME (MCDATE, MCSEC, MBDATE, NSTEP, DTIME, 0)
      CALL CALDYI (JYEAR, JMNTH, JDAY, CALDAY, CALSEC, ISTEPY,
     *             MCDATE, MCSEC, MBDATE, DTIME)
C
CL              SET HISTORY FLAGS: FAHIST(ACCUMULATION), FWHIST(WRITE),
CL              FSAVE(SAVE), FNEXF(NEXT FILE). TEST ON ISTEPY, THE
CL              TIMESTEP WITHIN CURRENT 365-DAY YEAR.
C
      FAHIST = .FALSE.
      FWHIST = .FALSE.
      FSAVE = .FALSE.
      FNEXF = .FALSE.
C
      DO 130 I=1,NHISTN
         IF (NHIST(I,1).LE.NHIST(I,2)) THEN
           IF (ISTEPY.GE.NHIST(I,1) .AND. ISTEPY.LE.NHIST(I,2))
     *        FAHIST = .TRUE.
         ELSE
           IF (ISTEPY.GE.NHIST(I,1) .OR.  ISTEPY.LE.NHIST(I,2))
     *        FAHIST = .TRUE.
         ENDIF
         IF (ISTEPY.EQ.NHIST(I,2)) FWHIST = .TRUE.
         IF (HISTLAST .AND. NSTEP.EQ.NSTOP) FWHIST = .TRUE.
  130 CONTINUE
C
      DO 132 I=1,NSAVEN
         IF (ISTEPY.EQ.NSAVE(I)) FSAVE = .TRUE.
  132 CONTINUE
C
      DO 134 I=1,NNEXFN
         IF (ISTEPY.EQ.NNEXF(I)) FNEXF = .TRUE.
  134 CONTINUE
C
CL              LAST TIMESTEP OR NEXT-FILE AUTOMATICALLY FORCES FSAVE
C
      IF (NSTEP.EQ.NSTOP .OR. FNEXF) FSAVE = .TRUE.
C
CL              INITIALIZE VARIOUS GRID ARRAYS FOR SIGETA/SLT SUBRS
C
      IF (NSTEP.EQ.NRSTRT) CALL SLTINI (ALONG,ALATG,SIG,SIGKMH,HYBTYPE)

c               Apply first-timestep perturbation test if requested

      if (nstep.eq.nrstrt .and. perturb.ne.0.) then

#if   defined ( aix )
        call random_seed()
#elif defined ( cray )
        zdum = ranset(12345)
#elif defined ( sgi )
        call srand(12345)
#elif defined ( sun )
        zdum = rand(12345)
#elif defined ( linux )
        zdum = 0
        call random_seed()
#else
        write (*,*) 'Error (linems, ran1): unrecognized machine type'
        call endrun (-1)
#endif
        write (nout,138) perturb
  138   format (/'*** Applying air-temperature perturbations',
     *          ' *(1 +/- ran*', e15.5,' )' )
   
        do jk=1,nlev
          do jj=1,norec
            do ji=1,nlon
#if   defined ( aix )
              call random_number(zran)
#elif defined ( cray )
              zran = ranf()
#elif defined ( sgi )
              zran = rand(0)
#elif defined ( sun )
              zran = rand(0)
#elif defined ( linux )
              call random_number(zran)
#else
              write(*,*) 'Error (linems,ran2):unrecognized machine type'
              call endrun (-1)
#endif
              zz = 1. + perturb*(2.*(zran-0.5))
c             zz = 1. + perturb * sin(2.*3.1415927*(ji+jj+jk)/20.)
              ytm1(ji,jj,jk) = ytm1(ji,jj,jk) * zz
              yt(ji,jj,jk) = yt(ji,jj,jk) * zz
            enddo
          enddo
        enddo

      endif
C
CL              TIME FILTER PART 1. STORE (1-2*EPS)*Y... + EPS*Y...M1
CL              IN TEMPORARY FILE TO BE READ BY TIME FILTER PART 2.
CL              DO BEFORE PHYS, WHICH CHANGES YQM1.
C
      CALL TIMESM (YUM1, YU,
#ifdef isotrac
     *             yqm1_o, yq_o, yqm1_h, yq_h, rtrqm1,rtrq,
#endif
     *             1)

C
CL              SET YPSM1 AND YPS. YPSM1(TIME N-1) IS NEEDED FOR
CL              PHYS (BUDGETS), SLT. YPS (TIME N) IS
CL              NEEDED IN SIGETA, PHYS(VDIF),OMCALC, DYN(STATS), SLT.
C
      DO 140 JJ = 1,NOREC
         DO 142 JI = 1,NLON
            YPS(JI,JJ) = EXP (YLPS(JI,JJ))
            YPSM1(JI,JJ) = EXP (YLPSM1(JI,JJ))
  142    CONTINUE
  140 CONTINUE
 
c       Calculate VOMEGA = d[p]/dt at time n on sigma coords
c       (for etadot and for spectral)
 
      CALL OMCALC
      CALL PUTHIS (VOMEGA, 35)
 
c       Calculate WSLT = d[ln(ps)]/dt at time n on sigma coords
 
      CALL LPSDOT (YU, YV, YD, YDPSL, YDPSM, DSIGMA, WSLT)
 
      call prinverv (wslt,   'd[ln(ps)]/dt * sigma(jk)')
      call prinverv (vomega, 'omega / ps')
 
c       Convert WSLT to d[eta]/dt (if hybrid coords used) or to
c       d[sig]/dt (if sigma coords used). In either case, WSLT is still
c       on sigma grid.
 
      CALL ETADOT (WSLT, VOMEGA, YPS, SIG, HYBTYPE)
      CALL PUTHIS (WSLT, 36)
 
      if (hybtype.eq.0) then
        call prinverv (wslt, 'sigdot')
      else
        call prinverv (wslt, 'etadot')
      endif
 
c       Interpolate U,V,T to eta grid for phys and slt, WSLT to eta grid
c       for slt, Q to sigma grid for vdif and spectral, all at time n.
c       (Y*_ETA,WSLT,YQ_SIG,YPM, YPTHIC* are local arrays in linems.)
c       (If hybrid-slt not used, SIGETA merely copies 1st array to 2nd.)
 
#if defined ( omp )
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared (yu,   yu_eta, yps, hybtype,
!$OMP*         yv,   yv_eta,
!$OMP*         yt,   yt_eta,
!$OMP*         wslt,
!$OMP*         yq,   yq_sig)
!$OMP* private (jj)
!$OMP* schedule(dynamic,2)
#endif
      DO 300 JJ=1,NOREC
        CALL SIGETA (YU,   YU_ETA, JJ, NOREC, YPS, JJ,  1, 1, HYBTYPE)
#ifdef omp
        CALL SIGETA (YV,   YV_ETA, JJ, NOREC, YPS, JJ,  1, 1, HYBTYPE)
        CALL SIGETA (YT,   YT_ETA, JJ, NOREC, YPS, JJ,  1, 1, HYBTYPE)
        CALL SIGETA (WSLT, WSLT,   JJ, NOREC, YPS, JJ,  1, 1, HYBTYPE)
#else
        CALL SIGETA (YV,   YV_ETA, JJ, NOREC, YPS, JJ,  1, 0, HYBTYPE)
        CALL SIGETA (YT,   YT_ETA, JJ, NOREC, YPS, JJ,  1, 0, HYBTYPE)
        CALL SIGETA (WSLT, WSLT,   JJ, NOREC, YPS, JJ,  1, 0, HYBTYPE)
#endif
        CALL SIGETA (YQ,   YQ_SIG, JJ, NOREC, YPS, JJ, -1, 1, HYBTYPE)
  300 CONTINUE
 
c       Calculate mid-pt pressures (YPM), etc, (including layer
c       pressure-thicknesses YPTHIC) on eta grid at time n, for phys
c       and slt. If hybrid coords not used, YPM to YPLCAP are just set
c       to appropriate sigma-grid values.
 
      CALL ETAPRES (YPS, SIG, SIGKMH, YPM, YPL, YPTHIC, YPMCAP, YPLCAP,
     *              CAPPA, HYBTYPE)
 
c     write(*,*)
c     do 8880 k=1,nlev
c       zmin = 1.e20
c       zmax = -1.e20
c       do 8881 j=1,norec
c       do 8881 i=1,nlon
c         zmin = min (zmin, ypm(i,j,k))
c         zmax = max (zmax, ypm(i,j,k))
c8881   continue
c       write(*,8882) k, zmin, zmax, zmax-zmin
c8882   format(' k=',i3,'  zmin=',f15.6,'  zmax=',f15.6,'  delt=',f15.6)
c8880 continue
c     write(*,*)
 
c       Also call stripped-down etapres_a with YPSM1, for YPTHICM
c       at time n-1 (used for budget and slt global conservation)
 
      CALL ETAPRES_A (YPSM1, SIG, SIGKMH, YPTHICM, HYBTYPE)
 
c        Call column and surface "physics"
 
      call pringlob (ypthicm, ypthic, 'prephys')
#ifdef isotrac
      call prout_i (ypthicm, ypthic, 'prePHYS   ')
#endif

      CALL PHYS (YUM1,   YVM1,   YTM1,   YQM1,
     *           YU,     YV,     YT,     YQ,
     *           YU_ETA, YV_ETA, YT_ETA, YQ_SIG,
     *           YLWC_C, YLWC_A, YLWC_S,
     *           FPLUMT, FPLUMQ, FPLUMU, FPLUMV,
     *           TRACER, TRAFLUX,YPM,    YPL,    YPTHIC, YPMCAP,
     *           YFI
#ifdef isotrac
     *          ,yqm1_o, yqm1_h, yq_o, yq_h, ylwc_c_o, ylwc_c_h,
     *           ylwc_a_o, ylwc_a_h, ylwc_s_o, ylwc_s_h,
     *           fplumq_o, fplumq_h,
     *           rtrqm1, rtrq, rtrfq
#endif
     *          )

#ifdef isotrac
      call prout_i (ypthic , ypthic, 'postPHYS  ')
#endif
      call pringlob (ypthicm, ypthic, 'postphys')
C
CL              SAVE PRE-ADVECTION FIELDS FOR LATER GLOBAL MASS
CL              ADJUSTMENT FOR TRACERS (TRAGLOB)
C
      CALL SCOPY (NLON*NOREC*NLEV*NTRACE, TRACER, 1, TROLD,     1)
      CALL SCOPY (NLON*NOREC*NLEV,        YPTHIC, 1, YPTHICOLD, 1)
#ifdef isotrac
      call scopy (nlon*norec*nlev, ylwc_c_o, 1, trold_o(1,1,1,1),1)
      call scopy (nlon*norec*nlev, ylwc_c_h, 1, trold_h(1,1,1,1),1)
      call scopy (nlon*norec*nlev, ylwc_a_o, 1, trold_o(1,1,1,2),1)
      call scopy (nlon*norec*nlev, ylwc_a_h, 1, trold_h(1,1,1,2),1)
      call scopy (nlon*norec*nlev, ylwc_s_o, 1, trold_o(1,1,1,3),1)
      call scopy (nlon*norec*nlev, ylwc_s_h, 1, trold_h(1,1,1,3),1)
#endif
C
CL              TWO DELTA-T LEAPFROG TIMESTEP AND HORIZONTAL DIFFUSION
CL              FROM N-1 TO N+1 (PUTS N+1 VALUES IN COML30 (Y*),
CL              DESTROYS COML10 (Y*M1))
C
      call pringlob (ypthicm, ypthic, 'prespec')

c     zke = 0.
c     zak = 0.
c     zps = 0.
c     zas = 0.
c     do jj=1,norec
c       do ji=1,nlon
c         do jk=1,nlev
c           zwei = sqrt(cs(jj)) * dsigma(jk)
c           zke = zke + (yu(ji,jj,jk)**2 + yv(ji,jj,jk)**2)*zwei
c           zak = zak + zwei
c         enddo
c         zwei = sqrt(cs(jj)) * dsigma(jk)
c         zps = zps + exp(ylps(ji,jj))*zwei
c         zas = zas + zwei
c       enddo
c     enddo
c     write (6,'(/a,i5,2f15.5)')
c    *  'luna: before dyn: nstep, zke, zps=',
c    *   nstep, sqrt(zke/zak), zps/zas

      CALL DYN (YQ_SIG)

c     zke = 0.
c     zak = 0.
c     zps = 0.
c     zas = 0.
c     do jj=1,norec
c       do ji=1,nlon
c         do jk=1,nlev
c           zwei = sqrt(cs(jj)) * dsigma(jk)
c           zke = zke + (yu(ji,jj,jk)**2 + yv(ji,jj,jk)**2)*zwei
c           zak = zak + zwei
c         enddo
c         zwei = sqrt(cs(jj)) * dsigma(jk)
c         zps = zps + exp(ylps(ji,jj))*zwei
c         zas = zas + zwei
c       enddo
c     enddo
c     write (6,'(a,i5,2f15.5)')
c    *  'luna: after  dyn: nstep, zke, zps=',
c    *   nstep, sqrt(zke/zak), zps/zas

c        Impose upper limit (namelist) on air temperatures 

      if (hotlimit.ne.-999.) then
        do k=1,nlev
          do j=1,norec
            do i=1,nlon
              yt(i,j,k) = min (yt(i,j,k), tmelt+hotlimit)
            enddo
          enddo
        enddo
      endif
C
CL              LOCALLY CORRECT SURFACE PRESSURE DUE TO THIS TIMESTEP'S
CL              NET WATER-MASS FLUX, RESET YLPS AND YPS, AND CALCULATE
CL              YPTHICP FOR TIME N+1 (FOR SLT GLOBAL CONSERVATION AND
CL              DRY MASS ADJUST)
C
      DO 400 JJ=1,NOREC
        DO 402 JI=1,NLON
          YPS(JI,JJ) = EXP (YLPS(JI,JJ))
          YPS(JI,JJ) = YPS(JI,JJ)
     *               - (  FRAINC(JI,JJ) + FSNOWC(JI,JJ)
     *                  + FRAINS(JI,JJ) + FSNOWS(JI,JJ) - YQFLX(JI,JJ) )
     *                 * GRAVIT * TWODT
          YLPS(JI,JJ) = LOG (YPS(JI,JJ))
  402   CONTINUE
  400 CONTINUE
 
      CALL ETAPRES_A (YPS, SIG, SIGKMH, YPTHICP, HYBTYPE)
C
CL              SLT ADVECTION FOR YQM1 -> YQ (LEAPFROG, TIME N-1 TO N+1)
CL              AND FOR TRACERS INCLUDING LIQUID WATER AMOUNTS AND
CL              FLUX RESERVOIRS (NON-LEAPFROG, TIME N-1 TO N)
C
C        Convert flux reservoirs from total-for-layer to per-kg
C        (fplum* are equivalenced to tracers). Only advect t-reservoir,
C        since advecting q-reservoir causes blowup, and physical
C        justification for advecting u,v-reservoirs is unclear.
C        (q,u,v reservoirs are tracer #s 5-7).
C
      DO 420 JK=1,NLEV
        DO 422 JJ=1,NOREC
          DO 424 JI=1,NLON
            ZM = YPTHIC(JI,JJ,JK)/GRAVIT
            FPLUMT(JI,JJ,JK) = FPLUMT(JI,JJ,JK) / ZM
c           FPLUMQ(JI,JJ,JK) = FPLUMQ(JI,JJ,JK) / ZM
c           FPLUMU(JI,JJ,JK) = FPLUMU(JI,JJ,JK) / ZM
c           FPLUMV(JI,JJ,JK) = FPLUMV(JI,JJ,JK) / ZM
  424     CONTINUE
  422   CONTINUE
  420 CONTINUE
c
      CALL SLTHYB (W, SIG, SIGKMH, DSIGMA,
     *             GRAVIT, REARTH, TWODT,
     *             YU_ETA, YV_ETA, WSLT, WSLTSAV,
     *             PRIESTQ, STRATTIME, STRATFLAG, SMOLAR,
     *             YPTHICM, YPTHICP,
     *             YQM1, YQ, TRACER, NTRACE, NTRACA,
     *             ALONSAV, ALATSAV, SIGSAV,
     *             TOTQO, TOTQN, TOTQA 
#ifdef isotrac
     *            ,yqm1_o, yq_o, yqm1_h, yq_h, ylwc_c_o, ylwc_c_h,
     *             ylwc_a_o, ylwc_a_h, ylwc_s_o, ylwc_s_h,
     *             rtrqm1, rtrq, nrtr
#endif
     *            )
C
C        Convert flux reservoirs back to layer totals
C
      DO 430 JK=1,NLEV
        DO 432 JJ=1,NOREC
          DO 434 JI=1,NLON
            ZM = YPTHICP(JI,JJ,JK)/GRAVIT
            FPLUMT(JI,JJ,JK) = FPLUMT(JI,JJ,JK) * ZM
c           FPLUMQ(JI,JJ,JK) = FPLUMQ(JI,JJ,JK) * ZM
c           FPLUMU(JI,JJ,JK) = FPLUMU(JI,JJ,JK) * ZM
c           FPLUMV(JI,JJ,JK) = FPLUMV(JI,JJ,JK) * ZM
  434     CONTINUE
  432   CONTINUE
  430 CONTINUE
 
c     write variables needed by slt to a file for offline driver
c     call sltoff (nstep, dtime, mcdate, calday, hybtype,
c    *             along, alatg, cosbud, sig, sigkmh, dsigma,
c    *             gravit, rearth, twodt,
c    *             yu_eta, yv_eta, wslt, wsltsav,
c    *             yps, ypm, ypl, ypthicm, ypthic, ypthicp,
c    *             albsud, albsui, solsud, solsui, cosq,
c    *             yt_eta, yqm1, yq,
c    *             ylwc_c, ylwc_a, ylwc_s)
C
C        Save slt advection rate to history file
C
      IF (TESTHIST(53,1)) THEN
        DO 500 JK = 1,NLEV
          DO 502 JJ = 1,NOREC
            DO 504 JI = 1,NLON
              HISTMP(JI,JJ,JK) = (YQ(JI,JJ,JK)-YQM1(JI,JJ,JK)) / TWODT
  504       CONTINUE
  502     CONTINUE
  500   CONTINUE
        CALL PUTHIS (HISTMP, 53)
      ENDIF
C
C        Increment budget slt conservative adjustment for q
C
      DWFIX  = DWFIX  + 0.5*   (TOTQA-TOTQN)
      DWFIXA = DWFIXA + 0.5*ABS(TOTQA-TOTQN)
c     write(*,520) nstep, totqo, totqn, totqa, (totqa-totqn)/totqn
c 520 format('nstep=',i8, '  totqo=',f14.8,'  totqn=',f14.8,
c    *       '  totqa=',f14.8,'  frac=',f14.8)
 
C
CL            DO DRY MASS ADJUSTMENT FOR TIME N+1 (RESETS YPS AND YLPS)
C
C     Compute global atmospheric mass (tmass,qmass are local vars)
      TMASS  = 0.
      DO 800 JJ=1,NOREC
        ZWEI = COSBUD(JJ)/NLON
        DO 802 JI=1,NLON
          TMASS  = TMASS  + ZWEI*YPS(JI,JJ)/GRAVIT
  802   CONTINUE
  800 CONTINUE
 
C     Compute global mass of water vapor (including reservoir fplumq)
c     and cloud liquid water
      QMASS = 0.
      DO 810 JK=1,NLEV
        DO 812 JJ=1,NOREC
          ZWEI  = COSBUD(JJ)/NLON
          DO 814 JI=1,NLON
            ZH2O = YQ(JI,JJ,JK)
     *           + YLWC_C(JI,JJ,JK) + YLWC_A(JI,JJ,JK) +YLWC_S(JI,JJ,JK)
            QMASS = QMASS + ZWEI*(YPTHICP(JI,JJ,JK)/GRAVIT)*ZH2O
     *                    + ZWEI*FPLUMQ(JI,JJ,JK)
  814     CONTINUE
  812   CONTINUE
  810 CONTINUE
 
C     Impose a globally uniform surface pressure correction to restore
C     dry mass to drymass (in Pa, set in topini), and reset ylps,ypthicp
c     for time n+1.
      ZADD = DRYMASS - (TMASS-QMASS)*GRAVIT
      DO 820 JJ=1,NOREC
        DO 822 JI=1,NLON
          YPS(JI,JJ) = YPS(JI,JJ) + ZADD
          YLPS(JI,JJ) = LOG (YPS(JI,JJ))
c         To avoid roundoff error in global mass fixer for tracers:
          YPS(JI,JJ) = EXP (YLPS(JI,JJ))
  822   CONTINUE
  820 CONTINUE
      CALL ETAPRES_A (YPS, SIG, SIGKMH, YPTHICP, HYBTYPE)
 
C     Accumulate dry-mass correction (in units of Pa) for budget
      DMASS  = DMASS  + 0.5*ZADD
      DMASSA = DMASSA + 0.5*ABS(ZADD)
c     write(6,901) nstep, tmass*gravit, qmass*gravit, zadd
c 901 format('nstep=',i4,' tmass=',f16.6,' qmass=',f16.6,' zadd=',f16.6)
 
C
CL              DO SECOND PART OF TIME FILTER, AND STORE COMPLETED
CL              FILTER FOR TIME N IN COML10
C
      CALL TIMESM (YUM1, YU,
#ifdef isotrac
     *             yqm1_o, yq_o, yqm1_h, yq_h, rtrqm1,rtrq,
#endif
     *             2)
C
CL              RESET YPSM1,YPTHIC FOR TIME N (FOR PUTHIS, BUDGET)
C
      DO 940 JJ = 1,NOREC
         DO 942 JI = 1,NLON
            YPSM1(JI,JJ) = EXP (YLPSM1(JI,JJ))
  942    CONTINUE
  940 CONTINUE
      CALL ETAPRES_A (YPSM1, SIG, SIGKMH, YPTHIC, HYBTYPE)
 
      call pringlob (ypthic, ypthicp, 'postfilt')
C
CL              GLOBAL MASS ADJUSTMENT FOR TRACERS, "+/-" METHOD.
CL              (1 TO 3 ARE CLOUD LIQUID WATER, 5-7 ARE FLUX RESERVOIRS
CL              FOR T,Q,U,V; NTRACA+1 TO NTRACE ARE PASSIVE
CL              TRACERS). DON'T DO FOR Q,U,V RESERVOIRS (SINCE NOT
CL              TRANSPORTED...SEE CALL TO SLTHYB ABOVE).
CL              FOR EXACT CONSERVATION, RELATE TRACER TO YPTHICP (YLPS)
CL              AT TIME N+1, WHICH WILL BE YPTHIC (YLPS) AT TIME N FOR
CL              THE NEXT TIMESTEP (UNTIL DYN),COPIED ABOVE TO YPTHICOLD.
C
      DO N=1,NTRACE
        IF (N.LE.4 .OR. N.GE.NTRACA+1)
     *    CALL TRAGLOB (TRACER(1,1,1,N), TROLD(1,1,1,N),
     *                  YPTHICP, YPTHICOLD,
     *                  TOTTO(N), TOTTN(N), TOTTA(N), W, GRAVIT, N)
      ENDDO
c     nz = min (ntraca+1, ntrace)
c     write(*,950) nstep, totto(nz), tottn(nz), totta(nz),
c    *             (totta(nz)-tottn(nz))/tottn(nz)
c 950 format('nstep=',i8,'  totto=',f14.8,'  tottn=',f14.8,
c    *       '  totta=',f14.8,'  frac=',f14.8)

#ifdef isotrac
c       Global mass adjustment for isotope cloud liquid water by
c       "+/-" method (traglob) and ratios extrema fixing (traglo_i)
      call traglob (ylwc_c_o, trold_o(1,1,1,1),
     *              ypthicp, ypthicold, dum1, dum2, tot, w, gravit, 0)
      call traglo_i(ylwc_c, trold(1,1,1,1), ylwc_c_o, trold_o(1,1,1,1),
     *              ypthicp, w, gravit, tot)
      call traglob (ylwc_c_h, trold_h(1,1,1,1),
     *              ypthicp, ypthicold, dum1, dum2, tot, w, gravit, 0)
      call traglo_i(ylwc_c, trold(1,1,1,1), ylwc_c_h, trold_h(1,1,1,1),
     *              ypthicp, w, gravit, tot)
      call traglob (ylwc_a_o, trold_o(1,1,1,2),
     *              ypthicp, ypthicold, dum1, dum2, tot, w, gravit, 0)
      call traglo_i(ylwc_a, trold(1,1,1,2), ylwc_a_o, trold_o(1,1,1,2),
     *              ypthicp, w, gravit, tot)
      call traglob (ylwc_a_h, trold_h(1,1,1,2),
     *              ypthicp, ypthicold, dum1, dum2, tot, w, gravit, 0)
      call traglo_i(ylwc_a, trold(1,1,1,2), ylwc_a_h, trold_h(1,1,1,2),
     *              ypthicp, w, gravit, tot)
      call traglob (ylwc_s_o, trold_o(1,1,1,3),
     *              ypthicp, ypthicold, dum1, dum2, tot, w, gravit, 0)
      call traglo_i(ylwc_s, trold(1,1,1,3), ylwc_s_o, trold_o(1,1,1,3),
     *              ypthicp, w, gravit, tot)
      call traglob (ylwc_s_h, trold_h(1,1,1,3),
     *              ypthicp, ypthicold, dum1, dum2, tot, w, gravit, 0)
      call traglo_i(ylwc_s, trold(1,1,1,3), ylwc_s_h, trold_h(1,1,1,3),
     *              ypthicp, w, gravit, tot)
      call prout_i (ypthic, ypthicp, 'FINAL     ')
#endif
C
CL              SAVE SURFACE TOPOGRAPHY: TRUNCATED (YFI/GRAVIT, COML10),
CL              AND DELTA TRUNC_MINUS_UNTRUNC (USING TOPOGUN, COMMUN)
C
      IF (TESTHIST(1,1) .OR. TESTHIST(10,1)) THEN
        DO 990 JJ=1,NOREC
          DO 992 JI=1,NLON
            HISTMP(JI,JJ,1) = YFI(JI,JJ)/GRAVIT
            HISTMP(JI,JJ,2) = TOPOGUN(JI,JJ)
  992     CONTINUE
  990   CONTINUE
        CALL PUTHIS (HISTMP(1,1,1),  1)
        CALL PUTHIS (HISTMP(1,1,2), 10)
      ENDIF
C
CL              SAVE 200,500,850 mb GEOPOTENTIALS FOR TIME N (COML10)
C
      IF (TESTHIST(20,1) .OR. TESTHIST(21,1) .OR. TESTHIST(22,1)) THEN
 
c         Compute 3-D geopotential on sigma surfaces (see Mod. Proc.
c         User's Guide pg.I-72, and CCM1 Descrip.eq.(3.b.11). Matrix G
c         in commap is set in setcom, already multiplied by RAIR.
 
        DO 1000 JK=1,NLEV
 
          DO 1002 JJ=1,NOREC
            DO 1004 JI=1,NLON
              HISTMP(JI,JJ,JK) = YFI(JI,JJ)/GRAVIT
 1004       CONTINUE
 1002     CONTINUE
 
          KOFF = NLEV*(JK-1)
          DO 1006 JL=1,NLEV
            ZHYD = G(JL+KOFF) / GRAVIT
            DO 1007 JJ=1,NOREC
              DO 1008 JI=1,NLON
c               Use t on sigma, q on eta coords (spatial error!)
                ZTV = YTM1(JI,JJ,JL)*(1+ZVIR*YQM1(JI,JJ,JL))
c               Or use unfiltered yq_sig (time error!)
c               ZTV = YTM1(JI,JJ,JL)*(1+ZVIR*YQ_SIG(JI,JJ,JL))
                HISTMP(JI,JJ,JK) = HISTMP(JI,JJ,JK) + ZHYD*ZTV
 1008         CONTINUE
 1007       CONTINUE
 1006     CONTINUE
 
 1000   CONTINUE
 
c          Vertically interpolate to requested pressure surfaces
c          and write to history file
 
        DO 1010 N=20,22
          IF (TESTHIST(N,1)) THEN
            IF (N.EQ.20) ZGEO = 200*1.E2
            IF (N.EQ.21) ZGEO = 500*1.E2
            IF (N.EQ.22) ZGEO = 850*1.E2
 
            ZFAC = LOG(SIG(NLEV)) / LOG(SIG(NLEV-1)/SIG(NLEV))
            ZLAPSE = 6.5E-3
c           ZLAPSE = 10.E-3
            ZPOW = GRAVIT / (RAIR*ZLAPSE)
 
            DO 1012 JJ=1,NOREC
              DO 1014 JI=1,NLON
 
                ZSIG = ZGEO/YPSM1(JI,JJ)
                DO 1016 JK=2,NLEV
                  IF (SIG(JK).GE.ZSIG) THEN
                    JKM = JK-1
                    JKP = JK
                    ZWEI = (SIG(JKP)-ZSIG) / (SIG(JKP)-SIG(JKM))
                    WORKH(JI,JJ) =     ZWEI *HISTMP(JI,JJ,JKM)
     *                           + (1.-ZWEI)*HISTMP(JI,JJ,JKP)
                    GOTO 1014
                  ENDIF
 1016           CONTINUE
 
c               Requested pressure level is below lowest model level,
c               so interpolate between there and sea-level-pressure
c               sigma value (where geopotential=0). For slp, use same
c               extrapolation formula as used for slp history field.
                ZTS = YTM1(JI,JJ,NLEV)
     *              + ZFAC*(YTM1(JI,JJ,NLEV)-YTM1(JI,JJ,NLEV-1))
                ZTS = MAX (YTM1(JI,JJ,NLEV-1), ZTS)
                ZSLS = (1. + ZLAPSE*YFI(JI,JJ)/(ZTS*GRAVIT))**ZPOW
 
                ZWEI = (ZSLS-ZSIG) / (ZSLS-SIG(NLEV))
                WORKH(JI,JJ) = ZWEI * HISTMP(JI,JJ,NLEV)
 
 1014         CONTINUE
 1012       CONTINUE
 
            CALL PUTHIS (WORKH, N)
          ENDIF
 1010   CONTINUE
 
      ENDIF
C
CL              SAVE FILTERED FIELDS FOR TIME N (COML10)
C
      CALL PUTHIS (YPSM1, 2)
      CALL PUTHIS (YTM1, 3)
      CALL PUTHIS (YTM1, 66)
      CALL PUTHIS (YTM1(1,1,NLEV), 43)
      CALL PUTHIS (YUM1, 4)
      CALL PUTHIS (YVM1, 5)
      CALL PUTHIS (YQM1, 6)
#ifdef isotrac
      call puthis (yqm1_o, 69)
      call puthis (yqm1_h, 70)
#ifdef isotracreg
      do it=1,nrtr
        call puthis (rtrq(1,1,1,it), 100 + it)
      enddo
#endif
#endif
C
CL              SAVE 3-D WIND SPEED FOR TIME N 
C
      IF (TESTHIST(60,1)) THEN
        DO 1020 JJ=1,NOREC
          DO 1022 JI=1,NLON
            HISTMP(JI,JJ,1) = SQRT (  YUM1(JI,JJ,NLEV)**2
     *                              + YVM1(JI,JJ,NLEV)**2)
 1022     CONTINUE
 1020   CONTINUE
        CALL PUTHIS (HISTMP, 60)
      ENDIF
C
CL              SAVE VARIOUS 3-D DYNAMICAL FLUXES FOR TIME N 
CL              (U*V, W*U, W*T, W*VORTICITY, W*POTENTIAL VORTICITY)
C
      if (testhist(90,1)) then
        do jk=1,nlev
          do jj=1,norec
            do ji=1,nlon
              histmp(ji,jj,jk) = yum1(ji,jj,jk) * yvm1(ji,jj,jk)
            enddo
          enddo
        enddo
        call puthis (histmp, 90)
      endif

      if (testhist(91,1)) then
        do jk=1,nlev
          do jj=1,norec
            do ji=1,nlon
              histmp(ji,jj,jk) = vomega(ji,jj,jk) * yum1(ji,jj,jk)
            enddo
          enddo
        enddo
        call puthis (histmp, 91)
      endif

      if (testhist(92,1)) then
        do jk=1,nlev
          do jj=1,norec
            do ji=1,nlon
              histmp(ji,jj,jk) = vomega(ji,jj,jk) * ytm1(ji,jj,jk)
            enddo
          enddo
        enddo
        call puthis (histmp, 92)
      endif

      if (testhist(93,1) .or. testhist(94,1)) then
        do jk=1,nlev
          jkp = min (jk+1,nlev)
          jkm = max (jk-1,1)
          do jj=1,norec
            jjp = min (jj+1,norec)
            jjm = max (jj-1,1)
            zslat  = sin(alatg(jj))
            zslatp = sin(alatg(jjp))
            zslatm = sin(alatg(jjm))
            zclat  = cos(alatg(jj))
            zclatp = cos(alatg(jjp))
            zclatm = cos(alatg(jjm))
            zdx = zclat*rearth*(2.*pi/nlon)

            do ji=1,nlon
              jip = mod (ji,       nlon) + 1
              jim = mod (ji+nlon-2,nlon) + 1

c             vomega * vorticity (planetary + relative):
              zvo = 2.*omega*zslat
     *            + (yvm1(jip,jj,jk) - yvm1(jim,jj,jk)) 
     *              / (2.*zdx)
     *            - (yum1(ji,jjp,jk)*zclatp - yum1(ji,jjm,jk)*zclatm)
     *              / (rearth*(zslatp-zslatm))
              histmp(ji,jj,jk) = vomega(ji,jj,jk) * zvo

              if (testhist(94,1)) then
c               vomega * potential vorticity (=g*vorticity*|dtheta/dp|):
                zpm    = sig(jkm)*ypsm1(ji,jj)
                zpp    = sig(jkp)*ypsm1(ji,jj)
                zthetm = ytm1(ji,jj,jkm) * ((1.e5/zpm)**(rair/cpair))
                zthetp = ytm1(ji,jj,jkp) * ((1.e5/zpp)**(rair/cpair))
                zdtdp = abs ((zthetp-zthetm)/(zpp-zpm))
                histmp2(ji,jj,jk) = vomega(ji,jj,jk) * gravit*zvo*zdtdp
              endif

            enddo
          enddo
        enddo
        call puthis (histmp, 93)
        if (testhist(94,1)) call puthis (histmp2, 94)
      endif
C
CL              COMPUTE AND SAVE HISTORY SEA-LEVEL PRESSURE
CL              (SEE MODULAR PROCESSOR USER'S GUIDE PROC02, PG.I-76).
C
      IF (TESTHIST(58,1)) THEN
        ZFAC = LOG(SIG(NLEV)) / LOG(SIG(NLEV-1)/SIG(NLEV))
        ZLAPSE = 6.5E-3
c       ZLAPSE = 10.E-3
        ZPOW = GRAVIT / (RAIR*ZLAPSE)
        DO 1030 JJ=1,NOREC
           DO 1032 JI=1,NLON
              ZTS = YTM1(JI,JJ,NLEV)
     *            + ZFAC*(YTM1(JI,JJ,NLEV)-YTM1(JI,JJ,NLEV-1))
              ZTS = MAX (YTM1(JI,JJ,NLEV-1), ZTS)
              HISTMP(JI,JJ,1)= YPSM1(JI,JJ)
     *                     * (1. + ZLAPSE*YFI(JI,JJ)/(ZTS*GRAVIT))**ZPOW
 1032     CONTINUE
 1030   CONTINUE
        CALL PUTHIS (HISTMP, 58)
      ENDIF
C
CL              SAVE 3-D AND VERTICALLY INTEGRATED TRACER FIELDS
CL              (NOT FOR #1-3, WHICH ARE LIQUID WATER CONTENTS,
CL              OR #4-7, WHICH ARE FLUX RESERVOIRS;
CL              NTRACA=7 SET IN COMMONS).
C
c   3-D tracer history: 
c   cosmogenic isotopes history saves done in subr cosmog
#ifndef cosmog
      IF (NTRACA.LT.NTRACE) THEN
        DO 1040 N=1,MIN(NTRACE-NTRACA,40)
          CALL PUTHIS (TRACER(1,1,1,NTRACA+N), IHTRACE+N-1)
          IF (TESTHIST(IHTRACE+40+N-1, 1)) THEN
            CALL ZERO (HISTMP(1,1,1), NLON*NOREC)
            DO 1042 JK=1,NLEV
              DO 1044 JJ=1,NOREC
                DO 1046 JI=1,NLON
                  HISTMP(JI,JJ,1)= HISTMP(JI,JJ,1)
     *              + TRACER(JI,JJ,JK,NTRACA+N)*YPTHICP(JI,JJ,JK)/GRAVIT
 1046           CONTINUE
 1044         CONTINUE
 1042       CONTINUE
            CALL PUTHIS (HISTMP(1,1,1), IHTRACE+40+N-1)
          ENDIF
 
 1040   CONTINUE
      ENDIF
#endif
C
CL              COMPUTE MIN AND MAX TEMPERATURES FOR PRINTOUT
C
      CALL THILO
C
CL              INITIALIAZE BUDGET IF PREVIOUS BUDGET NOT AVAILABLE
CL              FROM RESTART FILE (EG, COLD/INITIAL STARTS, EARLIER
CL              RESTART-FILE VERSION, CHANGED LATITUDE RESOL,...),
CL              OR FOR RESTARTS WITH A BUDGET DONE AT END OF PREVIOUS
CL              RUN (IN CASE MODEL CHANGED BY NAMELIST).
CL              DO ALSO AFTER SECOND TIMESTEP OF ANY RUN IN WHICH THE
CL              TOTAL ATMOS DRY MASS MAY HAVE BEEN CHANGED FROM THE
CL              PREVIOUS RUN DUE TO DIFFERENT TOPOGRAPHY IN TOPINI
CL              (IE, AFTER BOTH LEAPFROG YPS's HAVE HAD DRY-MASS
CL              CORRECTION).
 
CL              ELSE COMPUTE CURRENT BUDGET, PRINTOUT AND RESET, EVERY
CL              BUDFREQ SECS (FOR FILTERED TIME N) OR LAST STEP OF RUN.
C
      IF( (NSTEP.EQ.NRSTRT .AND. (NBUDGPRE.EQ.0 .OR.
     *                            ABS((NSTEP-1)*DTIME-BTIME0).LT.DTTOL))
     *    .OR.
     *    (NSTEP.EQ.NRSTRT+1 .AND. ABS(DRYMASS-DRYMASSIN).GT.0.1)
     *  ) THEN
c       Compute current state
        CALL BUDCLC (NSTEP*DTIME, YPTHIC, YPTHICP)
c       Zero accumulators and set old budget state vars to current
        CALL BUDINA
      ELSE
        N =  MAX ( NINT(BUDFREQ/DTIME), 1 )
        IF (MOD(NSTEP,N).EQ.0 .OR. NSTEP.EQ.NSTOP) THEN
c         Compute current state
          CALL BUDCLC (NSTEP*DTIME, YPTHIC, YPTHICP)
c         Print budget
          CALL BUDPRN (NOUT)
c         Zero accumulators and set old budget state vars to current
          CALL BUDINA
        ENDIF
      ENDIF
C
CL              CALL USER-DEFINED SUBROUTINE(S)
C
C     CALL USERDEF0 (NSTEP, NRSTRT, NSTOP, ISTEPY, DTIME, NUSER0)
C     CALL USERDEF1 (NSTEP, NRSTRT, NSTOP, ISTEPY, DTIME, NUSER1)
C     CALL USERDEF2 (NSTEP, NUSER2)
C
CL              PRINT OUT CURRENT TIME
C
      CALL SHOWTIME
C
CL              WRITE AGCM FIELDS TO REGIONAL-MODEL FILES
C
      CALL REGWRITE
C
CL              WRITE AGCM FIELDS TO STORM-TRACK HISTORY FILES
C
      CALL STORMWRITE
C
CL              PROCESSING FOR HISTORY AND RESTART FILES.
C
CL              INITIALIZE NSTEPH TO "UNSET", THEN SET NSTEPH TO MARK
CL              BEGINNING OF CURRENT HISTORY ACCUMULATION PERIOD
C
      IF (NSTEP.EQ.NRSTRT) NSTEPH = -1
      IF (FAHIST.AND.NSTEPH.EQ.-1) NSTEPH = NSTEP
C
CL              INCREMENT HISTORY WRITE COUNTER, WRITE TO HISTORY FILE,
CL              RESET RECORD NUMBER OFFSET, AND UNSET NSTEPH, IF ANY
CL              DATA RECORDS WRITTEN THIS RUN.
C
      IF (FWHIST .AND. NULDAT.GT.0) THEN
         CALL WRITH
         CALL HISSET (0)
         NSTEPH = -1
      ENDIF
C
CL              SAVE AND NEXT-FILE FOR HISTORY AND RESTART FILES.
C
C--------------------
      IF (FSAVE) THEN
C--------------------
C
CL              FLUSH OUTPUT BUFFER (BEFORE LARGE FILE-WRITE DELAYS)
C
         IER = FLUSH (NOUT)
C
C               SAVE HISTORY FILE TO MSS, IF NON-EMPTY AND ANY DATA
CL              RECORDS REQUESTED THIS RUN.
C
         IF (NUOFF.GT.0 .AND. NULDAT.GT.0) THEN
C
            CALL CLOSEF (NDATA)
            CALL DISPOSEF (MSNAMHIS(IFILE), MSPATH)
            IF (NSTEP.LT.NSTOP)
     *        CALL OPENF (NDATA, MSNAMHIS(IFILE), MSPATH, -1, 0)

            cname = msnamhis(ifile)
            if (hisnet) then
c             append '.nc' to name if not already there
              lnl = lenchr(cname)
              if (cname(lnl-2:lnl).ne.'.nc' .and.
     *            cname(lnl-2:lnl).ne.'.NC')  cname(lnl+1:) = '.nc'
            endif

            IF (.NOT.HISNET) WRITE (NOUT,1450)
            IF (     HISNET) WRITE (NOUT,1451)
            WRITE (NOUT,1452)
     *         NSTEP, ISTEPY*DTIME/86400., NSTEP*DTIME/86400.,
     *         MSPATH(1:LENCHR(MSPATH)), CNAME(1:LENCHR(CNAME)), NWRIT
 1450       FORMAT(/' ----------------------------------------------'
     *             /' ********** AGCM HISTORY FILE OUTPUT **********')
 1451       FORMAT(/' ----------------------------------------------'
     *             /' ***** AGCM HISTORY FILE OUTPUT (NetCDF) ******')
 1452       FORMAT( ' ----------------------------------------------'
     *             /' timestep                          = ',i10
     *             /' calendar days                     = ',f10.3
     *             /' elapsed days                      = ',f10.3
     *             /' agcm history file                 = ',a,a
     *             /' agcm history writes               = ',i10
     *             /' ----------------------------------------------')
#ifdef mss
            IF (MSPATH.EQ.' ') WRITE (NOUT,1453)
 1453       FORMAT( ' *  File NOT saved to MSS (MSPATH is blank!)  *'
     *             /' ----------------------------------------------')
#endif
            WRITE (NOUT,*)
C
         ENDIF
C
C               IF NEXT-FILE, START NEW HISTORY FILE AND INCREMENT IFILE
C               NB: IF LAST TIMESTEP IS A NEXT-FILE, IFILEPRE POINTS TO
C               THE FILES JUST SAVED, IFILE POINTS TO NEXT FILES, AND
C               NUOFF IS 0. THESE VARIABLES ARE WRITTEN TO RESTART FILE.
C
         IFILEPRE = IFILE
C
         IFENOUGH = 1
         IF (FNEXF) THEN
            IFILE = IFILE + 1
            CALL HISSET (1)
C           Open new local history file, unless last timestep
            IF (NSTEP.LT.NSTOP) THEN
               CALL CLOSEF (NDATA)
               CALL OPENF (NDATA, MSNAMHIS(IFILE), MSPATH, 1, NUSIZ)
c              If "not enough filenames entered" , fatal error, but
c              wait till restart file written below before stopping
               IF (IFILE.GT.NFILE) IFENOUGH = 0
            ENDIF
         ENDIF
C
C               WRITE RESTART FILE AND DISPOSE TO MSS
C
         CALL WSDS (MSNAMRES(IFILEPRE), MSPATH)
C
C               DISPOSE SECONDARY RESTART-FILE COPY IF REQUESTED
C
         IF (MSNAMRE2.NE.' ') THEN
c#####
c Either:
            CALL WSDS (MSNAMRE2, MSPATH)
c#####
c Or:
#ifdef mss
c           cname = msnamres(ifilepre)
c           lenc = lenchr(cname)
c           dname = msnamre2
c           lend = lenchr(dname)
#else
c           cname = concat (mspath, msnamres(ifilepre), lenc)
c           dname = concat (mspath, msnamre2,           lend)
#endif
c           cmd = 'cp '// cname(1:lenc) //' '// dname(1:lend)
c           ier = ishell (cmd(1:lenchr(cmd)))
c           if (ier.ne.0) then
c             write(nout,1458) cname(1:lenc), dname(1:lend), ier
c1458         format
c    *        (/' *** Error in LINEMS (copying secondary restart file):'
c    *         /'     Primary   restart file name  = ',a
c    *         /'     Secondary restart file name  = ',a
c    *         /'     ishell error code            = ',i8)
c             call endrun (-1)
c           endif
c
c           call disposef (dname, mspath)
c#####
         ENDIF
 
         WRITE (NOUT,1460)
     *      NSTEP, ISTEPY*DTIME/86400., NSTEP*DTIME/86400.,
     *      MSPATH(1:LENCHR(MSPATH)), MSNAMRES(IFILEPRE),
     *      VERSRES, IFILEPRE, FNEXF
 1460    FORMAT(
     *          /' ----------------------------------------------'
     *          /' ************* RESTART FILE OUTPUT ************'
     *          /' ----------------------------------------------'
     *          /' timestep                          = ',i10
     *          /' calendar days                     = ',f10.3
     *          /' elapsed days                      = ',f10.3
     *          /' restart file                      = ',a,a
     *          /' restart file version number       = ',f10.2
     *          /' filename list #                   = ',i10
     *          /' filename list # incremented (T/F) = ',l10
     *          /' ----------------------------------------------')
         IF (MSNAMRE2.NE.' ')
     *     WRITE (NOUT,1461) MSPATH(1:LENCHR(MSPATH)), MSNAMRE2
 1461      FORMAT(
     *           ' secondary copy                    = ',a,a
     *          /' ----------------------------------------------')
C
#ifdef mss
         IF (MSPATH.EQ.' ') WRITE (NOUT,1462)
 1462    FORMAT( ' *  File NOT saved to MSS (MSPATH is blank!)  *'
     *          /' ----------------------------------------------')
#endif
C
         WRITE (NOUT,*)
C
C        Stop if "not enough filenames" fatal error detected above
         IF (IFENOUGH.EQ.0) THEN
            WRITE(NOUT,*) '*** Not enough filenames input to',
     *                    ' complete this run. Stop in LINEMS.'
            CALL ENDRUN (-1)
         ENDIF
C
C----------
      ENDIF
C----------
C
CL              FLUSH OUTPUT BUFFER AT END OF EACH DAY
C
      N = MAX (NINT(1.*86400./DTIME), 1)
      IF (MOD(NSTEP,N).EQ.0) IER = FLUSH (NOUT)
C
      RETURN
      END
C
C ***************************** PHYS **********************************
C
      SUBROUTINE PHYS (YUM1,   YVM1,   YTM1,   YQM1,
     *                 YU,     YV,     YT,     YQ,
     *                 YU_ETA, YV_ETA, YT_ETA, YQ_SIG,
     *                 YLWC_C, YLWC_A, YLWC_S,
     *                 FPLUMT, FPLUMQ, FPLUMU, FPLUMV,
     *                 TRACER, TRAFLUX,YPM,    YPL,    YPTHIC, YPMCAP,
     *                 YFI 
#ifdef isotrac
     *                ,yqm1_o, yqm1_h, yq_o, yq_h, ylwc_c_o, ylwc_c_h,
     *                 ylwc_a_o, ylwc_a_h, ylwc_s_o, ylwc_s_h,
     *                 fplumq_o, fplumq_h,
     *                 rtrqm1, rtrq, rtrfq
#endif
     *                )
C----------------------------------------------------------------------
C
CL    DRIVER FOR "VERTICAL PHYSICS" ROUTINES. Y* ARE PASSED FROM
CL    LINEMS (NOT VIA COML30) TO ALLOW USE OF MIXED SIGMA-HYBRID COORDS)
CL
CL      Y[U,V,T]M1 ARE ON SIGMA GRID (TIME N-1),
CL      YQM1 IS ON HYBRID GRID (IF USED) (TIME N-1),
CL      Y[U,V,T], YQ_SIG ARE ON SIGMA GRID (TIME N),
CL      Y[U,V,T]_ETA, YQ, YP* ARE ON HYBRID GRID (IF USED) (TIME N)
CL      YLWC_[C,A,S],TRACER ARE ON HYBRID GRID (IF USED) (NON-LEAPFROG)
CL      FPLUM* ARE ON HYBRID GRID (IF USED) (NON-LEAPFROG)
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMBUD.h>
C---------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMDIF.h>
C---------------------------------------------------------------------
#include <COMDYR.h>
C---------------------------------------------------------------------
#include <COMFLG.h>
C---------------------------------------------------------------------
#include <COMHED.h>
C---------------------------------------------------------------------
#include <COMHIS.h>
C---------------------------------------------------------------------
#include <COMHIX.h>
C---------------------------------------------------------------------
#include <COMHIY.h>
C---------------------------------------------------------------------
#include <COMIMP.h>
C---------------------------------------------------------------------
#include <COMINI.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C---------------------------------------------------------------------
#include <COMMAP.h>
C---------------------------------------------------------------------
#include <COMMP2.h>
C ---------------------------------------------------------------------
#include <COMMUN.h>
C---------------------------------------------------------------------
#include <COMSDT.h>
C---------------------------------------------------------------------
#include <COMTIM.h>
C---------------------------------------------------------------------
#include <CRDCTL.h>
C---------------------------------------------------------------------
#include <LGASES.h>
C---------------------------------------------------------------------
      DIMENSION
     *  YUM1(NLON,NOREC,NLEV),      YVM1(NLON,NOREC,NLEV),
     *  YTM1(NLON,NOREC,NLEV),      YQM1(NLON,NOREC,NLEV),
     *  YU(NLON,NOREC,NLEV),        YV(NLON,NOREC,NLEV),
     *  YT(NLON,NOREC,NLEV),        YQ(NLON,NOREC,NLEV),
     *  YU_ETA(NLON,NOREC,NLEV),    YV_ETA(NLON,NOREC,NLEV),
     *  YT_ETA(NLON,NOREC,NLEV),    YQ_SIG(NLON,NOREC,NLEV),
     *  YLWC_C(NLON,NOREC,NLEV),    YLWC_A(NLON,NOREC,NLEV),
     *  YLWC_S(NLON,NOREC,NLEV),
     *  FPLUMT(NLON,NOREC,NLEV),    FPLUMQ(NLON,NOREC,NLEV),
     *  FPLUMU(NLON,NOREC,NLEV),    FPLUMV(NLON,NOREC,NLEV),
     *  TRACER(NLON,NOREC,NLEV,NTRACE), TRAFLUX(NLON,NOREC,NTRACB),
     *  YPM(NLON,NOREC,NLEV),       YPL(NLON,NOREC,NLEVP),
     *  YPTHIC(NLON,NOREC,NLEV),
     *  YPMCAP(NLON,NOREC,NLEV),
     *  YFI(NLON,NOREC)
#ifdef isotrac
      dimension
     *  yqm1_o(nlon,norec,nlev),    yqm1_h(nlon,norec,nlev),
     *  yq_o(nlon,norec,nlev),      yq_h(nlon,norec,nlev),
     *  ylwc_c_o(nlon,norec,nlev),  ylwc_c_h(nlon,norec,nlev),
     *  ylwc_a_o(nlon,norec,nlev),  ylwc_a_h(nlon,norec,nlev),
     *  ylwc_s_o(nlon,norec,nlev),  ylwc_s_h(nlon,norec,nlev),
     *  fplumq_o(nlon,norec,nlev),  fplumq_h(nlon,norec,nlev),
     *  rtrqm1(nlon,norec,nlev,nrtr), rtrq(nlon,norec,nlev,nrtr),
     *  rtrfq (nlon,norec,nlev,nrtr)
#endif
 
      DIMENSION
     *  YCLD(NLON,NOREC,NLEV),      YCLWP(NLON,NOREC,NLEV),
     *  TOTCLD(NLON,NOREC),         TOTCLDC(NLON,NOREC),
     *  TOTCLDA(NLON,NOREC),        TOTCLDS(NLON,NOREC),
     *  TOTLWC(NLON,NOREC),
     *  ZR(NLON), ZS(NLON), ZL(NLON), ZE(NLON), ZP(NLON),
     *  YRELHUM(NLON,NOREC,NLEV)
 
      DIMENSION
     *  PLUMEPREC(NLON,NOREC,NLEV),
     *  HEIMIN(NLON,NOREC),       HEIMAX(NLON,NOREC)
c     save since (i) convad-plume not called every timestep, and
c     (ii) history pbl heights are local 03:00,15:00 (was daily min/max)
      SAVE PLUMEPREC, HEIMIN, HEIMAX
#ifdef isotrac
      dimension
     *  plumepr_o(nlon,norec,nlev) ,plumepr_h(nlon,norec,nlev),
     *  plumepr_t(nlon,norec,nlev,nrtr)
      save plumepr_o ,plumepr_h, plumepr_t
#endif

#if defined (tropicheat)
      logical firsttrop
      data firsttrop /.true./
      dimension weitrop(nlev)
      save weitrop, firsttrop
#endif
 
      DIMENSION WORKH(NLON,NOREC)
      LOGICAL TESTHIST, FSAVE_LEAF
C---------------------------------------------------------------------
C
CL                 SET FLAGS FOR WHETHER TO DO SOLAR AND IR COMPUTATIONS
CL                 (IRAD*,FRAD* ARE IN CRDCTL)
C
      FRADLW   = MOD(NSTEP-1,IRADLW ).EQ.0 .OR. NSTEP.EQ.NRSTRT
      FRADLWA  = MOD(NSTEP-1,IRADLWA).EQ.0 .OR. NSTEP.EQ.NRSTRT
      FRADLW   = FRADLW .OR. FRADLWA
      FRADSW   = MOD(NSTEP-1,IRADSW ).EQ.0 .OR. NSTEP.EQ.NRSTRT
C
CL                 INITIALIZE VEG/SNOW/SOIL/OCEAN AT (RE)START.
CL                 (MUST CALL INISURF AND VEGDAT BEFORE ALBLSX CALLED
CL                 FROM RADCTL.) ALSO CALL SIZECHK TO CHECK SIZES OF
CL                 HISTORY FILES VS MSS LIMIT.
C
      IF (NSTEP.EQ.NRSTRT) THEN
         CALL INIFUD (NIN, NOUT, TIMOLSX, BUDFREQ, 1)
         CALL INISURF
     *        (NLON, NOREC, ALATG, W, 2, DTIME,
     *         YTSIN,    YOROOCN,  YOROLAND, YOROICE,
     *         DELEVTYPE,YFI,      VEGTYPE,  EVEDELTA, DSBIOME,
     *         WISCVEG,  PHOTOTYPE,TAUBIOME, DYNAIBIS, ASYNIBIS,
     *         SOILTEX,  DEPTHML,  DEPTHLAKE,DEPTHICE,DEPTHSNO,
     *         OCEANTYPE,DYNAMICE, PRESCOUV, SHOWMAPS,
     *         LONOLSX,  LATOLSX,  NSREST,   NSREST_ISOT,
     *         NNEXFN,   NRSTRT,   NSTOP,
     *         NHISI,    NDHIS,    HISTLAST)
         CALL SIZECHK
      ENDIF
C
CL                 SET VEGETATION DATA ON 1st TIMESTEP OF (RE)START AND
CL                 ONCE A DAY THEREAFTER.  MOD(NSTEP,N).EQ.1 COINCIDES
CL                 WITH NSTEP.EQ.NRSTRT AFTER A RESTART IF PREVIOUS RUN
CL                 LENGTH WAS A WHOLE NUMBER OF DAYS, SO THAT RESTART
CL                 OR NOT WON'T AFFECT RESULTS.
CL
CL                 ISTEPY_LEAF, FSAVE_LEAF ARE LOCAL VARS; FSAVE_LEAF
CL                 (= IF GENESIS SAVE THIS TIMESTEP) IS PASSED TO
CL                 VEGDAT_EVE TO CONTROL EVE/LEAF RESTART-FILE USE.
C
      N = NINT(86400./DTIME)
      IF (NSTEP.EQ.NRSTRT .OR. MOD(NSTEP,N).EQ.1) THEN
C
         ISTEPY_LEAF = NINT((CALSEC + 86400. - DTIME)/DTIME)
         FSAVE_LEAF = .FALSE.
         DO 100 I=1,NSAVEN
           IF (ISTEPY_LEAF.EQ.NSAVE(I)) FSAVE_LEAF = .TRUE.
  100    CONTINUE
         DO 110 I=1,NNEXFN
            IF (ISTEPY_LEAF.EQ.NNEXF(I)) FSAVE_LEAF = .TRUE.
  110    CONTINUE
         IF (NSTEP+N.GT.NSTOP) FSAVE_LEAF = .TRUE.
C
         CALL VEGDAT (JYEAR, JMNTH, JDAY, FSAVE_LEAF, TRIASSIC,
     *                NSREST_IBIS, HORRESLSX(1),HORRESLSX(2), TYPERESIN,
     *                SHOWMAPS, CO2PPV*1.E6)
C
      ENDIF
C
CL                 COMPUTE CLOUD AMOUNTS (YCLD,TOTCLD) AND LIQUID WATER
CL                 PATH (YCLWP) FOR RADCTL,CONVAD,SURFCTL.
C
#if defined ( omp )
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared (yt_eta, yq, ylwc_c, ylwc_a, ylwc_s,
!$OMP*         totcld, totcldc, totclda, totclds, totlwc,
!$OMP*         wplume0, wplume1, ypm, ypthic, ycld, yclwp, yrelhum)
!$OMP* private (jj)
!$OMP* schedule(dynamic,2)
#endif
C
      DO 150 JJ=1,NOREC
        CALL CLDCMP (YT_ETA, YQ, YLWC_C, YLWC_A, YLWC_S,
     *               TOTCLD, TOTCLDC, TOTCLDA, TOTCLDS, TOTLWC,
     *               WPLUME0, WPLUME1, YPM, YPTHIC, YCLD, YCLWP,YRELHUM,
     *               JJ)
  150 CONTINUE
C
CL        WRITE CLOUD DIAGNOSTICS
C
      IF (NSTEP.GT.NRSTRT) CALL CLDDIAGO (NOUT)
C
CL        SAVE CLOUD HISTORY FIELDS
 
      CALL PUTHIS (YCLD,    29)
      CALL PUTHIS (YRELHUM, 61)
      CALL PUTHIS (TOTCLD,  44)
      CALL PUTHIS (TOTCLDC, 45)
      CALL PUTHIS (TOTCLDA, 46)
      CALL PUTHIS (TOTCLDS, 47)
      CALL PUTHIS (TOTLWC,  62)
C
CL                 DO RADIATION COMPUTATIONS
C
      if (fradlw .or. fradsw) then
        call radctl (
c            input:
     *       yt_eta, yq, yrelhum, ypm, ypl, ypthic,
     *       ycld, yclwp, yps,
     *       cosq, fraq, albsud, albsui, firems, 
     *       yoroocn, yoroland, yoroice,
c            output:
     *       yqrs, hfrsa, hsabtp, hsabat, hclrss, hclrst, 
     *       hsolin, solsud, solsui, hfrsac, 
     *       yqrl, firabs, sldflx, hfrla, hfirtp, hclrls, hclrlt
#if defined (LKAEROSOL)
     *      ,sfcirrad
#endif
     *       )

      endif
C
CL                 STORE NET RADIATIVE HEATING RATE IN YT2 (HERE IN
CL                 UNITS OF FLUX CONVERGENCE PER LAYER, W M-2...
CL                 CONVERTED TO DEG K/S AT START OF DYN). MUST BE
CL                 DONE EVERY TIMESTEP (NOT IN RADCTL) SINCE YT2
CL                 IS CHANGED IN VDIF AND SPECTRAL SUBRS.
C

#if defined (tropicheat)
c     .005, .013, .033, .064, .099, .139,
c     .189, .251, .325, .409, .501, .598,
c     .695, .787, .866, .929, .970, .993 
      if (firsttrop) then
        zsiga = .15
        zsigb = .5
        zsigc = .75
        zweitot = 0.
        do k=1,nlev
          if (sig(k).lt.zsigb) then
            zz = (sig(k)-zsigb)/(zsigb-zsiga)
          else 
            zz = (sig(k)-zsigb)/(zsigc-zsigb)
          endif
          weitrop(k) = max (0., 1. - zz**2) * dsigma(k)
          zweitot = zweitot + weitrop(k)
        enddo
        do k=1,nlev
          weitrop(k) = weitrop(k) / zweitot
        enddo

        zweitot = 0.
        do k=1,nlev
          write (6,*) 'k, weitrop= ',k, weitrop(k) 
          zweitot = zweitot + weitrop(k)
        enddo
        write (6,*)   '      sum = ',zweitot
      endif

      ztropicheatp = 75.                     ! W/m2 total column
      ztropicheatm = ztropicheatp*90./270.   ! compensate,  0 zonal mean
      zghp = 0.
      zghm = 0.
      zgap = 0.
      zgam = 0.
#endif

      DO 200 JK=1,NLEV
        DO 202 JJ=1,NOREC
          DO 204 JI=1,NLON
            YT2(JI,JJ,JK) = (YQRS(JI,JJ,JK) - YQRS(JI,JJ,JK+1))
     *                    + (YQRL(JI,JJ,JK) - YQRL(JI,JJ,JK+1))
#if defined (tropicheat)
            zwei = 0.5*w(jj)/nlon
            zlat = alatg(jj)* 180./3.1415927
            zlon = along(ji)* 180./3.1415927
            if (zlon.gt.180.) zlon = zlon - 360.

            if (zlat.ge.-10.   .and. zlat.le.10.) then 
c             if (zlon.ge. 170.  .or.  zlon.le.-100.) then         ! 90d
              if (zlon.ge. -90.  .and. zlon.le.   0.) then         ! 90e
                yt2(ji,jj,jk) = yt2(ji,jj,jk) + ztropicheatp*weitrop(jk)
                zghp = zghp + ztropicheatp*weitrop(jk)*zwei
                zgap = zgap +              weitrop(jk)*zwei
              else 
                yt2(ji,jj,jk) = yt2(ji,jj,jk) - ztropicheatm*weitrop(jk)
                zghm = zghm - ztropicheatm*weitrop(jk)*zwei
                zgam = zgam +              weitrop(jk)*zwei
              endif
            endif
#endif
  204     CONTINUE
  202   CONTINUE
  200 CONTINUE

#if defined (tropicheat) 
      if (firsttrop) then
        write (nout,'(3(a,f9.5))')
     *    ' tropicheat: zghp=',zghp/zgap,
     *    ' zghm=',zghm/zgam,
     *    ' zght=',(zghp+zghm)/(zgap+zgam)
        ier = flush(nout)
        firsttrop = .false.
      endif
#endif

#ifdef isotrac
c     Zero isotopic liquid water contents corresponding to liquid water
c     contents zeroed above in CLDCMP when lower than 1e-10
      do jk=1,nlev
        do jj=1,norec
          do ji=1,nlon
            if (ylwc_c(ji,jj,jk).eq.0.) then
              ylwc_c_o(ji,jj,jk) =  0.
              ylwc_c_h(ji,jj,jk) =  0.
            endif
            if (ylwc_a(ji,jj,jk).eq.0.) then
              ylwc_a_o(ji,jj,jk) =  0.
              ylwc_a_h(ji,jj,jk) =  0.
            endif
            if (ylwc_s(ji,jj,jk).eq.0.) then
              ylwc_s_o(ji,jj,jk) =  0.
              ylwc_s_h(ji,jj,jk) =  0.
            endif
          enddo
        enddo
      enddo
#endif
C
CL                 STEP VEG/SNOW/SOIL/OCEAN MODELS THRU ONE AGCM STEP
CL                 (USING LAST TIMESTEP'S PRECIP)
C
      IACCUM = 0
      CALL SURFCTL
     *   (NSTEP,    DTIME,    FRADSW,
     *    COSQ,     SOLSUD,   SOLSUI,   FIRABS,
     *    YU_ETA(1,1,NLEV),   YV_ETA(1,1,NLEV),
     *    YT_ETA(1,1,NLEV),   YQ(1,1,NLEV),
     *    FRAINC,   FSNOWC,   FRAINS,   FSNOWS,  STOCHPR, YPS,  TOTCLD,
     *    FIREMS,   YTFLX,    YQFLX,    YUFLX,   YVFLX,
     *    SIG(NLEV),SIGKMH(NLEV),
     *    OMEGA,    SEQSST,
     *    DYNAMICE, DYNAMHEM,
     *    QFACTOR,  QNORWEG,    QICES,    QICEN,   QLIMIT,
     *    CO2PPV*1.E6, HOTLIMIT,
     *    NSREST,   NSREST_IBIS,VERSRESIN,FSAVE,   FNEXF,   ISTEPY,
     *    NRSTRT,   NSTOP,
     *    MBDATE,   MCDATE,     MCSEC,    JYEAR,   JMNTH,   JDAY,
     *    RUNTITLE, SHOWMAPS,
     *    BASEDATE, ENDDATE,    NSAVE,    NHISI,   NDHIS 
#ifdef isotrac
     *   ,yq_o(1,1,nlev), yq_h(1,1,nlev), frainc_o, frainc_h,
     *    fsnowc_o, fsnowc_h, frains_o, frains_h, fsnows_o, fsnows_h,
     *    yqflx_o, yqflx_h
#endif
     *   )
C
CL                 CONVERT SENSIB, VAP AND WIND STRESS FLUXES FROM LSX'S
CL                 DOWNWARD FLUXES TO UPWARD
C
      DO 320 JJ=1,NOREC
        DO 322 JI=1,NLON
          YTFLX(JI,JJ) = -YTFLX(JI,JJ)
          YQFLX(JI,JJ) = -YQFLX(JI,JJ)
          YUFLX(JI,JJ) = -YUFLX(JI,JJ)
          YVFLX(JI,JJ) = -YVFLX(JI,JJ)
#ifdef isotrac
          yqflx_o(ji,jj) = -yqflx_o(ji,jj)
          yqflx_h(ji,jj) = -yqflx_h(ji,jj)
#endif
  322   CONTINUE
  320 CONTINUE

#ifdef isotracreg
      do it=1,nrtr
        do jj=1,norec
          do ji=1,nlon
            rtrfv(ji,jj,it) = yqflx(ji,jj)*
     *                 cvmgt( cvmgt(1.,0.,nint(sourc_i(ji,jj)).eq.it),
     *                        rtrq(ji,jj,nlev,it)/yq(ji,jj,nlev),
     *                        yqflx(ji,jj).ge.0. )
          enddo
        enddo
      enddo
#endif
 
C
CL                 OVERALL LOOP OVER LATITUDE FOR CONVECTION AND
CL                 VERTICAL DIFFUSION. (BETTER TO COMBINE INTO ONE LOOP
CL                 FOR PARALLELIZATION).
C
CL                 INITIALIZE HISTORY PBL HEIGHTS
C
C     For daily min,max pbl heights:
C     IF (NSTEP.EQ.NRSTRT .OR. MOD(NSTEP,NINT(86400./DTIME)).EQ.1) THEN
C       CALL RESETR (HEIMIN, NLON*NOREC,  1.E20)
C       CALL RESETR (HEIMAX, NLON*NOREC, -1.E20)
C     ENDIF
C     For 03:00, 15:00 pbl heights:
      IF (NSTEP.EQ.NRSTRT) THEN
        CALL ZERO (HEIMIN, NLON*NOREC)
        CALL ZERO (HEIMAX, NLON*NOREC)
      ENDIF
C
CL              ZERO SURFACE PRECIP RATES, AND PRECIP FALL RATE 
CL              FOR COSMOGENIC ISOTOPES (ACCUMULATED BY CONVAD-REEVAP)
C
      CALL ZERO (FRAINS, NLON*NOREC)
      CALL ZERO (FSNOWS, NLON*NOREC)
      CALL ZERO (FRAINC, NLON*NOREC)
      CALL ZERO (FSNOWC, NLON*NOREC)
#ifdef isotrac
      call zero (frains_o, nlon*norec)
      call zero (frains_h, nlon*norec)
      call zero (fsnows_o, nlon*norec)
      call zero (fsnows_h, nlon*norec)
      call zero (frainc_o, nlon*norec)
      call zero (frainc_h, nlon*norec)
      call zero (fsnowc_o, nlon*norec)
      call zero (fsnowc_h, nlon*norec)
#endif
#ifdef isotracreg
      call zero (rtrprs,   nlon*norec*nrtr)
      call zero (rtrpss,   nlon*norec*nrtr)
      call zero (rtrprc,   nlon*norec*nrtr)
      call zero (rtrpsc,   nlon*norec*nrtr)
#endif

#ifdef cosmog
      CALL ZERO (PRECFLX, NLON*NOREC*NLEVP)
#endif
C
c     Horizontally diffuse heat reservoir fplumt slightly, in an
c     attempt to reduce intense stationary blobs of precip in ITCZ.
c     Second param is diffusion coeff (m s-2; 1.e4 ok, 1.e6 blows)
c     call difgrid (fplumt, 1.e4, dtime)
C
#if defined ( omp )
#ifndef isotrac
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared (yu_eta, yv_eta, yt_eta, yq,
!$OMP*         ylwc_c, ylwc_a, ylwc_s, ycld,
!$OMP*         fplumt, fplumq, fplumu, fplumv,
!$OMP*         tracer, traflux, ypm,  ypthic, ypmcap,
!$OMP*         plumeprec, heimin, heimax,
!$OMP*         yum1, yu, yvm1,   yv,     ytm1,   yt,  yq_sig,
!$OMP*         yqm1,
!$OMP*         yps, yfu, yfv, yt2)
!$OMP* private (jj)
!$OMP* schedule(dynamic,2)
#else
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP* shared (yu_eta, yv_eta, yt_eta, yq,
!$OMP*         ylwc_c, ylwc_a, ylwc_s, ycld,
!$OMP*         fplumt, fplumq, fplumu, fplumv,
!$OMP*         tracer, traflux, ypm,  ypthic, ypmcap,
!$OMP*         plumeprec, heimin, heimax,
!$OMP*         yum1, yu, yvm1,   yv,     ytm1,   yt,  yq_sig,
!$OMP*         yqm1,
!$OMP*         yps, yfu, yfv, yt2,
!$OMP*         yq_o, yq_h, ylwc_c_o, ylwc_c_h, ylwc_a_o,
!$OMP*         ylwc_a_h, ylwc_s_o, ylwc_s_h, fplumq_o, fplumq_h,
!$OMP*         plumepr_o, plumepr_h,
!$OMP*         rtrq, rtrfq, plumepr_t,
!$OMP*         yqm1_o, yqm1_h,
!$OMP*         rtrqm1)
!$OMP* private (jj)
!$OMP* schedule(dynamic,2)
#endif
#endif
C
C======================
      DO 400 JJ=1,NOREC
C======================
C
CL                 LOCAL SATURATION ADJUSTMENT, VERTICAL CONVECTION
CL                 (PBL AND FREE), AND REEVAPORATION OF PRECIP. ADDS
CL                 FLUX CONVERGENCES TO RESERVOIRS FPLUM[U,V,T,Q,TRA].
CL                 STEPS CLOUD LIQUID WATER CONTENTS THRU THIS TIMESTEP,
CL                 AND SAVES PRECIP FOR NEXT TIMESTEP'S CALL TO SURFCTL.
C
        CALL CONVAD (YU_ETA, YV_ETA, YT_ETA, YQ,
     *               YLWC_C, YLWC_A, YLWC_S, YCLD,
     *               FPLUMT, FPLUMQ, FPLUMU, FPLUMV,
     *               TRACER, TRAFLUX, YPM,  YPTHIC, YPMCAP,
     *               PLUMEPREC, HEIMIN, HEIMAX, JJ 
#ifdef isotrac
     *              ,yq_o, yq_h, ylwc_c_o, ylwc_c_h, ylwc_a_o,
     *               ylwc_a_h, ylwc_s_o, ylwc_s_h, fplumq_o, fplumq_h,
     *               plumepr_o, plumepr_h,
     *               rtrq, rtrfq, plumepr_t
#endif
     *              )
C
CL                 "BACKGROUND" VERTICAL DIFFUSION. USES UP SOME OF
CL                 THE REVERVOIRS FPLUM[U,V,T,Q] (CONVERTING THE
CL                 USAGE FOR U,V,T FROM ETA TO SIGMA COORDS). STORES
CL                 U,V,T TENDENCIES IN YT2,YFU,YFV, AND CHANGES YQM1
CL                 AND TRACER FIELDS. VDIF_SIG DOES U,V,T ON SIGMA
CL                 COORDS AND VDIF_ETA DOES Q,TRACERS ON ETA COORDS.
C
        CALL VDIF_SIG (YUM1, YU, YVM1,   YV,     YTM1,   YT,  YQ_SIG,
     *                 FPLUMU, FPLUMV, FPLUMT, JJ)
C
        CALL VDIF_ETA (YQM1, YQ, YU_ETA, YV_ETA, YT_ETA, FPLUMQ, TRACER,
     *                 YPM, YPTHIC, YPMCAP, JJ
#ifdef isotrac
     *                ,yqm1_o, yqm1_h, yq_o, yq_h, fplumq_o, fplumq_h,
     *                 ylwc_c_o, ylwc_c_h, ylwc_a_o, ylwc_a_h,
     *                 ylwc_s_o, ylwc_s_h,
     *                 rtrqm1, rtrq, rtrfq
#endif
     *                )
C
CL                 GRAVITY WAVE DRAG (FOR U AND V)
C
        CALL GWD (YU, YV, YT, YPS, YFU, YFV, YT2, JJ)
C
C=============
  400 CONTINUE
C=============
C
CL              STEP TRACER SOURCES AND SINKS FOR #NTRACA+1 TO #NTRACE,
CL              AND PRINT TRACER DIAGNOSTICS AT TIME N.
CL              (NOT TRACERS #1-3, WHICH ARE LIQUID WATER CONTENTS,
CL              OR TRACERS #4-7, WHICH ARE FLUX RESERVOIRS;
CL              NTRACA=7 SET IN COMMONS)
C
#ifdef cosmog
      CALL COSMOG (TRACER, YT_ETA, YQ, YPM, YPTHIC, PRECFLX)
#else
      CALL TRABUDG (TRACER, YPTHIC)
#endif
C
C                  SAVE 03:00 AND 15:00 (WAS DAILY MIN,MAX) PBL HEIGHTS
C                  ONCE EACH DAY (SET BY CONVAD-PLUME)
C
C     IF (MOD(NSTEP,NINT(86400./DTIME)).EQ.0) THEN
        CALL PUTHIS (HEIMIN, 63)
        CALL PUTHIS (HEIMAX, 64)
C     ENDIF
C
C                  WRITE GRAVITY WAVE DRAG DIAGNOSTICS
C
      CALL GWDO (NOUT)
C
C                  SAVE LAND-OCEAN ORO-LIKE MASK AVERAGED TO AGCM GRID
C                  (0=OCEAN/SEAICE, 1=LAND, 2=ICESHEET),ROUNDED TO INTS
C
      DO 8000 JJ=1,NOREC
        DO 8002 JI=1,NLON
          IF (YOROOCN(JI,JJ).GE.0.5) THEN
            HISTMP(JI,JJ,1) = 0.
          ELSE 
            IF (YOROLAND(JI,JJ).GE.YOROICE(JI,JJ)) THEN
              HISTMP(JI,JJ,1) = 1.
            ELSE 
              HISTMP(JI,JJ,1) = 2.
            ENDIF
          ENDIF
 8002   CONTINUE
 8000 CONTINUE
      CALL PUTHIS (HISTMP, 11)
C
CL                 SAVE NET SW+LW DOWNWARD FLUX AT SURFACE (W M-2),
CL                 3-D SW,LW,DOWN-UP FLUXES AT INTERFACES (W M-2),
CL                 3-D CLOUD LIQUID WATER CONTENTS
C
      CALL PUTHIS (SLDFLX, 39)
      CALL PUTHIS (YQRS, 40)
      CALL PUTHIS (YQRL, 41)
C
CL                 INCREMENT BUDGET NET RADIATIVE HEATING (W M-2)
C
      DO 8020 JK=1,NLEV
        DO 8022 JJ=1,NOREC
          ZWEI = DTIME/NLON
          DO 8024 JI=1,NLON
            DHRADS(JJ) = DHRADS(JJ) 
     *                 + (YQRS(JI,JJ,JK)-YQRS(JI,JJ,JK+1))*ZWEI
            DHRADL(JJ) = DHRADL(JJ)
     *                 + (YQRL(JI,JJ,JK)-YQRL(JI,JJ,JK+1))*ZWEI
 8024     CONTINUE
 8022   CONTINUE
 8020 CONTINUE
C
CL                 SAVE 3-D CLOUD LIQUID WATER CONTENTS
C
      IF (TESTHIST(59,1)) THEN
        DO 8030 JK=1,NLEV
          DO 8032 JJ=1,NOREC
            DO 8034 JI=1,NLON
              HISTMP(JI,JJ,JK) =
     *          YLWC_C(JI,JJ,JK) + YLWC_A(JI,JJ,JK) + YLWC_S(JI,JJ,JK)
 8034       CONTINUE
 8032     CONTINUE
 8030   CONTINUE
        CALL PUTHIS (HISTMP, 59)
      ENDIF
C
      CALL PUTHIS (YLWC_C, 54)
      CALL PUTHIS (YLWC_A, 55)
      CALL PUTHIS (YLWC_S, 57)
#ifdef isotrac
      call puthis (ylwc_c_o, 71)
      call puthis (ylwc_c_h, 72)
      call puthis (ylwc_a_o, 73)
      call puthis (ylwc_a_h, 74)
      call puthis (ylwc_s_o, 75)
      call puthis (ylwc_s_h, 76)
#endif
C
CL                 SAVE VISIBLE AND NEAR-IR SOLAR FLUXES AT SURFACE.
C
      IF (TESTHIST(50,1) .OR. TESTHIST(51,1) .OR. TESTHIST(52,1) .OR.
     *    TESTHIST(56,1)) THEN
        DO 810 JJ=1,NOREC
           DO 812 JI=1,NLON
C             Visible absorbed (direct and diffuse combined)
              ZAV = (1.-ALBSUD(JI,JJ,1))*SOLSUD(JI,JJ,1)
     *            + (1.-ALBSUI(JI,JJ,1))*SOLSUI(JI,JJ,1)
C             Visible incident (direct and diffuse combined)
              ZIV = SOLSUD(JI,JJ,1) + SOLSUI(JI,JJ,1)
C             Near-ir absorbed (direct and diffuse combined)
              ZAN = (1.-ALBSUD(JI,JJ,2))*SOLSUD(JI,JJ,2)
     *            + (1.-ALBSUI(JI,JJ,2))*SOLSUI(JI,JJ,2)
C             Near-ir incident (direct and diffuse combined)
              ZIN = SOLSUD(JI,JJ,2) + SOLSUI(JI,JJ,2)
              HISTMP(JI,JJ,1) = ZAV
              HISTMP(JI,JJ,2) = ZIV
              HISTMP(JI,JJ,3) = ZAN
              HISTMP(JI,JJ,4) = ZIN
  812      CONTINUE
  810   CONTINUE
        CALL PUTHIS (HISTMP(1,1,1), 50)
        CALL PUTHIS (HISTMP(1,1,2), 51)
        CALL PUTHIS (HISTMP(1,1,3), 52)
        CALL PUTHIS (HISTMP(1,1,4), 56)
      ENDIF

#if defined (LKAEROSOL)
      call puthis (sfcirrad, 67)
#endif
C
CL                 SAVE OTHER RADIATION QUANTITIES, STORED IN H* ARRAYS
C
      CALL PUTHIS (HFRSA, 24)
      CALL PUTHIS (HSABTP, 27)
      CALL PUTHIS (HSABAT, 28)
      CALL PUTHIS (HCLRSS, 34)
      CALL PUTHIS (HCLRST, 32)
      CALL PUTHIS (HSOLIN, 42)
C
      CALL PUTHIS (FIRABS, 23)
      CALL PUTHIS (HFRLA,  25)
      CALL PUTHIS (HFIRTP, 26)
      CALL PUTHIS (HCLRLS, 33)
      CALL PUTHIS (HCLRLT, 31)
C
      CALL PUTHIS (HFRSAC, 65)
 
c----
c     Dump long-wave heating rates, toa and sfc fluxes (used to
c     diagnose new ir code with other trace gases)
c     write(6,8880) 'yqrl',
c    *((jj,jk,(yqrl(ji,jj,jk)*86400,ji=1,nlon,4),jj=1,norec),jk=1,nlevp)
c     write(6,8880) 'hfirtp',
c    * (jj,0,(hfirtp(ji,jj),ji=1,nlon,4),jj=1,norec)
c     write(6,8880) 'hfrla',
c    * (jj,0,(hfrla(ji,jj),ji=1,nlon,4),jj=1,norec)
c8880 format(//a,':'/(2i4,19f10.6))
c
c     write(6,8881)
c8881 format(//'Zonal: yqrl:')
c     do 8882 jk=1,nlevp
c     do 8882 jj=1,norec
c       zzyqrl = 0.
c       do 8883 ji=1,nlon
c         zzyqrl = zzyqrl + yqrl(ji,jj,jk)*86400.
c8883   continue
c       zzyqrl = zzyqrl/nlon
c       write(6,8884) jj,jk,zzyqrl
c8884   format(2i4,f20.10)
c8882 continue
c
c     write(6,8885)
c8885 format(//'Zonal: hfirtp frla:')
c     do 8886 jj=1,norec
c       zzfirtp = 0.
c       zzfrla = 0.
c       do 8887 ji=1,nlon
c         zzfirtp = zzfirtp + hfirtp(ji,jj)
c         zzfrla  = zzfrla  + hfrla(ji,jj)
c8887   continue
c       zzfirtp = zzfirtp/nlon
c       zzfrla = zzfrla/nlon
c       write(6,8888) jj,zzfirtp, zzfrla
c8888   format(i4,2f20.10)
c8886 continue
c----
C
CL                 ACCUMULATE RADIATION AT TOP OF MODEL FOR BUDGETS.
CL                 ACCOUNT FOR "ABOVE-MODEL" ABSORBTION HSABAT IN
CL                 DSOLRF (FOR PLANETARY ALBEDO), BUT NOT IN DRADTP
CL                 (FOR MODEL ENERGY BUDGET). ALSO ACCUMULATE BUDGET
CL                 LOWEST-LEVEL TEMPS AND TOTAL CLOUDINESS.
C
      DO 830 JJ=1,NOREC
        ZZ = DTIME/NLON
        DO 832 JI=1,NLON
          DSOLIN(JJ) = DSOLIN(JJ) + HSOLIN(JI,JJ)*ZZ
          DSOLRF(JJ) = DSOLRF(JJ)
     *                 +(HSOLIN(JI,JJ)-HSABAT(JI,JJ)-HSABTP(JI,JJ))*ZZ
          DFIRTP(JJ) = DFIRTP(JJ) + HFIRTP(JI,JJ)*ZZ
          DRADTP(JJ) = DRADTP(JJ) + (HSABTP(JI,JJ)-HFIRTP(JI,JJ))*ZZ
          DTAIR(JJ)  = DTAIR(JJ)   + YT(JI,JJ,NLEV)*ZZ
          DTOTCLD(JJ)= DTOTCLD(JJ) + TOTCLD(JI,JJ)*ZZ
  832   CONTINUE
  830 CONTINUE
C
CL                 SAVE SURFACE-FLUX HISTORY VARIABLES (STORE AS:
CL                 SENSIBLE HEAT UPWARD, EVAPORATION UPWARD AS ~W/M2,
CL                 WIND STRESSES DOWNWARD (SAME SIGN AS WIND SPEED))
C
      CALL PUTHIS (YTFLX, 16)
      IF (TESTHIST(17,1) .OR. TESTHIST(18,1) .OR. TESTHIST(19,1) .OR.
     *    TESTHIST(37,1) .OR. TESTHIST(38,1)) THEN
        DO 900 JJ=1,NOREC
          DO 902 JI=1,NLON
            HISTMP(JI,JJ,1) =  YQFLX(JI,JJ) * LATVAP
            HISTMP(JI,JJ,2) = -YU(JI,JJ,NLEV)*YUFLX(JI,JJ)
            HISTMP(JI,JJ,3) = -YV(JI,JJ,NLEV)*YVFLX(JI,JJ)
            HISTMP(JI,JJ,4) = -YUFLX(JI,JJ)
            HISTMP(JI,JJ,5) = -YVFLX(JI,JJ)
  902     CONTINUE
  900   CONTINUE
        CALL PUTHIS (HISTMP(1,1,1), 17)
        CALL PUTHIS (HISTMP(1,1,2), 18)
        CALL PUTHIS (HISTMP(1,1,3), 19)
        CALL PUTHIS (HISTMP(1,1,4), 37)
        CALL PUTHIS (HISTMP(1,1,5), 38)
      ENDIF
#ifdef isotrac
      if (testhist(77,1) .or. testhist(78,1)) then
        do jj=1,norec
          do ji=1,nlon
            histmp(ji,jj,1) = yqflx_o(ji,jj) * latvap
            histmp(ji,jj,2) = yqflx_h(ji,jj) * latvap
          enddo
        enddo
        call puthis (histmp(1,1,1), 77)
        call puthis (histmp(1,1,2), 78)
      endif
#endif
C
CL                 SAVE NEW (TIME N) STRATIFORM AND CONVECTIVE PRECIP,
CL                 RAINFALL, SNOWFALL, PRECIP. CONVERT TO M/S.
C
      IF (TESTHIST(12,1) .OR. TESTHIST(13,1) .OR. TESTHIST(14,1) .OR.
     *    TESTHIST(15,1) .OR. TESTHIST(30,1)) THEN
        DO 1000 JJ=1,NOREC
          DO 1002 JI=1,NLON
             HISTMP(JI,JJ,1) = (FRAINS(JI,JJ)+FSNOWS(JI,JJ))/RHOH2O
             HISTMP(JI,JJ,2) = (FRAINC(JI,JJ)+FSNOWC(JI,JJ))/RHOH2O
             HISTMP(JI,JJ,3) = (FRAINC(JI,JJ)+FRAINS(JI,JJ))/RHOH2O
             HISTMP(JI,JJ,4) = (FSNOWC(JI,JJ)+FSNOWS(JI,JJ))/RHOH2O
             HISTMP(JI,JJ,5) = HISTMP(JI,JJ,1) + HISTMP(JI,JJ,2)
 1002     CONTINUE
 1000   CONTINUE
        CALL PUTHIS (HISTMP(1,1,1), 12)
        CALL PUTHIS (HISTMP(1,1,2), 13)
        CALL PUTHIS (HISTMP(1,1,3), 14)
        CALL PUTHIS (HISTMP(1,1,4), 15)
        CALL PUTHIS (HISTMP(1,1,5), 30)
      ENDIF
#ifdef isotrac
      if (testhist(79,1) .or. testhist(80,1) .or. testhist(81,1) .or.
     *    testhist(82,1) .or. testhist(83,1) .or. testhist(84,1) .or.
     *    testhist(85,1) .or. testhist(86,1) 
     *   ) then
        do jj=1,norec
          do ji=1,nlon
             histmp(ji,jj,1) = (frains_o(ji,jj)+fsnows_o(ji,jj))/rhoh2o
             histmp(ji,jj,2) = (frains_h(ji,jj)+fsnows_h(ji,jj))/rhoh2o
             histmp(ji,jj,3) = (frainc_o(ji,jj)+fsnowc_o(ji,jj))/rhoh2o
             histmp(ji,jj,4) = (frainc_h(ji,jj)+fsnowc_h(ji,jj))/rhoh2o
             histmp(ji,jj,5) = (frainc_o(ji,jj)+frains_o(ji,jj))/rhoh2o
             histmp(ji,jj,6) = (frainc_h(ji,jj)+frains_h(ji,jj))/rhoh2o
             histmp(ji,jj,7) = (fsnowc_o(ji,jj)+fsnows_o(ji,jj))/rhoh2o
             histmp(ji,jj,8) = (fsnowc_h(ji,jj)+fsnows_h(ji,jj))/rhoh2o
           enddo
         enddo
        call puthis (histmp(1,1,1), 79)
        call puthis (histmp(1,1,2), 80)
        call puthis (histmp(1,1,3), 81)
        call puthis (histmp(1,1,4), 82)
        call puthis (histmp(1,1,5), 83)
        call puthis (histmp(1,1,6), 84)
        call puthis (histmp(1,1,7), 85)
        call puthis (histmp(1,1,8), 86)
      endif
#endif

#ifdef isotracreg
      do it=1,nrtr
        if (testhist(100+nrtr+it,1)) then
          do jj=1,norec
            do ji=1,nlon
               histmp(ji,jj,1) = (rtrprc(ji,jj,it)+rtrprs(ji,jj,it)+
     *                       rtrpsc(ji,jj,it)+rtrpss(ji,jj,it))/rhoh2o
            enddo
          enddo
          call puthis (histmp(1,1,1), 100 + nrtr + it)
        endif
      enddo
#endif

CL                 COMPLETE BUDGET CALCULATIONS

      DO 1200 JJ=1,NOREC
        DO 1202 JI=1,NLON
          ZR(JI) = - (1.-ALBSUD(JI,JJ,1))*SOLSUD(JI,JJ,1)
     *             - (1.-ALBSUI(JI,JJ,1))*SOLSUI(JI,JJ,1)
     *             - (1.-ALBSUD(JI,JJ,2))*SOLSUD(JI,JJ,2)
     *             - (1.-ALBSUI(JI,JJ,2))*SOLSUI(JI,JJ,2)
     *             +  HFRLA(JI,JJ)
C         ZR(JI) = - HFRSA(JI,JJ)
C    *             - FIRABS(JI,JJ) + FIREMS(JI,JJ)
          ZS(JI) = YTFLX(JI,JJ)
          ZL(JI) = YQFLX(JI,JJ)
     *             *((YT_ETA(JI,JJ,NLEV)-TMELT)*CPWV + LATSUB)
     *           - (FRAINC(JI,JJ)+FRAINS(JI,JJ))
     *             *((MAX(YT_ETA(JI,JJ,NLEV),TMELT)-TMELT)*CH2O+LATICE)
     *           - (FSNOWC(JI,JJ)+FSNOWS(JI,JJ))
     *             *((MIN(YT_ETA(JI,JJ,NLEV),TMELT)-TMELT)*CICE)
          ZE(JI) = YQFLX(JI,JJ)
          ZP(JI) = FRAINC(JI,JJ) + FSNOWC(JI,JJ)
     *           + FRAINS(JI,JJ) + FSNOWS(JI,JJ)
 1202   CONTINUE
        DRFLUX(JJ) = DRFLUX(JJ) + SSUM(NLON,ZR,1)*DTIME/NLON
        DSFLUX(JJ) = DSFLUX(JJ) + SSUM(NLON,ZS,1)*DTIME/NLON
        DLFLUX(JJ) = DLFLUX(JJ) + SSUM(NLON,ZL,1)*DTIME/NLON
        DEFLUX(JJ) = DEFLUX(JJ) + SSUM(NLON,ZE,1)*DTIME/NLON
        DPFLUX(JJ) = DPFLUX(JJ) + SSUM(NLON,ZP,1)*DTIME/NLON
 1200 CONTINUE
C
      RETURN
      END
C
C ***************************** TIMESM ********************************
C

      subroutine timesm (arrm1, arr,
#ifdef isotrac
     *                   arrm1_o, arr_o, arrm1_h, arr_h, arrm1_t, arr_t,
#endif
     *                   iflag)

c----------------------------------------------------------------------
c
cl           If iflag eq 1, stores first part of time filter
cl           (1-2*eps)*arr (time n) + eps*arrm1 (time n-1) in
cl           storage, i.e., array filtarr 
cl
cl           If iflag ne 1, retrieves first part, adds eps*arr(time n+1)
cl           and stores in arrm1.
cl
cl           (arrm1 corresponds to yum1 in coml10, arr to yu in coml30.)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
#include <COMSDT.h>
c----------------------------------------------------------------------
#include <COMTIM.h>
c----------------------------------------------------------------------
      parameter (ntimfilt = (6*nlev + 1)*nlon*norec)
      dimension arrm1(ntimfilt), arr(ntimfilt)
#ifdef isotrac
      parameter (n_i = nlev*nlon*norec)
      parameter (n_t = nlev*nlon*norec*nrtr)
      dimension arrm1_o(n_i), arr_o(n_i), arrm1_h(n_i), arr_h(n_i),
     *          arrm1_t(n_t), arr_t(n_t)
#endif
 
c     filtarr is storage array.
c     Put in a common, not local saved array, to avoid
c     enlarging the size of the executable file.
      common /cfiltarc/ filtarr(ntimfilt)
#ifdef isotrac
      common /cfilt_i/  filt_o(n_i), filt_h(n_i)
      common /cfilt_t/  filt_t(n_t)
#endif

c----------------------------------------------------------------------
 
c-------------------------
      if (iflag.eq.1) then
c-------------------------
 
c        First part of time filter: put in storage
 
         do 100 i=1,ntimfilt
           filtarr(i) = (1.-2.*eps)*arr(i) + eps*arrm1(i)
  100    continue
#ifdef isotrac
         do 110 i=1,n_i
           filt_o(i)  = (1.-2.*eps)*arr_o(i)+eps*arrm1_o(i)
           filt_h(i)  = (1.-2.*eps)*arr_h(i)+eps*arrm1_h(i)
  110    continue
         do i=1,n_t
           filt_t(i)  = (1.-2.*eps)*arr_t(i)+eps*arrm1_t(i)
         enddo
#endif
 
c---------
      else
c---------
 
c        Second part of time filter:get from storage, put into coml10
 
         do 200 i=1,ntimfilt
           arrm1(i) = filtarr(i) + eps*arr(i)
  200    continue
#ifdef isotrac
         do 210 i=1,n_i
           arrm1_o(i)=filt_o(i)  + eps*arr_o(i)
           arrm1_h(i)=filt_h(i)  + eps*arr_h(i)
  210    continue
         do i=1,n_t
           arrm1_t(i) = filt_t(i)  + eps*arr_t(i)
         enddo
#endif
c----------
      endif
c----------
 
      return
      end
C
C ***************************** TRABUDG *******************************
C
      SUBROUTINE TRABUDG (TRACER, YPTHIC)
C----------------------------------------------------------------------
C
CL            (i)  STEP TRACER FIELDS DUE TO SOURCES AND SINKS
CL            (ii) CALCULATE AND PRINT TRACER DIAGNOSTICS
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      dimension
     *   tracer(nlon,norec,nlev,ntrace),
     *   ypthic(nlon,norec,nlev)
C----------------------------------------------------------------------
      dimension
     *   gz(nlev,ntrace), gtrop(ntrace), gstrat(ntrace), gtot(ntrace),
     *   zmasz(nlev),
     *   gmax(nlev,ntrace)
      logical firstwrite, printracer
      save firstwrite, printracer
      data firstwrite, printracer /.true., .true./
 
c     sigpaus = tropopause for purposes of tropos/stratos diagnostics
      dimension sigpaus(norec)
c     sigsink = sigma level below which tracers are zeroed
      save sigsink
c     data sigsink /.400/
      data sigsink /99.999/
C----------------------------------------------------------------------
 
c        Only print tracer diagnostics if # true tracers (ntracb) >= 2
c        (Even if want no true tracers, ntracb must = 1 for code) 

      if (ntracb.lt.2) return
 
c        Zero all tracers below sigma = sigsink ("tropospheric sink")
 
      do 50 n=ntraca+1,ntrace
        do 52 jk=1,nlev
          if(sig(jk).ge.sigsink) call zero (tracer(1,1,jk,n),nlon*norec)
   52   continue
   50 continue
 
c---------------------------
c        Diagnostics section:
c---------------------------
 
      if (.not. printracer) return
 
c        Only compute and write diagnostics at day boundaries
 
      if ( mod(nstep,nint(1.*86400./dtime)).ne.0 .and.
     *     nstep.ne.nrstrt ) return
 
c        Calculate atmos-mass amounts
 
      call zero (zmasz, nlev)
      zmastot  = 0.
      zmastrop = 0.
      zmastrat = 0.
 
      do 150 jk=1,nlev
        do 152 jj=1,norec
          sigpaus(jj) = .100*cs(jj) + .280*(1.-cs(jj))
          zwei = cosbud(jj)/nlon
          do 154 ji=1,nlon
            zdm = (ypthic(ji,jj,jk)/gravit) * zwei
            zmasz(jk) = zmasz(jk) + zdm
            zmastot   = zmastot   + zdm
            if (sig(jk).le.sigpaus(jj)) then
              zmastrat = zmastrat + zdm
            else
              zmastrop = zmastrop + zdm
            endif
  154     continue
  152   continue
  150 continue
 
c        Calculate tracer amounts
 
      do 200 n=ntraca+1,ntrace
 
        call zero (gz(1,n), nlev)
        gtot(n)   = 0.
        gtrop(n)  = 0.
        gstrat(n) = 0.
 
        do 210 jk=1,nlev
          do 212 jj=1,norec
            zwei = cosbud(jj)/nlon
            do 214 ji=1,nlon
              zdtrac = tracer(ji,jj,jk,n)*(ypthic(ji,jj,jk)/gravit)*zwei
              gz(jk,n) = gz(jk,n) + zdtrac
              gtot(n)  = gtot(n)  + zdtrac
              if (sig(jk).ge.0.200) then
                gtrop(n)  = gtrop(n)  + zdtrac
              else
                gstrat(n) = gstrat(n) + zdtrac
              endif
  214       continue
  212     continue
  210   continue
 
        do 220 jk=1,nlev
          gz(jk,n) = gz(jk,n) / zmasz(jk)
  220   continue
        gtot(n)   = gtot(n)   / zmastot
        gtrop(n)  = gtrop(n)  / zmastrop
        gstrat(n) = gstrat(n) / zmastrat
 
  200 continue
 
c        Write diagnostics (up to 3 tracers only) to files
c        fort.71,72,73, and max values to files fort.74,75,76.

      k1 = 9
      k2 = 1
      k3 = -1
 
      do 700 nt = 1, min (3, ntrace-ntraca)
        n = ntraca + nt
 
c          Layer-mean concentrations:
 
        iu  = 70 + nt
        if (firstwrite)
     *     write(iu,702) sigsink, sigpaus(1), sigpaus(norec/2),
     *                   nt, (sig(jk),jk=k1,k2,k3)
  702   format('tracer',5x,'sigsink=',f6.3,5x,'sigpaus(po/eq)=',2f7.3
     *         /i6
     *         /' calday',' runday',
     *          '      trop','     strat','     total',9f9.3/)
        write(iu,704) istepy*dtime/86400., nstep*dtime/86400.,
     *                gtrop(n),gstrat(n),gtot(n),(gz(jk,n),jk=k1,k2,k3)
  704   format(2f7.2, 3f10.7, 9f9.5)
        ier = flush(iu)
 
c          Layer-maximum concentrations:
 
        ium = 73 + nt
        do 712 jk=1,nlev
          gmax(jk,n) = -999.
          do 714 jj=1,norec
            do 716 ji=1,nlon
              gmax(jk,n) = max (gmax(jk,n), tracer(ji,jj,jk,n))
  716       continue
  714     continue
  712   continue
        if (firstwrite)
     *    write(ium,718) sigsink, sigpaus(1), sigpaus(norec/2),
     *                   nt,(sig(jk),jk=k1,k2,k3)
  718   format('tracer (max)',5x,'sigsink=',f6.3,5x,'sigpaus(po/eq)=',
     *         2f7.3
     *        /i6
     *        /' calday',' runday','      tottra','      totatm',
     *         9f9.3/)
        write(ium,719) istepy*dtime/86400., nstep*dtime/86400.,
     *                 gtot(n)*zmastot, zmastot,(gmax(jk,n),jk=k1,k2,k3)
  719   format(2f7.2, f12.6, f12.5, 9f9.5)
        ier = flush(ium)
 
  700 continue
      firstwrite = .false.
 
      return
      end
c
c **********************************************************************
c
#ifdef cosmog
      subroutine cosmog (tracer, yt, yq, ypm, ypthic, precflx)
c----------------------------------------------------------------------
c
c     If first timestep:
c       Reads 10-Be production rates (starbe) vs. latitude (alatbe)
c       and vertical column (azbe, g/cm2), for 6 different combinations
c       of geomagnetic field and solar force field, from file allbe.dat.
c       Vertical coordinate is mass of atmos column from TOA.
c       Units of starbe on this file are "stars per g.sec".
c       The 6 combinations, in order on the file, are:
c       (0.5,300), (0.5,600), (1,300), (1,600), (5,300), (5,600),
c       where the 1st value is geomagnetic field strength relative to
c       today, and the 2nd value is the solar force field parameter
c       (solar min/max).
c
c     Steps 10-Be tracer amounts (# 10-Be atoms/kg of air, ntraca+1 to
c     ntraca+6) through one timestep. Does same for 7-Be tracers
c     (ntraca+7 to ntraca+12), which have 2* yield (atoms per star)
c     and radioactive half life of 53.0 days.
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMBUD.h>
c----------------------------------------------------------------------
#include <COMCON.h>
c----------------------------------------------------------------------
#include <COMHIS.h>
c----------------------------------------------------------------------
#include <COMHIX.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
#include <COMMAP.h>
c----------------------------------------------------------------------
#include <COMTIM.h>
c----------------------------------------------------------------------
 
      dimension
     *   tracer(nlon,norec,nlev,ntrace),
     *   yt(nlon,norec,nlev),  yq(nlon,norec,nlev),
     *   ypm(nlon,norec,nlev), ypthic(nlon,norec,nlev),
     *   precflx(nlon,norec,nlevp)
 
c     For input from file allbe.dat (nb: ncasebe must = ntracb/2):
      parameter (nlatbe=11, nlevbe=16, ncasebe=6)
      dimension
     *  starbe(nlatbe,nlevbe,ncasebe), alatbe(nlatbe), azbe(nlevbe)
      character*80 chead
 
c     Production rates on agcm latitude grid and fine vertical grid
      parameter (nfine=1200, zfine=1.2e5)
      dimension jbe(norec), wjbe(norec), kbe(nfine), wkbe(nfine)
      save jbe, wjbe, kbe, wkbe, starbe
 
      dimension
     *  zrho(nlon), zdz(nlon), kfine(nlon)
 
      dimension
     *  gz(nlev,ntrace), gtrop(ntrace), gstrat(ntrace), gtot(ntrace),
     *  gflux(ntrace), zmasz(nlev)
      logical first, firstprint
      save first, firstprint, zareatot
      data first, firstprint /.true., .true./
 
c     sigpaus = tropopause for purposes of tropos/stratos diagnostics
      dimension sigpaus(norec)
 
      parameter (pi = 3.14159265358979)
 
c--------------------
      if (first) then
c--------------------
 
        if (ncasebe .ne. ntracb/2) then
          write (nout,100) ncasebe, ntracb
  100     format(/'*** Error (cosmog): ncasebe must equal 1/2 * ntracb'
     *           /'    ncasebe=',i4,'   ntracb=',i4)
          call endrun (-1)
        endif
 
c          Read starbe, alatbe and azbe from data file allbe.dat
 
        write (nout,200)
  200   format(/'Reading 10-Be production rate data file allbe.dat...')
        iu = 90
        open (iu, file='allbe.dat', status='old')
        do nc=1,ncasebe
          read (iu,205) chead, (alatbe(j),j=1,nlatbe)
  205     format(/a//9x,11(1x,f11.1))
          write (nout,210) chead(1:lenchr(chead))
  210     format('   ',a)
          do k=1,nlevbe
            read (iu,215) azbe(k), (starbe(j,k,nc),j=1,nlatbe)
  215       format(f7.2,2x,11(1x,e11.5))
          enddo
        enddo
        close (iu)
 
c          Calculate index and weights from agcm's alatg into alatbe
 
        do j=1,norec
          zlat = alatg(j)*180./pi
          zlat = max (alatbe(1)+.0001, min (alatbe(nlatbe)-.0001, zlat))
          do jj=2,nlatbe
            if (zlat.le.alatbe(jj)) then
              jbe(j)  = jj-1
              wjbe(j) = (alatbe(jj)-zlat)/(alatbe(jj)-alatbe(jj-1))
              goto 300
            endif
          enddo
  300     continue
        enddo
 
c          Calculate index and weights from nfine evenly spaced values
c          from 0 to zfine (Pa) into azbe (g/cm2). At each timestep,
c          index, from current levels to the nfine grid will
c          be set below (kfine), then kbe(kfine) will be used into azbe.
 
        do k=1,nfine
          zlev = ((k-.5)/nfine)*zfine * 0.1/gravit         ! Pa to g/cm2
          zlev = max (azbe(1)+.0001, min (azbe(nlevbe)-.0001, zlev))
          do kk=2,nlevbe
            if (zlev.le.azbe(kk)) then
              kbe(k)  = kk-1
              wkbe(k) = (azbe(kk)-zlev)/(azbe(kk)-azbe(kk-1))
              goto 400
            endif
          enddo
  400     continue
        enddo
 
        first = .false.
c----------
      endif
c----------
 
 
c        Overall loops over level and latitude
 
c==================
      do j=1,norec
        do k=1,nlev
c==================
          km1 = max (k-1,1)
 
          do i=1,nlon
            zrho(i)  = ypm(i,j,k) / (rair*(1.+zvir*yq(i,j,k))*yt(i,j,k))
            zdz(i)   = ypthic(i,j,k) / (zrho(i)*gravit)
            kfine(i) = nfine*(ypm(i,j,k)/zfine) + 1.0001
            kfine(i) = max (1, min (nfine, kfine(i)))
          enddo
 
c            Loop over 10-Be and 7-Be tracers
c            (combined into one loop since starbe's are same for each)
 
c-------------------------
          do nc=1,ntracb/2
c-------------------------
            n10 = ntraca + nc                                    ! 10-Be
            n7  = ntraca + ntracb/2 + nc                         !  7-Be
 
c              Step tracer amounts due to production (zprod) and wet
c              deposition (precflx). Change units of starbe from
c              "stars per g.sec" to "10-Be atoms/kg/s", where
c              1000 = g/kg and .0075 = # 10-Be atoms per star.
c              Multiply precflx terms by wet scavenging ratio = 100
c              = [10-Be concentration in precip / 10-Be conc. in air]
c              Do both for 10-Be, and for 7-Be with 2 * yield factor
c              and a half life of 53.0 days.
 
            zfac7 = 1. - exp(-log(2.)*dtime/(53.0*86400.))
 
            do i=1,nlon
              kf = kfine(i)
              zstar =
     *              wjbe(j) *    wkbe(kf) *starbe(jbe(j)  ,kbe(kf)  ,nc)
     *        + (1.-wjbe(j))*    wkbe(kf) *starbe(jbe(j)+1,kbe(kf)  ,nc)
     *        +     wjbe(j) *(1.-wkbe(kf))*starbe(jbe(j)  ,kbe(kf)+1,nc)
     *        + (1.-wjbe(j))*(1.-wkbe(kf))*starbe(jbe(j)+1,kbe(kf)+1,nc)
ctest         zstar = 0.
 
c10-Be:
              tracer(i,j,k,n10) = tracer(i,j,k,n10)
     *              + dtime * ( zstar * 1000. * .0075
     *                          + (   tracer(i,j,km1,n10)*precflx(i,j,k)
     *                              - tracer(i,j,k,n10)*precflx(i,j,k+1)
     *                            ) * 100. / (ypthic(i,j,k)/gravit)
     *                        )
c 7-Be:
              tracer(i,j,k,n7) = tracer(i,j,k,n7)
     *              + dtime * ( zstar * 1000. * .0075 * 2.
     *                          + (   tracer(i,j,km1,n7)*precflx(i,j,k)
     *                              - tracer(i,j,k, n7)*precflx(i,j,k+1)
     *                            ) * 100. / (ypthic(i,j,k)/gravit)
     *                        )
     *              - zfac7 * tracer(i,j,k,n7)
            enddo
 
c              Step tracer amounts due to dry deposition (lowest level
c              only), assuming dry deposition velocity = .002 m/s.
c              Also store net surface flux (# atoms/m2/s) for history.
 
            if (k.eq.nlev) then
              do i=1,nlon
                tracer(i,j,k,n10) = tracer(i,j,k,n10)
     *                          - dtime * tracer(i,j,k,n10) *.002/zdz(i)
                tracer(i,j,k,n7) = tracer(i,j,k,n7)
     *                          - dtime * tracer(i,j,k,n7)  *.002/zdz(i)
 
                histmp(i,j,nc) =
     *            tracer(i,j,k,n10) * ( precflx(i,j,k+1)*100.
     *                                  + zrho(i)*.002 )
                histmp(i,j,nc+ntracb/2)=
     *            tracer(i,j,k,n7)  * ( precflx(i,j,k+1)*100.
     *                                  + zrho(i)*.002 )
              enddo
            endif
 
c--------------
          enddo
c--------------
 
c============
        enddo
      enddo
c============
 
c        Save history fields
 
      do nc = 1,ntracb
        n = ntraca + nc
        call puthis (tracer(1,1,1,n), ihtrace+nc-1)
        call puthis (histmp(1,1,nc),  ihtrace+40+nc-1)
      enddo
 
c        Remainder is for writing diagnostics to scratch files.
c        Only compute and write diagnostics at day boundaries
 
      if ( mod(nstep,nint(1.*86400./dtime)).ne.0 .and.
     *     nstep.ne.nrstrt ) return
 
c        Calculate atmos-mass amounts
 
      call zero (zmasz, nlev)
      zmastot  = 0.
      zmastrop = 0.
      zmastrat = 0.
      zareatot = 0.
 
      do jj=1,norec
        sigpaus(jj) = .100*cs(jj) + .280*(1.-cs(jj))
        zwei = cosbud(jj)/nlon
        do jk=1,nlev
          do ji=1,nlon
            zdm = (ypthic(ji,jj,jk)/gravit) * zwei
            zmasz(jk) = zmasz(jk) + zdm
            zmastot   = zmastot   + zdm
            if (sig(jk).le.sigpaus(jj)) then
              zmastrat = zmastrat + zdm
            else
              zmastrop = zmastrop + zdm
            endif
          enddo
        enddo
        zareatot = zareatot + cosbud(jj)
      enddo
 
c        Calculate tracer amounts and surface tracer flux
 
      do nc = 1, ntracb
        n = ntraca + nc
 
        call zero (gz(1,n), nlev)
        gtot(n)   = 0.
        gtrop(n)  = 0.
        gstrat(n) = 0.
        gflux(n)  = 0.
 
        do jj=1,norec
          zwei = cosbud(jj)/nlon
          do jk=1,nlev
            do ji=1,nlon
              zdtrac = tracer(ji,jj,jk,n)*(ypthic(ji,jj,jk)/gravit)*zwei
              gz(jk,n) = gz(jk,n) + zdtrac
              gtot(n)  = gtot(n)  + zdtrac
              if (sig(jk).ge.0.200) then
                gtrop(n)  = gtrop(n)  + zdtrac
              else
                gstrat(n) = gstrat(n) + zdtrac
              endif
            enddo
          enddo
          do ji=1,nlon
            gflux(n) = gflux(n) + histmp(ji,jj,nc)*zwei
          enddo
        enddo
 
        do jk=1,nlev
          gz(jk,n) = gz(jk,n) / zmasz(jk)
        enddo
        gtot(n)   = gtot(n)   / zmastot
        gtrop(n)  = gtrop(n)  / zmastrop
        gstrat(n) = gstrat(n) / zmastrat
        gflux(n)  = gflux(n)  / zareatot
 
      enddo
 
c        Write diagnostics to files to files fort.[71-76] for 10-Be,
c        and to fort.[81-86} for 7-Be
c
 
      k1 = 9
      k2 = 1
      k3 = -1
 
c-------------------
      do iloop = 1,2                    ! loop #1 for 10-Be, #2 for 7-Be
c-------------------
      noffn  = (iloop-1)*(ntracb/2)
      noffiu = (iloop-1)*10
 
      do nh = 1, ntracb/2
        nc = noffn + nh
        n  = ntraca + noffn + nh
        iu = 70 + noffiu + nh
 
        if (firstprint) then
          ih = ihtrace + nc - 1
          write(iu,702)
     *      nc, cindex(1,ih), cindex(2,ih), sigpaus(1),sigpaus(norec/2),
     *      (sig(jk),jk=k1,k2,k3)
  702     format(
     *      ' tracer=',i3,'  (',a,2x,a,')   sigpaus(po/eq)=',2f7.3
     *     /' calday',' runday',
     *      '     trop','    strat','    total','   sfcflx', 9f9.3/)
        endif
 
        write(iu,704) istepy*dtime/86400., nstep*dtime/86400.,
     *                nint(gtrop(n)), nint(gstrat(n)), nint(gtot(n)),
     *                gflux(n), (nint(gz(jk,n)),jk=k1,k2,k3)
  704   format(2f7.2, 3i9, f9.2, 9i9)
 
        ier = flush(iu)
      enddo
 
c----------
      enddo
c----------
 
      firstprint = .false.
 
      return
      end
#endif
C
C **************************** PRINGLOB ********************************
C
      subroutine pringlob (ypthicm, ypthic, clabel)
 
c     Computes and prints out global mean atmospheric budget quantites,
c     similar to BUDCLC (but intended for diagnostics within timesteps).
 
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMDYR.h>
C----------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COML10.h>
C----------------------------------------------------------------------
#include <COML30.h>
C----------------------------------------------------------------------
#include <COML50.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      dimension ypthicm(nlon,norec,nlev), ypthic(nlon,norec,nlev)
C----------------------------------------------------------------------
      character clabel*(*), cerr*8
      logical zprint
      data zvapre1,zvapre,zheatpre1,ztotpre1,zplutpre,zpluqpre,
     *     zyt2pre,instepre /7*0., 0/
      data zprint /.false./
C----------------------------------------------------------------------
 
      if (.not. zprint) return
 
      zmass1= 0.
      zvap1 = 0.
      zheat1= 0.
      zke1= 0.
      zpot1= 0.
 
      zplut = 0.
      zpluq = 0.
 
      zmass = 0.
      zvap  = 0.
      zheat = 0.
      zplut = 0.
      zke = 0.
      zpot = 0.
 
      zyt2  = 0.
      zprec = 0.
      zevap = 0.
      zdiab = 0.
 
      do 100 jk = 1,nlev
        do 102 jj = 1,norec
          zwei = cosbud(jj)/nlon
          do 104 ji = 1,nlon
 
            zlwc = ylwc_c(ji,jj,jk)+ylwc_a(ji,jj,jk)+ylwc_s(ji,jj,jk)
 
            zm1= exp(ylpsm1(ji,jj))*dsigma(jk)/gravit
            zq1= yqm1(ji,jj,jk)
            zt1= ytm1(ji,jj,jk)
            zu1= yum1(ji,jj,jk) + fplumu(ji,jj,jk)/zm1
            zv1= yvm1(ji,jj,jk) + fplumv(ji,jj,jk)/zm1
 
            zm = exp(ylps(ji,jj))*dsigma(jk)/gravit
            zq = yq(ji,jj,jk)
            zt = yt(ji,jj,jk)
            zu = yu(ji,jj,jk) + fplumu(ji,jj,jk)/zm
            zv = yv(ji,jj,jk) + fplumv(ji,jj,jk)/zm
 
            if (hybtype.eq.0) then
              zme1= zm1
              zme = zm
            else
              zme1= ypthicm(ji,jj,jk)/gravit
              zme = ypthic(ji,jj,jk)/gravit
            endif
 
            zmass1= zmass1+ zwei*zm1
            zvap1 = zvap1 + zwei*(zq1+zlwc)*zme1
            zheat1= zheat1 + zwei
     *              * (   (1-zq1)*cpair*zt1*zm1
     *                  + zq1*(cpwv*(zt1-tmelt)+latsub)*zme1
     *                  + zlwc*latice*zme1 )
            zke1  = zke1 + zwei*0.5*(zu1**2+zv1**2)*zm1
 
            zpluq = zpluq + zwei*fplumq(ji,jj,jk)
            zplut = zplut + zwei
     *              * (  fplumt(ji,jj,jk)
     *                 + fplumq(ji,jj,jk)*(cpwv*(zt1-tmelt)+latsub) )
 
            zmass = zmass + zwei*zm
            zvap  = zvap  + zwei*(zq+zlwc)*zme
            zheat = zheat + zwei
     *              * (   (1.-zq)*cpair*zt*zm
     *                  + zq*(cpwv*(zt -tmelt)+latsub)*zme
     *                  + zlwc*latice*zme )
            zke= zke + zwei*0.5*(zu**2+zv**2)*zm
 
            if (nstep.gt.nrstrt) then
c             Subtract frictional heating from yt2: this was added to
c             yt2 in vdif, but cancels in *total* energy budget
              zyt2  = zyt2  + zwei * yt2(ji,jj,jk)
     *                      - zwei * zm *( -yfu(ji,jj,jk)*yu(ji,jj,jk)
     *                                     -yfv(ji,jj,jk)*yv(ji,jj,jk) )
              zdiab = zdiab + zwei*( yqrs(ji,jj,jk)-yqrs(ji,jj,jk+1)
     *                              +yqrl(ji,jj,jk)-yqrl(ji,jj,jk+1) )
            endif
  104     continue
  102   continue
  100 continue
 
      do 110 jj = 1,norec
        zwei = cosbud(jj)/nlon
        do 112 ji=1,nlon
          zpot1= zpot1+ zwei * exp(ylpsm1(ji,jj))*yfi(ji,jj) / gravit
          zpot = zpot + zwei * exp(ylps(ji,jj))  *yfi(ji,jj) / gravit
          if (nstep.gt.nrstrt) then
            zprec = zprec + zwei * ( frainc(ji,jj)+fsnowc(ji,jj)
     *                              +frains(ji,jj)+fsnows(ji,jj))
            zevap = zevap + zwei * yqflx(ji,jj)
            zdiab = zdiab + zwei *
     *            (   ytflx(ji,jj)
     *              + yqflx(ji,jj)
     *                *((yt(ji,jj,nlev)-tmelt)*cpwv + latsub)
     *              - (frainc(ji,jj)+frains(ji,jj))
     *                *((max(yt(ji,jj,nlev),tmelt)-tmelt)*ch2o+latice)
     *              - (fsnowc(ji,jj)+fsnows(ji,jj))
     *                *((min(yt(ji,jj,nlev),tmelt)-tmelt)*cice)
     *            )
          endif
  112   continue
  110 continue
 
      z1= zmass1*cpair
      z = zmass *cpair
      ztot1= zheat1+ zke1+ zpot1
      ztot = zheat + zke + zpot
 
      if (clabel.eq.'postphys') then
        write(cerr,'(f8.4)')
     *  ((zvap1-zvapre1)/twodt + (zpluq-zpluqpre)/dtime + (zprec-zevap))
     *  * 86400.
c       write(cerr,'(f8.4)')
c    *  (zheat1-zheatpre1)/twodt + (zplut-zplutpre)/dtime + zyt2 - zdiab
 
      else if (clabel.eq.'postspec') then
        write(cerr,'(f8.4)') (ztot-ztotpre1)/twodt - zyt2pre
 
      else
        cerr = ' '
      endif
 
      if (instepre.ne.nstep) write (*,8000) nstep
      write(*,8002)
     *  clabel(1:8),
     *  zvap1, zvap, (zheat1+zpot1)/z1, (zheat+zpot)/z,
     *  zke1/z1, zke/z, ztot1/z1, ztot/z,
     *  zyt2, zprec*dtime, zevap*dtime,
     *  cerr
 
 8000 format(/'ns=',i6,
     *        '    zvap1     zvap   zheat1    zheat',
     *        '     zke1      zke    ztot1     ztot',
     *        '     zyt2    zprec    zevap',
     *        '     err')
 8002 format(a8,':', 11f9.4, a)
 
      zvapre1 = zvap1
      zvapre = zvap
      zheatpre1 = zheat1
      zpluqpre = zpluq
      zplutpre = zplut
      zheatpre1 = zheat1
      ztotpre1 = ztot1
      zyt2pre = zyt2
      instepre = nstep
 
      return
      end
C
C **************************** PRINVERV ********************************
C
      subroutine prinverv (field, clabel)
 
c     Prints out mean and standard deviation of a 3-D field at
c     individual levels.
 
c     field  = 3-D field
c     clabel = field name
 
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMBUD.h>
C----------------------------------------------------------------------
#include <COMCON.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMMAP.h>
C----------------------------------------------------------------------
#include <COMMUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      dimension field(nlon,norec,nlev)
      character clabel*(*)
      dimension zav(nlev), zsd(nlev), zfac(nlev)
      logical zprint
      data zprint, izstepre /.false., 0/
C----------------------------------------------------------------------
 
      if (.not. zprint) return
 
      n = max (nint(0.25*86400./dtime), 1)
      if (mod(nstep,n).ne.0) return
 
c        Compute means
 
      do 100 jk = 1,nlev
        zav(jk) = 0.
        do 102 jj = 1,norec
          zwei = cosbud(jj)/nlon
          do 104 ji = 1,nlon
            zav(jk) = zav(jk) + field(ji,jj,jk)*zwei
  104     continue
  102   continue
 
c        Compute standard deviations
 
        zsd(jk) = 0.
        do 112 jj = 1,norec
          zwei = cosbud(jj)/nlon
          do 114 ji = 1,nlon
            zsd(jk) = zsd(jk) + ((field(ji,jj,jk)-zav(jk))**2)*zwei
  114     continue
  112   continue
        zsd(jk) = sqrt(zsd(jk))
  100 continue
 
c       Printout
 
      do 700 jk=1,nlev
        if (clabel(1:12).eq.'d[ln(ps)]/dt') then
          zfac(jk) = .864e5 * sig(jk)
        else if (clabel(1:5).eq.'omega') then
          zfac(jk) = .864e5 / 1.e5
        else
          zfac(jk) = .864e5
        endif
  700 continue
 
      if (nstep.ne.izstepre) then
        write(nout,705) nstep, (sig(jk),jk=1,12)
  705   format(80('-')/'nstep=',i6/12f11.3)
        izstepre = nstep
      endif
 
      write(nout,710)
     *  clabel,
     *  (zfac(jk)*zav(jk), jk=1,12),
     *  (zfac(jk)*zsd(jk), jk=1,12),
     *  (dsigma(jk)/(zfac(jk)*zsd(jk)), jk=1,12)
 
  710 format(
     *   'field=',a,'    mean,sd,dsigma/sd for k=1-12:'
     *  /12f11.5
     *  /12f11.5
     *  /12f11.5 )
 
      return
      end
c
c ***************************** equilr_i ******************************
c
c#ifdef isotrac
cc     subroutine equilr_i (T, qv, qv_o, qv_h, R_o, R_h)
cc
cc   Used to get isotopic ratios in isotopic equilibrium with the vapor
cc      Input:  T               temperature
cc              qv, qv_o, qv_h  water vapor specific humidities
cc      Output: R_o, R_h        equilibrium isotopic ratio
cc---------
c#include <COMCON.h>  ! smoj_*
c#include <COMFLG.h>  ! t_r_i
c#include <PARA_I.h>  ! av*_*()
cc---------
cc
c        if (T.ge.t_r_i) then
c          a_o = avl_o(T)
c          a_h = avl_h(T)
c        else
c          a_o = avi_o(T)
c          a_h = avi_h(T)
c        endif
cc
c        if (qv.ne.0.) then
c          R_o = a_o*qv_o/qv
c          R_h = a_h*qv_h/qv
cc       else
cc         R_o = 0.96*smoj_o
cc         R_h = 0.69*smoj_h
c        endif
cc
c      return
c      end
c#endif
c
c ***************************** prout_i *******************************
c
c>>>>>>>>>>>>> 
#ifdef isotrac
c>>>>>>>>>>>>> 

      subroutine prout_i (ypthicm, ypthic, titre)

#include <PARAMS.h>
#include <COMBUD.h>
#include <COMCON.h>
#include <COMLUN.h>
#include <COML10.h>
#include <COML30.h>
#include <COML50.h>
#include <COMMUN.h>
#include <COMTIM.h>
 
      dimension ypthicm(nlon,norec,nlev), ypthic(nlon,norec,nlev),
     ~          c(nrtr), zvap1_t(nrtr), zvap_t(nrtr), zpluq_t(nrtr),
     ~          zprec_t(nrtr), zevap_t(nrtr)
      character titre*10

      if (pii_i.eq.0 .or. pjj_i.eq.0) return

      iup = 80        ! logical unit nb for isotopes point print out
      iub = 81        ! logical unit nb for isotopes budget print out
      iut = 84        ! logical unit nb for tracers point print out
      iug = 98        ! logical unit nb for tracers budget print out
      ii  = pii_i     ! longitude indice
      jj  = pjj_i     ! latitude  indice

c           pdn_i     ! output every pdn_i timesteps
 
c     if (nstep.eq.nrstrt) then
      if (nstep.eq.nrstrt .and. titre.eq.'prePHYS   ') then
        open (iup, form='formatted', status='unknown', file='pnt_iso')
        open (iut, form='formatted', status='unknown', file='pnt_tra')
        open (iub, form='formatted', status='unknown', file='bud_iso')
        open (iug, form='formatted', status='unknown', file='bud_tra')
        write (iub,'(a)')  '^^Step ^Title^^^^ ^^^M-1(vap+liq) '//
     ~ '^^^^^M(vap+liq) ^^^^^^^^^^^^^Fv ^^^^^^^^^^^M(P) ^^^^^^^^^^^M(E)'
        write (iug,'(a)')  '^^Step ^Title^^^^ ^^^M-1(vap+liq) '//
     ~ '^^^^^M(vap+liq) ^^^^^^^^^^^^^Fv ^^^^^^^^^^^M(P) ^^^^^^^^^^^M(E)'
      endif
 1    format('### [',i2,',',i2,'] # ',a10,' # Step ',i6,' ###')
 2    format (i6,1x,a10,             5(1x,g15.8))
 3    format (6x,1x,a10,             5(1x,g15.8))
 4    format (6x,1x,'      t',i2,':',5(1x,g15.8))
 
c YQ<0 or YQM1 or TRACER>=1.01 or TRACER<=0.01 TRAP
c5    format (i6,a10,'>yq',a2,'(',i2,2(',',i2),')=',g15.8,
c    ~                                      '>c(',i2,')=',f7.2,'%')
c     write (85,*) 'step ',nstep
c#ifdef isotracreg
c     do i=1,nlon
c      do j=1,norec
c       do k=1,nlev
c         if (yq  (i,j,k).ne.0.) then
c           do it=1,nrtr
c             ct = rtrq  (i,j,k,it)/yq  (i,j,k)
c             if (yq  (i,j,k).lt.0. .or. ct.le.-0.01 .or. ct.ge.1.01)
c    ~        write (85,5) nstep,titre,'  ',i,j,k,yq  (i,j,k),it,ct*100.
c           enddo
c         endif
c         if (yqm1(i,j,k).ne.0.) then
c           do it=1,nrtr
c             ct = rtrqm1(i,j,k,it)/yqm1(i,j,k)
c             if (yqm1(i,j,k).lt.0. .or. ct.le.-0.01 .or. ct.ge.1.01)
c    ~        write (85,5) nstep,titre,'m1',i,j,k,yqm1(i,j,k),it,ct*100.
c           enddo
c         endif
c       enddo
c      enddo
c     enddo
c#endif
 
      if (nstep.ne.nrstrt.and.mod(nstep,pdn_i).ne.0) return
 
c ISOTOPES POINT PRINT OUT
 
      write (iup,1) ii,jj,titre,nstep
      write (iup,'(a194)')
     ~ '^k ^^^^^^^^^^^^^yq ^^^d18O ^^^^d2H ^^^^^^^^^^^yqm1 ^^^d18O '//
     ~ '^^^^d2H ^^^^^^^^^fplumq ^^^d18O ^^^^d2H ^^^^^^^^^ylwc_c '//
     ~ '^^^d18O ^^^^d2H ^^^^^^^^^ylwc_a ^^^d18O ^^^^d2H '//
     ~ '^^^^^^^^^ylwc_s ^^^d18O ^^^^d2H'
 
      do 10 k=1,nlev
        if (yq    (ii,jj,k).ne.0.) then
          doyq = 1e3*(yq    _o(ii,jj,k)/yq    (ii,jj,k)/smoj_o-1.)
          ddyq = 1e3*(yq    _h(ii,jj,k)/yq    (ii,jj,k)/smoj_h-1.)
        else
          doyq = -999.99
          ddyq = -999.99
        endif
        if (yqm1  (ii,jj,k).ne.0.) then
          doy1 = 1e3*(yqm1  _o(ii,jj,k)/yqm1  (ii,jj,k)/smoj_o-1.)
          ddy1 = 1e3*(yqm1  _h(ii,jj,k)/yqm1  (ii,jj,k)/smoj_h-1.)
        else
          doy1 = -999.99
          ddy1 = -999.99
        endif
        if (fplumq(ii,jj,k).ne.0.) then
          dofq = 1e3*(fplumq_o(ii,jj,k)/fplumq(ii,jj,k)/smoj_o-1.)
          ddfq = 1e3*(fplumq_h(ii,jj,k)/fplumq(ii,jj,k)/smoj_h-1.)
        else
          dofq = -999.99
          ddfq = -999.99
        endif
        if (ylwc_c(ii,jj,k).ne.0.) then
          dolc = 1e3*(ylwc_c_o(ii,jj,k)/ylwc_c(ii,jj,k)/smoj_o-1.)
          ddlc = 1e3*(ylwc_c_h(ii,jj,k)/ylwc_c(ii,jj,k)/smoj_h-1.)
        else
          dolc = -999.99
          ddlc = -999.99
        endif
        if (ylwc_a(ii,jj,k).ne.0.) then
          dola = 1e3*(ylwc_a_o(ii,jj,k)/ylwc_a(ii,jj,k)/smoj_o-1.)
          ddla = 1e3*(ylwc_a_h(ii,jj,k)/ylwc_a(ii,jj,k)/smoj_h-1.)
        else
          dola = -999.99
          ddla = -999.99
        endif
        if (ylwc_s(ii,jj,k).ne.0.) then
          dols = 1e3*(ylwc_s_o(ii,jj,k)/ylwc_s(ii,jj,k)/smoj_o-1.)
          ddls = 1e3*(ylwc_s_h(ii,jj,k)/ylwc_s(ii,jj,k)/smoj_h-1.)
        else
          dols = -999.99
          ddls = -999.99
        endif
       write (iup,'(i2,6(1x,g15.8,2(1x,f7.2)))')
     ~  k, yq    (ii,jj,k), doyq, ddyq, yqm1  (ii,jj,k), doy1, ddy1,
     ~      fplumq(ii,jj,k), dofq, ddfq, ylwc_c(ii,jj,k), dolc, ddlc,
     ~      ylwc_a(ii,jj,k), dola, ddla, ylwc_s(ii,jj,k), dols, ddls
  10  continue
 
c TRACERS POINT PRINT OUT
 
#ifdef isotracreg
      write (iut,1) ii,jj,titre,nstep
      write (iut,'(a13,$)') '^k ^^^^^^^^yq'
      do it=1,nrtr
        write (iut,'(a3,i2,a2,$)') ' c(',it,')%'
      enddo
      write (iut,'(a11,$)') ' ^^^^fplumq'
      do it=1,nrtr
        write (iut,'(a3,i2,a2,$)') ' c(',it,')%'
      enddo
      write (iut,'()')
 
      do k=1,nlev
        write (iut,'(i2,1x,g10.4,$)')    k, yq    (ii,jj,k)
        do it=1,nrtr
          if (yq    (ii,jj,k).ne.0.) then
            c(it) = 100.*rtrq (ii,jj,k,it)/yq    (ii,jj,k)
          else
            c(it) = -99.99
          endif
          write (iut,'(1x,f6.2,$)') c(it)
        enddo
        write (iut,'(1x,g10.4,$)')           fplumq(ii,jj,k)
        do it=1,nrtr
          if (fplumq(ii,jj,k).ne.0.) then
            c(it) = 100.*rtrfq(ii,jj,k,it)/fplumq(ii,jj,k)
          else
            c(it) = -99.99
          endif
          write (iut,'(1x,f6.2,$)') c(it)
        enddo
        write (iut,'()')
      enddo
 
      if (titre.eq.'postPHYS  ') then
        write (iut,'(a13,$)') '   ^^^^^^^^^P'
        do it=1,nrtr
          write (iut,'(a3,i2,a2,$)') ' c(',it,')%'
        enddo
        write (iut,'()')
 
        p = frainc(ii,jj)+frains(ii,jj)+fsnowc(ii,jj)+fsnows(ii,jj)
        write (iut,'(3x,g10.4,$)') p
        do it=1,nrtr
          if (p.ne.0.) then
            c(it) = 100.*(rtrprc(ii,jj,it)+rtrprs(ii,jj,it)+
     ~                    rtrpsc(ii,jj,it)+rtrpss(ii,jj,it))/p
          else
            c(it) = -99.99
          endif
          write (iut,'(1x,f6.2,$)') c(it)
        enddo
        write (iut,'()')
      endif
#endif
 
c ISOTOPES & TRACERS BUDGET PRINT OUT
 
      zvap1   = 0.
      zvap1_o = 0.
      zvap1_h = 0.
      zvap    = 0.
      zvap_o  = 0.
      zvap_h  = 0.
      zpluq   = 0.
      zpluq_o = 0.
      zpluq_h = 0.
      zprec   = 0.
      zprec_o = 0.
      zprec_h = 0.
      zevap   = 0.
      zevap_o = 0.
      zevap_h = 0.
#ifdef isotracreg
      do it=1,nrtr
        zvap1_t(it) = 0.
        zvap_t(it)  = 0.
        zpluq_t(it) = 0.
        zprec_t(it) = 0.
        zevap_t(it) = 0.
      enddo
#endif
 
      do 100 k = 1,nlev
        do 102 j = 1,norec
          zwei = cosbud(j)/nlon
          do 104 i = 1,nlon
            zme1   = ypthicm(i,j,k)/gravit
            zme    = ypthic(i,j,k) /gravit
            zlwc   = ylwc_c(i,j,k)  +ylwc_a(i,j,k)+
     ~               ylwc_s(i,j,k)
            zlwc_o = ylwc_c_o(i,j,k)+ylwc_a_o(i,j,k)+
     ~               ylwc_s_o(i,j,k)
            zlwc_h = ylwc_c_h(i,j,k)+ylwc_a_h(i,j,k)+
     ~               ylwc_s_h(i,j,k)
            zvap1  = zvap1  + zwei*(yqm1(i,j,k)  +zlwc)  *zme1
            zvap1_o= zvap1_o+ zwei*(yqm1_o(i,j,k)+zlwc_o)*zme1
            zvap1_h= zvap1_h+ zwei*(yqm1_h(i,j,k)+zlwc_h)*zme1
            zvap   = zvap   + zwei*(yq(i,j,k)    +zlwc)  *zme
            zvap_o = zvap_o + zwei*(yq_o(i,j,k)  +zlwc_o)*zme
            zvap_h = zvap_h + zwei*(yq_h(i,j,k)  +zlwc_h)*zme
            zpluq  = zpluq  + zwei*fplumq(i,j,k)
            zpluq_o= zpluq_o+ zwei*fplumq_o(i,j,k)
            zpluq_h= zpluq_h+ zwei*fplumq_h(i,j,k)
            if ((k.eq.nlev).and.(nstep.gt.nrstrt)) then
              zprec  = zprec  + zwei*(frainc(i,j)  +fsnowc(i,j)
     ~                               +frains(i,j)  +fsnows(i,j))
              zprec_o= zprec_o+ zwei*(frainc_o(i,j)+fsnowc_o(i,j)
     ~                               +frains_o(i,j)+fsnows_o(i,j))
              zprec_h= zprec_h+ zwei*(frainc_h(i,j)+fsnowc_h(i,j)
     ~                               +frains_h(i,j)+fsnows_h(i,j))
              zevap  = zevap  + zwei*yqflx(i,j)
              zevap_o= zevap_o+ zwei*yqflx_o(i,j)
              zevap_h= zevap_h+ zwei*yqflx_h(i,j)
            endif
  104     continue
  102   continue
  100 continue

#ifdef isotracreg
      do it=1,nrtr
        do k = 1,nlev
          do j = 1,norec
            zwei = cosbud(j)/nlon
            do i = 1,nlon
              zme1   = ypthicm(i,j,k)/gravit
              zme    = ypthic(i,j,k) /gravit
              zlwc   = ylwc_c(i,j,k)  +ylwc_a(i,j,k)+
     ~                 ylwc_s(i,j,k)
              if (yqm1(i,j,k).ne.0.) zvap1_t(it) = zvap1_t(it)+
     ~            zwei*rtrqm1(i,j,k,it)*(1.+zlwc/yqm1(i,j,k))*zme1
              if (yq(i,j,k)  .ne.0.) zvap_t(it)  = zvap_t(it) +
     ~            zwei*rtrq(i,j,k,it)  *(1.+zlwc/yq(i,j,k))  *zme
              zpluq_t(it) = zpluq_t(it)+zwei*rtrfq(i,j,k,it)
              if ((k.eq.nlev).and.(nstep.gt.nrstrt)) then
                zprec_t(it) = zprec_t(it)+zwei*
     ~                        (rtrprs(i,j,it)+rtrprc(i,j,it)
     ~                        +rtrpss(i,j,it)+rtrpsc(i,j,it))
                zevap_t(it) = zevap_t(it)+zwei*rtrfv(i,j,it)
              endif
            enddo
          enddo
        enddo
      enddo
#endif
 
      write (iub,2)   nstep,titre,
     ~       zvap1,      zvap,      zpluq,      zprec,      zevap
      write (iub,3)   '     18-O:',
     ~       zvap1_o,    zvap_o,    zpluq_o,    zprec_o,    zevap_o
      write (iub,3)   '      2-H:',
     ~       zvap1_h,    zvap_h,    zpluq_h,    zprec_h,    zevap_h
 
      write (iug,2)   nstep,titre,
     ~       zvap1,      zvap,      zpluq,      zprec,      zevap
#ifdef isotracreg
      do it=1,nrtr
        write (iug,4)          it,
     ~       zvap1_t(it),zvap_t(it),zpluq_t(it),zprec_t(it),zevap_t(it)
      enddo
#endif
 
c     ii = flush(iup)            !   ### for debug ... ###
c     ii = flush(iub)            !   ### for debug ... ###
c     ii = flush(iut)            !   ### for debug ... ###
c     ii = flush(iug)            !   ### for debug ... ###
c     ii = flush(85 )            !   ### for debug ... ###

      return
      end

c>>>>>
#endif
c>>>>>
c
c ***************************** hand1,2,3 *****************************
c
#ifdef sun
      integer function hand1()
      write (6,*) '*** STOP (exception: division) in hand1'
      hand1 = 0.
      call endrun (-1)
      end

      integer function hand2()
      write (6,*) '*** STOP (exception: overflow) in hand2'
      hand2 = 0.
      call endrun (-1)
      end

      integer function hand3()
      write (6,*) '*** STOP (exception: invalid) in hand3'
      hand3 = 0.
      call endrun (-1)
      end
#endif

