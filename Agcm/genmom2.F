c for Genesis v3.0

c nb: Test on Cray with indefs
c nb: Does mom use multiprocs? 

c*********************
#if defined (coupogcm)
c*********************

      program genmom2

c     This is the main driver for coupled A/OGCM.
c     It manages time stepping and coupling between GENESIS and MOM2.

c----------------------
c GENESIS common blocks
c----------------------
#include <compar.h>
#include <COMDYR.h>
#include <COMHIS.h>
#include <COMLUN.h>
#include <COMTIM.h>
#include <comicexy.h>
#include <commapl.h>
#include <commugxy.h>
#include <comocnxy.h>

c compar.h    nlon, nlat
c COMDYR.h    calday, jyear
c COMLUN.h    nout, iusurt, locname
c COMTIM.h    dtime, nrstrt, nstep, nstop, istepy
c commapl.h   alon, alat
c commugxy.h  lmaskxy, and for isot: rainaxy*, snowaxy*, avapaxy*
c comocnxy.h  fuocnet,fvocnet,ftocnxy,fwocnxy,fsolocnxy
c             uocnxy,  vocnxy, tocnxy, socnxy,  
c             and for isot: ro_*, ri_*
c comicexy.h  for isot: ficexy

c------------------
c MOM common blocks 
c------------------

#include <size.h>
#include <csbc.h>
#include <levind.h>

c size.h      imt, jmt, km, nt (needed only for csbc.h) 
c csbc.h      sbcatm, aland (nb: aland is integer) 
c             (includes aparam.h, for ima, imap2, jma)
c levind.h    kmt, kmu

c---------------

c  GENESIS SBCs (in comocnxy.h, except foocnxy,fdocnxy are local)
c     All net, *below* sea ice)
c     fuocnet  = lsx eastward wind stress on ocn, incl under ice(N/m**2)
c     fvocnet  = lsx northwrd wind stress on ocn, incl under ice(N/m**2)
c                (f[u,v]ocnet are modified from f[u,v]ocnxy in dynice,
c                 and shifted to NE corner rel to h-grid in surfctl)
c     ftocnxy  = lsx net heat flux into ocean, incl under ice (W/m**2)
c     fwocnxy  = lsx net fresh water flux P-E-d[ice], " " (kg/m**2/sec)
c     fsolocnxy= lsx solar flux absorbed by ocean (W/m**2)
c     foocnxy  = lsx net dO18 flux (permil * kg/m**2/sec)
c     fdocnxy  = lsx net dD   flux (permil * kg/m**2/sec)
c     tocnxy   = surface ocean temperature (K)
c     socnxy   = surface ocean salinity (fraction)
c     uocnxy   = eastward surface ocean velocity (m/sec)
c     vocnxy   = westward surface ocean velocity (m/sec)
c     ro_o      = surface ocean dO18/O16 (kg/kg)
c     ro_h      = surface ocean dD/H     (kg/kg)

c  MOM SBCs (on imap2,jma grid. imap2=ima+2,ima=nlon,jma=nlat) (csbc.h)
c  For 2x2, ima=180, jma=90, and sbcatm(imap2,jma,*), aland(imap2,jma), 
c  running from -1 to 361E, -89 to 89N.
c     sbcatm(1)  = taux (dynes/cm**2), NE corner rel to h grid
c     sbcatm(2)  = tauy (dynes/cm**2), NE corner rel to h grid
c     sbcatm(3)  = net heat flux (cal/cm**2/sec) (+ve if ocn warming) 
c     sbcatm(4)  = salt flux (frac*gm/cm**2/sec) (+ve if ocn more salty)
c     sbcatm(5)  = ocean absorbed shortwave (cal/cm**2/sec) 
c     sbcatm(6)  = dO18 flux (permil*g/cm**2/sec) 
c     sbcatm(7)  = dD flux (permil*g/cm**2/sec) 

c     sbcatm(8)  = sst (deg C)
c     sbcatm(9)  = sss (salinity, fraction  - .035)
c     sbcatm(10) = ss dO18 (permil)
c     sbcatm(11) = ss dD   (permil)
c     sbcatm(12) = surface u (cm/sec), NE corner rel to h grid
c     sbcatm(13) = surface v (cm/sec), NE corner rel to h grid
c     aland()    = land/ocean map (0=land, 1=ocean) on mom2 imap2,jma

c======
      parameter (pie=3.14159265358979)  ! pi in ndcon.h, from scalar.h
      dimension workl(nlon,nlat)

      logical ifexist, testhisl

      dimension 
     *  fuocnav(nlon,nlat), fvocnav(nlon,nlat),
     *  ftocnav(nlon,nlat), fwocnav(nlon,nlat), fsolocnav(nlon,nlat)
      save fuocnav, fvocnav, ftocnav, fwocnav, fsolocnav 
#ifdef isotrac
      dimension
     *  foocnxy(nlon,nlat), fdocnxy(nlon,nlat),
     *  foocnav(nlon,nlat), fdocnav(nlon,nlat)
      save foocnav, fdocnav
#endif

      dimension fwocnadj(nlon,nlat), runofadj(nlon,nlat),
     *          calvadj(nlon,nlat)

      dimension maskogcmh(nlon,nlat), maskogcmu(nlon,nlat),
     *          maskbothh(nlon,nlat), maskbothu(nlon,nlat)
      save maskogcmh, maskogcmu, maskbothh, maskbothu, odepth

c       Time scale, days, of exponential filters for ogcm heat,
c       h2o fluxes and wind stresses. Used in ocnaccum and below.
c       To not use any filter, set tocfilt to 0.
c     parameter (tocfilt = 0.)
      parameter (tocfilt = 5.)

#if defined (asyncmode)
      character*80 async_path, async_startdate, async_enddate
      character*96 async_surftype_file
      character*8  versasync
#if   defined (aix)
      character*8 idate,date
      character*8 itime,clock_
#elif defined (sgi)
      character*9 idate
      character*8 itime
#elif defined (sun)
      character*24 jdate
#elif defined (linux)
      character*8  kdate
      character*10 ktime
#endif
#endif

c     calendar days for ogcm history and restart writes, 
c     used in setting fwhisto, fsaveo below.

      parameter (nocmax=1000)
      dimension nochis(nocmax), nocsav(nocmax)

c     data (nochis(i),i=1,4) /1, 2, 3, 4/
c     data (nocsav(i),i=1,2) /2, 4/

c     data (nochis(i),i=1,2) 
c    *  /                      181,                          365/
      data (nochis(i),i=1,4)
     *  /        90,           181,           273,           365/
c     data (nochis(i),i=1,12)
c    *  /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/

      data (nocsav(i),i=1,1) /365/

c     only do ogcm history/restart writes once every nocyear years

c     data nocyear /1/  
      data nocyear /10/

      logical presc_rivers , presc_calvant
      data presc_rivers, presc_calvant /.true., .true./

      logical fwhisto, fsaveo
      save fwhisto, fsaveo, nochis, nochisn, nocsav, nocsavn, nocyear

      namelist /gencoup/ 
     *  nochis, nocsav, nocyear, presc_rivers, presc_calvant, 
#if defined (asyncmode)
     *  async_path, async_startdate, async_enddate, async_surftype_file,
#endif
     *  o18init

#if defined (asyncmode)
      dimension alatt31(48)
      data alatt31 /
     *  -1.5212131731, -1.4569823009, -1.3923721142, -1.3276762287,
     *  -1.2629478409, -1.1982037311, -1.1334508412, -1.0686925583,
     *  -1.0039307308, -0.9391664533, -0.8744004155, -0.8096330743,
     *  -0.7448647446, -0.6800956510, -0.6153259588, -0.5505557923,
     *  -0.4857852479, -0.4210144016, -0.3562433151, -0.2914720392,
     *  -0.2267006173, -0.1619290868, -0.0971574813, -0.0323858319,
     *   0.0323858319,  0.0971574813,  0.1619290868,  0.2267006173,
     *   0.2914720392,  0.3562433151,  0.4210144016,  0.4857852479,
     *   0.5505557923,  0.6153259588,  0.6800956510,  0.7448647446,
     *   0.8096330743,  0.8744004155,  0.9391664533,  1.0039307308,
     *   1.0686925583,  1.1334508412,  1.1982037311,  1.2629478409,
     *   1.3276762287,  1.3923721142,  1.4569823009,  1.5212131731 /
#endif

#ifdef isotrac
      save smoj_o, smoj_h, 
     *     o18std, d18std, o18riv, d18riv,
     *     ofstd,  dfstd,  ofriv,  dfriv
#endif

c----------------------------------------------------------------
c     G E N E S I S   A N D   M O M   I N I T I A L I Z A T I O N
c----------------------------------------------------------------

c         Namelist defaults (if not in data stmts above)

      do i=1,nocmax
        nochis(i) = 0
        nocsav(i) = 0
      enddo
      nocyear = 1
#if defined (asyncmode)
      async_path      = ' '
      async_startdate = ' '
      async_enddate   = ' '
      async_surftype_file = 'SURFTYP'
      nout = 6
#else
      nout = 6
c     nout = 7
#endif

#ifdef isotrac
        smoj_o = 2.2149e-3
        smoj_h = 3.2883e-4
c       smoj_o = 1.9934e-3 * (20./18.)
c       smoj_h = 0.1557e-3 * (19./18.) * 2.

        o18std  = 0.        ! std d18O ocean value for flux calcs (smow)
        d18std  = 7.*o18std ! (as in Mathieu et al.,2002)

        o18riv  = o18std ! riverine input (fn of lat, also for calvant?)
        d18riv  = d18std ! riverine input (fn of lat, also for calvant?)

        ofstd = (o18std/1000. + 1.) * smoj_o    ! smow to fraction
        dfstd = (d18std/1000. + 1.) * smoj_h    ! snow to fraction
        ofriv = (o18riv/1000. + 1.) * smoj_o    ! smow to fraction
        dfriv = (d18riv/1000. + 1.) * smoj_h    ! snow to fraction
#endif
c       Namelist: initial d18O ogcm (3-D) value (smow).
c       eg, 0 for modern, -1.2 for Cretaceous (default=-999., not set)
c       Used in setocn.F (init=.T.), ot read_restart (init=.F.)
        o18init = -999.

c        Read namelist ogcm history/restart write times, 
c        and if async: start and end dates, path to driver files,
c        lsx surface data file name

      write (nout,'(/a/)')
     *  'Reading namelist file interm_coup (gencoup)...'
      iunamel = 71
      open (iunamel, file='interm_coup', status='old',err=1000)
      read (iunamel,gencoup)
      close (iunamel)
      write (nout,*)

c        Calculate numbers of ogcm history, restart times

      nochisn = 0
      do i=1,nocmax
        if (nochis(i).ne.0) nochisn = i
      enddo
      nocsavn = 0
      do i=1,nocmax
        if (nocsav(i).ne.0) nocsavn = i
      enddo

c        If just one value entered,  convert to multiples of that
c        time from beginning of year

      if (nochisn.eq.1) then
        nochisn = 365/nochis(1) 
        if (nochisn.ge.2) then 
          do i=2,nochisn
            nochis(i) = nochis(1)*i
          enddo
        endif
      endif

      if (nocsavn.eq.1) then
        nocsavn = 365/nocsav(1) 
        if (nocsavn.ge.2) then 
          do i=2,nocsavn
            nocsav(i) = nocsav(1)*i
          enddo
        endif
      endif

c~~~~~~~~~~~~~~~~~~~~~~
#if defined (asyncmode)
c~~~~~~~~~~~~~~~~~~~~~~
      versasync = '1.0'

c        Write banner (like subr banner)
 
#if   defined (aix)
      idate = date()
      itime = clock_()
#elif defined (cray)
      call date  (idate)
      call clock (itime)
#elif defined (sgi)
      call date (idate)
      call time (itime)
#elif defined (sun)
      call fdate (jdate)
#elif defined (linux)
      call date_and_time (kdate, ktime)
#endif


#if defined (cray) || defined (sgi) || defined (aix)
      write(nout,100) versasync, idate, itime
 100  format(
     *  /17x,'=============================================='
     *  /17x,'| Asynchronous Gen-Mom2    version: ',  a8,' |'
     *  /17x,'| date: ',a9,   '          time:    ',  a8,' |'
     *  )
#elif defined (sun) 
      write(nout,100) versasync, jdate
 100  format(
     *  /17x,'=============================================='
     *  /17x,'| Asynchronous Gen-Mom2    version: ',  a8,' |'
     *  /17x,'|          ',                 a24,'          |'
     *  )
#elif defined (linux)
        write (nout,100) versasync, kdate, ktime
 100    format(
     *    /17X,'=============================================='
     *    /17x,'| Asynchronous Gen-Mom2    version: ',  a8,' |'
     *    /17x,'| date: ',a8,  '           time:  ',   a10,' |'
     *    )
#else
        write (nout,100) versasync
 100    format(
     *    /17X,'=============================================='
     *    /17x,'| Asynchronous Gen-Mom2    version: ',  a8,' |'
     *  )
#endif
      write(nout,200) 
 200  format(
     *   17x,'=============================================='
     *  /)

c        Set lsx latitude grid

      do i=1,nlon
       alon(i) = 2.*pie*(i-.5)/nlon
      enddo
      if (nlat.eq.48) then
        call scopy (nlat, alatt31, 1, alat, 1)
      else
        do j=1,nlat
          alat(j)  = -pie/2. + pie*(j-.5)/nlat
        enddo
      endif

c        MOM2 initialization

      write (nout,'(/a/)') ' Calling init_driver (mom2)'

c     1st arg: agcm timestep, sec, (if sync), passed to mom2 (not used).
c     2nd arg: ogcm (stp_driver outer) timestep, days, passed from mom2.
c     3rd arg: ogcm upper layer thickness, m,  passed from mom2.
c     4th arg: ogcm first timestep (segtim) counter.

      call init_driver (0., segtim, laststep, odepth, o18init)

c        Use segtim (passed from init_driver) as overall time step

      dtime = segtim*86400.    

c        Calculate nrstrt, nstop from Namelist dates

      if (async_startdate.eq.' ' .or. async_enddate.eq.' ') then
        write (nout,*) 
     *             '*** Must specify ASYNC_STARTDATE and ASYNC_ENDDATE',
     *             ' for this asynchronous run (Namelist, interm_coup)'
        call endrun (-1)
      endif

c        Convert dates to YYMMDD (or leave as day numbers if < 100000) 
c        (mbdate, medate), and  corresponding timestep #s (nrstrt,nstop)
c        Input dates (async_startdate, async_enddate) can either be 
c        'YY/MM/DD', YYMMDD or day numbers (< 100000).
c        Nb:  nrstrt is always set to 1.

      call dateconv2 (mbdate, async_startdate, 'STARTDATE')
      call dateconv2 (medate, async_enddate,   'ENDDATE')

      if (iabs(mbdate).ge.100000 .and. iabs(medate).ge.100000) then
        call datedif (ndays, mbdate, medate)
        nstop = nint (ndays*86400./dtime)
      else
        nstop = nint ((medate-mbdate+1)*86400./dtime)
      endif

      nrstrt = nint (laststep*segtim*86400./dtime) + 1

c~~~~
#else
c~~~~

c        Genesis initialization     

      call basic

c        MOM2 initialization

      write (nout,'(/a)') ' Calling MOM2 init_driver'
      call init_driver (dtime, segtim, laststep, odepth, o18init)

      if (nrstrt .ne. nint(laststep*segtim*86400./dtime) + 1) then
        write (nout,'(/a/2(a,i10,a,f10.3/))')
     *                 '*** Error: a/ogcm (re)start time mismatch:',
     *                 '    nrstrt  =', nrstrt,  '   dtime =',dtime,
     *                 '    laststep=', laststep,'   segtim=',segtim
        stop
      endif

c~~~~~
#endif
c~~~~~

      write (nout,'(/a,2f10.2/)') 'MOM2: segtim=',segtim
      ier = flush (nout)

c        Check a/ogcm horizontal dimensions match

      if (nlon.ne.imap2-2 .or. nlat.ne.jma) then
        write (nout,300) nlon, nlat, imap2-2, jma
  300   format(/'*** Error: LSX/OGCM grid mismatch:'
     *         /'    LSX        nlon,    nlat =',2i6
     *         /'    OGCM       imap2-2, jma  =',2i6)
        call endrun (-1)
      endif

c        Read lsx surface-type map. Need to do so here, can't wait for 
c        synchronous read in surfctl/inisurf/inilsx), since need 
c        lmaskxy to set first state vars below

#if defined (asyncmode)       
      iu = 17
      open (iu, file=async_surftype_file, status='old', err=2000)
#else
      call findiu (iusurt, loc)
      call openf (iusurt, locname(loc), ' ', 0, 0)
      iu = iusurt
#endif
      call readdata (iu, workl, nlon, nlat, 'SURT    ', 1, 1)
      close (iu)
      do j=1,nlat
        do i=1,nlon
          lmaskxy(i,j) = nint(workl(i,j))
          if (lmaskxy(i,j).eq.0. or. lmaskxy(i,j).eq.3)
     *      lmaskxy(i,j) = 3
        enddo
      enddo

c        Set aland(imap2,jma), used by mom2 (csbc.h, gosbc, gasbc).
c        Actually not used by mom2 after dp mods (except for diagnostic
c        averages in gasbc.F).
c        Now only used only to display sbcatm in printmapocn below. 

      do j=1,jma
        do i=1,imap2
          aland(i,j) = min (kmt(i,j+1),1)
        enddo
      enddo

c        Set maskogcm[h,u], equiv to mom's km[t,u], but on lsx grids.
c        maskogcmh is used in river/calv/pme*, and maskogcm[h,u] for 
c        printmapocn and history fields. 
c
c        Also set maskbothh and maskbothu if ocean in *both* lsx 
c        and ogcm, on lsx grids. Used in setting [t,s,o,u,v]ocnxy 
c        from sbcatm, masks in async_mom and subsequent fill_near of
c        fluxes, and fill_near for [u,v]ocnxy.

      ndiff = 0
      do j=1,nlat
        jp = min (j+1,nlat)
        do i=1,nlon
          maskogcmh(i,j) = min (kmt(i+1,j+1), 1)
          maskogcmu(i,j) = min (kmu(i+1,j+1), 1)
#if defined  (identical_grids)
          if ( (lmaskxy(i,j).eq.3 .and. maskogcmh(i,j).ne.1) .or.
     *         (lmaskxy(i,j).ne.3 .and. maskogcmh(i,j).eq.1) ) then
            if (ndiff.eq.0) write(nout,'(a)')
     *                             '*** "Identical" Grids Error'
            write (nout,'(a,4i4)') '    i,j,lmaskxy,maskogcmh=',
     *                                  i,j,lmaskxy(i,j),maskogcmh(i,j)
            ndiff = ndiff + 1
          endif
#endif
          if (lmaskxy(i,j).eq.3 .and. maskogcmh(i,j).eq.1) then
            maskbothh(i,j) = 1
          else
            maskbothh(i,j) = 0
          endif

c         u grid:
          ip = mod (i,nlon) + 1
          if (lmaskxy(i,j)  .eq.3 .and. lmaskxy(ip,j) .eq.3 .and.
     *        lmaskxy(i,jp) .eq.3 .and. lmaskxy(ip,jp).eq.3 .and.
     *        maskogcmu(i,j).eq.1) then
            maskbothu(i,j) = 1
          else
            maskbothu(i,j) = 0
          endif
        enddo
      enddo
#if ! defined (asyncmode)
      if (ndiff.gt.0) stop
#endif

c     set genesis's state vars for first calc of sfc fluxes.
c     (sbcatm* set by ogcm (cold/restart) by init_driver).
c     nb: u,v velocs (sbc[atm,ocn](8,9)) are on mom2's u-grid
c     (NE corner of h grid...dp mod in subr clinic). 
c     Set only for pts that are ocean in *both* lsx and mom,
c     since (i) tocn, etc, over lsx land (lake) is indep of mom, and 
c     (ii) sbcatm values not meaningful over mom's land. 

      do j=1,nlat
        do i=1,nlon
          if (maskbothh(i,j).eq.1) then
            tocnxy(i,j) = sbcatm(i+1,j,8) + 273.16
            socnxy(i,j) = sbcatm(i+1,j,9) + .035 
#ifdef isotrac
            ro_o(i,j)   = (sbcatm(i+1,j,10)/1000. + 1.) *smoj_o
            ro_h(i,j)   = (sbcatm(i+1,j,11)/1000. + 1.) *smoj_h
#endif
#if defined (asyncmode)       
          else  ! just for niceity:
            tocnxy(i,j) = 0.
            socnxy(i,j) = 0.
#ifdef isotrac
            ro_o(i,j)   = 0.
            ro_h(i,j)   = 0.
#endif
#endif
          endif
          if (maskbothu(i,j).eq.1) then
            uocnxy(i,j) = sbcatm(i+1,j,12) / 100.
            vocnxy(i,j) = sbcatm(i+1,j,13) / 100.
#if defined (asyncmode)       
          else! just for niceity:
            uocnxy(i,j) = 0.
            vocnxy(i,j) = 0.
#endif
          endif
        enddo
      enddo

c        Set surface flux accumulators initially to zero

      call zero (fuocnav, nlon*nlat)
      call zero (fvocnav, nlon*nlat)
      call zero (ftocnav, nlon*nlat)
      call zero (fwocnav, nlon*nlat)
      call zero (fsolocnav, nlon*nlat)
#ifdef isotrac
      call zero (foocnav, nlon*nlat)
      call zero (fdocnav, nlon*nlat)
#endif

c        Set ogcm history, restart write flags initially to false

      fwhisto = .false.
      fsaveo = .false.

c------------------------------------
c     M A I N   L O O P   S T A R T S
c------------------------------------
c
c     mom2 timesteps dtts, dtuv are set in namelist tsteps. 
c     dtts = temp, sal timestep (basic mom timestep...dtuv,dtsf finer)
c     Two loops in stp_driver: 
c     Outer loop (n=1,numseg) includes updates to/from sfc bcs 
c        (arrays, sbc[atm,ocn], subrs gosbc, gasbc), and integrates
c        (i.e., inner loop iterates) through segtim (days, namelist), 
c        Outer loop is 1 to numseg, for total of runlen (days, 
c        namelist). For gen-mom, segtim=runlen=1 day, numseg=1.
c     Inner loop (loop=1,ntspos) steps subr mom once each iter,
c        which integrates through one dtts (sec, namelist) timestep.
c        So ntspos is # of dtts in one segtim.

      napero = nint (segtim*86400./dtime) ! # of lsx timesteps between
                                           ! ogcm (stp_driver) calls

c+++++++++++++++++++++++
      do is=nrstrt,nstop
c+++++++++++++++++++++++
        nstep=is

c----------------------------
c     C A L L   G E N E S I S
c----------------------------

c~~~~~~~~~~~~~~~~~~~~~~
#if defined (asyncmode)
c~~~~~~~~~~~~~~~~~~~~~~
 
c          Calculate current date and seconds, etc (for calday)

        call nutime (mcdate, mcsec, mbdate, nstep, dtime, 0)
        call caldyi (jyear, jmnth, jday, calday, calsec, istepy,
     *               mcdate, mcsec, mbdate, dtime, .false., .false.)

c          Get surface forcing by reading stored Genesis data files

        call async_mom (fuocnet, fvocnet, ftocnxy, fwocnxy, fsolocnxy,
#ifdef isotrac
     *                  foocnxy, fdocnxy,
#endif
     *                  maskbothh,  maskbothu, tocnxy,  socnxy,  odepth,
     *                  async_path, dtime,    segtim*86400., calday,
     *                  nstep,      nstop)

        istepy = nint (calsec/dtime)

c~~~~
#else
c~~~~
        call linems

#ifdef isotrac
c          Calculate isotopic fluxes at ocean sfc (foocnxy, fdocnxy).
c          Use local, not glob uniform d18O,dD values for surface sea 
c          water (ofstd,dfstd) to account for non-zero local P-E 
c          (unlike salinity, glob unif value doesn't ensure zero glob
c          mean flux if glob mean P-E = 0). Instead, do simple 
c          global correction here to ensure glob mean flux = 0. 
c          This makes no attempt to account for runofadj. For now, 
c          without riverperry (presc_rivers = .F.), it is consistent
c          with the global uniform salinity flux (fwocnadj) correction 
c          in pmeadj. But if we add explicit river input, would need
c          to account for d18O,dD of river input and not use the simple
c          correction below.
c
c          Calculate using fractions, then convert back to delta's 
c          (smow) for combined fluxes at end. 
c          rainaxy*,snowaxy*,avapaxy* are in commugxy.h, 
c          ri_*,ficexy in comicexy.h,c, ro_* in comocnxy.h.

        zfoocntot = 0.
        zfdocntot = 0.
        zareo     = 0.

        do j=1,nlat
          zwei = cos(alat(j))*(pie/nlat)*(2.*pie/nlon)
          do i=1,nlon

            if (maskogcmh(i,j).eq.1) then

c               open ocean (neglect runoff contribs of rivers, antarc):

              zp =  rainaxy(i,j) + snowaxy(i,j)
              ze = -avapaxy(i,j)    ! convert lsx's +ve dnwd to +ve upwd
              zpo=  rainaxy_o(i,j) + snowaxy_o(i,j)
              zeo= -avapaxy_o(i,j)
              zph=  rainaxy_h(i,j) + snowaxy_h(i,j)
              zeh= -avapaxy_h(i,j)

c             zfoo = 1000. * (zpo - zeo - (zp-ze)*ofstd) / smoj_o
c             zfdo = 1000. * (zph - zeh - (zp-ze)*dfstd) / smoj_h
              zfoo = 1000. * (zpo - zeo - (zp-ze)*ro_o(i,j)) / smoj_o
              zfdo = 1000. * (zph - zeh - (zp-ze)*ro_h(i,j)) / smoj_h

c               crude estimate under sea ice (assuming +ve fwocnxy has 
c               sea ice isotopic composition):

              zp = max (fwocnxy(i,j), 0.)
              ze = min (fwocnxy(i,j), 0.)
              zpo = zp * ri_o(i,j)
              zeo = ze * ro_o(i,j)
              zph = zp * ri_h(i,j)
              zeh = ze * ro_h(i,j)

c             zfoo_i = 1000. * (zpo - zeo - (zp-ze)*ofstd) / smoj_o
c             zfdo_i = 1000. * (zph - zeh - (zp-ze)*dfstd) / smoj_h
              zfoo_i = 1000. * (zpo - zeo - (zp-ze)*ro_o(i,j)) / smoj_o
              zfdo_i = 1000. * (zph - zeh - (zp-ze)*ro_h(i,j)) / smoj_h

              foocnxy(i,j) = (1.-ficexy(i,j))*zfoo + ficexy(i,j)*zfoo_i
              fdocnxy(i,j) = (1.-ficexy(i,j))*zfdo + ficexy(i,j)*zfdo_i

              zfoocntot = zfoocntot + foocnxy(i,j)*zwei
              zfdocntot = zfdocntot + fdocnxy(i,j)*zwei
              zareo     = zareo     +              zwei

            else

              foocnxy(i,j) = 0.
              fdocnxy(i,j) = 0.

            endif

          enddo
        enddo

c         Correction to make global occean mean fluxes f[o,d]ocnxy zero
c         (so ocean mean d18O, dDH will be constant through run...
c         see comments above)

        do j=1,nlat
          zwei = cos(alat(j))*(pi/nlat)*(2.*pie/nlon)
          do i=1,nlon
            if (maskogcmh(i,j).eq.1) then
              foocnxy(i,j) = foocnxy(i,j) - zfoocntot/zareo
              fdocnxy(i,j) = fdocnxy(i,j) - zfdocntot/zareo
            endif
          enddo
        enddo
#endif

c~~~~~
#endif
c~~~~~

c          Set history, restart write flags for ogcm. Done at all 
c          atmos timesteps, so will be true for next call to ogcm
c          if any preceding atmos timestep was flagged. Then reset 
c          to false after each ogcm call below.

       if (mod (jyear,nocyear).eq.0) then
         do i=1,nochisn
           if (istepy.eq. nint(nochis(i)*86400./dtime)) fwhisto = .true.
         enddo

         do i=1,nocsavn
           if (istepy.eq. nint(nocsav(i)*86400./dtime)) fsaveo = .true.
         enddo
       endif

       if (nstep.eq.nstop) fsaveo = .true.

c      needed by mom2 history file closing/opening logic:
       if (fsaveo) fwhisto = .true.

c          Accumulate ocean SBCs from GENESIS

        call ocnaccum (fuocnav, fvocnav, ftocnav, fwocnav, fsolocnav,
     *                 fuocnet, fvocnet, ftocnxy, fwocnxy, fsolocnxy,
#ifdef isotrac
     *                 foocnav, fdocnav,
     *                 foocnxy, fdocnxy,
#endif
     *                 segtim*86400., dtime, nstep, nrstrt, tocfilt)

c--------------------
c     C A L L   M O M
c--------------------

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if (mod(nstep,napero).eq.0) then 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c            From here on, use temporary copy of p-e fwocnadj, 
c            otherwise adjustments added here to fwocnav would 
c            be accumulated in ocnaccum if using exp.filter

          call scopy (nlon*nlat, fwocnav, 1, fwocnadj, 1)

c            Fill out values to entire globe (replaces mom2's extrap),
c            so we can do river/calv/glob adjusts to fwocnadj here    

          call fill_near (fuocnav,   maskbothu, 2, 1)
          call fill_near (fvocnav,   maskbothu, 2, 1)
          call fill_near (ftocnav,   maskbothh, 1, 1)
          call fill_near (fwocnadj,  maskbothh, 1, 1)
          call fill_near (fsolocnav, maskbothh, 1, 1)
#ifdef isotrac
          call fill_near (foocnav,   maskbothh, 1, 1)
          call fill_near (fdocnav,   maskbothh, 1, 1)
#endif

c            Set Perry et al (1996) river and/or Arctic river runoff 

          call zero (runofadj, nlon*nlat)
          if (presc_rivers) then 
            call riverperry (runofadj,maskogcmh,alon,alat,nstep,dtime)
          endif
c777      call riverarctic(runofadj, maskogcmh, alon, alat)

c            Set Antarctic iceberg discharge

          call zero (calvadj, nlon*nlat)
          if (presc_calvant) then
            call calvantarc (calvadj, maskogcmh, alat)
          endif

c            Global uniform adjust to fwocnadj to make global
c            ocean mean  <s*(e-p-runof-calv)> = 0

          call pmeadj (fwocnadj, runofadj, calvadj, maskogcmh, socnxy, 
     *                 alon, alat, nstep, dtime)

c            Print ascii output maps   
          
         mday  = mod (nstep*dtime/86400. - .001, 365.) + 1 
c        if ( mday.eq.1. .or. mday.eq.181 ) then
c        if ( mday.eq.181 .or. mday.eq.365 ) then
c        if ( mod (nstep, max(1,nint(86400./dtime))*30) .eq. 0 ) then
c        if ( mod (nstep, max(1,nint(86400./dtime))*1 ) .eq. 0 ) then
         if (.false.) then
c        if (.true.) then
           iu = 92
           call printmapocn ('uocnxy (m/s):', uocnxy, .02,   0., 
     *                       maskogcmu, nlon, nlat, nstep, iu)
           call printmapocn ('vocnxy (m/s):', vocnxy, .02,   0.,
     *                       maskogcmu, nlon, nlat, nstep, iu)
           call printmapocn ('tocnxy (C):', tocnxy,  2.,273.15,
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('socnxy (f-.035):',  socnxy, .0005, .035,
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('fuocnav (N/m2):', fuocnav, .01,   0., 
     *                       maskogcmu, nlon, nlat, nstep, iu)
           call printmapocn ('fvocnav (N/m2):', fvocnav, .01,   0.,
     *                       maskogcmu, nlon, nlat, nstep, iu)
           call printmapocn ('ftocnav (W/m2):', ftocnav,  20.,  0.,
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('fsolocnav (W/m2):',fsolocnav, 5., 0., 
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('fwocn (mm/d):',fwocnadj,1./86400.,0.,
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('runof (mm/d):',runofadj,1./86400.,0.,
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('calv (mm/d):' ,calvadj, 1./86400.,0.,
     *                       maskogcmh, nlon, nlat, nstep, iu)
#ifdef isotrac
           call printmapocn ('foocnav (pm*mm/d):',foocnav, 1./86400.,0.,
     *                       maskogcmh, nlon, nlat, nstep, iu)
           call printmapocn ('fdocnav (pm*mm/d):',fdocnav, 1./86400.,0.,
     *                       maskogcmh, nlon, nlat, nstep, iu)
#endif
 
           call printmapocn ('mom U (m/s):', sbcatm(1,1,12), .02*100,0.,
     *                       aland, imap2, jma, nstep, iu)
           call printmapocn ('mom V (m/s):', sbcatm(1,1,13), .02*100,0.,
     *                       aland, imap2, jma, nstep, iu)
           call printmapocn ('mom T (C):', sbcatm(1,1,8), 2., 0.,
     *                       aland, imap2, jma, nstep, iu)
           call printmapocn ('mom S (f-.035):', sbcatm(1,1,9),.0005, 0.,
     *                       aland, imap2, jma, nstep, iu)
#ifdef isotrac
           call printmapocn ('mom dO18 (pm):', sbcatm(1,1,10), .25, 0.,
     *                       aland, imap2, jma, nstep, iu)
           call printmapocn ('mom dD   (pm):', sbcatm(1,1,11), .25, 0.,
     *                       aland, imap2, jma, nstep, iu)
#endif
         endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            Genesis history fields:
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          hvap2 = 2.5104e6
          hfus2 = 0.3336e6
          rhow2 = 1.e3

c            eastward wind stress (N/m2)

          if (testhisl(113)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) = fuocnav(i,j) * maskogcmu(i,j)
              enddo
            enddo
            call puthilsx (workl, 113)
          endif   

c            northward wind stress (N/m2)

          if (testhisl(114)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) = fvocnav(i,j) * maskogcmu(i,j)
              enddo
            enddo
            call puthilsx (workl, 114)
          endif   
   
c            net heat flux into ocean (below sea ice) (W/m2)
  
          if (testhisl(115)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) = ftocnav(i,j) * maskogcmh(i,j)
              enddo
            enddo
            call puthilsx (workl, 115)
          endif   
   
c           net P-E (+runof+calv) into ocean (below sea ice) (m/s)
   
          if (testhisl(116)) then
             do j=1,nlat
               do i=1,nlon
                 workl(i,j) = ( ( fwocnadj(i,j)
     *                            + runofadj(i,j)+calvadj(i,j)
     *                          ) / rhow2
     *                        ) * maskogcmh(i,j)
               enddo
             enddo
             call puthilsx (workl, 116)
          endif  
   
c            river runoff into ocean (m/s)
   
          if (testhisl(117)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) = (runofadj(i,j)/rhow2) * maskogcmh(i,j)
              enddo
            enddo
            call puthilsx (workl, 117)
          endif  
   
c            calving flux into ocean (m/s)
   
          if (testhisl(118)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) =  (calvadj(i,j)/rhow2) * maskogcmh(i,j)
              enddo
            enddo
            call puthilsx (workl, 118)
          endif  
   
#ifdef isotrac
c            isotopic d18O flux into ocean (delta_O18 * m/s)
   
          if (testhisl(119)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) = ( ( foocnav(i,j) 
ccriv*                   + runofadj(i,j)*(ofriv-ofstd    )*1000./smoj_o
ccriv*                   + calvadj (i,j)*(ofriv-ofstd    )*1000./smoj_o
criv *                   + runofadj(i,j)*(ofriv-ro_o(i,j))*1000./smoj_o
criv *                   + calvadj (i,j)*(ofriv-ro_o(i,j))*1000./smoj_o
     *                         ) / rhow2
     *                       ) * maskogcmh(i,j)
              enddo
            enddo
            call puthilsx (workl, 119)
          endif  

c            isotopic dD flux into ocean (delta_D * m/s)
   
          if (testhisl(120)) then
            do j=1,nlat
              do i=1,nlon
                workl(i,j) = ( ( fdocnav(i,j) 
ccriv*                   + runofadj(i,j)*(dfriv-dfstd    )*1000./smoj_h
ccriv*                   + calvadj (i,j)*(dfriv-dfstd    )*1000./smoj_h
criv *                   + runofadj(i,j)*(dfriv-ro_h(i,j))*1000./smoj_h
criv *                   + calvadj (i,j)*(dfriv-ro_h(i,j))*1000./smoj_h
     *                         ) / rhow2
     *                       ) * maskogcmh(i,j)
              enddo
            enddo
            call puthilsx (workl, 120)
          endif  
#endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         end genesis history fields
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c         Convert units for mom2:  
c            dynes/cm2    -->  0.1 N/m2
c            cal/cm2/s    -->  41860. J/m2/s 
c            g/cm2/s      -->  10 kg/m2/s
c            *-1*salinity --> salt flux

          do j=1,jma
            do i=2,imap2-1 
              sbcatm(i,j,1) = fuocnav (i-1,j) * 10.
              sbcatm(i,j,2) = fvocnav (i-1,j) * 10.
              sbcatm(i,j,3) = ftocnav (i-1,j) / 41860.
              sbcatm(i,j,4) = - ((  fwocnadj(i-1,j)
     *                             +runofadj(i-1,j)
     *                             +calvadj(i-1,j)) /10.)
     *                        * .035 
c    *                        * socnxy(i,j)
c             conservation test: global uniform salinity flux:
c             sbcatm(i,j,4) = - ((30./86400.) / 10.) * .035
              sbcatm(i,j,5) = fsolocnav(i-1,j) / 41860.
#ifdef isotrac
              sbcatm(i,j,6) = foocnav(i-1,j) / 10.
ccriv*        + (runofadj(i-1,j)*(ofriv-ofstd      )*1000./smoj_o) / 10.
ccriv*        + (calvadj (i-1,j)*(ofriv-ofstd      )*1000./smoj_o) / 10.
criv *        + (runofadj(i-1,j)*(ofriv-ro_o(i-1,j))*1000./smoj_o) / 10.
criv *        + (calvadj (i-1,j)*(ofriv-ro_o(i-1,j))*1000./smoj_o) / 10.

              sbcatm(i,j,7) = fdocnav(i-1,j) / 10.
ccriv*        + (runofadj(i-1,j)*(dfriv-dfstd      )*1000./smoj_h) / 10.
ccriv*        + (calvadj (i-1,j)*(dfriv-dfstd      )*1000./smoj_h) / 10.
criv *        + (runofadj(i-1,j)*(dfriv-ro_h(i-1,j))*1000./smoj_h) / 10.
criv *        + (calvadj (i-1,j)*(dfriv-ro_h(i-1,j))*1000./smoj_h) / 10.

#else
              sbcatm(i,j,6) = 0.
              sbcatm(i,j,7) = 0.
#endif
            enddo
          enddo

c+++++++++++++++++++++++++++++++++
c         diagnostic global means:
c+++++++++++++++++++++++++++++++++
          if (.false.) then
c         if (.true.) then
            iu = 187
            ztot = 0.
            ztoto = 0.
            ztotoa = 0.
            ztotd = 0.
            ztotda = 0.
            zare = 0.
            do j=1,nlat
              zwei = cos(alat(j))
              do i=1,nlon
                if (maskogcmh(i,j).eq.1) then
                  ztot = ztot - zwei*((fwocnadj(i,j)
     *                                 +runofadj(i,j)
     *                                 +calvadj(i,j)) ) * 86400.
#ifdef isotrac
                  ztoto = ztoto   +  zwei*foocnav(i,j) *86400.
                  ztotoa = ztotoa +  zwei*abs(foocnav(i,j))  * 86400.
                  ztotd = ztotd   +  zwei*fdocnav(i,j) *86400.
                  ztotda = ztotda +  zwei*abs(fdocnav(i,j))  * 86400.
#endif
                  zare = zare + zwei
                endif
              enddo
            enddo
            ztot = ztot /zare
            ztoto = ztoto /zare
            ztotoa = ztotoa /zare
            ztotd = ztotd /zare
            ztotda = ztotda /zare

            ztos = 0.
            ztoo = 0.
            ztod = 0.
            zars = 0.
            do j=1,jma
              zwei = cos(alat(j))
              do i=2,imap2-1
                if (kmt(i,j+1).ge.1) then
                  ztos = ztos + zwei*sbcatm(i,j,4)
                  ztoo = ztoo + zwei*sbcatm(i,j,6)
                  ztod = ztod + zwei*sbcatm(i,j,7)
                  zars = zars + zwei
                endif
              enddo
            enddo
            ztos = ztos /zars
            ztoo = ztoo /zars
            ztod = ztod /zars

            write (iu,'(a,i10)') 'nstep=',nstep
            write (iu,'(2(a,e15.8)/3(a,e15.8)/3(a,e15.8)/2(a,e15.8))')
     *        '  ztot  =',ztot, '  ztos  =',ztos,
     *        '  ztoto =',ztoto,'  ztotoa=',ztotoa,'  ztoo=',ztoo,
     *        '  ztotd =',ztotd,'  ztotda=',ztotda,'  ztod=',ztod,
     *        '  zare  =',zare, '  zars  =',zars
          endif
c++++++++++++++++++++++++++++++++++++
c         end diagnostic global means
c++++++++++++++++++++++++++++++++++++

c            Mom2 time step

          write (nout,'(/a)') ' Calling MOM2 stp_driver'

          call iomngr_resume

          call stp_driver (dtime, fwhisto, fsaveo)

          call release_all

c            Set genesis's sst,etc, from mom arrays, convert units.
c            Set only on ocn pts in both lsx and mom.
c            Don't change lsx's values for pts that are ocean in lsx,
c            lamd in mom (as comments above).

          do j=1,nlat
            do i=1,nlon
              if (maskbothh(i,j).eq.1) then
                tocnxy(i,j) = sbcatm(i+1,j,8) + 273.16
                socnxy(i,j) = sbcatm(i+1,j,9) + .035 
#ifdef isotrac
                ro_o(i,j)   = (sbcatm(i+1,j,10)/1000. + 1.) *smoj_o
                ro_h(i,j)   = (sbcatm(i+1,j,11)/1000. + 1.) *smoj_h
#endif
               endif
            enddo
          enddo

c           For u and v (used in surfctl/dynice), fill globally to
c           include pts that are ocean in lsx, land in mom.

          do j=1,nlat
            do i=1,nlon
              if (maskbothu(i,j).eq.1) then
                uocnxy(i,j) = sbcatm(i+1,j,12) / 100.
                vocnxy(i,j) = sbcatm(i+1,j,13) / 100.
              endif
            enddo
          enddo
          call fill_near (uocnxy, maskbothu, 2, 1)
          call fill_near (vocnxy, maskbothu, 2, 1)
c         zero where not needed by dynice (just for niceity):
          do j=1,nlat
            jp = min (j+1,nlat)
            do i=1,nlon
              ip = mod (i,nlon) + 1
              if ( .not. 
     *              (lmaskxy(i,j)  .eq.3 .and. lmaskxy(ip,j) .eq.3 .and.
     *               lmaskxy(i,jp) .eq.3 .and. lmaskxy(ip,jp).eq.3)
     *           ) then
                uocnxy(i,j) = 0.
                vocnxy(i,j) = 0.
              endif
            enddo
          enddo

c           Zero accumulated ocean fluxes, if not using exp. time filter

          if (tocfilt.eq.0.) then
            call zero (fuocnav, nlon*nlat)
            call zero (fvocnav, nlon*nlat)
            call zero (ftocnav, nlon*nlat)
            call zero (fwocnav, nlon*nlat)
            call zero (fsolocnav, nlon*nlat)
#ifdef isotrac
            call zero (foocnav, nlon*nlat)
            call zero (fdocnav, nlon*nlat)
#endif
          endif

c            Reset ogcm history, restart write flags to false

          fwhisto = .false.
          fsaveo = .false.

c>>>>>>>>>>>>
        endif
c>>>>>>>>>>>>
 
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

c          Flush stdout at end of each day (done in linems if synch,
c          but still do here to get ogcm output)
c#if defined (asyncmode)
        n = max (nint(1.*86400./dtime), 1)
        if (mod(nstep,n).eq.0) ier = flush (nout)
c#endif
 
c++++++++++
      enddo
c++++++++++
c--------------------------------
c     M A I N   L O O P   E N D S
c--------------------------------

      write (nout,500)
  500 format(/' --------------------------'
     *       /' ******* END OF RUN *******'
     *       /' --------------------------'/)
      call endrun (0)

c        Fatal error messages

 1000 write (nout,1002)
 1002 format (
     *  /'*** Error opening Genesis namelist input file interm_coup'/)
      call endrun (-1)

#if defined (asyncmode)
 2000 write (nout,2002)
     *  async_surftype_file(1:lenchr(async_surftype_file))
 2002 format (
     *  /'*** Error opening Genesis data file (async_surftype_file):'
     *  /'    ',a/)
      call endrun (-1)
#endif

      end
c
c-----------------------------------------------------------------------
c
      subroutine ocnaccum (fuocnav, fvocnav, ftocnav, fwocnav,fsolocnav,
     *                     fuocnet, fvocnet, ftocnxy, fwocnxy,fsolocnxy,
#ifdef isotrac
     *                     foocnav, fdocnav,
     *                     foocnxy, fdocnxy,
#endif
     *                     dtocn,   dtime,   nstep,   nrstrt, tocfilt)

c       Updates ocean fluxes f*ocnav for current f*ocnxy, using    
c       time-exponential filter with time scale of tocfilt days.
c       For now the filtered fluxes are not saved on restart file,
c       and start from zero for all runs. If tocfilt is zero, 
c       don't use any filter, and just average between each
c       ogcm timestep (and f*ocnav are rezero-ed in main).

#include <compar.h>

      dimension 
     *  fuocnav(nlon,nlat), fvocnav(nlon,nlat),
     *  ftocnav(nlon,nlat), fwocnav(nlon,nlat), fsolocnav(nlon,nlat),
     *  fuocnet(nlon,nlat), fvocnet(nlon,nlat),
     *  ftocnxy(nlon,nlat), fwocnxy(nlon,nlat), fsolocnxy(nlon,nlat)
#ifdef isotrac
     * ,foocnav(nlon,nlat), fdocnav(nlon,nlat),
     *  foocnxy(nlon,nlat), fdocnxy(nlon,nlat)
#endif

      if (tocfilt.eq.0.) then
        ztf = 1.
        ztg = 1. / max (1,nint(dtocn/dtime))
      else
c       For first few timesteps, ramp up from no-filter
        zelaps = (nstep-nrstrt)*dtime
        if (zelaps.eq.0.) then
          ztf = 0.
        else
          ztf = exp ( -dtime / min(zelaps,tocfilt*86400.) )
        endif
c       Or: start filtering from initial (zero) values
c       ztf = exp ( -dtime / (tocfilt*86400.) )
        ztg = 1.-ztf
      endif

      do j=1,nlat
        do i=1,nlon
          fuocnav(i,j)= ztf*fuocnav(i,j) + ztg*fuocnet(i,j)
          fvocnav(i,j)= ztf*fvocnav(i,j) + ztg*fvocnet(i,j)
          ftocnav(i,j)= ztf*ftocnav(i,j) + ztg*ftocnxy(i,j)
          fwocnav(i,j)= ztf*fwocnav(i,j) + ztg*fwocnxy(i,j)
          fsolocnav(i,j)= ztf*fsolocnav(i,j) + ztg*fsolocnxy(i,j)
#ifdef isotrac
          foocnav(i,j)= ztf*foocnav(i,j) + ztg*foocnxy(i,j)
          fdocnav(i,j)= ztf*fdocnav(i,j) + ztg*fdocnxy(i,j)
#endif
        enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c 
      subroutine riverperry (runofxy, maskxy, alon, alat, nstep, dtime)

c     Perry et al.(1996) rivers:
c     Adjusts runofxy (ocean P-E) for discharge at river mouths,
c     including compensation at river mouths to make global P-E = 0.
 
c     runofxy=  net fresh water flx from river runoff (kg/m2/s)
c     nstep   = lsx timestep counter from basedate of run
 
#include <compar.h>
#include <COMLUN.h>
 
      dimension runofxy(nlon,nlat), maskxy(nlon,nlat),
     *          alon(nlon),  alat(nlat)

      parameter (pi=3.14159265358979)

c     Adds observed discharges from present-day rivers (Perry et al.),
c     to runofxy, as "p-e" distributed over the grid box area of each
c     "coastal" ogcm h-grid point closest to each river mouth.
 
      parameter (nriver=981)
      dimension alonriv(nriver), alatriv(nriver), flowriv(nriver),
     *          iriv(nriver),    jriv(nriver)
      character cnriv*30, ccriv*16
 
      dimension icoast(nlon*nlat), jcoast(nlon*nlat),
     *          xcoast(nlon*nlat), ycoast(nlon*nlat),
     *          zcoast(nlon*nlat)
 
      save iriv, jriv, flowriv
      logical first
      save first
      data first /.true./
 
c====================
      if (first) then
c====================
 
c          Read global river discharge data (from Perry et al,1996,JGR,
c          101,D16,21339-21349; anon ftp to kosmos.agu.org, cd apend,
c          cd 96jd00932.dir). Convert lats and lons to radians, and
c          change flow units to kg/s.
 
        write (nout,*) ' Reading riverperry.dat...'
        iu = 90
        open (iu, file='riverperry.dat',status='old')
        do n=1,nriver
          read (iu,100)
     *      m, cnriv, ccriv, alatriv(n), alonriv(n), nsoriv, flowriv(n),
     *      sdevriv, tflowriv
c         write(nout,100)
c    *      m, cnriv, ccriv, alatriv(n), alonriv(n), nsoriv, flowriv(n),
c    *      sdevriv, tflowriv
          alonriv(n) = alonriv(n)*pi/180.
          alatriv(n) = alatriv(n)*pi/180.
          flowriv(n) = flowriv(n) *1.e3    ! m3/s to kg/s
        enddo
        close (iu)
  100   format(i5,a30,a16,f9.1,f9.1,i9,f15.1,f10.1,f13.1)
 
c          Make list of indices (icoast,jcoast) of ogcm h-grid points
c          neighboring land points
 
        ncoast = 0
        do j=1,nlat
          jp = min (j+1,nlat)
          jm = max (j-1,1)
          do i=1,nlon
            ip = mod (i+nlon  ,nlon) + 1
            im = mod (i+nlon-2,nlon) + 1
            if ( maskxy(i,j).eq.1 .and.
     *           (      maskxy(im,j ).ne.1 .or. maskxy(ip,j ).ne.1
     *             .or. maskxy(i, jm).ne.1 .or. maskxy(i ,jp).ne.1
cdiag*             .or. maskxy(im,jp).ne.1 .or. maskxy(ip,jp).ne.1
cdiag*             .or. maskxy(im,jm).ne.1 .or. maskxy(ip,jm).ne.1
     *           )
     *         ) then
              ncoast = ncoast + 1
              icoast(ncoast) = i
              jcoast(ncoast) = j
              xcoast(ncoast) = cos(alat(j))*cos(alon(i))
              ycoast(ncoast) = cos(alat(j))*sin(alon(i))
              zcoast(ncoast) = sin(alat(j))
            endif
          enddo
        enddo
 
c          For each river, find gcm indices (iriv,jriv) of
c          nearest coastal point
 
        zflowtot = 0.
        zflowexc = 0.
        do ir=1,nriver
          xriv = cos(alatriv(ir))*cos(alonriv(ir))
          yriv = cos(alatriv(ir))*sin(alonriv(ir))
          zriv = sin(alatriv(ir))
          zmin = 1.e20
          do n=1,ncoast
c           straight-line distance, squared:
            zz = (xriv-xcoast(n))**2 + (yriv-ycoast(n))**2
     *         + (zriv-zcoast(n))**2
c           great-circle distance (more cpu time):
c           zz = acos(xriv*xcoast(n) + yriv*ycoast(n) + zriv*zcoast(n))
 
            if (zz.lt.zmin) then
              iriv(ir) = icoast(n)
              jriv(ir) = jcoast(n)
              zmin  = zz
            endif
          enddo
c         Diagnostic:
c         i = iriv(ir)
c         j = jriv(ir)
c         write (nout,200) ir, alonriv(ir)*180./pi, alon(i)*180./pi,
c    *                         alatriv(ir)*180./pi, alat(j)*180./pi,
c    *                         maskxy(iriv(ir),jriv(ir))
c 200     format('river=',i3,'  lonriv,loncoas=',2f8.2,
c    *                       '  latriv,latcoas=',2f8.2,'  maskxy=',i2)

          zdist = sqrt(zmin)*6371.22e3 

c            Exclude rivers if nearest gcm coastal point is > 1100 km
c            from actual mouth

          if ( zdist .gt. 1100.e3
c    *         .or. alatriv(ir)*180./pi.gt.65.) then               ! 777
     *       ) then
            iriv(ir) = 0
            jriv(ir) = 0
            write (nout,250) ir,
     *        nint(alonriv(ir)*180./pi), nint(alatriv(ir)*180./pi),
     *        nint(flowriv(ir)*1.e-3), nint(zdist*1.e-3)
  250       format('*** nb: Excluded riverperry #',i4,
     *             ' lon,lat=',2i4,'  flow(m3/s)=',i6,'  km=',i4)
            zflowexc = zflowexc + flowriv(ir)
          endif
          zflowtot = zflowtot + flowriv(ir)

        enddo

        write (nout, 260) zflowexc/zflowtot
  260   format ('*** nb: fraction of total flow excluded =',f7.4/)
 
        first = .false.
 
c==========
      endif
c==========
 
c        Add each river's flow to runofxy. zdar*cos(lat) = grid box
c        area in m**2, and flowriv is in kg/s.
 
      zdar = (6371.22e3)**2 * (2.*pi/nlon) * (pi/nlat)
      zflo = 0.
      do ir=1,nriver
        if (iriv(ir).ne.0) then
          i = iriv(ir)
          j = jriv(ir)
          runofxy(i,j) = runofxy(i,j) + flowriv(ir)/(zdar*cos(alat(j)))
          zflo = zflo + flowriv(ir)
        endif
      enddo
 
c        Diagnostics
 
      if (nstep.eq.1 .or. mod (nstep, nint(30.*86400/dtime)).eq.0) then
        zpe2 = 0.
        zareo = 0.
        do j=1,nlat
          do i=1,nlon
            if (maskxy(i,j).eq.1) then
              zpe2  = zpe2  + runofxy(i,j)*zdar*cos(alat(j))
              zareo = zareo +            zdar*cos(alat(j))
            endif
          enddo
        enddo
        zareo = max (zareo, 1.e-3)
        zz = zareo / 86400.                          ! kg/m2/s to mm/day
        write (nout,300) 'total perry river flow ',zflo*1.e-3,zflo/zz
        write (nout,300) 'total p-e after  adjust',zpe2*1.e-3,zpe2/zz
  300   format(a, ' (m3/s, mm/day ocn) =', f10.1, f8.4)
      endif
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine riverarctic (runofxy, maskxy, alon, alat)

c     Adjusts runofxy (ocean P-E) for Arctic river discharge
 
c     runofxy  = net fresh water flx from river runoff (kg/m2/s)
 
#include <compar.h>
 
      dimension runofxy(nlon,nlat), maskxy(nlon,nlat), 
     *          alon(nlon), alat(nlat)
      parameter (pi=3.14159265358979)

c         Add Arctic P-E due to river influx poleward of 70N,
c         about 35 cm/yr (Coachman+Aagard,1974,pgs.1 and 6;
c         Hibler+Bryan,1987,p.993), converting from cm/yr to kg/m2/s.
 
      do j=1,nlat
        zlat = alat(j)*180./pi
        do i=1,nlon
          zlon = alon(i)*180./pi
          if ( ( (zlat.ge.70. and. (zlon.ge.50..and.zlon.le.300.))
     *           .or.
     *           (zlat.ge.84.)
     *         )
     *         .and. maskxy(i,j).eq.1
     *       ) then
            runofxy(i,j) = runofxy(i,j) + 35.*(10./(365.*.864e5))
c           runofxy(i,j) = runofxy(i,j) + 20.*(10./(365.*.864e5)) !low 
          endif
        enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine calvantarc (calvxy, maskxy, alat)

c     Adjusts calvxy for discharge of ice from modern Antarctica
 
c     calvxy  = net fresh water flx from iceberg/shelf calving (kg/m2/s)
 
#include <compar.h>
#include <COMLUN.h>
 
      dimension calvxy(nlon,nlat), maskxy(nlon,nlat), alat(nlat)

      dimension maskh_ant(nlon,nlat)
      save maskh_ant, area_ant

      parameter (pi=3.14159265358979)
      logical first
      save first
      data first /.true./
 
c--------------------
      if (first) then
c--------------------
 
c          Set mask for ogcm h-grid points to receive prescribed
c          runoff (extra P-E) from Antarctic discharge
 
        area_ant = 0.
        do j=1,nlat
          zlat = alat(j)*180./pi
          zdarea = cos(alat(j))*(2.*pi/nlon)*(pi/nlat)*(6371.22**2)
          if (zlat.lt.-60.) then
            jp = min (j+1,nlat)
            jm = max (j-1,1)
            do i=1,nlon
 
c             Either set points neighboring coast:
c             ip = mod (i+nlon  ,nlon) + 1
c             im = mod (i+nlon-2,nlon) + 1
c             if( maskxy(i,j).eq.1 .and.
c    *            (maskxy(im,jp).ne.1 .or. maskxy(ip,jp).ne.1 .or.
c    *             maskxy(im,j ).ne.1 .or. maskxy(ip,j ).ne.1 .or.
c    *             maskxy(im,jm).ne.1 .or. maskxy(ip,jm).ne.1 .or.
c    *             maskxy(i, jm).ne.1 .or. maskxy(i ,jp).ne.1 )
c    *          ) then
 
c             Or just set all ocean points poleward of a latitude:
              if (zlat.lt.-64. .and. maskxy(i,j).eq.1) then
c             if (zlat.lt.-60. .and. maskxy(i,j).eq.1) then
 
                area_ant = area_ant + zdarea
                maskh_ant(i,j) = 1
              else
                maskh_ant(i,j) = 0
              endif
            enddo
          else
            do i=1,nlon
              maskh_ant(i,j) = 0
            enddo
          endif
        enddo
c       write (nout,*) 'maskh_ant:'
c       do j=nlat,1,-1
c         write (nout,'(180i1)') (maskh_ant(i,j),i=1,nlon)
c       enddo
        write (nout,*) 'area_ant (km**2)=', area_ant
 
        first = .false.
 
c----------
      endif
c----------
 
c        Modify calvxy for Antarctic discharge. Apply to ocean pts
c        given by maskh_ant=1, set above, with total area = area_ant 
c        (km2), so that total is equal to the mean annual Antarctic mass
c        balance, 10 to 20 cm/yr * all-Antarctic area of 14.39e6 km2
c        (Genesis v2 value). Convert from cm/yr to kg/m2/s.
 
      if (area_ant.gt.0.) then
c       zz = 10. * (10./(365.*.864e5)) * (14.39e6/area_ant)  !low  value
        zz = 15. * (10./(365.*.864e5)) * (14.39e6/area_ant)  !med  value
c       zz = 20. * (10./(365.*.864e5)) * (14.39e6/area_ant)  !hi   value
      else
        zz = 0.
      endif
 
      do j=1,nlat
        do i=1,nlon
          calvxy(i,j) = calvxy(i,j) + maskh_ant(i,j)*zz
        enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pmeadj (fwocn, runof, calv, maskxy, socn, 
     *                   alon, alat, nstep, dtime)

c        Adjust fwocn (P-E) to make global ocean mean of
c        salinity*(E-P-runof-calv) = 0
 
c     fwocn  = lsx net fresh water flx P-E (kg/m2/s) (modified)
c     runof  = river runoff (kg/m2/s) (supplied)
c     calv   = iceberg calving (kg/m2/s) (supplied)
c     maskxy = ogcm mask (lsx array)  (supplied)
c     socn   = surface ocean salinity (fraction) (supplied)
c     nstep  = lsx timestep counter from basedate of run (supplied)
c     dtime  = lsx timstep (s) (supplied)
  
#include <compar.h>
#include <COMLUN.h>
 
      dimension
     *  fwocn(nlon,nlat), runof(nlon,nlat), calv(nlon,nlat),
     *  maskxy(nlon,nlat), socn(nlon,nlat), alon(nlon), alat(nlat) 

      dimension faclat(nlat)
      parameter (pi=3.14159265358979)
      logical firstsocn, firstdiag
      save firstsocn, firstdiag
      data firstsocn, firstdiag /.true., .true./

c        Compute global oceanic means of salinity*(E-P) and salinity,
c        and correct P-E globally by <s(E-P)> / <s>.

c     do j=1,nlat
c       zlat = alat(j)*180./pi
c       do i=1,nlon
c         zlon = alon(i)*180./pi
c         if (zlon.gt.180.) zlon = zlon - 360.
c         if (maskxy(i,j).eq.1) then 
c           if (zlat.gt. 62. .and. zlat.lt.80. .and.
c    *          zlon.gt.-18. .and. zlon.lt.50.) then 
c             fwocn(i,j) = -0.3/86400.
c           endif
c         endif
c       enddo
c     enddo

c     fudge to not apply P-E correction in high N. lats (N.Atl, GIN)
      do j=1,nlat
        zlat = alat(j)*180./pi
c       faclat(j) = max (0., min (1., (60.-zlat)/(60-56.)))        ! 777
        faclat(j) = 1.                                             ! 777
      enddo
 
      zgasp = 0.
      zgas  = 0.
      do j=1,nlat
        zwei = cos(alat(j))*(pi/nlat)*(2.*pi/nlon)
        do i=1,nlon
          if (maskxy(i,j).eq.1) then
            if (socn(i,j).eq.0. .and. firstsocn) then
              zsocn = .035   ! socn is init. to 0 if not coupled restart
            else
c             zsocn = socn(i,j)
              zsocn = .035
            endif
            zgasp = zgasp - zwei*zsocn*(fwocn(i,j)+runof(i,j)+calv(i,j))
            zgas  = zgas  + zwei*zsocn*faclat(j)
          endif
        enddo
      enddo
      firstsocn = .false.
 
      do j=1,nlat
        do i=1,nlon
          if (maskxy(i,j).eq.1) then
            fwocn(i,j) = fwocn(i,j) + faclat(j)*zgasp/zgas
          endif
        enddo
      enddo

c        diagnostics
 
      if (firstdiag .or. mod (nstep,nint(30.*86400/dtime)).eq.0) then
        zpeo = 0.
        zareo = 0.
        do j=1,nlat
          zwei = cos(alat(j))*(pi/nlat)*(2.*pi/nlon)
          do i=1,nlon
            if (maskxy(i,j).eq.1) then
              zpeo  = zpeo  + zwei*(fwocn(i,j)+runof(i,j)+calv(i,j))
              zareo = zareo + zwei
            endif
          enddo
        enddo
        zareo = max (zareo, 1.e-10)
 
        write(nout,180) '<s(E-P-R-C)>/<s>', (zgasp/zgas)*86400.
        write(nout,180) 'final <P-E+R+C> ', (zpeo/zareo)*86400.
  180   format(a,' (mm/day ocn) =',f12.6)
        firstdiag = .false.
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fill_near (arr, mask, igrid, iocean)
 
c        Sets all non-ocean points in arr to its value at the 
c        nearest ocean point (mask=iocean).
c        igrid is mask index; if first call for each igrid,
c        the lookup indices inear,jnear are calced (saved).

#include <compar.h>
#include <COMLUN.h>

      dimension arr(nlon,nlat), mask(nlon,nlat)

      parameter (maxgrid=3)
      dimension inear(nlon,nlat,maxgrid), jnear(nlon,nlat,maxgrid), 
     *          iffirsgrid(maxgrid)
 
      save inear, jnear, iffirsgrid
      data iffirsgrid /maxgrid*0/

      parameter (pi=3.14159265358979) 

#if defined  (identical_grids)
c     do nothing if asserted that lsx and ogcm ocean masks are identical
      return
#endif

      if (igrid.gt.maxgrid) then
        write (nout,*) 
     *    '*** Error (fill_near): igrid=',igrid,'  maxgrid=',maxgrid 
        call endrun(-1)
      endif
 
c-------------------------------------
      if (iffirsgrid(igrid).eq.0) then
c-------------------------------------

c          If first call for this grid, set lookup indices 
c          (inear, jnear) at each land (non-ocean) pt pointing 
c          to closest ocean point

c          (iland,jland) is the land (non-ocean) point.
 
        do jland=1,nlat
c         if (igrid.eq.1) then                  ! h-grid
c           zlat = -90. + ((jland-.5)/nlat)*180.
c         else                                  ! u-grid (NE corner)
c           zlat = -90. + ((jland   )/nlat)*180. 
c         endif
          zlat = -90. + ((jland-.5)/nlat)*180.  ! use for all, no zrat=0
 
          zdx = 2.*pi*cos(zlat*pi/180.) / nlon
          zdy = pi / nlat
          zrat = zdx/zdy
 
          do iland=1,nlon
 
            if (mask(iland,jland).ne.iocean) then
 
c                Search for nearest ocean pt(mask=iocean), on perimeters
c                of widening squares centered on (iland,jland),
c                with sides 2*jdist + 1 in latitude and 2*idist + 1 in
c                longitude, with idist and jdist adjusted (using zrat)
c                so that sides are approx equal lengths on a sphere.
c                  (j1,j2) are j-limits of current square.
c                  (i1,i2) are i-limits of current square.
c                  Wrap around in longitude and "reflect" in latitude.
 
              do jdist=0,nlat
                idist = min (nint(jdist/zrat), nlon)
 
                j1 = jland - jdist
                j2 = jland + jdist
                i1 = iland - idist
                i2 = iland + idist
 
                do jx = j1,j2
 
                  if (jx.eq.j1 .or. jx.eq.j2) then
                    i3 = 1                 ! search along whole E-W side
                  else
                    i3 = max(i2-i1,1)      ! search only the 2 E-W edges
                  endif
 
                  do ix = i1,i2,i3
 
c                   wraparound in longitude:
                    i = mod (ix + 2*nlon - 1, nlon) + 1
 
c                   "reflect" in latitude:
                    if (jx.gt.nlat) then
                      j = max (1, min (nlat, 2*nlat - jx + 1))
                      i = mod (i + nlon/2 - 1, nlon) + 1
                    else if (jx.lt.1) then
                      j = max (1, min (nlat,           - jx + 1))
                      i = mod (i + nlon/2 - 1, nlon) + 1
                    else
                      j = jx
                    endif
 
                    if (mask(i,j).eq.iocean) then
                      inear(iland,jland,igrid) = i
                      jnear(iland,jland,igrid) = j
                      goto 40
                    endif
 
                  enddo   ! ix
                enddo     ! jx
 
              enddo       ! jdist
              write(nout,*) '*** Search fail (fill_near)',
     *                   ' iland=',iland,' jland=',jland,' igrid=',igrid
              call endrun (-1)
   40         continue
 
            else          ! mask eq iocean
 
              inear(iland,jland,igrid) = 0
              jnear(iland,jland,igrid) = 0
 
            endif

          enddo           ! iland
        enddo             ! jland

 
        iffirsgrid(igrid) = 1

c----------
      endif
c----------

c       Copy data from nearest ocean point
 
      do j=1,nlat
        do i=1,nlon
          if (inear(i,j,igrid).ne.0) then
            arr(i,j) = arr(inear(i,j,igrid),jnear(i,j,igrid))
          endif
        enddo
      enddo
 
      return
      end
c
c-----------------------------------------------------------------------
c
c~~~~~~~~~~~~~~~~~~~~~~
#if defined (asyncmode)
c~~~~~~~~~~~~~~~~~~~~~~

      subroutine async_mom (fuocn,   fvocn,   ftocn,  fwocn, fsolocn,
#ifdef isotrac
     *                      foocn,   fdocn,
#endif
     *                      maskbothh, maskbothu, tocn, socn, odepth,
     *                      cpath,   dtime,   dtcoup,  calday,
     *                      nstep,   nstop)

c     Reads mid-monthly surface forcing fields from Genesis data files 
c     (generated by debplotdriveogcm,  which in turn has read Genesis 
c     history files written during previous synchronous-run period),
c     interpolates from mid-monthly to current time, and calculates
c     current fluxes to ocean.

#include <compar.h>
#include <COMLUN.h>

c     returned:
c     *ocn   = fluxes at ocean surface, below sea ice (corresponding to
c              *av arguments and *xy arrays in genmom2)
c     supplied:
c     maskbothh= 1 if lsx and mom2 ocean, else 0 (lsx h grid)
c     maskbothu= 1 if lsx and mom2 ocean, else 0 (lsx u-grid)
c     tocn   = top-layer ocean temperature (K)
c     socn   = top-layer ocean salinity (fraction)
c     cpath  = path to Genesis data files
c     dtime  = agcm/lsx timestep (s)
c     dtcoup = coupling timstep(how long ogcm is run for const atmos)(s)
c     calday = calendar day within year (1. to 366.)
c     nstep  = agcm/lsx timestep number
c     nstop  = lsx final timestep number for this run

      dimension
     *  fuocn(nlon,nlat),       fvocn(nlon,nlat), 
     *  ftocn(nlon,nlat),       fwocn(nlon,nlat),    fsolocn(nlon,nlat),
#ifdef isotrac
     *  foocn(nlon,nlat),       fdocn(nlon,nlat),
#endif
     *  maskbothh(nlon,nlat),   maskbothu(nlon,nlat),
     *  tocn(nlon,nlat),        socn(nlon,nlat)
      character*80 cpath

      parameter (pi=3.14159265358979)

c       Surface forcing arrays (read from files by debplotdriveogcm)
c         tx* = eastward  wind stress (N/m2)
c         ty* = northward wind stress (N/m2)
c         rs* = net downward (absorbed) solar surface flux (W/m2)
c         ri* = net downward infrared surface flux (W/m2)
c         ta* = air temperature at lowest agcm layer (C)
c         qa* = specific humidity at lowest agcm layer (Kg/Kg)
c         ps* = surface pressure (mb)
c         va* = scalar wind speed at lowest agcm layer (m/s)
c         pr* = precip from agcm to surface (mm/day)
c         av* = evaporation from surface to agcm (mm/day)
c         hi* = sea-ice thickness (m)
c         fi* = sea-ice fraction (0-1)
c         to* = ocean (mixed-layer) temperature (C)
c         po* = pme at ocean surface (under sea ice) (mm/day)

      dimension
     *  tx(nlon,nlat),   txprev(nlon,nlat), txnext(nlon,nlat),
     *  ty(nlon,nlat),   typrev(nlon,nlat), tynext(nlon,nlat),
     *  rs(nlon,nlat),   rsprev(nlon,nlat), rsnext(nlon,nlat),
     *  ri(nlon,nlat),   riprev(nlon,nlat), rinext(nlon,nlat),
     *  ta(nlon,nlat),   taprev(nlon,nlat), tanext(nlon,nlat),
     *  qa(nlon,nlat),   qaprev(nlon,nlat), qanext(nlon,nlat),
     *  ps(nlon,nlat),   psprev(nlon,nlat), psnext(nlon,nlat),
     *  va(nlon,nlat),   vaprev(nlon,nlat), vanext(nlon,nlat),
     *  pr(nlon,nlat),   prprev(nlon,nlat), prnext(nlon,nlat),
     *  av(nlon,nlat),   avprev(nlon,nlat), avnext(nlon,nlat),
     *  hi(nlon,nlat),   hiprev(nlon,nlat), hinext(nlon,nlat),
     *  fi(nlon,nlat),   fiprev(nlon,nlat), finext(nlon,nlat),
     *  to(nlon,nlat),   toprev(nlon,nlat), tonext(nlon,nlat),
     *  po(nlon,nlat),   poprev(nlon,nlat), ponext(nlon,nlat)
#ifdef isotrac
     * ,fo(nlon,nlat),   foprev(nlon,nlat), fonext(nlon,nlat),
     *  fd(nlon,nlat),   fdprev(nlon,nlat), fdnext(nlon,nlat)
#endif

      save txprev, txnext,
     *     typrev, tynext,
     *     rsprev, rsnext,
     *     riprev, rinext,
     *     taprev, tanext,
     *     qaprev, qanext,
     *     psprev, psnext,
     *     vaprev, vanext,
     *     prprev, prnext,
     *     avprev, avnext,
     *     hiprev, hinext,
     *     fiprev, finext,
     *     toprev, tonext,
     *     poprev, ponext
#ifdef isotrac
     *    ,foprev, fonext,
     *     fdprev, fdnext
#endif

c        Surface forcing unit numbers and scaling factors

      data
     *  iutx,iuty,iurs,iuri,iuta,iuqa,iups,
     *  iuva,iupr,iuav,iuhi,iufi,iuto,iupo,
     *  iufo,iufd
     *  / 71,  72,  73,  74,  75,  76,  77,
     *    78,  79,  80,  81,  82,  83,  84,
     *    85,  86/
      save
     *  iutx,iuty,iurs,iuri,iuta,iuqa,iups,
     *  iuva,iupr,iuav,iuhi,iufi,iuto,iupo,
     *  iufo,iufd
     *  scaletx, scalety, scalers, scaleri, scaleta, scaleqa, scaleps,
     *  scaleva, scalepr, scaleav, scalehi, scalefi, scaleto, scalepo,
     *  scalefo, scalefd

      parameter (nmonmid=12)
      dimension monmid(nmonmid), monlen(nmonmid)
      data monmid /  16,  45,  75, 105, 136, 166,
     *              197, 228, 258, 289, 319, 350 /
      data monlen /  31,  28,  31,  30,  31,  30,
     *               31,  31,  30,  31,  30,  31 /
      save monmid, monlen

      logical first, firstcalc
      data first, firstcalc /.true., .true./
      save first, firstcalc

      save stef,vonk,grav,tmelt,hvap,hfus,ch2o,cair,cvap,rair,rvap,
     *     cappa,rhow,rhoi,basefac,z0ocnh,emiocn,siga

      character*240 cfile

#include "comsat.h"

c--------------------
      if (first) then
c--------------------

c       If using table svp as in agcm (and not Lowe polys in comsat.h),
c       initialize table

#ifdef tablesvp
        call svpinit
#endif

        stef = 5.66961e-8
        vonk = 0.4
        grav = 9.80616
        tmelt = 273.16
        hvap = 2.5104e6
        hfus = 0.3336e6
        ch2o = 4.218e3
        cair = 1.00464e3
        cvap = 1.81e3
        rair = 287.04
        rvap = 461.
        cappa = rair / cair
        rhow = 1.e3
        rhoi = 0.9e3
        basefac = 70.   ! ocn<->sea ice base transfer coeff (W/m2/K,lsx)
        z0ocnh  = .0001 ! open ocean roughness length (m, lsx)
        emiocn = 1.     ! ir emissivity of open ocean (lsx)
        siga = .993     ! sigma (p/ps) value of lowest agcm layer
c       odepth = 25.    ! thickness of ogcm top layer (m) (now passed)

        m = lenchr(cpath)
        if (cpath.ne.' ' .and. cpath(m:m).ne.'/') then
           cpath(m+1:m+1) = '/'
        endif

c          Open surface-forcing files

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_TX'
        open(iutx, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_TY'
        open(iuty, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_RS'
        open(iurs, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_RI'
        open(iuri, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_TA'
        open(iuta, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_QA'
        open(iuqa, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_PS'
        open(iups, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_VA'
        open(iuva, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_PR'
        open(iupr, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_AV'
        open(iuav, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_HI'
        open(iuhi, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_FI'
        open(iufi, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_TO'
        open(iuto, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_PO'
        open(iupo, file=cfile, status='old', err=2000)

#ifdef isotrac
        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_FO'
        open(iufo, file=cfile, status='old', err=2000)

        cfile = cpath(1:lenchr(cpath)) // 'data_ogcmdrive_FD'
        open(iufd, file=cfile, status='old', err=2000)
#endif

        first = .false.

c----------
      endif
c----------

c        Main timestep loop

c     frequency for updating sfc-forcing time-interp
      nspf   = max (1, nint((86400.*1)/dtime)) 
c     nspf   = 1

c        Set surface forcing

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (nstep.eq.1 .or. mod(nstep,nspf).eq.0 .or. firstcalc) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c        Update time-interpolation, read next month's fields if time

      call dopreso (calday, iutx, tx, txprev, txnext, 'TX      ',
     *              scaletx, nlon, nlat)
      call dopreso (calday, iuty, ty, typrev, tynext, 'TY      ',
     *              scalety, nlon, nlat)
      call dopreso (calday, iurs, rs, rsprev, rsnext, 'RS      ',
     *              scalers, nlon, nlat)
      call dopreso (calday, iuri, ri, riprev, rinext, 'RI      ',
     *              scaleri, nlon, nlat)
      call dopreso (calday, iuta, ta, taprev, tanext, 'TA      ',
     *              scaleta, nlon, nlat)
      call dopreso (calday, iuqa, qa, qaprev, qanext, 'QA      ',
     *              scaleqa, nlon, nlat)
      call dopreso (calday, iups, ps, psprev, psnext, 'PS      ',
     *              scaleps, nlon, nlat)
      call dopreso (calday, iuva, va, vaprev, vanext, 'VA      ',
     *              scaleva, nlon, nlat)
      call dopreso (calday, iupr, pr, prprev, prnext, 'PR      ', 
     *              scalepr, nlon, nlat)
      call dopreso (calday, iuav, av, avprev, avnext, 'AV      ', 
     *              scaleav, nlon, nlat)
      call dopreso (calday, iuhi, hi, hiprev, hinext, 'HI      ', 
     *              scalehi, nlon, nlat)
      call dopreso (calday, iufi, fi, fiprev, finext, 'FI      ',
     *              scalefi, nlon, nlat)
      call dopreso (calday, iuto, to, toprev, tonext, 'TO      ',
     *              scaleto, nlon, nlat)
      call dopreso (calday, iupo, po, poprev, ponext, 'PO      ',
     *              scalepo, nlon, nlat)
#ifdef isotrac
      call dopreso (calday, iufo, fo, foprev, fonext, 'FO      ',
     *              scalefo, nlon, nlat)
      call dopreso (calday, iufd, fd, fdprev, fdnext, 'FD      ',
     *              scalefd, nlon, nlat)
#endif

c        Change various units from those on data files to MKS

      do j=1,nlat
        do i=1,nlon
c         deg C to deg K:
          ta(i,j) = ta(i,j) + tmelt
          to(i,j) = to(i,j) + tmelt

c         mb to N/m2:
          ps(i,j) = 100.*ps(i,j)

c         mm/day to kg/m2/s:
          pr(i,j) = pr(i,j)/86400.
          av(i,j) = av(i,j)/86400.
          po(i,j) = po(i,j)/86400.
#ifdef isotrac
          fo(i,j) = fo(i,j)/86400.
          fd(i,j) = fd(i,j)/86400.
#endif

c         increase wind speeds to compensate for using monthly means:
c         va(i,j) = 1.3*va(i,j)                                    ! 777
          va(i,j) = max (va(i,j), 1.)
        enddo
      enddo

c        Calculate net ocean surface heat flux (like ocnflx in lsx)

      tfac   = 1./(siga**cappa)
      cdhfac = 1.0
      cdefac = 1.35
      cdhfac = 1.0
      zlogh  = log (z0ocnh)

      do j=1,nlat
        do i=1,nlon

c--------------------------------------
          if (maskbothh(i,j).eq.1) then
c--------------------------------------

c           Atmos conditions at za, etc (like turini):
            psurf = ps(i,j)
            pa = psurf * siga
            rhoa = pa / (rair * ta(i,j) * (1.+(rvap/rair-1.)*qa(i,j)))
            cp = cair * (1. + (cvap/cair-1.)*qa(i,j))
            za = (psurf - pa) / (rhoa * grav)
            aloga = log (za)
            aloggh = zlogh

c           Ocean surface specific humidity:
            es = esat (tocn(i,j))
            qocn = qsat (es, psurf)

c           Calculate stratification factor for heat/vapor for ocean
c           surface to za, as in Louis's system (vi), 1982, ecmwf
c           workshop on planetary boundary layer parameterizations,
c           november 1981,59-79 (like fstrat):

            zht = ta(i,j)  *tfac*(1.+.622*qa(i,j))
            zhb = tocn(i,j)*     (1.+.622*qocn)

            rich = grav * max (za, 0.) * (zht-zhb)
     *           / ( 0.5*(zht+zhb) * (max(va(i,j),.01))**2 )
            rich = max (-2., min (1., rich))
            if (rich.le.0.) then
              xh = max (aloga-aloggh, .5)
              yh = (vonk/xh)**2 * exp (0.5*xh) * sqrt(-rich)
              strahu = 1. - 3*5*rich / (1. + 75*yh)
            else
              z = sqrt(1.+5*rich)
              strahu = 1. / (1.+ 3*5*rich*z)
            endif

c           Drag coefficients:
            cdh = cdhfac * va(i,j)*(vonk/(aloga-aloggh))**2 * strahu
            cde = cdefac * va(i,j)*(vonk/(aloga-aloggh))**2 * strahu
c           Crudely account for enhanced evap of ocean spray at high u
c           (currently commented out of lsx)
c           (Bortkovskii, in Smith,BLM,1989,47,277-293,Fig.1,Eq.13):
c           cde = cde + 1.5*1.e-3 * sqrt (max(va(i,j)-15.,0.)/5. )

c           Empirical adjustments to agree with synchronous ann means:
c           cdh = 0.65*cdh                                         ! 777
c           cde = 0.65*cde                                         ! 777
c           cdh = 0.75*cdh                                         ! 777
c           cde = 0.75*cde                                         ! 777

c           Upward sensible heat and water vapor fluxes 

c           zfseno = rhoa * cdh * cp * (tocn(i,j)-ta(i,j)*tfac)
c           zfvapo = rhoa * cde      * (qocn     -qa(i,j))

c           Or partially time-implicit in ta,qa,ux,uy (as lsx ocean.F):
            zthica = psurf*(1.-.985) / grav
            zmt = zthica/1800.
            zfach = zmt * (1.-exp(-rhoa*cdh*tfac/zmt))
            zface = zmt * (1.-exp(-rhoa*cde/zmt))
            zfseno = zfach * cp * (tocn(i,j)/tfac - ta(i,j))
            zfvapo = zface *      (qocn           - qa(i,j))

c           Net downward heat flux
            ztflx = rs(i,j)
     *            + ri(i,j)
     *            - zfseno
c    *            - hvap*zfvapo    ! evap calculated
     *            - hvap*av(i,j)   ! evap from gcm

            if (to(i,j).gt.tmelt-5.) then        ! not if to=null (land)
c             correct for change in ocean-surface emitted ir: 
              ztflx = ztflx - emiocn*stef*(tocn(i,j)**4-to(i,j)**4)
            endif

c           Under sea ice, linear exchange (as in mom4...none in mom2 ?)
            tfreez = tmelt - .054 * 1000.*socn(i,j)
            ztflxice = - basefac * max(tocn(i,j)-tfreez,0.)
c           ztflxice = 0.

            ftocn(i,j) = ztflx*(1.-fi(i,j)) + ztflxice*fi(i,j)

c           Limit to freezing top ogcm layer in next coupled timestep:
c           ztf =  max (tocn(i,j)-tfreez, 0.) 
            ztf =       tocn(i,j)-tfreez
            ftocn(i,j) = max (ftocn(i,j), -rhow*ch2o*odepth*ztf/dtcoup)

c           fwocn(i,j) = (pr(i,j)-zfvapo)*(1.-fi(i,j))  ! evap calced
c           fwocn(i,j) = (pr(i,j)-av(i,j))*(1.-fi(i,j)) ! gcm above s.i.
            fwocn(i,j) =  po(i,j)                       ! gcm below s.i.

            fsolocn(i,j) = rs(i,j)*(1.-fi(i,j))

#ifdef isotrac
            foocn(i,j) = fo(i,j)                        ! below s.i. 
            fdocn(i,j) = fd(i,j)                        ! below s.i.
#endif

c           Empirical adjustments to agree with synchronous ann means:
c           if (fi(i,j).lt.0.1) ftocn(i,j) = ftocn(i,j) + 2.0      ! 777
c           fwocn(i,j) = fwocn(i,j)*1.03                           ! 777
c           fsolocn(i,j) = 1.10 * fsolocn(i,j)                     ! 777
c           ftocn(i,j) = ftocn(i,j) + 2.                           ! 777

c-------------
          else
c-------------
            ftocn(i,j) = 0.
            fwocn(i,j) = 0.
            fsolocn(i,j) = 0.
#ifdef isotrac
            foocn(i,j) = 0.
            fdocn(i,j) = 0.
#endif
c--------------
          endif
c--------------

        enddo
      enddo

c        Set wind stresses (under ice), shifted to NE corner of 
c        grid boxes as expected by mom2. 

      do j=1,nlat
        jp = min (j+1, nlat)
        do i=1,nlon
          ip = mod (i,nlon) + 1
          if (maskbothu(i,j).eq.1) then
            fuocn(i,j) = 0.25 * (tx(i,j)+tx(ip,j)+tx(i,jp)+tx(ip,jp))
            fvocn(i,j) = 0.25 * (ty(i,j)+ty(ip,j)+ty(i,jp)+ty(ip,jp))
          else
            fuocn(i,j) = 0.
            fvocn(i,j) = 0.
          endif
        enddo
      enddo

      firstcalc = .false.

c>>>>>>>>>>
      endif
c>>>>>>>>>>

c        Diagnostic printout

c     iu = 6
c     zs = 86400.
c     call gcheck ('TX',tx,maskbothu,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('TY',ty,maskbothu,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('RS',rs,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('RI',ri,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('TA',ta,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('QA',qa,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('PS',ps,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('VA',va,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('PR',pr,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('AV',av,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('HI',hi,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('FI',fi,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('TO',to,maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('ftocn',ftocn, maskbothh,alat,nlon,nlat,1.,nstep,iu)
c     call gcheck ('fwocn',fwocn, maskbothh,alat,nlon,nlat,zs,nstep,iu)
c     call gcheck ('TOCN', tocn,  maskbothh,alat,nlon,nlat,1.,nstep,iu)

c       Clean up if last timestep of run

      if (nstep.eq.nstop) then
        close (iutx)
        close (iuty)
        close (iurs)
        close (iuri)
        close (iuta)
        close (iuqa)
        close (iups)
        close (iuva)
        close (iupr)
        close (iuav)
        close (iuhi)
        close (iufi)
        close (iuto)
      endif

      return

 2000 write (nout,2002) cfile(1:lenchr(cfile))
 2002 format
     *   (/'*** Error opening an ogcm surface-forcing driver file:'
     *    /'    ',a
     *    /)
      call endrun (-1)

      end
c
c-----------------------------------------------------------------------
c
      subroutine dopreso (calday, iu, field, fieldprev, fieldnext, ckey,
     *                    scale, nlon, nlat)

c    When necessary, reads mid-month driving fields from Genesis-format
c    data input file into arrays fieldprev and fieldnext, 
c    and lineary interpolates field[prev,next] in time between the
c    midmonth points to current array (field).

c calday   = day number within current year (1. to 366.) (supplied)
c iu       = unit number of already-open data input file  (supplied)
c field    = 2D current-time array (returned)
c fieldprev= 2D field at previous mid month (modified)
c fieldnext= 2D field at next mid month (modified)
c ckey     = alpha key for this data file (for readdata) (supplied)
c scale    = scaling factor for field (modified)
c nlon    = longitude index of field* (supplied)
c nlat    = latitude  index of field* (supplied)

#include <COMLUN.h>

      dimension field(nlon,nlat),
     *          fieldprev(nlon,nlat), fieldnext(nlon,nlat)

      character*(*) ckey

      parameter (nmon=12)
      dimension iday(0:nmon+1)
      logical first(200)
      dimension im2sav(200)
      save first, im2sav, iday
      data first /200*.true./
      data iday / -15.,  16.,  45.,  75., 105., 136., 166.,
     *                  197., 228., 258., 289., 319., 350.,381. /

c        Find interpolation indices for iday (nb: calday runs from
c        1. to 366. from 00:00 GMT Jan 1st to 24:00 GMT Dec 31st)

      do im=1,nmon+1
        if (calday.le.iday(im)) then
          im1 = im-1
          im2 = im
          wm1 = (iday(im2)-calday) / (iday(im2)-iday(im1))
          if (im1.eq.0) im1 = nmon
          if (im2.eq.nmon+1) im2 = 1
          goto 102
        endif
      enddo
      write (nout,*) '*** Error (dopreso): Bad calday=',calday
      call endrun (-1)
  102 continue

c        If first call, read through file to "prev" position, then read
c        in "next"

      if (first(iu)) then
        read (iu,200) scale
  200   format (7x,f12.5)
        do imon=1,im1
          call readdata (iu, fieldprev, nlon, nlat, ckey, 2, imon)
        enddo
        if (im2.eq.1) then
          rewind (iu)
          read (iu,200) scale
        endif
        call readdata (iu, fieldnext, nlon, nlat, ckey, 2, im2)
        im2sav(iu) = im2
        first(iu)=.false.
c       Divide by scaling factor:
        do j=1,nlat
          do i=1,nlon
            fieldprev(i,j) = fieldprev(i,j) / scale
            fieldnext(i,j) = fieldnext(i,j) / scale
          enddo
        enddo
      endif

c        If next-month time, copy "next" to "prev", read new "next"

      if (im2.ne.im2sav(iu)) then

        call scopy (nlon*nlat, fieldnext, 1, fieldprev, 1)

        if (im2.eq.1) then
          rewind (iu)
          read (iu,200) scale
        endif
        call readdata (iu, fieldnext, nlon, nlat, ckey, 2, im2)

        im2sav(iu) = im2

c          Divide by scaling factor

        do j=1,nlat
          do i=1,nlon
            fieldnext(i,j) = fieldnext(i,j) / scale
          enddo
        enddo

      endif

c        Interpolate to current time

      do j=1,nlat
        do i=1,nlon
          field(i,j) = (wm1*fieldprev(i,j) + (1.-wm1)*fieldnext(i,j))
        enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine gcheck (ckeyin, field, mask, alat, nlon, nlat, factor,
     *                   nstep, nspd, nspy)

#include <COMLUN.h>

      character*(*) ckeyin, ckey*8
      dimension field(nlon,nlat), mask(nlon,nlat), alat(nlat)
      parameter (pi=3.14159265358979)

      save npe,nhflx,petot,hflxtot
      data npe,nhflx,petot,hflxtot /0,0, 0., 0./

c        Prints diagnostics for global array field

      ckey = ckeyin

      gav   = 0.
      garea = 0.

      do j=1,nlat
        zwei = cos(alat(j))*(pi/nlat)*(2.*pi/nlon)
        do i=1,nlon
          if (mask(i,j).eq.1) then
c           gav   = gav   + zwei*abs(field(i,j))*factor
            gav   = gav   + zwei*    field(i,j) *factor
            garea = garea + zwei
          endif
        enddo
      enddo
      if (garea.gt.0) gav = gav/garea

c        Report averages up to this point in the run for PE and HFLX

      if (ckey.eq.'PE') then
        npe = npe + 1
        petot = petot + gav
        gav = petot/npe   
        nav = npe
      else if (ckey.eq.'HFLX') then
        nhflx = nhflx + 1
        hflxtot = hflxtot + gav
        gav = hflxtot/nhflx
        nav = nhflx
      else
        nav = 1
      endif

      if ( nstep.eq.1 .or. mod(nstep,nspd*30).eq.0
     *                .or. mod(nstep,nspy)   .eq.0 ) then
        write(nout,200) ckey, garea/(4.*pi), nav, gav
  200   format
     *    ('gcheck : ',a,'  fracoc=',f7.4,'  nav=',i5,'  oc.av.=',f13.5)
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      logical function testhisl (iv)
      testhisl = .false.
      return
      end

      subroutine puthilsx (dum, idum)
      return
      end
c
c-----------------------------------------------------------------------
c
#ifdef tablesvp

      subroutine svpinit

c  Table of saturation vapor pressure es from 173.16 to 376.16 K
c  in one degree increments. As in convec.F. Alternative to Lowe
c  polynomials in comsat.h. (namelist tablesvp).

      parameter (ntable=203, tmin=173.16)
      common /csatc/ table(0:ntable), dtable(0:ntable)

      data table     /.001403,.001719,.002101,.002561,.003117,.003784,
     1.004584,.005542,.006685,.008049,.009672, .01160, .01388, .01658,
     2 .01977, .02353, .02796, .03316, .03925, .04638, .05472, .06444,
     3 .07577, .08894,  .1042,  .1220,  .1425,  .1662,  .1936,  .2252,
     4  .2615,  .3032,  .3511,  .4060,  .4688,  .5406,  .6225,  .7159,
     5  .8223,  .9432,  1.080,  1.236,  1.413,  1.612,  1.838,  2.092,
     6  2.380,  2.703,  3.067,  3.476,  3.935,  4.449,  5.026,  5.671,
     7  6.393,  7.198,  8.097,  9.098,  10.21,  11.45,  12.83,  14.36,
     8  16.06,  17.94,  20.02,  22.33,  24.88,  27.69,  30.79,  34.21,
     9  37.98,  42.13,  46.69,  51.70,  57.20,  63.23,  69.85,  77.09,
     a  85.02,  93.70, 103.20, 114.66, 127.20, 140.81, 155.67, 171.69,
     1 189.03, 207.76, 227.96, 249.67, 272.98, 298.00, 324.78, 353.41,
     2 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78, 656.62,
     3 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2, 1147.4,
     4 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3, 1936.7,
     5 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1, 3167.1,
     6 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1, 5030.7,
     7 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7, 7780.2,
     8 8201.5, 8642.3, 9103.4, 9585.5, 10089., 10616., 11166., 11740.,
     9 12340., 12965., 13617., 14298., 15007., 15746., 16516., 17318.,
     a 18153., 19022., 19926., 20867., 21845., 22861., 23918., 25016.,
     1 26156., 27340., 28570., 29845., 31169., 32542., 33965., 35441.,
     2 36971., 38556., 40198., 41898., 43659., 45481., 47367., 49317.,
     3 51335., 53422., 55580., 57809., 60113., 62494., 64953., 67492.,
     4 70113., 72819., 75611., 78492., 81463., 84528., 87688., 90945.,
     5 94302., 97761.,101325.,104994.,108774.,108774./

      do i=1,ntable-1
        dtable(i) = (table(i+1)-table(i-1)) / 2.
      enddo
      dtable(0) = dtable(1)
      dtable(ntable) = dtable(ntable-1)

      return
      end
 
c-----------------------------------------------------------------------

      function esat(t)

c     returns svp (N/m2) as a function of temperature t (K)

      parameter (ntable=203, tmin=173.16)
      common /csatc/ table(0:ntable), dtable(0:ntable)
      it1 = int(t-tmin)
      zw  = (t-tmin) - it1
      it1 = max (0, min (ntable-1,it1))
      esat = (1.-zw)* table(it1) + zw* table(it1+1)
      return
      end
 
c-----------------------------------------------------------------------

      function qsat(es,p)

c     returns saturation specific humidity (kg/kg) as a function 
c     of svp es (N/m2) and ambient pressure (N/m2). Impose an upper
c     limit of 1 to avoid spurious values for very high es and/or
c     small p.

      qsat = .622 * es / max (p-(1.-.622)*es, .622*es)
      return
      end
 
c-----------------------------------------------------------------------

      function dsat(t,p)
c
c     returns d(sat. spec. humidity)/dt as a function of t (K)
c     and ambient pressure p (N/m2). Neglects q in denominator of the
c     q(vp) reln.

      parameter (ntable=203, tmin=173.16)
      common /csatc/ table(0:ntable), dtable(0:ntable)
      it1 = int(t-tmin)
      zw  = (t-tmin) - it1
      it1 = max (0, min (ntable-1,it1))
      es = (1.-zw)* table(it1) + zw* table(it1+1)
      ds = (1.-zw)*dtable(it1) + zw*dtable(it1+1)
      dsat = .622 * ds / max (p-(1.-.622)*es, .622*es)
c     write (166,*) 't=',t,'  p=',p,'  it1=',it1,' dtable=',dtable(it1)
      return
      end
 
#endif

c~~~~~
#endif
c~~~~~

c-----------------------------------------------------------------------

      subroutine printmapocn (ctitle, arr, scale, shift, 
     *                        masko, nx, ny, nstep, iuin)

c         Prints crude ascii map of arr(nx,ny) to unit |iuin|, 
c         masked by masko (0=land,else ocean) if iuin >= 0.
c         If iuin < 0, ignore masko, write all pts.
 
       character*(*) ctitle
       dimension arr(nx,ny), masko(nx,ny)

       integer i,j
       real zval
       character*1 cout(182,92)  ! so can be called from mom2 (imt,jmt)

       iu = iabs(iuin)

       do j = 1,ny
         do i = 1,nx
           zval = (arr(i,j)-shift) / scale
           ival = nint(abs(zval))
           if (masko(i,j).eq.0 .and. iuin.ge.0) then !ignore mask if -iu
             cout(i,j) = ' '   ! not ocean
           else
             if (arr(i,j).eq.0.1e30) then
               cout(i,j) = '?' ! should never be (mpp non-domain)
             else if (arr(i,j).eq.0.1e39) then
               cout(i,j) = '!' ! null vel in mom (?)
             else if (zval.ge.0. .and. ival.le.25) then
               cout(i,j) = char (ichar('A') + ival)
             else if (zval.ge.0. .and. ival.gt.25) then
               cout(i,j) = '*'
             else if (zval.lt.0. .and. ival.le.25) then
               cout(i,j) = char (ichar('a') + ival)
             else if (zval.lt.0. .and. ival.gt.25) then
               cout(i,j) = '#'
             else 
               cout(i,j) = '@'
             endif
           endif
         enddo
       enddo
       write (iu,'(/2a,f20.10,a,i8)')
     *    ctitle, '   scale=',scale, '   nstep=',nstep
       do j=ny,1,-1
         write (iu,'(182a1)') (cout(i,j),i=1,nx)
       enddo

c      write (iu,*)
c      write (iu,*) 'field=',ctitle, '  nstep=',nstep
c      do j=ny,1,-1
c        write (iu,'(182e15.5)') (arr(i,j),i=1,nx)
c      enddo

       return
       end 

c*****
#endif
c*****
