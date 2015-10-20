c -- $Id: histlsx.F,v 1.2 2001/07/03 21:29:28 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/histlsx.F,v $
c -- $Name:  $

c--------
c histlsx
c--------
c
c-----------------------------------------------------------------------
c
      subroutine datgelsx (dtgcm, nsrest, nrstrt, nstop,
     *                     nhisi, ndhis,  histlast)
c
c        Reads and processes lsx namelist inputlsx.
c        Similar to ccm subr datget in init.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
#include <comice.h>
#ifdef gcm
#include <COMLUN.h>
#endif
#include <comocn.h>
#ifdef isotrac
#include <comsoi.h>
#include <comsoixy.h>
#include <comocnxy.h>
#include <comicexy.h>
#endif
      dimension nhisi(nmark)
      logical histlast
c
      character*8 ctmp
      character cmon(12)*3
      save cmon
      data cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     *           'JUL','AUG','SEP','OCT','NOV','DEC' /
c
      namelist /inputlsx/
     *  msnam_a,   ndens_a,   nhisi_a,   ndhis_a,   histlast_a,
     *  cushnam_a, stomfac,   delsst,
     *  z0ocnh,    ficemaxn,  ficemaxs,  ticrit,    hiceexp,
     *  cdhfac,    cdefac,    cdmfac,
     *  dslons,    dslats,    dsvegs,    dssois,
     *  dslonn,    dslatn,    dsvegn,    dssoin,
     *  imagefile, wsoi_init, penfr,     fbrinmax,  fbrinsto
#ifdef isotrac
     * ,hctm_i,    him_i,     hldm_i,    hlsm_i,
     *  ietest,    oceano18
#endif
c
c        Set optionally-entered array variables to null
c        (since trailing words are left alone by namelist,
c        cannot set defaults till after namelist is read)
c
      do 120 i=1,pfile
        msnam_a(i) = ' '
  120 continue
#ifdef rcm
        msnam_a(1) = 'LSX'
#endif

      do 122 i=1,nmark
        nhisi_a(i) = 0
  122 continue

      do 124 i=1,nushmx
        cushnam_a(i) = ' '
  124 continue
c
c        Set optionally-entered single variables to defaults
c
      ndens_a = 1
      ndhis_a = 999999
      histlast_a = histlast                  ! default = agcm's histlast
c
      stomfac = 1.
      delsst = 0.
      z0ocnh = 0.
      ficemaxn = 0.
      ficemaxs = 0.
      ticrit = -999.
      hiceexp = 0.5 ! was 3. in v2.0.a.1
c
c v2:
      cdhfac = 1.0
      cdefac = 1.35
      cdmfac = 1.0
c v3:
c     cdhfac = 0.6
c     cdefac = 0.6
c     cdmfac = 0.6
c
      dslons(1) = -181.
      dslons(2) = -181.
      dslats(1) = -91.
      dslats(2) = -91.
      dsvegs    = 0
      dssois    = 0
c
      dslonn(1) = -181.
      dslonn(2) = -181.
      dslatn(1) = -91.
      dslatn(2) = -91.
      dsvegn    = 0
      dssoin    = 0
c
      imagefile = ' '
      wsoi_init = 0.6
      penfr = 0.17
      fbrinmax = 0.25
      fbrinsto = 0.75
#ifdef rcm
c     no brine for rcm...sea ice prescribed, and brine -> Nan
      fbrinmax = 0.
      fbrinsto = 0.
#endif

c
#ifdef isotrac
      hctm_i = 150.
      him_i  = 100.
      hldm_i = 1000.
      hlsm_i = 250.
      ietest = -1   ! gcm-strat calc (Mathieu et al., App. test "Ez")
      oceano18 = -999.
#endif
c
c        Set some quantities in COMLUN (set in gcm), used in lsx
c        (some just for writing to LSX history file headers in writhlsx)
c
#ifdef rcm
      mspathin = ' '
      msnamin  = ' '
      mspath   = ' '
      do i=1,pfile
        msnamhis(i) = ' '
        msnamres(i) = ' '
      enddo
      reghis_a = ' '
      fatalwms = .false.
      hisnet = .true.
#endif
c
c        Read namelist inputlsx for lsx variables. 
c        Namelist input file has been opened in datget. 
c
      write(ioterm2,2000)
 2000 format(/' Reading namelist input (group name: INPUTLSX)...'/)
#ifdef rcm
      read inputlsx
#else
      read (iunamel,inputlsx)
      call closef (iunamel) 
#endif
c
c        Process lsx filename lists
c
c---------
#ifdef rcm
c---------
      nfile = 0
      do i=1,pfile
        if (msnam_a(i).eq.' ') goto 200
        nfile = i
      enddo
  200 continue

c        If only one name entered, use for all (default set above 
c        to 'LSX'...see histpost, where files are renamed based on date)

      if (nfile.eq.1) then
        do i=2,pfile
          msnam_a(i) = msnam_a(1) 
        enddo
        nfile = pfile
      endif

      if (nsrest.eq.1) then
        if (nuoff_a.eq.0) ifile = ifile + 1
      else 
        ifile = 1
      endif

      if (ifile.gt.nfile) then
         write (6,'(/a,i6/a,i6)')
     *     '*** Error: current LSX history file #    (ifile) = ',ifile,
     *     '    exceeds number of filenames (msnam_a, nfile) = ',nfile
         call endrun (-1)
      endif
c----
#else
c----
      call datnalsx
c-----
#endif
c-----
c
c        Process lsx history write time lists
c
      call datimlsx (dtgcm, nhisi, ndhis)
c
c        Process lsx history field lists
c
      call datfdlsx
c
c        Estimate size of each history file in words nusiz_a (for openf)
c
c     iwpf = estimated # of history writes on each history file
c     (reduced if the length of run is less than 1 year)
c     iwpf = (nhistn_a/max(nnexfn,1))  ! nnexfn not avail in lsx
      iwpf = (nhistn_a)
     *       * min (1., (nstop-nrstrt+1)*dtgcm/(86400.*365.))
      if (nuldat_a.gt.0) then
c       add 1000 for an (over)estimate of one set of header records
        nusiz_a = (1000 + nurecl_a*nuldat_a) * max(iwpf,1)
      else
        nusiz_a = 0
      endif
c
c        Open lsx history file.
c
      if (nsrest.le.0) then
c
c           Open new history file for cold start or initial start.
c
         call hisselsx (1)
         call openf (ndata_a, msnam_a(ifile), mspath, 1, nusiz_a)
c
      else
c
         if (nuoff_a.eq.0) then
c
c           Open new history file for true restart
c
            call hisselsx (1)
            call openf (ndata_a, msnam_a(ifile), mspath, 1, nusiz_a)
c
         else
c
c           Open old history file to append to for restart
c
#ifndef rcm
            if (msnam_ax.ne.msnam_a(ifile)) write(ioterm2,4000)
     *         ifile, ifile, msnam_a(ifile), msnam_ax
 4000       format(
     *        /' *** Warning ***'
     *        /'     The current lsx history file to be appended to',
     *               ' by this run'
     *        /'     does not agree with its previously recorded name:'
     *        /'     Current position in lsx hist filename list = ',i8
     *        /'     Current lsx history file MSNAM_A (',i4,')    = ',a
     *        /'     Previous filename in that position         = ',a/)
#endif

            if (.not.hisnet) then
              call acquiref (msnam_a(ifile), mspath)
              call openf (ndata_a, msnam_a(ifile), mspath, 0, 0)
c
              ctmp = 'xxxxxxxx'
              read (ndata_a, err=4100, iostat=ios) ctmp
 4100         if (ctmp(5:8).eq.' ') ctmp(5:8) = 'CRAY'
              if (ctmp(5:8).ne.typehis_a) then
                write(ioterm2,4102) ifile, msnam_a(ifile),
     *                              ctmp(5:8), typehis_a
 4102           format(
     *            /' *** The current lsx history file to be',
     *                   ' appended to by this run'
     *            /'     has a different machine numerical format',
     *                   ' than the current program:'
     *            /'     MSNAM_A (',i4,')             = ',a
     *            /'     Machine format of file       = ',a
     *            /'     Machine format of program    = ',a/)
                call endrun (-1)
              endif
 
              if (ctmp(1:4).ne.'LSX') then
                write(ioterm2,4104) ifile, msnam_a(ifile), ctmp(1:4)
 4104           format(
     *            /' *** The current lsx history file to be',
     *                   ' appended to by this run'
     *            /'     has an incorrect first-4-byte label:'
     *            /'     MSNAM_A (',i4,')             = ',a
     *            /'     First word (should be LSX) = ',a /)
                call endrun (-1)
              endif

              call closef (ndata_a)
            endif
c
#ifdef rcm
            write (ioterm2,4106) 
     *        msnam_a(ifile)(1:lenchr(msnam_a(ifile))), nuoff_a
 4106       format (/'Opening previous LSX history file ',a,
     *               ' for append (nuoff_a=', i6, ')'/)
#endif
            call openf (ndata_a, msnam_a(ifile), mspath, -1, 0)
            call hisselsx (-1)
c
         endif
c
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine datnalsx
c
c        Set default lsx filenames, derived from
c        ccm history filenames msnamhis(1 to nfile).
c        Similar to ccm subr datnames in init.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
c
c        Possible "*" expansion for lsx file lists
c
      call starlist (msnam_a, pfile, 'LSX history')
c
c        Set default names from agcm history list
c
      do 140 i=1,nfile
        if (msnam_a(i).eq.' ') then
          do 142 l=len(msnamhis(i))-3,1,-1
            lenh = l
            if (msnamhis(i)(l:l).ne.' ') goto 144
  142     continue
  144     if (msnam_a(i).eq.' ') msnam_a(i) = msnamhis(i)(1:lenh)//'LSX'
        endif
  140 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine datimlsx (dtgcm, nhisi)
c
c        Processes user-entered lsx history-write times.
c        Similar to ccm subr dattimes in init.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
      dimension nhisi(nmark)
c
      dimension ndaydef(12)
      data ndaydef /  31,  59,  90, 120, 151, 181,
     *               212, 243, 273, 304, 334, 365  /
c
      ntpy = nint(365*86400./dtgcm)
c
c        Set default end times (nhisi_a) to agcm's nhisi; since nhisi is
c        in units of timesteps, set nhisi_a = -nhisi to handle ok below.
c
      if (nhisi_a(1).eq.0) then
c        do 100 i=1,12
c           nhisi_a(i) = ndaydef(i)
c 100    continue
         do 100 i=1,nmark
            nhisi_a(i) = -nhisi(i)
  100    continue
      endif
c
c        Count number of periods (nhistn_a), and change end times to
c        units of timesteps.
c
      nhistn_a = 0
      kount = 0
      do 110 i=1,nmark
         if (nhisi_a(i).eq.0) goto 112
c        convert any "mm/dd" input to day #
         call dateconv (nhisi_a(i), i, 'NHISI_A', kount)
         nhistn_a = i
         if (nhisi_a(i).gt.0) then
            nhisi_a(i) = max ( nint(nhisi_a(i)*86400./dtgcm), 1 )
         else if (nhisi_a(i).lt.0) then
            nhisi_a(i) = -nhisi_a(i)
         endif
  110 continue
  112 continue
c
c        If just one (non "mm/dd") write time has been entered,
c        convert to multiples of that time from beginning of year
c
      if (nhistn_a.eq.1 .and. kount.eq.0) then
         nhistn_a = ntpy / nhisi_a(1)
         nhistn_a = min (nmark, max (1, nhistn_a))
         if (nhistn_a.gt.1) then
            do 114 i=2,nhistn_a
               nhisi_a(i) = nhisi_a(1)*i
  114       continue
         endif
      endif
c
c        Set end times nhist_a(*,2) depending on nhisi_a,
c        constraining them to be within 1-year span
c
      do 120 i=1,nhistn_a
        nhist_a(i,2)= max (1, min (ntpy, nhisi_a(i)))
  120 continue
c
c        Set start times nhist_a(*,1) depending on ndhis_a. If ndhis_a
c        is not entered (=999999), set to agcm's ndhis. If ndhis_a is
c        zero,  start time = previous end time + 1 timestep,
c        ie, accumulated fields will be accumulating 100% of the time.
c        Convert entered user-entered ndhis_a (if positive) to units of
c        timesteps.
c
      if (ndhis_a.eq.999999) then
#ifdef rcm
         ndhis_a = 0
#else
         ndhis_a = ndhis
#endif
      else if (ndhis_a.gt.0) then
         ndhis_a = nint (ndhis_a*86400./dtgcm)
      else if (ndhis_a.lt.0) then
         ndhis_a = -ndhis_a
      endif
c
      if (ndhis_a.eq.0 .and. nhistn_a.gt.0) then
         nhist_a(1,1) = nhist_a(nhistn_a,2) + 1
         do 130 i=2,nhistn_a
            nhist_a(i,1) = nhist_a(i-1,2) + 1
  130    continue
      else
         do 132 i=1,nhistn_a
            nhist_a(i,1) = nhist_a(i,2) - ndhis_a + 1
  132    continue
      endif
 
c        Constraint start times to not overlap into previous periods
 
      if (nhist_a(1,1).le.0) nhist_a(1,1) = nhist_a(1,1) + ntpy
      if (nhist_a(1,1).gt.nhist_a(1,2) .and. nhistn_a.gt.0)
     *  nhist_a(1,1) = max ( nhist_a(nhistn_a,2) + 1, nhist_a(1,1) )
      do 134 i=2,nhistn_a
         nhist_a(i,1) = max ( nhist_a(i-1,2) + 1, nhist_a(i,1) )
  134 continue
c
c        Check reasonableness of history periods
c
      do 140 i=2,nhistn_a
        if (nhist_a(i,1).gt.nhist_a(i,2)) then
          write(ioterm2,1400) i-1, nhisi_a(i-1)*dtgcm/86400.,
     *                        i,   nhisi_a(i)  *dtgcm/86400.
 1400     format
     *      (/' *** lsx history-write times NHISI_A must increase',
     *              ' in time,'
     *       /'     and fall between 00:00 Jan 1 and 24:00 Dec 31:'
     *       /'     NHISI_A(',i3,') =',f8.3,'  NHISI_A(',i3,') =',f8.3)
           call endrun (-1)
        endif
  140 continue

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine datfdlsx
c
c        Processes user-entered lsx history field lists.
c        Similar to ccm subr datfield in init.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
c
c        Initialize hard-coded lsx history field list
c
      call fldlsx
c
c        Count number of history fields entered, and transfer
c        from input array cushnam_a to cushis_a, mushfl_a.
c        If no fields entered (nushis_a=0), rely on "no-trip"
c        do-loops below (elsewhere, handled by nuldat_a=0).
c
      nushis_a = 0
      do 200 i=1,nushmx
         if (cushnam_a(i).eq.' ') goto 202
         call upper (cushnam_a(i))
         nushis_a = i
         n = lenchr(cushnam_a(i))
         if (cushnam_a(i)(n:n).eq.'*') then
           mushfl_a(i) = 0                            ! instantaneous
           cushis_a(i) = cushnam_a(i)(1:n-1)
         else
           mushfl_a(i) = 1                            ! time-accumulated
           cushis_a(i) = cushnam_a(i)
         endif
  200 continue
  202 continue
c
c         Form index from user-entered list to program list,
c         and check for duplicate entries
c
      ier = 0
c
      do 300 i=1,nushis_a
         do 302 ip = 1,nindex
            if(cushis_a(i).eq.cindex(1,ip)) then
               mushis_a(i) = ip
               goto 300
            endif
  302    continue
         ier = ier + 1
         write(ioterm2,*) '*** Requested lsx history field ',
     *                    cushis_a(i),' is not a valid field label.'
  300 continue
c
      do 310 i=1,nushis_a-1
         do 312 j=i+1,nushis_a
            if (cushis_a(i).eq.cushis_a(j)) then
               ier = ier + 1
               write(ioterm2,*) '*** Duplicate lsx history',
     *           ' fields ',cushis_a(i),' at entry numbers ',i,' and ',j
            endif
  312    continue
  310 continue
c
      if (ier.gt.0) call endrun (-1)
c
c        Set "record" positions for each field (relative to start of
c        storage area) and the total records nuldat_a and nuldat_atot.
c        Requested fields are allocated forwards from the start of
c        the storage area (total=nuldat_a), and daily min/max fields are
c        allocated 2 extra records each (for holding min and max)
c        after those (combined total=nuldat_atot).
c
      nuldat_a = 0
      do 400 i = 1,nushis_a
         ip = mushis_a(i)
         mindex_a(3,ip) = nuldat_a + 1
         mindex_a(4,ip) = nuldat_a + mindex_a(2,ip)
         nuldat_a = mindex_a(4,ip)
         mindfl_a(ip) = mushfl_a(i)
  400 continue
c
c        Add any extra daily min/max records onto the end
c
      nuldat_atot = nuldat_a
      do 410 i = 1,nushis_a
         ip = mushis_a(i)
         if (mindex_a(5,ip).eq.1) then
           mindex_a(6,ip) = nuldat_atot + 1
           nuldat_atot = nuldat_atot + 2
         endif
  410 continue
 
c        Initialize trailing min/max recs (1st will hold min, 2nd max)
 
      do 412 i = 1,nushis_a
         ip = mushis_a(i)
         if (mindex_a(5,ip).eq.1) then
           call resetr (workl, nlon*nlat,  1.e20)
           call puthrlsx (mindex_a(6,ip),   workl, ip)
           call resetr (workl, nlon*nlat, -1.e20)
           call puthrlsx (mindex_a(6,ip)+1, workl, ip)
         endif
  412 continue
c
c        Set number of history header recs nulhdr_a, constrain data
c        packing density ndens_a, set hist-file record length
c        nurecl_a (determined by length of data records, = 1 label
c        word + unpacked/packed horizontal gridsize).
c
      nulhdr_a = 4
      if (nuldat_a.eq.0) nulhdr_a = 0
 
      iz = ndens_a
 
      ndens_a = max (1, min(4,ndens_a))
#ifdef nopack
      ndens_a = 1
#endif
 
      if (ndens_a.ne.iz) write (ioterm2,500) iz, ndens_a
  500 format(/' *** Warning: entered lsx history-file packing density',
     *       /'     NDENS_A =',i2,' has been reset to NDENS_A =',i2/)
 
      if (ndens_a.eq.1) then
        nurecl_a = 1 + nlon*nlat
      else
        nurecl_a = 1 + (nlon*nlat + ndens_a - 1)/ndens_a + 2
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fldlsx
c
c        Initializes hard-coded program list of all lsx history fields.
c        Similar to ccm subr fldlst in hisres.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
c
      do 5 j=1,nushmx
        mindex_a(1,j) = 0
        mindex_a(2,j) = 1
        mindex_a(3,j) = 0
        mindex_a(4,j) = 0
        mindex_a(5,j) = 0
        mindex_a(6,j) = 0
    5 continue
c
      ip = 0
c
c01            LMASK           lsx surface type(1=land,2=icesheet,2=ocn)
      ip = ip + 1
      cindex(1,ip) = 'LMASK'
      cindex(2,ip) = ' '
c02            TS              surface temperature (areally averaged)
      ip = ip + 1
      cindex(1,ip) = 'TS'
      cindex(2,ip) = 'deg K'
c03-08         T[1-6]          soil/icesheet layer temps, top to bottom
      do 3 m=1,6
        ip = ip + 1
        write( cindex(1,ip), "('T',i1)" ) m
        cindex(2,ip) = 'deg K'
    3 continue
c09-14         WET[1-6]        soil-layer liquid content, top to bottom
      do 9 m=1,6
        ip = ip + 1
        write( cindex(1,ip), "('WET',i1)" ) m
        cindex(2,ip) = 'fraction'
    9 continue
c15-20         WICE[1-6]       soil-layer ice content, top to bottom
      do 15 m=1,6
        ip = ip + 1
        write( cindex(1,ip), "('WICE',i1)" ) m
        cindex(2,ip) = 'fraction'
   15 continue
c21            SNOWH           snow thickness (liquid equivalent depth)
      ip = ip + 1
      cindex(1,ip) = 'SNOWH'
      cindex(2,ip) = 'm'
c22            SNOWF           fractional snow cover
      ip = ip + 1
      cindex(1,ip) = 'SNOWF'
      cindex(2,ip) = 'fraction'
c23            ICEH            seaice thickness
      ip = ip + 1
      cindex(1,ip) = 'ICEH'
      cindex(2,ip) = 'm'
c24            ICEF            fractional seaice cover
      ip = ip + 1
      cindex(1,ip) = 'ICEF'
      cindex(2,ip) = 'fraction'
c25            TA              lowest-level agcm temperature (from ccm)
      ip = ip + 1
      cindex(1,ip) = 'TA'
      cindex(2,ip) = 'deg K'
c26            TUL             upper leaf temperature
      ip = ip + 1
      cindex(1,ip) = 'TUL'
      cindex(2,ip) = 'deg K'
c27            TUS             upper stem temperature
      ip = ip + 1
      cindex(1,ip) = 'TUS'
      cindex(2,ip) = 'deg K'
c28            TLV             lower vegetation temperature
      ip = ip + 1
      cindex(1,ip) = 'TLV'
      cindex(2,ip) = 'deg K'
c29            RUNOFF          surface runoff
      ip = ip + 1
      cindex(1,ip) = 'RUNOFF'
      cindex(2,ip) = 'm/s'
c30            DRAIN           sub-surface drainage
      ip = ip + 1
      cindex(1,ip) = 'DRAIN'
      cindex(2,ip) = 'm/s'
c31            TRANSP          transpiration rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'TRANSP'
      cindex(2,ip) = 'm/s'
c32            INFIL           soil infiltration rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'INFIL'
      cindex(2,ip) = 'm/s'
c33            ROUGHL          roughness length (felt by agcm)
      ip = ip + 1
      cindex(1,ip) = 'ROUGHL'
      cindex(2,ip) = 'm'
c34            DSBIOME         Dorman-Sellers biome index
      ip = ip + 1
      cindex(1,ip) = 'DSBIOME'
      cindex(2,ip) = ' '
c35            TAUX            Eastward wind stress
      ip = ip + 1
      cindex(1,ip) = 'TAUX'
      cindex(2,ip) = 'N/m2'
c36            TAUY            Northward wind stress
      ip = ip + 1
      cindex(1,ip) = 'TAUY'
      cindex(2,ip) = 'N/m2'
c37            H2OINT          Net intercepted h2o (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'H2OINT'
      cindex(2,ip) = 'm'
c38            UICE            Eastward sea-ice velocity
      ip = ip + 1
      cindex(1,ip) = 'UICE'
      cindex(2,ip) = 'm/s'
c39            VICE            Northward sea-ice velocity
      ip = ip + 1
      cindex(1,ip) = 'VICE'
      cindex(2,ip) = 'm/s'
c40            AVAP            Net upward water vapor flux to agcm
      ip = ip + 1
      cindex(1,ip) = 'AVAP'
      cindex(2,ip) = 'm/s'
c41            QFLUX           Oceanic heat flux convergence
      ip = ip + 1
      cindex(1,ip) = 'QFLUX'
      cindex(2,ip) = 'W/m2'
c42            WISCVEG         Wisconsin vegetation category
      ip = ip + 1
      cindex(1,ip) = 'WISCVEG'
      cindex(2,ip) = ' '
c43            BRINE           Seaice brine resevoir/heat to melt column
      ip = ip + 1
      cindex(1,ip) = 'BRINE'
      cindex(2,ip) = 'fraction'
c44            FICEADV         Advective d(sea-ice fraction)/dt
      ip = ip + 1
      cindex(1,ip) = 'FICEADV'
      cindex(2,ip) = '1/s'
c45            HICEADV         Advective d(sea-ice thickness)/dt
      ip = ip + 1
      cindex(1,ip) = 'HICEADV'
      cindex(2,ip) = 'm/s'
c46            HOCEAN          Ocean mixed-layer depth (m)
      ip = ip + 1
      cindex(1,ip) = 'HOCEAN'
      cindex(2,ip) = 'm'
c47            TS_M            daily min/max of surface temp
      ip = ip + 1
      cindex(1,ip) = 'TS_M'
      cindex(2,ip) = 'deg K'
      mindex_a(2,ip) = 2
      mindex_a(5,ip) = 1
c48            TA_M            daily min/max of lowest agcm temp
      ip = ip + 1
      cindex(1,ip) = 'TA_M'
      cindex(2,ip) = 'deg K'
      mindex_a(2,ip) = 2
      mindex_a(5,ip) = 1
c49            TS2             2-m screen temperature
      ip = ip + 1
      cindex(1,ip) = 'TS2'
      cindex(2,ip) = 'deg K'
c50            TS2_M           daily min/max of 2-m screen temp
      ip = ip + 1
      cindex(1,ip) = 'TS2_M'
      cindex(2,ip) = 'deg K'
      mindex_a(2,ip) = 2
      mindex_a(5,ip) = 1
c51            RELHUM          relative humidity at lowest agcm level
      ip = ip + 1
      cindex(1,ip) = 'RELHUM'
      cindex(2,ip) = 'fraction'
c52            TOCEAN          ocean mixed-layer temperature
      ip = ip + 1
      cindex(1,ip) = 'TOCEAN'
      cindex(2,ip) = 'K'
c53            INTVAP          net evaporation from intercepted h2o
      ip = ip + 1
      cindex(1,ip) = 'INTVAP'
      cindex(2,ip) = 'm/s'
c54            SURVAP          net evaporation from soil/ice/snow/ocean
      ip = ip + 1
      cindex(1,ip) = 'SURVAP'
      cindex(2,ip) = 'm/s'
c55            FOGRAT          net rate of canopy-air fog formation
      ip = ip + 1
      cindex(1,ip) = 'FOGRAT'
      cindex(2,ip) = 'm/s'
c56            TAUXO            Eastward wind stress, ocean sfc(blw ice)
      ip = ip + 1
      cindex(1,ip) = 'TAUXO'
      cindex(2,ip) = 'N/m2'
c57            TAUYO            Northwd wind stress, ocean sfc(blw ice)
      ip = ip + 1
      cindex(1,ip) = 'TAUYO'
      cindex(2,ip) = 'N/m2'
c58            WFLOALL         Dwnd soil h2o flow (same as WFLO[1-7])
      ip = ip + 1
      cindex(1,ip) = 'WFLOALL'
      cindex(2,ip) = 'kg/m2/s' 
      mindex_a(2,ip) = nsoilay + 1
c59-61         TICE[1-3]       sea-ice layer temps, top to bottom
      do 59 m=1,3
        ip = ip + 1
        write( cindex(1,ip), "('TICE',i1)" ) m
        cindex(2,ip) = 'deg K'
   59 continue
c62            FWATER          prescribed open water fraction (0-1)
      ip = ip + 1
      cindex(1,ip) = 'FWATER'
      cindex(2,ip) = 'fraction'
c63            LAIU            Net one-sided LAI, upper story,*FUXY
      ip = ip + 1
      cindex(1,ip) = 'LAIU'
      cindex(2,ip) = ' '
c64            LAIL            Net one-sided LAI, lower story,*FLXY
      ip = ip + 1
      cindex(1,ip) = 'LAIL'
      cindex(2,ip) = ' '
c65            LAIT            Net one-sided LAI, both stories,*F[U,L]XY
      ip = ip + 1
      cindex(1,ip) = 'LAIT'
      cindex(2,ip) = ' '
c66            SAIU            Net one-sided SAI, upper story,*FUXY
      ip = ip + 1
      cindex(1,ip) = 'SAIU'
      cindex(2,ip) = ' '
c67            SAIL            Net one-sided SAI, lower story,*FLXY
      ip = ip + 1
      cindex(1,ip) = 'SAIL'
      cindex(2,ip) = ' '
c68            SAIT            Net one-sided SAI, both stories,*F[U,L]XY
      ip = ip + 1
      cindex(1,ip) = 'SAIT'
      cindex(2,ip) = ' '
c69            NPPU            IBIS npp of C per upper-story leaf area
      ip = ip + 1
      cindex(1,ip) = 'NPPU'
      cindex(2,ip) = 'kg/m2/yr'
c70            NPPL            IBIS npp of C per lower-story leaf area
      ip = ip + 1
      cindex(1,ip) = 'NPPL'
      cindex(2,ip) = 'kg/m2/yr'
c71            GPPU            IBIS gpp of C per upper-story leaf area
      ip = ip + 1
      cindex(1,ip) = 'GPPU'
      cindex(2,ip) = 'kg/m2/yr'
c72            GPPL            IBIS gpp of C per lower-story leaf area
      ip = ip + 1
      cindex(1,ip) = 'GPPL'
      cindex(2,ip) = 'kg/m2/yr'
c73            GPP             IBIS total gpp per overall area
      ip = ip + 1
      cindex(1,ip) = 'GPP'
      cindex(2,ip) = 'kg/m2/yr'
c74            STOMCU          IBIS stomatal conductance, upper leaves
      ip = ip + 1
      cindex(1,ip) = 'STOMCU'
      cindex(2,ip) = 'm/s'
c75            STOMCL          IBIS stomatal conductance, lower leaves
      ip = ip + 1
      cindex(1,ip) = 'STOMCL'
      cindex(2,ip) = 'm/s'
c76            STRESSU         soil moisture stress factor, upper veg
      ip = ip + 1
      cindex(1,ip) = 'STRESSU'
      cindex(2,ip) = 'fraction'
c77            STRESSL         soil moisture stress factor, lower veg
      ip = ip + 1
      cindex(1,ip) = 'STRESSL'
      cindex(2,ip) = 'fraction'
c78            RAINFA          rainfall rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'RAINFA'
      cindex(2,ip) = 'm/s'
c79            SNOWFA          snowfall rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'SNOWFA'
      cindex(2,ip) = 'm/s'
c80            PRECIP          precipitation rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'PRECIP'
      cindex(2,ip) = 'm/s'
c81            SOCEAN          ocean (ogcm) surface salinity (fraction)
      ip = ip + 1
      cindex(1,ip) = 'SOCEAN'
      cindex(2,ip) = 'fraction'
c82            H2OSOIL         total h2o in soil column(liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'H2OSOIL'
      cindex(2,ip) = 'm'
c83            PORESOIL        total pore space in soil column (m)
      ip = ip + 1
      cindex(1,ip) = 'PORESOIL'
      cindex(2,ip) = 'm'
c84            H2O30           total h2o in top 30 cm of soil (m)
      ip = ip + 1
      cindex(1,ip) = 'H2O30'
      cindex(2,ip) = ' '
c85            PORE30          total pore space in top 30 cm of soil (m)
      ip = ip + 1
      cindex(1,ip) = 'PORE30'
      cindex(2,ip) = ' '
c86            HFLXOCN         net heat flux into ocean (below sea ice)
      ip = ip + 1
      cindex(1,ip) = 'HFLXOCN'
      cindex(2,ip) = 'W/m2'
c87            PMEOCN          net P-E into ocean (below sea ice)
      ip = ip + 1
      cindex(1,ip) = 'PMEOCN'
      cindex(2,ip) = 'm/s'
c88-93         POROS[1-6]      soil layer porosity, top to bottom
      do 88 m=1,6
        ip = ip + 1
        write( cindex(1,ip), "('POROS',i1)" ) m
        cindex(2,ip) = ' '
   88 continue
c94            TALL            soil layer temperatutes (same as T[1-6])
      ip = ip + 1
      cindex(1,ip) = 'TSOIALL'
      cindex(2,ip) = ' deg K'
      mindex_a(2,ip) = nsoilay
c95            WETALL          soil layer moisture  (same as WET[1-6])
      ip = ip + 1
      cindex(1,ip) = 'WETALL'
      cindex(2,ip) = ' '
      mindex_a(2,ip) = nsoilay
c96            WICEALL         soil layer ice (same as WICE[1-6])
      ip = ip + 1
      cindex(1,ip) = 'WICEALL'
      cindex(2,ip) = ' '
      mindex_a(2,ip) = nsoilay
c97            POROSALL        soil layer porosities(same as POROS[1-6])
      ip = ip + 1
      cindex(1,ip) = 'POROSALL'
      cindex(2,ip) = ' '
      mindex_a(2,ip) = nsoilay

c---------------------
#if defined  (isotrac)
c---------------------
c98            HCICE           sea-ice surface reserv. (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'HCICE'
      cindex(2,ip) = 'kg/m2'
c99            CICE_O          18O ratio for sea-ice surface reservoir
      ip = ip + 1
      cindex(1,ip) = 'CICE_O'
      cindex(2,ip) = ' '
c100           CICE_H          2H ratio for sea-ice surface reservoir
      ip = ip + 1
      cindex(1,ip) = 'CICE_H'
      cindex(2,ip) = ' '
c101           HTOPSOIL        upper soil reservoir (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'HTOPSOIL'
      cindex(2,ip) = 'kg/m2'
c102           TOPSOI_O        18O ratio for upper soil reservoir
      ip = ip + 1
      cindex(1,ip) = 'TOPSOI_O'
      cindex(2,ip) = ' '
c103           TOPSOI_H        2H ratio for upper soil reservoir
      ip = ip + 1
      cindex(1,ip) = 'TOPSOI_H'
      cindex(2,ip) = ' '
c104           HBOTSOIL        lower soil reservoir (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'HBOTSOIL'
      cindex(2,ip) = 'kg/m2'
c105           BOTSOI_O        18O ratio for lower soil reservoir
      ip = ip + 1
      cindex(1,ip) = 'BOTSOI_O'
      cindex(2,ip) = ' '
c106           BOTSOI_H        2H ratio for lower soil reservoir
      ip = ip + 1
      cindex(1,ip) = 'BOTSOI_H'
      cindex(2,ip) = ' '
c107           HLAKE           lake reservoir (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'HLAKE'
      cindex(2,ip) = 'kg/m2'
c108           LAKE_O          18O ratio for lake reservoir & ocean
      ip = ip + 1
      cindex(1,ip) = 'LAKE_O'
      cindex(2,ip) = ' '
c109           LAKE_H          2H ratio for lake reservoir & ocean
      ip = ip + 1
      cindex(1,ip) = 'LAKE_H'
      cindex(2,ip) = ' '
c110           OCEAN_I         ocean type or isotope's
      ip = ip + 1
      cindex(1,ip) = 'OCEAN_I'
      cindex(2,ip) = ' '
c-----
#else
c-----
c98-110       spare
      do i=98,110
        ip = ip + 1
        cindex(1,ip) = 'spare'
        cindex(2,ip) = ' '
      enddo
c-----
#endif
c-----

c111           UOCEAN          ocean (ogcm) surface eastward velocity
      ip = ip + 1
      cindex(1,ip) = 'UOCEAN'
      cindex(2,ip) = 'm/s'
c112           VOCEAN          ocean (ogcm) surface northward velocity
      ip = ip + 1
      cindex(1,ip) = 'VOCEAN'
      cindex(2,ip) = 'm/s'
c113           TAUXOCN2        ocean (ogcm) eastward wind stress
      ip = ip + 1
      cindex(1,ip) = 'TAUXOCN2'
      cindex(2,ip) = 'N/m2'
c114           TAUYOCN2        ocean (ogcm) northward wind stress
      ip = ip + 1
      cindex(1,ip) = 'TAUYOCN2'
      cindex(2,ip) = 'N/m2'
c115           HFLXOCN2        ocean (ogcm) net surface heat flux
      ip = ip + 1
      cindex(1,ip) = 'HFLXOCN2'
      cindex(2,ip) = 'W/m2'
c116           PMEOCN2        ocean (ogcm) net p-e (+runof+calv)
      ip = ip + 1
      cindex(1,ip) = 'PMEOCN2'
      cindex(2,ip) = 'm/s'
c117           RUNOCN2        ocean (ogcm) runoff
      ip = ip + 1
      cindex(1,ip) = 'RUNOCN2'
      cindex(2,ip) = 'm/s'
c118           CALVOCN2       ocean (ogcm) calving
      ip = ip + 1
      cindex(1,ip) = 'CALVOCN2'
      cindex(2,ip) = 'm/s'
c119           FO18OCN2       ocean (ogcm) delta18O flux
      ip = ip + 1
      cindex(1,ip) = 'FO18OCN2'
      cindex(2,ip) = 'del*m/s'
c120           FDHOCN2        ocean (ogcm) deltaD flux
      ip = ip + 1
      cindex(1,ip) = 'FDHOCN2'
      cindex(2,ip) = 'del*m/s'
c121           RAINFA          rainfall rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'RAINFA'
      cindex(2,ip) = 'm/s'
c122           SNOWFA          snowfall rate (liq equiv depth)
      ip = ip + 1
      cindex(1,ip) = 'SNOWFA'
      cindex(2,ip) = 'm/s'
c123           BIOMEN         Biome3 or Biome4 biome index
      ip = ip + 1
      cindex(1,ip) = 'BIOMEN'
      cindex(2,ip) = ' '
c124-130        WFLO[1-7]     Downward liquid h2o across soil interfaces
      do 124 m=1,7
        ip = ip + 1
        write( cindex(1,ip), "('WFLO',i1)" ) m
        cindex(2,ip) = 'kg/m2/s'
  124 continue

c----------------
#if defined (rcm)
c----------------
c131           LON2D          longitude (deg E) on cross grid
      ip = ip + 1
      cindex(1,ip) = 'LON2D'
      cindex(2,ip) = 'deg E'
c132           LAT2D          latitude  (deg N) on cross grid
      ip = ip + 1
      cindex(1,ip) = 'LAT2D'
      cindex(2,ip) = 'deg N'
c133
      ip = ip + 1
      cindex(1,ip) = 'TS10'
      cindex(2,ip) = 'deg K'
c134
      ip = ip + 1
      cindex(1,ip) = 'QS10'
      cindex(2,ip) = 'kg/kg'
c135
      ip = ip + 1
      cindex(1,ip) = 'US10'
      cindex(2,ip) = 'm/s'
c136
      ip = ip + 1
      cindex(1,ip) = 'VS10'
      cindex(2,ip) = 'm/s'
c137
      ip = ip + 1
      cindex(1,ip) = 'PSURF'
      cindex(2,ip) = 'N/m2'
c138
      ip = ip + 1
      cindex(1,ip) = 'SOL_D_V'
      cindex(2,ip) = 'W/m2'
c139
      ip = ip + 1
      cindex(1,ip) = 'SOL_D_N'
      cindex(2,ip) = 'W/m2'
c140
      ip = ip + 1
      cindex(1,ip) = 'SOL_I_V'
      cindex(2,ip) = 'W/m2'
c141
      ip = ip + 1
      cindex(1,ip) = 'SOL_I_N'
      cindex(2,ip) = 'W/m2'
c142
      ip = ip + 1
      cindex(1,ip) = 'COSZEN'
      cindex(2,ip) = ' '
c143
      ip = ip + 1
      cindex(1,ip) = 'FIRA'
      cindex(2,ip) = ' '
c144
      ip = ip + 1
      cindex(1,ip) = 'TOTCLD'
      cindex(2,ip) = ' '
c145
      ip = ip + 1
      cindex(1,ip) = 'SOILH2O'
      cindex(2,ip) = ' '
c146
      ip = ip + 1
      cindex(1,ip) = 'SOILH2O3'
      cindex(2,ip) = ' '
c147
      ip = ip + 1
      cindex(1,ip) = 'SOILT'
      cindex(2,ip) = 'deg K'
c148
      ip = ip + 1
      cindex(1,ip) = 'SOILT3'
      cindex(2,ip) = 'deg K'
c----
#else
c----
c131-148     spare   
      do i=131,148
        ip = ip + 1
        cindex(1,ip) = 'spare' 
        cindex(2,ip) = ' '
      enddo   
c-----
#endif
c-----

c149
      ip = ip + 1
      cindex(1,ip) = 'TOPOGNOD'
      cindex(2,ip) = 'm'
c150
      ip = ip + 1
      cindex(1,ip) = 'TOPOGDEL'
      cindex(2,ip) = 'm'
c151
      ip = ip + 1
      cindex(1,ip) = 'DELEV'
      cindex(2,ip) = 'm'

      nindex = ip
      if (nindex.gt.nushmx) then
         write (ioterm2,9005) nindex, nushmx
 9005    format(' Not enough room in indices cindex, mindex for',
     *          ' program lsx history field list.'
     *         /' Program error in fldlsx. nindex,nushmx =',2i10)
         call endrun (-1)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine hisselsx (iflag)
c
c        Sets record-counter nuoff_a and write-counter nwrit_a,
c        and zero accumulation counters mindco_a, for lsx history
c        files. Similar to ccm subr hisset in hisres.f.
c
c        iflag = 1 for starting new file,
c                0 for ongoing file, increment for next write,
c               -1 for append to old file (nuoff_a, nwrit_a have
c                  been read from restart file.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
c
      if (iflag.eq.1) then
         nuoff_a = 0
         nwrit_a = 0
      else if (iflag.eq.0) then
         nuoff_a = nuoff_a + nulhdr_a + nuldat_a
         nwrit_a = nwrit_a + 1
      endif
c
      do 10 i=1,nushmx
        mindco_a(i) = 0
   10 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      logical function testhisl (iv)
 
c        Test if need to calculate a history field for this timestep.
c        This is a "shell" of puthilsx, with "testhisl = .true."
c        replacing any actual puthilsx operations. It is used just to
c        save cpu in histwrit.
 
c        Returns .true. if need to calculate this field (with program
c        index iv), or .false. if not.
 
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
 
      testhisl = .false.
 
c     not a requested field
      if (mindex_a(3,iv).eq.0) return
 
      if (mindex_a(5,iv).eq.1) then
c        daily min/max field and no accum/write in current day
         if (.not.frange) return
         testhisl = .true.
      endif
 
c     no accum/write this timestep
      if (.not.fahist_a) return
 
c     instantaneous field and no write this timestep
      if (mindfl_a(iv).eq.0 .and. .not.fwhist_a) return
 
      testhisl = .true.
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine puthilsx (arr, iv)
c
c        Writes one field in array arr to the lsx
c        history storage area. iv is the program's field index.
c        Each "record" in the storage area holds nlon*nlat values.
c        The storage area (accessed by gethrlsx,puthrlsx) holds data for
c        one history "write". Array arr is not changed. (Don't use
c        comgrd's workl since this is often used to pass arr.)
c        Similar to ccm subr puthis in hisres.f
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
c
      dimension arr(nlon,nlat,1)
      dimension tmpx(nlon,nlat)
c
c        Do nothing unless user has requested this field
c
      if (mindex_a(3,iv).eq.0) return
c
c        Do nothing if not in accumulation window
c
      if (.not.fahist_a) return
c
c        Do nothing for instantaneous field and not at write timestep
c
      if (mindfl_a(iv).eq.0 .and. .not.fwhist_a) return
c
c        Increment counter for accumulated fields
c
      if (mindfl_a(iv).eq.1) mindco_a(iv) = mindco_a(iv) + 1
 
c        Loop over "levels" (1, except 2 for min/max fields called
c        from putrange, and nsoilay for TSOIALL,WETALL,WICEALL,POROSALL)
 
      do 100 jk=1,mindex_a(2,iv)
c
c        Set storage-area "record" number
c
        irec = mindex_a(3,iv) + jk - 1
c
c          For non-first-time accumulated fields, read existing values
c          and add current values, else just copy current values
c
        if (mindfl_a(iv).eq.1 .and. mindco_a(iv).gt.1) then
          call gethrlsx (irec, tmpx, ivdum)
          do 60 jj=1,nlat
            do 62 ji=1,nlon
               tmpx(ji,jj) = tmpx(ji,jj) + arr(ji,jj,jk)
   62       continue
   60     continue
        else
          call scopy(nlon*nlat, arr(1,1,jk), 1, tmpx, 1)
        endif
c
c          Put values in storage area
c
        call puthrlsx (irec, tmpx, iv)
 
  100 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine putrange (arr, iv)
c
c        Called like puthilsx but for daily [min,max]
c        fields. Tracks current day's min and max using the
c        2 extra min,max storage area records for each
c        field. Then, if at 24:00 GMT or at a write timestep,
c        calls puthilsx to record the [min,max]
c        for the previous 24 hours, and re-initialize the 2 extra
c        min/max storage records to +/- 1.e20.
c
c        arr= current field values (supplied)
c        iv = program's field index (must be a daily min/max
c             field so that mindex_a(6,iv) has been set)
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
c
      dimension arr(nlon,nlat)
      dimension tmparr(nlon,nlat,2),
     *          tmpmin(nlon,nlat), tmpmax(nlon,nlat)
      equivalence (tmparr(1,1,1),tmpmin),  (tmparr(1,1,2),tmpmax)
c
c
c        Do nothing unless user has requested this field
c
      if (mindex_a(3,iv).eq.0) return
c
c        Check that field is a min/max type
c
      if (mindex_a(5,iv).ne.1) then
        write(ioterm2,5)
    5   format(/'*** Program Error ***',/
     *          '    putrange called for non-[min/max] field.',
     *          '    field index iv =',i6/)
        call endrun (-1)
      endif
 
c        Do nothing if current day does not contain "_a" accum/write
 
      if (.not.frange) return
 
c        Track current day's min and max (these trailing storage records
c        were first initialized in datfdlsx)
 
      call gethrlsx (mindex_a(6,iv),   tmpmin, ivdum)
      call gethrlsx (mindex_a(6,iv)+1, tmpmax, ivdum)
      do 10 jj=1,nlat
!DIR$ ivdep
        do 12 ji=1,nlon
          tmpmin(ji,jj) = min (tmpmin(ji,jj), arr(ji,jj))
          tmpmax(ji,jj) = max (tmpmax(ji,jj), arr(ji,jj))
   12   continue
   10 continue
      call puthrlsx (mindex_a(6,iv),   tmpmin, iv)
      call puthrlsx (mindex_a(6,iv)+1, tmpmax, iv)
 
c        If not at 24:00 GMT (frange24) or a write timestep (fwhist_a),
c        do nothing else.
 
      if (.not.(frange24.or.fwhist_a)) return
 
c        Else (remainder of subroutine), calls puthilsx
c        for min,max of previous 24 hrs. (Due to equivalences
c        above, tmparr(*,*,1) and (*,*,2) already hold the previous
c        24-hr min and max respectively).
 
      call puthilsx (tmparr, iv)
 
c        Re-initialize trailing min/max recs for next day's min,max
c        (use tmpmin/max just as temporary arrays)
 
      call resetr (tmpmin, nlon*nlat,  1.e20)
      call puthrlsx (mindex_a(6,iv),   tmpmin, iv)
      call resetr (tmpmax, nlon*nlat, -1.e20)
      call puthrlsx (mindex_a(6,iv)+1, tmpmax, iv)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine gethrlsx (irec, arr, iv)
 
c        Reads a field from the lsx history storage area.
c        Similar to ccm subr gethr in hisres.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
      dimension arr(nlon*nlat)
c     parameter (nhirec=50)
      parameter (nhirec=70)
      common /chirec/ hissav(nlon*nlat,nhirec), indivsav(nhirec)
 
      do 200 i=1,nlon*nlat
        arr(i) = hissav(i,irec)
  200 continue

      iv = indivsav(irec)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine puthrlsx (irec, arr, iv)
 
c        Writes a field to the lsx history storage area.
c        Similar to ccm subr gethr in hisres.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
      dimension arr(nlon*nlat)
      logical first
      data first /.true./
      save first
c     parameter (nhirec=50)
      parameter (nhirec=70)
      common /chirec/ hissav(nlon*nlat,nhirec), indivsav(nhirec)
 
      if (first) then
        if (nuldat_atot.gt.nhirec) then
          write(ioterm2,2000)
     *      nuldat_a, nuldat_atot-nuldat_a, nuldat_atot, nhirec
 2000     format(/' *** The number of requested lsx history fields',
     *           /'     is too large for the program memory buffer.'
     *           /'     # of requested lsx history fields     =',i5
     *           /'     # of buffers for daily min/max fields =',i5
     *           /'     Total                                 =',i5
     *           /'     Program limit                         =',i5
     *           /)
          call endrun (-1)
        endif

        do i=1,nhirec
          indivsav(i) = 0
        enddo
        first = .false.
      endif
 
      do 200 i=1,nlon*nlat
        hissav(i,irec) = arr(i)
  200 continue

      indivsav(irec) = iv
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine writhlsx (dtgcm,  nsrest, fsave,
     *                     mbdate, mcdate, mcsec, runtitle, nstep)
 
c If HISNET.EQ..F., (Namelist, in comlun), writes a history "write", 
c ie, one set of header recs and data recs, to the current lsx 
c history (binary stream) file.

c If HISNET.EQ..T., does same but to the current lsx *netCDF*
c history file. Uses variables (cusfielnet_a, numfielnet_a, nwrit_a) 
c in common (COMLUN.h). Adapted from convnet, August '00.
c These files are opened in openf, closed in closef (as for binary).

 
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commapl.h>
#include <comsoi.h>
c-----
      logical fsave
      character runtitle*80
      character crec*8, runtitlx*80
      character cdum*8, cdum16*16
      save cdum,cdum16
      data cdum,cdum16 /' ',' '/
      parameter (nlevdum=1)
      dimension sigdum(nlevdum), sigldum(nlevdum+1)
c-----
c for netCDF:
#ifdef NETCDF
#include <netcdf.inc>
#ifdef cray
      integer
#else
      integer*4 
#endif
     * lenattr, nlon4, nlat4, nsoilay4, nsoilayp4, ival4, idate4,
     * varid, rcode,                                  ! ncid_a in comlun
     * lonid, latid, zsoilid, zsoilpid, timid, minmaxid,
     * minmax(2), ndim, idim(4), start(4), count(4),
     * n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls. 
     * ncglobal4, ncunlim4,
     * nf_float4, nf_int4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_int4   = nf_int,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)

      save varid, rcode,
     *     lonid, latid, zsoilid, zsoilpid, timid, minmaxid,
     *     minmax, idim, start, count,
     *     n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15
      data minmax /1,2/
      data n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15
     *     /1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/
      character*8 choriz
      dimension ztmp(1000), zsod(nsoilay), zsodp(nsoilay+1)
#endif
c-----
 
c        If requested, put current file counter at end of runtitle
 
      runtitlx = runtitle
      i = lenchr(runtitle)
c     if (runtitle(i-2:i).eq.'...') write(runtitlx(i-2:i),'(i3.3)')ifile
      write(runtitlx(i+2:i+4),'(i3.3)') ifile

c==========================
      if (.not.hisnet) then
c==========================
 
      do 14 k=1,nlevdum
        sigdum(k)  = (k-.5)/nlevdum
        sigldum(k) = (k-1.)/nlevdum
   14 continue
      sigldum(nlevdum+1) = 1.
c     2x2 resolution corresponds to spectral T59 truncation
      ntrkdum = 59
      ntrmdum = 59
      ntrndum = 59

c        Regular binary files

c---------------------------------
c        Write "LSX" header record
c---------------------------------
 
      crec = 'LSX ' // typehis_a(1:4)
      nw = 16 + nlon + nlat + nlevdum + (nlevdum+1)
      write (ndata_a, err=900, iostat=ios)
     *  crec, nw, nulhdr_a, nuldat_a, ndens_a, nurecl_a,
     *  vershis_a, verslsx, mbdate, dtgcm,
     *  nlon, nlat, nlevdum, ntrkdum, ntrmdum, ntrndum,
     *  alon, alat, sigdum, sigldum
 
c----------------------------------
c        Write "FILE" header record
c----------------------------------
 
      crec = 'FILE'
      nw = 6
     *   + (  len(runtitlx)+len(mspathin)+len(msnamin)+len(mspath)
     *      + len(msnamhis(ifile)) + len(msnam_a(ifile))
     *      + len(cdum16) + len(msnamres(ifile))
     *     ) / 8
      write (ndata_a, err=900, iostat=ios)
     *  crec, nw, nfile, ifile, nsrest, cdum,
     *  runtitlx, mspathin, msnamin, mspath,
     *  msnamhis(ifile), msnam_a(ifile),
     *  cdum16, msnamres(ifile)
 
c-----------------------------------
c        Write "WRITE" header record
c-----------------------------------
 
      crec = 'WRITE'
      nw = 10
      write (ndata_a, err=900, iostat=ios)
     *  crec, nw, nwrit_a+1, ndhis_a, nsteph_a, nstep,
     *  nint(dtgcm)*(nsteph_a-1), nint(dtgcm)*nstep,
     *  mcdate, mcsec
 
c------------------------------------
c        Write "FIELDS" header record
c------------------------------------
 
      crec = 'FIELDS'
      nw = 3 + nushis_a*8
      write (ndata_a, err=900, iostat=ios)
     *  crec, nw, nushis_a,
     *  ( (cindex  (m,mushis_a(n)),m=1,2),
     *    (mindex_a(m,mushis_a(n)),m=1,4),
     *     mindfl_a(mushis_a(n)),
     *     mindco_a(mushis_a(n)),          n=1,nushis_a )
 
c--------------------------------------------------------------------
c        Write "DATALSX" data records. (If ndens_a > 1, record length
c        nurecl_a is < nlon*nlat, and the data is packed into workm).
c--------------------------------------------------------------------
 
      crec = 'DATALSX'
      do 100 irec = 1,nuldat_a
 
        call gethrlsx (irec, workl, iv)

c          For accumulated fields, divide by counter 
c          (all counters zero-ed by hisselsx after return to histpost)

        if (mindfl_a(iv).eq.1) then
           do jj=1,nlat
              do ji=1,nlon
                 workl(ji,jj) = workl(ji,jj) / mindco_a(iv)
              enddo
           enddo
        endif
 
        if(ndens_a.le.1) then
          write (ndata_a, err=910, iostat=ios) crec, workl
        else
          call packa (workl, nlon*nlat, workm1, ndens_a)
          write (ndata_a,err=910,iostat=ios) crec,
     *                                       (workm1(j),j=1,nurecl_a-1)
        endif
 
  100 continue
 
      return
 
c---------------
c Error messages
c---------------
 
  900 write(ioterm2,9000) crec, ios
 9000 format(/' *** Error in writing lsx history-file header record.'
     *       /'     crec   = ',a
     *       /'     iostat = ',i8
     *       /'     Stop in WRITHLSX.'/)
      call endrun (-1)
c
  910 write(ioterm2,9100) irec, ios, nurecl_a, ndens_a
 9100 format(/' *** Error in writing lsx history-file data record.'
     *       /'     irec     = ',i8
     *       /'     iostat   = ',i8
     *       /'     nurecl_a = ',i8
     *       /'     ndens_a  = ',i8
     *       /'     Stop in WRITHLSX.'/)
      call endrun (-1)

c=========
      else 
c=========
#ifdef NETCDF

c        netCDF files

c>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (nwrit_a.eq.0) then     ! # of previous writes on file (comlun)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>

c          Global attributes

        lenattr = lenchr(runtitlx)

        rcode = nf_put_att_text (ncid_a, ncglobal4, 'title',
     *                           lenattr, runtitlx)

        if (nlon.eq.32 .and. nlat.eq.16) then
          choriz = 'T10'
        else if (nlon.eq.48 .and. nlat.eq.40) then
          choriz = 'R15'
        else if (nlon.eq.96 .and. nlat.eq.48) then
          choriz = 'T31'
        else if (nlon.eq.128 .and. nlat.eq.64) then
          choriz = 'T42'
        else if (nlon.eq.180 .and. nlat.eq.90) then
          choriz = '2X2'
        else
          choriz = 'unknown'
        endif
        lenattr = lenchr(choriz)
        rcode = nf_put_att_text (ncid_a, ncglobal4, 'horiz_grid',
     *                           lenattr, choriz)

        rcode = nf_put_att_text (ncid_a, ncglobal4, 'projection',
     *                           n6, 'latlon')

        nlon4 = nlon
        nlat4 = nlat
        nsoilay4 = nsoilay
        nsoilayp4 = nsoilay + 1

c          Longitude dimension

        rcode = nf_def_dim (ncid_a, 'lon', nlon4, lonid)
        rcode = nf_def_var (ncid_a, 'lon', nf_float4, n1, lonid, varid)
        rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                           n9, 'longitude')
        rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                           n12, 'degrees east')
        rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                           n4, 'f8.3') 
        rcode = nf_enddef (ncid_a)
        do i=1,nlon
#if defined (rcm)
          ztmp(i) = alon2d(i,(nlat+1)/2)
#else
          ztmp(i) = alon(i)*180./pi
#endif
        enddo
#ifdef cray
        rcode = nf_put_vara_real (ncid_a, varid, n1, nlon4, ztmp)
#else
        rcode = nf_put_vara_double (ncid_a, varid, n1, nlon4, ztmp)
#endif
        rcode = nf_redef (ncid_a)

c          Latitude dimension

        rcode = nf_def_dim (ncid_a, 'lat', nlat4, latid)
        rcode = nf_def_var (ncid_a, 'lat', nf_float4, n1, latid, varid)
        rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                           n8, 'latitude')
        rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                           n13, 'degrees north')
        rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                           n4, 'f8.3') 
        rcode = nf_enddef (ncid_a)
        do j=1,nlat
#if defined (rcm)
          ztmp(j) = alat2d((nlon+1)/2,j)
#else
          ztmp(j) = alat(j)*180./pi
#endif
        enddo
#ifdef cray
        rcode = nf_put_vara_real (ncid_a, varid, n1, nlat4, ztmp)
#else
        rcode = nf_put_vara_double (ncid_a, varid, n1, nlat4, ztmp)
#endif
        rcode = nf_redef (ncid_a)

c          Soil mid-layer (zsod) and interfaces (zsodp) dimensions

        zsodp(1) = 0.
        do k=2,nsoilay+1
          zsodp(k) = zsodp(k-1) + hsod(k-1)
        enddo
        do k=1,nsoilay
          zsod(k) = 0.5*(zsodp(k) + zsodp(k+1))
        enddo

        rcode = nf_def_dim (ncid_a, 'zsoil', nsoilay4, zsoilid)
        rcode = nf_def_var (ncid_a, 'zsoil',nf_float4, n1,zsoilid,varid)
        rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                           n15, 'soil mid-depths')
        rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                           n1, 'm')
        rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                           n4, 'f6.3') 
        rcode = nf_enddef (ncid_a)
#ifdef cray
        rcode = nf_put_vara_real (ncid_a, varid, n1, nsoilay4, zsod)
#else
        rcode = nf_put_vara_double (ncid_a, varid, n1, nsoilay4, zsod)
#endif
        rcode = nf_redef (ncid_a)

        rcode = nf_def_dim (ncid_a,'zsoilp', nsoilayp4, zsoilpid)
        rcode = nf_def_var (ncid_a,'zsoilp',nf_float4,n1,zsoilpid,varid)
        rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                           n15, 'soil interfaces')
        rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                           n1, 'm')
        rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                           n4, 'f6.3') 
        rcode = nf_enddef (ncid_a)
#ifdef cray
        rcode = nf_put_vara_real (ncid_a, varid, n1, nsoilayp4, zsodp)
#else
        rcode = nf_put_vara_double (ncid_a, varid, n1, nsoilayp4, zsodp)
#endif
        rcode = nf_redef (ncid_a)

c          Min/Max dimension

        rcode = nf_def_dim (ncid_a, 'minmax', n2, minmaxid)
        rcode = nf_def_var (ncid_a, 'minmax', nf_int4, n1, minmaxid,
     *                      varid)
        rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                           n7, 'min/max')
        rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                           n6, '1 or 2')
        rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                           n2, 'i2') 
        rcode = nf_enddef (ncid_a)
        rcode = nf_put_vara_int (ncid_a, varid, n1, n2, minmax)
        rcode = nf_redef (ncid_a)

c          Time dimension (days since base date)

        rcode = nf_def_dim (ncid_a, 'time', ncunlim4, timid)
        rcode = nf_def_var (ncid_a, 'time', nf_float4, n1, timid, varid)
        rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                           n4, 'time')
        rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                           n4, 'days')
        rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                           n4, 'f9.3') 

c          Define "date" variable (not dimension) 
        ndim    = 1
        idim(1) = timid
        rcode = nf_def_var (ncid_a, 'date', nf_int4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid_a,varid,'long_name',n7, 'date+hr')
        rcode = nf_put_att_text (ncid_a,varid,'units',n10, 'yyyymmddhh')
        rcode = nf_put_att_text (ncid_a,varid,'FORTRAN_format',n3,'i10')

c>>>>>>>>>>
      endif
c>>>>>>>>>>

      rcode = nf_enddef (ncid_a)

c        Add current time (days since base date) to time dimension

      day = dtgcm * 0.5*((nsteph_a-1)+nstep) / 86400.
      rcode = nf_inq_varid (ncid_a, 'time', varid)
      ival4 = nwrit_a + 1
#ifdef cray
      rcode = nf_put_var1_real (ncid_a, varid, ival4, day)
#else
      rcode = nf_put_var1_double (ncid_a, varid, ival4, day)
#endif

c        Add current date (at end of averaging period) to date variable

      rcode = nf_inq_varid (ncid_a, 'date', varid)
      start(1) = nwrit_a + 1
      count(1) = 1
      idate4 = mcdate*100 + nint(mcsec/3600.)               ! yyyymmddhh
      rcode = nf_put_vara_int (ncid_a, varid, start, count, idate4)

      rcode = nf_redef (ncid_a)
 
c        Define variable dimensions and attributes, but only for first
c        encounter of each variable name for this file (if appending, 
c        call is ignored). (numfielnet_a zeroed in openf).

c~~~~~~~~~~~~~~~~~~~~
      do n=1,nushis_a
c~~~~~~~~~~~~~~~~~~~~

        ifdone = 0
        if (numfielnet_a.gt.0) then
          do m=1,numfielnet_a
            if (cusfielnet_a(m).eq.cindex(1,mushis_a(n))) ifdone = 1
          enddo
        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (ifdone.eq.0) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c         increment field number and name list for this file (in comlun)
          numfielnet_a = numfielnet_a + 1
          cusfielnet_a(numfielnet_a) = cindex(1,mushis_a(n))

          idim(1) = lonid
          idim(2) = latid
          if (mindex_a(2,mushis_a(n)).eq.1) then
            ndim    = 3
            idim(3) = timid
          else if (mindex_a(2,mushis_a(n)).eq.2) then
            ndim    = 4
            idim(3) = minmaxid
            idim(4) = timid
          else if (mindex_a(2,mushis_a(n)).eq.nsoilay) then
            ndim    = 4
            idim(3) = zsoilid
            idim(4) = timid
          else if (mindex_a(2,mushis_a(n)).eq.nsoilay+1) then
            ndim    = 4
            idim(3) = zsoilpid
            idim(4) = timid
          else
            write (ioterm2,*) '*** Error in writhlsx(netcdf):',
     *                        ' 3rd dim= ', mindex_a(2,mushis_a(n))
            call endrun (-1)
          endif
          rcode = nf_def_var (ncid_a, cindex(1,mushis_a(n)), nf_float4,
     *                        ndim, idim, varid)

          ival4 = lenchr(cindex(1,mushis_a(n)))
          rcode = nf_put_att_text (ncid_a, varid, 'long_name',
     *                             ival4, cindex(1,mushis_a(n)))

          ival4 = lenchr(cindex(2,mushis_a(n)))
          rcode = nf_put_att_text (ncid_a, varid, 'units',
     *                             ival4, cindex(2,mushis_a(n)))

          rcode = nf_put_att_text (ncid_a, varid, 'FORTRAN_format',
     *                             n5, 'e13.5')
c~~~~~~~~~~~~
        endif
c~~~~~~~~~~~~
      enddo
c~~~~~~~~~~

c        Leave define mode
  
      rcode = nf_enddef (ncid_a)

c        Assign values to variables

      do n=1,nushis_a

        rcode = nf_inq_varid (ncid_a, cindex(1,mushis_a(n)), varid)
        start(1) = 1
        count(1) = nlon
        start(2) = 1
        count(2) = nlat

        do j=1,mindex_a(2,mushis_a(n))
          if (mindex_a(2,mushis_a(n)).eq.1) then
            start(3) = nwrit_a + 1
            count(3) = 1
          else
            start(3) = j
            count(3) = 1
            start(4) = nwrit_a + 1
            count(4) = 1
          endif
          irec = mindex_a(3,mushis_a(n)) + j - 1
          call gethrlsx (irec, workl, iv)

c            For accumulated fields, divide by counter 
c            (all counters zero-ed by hisselsx after return to histpost)

          if (mindfl_a(iv).eq.1) then
             do jj=1,nlat
                do ji=1,nlon
                   workl(ji,jj) = workl(ji,jj) / mindco_a(iv)
                enddo
             enddo
          endif
#ifdef cray
          rcode = nf_put_vara_real (ncid_a, varid, start, count, workl)
#else
          rcode = nf_put_vara_double(ncid_a, varid, start, count, workl)
#endif
        enddo

      enddo

c        Return to define mode

      rcode = nf_redef (ncid_a)
      
      return

#endif
c==========
      endif
c==========
 
      end
c
c-----------------------------------------------------------------------
c
      subroutine histpre (dtgcm, fsave, fnexf, istepy, nstep,
     *                    nrstrt, nstop)
c
c        Sets lsx history flags, called from the
c        start of surfctl each timestep. Sets:
c           fahist_a (whether in accum period)
c           fwhist_a (whether at write timestep)
c           frange       (whether current day contains "_a" accum/write)
c           frange24     (whether end of timestep is at 24:00 GMT)
c        Tests on istepy, the timestep within current 365-day year.
c        Similar to section of ccm subr linems in ccm.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
C
      logical fsave, fnexf
 
      fahist_a = .false.
      fwhist_a = .false.
      frange   = .false.
      frange24 = .false.
 
      do 130 i=1,nhistn_a
        if (nhist_a(i,1).le.nhist_a(i,2)) then
           if (istepy.ge.nhist_a(i,1) .and. istepy.le.nhist_a(i,2))
     *        fahist_a = .true.
        else
           if (istepy.ge.nhist_a(i,1) .or.  istepy.le.nhist_a(i,2))
     *        fahist_a = .true.
        endif
        if (istepy.eq.nhist_a(i,2)) fwhist_a = .true.
        if (histlast_a .and. nstep.eq.nstop) fwhist_a = .true.
  130 continue

#ifdef rcm
c     Since rcm's timestep < lsx's, 1st call's istepy (must be called
c     on rcm's first timestep) from previous run's last istepy. So to
c     avoid spuroius LSX write, don't allow on first call.
      if (nstep.eq.nrstrt) fwhist_a = .false.
#endif
 
      npd = max (nint(86400./dtgcm), 1)
      iday = (istepy-1)/npd + 1
 
      do 132 i=1,nhistn_a
         iday1 = (nhist_a(i,1)-1)/npd + 1
         iday2 = (nhist_a(i,2)-1)/npd + 1
         if (iday1.le.iday2) then
           if (iday.ge.iday1 .and. iday.le.iday2) frange = .true.
         else
           if (iday.ge.iday1 .or.  iday.le.iday2) frange = .true.
         endif
  132 continue
 
      if ( mod(nstep,npd).eq.0 ) frange24 = .true.
      if (fwhist_a) frange24 = .true.

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine histpost (dtgcm,  nsrest, fsave, fnexf,
     *                     nstep,  nrstrt, nstop, istepy,
     *                     mbdate, mcdate, mcsec, runtitle
#ifdef rcm
     *                    ,idatefile
#endif
     *                     )
c
c        Do next-write, save and next-file processing for lsx
c        history files, called at the end of surfctl
c        each timestep. Similar to section of ccm subr linems in ccm.f.
c
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
c
      logical fsave, fnexf
      character*80 runtitle, cname, dname

      cname = msnam_a(ifile)
      if (hisnet) then
c       append '.nc' to name if not already there
        lnl = lenchr(cname)
        if (cname(lnl-2:lnl).ne.'.nc' .and.
     *      cname(lnl-2:lnl).ne.'.NC')  cname(lnl+1:) = '.nc'
      endif
c
c        Initialize nsteph_a to "unset", then set nsteph_a
c        to mark beginning of current accumulation period
c
      if (nstep.eq.nrstrt) nsteph_a = -1
      if (fahist_a.and.nsteph_a.eq.-1) nsteph_a = nstep
c
c        Increment write counter, write header and data, reset rec 
c        number offset, and unset nsteph_a (if any data recs written 
c        this run).
c
      if (fwhist_a .and. nuldat_a.gt.0) then
         call writhlsx (dtgcm,  nsrest, fsave,
     *                  mbdate, mcdate, mcsec, runtitle, nstep)
         call hisselsx(0)
         nsteph_a = -1
#ifdef rcm
         write (ioterm2,'(/a,a/)')
     *    'Writing to current LSX history file ', cname(1:lenchr(cname))
#endif
      endif

c     write (178,*)
c     write (178,'(a,i12)') 'nstep=',nstep
c     write (178,'(a,i12)') 'nstop=',nstop
c     write (178,'(a,l)')   'fsave=',fsave
c     write (178,'(a,l)')   'fnexf=',fnexf
c     write (178,'(a,i12)') 'idatefile=',idatefile
c     ier = flush(178)

c        Save and next-file for lsx history file
c
c>>>>>>>>>>>>>>>>>>>>
      if (fsave) then
c>>>>>>>>>>>>>>>>>>>>
c
c           Save history file to mss, if non-empty and any data
c           records requested this run
c
         if (nuoff_a.gt.0 .and. nuldat_a.gt.0) then
            call closef (ndata_a)
            call disposef (msnam_a(ifile), mspath)
            if (nstep.lt.nstop)
     *        call openf (ndata_a, msnam_a(ifile), mspath, -1, 0)

c----------
#ifndef rcm
c----------
            if (.not.hisnet) write(ioterm2,1450)
            if (     hisnet) write(ioterm2,1451)
            write(ioterm2,1452)
     *         nstep, istepy*dtgcm/86400., nstep*dtgcm/86400.,
     *         mspath(1:lenchr(mspath)), cname(1:lenchr(cname)), nwrit_a
 1450       format(/' ----------------------------------------------'
     *             /' ********** LSX HISTORY FILE OUTPUT ***********')
 1451       format(/' ----------------------------------------------'
     *             /' ****** LSX HISTORY FILE OUTPUT (NetCDF) ******')
 1452       format( ' ----------------------------------------------'
     *             /' timestep                          = ',i10
     *             /' calendar days                     = ',f10.3
     *             /' elapsed days                      = ',f10.3
     *             /' lsx history file                  = ',a,a
     *             /' lsx history writes                = ',i10
     *             /' ----------------------------------------------')
#ifdef mss
            if (mspath.eq.' ') write (ioterm2,1453)
 1453       format( ' *  File NOT saved to MSS (MSPATH is blank!)  *'
     *             /' ----------------------------------------------')
#endif
            write (ioterm2,*)
c-----
#endif
c-----
c
         endif
c
c           If next-file, start new history file
c
c~~~~~~~~~~~~~~~~~~~~~~~
         if (fnexf) then
c~~~~~~~~~~~~~~~~~~~~~~~
            call hisselsx(1)

#if defined (rcm)
            if (nstep.lt.nstop) call closef (ndata_a)

c              Make permament name using date (idatefile, passed from 
c              rcm) and rename just-finished current file to ..<date>.nc

            write (dname,'(a,".",i8.8,".nc")')
     *             msnam_a(ifile)(1:lenchr(msnam_a(ifile))), idatefile
            write (ioterm2,'(4a)') 'Renaming LSX history file ',
     *             cname(1:lenchr(cname)),' to ', dname(1:lenchr(dname))
            ier = ishell ('mv ' // cname(1:lenchr(cname)) //
     *                    ' '   // dname(1:lenchr(dname)) )
            if (ier.ne.0) then
              write (ioterm2,*)
     *           '*** Error (histpost, ishell mv), ier=',ier
              stop
            endif

            if (nstep.lt.nstop) then
              ifile = ifile + 1              ! (for gcm, done in genctl)
              if (ifile.le.nfile) then 
                call openf (ndata_a, msnam_a(ifile), mspath,1,nusiz_a)
              else
                write(ioterm2,*)
     *            '*** Not enough LSX history filenames input to',
     *            ' complete this run. Stop in HISTPOST.'
                call endrun (-1)
              endif
            endif
#else
            if (nstep.lt.nstop) then
               call closef (ndata_a)
               if (ifile+1.le.nfile)
     *           call openf (ndata_a, msnam_a(ifile+1),mspath,1,nusiz_a)
            endif
#endif

c~~~~~~~~~~~~~
         endif
c~~~~~~~~~~~~~
c
c>>>>>>>>>>
      endif
c>>>>>>>>>>
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine histwrit
c
c        For lsx history fields:
c        Assembles lsx program history fields if necessary, and calls
c        puthilsx to write to lsx history files if requested.
c
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
#include <commapl.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#include <com1d.h>
      logical testhisl
#include <comsat.h>
c
c
c        lsx land-ocean mask (1=land, 2=icesheet, 3=ocean)
c
      if (testhisl(1)) then
        do 10 jj=1,nlat
          do 12 ji=1,nlon
            workl(ji,jj) = lmaskxy(ji,jj)
   12     continue
   10   continue
        call puthilsx (workl, 1)
      endif
c
c        Areal mean surface temperature, and its daily min/max (K)
c
      if (testhisl(2) .or. testhisl(47)) then
        do 20 jj=1,nlat
          do 22 ji=1,nlon
            x = (1. - fixy(ji,jj))*tsoixy(ji,jj,1)
     *              + fixy(ji,jj) *tsnoxy(ji,jj,1)
            y = (1. - fi_icexy(ji,jj))*ticexy(ji,jj,1)
     *              + fi_icexy(ji,jj) *tsno_icexy(ji,jj,1)
            z = tocnxy(ji,jj)
            workl(ji,jj)=(1.-fwaterxy(ji,jj))*                        x
     *                 +     fwaterxy(ji,jj) * (      ficexy(ji,jj) * y
     *                                          + (1.-ficexy(ji,jj))* z)
   22     continue
   20   continue
        call puthilsx (workl, 2)
        call putrange (workl, 47)
      endif
c
c        Soil/icesheet layer temperatures (K)
c
      call puthilsx (tsoixy(1,1,1), 3)
      call puthilsx (tsoixy(1,1,2), 4)
      call puthilsx (tsoixy(1,1,3), 5)
      call puthilsx (tsoixy(1,1,4), 6)
      call puthilsx (tsoixy(1,1,5), 7)
      call puthilsx (tsoixy(1,1,6), 8)
c
c        Soil layer liquid content (0 to 1)
c
      call puthilsx (wsoixy(1,1,1), 9)
      call puthilsx (wsoixy(1,1,2), 10)
      call puthilsx (wsoixy(1,1,3), 11)
      call puthilsx (wsoixy(1,1,4), 12)
      call puthilsx (wsoixy(1,1,5), 13)
      call puthilsx (wsoixy(1,1,6), 14)
c
c        Soil layer ice content (0 to 1)
c
      call puthilsx (wisoixy(1,1,1), 15)
      call puthilsx (wisoixy(1,1,2), 16)
      call puthilsx (wisoixy(1,1,3), 17)
      call puthilsx (wisoixy(1,1,4), 18)
      call puthilsx (wisoixy(1,1,5), 19)
      call puthilsx (wisoixy(1,1,6), 20)
c
c        Snow thickness (m)
c
      if (testhisl(21)) then
        call zero (workl, nlon*nlat)
        call zero (workm, nlon*nlat)
        do 210 k=1,nsnolay
          do 212 jj=1,nlat
            do 214 ji=1,nlon
              workl(ji,jj) = workl(ji,jj) + hsnoxy(ji,jj,k)
              workm(ji,jj) = workm(ji,jj) + hsno_icexy(ji,jj,k)
  214       continue
  212     continue
  210   continue
        do 216 jj=1,nlat
          do 218 ji=1,nlon
            workl(ji,jj) = (1.-fwaterxy(ji,jj))*workl(ji,jj)
     *                        +fwaterxy(ji,jj) *workm(ji,jj)
c           For icesheet initiation points (don't use water fraction):
c           workl(ji,jj) = cvmgt (workl(ji,jj), workm(ji,jj),
c    *                            fwaterxy(ji,jj).lt.0.999)
  218     continue
  216   continue
        call puthilsx (workl, 21)
      endif
c
c        Snow fraction (0 to 1)
c
      if (testhisl(22)) then
        do 2290 jj=1,nlat
          do 2292 ji=1,nlon
            workl(ji,jj) = (1.-fwaterxy(ji,jj))*fixy(ji,jj)
     *                        +fwaterxy(ji,jj) *fi_icexy(ji,jj)
c           For icesheet initiation points (don't use water fraction):
c           workl(ji,jj) = cvmgt (fixy(ji,jj), fi_icexy(ji,jj),
c    *                            fwaterxy(ji,jj).lt.0.999)
 2292     continue
 2290   continue
        call puthilsx (workl, 22)
      endif
c
c        Seaice thickness (m)
c
      if (testhisl(23)) then
        call zero (workl, nlon*nlat)
        do 230 k=1,nicelay
          do 232 jj=1,nlat
            do 234 ji=1,nlon
              workl(ji,jj) = workl(ji,jj) + hicexy(ji,jj,k)
  234       continue
  232     continue
  230   continue
        call puthilsx (workl, 23)
      endif
c
c        Seaice fraction (0 to 1)
c
      call puthilsx (ficexy, 24)
c
c        Lowest-level agcm air temp, and its daily min/max (K)
c
      if (testhisl(25) .or. testhisl(48)) then
        call scopy (nlon*nlat, taxy, 1, workl, 1)
        call puthilsx (workl, 25)
        call putrange (workl, 48)
      endif
c
c        Upper-leaf temperature (K)
c
      call puthilsx (tuxy, 26)
c
c        Upper-stem temperature (K)
c
      call puthilsx (tsxy, 27)
c
c        Lower-veg temperature (K)
c
      call puthilsx (tlxy, 28)
c
c        Net surface runoff (m/s)
c
      if (testhisl(29)) then
        do 290 jj=1,nlat
          do 292 ji=1,nlon
            workl(ji,jj) = grunofxy(ji,jj) / rhow
  292     continue
  290   continue
        call puthilsx (workl, 29)
      endif
c
c        Net subsurface drainage (m/s)
c
      if (testhisl(30)) then
        do 300 jj=1,nlat
          do 302 ji=1,nlon
            workl(ji,jj) = gdrainxy(ji,jj) / rhow
  302     continue
  300   continue
        call puthilsx (workl, 30)
      endif
c
c        Net transpiration rate (m/s)
c
      if (testhisl(31)) then
        do 310 jj=1,nlat
          do 312 ji=1,nlon
            workl(ji,jj) = (gtransuxy(ji,jj)+gtranslxy(ji,jj)) / rhow
  312     continue
  310   continue
        call puthilsx (workl, 31)
      endif
c
c        Net transpiration rate (m/s)
c
      if (testhisl(32)) then
        do jj=1,nlat
          do ji=1,nlon
            workl(ji,jj) = ginfilxy(ji,jj) / rhow
          enddo
        enddo
        call puthilsx (workl, 32)
      endif
c
c        Roughness length (m)
c
      call puthilsx (z0lenxy, 33)
c
c        Dorman-Sellers/SDVM/IBIS biome/veg.-type index
c        (BIOME-model biome index is ibiomexy, see below)
c
      if (testhisl(34)) then
        do 340 jj=1,nlat
          do 342 ji=1,nlon
            if ( ivegxy(ji,jj).eq.nved .and. 
     *           (vegtype.eq.0.or.vegtype.eq.4.or.vegtype.eq.5.or.
     *            vegtype.eq.6.or.vegtype.eq.7) ) then
              workl(ji,jj) = 0.
            else
              workl(ji,jj) = ivegxy(ji,jj)
            endif
  342     continue
  340   continue
        call puthilsx (workl, 34)
      endif
c
c        Wind stress (N m-2)
c
      call puthilsx (atauxxy, 35)
      call puthilsx (atauyxy, 36)
c
c        Net intercepted water (m)
c
      if (testhisl(37)) then
        do 370 jj=1,nlat
          do 372 ji=1,nlon
            workl(ji,jj) = ( (wliquxy(ji,jj)+wsnouxy(ji,jj))
     *                          * laixy(ji,jj,2)*fuxy(ji,jj)
     *                     + (wliqsxy(ji,jj)+wsnosxy(ji,jj))
     *                          * saixy(ji,jj,2)*fuxy(ji,jj)
     *                     + (wliqlxy(ji,jj)+wsnolxy(ji,jj))
     *                          * (laixy(ji,jj,1)+saixy(ji,jj,1))
     *                          * flxy(ji,jj)*(1.-fixy(ji,jj))
     *                     ) / rhow
  372     continue
  370   continue
        call puthilsx (workl, 37)
      endif
c
c        #38 and #39 are for dynamic sea-ice u,v
c
c        Net upward water vapor flux to agcm (m/s)
c
      if (testhisl(40)) then
        do 400 jj=1,nlat
          do 402 ji=1,nlon
            workl(ji,jj) = -avapaxy(ji,jj) / rhow
  402     continue
  400   continue
        call puthilsx (workl, 40)
      endif
c
c        Oceanic heat transport convergence
c
      call puthilsx (qfluxxy, 41)
c
c        Wisconsin vegetation category index
c
      if (testhisl(42)) then
        do 420 jj=1,nlat
          do 422 ji=1,nlon
            workl(ji,jj) = ivegxy(ji,jj)
  422     continue
  420   continue
        call puthilsx (workl, 42)
      endif
c
c        Seaice brine reservoir as a fraction of heat needed to melt
c        entire column
c
      if (testhisl(43)) then
        call zero (workl, nlon*nlat)
        do 430 k=1,nicelay
          do 432 jj=1,nlat
            do 434 ji=1,nlon
              workl(ji,jj) = workl(ji,jj)
     *        + rhoi*hicexy(ji,jj,k)*(cice*(tmelt-ticexy(ji,jj,k))+hfus)
  434       continue
  432     continue
  430   continue
c
        do 436 jj=1,nlat
          do 438 ji=1,nlon
            if (ficexy(ji,jj).gt.0.) then
              workl(ji,jj) = brinexy(ji,jj) / workl(ji,jj)
            else
              workl(ji,jj) = 0.
            endif
  438     continue
  436   continue
        call puthilsx (workl, 43)
      endif
c
c        #44 and #45 are for dynamic sea-ice d(fice)/dt, d(hice)/dt
c
c        Ocean mixed-layer depth
c
      call puthilsx (hocnxy, 46)
c
c        #47 and #48 are for daily min/max of TS(#2) and TA(#25)
c
c        2-m  screen air temp and its daily min/max (#49,50)
c
      call puthilsx (ts2xy, 49)
      call putrange (ts2xy, 50)
c
c        Relative humidity at lowest agcm level
c
      if (testhisl(51)) then
        do 510 jj=1,nlat
          do 512 ji=1,nlon
            zes = esat (taxy(ji,jj))
            zp = psurfxy(ji,jj)*siga
            zmr = qaxy(ji,jj) / (1.-qaxy(ji,jj))
            workl(ji,jj) = min ((zmr/(.622+zmr))*(zp/zes), 1.)
  512     continue
  510   continue
        call puthilsx (workl, 51)
      endif
c
c        Ocean (mixed-layer or ogcm surface) temperature (K)
c
      if (testhisl(52)) then
        do jj=1,nlat
          do ji=1,nlon
            if (lmaskxy(ji,jj).eq.3) then
              workl(ji,jj) = tocnxy(ji,jj)
            else
              workl(ji,jj) = 0.
            endif
          enddo
        enddo
        call puthilsx (workl, 52)
      endif
c
c        Net evaporation from intercepted h2o (m/s)
c
      if (testhisl(53)) then
        do 530 jj=1,nlat
          do 532 ji=1,nlon
            workl(ji,jj) = ginvapxy(ji,jj) / rhow
  532     continue
  530   continue
        call puthilsx (workl, 53)
      endif
c
c        Net evaporation from soil/ice/snow/ocean (m/s)
c
      if (testhisl(54)) then
        do 540 jj=1,nlat
          do 542 ji=1,nlon
            workl(ji,jj) = gsuvapxy(ji,jj) / rhow
  542     continue
  540   continue
        call puthilsx (workl, 54)
      endif
c
c        Fog formation over land (m/s)
c
      if (testhisl(55)) then
        do 550 jj=1,nlat
          do 552 ji=1,nlon
            workl(ji,jj) = (fog12xy(ji,jj)+fog34xy(ji,jj)) / rhow
  552     continue
  550   continue
        call puthilsx (workl, 55)
      endif
c
c        Soil interfaces downward flow (multi-level field) (kg/m2/s)
c
      call puthilsx (wfloxy, 58)
c
c        Sea-ice layer temperatures (K)
c
      call puthilsx (ticexy(1,1,1), 59)
      call puthilsx (ticexy(1,1,2), 60)
      call puthilsx (ticexy(1,1,3), 61)
c
c        Prescribed open-water fraction (0-1)
c
      call puthilsx (fwaterxy, 62)
c
c        Net one-sided LAI, upper story, * FUXY
c
      if (testhisl(63)) then
        do 630 jj=1,nlat
          do 632 ji=1,nlon
            workl(ji,jj) = 0.5 * laixy(ji,jj,2)*fuxy(ji,jj)
  632     continue
  630   continue
        call puthilsx (workl, 63)
      endif
c
c        Net one-sided LAI, upper story, * FLXY
c
      if (testhisl(64)) then
        do 640 jj=1,nlat
          do 642 ji=1,nlon
            workl(ji,jj) = 0.5 * laixy(ji,jj,1)*flxy(ji,jj)
  642     continue
  640   continue
        call puthilsx (workl, 64)
      endif
c
c        Net one-sided LAI, summed over stories, * FUXY and FLXY
c
      if (testhisl(65)) then
        do 650 jj=1,nlat
          do 652 ji=1,nlon
            workl(ji,jj) = 0.5 * (  laixy(ji,jj,1)*flxy(ji,jj)
     *                            + laixy(ji,jj,2)*fuxy(ji,jj))
  652     continue
  650   continue
        call puthilsx (workl, 65)
      endif
c
c        Net one-sided SAI, upper story, * FUXY
c
      if (testhisl(66)) then
        do 660 jj=1,nlat
          do 662 ji=1,nlon
            workl(ji,jj) = 0.5 * saixy(ji,jj,2)*fuxy(ji,jj)
  662     continue
  660   continue
        call puthilsx (workl, 66)
      endif
c
c        Net one-sided SAI, lower story, * FLXY
c
      if (testhisl(67)) then
        do 670 jj=1,nlat
          do 672 ji=1,nlon
            workl(ji,jj) = 0.5 * saixy(ji,jj,1)*flxy(ji,jj)
  672     continue
  670   continue
        call puthilsx (workl, 67)
      endif
c
c        Net one-sided SAI, summed over stories, * FUXY and FLXY
c
      if (testhisl(68)) then
        do 680 jj=1,nlat
          do 682 ji=1,nlon
            workl(ji,jj) = 0.5 * (  saixy(ji,jj,1)*flxy(ji,jj)
     *                            + saixy(ji,jj,2)*fuxy(ji,jj))
  682     continue
  680   continue
        call puthilsx (workl, 68)
      endif
c
c----------
#ifdef ibis
c----------
c
c        For IBIS photosynthesis/stomatal conductance history fields:
c
c        IBIS net primary productivity per upper-story leaf area,
c        in kg_carbon/m2/yr. Convert from mole_co2/m2/s to
c        kg_carbon/m2/yr, so * 12 for atomic weight of carbon.
c
      if (testhisl(69)) then
        do 690 jj=1,nlat
          do 692 ji=1,nlon
            workl(ji,jj) = (  (1.-confracxy(ji,jj))*anudxy(ji,jj)
     *                       +    confracxy(ji,jj) *anucxy(ji,jj) )
     *                     * 12. * 1.e-3 * .864e5 * 365.
  692     continue
  690   continue
        call puthilsx (workl, 69)
      endif
c
c        IBIS net primary productivity per lower-story leaf area,
c        in kg_carbon/m2/yr. Convert from mole_co2/m2/s to
c        kg_carbon/m2/yr, so * 12 for atomic weight of carbon.
c
      if (testhisl(70)) then
        do 700 jj=1,nlat
          do 702 ji=1,nlon
            workl(ji,jj) = (  (1.-c4fracxy(ji,jj))*anl3xy(ji,jj)
     *                       +    c4fracxy(ji,jj) *anl4xy(ji,jj) )
     *                     * 12. * 1.e-3 * .864e5 * 365.
  702     continue
  700   continue
        call puthilsx (workl, 70)
      endif
c
c        IBIS gross primary productivity per upper-story leaf area,
c        in kg_carbon/m2/yr. Convert from mole_co2/m2/s to
c        kg_carbon/m2/yr, so * 12 for atomic weight of carbon.
c
      if (testhisl(71)) then
        do 710 jj=1,nlat
          do 712 ji=1,nlon
            workl(ji,jj) = (  (1.-confracxy(ji,jj))*agudxy(ji,jj)
     *                       +    confracxy(ji,jj) *agucxy(ji,jj) )
     *                     * 12. * 1.e-3 * .864e5 * 365.
  712     continue
  710   continue
        call puthilsx (workl, 71)
      endif
c
c        IBIS gross primary productivity per lower-story leaf area,
c        in kg_carbon/m2/yr. Convert from mole_co2/m2/s to
c        kg_carbon/m2/yr, so * 12 for atomic weight of carbon.
c
      if (testhisl(72)) then
        do 720 jj=1,nlat
          do 722 ji=1,nlon
            workl(ji,jj) = (  (1.-c4fracxy(ji,jj))*agl3xy(ji,jj)
     *                       +    c4fracxy(ji,jj) *agl4xy(ji,jj) )
     *                     * 12. * 1.e-3 * .864e5 * 365.
  722     continue
  720   continue
        call puthilsx (workl, 72)
      endif
c
c        IBIS gross primary productivity per overall surface area,
c        in kg_carbon/m2/yr. Convert from mole_co2/m2/s to
c        kg_carbon/m2/yr, so * 12 for atomic weight of carbon.
c
      if (testhisl(73)) then
        do 730 jj=1,nlat
          do 732 ji=1,nlon
            zzu =  (1.-confracxy(ji,jj))*agudxy(ji,jj)
     *            +    confracxy(ji,jj) *agucxy(ji,jj)
            zzl =  (1.-c4fracxy(ji,jj) )*agl3xy(ji,jj)
     *            +    c4fracxy(ji,jj)  *agl4xy(ji,jj)
            workl(ji,jj) =
     *       (1.-fwaterxy(ji,jj))
     *      *(                   fuxy(ji,jj)*laixy(ji,jj,2)*zzu
     *        + (1.-fixy(ji,jj))*flxy(ji,jj)*laixy(ji,jj,1)*zzl
     *       )
     *      * 12. * 1.e-3 * .864e5 * 365.
  732     continue
  730   continue
        call puthilsx (workl, 73)
      endif
c
c        IBIS stomatal conductance for upper-story leaves (m/s).
c        Conversion factor .029 = (rh2o/rair)*18*1.e-3 from
c        mole_h2o/m2/s per delta(svp/P) to Kg_h2o/m2/s per delta(q),
c        where 18 = molecular weight of h2o. Then divide by rhoa
c        to get to m/s = Kg_h2o/m2/s per delta(Kg_h2o/m3).
c
      if (testhisl(74)) then
        do 740 jj=1,nlat
          do 742 ji=1,nlon
            zrhoa = (psurfxy(ji,jj)*siga)/(rair*taxy(ji,jj))
            workl(ji,jj) = (  (1.-confracxy(ji,jj))*gsudxy(ji,jj)
     *                       +    confracxy(ji,jj) *gsucxy(ji,jj)
     *                     ) * .029 / zrhoa
  742     continue
  740   continue
        call puthilsx (workl, 74)
      endif
c
c        IBIS stomatal conductance for lower-story leaves (m/s)
c
      if (testhisl(75)) then
        do 750 jj=1,nlat
          do 752 ji=1,nlon
            zrhoa = (psurfxy(ji,jj)*siga)/(rair*taxy(ji,jj))
            workl(ji,jj) = (  (1.-c4fracxy(ji,jj))*gsl3xy(ji,jj)
     *                       +    c4fracxy(ji,jj) *gsl4xy(ji,jj)
     *                     ) * .029 / zrhoa
  752     continue
  750   continue
        call puthilsx (workl, 75)
      endif
 
c-----
#endif
c-----
 
c        soil moisture stress factor, upper and lower veg (see setsoi
c        and drystress...weighted by root fractions, also depends on
c        soil temperature for BATS stomatal resistance) (0-1)
 
      call puthilsx (uptotuxy, 76)
      call puthilsx (uptotlxy, 77)
 
c        Rainfall rate (m/s)
 
      if (testhisl(78)) then
        do 780 jj=1,nlat
          do 782 ji=1,nlon
            workl(ji,jj) = rainaxy(ji,jj) / rhow
  782     continue
  780   continue
        call puthilsx (workl, 78)
      endif
 
c        Snowfall rate (m/s)
 
      if (testhisl(79)) then
        do 790 jj=1,nlat
          do 792 ji=1,nlon
            workl(ji,jj) = snowaxy(ji,jj) / rhow
  792     continue
  790   continue
        call puthilsx (workl, 79)
      endif
 
c        Precipitation rate (m/s)
 
      if (testhisl(80)) then
        do 800 jj=1,nlat
          do 802 ji=1,nlon
            workl(ji,jj) = (rainaxy(ji,jj) + snowaxy(ji,jj)) / rhow
  802     continue
  800   continue
        call puthilsx (workl, 80)
      endif
 
c        Ocean (ogcm) surface salinity (fraction)
 
      if (testhisl(81)) then
        do jj=1,nlat
          do ji=1,nlon
            if (lmaskxy(ji,jj).eq.3) then
              workl(ji,jj) = socnxy(ji,jj)
            else
              workl(ji,jj) = 0.
            endif
          enddo
        enddo
        call puthilsx (workl, 81)
      endif
 
c        h2o and pore space in total soil column (m)
 
      if (testhisl(82) .or. testhisl(83)) then
        call zero (workl, nlon*nlat)
        call zero (workm, nlon*nlat)
        do 820 k=1,nsoilay
          do 822 jj=1,nlat
            do 824 ji=1,nlon
              zh = hsod(k)*porosxy(ji,jj,k)
              workl(ji,jj) = workl(ji,jj)
     *                     + zh*( (1.-wisoixy(ji,jj,k))*wsoixy(ji,jj,k)
     *                             + wisoixy(ji,jj,k) )
              workm(ji,jj) = workm(ji,jj) + zh
  824       continue
  822     continue
  820   continue
        call puthilsx (workl, 82)
        call puthilsx (workm, 83)
      endif
 
c        h2o and pore space in top 30 cm of soil column (m)
 
      if (testhisl(84) .or. testhisl(85)) then
        call zero (workl, nlon*nlat)
        call zero (workm, nlon*nlat)
        zcut = 0.30
        zha = 0.
        do 840 k=1,nsoilay
          zhb = zha + hsod(k)
          if (zha.lt.zcut) then
            zhsod = min(zhb,zcut) - zha
            do 842 jj=1,nlat
              do 844 ji=1,nlon
                zh = zhsod*porosxy(ji,jj,k)
                workl(ji,jj) = workl(ji,jj)
     *                       + zh*((1.-wisoixy(ji,jj,k))*wsoixy(ji,jj,k)
     *                               + wisoixy(ji,jj,k))
                workm(ji,jj) = workm(ji,jj) + zh
  844         continue
  842       continue
          endif
          zha = zhb
  840   continue
        call puthilsx (workl, 84)
        call puthilsx (workm, 85)
      endif

 
c          net heat flux into ocean (below sea ice) (W/m2)
 
      call puthilsx (ftocnxy, 86)
 
c          net P-E into ocean (below sea ice) (m/s)
 
      if (testhisl(87)) then
        do 870 jj=1,nlat
          do 872 ji=1,nlon
            workl(ji,jj) = fwocnxy(ji,jj) / rhow
  872     continue
  870   continue
        call puthilsx (workl, 87)
      endif   

c        Soil layer porosity (0 to 1)
 
      call puthilsx (porosxy(1,1,1), 88)
      call puthilsx (porosxy(1,1,2), 89)
      call puthilsx (porosxy(1,1,3), 90)
      call puthilsx (porosxy(1,1,4), 91)
      call puthilsx (porosxy(1,1,5), 92)
      call puthilsx (porosxy(1,1,6), 93)
 
c        Soil/icesheet layer temperatures (multi-level field) (K)
 
      call puthilsx (tsoixy, 94)
 
c        Soil layer liquid moisture contents (multi-level field) (0-1)
 
      call puthilsx (wsoixy, 95)
 
c        Soil layer ice moisture contents (multi-level field) (0-1)
 
      call puthilsx (wisoixy, 96)
 
c        Soil layer porosities (multi-level field) (0-1)
 
      call puthilsx (porosxy, 97)

c>>>>>>>>>>>>>>>>>>>>>
#if defined  (isotrac)
c>>>>>>>>>>>>>>>>>>>>>
c        Sea-ice surface reservoir

      call puthilsx (hi_i, 98)
      do 876 jj=1,nlat
        do 878 ji=1,nlon
          workl(ji,jj) = hi_i(ji,jj)*ri_o(ji,jj)
          workm(ji,jj) = hi_i(ji,jj)*ri_h(ji,jj)
  878   continue
  876 continue
      call puthilsx (workl, 99)
      call puthilsx (workm, 100)

c        Top soil reservoir

      call puthilsx (hct_i, 101)
      do 880 jj=1,nlat
        do 882 ji=1,nlon
          workl(ji,jj) = hct_i(ji,jj)*rct_o(ji,jj)
          workm(ji,jj) = hct_i(ji,jj)*rct_h(ji,jj)
  882   continue
  880 continue
      call puthilsx (workl, 102)
      call puthilsx (workm, 103)

c        Bottom soil reservoir

      call puthilsx (hcb_i, 104)
      do 884 jj=1,nlat
        do 886 ji=1,nlon
          workl(ji,jj) = hcb_i(ji,jj)*rcb_o(ji,jj)
          workm(ji,jj) = hcb_i(ji,jj)*rcb_h(ji,jj)
  886   continue
  884 continue
      call puthilsx (workl, 105)
      call puthilsx (workm, 106)

c        Lake reservoir and ocean/seas

      call puthilsx (hl_i, 107)
      do 887 jj=1,nlat
        do 888 ji=1,nlon
          workl(ji,jj) = hl_i(ji,jj)*ro_o(ji,jj)
          workm(ji,jj) = hl_i(ji,jj)*ro_h(ji,jj)
  888   continue
  887 continue
      call puthilsx (workl, 108)
      call puthilsx (workm, 109)

c        Ocean surface mask for isotopes

      do jj=1,nlat
        do ji=1,nlon
          workl(ji,jj) = mocean_i(ji,jj)
        enddo
      enddo
      call puthilsx (workl, 110)
c>>>>>
#endif
c>>>>>

c        Ocean (ogcm) surface eastward and northward velocities (m/s).
c        Shift lmaskxy since mom2,4's u-grid is at NE corner of h grid.

      if (testhisl(111) .or. testhisl(112)) then
        do jj=1,nlat
          jjp = min (jj,nlat)
          do ji=1,nlon
            jip = mod (ji,nlon) + 1
c           if (lmaskxy(ji,jj ).eq.3 .and. lmaskxy(jip,jj ).eq.3 .and.
c    *          lmaskxy(ji,jjp).eq.3 .and. lmaskxy(jip,jjp).eq.3 ) then
            if (lmaskxy(ji,jj ).eq.3 .or.  lmaskxy(jip,jj ).eq.3 .or. 
     *          lmaskxy(ji,jjp).eq.3 .or.  lmaskxy(jip,jjp).eq.3 ) then
              workl(ji,jj) = uocnxy(ji,jj)
              workm(ji,jj) = vocnxy(ji,jj)
            else
              workl(ji,jj) = 0.
              workm(ji,jj) = 0.
            endif
          enddo
        enddo
        call puthilsx (workl, 111)
        call puthilsx (workm, 112)
      endif

c        Rainfall rate (m/s)
 
      if (testhisl(121)) then
        do 1210 jj=1,nlat
          do 1212 ji=1,nlon
            workl(ji,jj) = rainaxy(ji,jj) / rhow
 1212     continue
 1210   continue
        call puthilsx (workl, 121)
      endif
 
c        Snowfall rate (m/s)
 
      if (testhisl(122)) then
        do 1220 jj=1,nlat
          do 1222 ji=1,nlon
            workl(ji,jj) = snowaxy(ji,jj) / rhow
 1222     continue
 1220   continue
        call puthilsx (workl, 122)
      endif
c
c        BIOME-model biome index
c
      if (testhisl(123)) then
        do 1230 jj=1,nlat
          do 1232 ji=1,nlon
            workl(ji,jj) = ibiomexy(ji,jj)
 1232     continue
 1230   continue
        call puthilsx (workl, 123)
      endif

c        Downward liquid h2o flow across soil interfaces (kg/m2/s)
 
      call puthilsx (wfloxy(1,1,1), 124)
      call puthilsx (wfloxy(1,1,2), 125)
      call puthilsx (wfloxy(1,1,3), 126)
      call puthilsx (wfloxy(1,1,4), 127)
      call puthilsx (wfloxy(1,1,5), 128)
      call puthilsx (wfloxy(1,1,6), 129)
      call puthilsx (wfloxy(1,1,7), 130)

c>>>>>>>>>>>>>>>>>
#if defined  (rcm)
c>>>>>>>>>>>>>>>>>
c        Longitudes and latitudes on cross grid

      call puthilsx (alon2d, 131)
      call puthilsx (alat2d, 132)

c        10-m air temp, specific humidity, u and v wind velocs

      call puthilsx (ts10xy, 133)
      call puthilsx (qs10xy, 134)
      do jj=1,nlat
        do ji=1,nlon
          zua = max ( sqrt (uxxy(ji,jj)**2 + uyxy(ji,jj)**2), 0.1 )
          workl(ji,jj) = us10xy(ji,jj) * uxxy(ji,jj)/zua
          workm(ji,jj) = us10xy(ji,jj) * uyxy(ji,jj)/zua
        enddo
      enddo
      call puthilsx (workl, 135)
      call puthilsx (workm, 136)

c       Surface pressure

      call puthilsx (psurfxy, 137)

c       Solar and ir radiation, downward above sfc 

      call puthilsx (soladxy(1,1,1), 138)  ! solar, direct, vis
      call puthilsx (soladxy(1,1,2), 139)  ! solar, direct, nir
      call puthilsx (solaixy(1,1,1), 140)  ! solar, direct, vis
      call puthilsx (solaixy(1,1,2), 141)  ! solar, direct, nir
      call puthilsx (coszenxy,       142)  ! cosine zenith angle
      call puthilsx (firaxy,         143)  ! downward ir

c        Total cloudiness (0-1)  

      call puthilsx (totcldxy, 144)

c        Soil h2o

      do jj=1,nlat
        do ji=1,nlon
          totpor = 0.
          toth2o = 0.
          totpor3 = 0.
          toth2o3 = 0.
          do k=1,nsoilay
            totpor = totpor + hsod(k)*porosxy(ji,jj,k)
            toth2o = toth2o + hsod(k)*porosxy(ji,jj,k)
     *      * (wsoixy(ji,jj,k)*(1.-wisoixy(ji,jj,k))+wisoixy(ji,jj,k))
            if (k.le.3) then
              totpor3 = totpor3 + hsod(k)*porosxy(ji,jj,k)
              toth2o3 = toth2o3 + hsod(k)*porosxy(ji,jj,k)
     *        *(wsoixy(ji,jj,k)*(1.-wisoixy(ji,jj,k))+wisoixy(ji,jj,k))
            endif
          enddo
          workl(ji,jj) = toth2o  / max (totpor, 1.e-6)
          workm(ji,jj) = toth2o3 / max (totpor3, 1.e-6)
        enddo
      enddo
      call puthilsx (workl, 145)
      call puthilsx (workm, 146)

c        Soil temperature

      do jj=1,nlat
        do ji=1,nlon
          toth = 0.
          tott = 0.
          toth3 = 0.
          tott3 = 0.
          do k=1,nsoilay
            toth = toth + hsod(k)
            tott = tott + hsod(k)*tsoixy(ji,jj,k)
            if (k.le.3) then
              toth3 = toth + hsod(k)
              tott3 = tott + hsod(k)*tsoixy(ji,jj,k)
            endif
          enddo
          workl(ji,jj) = tott  / max (toth,  1.e-6)
          workm(ji,jj) = tott3 / max (toth3, 1.e-6)
        enddo
      enddo
      call puthilsx (workl, 147)
      call puthilsx (workm, 148)
       
c>>>>>
#endif
c>>>>>

c       Uncorrected and corrected (delevtype) topography

      call puthilsx (topnodxy, 149)
      do jj=1,nlat
        do ji=1,nlon
          workl(ji,jj) = topnodxy(ji,jj) + delevxy(ji,jj)
        enddo
      enddo
      call puthilsx (workl, 150)
      call puthilsx (delevxy, 151)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine evewrite (dtgcm, nstep, nrstrt, fsave,
     *                     jyear, jmnth, jday)
 
c        Accumulate daily-mean sfc-air temp and precip for interactive
c        eve/leaf. Do nothing unless vegtype=3. Daily-mean
c        global arrays [T,P,R]EVE_D are in comeve, on restart file,
c        passed by vegdat_eve. Also accumulate monthly mean T,P and
c        rel.hum. in global arrays [T,P,R]EVE_M (comeve), written here
c        to file eve_*.cur files at end of each month, and then 
c        used to update eve_*.exp files at end of year (using running
c        mean filter), to be read by EVE (eco) program.
 
c        All arguments are supplied:
c        DTGCM    = length of time to run lsx this call (=agcm timestep)
c        NSTEP    = agcm timestep counter
c        NRSTRT   = first agcm timestep # for this run
c        FSAVE    = .TRUE. if Genesis restart file-save this timestep
c                   (not used)
c        JYEAR    = current year number (2001 AD,etc, or -21000 BP,etc)
c        JMNTH    = current month number (1-12)
c        JDAY     = current day number (1-31)
 
#include <compar.h>
#include <comeve.h>
#include <comgrd.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#include <com1d.h>
      logical fsave, curexist, expexist

      character*3 cyr
      dimension workcur(nlon,nlat,12), workexp(nlon,nlat,12)
      logical if_bod, if_eod, if_bom, if_eom, if_boy, if_eoy
      dimension ndaypm(12)
      save ndaypm, kount_d, kount_m
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
 
#include <comsat.h>
 
c        Only do for interactive Eve/Leaf
 
      if (vegtype.ne.3) return
 
c       Set beginning/end_of_day/month/year flags
 
      if_bod = .false.
      if_eod = .false.
      if_bom = .false.
      if_eom = .false.
      if_boy = .false.
      if_eoy = .false.
 
      n = nint(86400./dtgcm)
      if (nstep.eq.nrstrt .or. mod(nstep,n).eq.1) then
        if_bod = .true.
      else if ( mod(nstep,n).eq.0 ) then
        if_eod = .true.
      endif
 
      if (if_bod) then
        if (jday.eq.1) if_bom = .true.
      else if (if_eod) then
        if (jday.eq.ndaypm(jmnth)) if_eom = .true.
      endif
 
      if (if_bod .and. jday.eq. 1 .and. jmnth.eq. 1) if_boy = .true.
      if (if_eod .and. jday.eq.31 .and. jmnth.eq.12) if_eoy = .true.
 
c     write(*,8880)nstep, if_bod, if_eod, if_bom, if_eom, if_boy, if_eoy
c8880 format('evewrite: nstep=',i6,'  if[bod,eod,bom,eom,boy,eoy]=',6l3)
 
c     if (nstep.eq.nrstrt) then
c       write (*,8881) nstep, ((nint(teve_d(i,j)),i=1,nlon),j=1,nlat)
c8881   format(/'evewrite: nstep=', i6,'  FIRST TEVE_D :' / (180i3))
c       write (*,8882) nstep, ((nint(peve_d(i,j)),i=1,nlon),j=1,nlat)
c8882   format(/'evewrite: nstep=', i6,'  FIRST PEVE_D :' / (180i3))
c     endif
 
c       If beginning of day (or first timestep), zero *eve_d arrays
c       and counters
 
      if (if_bod .or. nstep.eq.nrstrt) then
        call zero (teve_d, nlon*nlat)
        call zero (peve_d, nlon*nlat)
        call zero (reve_d, nlon*nlat)
        kount_d = 0
      endif
 
c       If beginning of month (or first timestep), zero *eve_m arrays
c       and counters
 
      if (if_bom .or. nstep.eq.nrstrt) then
        call zero (teve_m(1,1), nlon*nlat)
        call zero (peve_m(1,1), nlon*nlat)
        call zero (reve_m(1,1), nlon*nlat)
        kount_m = 0
      endif
 
c        Accumulate current values in *eve_d, *eve_m  arrays. First calc
c        relative humidity at lowest agcm level (as in histwrit).
 
      do 240 jj=1,nlat
        do 242 ji=1,nlon
          zes = esat (taxy(ji,jj))
          zp = psurfxy(ji,jj)*siga
          zmr = qaxy(ji,jj) / (1.-qaxy(ji,jj))
          workl(ji,jj) = min ((zmr/(.622+zmr))*(zp/zes), 1.)
  242   continue
  240 continue
 
      do 280 jj=1,nlat
        do 282 ji=1,nlon
          teve_d(ji,jj) = teve_d(ji,jj) + ts2xy(ji,jj)
          peve_d(ji,jj) = peve_d(ji,jj) + rainaxy(ji,jj)+snowaxy(ji,jj)
          reve_d(ji,jj) = reve_d(ji,jj) + workl(ji,jj)
 
          teve_m(ji,jj) = teve_m(ji,jj) + ts2xy(ji,jj)
          peve_m(ji,jj) = peve_m(ji,jj) + rainaxy(ji,jj)+snowaxy(ji,jj)
          reve_m(ji,jj) = reve_m(ji,jj) + workl(ji,jj)
  282   continue
  280 continue
      kount_m = kount_m + 1
      kount_d = kount_d + 1
 
c        If end of day, normalize to daily means in *eve_d arrays,
c        change units to deg C, mm in prev 24 hrs, %.
 
      if (if_eod) then
        do 400 jj=1,nlat
          do 402 ji=1,nlon
            teve_d(ji,jj) = (teve_d(ji,jj)/kount_d) - 273.16
            peve_d(ji,jj) = peve_d(ji,jj) * dtgcm * (1000./rhow)
            reve_d(ji,jj) = 100.
     *                    * max (0., min (1.,reve_d(ji,jj)/kount_d))
  402     continue
  400   continue
      endif
 
c       If beginning of year, delete any existing eve_*.cur files,
c       and touch new ones (writes will be in append mode below)
 
      if (if_boy) then
        ier = ishell ('rm -f eve_temp.cur ' //
     *                      'eve_prec.cur ' //
     *                      'eve_relh.cur ' )
        ier = ishell ('touch eve_temp.cur ' //
     *                      'eve_prec.cur ' //
     *                      'eve_relh.cur ' )
      endif

c        If end of month, normalize *eve_m arrays to monthly means,
c        convert to units needed for eve_*.cur files (100ths of deg C,
c        10ths of mm/month, 100ths of %), and write(append) to those 
c        files.
 
      if (if_eom) then
 
        do 500 jj=1,nlat
          do 502 ji=1,nlon
            teve_m(ji,jj) = 100. * ((teve_m(ji,jj)/kount_m) - 273.16)
            peve_m(ji,jj) = 10.  * peve_m(ji,jj) * dtgcm*(1000./rhow)
            reve_m(ji,jj) = 10000.
     *                      * max (0., min (1., reve_m(ji,jj)/kount_m))
  502     continue
  500   continue
 
        call writevefile ('eve_temp.cur', teve_m, 1,1, 90, 1)
        call writevefile ('eve_prec.cur', peve_m, 1,1, 90, 1)
        call writevefile ('eve_relh.cur', reve_m, 1,1, 90, 1)
 
      endif

c        If end of year, use this year's data (*cur files) to update
c        multi-year running-mean data in *exp files. 

      if (if_eoy) then

        if (taubiome.gt.0.) then
          zw1 = exp (-1./taubiome)
        else
          zw1 = 0.
        endif
        zw2 = 1.-zw1

c       update temperature multi-year running-mean file:
        call readevefile ('eve_temp.cur', workcur, 12, 91, curexist)
        call readevefile ('eve_temp.exp', workexp, 12, 91, expexist)
        if (expexist) then
          do k=1,12
            do j=1,nlat
              do i=1,nlon
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              enddo
            enddo
          enddo
        else
          call scopy (nlon*nlat*12, workcur, 1, workexp, 1)
        endif
        call writevefile ('eve_temp.exp', workexp, 1,12, 91, 0)

c       update precip multi-year running-mean file:
        call readevefile ('eve_prec.cur', workcur, 12, 91, curexist)
        call readevefile ('eve_prec.exp', workexp, 12, 91, expexist)
        if (expexist) then
          do k=1,12
            do j=1,nlat
              do i=1,nlon
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              enddo
            enddo
          enddo
        else
          call scopy (nlon*nlat*12, workcur, 1, workexp, 1)
        endif
        call writevefile ('eve_prec.exp', workexp, 1,12, 91, 0)

c       update relative humidity multi-year running-mean file:
        call readevefile ('eve_relh.cur', workcur, 12, 91, curexist)
        call readevefile ('eve_relh.exp', workexp, 12, 91, expexist)
        if (expexist) then
          do k=1,12
            do j=1,nlat
              do i=1,nlon
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              enddo
            enddo
          enddo
        else
          call scopy (nlon*nlat*12, workcur, 1, workexp, 1)
        endif
        call writevefile ('eve_relh.exp', workexp, 1,12, 91, 0)

c          Rename *cur files to *cur.<year_number>
c 
c       write (cyr,'(i3.3)') mod(abs(jyear),100)
c       ier = ishell ('mv eve_temp.cur eve_temp.cur.' // cyr // ' ; ' //
c    *                'mv eve_prec.cur eve_prec.cur.' // cyr // ' ; ' //
c    *                'mv eve_relh.cur eve_relh.cur.' // cyr // ' ; ' )

c          Copy *exp files to *exp.<year_number>

        write (cyr,'(i3.3)') mod(abs(jyear),100)
        ier = ishell ('cp -p eve_temp.exp ' //
     *                     ' eve_temp.exp.' // cyr // ' ; ' //
     *                'cp -p eve_prec.exp ' //
     *                     ' eve_prec.exp.' // cyr // ' ; ' //
     *                'cp -p eve_relh.exp ' //
     *                     ' eve_relh.exp.' // cyr // ' ; ' )

      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine readevefile (cfile, arr, nval, iu, ifexist)
 
c        Reads a temp/prec/relh file for Eve model veg, in 180i6 format
c        (internal Genesis format, not final Eve format).
 
#include <compar.h>
#include <comgrd.h>
 
      character*(*) cfile
      dimension arr(nlon,nlat,nval)
      logical ifexist

c        If file doesn't exist, assume start up, zero array and return 

      inquire (file=cfile, exist=ifexist)
      if (.not.ifexist) then
        call zero (arr, nlon*nlat*nval)
        return
      endif
 
      write (ioterm2,*) '   reading file ',cfile
 
      open (iu, file=cfile, status='old')
 
      do k=1,nval
        do j=nlat,1,-1
          read (iu,'(180f6.0)') (arr(i,j,k), i=1,nlon)
        enddo
        read(iu,'()')
      enddo
      close (iu)
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine writevefile (cfile, arr, nvala, nvalb, iu, iflag)
 
c        Writes a temp/prec/relh file for Eve model veg, in 180i6 format
c        (internal Genesis format, not final Eve format).

c        iflag = 0 for new file (delete old),
c              = 1 for append to old file.
c        nvala,b = range of months to be written

#include <compar.h>
#include <comgrd.h>
 
      character*(*) cfile
      dimension arr(nlon,nlat,nvalb)
 
      write (ioterm2,*) '   writing file ',cfile
 
      if (iflag.eq.0) then 
        open (iu, file=cfile, status='unknown')
        close (iu, status='delete')
        open (iu, file=cfile, status='new', form='formatted')
      else
        open (iu, file=cfile, status='old', form='formatted',
#if defined ( cray ) || defined ( aix )
     *        position = 'append')
#else
     *        access   = 'append')
#endif
      endif
 
      do k=nvala,nvalb
        do j=nlat,1,-1
          write (iu,'(180i6)') (nint(arr(i,j,k)), i=1,nlon)
        enddo
        write(iu,'()')
      enddo
      close (iu)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine biowrite (dtgcm, nstep, nrstrt, fsave,
     *                     jyear, jmnth, jday)
 
c        Accumulate monthly-mean sfc-air temp, precip and cloudiness
c        for BIOME or SDVM veg model. Do nothing unless vegtype=5 or 7.
c        Use global arrays [T,P,C]BIOME (in comeve), written here to 
c        files biome_*.cur files at end of each month, to be read before
c        calling BIOME subr at beginning of each year.
 
c        All arguments are supplied:
c        DTGCM    = length of time to run lsx this call (=agcm timestep)
c        NSTEP    = agcm timestep counter
c        NRSTRT   = first agcm timestep # for this run
c        FSAVE    = .TRUE. if Genesis restart file-save this timestep
c                   (not used)
c        JYEAR    = current year number (2001, etc, or 1950-21000, etc)
c        JMNTH    = current month number (1-12)
c        JDAY     = current day number (1-31)
 
#include <compar.h>
#include <comeve.h>
#include <comgrd.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#include <com1d.h>
      logical fsave, curexist, expexist

      character*3 cyr
      dimension workcur(nlon,nlat,12), workexp(nlon,nlat,12)
      logical if_bod, if_eod, if_bom, if_eom, if_boy, if_eoy
      dimension ndaypm(12)
      save ndaypm, kount
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
 
c        Only do for interactive BIOME or interactive SDVM

      if (vegtype.ne.5 .and. vegtype.ne.7) return
 
c       Set beginning/end_of_day/month/year flags
 
      if_bod = .false.
      if_eod = .false.
      if_bom = .false.
      if_eom = .false.
      if_boy = .false.
      if_eoy = .false.
 
      n = nint(86400./dtgcm)
      if (nstep.eq.nrstrt .or. mod(nstep,n).eq.1) then
        if_bod = .true.
      else if ( mod(nstep,n).eq.0 ) then
        if_eod = .true.
      endif
 
      if (if_bod) then
        if (jday.eq.1) if_bom = .true.
      else if (if_eod) then
        if (jday.eq.ndaypm(jmnth)) if_eom = .true.
      endif
 
      if (if_bod .and. jday.eq. 1 .and. jmnth.eq. 1) if_boy = .true.
      if (if_eod .and. jday.eq.31 .and. jmnth.eq.12) if_eoy = .true.
 
c       If beginning of month (or first timestep), zero *biome arrays
c       and counters
 
      if (if_bom .or. nstep.eq.nrstrt) then
        call zero (tbiome(1,1), nlon*nlat)
        call zero (pbiome(1,1), nlon*nlat)
        call zero (cbiome(1,1), nlon*nlat)
        kount = 0
      endif
 
c        Accumulate current values in *biome  arrays.
 
      do jj=1,nlat
        do ji=1,nlon
          tbiome(ji,jj) = tbiome(ji,jj) + ts2xy(ji,jj)
          pbiome(ji,jj) = pbiome(ji,jj) + rainaxy(ji,jj)+snowaxy(ji,jj)
          cbiome(ji,jj) = cbiome(ji,jj) + totcldxy(ji,jj)
        enddo
      enddo
      kount = kount + 1
 
c       If beginning of year, delete any existing biome_* files,
c       and touch new ones (later open(old) below)
 
      if (if_boy) then
        ier = ishell ('rm -f biome_temp.cur ' //
     *                      'biome_prec.cur ' //
     *                      'biome_sunf.cur ' )
        ier = ishell ('touch biome_temp.cur ' //
     *                      'biome_prec.cur ' //
     *                      'biome_sunf.cur ' )
      endif
 
c        If end of month, normalize biome_* arrays to monthly means,
c        convert to units needed for biome_* files (deg C, mm/month,
c        sunshine %), and write to those files
 
      if (if_eom) then
 
        do jj=1,nlat
          do ji=1,nlon
            tbiome(ji,jj) = (tbiome(ji,jj)/kount) - 273.16
            pbiome(ji,jj) =  pbiome(ji,jj) * dtgcm * (1000./rhow)
            cbiome(ji,jj) = 100. * (1.-cbiome(ji,jj)/kount)
          enddo
        enddo
 
        open (90,file='biome_temp.cur',status='old',form='formatted',
#if defined ( cray ) || defined ( aix )
     *       position = 'append')
#else
     *       access   = 'append')
#endif
        do jj=1,nlat
           write (90,'(2i3,180f7.2)') jj,jmnth,(tbiome(ji,jj),ji=1,nlon)
        enddo
        close(90)
 
        open (90,file='biome_prec.cur',status='old',form='formatted',
#if defined ( cray ) || defined ( aix )
     *       position = 'append')
#else
     *       access   = 'append')
#endif
        do jj=1,nlat
           write (90,'(2i3,180f7.2)') jj,jmnth,(pbiome(ji,jj),ji=1,nlon)
        enddo
        close(90)
 
        open (90,file='biome_sunf.cur',status='old',form='formatted',
#if defined ( cray ) || defined ( aix )
     *       position = 'append')
#else
     *       access   = 'append')
#endif
        do jj=1,nlat
           write (90,'(2i3,180f7.2)') jj,jmnth,(cbiome(ji,jj),ji=1,nlon)
        enddo
        close(90)
 
      endif
 
c        If end of year, use this year's data (*cur files) to update
c        running-mean data in *exp files.
 
      if (if_eoy) then
 
        if (taubiome.gt.0.) then
          zw1 = exp (-1./taubiome)
        else
          zw1 = 0.
        endif
        zw2 = 1.-zw1
 
c       update temperature running-mean file:
        call readbiomefile ('biome_temp.cur', workcur, 12, 91, curexist)
        call readbiomefile ('biome_temp.exp', workexp, 12, 91, expexist)
        ier = ishell ('rm -f biome_temp.exp')
        open (91,file='biome_temp.exp',status='new',form='formatted')
        do k=1,12
          do j=1,nlat
            do i=1,nlon
              if (expexist) then
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              else
                workexp(i,j,k)= workcur(i,j,k)
              endif
            enddo
            write (91,'(2i3,180f7.2)') j,k,(workexp(i,j,k),i=1,nlon)
          enddo
        enddo
        close (91)
 
c       update precip running-mean file:
        call readbiomefile ('biome_prec.cur', workcur, 12, 91, curexist)
        call readbiomefile ('biome_prec.exp', workexp, 12, 91, expexist)
        ier = ishell ('rm -f biome_prec.exp')
        open (91,file='biome_prec.exp',status='new',form='formatted')
        do k=1,12
          do j=1,nlat
            do i=1,nlon
              if (expexist) then
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              else
                workexp(i,j,k)= workcur(i,j,k)
              endif
            enddo
            write (91,'(2i3,180f7.2)') j,k,(workexp(i,j,k),i=1,nlon)
          enddo
        enddo
        close (91)
 
c       update sunlight percent running-mean file:
        call readbiomefile ('biome_sunf.cur', workcur, 12, 91, curexist)
        call readbiomefile ('biome_sunf.exp', workexp, 12, 91, expexist)
        ier = ishell ('rm -f biome_sunf.exp')
        open (91,file='biome_sunf.exp',status='new',form='formatted')
        do k=1,12
          do j=1,nlat
            do i=1,nlon
              if (expexist) then
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              else
                workexp(i,j,k)= workcur(i,j,k)
              endif
            enddo
            write (91,'(2i3,180f7.2)') j,k,(workexp(i,j,k),i=1,nlon)
          enddo
        enddo
        close (91)
 
c          Rename *cur files to *cur.<year_number>
c 
c       write (cyr,'(i3.3)') mod(abs(jyear),100)
c       ier = ishell ('mv biome_temp.cur biome_temp.cur.' // cyr //
c    *                ' ; ' //
c    *                'mv biome_prec.cur biome_prec.cur.' // cyr //
c    *                ' ; ' //
c    *                'mv biome_sunf.cur biome_sunf.cur.' // cyr //
c    *                ' ; ' )

c          Copy *exp files to *exp.<year_number>
 
        write (cyr,'(i3.3)') mod(abs(jyear),100)
        ier = ishell ('cp -p biome_temp.exp biome_temp.exp.' // cyr //
     *                ' ; ' //
     *                'cp -p biome_prec.exp biome_prec.exp.' // cyr //
     *                ' ; ' //
     *                'cp -p biome_sunf.exp biome_sunf.exp.' // cyr //
     *                ' ; ' )
 
      endif
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine readbiomefile (cfile, arr, nval, iu, ifexist)
 
c        Reads a climate/soil data file for Biome model veg.
c        For present observed data, these have been written by
c        biomedata.f. For Genesis fields, they have been written by
c        subr biowrite in histlsx.f.
 
#include <compar.h>
#include <comgrd.h>
 
      character*(*) cfile
      dimension arr(nlon,nlat,nval)
      logical ifexist

c        If file doesn't exist, assume start up, zero array and return 

      inquire (file=cfile, exist=ifexist)
      if (.not.ifexist) then
        call zero (arr, nlon*nlat*nval)
        return
      endif
 
      write (ioterm2,*) '   reading file ',cfile
 
      open (iu, file=cfile, status='old')
 
      do k=1,nval
        do j=1,nlat
          read (iu,'(2i3,180f7.2)') jj, kk, (arr(i,j,k),i=1,nlon)
          if (jj.ne.j .or. kk.ne.k) then
            write (ioterm2,*) '*** Error (readbiomefile) file=',cfile,
     *                       ' jj=',jj,' kk=',kk,' j=',j,' k=',k
            call endrun (-1)
          endif
        enddo
      enddo
      close (iu)
 
      return
      end
c
c----------------------------------------------------------------------
c
c#######################
#if defined (WEATHERING)
c#######################
      subroutine weathwrite (dtgcm, nstep, nrstrt, fsave,
     *                       jyear, jmnth, jday)
 
c        Accumulate monthly-mean sfc-air temp, precip and total runoff
c        (sfc+drain) rates,for weathering model. Use arrays [T,P,R]WEATH
c        (local, saved), written to files weath_[temp,prec,runo].cur 
c        files at end of each month. At last timestep of year, these 
c        are combined with permanent files weath_*.exp, applying an 
c        exponential time filter to smooth out interannual variability 
c        (tauweath, years).
 
c        All arguments are supplied:
c        DTGCM    = length of time to run lsx this call (=agcm timestep)
c        NSTEP    = agcm timestep counter
c        NRSTRT   = first agcm timestep # for this run
c        FSAVE    = .TRUE. if Genesis restart file-save this timestep
c                   (not used)
c        JYEAR    = current year number (2001, etc, or 1950-21000, etc)
c        JMNTH    = current month number (1-12)
c        JDAY     = current day number (1-31)
 
#include <compar.h>
#include <comgrd.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#include <com1d.h>

      dimension tweath(nlon,nlat), pweath(nlon,nlat), rweath(nlon,nlat)
      save tweath, pweath, rweath

c     parameter (tauweath = 0.)   ! no filter
      parameter (tauweath = 2.)   ! exponential filter timescale (years)

      logical fsave, curexist, expexist

      character*3 cyr
      dimension workcur(nlon,nlat,12), workexp(nlon,nlat,12)
      logical if_bod, if_eod, if_bom, if_eom, if_boy, if_eoy
      dimension ndaypm(12)
      save ndaypm, kount
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
 
c       Set beginning/end_of_day/month/year flags
 
      if_bod = .false.
      if_eod = .false.
      if_bom = .false.
      if_eom = .false.
      if_boy = .false.
      if_eoy = .false.
 
      n = nint(86400./dtgcm)
      if (nstep.eq.nrstrt .or. mod(nstep,n).eq.1) then
        if_bod = .true.
      else if ( mod(nstep,n).eq.0 ) then
        if_eod = .true.
      endif
 
      if (if_bod) then
        if (jday.eq.1) if_bom = .true.
      else if (if_eod) then
        if (jday.eq.ndaypm(jmnth)) if_eom = .true.
      endif
 
      if (if_bod .and. jday.eq. 1 .and. jmnth.eq. 1) if_boy = .true.
      if (if_eod .and. jday.eq.31 .and. jmnth.eq.12) if_eoy = .true.
 
c       If beginning of month (or first timestep), zero *weath arrays
c       and counter
 
      if (if_bom .or. nstep.eq.nrstrt) then
        call zero (tweath(1,1), nlon*nlat)
        call zero (pweath(1,1), nlon*nlat)
        call zero (rweath(1,1), nlon*nlat)
        kount = 0
      endif
 
c        Accumulate current values in *weath  arrays
 
      do jj=1,nlat
        do ji=1,nlon
          tweath(ji,jj) = tweath(ji,jj) + ts2xy(ji,jj)
          pweath(ji,jj) = pweath(ji,jj) + rainaxy(ji,jj)+snowaxy(ji,jj)
          rweath(ji,jj) = rweath(ji,jj) +grunofxy(ji,jj)+gdrainxy(ji,jj)
        enddo
      enddo
      kount = kount + 1
 
c       If beginning of year, delete any existing weath_* files,
c       and touch new ones (later open(old) below)
 
      if (if_boy) then
        ier = ishell ('rm -f weath_temp.cur ' //
     *                      'weath_prec.cur ' //
     *                      'weath_runo.cur ' )
        ier = ishell ('touch weath_temp.cur ' //
     *                      'weath_prec.cur ' //
     *                      'weath_runo.cur ' )
      endif
 
c        If end of month, normalize weath_* arrays to monthly means,
c        convert to units needed for weath_* files (deg C, mm/day),
c        and write to those files
 
      if (if_eom) then
 
        do jj=1,nlat
          do ji=1,nlon
            tweath(ji,jj) = (tweath(ji,jj)/kount) - 273.16
c           mm/day:
            pweath(ji,jj) = (pweath(ji,jj)/kount) * (86400.* 1000./rhow)
            rweath(ji,jj) = (rweath(ji,jj)/kount) * (86400.* 1000./rhow)
c           monthly mm:
c           pweath(ji,jj) =  pweath(ji,jj) * dtgcm * (1000./rhow)
c           rweath(ji,jj) =  rweath(ji,jj) * dtgcm * (1000./rhow)
          enddo
        enddo
 
        open (90,file='weath_temp.cur',status='old',form='formatted',
#if defined ( cray ) || defined ( aix )
     *       position = 'append')
#else
     *       access   = 'append')
#endif
        do jj=1,nlat
           write (90,'(2i3,180f7.2)') jj,jmnth,(tweath(ji,jj),ji=1,nlon)
        enddo
        close(90)
 
        open (90,file='weath_prec.cur',status='old',form='formatted',
#if defined ( cray ) || defined ( aix )
     *       position = 'append')
#else
     *       access   = 'append')
#endif
        do jj=1,nlat
           write (90,'(2i3,180f7.2)') jj,jmnth,(pweath(ji,jj),ji=1,nlon)
        enddo
        close(90)
 
        open (90,file='weath_runo.cur',status='old',form='formatted',
#if defined ( cray ) || defined ( aix )
     *       position = 'append')
#else
     *       access   = 'append')
#endif
        do jj=1,nlat
           write (90,'(2i3,180f7.2)') jj,jmnth,(rweath(ji,jj),ji=1,nlon)
        enddo
        close(90)
 
      endif
 
c        If end of year, use this year's data (*cur files) to update
c        running-mean data in *exp files.
 
      if (if_eoy) then
 
        if (tauweath.gt.0.) then
          zw1 = exp (-1./tauweath)
        else
          zw1 = 0.
        endif
        zw2 = 1.-zw1
 
c       update temperature running-mean file:
        call readweathfile ('weath_temp.cur', workcur, 12, 91, curexist)
        call readweathfile ('weath_temp.exp', workexp, 12, 91, expexist)
        ier = ishell ('rm -f weath_temp.exp')
        open (91,file='weath_temp.exp',status='new',form='formatted')
        do k=1,12
          do j=1,nlat
            do i=1,nlon
              if (expexist) then
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              else
                workexp(i,j,k)= workcur(i,j,k)
              endif
            enddo
            write (91,'(2i3,180f7.2)') j,k,(workexp(i,j,k),i=1,nlon)
          enddo
        enddo
        close (91)

c       update precipitation running-mean file:
        call readweathfile ('weath_prec.cur', workcur, 12, 91, curexist)
        call readweathfile ('weath_prec.exp', workexp, 12, 91, expexist)
        ier = ishell ('rm -f weath_prec.exp')
        open (91,file='weath_prec.exp',status='new',form='formatted')
        do k=1,12
          do j=1,nlat
            do i=1,nlon
              if (expexist) then
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              else
                workexp(i,j,k)= workcur(i,j,k)
              endif
            enddo
            write (91,'(2i3,180f7.2)') j,k,(workexp(i,j,k),i=1,nlon)
          enddo
        enddo
        close (91)
 
c       update runoff running-mean file:
        call readweathfile ('weath_runo.cur', workcur, 12, 91, curexist)
        call readweathfile ('weath_runo.exp', workexp, 12, 91, expexist)
        ier = ishell ('rm -f weath_runo.exp')
        open (91,file='weath_runo.exp',status='new',form='formatted')
        do k=1,12
          do j=1,nlat
            do i=1,nlon
              if (expexist) then
                workexp(i,j,k) = zw1*workexp(i,j,k) + zw2*workcur(i,j,k)
              else
                workexp(i,j,k)= workcur(i,j,k)
              endif
            enddo
            write (91,'(2i3,180f7.2)') j,k,(workexp(i,j,k),i=1,nlon)
          enddo
        enddo
        close (91)
 
c          Rename *cur files to *cur.<year_number>
c 
c       write (cyr,'(i3.3)') mod(abs(jyear),100)
c       ier = ishell ('mv weath_temp.cur weath_temp.cur.' // cyr //
c    *                ' ; ' //
c    *                'mv weath_prec.cur weath_prec.cur.' // cyr //
c    *                ' ; ' //
c    *                'mv weath_runo.cur weath_runo.cur.' // cyr //
c    *                ' ; ' )

c          Copy *exp files to *exp.<year_number>

c       write (cyr,'(i3.3)') mod(abs(jyear),100)
c       ier = ishell ('cp -p weath_temp.exp weath_temp.exp.' // cyr //
c    *                ' ; ' //
c    *                'cp -p weath_prec.exp weath_prec.exp.' // cyr //
c    *                ' ; ' //
c    *                'cp -p weath_runo.exp weath_runo.exp.' // cyr //
c    *                ' ; ' )

      endif
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine readweathfile (cfile, arr, nval, iu, ifexist)
 
c        Reads a climate data file for weathering model.
c        They have been written by subr weathwrite in histlsx.F.
 
#include <compar.h>
#include <comgrd.h>
 
      character*(*) cfile
      dimension arr(nlon,nlat,nval)
      logical ifexist

c        If file doesn't exist, assume start up, zero array and return 

      inquire (file=cfile, exist=ifexist)
      if (.not.ifexist) then
        call zero (arr, nlon*nlat*nval)
        return
      endif
 
      write (ioterm2,*) '   reading file ',cfile
 
      open (iu, file=cfile, status='old')
 
      do k=1,nval
        do j=1,nlat
          read (iu,'(2i3,180f7.2)') jj, kk, (arr(i,j,k),i=1,nlon)
          if (jj.ne.j .or. kk.ne.k) then
            write (ioterm2,*) '*** Error (readweathfile) file=',cfile,
     *                       ' jj=',jj,' kk=',kk,' j=',j,' k=',k
            call endrun (-1)
          endif
        enddo
      enddo
      close (iu)
 
      return
      end
c#####
#endif
c#####
c
c----------------------------------------------------------------------
c
      subroutine regwrilsx (dtgcm, nstep, nrstrt, istepy, mbdate,
     *                      jyear, jmnth, jday, mcsec, runtitle)
 
c       Write instantaneous lsxm fields at specified interval (namelist
c       param REGFREQ) to a stream of files, one file per month, for
c       driving regional model RegCM2 (one-way nesting).
 
c       All arguments are supplied:
c       DTGCM    = length of time to run lsx this call (=agcm timestep)
c       NSTEP    = agcm timestep counter
c       NRSTRT   = first agcm timestep # for this run
c       ISTEPY   = timestep within current 365-day year
c       MBDATE   = base date for curr. sequence of runs (yymmdd)
c       JYEAR    = current year number (2001 AD,etc, or -21000 BP,etc)
c       JMNTH    = current month number (1-12)
c       JDAY     = current day number (1-31)
c       MCSEC    = seconds into current day (at end of timestep)
c       RUNTITLE = run title
 
#include <compar.h>
#include <comgrd.h>
#include <comhisl.h>
#include <comice.h>
#include <comicexy.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commudxy.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
 
      character*80 runtitle
 
      logical if_bod, if_eod, if_bom, if_eom
      dimension ndaypm(12)
      character cmon(12)*3, cyr*20, cname*200
      save ndaypm, cmon, cname
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data cmon /'JAN','FEB','MAR','APR','MAY','JUN',
     *           'JUL','AUG','SEP','OCT','NOV','DEC' /
 
      dimension
     *  regprecip(nlon,nlat), regtmin(nlon,nlat),  regtmax(nlon,nlat)
      logical first
      save regprecip, regtmin, regtmax, first, kount
      data first /.true./
 
      if (first) then
        kount = 0
        call zero (regprecip, nlon*nlat)
        call resetr (regtmin, nlon*nlat,  999.)
        call resetr (regtmax, nlon*nlat, -999.)
        first = .false.
      endif
 
c        If user has not specified permanent filename, do nothing
 
      if (reghis_a.eq.' ') return
 
c       Set beginning/end_of day/month flags
 
      if_bom = .false.
      if_eom = .false.
      if_bod = .false.
      if_eod = .false.
 
      n = nint(86400./dtgcm)
      if (nstep.eq.nrstrt .or. mod(nstep,n).eq.1) then
        if_bod = .true.
      else if ( mod(nstep,n).eq.0 ) then
        if_eod = .true.
      endif
 
      if (if_bod) then
        if (jday.eq.1) if_bom = .true.
      else if (if_eod) then
        if (jday.eq.ndaypm(jmnth)) if_eom = .true.
      endif
 
c     write (*,888) nstep, if_bod, if_eod, if_bom, if_eom
c888  format ('regwrilsx: nstep=',i6,'  if[bod,eod,bom,eom]=',4l3)
 
 
c        If beginning of month, open new file.
c        First assemble permanent name by appending "_year_month" to
c        requested prefix (namelist parameter REGHIS_A).
 
      if (if_bom) then
         write (cyr,'(i20)') jyear
         call lstrip (cyr)
         cname = reghis_a(1:lenchr(reghis_a)) // '_' //
     *          cyr(1:lenchr(cyr))       // '_' //
     *          cmon(jmnth)
         call openf (iureg_a, cname, regpath, 1, 0)
      endif
 
 
c        Accumulate time-mean precip, and track 24-hr max/min 2-m air
c        temperature (not needed by RegCM2, but included here)
 
      kount = kount + 1
      do 20 jj=1,nlat
        do 22 ji=1,nlon
          regprecip(ji,jj) = regprecip(ji,jj)
     *                     + ( rainaxy(ji,jj) + snowaxy(ji,jj) )
          regtmin(ji,jj) = min (regtmin(ji,jj), ts2xy(ji,jj))
          regtmax(ji,jj) = max (regtmax(ji,jj), ts2xy(ji,jj))
   22   continue
   20 continue
 
c        If at requested output time, write header and one set of
c        fields; normalize, write and reset precip accumulators.
 
      n = max (nint(regfreq/dtgcm), 1)
      if (mod(nstep,n).eq.0) then
 
 
c       Header info:
        nlevdum = 1
        write (iureg_a) typehis_a, runtitle, nlon, nlat, nlevdum,
     *                  jyear, jmnth, jday, mcsec
 
c       Land/icesheet/ocean mask (1=land, 2=icesheet, 3=ocean/seaice):
        do 80 jj=1,nlat
          do 82 ji=1,nlon
            workl(ji,jj) = lmaskxy(ji,jj)
   82     continue
   80   continue
        write (iureg_a) 'LMASK   ', workl
 
c       "Solid" top-layer temperature:
        do 100 jj=1,nlat
          do 102 ji=1,nlon
            workl(ji,jj) = 0.
            x = (1. - fixy(ji,jj))*tsoixy(ji,jj,1)
     *              + fixy(ji,jj) *tsnoxy(ji,jj,1)
            y = (1. - fi_icexy(ji,jj))*ticexy(ji,jj,1)
     *              + fi_icexy(ji,jj) *tsno_icexy(ji,jj,1)
            z = tocnxy(ji,jj)
            workl(ji,jj)=(1.-fwaterxy(ji,jj))*                        x
     *                 +     fwaterxy(ji,jj) * (      ficexy(ji,jj) * y
     *                                          + (1.-ficexy(ji,jj))* z)
  102     continue
  100   continue
        write (iureg_a) 'TS      ', workl
 
c       Snow fraction:
        do 110 jj=1,nlat
          do 112 ji=1,nlon
            zf = fwaterxy(ji,jj)
            workl(ji,jj) = (1.-zf)*fixy(ji,jj)
     *                   +     zf *fi_icexy(ji,jj)
  112     continue
  110   continue
        write (iureg_a) 'SNOWF   ', workl
 
c       Snow thickness:
        call zero (workl, nlon*nlat)
        do 120 k=1,nsnolay
          do 122 jj=1,nlat
            do 124 ji=1,nlon
              zf = fwaterxy(ji,jj)
              workl(ji,jj) = workl(ji,jj) + (1.-zf)*hsnoxy(ji,jj,k)
     *                                    +     zf *hsno_icexy(ji,jj,k)
  124       continue
  122     continue
  120   continue
        write (iureg_a) 'SNOWH   ', workl
 
c       Sea-ice fraction:
        write (iureg_a) 'ICEF    ', ficexy
 
c       Sea-ice thickness:
        call zero (workl, nlon*nlat)
        do 130 k=1,nicelay
          do 132 jj=1,nlat
            do 134 ji=1,nlon
              workl(ji,jj) = workl(ji,jj) + hicexy(ji,jj,k)
  134       continue
  132     continue
  130   continue
        write (iureg_a) 'ICEH    ', workl
 
c       Normalize, write, reset precip accumulators:
        if (kount.ge.1) then
          do 140 jj=1,nlat
            do 142 ji=1,nlon
              regprecip(ji,jj) = (regprecip(ji,jj)/kount)/rhow     ! m/s
  142       continue
  140     continue
        endif
        write (iureg_a) 'PRECIP  ', regprecip
        kount = 0
        call zero (regprecip, nlon*nlat)
 
c       24-hr min 2-m temperature:
        write (iureg_a) 'T2MIN   ', regtmin
 
c       24-hr max 2-m temperature:
        write (iureg_a) 'T2MAX   ', regtmax
 
      endif
 
c        If at end of a day (GMT), reset min/max temperature trackers
 
      if (if_eod) then
        call resetr (regtmin, nlon*nlat,  999.)
        call resetr (regtmax, nlon*nlat, -999.)
      endif
 
c        If end of month, dispose local file to permanent location.
 
      if (if_eom) then
        call closef (iureg_a)
        call disposef (cname, regpath)
 
        write (ioterm2,200)
     *    nstep, istepy*dtgcm/86400., nstep*dtgcm/86400.,
     *    regpath(1:lenchr(regpath)), cname(1:lenchr(cname))
  200   format(
     *         /' ----------------------------------------------'
     *         /' ******* LSX RegCM2-DRIVER FILE OUTPUT ********'
     *         /' ----------------------------------------------'
     *         /' timestep                          = ',i10
     *         /' calendar days                     = ',f10.3
     *         /' elapsed days                      = ',f10.3
     *         /' lsx regcm2-driver file name       = ',a,a
     *         /' ----------------------------------------------')
#ifdef mss
        if (regpath.eq.' ') write (ioterm2,202)
  202   format( ' *  File NOT saved to MSS (REGPATH is blank!) *'
     *         /' ----------------------------------------------')
#endif
        write (ioterm2,*)
 
      endif
 
      return
      end
