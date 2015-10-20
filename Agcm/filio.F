c -- $Id: filio.F,v 1.2 2001/07/03 21:29:18 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/filio.F,v $
c -- $Name:  $

c
c ***************************** OPENF *********************************
c
      subroutine openf (iu, cname, cpath, ifnew, ialloc)
c----------------------------------------------------------------------
c
c            Opens a file.
c
c            iu    = logical unit number
c            cname = local file name
c            cpath = local path name if not mss, ignored if mss
c            ifnew = -1 for existing file, position=append if history,
c                     0 for existing file, position=rewind,
c                     1 for new file (old is deleted)
c            ialloc = pre-allocated file size in words, or 0 if not
c                     needed (only used for new history files)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMHIS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
      character*(*) cname, cpath
      character concat*96, cmd*200, cstat*8, cpos*8, cform*16, lname*96
      character*200 ctmp
      character*8 cduma, cconv
c----------------------------------------------------------------------
c for netCDF:
#ifdef NETCDF
#include <netcdf.inc>
c     include '/usr/local/include/netcdf.inc'
#ifdef cray
      integer 
#else
      integer*4 
#endif
     * rcode, ncerr,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls. 
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
#endif
c-----------------------------------------------------------------------

      ier = 0
      ios = 0
 
#ifdef mss
      lname = cname
      lenl = lenchr(lname)
#else
      lname = concat (cpath, cname, lenl)
#endif
 
      if (ifnew.eq.-1) then
        cstat = 'old'
        if (iu.eq.ndata .or. iu.eq.ndata_a .or. iu.eq.iustorm) then
          cpos = 'append'
        else
          cpos = 'rewind'
        endif
      else if (ifnew.eq.0) then
        cstat = 'old'
        cpos = 'rewind'
      else
        cstat = 'new'
        cpos = 'rewind'
      endif
 
 
      if (ifnew.eq.1) call deletef (iu, lname)
 
c---------------------------------
c netcdf history files (agcm, lsx)
c---------------------------------
 
      if ( (iu.eq.ndata .or. iu.eq.ndata_a) .and. hisnet) then

#ifdef NETCDF

c       append '.nc' to name if not already there
        lnl = lenchr(lname)
        if  (lname(lnl-2:lnl).ne.'.nc'
     *  .and.lname(lnl-2:lnl).ne.'.NC') then
            lname(lnl+1:) = '.nc'
            lenl = lenchr(lname)
            if (ifnew.eq.1) call deletef (iu, lname)
        endif

c          netCDF history file:

        if (ifnew.eq.1) then

c         new netCDF file:
          if (iu.eq.ndata) then
            rcode = nf_create (lname, nf_clobber4, ncid)    !ncid,comlun
            numfielnet = 0             ! # of distinct fields, in comlun
          else if (iu.eq.ndata_a) then
            rcode = nf_create (lname, nf_clobber4,ncid_a) !ncid_a,comlun
            numfielnet_a = 0
          endif

#ifdef cray
c         pre-allocate disk space on Cray:
          if (ialloc.eq.0) return
          ia = (ialloc-1)/512 + 1
          ia = min (ia, 1000000)
          write(cmd,52) lname(1:lenl), ia, lname(1:lenl)
   52     format('assign -a ',a, ' -n ',i8, ' ', a)
          call assign (cmd, ier)
          if (ier.ne.0) goto 900
#endif

        else 

c         existing netCDF file (append):
          if (iu.eq.ndata) then
            rcode = nf_open(lname, nf_write4, ncid)
            numfielnet = 0
          else if (iu.eq.ndata_a) then
            rcode = nf_open(lname, nf_write4, ncid_a)
            numfielnet_a = 0
          endif

        endif

        if (rcode.ne.nf_noerr4) then
          if (iu.eq.ndata)   iz = ncid
          if (iu.eq.ndata_a) iz = ncid_a
          ctmp = nf_strerror (ncerr)
          write(nout,60) lname(1:lenl), cstat, iz, rcode,
     *                    ncerr, ctmp(1:lenchr(ctmp))
   60     format(/' *** Error in OPENF (opening a NetCDF file):'
     *           /'     filename = ',a
     *           /'     status   = ',a
     *           /'     ncid     = ',i8
     *           /'     rcode    = ',i8
     *           /'     ncerr    = ',i8
     *           /'     strerror = ',a)
          call endrun (-1)
        endif

        return

#endif
 
c----------------------------------------------
c binary history files (agcm, lsx, storm-track)
c----------------------------------------------

      else if (iu.eq.ndata .or. iu.eq.ndata_a .or. iu.eq.iustorm) then
 
c----------
#ifdef cray
c----------
c           If new file, set ia = preallocated disk space in 512-wd
c           blocks; (if ialloc=0, then no history fields requested this
c           run, so return without opening).
c           For new or old files, set ib = buffer size in 512-wd blocks.
 
         if (ifnew.eq.1) then
           if (ialloc.eq.0) return
           ia = (ialloc-1)/512 + 1
           ia = min (ia, 1000000)
         endif
 
         if (iu.eq.ndata) then
           ib = (2*nurecl-1)/512 + 1
         else if (iu.eq.ndata_a) then
           ib = (2*nurecl_a-1)/512 + 1
         else if (iu.eq.iustorm) then
           ib = (2*(1+nlon*norec)-1)/512 + 1
         endif
 
         if (ifnew.le.0) then
           write(cmd,100) lname(1:lenl), ib, lname(1:lenl)
  100      format('assign -a ',a, ' -b ',i8, ' -s bin ', a)
cz         write(cmd,100) lname(1:lenl),     lname(1:lenl)
cz100      format('assign -a ',a,            ' -s bin ', a)
         else
           write(cmd,102) lname(1:lenl), ib, ia, lname(1:lenl)
  102      format('assign -a ',a, ' -b ',i8, ' -n ',i8, ' -s bin ', a)
cz         write(cmd,102) lname(1:lenl),         lname(1:lenl)
cz102      format('assign -a ',a,                       ' -s bin ', a)
         endif
         call assign (cmd, ier)
         if (ier.ne.0) goto 900
 
         cform = 'unformatted'
c----
#else
c----
#ifdef sgi
         cform = 'system'
#else
         cform = 'unformatted'
#endif
c-----
#endif
c-----
 
#if defined ( sgi )
c position='append' doesn't work with Sgi "system" files. Use gotoeof:
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         form=cform, iostat=ios, err=900)
         if (ios.ne.0) goto 900
         if (cpos.eq.'append') then
           call gotoeof (lname(1:lenl), iu)
         endif
#elif defined ( cray ) || defined ( aix )
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         form=cform, position=cpos, iostat=ios, err=900)
         if (ios.ne.0) goto 900
#elif defined ( linux )
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         recordtype='stream', carriagecontrol='none',
     *         form=cform, position=cpos, iostat=ios, err=900)
         if (ios.ne.0) goto 900
#elif defined ( sun )
         if (cstat.eq.'old' .and. cpos.eq.'append') then
           open (unit=iu, file=lname(1:lenl), status=cstat,
     *           form=cform, access  =cpos, iostat=ios, err=900)   ! f77
c    *           form=cform, position=cpos, iostat=ios, err=900)   ! f90
         else
           open (unit=iu, file=lname(1:lenl), status=cstat,
     *           form=cform,                iostat=ios, err=900)
         endif
         if (ios.ne.0) goto 900
#else
         write (6,*) 'Error (openf): unknown machine type'
         call endrun (-1)
#endif

c--------------------------------------------------------------
c binary restart file (GENESIS or IBIS), or RegCM2 driver files
c--------------------------------------------------------------
 
      else if (iu.eq.nsre .or. iu.eq.iureg .or. iu.eq.iureg_a .or.
     *         iu.eq.nsreibis) then
 
c----------
#ifdef cray
c----------
         if (ifnew.le.0) then
           write(cmd,110) lname(1:lenl), lname(1:lenl)
  110      format('assign -a ',a, ' -b 256 ', ' -s bin ', a)
cz110      format('assign -a ',a,             ' -s bin ', a)
         else
c          ia is new-file preallocated disk space in 512-wd blks
           if (iu.eq.nsre) then
             ia = (19000000/8 - 1)/512 + 1         ! T31
           else
             ia = ( 3000000/8 - 1)/512 + 1         ! T31, one agcm write
           endif
           ia = min (ia, 1000000)
           write(cmd,112) lname(1:lenl), ia, lname(1:lenl)
  112      format('assign -a ',a, ' -b 256 ', ' -n ',i8, ' -s bin ', a)
cz         write(cmd,112) lname(1:lenl),     lname(1:lenl)
cz112      format('assign -a ',a,                        ' -s bin ', a)
         endif
         call assign (cmd, ier)
         if (ier.ne.0) goto 900
 
         cform = 'unformatted'
c----
#else
c----
#ifdef sgi
         cform = 'system'
#else
         cform = 'unformatted'
#endif
c-----
#endif
c-----

#if defined ( linux ) 
       if (cstat.eq.'old') then

c        check if Cray restart file...if so, use "convert='cray'":
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         recordtype='stream', carriagecontrol='none',
     *         form=cform, iostat=ios, err=900)
         read (iu, err=950, iostat=ios) cduma
         if (cduma(1:4).eq.'CRAY') then
           cconv = 'cray'
         else
           cconv = 'native'
         endif
         close (iu)
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         recordtype='stream', carriagecontrol='none',
     *         convert=cconv,
     *         form=cform, iostat=ios, err=900)
       else
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         recordtype='stream', carriagecontrol='none',
     *         form=cform, iostat=ios, err=900)
       endif
#else
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         form=cform, iostat=ios, err=900)
#endif
         if (ios.ne.0) goto 900
 
c-----------------------
c ascii data input files
c-----------------------
 
      else if ( iu.eq.ntopog .or. iu.eq.nozone   .or. iu.eq.ngwdrag .or.
     *          iu.eq.ntopog2.or. iu.eq.iuvegt   .or. iu.eq.iuwisc  .or.
     *          iu.eq.iusurt .or. iu.eq.iuwater  .or. iu.eq.nocean  .or.
     *          iu.eq.nsstice.or.
     *          iu.eq.iuocean_i.or. iu.eq.iusourc_i.or. 
     *          iu.eq.iunamel  .or.
     *          (iu.ge.iusoit(1).and.iu.le.iusoit(psoilun)) ) then
         cform = 'formatted'
         open (unit=iu, file=lname(1:lenl), status=cstat,
     *         form=cform, iostat=ios, err=900)
         if (ios.ne.0) goto 900
 
c----------
      endif
c----------
 
c        Normal return
 
      if (cpos.eq.'rewind') rewind (iu)
 
      return
 
c        Fatal error for assign (ier != 0) or open (ios != 0)
 
  900 if (ier.ne.0) then
        write(nout,902) lname(1:lenl), cmd(1:lenchr(cmd)), ier
  902   format(/' *** Error in OPENF (doing system command):'
     *         /'     filename       = ',a
     *         /'     system command = ',a
     *         /'     error code     = ',i8)
        iexpl = ier
      else
        write(nout,904) lname(1:lenl), iu, cstat, cform, cpos, ios
  904   format(/' *** Error in OPENF (opening a file):'
     *         /'     filename       = ',a
     *         /'     unit number    = ',i8
     *         /'     status         = ',a
     *         /'     form           = ',a
     *         /'     position/access= ',a
     *         /'     iostat code    = ',i8)
        iexpl = ios
      endif
 
#ifdef cray
      if (iexpl.le.900) then
        write (cmd,'("explain sys-",i8.8)') iexpl
      else
        write (cmd,'("explain lib-",i8.8)') iexpl
      endif
      write (nout,906) cmd(1:lenchr(cmd))
  906 format('     explanation (',a,'):')
      write (nout,'(40("="))')
      ier = flush(nout)
      jer = ishell(cmd(1:lenchr(cmd)))
      write (nout,'(40("="))')
#else
c perror done in endrun
c     write (nout,'(40("="))')
c     call perror (cmd)
c     write (nout,'("gerror: ",a)') cmd(1:lenchr(cmd))
c     write (nout,'(40("="))')
c     ier = flush(nout)
#endif
      call endrun (-1)

#if defined ( linux )
c        Fatal error for reading first byte of a restart file (linux)

  950 write(nout,952) lname(1:lenl), ios
  952 format(/' *** Error in reading first byte of restart file(OPENF):'
     *       /'     filename    = ',a
     *       /'     iostat code = ',i8)
      call endrun (-1)
#endif
 
      end
c
c****************************** CLOSEF ********************************
c
      subroutine closef (iu)
c----------------------------------------------------------------------
c
c              Closes a file.
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
c for netCDF:
#ifdef NETCDF
#include <netcdf.inc>
#ifdef cray
      integer 
#else
      integer*4 
#endif
     * rcode,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls. 
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
#endif
c-----------------------------------------------------------------------


      if (iu.eq.ndata .and. hisnet) then
#ifdef NETCDF
c       close netCDF agcm file
        rcode = nf_close (ncid)      ! ncid is in comlun

      else if (iu.eq.ndata_a .and. hisnet) then

c       close netCDF lsx file
        rcode = nf_close (ncid_a)    ! ncid_a is in comlun
#endif
      else

        close (iu, err=900)
  900   continue

      endif
 
      return
      end
c
c ***************************** ACQUIREF ******************************
c
      subroutine acquiref (cname, cpath)
c----------------------------------------------------------------------
c
c            Copies a file from Mass Store to local disk (if mss).
c            Does nothing if not mss, or if cpath = ' '.
c
c            cname = file name (supplied)
c            cpath = mss path name (supplied)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMHIS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
      character*(*) cname, cpath
      character concat*96, cmd*200, lname*96, mname*96, cpass*20
      dimension itmp(6)
c----------------------------------------------------------------------
 
 
#ifdef mss
 
      if (cpath.ne.' ') then
 
         ier = 0
 
         lname = cname
         lenl = lenchr(lname)
         mname = concat (cpath, cname, lenm)
 
         cpass = passrd
 
         call msread (ier, lname(1:lenl), mname(1:lenm),
     *                cpass, 'lcount=15, replace')
 
         if (ier.ne.0) then
           call mserror (cmd)
           write(nout,900) mname(1:lenm), lname(1:lenl), cpass, ier,
     *                     cmd(1:lenchr(cmd))
  900      format(/' *** Error in ACQUIREF (doing msread):'
     *            /'     MSS   file            = ',a
     *            /'     Local file            = ',a
     *            /'     read password         = ',a
     *            /'     msread error code     = ',i8
     *            /'     explanation (mserror) = ',a)
           call endrun (-1)
         else
            write(nout,950) lname(1:lenl), mname(1:lenm)
  950       format(/' Local file ',a,' acquired from MSS file ',a)
         endif
 
      endif
 
#endif
 
      return
      end
c
c ***************************** DISPOSEF ******************************
c
      subroutine disposef (cname, cpath)
c----------------------------------------------------------------------
c
c            Copies a file from local disk to Mass Store (if mss).
c            Does nothing if not mss, or if cpath = ' '.
c
c            cname = file name (supplied)
c            cpath = mss path name (supplied)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
      character*(*) cname, cpath
      character concat*96, cmd*200, lname*96, mname*96, cpass*20
c----------------------------------------------------------------------
 
 
#ifdef mss
 
      if (cpath.ne.' ') then
 
         ier = 0
 
         lname = cname
         lenl = lenchr(lname)
         mname = concat (cpath, cname, lenm)
 
         if (passrd.ne.' ' .and. passwd.ne.' ') then
           cpass =passrd(1:lenchr(passrd)) // ',' // passwd
         else if (passrd.ne.' ' .and. passwd.eq.' ') then
           cpass = passrd
         else if (passrd.eq.' ' .and. passwd.ne.' ') then
           cpass = ',' // passwd
         else
           cpass = ' '
         endif
 
         irtpd = nint(retpd)
 
         call mswrite (ier, lname(1:lenl), mname(1:lenm),
     *                 cpass, irtpd, 'lcount=15')
 
         if (ier.ne.0) then
           call mserror (cmd)
           write(nout,900) lname(1:lenl), mname(1:lenm), cpass,
     *                     nint(retpd), ier, cmd(1:lenchr(cmd))
  900      format(/' *** Error in DISPOSEF (doing mswrite):'
     *            /'     Local file            = ',a
     *            /'     MSS   file            = ',a
     *            /'     read,write passwords  = ',a
     *            /'     retention period      = ',i8
     *            /'     mswrite error code    = ',i8
     *            /'     explanation (mserror) = ',a)
c          if fatalwms=.f., let program continue even if dispose fails
           if (fatalwms) call endrun (-1)
         else
           write(nout,950) lname(1:lenl), mname(1:lenm)
  950      format(/' Local file ',a,' disposed to MSS file ',a)
         endif
 
       endif
 
#endif
 
      return
      end
c
c ***************************** DELETEF ********************************
c
      subroutine deletef (iu, cname)
c----------------------------------------------------------------------
c
c            Deletes a file. Called only from OPENF if opening
c            a new file (old is deleted).
c
c            iu = unit number (supplied)
c            cname = local file name (including path if not mss)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
      character*(*) cname
      character cmd*200, lname*96
c----------------------------------------------------------------------
 
 
      call closef (iu)
 
      lname = cname
      lenl = lenchr(lname)
 
#ifdef cray
      write(cmd,48) lname(1:lenl), lname(1:lenl)
   48 format('assign -a ',a,' ',a)
      call assign (cmd, ier)
#endif
      open  (iu, file=lname(1:lenl), status='unknown', err=50)
   50 close (iu, status='delete', err=52)
   52 continue
 
c     Alternatively:
c     ier = ishell ('rm -f ' // lname)
 
      return
      end
c
c ***************************** FINDIU ********************************
c
      subroutine findiu (iu, loc)
c----------------------------------------------------------------------
c
c            Finds index loc into file-info arrays for unit number iu.
c
c            iu = unit number (supplied)
c            loc = index (returned)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
 
 
      do 100 i=1,plun
        if (iu.eq.lociu(i)) then
          loc = i
          return
        endif
  100 continue
 
      write(nout,900) iu
  900 format(/' *** Error in FINDIU. Unknown unit number iu=',I6)
      call endrun (-1)
 
      end
c
c ***************************** CONCAT ********************************
c
      character*(*) function concat (ch1, ch2, len)
c----------------------------------------------------------------------
c
c            Concatanates ch1 and ch2 omitting trailing blanks/zeros,
c            and returns result. Also returns length of result in len.
c            (Used in forming permanent filenames in acquiref,disposef.)
c
c----------------------------------------------------------------------
      character*(*) ch1, ch2
c----------------------------------------------------------------------
 
      len1 = lenchr(ch1)
      len2 = lenchr(ch2)
 
      if (ch1(1:len1).eq.' ') then
         concat = ch2(1:len2)
         len = len2
      else if (ch2(1:len2).eq.' ') then
         concat = ch1(1:len1)
         len = len1
      else
         if (ch1(len1:len1).eq.'/') then
           concat = ch1(1:len1) // ch2(1:len2)
           len = len1 + len2
         else
           concat = ch1(1:len1) // '/' // ch2(1:len2)
           len = len1 + 1 + len2
         endif
      endif
 
      return
      end
c
c ***************************** SIZECHK *******************************
c
      subroutine sizechk
c----------------------------------------------------------------------
c
c            Checks that estimated sizes of all types of history files
c            do not exceed system maximum. Called from PHYS after lsx
c            initialiazation complete.
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
#include <COMLUN.h>
c----------------------------------------------------------------------
#include <COMTIM.h>
c----------------------------------------------------------------------
 
#ifdef mss
c     msslimit =  400000000   ! NCAR Mass Store file limit (bytes)
      msslimit = 6000000000   ! NCAR Mass Store file limit (bytes)
      if ( 8*nusiz  .gt.msslimit .or. 8*nusiz_a.gt.msslimit ) then
        write (nout,100) 8*nusiz, 8*nusiz_a, msslimit
  100   format(/' *** The estimated size of at least one type of',
     *          ' history file exceeds the'
     *         /'     Mass Store limit (i.e., the capacity',
     *          ' of one MSS tape cartridge).'
     *         /'     Estimated AGCM history file size (bytes):',i12
     *         /'     estimated LSX  history file size (bytes):',i12
     *         /'     MSS limit                        (bytes):',i12
     *         /)
        call endrun (-1)
      endif
#endif
 
      return
      end
c
c ***************************** GOTOEOF *******************************
c
c---------
#ifdef sgi
c---------
      subroutine gotoeof (cname, iu)
 
c        Positions local file cname at eof. File must be open with
c        unit number iu. Only called from openf for existing history
c        files, for SGI (since Fortran open with position='append'
c        doesn't work for SGI "system" files).

c        Used to not work on Cray...now does, with f90. 
c        Keep commented-out Cray code in case ever useful.
 
      character*(*) cname
      character*80  ctmp
 
c------------------
c#if defined (cray)
c------------------
c 
c     For Cray:
c     Uses stat and setpos (some returned words from stat, and
c     setpos itself, are Cray specific). stat gives size of file in
c     bytes, then setpos positions at eof.
c
c     integer stat, istat(200) ! make istat big enough to hold all
c                              ! returned values (53 on Unicos 6.1)
c
c     stat is used to find size of file. Other returned values:
c     istat(1)  - File mode
c     istat(2)  - I-node number
c     istat(3)  - Device
c     istat(4)  - Device ID
c     istat(5)  - Number of links
c     istat(6)  - User ID of file's owner
c     istat(7)  - Group ID of file's group
c     istat(8)  - Size in bytes
c     istat(9)  - Time last accessed (secs)
c     istat(10) - Time last modified
c     istat(11) - Time status last changed
c 
c        Find size of file using "stat" (must first copy cname to
c        local ctmp since cname*(*) can't be used in actual-arg exprns)
c
c     ctmp = cname(1:len(cname))
c     ier = stat (ctmp(1:lenchr(ctmp))//char(0), istat)
c     nbytes = istat(8)
c     if (ier.ne.0) then
c       write (*,900) ctmp(1:lenchr(ctmp)), iu, nbytes, ier
c 900   format (/' *** Error in GOTOEOF, calling STAT.'
c    *          /'     cname = ',a
c    *          /'     iu    = ',i8
c    *          /'     nbytes= ',i8
c    *          /'     ier   = ',i8 )
c       call endrun (-1)
c     endif
c
c        Use setpos to position to (end of) last word
c
c     iw = (nbytes-1)/8 + 1
c     call setpos (iu, 1, iw, ios)
c     if (ios.ne.0) then
c       write (*,910) ctmp(1:lenchr(ctmp)), iu, nbytes, ios
c 910   format (/' *** Error in GOTOEOF, calling SETPOS.'
c    *          /'     cname = ',a
c    *          /'     iu    = ',i8
c    *          /'     nbytes= ',i8
c    *          /'     ios   = ',i8 )
c       call endrun (-1)
c     endif
c
c----
c#else
c----
 
c     For Sgi: Uses fseek64 to position to eof.
 
      integer*8 fseek64, ioffset, ios
      integer*4 iuseek,  ifrom
 
      iuseek = iu
      ioffset = 0
      ifrom = 2    ! from eof
 
c#ifdef sgi
       ios = fseek64 (iuseek, ioffset, ifrom)
c#else
c      ios = fseek   (iuseek, ioffset, ifrom)    ! for Sun (not used)
c#endif
 
      if (ios.ne.0) then
        ctmp = cname(1:len(cname))
        write (*,910) ctmp(1:lenchr(ctmp)), iu, ios
  910   format (/' *** Error in GOTOEOF, calling FSEEK64.'
     *          /'     cname = ',a
     *          /'     iu    = ',i8
     *          /'     ios   = ',i8 )
        call endrun (-1)
      endif

c------
c#endif
c------
 
      return
      end
c-----
#endif
c-----
c
c ***************************** CONVERTYPE *****************************
c
      subroutine convertype (buf, nbuf, btypin, ftypin)
c
c        Converts real or integer words between CRAY and IEEE
c        binary machine representations.
c
c        buf   = buffer containing real or integer data,
c                always assumed real* or integer*8 (modified)
c        nbuf  = number of reals or integers (supplied)
c        btypin= 'r' for real, 'i' for integers (supplied)
c        ftype = 'CRAY' or 'IEEE', the machine-representation type
c                of buf values on input (supplied)
c
c----------------------------------------------------------------------
#include <PARAMS.h>
c----------------------------------------------------------------------
 
      dimension buf(nbuf)
      character*(*) btypin, ftypin
      character btype*1, ftype*4
 
      character cdum8*8, cdum4*4
      integer*4 idum4, nbuf4, ione4
      integer*8 idum8
      equivalence (cdum8,zdum8,idum8)

#if defined ( linux )
c do nothing for linux cluster machines (use convert='cray' in openf)
      return
#endif
 
      btype = btypin
      call upper (btype)
      if (btype.eq.'C') return
 
      ftype = ftypin
      call upper (ftype)
 
      nbuf4 = nbuf
      ione4 = 1
 
#ifdef cray
 
      if (ftype.eq.'IEEE') then
 
        if (btype.eq.'R') then
 
          do 5 i=1,nbuf
            ztmp = buf(i)
#if defined ( aix )
            call dptocf_ (ztmp, buf(i), ione4)
#else
            call dptocf (ztmp, buf(i), ione4)
#endif
    5     continue
 
        else if (btype.eq.'I') then
 
c         Must do the following since we only have conversion routines
c         for IEEE 4-byte integers
          do 10 i=1,nbuf
            zdum8 = buf(i)
            cdum4 = cdum8(5:8)
#if defined ( aix )
            call sptoci_ (cdum4, buf(i), ione4)
#else
            call sptoci (cdum4, buf(i), ione4)
#endif
   10     continue
 
        endif
 
      endif
 
#else
 
      if (ftype.eq.'CRAY') then
 
        if (btype.eq.'R') then
 
          do 15 i=1,nbuf
            ztmp = buf(i)
#if defined ( aix )
            call ctowsd_ (ztmp, buf(i), ione4)
#else
            call ctowsd (ztmp, buf(i), ione4)
#endif
   15     continue
c         call ctowsd (buf, buf, nbuf4)
 
        else if (btype.eq.'I') then
 
c         Must do the following since we only have conversion routines
c         for IEEE 4-byte integers
          do 20 i=1,nbuf
#if defined ( aix )
            call ctowsi_ (buf(i), idum4, ione4)
#else
            call ctowsi (buf(i), idum4, ione4)
#endif
            idum8 = idum4
            buf(i) = zdum8
   20     continue
 
        endif
 
      endif
 
#endif
 
      return
      end
