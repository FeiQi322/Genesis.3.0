c -- $Id: util.F,v 1.3 2001/07/03 21:29:21 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/util.F,v $
c -- $Name:  $

C
C ***************************** BANNER ********************************
C
      SUBROUTINE BANNER
C----------------------------------------------------------------------
C
CL            WRITES MODEL BANNER TO UNIT NOUT
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMHED.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------

#if   defined ( aix )
      character*8 idate,date
      character*8 itime,clock_
#elif defined ( sgi )
      CHARACTER*9 IDATE
      CHARACTER*8 ITIME
#elif defined ( sun )
      CHARACTER*24 JDATE
#elif defined ( linux )
      CHARACTER*8  KDATE
      CHARACTER*10 KTIME
#endif

C----------------------------------------------------------------------
 
c     WRITE(NOUT,'("1")')
 
#if   defined ( aix )
      idate = date()
      itime = clock_()
#elif defined ( cray )
      CALL DATE  (IDATE)
      CALL CLOCK (ITIME)
#elif defined ( sgi )
      CALL DATE (IDATE)
      CALL TIME (ITIME)
#elif defined ( sun )
      CALL FDATE (JDATE)
#elif defined ( linux )
      CALL DATE_AND_TIME (KDATE, KTIME)
#endif

#if defined ( cray ) || defined ( sgi ) || defined ( aix )
      WRITE(NOUT,100) VERSGEN, IDATE, ITIME
 100  FORMAT(
     *  /17X,'=============================================='
     *  /17X,'| GENESIS CLIMATE MODEL    VERSION: ',  A8,' |'
     *  /17X,'| DATE: ',A9,   '          TIME:    ',  A8,' |'
     *  )
#elif defined ( sun ) 
      WRITE(NOUT,100) VERSGEN, JDATE
 100  FORMAT(
     *  /17X,'=============================================='
     *  /17X,'| GENESIS CLIMATE MODEL    VERSION: ',  A8,' |'
     *  /17X,'|          ',                 A24,'          |'
     *  )
#elif defined ( linux )
      WRITE(NOUT,100) VERSGEN, KDATE, KTIME
 100  FORMAT(
     *  /17X,'=============================================='
     *  /17X,'| GENESIS CLIMATE MODEL    VERSION: ',  A8,' |'
     *  /17X,'| DATE: ',A8,  '           TIME:  ',   A10,' |'
     *  )
#else
      WRITE(NOUT,100) VERSGEN
 100  FORMAT(
     *  /17X,'=============================================='
     *  /17X,'| GENESIS CLIMATE MODEL    VERSION: ',  A8,' |'
     *  )
#endif

#ifdef omp 
      write (nout,"(
     *   17x,'| with omp multiprocessing                   |' )")
#endif

#ifdef isotrac
#ifdef isotracreg
      write (nout,"(
     *   17x,'| with water isotopes and source regions     |' )")
#else
      write (nout,"(
     *   17x,'| with water isotopes                        |' )")
#endif
#endif

#ifdef ibis 
#ifdef ibisint
      write (nout,"(
     *   17x,'| with interactive ibis vegetation           |' )")
#else
      write (nout,"(
     *   17x,'| with ibis vegetation                       |' )")
#endif
#endif

#ifdef cosmog 
      write (nout,"(
     *   17x,'| with cosmogenic tracers                    |' )")
#endif

#ifdef hydmodel 
      write (nout,"(
     *   17x,'| with hydrology model coupling              |' )")
#endif

#ifdef dotdot
      write (nout,"(
     *   17x,'| with dotdot monitoring                     |' )")
#endif

      WRITE(NOUT,200) 
 200  FORMAT(
     *   17X,'=============================================='
     *  /)
 
      RETURN
      END
C
C ***************************** CALDYI ********************************
C
      SUBROUTINE CALDYI (JYEAR, JMNTH, JDAY, CALDAY, CALSEC, ISTEPY,
     *                   MCDATE, MCSEC, MBDATE, DTIME)
C----------------------------------------------------------------------
C
CL            COMPUTE CALENDAR DAY INFO (JYEAR,JMNTH,JDAY,CALDAY,CALSEC,
CL            ISTEPY) FROM CURRENT DATE (MCDATE, IN FORM YYMMDD) AND
CL            SECONDS INTO CURRENT DAY (MCSEC). ALSO USES BASE DATE
CL            (MBDATE, IN FORM YYMMDD) AND TIMESTEP (DTIME, IN SECONDS).
c
C----------------------------------------------------------------------
#include <PARAMS.h>
C --------------------------------------------------------------------
      DIMENSION DAZ(12)
      DATA DAZ /  0.,  31.,  59.,  90., 120., 151.,
     *          181., 212., 243., 273., 304., 334. /
C --------------------------------------------------------------------
C
CL
CL            MCDATE IS IN FORM YYMMDD. IF >=0, YY MEASURED FROM 0 A.D.
CL            IF < 0, MEASURED FROM 1950 A.D. ("B.P.").
CL            ALLOW FOR NEGATIVE MCDATE BY USING ITS ABSOLUTE VALUE
CL            FOR JMNTH, JDAY.
C
      JYEAR = MCDATE/10000
      JMNTH = MOD(IABS(MCDATE)/100,100)
      JDAY  = MOD(IABS(MCDATE),100)
C
CL            CALCULATE SECONDS INTO CURRENT YEAR, FOR HISTORY TIMING
CL            (VIA ISTEPY)
C
      CALSEC = (DAZ(JMNTH) + JDAY - 1)*86400. + MCSEC
C
CL            CALCULATE DAY OF THE YEAR, FOR UPDATING PRESCRIBED
CL            BOUNDARY CONDS. CALDAY RUNS FROM 1.0 ON JAN 1 AT 12AM
CL            TO 366.0 ON DEC 31 AT 12PM. (IGNORES LEAP YEARS SINCE
CL            DATASETS DO TOO; 2/29 WOULD HAVE SAME CALDAY AS 3/1.)
C
      CALDAY = DAZ(JMNTH) + JDAY + MCSEC/86400.
C
C
CL            SET ISTEPY, THE TIMESTEP WITHIN CURRENT 365-DAY YEAR
C
      ISTEPY = NINT(CALSEC/DTIME)
 
c     If at end of a year, convert istepy to correspond to end of
c     day 365 rather than beginning of day 0, to be consistent with
c     conventions preventing 0 values for nhist,nsave,nnexf in dattimes
c     (also assumed in iceprint in dynice.f)
 
      IF (ISTEPY.EQ.0) ISTEPY = NINT(365*86400./DTIME)

      RETURN
      END
C
C ***************************** DATECONV *******************************
C
      subroutine dateconv (ival, index, clabel, kount)
c-----------------------------------------------------------------------
c
c             Checks if ival contains string "mm/dd" where mm is month #
c             (1-12) and dd is day # (1-31). If so, change ival to day #
c             within year (1-365, no leap years). First tests if any
c             byte holds a non-printing character, in which case ival
c             is assumed to be a binary integer and is unchanged.
c
c             index  = namelist index # of ival (for error message only)
c             clabel = namelist variable name (for error message only)
c             kount  = counter, incremented by 1 if ival contains string
c
c-----------------------------------------------------------------------
#include <PARAMS.h>
c-----------------------------------------------------------------------
      character clabel*(*), cval*8, cform*8
      dimension idaz(12)
      data idaz /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
c
c
      write(cval,'(a8)') ival
c
c        If any byte is a non-printing character, assume ival is
c        meant to be an integer value and return
c
      do 5 i=1,8
        if (ichar(cval(i:i)).lt.32 .or. ichar(cval(i:i)).gt.126) return
    5 continue
c
c        Presumably ival contains a string, so increment counter
c
      kount = kount + 1
c
c        is = locn of "/", il = locn of last non-blank
c
      is = 0
      il = 0
      do 10 i=1,8
        if (cval(i:i).eq.'/') is = i
        if (cval(i:i).ne.' ') il = i
   10 continue
      if (is.le.1 .or. is.eq.il) goto 900
c
c        read integer values before and after "/"
c
      write (cform,20) is-1
      read (cval(1:is-1), cform, err=900) imon
      write (cform,20) il-is
      read (cval(is+1:il), cform, err=900) iday
   20 format('(i',i1,')')
c
      if (imon.lt.1.or.imon.gt.12 .or. iday.lt.1.or.iday.gt.31) goto 900
c
c        set day # within year, and return
c
      ival = idaz(imon) + iday
      return
c
  900 write(6,9000) clabel, index, cval
 9000 format(/' *** Error in "month/day" input for parameter ',a
     *       /'     index = ',i4
     *       /'     input = ',a)
      call endrun (-1)
c
      end
C
C ***************************** DATECONV2 ******************************
C
      subroutine dateconv2 (ival, cval, clabel)
c-----------------------------------------------------------------------
c
c             Checks if cval contains string "yy/mm/dd" where yy is
c             year#, mm is month# (1-12) and dd is day# (1-31). If so,
c             sets ival to integer value YYMMDD. If not, sets ival to
c             integer value if cval contains a valid integer, or detects
c             error if not.
c
c             ival   = integer date or day # (returned)
c             cval   = character string (supplied)
c             clabel = namelist variable name (supplied)
c
c-----------------------------------------------------------------------
#include <PARAMS.h>
c-----------------------------------------------------------------------
      character cval*(*), clabel*(*), cform*8
      dimension idaz(12)
      data idaz /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
c
c        isa = locn of first "/", isb = locn of second "/",
c        iend = locn of last non-blank.
c
      isa = 0
      isb = 0
      iend= 0
      do 10 i=1,len(cval)
        if (cval(i:i).eq.'/') then
          if (isa.eq.0) then
            isa = i
          else if (isb.eq.0) then
            isb = i
          endif
        endif
        if (cval(i:i).ne.' ') iend = i
   10 continue
 
      if (iend.eq.0 .or. isa.eq.1 .or. isa.eq.iend .or. isb.eq.iend
     *    .or. (isa.gt.1 .and. isb.lt.isa+2) ) goto 900
 
c        If no "/", try reading as a single integer.
c        Else, read integer values before first "/", between first
c        two "/"s, and after the second "/".
 
      if (isa.eq.0) then
        write (cform,20) iend
   20   format('(i',i2.2,')')
        read (cval(1:iend), cform, err=900) ival
c       if integer (+days, -timesteps) entered for enddate, return
        if (clabel.eq.'ENDDATE' .and. iabs(ival).lt.100000) return
        iyear = ival/10000
        imon  = mod(iabs(iyear)/100,100)
        iday  = mod(iabs(iyear),100)
      else
        write (cform,20) isa-1
        read (cval(1:isa-1), cform, err=900) iyear
 
        write (cform,20) isb-isa-1
        read (cval(isa+1:isb-1), cform, err=900) imon
 
        write (cform,20) iend-isb
        read (cval(isb+1:iend), cform, err=900) iday
      endif
 
      if (imon.lt.1.or.imon.gt.12.or.iday.lt.1.or.iday.gt.31) goto 900
 
cv2.3 shift years (if positive) to 1900 A.D + iyear.
cv2.3 if (iyear.ge.0) iyear = 1900 + iyear
 
c        Set YYMMDD integer value and return
 
      if (iyear.ge.0) then
        ival = 10000*iyear + 100*imon + iday
      else
        ival = 10000*iyear - 100*imon - iday
      endif
c
      return
c
  900 write(6,9000) clabel, cval
 9000 format(/' *** Error in "year/month/day" or "day #" input for',
     *        ' parameter ',a
     *       /'     input = ',a)
      call endrun (-1)
c
      end
C
C ***************************** DATEDIF *******************************
C
      subroutine datedif (ndays, mdatea, mdateb)
c-----------------------------------------------------------------------
c
c             Computes number of days ndays between two (YYMMDD) dates
c             mdatea and mdateb, from 00:00 GMT on mdatea to 24:00 GMT
c             on mdateb.
c
C----------------------------------------------------------------------
      dimension idaz(12)
      save idaz
      data idaz /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
 
      jya = mdatea/10000
      jyb = mdateb/10000
      jma = mod(iabs(mdatea)/100,100)
      jmb = mod(iabs(mdateb)/100,100)
      jda = mod(iabs(mdatea),100)
      jdb = mod(iabs(mdateb),100)
 
      ndays =  365*jyb + idaz(jmb) + jdb
     *      - (365*jya + idaz(jma) + jda) + 1
 
      return
      end
C
C ***************************** ENDRUN ********************************
C
      SUBROUTINE ENDRUN (IFLAG)
C ---------------------------------------------------------------------
C
CL            CLOSE LOCAL FILES AND ABORT THE MODEL.
CL            ALL STOPPING POINTS GO THROUGH HERE...THERE ARE 
CL            NO OTHER STOP STATEMENST IN PROGRAM. IFLAG = 0 
CL            ONLY FOR NORMAL END FROM MAIN PROGRAM, SO $STATUS 
CL            (CALL EXIT) CAN BE USED BY RUN SCRIPT. 
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      character*200 cmd

C        Close always-open or possibly-open files

#if ! defined (asyncmode)
      CALL CLOSEF (NSSTICE)
      CALL CLOSEF (NDATA)
      CALL CLOSEF (NDATA_A)
      CALL CLOSEF (IUREG)
      CALL CLOSEF (IUREG_A)
      CALL CLOSEF (IUSTORM)
#endif

c-------------------------
      IF (IFLAG.NE.0) THEN
c-------------------------
#if defined ( sgi ) || defined ( sun )
        write (nout,'(40("="))')
        ier = flush(nout)
        call gerror (cmd)
c       call perror (cmd)  ! alternate, sun
        write (nout,'("gerror: ",a)') cmd(1:lenchr(cmd))
        write (nout,'(40("="))')
#endif

        WRITE(NOUT,900)
  900   FORMAT(/' *** PROGRAM ABORT ***'/)
 
c         Cray utility to print traceback to stdout
 
#ifdef cray
        CALL TRACEBK()
#endif

#if defined (linux)
        STOP 1
#else
c       CALL EXIT (IFLAG)
        CALL EXIT (-1)
#endif
c---------
      ELSE
c---------
        STOP
c----------
      ENDIF
c----------
 
      END
C
C ***************************** IVAR **********************************
C
      SUBROUTINE IVAR (KNAME,KVALUE)
C----------------------------------------------------------------------
C
CL            PRINT NAME AND VALUE OF INTEGER VARIABLE
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
C
      WRITE(NOUT,900) KNAME,KVALUE
  900 FORMAT(1H ,A8,2H =,I12)
C
      RETURN
      END
c
c ***************************** lenchr ********************************
c
      integer function lenchr (ch)
c     returns position of right-most non-blank,non-null character in ch.
      character*(*) ch
      do i = len(ch), 1, -1
        lenchr = i
        if (ch(i:i).ne.' ' .and. ch(i:i).ne.char(0)) return
      enddo
      lenchr = 1
      return
      end
C
C ***************************** LSTRIP ********************************
C
      subroutine lstrip (cin)
 
c        Strips leading blanks from character string cin
 
      character*(*) cin
      character*1000 ctmp
 
      do 10 i=1,len(cin)
         if (cin(i:i).ne.' ') then
            ctmp = cin(i:)
            cin = ctmp
            return
         endif
   10 continue
 
      return
      end
C
C ***************************** MAPOUT ********************************
C
      subroutine mapout (iu, iarr, nlon, nlat, alon, alat, maskval,
     *                   ctitle)
c
c        Prints 2-D global map in tabular form to unit iu. Depending on
c        range of iarr values, prints in '120a1' or in '24i5' format.
c        If can print as '120a1' and nlon le 120, prints one global map,
c        else prints as two separate hemispheres.
c
c        iu      = unit number for output
c        iarr    = data array
c        nlon    = longitude dimension (east from Greenwich)
c        nlat    = latitude dimension (south to north)
c        alon    = longitudes (radians)
c        alat    = latitudes (radians)
c        maskval = iarr value for which blanks are shown
c        ctitle  = title
 
c-----
c     On sgi, declare nint and int external, and supply own versions
c     in utilsx, to avoid sgi bug with -O3 and -i8.
#ifdef sgi
      external nint, int
#endif

#if defined (NOFLUSH) || defined (OWNFLUSH)
      external flush
      integer flush
#endif
c-----
 
      dimension iarr(nlon,nlat), alon(nlon), alat(nlat)
      character*(*) ctitle
      character ch(-30:35)*1
      character carr(360,360)*5    ! max resolution
      dimension alono(360)
      parameter (pi=3.14159265358979, dpr=(180./pi))
      save ch
      data ch/
     *        '#','*','=','-',
     *        'Z','Y','X','W','V','U','T','S','R','Q',
     *        'P','O','N','M','L','K','J','I','H','G',
     *        'F','E','D','C','B','A',
     *        '0','1','2','3','4','5','6','7','8','9',
     *        'a','b','c','d','e','f','g','h','i','j',
     *        'k','l','m','n','o','p','q','r','s','t',
     *        'u','v','w','x','y','z'/
 
 
c        Check hard-coded dimension limits
 
      if (nlon.gt.360 .or. nlat.gt.360) then
        write(iu,10) ctitle(1:lenchr(ctitle)), nlon, nlat, 360, 360
   10   format(/' *** Cannot print out map for field: ',a
     *         /'     Requested longitude, latitude =',2i6
     *         /'     Hard-coded lon,lat limits     =',2i6
     *         /'     To fix, increase sizes of carr and alono in',
     *          ' subroutine mapout.')
        return
      endif
 
c        Set output longitude scale (shifted so western hemis is first)
 
      do 20 i=1,nlon
        ii = mod (i+nlon/2-1, nlon) + 1
        alono(ii) = alon(i)
        if (alono(ii).gt.pi) alono(ii) = alono(ii) - 2.*pi
   20 continue
 
c        Find max and min data values
 
      minval =  100000000
      maxval = -100000000
      do 50 j=1,nlat
        do 52 i=1,nlon
           if (iarr(i,j).ne.maskval) then
              minval = min (iarr(i,j), minval)
              maxval = max (iarr(i,j), maxval)
           endif
   52   continue
   50 continue
 
c        Decide if can use single chars (nwidth=1), or must use ints
c        (nwidth=5). Also transfer iarr values to character array carr,
c        shifted so western hemis is first.
c
 
      if (minval.ge.-30 .and. maxval.le.35) then
 
        nwidth = 1
        do 100 j=1,nlat
          do 102 i=1,nlon
            ii = mod (i+nlon/2-1, nlon) + 1
            if (iarr(i,j).eq.maskval) then
              carr(ii,j) = ' '
            else
              carr(ii,j) = ch(iarr(i,j))
            endif
  102     continue
  100   continue
 
      else
 
        nwidth = 5
        do 110 j=1,nlat
          do 112 i=1,nlon
            ii = mod (i+nlon/2-1, nlon) + 1
            if (iarr(i,j).eq.maskval) then
              carr(ii,j) = ' '
            else
c             Offset used for land topography:
              if (iarr(i,j).gt.1000000/2) then
                write(carr(ii,j),'(i5)')
     *            min (9999, max(-9999, iarr(i,j)-1000000))
                do 114 m=1,5
                  if (carr(ii,j)(m:m).eq.' ') carr(ii,j)(m:m) = '.'
  114           continue
              else
                write(carr(ii,j),'(i5)')
     *            min (9999, max(-9999, iarr(i,j)))
              endif
            endif
  112     continue
  110   continue
 
      endif
 
c        Decide whether to print one global map (nhem=1), or as two
c        separate hemispheres (nhem=2)
 
      if (nwidth.eq.1 .and. nlon.le.120) then
        nhem = 1
      else
        nhem = 2
      endif
 
      do 200 ihem = 1,nhem
 
        if (nhem.eq.2  .and. ihem.eq.1) then
          i1 = 1
          i2 = nlon/2
          inc = (i2-i1) / (120/nwidth) + 1
          write(iu,210) char(12), ctitle(1:lenchr(ctitle)),
     *                  '(western hemisphere)',inc
  210     format(/a1,' ',a,'  ',a,'  (longitude increment =',i2,')')
        else if (nhem.eq.2 .and. ihem.eq.2) then
          i1 = nlon/2 + 1
          i2 = nlon
          inc = (i2-i1) / (120/nwidth) + 1
          write(iu,210) char(12), ctitle(1:lenchr(ctitle)),
     *                  '(eastern hemisphere)', inc
        else if (nhem.eq.1) then
          i1 = 1
          i2 = nlon
          inc = (i2-i1) / (120/nwidth) + 1
          write(iu,210) char(12), ctitle(1:lenchr(ctitle)),
     *                  ' ', inc
        endif
 
 
        if (nwidth.eq.1) then
 
          write(iu,220) (nint(alono(i)*dpr), i=i1,i2,inc*5)
  220     format(/5x,24i5)
          write(iu,*)
          do 225 j=nlat,1,-1
            write(iu,226) alat(j)*dpr, (carr(i,j)(1:1),i=i1,i2,inc)
  226       format(1x,f5.1,3x,120a1)
  225     continue
 
        else if (nwidth.eq.5) then
 
          write(iu,230) (nint(alono(i)*dpr), i=i1,i2,inc)
  230     format(/9x,24i5)
          write(iu,*)
          do 235 j=nlat,1,-1
            write(iu,236) alat(j)*dpr, (carr(i,j),i=i1,i2,inc)
  236       format(1x,f5.1,3x,24a5)
  235     continue
 
        endif
 
  200 continue
 
      write(iu,*)
 
c        flush file buffer
 
      ier = flush(iu)
 
      return
      end
C
C ***************************** MESAGE ********************************
C
      SUBROUTINE MESAGE (KMESS)
C----------------------------------------------------------------------
C
CL            PRINT 40-CHARACTER MESSAGE ON OUTPUT CHANNEL NOUT
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      DIMENSION KMESS(5)
C----------------------------------------------------------------------
C
      WRITE(NOUT,900) KMESS
  900 FORMAT(1X,5A8)
C
      RETURN
      END
C
C ***************************** NUTIME ********************************
C
      subroutine nutime (mdate, msec, mbdate, nstep, dtime, iflag)

c        Computes date mdate (yymmdd) and seconds msec
c        corresponding to seconds nstep*dtime from 00:00 gmt
c        on base date mbdate.
c        if iflag=0, returns prev day + 86400 sec if at midnight,
c        if iflag=1, returns next day +     0 sec if at midnight.

c        Use int*8 for local variable ntim (secs since basedate), 
c        otherwise it would exceed int*4 range (~2^31, ~few*10^9) 
c        for GCM runs longer than ~68 years from basedate. Also
c        need temporary variable n8 for calcs with ntim.


      integer*8 ntim, n8

      dimension imday(12)
      save imday
      data imday /31,28,31,30,31,30,31,31,30,31,30,31/

c     on sgi, declare nint and int external, and supply own versions
c     in utilsx, to avoid sgi bug with -o3 and -i8.
#ifdef sgi
      external nint, int
#endif

c        Calculate elapsed time in seconds

      n8 = nstep
      ntim = n8*nint(dtime)

c        Compute current date, by first initializing to base date
c        mbdate, and then stepping out month by month.
c        allow for negative mcdate by using its absolute value
c        and converting sign of mdate at end.

      iyr = mbdate/10000
      imo = mod(iabs(mbdate)/100,100)
      idy = mod(iabs(mbdate),100)

c     next two lines would allow for leap years (nb: 2000 is a leap
c     year by 400-yr rule) - but omit leap years from current model,
c     so comment out here and within 200 loop below (also see zencal).
cleap imday(2) = 28
cleap if (mod(iyr,4).eq.0) imday(2) = 29
c
c     if iflag = 0, subtract 1 sec so that midnight yields previous day
      if (iflag.eq.0) then
        idy = idy + (ntim-1)/86400
      else
        idy = idy + (ntim  )/86400
      endif
200   if (idy .gt. imday(imo)) then
         idy = idy - imday(imo)
         imo = imo + 1
210      if (imo .gt. 12) then
            iyr = iyr + 1
cleap       imday(2) = 28
cleap       if (mod(iyr,4).eq.0) imday(2) = 29
            imo = imo - 12
            go to 210
         endif
         go to 200
      endif
      mdate = isign  (iabs(iyr)*10000 + imo*100 + idy, iyr)

c         Compute seconds into current day. If iflag = 0, subtract 
c         1 sec so that midnight yields msec = 86400

      if (iflag.eq.0) then
c       msec = mod (ntim-1, 86400) + 1   !mod not generic w.int*8,IBM SP
        n8 = (ntim-1)/86400
        msec = (ntim-1) - n8*86400 + 1
      else
c       msec = mod (ntim, 86400)         !mod not generic w.int*8,IBM SP
        n8 = ntim/86400
        msec = ntim - n8*86400
      endif

      return
      end

C
C ***************************** SHOWTIME ******************************
C
      SUBROUTINE SHOWTIME
C----------------------------------------------------------------------
C
CL            PRINT OUT CURRENT TIME INFORMATION AT REGULAR INTERVALS
CL            (CALLED AT END OF TIMESTEP FROM LINEMS)
C
C----------------------------------------------------------------------
#include <PARAMS.h>
C----------------------------------------------------------------------
#include <COMDYR.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
#include <COMTIM.h>
C----------------------------------------------------------------------
      CHARACTER CTMP*5
C---------------------------------------------------------------------
 
c           Write out current-time info once a day
 
 888  N = MAX (NINT(1.*86400./DTIME), 1)
c888  N = 1
 
      IF (MOD(NSTEP,N).EQ.0) THEN
 
c           Display "A.D." if mbdate<0, "B.P." year if mbdate<0
 
         IF (MBDATE.GE.0) THEN
            IYR = JYEAR
            CTMP = ' A.D.'
         ELSE
            IYR = -JYEAR
            CTMP = ' B.P.'
         ENDIF
 
         WRITE(NOUT,802)
     *     NSTEP, NSTEP*DTIME/86400., ISTEPY*DTIME/86400.,
     *     IYR, JMNTH, JDAY, CTMP, MCSEC
  802    FORMAT(
     *   /' Timestep=',i8,
     *    '   Elapsed days=',f10.3,'   Calendar days=',f8.3,
     *    '   Date=',i8,'/',i2.2,'/',i2.2,a,' + ',i5,' secs')
 
c           If every timestep, presumably want to see immediately
 
         IF (N.EQ.1) IER = FLUSH(NOUT)
 
      ENDIF
 
      RETURN
      END
C
C ***************************** READDATA ******************************
C
      SUBROUTINE READDATA (IU, ARR, MLON, MLAT, CKEY, ICODE, INUM)
C ---------------------------------------------------------------------
C
C        Reads an agcm or lsx data input file, checking on label and
C        dimensions in header record. Data files are ascii, with 1
C        header record, 3 skip records (blank, longitude grid, blank),
C        and MLAT data records. Each data record contains a
C        latitude grid value, followed by MLON data values either
C        in A1 or F5.0 format (flagged by ICODE). The longitudes are
C        western hemisphere first, then eastern, ie, (I=MLON/2+1,MLON),
C        (I=1,MLON/2). The data records are in reverse latitude-grid
C        order (J=MLAT,1,-1). Data values are returned in *real* array
C        ARR.
C
C        Permissible characters in A1 format are 0-9, blank(=0),
C        A-Z (=10-35), or a-z(=10-35).
C
C        If file has more than one field, each set of data records is
C        preceded by the same 4 records as above. Each group is
C        separated by an extra blank line, between the last data record
C        of the preceding field and the header record of the next field.
C
C        IU    = unit number
C        ARR   = data array (returned)
C        MLON  = longitude dimension
C        MLAT  = latitude  dimension
C        CKEY  = expected 8-byte label in word # 1, rec # 1
C        ICODE = 1 for A1 format, 2 for F5.0 format.
C        INUM  = field number for this file
C
C ---------------------------------------------------------------------
#include <PARAMS.h>
C ---------------------------------------------------------------------
#include <COMLUN.h>
C ---------------------------------------------------------------------
      DIMENSION ARR(MLON,MLAT)
      CHARACTER CKEY*(*)
      CHARACTER CARR(360,360)*1, CKEYZ*8, COMMENT*100
C ---------------------------------------------------------------------
C
C        Find index in program's file list (just for messages)
C
#if ! defined (asyncmode)
      CALL FINDIU (IU, LOC)
#endif
C
C        Initialize record number (just for error messages)
C
      IREC = (INUM-1)*(4+MLAT) + MAX(INUM-2,0)
C
C        If past first field, skip blank line separator
C
 
      IF (INUM.GT.1) THEN
        IREC = IREC + 1
        READ(IU, *, ERR=900, IOSTAT=IOS)
      ENDIF
C
C        Read header record
C
      IREC = IREC + 1
      READ(IU, 10, ERR=900, IOSTAT=IOS) CKEYZ, MLONZ, MLATZ, COMMENT
   10 FORMAT(A8,2I8,8X,A)

      IF (CKEYZ.NE.CKEY .OR. MLONZ.NE.MLON .OR. MLATZ.NE.MLAT) THEN
         WRITE(NOUT,12) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC))),
     *                  COMMENT(1:LENCHR(COMMENT)),
#if defined (sun) || defined (linux)
     *                  CKEYZ,MLONZ,MLATZ, CKEY,MLON,MLAT
#else
     *                  CKEYZ,MLONZ,MLATZ, CKEY//'        ',MLON,MLAT
#endif

   12    FORMAT(/' *** Fatal Error.',
     *           ' Incorrect attributes for data file ',A
     *          /'     File comment                        = ',A,
     *          /'     File     label, longitude, latitude = ',A8, 2I6
     *          /'     Expected label, longitude, latitude = ',A8, 2I6/)
         CALL ENDRUN (-1)
      ELSE
         WRITE(NOUT,14) COMMENT(1:LENCHR(COMMENT)), MLONZ, MLATZ
   14    FORMAT('    File comment   : ',a
     *         /'    File resolution:',2I6)
      ENDIF
C
C        Skip 3 records (blank, longitude grid, blank)
C
      IREC = IREC + 3
      READ(IU, 50, ERR=900, IOSTAT=IOS)
   50 FORMAT(//)

c=========================
      IF (ICODE.EQ.1) THEN
c=========================

C          Read "A1-format" data records
C
        DO 100 J=MLAT,1,-1
          IREC = IREC + 1
          READ (IU,102,ERR=900,IOSTAT=IOS) (CARR(I,J),I=MLON/2+1,MLON),
     *                                     (CARR(I,J),I=1,MLON/2)
  102     FORMAT(5X,3X,360A1)
          DO 104 I=1,MLON
            IC = ICHAR(CARR(I,J)(1:1))
            IF (IC.EQ.ICHAR(' ')) THEN
              ARR(I,J) = 0
            ELSE IF (IC.GE.ICHAR('0') .AND. IC.LE.ICHAR('9')) THEN
              ARR(I,J) = IC - ICHAR('0')
            ELSE IF (IC.GE.ICHAR('A') .AND. IC.LE.ICHAR('Z')) THEN
              ARR(I,J) = IC - ICHAR('A') + 10
            ELSE IF (IC.GE.ICHAR('a') .AND. IC.LE.ICHAR('z')) THEN
              ARR(I,J) = IC - ICHAR('a') + 10
            ELSE
              GOTO 910
            ENDIF
  104     CONTINUE
  100   CONTINUE

#ifdef isotrac
c==============================
      ELSE IF (ICODE.EQ.3) THEN
c==============================

c          Read "i2-format" data records

        DO 120 J=MLAT,1,-1
          IREC = IREC + 1
          READ (IU,122,ERR=900,IOSTAT=IOS) (ARR(I,J),I=MLON/2+1,MLON),
     *                                     (ARR(I,J),I=1,MLON/2)
  122     FORMAT(5X,2X,360f2.0)
  120   CONTINUE
#endif

c=========
      ELSE
c=========

C          Read "F5.0-format" data records

        DO 150 J=MLAT,1,-1
          IREC = IREC + 1
          READ (IU,152,ERR=900,IOSTAT=IOS) (ARR(I,J),I=MLON/2+1,MLON),
     *                                     (ARR(I,J),I=1,MLON/2)
  152     FORMAT(5X,3X,360F5.0)
  150   CONTINUE

c==========
      ENDIF
c==========

      RETURN
C
C        Fatal errors
C
  900 WRITE(NOUT,902) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC))),
     *              COMMENT, IREC, IOS
  902 FORMAT(/' *** Error reading data file ',A
     *       /'     File comment = ',A
     *       /'     record number=',I6
     *       /'     Iostat code  =',I6/)
      CALL ENDRUN (-1)
C
  910 WRITE(NOUT,912) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC))),
     *                COMMENT, IREC, I, J, CARR(I,J)
  912 FORMAT(/' *** Bad character in data file ',A
     *       /'     File comment    = ',A
     *       /'     record number   =',I6
     *       /'     longitude index =',I6
     *       /'     latitude  index =',I6
     *       /'     bad character   = ',A1/)
      CALL ENDRUN (-1)
C
      END
C
C ***************************** STARLIST ******************************
C
      subroutine starlist (clist, nlist, cwhich)
c-----------------------------------------------------------------------
c
c        Expands '*' for a filename list, if just one filename
c        entered (must be first one). Replaces * (or *'s) with
c        '001','002',...'999' surrounded by any entered text.
c
c             clist  = array of filenames
c             nlist  = dimension of clist
c             cwhich = list identifier for message only
c
c-----------------------------------------------------------------------
#include <PARAMS.h>
c-----------------------------------------------------------------------
#include <COMLUN.h>
c-----------------------------------------------------------------------
      character clist(nlist)*(*), cwhich*(*)
      character corig*16, ctmp*4
c
c        Check if first entry is non-blank, rest blank. If not, return.
c
      do 10 n=nlist,1,-1
        if (clist(n).ne.' ') goto 12
   10 continue
      return
   12 if (n.ne.1) return
c
      corig = clist(1)
      lenc = len(corig)
c
c        Look for first '*'
c
      do 20 i=1,lenc
        ia = i
        if (corig(i:i).eq.'*') goto 22
   20 continue
      return
   22 continue
c
c        Look for last '*' (if last entered char, corig(ib+1:lenc)=' ')
c
      do 30 i=lenc,1,-1
        ib = i
        if (clist(1)(i:i).eq.'*') goto 32
   30 continue
      return
   32 continue
c
c        Expand list
c
      do 100 n=1,nlist
        write(ctmp,'(i4)') 1000 + n
        if (ia.gt.1) then
           clist(n) = corig(1:ia-1) // ctmp(2:4) // corig(ib+1:lenc)
        else
           clist(n) =                  ctmp(2:4) // corig(ib+1:lenc)
        endif
  100 continue
c
      nz = min (nlist,999)
      write(nout,900) cwhich, corig(1:lenchr(corig)),
     *               (clist(n )(1:lenchr(clist(n ))), n=1,3),
     *                clist(nz)(1:lenchr(clist(nz)))
  900 format(/' ',a,' filename list has a single entry ',a,
     *        ' containing * character(s).'
     *       /' List is expanded into:',3(2x,a),'  ...  ',a /)
c
      return
      end
C
C ***************************** UPPER **********************************
C
      subroutine upper (cin)
 
c        Upper-cases character string cin
 
      character*(*) cin
 
      ia = ichar('a')
      iz = ichar('z')
      ishif = ichar('A') - ichar('a')
 
      do 10 i=1,len(cin)
        if (ichar(cin(i:i)).ge.ia .and. ichar(cin(i:i)).le.iz)
     *    cin(i:i) = char (ichar(cin(i:i)) + ishif)
   10 continue
 
      return
      end
C
C ***************************** USERDEF0 *******************************
C
      subroutine userdef0 (nstep, nrstrt, nstop, istepy, dtime, iu)
 
c User-codeable routine, called from end of linems.
 
c nstep  = number of timesteps (from base date)
c nrstrt = value of nstep at first timestep of this run
c nstop  = value of nstep at last  timestep of this run
c istepy = number of timesteps within current calendar year
c dtime  = timestep (sec)
c iu     = logical unit number
 
c-----------------------------------------------------------------------
#include <PARAMS.h>
c-----------------------------------------------------------------------
#include <COMCON.h>
c-----------------------------------------------------------------------
#include <COMLUN.h>
c-----------------------------------------------------------------------
#include <COMMAP.h>
c-----------------------------------------------------------------------
#include <COMMUN.h>
c-----------------------------------------------------------------------
 
      return
      end
C
C ***************************** USERDEF1 *******************************
C
      subroutine userdef1 (nstep, nrstrt, nstop, istepy, dtime, iu)
 
c User-codeable routine, called from end of linems.
c Accumulates and prints daily mean precipitation at certain points.
 
c nstep  = number of timesteps (from base date)
c nrstrt = value of nstep at first timestep of this run
c nstop  = value of nstep at last  timestep of this run
c istepy = number of timesteps within current calendar year
c dtime  = timestep (sec)
c iu     = logical unit number
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMLUN.h>
c-----------------------------------------------------------------------
#include <COMMAP.h>
c-----------------------------------------------------------------------
#include <COMMUN.h>
c ---------------------------------------------------------------------
      dimension
     *   lonr15(13), latr15(13),
     *   lont31(17), latt31(17),
     *   lon(17),    lat(17),     val(17)
      save lon, lat, val, nval, kount
 
c     R15:
c     India:                   682 683 684 685 731 732 779
      data (lonr15(i),i=1,7)  / 10, 11, 12, 13, 11, 12, 11/
      data (latr15(i),i=1,7)  / 26, 26, 26, 26, 25, 25, 24/
c     N.America:               467 513 516 518 564 614
      data (lonr15(i),i=8,13) / 32, 33, 37, 38, 35, 36/
      data (latr15(i),i=8,13) / 30, 30, 30, 30, 29, 29/
 
c     T31:
c     India:
      data lont31
     *   / 19,20,21,22,23,  19,20,21,22,23, 20,21,22,23, 21,22, 21/
      data latt31
     *   / 32,32,32,32,32,  31,31,31,31,31, 30,30,30,30, 29,29, 28/
 
c ---------------------------------------------------------------------
 
c        Set coords, initialize accumulators, write header
 
      if (nstep.eq.nrstrt) then
 
         if (nlon.eq.48 .and. norec.eq.40) then
c          R15:
           nval = 13
           do 5 n=1,nval
             lon(n) = lonr15(n)
             lat(n) = latr15(n)
    5      continue
         else
c          T31:
           nval = 17
           do 6 n=1,nval
             lon(n) = lont31(n)
             lat(n) = latt31(n)
    6      continue
         endif
 
         call zero (val, nval)
         kount = 0
 
         zdr = 180./3.1415927
         write(iu,10)
     *      (lon(i)+(norec-lat(i))*nlon, i=1,nval)
         write(iu,11)
     *      (nint(along(lon(i))*zdr),nint(alatg(lat(i))*zdr), i=1,nval)
   10    format( /'precip (mm/day):'
     *          //'       ', 4x, 17i7)
   11    format(  '    day', 4x, 17(i4,',',i2) /)
 
      endif
 
c        Accumulate field values at particular points
 
      do 20 i=1,nval
        val(i) = val(i)
     *         + (  frainc(lon(i),lat(i)) + fsnowc(lon(i),lat(i))
     *            + frains(lon(i),lat(i)) + fsnows(lon(i),lat(i)) )
     *         * (8.64e7/rhoh2o)
   20 continue
      kount = kount + 1
 
c        If time, write out time-mean field values, zero accumulators
 
      if (mod(nstep,nint(86400./dtime)).eq.0 .and. nstep.ne.nrstrt) then
         day = nstep*dtime/86400.
         write(iu,30) day, ((val(i)/kount), i=1,nval)
   30    format(f7.2, 4x, 17f7.2)
         call zero (val, nval)
         kount = 0
c        flush file buffer
         ier = flush (iu)
      endif
 
      return
      end
C
C ***************************** USERDEF2 *******************************
C
      subroutine userdef2 (nstep, iu)
 
c User-codeable routine, called from end of linems.
c Prints instantaneous u,v,w,T,q for a particular agcm column.
 
c nstep  = number of timesteps (from base date)
c iu     = starting logical unit number
 
c ---------------------------------------------------------------------
#include <PARAMS.h>
c ---------------------------------------------------------------------
#include <COMCON.h>
c ---------------------------------------------------------------------
#include <COMLUN.h>
c-----------------------------------------------------------------------
#include <COML10.h>
c ---------------------------------------------------------------------
#include <COML30.h>
c ---------------------------------------------------------------------
#include <COMMAP.h>
c ---------------------------------------------------------------------
#include <COMMUN.h>
c ---------------------------------------------------------------------
      logical first
      save first, io, jo, ko1, ko2
      data first, io, jo, ko1, ko2 /.true., 49, 37, 1, 18/
c ---------------------------------------------------------------------
 
 
      if (first) then
        zlon = along(io)*180./3.1415927
        zlat = alatg(jo)*180./3.1415927
        write(iu  ,5) 'u (m/s)  ', zlon, zlat, (sig(k),k=ko1,ko2)
        write(iu+1,5) 'v (m/s)  ', zlon, zlat, (sig(k),k=ko1,ko2)
        write(iu+2,5) 'w (dp/dt)', zlon, zlat, (sig(k),k=ko1,ko2)
        write(iu+3,5) 't        ', zlon, zlat, (sig(k),k=ko1,ko2)
        write(iu+4,5) 'q (*1000)', zlon, zlat, (sig(k),k=ko1,ko2)
    5   format(a,'   lon=',f8.2,'   lat=',f8.2 / 6x,18f8.3 /)
        first = .false.
      endif
 
      write(iu  ,10) nstep, (yum1(io,jo,k), k=ko1,ko2)
      write(iu+1,10) nstep, (yvm1(io,jo,k), k=ko1,ko2)
      write(iu+2,10) nstep, (vomega(io,jo,k), k=ko1,ko2)
      write(iu+3,10) nstep, (ytm1(io,jo,k)-273.16, k=ko1,ko2)
      write(iu+4,10) nstep, (yqm1(io,jo,k)*1000., k=ko1,ko2)
   10 format(i4,2x,18f8.2)
 
      ier = flush (iu)
      ier = flush (iu+1)
      ier = flush (iu+2)
      ier = flush (iu+3)
      ier = flush (iu+4)
 
      return
      end
C
C ***************************** READER *********************************
C
c~~~~~~~~~~~~~~~~~~~~~~~~
#if ! defined (asyncmode) 
c~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE READER (IU, ARR, ILON,ILAT,ILEV, ICODEIN)
C----------------------------------------------------------------------
c
cl            Reads arr(ilon,ilat,ilev) from restart file on unit iu.
cl            If restart-file resolutions are different from model
cl            resolution, calls reader[h,v], else just read into arr.
cl
cl            icodein = 1 for agcm real fields, 2 for lsx real fields,
cl                     -1,-2 = same with no conversion in convertype.
cl
c
C----------------------------------------------------------------------
#include <PARAMS.h>
C---------------------------------------------------------------------
#include <COMFLG.h>
C----------------------------------------------------------------------
#include <COMHED.h>
C----------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      DIMENSION ARR (ILON,ILAT,ILEV)
c     20 is max restart-file vertical resol (see readerv):
      DIMENSION ARRES(NLON,NOREC,20)
      CHARACTER*1 BTYPE
C----------------------------------------------------------------------
 
      IF (ICODEIN.GE.0) THEN
        ICODE = ICODEIN
        BTYPE = 'R'
      ELSE
        ICODE = -ICODEIN
        BTYPE = 'C'
      ENDIF
 
c-------------------------
      IF (ICODE.EQ.1) THEN
c-------------------------
 
        if ( ilon.ne.nlon .or. ilat.ne.norec .or.
     *       (ilev.ne.nlev.and.ilev.ne.1) ) then
          write(*,8880) ilon,ilat,ilev,nlon,norec,nlev
 8880     format(/' *** Program error in reader (icode=1):',
     *            '     ilon,ilat,ilev,nlon,norec,nlev=',7i6/)
          call endrun (-1)
        endif
 
 
        IF (LEVRES.NE.NLEV .AND. ILEV.EQ.NLEV) THEN
 
c         Different agcm vertical resolution, 3-D field.
 
          DO 100 K=1,LEVRES
            IF (HORRES(1).NE.ILON .OR. HORRES(2).NE.ILAT) THEN
c             Also different agcm horizontal resolution
              CALL READERH
     *          (IU, ARRES(1,1,K), ILON, ILAT, HORRES(1), HORRES(2), 1,
     *           BTYPE)
            ELSE
c             Same horizontal resolution
              READ (IU, IOSTAT=IOS, ERR=20) (ARRES(I,1,K),I=1,ILON*ILAT)
              CALL CONVERTYPE (ARRES(1,1,K), ILON*ILAT, BTYPE,TYPERESIN)
            ENDIF
  100     CONTINUE
 
c         Interpolate to model vertical grid
          CALL READERV (ARR, ARRES, ILON, ILAT, ILEV, LEVRES)
 
 
        ELSE IF (HORRES(1).NE.ILON .OR. HORRES(2).NE.ILAT) THEN
 
c         Different agcm horizontal resolution, single-level field
c         or same vertical resolution.
 
          DO 150 K=1,ILEV
            CALL READERH
     *        (IU, ARR(1,1,K), ILON, ILAT, HORRES(1), HORRES(2), 1,
     *         BTYPE)
  150     CONTINUE
 
 
        ELSE
 
c         Same agcm resolutions
c         READ (IU, IOSTAT=IOS, ERR=20) ARR
c         One "record" per horiz slice, for machines without stream i/o
          DO K=1,ILEV
c           READ (IU, IOSTAT=IOS, ERR=20) (ARR(I,1,K),I=1,ILON*ILAT)
            READ (IU, IOSTAT=IOS, ERR=20) 
     *         ((ARR(I,J,K),I=1,ILON),J=1,ILAT)
          ENDDO
          CALL CONVERTYPE (ARR, ILON*ILAT*ILEV, BTYPE, TYPERESIN)
 
        ENDIF
 
c---------
      ELSE
c---------
 
        IF (HORRESLSX(1).NE.ILON .OR. HORRESLSX(2).NE.ILAT) THEN
 
c         Different lsx (horizontal) resolution
          DO 200 K=1,ILEV
            CALL READERH
     *        (IU, ARR(1,1,K), ILON, ILAT, HORRESLSX(1), HORRESLSX(2),2,
     *         BTYPE)
  200     CONTINUE
 
        ELSE
 
c         Same lsx (horizontal) resolution
c         READ (IU, IOSTAT=IOS, ERR=20) ARR
c         One "record" per horiz slice, for machines without stream i/o
          DO K=1,ILEV
c           READ (IU, IOSTAT=IOS, ERR=20) (ARR(I,1,K),I=1,ILON*ILAT)
            READ (IU, IOSTAT=IOS, ERR=20) 
     *        ((ARR(I,J,K),I=1,ILON),J=1,ILAT)
          ENDDO
          CALL CONVERTYPE (ARR, ILON*ILAT*ILEV, BTYPE, TYPERESIN)
 
        ENDIF
 
c----------
      ENDIF
c----------
 
      RETURN
 
   20 CALL FINDIU (IU, LOC)
      WRITE (NOUT,22) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC))),
     *                IU, ICODE, ILON, ILAT, ILEV, IOS
   22 FORMAT(//' *** Error or EOF reading restart file.'
     *        /'     local file name      =',a
     *        /'     unit number          =',i8
     *        /'     icode (1=agcm,2=sfc) =',i8
     *        /'     ilon,ilat,ilev       =',3i8
     *        /'     iostat error code    =',i8
     *        /'     Stop in READER.'/)
      CALL ENDRUN (-1)
 
      END
C
C ***************************** WRITER ********************************
C
       SUBROUTINE WRITER (IU, ARR, ILON,ILAT,ILEV, ICODE)
C----------------------------------------------------------------------
C
CL            Writes arr(ilon,ilat,ilev) to restart file (unit iu).
CL
CL            icode = 1 for agcm fields, = 2 for lsx fields
C
c---------------------------------------------------------------------
#include <PARAMS.h>
c---------------------------------------------------------------------
#include <COMLUN.h>
C----------------------------------------------------------------------
      DIMENSION ARR (ILON*ILAT*ILEV)
C----------------------------------------------------------------------
 
      IF (ILEV.GT.0) THEN
c         WRITE (IU, IOSTAT=IOS, ERR=20) ARR
c         One "record" per horiz slice, for machines without stream i/o
          DO K=1,ILEV
            WRITE (IU, IOSTAT=IOS, ERR=20)
     *        (ARR(I+(K-1)*ILON*ILAT),I=1,ILON*ILAT)
          ENDDO
      ENDIF
 
      RETURN
 
   20 CALL FINDIU (IU, LOC)
      WRITE (NOUT,22) LOCNAME(LOC)(1:LENCHR(LOCNAME(LOC))),
     *                IU, ICODE, ILON, ILAT, ILEV, IOS
   22 FORMAT(//' *** Error or EOF writing restart file.'
     *        /'     local file name      =',a
     *        /'     unit number          =',i8
     *        /'     icode (1=agcm,2=sfc) =',i8
     *        /'     ilon,ilat,ilev       =',3i8
     *        /'     iostat error code    =',i8
     *        /'     Stop in WRITER.'/)
      CALL ENDRUN (-1)
 
      END
c
c ***************************** READERH *******************************
c
      subroutine readerh (iu, arr, ilon, ilat, lonres, latres, icode,
     *                    btype)
 
c        Reads one horizontal slice from the restart file at
c        resolution (lonres,latres), interpolates to model grid at
c        resolution (ilon,ilat), and returns interpolated field in arr.
c        Interpolation is bilinear for agcm, nearest neighbor for
c        surface.
 
c        iu   = restart file unit number (supplied)
c        arr  = field on current model grid (returned)
c        ilon,ilat = dimensions of arr
c        lonres,latres = dimensions of restart file (supplied)
c        icode = 1 for agcm, 2 for surface (supplied).
c        btype = 'r' or 'c' for real or character type of arr
 
c---------------------------------------------------------------------
#include <PARAMS.h>
c---------------------------------------------------------------------
#include <COMHED.h>
c---------------------------------------------------------------------
      dimension arr(ilon,ilat)
      character*1 btype
c---------------------------------------------------------------------
      parameter (pi = 3.14159265358979)
c     360x360 is max restart-file horizontal resolution (new or old)
      dimension arres(360*360)
      dimension alon(360),    alat(360),
     *          alonres(360), alatres(360),
     *          indl(360),    indr(360),    weil(360),
     *          indb(360),    indt(360),   weib(360)
      save ilonprev, ilatprev, lonresprev, latresprev
      data ilonprev, ilatprev, lonresprev, latresprev /-1,-1,-1,-1/
      save indl, indr, weil, indb, indt, weib
c---------------------------------------------------------------------
 
c        If first call or resolutions changed since previous call
c        (eg, previous was for agcm, now for surface),
c        set model grid axes (alon,alat) and restart-file grid axes
c        (alonres,alatres), and calculate indices and weights to
c        interpolate from the restart grid to the model grid
c        (indl,indr,weil, indb,indt,weib). Ignore error in using
c        regularly spaced latitudes, not gaussian, for agcm.
 
c-----
      if (ilon.ne.ilonprev     .or. ilat.ne.ilatprev .or.
     *    lonres.ne.lonresprev .or. latres.ne.latresprev) then
c-----
 
        ilonprev = ilon
        ilatprev = ilat
        lonresprev = lonres
        latresprev = latres
 
        do 10 i=1,ilon
          alon(i) = ((i-.5)/ilon)*2.*pi
   10   continue
 
        do 12 j=1,ilat
          alat(j) = -.5*pi + ((j-.5)/ilat)*pi
   12   continue
 
        do 14 i=1,lonres
          alonres(i) = ((i-.5)/lonres)*2.*pi
   14   continue
 
        do 16 j=1,latres
          alatres(j) = -.5*pi + ((j-.5)/latres)*pi
   16   continue
 
        do 20 i=1,ilon
          val = alon(i)
          if (val.le.alonres(1)) then
            indl(i) = lonres
            indr(i) = 1
            weil(i) = (alonres(1)-val) / (2.*pi/lonres)
          else if (val.ge.alonres(lonres)) then
            indl(i) = lonres
            indr(i) = 1
            weil(i) = 1. - (val-alonres(lonres)) / (2.*pi/lonres)
          else
            do 22 io = 2,lonres
              if (val.lt.alonres(io)) then
                indl(i) = io-1
                indr(i) = io
                weil(i) = (alonres(io)-val) / (2.*pi/lonres)
                goto 24
              endif
   22       continue
   24       continue
          endif
   20   continue
 
        do 30 j=1,ilat
          val = alat(j)
          if (val.le.alatres(1)) then
            indb(j) = 1
            indt(j) = 1
            weib(j) = 1.
          else if (val.ge.alatres(latres)) then
            indb(j) = latres
            indt(j) = latres
            weib(j) = 1.
          else
            do 32 jo = 2,latres
              if (val.lt.alatres(jo)) then
                indb(j) = jo-1
                indt(j) = jo
                weib(j) =(alatres(jo)-val) / (alatres(jo)-alatres(jo-1))
                goto 34
              endif
   32       continue
   34       continue
          endif
   30   continue
 
c----------
      endif
c----------
 
c        Read one horizontal slice from restart file
 
      read (iu) (arres(m),m=1,lonres*latres)
      call convertype (arres, lonres*latres, btype, typeresin)
 
c        Do horizontal interpolation from restart-file grid to model
c        grid. For agcm (icode=1), do bilinear interp. For surface
c        (icode=2), take value of closest point, since for some lsx
c        fields, linear interp would cause problems...eg, with lmask,
c        blowups due to creating snow with very small thicknesses, etc.)
 
      if (icode.eq.1) then
 
        do 110 j=1,ilat
          do 112 i=1,ilon
            wlb =     weil(i) *    weib(j)
            wrb = (1.-weil(i))*    weib(j)
            wlt =     weil(i) *(1.-weib(j))
            wrt = (1.-weil(i))*(1.-weib(j))
            arr(i,j) = arres (indl(i)+(indb(j)-1)*lonres) * wlb
     *               + arres (indr(i)+(indb(j)-1)*lonres) * wrb
     *               + arres (indl(i)+(indt(j)-1)*lonres) * wlt
     *               + arres (indr(i)+(indt(j)-1)*lonres) * wrt
  112     continue
  110   continue
 
      else
 
        do 120 j=1,ilat
          do 122 i=1,ilon
            io = indl(i)
            if (weil(i).lt.0.5) io = indr(i)
            jo = indb(j)
            if (weib(j).lt.0.5) jo = indt(j)
            arr(i,j) = arres(io + (jo-1)*lonres)
  122     continue
  120   continue
 
      endif
 
      return
      end
c
c ***************************** READERV *******************************
c
       subroutine readerv (arr, arres, ilon, ilat, ilev, levres)
 
c        Called from subr reader to interpolate an agcm field for
c        a restart file with a different vertical resolution (ilev)
c        from the model's, by interpolating linearly in the vertical.
 
c        arr  = field on model grid (returned)
c        arres= field on restart-file's vertical grid (supplied)
c               (already on model's horizontal grid)
c        ilon,ilat,ilev = model dimensions of arr (=nlon,norec,nlev).
c        levres = restart-file's agcm vertical resolution
 
c---------------------------------------------------------------------
#include <PARAMS.h>
c---------------------------------------------------------------------
#include <COMLUN.h>
c---------------------------------------------------------------------
#include <COMMAP.h>
c---------------------------------------------------------------------
      dimension
     *  arr(ilon,ilat,ilev), arres(ilon,ilat,*)
c---------------------------------------------------------------------
      dimension
     *  sig6(6), sig12(12), sig16(16), sig18(18), sig20(20), sig36(36),
     *  sigres(100),
     *  jkm(nlev), jkp(nlev), wep(nlev)
 
      logical first
      save first, sigres, jkm, jkp, wep
      data first /.true./
c----------------------------------------------------------------------
c     6 layers:
      data sig6  /
     *  .020, .100, .250, .500, .850, .991 /
c----------------------------------------------------------------------
c     12 layers:
      data sig12 /
     *  .009, .025, .060, .110, .165, .245,
     *  .355, .500, .664, .811, .926, .991 /
c----------------------------------------------------------------------
c     16 layers (martin wild):
c     data sig16 /
c    *  .010, .030, .060, .110, .175, .255, .350, .450,
c    *  .550, .650, .745, .830, .900, .950, .980, .995 /
c----------------------------------------------------------------------
c     16 layers (med-res stratosphere):
      data sig16 /
     *  .010, .030, .060, .090, .120,
     *  .150, .180, .210, .250, .290,
     *  .355, .500, .664, .811, .926, .991 /
c----------------------------------------------------------------------
c     18 layers (ccm2):
      data sig18 /
     *  .005, .013, .033, .064, .099, .139,
     *  .189, .251, .325, .409, .501, .598,
     *  .695, .787, .866, .929, .970, .993 /
c----------------------------------------------------------------------
c     20 layers (hi-res stratosphere):
      data sig20 /
     *  .005, .010, .020, .040, .060, .080, .100,
     *  .120, .140, .160, .180, .200, .230, .280,
     *  .355, .500, .664, .811, .926, .991 /
c----------------------------------------------------------------------
c     36 layers (doubled sig18):
      data sig36 /
     *    .0050, .0090, .0130, .0230, .0330, .0485,
     *    .0640, .0815, .0990, .1190, .1390, .1640,
     *    .1890, .2200, .2510, .2880, .3250, .3670,
     *    .4090, .4550, .5010, .5495, .5980, .6465,
     *    .6950, .7410, .7870, .8265, .8660, .8975,
     *    .9290, .9495, .9700, .9815, .9930, .9965 /
c----------------------------------------------------------------------
 
c        If first call, set restart-file vertical grid (sigres) and
c        compute vertical-interpolation indices and weights into
c        model grid (sig)
 
c--------------------
      if (first) then
c--------------------
 
        if (levres.eq.6) then
          call scopy (levres, sig6, 1, sigres, 1)
        else if (levres.eq.12) then
          call scopy (levres, sig12, 1, sigres, 1)
        else if (levres.eq.16) then
          call scopy (levres, sig16, 1, sigres, 1)
        else if (levres.eq.18) then
          call scopy (levres, sig18, 1, sigres, 1)
        else if (levres.eq.20) then
          call scopy (levres, sig20, 1, sigres, 1)
        else if (levres.eq.36) then
          call scopy (levres, sig36, 1, sigres, 1)
        else
          write(nout,10) levres
   10     format
     *     (/' *** Error: Restart-file agcm vertical resolution LEVRES',
     *       ' =',i6
     *      /'     Must be one of the pre-set values (12,16,18,20)',
     *       ' in subr READERV.'/)
          call endrun (-1)
        endif
 
        do 100 k=1,ilev
          if (sigres(1).ge.sig(k)) then
            jkm(k) = 1
            jkp(k) = 1
            wep(k) = 1.
          else if (sigres(levres).le.sig(k)) then
            jkm(k) = levres
            jkp(k) = levres
            wep(k) = 1.
          else
            do 110 ko = 2,levres
              if (sigres(ko).ge.sig(k)) then
                jkm(k) = ko-1
                jkp(k) = ko
                wep(k) = (sig(k)-sigres(ko-1))/(sigres(ko)-sigres(ko-1))
                goto 112
              endif
  110       continue
  112       continue
          endif
  100   continue
 
        first = .false.
c----------
      endif
c----------
 
c         Vertically interpolate restart field arres to model field arr
 
      do 200 k=1,ilev
        do 210 j=1,ilat
          do 212 i=1,ilon
            arr(i,j,k) =  (1.-wep(k))*arres(i,j,jkm(k))
     *                 +      wep(k) *arres(i,j,jkp(k))
  212     continue
  210   continue
  200 continue
 
      return
      end

c~~~~~
#endif
c~~~~~
