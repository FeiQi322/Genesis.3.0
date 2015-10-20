c -- $Id: utilsx.F,v 1.3 2001/07/03 21:29:31 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/utilsx.F,v $
c -- $Name:  $

c-------
c utilsx
c-------
c
c----------------------------------------------------------------------
c
      subroutine setmap (alatg, wgaussg, lcent)
c
c        alatg  = gaussian latitude grid values, radians (from ccm)
c        wgaussg= gaussian weights for each latitude interval (from ccm)
c
c        sets index arrays indg,facg,indl,facl (in commapl) to be used
c        by ltog,gtol to map between ccm grid and lsx grid. 1-d
c        subscripts into the 2-d arrays are used for clarity, ie,
c        kg = ig+nlong*(jg-1) for ccm and kl = i+nlon*(j-1) for lsx.
c        the map from ccm to lsx is done by linear interpolation, and
c        indg(kl,1-4) = ccm index kg contributing a factor facg(kl,1-4)
c        times its field value at kg to that at kl. (this does not
c        conserve global integrals, so corrections are made in gtol.)
c        the map from lsx to ccm is done by simple "tiling", and
c        indl(kg,1-ntile) = lsx index kl contributing a factor
c        facl(kg,1-ntile) times its field value to that at kg. ntile
c        is the max number of lsx "tiles" that can be overlapped by
c        a ccm grid square (=20 for R15 and 2x2, set in compar).
c
c        also sets weig, weil (in commapl), areas of ccm and lsx boxes
c        relative to whole globe, used by gtol for global integrals.
c
c        the logic here and in vlocat is not general, but assumes that
c        lsx longitude box # 1 starts at 0, and ccm longitude box # 1
c        is centered on 0 (if lcent = true, or starts at 0 if false).
c        this avoids messy wraparound problems. if lcent = true, a
c        kludge is done for the lsx-to-ccm map by setting the edge
c        (rather than the center) of ccm box # 1 to 0 deg E, then
c        "duplicating" the contributions from lsx boxes in loop 222.
c
c        lcent is true  for ccm longitude box # 1 centered on 0 deg,
c        or false for that box's western edge at 0 deg. The lsx grid
c        is assumed to have box #1 left-hand edge at 0 deg.
c        If both grids are centered, fake by setting lcent=.false.,
c        as if both grids had edge at 0 deg. If both grids are edged,
c        then lcent=.false. Else lcent=.true. (Can't handle coarse grid
c        edged, fine grid centered).
c
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
      dimension alatg(nlatg), wgaussg(nlatg)
      logical lcent
c
      dimension blon(nlon+1),   blat(nlat+1),    wgauss(nlat),
     *          along(nlong),
     *          blong(nlong+1), blatg(nlatg+1),
     *          ig1(nlon),      ig2(nlon),       p(nlon),
     *          jg1(nlat),      jg2(nlat),       q(nlat),
     *          il1(nlong),     il2(nlong),
     *          jl1(nlatg),     jl2(nlatg)
c
      save zeroz
      data zeroz /0./
c
c
      twopi = 2.*pi
c
c        set lsx grids. alon, alat are centers of grid intervals, and
c        blon, blat are the edges. a sine scale is used for blat, since
c        this is used for the lsx-to-ccm tiling and the area between
c        two latitudes is proportional to d(sin(lat)).
c
      do i=1,nlon
        alon(i) = ((i-.5)/nlon)*twopi
      enddo
      do j=1,nlat
        alat(j) = -pi/2. + (j-0.5)*pi/nlat
      enddo
 
      do 10 i=1,nlon
        blon(i) = ((i-1.)/nlon)*twopi
   10 continue
      blon(nlon+1) = twopi
c
      if (nlat.ge.2) then
        do 14 j=2,nlat
          blat(j) = sin (0.5*(alat(j-1)+alat(j)))
   14   continue
      endif
      blat(1) = -1.
      blat(nlat+1) = 1.
 
      wtot = 0.
      do 16 j=1,nlat
        wgauss(j) = blat(j+1)-blat(j)
        wtot = wtot + wgauss(j)
   16 continue
      do 18 j=1,nlat
        wgauss(j) = wgauss(j) * (2./wtot)
   18 continue
c
c        set ccm grids. same comments as for lsx grids apply.
c        if lcent = true, ccm long box # 1 is centered on 0, but
c        set blong(1) to 0 for kludge in loop 222 below.
c        for latitude grid, use gaussian values alatg from ccm.
c
      if (lcent) then
        do 20 ig=1,nlong
          along(ig) = ((ig-1. )/nlong)*twopi
          blong(ig) = ((ig-1.5)/nlong)*twopi
   20   continue
        blong(1) = 0.
        blong(nlong+1) = ((nlong-0.5)/nlong)*twopi
      else
        do 22 ig=1,nlong
          along(ig) = ((ig-.5)/nlong)*twopi
          blong(ig) = ((ig-1.)/nlong)*twopi
   22   continue
        blong(nlong+1) = twopi
      endif
c
      if (nlatg.ge.2) then
        do 24 jg=2,nlatg
          blatg(jg) = sin (0.5*(alatg(jg-1)+alatg(jg)))
  24    continue
      endif
      blatg(1) = -1.
      blatg(nlatg+1) = 1.
c
c        *** same resolutions ***
c        if ccm and lsx resols are the same, enforce identical grids
c        (already the same for long, but not yet for gaussian lats.)
c

c#ifdef icedrive
      if (nlon.eq.nlong .and. nlat.eq.nlatg) then
        do 26 j=1,nlat
          jj = min (j,nlatg,nlat)   ! to avoid compiler warning
          alat(jj) = alatg(jj)
          wgauss(jj) = wgaussg(jj)
          blat(jj) = blatg(jj)
  26    continue
        jj = min (nlatg+1,nlat+1)   ! to avoid compiler warning
        blat(jj) = blatg(jj)
      endif
c#endif

c        set up for main ccm-to-lsx loop 110,112 below. ig1,ig2 are
c        the ccm center-longitudes that bracket each lsx longitude,
c        and p is the relative weight of ig1.
c
      do 100 i=1,nlon
        call vlocat(alon(i), along,nlong, twopi, ig1(i), ig2(i), p(i))
  100 continue
c
c        similarly, jg1,jg2 are the ccm center-latitudes that bracket
c        each lsx latitude, and q is the relative weight of jg1.
c
      do 102 j=1,nlat
        call vlocat(alat(j), alatg,nlatg, zeroz, jg1(j), jg2(j), q(j))
  102 continue
c
c        main loop to set ccm-to-lsx indices indg, facg. simple linear
c        interpolation is used between the 4 ccm center-points that
c        surround each lsx center-point.
c
      do 110 j=1,nlat
        do 112 i=1,nlon
          kl = i + nlon*(j-1)
c
          indg(kl,1) = ig1(i) + nlong*(jg1(j)-1)
          facg(kl,1) = p(i) * q(j)
c
          indg(kl,2) = ig2(i) + nlong*(jg1(j)-1)
          facg(kl,2) = (1.-p(i)) * q(j)
c
          indg(kl,3) = ig1(i) + nlong*(jg2(j)-1)
          facg(kl,3) = p(i) * (1.-q(j))
c
          indg(kl,4) = ig2(i) + nlong*(jg2(j)-1)
          facg(kl,4) = (1.-p(i)) * (1.-q(j))
  112   continue
  110 continue
c
c        set up for main lsx-to-ccm loop 210,212 below. jl1-to-jl2 is
c        the range of lsx latitude-intervals that are touched by each
c        ccm latitude-interval. this is done only for efficiency in
c        loop 210,212 below.
c
      do 200 jg=1,nlatg
       call vlocat(blatg(jg),   blat,nlat+1, zeroz, jl1(jg), m, f)
       call vlocat(blatg(jg+1), blat,nlat+1, zeroz, jl2(jg), m, f)
  200 continue
c
c        similarly, il1-to-il2 is the range of lsx longitude-intervals
c        that are touched by each ccm longitude-interval.
 
      do 202 ig=1,nlong
       call vlocat(blong(ig),   blon,nlon+1, twopi, il1(ig), m, f)
       call vlocat(blong(ig+1), blon,nlon+1, twopi, il2(ig), m, f)
  202 continue
c
c        main loop to set lsx-to-ccm indices indl, facl. for each
c        ccm box, the area of overlap with each of the lsx boxes it
c        touches is calculated. indl is the index of each such lsx box,
c        and its overlapped area divided by the ccm-box area is facl.
c
c        if lcent = true, for ccm long box # 1 (with fake blong(1)=0),
c        double the ccm box area areag, and have additional contrib in
c        loop 222 from lsx boxes "reflected" about longitude 0.
c
      do 210 jg=1,nlatg
        do 212 ig=1,nlong
          kg = ig + nlong*(jg-1)
          do 215 l=1,ntile
            indl(kg,l) = 1
            facl(kg,l) = 0.
  215     continue
c
          ng = 0
          areag = (blong(ig+1)-blong(ig)) * (blatg(jg+1)-blatg(jg))
          if (lcent .and. ig.eq.1) areag = 2.*areag
c
          do 220 j=jl1(jg),jl2(jg)
            dy = min(blat(j+1),blatg(jg+1)) - max(blat(j),blatg(jg))
            do 222 i=il1(ig),il2(ig)
              dx = min(blon(i+1),blong(ig+1)) - max(blon(i),blong(ig))
              area = dx*dy
              if (area.gt.0.) then
c
                 ng = ng + 1
                 indl(kg,ng) = i + nlon*(j-1)
                 facl(kg,ng) = area / areag
c
                 if (lcent .and. ig.eq.1) then
                   ng = ng + 1
                   indl(kg,ng) = (nlon+1-i) + nlon*(j-1)
                   facl(kg,ng) = area / areag
                 endif
c
                 if (ng.gt.ntile) then
                    write(ioterm2,900) ig,jg,ntile
                    call endrun (-1)
                 endif
c
              endif
  222       continue
  220     continue
c
  212   continue
  210 continue
c
c        set areal weights weil, weig. note sum of wgaussg should = 2.
c        don't forget to double weig for fake blong(1) if lcent = true.
c
      do 300 j=1,nlat
        do 302 i=1,nlon
          kl = i + nlon*(j-1)
          weil(kl) = (wgauss(j)/2.)*((blon(i+1)-blon(i))/twopi)
  302   continue
  300 continue
c
      do 310 jg=1,nlatg
        do 312 ig=1,nlong
          kg = ig + nlong*(jg-1)
          weig(kg) = (wgaussg(jg)/2.)* ((blong(ig+1)-blong(ig))/twopi)
          if (lcent .and. ig.eq.1) weig(kg) = 2.*weig(kg)
  312   continue
  310 continue
c
      return
c
  900 format(/' *** Fatal error ***'/
     *       /' number of lsx overlays for ccm surface grid box',2i4,
     *        ' exceeds limit',i4,' (setmap).'/)
c
      end
c
c----------------------------------------------------------------------
c
      subroutine vlocat (valin, arr, narr, wrapval, iv1, iv2, frac)
c
c        locates a given value in a monotonic array, and calculates
c        the appropriate linear interpolation factor. the array can
c        wrap around mod(wrapval), or not if wrapval = 0.
c
c        valin   = value to be located. If wrapval ne 0, valin is
c                  shifted mod wrapval (in copy to local variable val)
c                  to be ge 0 and le wrapval (supplied)
c        arr     = array of monotonically increasing, possibly
c                  irregularly spaced values (supplied)
c        narr    = size of arr (supplied)
c        wrapval = wraparound value for arr, or zero for no wraparound.
c                  If wrapval ne 0, arr(1) must be ge 0 and arr(narr)
c                  must be le wrapval. (supplied)
c        iv1,iv2 = indices such that arr(iv1) le val le arr(iv2)
c                  (returned)
c        frac    = linear-interp weighting factor for iv1 (returned)
c
c
#include <compar.h>
      dimension arr(narr)
 
 
c        Protection for single-point tests (local machine).
c        Code below requires narr >= 2.
 
      if (narr.eq.1) then
        iv1 = 1
        iv2 = 1
        frac = 1.
        return
      endif
c
      if (wrapval.gt.0. .and. (valin.lt.0. .or. valin.gt.wrapval)) then
        val = mod (valin, wrapval)
        if (val.lt.0.) val = val + wrapval
      else
        val = valin
      endif
 
      if (val.lt.arr(1)) then
         if (wrapval.eq.0.) then
            iv1 = 1
            iv2 = 2
            frac = 1.
         else
            iv1 = narr
            iv2 = 1
            frac = (arr(1)-val) / (arr(1)+wrapval-arr(narr))
         endif
c
      else if (val.gt.arr(narr)) then
         if (wrapval.eq.0.) then
            iv1 = narr-1
            iv2 = narr
            frac = 0.
         else
            iv1 = narr
            iv2 = 1
            frac = (arr(1)+wrapval-val) / (arr(1)+wrapval-arr(narr))
         endif
c
      else
         do 10 i=2,narr
            if (val.le.arr(i)) then
               iv1 = i-1
               iv2 = i
               frac = (arr(i)-val) / (arr(i)-arr(i-1))
               return
            endif
   10    continue
c
      endif
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine gtol (agcm, alsx, icode)
c
c        maps gcm global array agcm to lsx global array alsx. if icode
c        ne 0, corrects alsx to conserve global average either by
c        adding a constant amount (icode gt 0) or by multiplying by a
c        constant factor (icode lt 0). the latter is for fields
c        than cannot be negative (eg, rainfall) and/or whose zero
c        values should remain zero (eg, solar flux). for simplicity,
c        1d indices into 2d arrays are used (ie, lon + nlon*(lat-1)).
c        uses index arrays indg,facg and global integral weights
c        weig,weil set up by setmap.
c
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
      dimension agcm(nlong*nlatg), alsx(nlon*nlat)
c
c
c        *** same resolutions ***
c        if ccm and lsx resols are the same,just copy arrays and return
c
c#ifdef icedrive
      if (nlon.eq.nlong .and. nlat.eq.nlatg) then
        call scopy (nlon*nlat, agcm, 1, alsx, 1)
        return
      endif
c#endif
c
c        calculate new lsx field alsx from given ccm field agcm
c
      call zero (alsx, nlon*nlat)
c
      do 10 n=1,4
         do 12 kl=1,nlon*nlat
            alsx(kl) = alsx(kl) + facg(kl,n)*agcm(indg(kl,n))
   12    continue
   10 continue
c
c---------------------------
      if (icode.eq.0) return
c---------------------------
c
c        if icode ne 0, adjust new lsx fld so global avg is conserved
c
      avg = 0.
      do 20 kg=1,nlong*nlatg
        avg = avg + agcm(kg)*weig(kg)
   20 continue
c
      avl = 0.
      do 22 kl=1,nlon*nlat
        avl = avl + alsx(kl)*weil(kl)
   22 continue
c
c-------------------------
      if (icode.gt.0) then
c-------------------------
        dav = avg-avl
        do 30 kl=1,nlon*nlat
          alsx(kl) = alsx(kl) + dav
   30   continue
c------------------------------
      else if (icode.lt.0) then
c------------------------------
        if (abs(avl).gt.0.) then
          dav = avg/avl
          do 32 kl=1,nlon*nlat
            alsx(kl) = alsx(kl) * dav
   32     continue
        endif
c----------
      endif
c----------
c
c     write(6,900) icode,avg,avl
c 900 format(' gtol: icode, global means(g,l)=',i2,2e20.10)
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine gtolsol (agcm, alsx, cofx)
c
c        adjusts lsx field alsx, by multiplying by a constant factor,
c        so that global average of alsx*(1-cofx) = that of
c        agcm*(1-cofg), where cofg is the map of cofx to the gcm
c        grid. used to adjust sola[d,i]xy to make lsx global absorbed
c        solar flux consistent with that implied by gcm rad calcs
c        (called from surfctl for each waveband and direct/diffuse
c        separately, with cofx = lsx surface albedos). note that 1d
c        indices into 2d arrays are used (ie, lon + nlon*(lat-1)).
c        uses index arrays indg,facg and global weights weig,weil
c        set up by setmap.
c
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
      dimension agcm(nlong*nlatg), alsx(nlon*nlat), cofx(nlon*nlat)
c
c
c        interpolate lsx cofx (sfc albedo) to ccm grid, as was done
c        at end of alblsx. (use temporary workg1 in comgrd.)
c
      call ltog (cofx, workg1, 0)
c
c        compute ccm and lsx global mean absorbed solar fluxes
c
      avg = 0.
      do 20 kg=1,nlong*nlatg
        avg = avg + (1.-workg1(kg))*agcm(kg)*weig(kg)
   20 continue
c
      avl = 0.
      do 22 kl=1,nlon*nlat
        avl = avl + (1.-cofx(kl))*alsx(kl)*weil(kl)
   22 continue
c
c        adjust lsx solar incident flux to equate global absorbed flxs
c
      if (abs(avl).gt.0.) then
        dav = avg/avl
        do 32 kl=1,nlon*nlat
          alsx(kl) = alsx(kl) * dav
   32   continue
      endif
c
c     write(6,900) avg,avl
c 900 format(' gtolsol. avg,avl: ',2e20.10)
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine gtolnear (agcm, alsx)
c
c        maps gcm global array agcm to lsx global array alsx, simply
c        by taking the nearest ccm neighbor value to each lsx point.
c        used by other dataset generation programs using bats r15 soil
c        texture and color maps. makes use of index arrays indg,facg
c        set up by setmap.
c
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
      dimension agcm(nlong*nlatg), alsx(nlon*nlat)
c
c
c        *** same resolutions ***
c        if ccm and lsx resols are the same,just copy arrays and return
c
c#ifdef icedrive
      if (nlon.eq.nlong .and. nlat.eq.nlatg) then
        do 2 k=1,nlon*nlat
          alsx(k) = agcm(k)
    2   continue
        return
      endif
c#endif
c
      do 10 kl=1,nlon*nlat
        fmax = -1.e20
        do 12 n=1,4
          if (facg(kl,n).gt.fmax) then
            nmax = n
            fmax = facg(kl,n)
          endif
   12   continue
        alsx(kl) = agcm(indg(kl,nmax))
   10 continue
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine ltog (alsx, agcm, icode)
c
c        maps lsx global array alsx to gcm global array agcm.
c        similar comments as for gtol apply here, including correction
c        of alsx to conserve global average, either by adding
c        a constant amount (if icode gt 0) or by multiplying by a
c        constant factor (if icode lt 0). the indices indl, facl
c        used here represent "tiling" (set in setmap), but the
c        global integral is not conserved since the ccm latitudinal
c        integration is by gaussian quadrature, not by adding up the
c        squares.
c
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
      dimension alsx(nlon*nlat), agcm(nlong*nlatg)
c
c
c        *** same resolutions ***
c        if ccm and lsx resols are the same,just copy arrays and return
c
c#ifdef icedrive
      if (nlon.eq.nlong .and. nlat.eq.nlatg) then
        call scopy (nlon*nlat, alsx, 1, agcm, 1)
        return
      endif
c#endif
c
c        calculate new ccm field agcm from given lsx field alsx
c
      call zero (agcm, nlong*nlatg)
c
      do 10 n=1,ntile
         do 12 kg=1,nlong*nlatg
            agcm(kg) = agcm(kg) + facl(kg,n)*alsx(indl(kg,n))
   12    continue
   10 continue
c
c---------------------------
      if (icode.eq.0) return
c---------------------------
c
c        if icode ne 0, adjust new ccm fld so global avg is conserved
c
      avg = 0.
      do 20 kg=1,nlong*nlatg
        avg = avg + agcm(kg)*weig(kg)
   20 continue
c
      avl = 0.
      do 22 kl=1,nlon*nlat
        avl = avl + alsx(kl)*weil(kl)
   22 continue
c
c-------------------------
      if (icode.gt.0) then
c-------------------------
        dav = avl-avg
        do 30 kg=1,nlong*nlatg
          agcm(kg) = agcm(kg) + dav
   30   continue
c------------------------------
      else if (icode.lt.0) then
c------------------------------
        if (abs(avg).gt.0.) then
          dav = avl/avg
          do 32 kg=1,nlong*nlatg
            agcm(kg) = agcm(kg) * dav
   32     continue
        endif
c----------
      endif
c----------
c
c     write(6,900) icode,avl,avg
c 900 format(' ltog. icode, global means(g,l)=',i2,2e20.10)
c
      return
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE SETIND
c
c        Sets arrays of indices INDALL_[L,I,O] to the lon/lat locations
c        of land/icesheet/ocean points on the surface grid, to be
c        processed in 1D strips. Also sets NPALL_[L,I,O], the number of
c        points in each strip (= NVEC except for the last strip),
c        and NSTRIP_[L,I,O], the number of strips.
c
c        Called once at the beginning of a run; thereafter, GETIND
c        will copy the ISTRIP-th strip to 1D strip arrays and set NPOI
c        to NPALL_[L,I,O](ISTRIP), and PUTIND will copy them back to
c        global arrays. (INDALL_*, NPALL_*, NSTRIP_* and NPOI are in
c        comgrd.)
c
#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
#include <comocnxy.h>
 
c----------------------------
c        Land (non-icesheet):
c----------------------------
 
      NSTRIP_L = 1
      I = 0
      DO JJ=1,NLAT
        DO JI=1,NLON
          IF (LMASKXY(JI,JJ).EQ.1) THEN
            I = I+1
            INDALL_L (I,NSTRIP_L) = JI + (JJ-1)*NLON
            IF (I.EQ.NVEC) THEN
               NPALL_L (NSTRIP_L) = I
               NSTRIP_L = NSTRIP_L+1
               I = 0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c
#ifdef icedrive
#ifdef lake
      I = 0
      NSTRIP_L = 1
#endif
#endif
c
c        If I=0, no pts between end of last strip and global end
c
      IF (I.EQ.0) THEN
        NSTRIP_L = NSTRIP_L - 1
      ELSE
        NPALL_L (NSTRIP_L) = I
      ENDIF
 
c------------------
c        Ice sheet:
c------------------
 
      NSTRIP_I = 1
      I = 0
      DO JJ=1,NLAT
        DO JI=1,NLON
          IF (LMASKXY(JI,JJ).EQ.2) THEN
            I = I+1
            INDALL_I (I,NSTRIP_I) = JI + (JJ-1)*NLON
            IF (I.EQ.NVEC) THEN
               NPALL_I (NSTRIP_I) = I
               NSTRIP_I = NSTRIP_I+1
               I = 0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c
#ifdef hydmodel
      I = 0
      NSTRIP_I = 1
#endif
c
c        If I=0, no pts between end of last strip and global end
c
      IF (I.EQ.0) THEN
        NSTRIP_I = NSTRIP_I - 1
      ELSE
        NPALL_I (NSTRIP_I) = I
      ENDIF
 
c----------------------
c        Ocean/sea-ice:
c----------------------
 
      NSTRIP_O = 1
      I = 0
      DO JJ=1,NLAT
        DO JI=1,NLON
          IF (FWATERXY(JI,JJ).GT.0.) THEN
            I = I+1
            INDALL_O (I,NSTRIP_O) = JI + (JJ-1)*NLON
            IF (I.EQ.NVEC) THEN
               NPALL_O (NSTRIP_O) = I
               NSTRIP_O = NSTRIP_O+1
               I = 0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c
#ifdef icedrive
#ifndef lake
      I = 0
      NSTRIP_O = 1
#endif
#endif
c
c        If I=0, no pts between end of last strip and global end
c
      IF (I.EQ.0) THEN
        NSTRIP_O = NSTRIP_O - 1
      ELSE
        NPALL_O (NSTRIP_O) = I
      ENDIF
 
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE GETIND (ISTRIP, INDALL, NPALL, ITYPE)
c
c        Copies current strip from global arrays to 1D arrays.
c        Also sets NPOI and NSTRIP (comgrd).
c        Uses info from INDALL_[L,I,O] and NPALL_[L,I,O] previously set
c        in SETIND, passed to local arrays INDALL,NPALL.
c
c        ISTRIP = current strip number
c        INDALL = lon/lat encoded indices for each strip
c        NPALL  = number of points in each strip
c        ITYPE  = 1 for land strips,
c                 2 for icesheet strips (don't need veg),
c                 3 for ocean/seaice strips (don't need veg).
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
#include <commud.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <commus.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#include <com1d.h>
      DIMENSION INDALL(NVEC,NALL), NPALL(NALL)
c
c
      IF (ITYPE.EQ.1) THEN
        NSTRIP = NSTRIP_L
      ELSE IF (ITYPE.EQ.2) THEN
        NSTRIP = NSTRIP_I
      ELSE IF (ITYPE.EQ.3) THEN
        NSTRIP = NSTRIP_O
      ENDIF
c
      NPOI = NPALL(ISTRIP)
c
c        long,lat indices of corresponding lsx 2-D grid
c
      DO 50 I=1,NPOI
        IZ = INDALL(I,ISTRIP)
        LOCPOI(I,1) = MOD(IZ-1,NLON) + 1
        LOCPOI(I,2) = (IZ-1)/NLON + 1
   50 CONTINUE
c
c-------------------------
      IF (ITYPE.EQ.1) THEN
c-------------------------
c
c          Arrays (non-layer,...) in comveg[xy]
c
!DIR$ IVDEP
        DO 100 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
c         comveg:
          TU(I) = TUXY(JI,JJ)
          TS(I) = TSXY(JI,JJ)
          TL(I) = TLXY(JI,JJ)
          TLSUB(I) = TLSUBXY(JI,JJ)
          WLIQU(I) = WLIQUXY(JI,JJ)
          WLIQS(I) = WLIQSXY(JI,JJ)
          WLIQL(I) = WLIQLXY(JI,JJ)
          WSNOU(I) = WSNOUXY(JI,JJ)
          WSNOS(I) = WSNOSXY(JI,JJ)
          WSNOL(I) = WSNOLXY(JI,JJ)
          T12(I)   = T12XY(JI,JJ)
          T34(I)   = T34XY(JI,JJ)
          Q12(I)   = Q12XY(JI,JJ)
          Q34(I)   = Q34XY(JI,JJ)
          FOG12(I) = FOG12XY(JI,JJ)
          FOG34(I) = FOG34XY(JI,JJ)
          HFOG12(I) = HFOG12XY(JI,JJ)
          HFOG34(I) = HFOG34XY(JI,JJ)
#ifdef ibis
          VMAXSUNU(I,1) = VMAXSUNUXY(JI,JJ,1)
          VMAXSUNU(I,2) = VMAXSUNUXY(JI,JJ,2)
          VMAXSUNL(I,1) = VMAXSUNLXY(JI,JJ,1)
          VMAXSUNL(I,2) = VMAXSUNLXY(JI,JJ,2)
          VMAXSHAU(I,1) = VMAXSHAUXY(JI,JJ,1)
          VMAXSHAU(I,2) = VMAXSHAUXY(JI,JJ,2)
          VMAXSHAL(I,1) = VMAXSHALXY(JI,JJ,1)
          VMAXSHAL(I,2) = VMAXSHALXY(JI,JJ,2)
 
          CISUNUD(I) = CISUNUDXY(JI,JJ)
          CISUNUC(I) = CISUNUCXY(JI,JJ)
          CISUNL3(I) = CISUNL3XY(JI,JJ)
          CISUNL4(I) = CISUNL4XY(JI,JJ)
          CISHAUD(I) = CISHAUDXY(JI,JJ)
          CISHAUC(I) = CISHAUCXY(JI,JJ)
          CISHAL3(I) = CISHAL3XY(JI,JJ)
          CISHAL4(I) = CISHAL4XY(JI,JJ)
 
          CSSUNUD(I) = CSSUNUDXY(JI,JJ)
          CSSUNUC(I) = CSSUNUCXY(JI,JJ)
          CSSUNL3(I) = CSSUNL3XY(JI,JJ)
          CSSUNL4(I) = CSSUNL4XY(JI,JJ)
          CSSHAUD(I) = CSSHAUDXY(JI,JJ)
          CSSHAUC(I) = CSSHAUCXY(JI,JJ)
          CSSHAL3(I) = CSSHAL3XY(JI,JJ)
          CSSHAL4(I) = CSSHAL4XY(JI,JJ)
 
          CONFRAC(I) = CONFRACXY(JI,JJ)
          C4FRAC (I) = C4FRACXY (JI,JJ)
 
          ANUD(I) = ANUDXY(JI,JJ)
          ANUC(I) = ANUCXY(JI,JJ)
          ANL3(I) = ANL3XY(JI,JJ)
          ANL4(I) = ANL4XY(JI,JJ)
 
          AGUD(I) = AGUDXY(JI,JJ)
          AGUC(I) = AGUCXY(JI,JJ)
          AGL3(I) = AGL3XY(JI,JJ)
          AGL4(I) = AGL4XY(JI,JJ)
 
          GSUD(I) = GSUDXY(JI,JJ)
          GSUC(I) = GSUCXY(JI,JJ)
          GSL3(I) = GSL3XY(JI,JJ)
          GSL4(I) = GSL4XY(JI,JJ)
#endif
  100   CONTINUE
c
c          Arrays (stories) in comveg[xy]
c
        DO 110 IV=1,2
!DIR$ IVDEP
          DO 112 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            LAI(I,IV) = LAIXY(JI,JJ,IV)
            SAI(I,IV) = SAIXY(JI,JJ,IV)
            ORIEV(I,IV) = ORIEVXY(JI,JJ,IV)
            ORIEH(I,IV) = ORIEHXY(JI,JJ,IV)
            ZTOP(I,IV) = ZTOPXY(JI,JJ,IV)
            ZBOT(I,IV) = ZBOTXY(JI,JJ,IV)
            DLEAF(I,IV) = DLEAFXY(JI,JJ,IV)
            DSTEM(I,IV) = DSTEMXY(JI,JJ,IV)
            FSTOMA(I,IV) = FSTOMAXY(JI,JJ,IV)
            RSMIN(I,IV) = RSMINXY(JI,JJ,IV)
            RSMAX(I,IV) = RSMAXXY(JI,JJ,IV)
            PARCON(I,IV) = PARCONXY(JI,JJ,IV)
            VPDCON(I,IV) = VPDCONXY(JI,JJ,IV)
            TEMCON(I,IV) = TEMCONXY(JI,JJ,IV)
            TEMOPT(I,IV) = TEMOPTXY(JI,JJ,IV)
  112     CONTINUE
  110   CONTINUE
c
c          Arrays (stories, wavebands) in comveg[xy]
c
        DO 120 IV=1,2
          DO 122 IB=1,NBAND
!DIR$ IVDEP
            DO 124 I=1,NPOI
              JI = LOCPOI(I,1)
              JJ = LOCPOI(I,2)
              RHOVEG(I,IB,IV) = RHOVEGXY(JI,JJ,IB,IV)
              TAUVEG(I,IB,IV) = TAUVEGXY(JI,JJ,IB,IV)
  124       CONTINUE
  122     CONTINUE
  120   CONTINUE
c
c          Arrays (stories, soil layers) in comveg[xy]
c
        DO 130 IV=1,2
          DO 132 K=1,NSOILAY
!DIR$ IVDEP
            DO 134 I=1,NPOI
              JI = LOCPOI(I,1)
              JJ = LOCPOI(I,2)
              FROOT(I,K,IV) = FROOTXY(JI,JJ,K,IV)
  134       CONTINUE
  132     CONTINUE
  130   CONTINUE
c
c----------
      ENDIF
c----------
 
c-----------------------------------------
      IF (ITYPE.EQ.1 .OR. ITYPE.EQ.2) THEN
c-----------------------------------------
c
c          Arrays (non-layer,...) in comsno[xy], comsoi[x,y]
c
!DIR$ IVDEP
        DO 150 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
c         comsno:
          FI(I) = FIXY(JI,JJ)
          TI(I) = TIXY(JI,JJ)
c         comsoi:
          WPUD(I)   = WPUDXY(JI,JJ)
          WIPUD(I)  = WIPUDXY(JI,JJ)
          Z0SOI(I)  = Z0SOIXY(JI,JJ)
          ALBSAV(I) = ALBSAVXY(JI,JJ)
          ALBSAN(I) = ALBSANXY(JI,JJ)
          TG(I)     = TGXY(JI,JJ)
  150   CONTINUE
c
c          Arrays (soil layers) in comsoi[xy]
c
        DO 160 K=1,NSOILAY
!DIR$ IVDEP
          DO 162 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            ITEX(I,K) = ITEXXY(JI,JJ,K)
            TSOI(I,K) = TSOIXY(JI,JJ,K)
            WSOI(I,K) = WSOIXY(JI,JJ,K)
            WISOI(I,K)= WISOIXY(JI,JJ,K)
  162     CONTINUE
  160   CONTINUE
c
c          Arrays (snow layers) in comsno[xy]
c
        DO 170 K=1,NSNOLAY
!DIR$ IVDEP
          DO 172 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            TSNO(I,K) = TSNOXY(JI,JJ,K)
            HSNO(I,K) = HSNOXY(JI,JJ,K)
  172     CONTINUE
  170   CONTINUE
c
c---------
      ENDIF
c---------
 
c-------------------------
      IF (ITYPE.EQ.3) THEN
c-------------------------
 
c          Arrays (non-layer) in comocn[xy], comice[xy]
 
!DIR$ IVDEP
        DO 200 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
c         comocn:
          TOCN(I)  = TOCNXY(JI,JJ)
          SOCN(I)  = SOCNXY(JI,JJ)
          TFREEZ(I)= TFREEZXY(JI,JJ)
          HOCN(I)  = HOCNXY(JI,JJ)
c         for mom2,4:
          FHOCN(I)   = FHOCNXY(JI,JJ)
          FQOCN(I)   = FQOCNXY(JI,JJ)
          FSOLOCN(I) = FSOLOCNXY(JI,JJ)
          FIROCN(I)  = FIROCNXY(JI,JJ)
          FPRLOCN(I) = FPRLOCNXY(JI,JJ)
          FPRFOCN(I) = FPRFOCNXY(JI,JJ)
          FSALOCN(I) = FSALOCNXY(JI,JJ)
          PRESOCN(I) = PRESOCNXY(JI,JJ)

          FTOCN(I) = FTOCNXY(JI,JJ)
          FWOCN(I) = FWOCNXY(JI,JJ)
          FUOCN(I) = FUOCNXY(JI,JJ)
          FVOCN(I) = FVOCNXY(JI,JJ)
          QFLUX(I) = QFLUXXY(JI,JJ)
          FUICE(I) = FUICEXY(JI,JJ)
          FVICE(I) = FVICEXY(JI,JJ)
c         comice:
          FICE(I)  = FICEXY(JI,JJ)
          BRINE(I) = BRINEXY(JI,JJ)
          TG_ICE(I)= TG_ICEXY(JI,JJ)
c         comsno:
          FI(I)    = FI_ICEXY(JI,JJ)
          TI(I)    = TI_ICEXY(JI,JJ)
#ifdef isotrac
          xo_o(i)  = xoxy_o(ji,jj)
          xo_h(i)  = xoxy_h(ji,jj)
          qo_i(i)  = qoxy_i(ji,jj)
          zf_i(i)  = zfxy_i(ji,jj)
          fice_i(i)= ficexy_i(ji,jj)
          fvapo(i) = 0.
          fvi_i(i) = fvixy_i(ji,jj)
#endif
  200   CONTINUE
c
c          Arrays (ice layers) in comice[xy]
c
        DO 240 K=1,NICELAY
!DIR$ IVDEP
          DO 242 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            HICE(I,K) = HICEXY(JI,JJ,K)
            TICE(I,K) = TICEXY(JI,JJ,K)
  242     CONTINUE
  240   CONTINUE
c
c          Arrays (snow layers) in comice[xy] -> comsno
c
        DO 250 K=1,NSNOLAY
!DIR$ IVDEP
          DO 252 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            TSNO(I,K) = TSNO_ICEXY(JI,JJ,K)
            HSNO(I,K) = HSNO_ICEXY(JI,JJ,K)
  252     CONTINUE
  250   CONTINUE
 
c----------
      ENDIF
c----------
 
c
!DIR$ IVDEP
      DO 300 I=1,NPOI
        JI = LOCPOI(I,1)
        JJ = LOCPOI(I,2)
 
c       Arrays in comveg[xy]
 
        IVEG(I) = IVEGXY(JI,JJ)
        FU(I) = FUXY(JI,JJ)
        FL(I) = FLXY(JI,JJ)
 
c       Arrays in commug[xy]
 
        COSZEN(I) = COSZENXY(JI,JJ)
        FIRA(I) = FIRAXY(JI,JJ)
        UX(I) = UXXY(JI,JJ)
        UY(I) = UYXY(JI,JJ)
        TA(I) = TAXY(JI,JJ)
        QA(I) = QAXY(JI,JJ)
        RAINA(I) = RAINAXY(JI,JJ)
        SNOWA(I) = SNOWAXY(JI,JJ)
        PSURF(I) = PSURFXY(JI,JJ)
        AFIRE(I) = AFIREXY(JI,JJ)
        ASENA(I) = ASENAXY(JI,JJ)
        AVAPA(I) = AVAPAXY(JI,JJ)
        ATAUX(I) = ATAUXXY(JI,JJ)
        ATAUY(I) = ATAUYXY(JI,JJ)
        LMASK(I) = LMASKXY(JI,JJ)
        FWATER(I)= FWATERXY(JI,JJ)
 
c       Arrays in commud[xy]
 
        GINVAP(I) = GINVAPXY(JI,JJ)
        GSUVAP(I) = GSUVAPXY(JI,JJ)
        GTRANSU(I)= GTRANSUXY(JI,JJ)
        GTRANSL(I)= GTRANSLXY(JI,JJ)
        GRUNOF(I) = GRUNOFXY(JI,JJ)
        GDRAIN(I) = GDRAINXY(JI,JJ)
        GDRAINH(I) = GDRAINHXY(JI,JJ)
        DO K=1,NSOILAY+1
          WFLO(I,K) = WFLOXY(JI,JJ,K)
        ENDDO
        GABLAT(I) = GABLATXY(JI,JJ)
        GADJUST(I)= GADJUSTXY(JI,JJ)
        TS2(I)    = TS2XY(JI,JJ)
#if defined (rcm) 
        ts10(i) = ts10xy(ji,jj)
        qs10(i) = qs10xy(ji,jj)
        us10(i) = us10xy(ji,jj)
#endif
        Z0LEN(I)  = Z0LENXY(JI,JJ)
        GINFIL(I)  = GINFILXY(JI,JJ)
        GIRRIG(I)  = GIRRIGXY(JI,JJ)
 
  300 CONTINUE
c
c        Arrays (wavebands) in commug[xy]
c
      DO 310 IB=1,NBAND
!DIR$ IVDEP
        DO 312 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
          SOLAD(I,IB) = SOLADXY(JI,JJ,IB)
          SOLAI(I,IB) = SOLAIXY(JI,JJ,IB)
          ASURD(I,IB) = ASURDXY(JI,JJ,IB)
          ASURI(I,IB) = ASURIXY(JI,JJ,IB)
  312   CONTINUE
  310 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE PUTIND (ITYPE)
c
c        Copies current strip from 1D arrays to global arrays.
c        Uses info in LOCPOI(I,1:2) in comgrd, set previously in SETIND.
c        If ITYPE = 3 (ocean/seaice), weight fluxes to be returned to
c        agcm (A*), and overall surface albedos (ASUR*), with those
c        previously calculated for land or icesheet, using fraction of
c        open water FWATER (copied to local ZW).
c
c        ITYPE  = 1 for land processing,
c                 2 for icesheet processing (don't need veg),
c                 3 for ocean/seaice processing (don't need veg).
c
#include <compar.h>
#include <comgrd.h>
#include <comice.h>
#include <comicexy.h>
#include <commud.h>
#include <commudxy.h>
#include <commug.h>
#include <commugxy.h>
#include <commus.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
#ifdef isotrac
#include <com1d.h>
#endif
      DIMENSION ZW(NVEC)
 
c        Set flux and albedo weighting for ocean/seaice calls.
c
      IF (ITYPE.EQ.1 .OR. ITYPE.EQ.2) THEN
        DO 40 I=1,NPOI
          ZW(I) = 1.
   40   CONTINUE
      ELSE
        DO 50 I=1,NPOI
          ZW(I) = FWATER(I)
   50   CONTINUE
      ENDIF
c
c-------------------------
      IF (ITYPE.EQ.1) THEN
c-------------------------
c
c          Arrays (non-layer,...) in comveg[xy]
c
!DIR$ IVDEP
        DO 100 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
c         comveg:
          TUXY(JI,JJ) = TU(I)
          TSXY(JI,JJ) = TS(I)
          TLXY(JI,JJ) = TL(I)
          TLSUBXY(JI,JJ) = TLSUB(I)
          WLIQUXY(JI,JJ) = WLIQU(I)
          WLIQSXY(JI,JJ) = WLIQS(I)
          WLIQLXY(JI,JJ) = WLIQL(I)
          WSNOUXY(JI,JJ) = WSNOU(I)
          WSNOSXY(JI,JJ) = WSNOS(I)
          WSNOLXY(JI,JJ) = WSNOL(I)
          T12XY(JI,JJ)   = T12(I)
          T34XY(JI,JJ)   = T34(I)
          Q12XY(JI,JJ)   = Q12(I)
          Q34XY(JI,JJ)   = Q34(I)
          FOG12XY(JI,JJ) = FOG12(I)
          FOG34XY(JI,JJ) = FOG34(I)
          HFOG12XY(JI,JJ) = HFOG12(I)
          HFOG34XY(JI,JJ) = HFOG34(I)
#ifdef ibis
          VMAXSUNUXY(JI,JJ,1) = VMAXSUNU(I,1)
          VMAXSUNUXY(JI,JJ,2) = VMAXSUNU(I,2)
          VMAXSUNLXY(JI,JJ,1) = VMAXSUNL(I,1)
          VMAXSUNLXY(JI,JJ,2) = VMAXSUNL(I,2)
          VMAXSHAUXY(JI,JJ,1) = VMAXSHAU(I,1)
          VMAXSHAUXY(JI,JJ,2) = VMAXSHAU(I,2)
          VMAXSHALXY(JI,JJ,1) = VMAXSHAL(I,1)
          VMAXSHALXY(JI,JJ,2) = VMAXSHAL(I,2)
 
          CISUNUDXY(JI,JJ) = CISUNUD(I)
          CISUNUCXY(JI,JJ) = CISUNUC(I)
          CISUNL3XY(JI,JJ) = CISUNL3(I)
          CISUNL4XY(JI,JJ) = CISUNL4(I)
          CISHAUDXY(JI,JJ) = CISHAUD(I)
          CISHAUCXY(JI,JJ) = CISHAUC(I)
          CISHAL3XY(JI,JJ) = CISHAL3(I)
          CISHAL4XY(JI,JJ) = CISHAL4(I)
 
          CSSUNUDXY(JI,JJ) = CSSUNUD(I)
          CSSUNUCXY(JI,JJ) = CSSUNUC(I)
          CSSUNL3XY(JI,JJ) = CSSUNL3(I)
          CSSUNL4XY(JI,JJ) = CSSUNL4(I)
          CSSHAUDXY(JI,JJ) = CSSHAUD(I)
          CSSHAUCXY(JI,JJ) = CSSHAUC(I)
          CSSHAL3XY(JI,JJ) = CSSHAL3(I)
          CSSHAL4XY(JI,JJ) = CSSHAL4(I)
 
          CONFRACXY(JI,JJ) = CONFRAC(I)
          C4FRACXY (JI,JJ) = C4FRAC (I)
 
          ANUDXY(JI,JJ) = ANUD(I)
          ANUCXY(JI,JJ) = ANUC(I)
          ANL3XY(JI,JJ) = ANL3(I)
          ANL4XY(JI,JJ) = ANL4(I)
 
          AGUDXY(JI,JJ) = AGUD(I)
          AGUCXY(JI,JJ) = AGUC(I)
          AGL3XY(JI,JJ) = AGL3(I)
          AGL4XY(JI,JJ) = AGL4(I)
 
          GSUDXY(JI,JJ) = GSUD(I)
          GSUCXY(JI,JJ) = GSUC(I)
          GSL3XY(JI,JJ) = GSL3(I)
          GSL4XY(JI,JJ) = GSL4(I)
#endif
  100   CONTINUE
c
c          Arrays (stories) in comveg[xy]
c
        DO 110 IV=1,2
!DIR$ IVDEP
          DO 112 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            LAIXY(JI,JJ,IV) = LAI(I,IV)
            SAIXY(JI,JJ,IV) = SAI(I,IV)
            ORIEVXY(JI,JJ,IV) = ORIEV(I,IV)
            ORIEHXY(JI,JJ,IV) = ORIEH(I,IV)
            ZTOPXY(JI,JJ,IV) = ZTOP(I,IV)
            ZBOTXY(JI,JJ,IV) = ZBOT(I,IV)
            DLEAFXY(JI,JJ,IV) = DLEAF(I,IV)
            DSTEMXY(JI,JJ,IV) = DSTEM(I,IV)
            FSTOMAXY(JI,JJ,IV) = FSTOMA(I,IV)
            RSMINXY(JI,JJ,IV) = RSMIN(I,IV)
            RSMAXXY(JI,JJ,IV) = RSMAX(I,IV)
            PARCONXY(JI,JJ,IV) = PARCON(I,IV)
            VPDCONXY(JI,JJ,IV) = VPDCON(I,IV)
            TEMCONXY(JI,JJ,IV) = TEMCON(I,IV)
            TEMOPTXY(JI,JJ,IV) = TEMOPT(I,IV)
  112     CONTINUE
  110   CONTINUE
c
c          Arrays (stories, wavebands) in comveg[xy]
c
        DO 120 IV=1,2
          DO 122 IB=1,NBAND
!DIR$ IVDEP
            DO 124 I=1,NPOI
              JI = LOCPOI(I,1)
              JJ = LOCPOI(I,2)
              RHOVEGXY(JI,JJ,IB,IV) = RHOVEG(I,IB,IV)
              TAUVEGXY(JI,JJ,IB,IV) = TAUVEG(I,IB,IV)
  124       CONTINUE
  122     CONTINUE
  120   CONTINUE
c
c          Arrays (stories, soil layers) in comveg[xy]
c
        DO 130 IV=1,2
          DO 132 K=1,NSOILAY
!DIR$ IVDEP
            DO 134 I=1,NPOI
              JI = LOCPOI(I,1)
              JJ = LOCPOI(I,2)
              FROOTXY(JI,JJ,K,IV) = FROOT(I,K,IV)
  134       CONTINUE
  132     CONTINUE
  130   CONTINUE
c
c----------
      ENDIF
c----------
 
c-----------------------------------------
      IF (ITYPE.EQ.1 .OR. ITYPE.EQ.2) THEN
c-----------------------------------------
c
c          Arrays (non-layer,...) in comsno[xy], comsoi[x,y], commudxy
c
!DIR$ IVDEP
        DO 150 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
#ifdef isotrac
          fvcxy_i(ji,jj) = avapa(i)
#endif
c         comsno:
          FIXY(JI,JJ) = FI(I)
          TIXY(JI,JJ) = TI(I)
c         comsoi:
          WPUDXY(JI,JJ)   = WPUD(I)
          WIPUDXY(JI,JJ)  = WIPUD(I)
          Z0SOIXY(JI,JJ)  = Z0SOI(I)
          ALBSAVXY(JI,JJ) = ALBSAV(I)
          ALBSANXY(JI,JJ) = ALBSAN(I)
          TGXY(JI,JJ)     = TG(I)
c         commudxy:
          UPTOTUXY(JI,JJ) = UPTOTU(I)
          UPTOTLXY(JI,JJ) = UPTOTL(I)
  150   CONTINUE
c
c          Arrays (soil layers) in comsoi[xy] and commudxy
c
        DO 160 K=1,NSOILAY
!DIR$ IVDEP
          DO 162 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
c           comsoixy:
            ITEXXY(JI,JJ,K) = ITEX(I,K)
            TSOIXY(JI,JJ,K) = TSOI(I,K)
            WSOIXY(JI,JJ,K) = WSOI(I,K)
            WISOIXY(JI,JJ,K) = WISOI(I,K)
c           commudxy:
            POROSXY(JI,JJ,K)= POROS(I,K)
  162     CONTINUE
  160   CONTINUE
c
c          Arrays (snow layers) in comsno[xy]
c
        DO 170 K=1,NSNOLAY
!DIR$ IVDEP
          DO 172 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            TSNOXY(JI,JJ,K) = TSNO(I,K)
            HSNOXY(JI,JJ,K) = HSNO(I,K)
  172     CONTINUE
  170   CONTINUE
c
c---------
      ENDIF
c---------
 
c-------------------------
      IF (ITYPE.EQ.3) THEN
c-------------------------
 
c          Arrays (non-layer) in comocn[xy], comice[xy], comsno[xy]
 
!DIR$ IVDEP
        DO 200 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
c         comocn:
          TOCNXY(JI,JJ)  = TOCN(I)
          SOCNXY(JI,JJ)  = SOCN(I)
          TFREEZXY(JI,JJ)= TFREEZ(I)
          HOCNXY(JI,JJ)  = HOCN(I)
c         for mom2,4:
          FHOCNXY(JI,JJ)   = FHOCN(I)
          FQOCNXY(JI,JJ)   = FQOCN(I)
          FSOLOCNXY(JI,JJ) = FSOLOCN(I)
          FIROCNXY(JI,JJ)  = FIROCN(I)
          FPRLOCNXY(JI,JJ) = FPRLOCN(I)
          FPRFOCNXY(JI,JJ) = FPRFOCN(I)
          FSALOCNXY(JI,JJ) = FSALOCN(I)
          PRESOCNXY(JI,JJ) = PRESOCN(I) 

          FTOCNXY(JI,JJ) = FTOCN(I)
          FWOCNXY(JI,JJ) = FWOCN(I)
          FUOCNXY(JI,JJ) = FUOCN(I)
          FVOCNXY(JI,JJ) = FVOCN(I)
          QFLUXXY(JI,JJ) = QFLUX(I)
          FUICEXY(JI,JJ) = FUICE(I)
          FVICEXY(JI,JJ) = FVICE(I)
c         comice:
          FICEXY(JI,JJ)  = FICE(I)
          BRINEXY(JI,JJ) = BRINE(I)
          TG_ICEXY(JI,JJ)= TG_ICE(I)
c         comsno:
          FI_ICEXY(JI,JJ)= FI(I)
          TI_ICEXY(JI,JJ)= TI(I)
#ifdef isotrac
          xoxy_o(ji,jj)  = xo_o(i)
          xoxy_h(ji,jj)  = xo_h(i)
          qoxy_i(ji,jj)  = qo_i(i)
          zfxy_i(ji,jj)  = zf_i(i)
          ficexy_i(ji,jj)= fice_i(i)
          fvoxy_i(ji,jj) = -fvapo(i)
          fvixy_i(ji,jj) = fvi_i(i)
#endif
  200   CONTINUE
c
c          Arrays (sea-ice layers) in comice[xy]
c
        DO 240 K=1,NICELAY
!DIR$ IVDEP
          DO 242 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            HICEXY(JI,JJ,K) = HICE(I,K)
            TICEXY(JI,JJ,K) = TICE(I,K)
  242     CONTINUE
  240   CONTINUE
c
c          Arrays (snow layers) in comsno[xy]
c
        DO 250 K=1,NSNOLAY
!DIR$ IVDEP
          DO 252 I=1,NPOI
            JI = LOCPOI(I,1)
            JJ = LOCPOI(I,2)
            TSNO_ICEXY(JI,JJ,K) = TSNO(I,K)
            HSNO_ICEXY(JI,JJ,K) = HSNO(I,K)
  252     CONTINUE
  250   CONTINUE
 
c----------
      ENDIF
c----------
 
c
!DIR$ IVDEP
      DO 300 I=1,NPOI
        JI = LOCPOI(I,1)
        JJ = LOCPOI(I,2)
 
c       Arrays in comveg[xy]
 
        IVEGXY(JI,JJ) = IVEG(I)
        FUXY(JI,JJ) = FU(I)
        FLXY(JI,JJ) = FL(I)
 
c       Arrays in commug[xy]
 
        COSZENXY(JI,JJ) = COSZEN(I)
        FIRAXY(JI,JJ) = FIRA(I)
        UXXY(JI,JJ) = UX(I)
        UYXY(JI,JJ) = UY(I)
        TAXY(JI,JJ) = TA(I)
        QAXY(JI,JJ) = QA(I)
        RAINAXY(JI,JJ) = RAINA(I)
        SNOWAXY(JI,JJ) = SNOWA(I)
        PSURFXY(JI,JJ) = PSURF(I)
        AFIREXY(JI,JJ) = AFIRE(I)*ZW(I) + AFIREXY(JI,JJ)*(1.-ZW(I))
        ASENAXY(JI,JJ) = ASENA(I)*ZW(I) + ASENAXY(JI,JJ)*(1.-ZW(I))
        AVAPAXY(JI,JJ) = AVAPA(I)*ZW(I) + AVAPAXY(JI,JJ)*(1.-ZW(I))
        ATAUXXY(JI,JJ) = ATAUX(I)*ZW(I) + ATAUXXY(JI,JJ)*(1.-ZW(I))
        ATAUYXY(JI,JJ) = ATAUY(I)*ZW(I) + ATAUYXY(JI,JJ)*(1.-ZW(I))
        LMASKXY(JI,JJ) = LMASK(I)
        FWATERXY(JI,JJ)= FWATER(I)
 
c       Arrays in commud[xy]
 
        GINVAPXY(JI,JJ) = GINVAP(I)*ZW(I) + GINVAPXY(JI,JJ)*(1.-ZW(I))
        GSUVAPXY(JI,JJ) = GSUVAP(I)*ZW(I) + GSUVAPXY(JI,JJ)*(1.-ZW(I))
        GTRANSUXY(JI,JJ)= GTRANSU(I)*ZW(I)+ GTRANSUXY(JI,JJ)*(1.-ZW(I))
        GTRANSLXY(JI,JJ)= GTRANSL(I)*ZW(I)+ GTRANSLXY(JI,JJ)*(1.-ZW(I))
        GRUNOFXY(JI,JJ) = GRUNOF(I)*ZW(I) + GRUNOFXY(JI,JJ)*(1.-ZW(I))
        GDRAINXY(JI,JJ) = GDRAIN(I)*ZW(I) + GDRAINXY(JI,JJ)*(1.-ZW(I))
        GDRAINHXY(JI,JJ)= GDRAINH(I)*ZW(I) + GDRAINHXY(JI,JJ)*(1.-ZW(I))
        DO K=1,NSOILAY+1
          WFLOXY(JI,JJ,K) = WFLO(I,K)*ZW(I)+ WFLOXY(JI,JJ,K)*(1.-ZW(I))
        ENDDO
        GABLATXY(JI,JJ) = GABLAT(I)*ZW(I) + GABLATXY(JI,JJ)*(1.-ZW(I))
        GADJUSTXY(JI,JJ)= GADJUST(I)*ZW(I)+ GADJUSTXY(JI,JJ)*(1.-ZW(I))
        TS2XY(JI,JJ)    = TS2(I)   *ZW(I) + TS2XY(JI,JJ)   *(1.-ZW(I))
#if defined (rcm) 
        ts10xy(ji,jj)   = ts10(i)  *zw(i) + ts10xy(ji,jj)  *(1.-zw(i))
        qs10xy(ji,jj)   = qs10(i)  *zw(i) + qs10xy(ji,jj)  *(1.-zw(i))
        us10xy(ji,jj)   = us10(i)  *zw(i) + us10xy(ji,jj)  *(1.-zw(i))
#endif
        Z0LENXY(JI,JJ)  = Z0LEN(I) *ZW(I) + Z0LENXY(JI,JJ) *(1.-ZW(I))
        GINFILXY(JI,JJ) = GINFIL(I)*ZW(I) + GINFILXY(JI,JJ)*(1.-ZW(I))
        GIRRIGXY(JI,JJ) = GIRRIG(I)*ZW(I) + GIRRIGXY(JI,JJ)*(1.-ZW(I))
 
  300 CONTINUE
c
c        Arrays (wavebands) in commug[xy]
c
      DO 310 IB=1,NBAND
!DIR$ IVDEP
        DO 312 I=1,NPOI
          JI = LOCPOI(I,1)
          JJ = LOCPOI(I,2)
          SOLADXY(JI,JJ,IB) = SOLAD(I,IB)
          SOLAIXY(JI,JJ,IB) = SOLAI(I,IB)
          ASURDXY(JI,JJ,IB) = ASURD(I,IB)      *    ZW(I)
     *                      + ASURDXY(JI,JJ,IB)*(1.-ZW(I))
          ASURIXY(JI,JJ,IB) = ASURI(I,IB)      *    ZW(I)
     *                      + ASURIXY(JI,JJ,IB)*(1.-ZW(I))
  312   CONTINUE
  310 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
#if defined (rcm) || defined (hydmodel)
      SUBROUTINE READLSX (IU, IVERSRCMIN)
#else
      SUBROUTINE READLSX (IU, IUOUT, VERSRESIN, TYPERESIN,
     *                    LONRESLSX, LATRESLSX, NLONLSX, NLATLSX,
     *                    VERSRESIN_ISOT, NSREST_ISOT, cduma)
#endif
c
c        Reads LSX info from ccm restart file (called from ccm RSDS).
c        Use params (NLON,etc) for array sizes (instead of NLON,etc),
c        and pass output unit#, so READLSX can be called before INISURF.
c
c        IU       = unit number of ccm restart file           (supplied)
c        IUOUT    = unit number for messages                  (supplied)
c        VERSRESIN= version number of restart file being read (supplied)
c        TYPERESIN= machine-representation flag (CRAY/IEEE)   (supplied)
c        LONRESLSX = restart-file longitude resolution if user-specified
c                    (or 0 if not, and set to lsx NLON here)  (modified)
c        LATRESLSX = as lonreslsx, for restart file lat resol (modified)
c        NLONLSX   = surface model longitude resolution (=NLON,returned)
c        NLATLSX   = surface model latitude  resolution (=NLAT,returned)
c
#include <compar.h>
#include <combudl.h>
#include <comeve.h>
#include <comgrd.h>
#include <comice.h>
#if defined gcm && defined rcm
#include <COMLUN.h>
#endif
#include <comicexy.h>
#include <commudxy.h>
#include <commugxy.h>
#include <comocn.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comvegxy.h>

      character*8 cduma
      integer*8 i8a, i8b
 
#if defined (rcm) || defined (hydmodel)
      CHARACTER*8 TYPERESIN
#else
      CHARACTER*(*) TYPERESIN
#endif
      DIMENSION ZHTOT(NLON,NLAT)

#if defined (rcm) || defined (hydmodel)
      iuout = ioterm
      read (iu) versresin   ! versres for rcm
      write (iuout,'(/a/a,f10.3)')
     *   ' *** Reading LSX restart file section',
     *   '     versresin: ',versresin
      lonreslsx = 0
      latreslsx = 0
#endif
 
c        Set expected restart-file resols (namelist param horreslsx[])
c        to model resols if not entered. Must be done before first
c        lsx call to reader.
 
      IF (LONRESLSX.EQ.0) LONRESLSX = NLON
      IF (LATRESLSX.EQ.0) LATRESLSX = NLAT
 
c        Return surface model resolutions to agcm (used for info only)
 
      NLONLSX = NLON
      NLATLSX = NLAT
 
c        lsx header rec (moved to before surface fields to allow
c        for reading different resolutions without using horreslsx
c        in the future)
 
      IF (VERSRESIN.GE.2.01) THEN
#if defined (gcm) || defined (hydmodel)
        if (cduma(1:4).eq.'CRAY') then
          READ (IU, ERR=6000, IOSTAT=IOS) i8a, i8b
          NLONZ = i8a
          NLATZ = i8b
        else
          READ (IU, ERR=6000, IOSTAT=IOS) NLONZ,NLATZ
        endif
        CALL CONVERTYPE (NLONZ, 1, 'I', TYPERESIN)
        CALL CONVERTYPE (NLATZ, 1, 'I', TYPERESIN)
        WRITE(IUOUT,5010) NLONZ,NLATZ
 5010   FORMAT(9x,'surface lsx resolutions  :',2i6)
#else
        READ (IU,'(1X,I25)', ERR=6100, IOSTAT=IOS) NLONZ,NLATZ
#endif
        IF (NLONZ.NE.LONRESLSX .OR. NLATZ.NE.LATRESLSX) THEN
          WRITE(IUOUT,5020) NLONZ,LONRESLSX,NLATZ,LATRESLSX
 5020     FORMAT (
     *    //' Surface (LSX) resolutions of initial/restart file do not',
     *      ' match the expected'
     *    / ' values (i.e., the current model''s, or specified by',
     *      ' namelist horreslsx[])'
     *    /4x,'Initial/restart file NLON:',i8,'   expected NLON:',i8
     *    /4x,'Initial/restart file NLAT:',i8,'   expected NLAT:',i8
     *    /4x,'Stop in readlsx.'/ )
          CALL ENDRUN (-1)
        ENDIF
      ENDIF
 
c
c        combudl (lsx global budget, subr budall)
c
      IF (VERSRESIN.GE.2.00) THEN
#if defined (gcm) || defined (hydmodel)
        READ (IU, ERR=6100, IOSTAT=IOS) TMPBUDL
        CALL CONVERTYPE (TMPBUDL, NBUDGL, 'R', TYPERESIN)
#else
        READ (IU,'(1X,E25.15)', ERR=6100, IOSTAT=IOS) TMPBUDL
#endif
      ENDIF
c
c        ifile for rcm history file list
c
#if defined (rcm)
      read (iu) ifile
#endif
c
c        commug
c
      IF (VERSRESIN.EQ.1.00) THEN
        DO 100 IDUM = 1,1+2*NBAND             ! coszenxy,soladxy,solaixy
          CALL READER (IU, WORKL,  NLON,NLAT,1, 2)
  100   CONTINUE
      ENDIF
 
      IF (VERSRESIN.GE.3.00) THEN
        CALL READER (IU, WORKL, NLON,NLAT,1, 2)
        DO 104 JJ=1,NLAT
          DO 106 JI=1,NLON
            LMASKXY(JI,JJ) = NINT (WORKL(JI,JJ))
  106     CONTINUE
  104   CONTINUE
      ELSE
c       Fudge to read old Cray restart files with integer lmaskxy.
c       Relies on any restart version < 3.00 => CRAY type, so
c       can read directly with no conversion (...-2), then on Sgi
c       convert from CRAY to IEEE integers.
 
        CALL READER (IU, LMASKXY, NLON,NLAT,1, -2)
#ifdef gcm
#ifdef sgi
        CALL CONVERTYPE (LMASKXY, NLON*NLAT, 'I', TYPERESIN)
#endif
#endif
      ENDIF
 
      IF (VERSRESIN.GE.2.04) THEN
        CALL READER (IU, FWATERXY, NLON,NLAT,1, 2)
      ELSE
c       For pre-lake versions, set open water fraction to 0 or 1
c       depending on LMASK, and change any lmask=4 to 3
        DO 110 JJ=1,NLAT
          DO 112 JI=1,NLON
            IF (LMASKXY(JI,JJ).LE.2) THEN
              FWATERXY(JI,JJ) = 0.
            ELSE
              LMASKXY(JI,JJ) = 3
              FWATERXY(JI,JJ) = 1.
            ENDIF
  112     CONTINUE
  110   CONTINUE
      ENDIF
c
c        comocn
c
      CALL READER (IU, TOCNXY,  NLON,NLAT,1, 2)
      IF (VERSRESIN.GE.3.10) THEN
        CALL READER (IU, SOCNXY,  NLON,NLAT,1, 2)
      ELSE
        CALL RESETR (SOCNXY,    NLON*NLAT,  .035)
      ENDIF
      CALL READER (IU, HOCNXY,  NLON,NLAT,1, 2)
      IF (VERSRESIN.EQ.1.00) THEN
        DO 150 IDUM = 1,8  ! solo1xy,f[t,w,u,v]ocnxy,qfluxxy,f[u,v]icexy
          CALL READER (IU, WORKL,  NLON,NLAT,1, 2)
  150   CONTINUE
      ENDIF
      IF (VERSRESIN.LE.2.03) THEN !will be read below(comice) for > 2.03
        CALL READER (IU, FICEXY,  NLON,NLAT,1, 2)
        CALL READER (IU, BRINEXY, NLON,NLAT,1, 2)
      ENDIF
      IF (VERSRESIN.GE.2.00) THEN
        CALL READER (IU, FUICEXYM, NLON,NLAT,1, 2)
        CALL READER (IU, FVICEXYM, NLON,NLAT,1, 2)
        CALL READER (IU, FUOCNXYM, NLON,NLAT,1, 2)
        CALL READER (IU, FVOCNXYM, NLON,NLAT,1, 2)
      ELSE
        CALL ZERO (FUICEXYM, NLON*NLAT)
        CALL ZERO (FVICEXYM, NLON*NLAT)
        CALL ZERO (FUOCNXYM, NLON*NLAT)
        CALL ZERO (FVOCNXYM, NLON*NLAT)
      ENDIF
c
c        comsno
c
      CALL READER (IU, TSNOXY, NLON,NLAT,NSNOLAY, 2)
      CALL READER (IU, HSNOXY, NLON,NLAT,NSNOLAY, 2)
      CALL READER (IU, FIXY,   NLON,NLAT,1,       2)
      IF (VERSRESIN.GE.2.04) THEN
        CALL READER (IU, TIXY,       NLON,NLAT,1,       2)
        CALL READER (IU, TSNO_ICEXY, NLON,NLAT,NSNOLAY, 2)
        CALL READER (IU, HSNO_ICEXY, NLON,NLAT,NSNOLAY, 2)
        CALL READER (IU, FI_ICEXY,   NLON,NLAT,1,       2)
        CALL READER (IU, TI_ICEXY,   NLON,NLAT,1,       2)
      ELSE
        DO 180 JJ=1,NLAT
          DO 182 JI=1,NLON
 
            IF (FWATERXY(JI,JJ).EQ.0.) THEN
              DO 184 JK=1,NSNOLAY
                TSNO_ICEXY(JI,JJ,JK) = TMELT
                HSNO_ICEXY(JI,JJ,JK) = 0.
  184         CONTINUE
              FI_ICEXY(JI,JJ) = 0.
            ELSE
              DO 186 JK=1,NSNOLAY
                TSNO_ICEXY(JI,JJ,JK) = TSNOXY(JI,JJ,JK)
                HSNO_ICEXY(JI,JJ,JK) = HSNOXY(JI,JJ,JK)
                TSNOXY(JI,JJ,JK) = 0.
                HSNOXY(JI,JJ,JK) = 0.
  186         CONTINUE
              FI_ICEXY(JI,JJ) = FIXY(JI,JJ)
              FIXY(JI,JJ) = 0.
            ENDIF
 
            TIXY(JI,JJ) = TSNOXY(JI,JJ,1)
            TI_ICEXY(JI,JJ) = TSNO_ICEXY(JI,JJ,1)
 
  182     CONTINUE
  180   CONTINUE
      ENDIF
c
c        comsoi
c
      IF (VERSRESIN.LE.3.02) THEN
        CALL ZERO (ZHTOT, NLON*NLAT)
        DO 187 K=1,NSOILAY
          CALL READER (IU, WORKL, NLON,NLAT,1, 2)
          DO 188 JJ=1,NLAT
            DO 189 JI=1,NLON
              ZHTOT(JI,JJ) = ZHTOT(JI,JJ) + WORKL(JI,JJ)
  189       CONTINUE
  188     CONTINUE
  187   CONTINUE
      ENDIF
      CALL READER (IU, TSOIXY,  NLON,NLAT,NSOILAY, 2)
      CALL READER (IU, WSOIXY,  NLON,NLAT,NSOILAY, 2)
      CALL READER (IU, WISOIXY, NLON,NLAT,NSOILAY, 2)
      IF (VERSRESIN.GE.1.08) THEN
        CALL READER (IU, WPUDXY,  NLON,NLAT,1, 2)
        CALL READER (IU, WIPUDXY, NLON,NLAT,1, 2)
      ELSE
        CALL ZERO (WPUDXY, NLON*NLAT)
        CALL ZERO (WIPUDXY, NLON*NLAT)
      ENDIF
 
      IF (VERSRESIN.GE.2.04) THEN
        CALL READER (IU, TGXY,  NLON,NLAT,1,       2)
      ELSE
        DO 190 JJ=1,NLAT
          DO 192 JI=1,NLON
            TGXY(JI,JJ) = TSOIXY(JI,JJ,1)
 192      CONTINUE
 190    CONTINUE
      ENDIF
c
c        comice
c
      IF (VERSRESIN.GE.2.04) THEN
        CALL READER (IU, HICEXY,  NLON,NLAT,NICELAY, 2)
        CALL READER (IU, TICEXY,  NLON,NLAT,NICELAY, 2)
        CALL READER (IU, FICEXY,  NLON,NLAT,1, 2)
        CALL READER (IU, BRINEXY, NLON,NLAT,1, 2)
#ifdef rcm
        call zero (brinexy, nlon*nlat)  ! Nan problem in rcm rstrt files
                                        ! fbrinmax,fbrinsto=0 (no brine)
#endif

        CALL READER (IU, TG_ICEXY,NLON,NLAT,1, 2)
      ELSE
 
c         If no open water, ensure hocnxy = 0 (may be overridden
c         by data input files and iniocn). If some sea ice,
c         adapt hice,tice to new vertical resolution (nsoilay->nicelay)
c         and set tsoixy,etc, to null.
 
        DO 200 JJ=1,NLAT
          DO 202 JI=1,NLON
 
            IF (FWATERXY(JI,JJ).EQ.0.) THEN
              HOCNXY(JI,JJ) = 0.
              TOCNXY(JI,JJ) = TMELT
            ENDIF
 
            IF (FWATERXY(JI,JJ).EQ.0. .OR. FICEXY(JI,JJ).EQ.0.) THEN
 
              DO 210 JK=1,NICELAY
                HICEXY(JI,JJ,JK) = 0.
                TICEXY(JI,JJ,JK) = TMELT
  210         CONTINUE
              FICEXY(JI,JJ) = 0.
              BRINEXY(JI,JJ) = 0.
 
            ELSE
 
              HICEXY(JI,JJ,1) = MIN (.05, 0.5*ZHTOT(JI,JJ))
              TICEXY(JI,JJ,1) = TSOIXY(JI,JJ,1)
              DO 222 JK=2,NICELAY
                HICEXY(JI,JJ,JK) = (ZHTOT(JI,JJ)-HICEXY(JI,JJ,1))
     *                             / (NICELAY-1)
                JK1 = MAX (1, NINT(NSOILAY*FLOAT(JK-1)/FLOAT(NICELAY)))
                JK2 = MAX (1, NINT(NSOILAY*FLOAT(JK  )/FLOAT(NICELAY)))
                TICEXY(JI,JJ,JK) = 0.
                DO 224 JKK = JK1,JK2
                  TICEXY(JI,JJ,JK) = TICEXY(JI,JJ,JK)
     *                             + TSOIXY(JI,JJ,JKK) / (JK2-JK1+1)
  224           CONTINUE
  222         CONTINUE
 
              DO 226 JK=1,NSOILAY
                TSOIXY(JI,JJ,JK) = 0.
                WSOIXY(JI,JJ,JK) = 0.
                WISOIXY(JI,JJ,JK) = 0.
  226         CONTINUE
              WPUDXY(JI,JJ) = 0.
              WIPUDXY(JI,JJ) = 0.
 
            ENDIF
 
            TG_ICEXY(JI,JJ) = TICEXY(JI,JJ,1)
 
  202     CONTINUE
  200   CONTINUE
 
      ENDIF
c
c        comveg
c
      CALL READER (IU, TUXY,    NLON,NLAT,1, 2)
      CALL READER (IU, TSXY,    NLON,NLAT,1, 2)
      CALL READER (IU, TLXY,    NLON,NLAT,1, 2)
      CALL READER (IU, TLSUBXY, NLON,NLAT,1, 2)
      CALL READER (IU, WLIQUXY, NLON,NLAT,1, 2)
      CALL READER (IU, WLIQSXY, NLON,NLAT,1, 2)
      CALL READER (IU, WLIQLXY, NLON,NLAT,1, 2)
      CALL READER (IU, WSNOUXY, NLON,NLAT,1, 2)
      CALL READER (IU, WSNOSXY, NLON,NLAT,1, 2)
      CALL READER (IU, WSNOLXY, NLON,NLAT,1, 2)
      CALL READER (IU, T12XY,   NLON,NLAT,1, 2)
      CALL READER (IU, T34XY,   NLON,NLAT,1, 2)
      CALL READER (IU, Q12XY,   NLON,NLAT,1, 2)
      CALL READER (IU, Q34XY,   NLON,NLAT,1, 2)
      CALL READER (IU, FOG12XY, NLON,NLAT,1, 2)
      CALL READER (IU, FOG34XY, NLON,NLAT,1, 2)
      CALL READER (IU, HFOG12XY,NLON,NLAT,1, 2)
      CALL READER (IU, HFOG34XY,NLON,NLAT,1, 2)
      IF (VERSRESIN.EQ.1.00) THEN
        DO 300 IDUM = 1,6                      ! sol[u,s,l,g,i]1xy,fi1xy
          CALL READER (IU, WORKL, NLON,NLAT,1, 2)
  300   CONTINUE
      ENDIF
#ifdef ibis
      IF (VERSRESIN.GE.3.05) THEN
        CALL READER (IU, VMAXSUNUXY,    NLON,NLAT,2, 2)
        CALL READER (IU, VMAXSUNLXY,    NLON,NLAT,2, 2)
        CALL READER (IU, VMAXSHAUXY,    NLON,NLAT,2, 2)
        CALL READER (IU, VMAXSHALXY,    NLON,NLAT,2, 2)
        IF (VERSRESIN.GE.3.06) THEN
          CALL READER (IU, UPTOTUXY,    NLON,NLAT,1, 2)
          CALL READER (IU, UPTOTLXY,    NLON,NLAT,1, 2)
          CALL READER (IU, T15RUNXY,    NLON,NLAT,1, 2)
          CALL READER (IU, AGDD0XY,     NLON,NLAT,1, 2)
          CALL READER (IU, AGDD5XY,     NLON,NLAT,1, 2)
        ELSE
          CALL RESETR (UPTOTUXY,    NLON*NLAT,   1.E20)
          DO IDUM=1,3                    ! tdailxy,tcoldprevxy,pmonminxy
            CALL READER (IU, WORKL,     NLON,NLAT,1, 2)
          ENDDO
        ENDIF
      ELSE
c       Flag for later detection of no IBIS info (see writlsx, inilsx)
        CALL RESETR (VMAXSUNUXY,    NLON*NLAT*2, 1.E20)
        CALL RESETR (UPTOTUXY,      NLON*NLAT,   1.E20)
        IF (VERSRESIN.GE.3.01) THEN
c         v2.0.a vmax[sun,sha]xy
          DO 340 IDUM = 1,4
            CALL READER (IU, WORKL, NLON,NLAT,1, 2)
  340     CONTINUE
        ENDIF
      ENDIF
#else
      IF (VERSRESIN.GE.3.06) THEN
        DO IDUM = 1,13
          CALL READER (IU, WORKL, NLON,NLAT,1, 2)
        ENDDO
      ELSE IF (VERSRESIN.GE.3.05) THEN
        DO 350 IDUM = 1,11
          CALL READER (IU, WORKL, NLON,NLAT,1, 2)
  350   CONTINUE
      ELSE IF (VERSRESIN.GE.3.01) THEN
c       v2.0.a vmax[sun,sha]xy
        DO 352 IDUM = 1,4
          CALL READER (IU, WORKL, NLON,NLAT,1, 2)
  352   CONTINUE
      ENDIF
#endif
c
c        comaccxy
c
      IF (VERSRESIN.LE.1.02) THEN
        DO 400 IDUM = 1,2+2*3              ! tocaxy,hicaxy,wsoaxy,wisaxy
          CALL READER (IU, WORKL, NLON,NLAT,1, 2)
  400   CONTINUE
      ENDIF
c
c        comeve
c
      IF (VERSRESIN.GE.2.05) THEN
        CALL READER (IU, TEVE_D, NLON,NLAT,1, 2)
        CALL READER (IU, PEVE_D, NLON,NLAT,1, 2)
        CALL READER (IU, REVE_D, NLON,NLAT,1, 2)
      ELSE
        DO 500 JJ=1,NLAT
          DO 502 JI=1,NLON
            TEVE_D(JI,JJ) = 15.       ! deg C
            PEVE_D(JI,JJ) = 0.        ! mm per day
            REVE_D(JI,JJ) = 50.       ! %
  502     CONTINUE
  500   CONTINUE
      ENDIF
c
c        commud
c
      IF (VERSRESIN.GE.3.09) THEN
        CALL READER (IU, GIRRIGXY,  NLON,NLAT,1, 2)
      ELSE
        CALL ZERO (GIRRIGXY, NLON*NLAT)
      ENDIF

c
c     write (*,8881) ((nint(teve_d(i,j)),i=1,nlon),j=1,nlat)
c8881 format(/'readlsx:  FIRST TEVE_D :' / (180i3))
c     write (*,8882) ((nint(peve_d(i,j)),i=1,nlon),j=1,nlat)
c8882 format(/'readlsx:  FIRST PEVE_D :' / (180i3))
c
c     v2.0.a wisconsin twisc
      IF (VERSRESIN.GE.3.02 .AND. VERSRESIN.LT.3.05) THEN
        CALL READER (IU, WORKL, NLON,NLAT,1, 2)
      ENDIF
c
c     v2.0.a wisconsin total soil stress
      IF (VERSRESIN.GE.3.03 .AND. VERSRESIN.LT. 3.05) THEN
        CALL READER (IU, WORKL, NLON,NLAT,1, 2)
        CALL READER (IU, WORKL, NLON,NLAT,1, 2)
      ENDIF

c>>>>>>>>>>>>>
#ifdef isotrac
c>>>>>>>>>>>>>
      if (versresin_isot.ge.1.00) then
        write(6,6500)
 6500   format (/4x,'Reading surface isotopic fields...')
        call reader (iu, hi_i,   nlon,nlat,1, 2)
        call reader (iu, ri_o,   nlon,nlat,1, 2)
        call reader (iu, ri_h,   nlon,nlat,1, 2)
        call reader (iu, hct_i,  nlon,nlat,1, 2)
        call reader (iu, rct_o,  nlon,nlat,1, 2)
        call reader (iu, rct_h,  nlon,nlat,1, 2)
        call reader (iu, hcb_i,  nlon,nlat,1, 2)
        call reader (iu, rcb_o,  nlon,nlat,1, 2)
        call reader (iu, rcb_h,  nlon,nlat,1, 2)
        call reader (iu, hl_i,   nlon,nlat,1, 2)
        call reader (iu, ro_o,   nlon,nlat,1, 2)
        call reader (iu, ro_h,   nlon,nlat,1, 2)
      endif
c>>>>>
#else
c>>>>>
      if (versresin_isot.ge.1.00) then
        write(6,6504)
 6504   format (/4x,'Skipping surface isotopic  fields...')
        do iskip=1,12
          call reader (iu, workl, nlon, nlat, 1, 2)
        enddo
      endif
c>>>>>
#endif
c>>>>>

c        lsx trailer rec

#ifdef gcm
#if defined (rcm) || defined (hydmodel)

c       to detect need-to-append to history file. For gcm, nwrit_a and
c       nuoff_a are written in wsds (agcm). So for rcm, need them here.
        if (iversrcmin.ge.1) then
          read (iu) nwrit_a, nuoff_a 
        else
          nwrit_a = 0
          nuoff_a = 0
        endif

c       Last word check for rcm or hydmodel (done in rsds,wsds for gcm) 
        read (iu) versresinx
        if (versresinx.ne.versresin) then
          write (iuout,'(/a/a,f10.3/a,f10.3)')
     *      '*** Error reading restart file, lsx section:',
     *      '    first-word versresin = ', versresin,
     *      '    last -word versresin = ', versresinx
          call endrun (-1)
        else
          write (iuout,'(a,f10.3,a/)')
     *        '     versresin: ',versresinx,' (last word check)'
        endif
#else
      IF (VERSRESIN.LE.2.00) THEN
        READ (IU, ERR=6200, IOSTAT=IOS) NLONZ, NLATZ
        CALL CONVERTYPE (NLONZ, 1, 'I', TYPERESIN)
        CALL CONVERTYPE (NLATZ, 1, 'I', TYPERESIN)
        WRITE(IUOUT,5010) NLONZ,NLATZ
        IF (NLONZ.NE.LONRESLSX .OR. NLATZ.NE.LATRESLSX) THEN
          WRITE(IUOUT,5020) NLONZ,LONRESLSX,NLATZ,LATRESLSX
          CALL ENDRUN (-1)
        ENDIF
      ENDIF
#endif
#endif
c
      RETURN
c
c        Error messages
c
 6000 WRITE(IUOUT,6002) IOS
 6002 FORMAT(/' *** Error or EOF reading surface-model',
     *        ' header (resolution) record of initial/restart file.'
     *       /'     IOSTAT code = ',i8
     *       /'     Stop in READLSX.'/)
      CALL ENDRUN (-1)
c
 6100 WRITE(IUOUT,6102) IOS
 6102 FORMAT(/' *** Errr or EOF reading surface-model',
     *        ' budget record of initial/restart file.'
     *       /'     IOSTAT code = ',i8
     *       /'     Stop in READLSX.'/)
      CALL ENDRUN (-1)
c
 6150 WRITE(IUOUT,6152) IOS
 6152 FORMAT(/' *** Errr or EOF reading surface-model',
     *        ' accumulation number record of initial/restart file.'
     *       /'     IOSTAT code = ',i8
     *       /'     Stop in READLSX.'/)
      CALL ENDRUN (-1)
c
 6200 WRITE(IUOUT,6202) IOS
 6202 FORMAT(/' *** Error or EOF reading surface-model',
     *        ' trailer (resolution) record of initial/restart file.'
     *       /'     IOSTAT code = ',i8
     *       /'     Stop in READLSX.'/)
      CALL ENDRUN (-1)
c
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE WRITLSX (IU)
c
c        Writes LSX info to restart file (called from ccm WSDS)
c
c        IU        = unit number of ccm restart file
c
#include <compar.h>
#include <combudl.h>
#include <comeve.h>
#include <comgrd.h>
#include <comice.h>
#if defined (gcm) && defined (rcm)
#include <COMLUN.h>
#endif
#include <comicexy.h>
#include <commudxy.h>
#include <commugxy.h>
#include <comocnxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comvegxy.h>

#if defined (rcm)
      versreszz = 3.10   ! versres for rcm (should be current gcm value)
      write (iu) versreszz 
#elif defined (hydmodel)
      versreszz = 3.10   ! versres for hydmodel (" " " " " " " " " " ")
      write (iu) versreszz 
#endif

c        header (lsx resolutions) record
c
#if defined (gcm) || defined (hydmodel)
      WRITE (IU) NLON,NLAT
#else
      WRITE (IU,'(1X,I25)') NLON,NLAT
#endif
c
c        combudl (lsx global budget, subr budall)
c
#if defined (gcm) || defined (hydmodel)
      WRITE (IU) TMPBUDL
#else
      WRITE (IU,'(1X,E25.15)') TMPBUDL
#endif
c
c        ifile for rcm history file list
c
#if defined (rcm)
      write (iu) ifile
#endif
c
c        commug
c
      DO 100 JJ=1,NLAT
        DO 102 JI=1,NLON
          WORKL(JI,JJ) = LMASKXY(JI,JJ)
  102   CONTINUE
  100 CONTINUE
      CALL WRITER (IU, WORKL,   NLON,NLAT,1, 2)
      CALL WRITER (IU, FWATERXY,NLON,NLAT,1, 2)
c
c        comocn
c
      CALL WRITER (IU, TOCNXY,  NLON,NLAT,1, 2)
      CALL WRITER (IU, SOCNXY,  NLON,NLAT,1, 2)
      CALL WRITER (IU, HOCNXY,  NLON,NLAT,1, 2)
      CALL WRITER (IU, FUICEXYM,NLON,NLAT,1, 2)
      CALL WRITER (IU, FVICEXYM,NLON,NLAT,1, 2)
      CALL WRITER (IU, FUOCNXYM,NLON,NLAT,1, 2)
      CALL WRITER (IU, FVOCNXYM,NLON,NLAT,1, 2)
c
c        comsno
c
      CALL WRITER (IU, TSNOXY,     NLON,NLAT,NSNOLAY, 2)
      CALL WRITER (IU, HSNOXY,     NLON,NLAT,NSNOLAY, 2)
      CALL WRITER (IU, FIXY,       NLON,NLAT,1,       2)
      CALL WRITER (IU, TIXY,       NLON,NLAT,1,       2)
      CALL WRITER (IU, TSNO_ICEXY, NLON,NLAT,NSNOLAY, 2)
      CALL WRITER (IU, HSNO_ICEXY, NLON,NLAT,NSNOLAY, 2)
      CALL WRITER (IU, FI_ICEXY,   NLON,NLAT,1,       2)
      CALL WRITER (IU, TI_ICEXY,   NLON,NLAT,1,       2)
c
c        comsoi
c
      CALL WRITER (IU, TSOIXY,  NLON,NLAT,NSOILAY, 2)
      CALL WRITER (IU, WSOIXY,  NLON,NLAT,NSOILAY, 2)
      CALL WRITER (IU, WISOIXY, NLON,NLAT,NSOILAY, 2)
      CALL WRITER (IU, WPUDXY,  NLON,NLAT,1,       2)
      CALL WRITER (IU, WIPUDXY, NLON,NLAT,1,       2)
      CALL WRITER (IU, TGXY,    NLON,NLAT,1,       2)
c
c        comice
c
      CALL WRITER (IU, HICEXY,  NLON,NLAT,NICELAY, 2)
      CALL WRITER (IU, TICEXY,  NLON,NLAT,NICELAY, 2)
      CALL WRITER (IU, FICEXY,  NLON,NLAT,1, 2)
      CALL WRITER (IU, BRINEXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, TG_ICEXY,NLON,NLAT,1, 2)
c
c        comveg
c
      CALL WRITER (IU, TUXY,    NLON,NLAT,1, 2)
      CALL WRITER (IU, TSXY,    NLON,NLAT,1, 2)
      CALL WRITER (IU, TLXY,    NLON,NLAT,1, 2)
      CALL WRITER (IU, TLSUBXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, WLIQUXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, WLIQSXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, WLIQLXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, WSNOUXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, WSNOSXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, WSNOLXY, NLON,NLAT,1, 2)
      CALL WRITER (IU, T12XY,   NLON,NLAT,1, 2)
      CALL WRITER (IU, T34XY,   NLON,NLAT,1, 2)
      CALL WRITER (IU, Q12XY,   NLON,NLAT,1, 2)
      CALL WRITER (IU, Q34XY,   NLON,NLAT,1, 2)
      CALL WRITER (IU, FOG12XY, NLON,NLAT,1, 2)
      CALL WRITER (IU, FOG34XY, NLON,NLAT,1, 2)
      CALL WRITER (IU, HFOG12XY,NLON,NLAT,1, 2)
      CALL WRITER (IU, HFOG34XY,NLON,NLAT,1, 2)
#ifdef ibis
      IF (PHOTOTYPE.EQ.0) THEN
        CALL RESETR (WORKL, NLON*NLAT, 1.E20)    ! flag for no IBIS info
        DO IDUM=1,13
          CALL WRITER (IU, WORKL,      NLON,NLAT,1, 2)
        ENDDO
      ELSE
        CALL WRITER (IU, VMAXSUNUXY,   NLON,NLAT,2, 2)
        CALL WRITER (IU, VMAXSUNLXY,   NLON,NLAT,2, 2)
        CALL WRITER (IU, VMAXSHAUXY,   NLON,NLAT,2, 2)
        CALL WRITER (IU, VMAXSHALXY,   NLON,NLAT,2, 2)
        CALL WRITER (IU, UPTOTUXY,     NLON,NLAT,1, 2)
        CALL WRITER (IU, UPTOTLXY,     NLON,NLAT,1, 2)
        CALL WRITER (IU, T15RUNXY,     NLON,NLAT,1, 2)
        CALL WRITER (IU, AGDD0XY,      NLON,NLAT,1, 2)
        CALL WRITER (IU, AGDD5XY,      NLON,NLAT,1, 2)
      ENDIF
#else
      CALL RESETR (WORKL, NLON*NLAT, 1.E20)
      DO 400 IDUM=1,13
        CALL WRITER (IU, WORKL, NLON,NLAT,1, 2)
  400 CONTINUE
#endif
c
c        comeve
c
      CALL WRITER (IU, TEVE_D,  NLON,NLAT,1, 2)
      CALL WRITER (IU, PEVE_D,  NLON,NLAT,1, 2)
      CALL WRITER (IU, REVE_D,  NLON,NLAT,1, 2)
c
c        commud
c
      CALL WRITER (IU, GIRRIGXY,  NLON,NLAT,1, 2)

#ifdef isotrac
        call writer (iu, hi_i,   nlon,nlat,1, 2)
        call writer (iu, ri_o,   nlon,nlat,1, 2)
        call writer (iu, ri_h,   nlon,nlat,1, 2)
        call writer (iu, hct_i,  nlon,nlat,1, 2)
        call writer (iu, rct_o,  nlon,nlat,1, 2)
        call writer (iu, rct_h,  nlon,nlat,1, 2)
        call writer (iu, hcb_i,  nlon,nlat,1, 2)
        call writer (iu, rcb_o,  nlon,nlat,1, 2)
        call writer (iu, rcb_h,  nlon,nlat,1, 2)
        call writer (iu, hl_i,   nlon,nlat,1, 2)
        call writer (iu, ro_o,   nlon,nlat,1, 2)
        call writer (iu, ro_h,   nlon,nlat,1, 2)
#endif

c       Last word check for rcm or hydmodel (done in rsds,wsds for gcm) 

#if defined (rcm) || defined (hydmodel)
      write (iu) nwrit_a, nuoff_a
      write (iu) versreszz
#endif

      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE LINSOLVE (ARR, RHS, VEC, MPLATE, NPOI, ND, IUN)
c
c        Solves multiple linear systems of equations, vectorizing
c        over the number of systems. Basic Gaussian elimination is
c        used, with no pivoting (relies on all diagonal elements
c        being and staying significantly non-zero). For efficiency
c        a template array MPLATE is used to detect when an operation
c        is not necessary (element already zero or would add zeros),
c        assuming that every system has the same pattern of zero
c        elements. This template is first copied to MPLATEX since it
c        must be updated during the procedure in case an original-zero
c        pattern location becomes non-zero.
c
c        The first subscript in ARR, RHS, VEC is over the multiple
c        systems, and the others are the usual row, column subscripts.
c
c        ARR    = equation coefficients (supplied, overwritten)
c        RHS    = equation right-hand sides (supplied, overwritten)
c        VEC    = solution (returned)
c        MPLATE = pattern of zero elements of ARR (supplied)
c        NPOI   = number of systems, le NVEC (supplied)
c        ND     = number of equations (supplied)
c        IUN    = logical unit number for error messages (supplied)
c
#include <compar.h>
      PARAMETER (NDX=11)
      DIMENSION ARR(NVEC,ND,ND), RHS(NVEC,ND), VEC(NVEC,ND),
     *          MPLATE(ND,ND)
      DIMENSION F(NVEC), MPLATEX(NDX,NDX)
c
c
      IF (ND.GT.NDX) THEN
         WRITE(IUN,900) ND, NDX
  900    FORMAT(/' *** Fatal Error ***'/
     *          /' Number of LINSOLVE eqns',I4,' exceeds limit',I4)
         CALL ENDRUN (-1)
      ENDIF
c
c        Copy the zero template so it can be changed below
c
      DO 5 I=1,ND
         DO 6 J=1,ND
            MPLATEX(I,J) = MPLATE(I,J)
    6    CONTINUE
    5 CONTINUE
c
c        Zero all array elements below the diagonal, proceeding from
c        the first row to the last. Note that MPLATEX is set non-zero
c        for changed (I,J) locations, in loop 20.
c
      DO 10 ID=1, ND-1
         DO 12 I=ID+1,ND
c
            IF (MPLATEX(I,ID).NE.0) THEN
               DO 14 M=1,NPOI
                  F(M) = ARR(M,I,ID) / ARR(M,ID,ID)
   14          CONTINUE
c
               DO 20 J=ID,ND
                  IF (MPLATEX(ID,J).NE.0) THEN
                     DO 22 M=1,NPOI
                        ARR(M,I,J) = ARR(M,I,J) - F(M)*ARR(M,ID,J)
   22                CONTINUE
                     MPLATEX(I,J) = 1
                  ENDIF
   20          CONTINUE
c
               DO 30 M=1,NPOI
                  RHS(M,I) = RHS(M,I) - F(M)*RHS(M,ID)
   30          CONTINUE
            ENDIF
c
   12    CONTINUE
   10 CONTINUE
c
c        All array elements below the diagonal are zero, so can
c        immediately solve the equations in reverse order.
c
      DO 50 ID=ND,1,-1
c
         CALL ZERO (F, NPOI)
         IF (ID.LT.ND) THEN
            DO 52 J=ID+1,ND
               IF (MPLATEX(ID,J).NE.0) THEN
                  DO 54 M=1,NPOI
                     F(M) = F(M) + ARR(M,ID,J)*VEC(M,J)
   54             CONTINUE
               ENDIF
   52       CONTINUE
         ENDIF
c
         DO 56 M=1,NPOI
            VEC(M,ID) = (RHS(M,ID) - F(M)) / ARR(M,ID,ID)
   56    CONTINUE
c
   50 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE TRIDIA (NS, ND, NE, A, B, C, Y, X, ALPHA, GAMMA)
c
#include <compar.h>
c***
c***     Purpose:
c***     To compute the solution of many tridiagonal linear systems.
c***
c**      Arguments:
c**
c**      NS ..... The number of systems to be solved.
c**
c**      ND ..... First dimension of arrays (ge NS).
c**
c**      NE ..... The number of unknowns in each system.
c**               This must be > 2. Second dimension of arrays.
c**
c**      A ...... The subdiagonals of the matrices are stored
c**               in locations A(j,2) through A(j,NE).
c**
c**      B ...... The main diagonals of the matrices are stored
c**               in locations B(j,1) through B(j,NE).
c**
c**      C ...... The super-diagonals of the matrices are stored in
c**               locations C(j,1) through C(j,NE-1).
c**
c**      Y ...... The right hand side of the equations is stored in
c**               Y(j,1) through Y(j,NE).
c**
c**      X ...... The solutions of the systems are returned in
c**               locations X(j,1) through X(j,NE).
c**
c**      ALPHA .. Work array dimensioned ALPHA(ND,NE)
c**
c**      GAMMA .. Work array dimensioned GAMMA(ND,NE)
c**
c*       History:  Based on a streamlined version of the old NCAR
c*                 ULIB subr TRDI used in the PHOENIX climate
c*                 model of Schneider and Thompson (J.G.R., 1981).
c*                 Revised by Starley Thompson to solve multiple
c*                 systems and vectorize well on the CRAY-1.
c*                 Later revised to include a PARAMETER statement
c*                 to define loop limits and thus enable Cray short
c*                 vector loops.
c*
c*       Algorithm:  LU decomposition followed by solution.
c*                   NOTE: This subr executes satisfactorily
c*                   if the input matrix is diagonally dominant
c*                   and non-singular.  The diagonal elements are
c*                   used to pivot, and no tests are made to determine
c*                   singularity. If a singular or numerically singular
c*                   matrix is used as input a divide by zero or
c*                   floating point overflow will result.
c*
c*       Last Revision Date:      4 February 1988
c*
      DIMENSION       A(ND,NE), B(ND,NE), C(ND,NE), Y(ND,NE),
     *                X(ND,NE), ALPHA(ND,NE), GAMMA(ND,NE)
c
      NM1 = NE-1
c
c        Obtain the LU decompositions.
c
      DO 10 J=1,NS
         ALPHA(J,1) = 1./B(J,1)
         GAMMA(J,1) = C(J,1)*ALPHA(J,1)
   10 CONTINUE
      DO 11 I=2,NM1
         DO 12 J=1,NS
            ALPHA(J,I) = 1./(B(J,I)-A(J,I)*GAMMA(J,I-1))
            GAMMA(J,I) = C(J,I)*ALPHA(J,I)
   12    CONTINUE
   11 CONTINUE
c
c        Solve
c
      DO 20 J=1,NS
         X(J,1) = Y(J,1)*ALPHA(J,1)
   20 CONTINUE
      DO 21 I=2,NM1
         DO 22 J=1,NS
            X(J,I) = (Y(J,I)-A(J,I)*X(J,I-1))*ALPHA(J,I)
   22    CONTINUE
   21 CONTINUE
      DO 23 J=1,NS
         X(J,NE) = (Y(J,NE)-A(J,NE)*X(J,NM1))/
     *             (B(J,NE)-A(J,NE)*GAMMA(J,NM1))
   23 CONTINUE
      DO 24 I=1,NM1
         IB = NE-I
         DO 25 J=1,NS
            X(J,IB) = X(J,IB)-GAMMA(J,IB)*X(J,IB+1)
   25    CONTINUE
   24 CONTINUE
c
      RETURN
      END
c
c----------------------------------------------------------------------
c
#ifndef cray
 
      function cvmgt (x,y,l)
 
c        Duplicates Cray-Vector-Merge-GT function.
c        Only safe when first two arguments are reals (especially
c        if integers are 32 bits and reals are 64 bits).
 
      logical l
      if (l) then
        cvmgt = x
      else
        cvmgt = y
      endif
      return
      end
 
#endif
c
c----------------------------------------------------------------------
c
#ifndef cray
 
      integer function ishell (a)
 
c        Duplicates Cray ISHELL function
 
      character*(*) a

#if defined (aix)
      integer*4 result

      ! - result is the return status of a waitpid(2) call and hence
      ! - if status is non-zero it will be different from the return
      ! - status if the same command were exectued from command line
      ! - as waitpid(2) encodes additional information -- see man pg.
      call system(a,result)
      ishell = result
      return
      end
#endif

#if defined  (sun)
      integer*4 system
#elif ! defined (sgi) 
      integer system
#endif
      ishell =  system (a)
      return
      end
 
#endif
c
c----------------------------------------------------------------------
c
#if ! defined (cray) && ! defined (rcm)
 
c        Duplicates Cray WHEN[EQ,NE] (ints) and
c        WHEN[FLT,FLE,FGT,FGE] (reals) subrs
 
      subroutine wheneq (n, arr, inc, target, index, nval)
      integer arr(*), index(*), target
      nval = 0
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).eq.target) then
          nval = nval + 1
          index(nval) = i
        endif
        ia = ia + inc
  100 continue
      return
      end
 
      subroutine whenne (n, arr, inc, target, index, nval)
      integer arr(*), index(*), target
      nval = 0
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).ne.target) then
          nval = nval + 1
          index(nval) = i
        endif
        ia = ia + inc
  100 continue
      return
      end
 
      subroutine whenflt (n, arr, inc, target, index, nval)
      dimension arr(*), index(*)
      nval = 0
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).lt.target) then
          nval = nval + 1
          index(nval) = i
        endif
        ia = ia + inc
  100 continue
      return
      end
 
      subroutine whenfle (n, arr, inc, target, index, nval)
      dimension arr(*), index(*)
      nval = 0
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).le.target) then
          nval = nval + 1
          index(nval) = i
        endif
        ia = ia + inc
  100 continue
      return
      end
 
      subroutine whenfgt (n, arr, inc, target, index, nval)
      dimension arr(*), index(*)
      nval = 0
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).gt.target) then
          nval = nval + 1
          index(nval) = i
        endif
        ia = ia + inc
  100 continue
      return
      end
 
      subroutine whenfge (n, arr, inc, target, index, nval)
      dimension arr(*), index(*)
      nval = 0
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).ge.target) then
          nval = nval + 1
          index(nval) = i
        endif
        ia = ia + inc
  100 continue
      return
      end
 
#endif
c
c----------------------------------------------------------------------
c
#if ! defined  (cray) && ! defined (rcm)
 
c        Duplicates Cray ISRCH[EQ,LT,LE,GT,GE] functions
 
      function isrcheq (n, arr, inc, target)
      dimension arr(*)
      if (n.le.0) then
        isrcheq = 0
        return
      endif
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).eq.target) then
          isrcheq = i
          return
        endif
        ia = ia + inc
  100 continue
      isrcheq = n+1
      return
      end
 
      function isrchflt (n, arr, inc, target)
      dimension arr(*)
      if (n.le.0) then
        isrchflt = 0
        return
      endif
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).lt.target) then
          isrchflt = i
          return
        endif
        ia = ia + inc
  100 continue
      isrchflt = n+1
      return
      end
 
      function isrchfle (n, arr, inc, target)
      dimension arr(*)
      if (n.le.0) then
        isrchfle = 0
        return
      endif
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).le.target) then
          isrchfle = i
          return
        endif
        ia = ia + inc
  100 continue
      isrchfle = n+1
      return
      end
 
      function isrchfgt (n, arr, inc, target)
      dimension arr(*)
      if (n.le.0) then
        isrchfgt = 0
        return
      endif
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).gt.target) then
          isrchfgt = i
          return
        endif
        ia = ia + inc
  100 continue
      isrchfgt = n+1
      return
      end
 
      function isrchfge (n, arr, inc, target)
      dimension arr(*)
      if (n.le.0) then
        isrchfge = 0
        return
      endif
      ia = 1
      if (inc.lt.0) ia = 1 - (n-1)*inc
      do 100 i=1,n
        if (arr(ia).ge.target) then
          isrchfge = i
          return
        endif
        ia = ia + inc
  100 continue
      isrchfge = n+1
      return
      end
 
#endif
c
c----------------------------------------------------------------------
c
#if ! defined  (cray) && ! defined (rcm)
 
c        Duplicates Cray ISMAX, ISMIN, INTMAX functions
 
      function ismax (n, arr, inc)
      dimension arr(*)
      zmax = -1.e25
      ia = 1 + (n-1)*inc
      if (inc.lt.0) ia = 1
      do 100 i=1,n
        if (arr(ia).ge.zmax) then
          zmax = arr(ia)
          ismax = ia
        endif
        ia = ia - inc
  100 continue
      return
      end
 
      function ismin (n, arr, inc)
      dimension arr(*)
      zmin = 1.e25
      ia = 1 + (n-1)*inc
      if (inc.lt.0) ia = 1
      do 100 i=1,n
        if (arr(ia).le.zmin) then
          zmin = arr(ia)
          ismin = ia
        endif
        ia = ia - inc
  100 continue
      return
      end

      function intmax (n, ix, inc)
      dimension ix(*)
      mx = ix(1)
      intmax = 1
      do i = 1+inc, inc*n, inc
        if (ix(i).gt.mx) then
          mx = ix(i)
          intmax = i
       endif
      enddo
      return
      end
 
#endif
c
c----------------------------------------------------------------------
c
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if ! defined  (cray) && ! defined (rcm)
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
c        Duplicates Cray SGEDI,SGEFA and lower subrs (for matrix
c        inversion, called from SETTAU in spect.f, replacing old call
c        to MINV). SGEDI... code obtained by SCD's "dsl" utility
c        from the LINPACK library. On Crays, these are in /lib/libsci.a
c        and /lib/libu.a.
 
c-----
 
      SUBROUTINE SGEDI (A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(*),JOB
cdp   REAL A(LDA,1),DET(2),WORK(*)
      REAL A(LDA,LDA),DET(2),WORK(*)  ! dp
C
C     SGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        WORK    REAL(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     REAL(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF SGECO HAS SET RCOND .GT. 0.0 OR SGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,SSWAP
C     FORTRAN ABS,MOD
C
C     INTERNAL VARIABLES
C
      REAL T
      REAL TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         TEN = 10.0E0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0E0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0E0/A(K,K)
            T = -A(K,K)
            CALL SSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0E0
               CALL SAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0E0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL SAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL SSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
 
c-----
 
      SUBROUTINE SGEFA (A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
cdp   REAL A(LDA,1)
      REAL A(LDA,LDA)  ! dp
C
C     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,ISAMAX
C
C     INTERNAL VARIABLES
C
      REAL T
      INTEGER ISAMAX,J,K,KP1,L,NM1
C
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ISAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0E0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0E0/A(K,K)
            CALL SSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
 
c-----
 
      SUBROUTINE SAXPY (N,SA,SX,INCX,SY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL SX(*),SY(*),SA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (SA .EQ. 0.0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
      END
 
c-----
 
      SUBROUTINE SSCAL (N,SA,SX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL SA,SX(*)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END
 
c-----
 
      SUBROUTINE SSWAP (N,SX,INCX,SY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL SX(*),SY(*),STEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        STEMP = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        STEMP = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP
        STEMP = SX(I + 1)
        SX(I + 1) = SY(I + 1)
        SY(I + 1) = STEMP
        STEMP = SX(I + 2)
        SX(I + 2) = SY(I + 2)
        SY(I + 2) = STEMP
   50 CONTINUE
      RETURN
      END
 
c-----------------------------------------------------------------------
 
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL SX(*),SMAX
      INTEGER I,INCX,IX,N
C
      ISAMAX = 0
      IF( N .LT. 1 ) RETURN
      ISAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(SX(IX)).LE.SMAX) GO TO 5
         ISAMAX = I
         SMAX = ABS(SX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         IF(ABS(SX(I)).LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END
 
c<<<<<
#endif
c<<<<<
c
c-----------------------------------------------------------------------
c
c     On sgi, nint and int are declared external throughout (in compar
c     and params), in order to use own versions and avoid sgi bug
c     with -O3 and -i8.
c
c     Also the system's flush doesn't do anything. The flush below
c     assumes it is a "scratch" file with default name fort.<nn>.
 
#if defined (sgi)
      integer function nint(z)
      if (z.ge.0) then
        nint = z + 0.5
      else
        nint = z - 0.5
      endif
      return
      end
 
      integer function int(z)
      int = z
      return
      end
#endif
 
#if defined (NOFLUSH)
c     with coupled mom2,4, can't flush any unit number 
c     or else mom's get_unit() goes haywire
      integer function flush(iu)
      flush = 0
      return
      end
#elif defined (OWNFLUSH)
      integer function flush(iu)
      if (iu.eq.0 .or. iu.eq.5 .or. iu.eq.6) then  ! not for std[in,out]
c     if (iu.eq.0 .or. iu.eq.5             ) then  ! ok for fort.6
         flush = 0
         return
      endif
      close (iu, err = 1000, iostat = ios1)
 1000 open  (iu, status='old', access='append', err=2000, iostat=ios2)
 2000 if (ios1.eq.0 .and. ios2.eq.0) then
        flush = 0
      else
        flush = -1
      endif
      return
      end
#endif

c----------------------------------------------------------------------

#if defined ( aix )
      integer function flush(lun)
      ! - call XLF general service and utility procedure:
      integer*4 lun4
      lun4 = lun
      call flush_(lun4)
      flush = 0
      return
      end
#endif
