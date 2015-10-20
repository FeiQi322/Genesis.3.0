c -- $Id: vegdat.F,v 1.3 2001/07/03 21:29:31 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/vegdat.F,v $
c -- $Name:  $

c-------
c vegdat
c-------
c
c----------------------------------------------------------------------
c
      subroutine vegdat (jyear, jmonth, jday, fsave_leaf, triassic,
     *                   nsrest_ibis, lonreslsx, latreslsx, typeresin,
     *                   showmaps, co2ppm)
 
c        Calls appropriate vegdat_..., depending on vegtype (comgrd),
c        and vegdat_all for settings common to all.
c        See vegdat_... for argument descriptions.

 
#include <compar.h>
#include <comgrd.h>
#include <comvegxy.h>
      logical fsave_leaf, triassic, showmaps
      character*(*) typeresin
      dimension rootdepth (nlon,nlat,2)
 
 
c        Set D+S, Wisconsin, Eve/Leaf,  BIOME or SDVM vegetation,
c        depending on VEGTYPE
 
c---------------------------
      if (vegtype.eq.0) then
c---------------------------
 
        call vegdat_dorm (jyear, jmonth, jday, triassic, showmaps,
     *                    rootdepth, 0)
 
c--------------------------------
      else if (vegtype.eq.1) then
c--------------------------------
 
        call vegdat_wisc (showmaps, jyear, jmonth, nsrest_ibis,
     *                    lonreslsx, latreslsx, typeresin)
 
c--------------------------------------------------
      else if (vegtype.eq.2 .or. vegtype.eq.3) then
c--------------------------------------------------
 
        call vegdat_eve (jyear, jmonth, jday, fsave_leaf, showmaps,
     *                   rootdepth)
 
c          If requested, overwrite eve veg in specified rectangles
c          with a particular D+S veg type
 
        if (dsvegs.ne.0 .or. dsvegn.ne.0) then
 
          call vegdat_dorm (jyear, jmonth, jday, triassic, showmaps,
     *                      rootdepth, 1)
        else
 
c         zero D+S/Wisc vegetation index (to protect diagnostics):
          call izero (ivegxy, nlon*nlat)
 
        endif
 
c--------------------------------------------------
      else if (vegtype.eq.4 .or. vegtype.eq.5) then
c--------------------------------------------------
 
        call vegdat_biome (jyear, jmonth, jday, showmaps, rootdepth,
     *                     co2ppm)

c--------------------------------------------------
      else if (vegtype.eq.6 .or. vegtype.eq.7) then
c--------------------------------------------------

        call vegdat_sdvm (jyear, jmonth, jday, showmaps, rootdepth)
 
c----------
      endif
c----------
 
c        Set veg parameters that are common to both D+S and Eve/Leaf
 
      call vegdat_all (rootdepth)
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vegdat_dorm (jyear, jmonth, jday, triassic, showmaps,
     *                        rootdepth, ifoverwrite)
c
c        Sets current state of vegetation on lsx 2-D grid
c        using Dorman and Sellers (1989) JAM,28,833-855 (henceforth DS).
c
c        jyear  = year number (eg, 1990) (supplied)
c        jmonth = month number (1-12) (supplied)
c        jday   = day number (1-31) nb:agcm has no leap years (supplied)
c        triassic = .true. for Triassic lai,sai,fu,fl settings(supplied)
c        showmaps = .true. to print 2-D maps of prescribed fields (supp)
c        rootdepth= maximum root depth for each story (returned)
c        ifoverwrite = 0 for normal call,
c                    = 1 for after a call to vegdat_eve, to overwrite
c                      with specified D+S type in certain rectangles
c                      (requested by Namelist DS[LON,LAT,VEG][S,N})
c                    = 2 for call from vegdat_biome, with ivegxy already
c                      translated to closest D+S types for each biome #,
c                      to set physical attributes
c
#include <compar.h>
#include <comfud.h>
#include <comgrd.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commapl.h>
#include <commugxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>

      logical triassic, showmaps
      dimension rootdepth (nlon,nlat,2)
 
      dimension ndaypm(12), saiuzav(nved)
      logical lfirs, lfirsa, lfirsb, lfirsc
      save ndaypm, saiuzav, lfirs, lfirsa, lfirsb, lfirsc
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data lfirs, lfirsa, lfirsb, lfirsc /.true.,.true.,.true.,.true./

      character*1 cimage(nlon,nlat)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c comved.h contains D-S specific arrays, initialized by data
c statements below. Declared as local (dimension...) and saved,
c not in common, to allow different data settings in subr vegdat_sdvm.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <comved.h>    

c     total lai+sai, one-sided, lower story
      data tailz/
     *   0,   0,   0,   0,   0, 1.8, 0.7, 2.0,   0, 0.5,   0, 1.3,   0,
     *   0,   0,   0,   0,   0, 1.6, 0.8, 2.0,   0, 0.5,   0, 1.5,   0,
     *   0,   0,   0,   0,   0, 1.5, 0.9, 1.7,   0, 0.5,   0, 1.7,   0,
     *   0,   0,   0,   0,   0, 1.4, 1.0, 1.2,   0, 0.5,   0, 1.9,   0,
     *   0,   0,   0,   0,   0, 1.3, 1.5, 0.8,   0, 0.5,   0, 3.5,   0,
     *   0,   0,   0,   0,   0, 1.5, 3.4, 0.8,   0, 0.6,   0, 4.7,   0,
     *   0,   0,   0,   0,   0, 1.7, 4.3, 0.8,   0, 0.6,   0, 0.3,   0,
     *   0,   0,   0,   0,   0, 4.8, 3.8, 2.4,   0, 0.7,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 4.2, 1.8, 4.3,   0, 0.5,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 3.5, 1.0, 3.7,   0, 0.5,   0, 0.6,   0,
     *   0,   0,   0,   0,   0, 2.8, 0.9, 2.7,   0, 0.5,   0, 0.8,   0,
     *   0,   0,   0,   0,   0, 2.0, 0.8, 2.0,   0, 0.5,   0, 1.1,   0/
c
c     green lai, one-sided, lower story
      data gailz/
     *   0,   0,   0,   0,   0, 0.8, 0.4, 1.3,   0, 0.4,   0, 1.1,   0,
     *   0,   0,   0,   0,   0, 0.7, 0.5, 1.0,   0, 0.4,   0, 1.3,   0,
     *   0,   0,   0,   0,   0, 0.6, 0.6, 0.5,   0, 0.4,   0, 1.5,   0,
     *   0,   0,   0,   0,   0, 0.5, 0.7, 0.0,   0, 0.4,   0, 1.7,   0,
     *   0,   0,   0,   0,   0, 0.5, 1.2, 0.0,   0, 0.4,   0, 3.1,   0,
     *   0,   0,   0,   0,   0, 0.7, 3.0, 0.0,   0, 0.5,   0, 2.8,   0,
     *   0,   0,   0,   0,   0, 0.9, 3.5, 0.0,   0, 0.5,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 3.5, 1.5, 1.6,   0, 0.4,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 2.7, 0.7, 3.0,   0, 0.4,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 2.0, 0.6, 2.0,   0, 0.4,   0, 0.5,   0,
     *   0,   0,   0,   0,   0, 1.1, 0.5, 1.3,   0, 0.4,   0, 0.7,   0,
     *   0,   0,   0,   0,   0, 1.0, 0.4, 1.3,   0, 0.4,   0, 1.0,   0/
c
c     fractional cover, lower story
      data flz/
     *   0,   0,   0,   0,   0, .83, .45, .87,   0, .40,   0, .65,   0,
     *   0,   0,   0,   0,   0, .80, .50, .87,   0, .40,   0, .70,   0,
     *   0,   0,   0,   0,   0, .78, .54, .83,   0, .40,   0, .73,   0,
     *   0,   0,   0,   0,   0, .76, .57, .67,   0, .40,   0, .77,   0,
     *   0,   0,   0,   0,   0, .74, .71, .51,   0, .40,   0, .88,   0,
     *   0,   0,   0,   0,   0, .78, .87, .52,   0, .46,   0, .90,   0,
     *   0,   0,   0,   0,   0, .83, .89, .51,   0, .48,   0, .20,   0,
     *   0,   0,   0,   0,   0, .99, .88, .90,   0, .49,   0, .00,   0,
     *   0,   0,   0,   0,   0, .98, .75, .98,   0, .40,   0, .00,   0,
     *   0,   0,   0,   0,   0, .97, .58, .97,   0, .40,   0, .41,   0,
     *   0,   0,   0,   0,   0, .94, .54, .93,   0, .40,   0, .50,   0,
     *   0,   0,   0,   0,   0, .86, .50, .87,   0, .40,   0, .58,   0/
c
c     total lai+sai, one-sided, upper story
      data taiuz/
     * 5.0, 0.4, 3.5, 6.6, 0.3, 1.2,   0, 0.4, 1.0, 0.1,   0, 0.0,   0,
     * 5.0, 0.4, 3.6, 6.8, 0.3, 1.0,   0, 0.4, 1.0, 0.1,   0, 0.0,   0,
     * 5.0, 0.7, 3.9, 7.3, 0.3, 0.9,   0, 0.4, 0.8, 0.1,   0, 0.0,   0,
     * 5.0, 1.6, 4.9, 7.5, 1.0, 0.8,   0, 0.2, 0.3, 0.1,   0, 0.2,   0,
     * 5.0, 3.5, 5.5, 7.7, 1.6, 0.8,   0, 0.1, 0.6, 0.1,   0, 0.3,   0,
     * 5.0, 5.1, 6.2, 8.0, 2.4, 1.0,   0, 0.1, 0.0, 0.3,   0, 0.5,   0,
     * 5.0, 5.4, 6.4, 7.8, 4.3, 2.0,   0, 0.0, 0.1, 1.5,   0, 0.5,   0,
     * 5.0, 4.8, 6.5, 7.7, 2.9, 3.7,   0, 0.6, 0.3, 1.7,   0, 0.5,   0,
     * 5.0, 3.8, 5.7, 7.6, 2.0, 3.2,   0, 1.0, 0.5, 1.4,   0, 0.4,   0,
     * 5.0, 1.7, 4.6, 7.2, 1.3, 2.7,   0, 0.8, 0.6, 0.1,   0, 0.2,   0,
     * 5.0, 0.6, 3.2, 6.7, 0.8, 1.9,   0, 0.5, 0.7, 0.1,   0, 0.1,   0,
     * 5.0, 0.4, 3.2, 6.5, 0.5, 1.2,   0, 0.4, 0.9, 0.1,   0, 0.0,   0/
c
c     green lai, one-sided, upper story
      data gaiuz/
     * 4.5, 0.0, 3.0, 6.0, 0.0, 0.8,   0, 0.3, 0.9, 0.0,   0, 0.0,   0,
     * 4.5, 0.0, 3.1, 6.2, 0.0, 0.7,   0, 0.2, 0.8, 0.0,   0, 0.0,   0,
     * 4.5, 0.3, 3.4, 6.7, 0.0, 0.4,   0, 0.1, 0.2, 0.0,   0, 0.0,   0,
     * 4.5, 1.2, 4.3, 7.0, 0.6, 0.5,   0, 0.0, 0.2, 0.0,   0, 0.2,   0,
     * 4.5, 3.0, 5.0, 7.2, 1.2, 0.5,   0, 0.0, 0.0, 0.0,   0, 0.3,   0,
     * 4.5, 4.7, 5.7, 7.3, 2.0, 0.7,   0, 0.0, 0.0, 0.2,   0, 0.5,   0,
     * 4.5, 4.5, 6.0, 7.0, 2.6, 1.7,   0, 0.0, 0.0, 1.4,   0, 0.4,   0,
     * 4.5, 3.4, 5.3, 6.9, 1.7, 3.0,   0, 0.5, 0.2, 1.2,   0, 0.3,   0,
     * 4.5, 1.2, 4.0, 6.7, 1.0, 2.5,   0, 0.7, 0.4, 0.0,   0, 0.1,   0,
     * 4.5, 0.3, 2.8, 6.2, 0.5, 1.6,   0, 0.3, 0.5, 0.0,   0, 0.0,   0,
     * 4.5, 0.0, 2.8, 5.9, 0.2, 1.0,   0, 0.3, 0.6, 0.0,   0, 0.0,   0,
     * 4.5, 0.0, 2.8, 5.9, 0.0, 1.0,   0, 0.3, 0.8, 0.0,   0, 0.0,   0/
c
c     fractional cover, upper story
      data fuz/
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0,
     * .98, .75, .75, .75, .50, .30,   0, .10, .10, .30,   0, .08,   0/
c
c *** Transmittances ***
c
c     VIS leaf transmittance, lower story (1st row = live, 2nd = dead)
      data tauvlz/
     * .07, .07, .07, .07, .07, .07, .07, .07, .07, .07, .07, .07, .07,
     * .22, .22, .22, .22, .22, .22, .22, .22, .22, .22, .22, .22, .22/
c
c     VIS leaf transmittance, upper story (1st row = live, 2nd = dead)
      data tauvuz/
     * .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05,
     *.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001/
c
c     NIR leaf transmittance, lower story (1st row = live, 2nd = dead)
      data taunlz/
     * .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25,
     * .38, .38, .38, .38, .38, .38, .38, .38, .38, .38, .38, .38, .38/
c
c     NIR leaf transmittance, upper story (1st row = live, 2nd = dead)
      data taunuz/
     * .25, .25, .15, .10, .10, .25, .25, .25, .25, .25, .25, .25, .25,
     *.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001/
c
c *** Reflectances ***
c
c     VIS leaf reflectance, lower story (1st row = live, 2nd = dead)
      data rhovlz/
     * .11, .11, .11, .11, .11, .11, .11, .11, .11, .11, .11, .11, .11,
     * .36, .36, .36, .36, .36, .36, .36, .36, .36, .36, .36, .36, .36/
c
c     VIS leaf reflectance, upper story (1st row = live, 2nd = dead)
      data rhovuz/
     * .10, .10, .07, .07, .07, .10, .10, .10, .10, .10, .10, .10, .10,
     * .16, .16, .16, .16, .16, .16, .16, .16, .16, .16, .16, .16, .16/
c
c     NIR leaf reflectance, lower story (1st row = live, 2nd = dead)
      data rhonlz/
     * .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58,
     * .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58/
c
c     NIR leaf reflectance, upper story (1st row = live, 2nd = dead)
      data rhonuz/
     * .45, .45, .40, .35, .35, .45, .45, .45, .45, .45, .45, .45, .45,
     * .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39/
c
c *** Vegetation angle orientation (-1=vert, 0=random, 1=horiz) ***
c
c     chi factor, lower story
      data chiflz/
     * .00, .00, .00, .00, .00,-.30,-.30,-.30, .00, .20, .00,-.02, .00/
c
c     chi factor, upper story
      data chifuz/
     * .10, .25, .13, .01, .01, .01, .00, .10, .01, .20, .00, .25, .00/
c
c *** Vegetation heights ***
c
c     vegetation base, lower story (m)
      data z1lz/
     * .02, .02, .02, .02, .02, .02, .02, .02, .02, .02, .02, .02, .02/
c
c     vegetation top, lower story (m)
      data z2lz/
     * .02, .02, .02, .02, .02, .60, .60, .06, .02, .06, .02, .60, .02/
c
c     vegetation base, upper story (m)
      data z1uz/
     *  1., 11., 10.,  8.,  7., 10., .61, .10, .10, .10, .03, 11., .03/
c
c     vegetation top, upper story (m)
      data z2uz/
     * 35., 20., 20., 17., 14., 18., .61, .50, .50, .60, .03, 20., .03/
c
c
c *** PAR coefficients ***
c
c     "a", lower story (J m-3)
      data aparlz/
     *  2582,  2582,  2582,  2582,  2582,  2582,
     *  2582,  2582,  2582,  2582,  2582,  7459,  2582/
c
c     "b", lower story (W m-2)
      data bparlz/
     *  1.09,  1.09,  1.09,  1.09,  1.09,  1.09,
     *  1.09,  1.09,  1.09,  1.09,  1.09,  5.70,  1.09/
c
c     "c", lower story (s m-1)
      data cparlz/
     *   110,   110,   110,   110,   110,   110,
     *   110,   110,   110,   110,   110,    25,   110/
c
c     "a", upper story (J m-3)
      data aparuz/
     *  2336,  9802,  6336,  2870,  2870, 56651,
     * 10000, 93989, 93989,  9802, 10000,  9802, 10000/
c
c     "b", upper story (W m-2)
      data bparuz/
     *  0.01, 10.60,  7.10,  3.70,  3.70, 10.80,
     *    10,  0.01,  0.01, 10.55,    10, 10.60,    10/
c
c     "c", upper story (s m-1)
      data cparuz/
     *   154,   180,   207,   233,   233,   165,
     *   100,   855,   855,   180,   100,   180,   100/
c
c *** vpd coefficient ***
c
c     "h5", lower story (mb-1)
      data h5lz/
     * .0238, .0238, .0238, .0238, .0238, .0238,
     * .0238, .0238, .0238, .0238, .0238, .0001, .0238/
c
c     "h5", upper story (mb-1)
      data h5uz/
     * .0273, .0357, .0340, .0310, .0310, .0357,
     * .0100, .0275, .0275, .0275, .0100, .0357,  .0100/
c
c *** winter wheat (from DS Table 5 and table5 program) ***
c
c     total lai+sai, one-sided, for winter wheat, vs (lat_zone, month)
      data taiwheat/
     *     1.21,    1.46,    4.33,    0.79,    0.00,    0.00,    3.27,
     *     1.41,    1.76,    3.84,    0.00,    0.00,    0.00,    0.00,
     *     1.60,    2.81,    1.00,    0.00,    0.00,    0.00,    0.25,
     *     1.80,    4.90,    0.00,    0.00,    0.46,    0.60,    0.63,
     *     1.99,    2.14,    0.00,    0.00,    1.06,    0.93,    0.82,
     *     5.97,    0.00,    0.00,    0.00,    1.64,    1.26,    1.02,
     *     3.24,    0.00,    0.00,    0.00,    4.23,    1.59,    1.21,
     *     0.00,    0.00,    0.00,    0.00,    3.82,    1.93,    1.41,
     *     0.24,    0.00,    0.00,    0.00,    1.08,    4.19,    1.60,
     *     0.63,    0.59,    0.47,    0.22,    0.00,    4.85,    1.80,
     *     0.82,    0.88,    1.06,    5.83,    0.00,    2.16,    1.99,
     *     1.02,    1.17,    1.65,    3.44,    0.00,    0.00,    5.92/
c
c     green lai, one-sided, for winter wheat, vs (lat_zone, month)
      data gaiwheat/
     *     1.09,    1.32,    3.54,    0.08,    0.00,    0.00,    0.98,
     *     1.27,    1.58,    1.78,    0.00,    0.00,    0.00,    0.00,
     *     1.44,    2.40,    0.10,    0.00,    0.00,    0.00,    0.23,
     *     1.62,    3.26,    0.00,    0.00,    0.42,    0.54,    0.57,
     *     1.79,    0.38,    0.00,    0.00,    0.95,    0.84,    0.74,
     *     4.76,    0.00,    0.00,    0.00,    1.48,    1.14,    0.92,
     *     0.93,    0.00,    0.00,    0.00,    3.47,    1.44,    1.09,
     *     0.00,    0.00,    0.00,    0.00,    1.75,    1.73,    1.27,
     *     0.21,    0.00,    0.00,    0.00,    0.12,    3.44,    1.44,
     *     0.57,    0.53,    0.43,    0.20,    0.00,    3.19,    1.62,
     *     0.74,    0.79,    0.96,    4.56,    0.00,    0.39,    1.79,
     *     0.92,    1.05,    1.48,    1.21,    0.00,    0.00,    4.69/
c
c     fractional cover for winter wheat, vs (lat_zone, month)
      data frawheat/
     *     0.59,    0.66,    0.86,    0.24,    0.00,    0.00,    0.55,
     *     0.64,    0.73,    0.67,    0.00,    0.00,    0.00,    0.00,
     *     0.69,    0.82,    0.30,    0.00,    0.00,    0.00,    0.20,
     *     0.75,    0.79,    0.00,    0.00,    0.37,    0.43,    0.43,
     *     0.80,    0.41,    0.00,    0.00,    0.55,    0.52,    0.49,
     *     0.90,    0.00,    0.00,    0.00,    0.70,    0.60,    0.54,
     *     0.55,    0.00,    0.00,    0.00,    0.86,    0.69,    0.59,
     *     0.00,    0.00,    0.00,    0.00,    0.66,    0.78,    0.64,
     *     0.19,    0.00,    0.00,    0.00,    0.31,    0.85,    0.69,
     *     0.43,    0.42,    0.38,    0.18,    0.00,    0.78,    0.75,
     *     0.49,    0.50,    0.55,    0.88,    0.00,    0.42,    0.80,
     *     0.54,    0.58,    0.71,    0.59,    0.00,    0.00,    0.89/
c
c     "chi" orientation factor for winter wheat, vs (lat_zone, month)
      data chiwheat/
     *    -0.05,   -0.06,   -0.13,   -0.24,    0.00,    0.00,   -0.28,
     *    -0.06,   -0.08,   -0.24,    0.00,    0.00,    0.00,    0.00,
     *    -0.07,   -0.10,   -0.30,    0.00,    0.00,    0.00,   -0.01,
     *    -0.08,   -0.20,    0.00,    0.00,   -0.02,   -0.02,   -0.03,
     *    -0.09,   -0.30,    0.00,    0.00,   -0.05,   -0.04,   -0.04,
     *    -0.16,    0.00,    0.00,    0.00,   -0.07,   -0.06,   -0.04,
     *    -0.28,    0.00,    0.00,    0.00,   -0.13,   -0.07,   -0.05,
     *     0.00,    0.00,    0.00,    0.00,   -0.24,   -0.09,   -0.06,
     *    -0.01,    0.00,    0.00,    0.00,   -0.30,   -0.13,   -0.07,
     *    -0.03,   -0.02,   -0.02,   -0.01,    0.00,   -0.20,   -0.08,
     *    -0.04,   -0.04,   -0.05,   -0.17,    0.00,   -0.30,   -0.09,
     *    -0.04,   -0.05,   -0.07,   -0.27,    0.00,    0.00,   -0.16/
c++++++++++++++++++++++++++++++++++++++++++++++ 
c end of data statements for comved.h variables
c++++++++++++++++++++++++++++++++++++++++++++++ 
 
c======
#ifdef rcm
c     lookup table from BATS landuse categories to closest D+S biometype
      dimension lookup_bats(0:19)
      data lookup_bats/13,
     *                 12,  7,  4,  5,  2,  1,  7, 11, 10, 12,  9, 13, 
     *                  6, 13, 13,  8,  8,  4,  6/

c     lookup table from MM4 landuse categories to closest D+S biome type
      dimension lookup_mm4(0:14)
      data lookup_mm4 /13,
     *                  9, 12,  7,  2,  4,  3, 13,  7, 11, 10, 13,  1,
     *                  6, 13/

c     lookup table from MM5 landuse categories to closest D+S biome type
      dimension lookup_mm5(0:25)
      data lookup_mm5 /13,
     *                 11, 12, 12, 12,  7,  6,  7,  8,  8,  8,  2,  5,
     *                  1,  4,  3, 13,  7,  6, 11, 10, 10, 10, 10, 13,
     *                 13 /
#endif
c======

 
c          Compute annual mean upper-story SAI, used below to avoid
c          non-conservation of energy due to addition/removal of
c          wood mass at different seasonal temps.
 
      if (lfirs) then
        do iv=1,nved
          saiuzav(iv) = 0.
          do it=1,ntimds
            saiuzav(iv) = saiuzav(iv) + (taiuz(iv,it)-gaiuz(iv,it))
          enddo
          saiuzav(iv) = saiuzav(iv) / ntimds
        enddo

        do jj=1,nlat
          do ji=1,nlon
            ifimage(ji,jj) = 0   ! in comgrd, used in vegdat_biome
          enddo
        enddo

        lfirs = .false.
      endif
 
c======
 
c        If first call of this run and normal call (ifoverwrite=0),
c        set constants, read veg-type data file derived from DS Fig.1b,
c        put into ivegxy, and set ivegxy values to nved for non-land
c        points and for erroneous values.
c        Note that nominal values for all lsx veg params are set for
c        non-land pts (ivegxy=nved, lai = sai = fractional cover = 0),
c        since radiative and aero calcs are done whether vegetation
c        exists or not, to allow Cray vectorization. Same applies for
c        land points with only one story, and for bare ground.
c
c        If rcm, instead use vegrcm array (in commugxy, set in 
c        inisurf passed from the rcm). These are MM4 or MM5 landuse
c        categories (depending on ISURFTYPE=0/1/2,in commugxy), 
c        converted to closest D+S biome using lookup_... array 
c        (set above).
c
c        If hydrol, instead use veghyd array (in commugxy, set in 
c        inisurf passed from hydrol code). 

      if (lfirsa .and. ifoverwrite.eq.0) then

c-----------------------------
        if (dsbiome.eq.0) then
c-----------------------------

#ifdef rcm
          do jj=1,nlat
            do ii=1,nlon
              if (isurftype.eq.0) then
                workl(ii,jj) = lookup_bats(nint(vegrcm(ii,jj)))
              else if (isurftype.eq.1) then
                workl(ii,jj) = lookup_mm4(nint(vegrcm(ii,jj)))
              else
                workl(ii,jj) = lookup_mm5(nint(vegrcm(ii,jj)))
              endif
            enddo
          enddo
#elif hydmodel
          do jj=1,nlat
            do ii=1,nlon
              workl(ii,jj) = iveghyd(ii,jj)
            enddo
          enddo

#else
          call findiu (iuvegt, loc)
#ifdef gcm
          write(ioterm2,40) locname(loc)(1:lenchr(locname(loc)))
   40     format(/' Reading Dorman-Sellers vegetation data file'/4x,a)
#endif
          call openf (iuvegt, locname(loc), ' ', 0, 0)
          call readdata (iuvegt, workl, nlon, nlat, 'VEGT', 1, 1)
          call closef(iuvegt)
#endif

          do 42 jj=1,nlat
            do 44 ji=1,nlon
              ivegxy(ji,jj) = nint(workl(ji,jj))
              if (lmaskxy(ji,jj).ne.1) then
                ivegxy(ji,jj) = nved
              else
                if (ivegxy(ji,jj).lt.1 .or. ivegxy(ji,jj).gt.nved)
     *            ivegxy(ji,jj) = nved
              endif
   44       continue
   42     continue

c-----------
        else
c-----------

c         Namelist DSBIOME entered for global uniform single type:
          write(ioterm2,50) dsbiome
   50     format(/' *** Using single Dorman-Sellers vegetation biome',
     *            i4,' for all land points')
          do 52 jj=1,nlat
            do 54 ji=1,nlon
              if (lmaskxy(ji,jj).eq.1) then
                ivegxy(ji,jj) = dsbiome
              else
                ivegxy(ji,jj) = nved
              endif
   54       continue
   52     continue

c------------
        endif
c------------

c          Print out veg-type map if requested
 
        if (showmaps)
     *    call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, nved,
     *                 'Vegetation Type')
 
        lfirsa = .false.
      endif
 
c======
 
c        If first call of this run and "overwrite" call following
c        a call to vegdat_eve (ifoverwrite=1),  then:
c        If (dsvegs or dsvegn = 100): read ivegxy from D-S file.
c        Else: set ivegxy to null (0) everywhere, except within the
c              specified rectangles to dsveg[s,n] if land.
c        Null values (ivegxy=0) will be skipped in all loops below.
 
      if (lfirsb .and. ifoverwrite.eq.1) then
 
        if (dsvegs.eq.100 .or. dsvegn.eq.100) then
 
c            Read from data file
 
          call findiu (iuvegt, loc)
#ifdef gcm
          write(ioterm2,55) locname(loc)(1:lenchr(locname(loc)))
   55     format(/' Reading Dorman-Sellers "Overwrite" vegetation',
     *            ' data file'/4x,a)
#endif
          call openf (iuvegt, locname(loc), ' ', 0, 0)
          call readdata (iuvegt, workl, nlon, nlat, 'VEGT', 1, 1)
          call closef(iuvegt)
          do jj=1,nlat
            do ji=1,nlon
              ivegxy(ji,jj) = nint(workl(ji,jj))
              if (ivegxy(ji,jj).lt.0) then
                ivegxy(ji,jj) = 0
              else if (lmaskxy(ji,jj).ne.1. or. ivegxy(ji,jj).gt.nved)
     *          then
                ivegxy(ji,jj) = nved
              endif
            enddo
          enddo
 
        else
 
c            Set specified rectangles
 
          do jj=1,nlat
            zlat = alat(jj)*180./pi
            do ji=1,nlon
              zlon = alon(ji)*180./pi
              if (zlon.gt.180.) zlon = zlon - 360.
              if      ( zlat.ge.dslats(1) .and. zlat.le.dslats(2) .and.
     *                  zlon.ge.dslons(1) .and. zlon.le.dslons(2) ) then
c               southern hemisphere rectangle:
                if (lmaskxy(ji,jj).eq.1) then
                  ivegxy(ji,jj) = dsvegs
                else
                  ivegxy(ji,jj) = nved
                endif
              else if ( zlat.ge.dslatn(1) .and. zlat.le.dslatn(2) .and.
     *                  zlon.ge.dslonn(1) .and. zlon.le.dslonn(2) ) then
c               northern hemisphere rectangle:
                if (lmaskxy(ji,jj).eq.1) then
                  ivegxy(ji,jj) = dsvegn
                else
                  ivegxy(ji,jj) = nved
                endif
              else
c               null (leave Eve veg alone):
                ivegxy(ji,jj) = 0
              endif
            enddo
          enddo
 
        endif
 
c          Print out veg-type map if requested
 
        if (showmaps)
     *    call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, nved,
     *                 'D+S Overwrite Vegetation Type')
 
        lfirsb = .false.
      endif
 
c======

c................................
      if (lfirsc) then
        if(imagefile.ne.' ') then
c................................

c          If requested, read in IMAGE-model file (lat,lon,1/2) 
c          to superimpose agriculture (image type 1 or 2).

c       initialize cimage map 
        do jj=1,nlat
          do ji=1,nlon
            if (lmaskxy(ji,jj).eq.0) then
              cimage(ji,jj) = ' '
            else
              cimage(ji,jj) = '.'
            endif
          enddo
        enddo

        write (ioterm2,58) imagefile(1:lenchr(imagefile))
   58   format(/' *** Reading IMAGE file ',a,' for agriculture')
        iu = 166
        open (iu,file=imagefile, status='old')

        dlat = 180./nlat  
        dlon = 360./nlon  
        imag1 = 0
        imag2 = 0
        do k = 1,100000
          read (iu,*,end=580) zlat, zlon, imv
          jj = nint((zlat+90.)/dlat + 0.5)
          ji = nint( mod(zlon+360.,360.)/dlon + 0.5)
          if (jj.lt.1 .or. jj.gt.nlat .or. ji.lt.1 .or. ji.gt.nlon
     *        .or. imv.lt.1 .or. imv.gt.2) then
            write (ioterm,*) ' file error: k,zlat,zlon,ji,jj,imv=',
     *                       k,zlat,zlon,ji,jj,imv
            call endrun (-1)
          endif
          if (imv.eq.1 .and. lmaskxy(ji,jj).eq.1) then
            ivegxy(ji,jj) = 12     ! IMAGE # 1 to D-S # 12
            ifimage(ji,jj) = 1
            imag1 = imag1 + 1
            cimage(ji,jj) = '1'
          else if (imv.eq.2 .and. lmaskxy(ji,jj).eq.1) then
            ivegxy(ji,jj) = 7      ! IMAGE # 2 to D-S # 7
            ifimage(ji,jj) = 1
            imag2 = imag2 + 1
            cimage(ji,jj) = '2'
          endif
        enddo

  580   close (iu)
        write (ioterm2,582) k-1, imag1, imag2
  582   format(5x,i6,' records read.  imag1=',i6,'  imag2=',i6)

c       write scratch (image) map 
        do j=nlat,1,-1
          write (ioterm2,'(180a1)') (cimage(i,j),i=nlon/2+1,nlon),
     *                              (cimage(i,j),i=1,nlon/2)
        enddo

c          Print out veg-type map if requested
 
        if (showmaps)
     *    call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, nved,
     *                 'D+S+Image Vegetation Type')

c...........
        else
c...........

        if (showmaps)
     *    call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, nved,
     *                 'D+S+NoImage Vegetation Type')

c............
        endif
        lfirsc = .false.
      endif
c..........

c======
 
c        Compute ndaynor, number of days in current year since Jan 1st,
c        ndaysou = ndaynor shifted by 6 months for SH, and ndayequ for
c        the equator, used below for time indices into Figs 2(a)-(o).
c
      nday = 0
      do 60 jm = 1,jmonth
        nday = nday + ndaypm(jm)
   60 continue
      nday = nday - ndaypm(jmonth) + jday
      ndaynor = nday
      ndaysou = mod (nday-1+365/2, 365) + 1
      ndayequ = nday
c
c---------------------
      do 100 jj=1,nlat
c---------------------
c
c          Set LAI, SAI, fractional cover and orientation weights,
c          changing DS one-sided values to LSX two-sided LAI,SAI.
c          These are done separately from other parameters since they
c          are different for winter wheat (ivegxy=12).
c          (Note many lsx veg arrays have last subscript = 1 for lower
c          story, = 2 for upper.)
c
c          First, compute seasonal indices it1,it2,wt1,wt2 for DS
c          Figs 2a-o, linearly weighting between nday[nor,sou] and
c          ndayequ for latitudes equatorward of 25 deg [N,S].
c
        alatd = alat(jj)*180./pi
c
        w = min (1., abs(alatd/25.))
        if (alatd.ge.0.) then
          d = w*ndaynor + (1.-w)*ndayequ
        else
          d = w*ndaysou + (1.-w)*ndayequ
        endif
c
        t = 12. * (d-0.5)/365.
        it1 = t + 0.5
        it2 = it1 + 1
        wt1 = (it1+.5) - t
        wt2 = 1.-wt1
        if (it1.lt. 1) it1 = 12
        if (it2.gt.12) it2 = 1
c
        do 102 ji=1,nlon
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (ivegxy(ji,jj).ne.0) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          iv = ivegxy(ji,jj)
c
          flxy(ji,jj) = wt1*flz(iv,it1) + wt2*flz(iv,it2)
          fuxy(ji,jj) = wt1*fuz(iv,it1) + wt2*fuz(iv,it2)
c
          laixy(ji,jj,1) = 2.*(wt1*gailz(iv,it1) + wt2*gailz(iv,it2))
          saixy(ji,jj,1) = 2.*(wt1*tailz(iv,it1) + wt2*tailz(iv,it2))
     *                     - laixy(ji,jj,1)
          laixy(ji,jj,2) = 2.*(wt1*gaiuz(iv,it1) + wt2*gaiuz(iv,it2))
          saixy(ji,jj,2) = 2.*saiuzav(iv)
 
c         laixy(ji,jj,1) = laixy(ji,jj,1)*0.5 ! if need single sided
c         laixy(ji,jj,2) = laixy(ji,jj,2)*0.5 ! if need single sided
 
#ifndef gcm
c8880     Schwartz experiments
c8880     if (nint(xdate).ge.1000101.) then
c           flxy(ji,jj) = 0.
c           fuxy(ji,jj) = ((nint(xdate)/10000)-100)*.01
c         endif
#endif
 
          orievxy(ji,jj,1) = max (-chiflz(iv), 0.)
          orievxy(ji,jj,2) = max (-chifuz(iv), 0.)
          oriehxy(ji,jj,1) = max ( chiflz(iv), 0.)
          oriehxy(ji,jj,2) = max ( chifuz(iv), 0.)
c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>
  102   continue
c
c          Do special settings for Triassic if requested
c
        if (triassic) then
          do 103 ji=1,nlon
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if (ivegxy(ji,jj).ne.0) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            iv = ivegxy(ji,jj)
            if (iv.eq.4) then
              flxy(ji,jj) = 0.
              fuxy(ji,jj) = 0.75
              laixy(ji,jj,1) = 0.
              saixy(ji,jj,1) = 0.
              laixy(ji,jj,2) = 2.*6.58
              saixy(ji,jj,2) = 2.*7.28 - laixy(ji,jj,2)
            else if (iv.eq.8) then
              flxy(ji,jj) = 0.79
              fuxy(ji,jj) = 0.10
              laixy(ji,jj,1) = 2.*1.00
              saixy(ji,jj,1) = 2.*2.03 - laixy(ji,jj,1)
              laixy(ji,jj,2) = 2.*0.23
              saixy(ji,jj,2) = 2.*0.41 - laixy(ji,jj,2)
            else if (iv.eq.10) then
              flxy(ji,jj) = 0.9
              fuxy(ji,jj) = 0.
              laixy(ji,jj,1) = 1.0
              saixy(ji,jj,1) = 3.0
              laixy(ji,jj,2) = 0.
              saixy(ji,jj,2) = 0.
            endif
c>>>>>>>>>>>>>>>>
            endif
c>>>>>>>>>>>>>>>>
  103     continue
        endif
c
c          If winter wheat, replace lower-story laixy,saixy,flxy,orievxy
c          and oriehxy with values derived from DS Table 5.
c
c          First calculate latitude index iwheat (n to s) for the DS
c          delta-15-deg zone, and seasonal indices it1w,it2w,wt1w,wt2w.
c          (The latter are indep of latitude, but are calculated here
c          inside loop 100 for clarity.)
c
        iwheat = 1 + (52.5-alatd)/15.
        iwheat = min (nwheat, max(1,iwheat))
c
        t = 12. * (nday-0.5)/365.
        it1w = t + 0.5
        it2w = it1w + 1
        wt1w = (it1w+.5) - t
        wt2w = 1.-wt1w
        if (it1w.lt. 1) it1w = 12
        if (it2w.gt.12) it2w = 1
c
        do 104 ji=1,nlon
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (ivegxy(ji,jj).ne.0) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          iv = ivegxy(ji,jj)
c
          if (iv.eq.12) then
            flxy(ji,jj) = wt1w*frawheat(iwheat,it1w)
     *                  + wt2w*frawheat(iwheat,it2w)
c
            laixy(ji,jj,1) = 2.*(  wt1w*gaiwheat(iwheat,it1w)
     *                           + wt2w*gaiwheat(iwheat,it2w) )
            saixy(ji,jj,1) = 2.*(  wt1w*taiwheat(iwheat,it1w)
     *                           + wt2w*taiwheat(iwheat,it2w) )
     *                       - laixy(ji,jj,1)
c
            zchi = wt1w*chiwheat(iwheat,it1w)
     *           + wt2w*chiwheat(iwheat,it2w)
            orievxy(ji,jj,1) = max (-zchi, 0.)
            oriehxy(ji,jj,1) = max ( zchi, 0.)
          endif
c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>
  104   continue
 
#ifndef gcm
c8881   Schwartz experiments
c8881   if (nint(xdate).ge.2000101.) then
c         ivegxy(1,1) = 2
c         iv = ivegxy(1,1)
c
c         flxy(1,1) = 1.
c         fuxy(1,1) = .1
c         fuxy(1,1) = .2
c
c         laixy(1,1,1) = .05
c         saixy(1,1,1) = 10.
c         laixy(1,1,2) = .05
c         saixy(1,1,2) = 1.
c         if (xdate.gt.1001231.) then
c           laixy(1,1,1) = 10.
c           saixy(1,1,1) = .05
c           laixy(1,1,2) = 5.
c           saixy(1,1,2) = 1.
c         endif
c
c         orievxy(1,1,1) = 0.5*(orievxy(1,1,1)+orievxy(1,1,2))
c         orievxy(1,1,2) = orievxy(1,1,1)
c         oriehxy(1,1,1) = 0.5*(oriehxy(1,1,1)+oriehxy(1,1,2))
c         oriehxy(1,1,2) = oriehxy(1,1,1)
c
c         z2lz(iv) = .25
c         z1lz(iv) = .02
c         z2uz(iv) = 10.
c         z1uz(iv) = 3.
c
c         tauvlz(iv,1) = 0.5*(tauvlz(iv,1)+tauvlz(iv,2))
c         tauvlz(iv,2) = tauvlz(iv,1)
c         taunlz(iv,1) = 0.5*(taunlz(iv,1)+taunlz(iv,2))
c         taunlz(iv,2) = taunlz(iv,1)
c         rhovlz(iv,1) = 0.5*(rhovlz(iv,1)+rhovlz(iv,2))
c         rhovlz(iv,2) = rhovlz(iv,1)
c         rhonlz(iv,1) = 0.5*(rhonlz(iv,1)+rhonlz(iv,2))
c         rhonlz(iv,2) = rhonlz(iv,1)
c       endif
#endif
c
c          Set remaining veg parameters using DS data
c
        do 106 ji=1,nlon
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (ivegxy(ji,jj).ne.0) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          iv = ivegxy(ji,jj)
c
          ztopxy(ji,jj,1) = z2lz(iv)
          zbotxy(ji,jj,1) = z1lz(iv)
          ztopxy(ji,jj,2) = z2uz(iv)
          zbotxy(ji,jj,2) = z1uz(iv)
c
c         set up live/dead weights for lower-story radiative trans,refl.
          tot = laixy(ji,jj,1) + saixy(ji,jj,1)
          if (tot.gt.0.) then
            wlivl = laixy(ji,jj,1) / tot
            wdedl = saixy(ji,jj,1) / tot
          else
            wlivl = 0.5
            wdedl = 0.5
          endif
c
c         set up live/dead weights for upper-story radiative trans,refl.
          tot = laixy(ji,jj,2) + saixy(ji,jj,2)
          if (tot.gt.0.) then
            wlivu = laixy(ji,jj,2) / tot
            wdedu = saixy(ji,jj,2) / tot
          else
            wlivu = 0.5
            wdedu = 0.5
          endif
c
c         transmittance and reflectance, visible waveband, lower story
          tauvegxy(ji,jj,1,1) = wlivl*tauvlz(iv,1) + wdedl*tauvlz(iv,2)
          rhovegxy(ji,jj,1,1) = wlivl*rhovlz(iv,1) + wdedl*rhovlz(iv,2)
c
c         transmittance and reflectance, visible waveband, upper story
          tauvegxy(ji,jj,1,2) = wlivu*tauvuz(iv,1) + wdedu*tauvuz(iv,2)
          rhovegxy(ji,jj,1,2) = wlivu*rhovuz(iv,1) + wdedu*rhovuz(iv,2)
c
c         transmittance and reflectance, near-ir waveband, lower story
          tauvegxy(ji,jj,2,1) = wlivl*taunlz(iv,1) + wdedl*taunlz(iv,2)
          rhovegxy(ji,jj,2,1) = wlivl*rhonlz(iv,1) + wdedl*rhonlz(iv,2)
c
c         transmittance and reflectance, near-ir waveband, upper story
          tauvegxy(ji,jj,2,2) = wlivu*taunuz(iv,1) + wdedu*taunuz(iv,2)
          rhovegxy(ji,jj,2,2) = wlivu*rhonuz(iv,1) + wdedu*rhonuz(iv,2)
 
c8880     adjust vis vs nir optical properties for .7 vs .9 microns.
c8880     ww = .75
c         xv =     ww *tauvegxy(ji,jj,1,1) + (1.-ww)*tauvegxy(ji,jj,2,1)
c         xn = (1.-ww)*tauvegxy(ji,jj,1,1) +     ww *tauvegxy(ji,jj,2,1)
c         tauvegxy(ji,jj,1,1) = xv
c         tauvegxy(ji,jj,2,1) = xn
c         xv =     ww *rhovegxy(ji,jj,1,1) + (1.-ww)*rhovegxy(ji,jj,2,1)
c         xn = (1.-ww)*rhovegxy(ji,jj,1,1) +     ww *rhovegxy(ji,jj,2,1)
c         rhovegxy(ji,jj,1,1) = xv
c         rhovegxy(ji,jj,2,1) = xn
c         xv =     ww *tauvegxy(ji,jj,1,2) + (1.-ww)*tauvegxy(ji,jj,2,2)
c         xn = (1.-ww)*tauvegxy(ji,jj,1,2) +     ww *tauvegxy(ji,jj,2,2)
c         tauvegxy(ji,jj,1,2) = xv
c         tauvegxy(ji,jj,2,2) = xn
c         xv =     ww *rhovegxy(ji,jj,1,2) + (1.-ww)*rhovegxy(ji,jj,2,2)
c         xn = (1.-ww)*rhovegxy(ji,jj,1,2) +     ww *rhovegxy(ji,jj,2,2)
c         rhovegxy(ji,jj,1,2) = xv
c         rhovegxy(ji,jj,2,2) = xn
 
c         PAR stomatal resistance dependence (parcomxy will be modified
c         in vegdat_all)
 
          parconxy(ji,jj,1) = aparlz(iv)/cparlz(iv) + bparlz(iv)
          rsmaxxy(ji,jj,1)  = aparlz(iv)/bparlz(iv) + cparlz(iv)
          rsminxy(ji,jj,1)  = cparlz(iv)
          parconxy(ji,jj,2) = aparuz(iv)/cparuz(iv) + bparuz(iv)
          rsmaxxy(ji,jj,2)  = aparuz(iv)/bparuz(iv) + cparuz(iv)
          rsminxy(ji,jj,2)  = cparuz(iv)
 
c         VPD stomatal resistance dependence (1 mb = 100 Pascals)
          vpdconxy(ji,jj,1) =  (1./h5lz(iv))*100.
          vpdconxy(ji,jj,2) =  (1./h5uz(iv))*100.
 
c         linear dimensions for aerodynamic flux parameterization
          dleafxy(ji,jj,1) = .01
          dstemxy(ji,jj,1) = .01
          dleafxy(ji,jj,2) = .05
          dstemxy(ji,jj,2) = .20
c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>
  106   continue
c
c-------------
  100 continue
c-------------
c
c        Set fraction of leaf area with stomates (fstomaxy).
c        Also set max root depths (rootdepth, passed to vegdat_all)
c        to 0.1 times the geometric canopy height. (This ratio should
c        really vary with veg type from about .01 to .1 for pine trees
c        to 1 or more for desert plants with tap roots.)
c
      do 200 m=1,2
        do 202 jj=1,nlat
          do 204 ji=1,nlon
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if (ivegxy(ji,jj).ne.0) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            fstomaxy(ji,jj,m) = 0.5
            rootdepth(ji,jj,m) = 0.1*ztopxy(ji,jj,m)
c>>>>>>>>>>>>>>>>
            endif
c>>>>>>>>>>>>>>>>
  204     continue
  202   continue
  200 continue
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vegdat_eve (jyear, jmonth, jday, fsave_leaf, showmaps,
     *                       rootdepth)
c
c        Sets current state of vegetation on lsx 2-D grid
c        using Eve/Leaf module.
c
c        jyear  = year number (eg, 1990) (supplied)
c        jmonth = month number (1-12) (supplied)
c        jday   = day number (1-31) nb:agcm has no leap years (supplied)
c        fsave_leaf = .true. if Genesis "save" this timestep (supplied)
c                     (not used with new Eve3)
c        showmaps   = .true. to print 2-D maps of prescribed flds (supp)
c        rootdepth  = maximum root depth for each story (returned)
c
#include <compar.h>
#include <combudl.h>
#include <comeve.h>
#include <comfud.h>
#include <comgrd.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commapl.h>
#include <commugxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
 
      logical fsave_leaf, showmaps
      dimension rootdepth(nlon,nlat,2)
 
      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm
 
      character concat*96, cmd*200, mname*96, mnamea*96
 
      logical lfirs
      save lfirs
      data lfirs /.true./
 
c     Solar properties for dead leaves:
      data rho_ded_vis, rho_ded_nir, tau_ded_vis, tau_ded_nir
     *     /.36, .58, .22, .38/
 
c     Solar properties for stems:
      data rho_ste_vis, rho_ste_nir, tau_ste_vis, tau_ste_nir
     *     /.16, .39,.001,.001/
 
      save rho_ded_vis, rho_ded_nir, tau_ded_vis, tau_ded_nir,
     *     rho_ste_vis, rho_ste_nir, tau_ste_vis, tau_ste_nir
 
      dimension
     *  laid_lsx(nlon,nlat,2),
     *  flb_lsx(nlon,nlat,2),
     *  xl_broad_lsx(nlon,nlat,2),
     *  xl_needl_lsx(nlon,nlat,2),
     *  xw_broad_lsx(nlon,nlat,2),
     *  xw_needl_lsx(nlon,nlat,2)
 
      real laid_lsx
 
#ifdef cray
      dimension maskeve(nlon,nlat)
#else
      integer*4 maskeve(nlon,nlat)
#endif
      save maskeve
c=====
c for calling eveinterp2 (also set in eveinterp2):
c      parameter (nloneve=180,nlateve=90)
c#ifdef cray
c     dimension maskeve(nloneve,nlateve)
celse
c     integer*4 maskeve(nloneve,nlateve)
c#endif
c     save maskeve
c     dimension
c    *  fallxy_eve(nloneve,nlateve,2),
c    *  ztopxy_eve(nloneve,nlateve,2),
c    *  zbotxy_eve(nloneve,nlateve,2),
c    *  saixy_eve(nloneve,nlateve,2),
c    *  dstemxy_eve(nloneve,nlateve,2),
c    *  rootdepth_eve(nloneve,nlateve,2),
c    *  laixy_eve(nloneve,nlateve,2),
c    *  laid_eve(nloneve,nlateve,2),
c    *  flb_eve(nloneve,nlateve,2),
c    *  xl_broad_eve(nloneve,nlateve,2),
c    *  xl_needl_eve(nloneve,nlateve,2),
c    *  xw_broad_eve(nloneve,nlateve,2),
c    *  xw_needl_eve(nloneve,nlateve,2),
c    *  orievxy_eve(nloneve,nlateve,2),
c    *  rhovegxy_eve(nloneve,nlateve,nband,2),
c    *  tauvegxy_eve(nloneve,nlateve,nband,2),
c    *  rsmaxxy_eve(nloneve,nlateve,2),
c    *  rsminxy_eve(nloneve,nlateve,2),
c    *  parconxy_eve(nloneve,nlateve,2),
c    *  vpdconxy_eve(nloneve,nlateve,2)
c     real laixy_eve, laid_eve
c=====
c     For sgi (and sun?), must pass integer args to C-code eve_leaf
c     as integer*4 ("int" in eve_leaf), to avoid -i8. Same for maskeve
c     above. Reals are "double" in eve_leaf, so arrays here are ok.
c     nb: Fortran integer*4 <-> C int,  Fortran real*4 <-> C float, 
c         Fortran integer*8 <-> C long, Fortran real*8 <-> C double.
#ifndef cray
      integer*4 
     *  ifinter_eve, ifdelta_eve, leafsave, iday_leaf, jyear_eve
#endif
 
c        Compute iday_leaf, number of days into current year for
c        Eve/Leaf (Jan 1st = 1, etc).
 
      iday_leaf = 0
      if (jmonth.gt.1) then
        do 50 jm = 1,jmonth-1
          iday_leaf = iday_leaf + ndaypm(jm)
   50   continue
      endif
      iday_leaf = iday_leaf + jday
 
 
c        Shift [t,p]eve_d by 180 deg in longitude and reverse latitudes,
c        for interactive eve/leaf.
 
      if (vegtype.eq.3) then
        do 70 j=1,nlat
          jj = nlat - j + 1
          do 72 i=1,nlon
            ii = mod (i + nlon/2 - 1, nlon) + 1
            workl(i,j) = teve_d(ii,jj)
            workm(i,j) = peve_d(ii,jj)
   72     continue
   70   continue
      endif
 
c     if (jday.le.2 .and. jmonth.eq.1) then
c       write(*,8881) jday, ((nint(workl(i,j)),i=1,nlon),j=1,nlat)
c8881   format(/'vegdat_eve: jday=',i3,'  TEVE_D:' / (180i3))
c       write(*,8882) jday, ((nint(workm(i,j)),i=1,nlon),j=1,nlat)
c8882   format(/'vegdat_eve: jday=',i3,'  PEVE_D:' / (180i3))
c     endif
 
c        Set leaf restart-file control flag (leafsave):
c        0 = Eve/Leaf should do normal calculations
c        1 = Eve/Leaf should read its restart file first thing
c            and then do normal calculations
 
      if (lfirs) then
        leafsave = 1
        lfirs = .false.
      else 
        leafsave = 0
      endif
 
c        Call Eve/Leaf module to return current day's veg fields.
c        Most lsx arrays are in common comvegxy, except rootdepth
c        (passed array) and laid_lsx, flb_lsx, x[l,w]_*_lsx (local
c        arrays) used later for setting comvegxy arrays. All veg arrays
c        *are* reset by eve_leaf every call even though some only
c        change at the beginning of a year.
 
c     write(*,8880) leafsave, iday_leaf, jyear
c8880 format(/'Vegdat: calling eve_leaf... leafsave,iday,jyear=',3i8)
 
      ifinter_eve = vegtype-2  ! non-int.(0)/interactive(1) flag(comgrd)
      ifdelta_eve = evedelta   ! no delta(0)/delta(1) flag (comgrd)
      jyear_eve = jyear
 
      if (ifinter_eve.eq.0 .and. ifdelta_eve.eq.1) then
        ifdelta_eve = 0
        write (ioterm2, 80)
   80   format('/ *** Error (vegdat_eve): cannot have non-interactive',
     *         ' EVE and climate deltas')
        call endrun (-1)
      endif
 
c     new (post 3/98) eve_leaf (Eve3): 
#if defined ( cray ) || defined ( sun ) || defined ( sgi ) || defined (linux)
      call eve_leaf (
#else
      call eve_leaf_ (
#endif
     *  ifinter_eve,     ! non-interective/interactive flag (0/1;suppl)
     *  ifdelta_eve,     ! no-delta/delta flag (0/1;suppl)
     *  leafsave,        ! leaf-file read flag (0/1; supplied)
     *  iday_leaf,       ! day number of current year, (1-365; supplied)
     *  jyear_eve,       ! year number (2001,...; supplied)
 
     *  workm,           ! precip,  total over prev 24 hrs, mm,   suppl.
     *  workl,           ! 2-m temp, mean over prev 24 hrs,deg C, suppl
 
     *  maskeve,         ! eve-leaf map (0=ocn,1=land,2=ice), returned
c===============
c not eveinterp2:
c===============
     *  fallxy,          ! frational cover of canopy
     *  ztopxy,          ! geometric height of canopy top (m)
     *  zbotxy,          ! geometric height of canopy bottom (m)
     *  saixy,           ! stem area index (m2/m2)
     *  dstemxy,         ! aero dimension of stems (m)
     *  rootdepth,       ! maximum root depth (m)
 
     *  laixy,           ! *green* lai, one-sided (m2/m2)
     *  laid_lsx,        ! *brown* (dead) lai, one-sided (m2/m2)
     *  flb_lsx,         ! fraction of total lai that is broadleaf
 
     *  xl_broad_lsx,    ! aero length dimension for broadleaf (m)
     *  xl_needl_lsx,    ! aero length dimension for needleleaf (m)
     *  xw_broad_lsx,    ! aero width  dimension for broadleaf (m)
     *  xw_needl_lsx,    ! aero width  dimension for needleleaf (m)
 
     *  orievxy,         ! leaf orientation (-1=vert, 0=random, 1=horiz)
 
     *  rhovegxy,        ! *live* leaf reflectivity,   vis and nir
     *  tauvegxy,        ! *live* leaf transmissivity, vis and nir
 
     *  rsmaxxy,         ! stomatal conductance, maximum (m/s)
     *  rsminxy,         ! stomatal conductance, minimum (m/s)
     *  parconxy,        ! stomatal PAR constant (W/m2)
     *  vpdconxy         ! stomatal VPD constant (N/m2)
c=======================
c for calling eveinterp2:
c=======================
c    *  fallxy_eve,      ! frational cover of canopy
c    *  ztopxy_eve,      ! geometric height of canopy top (m)
c    *  zbotxy_eve,      ! geometric height of canopy bottom (m)
c    *  saixy_eve,       ! stem area index (m2/m2)
c    *  dstemxy_eve,     ! aero dimension of stems (m)
c    *  rootdepth_eve,   ! maximum root depth (m)
c
c    *  laixy_eve,       ! *green* lai, one-sided (m2/m2)
c    *  laid_eve,        ! *brown* (dead) lai, one-sided (m2/m2)
c    *  flb_eve,         ! fraction of total lai that is broadleaf
c
c    *  xl_broad_eve,    ! aero length dimension for broadleaf (m)
c    *  xl_needl_eve,    ! aero length dimension for needleleaf (m)
c    *  xw_broad_eve,    ! aero width  dimension for broadleaf (m)
c    *  xw_needl_eve,    ! aero width  dimension for needleleaf (m)
c
c    *  orievxy_eve,     ! leaf orientation (-1=vert, 0=random, 1=horiz)
c
c    *  rhovegxy_eve,    ! *live* leaf reflectivity,   vis and nir
c    *  tauvegxy_eve,    ! *live* leaf transmissivity, vis and nir
c
c    *  rsmaxxy_eve,     ! stomatal conductance, maximum (m/s)
c    *  rsminxy_eve,     ! stomatal conductance, minimum (m/s)
c    *  parconxy_eve,    ! stomatal PAR constant (W/m2)
c    *  vpdconxy_eve     ! stomatal VPD constant (N/m2)
c===================
     *  )
 
c     write(*,8881)
c8881 format('Vegdat: returned from eve_leaf'/)
 
c        For old Eve rf/leaf_restart_file: not used with new Eve3:
c        If save time and using MSS, save rf/leaf_restart_file to MSS.
c        Then jobscript can mswrite leaf_restart_file back to rf/ for
c        non-checkpointed restarts, ie, if rf/leaf_restart_file doesn't
c        exist. (For non-MSS, Eve/Leaf does all this in rf/ itself).
 
c#ifdef gcm
c#ifdef mss
c      if (leafsave.eq.2 .and. mspath.ne.' ') then
cc       to generic name:
c        mname  = concat (mspath, 'leaf_restart_file', lenm)
c        call mswrite (ier, 'rf/leaf_restart_file', mname(1:lenm),
c     *                ' ', 2000, 'lcount=15')
cc       to individual (day-year) name:
cc       write (mnamea, "(a,'_',i3.3,'_',i4.4)")
cc    *    'leaf_restart', iday_leaf, max(-999,min(9999,jyear))
cc       mname = concat (mspath, mnamea, lenm)
cc       call mswrite (ier, 'rf/leaf_restart_file', mname(1:lenm),
cc    *                ' ', 2000, 'lcount=15')
c      endif
c#else
cc     if (leafsave.eq.2 .and. mspath.ne.' ' .and.mspath.ne.'./') then
cc       mname  = concat (mspath, 'leaf_restart_file', lenm)
cc       cmd = 'cp rf/leaf_restart_file ' // mname(1:lenm)
cc       ier = ishell (cmd(1:lenchr(cmd)))
cc     endif
c#endif
c#endif
 
c=====
c Transfer from eve-leaf map to current lsx map:
c     call eveinterp2 (fallxy_eve(1,1,1), flxy, maskeve)
c     call eveinterp2 (fallxy_eve(1,1,2), fuxy, maskeve)
c     do m=1,2
c       call eveinterp2 (ztopxy_eve(1,1,m),   ztopxy(1,1,m),    maskeve)
c       call eveinterp2 (zbotxy_eve(1,1,m),   zbotxy(1,1,m),    maskeve)
c       call eveinterp2 (saixy_eve(1,1,m),    saixy(1,1,m),     maskeve)
c       call eveinterp2 (dstemxy_eve(1,1,m),  dstemxy(1,1,m),   maskeve)
c       call eveinterp2 (rootdepth_eve(1,1,m),rootdepth(1,1,m), maskeve)
c       call eveinterp2 (laixy_eve(1,1,m),    laixy(1,1,m),     maskeve)
c       call eveinterp2 (laid_eve(1,1,m),     laid_lsx(1,1,m),  maskeve)
c       call eveinterp2 (flb_eve(1,1,m),      flb_lsx(1,1,m),   maskeve)
c       call eveinterp2(xl_broad_eve(1,1,m),xl_broad_lsx(1,1,m),maskeve)
c       call eveinterp2(xl_needl_eve(1,1,m),xl_needl_lsx(1,1,m),maskeve)
c       call eveinterp2(xw_broad_eve(1,1,m),xw_broad_lsx(1,1,m),maskeve)
c       call eveinterp2(xw_needl_eve(1,1,m),xw_needl_lsx(1,1,m),maskeve)
c       call eveinterp2 (orievxy_eve(1,1,m),  orievxy(1,1,m),   maskeve)
c       do ib=1,nband
c         call eveinterp2 (rhovegxy_eve(1,1,ib,m), rhovegxy(1,1,ib,m),
c    *                     maskeve)
c         call eveinterp2 (tauvegxy_eve(1,1,ib,m), tauvegxy(1,1,ib,m),
c    *                     maskeve)
c       enddo
c       call eveinterp2 (rsmaxxy_eve(1,1,m),  rsmaxxy(1,1,m),  maskeve)
c       call eveinterp2 (rsminxy_eve(1,1,m),  rsminxy(1,1,m),  maskeve)
c       call eveinterp2 (parconxy_eve(1,1,m), parconxy(1,1,m), maskeve)
c       call eveinterp2 (vpdconxy_eve(1,1,m), vpdconxy(1,1,m), maskeve)
c     enddo
c=====
 
c       Diagnostic printout of leaf quantities for the lsx output point
 
      if (.true. .and. lonbud.ne.-999 .and. latbud.ne.-999) then
c     if (.false.) then
        iuz = 57
        iz = lonbud
        jz = latbud
        zlon = ((iz-.5)/nlon)*360.
        if (zlon.gt.180.) zlon = zlon - 360.
        zlat = -90. + ((jz-.5)/nlat)*180.
        write(iuz,*)
        write(iuz,*) 'vegdat_eve: zlon=',zlon,'  zlat=',zlat
        write(iuz,8888) 'flxy,fuxy   ', flxy(iz,jz), fuxy(iz,jz)
        write(iuz,8888) 'ztopxy      ',(ztopxy(iz,jz,m), m=1,2)
        write(iuz,8888) 'zbotxy      ',(zbotxy(iz,jz,m), m=1,2)
        write(iuz,8888) 'saixyy      ',(saixy(iz,jz,m), m=1,2)
        write(iuz,8888) 'dstemxy     ',(dstemxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'rootdepth   ',(rootdepth(iz,jz,m),m=1,2)
        write(iuz,8888) 'laixy       ',(laixy(iz,jz,m),m=1,2)
        write(iuz,8888) 'laid_lsx    ',(laid_lsx(iz,jz,m),m=1,2)
        write(iuz,8888) 'flb_lsx     ',(flb_lsx(iz,jz,m),m=1,2)
        write(iuz,8888) 'xl_broad_lsx',(xl_broad_lsx(iz,jz,m),m=1,2)
        write(iuz,8888) 'xl_needl_lsx',(xl_needl_lsx(iz,jz,m),m=1,2)
        write(iuz,8888) 'xw_broad_lsx',(xw_broad_lsx(iz,jz,m),m=1,2)
        write(iuz,8888) 'xw_needl_lsx',(xw_needl_lsx(iz,jz,m),m=1,2)
        write(iuz,8888) 'orievxy     ',(orievxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'rhovegxy(1) ',(rhovegxy(iz,jz,1,m),m=1,2)
        write(iuz,8888) 'rhovegxy(2) ',(rhovegxy(iz,jz,2,m),m=1,2)
        write(iuz,8888) 'tauvegxy(1) ',(tauvegxy(iz,jz,1,m),m=1,2)
        write(iuz,8888) 'tauvegxy(2) ',(tauvegxy(iz,jz,2,m),m=1,2)
        write(iuz,8888) 'rsmaxxy     ',(rsmaxxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'rsminxy     ',(rsminxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'parconxy    ',(parconxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'vpdconxy    ',(vpdconxy(iz,jz,m),m=1,2)
 8888   format(a,':  ',2f20.10)
      endif
 
c       Ensure no vegetation over ocean or icesheet
 
      do 100 jj=1,nlat
        do 102 ji=1,nlon
          if (lmaskxy(ji,jj).ne.1) then
            flxy(ji,jj) = 0.
            fuxy(ji,jj) = 0.
          endif
  102   continue
  100 continue
 
c        Set harmless values for non-veg points or empty stories
c        (eve_leaf returns values as zero if fractional cover or lai/sai
c        are zero, but lsx needs "safe" values for land/landice pts
c        for both stories even if no veg (with lai = sai = fractional
c        cover = 0), since radiative and aero calcs are done whether
c        canopies exist or not, to allow Cray vectorization.
 
      do 150 m=1,2
        do 152 jj=1,nlat
          do 154 ji=1,nlon
 
            if ( (m.eq.1 .and. flxy(ji,jj).eq.0.) .or.
     *           (m.eq.2 .and. fuxy(ji,jj).eq.0.)      ) then
              ztopxy(ji,jj,m) = 0.
              zbotxy(ji,jj,m) = 0.
              dstemxy(ji,jj,m) = .10
              rootdepth(ji,jj,m) = .01
 
              flb_lsx(ji,jj,m) = 1.
              xl_broad_lsx(ji,jj,m) = .05
              xl_needl_lsx(ji,jj,m) = .001
              xw_broad_lsx(ji,jj,m) = .05
              xw_needl_lsx(ji,jj,m) = .001
              orievxy(ji,jj,m) = 0.
 
              rhovegxy(ji,jj,1,m) = .1
              rhovegxy(ji,jj,2,m) = .5
              tauvegxy(ji,jj,1,m) = .05
              tauvegxy(ji,jj,2,m) = .25
              rsmaxxy(ji,jj,m) = 1.
              rsminxy(ji,jj,m) = .01
              parconxy(ji,jj,m) = 100.
              vpdconxy(ji,jj,m) = 1.e4
            endif
 
c           If sai is zero, set stem dimension
            if (saixy(ji,jj,m).eq.0.) then
              dstemxy(ji,jj,m) = 1.
            endif
 
c           If live and dead lai are both zero, set "total" quantities
            if (laixy(ji,jj,m).eq.0. .and. laid_lsx(ji,jj,m).eq.0.) then
              flb_lsx(ji,jj,m) = 1.
              xl_broad_lsx(ji,jj,m) = .05
              xl_needl_lsx(ji,jj,m) = .001
              xw_broad_lsx(ji,jj,m) = .05
              xw_needl_lsx(ji,jj,m) = .001
              orievxy(ji,jj,m) = 0.
            endif
 
c           If live lai is zero, set "live" quantities
            if (laixy(ji,jj,m).eq.0.) then
              rhovegxy(ji,jj,1,m) = .1
              rhovegxy(ji,jj,2,m) = .5
              tauvegxy(ji,jj,1,m) = .05
              tauvegxy(ji,jj,2,m) = .25
              rsmaxxy(ji,jj,m) = 1.
              rsminxy(ji,jj,m) = .01
              parconxy(ji,jj,m) = 100.
              vpdconxy(ji,jj,m) = 1.e4
            endif
 
  154     continue
  152   continue
  150 continue
 
c        Constrain geometric canopy heights to avoid blowups in turlsx
c        subrs, and root depth to avoid blowup in frootxy calc below
 
      do 160 jj=1,nlat
        do 162 ji=1,nlon
          zbotxy(ji,jj,1) = max (zbotxy(ji,jj,1), .02)
          ztopxy(ji,jj,1) = max (ztopxy(ji,jj,1), zbotxy(ji,jj,1)+.01)
          zbotxy(ji,jj,2) = max (zbotxy(ji,jj,2), ztopxy(ji,jj,1)+.01)
          ztopxy(ji,jj,2) = max (ztopxy(ji,jj,2), zbotxy(ji,jj,2)+.01)
 
          rootdepth(ji,jj,1) = max (rootdepth(ji,jj,1), .001)
          rootdepth(ji,jj,2) = max (rootdepth(ji,jj,2), .001)
  162   continue
  160 continue
 
c         Convert various Eve/Leaf quantities into lsx forms
c         (m=1 for lower story, m=2 for upper story)
 
      do 170 m=1,2
        do 172 jj=1,nlat
          do 174 ji=1,nlon
 
c              Double all lai,sai's since lsx uses 2-sided values
 
            saixy(ji,jj,m) = 2.*saixy(ji,jj,m)
            laixy(ji,jj,m) = 2.*laixy(ji,jj,m)
            laid_lsx(ji,jj,m) = 2.*laid_lsx(ji,jj,m)
 
c              Reduce transpiration for dead leaves using fstomaxy
c              (0.5 is for live leaves having stoma on one side only)
 
            fstomaxy(ji,jj,m) = 0.5 * laixy(ji,jj,m)
     *                 / max (laixy(ji,jj,m)+laid_lsx(ji,jj,m), epsilon)
 
c              Lump dead leaves with stems for solar refl/trans weights
 
            ztot = max (laixy(ji,jj,m)+laid_lsx(ji,jj,m)+saixy(ji,jj,m),
     *                  epsilon)
            wliv= laixy(ji,jj,m) / ztot
            wded = laid_lsx(ji,jj,m) / ztot
            wste = 1.-wliv-wded
            rhovegxy(ji,jj,1,m) = wliv * rhovegxy(ji,jj,1,m)
     *                          + wded * rho_ded_vis
     *                          + wste * rho_ste_vis
            rhovegxy(ji,jj,2,m) = wliv * rhovegxy(ji,jj,2,m)
     *                          + wded * rho_ded_nir
     *                          + wste * rho_ste_nir
            tauvegxy(ji,jj,1,m) = wliv * tauvegxy(ji,jj,1,m)
     *                          + wded * tau_ded_vis
     *                          + wste * tau_ste_vis
            tauvegxy(ji,jj,2,m) = wliv * tauvegxy(ji,jj,2,m)
     *                          + wded * tau_ded_nir
     *                          + wste * tau_ste_nir
 
c              Now can sum live leaves and dead leaves into laixy
 
            laixy(ji,jj,m) = laixy(ji,jj,m) + laid_lsx(ji,jj,m)
 
c              Leaf orientation
 
            zchi = orievxy(ji,jj,m)
            orievxy(ji,jj,m) = max (-zchi, 0.)
            oriehxy(ji,jj,m) = max ( zchi, 0.)
 
c              For aerodynamic linear leaf dimension, use:
c              (.5*width + .5*length) for broadleaf,
c              (.9*width + .1*length) for needleleaf
 
            wbroad = flb_lsx(ji,jj,m)
            wneedl = 1.-wbroad
            dleafxy(ji,jj,m) =
     *       wbroad*(.5*xw_broad_lsx(ji,jj,m)+.5*xl_broad_lsx(ji,jj,m))
     *      +wneedl*(.9*xw_needl_lsx(ji,jj,m)+.1*xl_needl_lsx(ji,jj,m))
 
c              Convert from stomatal conductance to stomatal resistance
c              (and so switch max and min)
 
            zrsmax = rsmaxxy(ji,jj,m)
            rsmaxxy(ji,jj,m) = 1./rsminxy(ji,jj,m)
            rsminxy(ji,jj,m) = 1./zrsmax
 
  174     continue
  172   continue
  170 continue
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vegdat_biome (jyear, jmonth, jday, showmaps, rootdepth,
     *                         co2ppm)
c
c        Sets current state of vegetation on lsx 2-D grid, by
c        reading biome_out.cur file (biome#s and monthly lai's)
c        that was written by separate program biomedrive.exe.
c        The sets other veg attributes by calling vegdat_dorm for
c        closest corresponding D+S type.
c
c        jyear  = year number (eg, 1990) (supplied)
c        jmonth = month number (1-12) (supplied)
c        jday   = day number (1-31) nb:agcm has no leap years (supplied)
c        showmaps   = .true. to print 2-D maps of prescribed flds (supp)
c        rootdepth  = maximum root depth for each story (returned)
c        co2ppm = atmospheric CO2 concentration, ppmv
c
#include <compar.h>
#include <combudl.h>
#include <comeve.h>
#include <comfud.h>
#include <comgrd.h>
 
#ifdef gcm
#include <COMLUN.h>
#endif
 
#include <commapl.h>
#include <commugxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
 
      logical showmaps
      dimension rootdepth (nlon,nlat,2)
 
      parameter (nmon=12)
 
      dimension
     *  laibiome(nlon,nlat,nmon),  laibiome_max(nlon,nlat)
      real laibiome, laibiome_max
 
cc Biome3.5:
cc    Lookup table from Biome3.5 biome # to closest Dorman+Sellers 
cc    type #. First value of each pair is D+S index,
cc    second value of each pair is 1 for lower story, 2 for upper.
c     parameter (nbio=24)
c     dimension ibio2ds(2,nbio)
c     data ibio2ds /
cc        D+S   u/l         biome3.5
c    *      1,    2,        !  1
c    *      1,    2,        !  2
c    *      6,    2,        !  3
c    *      1,    2,        !  4
c    *      2,    2,        !  5
c    *      4,    2,        !  6
c    *      6,    2,        !  7
c    *      3,    2,        !  8
c    *      5,    2,        !  9
c    *      4,    2,        ! 10
c    *      5,    2,        ! 11
c    *      6,    2,        ! 12
c    *      6,    2,        ! 13
c    *      6,    2,        ! 14
c    *      7,    1,        ! 15
c    *      7,    1,        ! 16
c    *      9,    1,        ! 17
c    *     10,    1,        ! 18
c    *     10,    1,        ! 19
c    *     10,    1,        ! 20
c    *     10,    1,        ! 21
c    *     10,    1,        ! 22
c    *     11,    1,        ! 23
c    *     13,    1/        ! 24
 
cc Biome4:
cc    Lookup table from Biome4 biome # to closest Dorman+Sellers type #.
cc    First  value of each pair is D+S index,
cc    second value of each pair is 1 for lower story, 2 for upper.
      parameter (nbio=28)
      dimension ibio2ds(2,nbio)
      data ibio2ds /
cc        D+S   u/l         biome4
     *      1,    2,        !  1
     *      1,    2,        !  2
     *      2,    2,        !  3
     *      2,    2,        !  4
     *      4,    2,        !  5
     *      3,    2,        !  6
     *      3,    2,        !  7
     *      4,    2,        !  8
     *      3,    2,        !  9
     *      4,    2,        ! 10
     *      5,    2,        ! 11
     *      6,    2,        ! 12
     *      8,    2,        ! 13
     *      8,    2,        ! 14
     *      2,    2,        ! 15
     *      6,    2,        ! 16
     *      4,    2,        ! 17
     *     10,    1,        ! 18
     *      7,    1,        ! 19
     *      7,    1,        ! 20
     *     11,    1,        ! 21
     *     10,    1,        ! 22
     *     10,    1,        ! 23
     *     10,    1,        ! 24
     *     10,    1,        ! 25
     *     10,    1,        ! 26
     *     11,    1,        ! 27
     *     13,    1/        ! 28
      save laibiome, laibiome_max, ibio2ds
 
      dimension ndaypm(nmon), daymid(0:nmon+1)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daymid / -15.,  16.,  45.,  75., 105., 136., 166.,
     *               197., 228., 258., 289., 319., 350.,381. /
 
      logical first
      data first /.true./
 
      if (first) then
 
        write (ioterm2,*) 'Reading biomes and lais (biome_out.cur)'
        iu = 92
        open (iu, file='biome_out.cur', status='old')
 
        do j=1,nlat
          read (iu,'(2i3,180i4)') jj, kk, (ibiomexy(i,j),i=1,nlon)
          if (jj.ne.j .and. kk.ne.0) then
            write (ioterm2,'(a,3i4)') '*** Error: j,jj,kk=', j, jj, kk
            call endrun (-1)
          endif
          do i=1,nlon
            if (lmaskxy(i,j).eq.1) then
              if (ibiomexy(i,j).ge.1 .and. ibiomexy(i,j).le.nbio) then
                ivegxy(i,j) = ibio2ds(1,ibiomexy(i,j)) !closest D+S type
              else
                write (6,*) '*** Warning: biome=', ibiomexy(i,j),
     *                      '  at i=', i,'  j=',j, ' (land)'
                ibiomexy(i,j) = 0
                ivegxy(i,j) = nved
              endif
            else
              ibiomexy(i,j) = 0
              ivegxy(i,j) = nved
            endif
          enddo
        enddo

c          Impose specified D/S veg type in specified rectangles
c          (like dssoi[s,n] in soil.F)

        if (dsvegs.ne.0 .or. dsvegn.ne.0) then
          do j=1,nlat
            zlat = alat(j)*180./pi
            do i=1,nlon
              zlon = alon(i)*180./pi
              if (zlon.gt.180.) zlon = zlon - 360.
              if ( zlat.ge.dslats(1) .and. zlat.le.dslats(2) .and.
     *             zlon.ge.dslons(1) .and. zlon.le.dslons(2) .and.
     *             dsvegs.ne.0 .and. lmaskxy(i,j).eq.1 ) then
                ivegxy(i,j) = dsvegs
              endif
              if ( zlat.ge.dslatn(1) .and. zlat.le.dslatn(2) .and.
     *             zlon.ge.dslonn(1) .and. zlon.le.dslonn(2) .and.
     *             dsvegn.ne.0 .and. lmaskxy(i,j).eq.1 ) then
                ivegxy(i,j) = dsvegn
              endif
            enddo
          enddo
        endif

        do k=1,nmon
          do j=1,nlat
            read (iu,'(2i3,180f4.0)') jj, kk, (laibiome(i,j,k),i=1,nlon)
            if (jj.ne.j .and. kk.ne.k) then
              write (ioterm2,'(a,4i4)')
     *          '*** Error: j,jj,k,kk=', j, jj, k, kk
              call endrun (-1)
            endif
            do i=1,nlon
              if (lmaskxy(i,j).eq.1) then
                laibiome(i,j,k) = laibiome(i,j,k) / 100. !convert from %
              else
                laibiome(i,j,k) = 0.
              endif
            enddo
          enddo
        enddo
 
        close (iu)
 
c          Set seasonal max lai (just to constrain sai's below)
 
        do j=1,nlat
          do i=1,nlon
            laibiome_max(i,j) = -1.e20
            do k=1,nmon
             laibiome_max(i,j) = max (laibiome_max(i,j),laibiome(i,j,k))
            enddo
          enddo
        enddo
 
c          Print out veg-type maps if requested
 
c       if (showmaps) then
          call mapout (ioterm2, ibiomexy, nlon, nlat, alon, alat, 0,
     *                 'BIOME biomes')
          call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, nved,
     *                 'BIOME -> DS Vegetation Type')
c       endif
 
        first = .false.
 
      endif
 
c        Call Dorman and Sellers subr to set physical parameters
c        according to ivegxy (closest D+S type to biome #s)
 
      call vegdat_dorm (jyear, jmonth, jday, .false., showmaps,
     *                  rootdepth, 2)

c        Find time-interpolation indices (im1,im2,wm1) to mid months
      nday = 0
      do jm = 1,jmonth
        nday = nday + ndaypm(jm)
      enddo
      zday = nday - ndaypm(jmonth) + jday
 
      do im=1,nmon+1
        if (zday.le.daymid(im)) then
          im1 = im-1
          im2 = im
          wm1 = (daymid(im2)-zday) / (daymid(im2)-daymid(im1))
          wm2 = 1.-wm1
          if (im1.eq.0) im1 = nmon
          if (im2.eq.nmon+1) im2 = 1
          goto 100
        endif
      enddo
  100 continue
 
c        Overwrite D+S lai and fractional cover by biome values.
c        Convert biome's single-sided lai's to two sided, and
c        use biome's fractional cover = min (1., single-sided lai / 2.)

c        Don't do if ifimage != 0 (i.e.,agriculture imposed in 
c        vegdat_dorm), or if in a rectangle with specified D/S veg type
 
      do j=1,nlat
        do i=1,nlon
          if (lmaskxy(i,j).eq.1 .and. ifimage(i,j).eq.0) then
 
            ifspecveg = 0
            if (dsvegs.ne.0 .or. dsvegn.ne.0) then
              zlat = alat(j)*180./pi
              zlon = alon(i)*180./pi
              if (zlon.gt.180.) zlon = zlon - 360.
              if ( zlat.ge.dslats(1) .and. zlat.le.dslats(2) .and.
     *             zlon.ge.dslons(1) .and. zlon.le.dslons(2) .and.
     *             dsvegs.ne.0 ) ifspecveg = 1
              if ( zlat.ge.dslatn(1) .and. zlat.le.dslatn(2) .and.
     *             zlon.ge.dslonn(1) .and. zlon.le.dslonn(2) .and.
     *             dsvegn.ne.0 ) ifspecveg = 1
            endif

            if (ifspecveg.eq.0) then
              zlai  = wm1*laibiome(i,j,im1) + wm2*laibiome(i,j,im2)
              zfrac = min (1., zlai/2.)
              zsaim = 0.5*laibiome_max(i,j)
              if (ibio2ds(2,ibiomexy(i,j)).eq.1) then
c               all lower story:
                flxy(i,j)    = zfrac
                laixy(i,j,1) = 2.*zlai
                saixy(i,j,1) = min (saixy(i,j,1), 2.*zsaim)
                fuxy(i,j)    = 0.
                laixy(i,j,2) = 0.
                saixy(i,j,2) = 0.
              else
c               all upper story:
                fuxy(i,j)    = zfrac
                laixy(i,j,2) = 2.*zlai
                saixy(i,j,2) = min (saixy(i,j,2), 2.*zsaim)
                flxy(i,j)    = 0.
                laixy(i,j,1) = 0.
                saixy(i,j,1) = 0.
              endif
            endif
 
          endif
        enddo
      enddo
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vegdat_sdvm (jyear, jmonth, jday, showmaps, rootdepth)

c     jyear  = year number (eg, 1990) (supplied)
c     jmonth = month number (1-12) (supplied)
c     jday   = day number (1-31) nb:agcm has no leap years (supplied)
c     showmaps = .true. to print 2-D maps of prescribed fields (supp)
c     rootdepth= maximum root depth for each story (returned)

c Obtained from vegdat_dorm.04.F,
c in /home/bac/proc/presday/simpveg/vegdat.u/ver2/
c
c Sets current state of vegetation for SDVM model (Brian Cosgrove)
c on lsx 2-D grid, based on condensed groups of Dorman and Sellers.
c Vegetation values have been condensed down into 5 groups.
c This condensation was accomplished through weighting the
c original vegetation values according to their area coverage. These
c groups consist of groups 1 and 3,4 and 2,5 and 6,7,8,9,10 and 11.
c Perpetual ice is left unaltered in its own category.
c These five groups were constructed on the basis of deciduousness
c Winter wheat values have not been weighted since they do not appear
c on the base natural vegetation map and are not predicted by the
c dynamic vegetation model--hence wheat is not used at all and should 
c not be used as a biome.

#include <compar.h>
#include <comfud.h>
#include <comgrd.h>

#ifdef gcm
#include <COMLUN.h>
#endif

#include <commapl.h>
#include <commugxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>

      logical showmaps
      dimension rootdepth (nlon,nlat,2)

      dimension ndaypm(12), saiuzav(nved)
      logical lfirs
      save ndaypm, saiuzav, lfirs
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data lfirs /.true./

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c comved.h contains D-S specific arrays, initialized by data
c statements below. Declared as local (dimension...) and saved,
c not in common, to allow different data settings in subr vegdat_sdvm.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <comved.h>
c
c     total lai+sai, one-sided, lower story
      data tailz/
     *   0,   0,   0,   0,   0, 0.9, 0.9, 0.9, 0.9, 0.9,   0, 1.3,   0,
     *   0,   0,   0,   0,   0, 0.9, 0.9, 0.9, 0.9, 0.9,   0, 1.5,   0,
     *   0,   0,   0,   0,   0, 0.9, 0.9, 0.9, 0.9, 0.9,   0, 1.7,   0,
     *   0,   0,   0,   0,   0, 0.8, 0.8, 0.8, 0.8, 0.8,   0, 1.9,   0,
     *   0,   0,   0,   0,   0, 0.9, 0.9, 0.9, 0.9, 0.9,   0, 3.5,   0,
     *   0,   0,   0,   0,   0, 1.5, 1.5, 1.5, 1.5, 1.5,   0, 4.7,   0,
     *   0,   0,   0,   0,   0, 1.8, 1.8, 1.8, 1.8, 1.8,   0, 0.3,   0,
     *   0,   0,   0,   0,   0, 2.3, 2.3, 2.3, 2.3, 2.3,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 1.9, 1.9, 1.9, 1.9, 1.9,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 1.4, 1.4, 1.4, 1.4, 1.4,   0, 0.6,   0,
     *   0,   0,   0,   0,   0, 1.2, 1.2, 1.2, 1.2, 1.2,   0, 0.8,   0,
     *   0,   0,   0,   0,   0, 0.9, 0.9, 0.9, 0.9, 0.9,   0, 1.1,   0/
c
c     green lai, one-sided, lower story
      data gailz/
     *   0,   0,   0,   0,   0, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 1.1,   0,
     *   0,   0,   0,   0,   0, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 1.3,   0,
     *   0,   0,   0,   0,   0, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 1.5,   0,
     *   0,   0,   0,   0,   0, 0.4, 0.4, 0.4, 0.4, 0.4,   0, 1.7,   0,
     *   0,   0,   0,   0,   0, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 3.1,   0,
     *   0,   0,   0,   0,   0, 1.1, 1.1, 1.1, 1.1, 1.1,   0, 2.8,   0,
     *   0,   0,   0,   0,   0, 1.3, 1.3, 0.3, 0.3, 0.3,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 1.3, 1.3, 1.3, 1.3, 1.3,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 1.1, 1.1, 1.1, 1.1, 1.1,   0, 0.0,   0,
     *   0,   0,   0,   0,   0, 0.8, 1.8, 1.8, 1.8, 1.8,   0, 0.5,   0,
     *   0,   0,   0,   0,   0, 0.6, 0.6, 0.6, 0.6, 0.6,   0, 0.7,   0,
     *   0,   0,   0,   0,   0, 0.6, 0.6, 0.6, 0.6, 0.6,   0, 1.0,   0/
c
c     fractional cover, lower story
      data flz/
     *   0,   0,   0,   0,   0, .49, .49, .49, .49, .59,   0, .65,   0,
     *   0,   0,   0,   0,   0, .50, .50, .50, .50, .50,   0, .70,   0,
     *   0,   0,   0,   0,   0, .50, .50, .50, .50, .50,   0, .73,   0,
     *   0,   0,   0,   0,   0, .49, .49, .49, .49, .49,   0, .77,   0,
     *   0,   0,   0,   0,   0, .50, .50, .50, .50, .50,   0, .88,   0,
     *   0,   0,   0,   0,   0, .58, .58, .58, .58, .58,   0, .90,   0,
     *   0,   0,   0,   0,   0, .59, .59, .59, .59, .59,   0, .20,   0,
     *   0,   0,   0,   0,   0, .67, .67, .67, .67, .67,   0, .00,   0,
     *   0,   0,   0,   0,   0, .61, .61, .61, .61, .61,   0, .00,   0,
     *   0,   0,   0,   0,   0, .56, .56, .56, .56, .56,   0, .41,   0,
     *   0,   0,   0,   0,   0, .54, .54, .54, .54, .54,   0, .50,   0,
     *   0,   0,   0,   0,   0, .51, .51, .51, .51, .51,   0, .58,   0/
c
c     total lai+sai, one-sided, upper story
      data taiuz/
     * 5.0, 0.3, 5.5, 5.5, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4,   0, 0.0,   0,
     * 5.0, 0.3, 5.7, 5.7, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4,   0, 0.0,   0,
     * 5.0, 0.5, 6.1, 6.1, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3,   0, 0.0,   0,
     * 5.0, 1.2, 6.6, 6.6, 1.2, 0.2, 0.2, 0.2, 0.2, 0.2,   0, 0.2,   0,
     * 5.0, 2.3, 7.0, 7.0, 2.3, 0.2, 0.2, 0.2, 0.2, 0.2,   0, 0.3,   0,
     * 5.0, 3.4, 7.4, 7.4, 3.4, 0.2, 0.2, 0.2, 0.2, 0.2,   0, 0.5,   0,
     * 5.0, 4.7, 7.3, 7.3, 4.7, 0.7, 0.7, 0.7, 0.7, 0.7,   0, 0.5,   0,
     * 5.0, 3.6, 7.3, 7.3, 3.6, 1.1, 1.1, 1.1, 1.1, 1.1,   0, 0.5,   0,
     * 5.0, 2.7, 6.9, 6.9, 2.7, 1.1, 1.1, 1.1, 1.1, 1.1,   0, 0.4,   0,
     * 5.0, 1.5, 6.3, 6.3, 1.5, 0.6, 0.6, 0.6, 0.6, 0.6,   0, 0.2,   0,
     * 5.0, 0.7, 5.5, 5.5, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 0.1,   0,
     * 5.0, 0.5, 5.3, 5.3, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4,   0, 0.0,   0/
c
c     green lai, one-sided, upper story
      data gaiuz/
     * 4.5, 0.0, 4.9, 4.9, 0.0, 0.3, 0.3, 0.3, 0.3, 0.3,   0, 0.0,   0,
     * 4.5, 0.0, 5.1, 5.1, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2,   0, 0.0,   0,
     * 4.5, 0.1, 5.5, 5.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,   0, 0.0,   0,
     * 4.5, 0.8, 6.0, 6.0, 0.8, 0.1, 0.1, 0.1, 0.1, 0.1,   0, 0.2,   0,
     * 4.5, 1.9, 6.4, 6.4, 1.9, 0.1, 0.1, 0.1, 0.1, 0.1,   0, 0.3,   0,
     * 4.5, 3.0, 6.7, 6.7, 3.0, 0.2, 0.2, 0.2, 0.2, 0.2,   0, 0.5,   0,
     * 4.5, 3.3, 6.6, 6.6, 3.3, 0.7, 0.7, 0.7, 0.7, 0.7,   0, 0.4,   0,
     * 4.5, 2.3, 6.3, 6.3, 2.3, 1.9, 1.9, 1.9, 1.9, 1.9,   0, 0.3,   0,
     * 4.5, 1.1, 5.7, 5.7, 1.1, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 0.1,   0,
     * 4.5, 0.4, 5.0, 5.0, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3,   0, 0.0,   0,
     * 4.5, 0.1, 4.8, 4.8, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3,   0, 0.0,   0,
     * 4.5, 0.0, 4.8, 4.8, 0.0, 0.3, 0.3, 0.3, 0.3, 0.3,   0, 0.0,   0/
c
c     fractional cover, upper story
      data fuz/
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0,
     * .98, .60, .75, .75, .60, .16, .16, .16, .16, .16,   0, .08,   0/
c
c *** Transmittances ***
c
c     VIS leaf transmittance, lower story (1st row = live, 2nd = dead)
      data tauvlz/
     * .07, .07, .07, .07, .07, .07, .07, .07, .07, .07, .07, .07, .07,
     * .22, .22, .22, .22, .22, .22, .22, .22, .22, .22, .22, .22, .22/
c
c     VIS leaf transmittance, upper story (1st row = live, 2nd = dead)
      data tauvuz/
     * .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05,
     *.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001/
c
c     NIR leaf transmittance, lower story (1st row = live, 2nd = dead)
      data taunlz/
     * .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25,
     * .38, .38, .38, .38, .38, .38, .38, .38, .38, .38, .38, .38, .38/
c
c     NIR leaf transmittance, upper story (1st row = live, 2nd = dead)
      data taunuz/
     * .25, .16, .12, .12, .16, .25, .25, .25, .25, .25, .25, .25, .25,
     *.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001/
c
c *** Reflectances ***
c
c     VIS leaf reflectance, lower story (1st row = live, 2nd = dead)
      data rhovlz/
     * .11, .11, .11, .11, .11, .11, .11, .11, .11, .11, .11, .11, .11,
     * .36, .36, .36, .36, .36, .36, .36, .36, .36, .36, .36, .36, .36/
c
c     VIS leaf reflectance, upper story (1st row = live, 2nd = dead)
      data rhovuz/
     * .10, .08, .07, .07, .08, .10, .10, .10, .10, .10, .10, .10, .10,
     * .16, .16, .16, .16, .16, .16, .16, .16, .16, .16, .16, .16, .16/
c
c     NIR leaf reflectance, lower story (1st row = live, 2nd = dead)
      data rhonlz/
     * .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58,
     * .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58, .58/
c
c     NIR leaf reflectance, upper story (1st row = live, 2nd = dead)
      data rhonuz/
     * .45, .39, .37, .37, .39, .45, .45, .45, .45, .45, .45, .45, .45,
     * .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39/
c
c *** Vegetation angle orientation (-1=vert, 0=random, 1=horiz) ***
c
c     chi factor, lower story
      data chiflz/
     * .00, .00, .00, .00, .00,-.11,-.11,-.11,-.11,-.11, .00,-.02, .00/
c
c     chi factor, upper story
      data chifuz/
     * .10, .10, .05, .05, .10, .08, .08, .08, .08, .08, .00, .25, .00/
c
c *** Vegetation heights ***
c
c     vegetation base, lower story (m) 
      data z1lz/
     * .02, .02, .02, .02, .02, .02, .02, .02, .02, .02, .02, .02, .02/
c
c     vegetation top, lower story (m)
      data z2lz/
     * .02, .02, .02, .02, .02, .30, .30, .30, .30, .30, .02, .60, .02/
c
c     vegetation base, upper story (m) 
      data z1uz/
     *1.00,8.52,8.71,8.71,8.52,1.66,1.66,1.66,1.66,1.66, .03, 11., .03/
c
c     vegetation top, upper story (m) 
      data z2uz/
     *35., 16.3, 18.1, 18.1, 16.3, 3.1, 3.1,3.1,3.1,3.1,.03, 20. , .03/
c
c
c *** PAR coefficients ***
c
c     "a", lower story (J m-3)
      data aparlz/
     *  2582,  2582,  2582,  2582,  2582,  2582,
     *  2582,  2582,  2582,  2582,  2582,  7459,  2582/
c
c     "b", lower story (W m-2)
      data bparlz/
     *  1.09,  1.09,  1.09,  1.09,  1.09,  1.09,
     *  1.09,  1.09,  1.09,  1.09,  1.09,  5.70,  1.09/
c
c     "c", lower story (s m-1)
      data cparlz/
     *   110,   110,   110,   110,   110,   110,
     *   110,   110,   110,   110,   110,    25,   110/
c
c     "a", upper story (J m-3)
      data aparuz/
     *  2336,  5504,  4100,  4100,  5504, 38487,
     * 38487, 38487, 38487, 38487, 10000,  9802, 10000/
c
c     "b", upper story (W m-2)
      data bparuz/
     *  0.01,  6.32,  4.91,  4.91,  6.32,  7.67,
     *  7.67,  7.67,  7.67,  7.67,    10, 10.60,    10/
c
c     "c", upper story (s m-1)
      data cparuz/
     *   154,   213,   224,   224,   213,   330,
     *   330,   330,   330,   330,   100,   180,   100/
c
c *** vpd coefficient ***
c
c     "h5", lower story (mb-1)
      data h5lz/
     * .0238, .0238, .0238, .0238, .0238, .0238,
     * .0238, .0238, .0238, .0238, .0238, .0001, .0238/
c
c     "h5", upper story (mb-1)
      data h5uz/
     * .0273, .0328, .0321, .0321, .0328, .0234,
     * .0234, .0234, .0234, .0234, .0100, .0357, .0100/
c
c *** winter wheat (from DS Table 5 and table5 program) ***
c
c     total lai+sai, one-sided, for winter wheat, vs (lat_zone, month)
      data taiwheat/
     *     1.21,    1.46,    4.33,    0.79,    0.00,    0.00,    3.27,
     *     1.41,    1.76,    3.84,    0.00,    0.00,    0.00,    0.00,
     *     1.60,    2.81,    1.00,    0.00,    0.00,    0.00,    0.25,
     *     1.80,    4.90,    0.00,    0.00,    0.46,    0.60,    0.63,
     *     1.99,    2.14,    0.00,    0.00,    1.06,    0.93,    0.82,
     *     5.97,    0.00,    0.00,    0.00,    1.64,    1.26,    1.02,
     *     3.24,    0.00,    0.00,    0.00,    4.23,    1.59,    1.21,
     *     0.00,    0.00,    0.00,    0.00,    3.82,    1.93,    1.41,
     *     0.24,    0.00,    0.00,    0.00,    1.08,    4.19,    1.60,
     *     0.63,    0.59,    0.47,    0.22,    0.00,    4.85,    1.80,
     *     0.82,    0.88,    1.06,    5.83,    0.00,    2.16,    1.99,
     *     1.02,    1.17,    1.65,    3.44,    0.00,    0.00,    5.92/
c
c     green lai, one-sided, for winter wheat, vs (lat_zone, month)
      data gaiwheat/
     *     1.09,    1.32,    3.54,    0.08,    0.00,    0.00,    0.98,
     *     1.27,    1.58,    1.78,    0.00,    0.00,    0.00,    0.00,
     *     1.44,    2.40,    0.10,    0.00,    0.00,    0.00,    0.23,
     *     1.62,    3.26,    0.00,    0.00,    0.42,    0.54,    0.57,
     *     1.79,    0.38,    0.00,    0.00,    0.95,    0.84,    0.74,
     *     4.76,    0.00,    0.00,    0.00,    1.48,    1.14,    0.92,
     *     0.93,    0.00,    0.00,    0.00,    3.47,    1.44,    1.09,
     *     0.00,    0.00,    0.00,    0.00,    1.75,    1.73,    1.27,
     *     0.21,    0.00,    0.00,    0.00,    0.12,    3.44,    1.44,
     *     0.57,    0.53,    0.43,    0.20,    0.00,    3.19,    1.62,
     *     0.74,    0.79,    0.96,    4.56,    0.00,    0.39,    1.79,
     *     0.92,    1.05,    1.48,    1.21,    0.00,    0.00,    4.69/
c
c     fractional cover for winter wheat, vs (lat_zone, month)
      data frawheat/
     *     0.59,    0.66,    0.86,    0.24,    0.00,    0.00,    0.55,
     *     0.64,    0.73,    0.67,    0.00,    0.00,    0.00,    0.00,
     *     0.69,    0.82,    0.30,    0.00,    0.00,    0.00,    0.20,
     *     0.75,    0.79,    0.00,    0.00,    0.37,    0.43,    0.43,
     *     0.80,    0.41,    0.00,    0.00,    0.55,    0.52,    0.49,
     *     0.90,    0.00,    0.00,    0.00,    0.70,    0.60,    0.54,
     *     0.55,    0.00,    0.00,    0.00,    0.86,    0.69,    0.59,
     *     0.00,    0.00,    0.00,    0.00,    0.66,    0.78,    0.64,
     *     0.19,    0.00,    0.00,    0.00,    0.31,    0.85,    0.69,
     *     0.43,    0.42,    0.38,    0.18,    0.00,    0.78,    0.75,
     *     0.49,    0.50,    0.55,    0.88,    0.00,    0.42,    0.80,
     *     0.54,    0.58,    0.71,    0.59,    0.00,    0.00,    0.89/
c
c     "chi" orientation factor for winter wheat, vs (lat_zone, month)
      data chiwheat/
     *    -0.05,   -0.06,   -0.13,   -0.24,    0.00,    0.00,   -0.28,
     *    -0.06,   -0.08,   -0.24,    0.00,    0.00,    0.00,    0.00,
     *    -0.07,   -0.10,   -0.30,    0.00,    0.00,    0.00,   -0.01,
     *    -0.08,   -0.20,    0.00,    0.00,   -0.02,   -0.02,   -0.03,
     *    -0.09,   -0.30,    0.00,    0.00,   -0.05,   -0.04,   -0.04,
     *    -0.16,    0.00,    0.00,    0.00,   -0.07,   -0.06,   -0.04,
     *    -0.28,    0.00,    0.00,    0.00,   -0.13,   -0.07,   -0.05,
     *     0.00,    0.00,    0.00,    0.00,   -0.24,   -0.09,   -0.06,
     *    -0.01,    0.00,    0.00,    0.00,   -0.30,   -0.13,   -0.07,
     *    -0.03,   -0.02,   -0.02,   -0.01,    0.00,   -0.20,   -0.08,
     *    -0.04,   -0.04,   -0.05,   -0.17,    0.00,   -0.30,   -0.09,
     *    -0.04,   -0.05,   -0.07,   -0.27,    0.00,    0.00,   -0.16/
c
c
c        If first call of this run, set constants, read veg-type data
c        file derived from DS Fig. 1b, put into ivegxy, and set ivegxy
c        values to nved for non-land points and for erroneous values.
c        Note that nominal values for all lsx veg params are set for 
c        non-land pts (ivegxy=nved, lai = sai = fractional cover = 0),
c        since radiative and aero calcs are done whether vegetation
c        exists or not, to allow Cray vectorization. Same applies for
c        land points with only one story, and for bare ground.
c
      if (lfirs) then
c
c          Compute annual mean upper-story SAI, used below to avoid
c          non-conservation of energy due to addition/removal of
c          wood mass at different seasonal temps.
c
        do iv=1,nved
          saiuzav(iv) = 0.
          do it=1,ntimds
            saiuzav(iv) = saiuzav(iv) + (taiuz(iv,it)-gaiuz(iv,it))
          enddo
          saiuzav(iv) = saiuzav(iv) / ntimds
        enddo
c
        if (dsbiome.eq.0) then
          call findiu (iuvegt, loc)
#ifdef gcm
          write(ioterm2,40) locname(loc)
   40     format(/' Reading SDVM vegetation data file'/4x,a)
#endif
          call openf (iuvegt, locname(loc), ' ', 0, 0)
          call readdata (iuvegt, workl, nlon, nlat, 'VEGSDVM', 1, 1)
          call closef(iuvegt)
          do 42 jj=1,nlat
            do 44 ji=1,nlon
              ivegxy(ji,jj) = workl(ji,jj) + .001
              if (lmaskxy(ji,jj).ne.1) then
                ivegxy(ji,jj) = nved
              else
                if (ivegxy(ji,jj).lt.1 .or. ivegxy(ji,jj).gt.nved) 
     *            ivegxy(ji,jj) = nved
              endif
   44       continue
   42     continue
        else
c         Namelist DSBIOME entered for global uniform single type:
          write(ioterm2,50) dsbiome
   50     format(/' *** Using single SDVM vegetation biome',
     *            i4,' for all land points')
          do 52 jj=1,nlat
            do 54 ji=1,nlon
              if (lmaskxy(ji,jj).eq.1) then
                ivegxy(ji,jj) = dsbiome
              else
                ivegxy(ji,jj) = nved
              endif
   54       continue
   52     continue
        endif

c          Print out veg-type map if requested

        if (showmaps)
     *    call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, nved,
     *                 'Vegetation Type')
c
        lfirs = .false.
      endif
c
c        Compute ndaynor, number of days in current year since Jan 1st,
c        ndaysou = ndaynor shifted by 6 months for SH, and ndayequ for
c        the equator, used below for time indices into Figs 2(a)-(o).
c
      nday = 0
      do 60 jm = 1,jmonth
        nday = nday + ndaypm(jm)
   60 continue
      nday = nday - ndaypm(jmonth) + jday
      ndaynor = nday
      ndaysou = mod (nday-1+365/2, 365) + 1
      ndayequ = nday
c
c---------------------
      do 100 jj=1,nlat
c---------------------
c
c          Set LAI, SAI, fractional cover and orientation weights,
c          changing DS one-sided values to LSX two-sided LAI,SAI. 
c          These are done separately from other parameters since they
c          are different for winter wheat (ivegxy=12).
c          (Note many lsx veg arrays have last subscript = 1 for lower
c          story, = 2 for upper.)
c
c          First, compute seasonal indices it1,it2,wt1,wt2 for DS
c          Figs 2a-o, linearly weighting between nday[nor,sou] and 
c          ndayequ for latitudes equatorward of 25 deg [N,S].
c
        alatd = alat(jj)*180./pi
c       
        w = min (1., abs(alatd/25.))
        if (alatd.ge.0.) then
          d = w*ndaynor + (1.-w)*ndayequ
        else
          d = w*ndaysou + (1.-w)*ndayequ
        endif
c
        t = 12. * (d-0.5)/365.
        it1 = t + 0.5
        it2 = it1 + 1
        wt1 = (it1+.5) - t
        wt2 = 1.-wt1
        if (it1.lt. 1) it1 = 12
        if (it2.gt.12) it2 = 1
c
        do 102 ji=1,nlon
          iv = ivegxy(ji,jj)
c
          flxy(ji,jj) = wt1*flz(iv,it1) + wt2*flz(iv,it2)
          fuxy(ji,jj) = wt1*fuz(iv,it1) + wt2*fuz(iv,it2)
c
          laixy(ji,jj,1) = 2.*(wt1*gailz(iv,it1) + wt2*gailz(iv,it2))
          saixy(ji,jj,1) = 2.*(wt1*tailz(iv,it1) + wt2*tailz(iv,it2))
     *                     - laixy(ji,jj,1)
          laixy(ji,jj,2) = 2.*(wt1*gaiuz(iv,it1) + wt2*gaiuz(iv,it2))
          saixy(ji,jj,2) = 2.*saiuzav(iv)

#ifndef gcm
c8880     Schwartz experiments
c8880     if (nint(xdate).ge.1000101.) then
c           flxy(ji,jj) = 0.
c           fuxy(ji,jj) = ((nint(xdate)/10000)-100)*.01
c         endif
#endif

          orievxy(ji,jj,1) = max (-chiflz(iv), 0.)
          orievxy(ji,jj,2) = max (-chifuz(iv), 0.)
          oriehxy(ji,jj,1) = max ( chiflz(iv), 0.)
          oriehxy(ji,jj,2) = max ( chifuz(iv), 0.)
  102   continue
c
c          If winter wheat, replace lower-story laixy,saixy,flxy,orievxy
c          and oriehxy with values derived from DS Table 5.
c
c          First calculate latitude index iwheat (n to s) for the DS
c          delta-15-deg zone, and seasonal indices it1w,it2w,wt1w,wt2w.
c          (The latter are indep of latitude, but are calculated here
c          inside loop 100 for clarity.)
c
        iwheat = 1 + (52.5-alatd)/15.
        iwheat = min (nwheat, max(1,iwheat))
c
        t = 12. * (nday-0.5)/365.
        it1w = t + 0.5
        it2w = it1w + 1
        wt1w = (it1w+.5) - t
        wt2w = 1.-wt1w
        if (it1w.lt. 1) it1w = 12
        if (it2w.gt.12) it2w = 1
c
        do 104 ji=1,nlon
          iv = ivegxy(ji,jj)
          if (iv.eq.12) then
            flxy(ji,jj) = wt1w*frawheat(iwheat,it1w)
     *                  + wt2w*frawheat(iwheat,it2w)

            laixy(ji,jj,1) = 2.*(  wt1w*gaiwheat(iwheat,it1w)
     *                           + wt2w*gaiwheat(iwheat,it2w) )
            saixy(ji,jj,1) = 2.*(  wt1w*taiwheat(iwheat,it1w)
     *                           + wt2w*taiwheat(iwheat,it2w) )
     *                       - laixy(ji,jj,1)
 
            zchi = wt1w*chiwheat(iwheat,it1w) 
     *           + wt2w*chiwheat(iwheat,it2w)
            orievxy(ji,jj,1) = max (-zchi, 0.)
            oriehxy(ji,jj,1) = max ( zchi, 0.)
          endif
  104   continue
c
c          Set remaining veg parameters using DS data
c
        do 106 ji=1,nlon
          iv = ivegxy(ji,jj)
c
          ztopxy(ji,jj,1) = z2lz(iv)
          zbotxy(ji,jj,1) = z1lz(iv)
          ztopxy(ji,jj,2) = z2uz(iv)
          zbotxy(ji,jj,2) = z1uz(iv)
c
c         set up live/dead weights for lower-story radiative trans,refl.
          tot = laixy(ji,jj,1) + saixy(ji,jj,1)
          if (tot.gt.0.) then
            wlivl = laixy(ji,jj,1) / tot
            wdedl = saixy(ji,jj,1) / tot
          else
            wlivl = 0.5
            wdedl = 0.5
          endif
c
c         set up live/dead weights for upper-story radiative trans,refl.
          tot = laixy(ji,jj,2) + saixy(ji,jj,2)
          if (tot.gt.0.) then
            wlivu = laixy(ji,jj,2) / tot
            wdedu = saixy(ji,jj,2) / tot
          else
            wlivu = 0.5
            wdedu = 0.5
          endif
c
c         transmittance and reflectance, visible waveband, lower story
          tauvegxy(ji,jj,1,1) = wlivl*tauvlz(iv,1) + wdedl*tauvlz(iv,2)
          rhovegxy(ji,jj,1,1) = wlivl*rhovlz(iv,1) + wdedl*rhovlz(iv,2)
c
c         transmittance and reflectance, visible waveband, upper story
          tauvegxy(ji,jj,1,2) = wlivu*tauvuz(iv,1) + wdedu*tauvuz(iv,2)
          rhovegxy(ji,jj,1,2) = wlivu*rhovuz(iv,1) + wdedu*rhovuz(iv,2)
c
c         transmittance and reflectance, near-ir waveband, lower story
          tauvegxy(ji,jj,2,1) = wlivl*taunlz(iv,1) + wdedl*taunlz(iv,2)
          rhovegxy(ji,jj,2,1) = wlivl*rhonlz(iv,1) + wdedl*rhonlz(iv,2)
c
c         transmittance and reflectance, near-ir waveband, upper story
          tauvegxy(ji,jj,2,2) = wlivu*taunuz(iv,1) + wdedu*taunuz(iv,2)
          rhovegxy(ji,jj,2,2) = wlivu*rhonuz(iv,1) + wdedu*rhonuz(iv,2)
 
c8880     adjust vis vs nir optical properties for .7 vs .9 microns.
c8880     ww = .75
c         xv =     ww *tauvegxy(ji,jj,1,1) + (1.-ww)*tauvegxy(ji,jj,2,1)
c         xn = (1.-ww)*tauvegxy(ji,jj,1,1) +     ww *tauvegxy(ji,jj,2,1)
c         tauvegxy(ji,jj,1,1) = xv
c         tauvegxy(ji,jj,2,1) = xn
c         xv =     ww *rhovegxy(ji,jj,1,1) + (1.-ww)*rhovegxy(ji,jj,2,1)
c         xn = (1.-ww)*rhovegxy(ji,jj,1,1) +     ww *rhovegxy(ji,jj,2,1)
c         rhovegxy(ji,jj,1,1) = xv
c         rhovegxy(ji,jj,2,1) = xn
c         xv =     ww *tauvegxy(ji,jj,1,2) + (1.-ww)*tauvegxy(ji,jj,2,2)
c         xn = (1.-ww)*tauvegxy(ji,jj,1,2) +     ww *tauvegxy(ji,jj,2,2)
c         tauvegxy(ji,jj,1,2) = xv
c         tauvegxy(ji,jj,2,2) = xn
c         xv =     ww *rhovegxy(ji,jj,1,2) + (1.-ww)*rhovegxy(ji,jj,2,2)
c         xn = (1.-ww)*rhovegxy(ji,jj,1,2) +     ww *rhovegxy(ji,jj,2,2)
c         rhovegxy(ji,jj,1,2) = xv
c         rhovegxy(ji,jj,2,2) = xn
 
c         PAR stomatal resistance dependence (parcomxy will be modified
c         in vegdat_all)

          parconxy(ji,jj,1) = aparlz(iv)/cparlz(iv) + bparlz(iv)
          rsmaxxy(ji,jj,1)  = aparlz(iv)/bparlz(iv) + cparlz(iv)
          rsminxy(ji,jj,1)  = cparlz(iv)
          parconxy(ji,jj,2) = aparuz(iv)/cparuz(iv) + bparuz(iv)
          rsmaxxy(ji,jj,2)  = aparuz(iv)/bparuz(iv) + cparuz(iv)
          rsminxy(ji,jj,2)  = cparuz(iv)

c         VPD stomatal resistance dependence (1 mb = 100 Pascals)
          vpdconxy(ji,jj,1) =  (1./h5lz(iv))*100.
          vpdconxy(ji,jj,2) =  (1./h5uz(iv))*100.

c         linear dimensions for aerodynamic flux parameterization
          dleafxy(ji,jj,1) = .01
          dstemxy(ji,jj,1) = .01
          dleafxy(ji,jj,2) = .05
          dstemxy(ji,jj,2) = .20
c
  106   continue
c
c-------------
  100 continue
c-------------
c
c        Set fraction of leaf area with stomates (fstomaxy). 
c        Also set max root depths (rootdepth, passed to vegdat_all)
c        to 0.1 times the geometric canopy height. (This ratio should
c        really vary with veg type from about .01 to .1 for pine trees
c        to 1 or more for desert plants with tap roots.)
c
      do 200 m=1,2
        do 202 jj=1,nlat
          do 204 ji=1,nlon
            fstomaxy(ji,jj,m) = 0.5
            rootdepth(ji,jj,m) = 0.1*ztopxy(ji,jj,m)
  204     continue
  202   continue
  200 continue
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vegdat_all (rootdepth)
c
c        Sets vegetation parameters that are the same for
c        both Dorman+Sellers and Eve/Leaf.
c
c        rootdepth= maximum root depth for each story (supplied)
c
#include <compar.h>
#include <comgrd.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
c
      dimension rootdepth(nlon,nlat,2)
c
      dimension d0root(nlon),  dtroot(nlon),
     *          aroot(nlon),   broot(nlon),  croot(nlon),
     *          xcroot(nlon),  ycroot(nlon),
     *          v1root(nlon),  v2root(nlon), vtroot(nlon),
     *          za(nlon),      zb(nlon),
     *          va(nlon),      vb(nlon)
 
 
c        Reduce all veg solar rho,tau by arbitrary factor:
c
c     do 50 m=1,2
c       do 52 jj=1,nlat
c         do 54 ji=1,nlon
c           rhovegxy(ji,jj,1,m) = 0.7*rhovegxy(ji,jj,1,m)
c           rhovegxy(ji,jj,2,m) = 0.7*rhovegxy(ji,jj,2,m)
c           tauvegxy(ji,jj,1,m) = 0.7*tauvegxy(ji,jj,1,m)
c           tauvegxy(ji,jj,2,m) = 0.7*tauvegxy(ji,jj,2,m)
c  54     continue
c  52   continue
c  50 continue
 
c        Reduce parcon (PAR stomatal-resistance normalizing factor)
c        by [1 - vis. single-scattering albedo], to allow approximately
c        for use of absorbed PAR, not incident (see SOLARF and TURTRA).
c        (Photosynthesis model will need absorbed PAR.)
c        Do only for BATS-SiB stomatal resistance param
c        (not used for IBIS photosynthesis/stomatal conductance)
 
      if (phototype.eq.0) then
        do 100 m=1,2
          do 102 jj=1,nlat
            do 104 ji=1,nlon
              parconxy(ji,jj,m) = parconxy(ji,jj,m)
     *                     *(1.-rhovegxy(ji,jj,1,m)-tauvegxy(ji,jj,1,m))
  104       continue
  102     continue
  100   continue
      endif
 
c        Set frootxy = fraction of roots in each soil layer
c        (vertical sum = 1).
c
c        The derivation of the expression below assumes power-law
c        conical geometry and uniform root density. For a linear cone,
c        aroot = 3, but for fairly deep root depths (>.93 m), the
c        profile has a higher power in z.
c
c        The ratio of max root depth to plant height should
c        vary with veg type from about .01 to .1 for pine trees to 1
c        or more for desert plants with tap roots for instance.
c        For Dorman+Sellers it is assumed to be .1 for all veg types.
c        This has been set earlier, passed as max root depths(rootdepth)
 
c        Don't do this for Wisconsin veg (frootxy has been set in
c        in vegdat_wisc/wisc_veg)
 
c---------------------------
      if (vegtype.ne.1) then
c---------------------------
 
      do 200 m=1,2
        do 202 jj=1,nlat
 
c            Set root depths and various depth-independent quantities
 
          do 210 ji=1,nlon
            d0root(ji) = rootdepth(ji,jj,m)
            dtroot(ji) = d0root(ji) - min (.02, .1*d0root(ji))
            dcroot     = min (.3, .3*d0root(ji))
            xcroot(ji) = (d0root(ji)-dcroot)/dtroot(ji)
            zza = 1.-3.*xcroot(ji)
            zzb = 1.-xcroot(ji)
            zzc = zzb
            ycroot(ji) = (-zzb - sqrt(zzb**2-4.*zza*zzc)) / (2.*zza)
            aroot(ji) = ycroot(ji)/xcroot(ji)
            broot(ji) = (1.-ycroot(ji))/(1.-xcroot(ji))
            croot(ji) = 1.-broot(ji)
            v1root(ji) = aroot(ji)**2 * xcroot(ji)**3
            v2root(ji) = (1.-ycroot(ji)**3) / broot(ji)
            vtroot(ji) = v1root(ji) + v2root(ji)
            za(ji) = 0.
            va(ji) = 1.
  210     continue
 
c            Step downwards though layers setting root amounts frootxy
 
          do 220 k=1,nsoilay
            do 222 ji=1,nlon
              zb(ji) = za(ji) + hsod(k)
              xb = max (0., min (1., (d0root(ji)-zb(ji))/dtroot(ji) ))
 
              if (xb.le.xcroot(ji)) then
                 vb(ji) = aroot(ji)**2 * xb**3 / vtroot(ji)
              else
                 vb(ji) = (  v1root(ji)
     *                       + (  (broot(ji)*xb + croot(ji))**3
     *                           - ycroot(ji)**3
     *                         ) / broot(ji)
     *                     ) / vtroot(ji)
              endif
 
              vb(ji) = cvmgt ( vb(ji), 0., k.lt.nsoilay )
 
              frootxy(ji,jj,k,m) = (va(ji) - vb(ji))
              za(ji) = zb(ji)
              va(ji) = vb(ji)
  222       continue
  220     continue
 
  202   continue
  200 continue
 
c----------
      endif
c----------
 
c        Same values for all points (for now)
c        Set tmpopt,temcon only for BATS-SiB stomatal resistance param
c        (not used for IBIS photosynthesis/stomatal conductance)
c
      if (phototype.eq.0) then
        do 300 m=1,2
          do 302 jj=1,nlat
            do 304 ji=1,nlon
              temconxy(ji,jj,m) = .0016
              temoptxy(ji,jj,m) = 298.
  304       continue
  302     continue
  300   continue
      endif
c
c        Same values for all points (ch[u,s,l] = ch2o * mm of liq H2O)
c
      chu = ch2o * 1.0        ! was 0.25 mm in early v2
      chs = ch2o * 10.        ! was 50   mm in early v2, 200 mm in v1.02
      chl = ch2o * 1.0
      alaimu = 12. ! was 6 in v1.02
      alaiml = 12. ! was 6 in v1.02
      blaim  = 12. ! was 6 in v1.02
      cleaf = .01
      cstem = .01
      cgrass = .01
      wliqumax = .1
      wliqsmax = .2  ! .1
      wliqlmax = .1
      wsnoumax = 2.
      wsnosmax = 2.
      wsnolmax = 2.
c     Esentially no drip below max capacity, but keep tblow (snow):
      tdripu = 1.e8
      tdrips = 1.e8
      tdripl = 1.e8
      tblowu = 2.*.864e5
      tblows = 2.*.864e5
      tblowl = 2.*.864e5
c
c        Some values are different for Wisconsin veg
c
      if (vegtype.eq.1) then
        chu = ch2o *  2.0
        chs = ch2o * 20.0
        chl = ch2o *  2.0
        alaimu = 10.
        alaiml = 10.
        wliqsmax = .1
        tdripu = 3.*3600.
        tdrips = 3.*3600.
        tdripl = 3.*3600.
      endif
c
      return
      end

c----------------------------------------------------------------------

c     nb: if eveinterp2 called (uncommented calls) above, for 
c     nloneve,nlateve different from nlon,nlat, then choose either 
c     subr eveinterp2_near or subr eveinterp2_aggr to be named 
c     subr eveinterp2.

      subroutine eveinterp2_near (arr_eve, arr_lsx, maskeve)
c     subroutine eveinterp2      (arr_eve, arr_lsx, maskeve)
 
c        Sets lsx-land points in arr_lsx (lsx array) to the value of
c        arr_eve (eve/leaf array) at nearest 2x2 eve/leaf land point.
c        Called from vegdat_eve. Eve/leaf's land/ocean mask is
c        supplied in maskeve.
c
#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
 
      parameter (nloneve=180, nlateve=90)
      dimension arr_eve(nloneve,nlateve), arr_lsx(nlon,nlat)
#ifdef cray
      dimension maskeve(nloneve,nlateve)
#else
      integer*4 maskeve(nloneve,nlateve)
#endif
 
      dimension inear(nlon,nlat), jnear(nlon,nlat)
 
      logical first
      save first, inear, jnear
      data first /.true./
 
c--------------------
      if (first) then
c--------------------
 
        write(99,*) 'maskeve (eveinterp2):'
        do j = nlateve,1,-1
          write(99,'(360i1)') (maskeve(i,j),i=1,nloneve)
        enddo
 
c          Set lookup indices (inear, jnear) pointing to closest
c          eve-leaf land point for each lsx land point.
c            (ilsx,jlsx) is the lsx land point.
c            (ieve,jeve) is the closest eve-leaf grid pt to (ilsx,jlsx).
 
        do jlsx=1,nlat
          zlat = -90. + ((jlsx-.5)/nlat)*180.
          jeve = nint(nlateve*(zlat+90.)/180. + .5)
          jeve = max (1, min (nlateve, jeve))
 
          zdx = 2.*pi*cos(zlat*pi/180.) / nloneve
          zdy = pi / nlateve
          zrat = zdx/zdy
 
          do ilsx=1,nlon
 
            if (lmaskxy(ilsx,jlsx).eq.1) then
 
              zlon = ((ilsx-.5)/nlon)*360.
              ieve = nint(nloneve*zlon/360. + .5)
              ieve = max (1, min (nloneve, ieve))
 
c                Search for nearest eve/leaf land pt (maskeve = 1), on
c                perimeters of widening squares centered on (ieve,jeve),
c                with sides 2*jdist + 1 in latitude and 2*idist + 1 in
c                longitude, with idist and jdist adjusted (using zrat)
c                so that sides are approx equal lengths on a sphere.
c                  (j1,j2) are j-limits of current square.
c                  (i1,i2) are i-limits of current square.
c                  Wrap around in longitude and "reflect" in latitude.
 
              do jdist=0,nlateve
                idist = min (nint(jdist/zrat), nloneve)
 
                j1 = jeve - jdist
                j2 = jeve + jdist
                i1 = ieve - idist
                i2 = ieve + idist
 
                do jx = j1,j2
 
                  if (jx.eq.j1 .or. jx.eq.j2) then
                    i3 = 1                 ! search along whole E-W side
                  else
                    i3 = max(i2-i1,1)      ! search only the 2 E-W edges
                  endif
 
                  do ix = i1,i2,i3
 
c                   wraparound in longitude:
                    i = mod (ix + 2*nloneve - 1, nloneve) + 1
 
c                   "reflect" in latitude:
                    if (jx.gt.nlateve) then
                      j = max (1, min (nlateve, 2*nlateve - jx + 1))
                      i = mod (i + nloneve/2 - 1, nloneve) + 1
                    else if (jx.lt.1) then
                      j = max (1, min (nlateve,           - jx + 1))
                      i = mod (i + nloneve/2 - 1, nloneve) + 1
                    else
                      j = jx
                    endif
 
                    if (maskeve(i,j).eq.1) then
                      inear(ilsx,jlsx) = i
                      jnear(ilsx,jlsx) = j
                      goto 40
                    endif
 
                  enddo   ! ix
                enddo     ! jx
 
              enddo       ! jdist
              write(*,*)
     *          '*** Search fail (eveinterp2) ilsx=',ilsx,' jlsx=',jlsx
              call endrun (-1)
   40         continue
 
            else          ! lmaskxy ne 1
 
              inear(ilsx,jlsx) = 0
              jnear(ilsx,jlsx) = 0
 
            endif
 
          enddo           ! ilsx
        enddo             ! jlsx
 
        first = .false.
 
c----------
      endif
c----------
 
c       Copy data from eve-leaf array to land points of lsx array
 
      do 100 j=1,nlat
        do 102 i=1,nlon
          if (lmaskxy(i,j).eq.1) then
            arr_lsx(i,j) = arr_eve ( inear(i,j), jnear(i,j) )
          else
            arr_lsx(i,j) = 0.
          endif
  102   continue
  100 continue
 
      return
      end

c----------------------------------------------------------------------

c     nb: if eveinterp2 called (uncommented calls) above, for 
c     nloneve,nlateve different from nlon,nlat, then choose either 
c     subr eveinterp2_near or subr eveinterp2_aggr to be named 
c     subr eveinterp2.

c     subroutine eveinterp2_aggr (arr_eve, arr_lsx, maskeve)
      subroutine eveinterp2      (arr_eve, arr_lsx, maskeve)
 
c        Interpolates arr_eve (eve/leaf array) to lsx land points 
c        in arr_lsx (lsx array). Assumes eve/leaf dimensions are finer
c        than lsx, and divide exactly into them, both running w->e 
c        and s->n, same longitude cut and both with grid box edges
c        at SP and NP.

#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
 
      parameter (nloneve=180, nlateve=90)
      dimension arr_eve(nloneve,nlateve), arr_lsx(nlon,nlat)
#ifdef cray
      dimension maskeve(nloneve,nlateve)
#else
      integer*4 maskeve(nloneve,nlateve)
#endif

      parameter (lonrat = nloneve/nlon, latrat = nlateve/nlat)
 
c     logical first
c     save first
c     data first /.true./
c     if (first) then
c       write(99,*) 'maskeve (eveinterp2):'
c       do j = nlateve,1,-1
c         write(99,'(360i1)') (maskeve(i,j),i=1,nloneve)
c       enddo
c       first = .false.
c     endif
 
      do jlsx=1,nlat
        jevea = (jlsx-1)*latrat + 1
        jeveb = jevea + latrat - 1
        do ilsx=1,nlon
          ievea = (ilsx-1)*lonrat + 1
          ieveb = ievea + lonrat - 1

          if (lmaskxy(ilsx,jlsx).eq.1) then
            zav = 0.
            nav = 0
            do j = jevea,jeveb
              do i = ievea,ieveb
                if (maskeve(i,j).eq.1) then 
                  zav = zav + arr_eve(i,j)
                  nav = nav + 1
                endif
              enddo
            enddo
            if (nav.gt.0) then
              arr_lsx(ilsx,jlsx) = zav/nav
            else
              write (6,'(/2a,2i4)')
     *         '*** Error (eveinterp2): no eve_leaf land in lsx box',
     *         ' i,j=',ilsx,jlsx
            endif

          else
            arr_lsx(ilsx,jlsx) = 0.
          endif

        enddo
      enddo
 
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine table5
c     program    table5
c
c        This is a separate program, stored here for convenience.
c        It is not to be linked into the lsx program.
c
c        Linearly interpolates data in Dorman and Sellers (1989),
c        JAM,28,833, Table 5, to print out mid-monthly seasonal curves
c        of winter-wheat total LAI, green LAI, fractional coverage,
c        and "chi" orientation factor, for each latitude zone. The data
c        for fractional coverage, and the bounds of the periods with
c        LAI=coverage=0, have been estimated from DS Fig.2(o). The
c        values printed here can then be edited into data statements
c        for the lsx subr VEGDAT.
c
      parameter (klat=7, kvar=4, ktim=8, kmon=12)
c
      dimension tim(ktim,klat), val(ktim,kvar), valout(kmon,kvar,klat),
     *          tmon(kmon)
      character*20 ctitl(kvar)
      data ctitl /'total LAI', 'green LAI', 'fractional cover',
     *            'chi factor'/
c
      data tim/
     *  249., 269., 138., 167., 187., 201., 207., 227.,
     *  260., 280.,  71.,  96., 115., 133., 142., 162.,
     *  270., 290.,   3.,  24.,  42.,  65.,  76.,  96.,
     *  280., 300., 301., 318., 334., 362.,  11.,  31.,
     *   88., 108., 186., 207., 224., 247., 259., 279.,
     *   77.,  97., 235., 278., 297., 316., 324., 344.,
     *   66.,  86., 321., 349.,   5.,  19.,  24.,  44./
c
      data val/
     *    0.,  0.5,  2.0,  6.0,  4.0,  3.0,  1.0,   0.,
     *    0.,  .45,  1.8,  4.8,  2.0,  0.6,  0.1,   0.,
     *    0.,  0.4,  0.8,  0.9,  0.7,  0.5,  0.3,   0.,
     *    0., -.02, -.09, -.16, -.23, -.30, -.30,   0./
c
c
      year = 365.
      do 50 imon=1,kmon
         tmon(imon) = (imon-0.5)*(year/kmon)
   50 continue
c
      do 100 ilat=1,klat
      do 100 ivar=1,kvar
      do 100 imon=1,kmon
        t = tmon(imon)
c
        do 110 it=1,ktim
          itp1 = it + 1
          if (it.eq.ktim) itp1 = 1
c
          if (t.ge.tim(it,ilat) .and. t.le.tim(itp1,ilat)) then
            ita = it
            itb = itp1
            wta = (tim(itp1,ilat)-t) / (tim(itp1,ilat)-tim(it,ilat))
            goto 120
c
          else if (t.ge.tim(it,ilat) .and.
     *             tim(it,ilat).gt.tim(itp1,ilat)) then
c           t is near year end, after last tim
            ita = it
            itb = itp1
            wta = (tim(itp1,ilat)+year-t)
     *          / (tim(itp1,ilat)+year-tim(it,ilat))
            goto 120
c
          else if (t.le.tim(itp1,ilat) .and.
     *             tim(it,ilat).gt.tim(itp1,ilat)) then
c           t is near year start, before first tim
            ita = it
            itb = itp1
            wta = (tim(itp1,ilat)-t)
     *          / (tim(itp1,ilat)+year-tim(it,ilat))
            goto 120
          endif
c
  110   continue
        write(6,*) ' Error 110. ilat,ivar,imon: ',ilat,ivar,imon
        call endrun (-1)
c
  120   if (wta.gt.1.) write(6,*) ' Error 120. ilat,ivar,imon,wta: ',
     *                            ilat,ivar,imon,wta
c
        valout(imon,ivar,ilat) =     wta *val(ita,ivar)
     *                         + (1.-wta)*val(itb,ivar)
c
  100 continue
c
c        Print out interpolated curves, roughly in data stmt format,
c        with months (Jan to Dec) running vertically and latitudes
c        (45N to 45S) horizontally, in four blocks for total lai,
c        green lai, fractional coverage and chi factor in that order.
c
      do 200 ivar=1,kvar
        write(6,210) ctitl(ivar)
  210   format(/1x,'      ',a)
        do 202 imon=1,kmon
          write(6,212) (valout(imon,ivar,ilat), ilat=1,klat)
  212     format('     * ',7(f8.2,','))
  202  continue
  200 continue
c
      call endrun (0)
      end
c
c----------------------------------------------------------------------
c
#if ! defined (gcm) || defined (rcm) || defined (hydmodel)
#if defined ( cray ) || defined ( sun ) || defined ( sgi ) || defined (linux) 
      subroutine eve_leaf ()
#else
      subroutine eve_leaf_ ()
#endif
      write(*,*) '*** Error: called stub eve_leaf'
      call endrun (-1)
      end
#endif
 
#ifndef ibis
      subroutine vegdat_wisc (l1,i1,i2,i3,i4,i5,c1)
      write(*,*) '*** Error: called stub vegdat_wisc'
      call endrun (-1)
      end
c     subroutine stomata (x1,x2,i1,i2)
c     write(*,*) '*** Error: called stub stomata'
c     call endrun (-1)
c     end
      subroutine drystress (x1,x2,x3,x4)
      dimension x1(*), x2(*), x3(*), x4(*)
      write(*,*) '*** Error: called stub drystress'
      call endrun (-1)
      end
c     subroutine ibis_stat (x1,i1,i2,i3,i4,i5,i6,i7,i8,i9)
c     write(*,*) '*** Error: called stub ibis_stat'
c     call endrun (-1)
c     end
#endif
 
