c -- $Id: ibis.F,v 1.2 2001/07/03 21:29:28 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/ibis.F,v $
c -- $Name:  $

c--------
c ibis2.0
c--------

c//////////
#ifdef ibis
c//////////

c Entire file is surrounded by "ifdef ibis".
c There are also "ifdef ibisint" segments inside for dynamic veg.
c Contains all code for:
c   IBIS photosynthesis  (PHOTOTYPE=1),
c   prescribed Wisconsin vegetation (VEGTYPE=1),
c   dynamic IBIS vegetation (DYNAIBIS=1)
c
c-----------------------------------------------------------------------
c
      subroutine stomata (dtstom, co2ppm, iter, niter)

c New IBIS version 1.1 photosynthesis (from Jon Foley, Univ.Wisconsin).
c (Sets stomatal conductances gs[sun,sha][u,l], in comveg).

c     dtstom = time step for this call, secs (supplied)
c     co2ppm = atmospheric (canopy-air) co2 amount, ppmv (supplied)
c     iter   = current iteration number in turlsx/turvap (supplied)
c     niter  = final   iteration number in turlsx/turvap (supplied)
c
#include <compar.h>
#include <comgrd.h>
#include <commug.h>
#include <commus.h>
#include <comsoi.h>
#include <comveg.h>
#include <com1d.h>
c
c local variables
c
      dimension
     >  agsunud(nvec), agsunuc(nvec),  agsunl3(nvec), agsunl4(nvec),
     >  agshaud(nvec), agshauc(nvec),  agshal3(nvec), agshal4(nvec),
     >  ansunud(nvec), ansunuc(nvec),  ansunl3(nvec), ansunl4(nvec),
     >  anshaud(nvec), anshauc(nvec),  anshal3(nvec), anshal4(nvec)

      real kc, ko, kco2, kc15, ko15,
     >     je, jc, js, jp, ji
c
c *******************************************************  
c model parameters 
c *******************************************************  
c
c atmospheric oxygen concentration (mol/mol)
c
      data o2conc  /0.209000/
c
c intrinsic quantum efficiency for c3 and c4 plants (dimensionless)
c
      data alpha3 /0.065/
      data alpha4 /0.050/
c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
      data tau15 /4500.0/
c
c o2/co2 kinetic parameters (mol/mol)
c
      data kc15 /1.5e-04/
      data ko15 /2.5e-01/
c
c photosynthesis coupling coefficients (dimensionless)
c
      data theta3 /0.90/
      data beta3  /0.90/
c
      data theta4 /0.90/
      data beta4  /0.90/
c
#include <comsat.h>
c
      co2conc = co2ppm*1.e-6
c
c *******************************************************  
c begin stomatal physiology model - for all points
c *******************************************************  
c
      do 100 i = 1, npoi
c
c *******************************************************  
c calculate boundary layer parameters (mol/m**2/s)
c *******************************************************  
c
        gbco2u = min(10., max(.10, su(i) * 25.5)) !gbco2=su/0.029*1.35
        gbco2l = min(10., max(.10, sl(i) * 25.5))
c 
c
c
c *******************************************************************
c calculate photosynthesis rates for upper sunny leaves - broadleaf
c *******************************************************************
c
        tleaf = tu(i) - 273.16
c
c calculate the parameter values which are a fnction of temperature
c
        tempvm = func1(tleaf) *
     >           exp(3500.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
        tau    = tau15 * 
     >           exp(-5000. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        kc     = kc15 *
     >           exp(6000.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
        ko     = ko15 *
     >           exp(1400.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxsunu(i,1) * tempvm
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.030 * vmax * max (0.50, uptotu(i))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstarud = min (0.9 * cisunud(i), o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpu(i) * 4.59e-06 * alpha3 * (cisunud(i) - gamstarud) / 
     >       (cisunud(i) + 2. * gamstarud)
c
        jc = vmax * (cisunud(i) - gamstarud) / 
     >       (cisunud(i) + kc * (1. + o2conc / ko))

c solve for intermediate photosynth. rate

        duma = theta3
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c sucrose synthesis limited rate of photosynth. (mol/m**2/s)

        js = 3.0 * vmax / 8.2 * (1.0 - gamstarud / cisunud(i)) +
     >       jp * gamstarud / cisunud(i)

        duma = beta3
        dumb = jp + js
        dumc = jp * js

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15


        agsunud(i) = min(dumq/duma, dumc/dumq) * uptotu(i)

c net photosynthesis (mol/m**2/s)
c
        ansunud(i) = agsunud(i) - rdark
c
c ********************************************************************
c calculate photosynthesis rates for upper shaded leaves - broadleaf
c ********************************************************************
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxshau(i,1) * tempvm
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.030 * vmax * max (0.50, uptotu(i))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstarud = min (0.9 * cishaud(i), o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapu(i) * 4.59e-06 * alpha3 * (cishaud(i) - gamstarud) / 
     >       (cishaud(i) + 2. * gamstarud)
c
        jc = vmax * (cishaud(i) - gamstarud) / 
     >       (cishaud(i) + kc * (1. + o2conc / ko))

c solve for intermediate photosynth. rate

        duma = theta3
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c sucrose synthesis limited rate of photosynth

        js = 3.0 * vmax / 8.2 * (1.0 - gamstarud / cishaud(i)) +
     >       jp * gamstarud / cishaud(i)
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15


        agshaud(i) = min(dumq/duma, dumc/dumq) * uptotu(i)

c net photosynthesis (mol/m**2/s)
c
        anshaud(i) = agshaud(i) - rdark
c
c ******************************************************************
c calculate photosynthesis rates for upper sunny leaves - conifers
c ******************************************************************
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxsunu(i,2) * tempvm
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.038 * vmax * max (0.50, uptotu(i))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstaruc = min(0.9 * cisunuc(i), o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpu(i) * 4.59e-06 * alpha3 * (cisunuc(i) - gamstaruc) / 
     >       (cisunuc(i) + 2. * gamstaruc)
c
        jc = vmax * (cisunuc(i) - gamstaruc) / 
     >       (cisunuc(i) + kc * (1. + o2conc / ko))

c solve for intermediate photosynth. rate

        duma = theta3
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.

        js = 3.0 * vmax / 8.2 * (1.0 - gamstaruc / cisunuc(i)) +
     >       jp * gamstaruc / cisunuc(i)
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15


        agsunuc(i) = min(dumq/duma, dumc/dumq) * uptotu(i)

c net photosynthesis (mol/m**2/s)
c
        ansunuc(i) = agsunuc(i) - rdark
c
c ******************************************************************
c calculate photosynthesis rates for upper shaded leaves - conifer
c ******************************************************************
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxshau(i,2) * tempvm
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.038 * vmax * max (0.50, uptotu(i))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstaruc = min(0.9 * cishauc(i), o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapu(i) * 4.59e-06 * alpha3 * (cishauc(i) - gamstaruc) / 
     >       (cishauc(i) + 2. * gamstaruc)
c
        jc = vmax * (cishauc(i) - gamstaruc) / 
     >       (cishauc(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

        duma = theta3
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.

        js = 3.0 * vmax / 8.2 * (1.0 - gamstaruc / cishauc(i)) +
     >       jp * gamstaruc / cishauc(i)
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        agshauc(i) = min(dumq/duma, dumc/dumq) * uptotu(i)

c net photosynthesis (mol/m**2/s)
c
        anshauc(i) = agshauc(i) - rdark
c
c **********************************************************
c calculate photosynthesis rates for lower sunny leaves - c3
c **********************************************************
c
        tleaf = tl(i) - 273.16
c
c calculate the parameter values which are a fnction of temperature
c
        tempvm = func1(tleaf) *
     >           exp(3500.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
        tau    = tau15 * 
     >           exp(-5000. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        kc     = kc15 *
     >           exp(6000.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
        ko     = ko15 *
     >           exp(1400.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxsunl(i,1) * tempvm
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.030 * vmax * max (0.50, uptotl(i))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstarl3 = min(0.9 * cisunl3(i), o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpl(i) * 4.59e-06 * alpha3 * (cisunl3(i) - gamstarl3) / 
     >       (cisunl3(i) + 2. * gamstarl3)
c
        jc = vmax * (cisunl3(i) - gamstarl3) / 
     >       (cisunl3(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

        duma = theta3
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.
c
        js = 3.0 * vmax / 8.2 * (1.0 - gamstarl3 / cisunl3(i)) +
     >       jp * gamstarl3 / cisunl3(i)
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        agsunl3(i) = min(dumq/duma, dumc/dumq) * uptotl(i)

c net photosynthesis (mol/m**2/s)
c
        ansunl3(i) = agsunl3(i) - rdark
c
c ***********************************************************
c calculate photosynthesis rates for lower shaded leaves - c3
c ***********************************************************
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxshal(i,1) * tempvm
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.030 * vmax * max (0.50, uptotl(i))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstarl3 = min(0.9 * cishal3(i), o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapl(i) * 4.59e-06  * alpha3 * (cishal3(i) - gamstarl3) /
     >       (cishal3(i) + 2. * gamstarl3)
c
        jc = vmax * (cishal3(i) - gamstarl3) / 
     >       (cishal3(i) + kc * (1. + o2conc / ko))
c 
c solve for intermediate photosynth. rate

        duma = theta3
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.

        js = 3.0 * vmax / 8.2 * (1.0 - gamstarl3 / cishal3(i)) +
     >       jp * gamstarl3 / cishal3(i)
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        agshal3(i) = min(dumq/duma, dumc/dumq) * uptotl(i)

c net photosynthesis (mol/m**2/s)
c
        anshal3(i) = agshal3(i) - rdark
c
c **********************************************************
c calculate photosynthesis rates for lower sunny leaves - c4
c **********************************************************
c
c calculate the parameter values which are a fnction of temperature
c
        tempvm = func2(tleaf) *
     >           exp(3500. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        vmax = vmaxsunl(i,2) * tempvm
c
c co2 efficiency (mol/m**2/s)
c (initial slope of photosynthetic co2 response curve)
c
        kco2 = 18.0e+03 * vmax
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.040 * vmax * max (0.50, uptotl(i))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpl(i) * 4.59e-06 * alpha4
c
        jc = vmax
c
c solve for intermediate photosynth. rate

        duma = theta4
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c co2 limited rate of photosynth.
c
        ji = kco2 * cisunl4(i)
c
        duma = beta4
        dumb = jp + ji
        dumc = jp * ji

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        agsunl4(i) = min(dumq/duma, dumc/dumq) * uptotl(i)

c net photosynthesis (mol/m**2/s)
c
        ansunl4(i) = agsunl4(i) - rdark
c
c ***********************************************************
c calculate photosynthesis rates for lower shaded leaves - c4
c ***********************************************************
c
c calculate the realized vmax (mol/m**2/s)
c
        vmax = vmaxshal(i,2) * tempvm
c
c co2 efficiency (mol/m**2/s)
c (initial slope of photosynthetic co2 response curve)
c
        kco2 = 18.0e+03 * vmax
c
c dark respiration (mol/m**2/s)
c
        rdark = 0.040 * vmax * max (0.50, uptotl(i))
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapl(i) * 4.59e-06 * alpha4
c
        jc = vmax
c
c solve for intermediate photosynth. rate

        duma = theta4
        dumb = je + jc
        dumc = je * jc

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        jp = min(dumq/duma, dumc/dumq)

c co2 limited rate of photosynth.

        ji = kco2 * cishal4(i)
c
        duma = beta4
        dumb = jp + ji
        dumc = jp * ji

        dume = max(dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

        agshal4(i) = min(dumq/duma, dumc/dumq) * uptotl(i)

c net photosynthesis (mol/m**2/s)
c
        anshal4(i) = agshal4(i) - rdark
c
c **********************************************************************
c average photosynthesis over sunlit and shaded fractions
c **********************************************************************
c
        agud(i) = fstoma(i,2) * (1. - fwetu(i)) *
     >           (sunfu(i) * agsunud(i) + (1.0 - sunfu(i)) * agshaud(i))
c
        aguc(i) = fstoma(i,2) * (1. - fwetu(i)) *
     >           (sunfu(i) * agsunuc(i) + (1.0 - sunfu(i)) * agshauc(i))
c
        agl3(i) = fstoma(i,1) * (1. - fwetl(i)) *
     >           (sunfl(i) * agsunl3(i) + (1.0 - sunfl(i)) * agshal3(i))
c
        agl4(i) = fstoma(i,1) * (1. - fwetl(i)) *
     >           (sunfl(i) * agsunl4(i) + (1.0 - sunfl(i)) * agshal4(i))
c
        anud(i) = fstoma(i,2) * (1. - fwetu(i)) *
     >           (sunfu(i) * ansunud(i) + (1.0 - sunfu(i)) * anshaud(i))
c
        anuc(i) = fstoma(i,2) * (1. - fwetu(i)) *
     >           (sunfu(i) * ansunuc(i) + (1.0 - sunfu(i)) * anshauc(i))
c
        anl3(i) = fstoma(i,1) * (1. - fwetl(i)) *
     >           (sunfl(i) * ansunl3(i) + (1.0 - sunfl(i)) * anshal3(i))
c
        anl4(i) = fstoma(i,1) * (1. - fwetl(i)) *
     >           (sunfl(i) * ansunl4(i) + (1.0 - sunfl(i)) * anshal4(i))
c
c **********************************************************************
c calculate new values of cs
c **********************************************************************
c
        cssunud(i) = 0.5 * cssunud(i) +
     >               0.5 * (co2conc - ansunud(i) / gbco2u)
        csshaud(i) = 0.5 * csshaud(i) +
     >               0.5 * (co2conc - anshaud(i) / gbco2u)
        cssunuc(i) = 0.5 * cssunuc(i) +
     >               0.5 * (co2conc - ansunuc(i) / gbco2u)
        csshauc(i) = 0.5 * csshauc(i) +
     >               0.5 * (co2conc - anshauc(i) / gbco2u)
        cssunl3(i) = 0.5 * cssunl3(i) +
     >               0.5 * (co2conc - ansunl3(i) / gbco2l)
        csshal3(i) = 0.5 * csshal3(i) +
     >               0.5 * (co2conc - anshal3(i) / gbco2l)
        cssunl4(i) = 0.5 * cssunl4(i) +
     >               0.5 * (co2conc - ansunl4(i) / gbco2l)
        csshal4(i) = 0.5 * csshal4(i) +
     >               0.5 * (co2conc - anshal4(i) / gbco2l)
c
c define minimum cs (mol/mol)
c
        csminud = gamstarud
        csminuc = gamstaruc
        csminl3 = gamstarl3
        csminl4 = 30.0e-06
c
        cssunud(i) = max (csminud, cssunud(i))
        csshaud(i) = max (csminud, csshaud(i))
        cssunuc(i) = max (csminuc, cssunuc(i))
        csshauc(i) = max (csminuc, csshauc(i))
        cssunl3(i) = max (csminl3, cssunl3(i))
        csshal3(i) = max (csminl3, csshal3(i))
        cssunl4(i) = max (csminl4, cssunl4(i))
        csshal4(i) = max (csminl4, csshal4(i))
c
c **********************************************************************
c calculate stomatal conductance
c **********************************************************************
c
c calculate the relative humidity in the canopy air space
c
        esat12 = esat (t12(i))
        qsat12 = qsat (esat12, psurf(i))
        rh12   = max(0.20, q12(i) / qsat12)
c
        esat34 = esat (t34(i))
        qsat34 = qsat (esat34, psurf(i))
        rh34   = max(0.20, q34(i) / qsat34)
c
c determine the stomatal conductance (mol-h2o/m**2/s)
c
        gssunud(i) = 0.5 * gssunud(i) + 0.5 *
     >               (9. * ansunud(i) * rh12 / cssunud(i) + 0.010)
c
        gsshaud(i) = 0.5 * gsshaud(i) + 0.5 *
     >               (9. * anshaud(i) * rh12 / csshaud(i) + 0.010)
c
        gssunuc(i) = 0.5 * gssunuc(i) + 0.5 *
     >               (6. * ansunuc(i) * rh12 / cssunuc(i) + 0.01)
c
        gsshauc(i) = 0.5 * gsshauc(i) + 0.5 *
     >               (6. * anshauc(i) * rh12 / csshauc(i) + 0.01)
c
        gssunl3(i) = 0.5 * gssunl3(i) + 0.5 *
     >               (9. * ansunl3(i) * rh34 / cssunl3(i) + 0.010)
c
        gsshal3(i) = 0.5 * gsshal3(i) + 0.5 *
     >               (9. * anshal3(i) * rh34 / csshal3(i) + 0.010)
c
        gssunl4(i) = 0.5 * gssunl4(i) + 0.5 *
     >               (3. * ansunl4(i) * rh34 / cssunl4(i) + 0.040)
c
        gsshal4(i) = 0.5 * gsshal4(i) + 0.5 *
     >               (3. * anshal4(i) * rh34 / csshal4(i) + 0.040)
c
c set minimum stomatal conductance (mol/m**2/s)
c
        gssunud(i) = max(gssunud(i), 0.010) 
        gsshaud(i) = max(gsshaud(i), 0.010)
        gssunuc(i) = max(gssunuc(i), 0.01) 
        gsshauc(i) = max(gsshauc(i), 0.01)
        gssunl3(i) = max(gssunl3(i), 0.010) 
        gsshal3(i) = max(gsshal3(i), 0.010)
        gssunl4(i) = max(gssunl4(i), 0.040) 
        gsshal4(i) = max(gsshal4(i), 0.040)
c
c set maximum stomatal conductance (mol/m**2/s)
c
        gsmax = 2.0
c
        gssunud(i) = min(gssunud(i), gsmax) 
        gsshaud(i) = min(gsshaud(i), gsmax)
        gssunuc(i) = min(gssunuc(i), gsmax) 
        gsshauc(i) = min(gsshauc(i), gsmax)
        gssunl3(i) = min(gssunl3(i), gsmax) 
        gsshal3(i) = min(gsshal3(i), gsmax)
        gssunl4(i) = min(gssunl4(i), gsmax) 
        gsshal4(i) = min(gsshal4(i), gsmax)
c
c average stomatal conductance over sunlit and shaded fractions
c
        gsud(i) = fstoma(i,2) * (1.0 - fwetu(i)) *
     >           (sunfu(i) * gssunud(i) + (1.0 - sunfu(i)) * gsshaud(i))
c
        gsuc(i) = fstoma(i,2) * (1.0 - fwetu(i)) *
     >           (sunfu(i) * gssunuc(i) + (1.0 - sunfu(i)) * gsshauc(i))
c
        gsl3(i) = fstoma(i,1) * (1.0 - fwetl(i)) *
     >           (sunfl(i) * gssunl3(i) + (1.0 - sunfl(i)) * gsshal3(i))
c
        gsl4(i) = fstoma(i,1) * (1.0 - fwetl(i)) *
     >           (sunfl(i) * gssunl4(i) + (1.0 - sunfl(i)) * gsshal4(i))
c
c **********************************************************************
c calculate new values of ci
c **********************************************************************
c
        cisunud(i) = 0.5 * cisunud(i) +
     >               0.5 * (cssunud(i) - 1.6 * ansunud(i) / gssunud(i))
        cishaud(i) = 0.5 * cishaud(i) +
     >               0.5 * (csshaud(i) - 1.6 * anshaud(i) / gsshaud(i))
        cisunuc(i) = 0.5 * cisunuc(i) +
     >               0.5 * (cssunuc(i) - 1.6 * ansunuc(i) / gssunuc(i))
        cishauc(i) = 0.5 * cishauc(i) +
     >               0.5 * (csshauc(i) - 1.6 * anshauc(i) / gsshauc(i))
        cisunl3(i) = 0.5 * cisunl3(i) +
     >               0.5 * (cssunl3(i) - 1.6 * ansunl3(i) / gssunl3(i))
        cishal3(i) = 0.5 * cishal3(i) +
     >               0.5 * (csshal3(i) - 1.6 * anshal3(i) / gsshal3(i))
        cisunl4(i) = 0.5 * cisunl4(i) +
     >               0.5 * (cssunl4(i) - 1.6 * ansunl4(i) / gssunl4(i))
        cishal4(i) = 0.5 * cishal4(i) +
     >               0.5 * (csshal4(i) - 1.6 * anshal4(i) / gsshal4(i))
c
c define minimum ci (mol/mol)
c
        ciminud = max (gamstarud * 1.5, 50.e-06)
        ciminuc = max (gamstaruc * 1.5, 50.e-06)
        ciminl3 = max (gamstarl3 * 1.5, 50.e-06)
        ciminl4 = 30.e-06
c
        cisunud(i) = max (ciminud, cisunud(i))
        cishaud(i) = max (ciminud, cishaud(i))
        cisunuc(i) = max (ciminuc, cisunuc(i))
        cishauc(i) = max (ciminuc, cishauc(i))
        cisunl3(i) = max (ciminl3, cisunl3(i))
        cishal3(i) = max (ciminl3, cishal3(i))
        cisunl4(i) = max (ciminl4, cisunl4(i))
        cishal4(i) = max (ciminl4, cishal4(i))
c
 100  continue
c
c        If last turvap iteration, update optimum Rubisco capacities
c
      if (iter.eq.niter) call optophys (dtstom)
c     if (iter.eq.1) call optophys (dtstom)

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine optophys (dtime)

c Time steps IBIS-model running-mean equation for optimum Rubisco
c capacities for upper/lower stories, sunlit/shaded fractions
c (vmax[sun,sha](*,2) in comveg, corresponding to vmax*xy in comvegxy)

c     All arguments are supplied:
c     dtime = time step for this call, seconds 

#include <compar.h>
#include <comgrd.h>
#include <comveg.h>
c
c declaration of variables
c
      real kc, ko, kco2, kc15, ko15,
     >     je, jc, js, jp, ji
c
      dimension an(3)
c
c *******************************************************  
c parameters that might change
c *******************************************************  
c
c atmospheric oxygen concentration (mol/mol)
c
      data o2conc  /0.209000/
c
c *******************************************************  
c physiological parameters
c *******************************************************  
c
c quantum efficiency for c3 and c4 plants (dimensionless)
c
      data alpha3 /0.065/
      data alpha4 /0.050/
c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
      data tau15 /4500.0/
c
c o2/co2 kinetic parameters (mol/mol)
c
      data kc15 /1.5e-04/
      data ko15 /2.5e-01/
c
c collatz coupling coefficients (dimensionless)
c
      data theta3 /0.90/
      data beta3  /0.90/
c
      data theta4 /0.90/
      data beta4  /0.90/
c
c vmax adjustment period or 'memory' (days)
c
      vmem = 30.0 * 86400.0
c
c vmax perturbation during any timestep (mol/m**2/s)
c
      dvopt = 5.0e-06
c
c timestep factor
c
      dtfac = 3. * (50.e-6 / vmem) / (dvopt / dtime)
c
c *****************
c do for all points
c *****************
c
      do 100 i = 1, npoi

c ********************************
c upper sunny leaves - broadleaf
c ********************************
c
c define leaf temperature for the leaf
c
        tleaf = tu(i) - 273.16
c
c calculate the parameter values which depend on temperature
c
        tempvm = func1(tleaf) *
     >           exp(3500. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        tau = tau15 * 
     >        exp(-5000. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        kc  = kc15 *
     >        exp(6000.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
        ko  = ko15 *
     >        exp(1400.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cisunud(i),
     >                o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpu(i) * 4.59e-06 * alpha3 * (cisunud(i) - gamstar) / 
     >       (cisunud(i) + 2. * gamstar)
c
c trial values of vmax
c
        do 110 j = -1, 1
c
          vmax = max (0., vmaxsunu(i,1) + j * dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax * (cisunud(i) - gamstar) / 
     >         (cisunud(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

          duma = theta3
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.

          js = 3.0 * avmax / 8.2 * (1.0 - gamstar / cisunud(i)) +
     >         jp * gamstar / cisunud(i)
c
          duma = beta3
          dumb = jp + js
          dumc = jp * js

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.030 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 110    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxsunu(i,1) = vmaxsunu(i,1) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxsunu(i,1) = vmaxsunu(i,1) + dtfac * (anmax - an(2))
        endif
c
c *********************************
c upper shaded leaves - broadleaf
c *********************************
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cishaud(i), gamstar)
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapu(i) * 4.59e-06 * alpha3 * (cishaud(i) - gamstar) / 
     >       (cishaud(i) + 2. * gamstar)
c
c trial values of vmax
c
        do 120 j = -1, 1
c
          vmax = max (0., vmaxshau(i,1) + j * dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax * (cishaud(i) - gamstar) / 
     >         (cishaud(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

          duma = theta3
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.
c
          js = 3.0 * avmax / 8.2 * (1.0 - gamstar / cishaud(i)) +
     >         jp * gamstar / cishaud(i)
c
          duma = beta3
          dumb = jp + js
          dumc = jp * js

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.030 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 120    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxshau(i,1) = vmaxshau(i,1) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxshau(i,1) = vmaxshau(i,1) + dtfac * (anmax - an(2))
        endif
c
c ******************************
c upper sunny leaves - conifer
c ******************************
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cisunuc(i),
     >                o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpu(i) * 4.59e-06 * alpha3 * (cisunuc(i) - gamstar) / 
     >       (cisunuc(i) + 2. * gamstar)
c
c trial values of vmax
c
        do 130 j = -1, 1
c
          vmax = max (0., vmaxsunu(i,2) + j * dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax * (cisunuc(i) - gamstar) / 
     >         (cisunuc(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

          duma = theta3
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.
c
          js = 3.0 * avmax / 8.2 * (1.0 - gamstar / cisunuc(i)) +
     >         jp * gamstar / cisunuc(i)
c
          duma = beta3
          dumb = jp + js
          dumc = jp * js

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.038 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 130    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxsunu(i,2) = vmaxsunu(i,2) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxsunu(i,2) = vmaxsunu(i,2) + dtfac * (anmax - an(2))
        endif
c
c *****************************
c upper shaded leaves - conifer
c *****************************
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cishauc(i), gamstar)
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapu(i) * 4.59e-06 * alpha3 * (cishauc(i) - gamstar) / 
     >       (cishauc(i) + 2. * gamstar)
c
c trial values of vmax
c
        do 140 j = -1, 1
c
          vmax = max (0., vmaxshau(i,2) + j * dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax * (cishauc(i) - gamstar) / 
     >         (cishauc(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

          duma = theta3
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)

c sucrose synth. limited rate of photosynth.
c
          js = 3.0 * avmax / 8.2 * (1.0 - gamstar / cishauc(i)) +
     >         jp * gamstar / cishauc(i)
c
          duma = beta3
          dumb = jp + js
          dumc = jp * js

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.038 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 140    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxshau(i,2) = vmaxshau(i,2) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxshau(i,2) = vmaxshau(i,2) + dtfac * (anmax - an(2))
        endif
c
c ***********************
c lower sunny leaves - c3
c ***********************
c
c define leaf temperature for the leaf
c
        tleaf = tl(i) - 273.16
c
c calculate the parameter values which depend on temperature
c
        tempvm = func1(tleaf) *
     >           exp(3500. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        tau = tau15 * 
     >        exp(-5000. * (3.47e-03 - 1./(tleaf + 273.16)))
c
        kc  = kc15 *
     >        exp(6000.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
        ko  = ko15 *
     >        exp(1400.  * (3.47e-03 - 1./(tleaf + 273.16)))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cisunl3(i),
     >                o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpl(i) * 4.59e-06 * alpha3 * (cisunl3(i) - gamstar) / 
     >       (cisunl3(i) + 2. * gamstar)
c
c trial values of vmax
c
        do 150 j = -1, 1
c
          vmax = max (0., vmaxsunl(i,1) + j*dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax * (cisunl3(i) - gamstar) / 
     >         (cisunl3(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

          duma = theta3
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)
c
c sucrose synth. limited rate of photosynth.

          js = 3.0 * avmax / 8.2 * (1.0 - gamstar / cisunl3(i)) +
     >         jp * gamstar / cisunl3(i)
c
          duma = beta3
          dumb = jp + js
          dumc = jp * js

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.030 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 150    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxsunl(i,1) = vmaxsunl(i,1) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxsunl(i,1) = vmaxsunl(i,1) + dtfac * (anmax - an(2))
        endif
c
c ************************
c lower shaded leaves - c3
c ************************
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cishal3(i), gamstar)
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapl(i) * 4.59e-06 * alpha3 * (cishal3(i) - gamstar) / 
     >       (cishal3(i) + 2. * gamstar)
c
c trial values of vmax
c
        do 160 j = -1, 1
c
          vmax = max (0., vmaxshal(i,1) + j*dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax * (cishal3(i) - gamstar) / 
     >         (cishal3(i) + kc * (1. + o2conc / ko))
c
c solve for intermediate photosynth. rate

          duma = theta3
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)
c
c sucrose synth. limited rate of photosynth.
c
          js = 3.0 * avmax / 8.2 * (1.0 - gamstar / cishal3(i)) +
     >         jp * gamstar / cishal3(i)
c
          duma = beta3
          dumb = jp + js
          dumc = jp * js

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.030 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 160    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxshal(i,1) = vmaxshal(i,1) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxshal(i,1) = vmaxshal(i,1) + dtfac * (anmax - an(2))
        endif
c
c ***********************
c lower sunny leaves - c4
c ***********************
c
c calculate the parameter values which depend on temperature
c
        tempvm = func2(tleaf) *
     >           exp(3500. * (3.47e-03 - 1./(tleaf + 273.16)))
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cisunl4(i),
     >                o2conc / (2. * tau))
c
c gross photosynthesis (mol/m**2/s)
c
        je = sunpl(i) * 4.59e-06 * alpha4
c
c trial values of vmax
c
        do 170 j = -1, 1
c
          vmax = max (0., vmaxsunl(i,2) + j*dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax
c
c solve for intermediate photosynth. rate

          duma = theta4
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)
c
c co2 efficiency (mol/m**2/s)
c (initial slope of photosynthetic co2 response curve)
c
          kco2 = 18.e+03 * avmax
c
          ji = kco2 * cisunl4(i)
c
          duma = beta4
          dumb = jp + ji
          dumc = jp * ji

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.040 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 170    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxsunl(i,2) = vmaxsunl(i,2) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxsunl(i,2) = vmaxsunl(i,2) + dtfac * (anmax - an(2))
        endif
c
c ************************
c lower shaded leaves - c4
c ************************
c
c carbon dioxide concentrations (mol/mol)
c
        gamstar = min(0.9 * cishal4(i), gamstar)
c
c gross photosynthesis (mol/m**2/s)
c
        je = shapl(i) * 4.59e-06 * alpha4
c
c trial values of vmax
c
        do 180 j = -1, 1
c
          vmax = max (0., vmaxshal(i,2) + j*dvopt)
c
          avmax = vmax * tempvm
c
          jc = avmax
c
c solve for intermediate photosynth. rate

          duma = theta4
          dumb = je + jc
          dumc = je * jc

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          jp = min(dumq/duma, dumc/dumq)
c
c co2 efficiency (mol/m**2/s)
c (initial slope of photosynthetic co2 response curve)
c
          kco2 = 18.e+03 * avmax
c
          ji = kco2 * cishal4(i)
c
          duma = beta4
          dumb = jp + ji
          dumc = jp * ji

          dume = max(dumb**2 - 4. * duma * dumc, 0.)
          dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

          ag = min(dumq/duma, dumc/dumq)
c
c dark respiration (mol/m**2/s)
c
          rdark = 0.040 * avmax
c
c net photosynthesis (mol/m**2/s)
c
          an(j+2) = ag - rdark
c
 180    continue
c
c update running mean of optimum vmax
c
        anmax = max (an(1), an(2), an(3))
c
        if (an(1).eq.anmax) then
          vmaxshal(i,2) = vmaxshal(i,2) - dtfac * (anmax - an(2))
        else if (an(3).eq.anmax) then
          vmaxshal(i,2) = vmaxshal(i,2) + dtfac * (anmax - an(2))
        endif
c
c constrain the minimum and max values of vmax
c
        vmaxsunu(i,1) = min(200.e-6, max(5.0e-06, vmaxsunu(i,1)))
        vmaxsunu(i,2) = min(200.e-6, max(5.0e-06, vmaxsunu(i,2)))
        vmaxsunl(i,1) = min(200.e-6, max(5.0e-06, vmaxsunl(i,1)))
        vmaxsunl(i,2) = min(200.e-6, max(5.0e-06, vmaxsunl(i,2)))
        vmaxshau(i,1) = min(200.e-6, max(5.0e-06, vmaxshau(i,1)))
        vmaxshau(i,2) = min(200.e-6, max(5.0e-06, vmaxshau(i,2)))
        vmaxshal(i,1) = min(200.e-6, max(5.0e-06, vmaxshal(i,1)))
        vmaxshal(i,2) = min(200.e-6, max(5.0e-06, vmaxshal(i,2)))
c
 100  continue

      return
      end
c
c ---------------------------------------------------------------------
c
      subroutine drystress (stressu, stressl, stresstu, stresstl)

#include <compar.h>
#include <comgrd.h>
#include <comsoi.h>
#include <comveg.h>

      dimension stressu(nvec,nsoilay), stressl(nvec,nsoilay),
     *          stresstu(nvec),        stresstl(nvec)

      real stressfac, znorm, zw, zwilt, zup
c
c stressfac determines the 'strength' of the soil moisture stress
c on gross photosynthesis rates
c
      stressfac = -8.0
      znorm = (1.0 - exp(stressfac * 1.0))
c
      call zero (stresstl, nvec)
      call zero (stresstu, nvec)
c
      do 100 k = 1, nsoilay
        do 110 i = 1, npoi
c
          zw    = max (0.0,
     >                      min (1.0,
     >                                (wsoi(i,k) - swilt(i,k)) /
     >                                (1.0 - swilt(i,k))
     >                          )
     >                )
c
          zwilt = max (0.0,
     >                      min (1.0,
     >                           (1.0 - exp(stressfac * zw)) / znorm
     >                          )
     >                )
c
          zup   = (1.0 - min (wisoi(i,k), 1.0)) * zwilt
c
          stressl(i,k) = zup * froot(i,k,1)
          stressu(i,k) = zup * froot(i,k,2)
c
          stresstl(i) = stresstl(i) + stressl(i,k)
          stresstu(i) = stresstu(i) + stressu(i,k)
c
 110    continue
 100  continue
c
      return
      end
c
c----------------------------------------------------------------------
c
      real function func1(t)
c
      func1 = 1. / (
     >              (1.0 + exp(0.40 * (5.0 - t))) *
     >              (1.0 + exp(0.40 * (t - 40.0)))
     >             )
c
      return
      end
c
c----------------------------------------------------------------------
c
      real function func2(t)
c
      func2 = 1. / (
     >              (1.0 + exp(0.40 * (10.0 - t))) *
     >              (1.0 + exp(0.40 * (t - 45.0)))
     >             )
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vegdat_wisc (showmaps, jyear, jmonth, nsrest_ibis,
     *                        lonreslsx, latreslsx, typeresin)
 
c        Sets current state of vegetation on lsx 2-D grid
c        using `Wisconsin' vegetation categories.
c
c        showmaps  = .true. to print 2-D maps of prescribed flds (suppl)
c        jyear     = current calendar year # (eg, 2001) (supplied)
c        jmonth    = current month # (1-12) (supplied)
c        nsrest_ibis = -1/[0,1] flag for coldstart/restart (suppl)
c        lonreslsx = expected long. resol. of IBIS restart file (suppl)
c        latreslsx = expected lat.  resol. of IBIS restart file (suppl)
c        typeresin = CRAY/IEEE flag,machine type of restart files(modif)
c
#include <compar.h>
#include <combudl.h>
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
 
      character*(*) typeresin
 
      logical showmaps
 
c     Solar properties for stems:
      data rho_ste_vis, rho_ste_nir, tau_ste_vis, tau_ste_nir
     *     /.16, .39,.001,.001/
      save rho_ste_vis, rho_ste_nir, tau_ste_vis, tau_ste_nir
 
      logical first
      save first
      data first /.true./
 
 
c====================
      if (first) then
c====================
 
c          If first call and non-interactive veg, read in WISCVEG
c          data file into ivegxy(or set to uniform value if namelist
c          param wiscveg ne 0)
 
c------------------------------
        if (dynaibis.eq.0) then
c------------------------------
 
          if (wiscveg.eq.0) then
 
            call findiu (iuwisc, loc)
#ifdef gcm
            write(ioterm2,40) locname(loc)(1:lenchr(locname(loc)))
   40       format(/' Reading Wisconsin vegetation data file ',a)
#endif
            call openf (iuwisc, locname(loc), ' ', 0, 0)
            call readdata (iuwisc, workl, nlon, nlat, 'WISC', 2, 1)
            call closef(iuwisc)
            do 42 jj=1,nlat
              do 44 ji=1,nlon
                ivegxy(ji,jj) = workl(ji,jj) + .001
                if (lmaskxy(ji,jj).ne.1) then
                  ivegxy(ji,jj) = 0
                endif
   44         continue
   42       continue
 
          else
 
c           Namelist wiscveg entered for global uniform single type:
            write(ioterm2,50) wiscveg
   50       format(/' *** Using single Wisconsin vegetation category',
     *              i5,' for all land points')
            do 52 jj=1,nlat
              do 54 ji=1,nlon
                if (lmaskxy(ji,jj).eq.1) then
                  ivegxy(ji,jj) = wiscveg
                else
                  ivegxy(ji,jj) = 0
                endif
   54         continue
   52       continue
 
          endif
 
c            Print out vegetation-type map if requested
 
          if (showmaps)
     *      call mapout (ioterm2, ivegxy, nlon, nlat, alon, alat, 0,
     *                   'Wisconsin Vegetation Categories')
 
c-----------
        else
c-----------
 
c            If dynamic IBIS veg requested, read IBIS veg restart file
 
          call ibis_read (nsrest_ibis, nlon, nlat, typeresin)
 
c------------
        endif
c------------
 
        first = .false.
 
c==========
      endif
c==========
 
c        Call Wisconsin veg subr to return current day's veg fields.
c        lsx arrays are in common comvegxy. All veg arrays *are* reset
c        by wisc_veg every call even though some are constant in time.
 
      call wisc_veg (hsod, nsrest_ibis, jyear, jmonth)
 
c        Diagnostic printout of leaf quantities for the lsx output point
 
      if (.true. .and. lonbud.ne.-999 .and. latbud.ne.-999) then
c     if (.false.) then
        iuz = 57
        iz = lonbud
        jz = latbud
        zlon = ((iz-.5)/nlon)*360.
        if (zlon.gt.180.) zlon = zlon - 360.
        zlat = -90. + ((jz-.5)/nlat)*180.
        write(iuz,*)
        write(iuz,*) 'vegdat_wisc: zlon=',zlon,'  zlat=',zlat
        write(iuz,8888) 'flxy,fuxy   ', flxy(iz,jz), fuxy(iz,jz)
        write(iuz,8888) 'ztopxy      ',(ztopxy(iz,jz,m), m=1,2)
        write(iuz,8888) 'zbotxy      ',(zbotxy(iz,jz,m), m=1,2)
        write(iuz,8888) 'laixy       ',(laixy(iz,jz,m), m=1,2)
        write(iuz,8888) 'saixy       ',(saixy(iz,jz,m), m=1,2)
        write(iuz,8888) 'fstomaxy    ',(fstomaxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'dleafxy     ',(dleafxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'dstemxy     ',(dstemxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'frootxy(1)  ',(frootxy(iz,jz,1,m),m=1,2)
        write(iuz,8888) 'frootxy(2)  ',(frootxy(iz,jz,2,m),m=1,2)
        write(iuz,8888) 'frootxy(3)  ',(frootxy(iz,jz,3,m),m=1,2)
        write(iuz,8888) 'frootxy(4)  ',(frootxy(iz,jz,4,m),m=1,2)
        write(iuz,8888) 'frootxy(5)  ',(frootxy(iz,jz,5,m),m=1,2)
        write(iuz,8888) 'frootxy(6)  ',(frootxy(iz,jz,6,m),m=1,2)
        write(iuz,8888) 'confracxy   ', confracxy(iz,jz)
        write(iuz,8888) 'c4fracxy    ', c4fracxy(iz,jz)
        write(iuz,8888) 'orievxy     ',(orievxy(iz,jz,m),m=1,2)
        write(iuz,8888) 'rhovegxy(1) ',(rhovegxy(iz,jz,1,m),m=1,2)
        write(iuz,8888) 'rhovegxy(2) ',(rhovegxy(iz,jz,2,m),m=1,2)
        write(iuz,8888) 'tauvegxy(1) ',(tauvegxy(iz,jz,1,m),m=1,2)
        write(iuz,8888) 'tauvegxy(2) ',(tauvegxy(iz,jz,2,m),m=1,2)
c       only for BATS stomatal resistance (phototype=0):
c       write(iuz,8888) 'rsmaxxy     ',(rsmaxxy(iz,jz,m),m=1,2)
c       write(iuz,8888) 'rsminxy     ',(rsminxy(iz,jz,m),m=1,2)
c       write(iuz,8888) 'parconxy    ',(parconxy(iz,jz,m),m=1,2)
c       write(iuz,8888) 'vpdconxy    ',(vpdconxy(iz,jz,m),m=1,2)
 8888   format(a,':  ',2f20.10)
        ier = flush(iuz)
      endif
 
c       Ensure no vegetation over ocean or icesheet,
c       and ensure no tiny amounts of lai or sai
 
 
c     zlaimin = 0.01
      zlaimin = 0.05   ! ibis2.0
      zsaimin = 0.001
 
      do 100 jj=1,nlat
        do 102 ji=1,nlon
          if (lmaskxy(ji,jj).ne.1) then
            flxy(ji,jj) = 0.
            fuxy(ji,jj) = 0.
          endif
 
          if (laixy(ji,jj,1).lt.zlaimin) laixy(ji,jj,1) = 0.
          if (saixy(ji,jj,1).lt.zsaimin) saixy(ji,jj,1) = 0.
          if (laixy(ji,jj,1).eq.0. .and. saixy(ji,jj,1).eq.0.)
     *      flxy(ji,jj) = 0.
 
          if (laixy(ji,jj,2).lt.zlaimin) laixy(ji,jj,2) = 0.
          if (saixy(ji,jj,2).lt.zsaimin) saixy(ji,jj,2) = 0.
          if (laixy(ji,jj,2).eq.0. .and. saixy(ji,jj,2).eq.0.)
     *      fuxy(ji,jj) = 0.
  102   continue
  100 continue
 
c        Set harmless values for non-veg points or empty stories
c        (wisc_veg returns values as zero if fractional cover or lai/sai
c        are zero, but lsx needs "safe" values for land/landice pts
c        for both stories even if no veg (with lai = sai = fractional
c        cover = 0), since radiative and aero calcs are done whether
c        canopies exist or not, to allow Cray vectorization.
c        (nb: wisc_veg has set frootxy ok for all points).
 
      do 150 m=1,2
        do 152 jj=1,nlat
          do 154 ji=1,nlon
 
            if ( (m.eq.1 .and. flxy(ji,jj).eq.0.) .or.
     *           (m.eq.2 .and. fuxy(ji,jj).eq.0.)      ) then
              ztopxy(ji,jj,m) = 0.
              zbotxy(ji,jj,m) = 0.
              dstemxy(ji,jj,m) = .10
 
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
 
c           If lai is zero, set leaf quantities
            if (laixy(ji,jj,m).eq.0.) then
              dleafxy(ji,jj,m) = .05
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
c
c        Constrain geometric canopy heights to avoid blowups in turlsx
c       (Redundant because these limits are set in dynaveg)
 
c     do 160 jj=1,nlat
c       do 162 ji=1,nlon
c         zbotxy(ji,jj,1) = max (zbotxy(ji,jj,1), .02)
c         ztopxy(ji,jj,1) = max (ztopxy(ji,jj,1), zbotxy(ji,jj,1)+.01)
c         zbotxy(ji,jj,2) = max (zbotxy(ji,jj,2), ztopxy(ji,jj,1)+.01)
c         ztopxy(ji,jj,2) = max (ztopxy(ji,jj,2), zbotxy(ji,jj,2)+.01)
c 162   continue
c 160 continue
 
c         Convert various wisc_veg quantities into lsx forms
c         (m=1 for lower story, m=2 for upper story)
 
      do 170 m=1,2
        do 172 jj=1,nlat
          do 174 ji=1,nlon
 
c              Double all lai's since lsx uses 2-sided values
c              *** Don't *** do for sai, since (i) lsx is ambiguous
c              whether this is necessary, and (ii) saixy are only set
c              by wisc_veg at the first timestep, and at end of each
c              year (interactive), so can't double it here once a day!)
 
            laixy(ji,jj,m) = 2.*laixy(ji,jj,m)
 
c              Lump leaves with stems for solar refl/trans weights
 
            ztot = max (laixy(ji,jj,m)+saixy(ji,jj,m), epsilon)
            wliv= laixy(ji,jj,m) / ztot
            wste = 1.-wliv
            rhovegxy(ji,jj,1,m) = wliv * rhovegxy(ji,jj,1,m)
     *                          + wste * rho_ste_vis
            rhovegxy(ji,jj,2,m) = wliv * rhovegxy(ji,jj,2,m)
     *                          + wste * rho_ste_nir
            tauvegxy(ji,jj,1,m) = wliv * tauvegxy(ji,jj,1,m)
     *                          + wste * tau_ste_vis
            tauvegxy(ji,jj,2,m) = wliv * tauvegxy(ji,jj,2,m)
     *                          + wste * tau_ste_nir
 
c              Leaf orientation
 
            zchi = orievxy(ji,jj,m)
            orievxy(ji,jj,m) = max (-zchi, 0.)
            oriehxy(ji,jj,m) = max ( zchi, 0.)
 
c              Convert from stomatal conductance to stomatal resistance
c              (and so switch max and min)
 
            zrsmax = rsmaxxy(ji,jj,m)
            rsmaxxy(ji,jj,m) = 1./rsminxy(ji,jj,m)
            rsminxy(ji,jj,m) = 1./zrsmax
 
  174     continue
  172   continue
  170 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wisc_veg (hsod, nsrest_ibis, jyear, jmonth)
 
c     If DYNAIBIS = 1, sets vegetation attributes from dynamic IBIS
c                      vegetation state.
c     If DYNAIBIS = 0,prescribes `Wisconsin' vegetation attributes
c     from 4-digit codes in WISCVEG data input file.
 
c     supplied:
c     hsod = soil layer depths (assumed uniform globally)
c     nsrest_ibis = -1/[0,1] flag for coldstart/restart (suppl)
c     jyear = current calendar year (eg, 2001)
c     jmonth = current month # (1-12)
 
c------------------
#include <compar.h>
#include <comgrd.h>
#include <comvegxy.h>
#include <commudxy.h>
c------------------
      dimension hsod(nsoilay)
 
c The following variables are used or set by wisc_veg.
c They are all in comvegxy, except uptot[l,u]xy in commudxy.
c Since some are set by dynaveg (called from wiscveg for first timestep
c in a run), they are not passed as arguments as for other veg methods.
c
c supplied:
c    *  ivegxy,          ! Wisconsin veg categories from data file
c    *  uptotlxy,        ! soil moisture stress factor, lower story(0-1)
c    *  uptotuxy,        ! soil moisture stress factor, upper story(0-1)
c    *  t15runxy,        ! 15-day exponential running mean 2-ma temp (K)
c    *  agdd0xy,         ! degree days above 0 C based on t15runxy,
c    *  agdd5xy,         ! degree days above 5 C based on t15runxy ,
c                        ! reset to 0 if t15runxy drops below threshold
c
c returned:
c    *  flxy,            ! fractional cover of lower canopy
c                          (if not set already from IBIS restart file)
c    *  fuxy,            ! fractional cover of upper canopy
c                          (if not set already from IBIS restart file)
c    *  zbotxy,          ! geometric height of canopy bottom (m)
c    *  ztopxy,          ! geometric height of canopy top (m)
c    *  laixy,           ! leaf area index, one-sided (m2/m2)
c    *  saixy,           ! stem area index, one-sided (m2/m2)
c    *  fstomaxy,        ! fraction of lai that is transpiring
c    *  dleafxy,         ! aero dimension for leaves (m)
c    *  dstemxy,         ! aero dimension for stems (m)
c    *  frootxy,         ! fraction of roots in each soil layer (m)
c    *  confracxy,       ! fraction of conifers (vs deciduous), upper
c    *  c4fracxy,        ! fraction of c4 (vs c3), lower
c
c    *  orievxy,         ! leaf orientation (-1=vert, 0=random, 1=horiz)
c
c    *  rhovegxy,        ! *lai* leaf reflectivity,   vis and nir
c    *  tauvegxy,        ! *lai* leaf transmissivity, vis and nir
c
c only for BATS stomatal resistance (phototype=0):
c    *  rsmaxxy,         ! stomatal conductance, maximum (m/s)
c    *  rsminxy,         ! stomatal conductance, minimum (m/s)
c    *  parconxy,        ! stomatal PAR constant (W/m2)
c    *  vpdconxy         ! stomatal VPD constant (N/m2)
 
 
c internal variables
c
      integer tz, ct, ph, cp, temp
      dimension tempu(nlon,nlat),  templ(nlon,nlat),
     *          dropu(nlon,nlat),  dropl(nlon,nlat)
 
c     for non-interactive vegetation:
      dimension plail3(nlon,nlat), plail4(nlon,nlat),
     *          plaiud(nlon,nlat), plaiue(nlon,nlat)
      save plail3, plail4, plaiud, plaiue
 
      logical firstroot, firstdyna, firstfile
      save firstroot, firstdyna, firstfile
      data firstroot, firstdyna, firstfile /.true., .true., .true./
 
 
c        If first call, set soil-layer root fractions
 
c------------------------
      if (firstroot) then
c------------------------
 
c using new root calculations from ibis 1.4 as taken from
c Jackson et al., 1996: Oecologia, 108, 389-411.
 
        beta1 = 0.95 !grass
        beta2 = 0.97 !tree
 
        zdepb = 0.
 
        do 50 k=1,nsoilay
 
          zdepa = zdepb
          zdepb = zdepa + hsod(k) * 100.
 
          do 52 j=1,nlat
            do 54 i=1,nlon
 
              frootxy(i,j,k,1) = (1. - beta1 ** zdepb) -
     *                           (1. - beta1 ** zdepa)
              frootxy(i,j,k,2) = (1. - beta2 ** zdepb) -
     *                           (1. - beta2 ** zdepa)
 
  54        continue
  52      continue
  50    continue
 
        firstroot = .false.
 
c----------
      endif
c----------
 
c Phenology factors: same for dynamic and non-interactive veg
c -----------------
 
      do 70 j=1,nlat
        do 72 i=1,nlon
 
c cold deciduous tree phenology
 
          threshold = 278.16
          if (t15runxy(i,j).lt.threshold) then
            agdd5xy(i,j) = 0.
          else
            agdd5xy(i,j) = agdd5xy(i,j) + t15runxy(i,j) - threshold
          endif
c         tempu(i,j) = min (1., agdd5xy(i,j) / 100.)
          tempu(i,j) = min (1., agdd5xy(i,j) / 50.)    ! ibis2.0
 
c cold deciduous grass phenology
 
          threshold = 273.16
          if (t15runxy(i,j).lt.threshold) then
            agdd0xy(i,j) = 0.
          else
            agdd0xy(i,j) = agdd0xy(i,j) + t15runxy(i,j) - threshold
          endif
          templ(i,j) = min (1., agdd0xy(i,j) / 50.)
 
c drought deciduous phenology
 
          dropl(i,j) = 1.0
          dropu(i,j) = 1.0
          if (uptotlxy(i,j).lt.0.10) dropl(i,j) = 0.0
          if (uptotuxy(i,j).lt.0.25) dropu(i,j) = 0.0
 
   72   continue
   70 continue
 
 
c****************************
      if (dynaibis.eq.1) then
c****************************
 
c=============
#ifdef ibisint
c=============
 
c For **dynamic (interactive)** vegetation, set [l,s]aixy,
c f[l,u]xy, z[bot,top]xy, [con,c4]fracxy, frac[u,l]xy.
 
c If first call this run, dynaveg (0,...) is called to initialize
c diagnostic vegetation arrays from prognostic variables
c obtained from restart file (or crudely initialized for
c cold start) by ibis_read. The diagnostic variables set
c in this way that are needed here are saixy, plai[u,l]xy
c (which is used below for laixy), ztopxy, zbotxy,
c and also flxy,fuxy if not set already from IBIS restart file.
c
        if (firstdyna) then
          call dynaveg (0, '(re)start', jyear)
          firstdyna = .false.
        endif
 
        do 100 j=1,nlat
          do 102 i=1,nlon
c
c lower canopy one-sided leaf area index (area-weighted)
c
            laixy(i,j,1) = (plailxy(i,j,1) + plailxy(i,j,2))
     >                     * dropl(i,j) * templ(i,j)
c
c lower canopy fractions
c
            fraclxy(i,j,1) = plailxy(i,j,1) * dropl(i,j) * templ(i,j)
     >                       / max(laixy(i,j,1), epsilon)
c
            fraclxy(i,j,2) = plailxy(i,j,2) * dropl(i,j) * templ(i,j)
     >                       / max(laixy(i,j,1), epsilon)
c
            c4fracxy(i,j) = fraclxy(i,j,2)
c
c upper canopy one-sided leaf area index (area-weighted)
c
            laixy(i,j,2) = plaiuxy(i,j,1)              +
     >                     plaiuxy(i,j,2) * dropu(i,j) +
     >                     plaiuxy(i,j,3)              +
     >                     plaiuxy(i,j,4)              +
     >                     plaiuxy(i,j,5) * tempu(i,j) +
     >                     plaiuxy(i,j,6)              +
     >                     plaiuxy(i,j,7) * tempu(i,j)
c
c upper canopy fractions
c
            fracuxy(i,j,1) = plaiuxy(i,j,1)              /
     >                       max(laixy(i,j,2), epsilon)
c
            fracuxy(i,j,2) = plaiuxy(i,j,2) * dropu(i,j) /
     >                       max(laixy(i,j,2), epsilon)
c
            fracuxy(i,j,3) = plaiuxy(i,j,3)              /
     >                       max(laixy(i,j,2), epsilon)
c
            fracuxy(i,j,4) = plaiuxy(i,j,4)              /
     >                       max(laixy(i,j,2), epsilon)
c
            fracuxy(i,j,5) = plaiuxy(i,j,5) * tempu(i,j) /
     >                       max(laixy(i,j,2), epsilon)
c
            fracuxy(i,j,6) = plaiuxy(i,j,6)              /
     >                       max(laixy(i,j,2), epsilon)
c
            fracuxy(i,j,7) = plaiuxy(i,j,7) * tempu(i,j) /
     >                       max(laixy(i,j,2), epsilon)
c
            confracxy(i,j) = fracuxy(i,j,4) + fracuxy(i,j,6)
 
c adjust the LAI by fractional veg cover, ie, convert from
c area-weighted to non-area-weighted as used by LSX.
 
c           laixy(i,j,1) = max( 0.01, laixy(i,j,1)/max(0.02,flxy(i,j)) )
c           laixy(i,j,2) = max( 0.01, laixy(i,j,2)/max(0.02,fuxy(i,j)) )
cibis2.0:
            laixy(i,j,1) = max( 0.05, laixy(i,j,1)/max(0.10,flxy(i,j)) )
            laixy(i,j,2) = max( 0.05, laixy(i,j,2)/max(0.10,fuxy(i,j)) )
 
  102     continue
  100   continue
 
c=====
#endif
c=====
 
c*********
      else
c*********
 
c For **non-interactive** vegetation, set [l,s]aixy, f[l,u]xy,
c z[bot,top]xy, [con,c4]fracxy. Only need to read the Wisconsin veg
c data file once, and set non-varying veg arrays (f[l,u]xy,
c z[bot,top]xy...also plai[l3,l4,ue,ud] which are local saved arrays).
 
c++++++++++++++++++++++++++
        if (firstfile) then
c++++++++++++++++++++++++++
 
          do 200 j=1,nlat
            do 202 i=1,nlon
 
c Four-digit vegetation codes (in ivegxy from file WISCVEG) mean:
c
c       cover type, first digit of code
c                          1 = forest
c                          2 = woodland
c                          3 = savanna
c                          4 = shrubland
c                          5 = grassland
c                          6 = tundra
c                          7 = desert
c
c       temperature zone, second digit of code
c                          1 = tropical
c                          2 = temperate
c                          3 = boreal and subpolar
c
c       phenology, third digit of code
c                          0 = mostly evergreen upper canopy
c                          5 = mixed upper canopy
c                          9 = mostly deciduous upper canopy
c
c       C3/C4 carbon pathway, as a fraction of C4 plants
c                          0 = mostly C3 plants
c                          9 = mostly C4 plants
c
c     from this information, set
c
c       one-sided lai for upper canopy, evergreen and deciduous
c       one-sided lai for lower story, C3 and C4
c       top and bottom height in meters for upper and lower story
c
c     also set these (for now) constants
c
c       leaf orientation for upper and lower canopy
c       fractional soil-layer root amounts for upper and lower canopy
c       ... and others ...
c
c break 4-digit veg code into four parts:
c
c - cover type ct
c - temperature zone tz
c - phenology ph
c - carbon path cp
c
              temp = ivegxy(i,j)
c
              cp = mod(temp,10)
              temp = temp/10
              ph = mod(temp,10)
              temp = temp/10
              tz = mod(temp,10)
              ct = temp/10
c
c assign vegetation parameters here:
c
c ------
c forest
c ------
c
              if (ct.eq.1) then
c
c tropical
c
                if (tz.eq.1) then
                  xlaiu          =  6.00   ! total lai - upper
                  xlail          =  0.00   ! total lai - lower
                  saixy(i,j,1)   =  0.00   ! sai - lower
                  saixy(i,j,2)   =  0.50   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  0.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  5.00   ! vegetation base - upper
                  ztopxy(i,j,2)  = 25.00   ! vegetation top  - upper
c
c temperate
c
                else if (tz.eq.2) then
                  xlaiu          =  5.00   ! total lai - upper
                  xlail          =  0.50   ! total lai - lower
                  saixy(i,j,1)   =  0.05   ! sai - lower
                  saixy(i,j,2)   =  0.50   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  0.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  5.00   ! vegetation base - upper
                  ztopxy(i,j,2)  = 20.00   ! vegetation top  - upper
c
c boreal
c
                else if (tz.eq.3) then
                  xlaiu          =  4.00   ! total lai - upper
                  xlail          =  0.50   ! total lai - lower
                  saixy(i,j,1)   =  0.05   ! sai - lower
                  saixy(i,j,2)   =  0.50   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  0.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  5.00   ! vegetation base - upper
                  ztopxy(i,j,2)  = 20.00   ! vegetation top  - upper
                endif
c
c --------
c woodland
c --------
c
              else if (ct.eq.2) then
c
c tropical
c
                if (tz.eq.1) then
                  xlaiu          =  3.50   ! total lai - upper
                  xlail          =  1.00   ! total lai - lower
                  saixy(i,j,1)   =  0.05   ! sai - lower
                  saixy(i,j,2)   =  0.20   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  0.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  5.00   ! vegetation base - upper
                  ztopxy(i,j,2)  = 20.00   ! vegetation top  - upper
c
c temperate
c
                else if (tz.eq.2) then
                  xlaiu          =  2.00   ! total lai - upper
                  xlail          =  1.00   ! total lai - lower
                  saixy(i,j,1)   =  0.05   ! sai - lower
                  saixy(i,j,2)   =  0.20   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  0.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  5.00   ! vegetation base - upper
                  ztopxy(i,j,2)  = 10.00   ! vegetation top  - upper
c
c boreal
c
                else if (tz.eq.3) then
                  xlaiu          =  2.00   ! total lai - upper
                  xlail          =  1.00   ! total lai - lower
                  saixy(i,j,1)   =  0.05   ! sai - lower
                  saixy(i,j,2)   =  0.20   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  0.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  5.00   ! vegetation base - upper
                  ztopxy(i,j,2)  = 10.00   ! vegetation top  - upper
                endif
c
c -------
c savanna
c -------
c
              else if (ct.eq.3) then
c
c tropical
c
                if (tz.eq.1) then
                  xlaiu          =  2.00   ! total lai - upper
                  xlail          =  4.00   ! total lai - lower
                  saixy(i,j,1)   =  0.20   ! sai - lower
                  saixy(i,j,2)   =  0.10   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  1.00   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  3.00   ! vegetation base - upper
                  ztopxy(i,j,2)  =  7.00   ! vegetation top  - upper
c
c temperate
c
                else if (tz.eq.2) then
                  xlaiu          =  1.00   ! total lai - upper
                  xlail          =  3.00   ! total lai - lower
                  saixy(i,j,1)   =  0.20   ! sai - lower
                  saixy(i,j,2)   =  0.10   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  1.00   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  3.00   ! vegetation base - upper
                  ztopxy(i,j,2)  =  7.00   ! vegetation top  - upper
c
c boreal
c
                else if (tz.eq.3) then
                  xlaiu          =  1.00   ! total lai - upper
                  xlail          =  3.00   ! total lai - lower
                  saixy(i,j,1)   =  0.20   ! sai - lower
                  saixy(i,j,2)   =  0.10   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  1.00   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  3.00   ! vegetation base - upper
                  ztopxy(i,j,2)  =  5.00   ! vegetation top  - upper
                endif
c
c ---------
c shrubland
c ---------
c
              else if (ct.eq.4) then
c
                xlaiu          =  1.00   ! total lai - upper
                xlail          =  1.00   ! total lai - lower
                saixy(i,j,1)   =  0.20   ! sai - lower
                saixy(i,j,2)   =  0.10   ! sai - upper
                zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                ztopxy(i,j,1)  =  0.30   ! vegetation top  - lower
                zbotxy(i,j,2)  =  0.50   ! vegetation base - upper
                ztopxy(i,j,2)  =  2.00   ! vegetation top  - upper
c
c ---------
c grassland
c ---------
c
              else if (ct.eq.5) then
c
c tropical
c
                if (tz.eq.1) then
                  xlaiu          =  0.00   ! total lai - upper
                  xlail          =  5.00   ! total lai - lower
                  saixy(i,j,1)   =  0.00   ! sai - lower
                  saixy(i,j,2)   =  0.00   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  1.50   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  0.00   ! vegetation base - upper
                  ztopxy(i,j,2)  =  0.00   ! vegetation top  - upper
c
c temperate
c
                else if (tz.eq.2) then
                  xlaiu          =  0.00   ! total lai - upper
                  xlail          =  4.00   ! total lai - lower
                  saixy(i,j,1)   =  0.00   ! sai - lower
                  saixy(i,j,2)   =  0.00   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  1.00   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  0.00   ! vegetation base - upper
                  ztopxy(i,j,2)  =  0.00   ! vegetation top  - upper
c
c boreal
c
                else if (tz.eq.3) then
                  xlaiu          =  0.00   ! total lai - upper
                  xlail          =  3.00   ! total lai - lower
                  saixy(i,j,1)   =  0.00   ! sai - lower
                  saixy(i,j,2)   =  0.00   ! sai - upper
                  zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                  ztopxy(i,j,1)  =  1.00   ! vegetation top  - lower
                  zbotxy(i,j,2)  =  0.00   ! vegetation base - upper
                  ztopxy(i,j,2)  =  0.00   ! vegetation top  - upper
                endif
c
c ------
c tundra
c ------
c
              else if (ct.eq.6) then
c
                xlaiu          =  0.00   ! total lai - upper
                xlail          =  1.00   ! total lai - lower
                saixy(i,j,1)   =  0.10   ! sai - lower
                saixy(i,j,2)   =  0.00   ! sai - upper
                zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                ztopxy(i,j,1)  =  0.30   ! vegetation top  - lower
                zbotxy(i,j,2)  =  0.00   ! vegetation base - upper
                ztopxy(i,j,2)  =  0.00   ! vegetation top  - upper
c
c ------
c desert
c ------
c
              else if (ct.eq.7) then
c
                xlaiu          =  0.00   ! total lai - upper
                xlail          =  0.10   ! total lai - lower
                saixy(i,j,1)   =  0.01   ! sai - lower
                saixy(i,j,2)   =  0.00   ! sai - upper
                zbotxy(i,j,1)  =  0.02   ! vegetation base - lower
                ztopxy(i,j,1)  =  0.30   ! vegetation top  - lower
                zbotxy(i,j,2)  =  0.00   ! vegetation base - upper
                ztopxy(i,j,2)  =  0.00   ! vegetation top  - upper
c
c ---------
c ice/ocean
c ---------
c
              else
c
                xlaiu          =  0.00   ! total lai - upper
                xlail          =  0.00   ! total lai - lower
                saixy(i,j,1)   =  0.00   ! sai - lower
                saixy(i,j,2)   =  0.00   ! sai - upper
                zbotxy(i,j,1)  =  0.00   ! vegetation base - lower
                ztopxy(i,j,1)  =  0.00   ! vegetation top  - lower
                zbotxy(i,j,2)  =  0.00   ! vegetation base - upper
                ztopxy(i,j,2)  =  0.00   ! vegetation top  - upper
c
              endif
c
c Using simple phenology rules, determine actual canopy lai
 
c determine internal lai parameters
c
              plail4(i,j)  = xlail * (cp / 10.0 + 0.05) !c4 grass lai
              plail3(i,j)  = xlail - plail4(i,j)        !c3 grass lai
c
              plaiud(i,j) = xlaiu * (ph / 10.0 + 0.05) !decid.  tree lai
              plaiue(i,j) = xlaiu - plaiud(i,j)        !evergrn tree lai
c
c fractional cover
c
c if there is a continuous canopy cover (one-sided lai > 1)
c then assume the fractional cover is 1.0
c
              if (plail3(i,j)+plail4(i,j).ge.1.) then
                flxy(i,j) = 1.0
              else
                flxy(i,j) = plail3(i,j) + plail4(i,j)
              endif
 
              if (plaiue(i,j)+plaiud(i,j).ge.1.) then
                fuxy(i,j) = 1.0
              else
                fuxy(i,j) = plaiue(i,j) + plaiud(i,j)
              endif
 
  202       continue
  200     continue
 
          firstfile = .false.
c++++++++++++
        endif
c++++++++++++
 
        do 210 j=1,nlat
          do 212 i=1,nlon
c
            laixy(i,j,1) = (plail3(i,j) + plail4(i,j))
     *                     * dropl(i,j) * templ(i,j)
c
            laixy(i,j,2) = plaiue(i,j)
     *                   + plaiud(i,j) * dropu(i,j) * tempu(i,j)
c
            confracxy(i,j) = plaiue(i,j) / (laixy(i,j,2) + epsilon)
c           confracxy(i,j) = 0.0   ! if no conifer physiology difference
c
            c4fracxy(i,j)  = plail4(i,j) * dropl(i,j) * templ(i,j)
     *                       / (laixy(i,j,1) + epsilon)
c
c adjust lai for fractional cover (re-set to 1 for small lai's with no
c drop-leaves)
c
            laixy(i,j,1) = laixy(i,j,1) / (flxy(i,j) + epsilon)
            laixy(i,j,2) = laixy(i,j,2) / (fuxy(i,j) + epsilon)
 
  212     continue
  210   continue
 
c**********
      endif
c**********
 
c         Set remaining veg parameters
c         (same for interactive and non-interactive veg)
 
      do 500 j=1,nlat
        do 502 i=1,nlon
c
c set orientation: -1 = vertical, 0 = random, 1 = horizontal
c
          orievxy(i,j,1) = -0.5
          orievxy(i,j,2) =  0.0
c
c radiative characteristics:
c
c tau and rho, visible, lower story
c ---------------------------------
c
c adjusted from Sellers et al. (1996) to account for visible
c waveband in GENESIS going upto 0.9 microns:
c
c rhoveg = 46/65 rhoveg(vis) + 19/65 rhoveg(nir)
c tauveg = 46/65 tauveg(vis) + 19/65 tauveg(nir)
c
          tauvegxy(i,j,1,1) = 0.12
          rhovegxy(i,j,1,1) = 0.25
c
c non adjusted
c
c         tauvegxy(i,j,1,1) = 0.07
c         rhovegxy(i,j,1,1) = 0.10
c
c tau and rho, visible, upper story
c ---------------------------------
c
c adjusted from Sellers et al. (1996) to account for visible
c waveband in GENESIS going upto 0.9 microns:
c
c rhoveg = 46/65 rhoveg(vis) + 19/65 rhoveg(nir)
c tauveg = 46/65 tauveg(vis) + 19/65 tauveg(nir)
c
          tauvegxy(i,j,1,2) = 0.09
          rhovegxy(i,j,1,2) = 0.19
c
c non adjusted
c
c         tauvegxy(i,j,1,2) = 0.05
c         rhovegxy(i,j,1,2) = 0.10
c
c tau and rho, near ir, lower story
c ---------------------------------
c
          tauvegxy(i,j,2,1) = 0.25
          rhovegxy(i,j,2,1) = 0.60
c
c tau and rho, near ir, upper story
c ---------------------------------
c
          tauvegxy(i,j,2,2) = 0.20
          rhovegxy(i,j,2,2) = 0.40
c
c linear dimensions for aerodynamic flux parameterization
c
          dleafxy(i,j,1) = .10
          dstemxy(i,j,1) = .10
c
          dleafxy(i,j,2) = .10
          dstemxy(i,j,2) = .20
c
c fraction of lai that is transpiring (0.5 for one side only,
c since lsx will double lai values to use two-sided lai)
c
          fstomaxy(i,j,1) = 0.5
          fstomaxy(i,j,2) = 0.5
c
c Parameters for BATS stomatal conductance parameterization, in case
c requested (phototype = 0) with Wisconsin veg.
c Use D-S biome #1 values for upper story, #7 for lower story, and set
c lsx variables as in vegdat_dorm (except 1/stomatal_resistance).
c
c         Lower story:
          rsminxy(i,j,1)  = 1. / ((2582./1.09) + 110.)
          rsmaxxy(i,j,1)  = 1. / 110.
          parconxy(i,j,1) = (2582./110.) + 1.09
          vpdconxy(i,j,1) =  (1./.0238)*100.
c         Upper story:
          rsminxy(i,j,2)  = 1. / ((2336./.01) + 154.)
          rsmaxxy(i,j,2)  = 1. / 154.
          parconxy(i,j,2) = (2336./154.) + .01
          vpdconxy(i,j,2) = (1./.0273)*100.
c
 502    continue
 500  continue
 
      return
      end
c
c ---------------------------------------------------------------------
c
      subroutine ibis_stat (dtgcm, nsrest_ibis, nstep, nrstrt, nstop,
     *                      istepy, mbdate, jyear, jmnth, jday)
 
c        Accumulates daily, monthly and annual mean statistics needed
c        by IBIS (Wisconsin) vegetation: by subr wisc_veg for phenology,
c        and at end of each year by subrs carbon and dynaveg for
c        interactive veg only.
 
c        All arguments are supplied:
c        DTGCM    = timestep for run lsx this call (=agcm timestep) (s)
c        NSREST_IBIS = -1/[0,1] flag for coldstart/restart (suppl)
c        NSTEP    = agcm timestep counter
c        NRSTRT   = first agcm timestep # for this run
c        NSTOP    = value of NSTEP for last timestep this run
c        ISTEPY   = timestep within current 365-day year
c        MBDATE   = base date for curr. sequence of runs (yymmdd)
c        JYEAR    = current year number (2001 AD,etc, or -21000 BP,etc)
c        JMNTH    = current month number (1-12)
c        JDAY     = current day number (1-31)
 
#include <compar.h>
#include <comgrd.h>
#include <commudxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
 
      character cyr*20, cdate*20
      logical if_bod, if_eod, if_bom, if_eom, if_boy, if_eoy
      dimension ndaypm(12)
      save ndaypm, nytimes, nmtimes, ndtimes
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
 
c        Accumulate current values for non-dynamic or dynamic veg:
c        15-day exponential running mean:
 
      zw15 = exp(-dtgcm/(15.*86400.))
      do j=1,nlat
        do i=1,nlon
          t15runxy(i,j) = zw15*t15runxy(i,j) + (1.-zw15)*ts2xy(i,j)
        enddo
      enddo
 
c        Remainder is only needed for dynamic vegetation
 
c=============
#ifdef ibisint
c=============
c****************************
      if (dynaibis.eq.1) then
c****************************
 
c        Set beginning/end_of_day/month/year flags
 
      if_bod = .false.
      if_eod = .false.
      if_bom = .false.
      if_eom = .false.
      if_boy = .false.
      if_eoy = .false.
 
      n = nint(86400./dtgcm)
      if (mod(nstep,n).eq.1) then
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
 
c        For debugging of short runs, force eo[d,m,y] if last timestep
 
      if (nstep.eq.nstop) then
        if_eod = .true.
        if_eom = .true.
        if_eoy = .true.
      endif
 
c     write(*,8880)nstep, if_bod, if_eod, if_bom, if_eom, if_boy, if_eoy
c8880 format('ibis_stat: nstep=',i6,' if[bod,eod,bom,eom,boy,eoy]=',6l3)
 
c       If beginning of year (or first timestep), reset annual arrays
c       and counter
 
      if (if_boy .or. nstep.eq.nrstrt) then
        call zero   (anppuxy,      nlon*nlat*npftu)
        call zero   (anpplxy,      nlon*nlat*npftl)
        call zero   (agppxy,       nlon*nlat)
        call resetr (tcoldxy,      nlon*nlat,  1.e20)
        call resetr (twarmxy,      nlon*nlat, -1.e20)
        call zero   (gdd0,         nlon*nlat)
        call zero   (gdd5,         nlon*nlat)
        call zero   (avgfrespaxy,  nlon*nlat)
        call zero   (avgfrespbxy,  nlon*nlat)
        nytimes = 0
      endif
 
c       If beginning of month (or first timestep), reset monthly arrays
c       and counter
 
      if (if_bom .or. nstep.eq.nrstrt) then
        call zero (tmonxy, nlon*nlat)
        nmtimes = 0
      endif
 
c       If beginning of day (or first timestep), reset daily arrays
c       and counter (tdailxy currently used only as diagnostic)
 
      if (if_bod .or. nstep.eq.nrstrt) then
        call zero (tdailxy, nlon*nlat)
        ndtimes = 0
      endif
 
c        Increment counters
 
      nytimes = nytimes + 1
      nmtimes = nmtimes + 1
      ndtimes = ndtimes + 1
 
c        Accumulate current values for dynamic veg:
c        tdailxy     = 24-hr mean 2-m temperature (K) (diagnostic only)
c        tmonxy      = monthly mean 2-m temperature (K)
c        anpp[u,l]xy = annual npp, without yet adjusting for stem and
c                      root respiration, or growth respiration (kg/m2/s)
c        agppxy      = annual gpp, all veg combined (diag only)(kg/m2/s)
c        gdd[0,5] = degree days above 0 or 5 deg C (days*deg C)
c        avgfresp[a,b]xy = temperature-dependence term for respiration
c                          for upper stems and top ~1m of soil
 
      do j=1,nlat
        do i=1,nlon
c            Daily:
 
          tdailxy(i,j) = ( (ndtimes-1)*tdailxy(i,j) + ts2xy(i,j) )
     *                   / ndtimes
 
c            Monthly:
 
          tmonxy(i,j)  = ( (nmtimes-1)*tmonxy(i,j)  + ts2xy(i,j) )
     *                   / nmtimes
c            Annual:
 
          anppuxy(i,j,1) = (
     >      (nytimes-1) * anppuxy(i,j,1) +
     >       1.0e+06 * anudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,1) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anppuxy(i,j,2) = (
     >      (nytimes-1) * anppuxy(i,j,2) +
     >       1.0e+06 * anudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,2) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anppuxy(i,j,3) = (
     >      (nytimes-1) * anppuxy(i,j,3) +
     >       1.0e+06 * anudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,3) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anppuxy(i,j,4) = (
     >      (nytimes-1) * anppuxy(i,j,4) +
     >       1.0e+06 * anucxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,4) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anppuxy(i,j,5) = (
     >      (nytimes-1) * anppuxy(i,j,5) +
     >       1.0e+06 * anudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,5) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anppuxy(i,j,6) = (
     >      (nytimes-1) * anppuxy(i,j,6) +
     >       1.0e+06 * anucxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,6) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anppuxy(i,j,7) = (
     >      (nytimes-1) * anppuxy(i,j,7) +
     >       1.0e+06 * anudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,7) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anpplxy(i,j,1) = (
     >      (nytimes-1) * anpplxy(i,j,1) +
     >       1.0e+06 * anl3xy(i,j) *
     >       laixy(i,j,1) * flxy(i,j) * fraclxy(i,j,1) *
     >       (1. - fixy(i,j)) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          anpplxy(i,j,2) = (
     >      (nytimes-1) * anpplxy(i,j,2) +
     >       1.0e+06 * anl4xy(i,j) *
     >       laixy(i,j,1) * flxy(i,j) * fraclxy(i,j,2) *
     >       (1. - fixy(i,j)) *
     >       365.0 * 86400.0 * 12.0e-9
     >                     ) / nytimes
c
          agppxy(i,j) =  (
     >      (nytimes-1) * agppxy(i,j) +
c
     >       1.0e+06 * agudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,1) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,2) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,3) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agucxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,4) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,5) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agucxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,6) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agudxy(i,j) *
     >       laixy(i,j,2) * fuxy(i,j) * fracuxy(i,j,7) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agl3xy(i,j) *
     >       laixy(i,j,1) * flxy(i,j) * fraclxy(i,j,1) *
     >       (1. - fixy(i,j)) *
     >       365.0 * 86400.0 * 12.0e-9   +
c
     >       1.0e+06 * agl4xy(i,j) *
     >       laixy(i,j,1) * flxy(i,j) * fraclxy(i,j,2) *
     >       (1. - fixy(i,j)) *
     >       365.0 * 86400.0 * 12.0e-9
c
     >                     ) / nytimes
c
          gdd0(i,j) = gdd0(i,j) + max (ts2xy(i,j)-tmelt,    0.)
     *                            * dtgcm/86400.
          gdd5(i,j) = gdd5(i,j) + max (ts2xy(i,j)-tmelt-5., 0.)
     *                            * dtgcm/86400.
 
        enddo
      enddo
 
c Accumulate annual stem and soil respiration temperature-dependence
c terms on a 15 C base. First, compute mean temperature of first
c 4 layers of soil (~1 m). (Only needed for dynamic vegetation).
 
      zdep = hsod(1) + hsod(2) + hsod(3) + hsod(4)
      call zero (workl, nlon*nlat)
      do 250 k=1,4
        do 252 j=1,nlat
          do 254 i=1,nlon
            workl(i,j) =  workl(i,j) + tsoixy(i,j,k)*hsod(k)/zdep
 254      continue
 252    continue
 250  continue
 
      do 260 j=1,nlat
        do 262 i=1,nlon
          zfunca = exp(3000.0 *
     >             (1. / 288.16 - 1. / tsxy(i,j) ))    ! upper stem temp
c
          zfuncb = exp(3000.0 *
     >             (1. / 288.16 - 1. / workl(i,j) ))   ! mean soil temp
c
          avgfrespaxy(i,j) = ((nytimes-1) * avgfrespaxy(i,j) +
     >                         zfunca) / nytimes
c
          avgfrespbxy(i,j) = ((nytimes-1) * avgfrespbxy(i,j) +
     >                         zfuncb) / nytimes
 262    continue
 260  continue
 
c        If end of month,set monthly quantities(only needed for dyn.veg)
c
c        tcoldxy     = minimum monthly mean 2-m temperature (K)
c        twarmxy     = maximum monthly mean 2-m temp. (dyn veg only) (K)
 
      if (if_eom) then
 
        do 300 j=1,nlat
          do 302 i=1,nlon
            tcoldxy(i,j) = min (tcoldxy(i,j), tmonxy(i,j))
            twarmxy(i,j) = max (twarmxy(i,j), tmonxy(i,j))
  302     continue
  300   continue
 
      endif
 
c
c        If end of year, save tcoldprevxy for
c        next year's phenology (written to Genesis restart file),
c        and twarmprevxy, gdd[0,5]prev for dynamic veg.
c        If dynamic veg, also call carbon and dynaveg, and write
c        IBIS restart file.
 
c to avoid bug with palaeo calendar:

c     if (if_eoy) then
      if (if_eoy .and. nytimes.gt.360) then   ! ibis2.0

c          Apply a 5-year exponential running-mean time filter
 
        zwei = exp(-1./5.)
        do j=1,nlat
          do i=1,nlon
            tcoldprevxy(i,j) =     zwei *tcoldprevxy(i,j)
     *                       + (1.-zwei)*tcoldxy(i,j)
            twarmprevxy(i,j) =     zwei *twarmprevxy(i,j)
     *                       + (1.-zwei)*twarmxy(i,j)
            gdd0prev(i,j)    =     zwei *gdd0prev(i,j)
     *                       + (1.-zwei)*gdd0(i,j)
            gdd5prev(i,j)    =     zwei *gdd5prev(i,j)
     *                       + (1.-zwei)*gdd5(i,j)
          enddo
        enddo
 
c       Assemble current year and date (for diagnostics)
        write (cyr,'(i20)') jyear
        call lstrip (cyr)
 
        n = nint(86400./dtgcm)
        zhr = (mod(nstep-1,n)+1) * dtgcm/3600.
        write (cdate,"(a,'/',i2.2,'/',i2.2,' ',f5.2,' GMT')")
     *    cyr(1:lenchr(cyr)), jmnth, jday, zhr
 
        call carbon
        call dynaveg (1, cdate, jyear)
        call ibis_write (dtgcm, nstep, istepy, cyr)
 
      endif
 
c**********
      endif
c**********
c=====
#endif
c=====
 
      return
      end
c
c ---------------------------------------------------------------------
c
c=============
#ifdef ibisint
c=============
 
      subroutine ibis_write (dtgcm, nstep, istepy, cyr)
 
c        *** Should be called only at the end of the last time step
c        in a year, after carbon and dynaveg are called ***.
c        Writes IBIS restart file and disposes to its permanent
c        location. Assembles permanent name by appending
c        "_year" to requested prefix (namelist parameter IBISRES).
 
c        All arguments are supplied:
c        DTGCM    = timestep for run lsx this call (=agcm timestep) (s)
c        NSTEP    = agcm timestep counter
c        ISTEPY   = timestep within current 365-day year
c        CYR      = current year number (character)
 
#include <compar.h>
#include <comgrd.h>
#include <commudxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
 
#ifdef gcm
#include <COMLUN.h>
#else
      data nsreibis /95/
#endif
 
      character cyr*(*)
      character cname*200, clab*8, cver*8
      save cname
      data versibis /1.01/
 
 
c        Assemble filename, and open a new IBIS restart file.
 
#ifdef gcm
      if (ibisres.eq.' ') then
        cname = msnamhis(ifile)(1:lenchr(msnamhis(ifile)))
     *          // 'RESI'
      else
        cname = ibisres(1:lenchr(ibisres))
     *          // '_' // cyr(1:lenchr(cyr))
      endif
#endif
      call openf (nsreibis, cname, mspath, 1, 0)
 
c        Write header words
 
#ifdef cray
      clab = 'IBISCRAY'
#else
      clab = 'IBISIEEE'
#endif
      write (cver,'(f8.2)') versibis
      write (nsreibis) clab, cver, nlon, nlat
 
c        Write data fields
 
      call writer (nsreibis, cluxy,  nlon, nlat, npftu, 2)
      call writer (nsreibis, cwuxy,  nlon, nlat, npftu, 2)
      call writer (nsreibis, cruxy,  nlon, nlat, npftu, 2)
      call writer (nsreibis, cllxy,  nlon, nlat, npftl, 2)
      call writer (nsreibis, crlxy,  nlon, nlat, npftl, 2)
 
      call writer (nsreibis, existu, nlon, nlat, npftu, 2)
      call writer (nsreibis, existl, nlon, nlat, npftl, 2)
 
      call writer (nsreibis, tdailxy,       nlon, nlat, 1, 2)
      call writer (nsreibis, t15runxy,      nlon, nlat, 1, 2)
      call writer (nsreibis, tcoldprevxy,   nlon, nlat, 1, 2)
      call writer (nsreibis, twarmprevxy,   nlon, nlat, 1, 2)
      call writer (nsreibis, gdd0prev,      nlon, nlat, 1, 2)
      call writer (nsreibis, gdd5prev,      nlon, nlat, 1, 2)
      call writer (nsreibis, uptotlxy,      nlon, nlat, 1, 2)
      call writer (nsreibis, uptotuxy,      nlon, nlat, 1, 2)
      call writer (nsreibis, agdd0xy,       nlon, nlat, 1, 2)
      call writer (nsreibis, agdd5xy,       nlon, nlat, 1, 2)
      call writer (nsreibis, flxy,          nlon, nlat, 1, 2)
      call writer (nsreibis, fuxy,          nlon, nlat, 1, 2)
 
c        Write last-word check
 
      write (nsreibis) cver
 
      call closef (nsreibis)
 
c        Dispose IBIS restart file to permanent location
 
#ifdef gcm
      call closef (nsreibis)
      call disposef (cname, mspath)
 
      write (ioterm2,200)
     *  nstep, istepy*dtgcm/86400., nstep*dtgcm/86400.,
     *  mspath(1:lenchr(mspath)), cname(1:lenchr(cname))
  200 format(
     *       /' ----------------------------------------------'
     *       /' **** IBIS VEGETATION RESTART FILE OUTPUT *****'
     *       /' ----------------------------------------------'
     *       /' timestep                          = ',i10
     *       /' calendar days                     = ',f10.3
     *       /' elapsed days                      = ',f10.3
     *       /' IBIS vegetation restart file name = ',a,a
     *       /' ----------------------------------------------')
#ifdef mss
      if (mspath.eq.' ') write (ioterm2,202)
  202 format( ' *  File NOT saved to MSS (MSPATH is blank!)  *'
     *       /' ----------------------------------------------')
#endif
      write (ioterm2,*)
#endif
 
      return
      end
 
c=====
#endif
c=====
c
c ---------------------------------------------------------------------
c
c=============
#ifdef ibisint
c=============
 
      subroutine ibis_read (nsrest_ibis, lonreslsx,latreslsx, typeresin)
 
c        Reads IBIS restart/initial file, or sets simple initial
c        values (cold start), depending on NSREST_IBIS
 
c        Supplied:
c        NSREST_IBIS = -1/[0,1] flag for coldstart/restart (suppl)
c        LONRESLSX  = expected IBIS restart-file longitude resolution
c        LATRESLSX  = expected IBIS restart-file latitude  resolution
c        TYPERESIN  = label (CRAY or IEEE) for type conversions
 
c        Copies IBIS restart file from permanent location and reads it.
 
#include <compar.h>
#include <comgrd.h>
#include <commudxy.h>
#include <commugxy.h>
#include <comsno.h>
#include <comsnoxy.h>
#include <comsoi.h>
#include <comsoixy.h>
#include <comveg.h>
#include <comvegxy.h>
 
#ifdef gcm
#include <COMLUN.h>
#else
      character ibisresin*16, ibispathin*80
      data ibisresin, ibispathin /' ', ' '/
      data nsreibis /95/
#endif
 
      character typeresin*(*), typeresinsav*8, typeibisin*8
      character clab*8, cver*8, cpath*80, cfile*80
 
c     added for "initial" ascii file:
      dimension speclau(npftu), speclal(npftl)
 
c     version number for binary file:
      parameter (nverspre=2)
      dimension verspre(nverspre)
      save verspre
      data verspre / 1.0, 1.01 /
 
 
c         Acquire and read IBIS restart file (nsrest_ibis = 0,1)
c         or do IBIS cold start (nsrest_ibis= -1)
 
c-----------------------------------------------------
c     no distinction between "initial" and "re" start for ibis yet
      if (nsrest_ibis.eq.0 .or. nsrest_ibis.eq.1) then
c-----------------------------------------------------
 
        write (ioterm2,"(/1x,50('-'))")
 
#ifdef gcm
        call acquiref (ibisresin, ibispathin)
#endif
 
        write (ioterm2,1000) ibispathin(1:lenchr(ibispathin)), ibisresin
 1000   format (' Opening binary IBIS vegetation restart file ',a,'/',a)
 
        call openf (nsreibis, ibisresin, ibispathin, 0, 0)
 
        write(ioterm2,1002)
 1002   format (/4x,'Reading header records ...')
 
        read (nsreibis, err=8000, iostat=ios) clab, cver, nlonx, nlatx
 
        if (clab.eq.'IBISCRAY' .or. clab.eq.'IBISIEEE') then
          typeibisin = clab(5:8)
          versibisin = 0.
          read (cver,'(f8.2)',err=1004) versibisin
 1004     continue
        else
          write(ioterm2,1006) ibisresin(1:lenchr(ibisresin)), clab
 1006     format (//'*** Incorrect first-word label of IBIS',
     *              ' initial/restart file ',a,':'
     *             /4x,'First 8 bytes   : ', a
     *             /4x,'Expected labels : IBISCRAY or IBISIEEE'
     *             /4x,'Stop in ibis_read.'/)
          call endrun (-1)
        endif
 
        do 1010 i=1,nverspre
          if (abs(versibisin-verspre(i)).lt.epsilon) then
            versibisin = verspre(i)
            goto 1016
          endif
 1010   continue
        write(ioterm2,1012) ibisresin(1:lenchr(ibisresin)), versibisin,
     *                     (verspre(i),i=1,nverspre)
 1012   format (
     * //'*** Incompatible version number of IBIS initial/restart file:'
     *    /4x,'Version number of IBIS initial/restart file ',a,':'
     *    /4x, f11.7
     *    /4x,'Version numbers accepted by program:'
     *    /(4x,9(f6.2)) )
        write(ioterm2,1014)
 1014   format(
     *    /4x,'Stop in ibis_read.'/)
        call endrun (-1)
 1016   continue
 
#ifdef gcm
        call convertype (nlonx, 1, 'I', typeibisin)
        call convertype (nlatx, 1, 'I', typeibisin)
#endif
 
        write(ioterm2,1020)
     *    versibisin, typeibisin, nlonx, nlatx
 1020   format ( 7x,'file version number :',f11.2
     *          /7x,'file machine type   : ',a
     *          /7x,'file resolutions    :',2i6)
 
        if (nlonx.ne.lonreslsx .or. nlatx.ne.latreslsx) then
          write(ioterm2,1030) nlonx, lonreslsx, nlatx, latreslsx
 1030     format (
     *      //' Resolutions of IBIS initial/restart file do not',
     *        ' match the expected'
     *      / ' values (i.e., the current model''s, or specified by',
     *        ' namelist HORRESLSX[])'
     *      /4x,'Initial/restart file NLON:',i8,'   expected NLON:',i8
     *      /4x,'Initial/restart file NLAT:',i8,'   expected NLAT:',i8
     *      /4x,'Stop in ibis_read.'/ )
          call endrun (-1)
        endif
 
c          Read IBIS data arrays
 
        write(ioterm2,6000)
 6000   format (/4x,'Reading data   fields...')
 
c          Temporarily set typeresin for the IBIS file type
c          (typeresin is passed here from common comhed, used by subrs
c          convertype and reader in util.f for all machine-type
c          conversions). typeresin will be reset to passed value
c          before returning below.
 
        typeresinsav = typeresin
        typeresin = typeibisin
 
        call reader (nsreibis, cluxy,  nlon, nlat, npftu, 2)
        call reader (nsreibis, cwuxy,  nlon, nlat, npftu, 2)
        call reader (nsreibis, cruxy,  nlon, nlat, npftu, 2)
        call reader (nsreibis, cllxy,  nlon, nlat, npftl, 2)
        call reader (nsreibis, crlxy,  nlon, nlat, npftl, 2)
 
        call reader (nsreibis, existu, nlon, nlat, npftu, 2)
        call reader (nsreibis, existl, nlon, nlat, npftl, 2)
 
c         uptot[l,u]xy, agdd[0,5]xy and t15runxy have already
c         been initialized in readlsx/inilsx, since they are on the
c         Genesis restart file (since needed for non-dynamic phenology).
c         However, if they are on the Ibis restart file, they are
c         overwritten here.
 
        if (versibisin.ge.1.01) then
          call reader (nsreibis, tdailxy,       nlon, nlat, 1, 2)
          call reader (nsreibis, t15runxy,      nlon, nlat, 1, 2)
          call reader (nsreibis, tcoldprevxy,   nlon, nlat, 1, 2)
          call reader (nsreibis, twarmprevxy,   nlon, nlat, 1, 2)
          call reader (nsreibis, gdd0prev,      nlon, nlat, 1, 2)
          call reader (nsreibis, gdd5prev,      nlon, nlat, 1, 2)
          call reader (nsreibis, uptotlxy,    nlon, nlat, 1, 2)
          call reader (nsreibis, uptotuxy,    nlon, nlat, 1, 2)
          call reader (nsreibis, agdd0xy,     nlon, nlat, 1, 2)
          call reader (nsreibis, agdd5xy,     nlon, nlat, 1, 2)
          call reader (nsreibis, flxy,        nlon, nlat, 1, 2)
          call reader (nsreibis, fuxy,        nlon, nlat, 1, 2)
        else
c         gtransmaxuxy:
          call reader (nsreibis, workl,         nlon, nlat, 1, 2)
          call reader (nsreibis, tdailxy,       nlon, nlat, 1, 2)
          call reader (nsreibis, tcoldprevxy,   nlon, nlat, 1, 2)
c         pmonminprevxy:
          call reader (nsreibis, workl,         nlon, nlat, 1, 2)
          call reader (nsreibis, twarmprevxy,   nlon, nlat, 1, 2)
          call reader (nsreibis, gdd0prev,      nlon, nlat, 1, 2)
          call reader (nsreibis, gdd5prev,      nlon, nlat, 1, 2)
c         flag for later detection of no flxy,fuxy info, for dynaveg(0):
          call resetr (flxy, nlon*nlat, 1.e20)
        endif
 
c          Last-word check
 
        cver = ' '
        read (nsreibis, err=6502, iostat=ios) cver
        read (cver,'(f8.2)',err=6502) zver
        do 6500 i=1,nverspre
          if (abs(zver-verspre(i)).lt.epsilon) zver = verspre(i)
 6500   continue
 6502   if (zver.ne.versibisin) then
          write(ioterm2,6504) versibisin, zver, cver, ios
 6504     format(/' *** Error: incorrect last word on',
     *                  ' IBIS initial/restart file'
     *           /'     expected     = ',f8.3
     *           /'     found (real) = ',f25.15
     *           /'     found (char) = ',a8
     *           /'     iostat       = ',i8
     *           /'     Stop in ibis_read.'/)
          call endrun (-1)
        endif
 
c          Reset passed machine-type for restart files
 
        typeresin = typeresinsav
 
        call closef (nsreibis)
 
        write (ioterm2,"(1x,50('-'))")
 
c-------------------------------------
      else if (nsrest_ibis.eq.-1) then
c-------------------------------------
 
c          IBIS cold start. Initialize to small amounts of veg on all
c          land points. zlai[u,l] are initial total [upper,lower]
c          one-sided lai's.
 
        versibisin = 0.
 
        zlaiu = .35
        zlail = .10
 
        do 750 j=1,nlat
          do 752 i=1,nlon
 
            if (lmaskxy(i,j).eq.1) then
 
              do 760 n=1,npftu
                cluxy(i,j,n) = zlaiu/(21.*npftu) !avg speclau in dynaveg
                cwuxy(i,j,n) = 0.
                cruxy(i,j,n) = 0.
                existu(i,j,n)= 1.
  760         continue
              do 762 n=1,npftl
                cllxy(i,j,n) = zlail/(20.*npftl)     !speclal in dynaveg
                crlxy(i,j,n) = 0.
                existl(i,j,n)= 1.
  762         continue
 
            else
 
              do 770 n=1,npftu
                cluxy(i,j,n) = 0.
                cwuxy(i,j,n) = 0.
                cruxy(i,j,n) = 0.
                existu(i,j,n)= 0.
  770         continue
              do 772 n=1,npftl
                cllxy(i,j,n) = 0.
                crlxy(i,j,n) = 0.
                existl(i,j,n)= 0.
  772         continue
 
            endif
 
  752     continue
  750   continue
 
c----------
      endif
c----------
 
c        For cold starts, set remaining (diagnostic or accumulated)
c        variables
 
c--------------------------------
      if (nsrest_ibis.eq.-1) then
c--------------------------------
 
        call resetr (tdailxy,       nlon*nlat, 273.16 + 15.)
        call resetr (tcoldprevxy,   nlon*nlat, 273.16 -  5.)
        call resetr (twarmprevxy,   nlon*nlat, 273.16 + 25.)
        call resetr (gdd0prev,      nlon*nlat, 0.)
        call resetr (gdd5prev,      nlon*nlat, 0.)
 
c          uptot[l,u]xy, t15runxy and agdd[0,5]xy have already been
c          initialized in readlsx/inilsx, since they are on the Genesis
c          restart file (since needed for non-dynamic phenology).
 
c       call resetr (uptotlxy,      nlon*nlat, 0.)
c       call resetr (uptotuxy,      nlon*nlat, 0.)
c       call resetr (t15runxy,      nlon*nlat, 273.16 + 15.)
c       call resetr (agdd0xy,       nlon*nlat, 0.)
c       call resetr (agdd5xy,       nlon*nlat, 0.)
 
c          flag for later detection of no flxy,fuxy info, for dynaveg(0)
 
        call resetr (flxy, nlon*nlat, 1.e20)
 
c----------
      endif
c----------
 
      return
 
c=====
 
c        Fatal errors for reading restart-file header records
 
 8000 write(ioterm2,8002) ios
 8002 format(/' *** Error or EOF reading header record number',i3,
     *              ' of IBIS initial/restart file.'
     *       /'     iostat = ',i8
     *       /'     Stop in ibis_read.'/)
      call endrun (-1)
      end
 
c=====
#endif
c=====
c
c ---------------------------------------------------------------------
c
c=============
#ifdef ibisint
c=============
 
      subroutine carbon
 
c        Called from ibis_stat at the end of each year, before dynaveg.
c        Adjusts net primary productivity (anpp[u,l]xy) for each
c        plant type, due to sapwood respiration (upper story) and
c        root respiration (both stories). Then dynaveg will use
c        anpp[u,l]xy to step the dynamic carbon stores through one year.
c
c        Also sets diagnostic annual all-veg (adjusted) npp (anppxy).
 
#include <compar.h>
#include <comgrd.h>
#include <comvegxy.h>
 
c
c local variables
c
      rgrowth = 0.3333
c
c initialize diagnostic all-veg npp
c
      call zero (anppxy, nlon*nlat)
c
c--------------------------------
c woody plant types (upper story)
c--------------------------------
c
      do 100 n=1,npftu
        do 102 j=1,nlat
          do 104 i=1,nlon
c
c sapwood respiration
c
            arsuxy = 0.02 *cwuxy(i,j,n) *sapfrac(i,j) *avgfrespaxy(i,j)
c
c fine root respiration
c
            arruxy = 0.20 * cruxy(i,j,n) * avgfrespbxy(i,j)
c
c determine total npp
c
            anppuxy(i,j,n) = max(0.0, ((1.0 - rgrowth)  * (
     >                                 anppuxy(i,j,n) -
     >                                 arsuxy         -
     >                                 arruxy             )
     >                          )     )
c
            anppxy(i,j) = anppxy(i,j) + anppuxy(i,j,n)
c
 104      continue
 102    continue
 100  continue
c
c-------------------------------------
c herbaceous plant types (lower story)
c-------------------------------------
c
      do 200 n=1,npftl
        do 202 j=1,nlat
          do 204 i=1,nlon
c
c fine root respiration
c
            arrlxy = 0.200 * crlxy(i,j,n) * avgfrespbxy(i,j)
c
c determine total npp
c
            anpplxy(i,j,n) = max(0.0, ((1.0 - rgrowth)  * (
     >                                 anpplxy(i,j,n) -
     >                                 arrlxy             )
     >                          )     )
c
            anppxy(i,j) = anppxy(i,j) + anpplxy(i,j,n)
c
 204      continue
 202    continue
 200  continue
c
      return
      end
 
c=====
#endif
c=====
c
c ---------------------------------------------------------------------
c
c=============
#ifdef ibisint
c=============
 
      subroutine dynaveg (iflag, cdate, iyr)

c        IBIS dynamic vegetation, called at end of each year from subr
c        ibis_stat (IFLAG=1). Also called once at start of run from
c        wisc_veg (IFLAG=0) to set non-prognostic arrays.
c        If IFLAG=1, steps carbon stores c[l,w,r]uxy and c[l,r]lxy
c        for one calendar year, and sets climate range limits
c        exist[u,l]xy, using statistics accumulated through the
c        previous year by ibis_stat. Then sets non-prognostic arrays 
c        (plai[u,l]xy, f[l,u]xy, frac[u,l]xy, etc.
c        Also calls mapibis to write out sets of 2-D maps.
c
c        IFLAG = 0 just to initialize *diagnostic* arrays (only once
c                  at start of run, called from vegdat_wisc)
c                1 for normal call at end of each year (from ibis_stat)
c        CDATE = current date (for ascii IBIS output file header)
c
#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
#include <comvegxy.h>

#ifdef gcm
#include <COMLUN.h>
#endif

      character cdate*(*)
      character clabel*100
      integer iyr

c        Local arrays (not saved)...not used elsewhere in ibis
c        except temporarily in this subr and/or for diagnostic output

      dimension 
     *  biomassuxy(nlon,nlat,npftu),  biomasslxy(nlon,nlat,npftu),
     *  effhtuxy(nlon,nlat),          effhtlxy(nlon,nlat),
     *  woodmasstu(nlon,nlat),
     *  totlaiuxy(nlon,nlat),         totlailxy(nlon,nlat),
     *  totbiouxy(nlon,nlat),         totbiolxy(nlon,nlat),
     *  disturbf(nlon,nlat),          disturbo(nlon,nlat)

c        pft characteristics (saved)

      dimension aleafu(npftu),   aleafl(npftl),
     *          awoodu(npftu),
     *          arootu(npftu),   arootl(npftl),
     *          tauleafu(npftu), tauleafl(npftl),
     *          tauwoodu(npftu), tauwood(npftu),
     *          speclau(npftu),  speclal(npftl)
      save tauleafu,  tauleafl, tauwood, tauroot,
     *     speclau,   speclal,
     *     totlaiuxy, totlailxy

      logical first
      save first
      data first /.true./ 


      if (first) then

c ***************************
c specify pft characteristics
c ***************************
 
c turnover parameters:

        tauleafu(1) = 1.02
        tauleafu(2) = 1.00
        tauleafu(3) = 1.02
        tauleafu(4) = 2.00
        tauleafu(5) = 1.00
        tauleafu(6) = 3.25
        tauleafu(7) = 1.00
        tauleafl(1) = 1.50
        tauleafl(2) = 1.25

        tauwood(1) = 50.0
        tauwood(2) = 50.0
        tauwood(3) = 50.0
        tauwood(4) = 50.0
        tauwood(5) = 50.0
        tauwood(6) = 50.0
        tauwood(7) = 50.0

        tauroot = 1

c specific leaf area parameters:
c also defined elsewhere!!!

        speclau(1) = 25.0
        speclau(2) = 25.0
        speclau(3) = 25.0
        speclau(4) = 12.5
        speclau(5) = 25.0
        speclau(6) = 12.5
        speclau(7) = 25.0

        speclal(1) = 20.0
        speclal(2) = 20.0

        first = .false.

      endif

c **********************************************
c step dynamical carbon stores through one year,
c and set  climate-range limits exist[u,l]xy,
c for call at end of each year.
c **********************************************

c     kyear = 2005          ! used below in several places
      kyear = -999999999    ! always have dynamics

c>>>>>>>>>>>>>>>>>>>>>>>>>
      if (iflag.eq.1) then
c>>>>>>>>>>>>>>>>>>>>>>>>>

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (iyr.gt.kyear) then ! veg dynamics after yr kyear
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c Iterate multiple times if requested
c "periodically synchronous coupling"
c (default for asynibis is 1)

c-------------------------------
      do 100 iloop = 1, asynibis
c-------------------------------

       woodnorm = 4.5 !used with woodold and woodnew below
       shape = -3.
       fnorm = 1.0 - exp(shape)

       do 102 j=1,nlat
        do 104 i=1,nlon

        if (lmaskxy(i,j).eq.1) then

         wood = 0. !used for fuxy further down

c fixed disturbance regime

         disturbf(i,j) = 0.005
         disturbo(i,j) = 0.005

         do 106 k = 1, npftu

          aleafu(k) = 0.25 + 0.10 * exp(-0.5 * totlaiuxy(i,j))
          arootu(k) = 1.25 * aleafu(k)
          awoodu(k) = 1. - aleafu(k) - arootu(k)
          taufin = 5.
          tauwoodu(k) = tauwood(k) - (tauwood(k) - taufin) *
     >                  (1. - existu(i,j,k))

c update carbon in leaves
c -----------------------

          cluxy(i,j,k) = cluxy(i,j,k) * exp(-1./tauleafu(k)) +
     >                   aleafu(k) * tauleafu(k)
     >                   * max(0., anppuxy(i,j,k)) *
     >                   (1. - exp(-1./tauleafu(k)))
          cluxy(i,j,k) = max(cluxy(i,j,k) *
     >                       (1. - disturbf(i,j) - disturbo(i,j)),
     >                       existu(i,j,k) * 0.05 / speclau(k))
                      
          if (k.le.npftl) then

            aleafl(k) = 0.30 + 0.20 * exp(-0.5 * totlailxy(i,j))
            arootl(k) = 1.0 - aleafl(k)

            cllxy(i,j,k) = cllxy(i,j,k) * exp(-1./tauleafl(k)) +
     >                     aleafl(k) * tauleafl(k)
     >                     * max(0., anpplxy(i,j,k)) *
     >                     (1. - exp(-1./tauleafl(k)))
            cllxy(i,j,k) = max(cllxy(i,j,k) *
     >                         (1. - disturbf(i,j) - disturbo(i,j)),
     >                         existl(i,j,k) * 0.05 / speclal(k))

          endif

c update carbon in wood
c ---------------------

          cwuxy(i,j,k) = cwuxy(i,j,k) * exp(-1./tauwoodu(k)) +
     >                   awoodu(k) * tauwoodu(k) 
     >                   * max(0., anppuxy(i,j,k)) *
     >                   (1. - exp(-1./tauwoodu(k)))
          cwuxy(i,j,k) = max(0.0, cwuxy(i,j,k))
          wood = wood + cwuxy(i,j,k)
          cwuxy(i,j,k) = max(0.0, cwuxy(i,j,k) *
     >                       (1. - disturbf(i,j) - disturbo(i,j)))

c update carbon in fine roots
c ---------------------------

c recall: tauroot=1

          cruxy(i,j,k) = cruxy(i,j,k) * exp(-1./tauroot) +
     >                   arootu(k) * tauroot
     >                   * max(0., anppuxy(i,j,k)) *
     >                   (1. - exp(-1./tauroot))
          cruxy(i,j,k) = max(0.0, cruxy(i,j,k) *
     >                       (1. - disturbf(i,j) - disturbo(i,j)))

          if (k.le.npftl) then

            crlxy(i,j,k) = crlxy(i,j,k) * exp(-1./tauroot) +
     >                     arootl(k) * tauroot
     >                     * max(0., anpplxy(i,j,k)) *
     >                     (1. - exp(-1./tauroot))
            cruxy(i,j,k) = max(0.0, cruxy(i,j,k) *
     >                         (1. - disturbf(i,j) - disturbo(i,j)))

          endif

 106     continue  !npftu

c fractional cover of lower and upper canopy

         fuxy(i,j) = (1.0 - exp(shape * wood / woodnorm )) / fnorm
         flxy(i,j) = (1.0 - exp(shape * totlailxy(i,j) / 1.0)) / fnorm

         fuxy(i,j) = max(0.10, min(0.98,fuxy(i,j)))
         flxy(i,j) = max(0.10, min(0.98,flxy(i,j)))

c apply disturbances

         fuxy(i,j) = fuxy(i,j) * (1. - disturbf(i,j) - disturbo(i,j)) 
         flxy(i,j) = flxy(i,j) * (1. - disturbf(i,j) - disturbo(i,j)) 

c constrain the fractional cover bet 10% and 98%

         fuxy(i,j) = max(0.10, min(0.98,fuxy(i,j)))
         flxy(i,j) = max(0.10, min(0.98,flxy(i,j)))

        endif      !only land points

 104    continue   !nlon
 102   continue    !nlat

c-------------
 100  continue
c-------------

c determine which plant types can exist in a given grid cell
c according to climate-range criteria (and only over land points)

      call zero (existu, nlon*nlat*npftu)
      call zero (existl, nlon*nlat*npftl)

      tplus   = 3.00           ! temp on warm side tlimit when exist = 0
      tshape1 = 5.00           ! controls shape of cold side temp limit
      gshape1 = 0.05           ! controls shape of cold gdd limit
c                              ! larger shape1 => steeper curve
      shape2  = 0.50           ! controls shape of warm side temp limit
c                              ! shape2 = 1.0 => straight line
      do 200 j=1,nlat
        do 202 i=1,nlon

          if (lmaskxy(i,j).eq.1) then

            ztcold = tcoldprevxy(i,j) - 273.16
            ztwarm = twarmprevxy(i,j) - 273.16

c tropical evergreen woody plants

c           if (ztcold.gt.15.5)             existu(i,j,1) = 1.
        existu(i,j,1) = min (1., exp(tshape1*(ztcold - 15.5)))

c tropical raingreen woody plants

c           if (ztcold.gt.15.5)             existu(i,j,2) = 1.
        existu(i,j,2) = min (1., exp(tshape1*(ztcold - 15.5)))

c warm-temperate evergreen woody plants

c           if ((ztcold.gt.5.0).and.
c    >          (ztcold.lt.15.5))           existu(i,j,3) = 1.
       existu(i,j,3) = min (1., exp(tshape1*(ztcold - 5.0))) *
     > min (1.,
     > (min (0., ztcold - 15.5 - tplus) / (-tplus))**shape2)

c cool-temperate evergreen woody plants

c           if ((ztcold.gt.-19.0).and.
c    >          (ztcold.lt.5.0).and.
c    >          (gdd5prev(i,j).gt.1200.0))  existu(i,j,4) = 1.
      existu(i,j,4) = min (1., exp(tshape1*(ztcold + 19.0))) *
     >   min ( 1.,
     >  (min (0., ztcold - 5. - tplus) / (-tplus))**shape2) *
     >   min (1., exp(gshape1*(gdd5prev(i,j) - 1200.0)))

c temperate summergreen woody plants

c           if ((ztcold.gt.-15.0).and.
c    >          (ztcold.lt.15.5).and.
c    >          (gdd5prev(i,j).gt.1200.0))  existu(i,j,5) = 1.
      existu(i,j,5) = min (1., exp(tshape1*(ztcold + 15.0))) *
     >   min ( 1.,
     >  (min (0., ztcold - 15.5 - tplus) / (-tplus))**shape2) *
     >   min (1., exp(gshape1*(gdd5prev(i,j) - 1200.0)))

c boreal evergreen woody plants

c           if ((ztcold.gt.-35.0).and.
c    >          (ztcold.lt.-2.0).and.
c    >          (gdd5prev(i,j).gt.350.0))   existu(i,j,6) = 1.
      existu(i,j,6) = min (1., exp(tshape1*(ztcold + 35.0))) *
     >   min ( 1.,
     >  (min (0., ztcold + 2. - tplus) / (-tplus))**shape2) *
     >   min (1., exp(gshape1*(gdd5prev(i,j) - 350.0)))

c boreal summergreen woody plants

c           if ((ztcold.lt.5.0).and.
c    >          (gdd5prev(i,j).gt.350.0))   existu(i,j,7) = 1.
       existu(i,j,7) = min ( 1.,
     >    (min (0., ztcold - 5. - tplus) / (-tplus))**shape2) *
     >     min (1., exp(gshape1*(gdd5prev(i,j) - 350.0)))

c c3 herbaceous plants

c           if (gdd0prev(i,j).gt.100.0)     existl(i,j,1) = 1.
      existl(i,j,1) = min (1.,
     >   exp(gshape1*(gdd0prev(i,j) - 100.0)))

c c4 herbaceous plants

c           if ((gdd0prev(i,j).gt.100.0).and.
c    >          (ztwarm.gt.22.0))           existl(i,j,2) = 1.
       existl(i,j,2) = min (1., exp(tshape1*(ztwarm - 22.))) *
     >    min (1., exp(gshape1*(gdd0prev(i,j) - 100.0)))

          endif !do this only on land points

  202   continue
  200 continue

c~~~~~~~~~~
      endif
c~~~~~~~~~~

c>>>>>>>>>>
      endif
c>>>>>>>>>>

c *********************************************
c Impose climate-range limits on carbon stores,
c Set plai* from cl*,
c Set saixy,
c Set various temporary sums over plant types
c *********************************************

      call zero (saixy,      nlon*nlat*2)
      call zero (effhtuxy,   nlon*nlat)
      call zero (effhtlxy,   nlon*nlat)
      call zero (woodmasstu, nlon*nlat)
      call zero (totlaiuxy, nlon*nlat)
      call zero (totlailxy, nlon*nlat)
      call zero (totbiouxy, nlon*nlat)
      call zero (totbiolxy, nlon*nlat)

      do 300 n=1,npftu
        do 302 j=1,nlat
          do 304 i=1,nlon

#ifdef sgi
           if (lmaskxy(i,j).eq.1) then
#endif

             if (iyr.gt.kyear) then ! Begin veg dynamics after yr kyear
               cluxy(i,j,n) = existu(i,j,n) * cluxy(i,j,n)
               cwuxy(i,j,n) = existu(i,j,n) * cwuxy(i,j,n)
               cruxy(i,j,n) = existu(i,j,n) * cruxy(i,j,n)
             endif

c              Small initial amounts of a plant type are created here
c              for grid points where previously there was none 

            plaiuxy(i,j,n) = max ( .05 * existu(i,j,n),
     >                             cluxy(i,j,n)*speclau(n) )
            cluxy(i,j,n) = plaiuxy(i,j,n) / speclau(n)

            saixy(i,j,2) = saixy(i,j,2) + 0.1 * plaiuxy(i,j,n)

            biomassuxy(i,j,n)=cluxy(i,j,n) + cwuxy(i,j,n) + cruxy(i,j,n)

            effhtuxy(i,j) = effhtuxy(i,j) + 2.50*biomassuxy(i,j,n)*0.75

            woodmasstu(i,j) = woodmasstu(i,j) + cwuxy(i,j,n)

            totlaiuxy(i,j) = totlaiuxy(i,j) + plaiuxy(i,j,n)

            totbiouxy(i,j) = totbiouxy(i,j) + biomassuxy(i,j,n)

#ifdef sgi
           endif !only on land points
#endif

  304     continue
  302   continue
  300 continue

      do 310 n=1,npftl
        do 312 j=1,nlat
          do 314 i=1,nlon

#ifdef sgi
           if (lmaskxy(i,j).eq.1) then
#endif

             if (iyr.gt.kyear) then ! Begin veg dynamics after yr kyear
               cllxy(i,j,n) = existl(i,j,n) * cllxy(i,j,n)
               crlxy(i,j,n) = existl(i,j,n) * crlxy(i,j,n)
             endif

c              Small initial amounts of a plant type are created here
c              for grid points where previously there was none 

            plailxy(i,j,n) = max ( .05 * existl(i,j,n),
     >                             cllxy(i,j,n)*speclal(n) )
            cllxy(i,j,n) = plailxy(i,j,n) / speclal(n)

            saixy(i,j,1) = saixy(i,j,1) + 0.01 * plailxy(i,j,n)

            biomasslxy(i,j,n) = cllxy(i,j,n) + crlxy(i,j,n)

            effhtlxy(i,j) = effhtlxy(i,j) + plailxy(i,j,n) * 0.20

            totlailxy(i,j) = totlailxy(i,j) + plailxy(i,j,n)

            totbiolxy(i,j) = totbiolxy(i,j) + biomasslxy(i,j,n)

#ifdef sgi
           endif !only on land pts
#endif

  314     continue
  312   continue
  310 continue

c *****************************************************
c calculate various aggregate quantites for each canopy
c *****************************************************

c set grass and tree heights

      do 400 j=1,nlat
        do 402 i=1,nlon

#ifdef sgi
         if (lmaskxy(i,j).eq.1) then
#endif

          zbotxy(i,j,1) =  0.01
          ztopxy(i,j,1) =  max (0.10, effhtlxy(i,j))

          zbotxy(i,j,2) =  max (1.00, ztopxy(i,j,1) + 0.50)
          ztopxy(i,j,2) =  max (zbotxy(i,j,2)+1.0, effhtuxy(i,j))
c
c calculate sapwood biomass
c
          sapspeed  = 25.0                                  ! (m/day)
c                                                           ! / m**2 sa
          trans     = 0.003                                 ! (m/day) 
          saparea   = trans / sapspeed                      ! m**2
c         sapvolume = saparea * effhtuxy(i,j) * 0.75        ! m**3
          sapvolume = saparea * ztopxy(i,j,2) * 0.75        ! m**3
          denswood  = 400.0                                 ! kg/m**3
c
c determine the sapwood fraction (allowed to range from 0.1 to 0.9)
c
          sapfrac(i,j) = min ( 0.9, max ( 0.1,
     >                   sapvolume*denswood/(woodmasstu(i,j)+epsilon)
     >                   ) )

c only before veg. dynamics begin:
c if there is a continuous canopy cover (one-sided lai > 1)
c then assume the fractional cover is 1.0; else (one-sided lai < 1),
c reduce the fractional cover accordingly.
 
          if (iyr.le.kyear .or. iflag.eq.0) then

            flxy(i,j) = 1.0

            if (totlaiuxy(i,j).ge.0.98) then
              fuxy(i,j) = .98
            else
              fuxy(i,j) = max(0.1,totlaiuxy(i,j) / 1.0)
            endif

          endif

#ifdef sgi
         endif !only on land pts
#endif

  402   continue
  400 continue
 
c        Write 2-D maps to IBIS ascii output file 
c        (just for call at end of each year)
 
c<<<<<<<<<<<<<<<<<<<<<<<<<
      if (iflag.eq.1) then
c<<<<<<<<<<<<<<<<<<<<<<<<<

        call mapibis (totlaiuxy, 100., 'totlaiuxy (m2/m2)', cdate)
        call mapibis (totlailxy, 100., 'totlailxy (m2/m2)', cdate)

        do 520 n=1,npftu
          write (clabel, '(a,i1,2x,a)') 'plaiuxy_', n, '(m2/m2)'
          call mapibis (plaiuxy(1,1,n), 100., clabel, cdate)
  520   continue

        do 525 n=1,npftl
          write (clabel, '(a,i1,2x,a)') 'plailxy_', n, '(m2/m2)'
          call mapibis (plailxy(1,1,n), 100., clabel, cdate)
  525   continue

        call mapibis (totbiouxy, 100., 'totbiouxy (m2/m2)', cdate)
        call mapibis (totbiolxy, 100., 'totbiolxy (m2/m2)', cdate)

        do 540 n=1,npftu
          write (clabel, '(a,i1,2x,a)') 'biomassu_', n, '(kg_C/m2)'
          call mapibis (biomassuxy(1,1,n), 100., clabel, cdate)
  540   continue

        do 545 n=1,npftl
          write (clabel, '(a,i1,2x,a)') 'biomassl_', n, '(kg_C/m2)'
          call mapibis (biomasslxy(1,1,n), 100., clabel, cdate)
  545   continue

c---

        call mapibis (sapfrac,  100., 'sapfrac (fraction)', cdate)
        call mapibis (effhtuxy, 100., 'effhtuxy (m)',       cdate)
        call mapibis (effhtlxy, 100., 'effhtlxy (m)',       cdate)

c---

        call mapibis (agppxy, 100., 'agppxy (kg_C/m2/yr)', cdate)
        call mapibis (anppxy, 100., 'anppxy (kg_C/m2/yr)', cdate)

        do 560 n=1,npftu
          write (clabel, '(a,i1,1x,a)') 'anppuxy_', n, '(kg_C/m2/yr)'
          call mapibis (anppuxy(1,1,n), 100., clabel, cdate)
  560   continue

        do 565 n=1,npftl
          write (clabel, '(a,i1,1x,a)') 'anpplxy_', n, '(kg_C/m2/yr)'
          call mapibis (anpplxy(1,1,n), 100., clabel, cdate)
  565   continue

c---

c       call mapibis (gtransmaxuxy,100.,'gtransmaxuxy (mm/day)',  cdate)
        call mapibis (tcoldxy,   -tmelt,'tcoldxy (K)',            cdate)
        call mapibis (twarmxy,   -tmelt,'twarmxy (K)',            cdate)
        call mapibis (tdailxy,   -tmelt,'tdailxy (K)',            cdate)
        call mapibis (gdd0,          1.,'gdd0 (C*day)',           cdate)
        call mapibis (gdd5,          1.,'gdd5 (C*day)',           cdate)
        call mapibis (avgfrespaxy, 100.,'avgfrespaxy (non-dim)',  cdate)
        call mapibis (avgfrespbxy, 100.,'avgfrespbxy (non-dim)',  cdate)

c write fort.96 to mass storage

#ifdef gcm
        call disposef('fort.96',mspath)
#endif

c<<<<<<<<<<
      endif
c<<<<<<<<<<

c       Set ivegxy (history-file vegetation types, diagnostic only)

      call vegmap (totlaiuxy, totlailxy)

      return
      end
c=====
#endif
c=====
c
c ---------------------------------------------------------------------
c
c=============
#ifdef ibisint
c=============
 
      subroutine mapibis (arr, scale, clabel, cdate)
 
c        Writes one 2-D map to an ascii file for IBIS dynamical veg
c        output (called from dynaveg).
c        (Format statements have a hard-coded limit of nlon <= 180).
c
c        arr   = 2-D array to be displayed
c        scale = if > 0, scaling  factor for arr (output is arr*scale)
c                if <=0, additive shift  for arr (output is arr+scale)
c        clabel = field label
c        cdate  = date label
 
#include <compar.h>
#include <comgrd.h>
#include <commapl.h>
#include <commugxy.h>
 
      dimension arr(nlon,nlat)
      character*(*) clabel, cdate
 
      character carr(nlon)*5, cyr*20
      parameter (dpr=(180./pi))
      data iuout /96/
 
 
c        Write header line
 
      write(iuout,100)
     *  clabel(1:lenchr(clabel)), scale,  cdate(1:lenchr(cdate))
  100 format(/a,'   scale=',f12.6,'   date=',a)
 
c        Write longitude grid values (-180 to 180 deg E)
 
      write(iuout,200) (nint(alon(i)*dpr-360.), i=nlon/2+1,nlon),
     *                 (nint(alon(i)*dpr),      i=1,nlon/2)
  200 format(/1x,5x,3x,180i5)
 
c        Write 2-D data, with blanks for ocean/icesheet
 
      do 300 j=nlat,1,-1
        do 302 i=1,nlon
          if (lmaskxy(i,j).eq.1) then
            if (scale.gt.0.) then
              write (carr(i),'(i5)') nint(arr(i,j)*scale)
            else
              write (carr(i),'(i5)') nint(arr(i,j)+scale)
            endif
          else
            carr(i) = ' '
          endif
  302   continue
        write(iuout,304) alat(j)*dpr, (carr(i),i=nlon/2+1,nlon),
     *                                (carr(i),i=1,nlon/2)
  304   format(1x,f5.1,3x,180a5)
  300 continue
 
c        flush file buffer
 
      ier = flush(iuout)
 
      return
      end
 
c=====
#endif
c=====
c
c ---------------------------------------------------------------------
c
c=============
#ifdef ibisint
c=============

      subroutine vegmap (totlaiuxy, totlailxy)

c Assigns vegetation classes (types) to each grid point, based
c on PFTs and their LAIs. Adapted from ibis2.1. Sets ivegxy
c (in comvegxy) to these types, for LSX history file.  

#include <compar.h>
#include <comgrd.h>
#include <commugxy.h>
#include <comvegxy.h>

      dimension  totlaiuxy(nlon,nlat), totlailxy(nlon,nlat)
      real maxlai
      integer domtree

c ---------------
c PFTs (ibis2.0):
c u1: tropical evergreen woody plants
c u2: tropical raingreen woody plants
c u3: warm-temperate evergreen woody plants
c u4: cool-temperate evergreen woody plants
c u5: temperate summergreen woody plants
c u6: boreal evergreen woody plants
c u7: boreal summergreen woody plants
c l1: c3 herbaceous plants
c l2: c4 herbaceous plants

c standard ibis vegetation classes (from ibis2.1):
c  1: tropical evergreen forest/woodland
c  2: tropical deciduous forest/woodland
c  3: temperate evergreen broadleaf forest/woodland
c  4: temperate evergreen conifer forest/woodland
c  5: temperate deciduous forest/woodland
c  6: boreal evergreen forest/woodland
c  7: boreal deciduous forest/woodland
c  8: mixed forest/woodland
c  9: savanna
c 10: grassland/steppe 
c 11: dense shrubland
c 12: open shrubland
c 13: tundra
c 14: desert 
c 15: polar desert/rock/ice
c ---------------

c================
      do j=1,nlat
      do i=1,nlon
c================

        ivegxy(i,j) = 0

c----------------------------------
        if (lmaskxy(i,j).eq.1) then
c----------------------------------

c determine total lai and tree, shrub, and grass fractions

        treelai   = totlaiuxy(i,j) 
        shrublai  = -1.                           ! no shrubs in ibis2.0
        grasslai  = plailxy(i,j,1) + plailxy(i,j,2)
        totlai    = max (0.01, totlailxy(i,j) + totlaiuxy(i,j))

c determine dominant tree type by lai dominance

        domtree = 0
        maxlai = 0.0
        do k=1,npftu
          if (plaiuxy(i,j,k).gt.maxlai) then
            domtree = k
            maxlai = plaiuxy(i,j,k)
          endif
        enddo

c dominant type:  tropical broadleaf evergreen tree

        if (domtree.eq.1) then

          if (treelai.gt.2.5) ivegxy(i,j) =  1   ! trop.ever.forest/wood
          if (treelai.le.2.5) ivegxy(i,j) =  9   ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
            if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
          endif

c dominant type:  tropical broadleaf drought-deciduous tree

        else  if (domtree.eq.2) then

          if (treelai.gt.2.5) ivegxy(i,j) =  2    ! trop.dec.forest/wood
          if (treelai.le.2.5) ivegxy(i,j) =  9    ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
            if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
          endif

c dominant type:  warm-temperate broadleaf evergreen tree

        else if (domtree.eq.3) then

          if (treelai.gt.2.5) ivegxy(i,j) =  3    ! temp.ever.broad.f/w
          if (treelai.le.2.5) ivegxy(i,j) =  9    ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
            if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
          endif

c dominant type:  temperate conifer evergreen tree

        else if (domtree.eq.4) then

          if (treelai.gt.1.5) ivegxy(i,j) =  4    ! temp.ever.conif.f/w
          if (treelai.le.1.5) ivegxy(i,j) =  9    ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
            if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
          endif

c dominant type:  temperate broadleaf deciduous tree

        else if (domtree.eq.5) then

          if (treelai.gt.1.5) ivegxy(i,j) =  5    ! temp.decid.f/w
          if (treelai.le.1.5) ivegxy(i,j) =  9    ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
            if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
          endif

c dominant type:  boreal conifer evergreen tree

        else if (domtree.eq.6) then

                              ivegxy(i,j) =  6    ! boreal ever. f/w
c         if (treelai.gt.1.0) ivegxy(i,j) =  6    ! boreal ever. f/w
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
c           if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
c         endif


c dominant type:  boreal cold-deciduous tree

        else if (domtree.eq.7) then

                              ivegxy(i,j) =  7    ! boreal decid. f/w
c         if (treelai.gt.1.0) ivegxy(i,j) =  7    ! boreal decid. f/w
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
c           if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
c         endif

        endif

c temperate/boreal forest mixtures

        if ((domtree.ge.4).and.(domtree.le.7)) then
          ratio = (plaiuxy(i,j,5) + plaiuxy(i,j,7)) / 
     >            (plaiuxy(i,j,4) + plaiuxy(i,j,5) + plaiuxy(i,j,6) +
     >             plaiuxy(i,j,7))
          if (treelai.gt.1.0) then
            if ((ratio.gt.0.45).and.(ratio.lt.0.55)) ivegxy(i,j) = 7 
          endif
          if ((domtree.le.5).and.(treelai.le.1.0)) then
            if (grasslai.ge.shrublai) ivegxy(i,j) = 10    ! grassland
            if (shrublai.ge.grasslai) ivegxy(i,j) = 11    ! closed shrub
          endif
        endif

c no tree is dominant

        if (domtree.eq.0) then
          if (treelai.gt.1.0) ivegxy(i,j) =  9    ! savanna
          if (treelai.le.1.0) ivegxy(i,j) = 10    ! grassland
        endif

c overriding vegtation classifications

        if (totlai.lt.1.0)    ivegxy(i,j) = 12    ! open shrubland
        if (totlai.le.0.4)    ivegxy(i,j) = 14    ! desert

c overriding climatic rules

        if (gdd5prev(i,j).lt.350.0) then
          if (totlai.ge.0.4)    ivegxy(i,j) = 13    ! tundra
          if (totlai.lt.0.4)    ivegxy(i,j) = 15    ! polar desert
        endif

        if (gdd0prev(i,j).lt.100.0) ivegxy(i,j) = 15    ! polar desert

c------------
        endif
c------------

c==========
      enddo
      enddo
c==========

      return
      end
c=====
#endif
c=====
c
c ---------------------------------------------------------------------
c
c Stubs for compiling with -D ibis but not -D ibisint
 
c==============
#ifndef ibisint
c==============
      subroutine ibis_write ()
      write(*,*) '*** Error: called stub ibis_write'
      call endrun (-1)
      end
      subroutine ibis_read ()
      write(*,*) '*** Error: called stub ibis_read'
      call endrun (-1)
      end
      subroutine carbon ()
      write(*,*) '*** Error: called stub carbon'
      call endrun (-1)
      end
      subroutine dynaveg ()
      write(*,*) '*** Error: called stub dynaveg'
      call endrun (-1)
      end
c=====
#endif
c=====
 
c End of "ifdef ibis" for the whole file. (Stubs are on end of vegdat.f
c for compiling without -D ibis and without this file, to satisfy
c links from within Genesis.
 
 
c//////////
#else
c//////////
 
c     Just to satisfy linker (non-empty ibis.o)
      subroutine ibis_dum
      return
      end
 
c/////
#endif
c/////
