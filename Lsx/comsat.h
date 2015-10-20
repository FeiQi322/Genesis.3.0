c -- $Id: comsat.h,v 1.2 2001/07/03 21:29:26 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comsat.h,v $
c -- $Name:  $

c-------
c comsat 
c-------

c    If tablesvp defined, use functions and table in inisurf.F
c    instead of statement functions and Lowe polynomials as here

#ifndef tablesvp

c     Statement functions and associated parameters.
c
c----------------------------------------------------------------------
c
c        Polynomials for svp(t), d(svp)/dt over water and ice are from
c        Lowe(1977),JAM,16,101-103.
c
      PARAMETER (ASAT0 =  6.107799961,
     *           ASAT1 =  4.436518521E-1,
     *           ASAT2 =  1.428945805E-2,
     *           ASAT3 =  2.650648471E-4,
     *           ASAT4 =  3.031240396E-6,
     *           ASAT5 =  2.034080948E-8,
     *           ASAT6 =  6.136820929E-11 )
c
      PARAMETER (BSAT0 =  6.109177956,
     *           BSAT1 =  5.03469897E-1,
     *           BSAT2 =  1.886013408E-2,
     *           BSAT3 =  4.176223716E-4,
     *           BSAT4 =  5.824720280E-6,
     *           BSAT5 =  4.838803174E-8,
     *           BSAT6 =  1.838826904E-10 )
c
      PARAMETER (CSAT0 =  4.438099984E-1,
     *           CSAT1 =  2.857002636E-2,
     *           CSAT2 =  7.938054040E-4,
     *           CSAT3 =  1.215215065E-5,
     *           CSAT4 =  1.036561403E-7,
     *           CSAT5 =  3.532421810E-10,
     *           CSAT6 = -7.090244804E-13 )
c
      PARAMETER (DSAT0 =  5.030305237E-1,
     *           DSAT1 =  3.773255020E-2,
     *           DSAT2 =  1.267995369E-3,
     *           DSAT3 =  2.477563108E-5,
     *           DSAT4 =  3.005693132E-7,
     *           DSAT5 =  2.158542548E-9,
     *           DSAT6 =  7.131097725E-12 )
c
c----------
c
c         Statement functions TSATL,TSATI are used below so that Lowe's
c         polyomial for liquid is used if T gt TMELT, or for ice if 
c         T lt TMELT. Also impose range of validity for Lowe's polys.
c
       TSATL(T) = MIN (100., MAX (T-TMELT, 0.))
       TSATI(T) = MAX (-55., MIN (T-TMELT, 0.))
c
c----------
c
c         Statement function ESAT is svp in N/m**2, with T in deg K. 
c         (100 * Lowe's poly since 1 mb = 100 N/m**2.)
c
      ESAT (T) = 
     *  100.*(
     *    CVMGT (ASAT0, BSAT0, T.GE.TMELT)
C
     *    + TSATL(T)*(ASAT1 + TSATL(T)*(ASAT2 + TSATL(T)*(ASAT3
     *    + TSATL(T)*(ASAT4 + TSATL(T)*(ASAT5 + TSATL(T)* ASAT6)))))
C
     *    + TSATI(T)*(BSAT1 + TSATI(T)*(BSAT2 + TSATI(T)*(BSAT3
     *    + TSATI(T)*(BSAT4 + TSATI(T)*(BSAT5 + TSATI(T)* BSAT6)))))
     *  )
c
c----------
c
c         Statement function DSAT is d(sat. spec. humidity)/dt, with T 
c         in deg K, and neglecting q in denominator of the q(vp) reln.
c         (100 * Lowe's poly since 1 mb = 100 N/m**2.)
c
      DSAT (T, P) = 
     *  100.*.622*(
     *    CVMGT (CSAT0, DSAT0, T.GE.TMELT)
C
     *    + TSATL(T)*(CSAT1 + TSATL(T)*(CSAT2 + TSATL(T)*(CSAT3
     *    + TSATL(T)*(CSAT4 + TSATL(T)*(CSAT5 + TSATL(T)* CSAT6)))))
C
     *    + TSATI(T)*(DSAT1 + TSATI(T)*(DSAT2 + TSATI(T)*(DSAT3
     *    + TSATI(T)*(DSAT4 + TSATI(T)*(DSAT5 + TSATI(T)* DSAT6)))))
     *  ) / P
c
c----------
c
c        Statement function QSAT is saturation specific humidity,
c        with svp E and ambient pressure P in N/m**2. Impose an upper
c        limit of 1 to avoid spurious values for very high svp
c        and/or small P.
c
       QSAT (E, P) = .622*E / MAX ( P-(1.-.622)*E, .622*E )
c
c----------
c
c     p = 98222.
c     tmelt = 273.16
c     do 10 i=1,501
c       tc= -100. + (i-1) 
c       t = tc + tmelt
c       e = esat(t)
c       d = dsat(t,p)
c       q = qsat(e,p)
c       if (i.eq.1) write(*,900)
c 900   format('    temp',
c    *         '          esat',  
c    *         '    d(qsat)/dT',  
c    *         '          qsat',
c    *         '  [10*dq/dT]/qsat'
c    *         '   T*[dq/dT]/qsat'
c    *         /)
c       write(*,902) tc, e, d, q, d*10./q, d*t/q
c 902   format(f8.2, f14.4, f14.8, f14.5, f17.3, f17.3)
c  10 continue
c     stop
c     end
c
c     function cvmgt (x,y,l)
c     logical l
c     if (l) then
c       cvmgt = x
c     else
c       cvmgt = y
c     endif
c     return
c     end

#endif

c
c----------------------------------------------------------------------
c
c        Statement functions HVAPF, HSUBF correct the latent heats of
c        vaporization (liquid-vapor) and sublimation (ice-vapor) to
c        allow for the concept that the phase change takes place at
c        TMELT, and the various phases are cooled/heated to that 
c        temperature before/after the change. This concept is not
c        physical but is needed to balance the "black-box" energy 
c        budget. Similar correction is applied in CONVAD in the agcm
c        for precip. Needs common comgrd for the physical constants.
c        Argument T is the temp of the liquid or ice, and TAIR is the
c        temp of the delivered or received vapor.
c
      HVAPF(T,TAIR) = HVAP + CVAP*(TAIR-TMELT) - CH2O*(T-TMELT)
      HSUBF(T,TAIR) = HSUB + CVAP*(TAIR-TMELT) - CICE*(T-TMELT)
c
c ---------------------------------------------------------------------
c
c        Statement function VMTEMP is the dependence of Rubisco Capacity
c        (Vmax) on temperature (deg C) for IBIS model photosynthesis 
c
      vmtemp(t) = 1./( (1.+exp(0.35*(9.-t))) * (1.+exp(0.35*(t-41.))) )
c
c----------------------------------------------------------------------

