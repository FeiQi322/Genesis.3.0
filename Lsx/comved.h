c -- $Id: comved.h,v 1.2 2001/07/03 21:29:27 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comved.h,v $
c -- $Name:  $

c-------
c comved
c-------
c
c     Vegetation database from Dorman and Sellers(1989),JAM,28,833-855.
c
c     D-S specific arrays, initialized by data statements below in
c     subr vegdat_dorm and vegdat_sdvm. Declared as local 
c    (dimension...) and saved, not in common, to allow different data 
c     settings in these subroutines. (And only used in these 2 subrs).
c
c-----
c        NVED       = number of types in D-S vegetation database, +1 for
c                     non-land pts (=13). NVED is in compar.h, to allow
c                     use in histwrit and earlier in vegdat_dorm,sdvm.
c        NTIMDS     = number of time-dependent data points per year
c        NWHEAT     = number of winter-wheat latitude zones(DS Table 5)
c
      PARAMETER (NTIMDS=12, NWHEAT=7)
c----
c
c        TAI[L,U]Z  = total one-sided LAI+SAI for [lower,upper] story
c        GAI[L,U]Z  = green one-sided LAI
c        F[L,U]Z    = fractional areal cover
c
c        TAUV[L,U]Z = VIS leaf transmittance, (*,1)=live, (*,2)=dead
c        RHOV[L,U]Z = VIS leaf reflectance,   (*,1)=live, (*,2)=dead
c        TAUN[L,U]Z = NIR leaf transmittance, (*,1)=live, (*,2)=dead
c        RHON[L,U]Z = NIR leaf reflectance,   (*,1)=live, (*,2)=dead
c
c        CHIF[L,U]Z = "chi" orientation factor for veg. elements
c                     (-1 = all vert, 0 = all random, 1 = all horiz)
c
c        Z1[L,U]Z   = height of vegetation base above ground (m)
c        Z2[L,U]Z   = height of vegetation top  above ground (m)
c
c        APAR[L,U]Z = "a"  PAR coeff for stomatal resistance (J m-3)
c        BPAR[L,U]Z = "b"  PAR coeff for stomatal resistance (W m-2)
c        CPAR[L,U]Z = "c"  PAR coeff for stomatal resistance (s m-1)
c        H5[L,U]Z   = "h5" vpd coeff for stomatal resistance (mb-1)
c
c        TAIWHEAT   = total one-sided LAI+SAI for winter wheat(Table 5)
c        GAIWHEAT   = green one-sided LAI+SAI for winter wheat(Table 5)
c        FRAWHEAT   = fractional areal cover  for winter wheat(Table 5)
c        CHIWHEAT   = chi orientation factor  for winter wheat(Table 5)
c
      DIMENSION
c
     *  TAILZ(NVED,NTIMDS),  GAILZ(NVED,NTIMDS),  FLZ(NVED,NTIMDS),
     *  TAIUZ(NVED,NTIMDS),  GAIUZ(NVED,NTIMDS),  FUZ(NVED,NTIMDS),
c
     *  TAUVLZ(NVED,2),    RHOVLZ(NVED,2), 
     *  TAUVUZ(NVED,2),    RHOVUZ(NVED,2), 
     *  TAUNLZ(NVED,2),    RHONLZ(NVED,2), 
     *  TAUNUZ(NVED,2),    RHONUZ(NVED,2), 
c
     *  CHIFLZ(NVED),      CHIFUZ(NVED),
c
     *  Z1LZ(NVED),        Z2LZ(NVED),
     *  Z1UZ(NVED),        Z2UZ(NVED),
c
     *  APARLZ(NVED),      BPARLZ(NVED),      CPARLZ(NVED),
     *  APARUZ(NVED),      BPARUZ(NVED),      CPARUZ(NVED),
     *  H5LZ(NVED),        H5UZ(NVED),
c
     *  TAIWHEAT(NWHEAT,NTIMDS),  GAIWHEAT(NWHEAT,NTIMDS), 
     *  FRAWHEAT(NWHEAT,NTIMDS),  CHIWHEAT(NWHEAT,NTIMDS)
c
      save 
     *  tailz,  gailz,  flz,    taiuz,  gaiuz,    fuz, 
     *  tauvlz, tauvuz, taunlz, taunuz, rhovlz,   rhovuz, 
     *  rhonlz, rhonuz, chiflz, chifuz, z1lz,     z2lz,
     *  z1uz,   z2uz,   aparlz, bparlz, cparlz,   aparuz,
     *  bparuz, cparuz, h5lz,   h5uz,   taiwheat, gaiwheat,
     *  frawheat, chiwheat
c
