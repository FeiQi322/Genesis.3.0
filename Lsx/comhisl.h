c -- $Id: comhisl.h,v 1.2 2001/07/03 21:29:23 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/comhisl.h,v $
c -- $Name:  $

c--------
c comhisl
c--------
c
c    For lsx (*_a) history files.
c    Similar to ccm common comhis.
c
c-----------------------------------------------------------------------
c
c        verslsx     = lsx model version number
c        cindex(1,*) = program lsx field labels (1-8 chars)
c        cindex(2,*) = units for each program lsx field (1-8 chars)
c        nindex      = number of hard-coded program lsx  fields
c
      common /comhisl/
     *   verslsx,               cindex(2,nushmx),    nindex
      character verslsx*8, cindex*8
c
c-----------------------------------------------------------------------
c
c        vershis_a   = lsx history-file version number
c        cushnam_a   = namelist input for cushis_a and mushfl_a
c        cushis_a(*) = list of user-requested lsx hist fields(1-8 chars)
c        mushis_a(*) = position in cindex for each user field
c        mushfl_a(*) = 0 for instantaneous, 1 for accumulated 
c        nushis_a    = number of user-requested fields
c        mindex_a(1,*) = spare
c        mindex_a(2,*) = number of levels
c        mindex_a(3,*) = 1st  relative record # in "lsx" storage area 
c        mindex_a(4,*) = last relative record # in "lsx" storage area 
c        mindex_a(5,*) = 0, or 1 if a daily min/max field
c        mindex_a(6,*) = 1st relative record # in "min/max" storage area
c        mindfl_a(*)   = mushfl for this program field
c        mindco_a(*)   = counter for accumulated fields
c
c        nulhdr_a = number of header records for each lsx hist write
c        nuldat_a = number of "requested" data records for each lsx 
c                   history write
c        nuldat_atot = as nuldat_a but including extra pairs of recs
c                      on the end for daily min/max calculations
c
c        nhisi_a  = namelist input for nhist
c        ndhis_a  = time interval for each lsx hist accumulation period
c        histlast_a = .true. to force history write at last timestep
c        nhist_a  = start,end times for each lsx accum period (per yr)
c        nhistn_a = number of lsx acumulation intervals (per year)
c        fwhist_a = .true. if lsx history write occurs this timestep
c        fahist_a = .true. if within an lsx accum period this timestep
c        frange   = .true. if current day contains fahist_a or fwhist_a
c        frange24 = .true. if (end of) timestep is at 24:00 GMT 
c
c        nsteph_a = first timestep in current lsx accumulation period
c
      common /comhis_a/
c
     *   vershis_a,           typehis_a,
c
     *   cushnam_a(nushmx),   cushis_a(nushmx),
     *   mushis_a(nushmx),    mushfl_a(nushmx),  nushis_a,
     *   mindex_a(6,nushmx),  mindfl_a(nushmx),  mindco_a(nushmx),
c
     *   nulhdr_a,            nuldat_a,          nuldat_atot,
c
     *   nhisi_a(nmark),      ndhis_a,           histlast_a,
     *   nhist_a(nmark,2),    nhistn_a,
     *   fwhist_a,            fahist_a, 
     *   frange,              frange24,
     *   nsteph_a
c
      character  typehis_a*8, cushnam_a*8, cushis_a*8
      logical    histlast_a,   
     *           fwhist_a,    fahist_a,
     *           frange,      frange24
c
c-----------------------------------------------------------------------
c
