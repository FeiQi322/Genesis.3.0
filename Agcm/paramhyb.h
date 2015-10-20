c -- $Id: paramhyb.h,v 1.2 2001/07/03 21:29:19 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/paramhyb.h,v $
c -- $Name:  $

c paramhyb

c--------------------------------------------------------------
c Parameters and statement functions for slt hybrid coordinates
c--------------------------------------------------------------

c p0      = hybrid reference pressure (N m-2)
c nps     = number of fine-spaced surface-pressure (ps) values 
c           for second index of sigval,etc
c psa,psb = surface pressure range (N m-2) for index 1 to nps

c sigval  = sigma value for mid-point eta-grid and ps index
c sigvall = sigma value for interface eta-grid and ps index
c etaval  = eta value   for mid-point sigma-grid and ps index
c detadp  = d(eta)/d(p)  at const ps for given mid-pt sigma and ps index
c detadps = d(eta)/d(ps) at const p  for given mid-pt sigma and ps index

      parameter (p0=101300.)   ! Genesis max(ps)    ~= 106600
                               ! Simmons+Burridge p0 = 101300

      parameter (nps = 1001, psa = 30000., psb = 130000.)

      common /chybset/ 
     *  sigval(nlev,nps),  sigvall(nlev+1,nps), 
     *  etaval(nlev,nps),  detadp(nlev,nps),     detadps(nlev,nps)

