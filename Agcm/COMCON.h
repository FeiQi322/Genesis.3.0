c -- $Id: COMCON.h,v 1.2 2001/07/03 21:29:05 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMCON.h,v $
c -- $Name:  $

C
CL            PHYSICAL CONSTANTS
C
      COMMON /COMCON/
     +   CALUNI,   CAPPA,     CARMAN,   CPVIR,    CPWV,
     +   CPAIR,    DAYLENCOR, EZ,       OMEGA,    GRAVIT,    
     *   DAYLENRAD,SSLON,     LATVAP,   LATICE,   LATSUB,   
     +   RAIR,     REARTH,    RH2O,     RHOH2O,   RA,        RGA,
     +   STEBOL,   TMELT,     ZVIR,
     +   CH2O,     CICE
#ifdef isotrac
c    ratios of isotopic species and pure water diffusion coefficients
c    in the vapor phase / VSMOW standard ratios (time 2*19/18 for 2H
c    and time 20/18 for 18O)
     +  ,dr_o,     dr_h,     smoj_o,   smoj_h
#endif
      REAL  LATVAP,LATICE,LATSUB

