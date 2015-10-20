c -- $Id: commud.h,v 1.2 2001/07/03 21:29:24 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/commud.h,v $
c -- $Name:  $

c-------
c commud
c-------
c
c     1D arrays for various history-file and budget-related diagnostic
c     fields,  corresponding to global arrays in commudxy (except gevap)
c
c----------------------------------------------------------------------
c
c
c        GINVAP = total evaporation rate from all intercepted H2O
c        GSUVAP = total evaporation rate from sfc (snow/soil/ice/ocean)
c        GTRANSU= total transpiration rate from upper veg
c        GTRANSL= total transpiration rate from lower veg
c        GRUNOF = surface runoff rate
c        GDRAIN = drainage rate out of bottom of lowest soil layer
c        GDRAINH= drainage rate into top of lowest soil layer(hydmodel)
c        WFLO    =downward h2o flow across layer boundaries (kg/s/m2)
c        GABLAT = net icesheet ablation rate
c        GADJUST= h2o flux due to adjustments in subr wadjust
c        TS2    = 2-m surface-air temperature
c        Z0LEN  = roughness length felt by agcm
c
c        GINFIL = liquid h2o flux into top soil layer (infiltration)
c        GEVAP  = evaporation rate from soil surface
c        GIRRIG = irrigation reservoir (adds to rainfall) 
      COMMON  /DIALSX/
     *  GINVAP(NVEC),      GSUVAP(NVEC), 
     *  GTRANSU(NVEC),     GTRANSL(NVEC),
     *  GRUNOF(NVEC),      GDRAIN(NVEC),      GDRAINH(NVEC),
     *  WFLO(NVEC, NSOILAY+1),
     *  GABLAT(NVEC),
     *  GADJUST(NVEC),     TS2(NVEC),         Z0LEN(NVEC),
     *  GINFIL(NVEC),      GEVAP(NVEC),
     *  GIRRIG(NVEC)

#if defined (rcm)
     * ,ts10(nvec),        qs10(nvec),        us10(nvec)
#endif

!$OMP THREADPRIVATE ( /DIALSX/ )
c
