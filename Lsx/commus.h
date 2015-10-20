c -- $Id: commus.h,v 1.2 2001/07/03 21:29:25 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Lsx/commus.h,v $
c -- $Name:  $

c-------
c commus
c-------
c
c     1D strip arrays for communication between surface models
c     (with no corresponding global arrays)
c
c----------------------------------------------------------------------
c
c        HEATI  = net heat flux into snow surface
c        HEATG  = net heat flux into soil surface
c        UPSOIU = H2O uptake from each soil layer due to transpiration
c                 in upper story
c        UPSOIL = H2O uptake from each soil layer due to transpiration
c                 in lower story
c 
c----
c        For BATS-SiB stomatal method:
c        UPLAYU = max upper-story water uptake rate, each soil layer
c        UPLAYL = max lower-story water uptake rate, each soil layer
c        For IBIS stomatal method:
c        UPLAYU = upper-story soil-moisture stress term, each soil layer
c        UPLAYL = lower-story soil-moisture stress term, each soil layer
c----
c        UPTOT[U,L] = sum of UPLAY[U,L] over soil layers 
c        UPSIZ[U,L] = (LAI+SAI)/BLAIM for BATS/SiB stom.res method only
c        QGLIF  = fraction of soil evap (fvapg) from soil liquid
c        fudgbud[r,s] = adjusted fvapg->raing/snowg if condensation 
c        HFLO   = dwnwd heat  transport thru soil/landice layer bdries
c        FPEN   = dwnwd penetrat. solar flux thru landice layer bdries
c        CONSOI = thermal conductivity of each soil layer
c        HVASUG = latent heat of vap/subl, for soil or ice surface
c        CONSNO = thermal conductivity of each snow layer
c        HVASUI = latent heat of vap/subl, for snow surface
c        HFLOICE= dwnwd heat  transport thru seaice layer bdries
c        FPENICE= dwnwd penetrat. solar flux thru seaice layer bdries
c        SITOP  = rate of melt of top layer of seaice (mass flux)
c        SIBASE = rate of melt(+)/freeze(-) at seaice base (mass flux)
c
      COMMON  /SURLSX/
c
c          lsx to snow,soil (HEATG is areally averaged by SNOW for 
c          the soil model.)
c
     *  HEATI(NVEC),              HEATG(NVEC),
     *  UPSOIU(NVEC,NSOILAY),     UPSOIL(NVEC,NSOILAY), 
c
c          soil/icesheet to lsx and snow (set in SETSOI)
c
     *  UPLAYU(NVEC,NSOILAY),     UPLAYL(NVEC,NSOILAY),
     *  UPTOTU(NVEC),             UPTOTL(NVEC),
     *  UPSIZU(NVEC),             UPSIZL(NVEC),
     *  QGLIF(NVEC),              
     *  fudgbudr(nvec),           fudgbuds(nvec),
     *  HFLO(NVEC,NSOILAY+1),     FPEN(NVEC,NSOILAY+1),
     *  CONSOI(NVEC,NSOILAY),     HVASUG(NVEC),
c
c          snow/icesheet to lsx, seaice and snow (set in SETSNO)
c
     *  CONSNO(NVEC,NSNOLAY),     HVASUI(NVEC),
c
c          seaice1 to seaice2 (set in seaice1)
c
     *  HFLOICE(NVEC,NICELAY+1),  FPENICE(NVEC,NICELAY+1),
c
c          seaice to nowhere currently (set and used in SEAICE)
c
     *  SITOP(NVEC), SIBASE(NVEC)
!$OMP THREADPRIVATE ( /SURLSX/ )
c
