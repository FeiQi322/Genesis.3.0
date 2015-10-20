c--------
c comrich 
c--------

c----------------------------------------------------------------------

c     Statement function giving weight for Prandtl logarithmic
c     profile vs. linear profile in vertical interpolation
c     for TS2 (between lowest gcm level TA and surface-layer temps),
c     used in subrs tscreen (land), ocnflx (ocean/sea ice) and
c     sheetflx (ice sheet). 

c     For usual Ri, weiprand=1. For very stable stratif (Ri > ~ 0.5),
c     "distance to boundary" scaling for vertical mixing is not valid
c     (M.O.length < z), so better to use linear interp.  

      weiprand (rich) = max (0., min (1., (.50-rich)/.25 ))

c----------------------------------------------------------------------

