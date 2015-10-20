
#ifdef isotrac

c-----------
c   Equilibrium fractionation factor between liquid and vapour and
c   solid and vapour for 18O and 2H (T in Kelvin).
c-----
      avl_h(T) = exp(24844./T/T-76.248/T+0.052612)
      avl_o(T) = exp(1137./T/T-0.4156/T-0.0020667)
      avi_h(T) = exp(16288/T/T-0.0934)
      avi_o(T) = exp(11.839/T-0.028224)
c-----------
c   dv_i: coefficient of diffusion of water vapour in air
c         (T in K, p in Pa).
c-----
      dv_i(T,P)= 2.17e-5*(101325./P)*(T/273.15)**1.88

#endif

