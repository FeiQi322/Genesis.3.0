c -- $Id: COMLUN.h,v 1.2 2001/07/03 21:29:08 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMLUN.h,v $
c -- $Name:  $

C
c        Logical unit numbers, path names, file names, rec/file sizes,
c        write password, retention period.
c        This is the only common shared with lsx code, so that
c        writh and writhlsx can access both agcm and lsx filenames,
c        and openf can access file and record sizes set in agcm and lsx.
c
c        "MSS" is the Mass Store pathname if mss is used, or the local
c        pathname if not (see acquiref and disposef in filio.f). 
c
c        plun   = # of local disk-file luns, excluding nin,nout,nuser*
c        pfile  = maximum length of user-specified filename lists
c        nsre...= disk-file logical unit numbers
c        nin    = standard input  logical unit number
c        nout   = standard output logical unit number
c        nuser* = logical unit numbers for user-codeable output(userdef)
c        lociu  = identical to luns nsre..., used for lookup into
c                 locname(plun) (see filio.f)
c
c        locname  = local file names
c        mspath   = MSS pathname for hist/rest files written by this run
c        mspathin = MSS pathname for init/rest file to start this run
c        msnamin  = MSS filename of  init/rest file to start this run
c        
c        nfile    = length of restart/history filename lists
c        ifile    = current position in restart/history filename lists
c        ifilepre = ifile unincremented by next-file (for restart file)
c
c        msnamres = list of restart filenames written to MSS
c        msnamre2 = a single restart filename written to MSS at same
c                   time as msnamres entries
c
c        msnamhis = list of agcm history filenames written to MSS
c        msnamhisx= current agcm history filename (for restart file)
c        nuoff    = # of records on current history file
c        nwrit    = # of history writes on current history file
c        ndens    = packing density for agcm history file data recs
c        nurecl   = agcm history file record length in words
c        nusiz    = size of one agcm history file in words
c
c        msnam_a, etc = corresp variables for lsx history file
c
c        passrd   = read  password   for all MSS files accessed this run
c        passwd   = write password   for all MSS files written this run
c        retpd    = retention period for all MSS files written this run
c
c     psoilun is for the soil texture files (one per each soil layer),
c     so must be >= lsx's psoilay. But not so big that it collides with
c     lsx's iutab* (50,...). Since lociu(i) is set to 10+i in init.f,
c     this means plun must be <= 39.

      integer presoil, psoilun, postsoil, plun, pfile 
      parameter (presoil=13, psoilun=12, postsoil=7,
     *           plun=presoil+psoilun+postsoil, 
     *           pfile=2000)
c
      common /comlun/

c        data input files:
     *   ntopog,     nozone,     ngwdrag,    ntopog2,   
     *   iuvegt,     iuwisc,     iusurt,     iuwater,
     *   nocean,     nsstice, 
     *   iuocean_i,  iusourc_i,  iunamel,  

c        soil texture data input files: 
     *   iusoit(psoilun),

c        restart,history files:
     *   nsre,      ndata,     ndata_a,
     *   iureg,     iureg_a,   iustorm,
     *   nsreibis,

c        misc files (not spanned by lociu):
     *   nin,       nout,      nuser0,    nuser1,    nuser2

      dimension lociu(plun)
      equivalence (lociu,ntopog)
c
      common /comlunc/ 
     *   locname_nl(plun), locname(plun),
     *   mspath,           mspathin,       msnamin,
     *   nfile,            ifile,          ifilepre,
     *   msnamres(pfile),  msnamre2,
     *   msnamhis(pfile),  msnam_a(pfile),
     *   msnamhisx,        msnam_ax,
     *   nuoff,            nuoff_a,
     *   nwrit,            nwrit_a,
     *   ndens,            ndens_a,        ndens_storm,
     *   nurecl,           nurecl_a,
     *   nusiz,            nusiz_a,
     *   regpath,          reghis,         reghis_a,       regfreq,
     *   stormpath,        stormhis,                       stormfreq,
     *   ibispathin,       ibisresin,      ibisres,
     *   passrd,           passwd,         retpd, 
     *   fatalwms
#ifdef isotrac
c    (i,j) point output and frequency (timesteps); optional reset of
c    surface reservoirs to hrst_i.
     *  ,pii_i,     pjj_i,     pdn_i,     hrst_i
#endif
c
      character
     *   locname_nl*96,    locname*96,
     *   mspath*80,        mspathin*80,    msnamin*16,
     *   msnamres*16,      msnamre2*16,
     *   msnamhis*16,      msnam_a*16,
     *   msnamhisx*16,     msnam_ax*16,
     *   regpath*80,       reghis*16,      reghis_a*16,
     *   stormpath*80,     stormhis*16,
     *   ibispathin*80,    ibisresin*16,   ibisres*16,
     *   passrd*8,         passwd*8
c
      logical fatalwms

c     for netcdf agcm and lsx (_a)history files:
      common /chisnetc/
     *   hisnet,       ncid,         ncid_a,
     *   cusfielnet,   numfielnet,
     *   cusfielnet_a, numfielnet_a
      character*8 cusfielnet(2000), cusfielnet_a(2000)
      logical hisnet
#ifdef cray
      integer
#else
      integer*4
#endif
     *  ncid, ncid_a

#ifdef isotrac
      integer pii_i,       pjj_i,       pdn_i
#endif


