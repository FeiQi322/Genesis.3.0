C     AGCM dimension sizes (grid, and spectral truncation):

c========================
#if defined  (STANDALONE)
c========================

c     for standalone driverad (radiation test)
      integer nlon, norec, nlev
      parameter (nlon=1, norec=1, nlev=18)
c     parameter (nlon=1, norec=48, nlev=18)

      integer ntrm, ntrn, ntrk
      parameter (ntrm=1, ntrn=1, ntrk=1)

c====
#else
c====

C
CL             AGCM DIMENSIONS
C
      INTEGER 
     * NLON, NOREC, NLEV,
     * NTRM, NTRN, NTRK

      PARAMETER
     * (
#ifdef A_T10   
     *  NTRM=10, NTRN=10, NTRK=10, NLON= 32, NOREC=16, NLEV=6    !T10   
#endif
#ifdef A_R15_12
     *  NTRM=15, NTRN=15, NTRK=30, NLON= 48, NOREC=40, NLEV=12   !R15_12
#endif
#ifdef A_R15
     *  NTRM=15, NTRN=15, NTRK=30, NLON= 48, NOREC=40, NLEV=18   !R15
#endif
#ifdef A_R30
     *  NTRM=30, NTRN=30, NTRK=60, NLON= 96, NOREC=80, NLEV=18   !R30
#endif
#ifdef A_T31_12
     *  NTRM=31, NTRN=31, NTRK=31, NLON= 96, NOREC=48, NLEV=12   !T31_12
#endif
#ifdef A_T31
     *  NTRM=31, NTRN=31, NTRK=31, NLON= 96, NOREC=48, NLEV=18   !T31
#endif
#ifdef A_T31_36
     *  NTRM=31, NTRN=31, NTRK=31, NLON= 96, NOREC=48, NLEV=36   !T31_36
#endif
#ifdef A_T42_12
     *  NTRM=42, NTRN=42, NTRK=42, NLON=128, NOREC=64, NLEV=12   !T42_12
#endif
#ifdef A_T42
     *  NTRM=42, NTRN=42, NTRK=42, NLON=128, NOREC=64, NLEV=18   !T42
#endif
#ifdef A_T63_12
     *  NTRM=63, NTRN=63, NTRK=63, NLON=192, NOREC=96, NLEV=12   !T63_12
#endif
#ifdef A_T63
     *  NTRM=63, NTRN=63, NTRK=63, NLON=192, NOREC=96, NLEV=18   !T63
#endif
#ifdef A_T106_12
     *  NTRM=106,NTRN=106,NTRK=106,NLON=320, NOREC=160,NLEV=12  !T106_12
#endif
#ifdef A_T106
     *  NTRM=106,NTRN=106,NTRK=106,NLON=320, NOREC=160,NLEV=18   !T106
#endif
     * )

c=====
#endif
c=====
