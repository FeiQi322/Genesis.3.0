c -- $Id: PARAMS.h,v 1.2 2001/07/03 21:29:15 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/PARAMS.h,v $
c -- $Name:  $

C
CL             AGCM PARAMETER DEFINITIONS
C
c    General save statement for all subrs, to assign all local
c    variables to static storage (not stack) on Cray
#ifdef saveall
      SAVE                                           
#endif

c     On sgi, declare nint and int external, and supply own versions
c     in utilsx, to avoid sgi bug with -O3 and -i8.
#ifdef sgi
      external nint, int
#endif

#if defined (NOFLUSH) || defined (OWNFLUSH)
      external flush
      integer flush
#endif

C     AGCM dimension sizes (grid, and spectral truncation):

c=======================
#include <AGCMSIZE.h>
c=======================

      INTEGER 
     * NOREC2, NLEVM, NLEVP, NLEVP2, NLEV2, NLEV2K,
     * MMAX, NMAX, NPMAX, MMAX2, NCRAY,
     * NUSHMX, NMARK,
     * NFFT1,  NLOG1, NLOGIC,
     * NTRACE

C       MACHINE WORD SIZE (OPTIMUM FFT BLOCK SIZE)
      PARAMETER (NCRAY=64)

C       SECONDARY SPECTRAL AND GRID DIMENSION VALUES
      PARAMETER (
     *  MMAX=NTRM+1, NMAX=NTRN+1, NPMAX=NTRK+1, MMAX2 = 2*MMAX,
     *  NLEVM=NLEV-1, NLEVP=NLEV+1, NLEVP2=NLEV+2, NLEV2=NLEV*NLEV,
     *  NLEV2K=NLEV2*NPMAX, NOREC2=NOREC/2 )

C       MAX # OF HISTORY FIELDS (PROGRAM LIST),
C       MAX # OF HISTORY PERIODS/SAVES/NEXT-FILES PER YEAR
      PARAMETER (NUSHMX=200, NMARK=20000)

C       FFT PARAMETERS
      PARAMETER (
     *  NFFT1=3*NLON/2+1,
     *  NLOG1=(100+NLEV)/NLON, NLOGIC=NLON*(NLOG1+1) )

C       NTRACE = TOTAL NUMBER OF TRACER FIELDS (IN COML50).
C       FIELDS 1 TO NTRACA ARE EQUIVALENCED TO LIQUID WATER
C       CONTENTS AND FLUX RESERVOIRS; FIELDS NTRACA+1
C       TO NTRACA+NTRACB (NTRACE) ARE AVAILABLE FOR TRUE TRACERS.
      PARAMETER (NTRACA = 7)
#ifdef cosmog
      PARAMETER (NTRACB = 12)
#else
      PARAMETER (NTRACB = 1)
#endif
      PARAMETER (NTRACE = NTRACA + NTRACB)

C       Length of horizontal strips, and # of strips, for radiation code
      PARAMETER (NRAD    = NLON,
     *           NSTRIPR = (NLON*NOREC-1)/NRAD + 1 )
 
c#ifdef isotrac
c   number of regional tracers
      parameter (nrtr =  1)
c#endif

