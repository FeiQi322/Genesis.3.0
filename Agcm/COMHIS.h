c -- $Id: COMHIS.h,v 1.2 2001/07/03 21:29:07 bill Exp $
c -- $Source: /usr/local/portgcm/src/cvsroot/genesis.2.2/Agcm/COMHIS.h,v $
c -- $Name:  $

C
CL            INDICES FOR HISTORY FIELDS. CUSHIS,MUSHIS ARE FOR USER-
CL            REQUESTED LIST. CINDEX AND MINDEX ARE FOR HARD-CODED
CL            PROGRAM LIST.
C
C                CUSHNAM(*)  = namelist input for CUSHIS amd MUSHFL
C                CUSHIS(*)   = list of user-requested labels (1-8 chars)
C                CINDEX(1,*) = program field labels (1-8 chars)
C                CINDEX(2,*) = units for each program field (1-8 chars)
C                MUSHIS(*)   = position in CINDEX for each user label
C                MUSHFL(*)   = 0 for instantaneous, 1 for accumulated 
C                NUSHIS      = number of user-requested fields
C                MINDEX(1,*) = spare
C                MINDEX(2,*) = number of levels (1 or 12)
C                MINDEX(3,*) = 1st  rec # rel to last header record
C                MINDEX(4,*) = last rec # rel to last header record
C                MINDFL(*)   = MUSHFL for this program field
C                MINDCO(*)   = counter for accumulated fields
C                NINDEX      = number of hard-coded program fields
C                IHTRACE     = position in CINDEX of first true tracer
C                              field name
C
C                NULHDR = number of header records for each hist write
C                NULDAT = number of data records for each hist write
C
C                NHISI  = namelist input for NHIST
C                NDHIS  = time interval for each accumulation period
C                HISTLAST = .T. to force history write at last timestep
C                NHIST  = start,end times for each accum period (per yr)
C                NHISTN = number of acumulation intervals (per year)
C                FWHIST = .TRUE. if history write this timestep
C                FAHIST = .TRUE. if within an accum period this timestep
C
C                NSAVE  = times for history file saves (per year)
C                NSAVEN = number of save times (per year)
C                FSAVE  = .TRUE. if save this timestep
C
C                NNEXF  = times for history next-files (per year)
C                NNEXFN = number of next-file times (per year)
C                FNEXF  = .TRUE. if next-file this timestep
C
      COMMON /COMHIS/
     +   CUSHNAM(NUSHMX),   CUSHIS(NUSHMX),   CINDEX(2,NUSHMX),
     +   MUSHIS(NUSHMX),    MUSHFL(NUSHMX),   NUSHIS,
     +   MINDEX(4,NUSHMX),  MINDFL(NUSHMX),   MINDCO(NUSHMX),  NINDEX,
     +   IHTRACE,
C
     +   NULHDR,   NULDAT,
C
     +   NHISI(NMARK),   NDHIS,   HISTLAST,
     +   NHIST(NMARK,2), NHISTN,  FWHIST,  FAHIST,
     +   NSAVE(NMARK),   NSAVEN,  FSAVE,
     +   NNEXF(NMARK),   NNEXFN,  FNEXF
C
      CHARACTER*8 CUSHNAM, CUSHIS, CINDEX
      LOGICAL HISTLAST, FWHIST, FAHIST, FSAVE, FNEXF

