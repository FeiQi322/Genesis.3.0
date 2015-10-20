#include <stdio.h>
#include <stdlib.h>
 
/*--------------------------------------------------------------------*/
 
/*  ctowsf_(), ctowsd_() and ctowsi_(): routines to convert Cray Fortran
 *	single precision reals and integers (64 bits each) to IEEE format
 *	(used by most workstations);
 *	it is assumed that char == 8 bits on all machines.
 *  Written April 1987 by Chuck D'Ambra.
 *  Last modified 31 March 1994 by CD.
 *  Next line is for netsend:
# '#' @
 */
 
#define CBIAS		040000	/*  Cray single precision exponent bias  */
#define SBIAS		126	/*  IEEE single precision exponent bias  */
#define DBIAS		1022	/*  IEEE double precision exponent bias  */
#define CSIGNMASK	0200	/*  Mask to get 1st of 8 bits  */
#define LCDIF		sizeof(long) - sizeof(char)
#define MBIAS_A		1022	/*  IEEE f.p. exponent bias  */
#define MBIAS_B		126 	/*  IEEE f.p. exponent bias  */
 
void ctowsf_();
void ctowsd_();
void ctowsi_();
void checkrng();
 
/*  Next 3 "dummy" routines handle difference between f77 and CFT  */
int CTOWSF(in, out, n)
unsigned char in[], out[];
int *n;
{
  ctowsf_(in, out, n);
  return(0);
}
 
 
int CTOWSD(in, out, n)
unsigned char in[], out[];
int *n;
{
  ctowsd_(in, out, n);
  return(0);
}
 
 
int CTOWSI(in, out, n)
unsigned char in[], out[];
int *n;
{
  ctowsi_(in, out, n);
  return(0);
}
 
/*--------------------------------------------------------------------*/
 
void ctowsf_(in, out, n)
unsigned char in[], out[];
int *n;
{	/*  Cray single precision floats to IEEE single precision floats  */
  register int i, j, joff;
  register unsigned char *p, sign, *expp;
  unsigned char *manp[3];
  long exp;
  static unsigned char maxman = 0xff, minman = 0;
 
/*  Cray 64 bit float representation:
 *	bits 0-47:	mantissa
 *	bits 48-62:	exponent
 *	bit 63:		sign
 */
 
  for (i = 0; i < *n; i++)
  {
    p = in+i*8;
    joff = i * 4;
 
    sign = *(p) & CSIGNMASK;	/*  sign bit  */
    for (j = 0; j < 3; j++)
      manp[j] = p+j+2;		/*  point at mantissa  */
 
    exp = (*(p+1)) | ((*p & 0x7f) << 8);
    if (exp != 0)
      exp = exp - CBIAS + SBIAS;
 
/*  If input is outside range representable on IEEE, set to
    closest representable number.  */
    if (exp > 254)		/*  Too large  */
    {
      exp = 254;
      for (j = 0; j < 3; j++)	/*  Set all bits in mantissa  */
        manp[j] = &maxman;
    }
    else if (exp < 0)		/*  Too small  */
    {
      exp = 1;
      for (j = 0; j < 3; j++)	/*  Turn off all bits in mantissa  */
        manp[j] = &minman;
    }
 
/*  Pack it into 32 bit IEEE float representation:
 *	bits 0-22:	mantissa
 *	bits 23-30:	exponent
 *	bit 31:		sign
 */
    expp = ((unsigned char *)(&exp)) + LCDIF;
    out[joff] = sign | (*expp >> 1);
    out[joff+1] = (*manp[0] & (~CSIGNMASK)) | ((*expp & 1) << 7);
    out[joff+2] = *manp[1];
    out[joff+3] = *manp[2];
  }
}
 
/*--------------------------------------------------------------------*/
 
void ctowsd_(in, out, n)
unsigned char in[], out[];
int *n;
{	/*  Cray single precision floats to IEEE double precision floats  */
  register int i, j, joff = 0;
  register unsigned char *p, sign;
  unsigned char *manp[6];
  register short exp;
  long dbias = DBIAS - CBIAS;
  static unsigned char maxman = 0xff, minman = 0;
 
/*  Cray 64 bit float representation:
 *	bits 0-47:	mantissa
 *	bits 48-62:	exponent
 *	bit 63:		sign
 */
 
  p = in;
  for (i = 0; i < *n; i++)
  {
    sign = *p & CSIGNMASK;	/*  sign bit  */
 
    exp = ((*p & 0x7f) << 8) | *(p+1);	/*  exponent  */
    if (exp != 0)
      exp = exp + dbias;	/*  convert to IEEE biased exponent  */
 
    p += 2;
    for (j = 0; j < 6; j++)
      manp[j] = p++;		/*  point at mantissa bytes  */
 
/*  If input is outside range representable in double precision IEEE,
    set to closest representable number.  */
    if (exp > 2048)		/*  Too large  */
    {
      exp = 2048;
      for (j = 0; j < 6; j++)	/*  Set all bits in mantissa  */
        manp[j] = &maxman;
    }
    else if (exp < 0)		/*  Too small  */
    {
      exp = 1;
      for (j = 0; j < 6; j++)	/*  Turn off all bits in mantissa  */
        manp[j] = &minman;
    }
 
/*  Pack it into 64 bit IEEE float representation:
 *	bits 0-51:	mantissa
 *	bits 52-62:	exponent
 *	bit 63:		sign
 */
 
/*  1st byte:  sign and 1st 7 bits of exponent */
    out[joff++] = sign | (exp >> 4);
/*  2nd byte:  last 4 bits of exponent and 1st 4 bits of mantissa  */
    out[joff++] = ((exp & 0xf) << 4) | ((*manp[0] & 0x7f) >> 3);
/*  remaining bytes contain mantissa  */
    out[joff++] = (*manp[0] << 5) | (*manp[1] >> 3);
    out[joff++] = (*manp[1] << 5) | (*manp[2] >> 3);
    out[joff++] = (*manp[2] << 5) | (*manp[3] >> 3);
    out[joff++] = (*manp[3] << 5) | (*manp[4] >> 3);
    out[joff++] = (*manp[4] << 5) | (*manp[5] >> 3);
    out[joff++] = *manp[5] << 5;
  }
}
 
/*--------------------------------------------------------------------*/
 
void ctowsi_(in, out, n)
unsigned char in[], out[];
int *n;
{	/*  Cray integers to IEEE integers  */
  register int i, j, joff;
  register unsigned char *pin, sign;
  unsigned char *p[4];
 
  for (i = 0; i < *n; i++)
  {
    pin = in+i*8;
    joff = i * 4;
 
    sign = *(pin) & CSIGNMASK;	/*  sign bit  */
    for (j = 0; j < 4; j++)
      p[j] = pin+j+4;		/*  point to last 4 bytes  */
 
    checkrng(sign,pin,p);	/*  check if int w/i IEEE limits  */
 
/*  Pack it into 4 byte IEEE integer representation  */
 
    out[joff] = sign | (*p[0] & (~CSIGNMASK));
    out[joff+1] = *p[1];
    out[joff+2] = *p[2];
    out[joff+3] = *p[3];
  }
}
 
/*--------------------------------------------------------------------*/
 
void checkrng(sign,cp,p)
unsigned char sign, *cp, *p[4];
{
  static unsigned char maxint = 0xff, minint = 0;
 
  if (sign == 0)	/*  int >= 0; 1st 33 bits must be 0.  */
  {
    if ((*(cp+4) & CSIGNMASK) != 0  ||  *cp != 0  ||  *(cp+1) != 0
        || *(cp+2) != 0  ||  *(cp+3) != 0)
      p[0] = p[1] = p[2] = p[3] = &maxint;
  }
  else			/*  int < 0; 1st 33 bits must be 1.  */
  {
    if ((*(cp+4) & CSIGNMASK) < 0200 || *cp != 0377  ||  *(cp+1) != 0377
        || *(cp+2) != 0377  ||  *(cp+3) != 0377)
      p[0] = p[1] = p[2] = p[3] = &minint;
  }
}
 
/*--------------------------------------------------------------------*/
 
/* Copyright 1994 University Corporation for Atmospheric Research (UCAR).
**	All rights reserved.
**
** Permission to use, copy, and modify this software and its documentation
** for any non-commercial purpose is hereby granted without fee, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation. UCAR makes no representations about the suitability of
** this software for any purpose.   It is provided "as is" without express
** or implied warranty.
*/
 
/*
** dctocf():  routine to convert  double IEEE to CRAY single precision reals.
**
**	it is assumed that char == 8 bits on all machines.
**
**  Created from a hack of D'Ambra's  C->IEEE code
**  by Dan Anderson for rs6k 14 March 1992
** Code cleanup Peter Morreale January 1994.
*/
	
void dptocf_();
 
/*-------------------------------------------------------------------------*/
 
/* Dummy routine to handle difference between f77 and CFT */
 
int DPTOCF (in, out, n)
unsigned char in[], out[];
int *n;
{
  dptocf_ (in, out, n);
  return(0);
}
 
/*-------------------------------------------------------------------------*/
 
void dptocf_(in, out, n)
unsigned char in[], out[];
int *n;
{
	register int 		i, j, joff;
	unsigned char 		*p, sign, *manp[7];
	long 			exp;
	static unsigned char 	minman = 0, and;
	
	/*
	** Cray 64 bit float representation:
	**	bits 0-47:	mantissa
	**	bits 48-62:	exponent
	**	bit 63:		sign
	*/
	
	for (i = 0; i < *n; i++) {
	    p = in+i*8;
	    joff = i * 8;
	
	    sign = *(p) & CSIGNMASK;	/*  sign bit  */
	
	    manp[0] = p+1;	
	    manp[1] = p+2;	
	    manp[2] = p+3;	
	    manp[3] = p+4;	
	    manp[4] = p+5;	
	    manp[5] = p+6;	
	    manp[6] = p+7;	
	
 
	    /*
	    ** Get the exponet
	    */
	    exp = ((*(p+1))>>4) | ((*p & 0x7f) << 4);
	    if (exp != 0) exp = exp - MBIAS_A + CBIAS;
	
	    /*
	    ** If input is outside range representable on IEEE, set to
	    ** closest representable number.
	    */
	    and=0xff;
	    if (exp <= 0) {  /* too small */
	        exp = 0;
	        /*
	        ** Turn off all bits in mantissa
	        */
	        for (j = 0; j < 7; j++) manp[j] = &minman;
	        and=0x00;
	    }
	
	    /*
	    ** Pack it into 64 bit CRAY float representation:
	    **	bits 0-47:	mantissa
	    **	bits 48-62:	exponent
	    **	bit 63:		sign
	    */
	
	    out[joff] = sign | (exp >> 8);  /* sign+upper 7 bits of exponent */
 
	    out[joff+1] = ((exp & 0xff));
	    out[joff+2] = (0x80 & and)|((*manp[0] & 0xf)<<3)|((*manp[1])>>5);
	    out[joff+3] = ((*manp[1] & 0x1f) <<3) | ((*manp[2] ) >>5);
	    out[joff+4] = ((*manp[2] & 0x1f) <<3) | ((*manp[3] ) >>5);
	    out[joff+5] = ((*manp[3] & 0x1f) <<3) | ((*manp[4] ) >>5);
	    out[joff+6] = ((*manp[4] & 0x1f) <<3) | ((*manp[5] ) >>5);
	    out[joff+7] = ((*manp[5] & 0x1f) <<3) | ((*manp[6] ) >>5);
	}
}
	
/*-------------------------------------------------------------------------*/
 
/* Copyright 1994 University Corporation for Atmospheric Research (UCAR).
**	All rights reserved
**
** Permission to use, copy, and modify this software and its documentation
** for any non-commercial purpose is hereby granted without fee, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation. UCAR makes no representations about the suitability of
** this software for any purpose.   It is provided "as is" without express
** or implied warranty.
*/
 
/*  sptoci():  routines to convert IEEE  single
**  precision integers (32 bits each) to CRAY (64 bit) format;
**	 it is assumed that char == 8 bits on all machines.
**  Hack of D'Ambras ctospf code
**  Created by Dan Anderson in March of 92
*/
	
void sptoci_();
 
/*-------------------------------------------------------------------------*/
 
/* Dummy routine to handle difference between f77 and CFT */
 
int SPTOCI (in, out, n)
unsigned char in[], out[];
int *n;
{
  sptoci_ (in, out, n);
  return(0);
}
 
/*-------------------------------------------------------------------------*/
 
void sptoci_(in, out, n)
unsigned char in[], out[];
int *n;
{
	register int 		i, j;
	register unsigned char 	*pin;
	
	for (i = 0; i < *n; i++) {
	
	    pin = in + i*4;
	    j = i * 8;
	
	    out[j] = *(pin) & CSIGNMASK;	/*  sign bit  */
	
	    if( out[j] != 0x00) {
		out[j] = 0xff;
		out[j+1]  = 0xff;
		out[j+2]  = 0xff;
		out[j+3]  = 0xff;
	    }
	    else {
		out[j+1]  = 0x00;
		out[j+2]  = 0x00;
		out[j+3]  = 0x00;
	    }
	
	    out[j+4] = *pin;
	    out[j+5] = *(pin+1);
	    out[j+6] = *(pin+2);
	    out[j+7] = *(pin+3);
	}
}
 
/*-------------------------------------------------------------------------*/
