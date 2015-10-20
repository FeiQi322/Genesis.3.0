/* scf:  better fscanf and other library functions */



# include	<ctype.h>
# include	<math.h>
# include	<stdio.h>
# include	<stdlib.h>
# include	<string.h>



double	string_to_double ();


void		reverse  (),
			reverses ();



void	fscanfi (sfp, sfw, svp)

FILE	*sfp;		/* scan file pointer     */

int	sfw,		/* scan field width      */
		*svp;		/* scan variable pointer */

{

	char	*ncs,
			ncs_0 [30];

	int	i;

	ncs = ncs_0;

	for (i = 0; i < sfw; i++)
		*ncs++ = fgetc(sfp);

	*ncs = '\0';

	*svp = atoi(ncs_0);

}



void	fscanfd (sfp, sfw, svp)

FILE		*sfp;		/* scan file pointer     */

int		sfw;		/* scan field width      */

double	*svp;		/* scan variable pointer */

{

	char	*ncs,
			ncs_0 [30];

	int	i;

	ncs = ncs_0;

	for (i = 0; i < sfw; i++)
		*ncs++ = fgetc(sfp);

	*ncs = '\0';

	*svp = string_to_double(ncs_0);

}



void	fscanff (sfp, sfw, svp)

FILE		*sfp;		/* scan file pointer     */

int		sfw;		/* scan field width      */

float		*svp;		/* scan variable pointer */

{

	char	*ncs,
			ncs_0 [30];

	int	i;

	ncs = ncs_0;

	for (i = 0; i < sfw; i++)
		*ncs++ = fgetc(sfp);

	*ncs = '\0';

	*svp = (float) string_to_double(ncs_0);

}



void	fskipc (scfp, scn)

FILE		*scfp;	/* skip character file pointer */

int		scn;		/* skip character number       */

{

	char	ch0;

	int	i;

	for (i = 0; i < scn; i++)
		ch0 = fgetc(scfp);

}



void	fskips (ssfp, ssn)

FILE		*ssfp;	/* skip string file pointer */

int		ssn;		/* skip string number       */

{

	int	i;

	for (i = 0; i < ssn; i++)
		fscanf(ssfp, "%*s");

}



void	fskipl (slfp, sln)

FILE		*slfp;	/* skip line file pointer */

int		sln;		/* skip line number       */

{

	char	ch0;

	while (sln)
	{
		while ((ch0 = fgetc(slfp)) != '\n');

		--sln;
	}

}



double	string_to_double (str)

char	str[];

{

	double	val, power;
	int		i, sign;

	for (i = 0; isspace(str[i]); i++)
		;

	sign = (str[i] == '-') ? -1 : 1;

	if (str[i] == '+' || str[i] == '-')
		i++;

	for (val = 0.0; isdigit(str[i]); i++)
		val = 10.0 * val + (str[i] - '0');

	if (str[i] == '.')
		i++;

	for (power = 1.0; isdigit(str[i]); i++)
	{
		val = 10.0 * val + (str[i] - '0');
		power *= 10.0;
	}

	return (sign * val / power);

}



int	str2int (str, fw)

char	str[];

int	fw;		/* field width */

{

	int	i, j,
			sign,
			val;

	for (i = 0; (i < (fw-1)) && (isspace(str[i]) || (str[i] == '0')); i++);

	sign = (str[i] == '-') ? -1 : 1;

	if (str[i] == '+' || str[i] == '-')
		i++;

	for ( ; (i < (fw-1)) && (str[i] == '0'); i++);

	if ((i == fw) || !isdigit(str[i]))
	{
		printf("\n\n  No decimal-digit character found in str2int (scf.c) \n\n");

		printf("    FW     = %d \n", fw);

		for (j = 0; j < fw; j++)
			printf("    STR[%d] = %c \n", j+1, str[j]);

		printf("\n\n  Exiting ..... \n\n\n");

		exit(0);
	}

	for (val = 0; (i < fw) && isdigit(str[i]); i++)
		val = (10 * val) + (str[i] - '0');

	return (sign * val);

}



void	int2str (cs, n, fw, sf)

char	cs[];		/* char string */

int   n,			/* integer     */
		fw,		/* field width */
		sf;		/* sign flag   */

{

	char	s [70];

	int	i,
			nfw,		/* numeric field width */
			sign;

	if ((sign = n) < 0)		/* record sign */
	{
		if (sf)
		{
			n = -n;				/* make n positive */
		}
		else
		{
			printf("\n\n  Inconsistent args in int2str (scf.c): n<0, sf=0 \n\n");
			printf("  Exiting ..... \n\n\n");

			exit(0);
		}
	}

	i = 0;

	do
	{											/* generate digits in reverse order */
		s[i++] = (n % 10) + '0';		/* get next digit */
	}
	while ((n /= 10) > 0);				/* delete it */

	if (sf)
		nfw = fw - 1;
	else
		nfw = fw;

	while (i < nfw)
		s[i++] = '0';

	if (sf)
	{
		if (sign >= 0)
			s[i++] = '+';
		else
			s[i++] = '-';
	}

	s[i] = '\0';

	if (i <= fw)
	{
		reverses(s, i);
	}
	else
	{
		printf("\n\n  Numerical overflow in int2str (scf.c): strlen(s)>fw \n\n");

		printf("    S    = %s \n", s);
		printf("    N    = %d \n", n);
		printf("    SF   = %d \n", sf);
		printf("    SIGN = %d \n", sign);

		printf("    I    = %d \n", i);
		printf("    FW   = %d \n", fw);
		printf("    NFW  = %d \n", nfw);

		printf("\n\n  Exiting ..... \n\n\n");

		exit(0);
	}

	for (i = 0; i < fw; i++)
		cs[i] = s[i];

}



char  *int_to_ascii (n, s)

int   n;

char s[];

{

	int   i,
			sign;

	if ((sign = n) < 0)  /* record sign */
		n = -n;           /* make n positive */

	i = 0;

	do
	{                             /* generate digits in reverse order */
		s[i++] = n % 10 + '0';     /* get next digit */
	}
	while ((n /= 10) > 0);     /* delete it */

	if (sign < 0)
		s[i++] = '-';

	s[i] = '\0';

	reverses(s, i);

	return(s);

}



void  reverse (s)

char  s[];

{

	int   c, i, j;

	for (i = 0, j = strlen(s)-1; i < j; i++, j--)
	{
		c    = s[i];
		s[i] = s[j];
		s[j] = c;
	}

}



void  reverses (s, sl)

char  s[];

int	sl;

{

	int   c, i, j;

	for (i = 0, j = sl-1; i < j; i++, j--)
	{
		c    = s[i];
		s[i] = s[j];
		s[j] = c;
	}

}



void	report (fn)

char	*fn;

{

	printf("\n\n  Can't open file:  %s \n\n\n", fn);
	exit(0);

}



FILE	*open_file (fn)

char	*fn;

{

	FILE	*fp;

	if ((fp = fopen(fn, "r")) == (FILE *) NULL)
		report(fn);

	return(fp);

}



FILE	*open_write (fn)

char	*fn;

{

	FILE	*fp;

	if ((fp = fopen(fn, "w")) == (FILE *) NULL)
		report(fn);

	return(fp);

}



char	*lc_to_uc (cs0)

char	*cs0;

{

	char	*cs1;

	int	i;

	cs1 = (char *) calloc(strlen(cs0), 1);

	for (i = 0; i < strlen(cs0); i++)
		cs1[i] = toupper(cs0[i]);

	return(cs1);

}



char	*uc_to_lc (cs0)

char	*cs0;

{

	char	*cs1;

	int	i;

	cs1 = (char *) calloc(strlen(cs0), 1);

	for (i = 0; i < strlen(cs0); i++)
		cs1[i] = tolower(cs0[i]);

	return(cs1);

}



double   anint(x)

double   x;

{

	double   y;

	if (x >= 0.0)
		y = floor(x+0.5);
	else
		y = ceil(x-0.5);

	return(y);

}



int   nint(x)

double   x;

{

	int   y;

	if (x >= 0.0)
		y = (int) floor(x+0.5);
	else
		y = (int) ceil(x-0.5);

	return(y);

}

