/* LEAFF:  the Leaf Area Index (LAI) and canopy phenology model:
           the main function */


# include	"canopy_fps.h"		/* Canopy file pointers             */
# include	"ddf.h"				/* Deciduous-day functions          */
# include	"gv.h"				/* Global variables                 */
# include	"lf.h"				/* Life forms                       */
# include   "lfo.h"				/* Life-form overlaps               */
# include	"lfpcl.h"			/* Life-form plant community levels */
# include   "lfsi.h"				/* Life-form shade intolerance      */
# include	<math.h>
# include	<stdio.h>
# include	<string.h>


# define		daily       1
# define		image       0
# define		leafy       0
# define		leafy2      0
# define		trace       0



enum seasonal_cycle { winter = 0, spring, summer, autumn, winter0 };


extern
enum seasonal_cycle	lf_season [num_lf];


extern
char		earth [n_lats] [n_lons];


	/* Daily LAIs:      1 -> green leaves                                     */
	/*                  0 -> dead leaves                                      */

	/* Daily LAIs:          0 -> trees with (partial) vines and epiphytes     */
	/*                      1 -> trees                                        */
	/*                      2 -> vines                                        */
	/*                      3 -> epiphytes                                    */
	/*                      4 -> shrubs                                       */
	/*                      5 -> grasses                                      */
	/*                      6 -> total plant community within the VAR         */
	/*                      7 -> total plant community at the landscape level */

extern
double	canopy_lai,
			canopy_laid,
			seasonal_lai [2] [num_lf];


extern
int		leaf_fall_day0 [4],
			leaf_fall_day1 [4],

			winter_day;


/* EVE Plant Community Upper and Lower Canopy Parameters (version 1) */


typedef
double	dcc;


extern
struct constant_canopy
{
	dcc   fc,	 /* fractional cover, canopy, landscape level */
			fcc,	 /* fractional cover, canopy                  */

			z1,	 /* height of canopy, top    (m) */
			z0,	 /* height of canopy, bottom (m) */

			sai,	 /* stem (surface) area index, within canopy (m^2 / m^2) */
			sd,	 /* typical stem dimension                   (m)         */

			rd;	 /* root depth, conical form (m) */
};


extern
struct seasonal_canopy
{
	dcc	lai,	 /* green leaf (surf) area index, within canopy (m^2 / m^2) */
			laid,	 /*  dead leaf (surf) area index, within canopy (m^2 / m^2) */

			flb,	 /* fraction of leaves that are broadleaf */

			llb,	 /* leaf length, broadleaf  (m) */
			lln,	 /* leaf length, needleleaf (m) */
			lwb,	 /* leaf width, broadleaf   (m) */
			lwn,	 /* leaf width, needleleaf  (m) */

			lo,	 /* leaf orientation (-1:vertical, 0:random, 1:horizontal) */

			lrv,	 /* green leaf reflectivity, visible   */
			lri,	 /* green leaf reflectivity, near IR   */
			ltv,	 /* green leaf transmissivity, visible */
			lti,	 /* green leaf transmissivity, near IR */

			scmax, /* green leaf stomatal conductance, maximum     (m / s)   */
			scmin, /* green leaf stomatal conductance, minimum     (m / s)   */
			scpar, /* green leaf stomatal conductance PAR constant (W / m^2) */
			scvpd; /* green leaf stomatal conductance VPD constant (N / m^2) */
};


extern
struct canopies		/* level 1: upper canopy,  level 0: lower canopy */
{
	struct constant_canopy	constant [2];
	struct seasonal_canopy	seasonal [2];

}	canopy;



FILE		*cfp1,
			*cfp2,
			*openff (),
			*open_file ();


char		*filename (),
			*int_to_ascii ();


int		ldm[12] = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
			mmd[12] = {16, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350},

			month = 1,

			output_day;


void		animate_leaves (),
			animate_sc (),
			assign_ddf (),
			close_leaf_files (),
#ifdef cray
			EVE_LEAF (),
#else
			eve_leaf_ (),
#endif
			initialization (),
			monitor_leaf (),
			open_leaf_files (),
			record_cc (),
			record_pc2 (),
			set_pc_lfs (),
			spin_up (),
			zero_vectors ();


typedef
int		ged [n_lats] [n_lons];


typedef
double	ccd      [2][n_lats] [n_lons],
			ccd2 [2] [2][n_lats] [n_lons],
			lcd         [n_lats] [n_lons];



#ifdef cray
void	EVE_LEAF
#else
void	eve_leaf_
#endif

		(

int	*iflag,				/* interactive flag */
int	*fiflag,				/* filter flag      */

int	*first_day,			/* 0 = no action      */
								/* 1 = read LEAF file */

int	*date1,				/* Day number  */
int	*year,				/* Year number */

		/* Daily climate */

lcd	global_prec,
lcd	global_temp,

		/* Earth map */

ged	global_earth,

		/* Constant Canopy parameters */

ccd   cc_fc,		/* fractional cover, canopy, landscape level */

ccd	cc_z1,		/* height of canopy, top    (m) */
ccd	cc_z0,		/* height of canopy, bottom (m) */

ccd	cc_sai,		/* stem (surface) area index, within canopy (m^2 / m^2) */
ccd	cc_sd,		/* typical stem dimension                   (m)         */

ccd	cc_rd,		/* root depth, conical form (m) */

		/* Seasonal Canopy parameters */

ccd	sc_lai,		/* green leaf (surf) area index, within canopy (m^2 / m^2) */
ccd	sc_laid,		/*  dead leaf (surf) area index, within canopy (m^2 / m^2) */

ccd	sc_flb,		/* fraction of leaves that are broadleaf */

ccd	sc_llb,		/* leaf length, broadleaf  (m) */
ccd	sc_lln,		/* leaf length, needleleaf (m) */
ccd	sc_lwb,		/* leaf width, broadleaf   (m) */
ccd	sc_lwn,		/* leaf width, needleleaf  (m) */

ccd	sc_lo,		/* leaf orientation (-1:vertical, 0:random, 1:horizontal) */

ccd2	sc_lrvi,		/* green leaf reflectivity, visible and near IR   */
ccd2	sc_ltvi,		/* green leaf transmissivity, visible and near IR */

ccd	sc_scmax,	/* green leaf stomatal conductance, maximum     (m / s)   */
ccd	sc_scmin,	/* green leaf stomatal conductance, minimum     (m / s)   */
ccd	sc_scpar,	/* green leaf stomatal conductance PAR constant (W / m^2) */
ccd	sc_scvpd		/* green leaf stomatal conductance VPD constant (N / m^2) */

		)

{

	double	vfc;

	int		i, j,

				frequency,

				jj = 0;

	if (trace)
		global_monitor = 1;
	else
		global_monitor = 0;

	if (global_monitor)
		monitor_leaf(1, "eve_leaf");

	/* printf("\n  The Leaf model. \n\n");                     */
	/* printf("    first_day = %5d \n", *first_day);           */
	/* printf("    date1     = %5d \n", *date1);               */
	/* printf("    year      = %5d \n\n", *year);              */
	/* printf("    cc_fc     = %11.5lf \n\n", cc_fc[0][0][0]); */

	interactive = *iflag;
	filtered    = *fiflag;

	site_dimensions = 2;		/* global mode */

	deg_per_lat = 180.0 / ((double) n_lats);
	deg_per_lon = 360.0 / ((double) n_lons);

	if (*first_day == 1)
	{
	    printf("\n  The Leaf model. \n\n");
	    printf("    date1     = %5d \n", *date1);
	    printf("    year      = %5d \n\n", *year); 
		printf("    Calculating DDFs and reading LEAF file \n\n");

		assign_ddf();
		read_leaf_file();

		if (interactive && filtered)
		{
			printf("    Reading climate files and calculating deltas \n\n");

			read_deltas();
		}

		if (*date1 > 1)
		{
			printf("    Spinning up LEAF to day %d \n\n", *date1);

			spin_up(*date1);
		}

		for (i = 0; i < n_lats; i++)
		{
			latgc = (n_lats - 1) - i;

			for (j = 0; j < n_lons; j++)
			{
				longc = (j + (n_lons / 2)) % n_lons;

				global_earth[latgc][longc] = (int) (earth[i][j] - '0');
			}
		}

	}

	date = *date1;

	if (*year == 2001)
		frequency = 1;
	else
		frequency = 5;

	if (leafy && ((daily && ((date % frequency) == 0)) ||
					  (date == mmd[month-1])))
		output_day = 1;
	else
		output_day = 0;

	open_leaf_files(*year);

	for (lat = 0; lat < n_lats; lat++)
	{
/*
		printf("               Lat = %3d \n", lat+1);
*/
		latitude =
		 (90.0 - (deg_per_lat / 2.0)) - (((double) lat) * deg_per_lat);

		/* set starting point for DDs:  one month after the winter solstice */

		if (latitude >= 0.0)
			winter_day =  21;
		else
			winter_day = 203;

		latgc = (n_lats - 1) - lat;

		for (lon = 0; lon < n_lons; lon++)
		{
			longc = (lon + (n_lons / 2)) % n_lons;

			if (earth[lat][lon] == '1')
			{
				read_data(jj);
/*
				if (!jj)
					printf("       LFC (%5d: %3d, %3d) = %2d \n",
					 jj+1, lat+1, lon+1, leaf_lfc);
*/
				if (leaf_lfc)
				{
					set_pc_lfs();

					if (interactive && filtered)
						set_deltas(jj);

					if (!interactive || filtered)
					{
						calc_daily_climate();
					}
					else
					{
						leaf_prec = 0.0;
						leaf_temp = 0.0;
					}

					if (interactive)
						add_gcm_climate(global_prec, global_temp);

					calc_inertial_climate();

					green_leaf(jj);		/* calculate leaf phenology             */

					calc_lais();			/* leaf phenology:  1x1 landscape level */
				}
				else
				{
					canopy_lai  = 0.0;
					canopy_laid = 0.0;

					for (i = 0; i < num_lf; i++)
					{
						seasonal_lai[1][i] = 0.0;
						seasonal_lai[0][i] = 0.0;
					}
				}

				write_data(jj);

				++jj;
			}
			else
			{
				leaf_var = 0.0;

				canopy_lai  = 0.0;
				canopy_laid = 0.0;

				for (i = 0; i < num_lf; i++)
				{
					leaf_fc[i] = 0.0;

					seasonal_lai[1][i] = 0.0;
					seasonal_lai[0][i] = 0.0;
				}
			}

			if ((earth[lat][lon] == '1') && leaf_lfc)
			{
				create_2_canopies();		/* calculate 2-level canopy properties */

				for (i = 0; i < 2; i++)
				{
					/* Constant Canopy parameters */

					cc_fc [i][latgc][longc] = canopy.constant[i].fc;

					cc_z1 [i][latgc][longc] = canopy.constant[i].z1;
					cc_z0 [i][latgc][longc] = canopy.constant[i].z0;

					cc_sai[i][latgc][longc] = canopy.constant[i].sai;
					cc_sd [i][latgc][longc] = canopy.constant[i].sd;

					cc_rd [i][latgc][longc] = canopy.constant[i].rd;

					/* Seasonal Canopy parameters */

					sc_lai    [i][latgc][longc] = canopy.seasonal[i].lai;
					sc_laid   [i][latgc][longc] = canopy.seasonal[i].laid;

					sc_flb    [i][latgc][longc] = canopy.seasonal[i].flb;

					sc_llb    [i][latgc][longc] = canopy.seasonal[i].llb;
					sc_lln    [i][latgc][longc] = canopy.seasonal[i].lln;
					sc_lwb    [i][latgc][longc] = canopy.seasonal[i].lwb;
					sc_lwn    [i][latgc][longc] = canopy.seasonal[i].lwn;

					sc_lo     [i][latgc][longc] = canopy.seasonal[i].lo;

					sc_lrvi[i][0][latgc][longc] = canopy.seasonal[i].lrv;
					sc_lrvi[i][1][latgc][longc] = canopy.seasonal[i].lri;
					sc_ltvi[i][0][latgc][longc] = canopy.seasonal[i].ltv;
					sc_ltvi[i][1][latgc][longc] = canopy.seasonal[i].lti;

					sc_scmax  [i][latgc][longc] = canopy.seasonal[i].scmax;
					sc_scmin  [i][latgc][longc] = canopy.seasonal[i].scmin;
					sc_scpar  [i][latgc][longc] = canopy.seasonal[i].scpar;
					sc_scvpd  [i][latgc][longc] = canopy.seasonal[i].scvpd;
				}
			}
			else
			{
				zero_vectors();

				for (i = 0; i < 2; i++)
				{
					/* Constant Canopy parameters */

					cc_fc [i][latgc][longc] = 0.0;

					cc_z1 [i][latgc][longc] = 0.0;
					cc_z0 [i][latgc][longc] = 0.0;

					cc_sai[i][latgc][longc] = 0.0;
					cc_sd [i][latgc][longc] = 0.0;

					cc_rd [i][latgc][longc] = 0.0;

					/* Seasonal Canopy parameters */

					sc_lai    [i][latgc][longc] = 0.0;
					sc_laid   [i][latgc][longc] = 0.0;

					sc_flb    [i][latgc][longc] = 0.0;

					sc_llb    [i][latgc][longc] = 0.0;
					sc_lln    [i][latgc][longc] = 0.0;
					sc_lwb    [i][latgc][longc] = 0.0;
					sc_lwn    [i][latgc][longc] = 0.0;

					sc_lo     [i][latgc][longc] = 0.0;

					sc_lrvi[i][0][latgc][longc] = 0.0;
					sc_lrvi[i][1][latgc][longc] = 0.0;
					sc_ltvi[i][0][latgc][longc] = 0.0;
					sc_ltvi[i][1][latgc][longc] = 0.0;

					sc_scmax  [i][latgc][longc] = 0.0;
					sc_scmin  [i][latgc][longc] = 0.0;
					sc_scpar  [i][latgc][longc] = 0.0;
					sc_scvpd  [i][latgc][longc] = 0.0;
				}
			}

			if (leafy && output_day)
			{
				animate_leaves();

				if (leafy2)
					for (i = 0; i < 2; i++)
						record_pc2(cfp[i].lai, 3, canopy.seasonal[i].lai);
			}

			if (image)
			{
				if (date == 1)
					record_cc();

				if (output_day)
					animate_sc();
			}
		}
	}

	close_leaf_files();
/*
	if (*first_day == 2)
		write_leaf_file(*year);

	if (date == mmd[0])

	if (date == 1)
	{
		write_inertial_vars();
		write_lf_variables();
	}
*/
	if (date == ldm[month-1])
		++month;

	if (global_monitor)
		monitor_leaf(0, "eve_leaf");

}



void	spin_up (dayn)

int	dayn;							/* spin-up the leaf phenology model */
										/* from day 1 to day n              */
{

	int	jj = 0;

	if (global_monitor)
		monitor_leaf(1, "spin_up");
/*
	printf("\n    Leaf spin-up: \n\n");
*/
	for (lat = 0; lat < n_lats; lat++)
	{
/*
		printf("             Lat = %3d \n", lat+1);
*/
		fflush(stdout);

		latitude = (90.0 - (deg_per_lat / 2.0)) - (((double) lat) * deg_per_lat);

		/* set starting point for DDs:  one month after the winter solstice */

		if (latitude >= 0.0)
			winter_day =  21;
		else
			winter_day = 203;

		for (lon = 0; lon < n_lons; lon++)
			if (earth[lat][lon] == '1')
			{
				initialization();

				read_data(jj);
/*
				if (!jj)
				{
					printf("       LFC (%5d: %3d, %3d) = %3d \n",
					 jj+1, lat+1, lon+1, leaf_lfc);

					printf("       Leaf-Fall Day0[0]     = %3d \n",
					 leaf_fall_day0[0]);
					printf("       Leaf-Fall Day0[1]     = %3d \n",
					 leaf_fall_day0[1]);
					printf("       Leaf-Fall Day0[2]     = %3d \n",
					 leaf_fall_day0[2]);
					printf("       Leaf-Fall Day0[3]     = %3d \n\n",
					 leaf_fall_day0[3]);

					exit(0);
				}
*/
				if (leaf_lfc)
				{
					set_pc_lfs();

					for (date = 1; date < dayn; date++)
					{
						calc_daily_climate();
						calc_inertial_climate();

						green_leaf(jj);
					}

					write_data(jj);
				}

				++jj;
			}
	}

	if (global_monitor)
		monitor_leaf(0, "spin_up");

}



void	initialization ()	

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "initialization");

	for (i = 0; i < num_lf; i++)
	{
		seasonal_lai[1][i] = 0.0;
		seasonal_lai[0][i] = 0.0;
	}
/*
	for (i = 0; i < 4; i++)
	{
		leaf_fall_day0[i] = 0;
		leaf_fall_day1[i] = 0;
	}
*/
	if (latitude >= 0.0)
		for (i = 0; i < num_lf; i++)
			lf_season[i] = winter0;
	else
		for (i = 0; i < num_lf; i++)
			lf_season[i] = summer;

	if (global_monitor)
		monitor_leaf(0, "initialization");

}



void	set_pc_lfs ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "set_pc_lfs");

	/* set predicted LF vector */

	for (i = 0; i < num_lf; i++)
		if (leaf_fc[i] > 0.0)
			leaf_lf[i] = 1;
		else
			leaf_lf[i] = 0;

	if (global_monitor)
		monitor_leaf(0, "set_pc_lfs");

}



void	zero_vectors ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "zero_vectors");

	if (leafy2)
		for (i = 0; i < 2; i++)
			canopy.seasonal[i].lai = 0.0;

	if (image)
	{
		for (i = 0; i < 2; i++)
		{
			canopy.constant[i].fc    = 0.0;
			canopy.constant[i].fcc   = 0.0;

			canopy.constant[i].z1    = 0.0;
			canopy.constant[i].z0    = 0.0;

			canopy.constant[i].sai   = 0.0;
			canopy.constant[i].sd    = 0.0;

			canopy.constant[i].rd    = 0.0;

			canopy.seasonal[i].lai   = 0.0;
			canopy.seasonal[i].laid  = 0.0;

			canopy.seasonal[i].flb   = 0.0;

			canopy.seasonal[i].llb   = 0.0;
			canopy.seasonal[i].lln   = 0.0;
			canopy.seasonal[i].lwb   = 0.0;
			canopy.seasonal[i].lwn   = 0.0;

			canopy.seasonal[i].lo    = 0.0;

			canopy.seasonal[i].lrv   = 0.0;
			canopy.seasonal[i].lri   = 0.0;
			canopy.seasonal[i].ltv   = 0.0;
			canopy.seasonal[i].lti   = 0.0;

			canopy.seasonal[i].scmax = 0.0;
			canopy.seasonal[i].scmin = 0.0;
			canopy.seasonal[i].scpar = 0.0;
			canopy.seasonal[i].scvpd = 0.0;
		}
	}

	if (global_monitor)
		monitor_leaf(0, "zero_vectors");

}



void	open_leaf_files (int year1)

{

	char	fn   [70],
			fn11 [70],
			fn22 [70],
			fnc  [70],
			fns  [70],

			time [7];

	int	i;

	if (global_monitor)
		monitor_leaf(1, "open_leaf_files");

	if (leafy && output_day)
		if (daily)
		{
			strcpy(fn11, "leaf_365/leaf1.");
			strcpy(fn22, "leaf_365/leaf0.");

			strcat(fn11, int_to_ascii(year1, time));
			strcat(fn22, int_to_ascii(year1, time));

			strcat(fn11, ".");
			strcat(fn22, ".");

			if (date < 10)
			{
				strcat(fn11, "00");
				strcat(fn22, "00");
			}
			else
				if (date < 100)
				{
					strcat(fn11, "0");
					strcat(fn22, "0");
				}

			strcpy(fns, int_to_ascii(date, time));

			strcat(fn11, fns);
			strcat(fn22, fns);
/*
			if (global_monitor)
*/
				printf("  Opening file:  %s \n", fn11);

			cfp1 = fopen(fn11, "w");
/*
			if (global_monitor)
*/
				printf("  Opening file:  %s \n", fn22);

			cfp2 = fopen(fn22, "w");
		}
		else
			if (date == mmd[month-1])
			{
				strcpy(fn11, "leaf_012/leaf.");
				strcpy(fn22, "leaf_012/leaf_22.");

				if (month < 10)
				{
					strcat(fn11, "0");
					strcat(fn22, "0");
				}

				strcpy(fns, int_to_ascii(month, time));

				strcat(fn11, fns);
				strcat(fn22, fns);

				if (global_monitor)
					printf("  Opening file:  %s \n", fn11);

				cfp1 = fopen(fn11, "w");

				if (global_monitor)
					printf("  Opening file:  %s \n", fn22);

				cfp2 = fopen(fn22, "w");
			}

	if (image)
	{
		if (daily)
			strcpy(fnc, "leaf_365/cc.");
		else
			strcpy(fnc, "leaf_012/cc.");

		if (date == 1)
			for (i = 0; i < 2; i++)
			{
				cfp[i].fc  = openff(filename(0, i, fnc, "fc"));
				cfp[i].fcc = openff(filename(0, i, fnc, "fcc"));

				cfp[i].z1  = openff(filename(0, i, fnc, "z1"));
				cfp[i].z0  = openff(filename(0, i, fnc, "z0"));

				cfp[i].sai = openff(filename(0, i, fnc, "sai"));
				cfp[i].sd  = openff(filename(0, i, fnc, "sd"));

				cfp[i].rd  = openff(filename(0, i, fnc, "rd"));
			}

		if (daily)
			strcpy(fns, "leaf_365/sc.");
		else
			if (date == mmd[month-1])
				strcpy(fns, "leaf_012/sc.");

		if (daily)
			for (i = 0; i < 2; i++)
			{
				cfp[i].lai   = openff(filename(1, i, fns, "lai"));
				cfp[i].laid  = openff(filename(1, i, fns, "laid"));
			}
		else
			if (date == mmd[month-1])
				for (i = 0; i < 2; i++)
				{
					cfp[i].lai   = openff(filename(2, i, fns, "lai"));
					cfp[i].laid  = openff(filename(2, i, fns, "laid"));
				}

		if (date == mmd[month-1])
			for (i = 0; i < 2; i++)
			{
				cfp[i].flb   = openff(filename(2, i, fns, "flb"));

				cfp[i].llb   = openff(filename(2, i, fns, "llb"));
				cfp[i].lln   = openff(filename(2, i, fns, "lln"));
				cfp[i].lwb   = openff(filename(2, i, fns, "lwb"));
				cfp[i].lwn   = openff(filename(2, i, fns, "lwn"));

				cfp[i].lo    = openff(filename(2, i, fns, "lo"));

				cfp[i].lrv   = openff(filename(2, i, fns, "lrv"));
				cfp[i].lri   = openff(filename(2, i, fns, "lri"));
				cfp[i].ltv   = openff(filename(2, i, fns, "ltv"));
				cfp[i].lti   = openff(filename(2, i, fns, "lti"));

				cfp[i].scmax = openff(filename(2, i, fns, "scmax"));
				cfp[i].scmin = openff(filename(2, i, fns, "scmin"));
				cfp[i].scpar = openff(filename(2, i, fns, "scpar"));
				cfp[i].scvpd = openff(filename(2, i, fns, "scvpd"));
			}
	}

	if (global_monitor)
		monitor_leaf(0, "open_leaf_files");

}



FILE	*openff (fn0)

char	*fn0;

{

	FILE	*fp;

	char	fn [70];

	strcpy(fn, fn0);

	fp = fopen(fn, "w");

	return(fp);

}



char	*filename (sc, cl, fnp, fns)

int	sc,
		cl;

char	fnp[],
		fns[];

{

	char	cs  [7],
			fn [70];

	if (global_monitor)
		monitor_leaf(1, "filename");

	strcpy(fn, fnp);
	strcat(fn, int_to_ascii(cl, cs));
	strcat(fn, ".");
	strcat(fn, fns);

	switch (sc)
	{
		case 0:
			break;

		case 1:			/* day name */

			strcat(fn, ".");

			if (date < 10)
				strcat(fn, "00");
			else
				if (date < 100)
					strcat(fn, "0");

			strcat(fn, int_to_ascii(date, cs));

			break;

		case 2:			/* month name */

			strcat(fn, ".");

			if (month < 10)
				strcat(fn, "0");

			strcat(fn, int_to_ascii(month, cs));

			break;
	}

	if (global_monitor)
	{
		monitor_leaf(0, "filename");
		printf("  Opening file:  %s \n", fn);
	}

	return(fn);

}



void	close_leaf_files ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "close_leaf_files");

	if (leafy && output_day)
	{
		fclose(cfp1);
		fclose(cfp2);
	}

	if (leafy2)
		for (i = 0; i < 2; i++)
			fclose(cfp[i].lai);

	if (image)
	{
		if (date == 1)
			for (i = 0; i < 2; i++)
			{
				fclose(cfp[i].fc);
				fclose(cfp[i].fcc);

				fclose(cfp[i].z1);
				fclose(cfp[i].z0);

				fclose(cfp[i].sai);
				fclose(cfp[i].sd);

				fclose(cfp[i].rd);
			}

		if (daily || (date == mmd[month-1]))
			for (i = 0; i < 2; i++)
			{
				fclose(cfp[i].lai);
				fclose(cfp[i].laid);
			}

		if (date == mmd[month-1])
			for (i = 0; i < 2; i++)
			{
				fclose(cfp[i].flb);

				fclose(cfp[i].llb);
				fclose(cfp[i].lln);
				fclose(cfp[i].lwb);
				fclose(cfp[i].lwn);

				fclose(cfp[i].lo);

				fclose(cfp[i].lrv);
				fclose(cfp[i].lri);
				fclose(cfp[i].ltv);
				fclose(cfp[i].lti);

				fclose(cfp[i].scmax);
				fclose(cfp[i].scmin);
				fclose(cfp[i].scpar);
				fclose(cfp[i].scvpd);
			}
	}

	if (global_monitor)
		monitor_leaf(0, "close_leaf_files");

}



void	record_cc ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "record_cc");

	for (i = 0; i < 2; i++)
	{
		record_pc2(cfp[i].fc,  4, canopy.constant[i].fc);
		record_pc2(cfp[i].fcc, 4, canopy.constant[i].fcc);
		record_pc2(cfp[i].z1,  2, canopy.constant[i].z1);
		record_pc2(cfp[i].z0,  2, canopy.constant[i].z0);
		record_pc2(cfp[i].sai, 4, canopy.constant[i].sai);
		record_pc2(cfp[i].sd,  4, canopy.constant[i].sd);
		record_pc2(cfp[i].rd,  3, canopy.constant[i].rd);
	}

	if (global_monitor)
		monitor_leaf(0, "record_cc");

}



void	animate_leaves ()

{

	if (global_monitor)
		monitor_leaf(1, "animate_leaves");

	fprintf(cfp1, "%5d", nint(canopy_lai  * 1000.0));
	fprintf(cfp2, "%5d", nint(canopy_laid * 1000.0));

	if (lon == 359)
	{
		fprintf(cfp1, "\n");
		fprintf(cfp2, "\n");
	}

	if (global_monitor)
		monitor_leaf(0, "animate_leaves");

}



void	animate_sc ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "animate_sc");

	for (i = 0; i < 2; i++)
	{
		record_pc2(cfp[i].lai,   3, canopy.seasonal[i].lai);
		record_pc2(cfp[i].laid,  3, canopy.seasonal[i].laid);
	}

	if (date == mmd[month-1])
		for (i = 0; i < 2; i++)
		{
			record_pc2(cfp[i].flb,   4, canopy.seasonal[i].flb);
			record_pc2(cfp[i].llb,   4, canopy.seasonal[i].llb);
			record_pc2(cfp[i].lln,   4, canopy.seasonal[i].lln);
			record_pc2(cfp[i].lwb,   4, canopy.seasonal[i].lwb);
			record_pc2(cfp[i].lwn,   4, canopy.seasonal[i].lwn);
			record_pc2(cfp[i].lo,    3, canopy.seasonal[i].lo);
			record_pc2(cfp[i].lrv,   4, canopy.seasonal[i].lrv);
			record_pc2(cfp[i].lri,   4, canopy.seasonal[i].lri);
			record_pc2(cfp[i].ltv,   4, canopy.seasonal[i].ltv);
			record_pc2(cfp[i].lti,   4, canopy.seasonal[i].lti);
			record_pc2(cfp[i].scmax, 7, canopy.seasonal[i].scmax);
			record_pc2(cfp[i].scmin, 7, canopy.seasonal[i].scmin);
			record_pc2(cfp[i].scpar, 2, canopy.seasonal[i].scpar);
			record_pc2(cfp[i].scvpd, 0, canopy.seasonal[i].scvpd);
		}

	if (global_monitor)
		monitor_leaf(0, "animate_sc");

}



void	record_pc2 (fp, sf, pcv)

FILE		*fp;

int		sf;

double	pcv;

{

	if (global_monitor)
		monitor_leaf(1, "record_pc2");
/*
	printf("  PCV = %7.5lf   SF = %1d \n", pcv, sf);
	printf("  PCV-> %5d \n", nint(pcv * pow(10.0, (double) sf)));
*/
	fprintf(fp, "%5d", nint(pcv * pow(10.0, (double) sf)));

	if (lon == (n_lons - 1))
		fprintf(fp, "\n");

	if (global_monitor)
		monitor_leaf(0, "record_pc2");

}



void	assign_ddf ()

{

	double	c1, c2,
	
				ip = 0.0,

				latd;

	int		i, j,

				lat1,
				lat2;

	if (global_monitor)
		monitor_leaf(1, "assign_ddf");
/*
	printf("    Assigning DDF data. \n");
*/
	for (i = 0; i < 4; i++)
	{
/*
		printf("    Summergreen Class:  %1d \n", i+1);
*/
		for (j = 0; j < 180; j++)
		{
			lf_day00[i][j] = ddf_vector[i][j].lfday00;
			dd_ave0 [i][j] = ddf_vector[i][j].ddave0;
		}
	}

	if (resolution == 1)
	{
		for (i = 0; i < 4; i++)
			for (j = 0; j < n_lats; j++)
			{
				lf_day0[i][j] = lf_day00[i][j];
				dd_ave [i][j] = dd_ave0 [i][j];
			}
	}
	else		/* interpolate to current resolution */
	{
		for (i = 0; i < 4; i++)
			for (j = 0; j < n_lats; j++)
			{
				latitude =
				 (90.0 - (deg_per_lat / 2.0)) - (((double) j) * deg_per_lat);

				if (latitude > 89.5)
				{
					lf_day0[i][j] = lf_day00[i][0];
					dd_ave [i][j] = dd_ave0 [i][0];
				}
				else
					if (latitude < -89.5)
					{
						lf_day0[i][j] = lf_day00[i][179];
						dd_ave [i][j] = dd_ave0 [i][179];
					}
					else
					{
						if ((latitude - floor(latitude)) == 0.5)
						{
							lat1 = (int) (89.5 - latitude);

							lf_day0[i][j] = lf_day00[i][lat1];
							dd_ave [i][j] = dd_ave0 [i][lat1];
						}
						else
						{
							latd = 89.5 - latitude;

							c1 = modf(latd, &ip);
							c2 = 1.0 - c1;

							lat1 = (int) latd;
							lat2 = lat1 + 1;

							lf_day0[i][j] = nint((c1 * lf_day00[i][lat1]) +
														(c2 * lf_day00[i][lat2]));

							dd_ave[i][j] = (c1 * dd_ave0[i][lat1]) +
												(c2 * dd_ave0[i][lat2]);
						}
					}
			}
	}

	if (global_monitor)
		monitor_leaf(0, "assign_ddf");

}



void  monitor_leaf (direction, fname)

int   direction;

char  *fname;

{

	switch (direction)
	{
		case 1:
			printf("Entering:  %s \n", fname);
			break;

		case 0:
			printf("Leaving:   %s \n", fname);
			break;
	}

	fflush(stdout);

}

