/* Global Variables for LEAF */



# include	"vsf.h"			/* Vector Size File */
# include	<stdio.h>



/* Prototypes */

extern
FILE		*open_file (),
			*open_write ();

extern
char		*int_to_ascii ();

extern   
double	anint ();

extern
int		nint ();



extern
double	P [12],
			T [12],
			H [12],

			latitude,
			longitude,

			deg_per_lat,
			deg_per_lon,

			leaf_fc [num_lf],

			leaf_var,

			leaf_dd,
			leaf_dlt,
			leaf_iprec,
			leaf_itemp,
			leaf_prec,
			leaf_temp,

			temp_max,
			temp_min,

			dd_ave  [4] [n_lats],
			dd_ave0 [4] [180];


extern
int		date,

			leaf_lf [num_lf],
			leaf_lfc,

			lf_day0  [4] [n_lats],
			lf_day00 [4] [180],

			global_monitor,

			interactive,
			filtered,

			lat,   lon,
			latgc, longc,

			site_dimensions;


