/* EVE Plant Community Upper and Lower Canopy Parameters (version 1) */



# include   "vsf.h"



typedef
double	dcc;



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


struct canopies		/* level 1: upper canopy,  level 0: lower canopy */
{
	struct constant_canopy	constant [2];
	struct seasonal_canopy	seasonal [2];

}	canopy;
