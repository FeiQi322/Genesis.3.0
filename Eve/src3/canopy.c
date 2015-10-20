/* Canopy:  the 2-level Canopy Model */



# include   "canopy.h"		/* EVE plant canopy parameters file                */
# include	"gv.h"			/* Global Variables file                           */
# include   "lfg.h"			/* Life-Form Geometry file                         */
# include   "lfpcl2.h"		/* Life-Form Plant Community Levels file (2-level) */
# include   "lsc.h"			/* Leaf Spectral Coefficients                      */
# include   "scc.h"			/* Stomatal Conductance Coefficients               */
# include	<math.h>
# include	<stdio.h>
# include	<string.h>


# define		conservation		1



enum canopy_levels	{ tree = 1, vine, epiphyte, shrub, grass };

enum stem_condition	{ herb = 0, wood, ngherb };


extern
enum canopy_levels	lf_pcl [num_lf];

extern
enum stem_condition	lf_stem [num_lf];


extern
struct lai_parameters
{
	int		deciduous,		/* Deciduousness:  0 -> evergreen                */
				leaf_type;		/*                 1 -> raingreen                */
									/*                 2 -> summergreen              */
	double	ip_max,			/*                 3 -> raingreen or summergreen */
				dd_max;

	int		day_0;

	double	lat_0,
				d0_l0,
				tau,

				ip_min,
				dd_min;

}	leaf [num_lf];


extern
double	alpha,	/* (V and E) tree lai displacement coefficient */

			lf_lai               [num_lf],
			lf_overlap           [num_lf],
			lf_shade_intolerance [num_lf],

			seasonal_lai     [2] [num_lf];



double	canopy_fc [2],

			fsun [num_lf],

			lvcf;					/* leaf volume conservation factor */



void	calc_fsuns (),
		calc_fc2 (),
		constant_canopies (),
		seasonal_canopies ();


void	create_2_canopies ()

{

	constant_canopies();
	seasonal_canopies();

}



void	constant_canopies ()

{

	double	canopy_fc_sum [2],
				canopy_leaf_v [2] [2],

				fci,
				lvg;

	int		i,
				cl;		/* canopy level */

	for (cl = 0; cl < 2; cl++)
	{
		canopy.constant[cl].fc  = 0.0;
		canopy.constant[cl].fcc = 0.0;
		canopy.constant[cl].z1  = 0.0;
		canopy.constant[cl].z0  = 0.0;
		canopy.constant[cl].sai = 0.0;
		canopy.constant[cl].sd  = 0.0;
		canopy.constant[cl].rd  = 0.0;
	}

	if (leaf_lfc)
	{
		calc_fsuns();
		calc_fc2();

		for (cl = 0; cl < 2; cl++)
		{
			canopy.constant[cl].fcc = canopy_fc[cl];

			canopy_fc_sum[cl]    = 0.0;
			canopy_leaf_v[cl][1] = 0.0;
			canopy_leaf_v[cl][0] = 0.0;
		}

		for (i = 0; i < num_lf; i++)
			if (leaf_lf[i])
			{
				cl = lf_pcl_2[i];

				fci = leaf_fc[i];

				switch (lf_pcl[i])
				{
					case tree:
					case shrub:
					case grass:

						canopy_fc_sum[cl] += fci;

						canopy.constant[cl].z1 += fci * lfg[i].z1;
						canopy.constant[cl].z0 += fci * lfg[i].z0;

					case vine:

						lvg = fci * lf_lai[i];

						/* estimate SAI (stem area index) for LFs with wood stems */

						if (lf_stem[i] == 1)
						{
							canopy_leaf_v[cl][1] += lvg;

							canopy.constant[cl].sai += fci * (0.10 * lf_lai[i]);

							canopy.constant[cl].sd  += lvg * lfg[i].sd;
						}

						canopy_leaf_v[cl][0] += lvg;

						canopy.constant[cl].rd += lvg * lfg[i].rd;

						break;

					case epiphyte:
						break;
				}
			}

		/* renormalize non-fc properties */

		for (cl = 0; cl < 2; cl++)
		{
			if (canopy_fc_sum[cl] > 0.0)
			{
				canopy.constant[cl].z1 /= canopy_fc_sum[cl];
				canopy.constant[cl].z0 /= canopy_fc_sum[cl];
			}

			if (canopy.constant[cl].fcc > 0.0)
				canopy.constant[cl].sai /= canopy.constant[cl].fcc;

			if (canopy_leaf_v[cl][1] > 0.0)
				canopy.constant[cl].sd /= (100.0 * canopy_leaf_v[cl][1]);

			if (canopy_leaf_v[cl][0] > 0.0)
				canopy.constant[cl].rd /= (100.0 * canopy_leaf_v[cl][0]);
		}

		/* fine-tune 2-level canopy properties */

		if (conservation)
		{
			if (canopy.constant[1].fc > 0.0)
				lvcf = (leaf_var * canopy.constant[1].fcc) /
										 canopy.constant[1].fc;
			else
				if (canopy.constant[0].fc > 0.0)
					lvcf = (leaf_var * canopy.constant[0].fcc) /
											 canopy.constant[0].fc;
				else
					lvcf = 1.0;

			if (lvcf != 1.0)
				for (cl = 0; cl < 2; cl++)
					canopy.constant[cl].sai *= lvcf;
		}
/*
		for (cl = 0; cl < 2; cl++)
			if (canopy.constant[cl].fc > 1.0)
				canopy.constant[cl].fc = 1.0;

		if (canopy.constant[0].z0 < 0.01)
			canopy.constant[0].z0 = 0.01;

		if (canopy.constant[0].z1 < (canopy.constant[0].z0 + 0.01))
			canopy.constant[0].z1 = canopy.constant[0].z0 + 0.01;

		if (canopy.constant[1].z0 < (canopy.constant[0].z1 + 0.01))
			canopy.constant[1].z0 = canopy.constant[0].z1 + 0.01;

		if (canopy.constant[1].z1 < (canopy.constant[1].z0 + 0.01))
			canopy.constant[1].z1 = canopy.constant[1].z0 + 0.01;
*/
	}

}



void	calc_fsuns ()

{

	double	fci,

				sun_alpha,
				sun_beta,
				sun_gamma,

				sunlit_fc  [2],
				sunlit_lai [2],
				sunlit_ovf [2];

	int		i,
				cl;		/* canopy level */

	for (cl = 0; cl < 2; cl++)
	{
		sunlit_fc [cl] = 0.0;
		sunlit_lai[cl] = 0.0;
		sunlit_ovf[cl] = 0.0;
	}

	for (i = 0; i < num_lf; i++)
	{
		fsun[i] = 0.0;

		if (leaf_lf[i])
		{
			cl = lf_pcl_2[i];

			fci = leaf_fc[i];

			switch (lf_pcl[i])
			{
				case tree:

					fsun[i] = 1.0;

					sunlit_fc [cl] += fci;
					sunlit_lai[cl] += fci * lf_lai[i];
					sunlit_ovf[cl] += fci * lf_overlap[i];

					break;

				case vine:
				case epiphyte:
					break;

				case shrub:

					if (!cl)
					{
						fsun[i] = 1.0;

						sunlit_fc [cl] += fci;
						sunlit_lai[cl] += fci * lf_lai[i];
						sunlit_ovf[cl] += fci * lf_overlap[i];
					}

					break;

				case grass:
					break;
			}
		}
	}

	/* renormalize non-fc properties */

	for (cl = 0; cl < 2; cl++)
		if (sunlit_fc[cl] > 0.0)
		{
			sunlit_lai[cl] /= sunlit_fc[cl];
			sunlit_ovf[cl] /= sunlit_fc[cl];
		}

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
			switch (lf_pcl[i])
			{
				case tree:
				case shrub:
				case grass:

					if (fsun[i] == 0.0)
					{
						cl = lf_pcl_2[i];

						sun_alpha = 1.0 - sunlit_fc[cl];

						sun_beta  = sunlit_ovf[cl] * sunlit_fc[cl];

						sun_gamma = lf_shade_intolerance[i] * sunlit_lai[cl];

						fsun[i] =
						 sun_alpha / (sun_alpha + (sun_beta * exp(-sun_gamma)));
					}
					break;

				case vine:
				case epiphyte:
					break;
			}
/*
	printf("\n\n");

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
			printf("     CL = %1d     LF = %3d     FSUN = %7.5lf \n",
			 lf_pcl_2[i], i+1, fsun[i]);
*/

}



void	calc_fc2 ()

{

	double	fc_alpha,
				fc_beta;

	int		i,
				cl;		/* canopy level */

	for (cl = 0; cl < 2; cl++)
		canopy_fc[cl] = 0.0;

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
		{
			cl = lf_pcl_2[i];

			canopy_fc[cl] += leaf_fc[i] * fsun[i];
/*
			printf("\n  LF = %3d   FCi: %7.5lf * FSun: %7.5lf = FC: %7.5lf",
			 i+1, leaf_fc[i], fsun[i], (leaf_fc[i]*fsun[i]));
*/
		}
/*
	printf("\n\n  FC1 = %7.5lf   FC0 = %7.5lf \n", canopy_fc[1], canopy_fc[0]);
*/
	if (canopy_fc[1] == 0.0)
	{
		if (canopy_fc[0] == 0.0)
		{
			canopy.constant[1].fc = 0.0;
			canopy.constant[0].fc = 0.0;
		}
		else
		{
			canopy.constant[1].fc = 0.0;
			canopy.constant[0].fc = leaf_var;
		}
	}
	else
		if (canopy_fc[0] == 0.0)
		{
			canopy.constant[1].fc = leaf_var;
			canopy.constant[0].fc = 0.0;
		}
		else
		{
			fc_alpha = canopy_fc[1] / canopy_fc[0];

			fc_beta = fc_alpha + 1.0;

			canopy.constant[1].fc = (fc_beta -
			 sqrt(pow(fc_beta, 2.0) - (4.0 * fc_alpha * leaf_var))) / 2.0;

			canopy.constant[0].fc = canopy.constant[1].fc / fc_alpha;
		}

}



void	seasonal_canopies ()

{

	double	bf,
				nf,

				fci,				/* fractional cover,  LF i */
				lvd,				/*  dead leaf volume, LF i */
				lvg,				/* green leaf volume, LF i */
				lvi;				/* total leaf volume, LF i */

	int		cl,				/* canopy level */
				i;

	for (cl = 0; cl < 2; cl++)
	{
		canopy.seasonal[cl].lai   = 0.0;
		canopy.seasonal[cl].laid  = 0.0;

		canopy.seasonal[cl].flb   = 0.0;

		canopy.seasonal[cl].llb   = 0.0;
		canopy.seasonal[cl].lln   = 0.0;
		canopy.seasonal[cl].lwb   = 0.0;
		canopy.seasonal[cl].lwn   = 0.0;

		canopy.seasonal[cl].lo    = 0.0;

		canopy.seasonal[cl].lrv   = 0.0;
		canopy.seasonal[cl].lri   = 0.0;
		canopy.seasonal[cl].ltv   = 0.0;
		canopy.seasonal[cl].lti   = 0.0;

		canopy.seasonal[cl].scmax = 0.0;
		canopy.seasonal[cl].scmin = 0.0;
		canopy.seasonal[cl].scpar = 0.0;
		canopy.seasonal[cl].scvpd = 0.0;
	}

	if (leaf_lfc)
	{
		for (i = 0; i < num_lf; i++)
			if (leaf_lf[i])
			{
				cl = lf_pcl_2[i];

				fci = leaf_fc[i];

				if (lf_pcl[i] == tree)
				{
					lvg = fci * seasonal_lai[1][i] * alpha;
					lvd = fci * seasonal_lai[0][i] * alpha;
				}
				else
				{
					lvg = fci * seasonal_lai[1][i];
					lvd = fci * seasonal_lai[0][i];
				}

				lvi = lvg + lvd;

				canopy.seasonal[cl].lai  += lvg;
				canopy.seasonal[cl].laid += lvd;

				canopy.seasonal[cl].lo += lvi * lfg[i].lo;

				canopy.seasonal[cl].lrv += lvg * lf_lsc[i].alpha_vis;
				canopy.seasonal[cl].lri += lvg * lf_lsc[i].alpha_nir;
				canopy.seasonal[cl].ltv += lvg * lf_lsc[i].tau_vis;
				canopy.seasonal[cl].lti += lvg * lf_lsc[i].tau_nir;

				canopy.seasonal[cl].scmax += lvg * lf_scc[i].scmax;
				canopy.seasonal[cl].scmin += lvg * lf_scc[i].scmin;
				canopy.seasonal[cl].scpar += lvg * lf_scc[i].scpar;
				canopy.seasonal[cl].scvpd += lvg * lf_scc[i].scvpd;

				if (!leaf[i].leaf_type)
				{
					canopy.seasonal[cl].flb += lvi;

					canopy.seasonal[cl].llb += lvi * lfg[i].ll;
					canopy.seasonal[cl].lwb += lvi * lfg[i].lw;
				}
				else
				{
					canopy.seasonal[cl].lln += lvi * lfg[i].ll;
					canopy.seasonal[cl].lwn += lvi * lfg[i].lw;
				}
/*
				if (i == 6)
					printf("\n\n");

				printf("  CL = %1d LF = %3d lvg:%8.4lf scmax:%8.5lf scvpd:%8.4lf\n",
				 cl, i+1, lvg, lf_scc[i].scmax, lf_scc[i].scvpd);
*/
			}

		/* renormalize non-fc properties */

		for (cl = 0; cl < 2; cl++)
		{
			lvg = canopy.seasonal[cl].lai;
			lvd = canopy.seasonal[cl].laid;

			lvi = lvg + lvd;

			if (canopy.constant[cl].fcc > 0.0)
			{
				canopy.seasonal[cl].lai  /= canopy.constant[cl].fcc;
				canopy.seasonal[cl].laid /= canopy.constant[cl].fcc;
			}

			if (lvi > 0.0)
			{
				bf = canopy.seasonal[cl].flb;

				nf = lvi - bf;

				canopy.seasonal[cl].flb /= lvi;

				if (bf > 0.0)
				{
					canopy.seasonal[cl].llb /= (100.0 * bf);
					canopy.seasonal[cl].lwb /= (100.0 * bf);
				}

				if (nf > 0.0)
				{
					canopy.seasonal[cl].lln /= (100.0 * nf);
					canopy.seasonal[cl].lwn /= (100.0 * nf);
				}

				canopy.seasonal[cl].lo /= lvi;
			}

			if (lvg > 0.0)
			{
				canopy.seasonal[cl].lrv /= lvg;
				canopy.seasonal[cl].lri /= lvg;
				canopy.seasonal[cl].ltv /= lvg;
				canopy.seasonal[cl].lti /= lvg;

				canopy.seasonal[cl].scmax /= (1000.0 * lvg);
				canopy.seasonal[cl].scmin /= (1000.0 * lvg);
				canopy.seasonal[cl].scpar /= lvg;
				canopy.seasonal[cl].scvpd *= (100.0 / lvg);
			}
/*
			if (!cl)
				printf("\n");

			printf("  CL = %1d          LVG:%8.4lf SCMAX:%8.5lf SCVPD:%8.4lf\n",
			 cl, lvg, canopy.seasonal[cl].scmax * 1000.0,
						 canopy.seasonal[cl].scvpd / 100.0);
*/
		}

		if (conservation && (lvcf != 1.0))
			for (cl = 0; cl < 2; cl++)
			{
				canopy.seasonal[cl].lai  *= lvcf;
				canopy.seasonal[cl].laid *= lvcf;
			}
	}

}

