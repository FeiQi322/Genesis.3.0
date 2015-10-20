/* Leaf:  the Leaf Area Index (LAI) Model:  life-form level */



# include	"gv.h"				/* Global Variables file                   */
# include	"lai.h"				/* Leaf Area Indices file                  */
# include	"leaf.h"				/* Leaf Area Index parameters file         */
# include	"leaf_type.h"		/* Summergreen LFs Leaf-Fall type file     */
# include	"med.h"				/* Mediterranean (nonsummergreen) LFs file */
# include	"stem.h"				/* Life-Form Stem condition file           */
# include	"summergreen.h"	/* Summergreen LFs LAI parameters file     */
# include	<math.h>
# include	<stdio.h>


# define		test_day		500
# define		test_lf		500



enum canopy_levels	{ tree = 1, vine, epiphyte, shrub, grass };


extern
enum canopy_levels	lf_pcl [num_lf];


extern
int		winter_day;



enum
seasonal_cycle { winter = 0, spring,
									  summer,
									  autumn,
									  winter0 }	lf_season [num_lf];


char		lfd0 [iflpm] [4] [3];


double	seasonal_lai [2] [num_lf],
			
			ddpcd [4] = { 52.0, 39.0, 40.0, 75.0 };


int		leaf_fall_day0 [4],
			leaf_fall_day1 [4],

			lf_deciduous [num_lf];


void		calc_leaf_days ();



void	green_leaf (jj)

int	jj;

{

	double	ltc = 0.818731, /* Leaf Time Constant:  ltc = exp(-(t1-t0)/tau) */
									 /*  where (t1-t0) is the time step (1 day) and  */
									 /*  tau is the leaf residence time (5 days).    */

				htc = 0.818731, /* Herb Time Constant:  htc = exp(-(t1-t0)/tau)  */
									 /*  where (t1-t0) is the time step (1 day) and   */
									 /*  tau is the herb leaf response time (5 days). */

				rtc = 0.818731, /* RG Time Constant:  rtc = exp(-(t1-t0)/tau)  */
									 /*  where (t1-t0) is the time step (1 day) and */
									 /*  tau is the leaf residence time (5 days).   */

				stc = 0.716531, /* Spring Time Constant:  stc = exp(-(t1-t0)/tau) */
									 /*  where (t1-t0) is the time step (1 day) and    */
									 /*  tau is the spring leaf response time (3 days).*/

				dd_range,

				delta_lai,

				iprec_diff,

				itemp_diff,
				itemp_max = 15.0,
				itemp_min = 10.0,

				lai_0,
				lai_eq,

				lai_factor,
				lai_itemp_factor,

				lff,
				lf_tau,

				mediterranean_fc,
				raingreen_fc,
				summergreen_fc;

	int		i, j,

				deciduousness,
				lft;

	if (global_monitor)
		monitor_leaf(1, "green_leaf");

	/* Assign multi-modal LFs to raingreen, mediterranean raingreen, or */
	/* summergreen state.  If any mediterranean LFs are present, all    */
	/* multi-modal and raingreen LFs are mediterranean raingreen.       */
	/* Non-mediterranean raingreen LFs are also winter deciduous.       */

	/* Leaf file deciduousness:  0 -> evergreen                         */
	/*                           1 -> raingreen                         */
	/*                           2 -> summergreen                       */
	/*                           3 -> raingreen or summergreen          */

	mediterranean_fc = 0.0;
	raingreen_fc     = 0.0;
	summergreen_fc   = 0.0;

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
		{
			if (mediterranean[i])
				mediterranean_fc += leaf_fc[i];

			switch (leaf[i].deciduous)
			{
				case 0:		/* evergreen */

					break;

				case 1:		/* raingreen */

					raingreen_fc += leaf_fc[i];

					break;

				case 2:		/* summergreen */

					summergreen_fc += leaf_fc[i];

					break;

				case 3:		/* raingreen or summergreen */

					break;
			}
		}

	if (mediterranean_fc > 0.0)
	{
		deciduousness = 1;								/* raingreen   */
	}
	else
	{
		if ((raingreen_fc > 0.0) || (summergreen_fc > 0.0))
		{
			if (raingreen_fc >= summergreen_fc)
				deciduousness = 1;						/* raingreen   */
			else
				deciduousness = 2;						/* summergreen */
		}
		else
		{
			if (leaf_dlt < 10.0)
				deciduousness = 1;						/* raingreen   */
			else
				deciduousness = 2;						/* summergreen */
		}
	}

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
		{
			if (leaf[i].deciduous == 3)
				lf_deciduous[i] = deciduousness;
			else
				lf_deciduous[i] = leaf[i].deciduous;
		}

	/* calculate daily inertial temperature factors for raingreen LFs */

	itemp_diff = itemp_max - itemp_min;

	if (mediterranean_fc == 0.0)
	{
		if (leaf_itemp >= itemp_max)
			lai_itemp_factor = 1.0;
		else
			if (leaf_itemp <= itemp_min)
				lai_itemp_factor = 0.0;
			else
				lai_itemp_factor = (leaf_itemp - itemp_min) / itemp_diff;
	}
	else
		lai_itemp_factor = 1.0;

	/* Calculate daily LAIs for each life form */

	if (date == test_day)
	{
		printf("\n\n    LF %1d:  plf = %1d, lfd = %1d, ",
		 test_lf+1, leaf_lf[test_lf], lf_deciduous[test_lf]);

		printf("wd = %1d, lf_season = %3d \n", winter_day, lf_season[test_lf]);
	}

	if ((summergreen_fc > 0.0) || (deciduousness == 2))
		calc_leaf_days(jj);

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
		{
			switch (lf_deciduous[i])
			{
				case 0:		/* evergreen */

					seasonal_lai[1][i] = lf_lai[i];
					seasonal_lai[0][i] = 0.0;

					break;

				case 1:		/* raingreen */

					iprec_diff = leaf[i].ip_max - leaf[i].ip_min;

					if (leaf_iprec >= leaf[i].ip_max)
						lai_factor = lai_itemp_factor;
					else
						if (leaf_iprec <= leaf[i].ip_min)
							lai_factor = 0.0;
						else
						{
							lff = (leaf_iprec - leaf[i].ip_min) / iprec_diff;

							lai_factor = lai_itemp_factor * lff;

							if (lf_stem[i])
							{
								lf_tau = (20.0 * lff) + 5.0;

								rtc = exp(-1.0 / lf_tau);
							}
						}

					lai_eq = lf_lai[i] * lai_factor;

					lai_0 = seasonal_lai[1][i];

					seasonal_lai[1][i] = lai_eq - ((lai_eq - lai_0) * stc);

					if (seasonal_lai[1][i] < 0.00001)
						seasonal_lai[1][i] = 0.0;

					if (lf_stem[i])
					{
						delta_lai = lai_0 - seasonal_lai[1][i];

						if (delta_lai < 0.0)
						{
							delta_lai = 0.0;

							rtc = 0.818731;		/* lf_tau = 5 days */
						}

						seasonal_lai[0][i] *= rtc;
						seasonal_lai[0][i] += delta_lai;

						if (seasonal_lai[0][i] < 0.00001)
							seasonal_lai[0][i] = 0.0;
					}

					break;

				case 2:		/* summergreen */

					lft = lf_type[i] - 1;

					dd_range = leaf[i].dd_max - leaf[i].dd_min;

					/* calculate summergreen LAIs */

					if (date == winter_day)
						lf_season[i] = winter;

					if ((test_day <= 365) && (i == test_lf))
						printf("\n    Day %3d: lf_season = %3d ", date, lf_season[i]);

					if (lf_season[i] != winter0)
						switch (lf_season[i])
						{
							case winter:

								/* check for absence of growing season */

								if (leaf_fall_day0[lft] && (date > leaf_fall_day0[lft]))
									lf_season[i] = winter0;
								else
									if (leaf_dd > leaf[i].dd_min)
										lf_season[i] = spring;

								break;

							case spring:

								/* check for short growing season */

								if (leaf_fall_day0[lft] && (date > leaf_fall_day0[lft]))
									lf_season[i] = autumn;
								else
									if (leaf_dd >= leaf[i].dd_max)
										lf_season[i] = summer;

								break;

							case summer:

								if (leaf_fall_day0[lft] && (date > leaf_fall_day0[lft]))
									lf_season[i] = autumn;

								break;

							case autumn:

								if (date >= leaf_fall_day1[lft])
									lf_season[i] = winter0;

								break;
						}

					if ((test_day <= 365) && (i == test_lf))
						printf("-> %3d ", lf_season[i]);

					switch (lf_season[i])
					{
						case winter0:
						case winter:

							lai_factor = 0.0;

							break;

						case spring:

							lai_factor = (leaf_dd - leaf[i].dd_min) / dd_range;

							break;

						case summer:

							lai_factor = 1.0;

							break;

						case autumn:

							lai_factor =
							 ((double) (leaf_fall_day1[lft] - date)) / leaff[lft].tau;

							if (lf_stem[i])
							{
								lf_tau = (95.0 * lai_factor) + 5.0;

								ltc = exp(-1.0 / lf_tau);
							}

							break;
					}

					lai_eq = lf_lai[i] * lai_factor;

					lai_0 = seasonal_lai[1][i];
/*
					if ((lf_season[i] == autumn) && (lai_eq > lai_0))
						lai_eq = lai_0;
*/
					seasonal_lai[1][i] = lai_eq - ((lai_eq - lai_0) * stc);

					if (seasonal_lai[1][i] < 0.00001)
						seasonal_lai[1][i] = 0.0;

					if (date == test_day)
					{
						printf("\n\n    lf_season (day = %3d) = %3d \n",
						 test_day, lf_season[i]);
						printf("    LF (%3d) lai     = %9.5lf \n",
						 i+1, seasonal_lai[1][i]);

						printf("    tree lai_eq      = %9.5lf \n", lai_eq);
						printf("    tree lai_0       = %9.5lf \n", lai_0);
						printf("    stc              = %9.5lf \n", stc);
						printf("    dd_max           = %9.5lf \n", leaf[i].dd_max);
						printf("    dd               = %9.5lf \n", leaf_dd);
						printf("    dd_min           = %9.5lf \n", leaf[i].dd_min);
						printf("    leaf_fall_day1   = %3d \n", leaf_fall_day1[lft]);
						printf("    date             = %3d \n", date);
						printf("    leaf_fall_day0   = %3d \n", leaf_fall_day0[lft]);
						printf("    lai_factor       = %9.5lf \n", lai_factor);
						printf("    lai_itemp_factor = %9.5lf \n",
						 lai_itemp_factor);
					}

					if (lf_stem[i])
					{
						delta_lai = lai_0 - seasonal_lai[1][i];

						if (delta_lai < 0.0)
							delta_lai = 0.0;

						seasonal_lai[0][i] *= ltc;
						seasonal_lai[0][i] += delta_lai;

						if (seasonal_lai[0][i] < 0.00001)
							seasonal_lai[0][i] = 0.0;
					}

					if (date == test_day)
					{
						printf("\n    LF (%3d) laid    = %9.5lf \n",
						 i+1, seasonal_lai[0][i]);

						if (lf_stem[i])
						{
							printf("    delta_lai        = %9.5lf \n", delta_lai);
							printf("    LF laid_0        = %9.5lf \n", lai_0);
						}
					}

					break;

				case 3:		/* raingreen or summergreen */

					printf("\n\n  Undefined deciduous state. \n\n");
					exit(0);

					break;
			}

			/* Herbs have seasonal and background */
			/*  components of LAId (dead leaves). */

			if (!lf_stem[i])
			{
				seasonal_lai[0][i] = (0.5 * lf_lai[i]) - (0.4 * seasonal_lai[1][i]);
			}
		}
		else
		{
			seasonal_lai[1][i] = 0.0;
			seasonal_lai[0][i] = 0.0;
		}

	if (global_monitor)
		monitor_leaf(0, "green_leaf");

}



void	calc_leaf_days (jj)

int	jj;

{

	double	dd,
				dd0,
				dd_pcd;

	int		i,
				latj,
				lfdd,
				time_shift;

	if (global_monitor)
		monitor_leaf(1, "calc_leaf_days");

	if (site_dimensions < 2)
		latj = nint(((90.0 - (deg_per_lat / 2.0)) - latitude) / deg_per_lat);
	else
		latj = lat;

	if ((date == 1) || (date == winter_day))
	{
		for (i = 0; i < 4; i++)
			leaf_fall_day0[i] = 0;
	}
	else
	{
		for (i = 0; i < 4; i++)
		{
			leaf_fall_day0[i] = str2int(lfd0[jj][i], 3);

			if (!leaf_fall_day0[i])
			{
				lfdd = lf_day0[i][latj] - 15;

				if (date == lfdd)
				{
					dd  = leaf_dd;
					dd0 = dd_ave[i][latj];

					dd_pcd = 100.0 * (dd - dd0) / dd0;

					time_shift = nint(7.0 * (dd_pcd / ddpcd[i]));

					if (time_shift > 7)
						time_shift = 7;
					else
						if (time_shift < -7)
							time_shift = -7;

					leaf_fall_day0[i] = lf_day0[i][latj] + time_shift;
					leaf_fall_day1[i] = leaf_fall_day0[i] + (int) leaff[i].tau;
				}
			}
			else
				leaf_fall_day1[i] = leaf_fall_day0[i] + (int) leaff[i].tau;
		}
	}

	for (i = 0; i < 4; i++)
		int2str(lfd0[jj][i], leaf_fall_day0[i], 3, 0);

	if (global_monitor)
		monitor_leaf(0, "calc_leaf_days");

}

