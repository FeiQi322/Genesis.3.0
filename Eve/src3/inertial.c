/* Inertia:  the Leaf model's climate code */



# include	"gv.h"				/* Global Variables file */
# include	<math.h>
# include	<stdio.h>



double	ddays [12]  = { 31., 28., 31., 30., 31., 30.,
								 31., 31., 30., 31., 30., 31. },

			ctc = 0.935507; 	/* Climate Time Constant:  ctc = exp(-(t1-t0)/tau) */
									/*  where (t1-t0) is the time step (1 day) and     */
									/*  tau is the climate response time (15 days).    */
int		winter_day;


void		add_gcm_climate (),
			calc_inertial_climate ();



void	calc_daily_climate ()

{

	double	prec0, prec1,
				temp0, temp1,

				t0,    t1,
				t01,				/* t01 = t0 + t1 = mdm1 - mdm0 */

				ii,				/* ii = (double) i             */

				fdm,				/*  first day of month 1       */
				mdm0,				/* middle day of month 0       */
				mdm1;				/* middle day of month 1       */

	int		i,

				m0, m1;

	if (global_monitor)
		monitor_leaf(1, "calc_daily_climate");

	ctc = exp(-1.0 / 15.0);

	/* find daily climate (possibly deltas) */

	m0 = 11;
	m1 =  0;

	mdm0 = -15.0;
	fdm  =   1.0;
	mdm1 =  16.0;

	t01 = 31.0;

	ii = (double) date;

	while (ii > mdm1)
	{
		m0 = (m0 + 1) % 12;
		m1 = (m1 + 1) % 12;

		fdm += ddays[m0];

		mdm0 = fdm - ((ddays[m0] + 1.0) / 2.0);
		mdm1 = fdm + ((ddays[m1] - 1.0) / 2.0);

		t01 = mdm1 - mdm0;
	}

	if (ii == mdm1)
	{
		leaf_temp = T[m1];
		leaf_prec = P[m1] / ddays[m1];
	}
	else
	{
		t0 = ii - mdm0;
		t1 = mdm1 - ii;

		temp0 = T[m0];
		temp1 = T[m1];

		prec0 = P[m0] / ddays[m0];
		prec1 = P[m1] / ddays[m1];

		leaf_temp = ((t1 * temp0) + (t0 * temp1)) / t01;
		leaf_prec = ((t1 * prec0) + (t0 * prec1)) / t01;
	}

	if (global_monitor)
		monitor_leaf(0, "calc_daily_climate");

}



void	add_gcm_climate (global_prec, global_temp)

double	global_prec [n_lats] [n_lons],
			global_temp [n_lats] [n_lons];

{

	leaf_prec += global_prec[lat][lon];
	leaf_temp += global_temp[lat][lon];

	if (leaf_prec < 0.0)
		leaf_prec = 0.0;

}



void	calc_inertial_climate ()

{

	double	atc;

	atc = 1.0 - ctc;

	/* calculate degree-day sum */

	if (date == winter_day)
		leaf_dd = 0.0;

	if (leaf_temp > 0.0)
		leaf_dd += leaf_temp;

	/* calculate inertial temperature and precipitation */

	leaf_itemp = (atc * leaf_temp) + (ctc * leaf_itemp);
	leaf_iprec = (atc * leaf_prec) + (ctc * leaf_iprec);
/*
	leaf_itemp = leaf_temp - ((leaf_temp - leaf_itemp) * ctc);
	leaf_iprec = leaf_prec - ((leaf_prec - leaf_iprec) * ctc);
*/

}


