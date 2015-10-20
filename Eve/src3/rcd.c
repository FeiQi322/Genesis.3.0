/* RCD:  read climate deltas */



# include	"gv.h"
# include	<stdio.h>


# define		high_output		0



extern
char	earth [n_lats] [n_lons];



char	fnp [5] [19] = { "dat/dat_",
							  "ctl/ctl_",
							  "exp/exp_",
							  "climate/",
							  "cds/cds_" },

		stn [6] [19] = { "Data            ",
							  "Control         ",
							  "Experiment      ",
							  "EVE Climate     ",
							  "Adjusted Climate",
							  "Inertial Climate" },

		svn [3] [19] = { "prec.",
							  "temp.",
							  "relh." },

		svt [3] [19] = { "precipitation",
							  "temperature",
							  "relative humidity" },

		res [5] [19] = { "",
							  "1x1",
							  "2x2",
							  "30x30",
							  "t31" };


struct delta_data
{

	char	prec [12] [6],			/* Climate deltas = data - control */
			temp [12] [5];			/*  filter applied to GCM climate  */

}	delta [iflpm];


void	set_deltas ();



void	read_deltas ()

{

	FILE	*fp [4];

	char	fn [77];

	int	i, j, k,
	
			jj = 0,
			
			precd,			/* precipitation, data    */
			precc,			/* precipitation, control */
			tempd,			/* temperature, data      */
			tempc;			/* temperature, control   */

	if (global_monitor)
		monitor_leaf(1, "read_deltas");

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
		{
			strcpy(fn, fnp[j]);
			strcat(fn, svn[i]);
			strcat(fn, res[resolution]);

			if (high_output)
				printf("      Opening file:  %s \n", fn);

			fp[(2*i)+j] = open_file(fn);
		}
	}

	for (i = 0; i < n_lats; i++)			/* prec and temp data is stored at */
		for (j = 0; j < n_lons; j++)		/*  0.1 x precision of the control */
		{
			if (earth[i][j] == '1')			/*  _||_  */
			{										/*  \  /  */
				for (k = 0; k < 12; k++)	/*   \/   */
				{
					fscanf(fp[0], "%d", &precd);
					fscanf(fp[1], "%d", &precc);
					fscanf(fp[2], "%d", &tempd);
					fscanf(fp[3], "%d", &tempc);

/* dp mod 8/8/00 dat precision same as ctl */
				     /* int2str(delta[jj].prec[k], ((10 * precd) - precc), 6, 1); */
				     /* int2str(delta[jj].temp[k], ((10 * tempd) - tempc), 5, 1); */
					int2str(delta[jj].prec[k], (      precd  - precc), 6, 1);
					int2str(delta[jj].temp[k], (      tempd  - tempc), 5, 1); 
				}

				++jj;
			}
			else
			{
				fskipl(fp[0], 1);
				fskipl(fp[1], 1);
				fskipl(fp[2], 1);
				fskipl(fp[3], 1);
			}
		}

	fclose(fp[0]);
	fclose(fp[1]);
	fclose(fp[2]);
	fclose(fp[3]);

	if (global_monitor)
		monitor_leaf(0, "read_deltas");

}



void	set_deltas (jj)

int	jj;

{

	int	i;

	for (i = 0; i < 12; i++)
	{
		P[i] = 0.10 * (double) str2int(delta[jj].prec[i], 6);	/* 0.1 mm */
		T[i] = 0.01 * (double) str2int(delta[jj].temp[i], 5);	/* 0.01 C */
	}

}


