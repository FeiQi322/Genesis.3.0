/* IRF:  Initial Read Functions library */



# include	"gv.h"
# include	<malloc.h>
# include	<math.h>
# include	<stdio.h>
# include	<string.h>



enum seasonal_cycle { winter = 0, spring, summer, autumn, winter0 };


extern
enum seasonal_cycle	lf_season [num_lf];


extern
double	seasonal_lai [2] [num_lf];


extern
struct lf_data
{

	char		lfn    [3],			/* Life-form number        */
				fc     [5],			/* Fractional cover        */
				lai    [6],			/* Leaf area index         */
				laid   [6],			/* Leaf area index (brown) */
				season [1];			/* Summergreen season      */

	struct
	lf_data	*ptr;

}	*active [iflpm],
	*new;


extern
struct eco_data
{

	char		iflpn [6],			/* Ice-free land point number */

				lat [5],				/* Latitude in degrees and minutes  */
				lon [6],				/* Longitude in degrees and minutes */

				prec [12] [5],		/* Climate used to drive EVE in creating  */
				temp [12] [5],		/*  the plant community in the LEAF file. */

				iprec [5],			/* Inertial-precipitation */
				itemp [5],			/* Inertial-temperature   */
				degd  [7],			/* Degree-days            */

				var [5],				/* Vegetation-area-ratio */

				lfc [3];				/* Life-form count */

	struct
	lf_data	*ptr;

}	eco [iflpm];



FILE	*fpc,
		*fpp,
		*fpt,

		*open_file ();


char	*int_to_ascii (),

		earth [n_lats] [n_lons],
		
		earth_surf_file  [77],
		global_directory [77],
		global_prec_file [77],
		global_temp_file [77],

		nres [5] [19] = { "",
								"1x1",
								"2x2",
								"30x30",
								"t31" };


int	high_precision,

		iflp = 0;


void	read_climate_and_pc (),
		write_climate_and_pc (),

		read_esd (),
		write_esd (),

		open_files (),
		close_files (),

		read_climate (),
		write_climate (),

		read_pc (),
		write_pc (),

		read_lf (),
		write_lf (),

		initialize_pcv (),
		initializer (),
		set_constants (),

		calc_coords (),
		skip_data (),
		zero_data ();



void	read_and_compress_data ()
 
{
 
	int   i, j, k;
 
	if (global_monitor)
		monitor_leaf(1, "read_and_compress_data");
/*
	printf("\n\n  IRF:  Initial Read Functions program \n\n\n");
*/
	initializer();

	if (input_files)
		initialize_pcv();

	initialize(1);
	read_climate_and_pc();
/*
	initialize(0);
	write_climate_and_pc();
 
	printf("\n\n  It is done !!!  \n\n\n ");
*/
	if (global_monitor)
		monitor_leaf(0, "read_and_compress_data");
 
}
 
 

void	read_climate_and_pc ()

{

	char	fn  [77],
			fnp [77],
			lfs [7];

	int	i, j,

			jj = 0,

			lf;

	if (global_monitor)
		monitor_leaf(1, "read_climate_and_pc");

	read_esd();

	open_files(1);

	printf("    Reading:  Climate and Plant Community description. \n\n");

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			if (earth[i][j] == '1')
			{
				read_climate();

				read_pc();

				set_constants(jj);

				int2str(eco[jj].iflpn, (jj + 1), 6, 0);

				++jj;
			}
			else
				skip_data(1);
		}

		skip_data(0);
	}

	close_files(0);

	printf("    Reading:  LF FCs. \n");

	strcpy(fnp, global_directory);
	strcat(fnp, "/lf_");

	for (lf = 1; lf <= num_lf; lf++)
	{
		printf("              LF = %3d \n", lf);

		strcpy(fn, fnp);

		if (lf < 10)
			strcat(fn, "00");
		else
			if (lf < 100)
				strcat(fn, "0");

		int_to_ascii(lf, lfs);
		strcat(fn, lfs);

		read_lf(lf, fn);
	}

	printf("\n");

	if (global_monitor)
		monitor_leaf(0, "read_climate_and_pc");

}



void	write_climate_and_pc ()

{

	char	fn [70],
			lfs [7];

	int	i, j,
	
			jj = 0,
			
			lf;

	if (global_monitor)
		monitor_leaf(1, "write_climate_and_pc");

	write_esd();

	printf("    Writing:  Climate and Plant Community description. \n");

	open_files(0);

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			if (earth[i][j] == '1')
			{
				read_data(jj);

				write_climate();

				write_pc();

				++jj;
			}
			else
				zero_data(1);
		}

		zero_data(0);
	}

	close_files();

	printf("    Writing:  LF FCs. \n");

	for (lf = 1; lf <= num_lf; lf++)
	{
		printf("              LF = %3d \n", lf);

		strcpy(fn, "Test/eco_leaf/lf_");

		if (lf < 10)
			strcat(fn, "00");
		else
			if (lf < 100)
				strcat(fn, "0");

		int_to_ascii(lf, lfs);
		strcat(fn, lfs);

		write_lf(lf, fn);
	}

	if (global_monitor)
		monitor_leaf(0, "write_climate_and_pc");

}



void	read_esd ()

{

	FILE	*fp;

	char	fn [77];

	int	i, j,

			esp = 0;			/* Earth Surface Points */

	if (global_monitor)
		monitor_leaf(1, "read_esd");

	printf("\n    Reading:  Earth surface data \n");

	fp = open_file(earth_surf_file);

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			fscanf(fp, "%1c", &earth[i][j]);

			if (earth[i][j] == '1')
			{
				calc_coords(i, j, iflp);

				++iflp;
			}

			++esp;
		}

		fskipc(fp, 1);
	}

	fclose(fp);

	printf("    Number of ice-free land points (maximum) = %7d (%7d) \n",
	 iflp, iflpm);
	printf("    Number of earth surface points           = %7d \n\n", esp);
	printf("    Fraction of ice-free land points         = %15.7lf \n\n",
	 (double) iflp / (double) esp);

	if (iflp > iflpm)
	{
		printf("    Memory allocation exceeded. \n");
		printf("    Set constant iflpm in vsf.h to at least %7d. \n", iflp);
		printf("    Exiting program. \n\n");

		exit(0);
	}

	if (global_monitor)
		monitor_leaf(0, "read_esd");

}



void	write_esd ()

{

	FILE	*fp;

	int	i, j;

	if (global_monitor)
		monitor_leaf(1, "write_esd");

	printf("    Writing:  Earth surface data \n");

	fp = fopen("Test/data/earth_file", "w");

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			fprintf(fp, "%1c", earth[i][j]);
		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	if (global_monitor)
		monitor_leaf(0, "write_esd");

}



void	open_files (mode)

int	mode;

{

	char	fn [77];

	if (global_monitor)
		monitor_leaf(1, "open_files");
/*
	printf("    Opening:  precipitation, temperature, and VAR files. \n");
*/
	if (mode)
	{
		fpp = open_file(global_prec_file);
		fpt = open_file(global_temp_file);

		strcpy(fn, global_directory);
		strcat(fn, "/sv_var");

		fpc = open_file(fn);
	}
	else
	{
		fpp = fopen("Test/climate/prec_file", "w");
		fpt = fopen("Test/climate/temp_file", "w");

		fpc = fopen("Test/eco_leaf/sv_var", "w");
	}

	if (global_monitor)
		monitor_leaf(0, "open_files");

}



void	close_files ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "close_files");
/*
	printf("    Closing:  precipitation, temperature, and VAR files. \n");
*/
	fclose(fpp);
	fclose(fpt);
	fclose(fpc);

	if (global_monitor)
		monitor_leaf(0, "close_files");

}



void	read_climate ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "read_climate");

	if (high_precision)
		for (i = 0; i < 12; i++)
		{
			fscanfd(fpp, 6, &P[i]);
			fscanfd(fpt, 6, &T[i]);

			P[i] *= 0.1;
			T[i] *= 0.01;
		}
	else
		for (i = 0; i < 12; i++)
		{
			fscanfd(fpp, 5, &P[i]);
			fscanfd(fpt, 5, &T[i]);

			T[i] *= 0.1;
		}

	fskipc(fpp, 1);
	fskipc(fpt, 1);

	if (global_monitor)
		monitor_leaf(0, "read_climate");

}



void	write_climate ()

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "write_climate");

	if (high_precision)
		for (i = 0; i < 12; i++)
		{
			fprintf(fpp, "%6d", nint( 10.0 * P[i]));
			fprintf(fpt, "%6d", nint(100.0 * T[i]));
		}
	else
		for (i = 0; i < 12; i++)
		{
			fprintf(fpp, "%5d", nint(       P[i]));
			fprintf(fpt, "%5d", nint(10.0 * T[i]));
		}

	fprintf(fpp, "\n");
	fprintf(fpt, "\n");

	if (global_monitor)
		monitor_leaf(0, "write_climate");

}



void	read_pc ()

{

	if (global_monitor)
		monitor_leaf(1, "read_pc");

	fscanfd(fpc, 4, &leaf_var);

	leaf_var *= 0.001;

	if (global_monitor)
		monitor_leaf(0, "read_pc");

}



void	write_pc ()

{

	if (global_monitor)
		monitor_leaf(1, "write_pc");

	fprintf(fpc, "%4d", nint(1000.0 * leaf_var));

	if (global_monitor)
		monitor_leaf(0, "write_pc");

}



void	read_lf (lf, fn)

int	lf;

char	*fn;

{

	FILE	*fp;
	
	int	i, j,

			jj = 0,
				
			lf_count = 0,

			lf_cover;

	if (global_monitor)
		monitor_leaf(1, "read_lf");

   fp = open_file(fn);

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			if (earth[i][j] == '1')
			{
				fscanfi(fp, 4, &lf_cover);

				if (lf_cover != 0)
				{
					new = (struct lf_data *) calloc (1, sizeof(struct lf_data));

					if (new == 0)
					{
						printf("Not enough memory for LF %d in Vector %d \n", lf, jj);
						exit(0);
					}

					/* add life form to the LF list at IFLP jj */

					if (active[jj] != 0)
						active[jj]->ptr = new;
					else
						eco[jj].ptr = new;

					active[jj] = new;

					/* increment the life-form count */

					lf_count = str2int(eco[jj].lfc, 3);

					int2str(eco[jj].lfc, ++lf_count, 3, 0);

					/* set the LF number and fractional cover */

					int2str(active[jj]->lfn, lf, 3, 0);
					int2str(active[jj]->fc, (10 * lf_cover), 5, 0);

					/* initialize LAI fields to 0 */

					int2str(active[jj]->lai,    0, 6, 0);
					int2str(active[jj]->laid,   0, 6, 0);
					int2str(active[jj]->season, 0, 1, 0);

					/* set the new pointer at the end of the linked LF list to 0 */

					active[jj]->ptr = (struct lf_data *) 0;
				}

				++jj;
			}
			else
				fskipc(fp, 4);
		}

		fskipc(fp, 1);
	}

	fclose(fp);

	if (global_monitor)
		monitor_leaf(0, "read_lf");

}



void	write_lf (lf, fn)

int	lf;

char	*fn;

{

	FILE	*fp;

	int	alf,

			i, j,

			jj = 0,
			
			lf_cover;

	if (global_monitor)
		monitor_leaf(1, "write_lf");

	fp = fopen(fn, "w");

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			if (earth[i][j] == '1')
			{
				if (active[jj] != 0)
				{
					alf = str2int(active[jj]->lfn, 3);

					if (alf == lf)
					{
						lf_cover = str2int(active[jj]->fc, 5);

						fprintf(fp, "%4d", nint(0.1 * (double) lf_cover));

						active[jj] = active[jj]->ptr;
					}
					else
						fprintf(fp, "   0");
				}
				else
					fprintf(fp, "   0");

				++jj;
			}
			else
				fprintf(fp, "   0");
		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	if (global_monitor)
		monitor_leaf(0, "write_lf");

}



void	set_constants (jj)

int	jj;

{

	int		alf, i;

	struct
	lf_data	*lf_ptr;

	if (global_monitor)
		monitor_leaf(1, "set_constants");

	for (i = 0; i < 12; i++)
	{
		int2str(eco[jj].prec[i], nint( 10.0 * P[i]), 5, 0);		/* 0.1 mm */
		int2str(eco[jj].temp[i], nint(100.0 * T[i]), 5, 1);		/* 0.01 C */
	}

	int2str(eco[jj].var, nint(10000.0 * leaf_var), 5, 0);

	if (global_monitor)
		monitor_leaf(0, "set_constants");

}



void	skip_data (mode)

int	mode;

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "skip_data");

	if (mode)
	{
		fskipl(fpp, 1);
		fskipl(fpt, 1);
		fskipc(fpc, 4);
	}
	else
	{
		fskipc(fpp, 1);
		fskipc(fpt, 1);
		fskipc(fpc, 1);
	}

	if (global_monitor)
		monitor_leaf(0, "skip_data");

}



void	zero_data (mode)

int	mode;

{

	int	i;

	if (global_monitor)
		monitor_leaf(1, "zero_data");

	if (mode)
	{
		for (i = 0; i < 12; i++)
		{
			fprintf(fpp, "    0");
			fprintf(fpt, "    0");
		}

		fprintf(fpp, "\n");
		fprintf(fpt, "\n");

		fprintf(fpc, "   0");
	}
	else
	{
		fprintf(fpp, "\n");
		fprintf(fpt, "\n");
		fprintf(fpc, "\n");
	}

	if (global_monitor)
		monitor_leaf(0, "zero_data");

}



void	calc_coords(lat, lon, jj)

int	lat, lon,

		jj;					/* IFLP number */

{

	char	coordi,
			coordj,
			
			lat_str [7],
			lon_str [7];

	int	i = 0,
			j = 0,
	
			deg_lat,
			deg_lon,

			min_lat,
			min_lon;

	latitude = (90.0 - (deg_per_lat / 2.0)) - (((double) lat) * deg_per_lat);

	longitude = (((double) lon) * deg_per_lon) - (180.0 - (deg_per_lon / 2.0));

	if (latitude >= 0.0)
	{
		coordi = 'N';
	}
	else
	{
		coordi = 'S';

		latitude = -latitude;
	}

	if (longitude < 0.0)
	{
		longitude = -longitude;

		coordj = 'W';
	}
	else
	{
		coordj = 'E';
	}

	deg_lat = (int) latitude;
	deg_lon = (int) longitude;

	min_lat = nint(60.0 * (latitude - (double) deg_lat));
	min_lon = nint(60.0 * (longitude - (double) deg_lon));

	/* build latitude character string in reverse */

	do
	{
		lat_str[i++] = (min_lat % 10) + '0';
	}
	while ((min_lat /= 10) > 0);

	while (i < 2)
		lat_str[i++] = '0';

	lat_str[i++] = coordi;

	do
	{
		lat_str[i++] = (deg_lat % 10) + '0';
	}
	while ((deg_lat /= 10) > 0);

	while (i < 5)
		lat_str[i++] = ' ';

	lat_str[i] = '\0';

	reverses(lat_str, i);

	/* build longitude character string in reverse */

	do
	{
		lon_str[j++] = (min_lon % 10) + '0';
	}
	while ((min_lon /= 10) > 0);

	while (j < 2)
		lon_str[j++] = '0';

	lon_str[j++] = coordj;

	do
	{
		lon_str[j++] = (deg_lon % 10) + '0';
	}
	while ((deg_lon /= 10) > 0);

	while (j < 6)
		lon_str[j++] = ' ';

	lon_str[j] = '\0';

	reverses(lon_str, j);

	/* assign lat and lon strings to IFLP data array */

	for (i = 0; i < 5; i++)
		eco[jj].lat[i] = lat_str[i];

	for (j = 0; j < 6; j++)
		eco[jj].lon[j] = lon_str[j];

}



void  initialize_pcv ()

{

	FILE	*fp0;

	if (global_monitor)
		monitor_leaf(1, "initialize_pcv");

	if ((fp0 = fopen("data/pcv_0", "r")) != NULL)
	{
		fscanf(fp0, "%*s%*s%*s%*s%*s%*s");

		fscanf(fp0, "%*s %*d");							/* site_dimensions       */
		fscanf(fp0, "%*s %*d");							/* evolution_flag        */
		fscanf(fp0, "%*s %*d");							/* tsd_data_flag         */
		fscanf(fp0, "%*s %*d");							/* fire_flag             */

		fscanf(fp0, "%*s %*d");							/* global_monitor        */
		fscanf(fp0, "%*s %*d");							/* program_control_flag  */
		fscanf(fp0, "%*s %*d");							/* ecv_mapping_flag      */

		fscanf(fp0, "%*s %*d");							/* fc_var_flag           */
		fscanf(fp0, "%*s %d", &high_precision);
		fscanf(fp0, "%*s %*d");							/* predictor_output_flag */
		fscanf(fp0, "%*s %*d");							/* monitor_ai_flag       */
		fscanf(fp0, "%*s %*d");							/* monitor_wb_flag       */
		fscanf(fp0, "%*s %*d");							/* monitor_energy        */
		fscanf(fp0, "%*s %*d");							/* vlm_flag              */

		fscanf(fp0, "%*s %*d");							/* transect              */
		fscanf(fp0, "%*s %*d");							/* site                  */

		fscanf(fp0, "%*s %s", earth_surf_file);
		fscanf(fp0, "%*s %*s");							/* continents_file       */

		fscanf(fp0, "%*s %s", global_prec_file);
		fscanf(fp0, "%*s %s", global_temp_file);
		fscanf(fp0, "%*s %*s");							/* global_relh_file      */

		fscanf(fp0, "%*s %*s");							/* single_site_directory */
		fscanf(fp0, "%*s %*s");							/* site_loop_directory   */
		fscanf(fp0, "%*s %s", global_directory);
		fscanf(fp0, "%*s %*s");							/* ecv_directory         */
		fscanf(fp0, "%*s %*s");							/* climate_directory     */

		fscanf(fp0, "%*s %*d");							/* output_mode           */
		fscanf(fp0, "%*s %*d");							/* output_destination    */
		fscanf(fp0, "%*s %*d");							/* restart_mode          */
/*
		strcat(global_directory, "/");
*/
		fclose(fp0);
	}
	else
	{
		printf("\n\n    File data/pcv_0 not found.       Using defaults. \n");
	}

	if (global_monitor)
		monitor_leaf(0, "initialize_pcv");

}



void	initializer () /* initial default values for pcv_0 */

{

	if (global_monitor)
		monitor_leaf(1, "initializer");

	/* Program Control Variables */

	high_precision = 1 ;

	strcpy(earth_surf_file,  "data/earth.");
	strcpy(global_prec_file, "climate/prec.");
	strcpy(global_temp_file, "climate/temp.");

	strcat(earth_surf_file,  nres[resolution]);
	strcat(global_prec_file, nres[resolution]);
	strcat(global_temp_file, nres[resolution]);

	strcpy(global_directory, "eco_leaf");

	if (global_monitor)
		monitor_leaf(0, "initializer");

}


