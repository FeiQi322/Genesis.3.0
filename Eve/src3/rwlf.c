/* RWLF:  Read and Write LEAF File function library */



# include	"gv.h"
# include	<malloc.h>
# include	<math.h>
# include	<stdio.h>
# include	<string.h>


# define		version "1.00"



enum seasonal_cycle { winter = 0, spring, summer, autumn, winter0 };


extern
enum seasonal_cycle	lf_season [num_lf];


extern
double	seasonal_lai [2] [num_lf];



FILE		*open_file ();


char		*int_to_ascii (),

			surfs [3] = { '_', '/', '\\' },

			earth [n_lats] [n_lons];


int		iflp,

			nlats,
			nlons;


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


void	initialize (),

		read_data (),
		write_data (),

		read_earth (),
		write_earth (),

		read_leaf_file (),
		write_leaf_file ();



void	test_main ()

{

	int	i;

	read_leaf_file();

	for (i = 0; i < iflp; i++)
	{
		read_data(i);
		write_data(i);
	}

	write_leaf_file(2001);

	printf("\n\n    It is done !!! \n\n\n");

}



void	read_data (jj)

int	jj;

{

	int		alf, i;

	struct
	lf_data	*lf_ptr;

	if (global_monitor)
		monitor_leaf(1, "read_data");

	for (i = 0; i < 12; i++)
	{
		P[i] = 0.10 * (double) str2int(eco[jj].prec[i], 5);		/* 0.1 mm */
		T[i] = 0.01 * (double) str2int(eco[jj].temp[i], 5);		/* 0.01 C */
	}

	leaf_iprec = 0.1000 * (double) str2int(eco[jj].iprec, 5);	/* 0.1 mm     */
	leaf_itemp = 0.0100 * (double) str2int(eco[jj].itemp, 5);	/* 0.01 C     */
	leaf_dd    = 0.0100 * (double) str2int(eco[jj].degd,  7);	/* 0.1 C-days */
	leaf_var   = 0.0001 * (double) str2int(eco[jj].var,   5);
	leaf_lfc   =                   str2int(eco[jj].lfc,   3);

	for (i = 0; i < num_lf; i++)
	{
		leaf_fc        [i] = 0.0;
		seasonal_lai[1][i] = 0.0;
		seasonal_lai[0][i] = 0.0;
		lf_season      [i] = 0;
	}

	lf_ptr = eco[jj].ptr;

	while (lf_ptr != 0)
	{
		alf = str2int(lf_ptr->lfn, 3) - 1;

		leaf_fc        [alf] = 0.0001 * (double) str2int(lf_ptr->fc,     5);
		seasonal_lai[1][alf] = 0.0001 * (double) str2int(lf_ptr->lai,    6);
		seasonal_lai[0][alf] = 0.0001 * (double) str2int(lf_ptr->laid,   6);
		lf_season      [alf] =                   str2int(lf_ptr->season, 1);

		lf_ptr = lf_ptr->ptr;
	}

	/* find delta_T = T_max - T_min */
   
	temp_max = T[0];
 
	for (i = 1; i < 12; i++)
		if (temp_max < T[i])
			temp_max = T[i];
 
	temp_min = T[0];
 
	for (i = 1; i < 12; i++)
		if (temp_min > T[i])
			temp_min = T[i];
 
	leaf_dlt = temp_max - temp_min;

	if (global_monitor)
		monitor_leaf(0, "read_data");

}



void	write_data (jj)	/* only record the variables, not the constants */

int	jj;

{

	int		alf, i;

	struct
	lf_data	*lf_ptr;

	if (global_monitor)
		monitor_leaf(1, "write_data");
/*
	for (i = 0; i < 12; i++)
	{
		int2str(eco[jj].prec[i], nint( 10.0 * P[i]), 5, 0);
		int2str(eco[jj].temp[i], nint(100.0 * T[i]), 5, 1);
	}
*/
	int2str(eco[jj].iprec, nint(   10.0 * leaf_iprec), 5, 0);	/* 0.1 mm     */
	int2str(eco[jj].itemp, nint(  100.0 * leaf_itemp), 5, 1);	/* 0.01 C     */
	int2str(eco[jj].degd,  nint(  100.0 * leaf_dd),    7, 0);	/* 0.1 C-days */
/*
	int2str(eco[jj].var,   nint(10000.0 * leaf_var),   5, 0);
	int2str(eco[jj].lfc,                  leaf_lfc,    3, 0);
*/
	lf_ptr = eco[jj].ptr;

	while (lf_ptr != 0)
	{
		alf = str2int(lf_ptr->lfn, 3) - 1;
/*
		int2str(lf_ptr->fc,     nint(10000.0 * leaf_fc        [alf]), 5, 0);
*/
		int2str(lf_ptr->lai,    nint(10000.0 * seasonal_lai[1][alf]), 6, 0);
		int2str(lf_ptr->laid,   nint(10000.0 * seasonal_lai[0][alf]), 6, 0);
		int2str(lf_ptr->season, (int)          lf_season      [alf],  1, 0);

		lf_ptr = lf_ptr->ptr;
	}

	if (global_monitor)
		monitor_leaf(0, "write_data");

}



void	read_leaf_file ()

{

	FILE		*fp;

	char		cs [77],
				vn  [7];

	int		i, j, lfc;

	struct
	lf_data	*lf_ptr;

	if (global_monitor)
		monitor_leaf(1, "read_leaf_file");

	initialize(1);

	printf("\n\n    Reading:  LEAF file ");

	fp = open_file("leaff/leaf_file");

	/* read the header */

	fscanf(fp, "%19c", &cs);
	cs[19] = '\0';
	fscanf(fp, "%s", &vn);

	if (strcmp(cs, "LEAF_file:  version") || strcmp(vn, version))
	{
		printf("\n\n    Incompatible LEAF file version. \n\n");
		printf("        cs = %s (%2d) \n", cs, strlen(cs));
		printf("        vn = %s (%2d) \n", vn, strlen(vn));
		printf("\n      Exiting ....... \n\n\n");
		exit(0);
	}

	fskips(fp, 3);
	fscanf(fp, "%d", &nlats);
	fskips(fp, 2);
	fscanf(fp, "%d", &nlons);
	fskips(fp, 4);
	fscanf(fp, "%d", &iflp);

	printf("with %d ice-free land points \n", iflp);

	fskipl(fp, 9);

	/* read the Earth map */

	read_earth(fp);

	/* read the records for each ice-free land point */

	fskipl(fp, 5);

	for (i = 0; i < iflp; i++)
	{
		fscanf(fp, "%*s%*1c%6c", &eco[i].iflpn);
		fscanf(fp, "%*s%*3c%5c", &eco[i].lat);
		fscanf(fp, "%*s%*1c%6c", &eco[i].lon);

		fscanf(fp, "%*s%*1c");

		for (j = 0; j < 12; j++)
			fscanf(fp, "%*1c%5c", &eco[i].prec[j]);

		fscanf(fp, "%*s%*1c");

		for (j = 0; j < 12; j++)
			fscanf(fp, "%*1c%5c", &eco[i].temp[j]);

		fscanf(fp, "%*s%*4c%5c", &eco[i].iprec);
		fscanf(fp, "%*s%*4c%5c", &eco[i].itemp);
		fscanf(fp, "%*s%*2c%7c", &eco[i].degd);
		fscanf(fp, "%*s%*2c%5c", &eco[i].var);
		fscanf(fp, "%*s%*3c%3c", &eco[i].lfc);

		lfc = str2int(eco[i].lfc, 3);

		if (lfc)
		{
			fskipl(fp, 4);

			lf_ptr = (struct lf_data *) 0;

			for (j = 0; j < lfc; j++)
			{
				new = (struct lf_data *) calloc (1, sizeof(struct lf_data));

				if (new == 0)
				{
					printf("Not enough memory for LF %d in Vector %d \n", j, i);
					exit(0);
				}

				if (lf_ptr != 0)
					lf_ptr->ptr = new;
				else
					eco[i].ptr = new;

				lf_ptr = new;

				fscanf(fp, "%3c",     &lf_ptr->lfn);
				fscanf(fp, "%*1c%5c", &lf_ptr->fc);
				fscanf(fp, "%*1c%6c", &lf_ptr->lai);
				fscanf(fp, "%*1c%6c", &lf_ptr->laid);
				fscanf(fp, "%*1c%1c%*1c", &lf_ptr->season);
			}

			fskipl(fp, 2);
		}
		else
			fskipl(fp, 3);
	}

	fclose(fp);

	if (global_monitor)
		monitor_leaf(0, "read_leaf_file");

}



void	write_leaf_file (int year)

{

	FILE		*fp;

	char		command [77],
				time     [7];

	int		i, j,
				lfc;

	struct
	lf_data	*lf_ptr;

	if (global_monitor)
		monitor_leaf(1, "write_leaf_file");

	printf("\n    Writing:  LEAF file with %d ice-free land points \n", iflp);

	fp = fopen("leaff/leaf_file", "w");

	/* write the header */

	fprintf(fp, "LEAF_file:  version %s\n\n", version);
	fprintf(fp, "Resolution:  nlats = %3d\n", n_lats);
	fprintf(fp, "             nlons = %3d\n\n", n_lons);
	fprintf(fp, "Ice-free land points = %d\n\n", iflp);
	fprintf(fp, "Precipitation: 0.1 mm\n");
	fprintf(fp, "Temperature:   0.01 C\n");

	for (i = 0; i < 75; i++)
		fprintf(fp, "_");

	fprintf(fp, "\n\nEarth map:\n");

	for (i = 0; i < 75; i++)
		fprintf(fp, "_");

	fprintf(fp, "\n\n");

	/* write the Earth map */

	write_earth(fp);

	/* write the records for each ice-free land point */

	for (i = 0; i < 75; i++)
		fprintf(fp, "_");

	fprintf(fp, "\n\nLEAF data for %d ice-free land points:\n", iflp);

	for (i = 0; i < 75; i++)
		fprintf(fp, "_");

	fprintf(fp, "\n\n");

	for (i = 0; i < iflp; i++)
	{
/*
		printf("      %1d  IFLP %6.6s\n", i+1, eco[i].iflpn);
*/
		fprintf(fp, "IFLP %6.6s\n\n", eco[i].iflpn);
		fprintf(fp, "Latitude:   %5.5s\n", eco[i].lat);
		fprintf(fp, "Longitude: %6.6s\n\n", eco[i].lon);

		fprintf(fp, "P: ");

		for (j = 0; j < 12; j++)
			fprintf(fp, " %5.5s", eco[i].prec[j]);

		fprintf(fp, "\nT: ");

		for (j = 0; j < 12; j++)
			fprintf(fp, " %5.5s", eco[i].temp[j]);

		fprintf(fp, "\n\n");
		fprintf(fp, "P(inertial):    %5.5s\n", eco[i].iprec);
		fprintf(fp, "T(inertial):    %5.5s\n", eco[i].itemp);
		fprintf(fp, "Degree-Days:  %7.7s\n\n", eco[i].degd);
		fprintf(fp, "VAR:  %5.5s\n\n", eco[i].var);
		fprintf(fp, "LFC = %3.3s life forms:\n", eco[i].lfc);

		lfc = str2int(eco[i].lfc, 3);

		if (lfc)
		{
			fprintf(fp, "\nLF   FC    LAI    LAId  S\n\n");

			lf_ptr = eco[i].ptr;

			while (lf_ptr != 0)
			{
				fprintf(fp, "%3.3s",    lf_ptr->lfn);
				fprintf(fp, " %5.5s",   lf_ptr->fc);
				fprintf(fp, " %6.6s",   lf_ptr->lai);
				fprintf(fp, " %6.6s",   lf_ptr->laid);
				fprintf(fp, " %1.1s\n", lf_ptr->season);

				lf_ptr = lf_ptr->ptr;
			}
		}

		for (j = 0; j < 75; j++)
			fprintf(fp, "_");

		fprintf(fp, "\n\n");
	}

	fprintf(fp, "\n");

	fclose(fp);

	strcpy(command, "cp leaff/leaf_file leaff/leaf_file_");

	if (date < 10)
		strcat(command, "00");
	else
		if (date < 100)
			strcat(command, "0");

	strcat(command, int_to_ascii(date, time));
	strcat(command, "_");

	strcat(command, int_to_ascii(year, time));

	system(command);

	if (global_monitor)
		monitor_leaf(0, "write_leaf_file");

}
 
 

void	initialize (mode)
 
int	mode;
 
{
 
	int	i, j, k;
 
	if (global_monitor)
		monitor_leaf(1, "initialize");
 
	if (mode)
	{
		for (i = 0; i < iflpm; i++)
		{
			for (j = 0; j < 6; j++)
				eco[i].iflpn[j] = '0';

			for (j = 0; j < 5; j++)
				eco[i].lat[j] = '0';

			for (j = 0; j < 6; j++)
				eco[i].lon[j] = '0';

			for (j = 0; j < 12; j++)
				for (k = 0; k < 5; k++)
				{
					eco[i].prec[j][k] = '0';
					eco[i].temp[j][k] = '0';
				}

			for (j = 0; j < 5; j++)
			{
				eco[i].iprec[j] = '0';
				eco[i].itemp[j] = '0';
				eco[i].var  [j] = '0';
			}

			for (j = 0; j < 7; j++)
				eco[i].degd [j] = '0';

			for (j = 0; j < 3; j++)
				eco[i].lfc[j] = '0';

			eco[i].ptr = (struct lf_data *) 0;
			active[i]  = (struct lf_data *) 0;
		}
	}
	else
		for (i = 0; i < iflpm; i++)
			active[i] = eco[i].ptr;

	if (global_monitor)
		monitor_leaf(0, "initialize");

}
 


void	read_earth (fp)

FILE	*fp;

{

	int	i, j,

			esp = 0,
			lps = 0,
			ops = 0,
			ips = 0,

			*sp [3];

	if (global_monitor)
		monitor_leaf(1, "read_earth");
/*
	printf("\n    Reading:  Earth surface data \n");
*/
	sp[0] = &ops;
	sp[1] = &lps;
	sp[2] = &ips;

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			fscanf(fp, "%1c", &earth[i][j]);
 
			++(*sp[(int)(earth[i][j]-'0')]);
 
			++esp;
		}
 
		fskipc(fp, 1);
	}

	printf("    Number of ocean points           = %7d \n",   ops);
	printf("    Number of land points            = %7d \n",   lps);
	printf("    Number of land ice points        = %7d \n\n", ips);
	printf("    Number of earth surface points   = %7d \n\n", esp);
    /* Add flush for SP */
    fflush(stdout);

/*
	fclose(fp);
*/
	if (global_monitor)
		monitor_leaf(0, "read_earth");

}



void	write_earth (fp)

FILE	*fp;

{

	int	i, j,

			esp = 0,
			lps = 0,
			ops = 0,
			ips = 0,

			*sp [3];

	if (global_monitor)
		monitor_leaf(1, "write_earth");
/*
	printf("\n    Writing:  Earth surface data \n");
*/
	sp[0] = &ops;
	sp[1] = &lps;
	sp[2] = &ips;

	for (i = 0; i < n_lats; i++)
	{
		for (j = 0; j < n_lons; j++)
		{
			fprintf(fp, "%1c", earth[i][j]);
 
			++(*sp[(int)(earth[i][j]-'0')]);
 
			++esp;
		}
 
		fprintf(fp, "\n");
	}

	printf("    Number of ocean points           = %7d \n",   ops);
	printf("    Number of land points            = %7d \n",   lps);
	printf("    Number of land ice points        = %7d \n\n", ips);
	printf("    Number of earth surface points   = %7d \n\n", esp);
/*
	fclose(fp);
*/
	if (global_monitor)
		monitor_leaf(0, "write_earth");

}


