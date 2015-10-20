/* LAIs:  the Leaf Area Index (LAI) Model:  plant community level */



# include	"gv.h"				/* Global Variables file */



enum canopy_levels	{ tree = 1, vine, epiphyte, shrub, grass };


extern
enum canopy_levels	lf_pcl [num_lf];


extern
double	lf_lai [num_lf],

			seasonal_lai [2] [num_lf];



double	alpha,	/* (V and E) tree lai displacement coefficient */

			canopy_lai,
			canopy_laid;



void	calc_lais ()

{

	double	lai_fc,
				laid_fc,
				tree_lai_max,
				vine_epiphyte_lai;

	int		i;

	if (global_monitor)
		monitor_leaf(1, "calc_lais");

	canopy_lai        = 0.0;
	canopy_laid       = 0.0;
	tree_lai_max      = 0.0;
	vine_epiphyte_lai = 0.0;

	/* sum daily LAIs for each canopy level */

	for (i = 0; i < num_lf; i++)
		if (leaf_lf[i])
		{
			lai_fc  = leaf_fc[i] * seasonal_lai[1][i];
			laid_fc = leaf_fc[i] * seasonal_lai[0][i];

			switch (lf_pcl[i])
			{
				case tree:

					tree_lai_max += leaf_fc[i] * lf_lai[i];

				case shrub:
				case grass:

					canopy_lai  += lai_fc;
					canopy_laid += laid_fc;

					break;

				case vine:
				case epiphyte:

					vine_epiphyte_lai += lai_fc;

					break;
			}
		}

	canopy_lai  *= leaf_var;
	canopy_laid *= leaf_var;

	/* Adjust the tree LAIs to include vines and epiphytes */

	if (tree_lai_max > 0.0)
		alpha = 1.0 - (vine_epiphyte_lai / tree_lai_max);
	else
		alpha = 1.0;

	if (global_monitor)
		monitor_leaf(0, "calc_lais");

}

