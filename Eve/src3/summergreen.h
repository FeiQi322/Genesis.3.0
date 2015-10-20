/* Leaf Area Index file:  summergreen LFs seasonal LAI parameters */



struct summer_lai_parameters
{
	int		day_0;

	double	lat_0,
				d0_l0,
				tau;

}	leaff [4] =

	{

/*       Day_0    Lat_0    D0/L0     Tau
    ------------------------------------------------------------------------
*/
/*   1. Universal Temperate         
*/
          270,     45.,     2.3,     20.,

/*   2. Temperate 
*/
          273,     40.,     2.5,     25.,

/*   3. Boreal 
*/
          245,     50.,     0.8,     20.,

/*   4. Tundra    
*/
          234,     55.,     0.5,     15.

	};

