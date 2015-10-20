/* Leaf Area Index file:  life-form seasonal LAI parameters */



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

}	leaf [110] =

	{

/*    Decid Lftyp     IP        DD       Day_0   Lat_0    D0/L0     Tau
    ------------------------------------------------------------------------
*/
/*   1. Equatorial Rainforest Broadleaf Trees         
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   2. Equatorial-Montane Rainforest Broadleaf Trees 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   3. Tropical Evergreen Microphyll-Broadleaf Trees 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   4. Tropical-Monsoon Evergreen Broadleaf Trees    
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   5. Tropical Evergreen Sclerophyll-Broadleaf Trees
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   6. Mediterranean Evergreen Broadleaf Trees       
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   7. Warm-Temperate-Mesic Evergreen Broadleaf Trees
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   8. Warm-Temperate Evergreen Broadleaf Trees      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*   9. Tropical-Monsoon Raingreen Broadleaf Trees    
*/
        1,    0,     3.00,       0.,       0,      0.,     0.0,      0.,
                     1.00,       0.,

/*  10. Tropical-Montane Raingreen Broadleaf Trees    
*/
        1,    0,     2.50,       0.,       0,      0.,     0.0,      0.,
                     1.00,       0.,

/*  11. Tropical-Xeric Raingreen Broadleaf Trees      
*/
        1,    0,     2.50,       0.,       0,      0.,     0.0,      0.,
                     1.00,       0.,

/*  12. Equatorial-Xeric Raingreen Broadleaf Trees    
*/
        1,    0,     2.50,       0.,       0,      0.,     0.0,      0.,
                     1.00,       0.,

/*  13. Western-Temperate Summergreen Broadleaf Trees 
*/
        2,    0,     0.00,     800.,     273,     40.,     2.5,     25.,
                     0.00,     350.,

/*  14. Eastern-Temperate Summergreen Broadleaf Trees 
*/
        2,    0,     0.00,     800.,     273,     40.,     2.5,     25.,
                     0.00,     350.,

/*  15. Eastern-Temp-Xeric Summergreen Broadleaf Trees
*/
        2,    0,     0.00,     800.,     273,     40.,     2.5,     25.,
                     0.00,     350.,

/*  16. Cool-Temp/Boreal Summergreen Broadleaf Trees  
*/
        2,    0,     0.00,     450.,     245,     50.,     0.8,     20.,
                     0.00,     100.,

/*  17. Tropical Evergreen Linearleaf Trees           
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  18. Tropical-Xeric Evergreen Needleleaf Trees     
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  19. Mediterranean Evergreen Needleleaf Trees      
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  20. Heliophilic Evergreen Needleleaf Trees        
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  21. Warm-Temp-Montane Evergreen Needleleaf Trees  
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  22. Temperate-Mesic Evergreen Needleleaf Trees    
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  23. Western-Temperate Evergreen Needleleaf Trees  
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  24. Eastern-Temperate Evergreen Needleleaf Trees  
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  25. Cool-Temp/Boreal Evergreen Needleleaf Trees   
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  26. Cool-Temp/Boreal Summergreen Needleleaf Trees 
*/
        2,    1,     0.00,     400.,     245,     50.,     0.8,     20.,
                     0.00,      75.,

/*  27. Hydrophilic Summergreen Needleleaf Trees      
*/
        2,    1,     0.00,     650.,     270,     45.,     2.3,     20.,
                     0.00,     250.,

/*  28. Tropical Evergreen Broadleaf Small-Trees      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  29. Warm-Temperate Evergreen Broadleaf Small-Trees
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  30. Cool-Maritime Evergreen Broadleaf Small-Trees 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  31. Tropical Raingreen Broadleaf Small-Trees      
*/
        1,    0,     2.00,       0.,       0,      0.,     0.0,      0.,
                     1.00,       0.,

/*  32. Temperate Summergreen Broadleaf Small-Trees   
*/
        2,    0,     0.00,     700.,     273,     40.,     2.5,     25.,
                     0.00,     350.,

/*  33. Temperate Evergreen Needleleaf Small-Trees    
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  34. Tropical Evergreen Broadleaf Dwarf-Trees      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  35. Tropical-Montane Evergreen Dwarf-Trees        
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  36. Temp-Montane/Boreal Needleleaf Dwarf-Trees    
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  37. Palmiform Rosette-Trees                       
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  38. Palmiform Small-Rosette-Trees                 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  39. Fern Small-Rosette-Trees                      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  40. Tropical-Montane Small-Rosette-Trees          
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  41. Xeric Small-Rosette-Trees                     
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  42. Raingreen Broadleaf Arborescent-Shrubs        
*/
        1,    0,     2.00,       0.,       0,      0.,     0.0,      0.,
                     0.70,       0.,

/*  43. Evergreen Broadleaf Arborescent-Shrubs        
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  44. Summergreen Broadleaf Arborescent-Shrubs      
*/
        2,    0,     0.00,     600.,     270,     45.,     2.3,     20.,
                     0.00,     100.,

/*  45. Xeric Leafless Arborescent-Shrubs             
*/
        3,    0,     1.00,     250.,     270,     45.,     2.3,     20.,
                     0.50,     100.,

/*  46. Tropical Evergreen Broadleaf Shrubs           
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  47. Hot-Desert Evergreen Broadleaf Shrubs         
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  48. Desert Evergreen Leaf-Succulent Shrubs        
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  49. Mediterranean Evergreen Broadleaf Shrubs      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  50. Mesic-Montane Evergreen Broadleaf Shrubs      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  51. Warm-Temperate Evergreen Broadleaf Shrubs     
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  52. Cool-Temperate-Xeric Evergreen Shrubs         
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  53. Temperate Summergreen Broadleaf Shrubs        
*/
        2,    0,     0.00,     400.,     273,     40.,     2.5,     25.,
                     0.00,     250.,

/*  54. Cool-Temp/Boreal Summergreen Broadleaf Shrubs 
*/
        2,    0,     0.00,     300.,     245,     50.,     0.8,     20.,
                     0.00,      70.,

/*  55. Temperate-Xeric Summergreen Broadleaf Shrubs  
*/
        2,    0,     0.00,     700.,     273,     40.,     2.5,     25.,
                     0.00,     350.,

/*  56. Mediterranean Summergreen Broadleaf Shrubs    
*/
        2,    0,     0.00,     300.,     273,     40.,     2.5,     25.,
                     0.00,     100.,

/*  57. Mediterranean Evergreen Needleleaf Shrubs     
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  58. Evergreen Needleleaf Shrubs                   
*/
        0,    1,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  59. Xeric Dwarf-Shrubs                            
*/
        3,    0,     0.70,     500.,     270,     45.,     2.3,     20.,
                     0.30,     100.,

/*  60. Mediterranean Dwarf-Shrubs                    
*/
        1,    0,     2.00,       0.,       0,      0.,     0.0,      0.,
                     0.50,       0.,

/*  61. Cool-Maritime Evergreen Dwarf-Shrubs          
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  62. Cool-Temp/Boreal Evergreen Dwarf-Shrubs       
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  63. Tundra Evergreen Dwarf-Shrubs                 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  64. Tundra Summergreen Dwarf-Shrubs               
*/
        2,    0,     0.00,      75.,     234,     55.,     0.5,     15.,
                     0.00,      25.,

/*  65. Cool-Maritime/Tropical-Alpine Cushion-Shrubs  
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  66. Xeric Cushion-Shrubs                          
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  67. Mesic Palmiform Rosette-Shrubs                
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  68. Xeric Rosette-Shrubs                          
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  69. Arborescent Stem-Succulents                   
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  70. Branched Stem-Succulents                      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  71. Unbranched Stem-Succulents                    
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  72. Arborescent Grasses                           
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  73. Tall-Cane Grasses                             
*/
        3,    0,     3.50,     700.,     270,     45.,     2.3,     20.,
                     0.50,     125.,

/*  74. Sclerophyllous Grasses                        
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  75. Tall Grasses                                  
*/
        3,    0,     3.00,     700.,     270,     45.,     2.3,     20.,
                     0.50,     125.,

/*  76. Short-Sward Graminoids                        
*/
        3,    0,     2.00,     475.,     270,     45.,     2.3,     20.,
                     0.50,     125.,

/*  77. Short-Bunch Graminoids                        
*/
        3,    0,     2.00,     525.,     270,     45.,     2.3,     20.,
                     0.50,     125.,

/*  78. Desert Short-Bunch Graminoids                 
*/
        1,    0,     1.00,       0.,       0,      0.,     0.0,      0.,
                     0.50,       0.,

/*  79. Temperate Tall Tussock-Graminoids             
*/
        2,    0,     0.00,     700.,     273,     40.,     2.5,     25.,
                     0.00,     125.,

/*  80. Temperate Short Tussock-Graminoids            
*/
        2,    0,     0.00,     525.,     273,     40.,     2.5,     25.,
                     0.00,     125.,

/*  81. Tropical-Montane Short Tussock-Graminoids     
*/
        3,    0,     2.00,     350.,     270,     45.,     2.3,     20.,
                     0.50,     125.,

/*  82. Cool-Maritime Short Tussock-Graminoids        
*/
        2,    0,     0.00,     300.,     270,     45.,     2.3,     20.,
                     0.00,     125.,

/*  83. Tundra Short Tussock-Graminoids               
*/
        2,    0,     0.00,      25.,     234,     55.,     0.5,     15.,
                     0.00,       5.,

/*  84. Tropical Evergreen Forbs                      
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  85. Temperate Evergreen Forbs                     
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  86. Succulent Forbs                               
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  87. Raingreen Forbs                               
*/
        1,    0,     1.50,       0.,       0,      0.,     0.0,      0.,
                     0.80,       0.,

/*  88. Summergreen Forbs                             
*/
        2,    0,     0.00,     400.,     270,     45.,     2.3,     20.,
                     0.00,     125.,

/*  89. Desert Ephemeral Herbs                        
*/
        1,    0,     0.50,       0.,       0,      0.,     0.0,      0.,
                     0.10,       0.,

/*  90. Tropical-Alpine Raingreen Herbs               
*/
        1,    0,     2.00,       0.,       0,      0.,     0.0,      0.,
                     0.50,       0.,

/*  91. Tundra Summergreen Herbs                      
*/
        2,    0,     0.00,      25.,     234,     55.,     0.5,     15.,
                     0.00,       5.,

/*  92. Desert Cushion-Herbs                          
*/
        1,    0,     0.70,       0.,       0,      0.,     0.0,      0.,
                     0.20,       0.,

/*  93. Tropical-Montane Cushion-Herbs                
*/
        3,    0,     2.00,     300.,     270,     45.,     2.3,     20.,
                     0.50,     125.,

/*  94. Tundra Summergreen Cushion-Herbs              
*/
        2,    0,     0.00,      30.,     234,     55.,     0.5,     15.,
                     0.00,       5.,

/*  95. Tropical Liana Evergreen Vines                
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  96. Evergreen Vines                               
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/*  97. Raingreen Vines                               
*/
        1,    0,     2.00,       0.,       0,      0.,     0.0,      0.,
                     1.00,       0.,

/*  98. Summergreen Vines                             
*/
        2,    0,     0.00,     300.,     270,     45.,     2.3,     20.,
                     0.00,     150.,

/*  99. Tropical Broadleaf Epiphytes                  
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 100. Broadleaf Epiphytes                           
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 101. Narrowleaf Epiphytes                          
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 102. Evergreen Ferns                               
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 103. Summergreen Ferns                             
*/
        2,    0,     0.00,     500.,     270,     45.,     2.3,     20.,
                     0.00,     300.,

/* 104. Mesic Peat-Forming Bryophytes                 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 105. Mesic Bryophytes                              
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 106. Epiphytic Bryophytes                          
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 107. Mesic Lichens                                 
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 108. Boreal Lichens                                
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 109. Desert Cryptogams                             
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.,

/* 110. Polar Cryptogams                              
*/
        0,    0,     0.00,       0.,       0,      0.,     0.0,      0.,
                     0.00,       0.

	};

