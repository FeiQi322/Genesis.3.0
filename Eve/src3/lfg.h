/* Life-Form Geometry */



struct life_form_geometry
{
	double	z1,	/* height of canopy, top      (m) */
				z0,	/* height of canopy, bottom   (m) */
				ll,	/* leaf length               (cm) */
				lw,	/* leaf width                (cm) */
				sd,	/* stem dimension            (cm) */
				rd,	/* root depth (conical form) (cm) */

				lo;	/* leaf orientation:                          */
						/*  1 = horizontal, 0 = random, -1 = vertical */

}  lfg [110] =

	{

/*        Z1       Z0      LL      LW      SD      RD       LO
-----------------------------------------------------------------------------
*/
/*   1. Equatorial Rainforest Broadleaf Trees         
*/
         40.00,  28.00,   17.0,   7.00,   40.3,   100.0,   0.10,

/*   2. Equatorial-Montane Rainforest Broadleaf Trees 
*/
         30.00,  21.00,   10.0,   4.00,   24.1,   100.0,   0.10,

/*   3. Tropical Evergreen Microphyll-Broadleaf Trees 
*/
         30.00,  21.00,    2.5,   1.00,   24.1,   100.0,   0.10,

/*   4. Tropical-Monsoon Evergreen Broadleaf Trees    
*/
         30.00,  15.00,    7.0,   3.00,   24.1,   100.0,   0.10,

/*   5. Tropical Evergreen Sclerophyll-Broadleaf Trees
*/
         35.00,  20.00,   15.0,   4.00,   27.4,   425.0,   0.10,

/*   6. Mediterranean Evergreen Broadleaf Trees       
*/
         20.00,   5.00,    6.0,   3.00,   15.8,   397.5,   0.10,

/*   7. Warm-Temperate-Mesic Evergreen Broadleaf Trees
*/
         20.00,   5.00,   12.0,   7.00,   19.0,   397.5,   0.10,

/*   8. Warm-Temperate Evergreen Broadleaf Trees      
*/
         15.00,   3.00,    9.0,   4.50,   17.2,   358.7,   0.10,

/*   9. Tropical-Monsoon Raingreen Broadleaf Trees    
*/
         30.00,  21.00,    7.0,   3.00,   24.1,   250.0,   0.25,

/*  10. Tropical-Montane Raingreen Broadleaf Trees    
*/
         25.00,  15.00,    6.0,   3.00,   20.0,   250.0,   0.25,

/*  11. Tropical-Xeric Raingreen Broadleaf Trees      
*/
         15.00,   5.00,    2.5,   1.00,   11.7,   358.7,   0.25,

/*  12. Equatorial-Xeric Raingreen Broadleaf Trees    
*/
         15.00,   5.00,    2.5,   1.00,   11.7,   358.7,   0.25,

/*  13. Western-Temperate Summergreen Broadleaf Trees 
*/
         25.00,  10.00,   15.0,  10.00,   20.0,   418.5,   0.25,

/*  14. Eastern-Temperate Summergreen Broadleaf Trees 
*/
         25.00,  10.00,   15.0,   9.00,   23.9,   418.5,   0.25,

/*  15. Eastern-Temp-Xeric Summergreen Broadleaf Trees
*/
         25.00,  10.00,   15.0,  10.00,   23.9,   418.5,   0.25,

/*  16. Cool-Temp/Boreal Summergreen Broadleaf Trees  
*/
         20.00,   8.00,    5.0,   3.00,   14.3,   397.5,   0.25,

/*  17. Tropical Evergreen Linearleaf Trees           
*/
         40.00,  15.00,   15.0,   1.00,   40.3,   425.0,   0.10,

/*  18. Tropical-Xeric Evergreen Needleleaf Trees     
*/
         25.00,  10.00,   15.0,   0.10,   20.0,   418.5,   0.01,

/*  19. Mediterranean Evergreen Needleleaf Trees      
*/
         15.00,   5.00,   17.0,   0.10,   14.1,   358.7,   0.01,

/*  20. Heliophilic Evergreen Needleleaf Trees        
*/
         30.00,  15.00,   25.0,   0.10,   24.9,   425.0,   0.01,

/*  21. Warm-Temp-Montane Evergreen Needleleaf Trees  
*/
         30.00,   5.00,   15.0,   0.10,   26.5,   425.0,   0.01,

/*  22. Temperate-Mesic Evergreen Needleleaf Trees    
*/
         60.00,  20.00,    2.0,   0.10,   69.8,   425.0,   0.01,

/*  23. Western-Temperate Evergreen Needleleaf Trees  
*/
         30.00,   5.00,   15.0,   0.10,   28.1,   425.0,   0.01,

/*  24. Eastern-Temperate Evergreen Needleleaf Trees  
*/
         25.00,   5.00,   10.0,   0.10,   25.5,   418.5,   0.01,

/*  25. Cool-Temp/Boreal Evergreen Needleleaf Trees   
*/
         18.00,   5.00,    2.0,   0.10,   11.8,   384.3,   0.01,

/*  26. Cool-Temp/Boreal Summergreen Needleleaf Trees 
*/
         17.00,   5.00,    2.5,   0.10,   11.0,   376.6,   0.01,

/*  27. Hydrophilic Summergreen Needleleaf Trees      
*/
         35.00,  20.00,    1.6,   0.10,   36.2,   425.0,   0.01,

/*  28. Tropical Evergreen Broadleaf Small-Trees      
*/
          9.00,   2.00,   10.0,   4.00,    6.6,   100.0,   0.10,

/*  29. Warm-Temperate Evergreen Broadleaf Small-Trees
*/
          7.00,   1.50,    9.0,   4.50,    5.0,   244.8,   0.10,

/*  30. Cool-Maritime Evergreen Broadleaf Small-Trees 
*/
          7.00,   1.50,    9.0,   4.50,    5.0,   244.8,   0.10,

/*  31. Tropical Raingreen Broadleaf Small-Trees      
*/
          7.00,   1.50,    5.0,   2.50,    5.0,   244.8,   0.25,

/*  32. Temperate Summergreen Broadleaf Small-Trees   
*/
          7.00,   1.50,   15.0,  10.00,    5.0,   244.8,   0.25,

/*  33. Temperate Evergreen Needleleaf Small-Trees    
*/
          7.00,   1.00,    3.0,   0.10,    5.0,   244.8,   0.01,

/*  34. Tropical Evergreen Broadleaf Dwarf-Trees      
*/
          3.00,   0.50,    2.5,   1.00,    2.5,   147.7,   0.10,

/*  35. Tropical-Montane Evergreen Dwarf-Trees        
*/
          5.00,   1.00,    5.0,   2.00,    3.8,   100.0,   0.10,

/*  36. Temp-Montane/Boreal Needleleaf Dwarf-Trees    
*/
          3.00,   0.25,    2.0,   0.10,    8.2,   147.7,   0.01,

/*  37. Palmiform Rosette-Trees                       
*/
         20.00,  15.00,  100.0,  30.00,   15.8,   397.5,   0.10,

/*  38. Palmiform Small-Rosette-Trees                 
*/
          7.00,   2.00,  100.0,  30.00,   10.3,   244.8,   0.10,

/*  39. Fern Small-Rosette-Trees                      
*/
          7.00,   3.50,   30.0,   5.00,    7.3,   244.8,   0.10,

/*  40. Tropical-Montane Small-Rosette-Trees          
*/
          4.00,   3.00,   50.0,  15.00,    3.1,   176.0,   0.10,

/*  41. Xeric Small-Rosette-Trees                     
*/
          3.00,   1.50,   50.0,   3.00,    2.5,   147.7,   0.10,

/*  42. Raingreen Broadleaf Arborescent-Shrubs        
*/
          7.00,   1.00,   15.0,   3.00,    5.0,   244.8,   0.25,

/*  43. Evergreen Broadleaf Arborescent-Shrubs        
*/
         10.00,   1.00,   15.0,   3.00,    7.6,   296.9,   0.10,

/*  44. Summergreen Broadleaf Arborescent-Shrubs      
*/
          6.00,   1.00,   15.0,   3.00,    4.5,   224.1,   0.25,

/*  45. Xeric Leafless Arborescent-Shrubs             
*/
          6.00,   1.00,    1.0,   0.20,    4.5,   224.1,   0.10,

/*  46. Tropical Evergreen Broadleaf Shrubs           
*/
          3.00,   0.50,   10.0,   4.00,    2.4,   100.0,   0.10,

/*  47. Hot-Desert Evergreen Broadleaf Shrubs         
*/
          2.00,   0.00,    0.7,   0.30,    1.6,   115.2,   0.10,

/*  48. Desert Evergreen Leaf-Succulent Shrubs        
*/
          1.00,   0.00,    4.0,   3.00,    0.8,    75.9,   0.10,

/*  49. Mediterranean Evergreen Broadleaf Shrubs      
*/
          3.00,   0.25,    5.0,   2.50,    2.4,   300.0,   0.10,

/*  50. Mesic-Montane Evergreen Broadleaf Shrubs      
*/
          3.00,   1.50,   15.0,   4.00,    2.4,   147.7,   0.10,

/*  51. Warm-Temperate Evergreen Broadleaf Shrubs     
*/
          3.00,   1.50,    6.0,   2.50,    2.4,   147.7,   0.10,

/*  52. Cool-Temperate-Xeric Evergreen Shrubs         
*/
          0.75,   0.00,    2.5,   0.75,    0.6,    64.1,   0.10,

/*  53. Temperate Summergreen Broadleaf Shrubs        
*/
          3.00,   1.00,    8.0,   4.00,    2.4,   147.7,   0.25,

/*  54. Cool-Temp/Boreal Summergreen Broadleaf Shrubs 
*/
          2.50,   0.80,    6.0,   3.00,    2.0,   132.0,   0.25,

/*  55. Temperate-Xeric Summergreen Broadleaf Shrubs  
*/
          3.00,   0.25,    7.0,   4.00,    2.4,   147.7,   0.25,

/*  56. Mediterranean Summergreen Broadleaf Shrubs    
*/
          3.00,   0.25,    3.5,   2.50,    2.4,   300.0,   0.25,

/*  57. Mediterranean Evergreen Needleleaf Shrubs     
*/
          2.00,   0.25,    0.7,   0.10,    1.6,   300.0,   0.01,

/*  58. Evergreen Needleleaf Shrubs                   
*/
          2.00,   0.00,    1.0,   0.15,    1.6,   115.2,   0.01,

/*  59. Xeric Dwarf-Shrubs                            
*/
          0.60,   0.00,    1.0,   0.20,    0.5,    56.2,   0.10,

/*  60. Mediterranean Dwarf-Shrubs                    
*/
          0.70,   0.00,    2.5,   1.50,    0.6,    70.0,   0.10,

/*  61. Cool-Maritime Evergreen Dwarf-Shrubs          
*/
          0.50,   0.00,    2.0,   0.20,    0.4,    50.0,   0.10,

/*  62. Cool-Temp/Boreal Evergreen Dwarf-Shrubs       
*/
          0.50,   0.00,    2.0,   0.20,    0.4,    50.0,   0.10,

/*  63. Tundra Evergreen Dwarf-Shrubs                 
*/
          0.15,   0.00,    2.0,   0.50,    0.1,    15.0,   0.10,

/*  64. Tundra Summergreen Dwarf-Shrubs               
*/
          0.50,   0.00,    2.0,   1.00,    0.4,    50.0,   0.25,

/*  65. Cool-Maritime/Tropical-Alpine Cushion-Shrubs  
*/
          0.30,   0.00,    2.0,   1.00,    0.2,    30.0,   0.25,

/*  66. Xeric Cushion-Shrubs                          
*/
          0.20,   0.00,    1.0,   0.50,    0.2,    20.0,   0.25,

/*  67. Mesic Palmiform Rosette-Shrubs                
*/
          2.00,   0.50,   50.0,  50.00,    1.6,   115.2,   0.10,

/*  68. Xeric Rosette-Shrubs                          
*/
          0.45,   0.00,   45.0,   5.00,    0.4,    45.0,   0.10,

/*  69. Arborescent Stem-Succulents                   
*/
          7.00,   3.00,  100.0,  20.00,   14.7,   244.8,  -0.50,

/*  70. Branched Stem-Succulents                      
*/
          1.00,   0.15,   10.0,   5.00,    0.8,    75.9,   0.00,

/*  71. Unbranched Stem-Succulents                    
*/
          0.30,   0.00,   30.0,  15.00,    0.2,    30.0,  -0.50,

/*  72. Arborescent Grasses                           
*/
          7.00,   3.00,   15.0,   3.00,    5.0,   100.0,   0.25,

/*  73. Tall-Cane Grasses                             
*/
          3.00,   0.00,   60.0,   3.00,    0.0,   100.0,  -0.25,

/*  74. Sclerophyllous Grasses                        
*/
          0.50,   0.00,   50.0,   1.00,    0.0,    50.0,  -0.10,

/*  75. Tall Grasses                                  
*/
          1.50,   0.00,   50.0,   1.00,    0.0,   200.0,  -0.25,

/*  76. Short-Sward Graminoids                        
*/
          0.40,   0.00,   50.0,   1.00,    0.0,   100.0,  -0.30,

/*  77. Short-Bunch Graminoids                        
*/
          0.40,   0.00,   50.0,   1.00,    0.0,   100.0,  -0.10,

/*  78. Desert Short-Bunch Graminoids                 
*/
          0.20,   0.00,   20.0,   1.00,    0.0,    20.0,  -0.20,

/*  79. Temperate Tall Tussock-Graminoids             
*/
          2.00,   0.00,  100.0,   3.00,    0.0,   115.2,  -0.10,

/*  80. Temperate Short Tussock-Graminoids            
*/
          0.50,   0.00,   50.0,   1.00,    0.0,    50.0,  -0.20,

/*  81. Tropical-Montane Short Tussock-Graminoids     
*/
          0.30,   0.00,   30.0,   1.00,    0.0,    30.0,  -0.20,

/*  82. Cool-Maritime Short Tussock-Graminoids        
*/
          0.30,   0.00,   30.0,   1.00,    0.0,    30.0,  -0.20,

/*  83. Tundra Short Tussock-Graminoids               
*/
          0.30,   0.00,   30.0,   1.00,    0.0,    30.0,  -0.20,

/*  84. Tropical Evergreen Forbs                      
*/
          1.00,   0.00,   25.0,   5.00,    0.0,    75.9,   0.10,

/*  85. Temperate Evergreen Forbs                     
*/
          0.10,   0.00,    3.0,   1.00,    0.0,    10.0,   0.10,

/*  86. Succulent Forbs                               
*/
          0.10,   0.00,    1.0,   0.50,    0.0,    10.0,   0.05,

/*  87. Raingreen Forbs                               
*/
          0.50,   0.00,    7.0,   3.00,    0.0,   100.0,   0.10,

/*  88. Summergreen Forbs                             
*/
          0.40,   0.00,    7.0,   3.00,    0.0,    40.0,   0.10,

/*  89. Desert Ephemeral Herbs                        
*/
          0.30,   0.00,   10.0,   1.00,    0.0,    30.0,   0.00,

/*  90. Tropical-Alpine Raingreen Herbs               
*/
          0.20,   0.00,   20.0,   1.00,    0.0,    20.0,   0.10,

/*  91. Tundra Summergreen Herbs                      
*/
          0.20,   0.00,   20.0,   1.00,    0.0,    20.0,   0.10,

/*  92. Desert Cushion-Herbs                          
*/
          0.15,   0.00,    1.0,   0.50,    0.0,    15.0,   0.25,

/*  93. Tropical-Montane Cushion-Herbs                
*/
          0.15,   0.00,    1.0,   0.50,    0.0,    15.0,   0.25,

/*  94. Tundra Summergreen Cushion-Herbs              
*/
          0.15,   0.00,    1.0,   0.50,    0.0,    15.0,   0.25,

/*  95. Tropical Liana Evergreen Vines                
*/
          0.00,   0.00,   15.0,   7.00,    1.0,   100.0,   0.50,

/*  96. Evergreen Vines                               
*/
          0.00,   0.00,    7.0,   3.00,    1.0,   100.0,   0.50,

/*  97. Raingreen Vines                               
*/
          0.00,   0.00,    5.0,   2.00,    1.0,   300.0,   0.50,

/*  98. Summergreen Vines                             
*/
          0.00,   0.00,   15.0,  12.00,    1.0,   100.0,   0.50,

/*  99. Tropical Broadleaf Epiphytes                  
*/
          0.00,   0.00,   50.0,   7.00,    0.0,     0.0,   0.10,

/* 100. Broadleaf Epiphytes                           
*/
          0.00,   0.00,    3.0,   2.00,    0.0,     0.0,   0.10,

/* 101. Narrowleaf Epiphytes                          
*/
          0.00,   0.00,    6.0,   0.20,    0.0,     0.0,   0.10,

/* 102. Evergreen Ferns                               
*/
          1.00,   0.00,   50.0,  10.00,    0.0,    75.9,   0.10,

/* 103. Summergreen Ferns                             
*/
          0.50,   0.00,   25.0,   5.00,    0.0,    50.0,   0.10,

/* 104. Mesic Peat-Forming Bryophytes                 
*/
          0.15,   0.00,    0.3,   0.20,    0.0,    15.0,   0.10,

/* 105. Mesic Bryophytes                              
*/
          0.03,   0.00,    0.2,   0.10,    0.0,     3.0,   0.10,

/* 106. Epiphytic Bryophytes                          
*/
          0.00,   0.00,    0.2,   0.10,    0.0,     0.0,   0.10,

/* 107. Mesic Lichens                                 
*/
          0.01,   0.00,    0.5,   0.50,    0.0,     0.5,   0.70,

/* 108. Boreal Lichens                                
*/
          0.04,   0.00,    0.5,   0.20,    0.0,     0.5,  -0.50,

/* 109. Desert Cryptogams                             
*/
          0.002,  0.00,    0.1,   0.10,    0.0,     0.1,   0.50,

/* 110. Polar Cryptogams                              
*/
          0.01,   0.00,    0.2,   0.10,    0.0,     1.0,  -0.10

	};

