/* Stomatal Conductance Coefficients */



struct sc_coefficients
{
	double	scmax,			/* mm / s   */
				scmin,			/* mm / s   */
				scpar,			/*  W / m^2 */
				scvpd;			/* mb       */

}	lf_scc [110] =

	{

/*               SCMAX          SCMIN           SCPAR            SCVPD       */
/*         ----------------------------------------------------------------- */

/*   1. Equatorial Rainforest Broadleaf Trees         
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   2. Equatorial-Montane Rainforest Broadleaf Trees 
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   3. Tropical Evergreen Microphyll-Broadleaf Trees 
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   4. Tropical-Monsoon Evergreen Broadleaf Trees    
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   5. Tropical Evergreen Sclerophyll-Broadleaf Trees
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   6. Mediterranean Evergreen Broadleaf Trees       
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   7. Warm-Temperate-Mesic Evergreen Broadleaf Trees
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   8. Warm-Temperate Evergreen Broadleaf Trees      
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*   9. Tropical-Monsoon Raingreen Broadleaf Trees    
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  10. Tropical-Montane Raingreen Broadleaf Trees    
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  11. Tropical-Xeric Raingreen Broadleaf Trees      
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  12. Equatorial-Xeric Raingreen Broadleaf Trees    
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  13. Western-Temperate Summergreen Broadleaf Trees 
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  14. Eastern-Temperate Summergreen Broadleaf Trees 
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  15. Eastern-Temp-Xeric Summergreen Broadleaf Trees
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  16. Cool-Temp/Boreal Summergreen Broadleaf Trees  
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  17. Tropical Evergreen Linearleaf Trees           
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  18. Tropical-Xeric Evergreen Needleleaf Trees     
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  19. Mediterranean Evergreen Needleleaf Trees      
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  20. Heliophilic Evergreen Needleleaf Trees        
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  21. Warm-Temp-Montane Evergreen Needleleaf Trees  
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  22. Temperate-Mesic Evergreen Needleleaf Trees    
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  23. Western-Temperate Evergreen Needleleaf Trees  
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  24. Eastern-Temperate Evergreen Needleleaf Trees  
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  25. Cool-Temp/Boreal Evergreen Needleleaf Trees   
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  26. Cool-Temp/Boreal Summergreen Needleleaf Trees 
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  27. Hydrophilic Summergreen Needleleaf Trees      
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  28. Tropical Evergreen Broadleaf Small-Trees      
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  29. Warm-Temperate Evergreen Broadleaf Small-Trees
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  30. Cool-Maritime Evergreen Broadleaf Small-Trees 
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  31. Tropical Raingreen Broadleaf Small-Trees      
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  32. Temperate Summergreen Broadleaf Small-Trees   
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  33. Temperate Evergreen Needleleaf Small-Trees    
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  34. Tropical Evergreen Broadleaf Dwarf-Trees      
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  35. Tropical-Montane Evergreen Dwarf-Trees        
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  36. Temp-Montane/Boreal Needleleaf Dwarf-Trees    
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  37. Palmiform Rosette-Trees                       
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  38. Palmiform Small-Rosette-Trees                 
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  39. Fern Small-Rosette-Trees                      
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  40. Tropical-Montane Small-Rosette-Trees          
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  41. Xeric Small-Rosette-Trees                     
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  42. Raingreen Broadleaf Arborescent-Shrubs        
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  43. Evergreen Broadleaf Arborescent-Shrubs        
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  44. Summergreen Broadleaf Arborescent-Shrubs      
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  45. Xeric Leafless Arborescent-Shrubs             
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  46. Tropical Evergreen Broadleaf Shrubs           
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  47. Hot-Desert Evergreen Broadleaf Shrubs         
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  48. Desert Evergreen Leaf-Succulent Shrubs        
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  49. Mediterranean Evergreen Broadleaf Shrubs      
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  50. Mesic-Montane Evergreen Broadleaf Shrubs      
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  51. Warm-Temperate Evergreen Broadleaf Shrubs     
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  52. Cool-Temperate-Xeric Evergreen Shrubs         
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  53. Temperate Summergreen Broadleaf Shrubs        
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  54. Cool-Temp/Boreal Summergreen Broadleaf Shrubs 
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  55. Temperate-Xeric Summergreen Broadleaf Shrubs  
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  56. Mediterranean Summergreen Broadleaf Shrubs    
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  57. Mediterranean Evergreen Needleleaf Shrubs     
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  58. Evergreen Needleleaf Shrubs                   
*/
              4.291845494,   0.991398944,    16.017596567,   32.258064516,

/*  59. Xeric Dwarf-Shrubs                            
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  60. Mediterranean Dwarf-Shrubs                    
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  61. Cool-Maritime Evergreen Dwarf-Shrubs          
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  62. Cool-Temp/Boreal Evergreen Dwarf-Shrubs       
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  63. Tundra Evergreen Dwarf-Shrubs                 
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  64. Tundra Summergreen Dwarf-Shrubs               
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  65. Cool-Maritime/Tropical-Alpine Cushion-Shrubs  
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  66. Xeric Cushion-Shrubs                          
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  67. Mesic Palmiform Rosette-Shrubs                
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  68. Xeric Rosette-Shrubs                          
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  69. Arborescent Stem-Succulents                   
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  70. Branched Stem-Succulents                      
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  71. Unbranched Stem-Succulents                    
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  72. Arborescent Grasses                           
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  73. Tall-Cane Grasses                             
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  74. Sclerophyllous Grasses                        
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  75. Tall Grasses                                  
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  76. Short-Sward Graminoids                        
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  77. Short-Bunch Graminoids                        
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  78. Desert Short-Bunch Graminoids                 
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  79. Temperate Tall Tussock-Graminoids             
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  80. Temperate Short Tussock-Graminoids            
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  81. Tropical-Montane Short Tussock-Graminoids     
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  82. Cool-Maritime Short Tussock-Graminoids        
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  83. Tundra Short Tussock-Graminoids               
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  84. Tropical Evergreen Forbs                      
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  85. Temperate Evergreen Forbs                     
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  86. Succulent Forbs                               
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  87. Raingreen Forbs                               
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  88. Summergreen Forbs                             
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/*  89. Desert Ephemeral Herbs                        
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  90. Tropical-Alpine Raingreen Herbs               
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  91. Tundra Summergreen Herbs                      
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  92. Desert Cushion-Herbs                          
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/*  93. Tropical-Montane Cushion-Herbs                
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  94. Tundra Summergreen Cushion-Herbs              
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/*  95. Tropical Liana Evergreen Vines                
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  96. Evergreen Vines                               
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/*  97. Raingreen Vines                               
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  98. Summergreen Vines                             
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/*  99. Tropical Broadleaf Epiphytes                  
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/* 100. Broadleaf Epiphytes                           
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/* 101. Narrowleaf Epiphytes                          
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/* 102. Evergreen Ferns                               
*/
              6.493506494,   0.004278002,    15.178831169,   36.630036630,

/* 103. Summergreen Ferns                             
*/
              5.555555556,   0.905209223,    65.055555556,   28.011204482,

/* 104. Mesic Peat-Forming Bryophytes                 
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/* 105. Mesic Bryophytes                              
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/* 106. Epiphytic Bryophytes                          
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/* 107. Mesic Lichens                                 
*/
              9.090909091,   0.403419816,    24.562727273,   42.016806723,

/* 108. Boreal Lichens                                
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364,

/* 109. Desert Cryptogams                             
*/
              1.169590643,   0.000106386,   109.938654971,   36.363636364,

/* 110. Polar Cryptogams                              
*/
              5.555555556,   0.901632339,    65.005555556,   36.363636364

	};

