/* Life-Form Plant Community Levels */



# define		num_trees          29
# define		num_vines           4
# define		num_epiphytes       4
# define		num_shrubs         45
# define		num_grasses        28


enum canopy_levels	{ tree = 1, vine, epiphyte, shrub, grass };



enum canopy_levels

		lf_pcl [110] =
		{
         1,    /*    1. Equatorial Rainforest Broadleaf Trees           */
         1,    /*    2. Equatorial-Montane Rainforest Broadleaf Trees   */
         1,    /*    3. Tropical Evergreen Microphyll-Broadleaf Trees   */
         1,    /*    4. Tropical-Monsoon Evergreen Broadleaf Trees      */
         1,    /*    5. Tropical Evergreen Sclerophyll-Broadleaf Trees  */
         1,    /*    6. Mediterranean Evergreen Broadleaf Trees         */
         1,    /*    7. Warm-Temperate-Mesic Evergreen Broadleaf Trees  */
         1,    /*    8. Warm-Temperate Evergreen Broadleaf Trees        */
         1,    /*    9. Tropical-Monsoon Raingreen Broadleaf Trees      */
         1,    /*   10. Tropical-Montane Raingreen Broadleaf Trees      */
         1,    /*   11. Tropical-Xeric Raingreen Broadleaf Trees        */
         1,    /*   12. Equatorial-Xeric Raingreen Broadleaf Trees      */
         1,    /*   13. Western-Temperate Summergreen Broadleaf Trees   */
         1,    /*   14. Eastern-Temperate Summergreen Broadleaf Trees   */
         1,    /*   15. Eastern-Temp-Xeric Summergreen Broadleaf Trees  */
         1,    /*   16. Cool-Temp/Boreal Summergreen Broadleaf Trees    */
         1,    /*   17. Tropical Evergreen Linearleaf Trees             */
         1,    /*   18. Tropical-Xeric Evergreen Needleleaf Trees       */
         1,    /*   19. Mediterranean Evergreen Needleleaf Trees        */
         1,    /*   20. Heliophilic Evergreen Needleleaf Trees          */
         1,    /*   21. Warm-Temp-Montane Evergreen Needleleaf Trees    */
         1,    /*   22. Temperate-Mesic Evergreen Needleleaf Trees      */
         1,    /*   23. Western-Temperate Evergreen Needleleaf Trees    */
         1,    /*   24. Eastern-Temperate Evergreen Needleleaf Trees    */
         1,    /*   25. Cool-Temp/Boreal Evergreen Needleleaf Trees     */
         1,    /*   26. Cool-Temp/Boreal Summergreen Needleleaf Trees   */
         1,    /*   27. Hydrophilic Summergreen Needleleaf Trees        */
         4,    /*   28. Tropical Evergreen Broadleaf Small-Trees        */
         4,    /*   29. Warm-Temperate Evergreen Broadleaf Small-Trees  */
         4,    /*   30. Cool-Maritime Evergreen Broadleaf Small-Trees   */
         4,    /*   31. Tropical Raingreen Broadleaf Small-Trees        */
         4,    /*   32. Temperate Summergreen Broadleaf Small-Trees     */
         4,    /*   33. Temperate Evergreen Needleleaf Small-Trees      */
         4,    /*   34. Tropical Evergreen Broadleaf Dwarf-Trees        */
         4,    /*   35. Tropical-Montane Evergreen Dwarf-Trees          */
         1,    /*   36. Temp-Montane/Boreal Needleleaf Dwarf-Trees      */
         1,    /*   37. Palmiform Rosette-Trees                         */
         4,    /*   38. Palmiform Small-Rosette-Trees                   */
         4,    /*   39. Fern Small-Rosette-Trees                        */
         4,    /*   40. Tropical-Montane Small-Rosette-Trees            */
         4,    /*   41. Xeric Small-Rosette-Trees                       */
         4,    /*   42. Raingreen Broadleaf Arborescent-Shrubs          */
         4,    /*   43. Evergreen Broadleaf Arborescent-Shrubs          */
         4,    /*   44. Summergreen Broadleaf Arborescent-Shrubs        */
         4,    /*   45. Xeric Leafless Arborescent-Shrubs               */
         4,    /*   46. Tropical Evergreen Broadleaf Shrubs             */
         4,    /*   47. Hot-Desert Evergreen Broadleaf Shrubs           */
         4,    /*   48. Desert Evergreen Leaf-Succulent Shrubs          */
         4,    /*   49. Mediterranean Evergreen Broadleaf Shrubs        */
         4,    /*   50. Mesic-Montane Evergreen Broadleaf Shrubs        */
         4,    /*   51. Warm-Temperate Evergreen Broadleaf Shrubs       */
         4,    /*   52. Cool-Temperate-Xeric Evergreen Shrubs           */
         4,    /*   53. Temperate Summergreen Broadleaf Shrubs          */
         4,    /*   54. Cool-Temp/Boreal Summergreen Broadleaf Shrubs   */
         4,    /*   55. Temperate-Xeric Summergreen Broadleaf Shrubs    */
         4,    /*   56. Mediterranean Summergreen Broadleaf Shrubs      */
         4,    /*   57. Mediterranean Evergreen Needleleaf Shrubs       */
         4,    /*   58. Evergreen Needleleaf Shrubs                     */
         4,    /*   59. Xeric Dwarf-Shrubs                              */
         4,    /*   60. Mediterranean Dwarf-Shrubs                      */
         4,    /*   61. Cool-Maritime Evergreen Dwarf-Shrubs            */
         4,    /*   62. Cool-Temp/Boreal Evergreen Dwarf-Shrubs         */
         4,    /*   63. Tundra Evergreen Dwarf-Shrubs                   */
         4,    /*   64. Tundra Summergreen Dwarf-Shrubs                 */
         4,    /*   65. Cool-Maritime/Tropical-Alpine Cushion-Shrubs    */
         4,    /*   66. Xeric Cushion-Shrubs                            */
         4,    /*   67. Mesic Palmiform Rosette-Shrubs                  */
         4,    /*   68. Xeric Rosette-Shrubs                            */
         4,    /*   69. Arborescent Stem-Succulents                     */
         4,    /*   70. Branched Stem-Succulents                        */
         4,    /*   71. Unbranched Stem-Succulents                      */
         4,    /*   72. Arborescent Grasses                             */
         4,    /*   73. Tall-Cane Grasses                               */
         5,    /*   74. Sclerophyllous Grasses                          */
         4,    /*   75. Tall Grasses                                    */
         5,    /*   76. Short-Sward Graminoids                          */
         5,    /*   77. Short-Bunch Graminoids                          */
         5,    /*   78. Desert Short-Bunch Graminoids                   */
         5,    /*   79. Temperate Tall Tussock-Graminoids               */
         5,    /*   80. Temperate Short Tussock-Graminoids              */
         5,    /*   81. Tropical-Montane Short Tussock-Graminoids       */
         5,    /*   82. Cool-Maritime Short Tussock-Graminoids          */
         5,    /*   83. Tundra Short Tussock-Graminoids                 */
         5,    /*   84. Tropical Evergreen Forbs                        */
         5,    /*   85. Temperate Evergreen Forbs                       */
         5,    /*   86. Succulent Forbs                                 */
         5,    /*   87. Raingreen Forbs                                 */
         5,    /*   88. Summergreen Forbs                               */
         5,    /*   89. Desert Ephemeral Herbs                          */
         5,    /*   90. Tropical-Alpine Raingreen Herbs                 */
         5,    /*   91. Tundra Summergreen Herbs                        */
         5,    /*   92. Desert Cushion-Herbs                            */
         5,    /*   93. Tropical-Montane Cushion-Herbs                  */
         5,    /*   94. Tundra Summergreen Cushion-Herbs                */
         2,    /*   95. Tropical Liana Evergreen Vines                  */
         2,    /*   96. Evergreen Vines                                 */
         2,    /*   97. Raingreen Vines                                 */
         2,    /*   98. Summergreen Vines                               */
         3,    /*   99. Tropical Broadleaf Epiphytes                    */
         3,    /*  100. Broadleaf Epiphytes                             */
         3,    /*  101. Narrowleaf Epiphytes                            */
         4,    /*  102. Evergreen Ferns                                 */
         4,    /*  103. Summergreen Ferns                               */
         5,    /*  104. Mesic Peat-Forming Bryophytes                   */
         5,    /*  105. Mesic Bryophytes                                */
         3,    /*  106. Epiphytic Bryophytes                            */
         5,    /*  107. Mesic Lichens                                   */
         5,    /*  108. Boreal Lichens                                  */
         5,    /*  109. Desert Cryptogams                               */
         5     /*  110. Polar Cryptogams                                */

		};

