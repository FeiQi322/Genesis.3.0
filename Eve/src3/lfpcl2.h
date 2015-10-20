/* Life-Form Plant Community Levels:  2 - level scheme */



# define		num_uc             55
# define		num_lc             55


/* upper canopy = 1 */
/* lower canopy = 0 */


enum canopy_levels_2		{ lower = 0, upper };



enum canopy_levels_2

		lf_pcl_2 [110] =
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
         1,    /*   28. Tropical Evergreen Broadleaf Small-Trees        */
         1,    /*   29. Warm-Temperate Evergreen Broadleaf Small-Trees  */
         1,    /*   30. Cool-Maritime Evergreen Broadleaf Small-Trees   */
         1,    /*   31. Tropical Raingreen Broadleaf Small-Trees        */
         1,    /*   32. Temperate Summergreen Broadleaf Small-Trees     */
         1,    /*   33. Temperate Evergreen Needleleaf Small-Trees      */
         1,    /*   34. Tropical Evergreen Broadleaf Dwarf-Trees        */
         1,    /*   35. Tropical-Montane Evergreen Dwarf-Trees          */
         1,    /*   36. Temp-Montane/Boreal Needleleaf Dwarf-Trees      */
         1,    /*   37. Palmiform Rosette-Trees                         */
         1,    /*   38. Palmiform Small-Rosette-Trees                   */
         1,    /*   39. Fern Small-Rosette-Trees                        */
         1,    /*   40. Tropical-Montane Small-Rosette-Trees            */
         1,    /*   41. Xeric Small-Rosette-Trees                       */
         1,    /*   42. Raingreen Broadleaf Arborescent-Shrubs          */
         1,    /*   43. Evergreen Broadleaf Arborescent-Shrubs          */
         1,    /*   44. Summergreen Broadleaf Arborescent-Shrubs        */
         1,    /*   45. Xeric Leafless Arborescent-Shrubs               */
         0,    /*   46. Tropical Evergreen Broadleaf Shrubs             */
         0,    /*   47. Hot-Desert Evergreen Broadleaf Shrubs           */
         0,    /*   48. Desert Evergreen Leaf-Succulent Shrubs          */
         0,    /*   49. Mediterranean Evergreen Broadleaf Shrubs        */
         0,    /*   50. Mesic-Montane Evergreen Broadleaf Shrubs        */
         0,    /*   51. Warm-Temperate Evergreen Broadleaf Shrubs       */
         0,    /*   52. Cool-Temperate-Xeric Evergreen Shrubs           */
         0,    /*   53. Temperate Summergreen Broadleaf Shrubs          */
         0,    /*   54. Cool-Temp/Boreal Summergreen Broadleaf Shrubs   */
         0,    /*   55. Temperate-Xeric Summergreen Broadleaf Shrubs    */
         0,    /*   56. Mediterranean Summergreen Broadleaf Shrubs      */
         0,    /*   57. Mediterranean Evergreen Needleleaf Shrubs       */
         0,    /*   58. Evergreen Needleleaf Shrubs                     */
         0,    /*   59. Xeric Dwarf-Shrubs                              */
         0,    /*   60. Mediterranean Dwarf-Shrubs                      */
         0,    /*   61. Cool-Maritime Evergreen Dwarf-Shrubs            */
         0,    /*   62. Cool-Temp/Boreal Evergreen Dwarf-Shrubs         */
         0,    /*   63. Tundra Evergreen Dwarf-Shrubs                   */
         0,    /*   64. Tundra Summergreen Dwarf-Shrubs                 */
         0,    /*   65. Cool-Maritime/Tropical-Alpine Cushion-Shrubs    */
         0,    /*   66. Xeric Cushion-Shrubs                            */
         0,    /*   67. Mesic Palmiform Rosette-Shrubs                  */
         0,    /*   68. Xeric Rosette-Shrubs                            */
         1,    /*   69. Arborescent Stem-Succulents                     */
         0,    /*   70. Branched Stem-Succulents                        */
         0,    /*   71. Unbranched Stem-Succulents                      */
         1,    /*   72. Arborescent Grasses                             */
         0,    /*   73. Tall-Cane Grasses                               */
         0,    /*   74. Sclerophyllous Grasses                          */
         0,    /*   75. Tall Grasses                                    */
         0,    /*   76. Short-Sward Graminoids                          */
         0,    /*   77. Short-Bunch Graminoids                          */
         0,    /*   78. Desert Short-Bunch Graminoids                   */
         0,    /*   79. Temperate Tall Tussock-Graminoids               */
         0,    /*   80. Temperate Short Tussock-Graminoids              */
         0,    /*   81. Tropical-Montane Short Tussock-Graminoids       */
         0,    /*   82. Cool-Maritime Short Tussock-Graminoids          */
         0,    /*   83. Tundra Short Tussock-Graminoids                 */
         0,    /*   84. Tropical Evergreen Forbs                        */
         0,    /*   85. Temperate Evergreen Forbs                       */
         0,    /*   86. Succulent Forbs                                 */
         0,    /*   87. Raingreen Forbs                                 */
         0,    /*   88. Summergreen Forbs                               */
         0,    /*   89. Desert Ephemeral Herbs                          */
         0,    /*   90. Tropical-Alpine Raingreen Herbs                 */
         0,    /*   91. Tundra Summergreen Herbs                        */
         0,    /*   92. Desert Cushion-Herbs                            */
         0,    /*   93. Tropical-Montane Cushion-Herbs                  */
         0,    /*   94. Tundra Summergreen Cushion-Herbs                */
         1,    /*   95. Tropical Liana Evergreen Vines                  */
         1,    /*   96. Evergreen Vines                                 */
         1,    /*   97. Raingreen Vines                                 */
         1,    /*   98. Summergreen Vines                               */
         1,    /*   99. Tropical Broadleaf Epiphytes                    */
         1,    /*  100. Broadleaf Epiphytes                             */
         1,    /*  101. Narrowleaf Epiphytes                            */
         0,    /*  102. Evergreen Ferns                                 */
         0,    /*  103. Summergreen Ferns                               */
         0,    /*  104. Mesic Peat-Forming Bryophytes                   */
         0,    /*  105. Mesic Bryophytes                                */
         1,    /*  106. Epiphytic Bryophytes                            */
         0,    /*  107. Mesic Lichens                                   */
         0,    /*  108. Boreal Lichens                                  */
         0,    /*  109. Desert Cryptogams                               */
         0     /*  110. Polar Cryptogams                                */

		};

