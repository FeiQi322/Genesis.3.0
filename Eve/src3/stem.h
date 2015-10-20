/* Life-Form Stem condition */



# define		num_ws             76
# define		num_hs             22
# define		num_nghs           12


/* wood stem       = 1                     */
/* herbaceous stem = 0, 2 (nongrassy herb) */


enum stem_condition		{ herb = 0, wood, ngherb };



enum stem_condition

		lf_stem [110] =
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
         1,    /*   46. Tropical Evergreen Broadleaf Shrubs             */
         1,    /*   47. Hot-Desert Evergreen Broadleaf Shrubs           */
         1,    /*   48. Desert Evergreen Leaf-Succulent Shrubs          */
         1,    /*   49. Mediterranean Evergreen Broadleaf Shrubs        */
         1,    /*   50. Mesic-Montane Evergreen Broadleaf Shrubs        */
         1,    /*   51. Warm-Temperate Evergreen Broadleaf Shrubs       */
         1,    /*   52. Cool-Temperate-Xeric Evergreen Shrubs           */
         1,    /*   53. Temperate Summergreen Broadleaf Shrubs          */
         1,    /*   54. Cool-Temp/Boreal Summergreen Broadleaf Shrubs   */
         1,    /*   55. Temperate-Xeric Summergreen Broadleaf Shrubs    */
         1,    /*   56. Mediterranean Summergreen Broadleaf Shrubs      */
         1,    /*   57. Mediterranean Evergreen Needleleaf Shrubs       */
         1,    /*   58. Evergreen Needleleaf Shrubs                     */
         1,    /*   59. Xeric Dwarf-Shrubs                              */
         1,    /*   60. Mediterranean Dwarf-Shrubs                      */
         1,    /*   61. Cool-Maritime Evergreen Dwarf-Shrubs            */
         1,    /*   62. Cool-Temp/Boreal Evergreen Dwarf-Shrubs         */
         1,    /*   63. Tundra Evergreen Dwarf-Shrubs                   */
         1,    /*   64. Tundra Summergreen Dwarf-Shrubs                 */
         1,    /*   65. Cool-Maritime/Tropical-Alpine Cushion-Shrubs    */
         1,    /*   66. Xeric Cushion-Shrubs                            */
         1,    /*   67. Mesic Palmiform Rosette-Shrubs                  */
         1,    /*   68. Xeric Rosette-Shrubs                            */
         1,    /*   69. Arborescent Stem-Succulents                     */
         1,    /*   70. Branched Stem-Succulents                        */
         1,    /*   71. Unbranched Stem-Succulents                      */
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
         2,    /*   99. Tropical Broadleaf Epiphytes                    */
         2,    /*  100. Broadleaf Epiphytes                             */
         2,    /*  101. Narrowleaf Epiphytes                            */
         2,    /*  102. Evergreen Ferns                                 */
         2,    /*  103. Summergreen Ferns                               */
         2,    /*  104. Mesic Peat-Forming Bryophytes                   */
         2,    /*  105. Mesic Bryophytes                                */
         2,    /*  106. Epiphytic Bryophytes                            */
         2,    /*  107. Mesic Lichens                                   */
         2,    /*  108. Boreal Lichens                                  */
         2,    /*  109. Desert Cryptogams                               */
         2     /*  110. Polar Cryptogams                                */

		};

