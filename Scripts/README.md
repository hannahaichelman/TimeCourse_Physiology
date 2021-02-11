# Scripts

TimeCoursePhysiologyAnalysis_github.R is used to analyze changes in individual coral physiology metrics through time as well as correlations of physiology metrics through time. This script was used to make Figures 2, 3, 4, 6, S3 and S4.

PCAs_github.R is used to explore and plot PCAs that consider how holobiont physiology changes through time. This script was used to make Figure 5 and Figure S2. 

BW_DW_regression_github.R is used to check the correlations between buoyant weight and dry weight for the coral fragments, in order to later calculate calcification rate. 


R Session Info for most recent analysis:
``` r
> sessionInfo()
R version 4.0.1 (2020-06-06)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] parallel  splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vegan_2.5-6                permute_0.9-5              PerformanceAnalytics_2.0.4 xts_0.12.1                
 [5] zoo_1.8-8                  factoextra_1.0.7           FactoMineR_2.3             cluster_2.1.0             
 [9] ggfortify_0.4.11           ggpubr_0.4.0               readxl_1.3.1               gamlss_5.2-0              
[13] nlme_3.1-148               gamlss.dist_5.1-7          gamlss.data_5.1-4          goft_1.3.6                
[17] sn_1.6-2                   fitdistrplus_1.1-3         MuMIn_1.43.17              corrplot_0.84             
[21] lmerTest_3.1-3             multcomp_1.4-13            TH.data_1.0-10             mvtnorm_1.1-1             
[25] car_3.0-8                  carData_3.0-4              lme4_1.1-23                Matrix_1.2-18             
[29] cowplot_1.0.0              plotly_4.9.2.1             forcats_0.5.0              stringr_1.4.0             
[33] dplyr_1.0.0                purrr_0.3.4                readr_1.3.1                tidyr_1.1.0               
[37] tibble_3.0.1               tidyverse_1.3.0            Rmisc_1.5                  plyr_1.8.6                
[41] MASS_7.3-51.6              lsmeans_2.30-0             emmeans_1.5.2-1            Hmisc_4.4-0               
[45] Formula_1.2-3              survival_3.2-3             lattice_0.20-41            chron_2.3-56              
[49] lubridate_1.7.9            plotrix_3.7-8              segmented_1.2-0            ggplot2_3.3.2             

loaded via a namespace (and not attached):
 [1] backports_1.1.8      lazyeval_0.2.2       sp_1.4-2             digest_0.6.25        htmltools_0.5.0     
 [6] fansi_0.4.1          magrittr_1.5         checkmate_2.0.0      openxlsx_4.1.5       modelr_0.1.8        
[11] sandwich_2.5-1       jpeg_0.1-8.1         colorspace_1.4-1     blob_1.2.1           rvest_0.3.5         
[16] ggrepel_0.8.2        haven_2.3.1          xfun_0.15            crayon_1.3.4         jsonlite_1.7.1      
[21] glue_1.4.1           gtable_0.3.0         abind_1.4-5          scales_1.1.1         DBI_1.1.0           
[26] rstatix_0.6.0        Rcpp_1.0.4.6         viridisLite_0.3.0    xtable_1.8-4         htmlTable_2.0.1     
[31] tmvnsim_1.0-2        flashClust_1.01-2    foreign_0.8-80       htmlwidgets_1.5.1    httr_1.4.2          
[36] RColorBrewer_1.1-2   acepack_1.4.1        ellipsis_0.3.1       pkgconfig_2.0.3      nnet_7.3-14         
[41] dbplyr_1.4.4         utf8_1.1.4           tidyselect_1.1.0     rlang_0.4.8          munsell_0.5.0       
[46] cellranger_1.1.0     tools_4.0.1          cli_2.0.2            generics_0.0.2       broom_0.5.6         
[51] yaml_2.2.1           knitr_1.29           fs_1.4.1             zip_2.0.4            leaps_3.1           
[56] xml2_1.3.2           compiler_4.0.1       rstudioapi_0.11      curl_4.3             png_0.1-7           
[61] ggsignif_0.6.0       reprex_0.3.0         statmod_1.4.34       stringi_1.4.6        nloptr_1.2.2.2      
[66] vctrs_0.3.1          pillar_1.4.4         lifecycle_0.2.0      estimability_1.3     data.table_1.12.8   
[71] R6_2.4.1             latticeExtra_0.6-29  gridExtra_2.3        rio_0.5.16           codetools_0.2-16    
[76] boot_1.3-25          assertthat_0.2.1     withr_2.2.0          mnormt_2.0.2         mgcv_1.8-31         
[81] hms_0.5.3            quadprog_1.5-8       grid_4.0.1           rpart_4.1-15         coda_0.19-3         
[86] minqa_1.2.4          numDeriv_2016.8-1.1  scatterplot3d_0.3-41 base64enc_0.1-3      tinytex_0.24  
```
