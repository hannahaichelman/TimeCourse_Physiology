# TimeCourse_Physiology
Data and scripts associated with Aichelman et al. (in revision, Limnology and Oceanography). 

These scripts are used to analyzing species and reef-zone specific responses to global change stressors through time. We used a 95-day laboratory experiment to investigate how treatment duration influenced physiological responses of two Caribbean reef-building coral species (Siderastrea siderea and Pseudodiploria strigosa) from two reef zones on the Belize Mesoamerican Barrier Reef System under ocean warming (28, 31°C), acidification (~pCO2 400–2800 µatm), and their interactions. Calcification rate, total host protein and carbohydrate, chlorophyll a pigment, and symbiont cell density were quantified every 30 days from fragments of the same coral colony to characterize acclimatory responses of each coral.

We analyzed the data both considering the response of individual physiology parameters as well as a Principal Components Analysis (PCA) approach to consider shifts in holobiont physiology through time.

Contact Info: hannahaichelman@gmail.com


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
 [1] gamlss_5.2-0       nlme_3.1-148       gamlss.dist_5.1-7  gamlss.data_5.1-4  goft_1.3.6         sn_1.6-2          
 [7] MuMIn_1.43.17      corrplot_0.84      lmerTest_3.1-3     multcomp_1.4-13    TH.data_1.0-10     mvtnorm_1.1-1     
[13] car_3.0-8          carData_3.0-4      lme4_1.1-23        Matrix_1.2-18      cowplot_1.0.0      plotly_4.9.2.1    
[19] forcats_0.5.0      stringr_1.4.0      dplyr_1.0.0        purrr_0.3.4        readr_1.3.1        tidyr_1.1.0       
[25] tibble_3.0.1       tidyverse_1.3.0    Rmisc_1.5          plyr_1.8.6         lsmeans_2.30-0     emmeans_1.5.2-1   
[31] Hmisc_4.4-0        Formula_1.2-3      lattice_0.20-41    chron_2.3-56       lubridate_1.7.9    plotrix_3.7-8     
[37] segmented_1.2-0    ggplot2_3.3.2      fitdistrplus_1.1-3 survival_3.2-3     MASS_7.3-51.6     

loaded via a namespace (and not attached):
 [1] minqa_1.2.4         colorspace_1.4-1    ellipsis_0.3.1      rio_0.5.16          estimability_1.3   
 [6] htmlTable_2.0.1     base64enc_0.1-3     fs_1.4.1            rstudioapi_0.11     fansi_0.4.1        
[11] xml2_1.3.2          codetools_0.2-16    mnormt_2.0.2        knitr_1.29          jsonlite_1.7.1     
[16] nloptr_1.2.2.2      broom_0.5.6         cluster_2.1.0       dbplyr_1.4.4        png_0.1-7          
[21] compiler_4.0.1      httr_1.4.2          backports_1.1.8     assertthat_0.2.1    lazyeval_0.2.2     
[26] cli_2.0.2           acepack_1.4.1       htmltools_0.5.0     tools_4.0.1         coda_0.19-3        
[31] gtable_0.3.0        glue_1.4.1          tinytex_0.24        Rcpp_1.0.4.6        cellranger_1.1.0   
[36] vctrs_0.3.1         xfun_0.15           openxlsx_4.1.5      rvest_0.3.5         lifecycle_0.2.0    
[41] statmod_1.4.34      zoo_1.8-8           scales_1.1.1        hms_0.5.3           sandwich_2.5-1     
[46] RColorBrewer_1.1-2  yaml_2.2.1          curl_4.3            gridExtra_2.3       rpart_4.1-15       
[51] latticeExtra_0.6-29 stringi_1.4.6       checkmate_2.0.0     boot_1.3-25         zip_2.0.4          
[56] rlang_0.4.8         pkgconfig_2.0.3     htmlwidgets_1.5.1   tidyselect_1.1.0    magrittr_1.5       
[61] R6_2.4.1            generics_0.0.2      DBI_1.1.0           pillar_1.4.4        haven_2.3.1        
[66] foreign_0.8-80      withr_2.2.0         abind_1.4-5         sp_1.4-2            nnet_7.3-14        
[71] modelr_0.1.8        crayon_1.3.4        utf8_1.1.4          tmvnsim_1.0-2       jpeg_0.1-8.1       
[76] grid_4.0.1          readxl_1.3.1        data.table_1.12.8   blob_1.2.1          reprex_0.3.0       
[81] digest_0.6.25       xtable_1.8-4        numDeriv_2016.8-1.1 munsell_0.5.0       viridisLite_0.3.0 
```
