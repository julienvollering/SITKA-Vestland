# Tuning hyperparameters by spatially-stratified CV

library(tidyverse)
library(terra)
library(sf)

load(file = "SDM-3_tuningdata.RData")
filterfactors <- 2^(seq.int(0,7))

## Hyperparameter design:

design <- tidyr::expand_grid(sp = c('ps', 'th'),
                             filtering = filterfactors,
                             regularization = 2^(seq.int(0,3)),
                             fc = c("lqph", "lqp", "lq"),
                             testfold = 1:5)

## Define one iteration of hyperparameter testing

#sp <- "th"; filtering <- 1; regularization <- 1; fc <- 'lqph'; testfold <- 1; df <- df
calc_fold_auc <- function(sp, filtering, regularization, fc, testfold, df) {
  
  dftrain <- df %>% 
    filter(species == sp, fold != testfold)
  bkgrows <- which(pull(dftrain, paste0("presmesh", filtering)) == 0)
  set.seed(42)
  dftrain <- dftrain[-sample(bkgrows, length(bkgrows)-1e4),]

  p <- dftrain %>% 
    pull(paste0("presmesh", filtering))
  data <- dftrain %>% 
    select(starts_with("CHELSA"))
  mod <- maxnet::maxnet(p = p, data = data, regmult = regularization, 
                        maxnet::maxnet.formula(p, data, classes = fc))
  
  dftest <- df %>% 
    filter(species == sp, fold == testfold)
  raw <- predict(mod, dftest, clamp = FALSE, type = "exponential")[,1]
  cont <- as.matrix(table(raw, pull(dftest, presmesh1)))
  cont <- cont[order(as.numeric(rownames(cont)), decreasing = T), ]
  falspos <- c(0, unname(cumsum(cont[, "0"])))
  truepos <- c(0, unname(cumsum(cont[, "1"])))
  fpr <- falspos/sum(cont[, "0"])
  tpr <- truepos/sum(cont[, "1"])
  
  hgtl <- tpr[-length(tpr)]
  hgtr <- tpr[-1]
  wdth <- diff(fpr)
  auc <- sum(((hgtl + hgtr)/2) * wdth)
  return(auc)
}
# profvis::profvis({
#   calc_fold_auc(sp = "ps", filtering = 4, regularization = 4, fc = "lq", testfold = 1, df = df)
# })

## Iterate

possibly_calc_fold_auc <- possibly(calc_fold_auc, otherwise = NA)

write_fold_auc <- function(sp, filtering, regularization, fc, testfold, df, 
                           filename) {
  auc <- possibly_calc_fold_auc(sp, filtering, regularization, fc, testfold, df)
  resrow <- data.frame(sp = sp, filtering = filtering, 
                       regularization = regularization, fc = fc, 
                       testfold = testfold, auc = auc)
  write.table(resrow, file = filename, append = TRUE, sep = ',', 
              row.names = FALSE, col.names = !file.exists(filename))
}

library(furrr)
plan(multisession, workers = 4)
future_pwalk(design, write_fold_auc, df = df, filename = "SDM-4_tuningmaxnet.csv")

sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
# [3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] sf_1.0-8        terra_1.6-17    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.10    purrr_0.3.4     readr_2.1.2    
# [8] tidyr_1.2.0     tibble_3.1.8    ggplot2_3.3.6   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.2   haven_2.5.0        colorspace_2.0-3   vctrs_0.4.2        generics_0.1.3    
# [6] utf8_1.2.2         rlang_1.0.6        e1071_1.7-11       pillar_1.8.1       glue_1.6.2        
# [11] withr_2.5.0        DBI_1.1.3          dbplyr_2.1.1       modelr_0.1.8       readxl_1.4.0      
# [16] lifecycle_1.0.3    munsell_0.5.0      gtable_0.3.1       cellranger_1.1.0   rvest_1.0.2       
# [21] codetools_0.2-18   maxnet_0.1.4       tzdb_0.3.0         class_7.3-20       fansi_1.0.3       
# [26] broom_0.8.0        Rcpp_1.0.9         KernSmooth_2.23-20 scales_1.2.1       backports_1.4.1   
# [31] classInt_0.4-7     jsonlite_1.8.2     fs_1.5.2           hms_1.1.1          stringi_1.7.6     
# [36] grid_4.1.3         cli_3.4.1          tools_4.1.3        magrittr_2.0.3     proxy_0.4-27      
# [41] crayon_1.5.1       pkgconfig_2.0.3    ellipsis_0.3.2     xml2_1.3.3         reprex_2.0.1      
# [46] lubridate_1.8.0    assertthat_0.2.1   httr_1.4.3         rstudioapi_0.13    R6_2.5.1          
# [51] units_0.8-0        compiler_4.1.3    
