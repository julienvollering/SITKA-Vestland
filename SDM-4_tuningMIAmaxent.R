# Tuning hyperparameters by spatially-stratified CV

library(tidyverse)
library(terra)
library(sf)

load(file = "SDM-3_tuningdata.RData")
filterfactors <- 2^(seq.int(0,7))

## Hyperparameter design:

design <- tidyr::expand_grid(sp = c('ps', 'th'),
                             filtering = filterfactors,
                             transformtypes = c("LDHfHr", "LD"),
                             interact = c(TRUE, FALSE),
                             a = 10^(c(-3,-6,-9,-12)),
                             testfold = 1:5)
designdone <- read_csv("SDM-4_tuningMIAmaxent.csv")
design <- bind_rows(designdone, design) %>% 
  distinct(sp, filtering, transformtypes, interact, a, testfold) %>% 
  slice(-(1:nrow(designdone)))


## Define one iteration of hyperparameter testing

#sp <- "ps"; filtering <- 128; transformtypes <- "LDHfHr"; interact = TRUE; a <- 1e-12; testfold <- 4
calc_fold_auc <- function(sp, filtering, transformtypes, interact, a, testfold, df) {
  
  dftrain <- df %>% 
    filter(species == sp, fold != testfold)
  bkgrows <- which(pull(dftrain, paste0("presmesh", filtering)) == 0)
  set.seed(42)
  dftrain <- dftrain[-sample(bkgrows, length(bkgrows)-1e4),]
  dftrain <- select(dftrain, paste0("presmesh", filtering), starts_with("CHELSA"))
  
  if (transformtypes == "LD") {
    tt <- c("L", "D")
  } else if (transformtypes == "LDHfHr") {
    tt <- c("L", "D", "HF", "HR")
  }
  dv <- MIAmaxent::deriveVars(dftrain, transformtype = tt, quiet = TRUE)
  seldv <- MIAmaxent::selectDVforEV(dv$dvdata, alpha = a, quiet = TRUE)
  selev <- MIAmaxent::selectEV(seldv$dvdata, alpha = a, interaction = interact, quiet = TRUE)
  
  dftest <- df %>% 
    filter(species == sp, fold == testfold)
  PRO <- MIAmaxent::projectModel(selev$selectedmodel, dv$transformations, dftest)$output[,1]
  cont <- as.matrix(table(PRO, pull(dftest, presmesh1)))
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
# calc_fold_auc(sp = "th", filtering = 32, transformtypes = "LD", interact=TRUE, a = 1e-12, testfold = 2, df)

## Iterate

possibly_calc_fold_auc <- possibly(calc_fold_auc, otherwise = NA)

write_fold_auc <- function(sp, filtering, transformtypes, interact, a, testfold, df, 
                           filename) {
  auc <- possibly_calc_fold_auc(sp, filtering, transformtypes, interact, a, testfold, df)
  resrow <- data.frame(sp = sp, filtering = filtering, 
                       transformtypes = transformtypes, interact = interact, a = a, 
                       testfold = testfold, auc = auc)
  write.table(resrow, file = filename, append = TRUE, sep = ',', 
              row.names = FALSE, col.names = !file.exists(filename))
}

library(furrr)
plan(multisession, workers = 4)
future_pwalk(design, write_fold_auc, df = df, filename = "SDM-4_tuningMIAmaxent.csv")

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
#   [1] tidyselect_1.1.2     haven_2.5.0          colorspace_2.0-3     vctrs_0.4.2          generics_0.1.3      
# [6] utf8_1.2.2           rlang_1.0.6          e1071_1.7-11         pillar_1.8.1         glue_1.6.2          
# [11] withr_2.5.0          DBI_1.1.3            bit64_4.0.5          dbplyr_2.1.1         modelr_0.1.8        
# [16] readxl_1.4.0         lifecycle_1.0.3      munsell_0.5.0        gtable_0.3.1         cellranger_1.1.0    
# [21] rvest_1.0.2          codetools_0.2-18     MIAmaxent_1.2.0.9000 maxnet_0.1.4         tzdb_0.3.0          
# [26] parallel_4.1.3       class_7.3-20         fansi_1.0.3          broom_0.8.0          Rcpp_1.0.9          
# [31] KernSmooth_2.23-20   scales_1.2.1         backports_1.4.1      classInt_0.4-7       vroom_1.5.7         
# [36] jsonlite_1.8.2       bit_4.0.4            fs_1.5.2             hms_1.1.1            stringi_1.7.6       
# [41] grid_4.1.3           cli_3.4.1            tools_4.1.3          magrittr_2.0.3       proxy_0.4-27        
# [46] crayon_1.5.1         pkgconfig_2.0.3      ellipsis_0.3.2       xml2_1.3.3           reprex_2.0.1        
# [51] lubridate_1.8.0      assertthat_0.2.1     httr_1.4.3           rstudioapi_0.13      R6_2.5.1            
# [56] units_0.8-0          compiler_4.1.3    