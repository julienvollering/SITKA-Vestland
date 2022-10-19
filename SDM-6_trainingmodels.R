# Training final models

library(tidyverse)

load(file = "SDM-5_selectingmodels.RData")
str(df)
df <- df %>% 
  rename_with(~ gsub(".1981.2010.V.2.1$", "", .x), ends_with(".1981.2010.V.2.1")) #compatible with names in projection

train <- function(sp, filtering, regularization, fc, transformtypes, interact, a, 
                  df, ...) {
  
  dftrain <- df %>% 
    filter(species == sp)
  bkgrows <- which(pull(dftrain, paste0("presmesh", filtering)) == 0)
  set.seed(42)
  dftrain <- dftrain[-sample(bkgrows, length(bkgrows)-1e4),]
  
  if (!is.na(regularization)) {
    p <- dftrain %>% 
      pull(paste0("presmesh", filtering))
    data <- dftrain %>% 
      select(starts_with("CHELSA"))
    mod <- maxnet::maxnet(p = p, data = data, regmult = regularization, 
                          maxnet::maxnet.formula(p, data, classes = fc))
    
  } else {
    dftrain <- select(dftrain, paste0("presmesh", filtering), starts_with("CHELSA"))
    if (transformtypes == "LD") {
      tt <- c("L", "D")
    } else if (transformtypes == "LDHfHr") {
      tt <- c("L", "D", "HF", "HR")
    }
    dv <- MIAmaxent::deriveVars(dftrain, transformtype = tt, quiet = TRUE)
    seldv <- MIAmaxent::selectDVforEV(dv$dvdata, alpha = a, quiet = TRUE)
    selev <- MIAmaxent::selectEV(seldv$dvdata, alpha = a, interaction = interact, quiet = TRUE)
    mod <- selev$selectedmodel
    mod$transformations <- dv$transformations
  }
  return(mod)
}

hyper <-  hyper %>%
  mutate(model = pmap(., train, df))

## Find threshold values

#hyper <- readRDS("hyper.rds")
predict_maxnet_MIAmaxent <- function(model, newdata) {
  if ("maxnet" %in% class(model)) {
    out <- predict(model, newdata, clamp = FALSE, type = "exponential")
    return(out)
  } else if ("iwlr" %in% class(model)) {
    out <- MIAmaxent::projectModel(model, model$transformations, newdata, raw = TRUE, rescale = FALSE)
    return(out$output[,1])
  }
}

find_threshold <- function(preds, pres, sensitivity = 0.9) {
  cont <- as.matrix(table(preds, pres)) 
  cont <- cont[order(as.numeric(rownames(cont)), decreasing = T), ]
  truepos <- c(0, unname(cumsum(cont[, "1"])))
  tpr <- truepos/sum(cont[, "1"])
  threshold <- as.numeric(row.names(cont)[which(tpr >= sensitivity)[1]])
  return(threshold)
}

### P.s.

dfps <- df %>% 
  filter(species == "ps")

trainpreds <- hyper %>% 
  filter(sp == "ps") %>% 
  pull(model) %>% 
  map_dfc(~ predict_maxnet_MIAmaxent(model = .x, newdata = dfps))
scalersps <- colSums(trainpreds)
trainpreds <- trainpreds %>% 
  sweep(2, scalersps, `/`) %>% 
  as_tibble()
colSums(trainpreds)

thresholdsps <- map_dbl(trainpreds, find_threshold, pres = dfps$presmesh1, sensitivity = 0.95)

trainpreds <- trainpreds %>% 
  mutate(ensemble = rowSums(trainpreds))
thresholdpsensemble <- find_threshold(trainpreds$ensemble, dfps$presmesh1, sensitivity = 0.95)

### T.h.

dfth <- df %>% 
  filter(species == "th")

trainpreds <- hyper %>% 
  filter(sp == "th") %>% 
  pull(model) %>% 
  map_dfc(~ predict_maxnet_MIAmaxent(model = .x, newdata = dfth))
scalersth <- colSums(trainpreds)
trainpreds <- trainpreds %>% 
  sweep(2, scalersth, `/`) %>% 
  as_tibble()
colSums(trainpreds)

thresholdsth <- map_dbl(trainpreds, find_threshold, pres = dfth$presmesh1, sensitivity = 0.95)

trainpreds <- trainpreds %>% 
  mutate(ensemble = rowSums(trainpreds))
thresholdthensemble <- find_threshold(trainpreds$ensemble, dfth$presmesh1, sensitivity = 0.95)

hyper <- hyper %>% 
  mutate(
    scaler = c(scalersps, scalersth),
    threshold = c(thresholdsps, thresholdsth))
save(list = c('hyper', 'thresholdpsensemble', 'thresholdthensemble'), 
     file = "SDM-6_trainingmodels.RData")

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
#   [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.10    purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.8   
# [8] ggplot2_3.3.6   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9           cellranger_1.1.0     pillar_1.8.1         compiler_4.1.3       dbplyr_2.1.1        
# [6] tools_4.1.3          lubridate_1.8.0      jsonlite_1.8.2       lifecycle_1.0.3      gtable_0.3.1        
# [11] pkgconfig_2.0.3      rlang_1.0.6          reprex_2.0.1         DBI_1.1.3            cli_3.4.1           
# [16] rstudioapi_0.13      haven_2.5.0          xml2_1.3.3           withr_2.5.0          httr_1.4.3          
# [21] fs_1.5.2             generics_0.1.3       vctrs_0.4.2          hms_1.1.1            maxnet_0.1.4        
# [26] MIAmaxent_1.2.0.9000 grid_4.1.3           tidyselect_1.1.2     glue_1.6.2           R6_2.5.1            
# [31] fansi_1.0.3          readxl_1.4.0         tzdb_0.3.0           modelr_0.1.8         magrittr_2.0.3      
# [36] backports_1.4.1      scales_1.2.1         ellipsis_0.3.2       units_0.8-0          rvest_1.0.2         
# [41] assertthat_0.2.1     colorspace_2.0-3     utf8_1.2.2           stringi_1.7.6        munsell_0.5.0       
# [46] broom_0.8.0          crayon_1.5.1        