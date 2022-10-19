# Make data frame for tuning hyperparameters by spatially-stratified CV

library(tidyverse)
library(terra)
library(sf)

pastEAland <- rast("./data/climate/pastEAland.tif")
names(pastEAland) <- make.names(names(pastEAland), allow_ = FALSE)
load("SDM-2_CVfolds.RData")

## Full data frame for P.s.

saPS <- mask(pastEAland, vect(st_union(st_buffer(psEA, dist = 200e3))))
PS <- as.data.frame(saPS, xy = TRUE, cells = TRUE)
PS <- PS %>% 
  mutate(fold = case_when(x > quintilexPS[1] & x < quintilexPS[2] ~ 1,
                          x > quintilexPS[2] & x < quintilexPS[3] ~ 2,
                          x > quintilexPS[3] & x < quintilexPS[4] ~ 3,
                          x > quintilexPS[4] & x < quintilexPS[5] ~ 4,
                          x > quintilexPS[5] & x < quintilexPS[6] ~ 5))

filterfactors <- 2^(seq.int(0,7))
meshes <- filterfactors %>% 
  map(~ aggregate(pastEAland$CHELSA.bio1.1981.2010.V.2.1, fact = .)) %>% 
  set_names(filterfactors)
psEA <- psEA %>% 
  mutate(cellnr = cells(saPS, vect(psEA))[,2],
         mesh1nr = cells(meshes$`1`, vect(psEA))[,2],
         mesh2nr = cells(meshes$`2`, vect(psEA))[,2],
         mesh4nr = cells(meshes$`4`, vect(psEA))[,2],
         mesh8nr = cells(meshes$`8`, vect(psEA))[,2],
         mesh16nr = cells(meshes$`16`, vect(psEA))[,2],
         mesh32nr = cells(meshes$`32`, vect(psEA))[,2],
         mesh64nr = cells(meshes$`64`, vect(psEA))[,2],
         mesh128nr = cells(meshes$`128`, vect(psEA))[,2])

psEA <- group_by(psEA, fold)
presmesh1cells <- group_map(psEA, function(x, ...) {
    distinct(x, mesh1nr, .keep_all = TRUE) %>% 
      pull(cellnr)
  } ) %>% 
  unlist()
presmesh2cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh2nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh4cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh4nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh8cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh8nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh16cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh16nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh32cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh32nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh64cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh64nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh128cells <- group_map(psEA, function(x, ...) {
  distinct(x, mesh128nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
PS <- mutate(PS, 
       presmesh1 = if_else(cell %in% presmesh1cells, 1, 0),
       presmesh2 = if_else(cell %in% presmesh2cells, 1, 0),
       presmesh4 = if_else(cell %in% presmesh4cells, 1, 0),
       presmesh8 = if_else(cell %in% presmesh8cells, 1, 0),
       presmesh16 = if_else(cell %in% presmesh16cells, 1, 0),
       presmesh32 = if_else(cell %in% presmesh32cells, 1, 0),
       presmesh64 = if_else(cell %in% presmesh64cells, 1, 0),
       presmesh128 = if_else(cell %in% presmesh128cells, 1, 0)) 
summary(PS)

PS <- select(PS, -cell, -x, -y)

## Full data frame for T.h.

saTH <- mask(pastEAland, vect(st_union(st_buffer(thEA, dist = 200e3))))
TH <- as.data.frame(saTH, xy = TRUE, cells = TRUE)
TH <- TH %>% 
  mutate(fold = case_when(x > quintilexTH[1] & x < quintilexTH[2] ~ 1,
                          x > quintilexTH[2] & x < quintilexTH[3] ~ 2,
                          x > quintilexTH[3] & x < quintilexTH[4] ~ 3,
                          x > quintilexTH[4] & x < quintilexTH[5] ~ 4,
                          x > quintilexTH[5] & x < quintilexTH[6] ~ 5))

filterfactors <- 2^(seq.int(0,7))
meshes <- filterfactors %>% 
  map(~ aggregate(pastEAland$CHELSA.bio1.1981.2010.V.2.1, fact = .)) %>% 
  set_names(filterfactors)
thEA <- thEA %>% 
  mutate(cellnr = cells(saTH, vect(thEA))[,2],
         mesh1nr = cells(meshes$`1`, vect(thEA))[,2],
         mesh2nr = cells(meshes$`2`, vect(thEA))[,2],
         mesh4nr = cells(meshes$`4`, vect(thEA))[,2],
         mesh8nr = cells(meshes$`8`, vect(thEA))[,2],
         mesh16nr = cells(meshes$`16`, vect(thEA))[,2],
         mesh32nr = cells(meshes$`32`, vect(thEA))[,2],
         mesh64nr = cells(meshes$`64`, vect(thEA))[,2],
         mesh128nr = cells(meshes$`128`, vect(thEA))[,2])

thEA <- group_by(thEA, fold)
presmesh1cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh1nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh2cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh2nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh4cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh4nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh8cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh8nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh16cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh16nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh32cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh32nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh64cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh64nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
presmesh128cells <- group_map(thEA, function(x, ...) {
  distinct(x, mesh128nr, .keep_all = TRUE) %>% 
    pull(cellnr)
} ) %>% 
  unlist()
TH <- mutate(TH, 
             presmesh1 = if_else(cell %in% presmesh1cells, 1, 0),
             presmesh2 = if_else(cell %in% presmesh2cells, 1, 0),
             presmesh4 = if_else(cell %in% presmesh4cells, 1, 0),
             presmesh8 = if_else(cell %in% presmesh8cells, 1, 0),
             presmesh16 = if_else(cell %in% presmesh16cells, 1, 0),
             presmesh32 = if_else(cell %in% presmesh32cells, 1, 0),
             presmesh64 = if_else(cell %in% presmesh64cells, 1, 0),
             presmesh128 = if_else(cell %in% presmesh128cells, 1, 0)) 
summary(TH)

TH <- select(TH, -cell, -x, -y)

## Combine species

df <- bind_rows(mutate(PS, species = 'ps'), mutate(TH, species = 'th'))
save(df, file = "SDM-3_tuningdata.RData")

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
# [21] codetools_0.2-18   tzdb_0.3.0         class_7.3-20       fansi_1.0.3        broom_0.8.0       
# [26] Rcpp_1.0.9         KernSmooth_2.23-20 scales_1.2.1       backports_1.4.1    classInt_0.4-7    
# [31] jsonlite_1.8.2     fs_1.5.2           hms_1.1.1          stringi_1.7.6      grid_4.1.3        
# [36] cli_3.4.1          tools_4.1.3        magrittr_2.0.3     proxy_0.4-27       crayon_1.5.1      
# [41] pkgconfig_2.0.3    ellipsis_0.3.2     xml2_1.3.3         reprex_2.0.1       lubridate_1.8.0   
# [46] assertthat_0.2.1   httr_1.4.3         rstudioapi_0.13    R6_2.5.1           units_0.8-0       
# [51] compiler_4.1.3