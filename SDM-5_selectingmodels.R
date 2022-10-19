# Selecting hyperparameters for final models

library(tidyverse)
library(terra)
library(sf)

tuningmaxnet <- read_csv("SDM-4_tuningmaxnet.csv")
filter(tuningmaxnet, is.na(auc)) %>% 
  nrow()
# tuningmaxnet %>% 
#   select(-auc) %>% 
#   distinct() %>% 
#   nrow()
tuningMIAmaxent <- read_csv("SDM-4_tuningMIAmaxent.csv")
filter(tuningMIAmaxent, is.na(auc)) %>% 
  nrow()
# tuningMIAmaxent %>% 
#   select(-auc) %>% 
#   distinct() %>% 
#   nrow()

## Hyperparameters for P.s.

ggplot(data = filter(tuningmaxnet, sp == "ps")) +
  geom_jitter(mapping = aes(x=filtering, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "ps")) +
  geom_jitter(mapping = aes(x=regularization, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "ps")) +
  geom_jitter(mapping = aes(x=fc, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "ps")) +
  geom_jitter(mapping = aes(x=testfold, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "ps", testfold == 5)) +
  geom_jitter(mapping = aes(x=filtering, y = auc))

tuningmaxnetcv <- group_by(tuningmaxnet, sp, filtering, regularization, fc) %>% 
  summarise(auc = mean(auc))

ggplot(data = filter(tuningMIAmaxent, sp == "ps")) +
  geom_jitter(mapping = aes(x=filtering, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "ps")) +
  geom_jitter(mapping = aes(x=transformtypes, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "ps")) +
  geom_jitter(mapping = aes(x=interact, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "ps")) +
  geom_jitter(mapping = aes(x=factor(a), y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "ps")) +
  geom_jitter(mapping = aes(x=testfold, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "ps", testfold == 5)) +
  geom_jitter(mapping = aes(x=filtering, y = auc))

tuningMIAmaxentcv <- group_by(tuningMIAmaxent, sp, filtering, transformtypes, interact, a) %>% 
  summarise(auc = mean(auc))

tuningcv <- bind_rows(tuningmaxnetcv, tuningMIAmaxentcv)
cutoff <- tuningcv %>% 
  filter(sp == "ps") %>% 
  pull(auc) %>% 
  max(na.rm = TRUE) %>% 
  `*`(0.99)
hyperps <- tuningcv %>% 
  filter(sp == "ps") %>% 
  filter(auc >= cutoff) %>% 
  ungroup() %>% 
  arrange(desc(auc), desc(transformtypes), desc(interact), a) %>% 
  distinct(regularization, auc, .keep_all = TRUE)

## Hyperparameters for T.h.

### maxnet

ggplot(data = filter(tuningmaxnet, sp == "th")) +
  geom_jitter(mapping = aes(x=filtering, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "th")) +
  geom_jitter(mapping = aes(x=regularization, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "th")) +
  geom_jitter(mapping = aes(x=fc, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "th")) +
  geom_jitter(mapping = aes(x=testfold, y = auc))
ggplot(data = filter(tuningmaxnet, sp == "th", testfold == 5)) +
  geom_jitter(mapping = aes(x=filtering, y = auc))

ggplot(data = filter(tuningMIAmaxent, sp == "th")) +
  geom_jitter(mapping = aes(x=filtering, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "th")) +
  geom_jitter(mapping = aes(x=transformtypes, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "th")) +
  geom_jitter(mapping = aes(x=interact, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "th")) +
  geom_jitter(mapping = aes(x=factor(a), y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "th")) +
  geom_jitter(mapping = aes(x=testfold, y = auc))
ggplot(data = filter(tuningMIAmaxent, sp == "th", testfold == 5)) +
  geom_jitter(mapping = aes(x=filtering, y = auc))

cutoff <- tuningcv %>% 
  filter(sp == "th") %>% 
  pull(auc) %>% 
  max(na.rm = TRUE) %>% 
  `*`(0.99)
hyperth <- tuningcv %>% 
  filter(sp == "th") %>% 
  filter(auc >= cutoff) %>% 
  ungroup() %>% 
  arrange(desc(auc), desc(transformtypes), desc(interact), a) %>% 
  distinct(regularization, auc, .keep_all = TRUE)

hyper <- bind_rows(hyperps, hyperth)

## Make data frame for training final models 

pastEAland <- rast("./climate/pastEAland.tif")
names(pastEAland) <- make.names(names(pastEAland), allow_ = FALSE)
load("SDM-2_CVfolds.RData")

### Full data frame for P.s.

saPS <- mask(pastEAland, vect(st_union(st_buffer(psEA, dist = 200e3))))
PS <- as.data.frame(saPS, xy = TRUE, cells = TRUE)

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

presmesh1cells <- distinct(psEA, mesh1nr, .keep_all = TRUE) %>% 
    pull(cellnr)
presmesh2cells <- distinct(psEA, mesh2nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh4cells <- distinct(psEA, mesh4nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh8cells <- distinct(psEA, mesh8nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh16cells <- distinct(psEA, mesh16nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh32cells <- distinct(psEA, mesh32nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh64cells <- distinct(psEA, mesh64nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh128cells <- distinct(psEA, mesh128nr, .keep_all = TRUE) %>% 
  pull(cellnr)

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

### Full data frame for T.h.

saTH <- mask(pastEAland, vect(st_union(st_buffer(thEA, dist = 200e3))))
TH <- as.data.frame(saTH, xy = TRUE, cells = TRUE)

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

presmesh1cells <- distinct(thEA, mesh1nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh2cells <- distinct(thEA, mesh2nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh4cells <- distinct(thEA, mesh4nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh8cells <- distinct(thEA, mesh8nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh16cells <- distinct(thEA, mesh16nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh32cells <- distinct(thEA, mesh32nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh64cells <- distinct(thEA, mesh64nr, .keep_all = TRUE) %>% 
  pull(cellnr)
presmesh128cells <- distinct(thEA, mesh128nr, .keep_all = TRUE) %>% 
  pull(cellnr)

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

### Combine species

df <- bind_rows(mutate(PS, species = 'ps'), mutate(TH, species = 'th'))

## Save data

rm(list=setdiff(ls(), c("hyper", "df")))
save.image("SDM-5_selectingmodels.RData")

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
# [11] withr_2.5.0        DBI_1.1.3          bit64_4.0.5        dbplyr_2.1.1       modelr_0.1.8      
# [16] readxl_1.4.0       lifecycle_1.0.3    munsell_0.5.0      gtable_0.3.1       cellranger_1.1.0  
# [21] rvest_1.0.2        codetools_0.2-18   tzdb_0.3.0         parallel_4.1.3     class_7.3-20      
# [26] fansi_1.0.3        broom_0.8.0        Rcpp_1.0.9         KernSmooth_2.23-20 scales_1.2.1      
# [31] backports_1.4.1    classInt_0.4-7     vroom_1.5.7        jsonlite_1.8.2     bit_4.0.4         
# [36] fs_1.5.2           hms_1.1.1          stringi_1.7.6      grid_4.1.3         cli_3.4.1         
# [41] tools_4.1.3        magrittr_2.0.3     proxy_0.4-27       crayon_1.5.1       pkgconfig_2.0.3   
# [46] ellipsis_0.3.2     xml2_1.3.3         reprex_2.0.1       lubridate_1.8.0    assertthat_0.2.1  
# [51] httr_1.4.3         rstudioapi_0.13    R6_2.5.1           units_0.8-0        compiler_4.1.3  