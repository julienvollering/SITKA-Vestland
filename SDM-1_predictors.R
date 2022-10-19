# Preparing past climate data

library(tidyverse)
library(terra)
library(sf)

chelsafiles <- list.files(path = "./data/climate/past", full.names = TRUE)
past <- rast(chelsafiles)
plot(past$`CHELSA_bio1_1981-2010_V.2.1`)

psfiltrert <- st_read("./data/0025471-220831081235567/PSfiltrert.shp")
ps <- filter(psfiltrert, intHarris == 1 | retained == 1)
psbbox <- ps %>% 
  st_transform(crs = st_crs("ESRI:102017")) %>% 
  st_buffer(dist = 300e3) %>% 
  st_bbox()

thfiltrert <- st_read("./data/0026583-220831081235567/THfiltrert.shp")
th <- filter(thfiltrert, intPackee == 1 | retained == 1)
thbbox <- th %>% 
  st_transform(crs = st_crs("ESRI:102017")) %>% 
  st_buffer(dist = 300e3) %>% 
  st_bbox()

psbbox
thbbox
b <- rast(xmin=-4914757, xmax=-1272575, ymin=1498388, ymax=3444981,
          resolution = 1000, crs="ESRI:102017")
pastEA <- project(past, b, method = "bilinear")

plot(pastEA$`CHELSA_bio1_1981-2010_V.2.1`)
th %>% 
  st_transform(crs = st_crs("ESRI:102017")) %>% 
  st_geometry() %>% 
  plot(add=TRUE)

template <- pastEA$`CHELSA_bio1_1981-2010_V.2.1`
values(template) <- NA
tm03 <- st_read("./data/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp")
tm03EA <- tm03 %>% 
  filter(ISO3 %in% c("USA", "CAN")) %>% 
  st_transform(st_crs("ESRI:102017"))
landmask <- rasterize(tm03EA, template, field = 1)
landmask <- bind_rows(ps, th) %>% 
  st_transform(st_crs("ESRI:102017")) %>% 
  rasterize(landmask, field = 1, update = TRUE)
pastEAland <- mask(pastEA, landmask)
plot(pastEAland$`CHELSA_bio1_1981-2010_V.2.1`)

writeRaster(pastEAland, "./data/climate/pastEAland.tif")

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