# Cartography with results

library(tidyverse)
library(terra)
library(sf)

# P.s. 

## EPSG:3857

files <- list.files("./resultater/", pattern = "wm", full.names = TRUE)[1:3]

r <- map(files, rast)
scaler <- r[[3]]$ensemble %>% 
  values() %>% 
  mean(na.rm = TRUE)

carto <- r[[3]]$ensemble/scaler
names(carto) <- "ssp126_2011.2040"
carto$ssp126_2011.2040_bin <- r[[3]]$ensemblebin
carto$ssp126_2071.2100 <- r[[1]]$ensemble/scaler
carto$ssp126_2071.2100_bin <- r[[1]]$ensemblebin
carto$ssp370_2071.2100 <- r[[2]]$ensemble/scaler
carto$ssp370_2071.2100_bin <- r[[2]]$ensemblebin

plot(r[[3]]$dissent)
m <- matrix(c(0, 1, 
              1, 1, 
              2, 2, 
              3, 3), ncol = 2, byrow = TRUE)
carto$ssp126_2011.2040_uncertainty <- classify(r[[3]]$dissent, m)
writeRaster(carto, "./resultater/Picea_sitchensis-epsg3857.tif", overwrite = TRUE)

## EPSG:25832

files <- list.files("./resultater/", pattern = "utm", full.names = TRUE)[1:3]

r <- map(files, rast)
scaler <- r[[3]]$ensemble %>% 
  values() %>% 
  mean(na.rm = TRUE)

carto <- log10(r[[3]]$ensemble/scaler)
names(carto) <- "ssp126_2011.2040"
carto$ssp126_2011.2040_bin <- r[[3]]$ensemblebin
carto$ssp126_2071.2100 <- log10(r[[1]]$ensemble/scaler)
carto$ssp126_2071.2100_bin <- r[[1]]$ensemblebin
carto$ssp370_2071.2100 <- log10(r[[2]]$ensemble/scaler)
carto$ssp370_2071.2100_bin <- r[[2]]$ensemblebin

plot(r[[3]]$dissent)
carto$ssp126_2011.2040_uncertainty <- classify(r[[3]]$dissent, m)
range(c(values(carto$ssp126_2011.2040), 
        values(carto$ssp126_2071.2100), 
        values(carto$ssp370_2071.2100)), na.rm = TRUE)
writeRaster(carto, "./resultater/Picea_sitchensis-epsg25832.tif", overwrite = TRUE)

r2 <- r[[3]]$ensemblebin
values(r2) <- values(r2) #force memory storage for a.polygons values to work
r2 <- r2 %>% 
  as.polygons(values = TRUE) %>% 
  st_as_sf() %>% 
  st_cast(to = "MULTILINESTRING")
s1 <- st_intersection(filter(r2, ensemblebin == 0), filter(r2, ensemblebin == 1)) %>% 
  st_collection_extract(type = "LINESTRING") %>% 
  mutate(scene = "ssp126_2011.2040") %>% 
  select(scene)

r2 <- r[[1]]$ensemblebin
values(r2) <- values(r2) #force memory storage for as.polygons values to work
r2 <- r2 %>% 
  as.polygons(values = TRUE) %>% 
  st_as_sf() %>% 
  st_cast(to = "MULTILINESTRING")
s2 <- st_intersection(filter(r2, ensemblebin == 0), filter(r2, ensemblebin == 1)) %>% 
  st_collection_extract(type = "LINESTRING") %>% 
  mutate(scene = "ssp126_2071.2100") %>% 
  select(scene)

r2 <- r[[2]]$ensemblebin
values(r2) <- values(r2) #force memory storage for as.polygons values to work
r2 <- r2 %>% 
  as.polygons(values = TRUE) %>% 
  st_as_sf() %>% 
  st_cast(to = "MULTILINESTRING")
s3 <- st_intersection(filter(r2, ensemblebin == 0), filter(r2, ensemblebin == 1)) %>% 
  st_collection_extract(type = "LINESTRING") %>% 
  mutate(scene = "ssp370_2071.2100") %>% 
  select(scene)

bind_rows(s1, s2, s3) %>% 
  st_write("./resultater/Picea_sitchensis-binary-epsg25832.shp")

# T.h. 

## EPSG:3857

files <- list.files("./resultater/", pattern = "wm", full.names = TRUE)[4:6]

r <- map(files, rast)
scaler <- r[[3]]$ensemble %>% 
  values() %>% 
  mean(na.rm = TRUE)

carto <- r[[3]]$ensemble/scaler
names(carto) <- "ssp126_2011.2040"
carto$ssp126_2011.2040_bin <- r[[3]]$ensemblebin
carto$ssp126_2071.2100 <- r[[1]]$ensemble/scaler
carto$ssp126_2071.2100_bin <- r[[1]]$ensemblebin
carto$ssp370_2071.2100 <- r[[2]]$ensemble/scaler
carto$ssp370_2071.2100_bin <- r[[2]]$ensemblebin

plot(r[[3]]$dissent)
m <- matrix(c(0, 1, 
              1, 1, 
              2, 1,
              3, 2,
              4, 2,
              5, 3,
              6, 3), ncol = 2, byrow = TRUE)
carto$ssp126_2011.2040_uncertainty <- classify(r[[3]]$dissent, m)
writeRaster(carto, "./resultater/Tsuga_heterophylla-epsg3857.tif", overwrite = TRUE)

## EPSG:25832

files <- list.files("./resultater/", pattern = "utm", full.names = TRUE)[4:6]

r <- map(files, rast)
scaler <- r[[3]]$ensemble %>% 
  values() %>% 
  mean(na.rm = TRUE)

carto <- log10(r[[3]]$ensemble/scaler)
names(carto) <- "ssp126_2011.2040"
carto$ssp126_2011.2040_bin <- r[[3]]$ensemblebin
carto$ssp126_2071.2100 <- log10(r[[1]]$ensemble/scaler)
carto$ssp126_2071.2100_bin <- r[[1]]$ensemblebin
carto$ssp370_2071.2100 <- log10(r[[2]]$ensemble/scaler)
carto$ssp370_2071.2100_bin <- r[[2]]$ensemblebin

plot(r[[3]]$dissent)
carto$ssp126_2011.2040_uncertainty <- classify(r[[3]]$dissent, m)
range(c(values(carto$ssp126_2011.2040), 
        values(carto$ssp126_2071.2100), 
        values(carto$ssp370_2071.2100)), na.rm = TRUE)
writeRaster(carto, "./resultater/Tsuga_heterophylla-epsg25832.tif", overwrite = TRUE)

r2 <- r[[3]]$ensemblebin
values(r2) <- values(r2) #force memory storage for a.polygons values to work
r2 <- r2 %>% 
  as.polygons(values = TRUE) %>% 
  st_as_sf() %>% 
  st_cast(to = "MULTILINESTRING")
s1 <- st_intersection(filter(r2, ensemblebin == 0), filter(r2, ensemblebin == 1)) %>% 
  st_collection_extract(type = "LINESTRING") %>% 
  mutate(scene = "ssp126_2011.2040") %>% 
  select(scene)

r2 <- r[[1]]$ensemblebin
values(r2) <- values(r2) #force memory storage for as.polygons values to work
r2 <- r2 %>% 
  as.polygons(values = TRUE) %>% 
  st_as_sf() %>% 
  st_cast(to = "MULTILINESTRING")
s2 <- st_intersection(filter(r2, ensemblebin == 0), filter(r2, ensemblebin == 1)) %>% 
  st_collection_extract(type = "LINESTRING") %>% 
  mutate(scene = "ssp126_2071.2100") %>% 
  select(scene)

r2 <- r[[2]]$ensemblebin
values(r2) <- values(r2) #force memory storage for as.polygons values to work
r2 <- r2 %>% 
  as.polygons(values = TRUE) %>% 
  st_as_sf() %>% 
  st_cast(to = "MULTILINESTRING")
s3 <- st_intersection(filter(r2, ensemblebin == 0), filter(r2, ensemblebin == 1)) %>% 
  st_collection_extract(type = "LINESTRING") %>% 
  mutate(scene = "ssp370_2071.2100") %>% 
  select(scene)

bind_rows(s1, s2, s3) %>% 
  st_write("./resultater/Tsuga_heterophylla-binary-epsg25832.shp")


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