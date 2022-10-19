# Projecting models

library(tidyverse)
library(terra)
library(sf)
library(maxnet)
library(MIAmaxent)

load(file = "SDM-6_trainingmodels.RData")

## Preparing Norway climate data

present <- rast(list.files(path = "./data/climate/present", full.names = TRUE))
future126 <- rast(list.files(path = "./data/climate/future", full.names = TRUE, pattern = "ssp126"))
future370 <- rast(list.files(path = "./data/climate/future", full.names = TRUE, pattern = "ssp370"))

### UTM

tm03 <- st_read("./data/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp")
tm03utm <- tm03 %>% 
  filter(ISO3 == "NOR") %>% 
  st_transform(st_crs("EPSG:25832"))

st_bbox(tm03utm)
b <- rast(xmin=250e3, xmax=1350e3, ymin=6410e3, ymax=8000e3,
          resolution = 1000, crs="EPSG:25832")
presentutm <- project(present, b, method = "bilinear")
future126utm <- project(future126, b, method = "bilinear")
future370utm <- project(future370, b, method = "bilinear")

templateutm <- presentutm$`CHELSA_bio1_2011-2040_gfdl-esm4_ssp126_V.2.1`
values(templateutm) <- NA
names(templateutm) <- "template"

st_layers("./data/Basisdata_0000_Norge_25833_N250Kartdata_GML/Basisdata_0000_Norge_25833_N250Hoyde_GML.gml")
N250 <- st_read("./data/Basisdata_0000_Norge_25833_N250Kartdata_GML/Basisdata_0000_Norge_25833_N250Hoyde_GML.gml",
                "Høydelag")
N250utm <- st_transform(N250, "EPSG:25832")

landmask <- rasterize(N250utm, templateutm, field = 1, touches = TRUE)
presentutm <- mask(presentutm, landmask)
future126utm <- mask(future126utm, landmask)
future370utm <- mask(future370utm, landmask)

presentutm <- as.data.frame(presentutm, cells = TRUE, na.rm = TRUE)
names(presentutm) <- make.names(names(presentutm), allow_ = FALSE)
presentutm <- presentutm %>% 
  rename_with(~ gsub(".2011.2040.gfdl.esm4.ssp126.V.2.1$", "", .x), ends_with(".2011.2040.gfdl.esm4.ssp126.V.2.1")) #compatible with names in projection
future126utm <- as.data.frame(future126utm, cells = TRUE, na.rm = TRUE)
names(future126utm) <- make.names(names(future126utm), allow_ = FALSE)
future126utm <- future126utm %>% 
  rename_with(~ gsub(".2071.2100.gfdl.esm4.ssp126.V.2.1$", "", .x), ends_with(".2071.2100.gfdl.esm4.ssp126.V.2.1")) #compatible with names in projection
future370utm <- as.data.frame(future370utm, cells = TRUE, na.rm = TRUE)
names(future370utm) <- make.names(names(future370utm), allow_ = FALSE)
future370utm <- future370utm %>% 
  rename_with(~ gsub(".2071.2100.gfdl.esm4.ssp370.V.2.1$", "", .x), ends_with(".2071.2100.gfdl.esm4.ssp370.V.2.1")) #compatible with names in projection

### Web Mercator

tm03wm <- tm03 %>% 
  filter(ISO3 == "NOR") %>% 
  st_transform(st_crs("EPSG:3857"))

st_bbox(tm03wm)
b <- rast(xmin=500e3, xmax=3475e3, ymin=7950e3, ymax=11475e3,
          resolution = 1000, crs="EPSG:3857")
presentwm <- project(present, b, method = "bilinear")
future126wm <- project(future126, b, method = "bilinear")
future370wm <- project(future370, b, method = "bilinear")

templatewm <- presentwm$`CHELSA_bio1_2011-2040_gfdl-esm4_ssp126_V.2.1`
values(templatewm) <- NA
names(templatewm) <- "template"

N250wm <- st_transform(N250, "EPSG:3857")

landmask <- rasterize(N250wm, templatewm, field = 1, touches = TRUE)
presentwm <- mask(presentwm, landmask)
future126wm <- mask(future126wm, landmask)
future370wm <- mask(future370wm, landmask)

presentwm <- as.data.frame(presentwm, cells = TRUE, na.rm = TRUE)
names(presentwm) <- make.names(names(presentwm), allow_ = FALSE)
presentwm <- presentwm %>% 
  rename_with(~ gsub(".2011.2040.gfdl.esm4.ssp126.V.2.1$", "", .x), ends_with(".2011.2040.gfdl.esm4.ssp126.V.2.1")) #compatible with names in projection
future126wm <- as.data.frame(future126wm, cells = TRUE, na.rm = TRUE)
names(future126wm) <- make.names(names(future126wm), allow_ = FALSE)
future126wm <- future126wm %>% 
  rename_with(~ gsub(".2071.2100.gfdl.esm4.ssp126.V.2.1$", "", .x), ends_with(".2071.2100.gfdl.esm4.ssp126.V.2.1")) #compatible with names in projection
future370wm <- as.data.frame(future370wm, cells = TRUE, na.rm = TRUE)
names(future370wm) <- make.names(names(future370wm), allow_ = FALSE)
future370wm <- future370wm %>% 
  rename_with(~ gsub(".2071.2100.gfdl.esm4.ssp370.V.2.1$", "", .x), ends_with(".2071.2100.gfdl.esm4.ssp370.V.2.1")) #compatible with names in projection

dflist <- list(presentutm = presentutm, 
               future126utm = future126utm, 
               future370utm = future370utm,
               presentwm = presentwm, 
               future126wm = future126wm, 
               future370wm = future370wm)
templatelist <- list(presentutm = templateutm, 
                     future126utm = templateutm, 
                     future370utm = templateutm,
                     presentwm = templatewm, 
                     future126wm = templatewm, 
                     future370wm = templatewm)

## Projecting

hyper
predict_maxnet_MIAmaxent <- function(model, newdata) {
  if ("maxnet" %in% class(model)) {
    out <- predict(model, newdata, clamp = FALSE, type = "exponential")
    return(out)
  } else if ("iwlr" %in% class(model)) {
    out <- MIAmaxent::projectModel(model, model$transformations, newdata, raw = TRUE, rescale = FALSE)
    return(out$output[,1])
  }
}

thresholdlist <- list(ps = thresholdpsensemble, th = thresholdthensemble)

#j <- "ps"; i <- names(dflist)[1]
for (j in c("ps", "th")) {
  for (i in names(dflist)) {
    projpreds <- hyper %>% 
      filter(sp == j) %>% 
      pull(model) %>% 
      map_dfc(~ predict_maxnet_MIAmaxent(model = .x, newdata = dflist[[i]]))
    scalers <- hyper %>% 
      filter(sp == j) %>% 
      pull(scaler)
    projpreds <- projpreds %>% 
      sweep(2, scalers, `/`) %>% 
      as_tibble()
    
    thresholds <- hyper %>% 
      filter(sp == j) %>% 
      pull(threshold) 
    projpredsbin <- map2_dfc(projpreds, thresholds, ~ .x >= .y)
    projpredsbin <- projpredsbin %>% 
      mutate(rowSums = rowSums(projpredsbin),
             dissent = pmin(abs(0-rowSums), abs(ncol(projpreds)-rowSums))) %>% 
      select(-rowSums)
    
    projpredsens <- projpreds %>% 
      mutate(ensemble = rowSums(projpreds))
    head(projpredsens)
    
    #### Transfer values to rasters
    r <- templatelist[[i]]
    r[dflist[[i]]$cell] <- projpredsens$ensemble
    names(r) <- "ensemble"
    r$ensemblebin <- r$ensemble >= thresholdlist[[j]]
    r$dissent <- NA
    r$dissent[dflist[[i]]$cell] <- projpredsbin$dissent
    
    writeRaster(r, paste0("./resultater/", j, "-", i, ".tif"), overwrite = TRUE)
  }
}

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
#   [1] MIAmaxent_1.2.0.9000 maxnet_0.1.4         sf_1.0-8             terra_1.6-17         forcats_0.5.1       
# [6] stringr_1.4.0        dplyr_1.0.10         purrr_0.3.4          readr_2.1.2          tidyr_1.2.0         
# [11] tibble_3.1.8         ggplot2_3.3.6        tidyverse_1.3.1     
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