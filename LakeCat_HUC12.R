library(sf)
library(dplyr)
library(readr)
gdb <- 'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Resource/PHYSICAL/HYDROLOGY/NHDPlusV21/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb'
st_layers(gdb)
lakes <- read_sf(dsn=gdb, layer='NHDWaterbody')
lakes <- lakes %>% 
  dplyr::filter(FTYPE %in% c('LakePond','Reservoir')) %>% 
  dplyr::select(COMID,GNIS_NAME,REACHCODE,FTYPE,ONOFFNET,LakeArea) %>% 
  st_centroid()

huc12 <- read_sf(dsn=gdb, layer='HUC12')
lakes <- st_transform(lakes, 5070)
huc12 <- st_transform(huc12, 5070)
lakes <- st_join(lakes, huc12)
any(is.na(lakes$HUC_12))
missing <- lakes[is.na(lakes$HUC_12),]
missing <- st_join(missing, huc12,join=st_nearest_feature)
lakes <- lakes[!is.na(lakes$HUC_12),]
lakes <- lakes %>% 
  dplyr::select(COMID,GNIS_NAME,REACHCODE,FTYPE,ONOFFNET,
                LakeArea,HUC_12,HU_12_NAME)
missing <- missing %>% 
  dplyr::select(COMID,GNIS_NAME,REACHCODE,FTYPE,ONOFFNET,
                LakeArea,HUC_12=HUC_12.y,HU_12_NAME=HU_12_NAME.y)
lakes <- rbind(lakes, missing)
write_csv(st_drop_geometry(lakes),'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/LakeCat_HUC12_lookup.csv')
