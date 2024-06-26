# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:37:20 2024

@author: mweber
"""

import geopandas as gpd
import pandas as pd

wbdy = gpd.read_file("G:/NHDPlusV21/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb",
                     layer="NHDWaterbody",driver = 'FileGDB').to_crs(epsg="5070")
wbdy = wbdy.loc[wbdy["FTYPE"].isin(["LakePond", "Reservoir"])]

cnty = gpd.read_file('O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Resource/DEMOGRAPHIC/BOUNDARIES/NATIONAL/tl_2018_us_county.shp').to_crs(epsg="5070")

wbdy_cnty = gpd.sjoin(left_df=wbdy[['COMID','GNIS_ID','GNIS_NAME','REACHCODE','geometry']], right_df=cnty[['STATEFP','COUNTYFP','GEOID','geometry']], how="left", predicate="intersects")

hydrgn = gpd.read_file("G:/NHDPlusV21/NHDPlusGlobalData/BoundaryUnit.shp").to_crs(epsg="5070")
hydrgn = hydrgn.loc[hydrgn["UnitType"]=='VPU']

wbdy_hydrgn = gpd.sjoin(left_df=wbdy[['COMID','GNIS_ID','GNIS_NAME','REACHCODE','geometry']], right_df=hydrgn[['DrainageID','UnitID','UnitName','geometry']], how="left", predicate="intersects") 

wbdy_cnty = pd.DataFrame(wbdy_cnty.drop(columns='geometry'))
wbdy_hydrgn = pd.DataFrame(wbdy_hydrgn.drop(columns='geometry'))  

wbdy_cnty.to_csv('O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/LakeCat_County_State.csv', index=False)  
wbdy_hydrgn.to_csv('O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/LakeCat_HydroRegion.csv', index=False)  
