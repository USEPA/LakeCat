# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:25:07 2016

@author: mweber
"""
import geopandas as gpd
from geopandas.tools import overlay
import pandas as pd

states = gpd.GeoDataFrame.from_file('L:/Priv/CORFiles/Geospatial_Library/Data/RESOURCE/POLITICAL/BOUNDARIES/NATIONAL/TIGER_2010_State_Boundaries.shp')
states.plot()
states.head()
# only keep CONUS states
states.NAME10.unique()
states = states.loc[~states['NAME10'].isin(['Puerto Rico','Alaska','Hawaii'])]
states.plot()

eco9 = gpd.GeoDataFrame.from_file('L:/Priv/ARM Data/Omernik_Ecoregions_And_Aggregated_Ecoregions/Aggr_Ecoregions9_20121005.shp')
eco9.head()
lakes = gpd.GeoDataFrame.from_file('L:/Priv/CORFiles/Geospatial_Library/Data/Project/LakeCat/NHDPlusV21_Lakes.shp')
lakes.head()

# Run identity, get area sums, calculate dominant ecoregion for border lakes
lakes.crs
eco9.crs
lakes = lakes.to_crs(eco9.crs)
lakes_eco = overlay(lakes, eco9, how="identity")
lakes_eco.head()
lakes_eco['AREASKM'] = lakes_eco.area * 1e-6
area_sum = lakes_eco.groupby(by=['COMID'])['AREASKM'].sum().to_frame()
area_sum.head()
area_sum['COMID'] = area_sum.index
area_sum.rename(columns={'AREASKM':'SUMLAKEAREA_SKM'}, inplace=True)
lakes_sum = pd.merge(lakes_eco, area_sum, on='COMID', how='left')
lakes_sum.head()

lakes_sum.loc[lakes_sum['COMID'] == 120052975]
lakes_eco.loc[lakes_eco['COMID'] == 120052975]