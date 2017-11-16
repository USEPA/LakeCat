# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 15:29:55 2017

@author: RHill04
"""

import geopandas as gpd
import pandas as pd
import numpy as np

def watershed(focal_com, comids, lengths, up, basins):
    try:
        l = np.asscalar(lengths[np.in1d(comids, focal_com)])
        start = np.sum(lengths[:np.asscalar(np.where(np.in1d(comids, focal_com))[0])])
        uplist = up[start:start+l]
    except: 
        uplist = [focal_com]
    tmp_basin = basins.loc[basins['COMID'].isin(uplist)].dissolve(by='dummy')
    return tmp_basin

#Read in numpy files
wd = 'L:/Priv/CORFiles/Geospatial_Library/Data/Project/LakeCat/LakeCat_Framework/LakeCat_npy/'
off_comids = np.load(wd + '/chiildren/accum.npz')['comids']
on_comids = np.load(wd + '/onNet_LakeCat.npz')['vpus'].item()['17']
lengths = np.load(wd + '/accum.npz')['lengths']
up = np.load(wd + '/accum.npz')['upstream']
    #Read in lakes shapefile (398 lakes) and get their COMIDs
lk_dir = r'D:\Scratch\JohnIiames'
lakes = gpd.read_file(lk_dir + '/lakes_900m_GLB_noOverlap.shp')
lakes = np.array(lakes.COMID).astype(int)
    #Read in off-network basins shapefile (big)
basins = gpd.read_file(lk_dir + '/allBasins.shp')
basins['dummy'] = 1 #Add dummy column for dissolving
    #Use a LakeCat table to find off-network lakes
lakecat = pd.read_csv(lk_dir + '/BFI_Final.csv')
offnet = np.array(lakecat.loc[lakecat['inStreamCat'] == 0, 'COMID'])

lakesLC = lakes[np.in1d(lakes, offnet)]

for lake in lakesLC:
    print lake
    out_ws = watershed(lake, comids, lengths, up, basins)
    out_ws.crs = basins.crs
    out_ws.to_file(filename = lk_dir + '/out_shps/ws_' + str(lake) + '.shp', driver = 'ESRI Shapefile')
    #world.to_file(filename=temp_shp,driver='ESRI Shapefile',crs_wkt=prj)
    print '------------------------------------'


