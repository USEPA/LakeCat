# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 18:53:07 2021

@author: Rdebbout
"""

import os
import sys
import geopandas as gpd

from LakeCat_functions import NHDtblMerge, makeBasins, makeNParrays
from lake_cat_config import NHD_DIR as nhd
out = "framework"
if not os.path.exists(f"{out}/rasters"):
    if not os.path.exists(out):
        os.mkdir(out)
    os.mkdir(f"{out}/rasters")
    os.mkdir(f"{out}/rasters/lakes")
    os.mkdir(f"{out}/rasters/lakes/scratchArc")
    os.mkdir(f"{out}/rasters/wsheds")
    os.mkdir(f"{out}/shps")
    os.mkdir(f"{out}/joinTables")
    os.mkdir(f"{out}/LakeCat_npy")
    
# print(os.path.exists(f"{nhd}/NHDPlusGlobalData/BoundaryUnit.shp"))
NHDbounds = gpd.read_file(f"{nhd}/NHDPlusGlobalData/BoundaryUnit.shp").to_crs(epsg="5070")
NHDbounds.drop(['AreaSqKM','DrainageID','Shape_Area','Shape_Leng','UnitName'], axis=1, inplace=True)
# print(NHDbounds)

if not os.path.exists("framework/Lake_QA.csv"):
    NHDtblMerge(nhd, NHDbounds, "framework")
# makeBasins(nhd, NHDbounds, "framework")
makeNParrays("framework")

# old = "framework_first_run"
# new = "framework"
# for zone in inputs:
#     print(zone)
#     # net1 = gpd.read_file(f"{old}/off_net_{zone}.shp")
#     # net2 = gpd.read_file(f"{new}/off_net_{zone}.shp")
#     # # break
#     # assert len(net1) == len(net2)
#     # assert all(net1.columns == net2.columns)
#     # assert all(net1.sort_values("COMID").COMID == net1.sort_values("COMID").COMID)
#     j1 = pd.read_csv(f"{old}/joinTables/join_{zone}.csv")
#     j2 = pd.read_csv(f"{new}/joinTables/join_{zone}.csv")
#     assert len(j1) == len(j2)
#     assert all(j1.columns == j2.columns)
#     assert j1.equals(j2)
#     break