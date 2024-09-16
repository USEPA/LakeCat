# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 09:29:27 2024

@author: mweber
"""

import pandas as pd
import numpy as np
import os
import sys
import glob

comid_wbcomid_lookups = r'O:\PRIV\CPHEA\PESD\COR\CORFILES\Geospatial_Library_Projects\LakeCat\LakeCat_Framework\joinTables' 
streamcat_results = 'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/AmaliaHandler/Allocation_and_Accumulation/Final_StreamCat'
lookup_files = glob.glob(os.path.join(comid_wbcomid_lookups, "*.csv"))     # advisable to use os.path.join as this makes concatenation OS independent
lakecat_results = 'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/AmaliaHandler/Allocation_and_Accumulation/Final_LakeCat_OnNetwork'

df_from_each_file = (pd.read_csv(f) for f in lookup_files)
lookup   = pd.concat(df_from_each_file, ignore_index=True)
lookup = lookup.rename(columns={'catCOMID': 'COMID', 'wbCOMID': 'WBCOMID'})

streamcat_files = glob.glob(os.path.join(streamcat_results, "*.csv"))

for sfile in streamcat_files:
    filename = sfile.rpartition('\\')[-1]
    filename = filename.replace('StreamCat', 'LakeCat')
    lake_file = pd.read_csv(sfile)
    lake_file = pd.merge(lake_file, lookup, how='inner', on='COMID')
    lake_file = lake_file.drop(['COMID', 'CatAreaSqKm_y'], axis=1)
    lake_file = lake_file.rename(columns={'WBCOMID': 'COMID','CatAreaSqKm_x': 'CatAreaSqKm'})
    cols = lake_file.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    lake_file = lake_file[cols] 
    lake_file[['COMID']] = lake_file[['COMID']].astype(int)
    lake_file.to_csv(lakecat_results + '/' + filename)
