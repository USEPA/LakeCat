# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 12:57:40 2016

@author: Rdebbout
"""

import os
import sys
import arcpy
import numpy as np
import pandas as pd
arcpy.CheckOutExtension("spatial")
from arcpy.sa import TabulateArea
sys.path.append('D:/Projects/StreamCat')
from StreamCat_functions import dbf2DF
from datetime import datetime as dt
arcpy.env.cellSize = "30"
NHD_dir = 'D:/NHDPlusV21'  
catRas = 'D:/Projects/lakesAnalysis/landscapeLayers/nlcd2006.tif'
inputs = np.load('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir).item()
print 'beginning....'
start = dt.now()
count = 0
for zone in inputs:
    for rpu in inputs[zone]:
        print rpu
        out = 'D:/Projects/LakeCat/Allocation_Accumulation/TabRpus/zonal_%srpu.dbf' % rpu
        if not os.path.exists(out):        
            raster = 'D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_1/reg%s_wtshds.tif' % rpu
            TabulateArea(raster, "VALUE", catRas, "Value", out, "30")
        if not os.path.exists('D:/Projects/LakeCat/Allocation_Accumulation/TabRpus/zonal_ALL3.csv'):
            tbl = dbf2DF(out)
            if count == 0:
                final = tbl
            if count > 0:
                final = pd.concat([final, tbl])
            count += 1
# need to merge with table of complete off-net so we don't lose records fromm 
# zstats where there is no coverage
allBasins = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_1/ALL_BASINS.tif.vat.dbf')[['VALUE']]
final = pd.merge(allBasins, final, on='VALUE', how='left')
final.to_csv('D:/Projects/LakeCat/Allocation_Accumulation/TabRpus/zonal_ALL.csv', index=False)
print 'Total time: %s' % str(dt.now() - start)


#            if accum_type == 'Categorical':
#                print 'running tabArea'
#                arcpy.AddField_management(inZoneData, "CALC_FIELD", "LONG", "", "")
#                arcpy.CalculateField_management(inZoneData, "CALC_FIELD", '!OID!', "PYTHON")
#                TabulateArea(inZoneData, "CALC_FIELD", LandscapeLayer, "Value", outTable, "30")
#                arcpy.JoinField_management(outTable,"CALC_FIELD",inZoneData,"CALC_FIELD",fields="Value")
#                arcpy.DeleteField_management(inZoneData, "CALC_FIELD")
#                arcpy.DeleteField_management(outTable, "CALC_FIELD")
#            if accum_type == 'Continuous':
#                print 'running ZstatsTable'
#                ZonalStatisticsAsTable(inZoneData, "Value", LandscapeLayer, outTable, "DATA", "ALL")