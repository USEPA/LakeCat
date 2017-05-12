# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:43:33 2016

@author: Rdebbout
"""

##############################################################################
# Find lake basins that are larger than CatAreaSqKM by using COUNT in reg_wtshds.dbf
# use checks list to determine if any raster counts are bigger than the catgrid counts

import sys
import os
import pandas as pd
import numpy as np
sys.path.append('D:/Projects/StreamCat') # change to your streamcat directory
from StreamCat_functions import makeRPUdict, dbf2DF

NHD_dir ='D:/NHDPlusV21'
if not os.path.exists('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir):
    rpuDict = makeVPUdict(NHD_dir)
else:
    rpuDict = np.load('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir).item()

#checkDict = {}
count = 0 
for zone in rpuDict:  
    print zone
    cats = pd.read_csv('D:/Projects/lakesAnalysis/710_basin_VIS/LakeCatJoin/lakeJoin%s.csv' % zone)
    cats.GRIDCODE = cats.GRIDCODE.values.astype(int)
    for rpu in rpuDict[zone]:
        print rpu        
        tbl = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_1/reg%s_wtshds.tif.vat.dbf' % rpu)
        both = pd.merge(cats, tbl, how='inner', left_on='GRIDCODE', right_on='VALUE')
        print 'Length both: ' + str(len(both))
        both['AreaSqKM_basin'] = (both.COUNT * 900) * 1e-6
        probs = both.ix[both.AreaSqKM < both.AreaSqKM_basin]
        probs = probs.drop(probs.ix[abs(probs.AreaSqKM - probs.AreaSqKM_basin) < .01].index, axis=0).copy()
        probs['VPU'] = zone
        probs['RPU'] = rpu
        probs.head()
        print 'Length probs: ' + str(len(probs))
        if count == 0:
            final = probs
        if count > 0:
            final = pd.concat([final,probs])
        count += 1
        
final = final.drop(['VALUE'], axis=1)
final = final[['COMID','GRIDCODE', 'AreaSqKM', 'FEATUREID', 'COUNT', 'AreaSqKM_basin', 'VPU', 'RPU']]
final.columns = ['WBCOMID', 'GRIDCODE', 'AreaSqKM', 'CATCOMID', 'ras_COUNT', 'AreaSqKM_basin', 'VPU', 'RPU']    
final.to_csv('C:/Users/Rdebbout/Desktop/ProblemLakeBasins_run6.csv',index=False)
final = final.drop(['CALC_FIELD'], axis=1)        



        if len(tbl) != len(cats):
            outs = []
            if len(tbl) > len(cats):
                for x in tbl.VALUE.values:
                    if x not in cats.COMID.values:
                        outs.append(x)
            if len(tbl) < len(cats):
                for x in cats.COMID.values:
                    if x not in tbl.VALUE.values:
                        outs.append(x)
            checkDict[rpu] = outs
##############################################################################
count = 0
for zone in rpuDict:  
    #print zone
    for rpu in rpuDict[zone]:
        #print rpu        
        tbl = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_1/reg%s_wtshds.tif.vat.dbf' % rpu)
        
        if not len([x for x in tbl.VALUE.values if x > 255793]) > 0:
            print rpu