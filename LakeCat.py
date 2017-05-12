# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:24:06 2016

@author: Rdebbout
"""

import sys
import os
import numpy as np
import pandas as pd
import geopandas as gpd
ctl = pd.read_csv(sys.argv[1])
#ctl = pd.read_csv('D:/Projects/LakeCatOutput/ControlTable_LakeCat_RD.csv')
from datetime import datetime as dt
import arcpy
arcpy.CheckOutExtension("spatial")
from arcpy.sa import ZonalStatisticsAsTable, TabulateArea
sys.path.append(ctl.DirectoryLocations.values[5])  #  sys.path.append('D:/Projects/lakesAnalysis/Scripts')
from LakeCat_functions import dbf2DF, createAccumTable, getOnNetLakes2, chkColumnLength, PointInPoly

arcpy.env.cellSize = "30"
ingrid_dir = ctl.DirectoryLocations.values[0]
inZoneData = ctl.DirectoryLocations.values[1]
out_dir = ctl.DirectoryLocations.values[2]
npy_files = ctl.DirectoryLocations.values[3]
pct_full_file = ctl.DirectoryLocations.values[4]
StreamCat = ctl.DirectoryLocations.values[6]
LakeComs = ctl.DirectoryLocations.values[7]
NHD_dir = ctl.DirectoryLocations.values[8]
onNet_npy = ctl.DirectoryLocations.values[9]
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    os.mkdir('%s/ZStats' % out_dir)
for line in range(len(ctl.values)):  # loop through each FullTableName in control table
    if ctl.run[line] == 1:
        break
        print 'running....%s' % ctl.LandscapeLayer[line]
        accum_type = ctl.accum_type[line]
        LandscapeLayer = '%s/%s' % (ingrid_dir, ctl.LandscapeLayer[line])
        metric = ctl.MetricName[line]
        if not os.path.exists('%s/ZStats/%s' % (out_dir,metric)):
            os.mkdir('%s/ZStats/%s' % (out_dir,metric))
        name = ctl.FullTableName[line]
        summaryfield = None
        if type(ctl.summaryfield[line]) == str:
            summaryfield = ctl.summaryfield[line].split(';')
        start = dt.now()

        if accum_type != 'Point':            
            inputs = np.load('%s/StreamCat_npy/zoneInputs.npy' % NHD_dir).item()
            rpus = np.load('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir).item()
            if metric == 'Elev':
                outTable = "%s/%s_temp.csv" % (out_dir, name)  # made and later erased after read in
            else:   
                outTable = "%s/%s.csv" % (out_dir, name)
            count = 0
            for zone in rpus.keys():
                for rpu in rpus[zone]:
                    if metric == 'Elev':
                        LandscapeLayer = '%s/NHDPlus%s/NHDPlus%s/NEDSnapshot/Ned%s/%s' % (NHD_dir, inputs[zone], zone, rpu, ctl.LandscapeLayer[line])
                    out = '{1}/ZStats/{1}/{1}_{2}.dbf'.format(out_dir, metric, rpu)
                    if not os.path.exists(out):
                        raster = '%s/wtshds_%s.tif' % (inZoneData, rpu)
                        if metric == 'Elev':
                            ZonalStatisticsAsTable(raster, "Value", LandscapeLayer, out, "DATA", "ALL")
                        else:
                            if accum_type == 'Categorical':
                                TabulateArea(raster, "VALUE", LandscapeLayer, "Value", out, "30")
                            if accum_type == 'Continuous':
                                ZonalStatisticsAsTable(raster, "Value", LandscapeLayer, out, "DATA", "SUM")
                    tbl = dbf2DF(out)
                    if count == 0:
                        final = tbl
                    if count > 0:
                        final = pd.concat([final, tbl])
                    count += 1
            final.to_csv(outTable, index=False)
            
#        if accum_type == 'Continuous' and metric != 'Elev':
#            
#            if LandscapeLayer.count('.tif') or LandscapeLayer.count('.img'):
#                outTable ="%s/zstats_%s.dbf" % (out_dir,LandscapeLayer.split("/")[-1].split(".")[0])
#            else:
#                outTable ="%s/zstats_%s.dbf" % (out_dir,LandscapeLayer.split("/")[-1])
#            if not os.path.exists(outTable):
#                print 'running ZstatsTable'
#                ZonalStatisticsAsTable(inZoneData + '/ALL_BASINS.tif', "Value", LandscapeLayer, outTable, "DATA", "ALL")
#                
        if accum_type == 'Point':
            
            pct_full = pd.read_csv(pct_full_file)
            points = gpd.GeoDataFrame.from_file(LandscapeLayer)
            basins = '%s/all_basins.shp' % (inZoneData)
            join = PointInPoly(points, basins, pct_full, summaryfield)
            join.rename(columns={'GRIDCODE': 'COMID'}, inplace=True)
## CAT STATS            
        print 'ZonalStats Results Complete in : ' + str(dt.now() - start)
        #if not os.path.exists('%s/%s.csv' % (out_dir, metric)):

        start2 = dt.now()
        
        if accum_type != 'Point':
            basin_tbl = pd.DataFrame()
            for zone in rpus.keys():
                for rpu in rpus[zone]:
                    b_tbl = dbf2DF('%s/wtshds_%s.tif.vat.dbf' % (inZoneData, rpu))
                    b_tbl['BSNAREASQKM'] = (b_tbl.COUNT * 900) * 1e-6
                    b_tbl = b_tbl[['VALUE','BSNAREASQKM','COUNT']]
                    b_tbl.columns = ['COMID', 'AreaSqKm', 'COUNT']
                    basin_tbl = pd.concat([basin_tbl,b_tbl])
            
            
tit = dbf2DF('D:/Projects/LakeCat/offnet_LakeCat_framework/off-network.dbf')
tot = []
for it in tit.COMID.values:
    if not it in basin_tbl.COMID.values:
        tot.append(it)


            
        if accum_type == 'Categorical':
            
            tbl = pd.read_csv(outTable)
            tbl = chkColumnLength(tbl,LandscapeLayer)
            cols = tbl.columns.tolist()[1:]
            tbl['AREA'] = tbl[tbl.columns.tolist()[1:]].sum(axis=1)
            join = pd.merge(b_tbl, tbl, how='left', left_on='COMID', right_on='VALUE')
            join['PctFull'] = (((join.AREA * 1e-6) / join.AreaSqKm) * 100)
            join = join[['COMID', 'AreaSqKm'] + cols + ['PctFull']]
            cols = join.columns[1:]
            join.columns = np.append('COMID', 'Cat' + cols.values)
            join = join.fillna(0)
            
        if accum_type == 'Continuous':
            if name == 'Elev': # done by rpu and written to csv already
                tbl = pd.read_csv(outTable)[['VALUE', 'AREA', 'COUNT','SUM']]
                os.remove(outTable)
            else:
                tbl = pd.read_csv(outTable)[['VALUE', 'AREA', 'COUNT','SUM']]
            join = pd.merge(basin_tbl, tbl, how='left', left_on='COMID', right_on='VALUE')
            join['CatPctFull'] = ((join.COUNT_x / join.COUNT_y) * 100)
            join = join[['COMID','AreaSqKm','COUNT_x','SUM', 'CatPctFull']]
            join.columns = ['COMID', 'CatAreaSqKm', 'CatCount', 'CatSum', 'CatPctFull']
            join.CatPctFull = join.CatPctFull.fillna(0)                

        accum = createAccumTable(join, npy_files, 'UpCat')
        accum2 = createAccumTable(join, npy_files, 'Ws')
        upFinal = pd.merge(join, accum, on='COMID')
        final = pd.merge(upFinal,accum2, on='COMID')
        cols = final.columns[1:].tolist()
        
        # The next 5 lines take the COMID that is in the final DF and translates it to Waterbody COMID
        # til this point we had been using the GRIDCODE under the column title COMID to use in the 
        # Accumulation function without having to re-write it
        
        shptbl = dbf2DF('D:/Projects/lakesAnalysis/new/Isolated_AlbersOut1_Arc.dbf')[['COMID','GRIDCODE']]
        final = pd.merge(final, shptbl, left_on='COMID', right_on='GRIDCODE', how='left')
        final = final.drop(['COMID_x','GRIDCODE'],axis=1)
        final = final[['COMID_y'] + cols] 
        final.columns = ['COMID'] + cols
        final['inStreamCat'] = 1
        print "Length of final: " + str(len(final))
        
        # goto StreamCat to get On-Net-work lake results from assoc. COMIDs
        
        onNet = getOnNetLakes2(name, StreamCat, LakeComs, onNet_npy)
        onNet['inStreamCat'] = 0
        print "Length of on_Net: " + str(len(onNet))
        out = pd.concat([final, onNet])
        out.to_csv('%s/%s.csv' % (out_dir, name), index=False)
        print 'Accumulation Results Complete in : ' + str(dt.now() - start2)
