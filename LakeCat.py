# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:24:06 2016


#ctl = pd.read_csv('D:/Projects/LakeCat_scrap/ControlTable_LakeCat_RD.csv')

@author: Rdebbout
"""

import os
import sys
import arcpy
import numpy as np
import pandas as pd
import geopandas as gpd
#ctl = pd.read_csv(sys.argv[1])
ctl = pd.read_csv(r'F:/GitProjects/LakeCat/ControlTable_LakeCat.csv')
ctl = pd.read_csv('F:/GitProjects/NARS/Landscape Metrics/ControlTable_LakeCat_NLA17.csv')
arcpy.CheckOutExtension("spatial")
from datetime import datetime as dt
#sys.path.append(ctl.DirectoryLocations.values[3])  
sys.path.append('F:/GitProjects/LakeCat')
from arcpy.sa import ZonalStatisticsAsTable, TabulateArea
from LakeCat_functions import inputs, rpus, dbf2DF, getOnNetLakes2, chkColumnLength, PointInPoly, Accumulation

arcpy.env.cellSize = "30"
ingrid_dir = ctl.DirectoryLocations.values[0]

frame = ctl.DirectoryLocations.values[1]

out_dir = ctl.DirectoryLocations.values[2]
StreamCat = ctl.DirectoryLocations.values[4]
NHD_dir = ctl.DirectoryLocations.values[5]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    os.mkdir('%s/ZStats' % out_dir)
    
lls = [line for line in ctl.index if ctl.run[line] == 1]

for ll in lls:  # loop through each FullTableName in control table
    print 'running....%s' % ctl.LandscapeLayer[ll]
    accum_type = ctl.accum_type[ll]
    LLyr = '%s/%s' % (ingrid_dir, ctl.LandscapeLayer[ll])
    metric = ctl.MetricName[ll]
    name = ctl.FullTableName[ll]
    ddir = '%s/ZStats/%s' % (out_dir,name)
    if not os.path.exists(ddir) and accum_type != 'Point':
        os.mkdir(ddir)
    summaryfield = None
    if type(ctl.summaryfield[ll]) == str:
        summaryfield = ctl.summaryfield[ll].split(';')
    start = dt.now()

    if accum_type != 'Point':             
        csv = "%s/%s.csv" % (out_dir, name)
        stats = pd.DataFrame()
        for zone in rpus.keys():
            pre = '%s/NHDPlus%s/NHDPlus%s' % (NHD_dir, inputs[zone], zone)
            for rpu in rpus[zone]:
                if metric == 'Elev':
                    LLyr = '%s/NEDSnapshot/Ned%s/%s' % (pre, rpu, 
                                                        ctl.LandscapeLayer[ll])
                out = '{0}/ZStats/{1}/{1}_{2}.dbf'.format(out_dir, name, rpu)
                if not os.path.exists(out):
                    raster = '%s/rasters/wsheds/wtshds_%s.tif' % (frame, rpu)
                    if accum_type == 'Categorical':
                        TabulateArea(raster, "Value", LLyr, "Value", 
                                     out, "30")
                    if accum_type == 'Continuous':
                        ZonalStatisticsAsTable(raster, "Value", LLyr, 
                                               out, "DATA", "ALL")
                tbl = dbf2DF(out)
                tbl.rename(columns={'VALUE':'UID'},inplace=True)
                stats = pd.concat([stats, tbl])
        stats.to_csv(csv, index=False)
        
    if accum_type == 'Point':
        
        pct_full = pd.read_csv('%s/border/pct_full.csv' % frame)
        points = gpd.GeoDataFrame.from_file(LLyr)
        basins = '%s/shps/allBasins.shp' % (frame)
        stats = PointInPoly(points, basins, pct_full, 'UID', summaryfield)
         
    print 'ZonalStats Results Complete in : ' + str(dt.now() - start)
    
    if accum_type != 'Point':
        b = pd.DataFrame()
        for zone in rpus.keys():
            for rpu in rpus[zone]:
                b_ = dbf2DF('%s/rasters/wsheds/wtshds_%s.tif.vat.dbf' % (frame, 
                                                                         rpu))
                b_['BSNAREASQKM'] = (b_.COUNT * 900) * 1e-6
                b_ = b_[['VALUE', 'BSNAREASQKM', 'COUNT']]
                b_.columns = ['UID', 'AreaSqKm', 'COUNT']
                b = pd.concat([b,b_])
        
        
    if accum_type == 'Categorical':
        stats = chkColumnLength(stats,LLyr)
        cols = stats.columns.tolist()[1:]
        stats['AREA'] = stats[stats.columns.tolist()[1:]].sum(axis=1)
        stats = pd.merge(b, stats, how='left', on='UID')
        stats['PctFull'] = (((stats.AREA * 1e-6) / stats.AreaSqKm) * 100)
        stats = stats[['UID', 'AreaSqKm'] + cols + ['PctFull']]
        cols = stats.columns[1:]
        stats.columns = np.append('UID', 'Cat' + cols.values)
        stats = stats.fillna(0)
        
    if accum_type == 'Continuous':
        stats = pd.merge(b, stats, how='left', on='UID')
        stats['CatPctFull'] = ((stats.COUNT_y / stats.COUNT_x) * 100)
        if name == 'Elev':
            stats = stats[['UID','AreaSqKm','COUNT_x','SUM', 
                           'MAX', 'MIN', 'CatPctFull']]
            stats.columns = ['UID', 'CatAreaSqKm', 'CatCount', 'CatSum', 
                             'CatMax', 'CatMin', 'CatPctFull']
        else: 
            stats = stats[['UID','AreaSqKm','COUNT_x','SUM', 'CatPctFull']]
            stats.columns = ['UID', 'CatAreaSqKm', 'CatCount', 
                             'CatSum', 'CatPctFull']
        stats.CatPctFull = stats.CatPctFull.fillna(0)                
    start2 = dt.now()
    npy = '%s/LakeCat_npy' % frame
    accum = np.load('%s/bastards/accum.npz' % npy)
    up = Accumulation(stats, accum['comids'], 
                           accum['lengths'], 
                           accum['upstream'], 
                           'UpCat','UID')
    accum = np.load('%s/children/accum.npz' % npy)
    ws = Accumulation(stats, accum['comids'], 
                           accum['lengths'], 
                           accum['upstream'], 
                           'Ws','UID')    
    stats = pd.merge(stats, up, on='UID')
    stats = pd.merge(stats, ws, on='UID')
    cols = stats.columns[1:].tolist()
    # goto StreamCat to get On-Net-work lake results from assoc. COMIDs
    stats['inStreamCat'] = 0
    # Join UID to COMID for final deliverable..
    lks = dbf2DF('%s/off-network.dbf' % frame)[['COMID','UID']]
    off = pd.merge(lks,stats,on='UID',how='right')
    off.drop('UID',axis=1,inplace=True)
    on = getOnNetLakes2(name, StreamCat, 
                           '%s/joinTables' % frame , 
                           '%s/onNet_LakeCat.npz' % npy, 
                           NHD_dir)
    on['inStreamCat'] = 1
    print "Length of on_Net: " + str(len(on))
    tot = pd.concat([off, on])
    tot.to_csv('%s/%s.csv' % (out_dir, name), index=False)
    print 'Accumulation Results Complete in : ' + str(dt.now() - start2)
