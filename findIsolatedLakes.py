# -*- coding: utf-8 -*-
"""
Created on Fri May 06 10:30:14 2016

@author: Rdebbout
"""

import os
import sys
import arcpy
import numpy as np
import pandas as pd
import rasterio as rs
import geopandas as gpd
from rasterio import features
from arcpy.sa import Watershed
from geopandas.tools import sjoin
from datetime import datetime as dt
sys.path.append('D:/Projects/LakeCat')
from LakeCat_functions import dbf2DF, NHD_Dict, DF2dbf, makeRat, purge


arcpy.CheckOutExtension("spatial")
arcpy.env.outputCoordinateSystem = 'PROJCS["NAD_1983_Contiguous_USA_Albers",'\
'GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",'\
'SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],'\
'UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],'\
'PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],'\
'PARAMETER["central_meridian",-96.0],PARAMETER["standard_parallel_1",29.5],'\
'PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_origin",23.0],'\
'UNIT["Meter",1.0]]'

NHD_dir = 'D:/NHDPlusV21'                    
inputs = NHD_Dict(NHD_dir)  # dictionaries to iterate thru NHD folder structure
rasterUnits = NHD_Dict(NHD_dir, unit='RPU')
outdir = 'D:/Projects/LakeCat/play'

# Shapefiles to store on and off-network lakes
outOn = '{}/NetworkLakes.shp'.format(outdir)
outOff = '{}/IsolatedLakes.shp'.format(outdir)
# create directories to hold intermediate data
if not os.path.exists(outdir):
    os.mkdir(outdir)
    os.mkdir("%s/rasters" % outdir)
    os.mkdir("%s/rasters/scratchArc" % outdir)
    os.mkdir("%s/joinTables" % outdir)
arcpy.env.workspace = "%s/rasters/scratchArc" % outdir   
boundShp = gpd.read_file(
            "%s/NHDPlusGlobalData/BoundaryUnit.shp" % NHD_dir).drop(
            ['AreaSqKM','DrainageID','Shape_Area',
            'Shape_Leng','UnitName'], axis=1)

vpus = boundShp.query("UnitType == 'VPU'")
rpus = boundShp.query("UnitType == 'RPU'")
count = 0
upOn = {}  # initialize dict to hold assoc. catchments for lakes lowest HYDROSEQ
for zone in inputs:
    #break
    hr = inputs[zone]
    pre = "%s/NHDPlus%s/NHDPlus%s" % (NHD_dir, hr, zone)    
    # read in lakes and handle caps for column titles, zone 14 improper
    wbShp = gpd.read_file("%s/NHDSnapshot/Hydrography/NHDWaterbody.shp"%(pre))
    wbShp.columns = wbShp.columns[:-1].str.upper().tolist() + ['geometry'] 
    wbShp.drop(['ELEVATION','FCODE','FDATE','GNIS_ID','GNIS_NAME',
                        'REACHCODE','RESOLUTION','SHAPE_AREA','SHAPE_LENG'],
                        axis=1, inplace=True)
    wbShp = wbShp.loc[wbShp['FTYPE'].isin(['LakePond','Reservoir'])]
    vpu = vpus.query("UnitID == '%s'" % zone)
    lakes = sjoin(wbShp, vpu, op='within').drop(['Hydroseq','UnitType','index_right'], axis=1)
    lakes.columns = ['AREASQKM', 'COMID', 'FTYPE', 'geometry', 'VPU']
    fl = dbf2DF("%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(pre))
    cat = gpd.read_file('%s/NHDPlusCatchment/Catchment.dbf'%(pre))
    vaa = dbf2DF('%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(pre))
    #-----------------     allTbls.columns.tolist()  vaa.columns.tolist()
    #  Merge them all together
    allTbls = pd.merge(cat,fl,left_on='FEATUREID',
                       right_on='COMID', how='inner') 
    allTbls = pd.merge(lakes.drop('geometry', axis=1),allTbls,left_on='COMID',
                       right_on='WBAREACOMI', how='left',
                       suffixes = ('_wb','_fl'))
    allTbls = pd.merge(allTbls,vaa,left_on='COMID_fl',
                       right_on='COMID', how='left',
                       suffixes = ('_cat','_vaa'))
    onNetDF = pd.DataFrame(columns=('catCOMID','CatAreaSqKm', 'wbCOMID'))  # create data frame to link catchment COMID to waterbody COMID
    # iterate through table while grouping by the waterbody COMID to select out associated catchments
    # hold on to AREASQKM for comparing with the size of off-network basin creation
    catDict = {}
    for name, group in allTbls.groupby('COMID_wb'):
        if not pd.isnull(group.FEATUREID).any():
            base = group.ix[group.HYDROSEQ.idxmin()]
            onNetDF = onNetDF.append(pd.Series([int(base.COMID_fl),  #there is a chance here that some zones don't cp. all letters!
                                                int(base.COMID_wb),
                                                base.AREASQKM_cat],                                                
                                                index=['catCOMID',
                                                'wbCOMID',
                                                'CatAreaSqKm']),
                                         ignore_index=True)
            catDict[int(base.COMID_fl)] = group.FEATUREID.astype(int).tolist()
    # create numpy arrays for connected catchments to waterbody            
    allLinkd = map(lambda x: catDict[x], catDict.keys())    
    upOn[zone] = {'comids':np.array(catDict.keys()),
                    'lengths':np.array([len(v) for v in allLinkd]),
                    'allLinkd':np.int32(np.hstack(np.array(allLinkd)))}
    onNetDF.to_csv("%s/joinTables/join_%s.csv" % (outdir, zone))
    offLks = lakes.ix[~lakes.COMID.isin(onNetDF.wbCOMID)].copy()  #.drop('index_right', axis=1)
    # make lake and watershed rasters
    for rpu in rasterUnits[zone]:
        #break
        tryLks = offLks.copy()
        if len(rasterUnits[zone]) > 1:
            rpuShp = rpus.query("UnitID == '%s'" % rpu).drop(['Hydroseq','UnitType'], axis=1)
            tryLks = sjoin(tryLks, rpuShp, op='within').drop('index_right', axis=1)
            tryLks.rename(columns={'UnitID': 'RPU'}, inplace=True)
        if len(rasterUnits[zone]) == 1:
            tryLks['RPU'] = rpu
        start = dt.now()
        count_in = 0
        probs = [47]
        while len(probs) > 0:       
            with rs.open("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu)) as rst:
                if rst.crs != tryLks.crs:
                    tryLks.to_crs(rst.crs, inplace=True)
                meta = rst.meta.copy()
                meta.update(compress='lzw')
                meta.update(nodata=0,
                            dtype=rs.uint32,
                            driver='GTiff',
                            crs={'init': u'epsg:5070'})
                with rs.open("%s/rasters/lakes_%s.tif" % (outdir, rpu),
                             'w', **meta) as out:
                    out_arr = out.read(1)
                    shapes = ((g,v) for g,v in zip(tryLks.geometry,tryLks.COMID))
                    burned = features.rasterize(shapes=shapes, fill=0,
                                                out=out_arr,
                                                out_shape=out_arr.shape,
                                                transform=out.transform)
                    out.write(burned.astype(np.uint32), indexes=1)     
            outWshed = Watershed("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu),
                      "%s/rasters/lakes_%s.tif" % (outdir, rpu),
                      "VALUE")
            outWshed.save("%s/rasters/wtshds_%s.tif"%(outdir,rpu))           
            # compare basin sizes ---------------------------------------------
            rat = makeRat("%s/rasters/wtshds_%s.tif"%(outdir,rpu))
            # Write out table for potential use in Arc ------------------------
            DF2dbf(rat, "%s/rasters/wtshds_%s.tif.vat.dbf"%(outdir,rpu))
            centroids = tryLks.to_crs({'init': u'epsg:4269'})  #.copy()
            centroids.geometry = centroids.centroid
            lkCat = sjoin(centroids, cat, op='within')      
            both = pd.merge(lkCat[['COMID','FEATUREID','AreaSqKM']], rat, how='inner', left_on='COMID', right_on='VALUE')
            both['AreaSqKM_basin'] = (both.COUNT * 900) * 1e-6
            probs = both.ix[both.AreaSqKM < both.AreaSqKM_basin].copy()
            probs['diff'] = abs(probs.AreaSqKM - probs.AreaSqKM_basin)
            probs = probs.drop(probs.ix[abs(probs.AreaSqKM - probs.AreaSqKM_basin) < .01].index, axis=0).copy()  #this was in the old script and not sur if it should be thrown out!
            probs['VPU'] = zone
            probs['RPU'] = rpu        
            tryLks.ix[tryLks.COMID.isin(probs.COMID)].index
            tryLks.drop(tryLks.ix[tryLks.COMID.isin(probs.COMID)].index, inplace=True)
            if len(probs) > 0:  #delete files for another iteration of while
                outWshed = None
                rat = None
                purge('%s/rasters' % outdir, rpu)        
                if count_in == 0:
                    problems = probs.copy()
                if count_in > 0:
                    problems = pd.concat([problems, probs])                
            count_in +=1
        print '%s watershed raster built in : %s' % (rpu, str(dt.now() - start))
    if count == 0:
        allOffLks = tryLks.to_crs({'init': u'epsg:5070'}).copy()
        probDF = problems.copy()
    if count > 0:
        allOffLks = pd.concat([allOffLks, tryLks.to_crs({'init': u'epsg:5070'}).copy()])
        probDF = pd.concat([probDF,problems.copy()])
    count += 1
allOffLks.to_file("%s/off_network.shp"%(outdir))
probDF.to_csv()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
np.savez_compressed('%s/onNet_LakeCat.npz' % (outdir), Connect_arrays=upOn)

        
        
#        allTbls.columns.tolist()
#        lakes.columns.tolist()
#                                offLks.to_file("%s/off_network.shp"%(outdir))
#        offLks.head()

#        warp.reproject("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu),
#                       "%s/NHDPlusFdrFac%s/fdr_5070.tif" % (pre, rpu),
#                        dst_crs={'init': u'epsg:5070'})


#writeable_ds = gdal.Open("%s/rasters/lakes_%s.tif" % (outdir, rpu), gdal.GA_Update)                
#band = writeable_ds.GetRasterBand(1)                
#band.Get                
#rat = gdal.RasterAttributeTable()
#rat.CreateColumn("Value", gdalconst.GFT_Integer, gdalconst.GFU_MinMax)
#
#rat_data = band.ReadAsArray(0, 0, writeable_ds.RasterXSize, writeable_ds.RasterYSize)
#vals = np.unique(rat_data).tolist()
#vals.remove(0)
#for i in range(len(vals)):
#        rat.SetValueAsInt(i, 0, int(vals[i]))
#band.WriteArray(rat_data)
#  
#cols = burned.shape[1]
#rows = burned.shape[0]
#driver = gdal.GetDriverByName('GTiff')
#outRaster = driver.Create("%s/rasters/pancakes_%s.tif" % (outdir, rpu), cols, rows, 1, gdal.GDT_UInt32)
#out.transform  
#originX = out.transform[2]
#originY = out.transform[5] 
#pixelWidth = out.transform[0]  
#pixelHeight = out.transform[4]   
#outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
#outband = outRaster.GetRasterBand(1)
#outband.WriteArray(burned)
#outband.FlushCache()


#            if not os.path.exists("%s/NHDPlusFdrFac%s/fdr_5070.tif" % (pre, rpu)):
#                data = rst.read(1)
#                newProj = rst.meta.copy()
#                newProj.update(crs={'init': u'epsg:5070'},
#                                                  driver='GTiff')
#                with rs.open("%s/NHDPlusFdrFac%s/fdr_5070.tif" % (pre, rpu),
#                             'w', **newProj) as out:
#                    out.write(data, indexes=1)



#167679192 in offLks.COMID.values
#
#0: (307, 22, 82, 255)
#1: (124, 83, 97, 255)
#2: (261, 89, 85, 255)
#4: (50, 99, 95, 255)
#8: (144, 80, 36, 255)
#16: (8, 98, 76, 255)
#32: (226, 74, 30, 255)
#64: (185, 65, 92, 255)    
#128: (217, 38, 128, 255)        
#offLks.COMID.min()
#offLks.COMID.max()