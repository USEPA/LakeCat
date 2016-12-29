# -*- coding: utf-8 -*-
"""
Created on Fri May 06 10:30:14 2016

@author: Rdebbout
"""

import os
import sys
import geopandas as gpd
sys.path.append('D:/Projects/LakeCat')
from LakeCat_functions import NHDtblMerge, makeBasins, makeFlowTbls
                 
def main(NHDdir, outdir):
    
    # create directories to hold intermediate data
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        os.mkdir("%s/rasters" % outdir)
        os.mkdir("%s/rasters/scratchArc" % outdir)
        os.mkdir("%s/joinTables" % outdir)
        os.mkdir("%s/flowTables" % outdir)
    
    
    NHDbounds = gpd.read_file(
                            "%s/NHDPlusGlobalData/BoundaryUnit.shp" % NHDdir).drop(
                                            ['AreaSqKM','DrainageID','Shape_Area',
                                             'Shape_Leng','UnitName'], axis=1)
    
    NHDtblMerge(NHDdir, NHDbounds, outdir)
    makeBasins(NHDdir, NHDbounds, outdir)
    makeFlowTbls(NHDdir, outdir)
    
if __name__ == '__main__':
    

outdir = 'D:/Projects/LakeCat/play2'
###############################################################################
#out ='D:/Projects/LakeCat/play62post'
#nhd = 'D:/NHDPlusv21'
#def makeBasins (nhd, out):
#    problems = gpd.GeoDataFrame()  # holding for overdrawn basin delineations
#    inputs = NHDdict(NHDdir)
#    rasterUnits = NHDdict(NHDdir, unit='RPU')
#    rpus = boundShp.query("UnitType == 'RPU'").copy()
#    Obounds = gpd.read_file("%s/out_of_bounds.shp" % out)
#    for zone in inputs:
#        # sjoin and add back-in lakes that are in other zones 
#        #break
#        print zone
#        hr = inputs[zone]
#        pre = "%s/NHDPlus%s/NHDPlus%s" % (nhd, hr, zone)
#        addLks = Obounds.ix[Obounds.UnitID == zone].copy()
#        offLks = gpd.read_file("%s/off_net_%s.shp" % (out, zone))
#        offLks = pd.concat([offLks,addLks]).reset_index().drop('index',axis=1)
#        # make lake and watershed rasters
#        r_count = 0
#        for rpu in rasterUnits[zone]:
#            #break
#            lakes = offLks.copy()
#            if len(rasterUnits[zone]) > 1:
#                rpuShp = rpus.query("UnitID == '%s'" % rpu).drop(['Hydroseq',
#                                                                'UnitType'],axis=1)
#                lakes = sjoin(lakes, rpuShp, op='within').drop('index_right',
#                                                                axis=1)
#                lakes.rename(columns={'UnitID': 'RPU'}, inplace=True)
#            if len(rasterUnits[zone]) == 1:
#                lakes['RPU'] = rpu
#            start = dt.now()
#            count_in = 0
#            bigs = [47]
#            # while len(bigs) > 0:       
#            with rs.open("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu)) as fdr:
#                if fdr.crs != lakes.crs:
#                    lakes.to_crs(fdr.crs, inplace=True)
#                meta = fdr.meta.copy()
#                meta.update(compress='lzw')
#                meta.update(nodata=0,
#                            dtype=rs.uint32,
#                            driver='GTiff',
#                            crs={'init': u'epsg:5070'})
#                with rs.open("%s/rasters/lakes_%s.tif" % (outdir, rpu),
#                             'w', **meta) as lksRas:
#                    lksArray = lksRas.read(1)
#                    shapes = ((g,v) for g,v in zip(lakes.geometry,lakes.COMID))
#                    burned = features.rasterize(shapes=shapes, fill=0,
#                                                out=lksArray,
#                                                out_shape=lksArray.shape,
#                                                transform=lksRas.transform)
#                    lksRas.write(burned.astype(np.uint32), indexes=1)    
#            outWshed = Watershed("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu),
#                      "%s/rasters/lakes_%s.tif" % (outdir, rpu),
#                      "VALUE")
#            outWshed.save("%s/rasters/wtshds_%s.tif"%(outdir,rpu))           
#            # compare basin sizes ---------------------------------------------
#            rat = makeRat("%s/rasters/wtshds_%s.tif"%(outdir,rpu))
#            # Write out table for potential use in Arc ------------------------
#            DF2dbf(rat, "%s/rasters/wtshds_%s.tif.vat.dbf"%(outdir,rpu))
#            centroids = lakes.to_crs({'init': u'epsg:4269'}).copy()
#            centroids.geometry = centroids.centroid
#            cat = gpd.read_file('%s/NHDPlusCatchment/Catchment.shp'%(pre))
#            lkCat = sjoin(centroids, cat, op='within', how='left')      
#            both = pd.merge(lkCat[['COMID','FEATUREID','AreaSqKM']], rat,
#                            how='inner', left_on='COMID', right_on='VALUE')
#            both['AreaSqKM_basin'] = (both.COUNT * 900) * 1e-6
#            bigs = both.ix[both.AreaSqKM < both.AreaSqKM_basin].copy()
#            bigs['diff'] = abs(bigs.AreaSqKM - bigs.AreaSqKM_basin)
#            #bigs = bigs.drop(bigs.ix[abs(bigs.AreaSqKM - bigs.AreaSqKM_basin) < .01].index, axis=0).copy()  #this was in the old script and not sur if it should be thrown out!
#            bigs['VPU'] = zone
#            bigs['RPU'] = rpu
#            problems = pd.concat([problems,bigs])
#            print problems
#    problems.to_csv(,index=False)
#    
#    
#    
#            
#            lakes.ix[lakes.COMID.isin(bigs.COMID)].index
#            lakes.drop(lakes.ix[lakes.COMID.isin(bigs.COMID)].index,
#                                  inplace=True)
#            if len(bigs) > 0:  #delete files for another iteration of while
#                outWshed = None
#                rat = None
#                purge('%s/rasters' % outdir, rpu)        
#                if count_in == 0:
#                    problems = bigs.copy()
#                if count_in > 0:
#                    problems = pd.concat([problems, bigs])                
#            count_in +=1
#                
#                
#                
#                
#                
#            print '%s watershed raster built in : %s' % (rpu, str(dt.now() - start))
#            if r_count == 0:
#                p = problems.copy()
#            if r_count > 0:
#                p = pd.concat([p, problems.copy()])
#            r_count += 1
#        if count == 0:
#            outOff = lakes.to_crs({'init': u'epsg:5070'}).copy()
#            probDF = problems.copy()
#        if count > 0:
#            outOff = pd.concat([outOff, lakes.to_crs({'init': u'epsg:5070'}).copy()])
#            probDF = pd.concat([probDF,problems.copy()])
#        count += 1
#    outOff.to_file("%s/off_network.shp" % outdir)
#    probDF.sort_values('diff').to_csv("%s/problemLakes.csv" % outdir, index=False)
#    print "You're ready to start processing Landscape Layers with LakeCat!"
###############################################################################        
        
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