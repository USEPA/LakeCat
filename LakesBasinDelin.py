#--------------------------------------------------------
# Name: LakesBasinZonalStats.py
# Purpose: Select NHDPlus lakes that are 'off-network',
#          grid the lakes, build watersheds by raster
#          processing unit, and run zonal statistics
#          on basins for selected landscape rasters
#
# Author: Marc Weber
# Created 10/21/2014
# ArcGIS Version:  10.2.1
# Python Version:  2.7
#--------------------------------------------------------

import os
import geopandas as gpd
from arcpy import env, PolygonToRaster_conversion
from arcpy.sa import Watershed 

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Set Geoprocessing environments

crs = 'PROJCS["NAD_1983_Contiguous_USA_Albers",'\
'GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",'\
'SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],'\
'UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],'\
'PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],'\
'PARAMETER["central_meridian",-96.0],PARAMETER["standard_parallel_1",29.5],'\
'PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_origin",23.0],'\
'UNIT["Meter",1.0]]'

env.outputCoordinateSystem = crs

NHD_dir = 'D:/NHDPlusV21'                   
outdir = 'D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_1'
if not os.path.exists(outdir):
    os.mkdir(outdir)
wbs = "D:/Projects/lakesAnalysis/new/Isolated_AlbersOut1_Arc.shp"  
   
#RPU = 'D:/Projects/lakesAnalysis/RPUs_Global_ALBERS.shp'  
#arcpy.MakeFeatureLayer_management(RPU, "RPUs") 
rpu = gpd.GeoDataFrame.from_file('{}/NHDPlusGlobalData/BoundaryUnit.shp'.format(NHD_dir))
wb = gpd.GeoDataFrame.from_file(wbs)

for x in wb.VPU.unique():
    print x
# Both VPU and RPU should be put into the shp of Isolated so we can lookup from there
# else we are just copying what has already happened in the last script, also, this
# might be able to be modularized to be used in while loop for removing lakes we can't use
rows = arcpy.SearchCursor(RPU)
row = rows.next()
while row:
    unit = row.getValue("UnitID")
    hydroregion = row.getValue("DrainageID")
    lakesras = '%s/test_reg%s_lakes.tif'%(outdir,str(unit))
    if not os.path.exists(lakesras):
        print unit
        for root, dirs, files in os.walk(NHD_dir):
            for name in dirs:
                if unit in name and 'FdrFac' in name:
                    fdr = (os.path.join(root, name) + '/fdr')
        env.snapRaster = fdr
        desc = arcpy.Describe(fdr)
        env.extent = desc.extent   
        arcpy.MakeFeatureLayer_management("RPUs", "cursor_layer", "UnitID = '%s'" % str(unit) )
        arcpy.SelectLayerByLocation_management("NHDWaterbody", "HAVE_THEIR_CENTER_IN", "cursor_layer", "", "NEW_SELECTION")
        print 'selected'
        #arcpy.CopyFeatures_management("NHDWaterbody", "%s/lks%s.shp" % (outdir, unit))    
        PolygonToRaster_conversion("NHDWaterbody", "GRIDCODE", lakesras, "CELL_CENTER", "NONE", "30")
        arcpy.Delete_management("cursor_layer")
        print lakesras
        print 'done'
    outwatshd = '%s/reg%s_wtshds.tif'%(outdir,unit)
    if not os.path.exists(outwatshd):
        print 'making watershed'
        Watershed(fdr, lakesras, outwatshd, "VALUE")    
    row = rows.next()
##########################################################################################################   