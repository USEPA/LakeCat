#--------------------------------
# Name: Lakes_zonal_stats
# Purpose: run zonal stats on
#          gridded watersheds for
#          'off-network' NHDPlus lakes
# Author: Marc Weber
# Created 7/25/2013
# Updated 10/1/2014
# ArcGIS Version:  10.2
# Python Version:  2.7
#--------------------------------
import os, sys, arcpy
from arcpy import env
from arcpy.sa import *
from datetime import datetime as dt
import winsound, pymsgbox
#####################################################################################################################
def selectHR(filestr):
    hr = filestr[3:5]
    if filestr[3:6] in ['03a','03b']:
        hr = '03N'
    if filestr[3:6] in ['03c','03d']:
        hr = '03S'
    if filestr[3:6] in ['03e','03f']:
        hr = '03W'  
    if filestr[3:6] == '03g':
        hr = '08'
    if filestr[3:6] in ['10a','10b','10c','10d']:
        hr = '10L'
    if filestr[3:6] in ['10e','10f','10g','10h','10i']:
        hr = '10U'
    return hr
#####################################################################################################################
# Settings
# Check out Spatial Analyst extension license
arcpy.CheckOutExtension("spatial")
arcpy.OverWriteOutput = 1
# Set Geoprocessing environments
# arcpy.env.resamplingMethod = "NEAREST"
arcpy.env.outputCoordinateSystem = "PROJCS['NAD_1983_Albers',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-96.0],PARAMETER['Standard_Parallel_1',29.5],PARAMETER['Standard_Parallel_2',45.5],PARAMETER['Latitude_Of_Origin',23.0],UNIT['Meter',1.0]]"
# arcpy.env.parallelProcessingFactor = "60%"
arcpy.env.pyramid = "NONE"
#processing cell size
arcpy.env.cellSize = "30"

#####################################################################################################################
#landscape_layer = ['clay','huden2010','imp2006','om','perm','popden2010','precip','rckdep','roadstrm','runoff', 'sand', 'wtdep']
landscape_layer = ['pestic']

# Variables
raster_dir = "L:/Priv/CORFiles/Geospatial_Library/Data/Project/SSWR1.1B/LandscapeRasters/QAComplete"
lake_dir = "D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters"

inputs = {'CA':['18'],'CO':['14','15'],'GB':['16'],'GL':['04'],'MA':['02'],'MS':['05','06','07','08','10L','10U','11'],'NE':['01'],'PN':['17'],\
          'RG':['13'],'SA':['03N','03S','03W'],'SR':['09'],'TX':['12']}

mask = 'False'
categorical = 'False'      
#####################################################################################################################
startTime = dt.now()
for l in landscape_layer:
    if not l.count('elev'):
        ingrid = "%s/%s.tif"%(raster_dir,l)
    out_dir = "D:/Projects/lakesAnalysis/Output/Zonal_Stats_Output2/%s"%(l)
    if mask=='True':
        out_dir = "D:/Projects/lakesAnalysis/Output/%s_MidSlp"%(l)
        #os.mkdir(out_dir)
    print out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for files in os.listdir(lake_dir):
        if files.count('.tif') and files.count('wtshds') and not files.count('.vat') and not files.count('.xml'):
            print files
			hr = selectHR(files)
            #get NHD elevation raster if landscape layer is elevation              
            if l.count('elev'):
                for regions in inputs.keys():
                    if hr in inputs[regions]:
                        rgn = regions                                    
                ingrid = "C:/Users/Rdebbout/temp/NHDPlusV21/NHDPlus%s/NHDPlus%s/NEDSnapshot/Ned%s/elev_cm"%(rgn, hr,files[3:6])               
            # Set variables for zonal stats
            inZoneData = "%s/%s"%(lake_dir,files)
            arcpy.env.snapRaster = inZoneData
            arcpy.env.resamplingMethod = "NEAREST"
            #use mask just for riparian forested
            if l.count('nlcd'):
                arcpy.AddField_management(inZoneData, "CALC_FIELD", "LONG", "", "")
                arcpy.CalculateField_management(inZoneData, "CALC_FIELD", '!OID!', "PYTHON")            
            if mask=='True' and categorical=='True':         
                if ingrid.count('.tif') or ingrid.count('.img'):
                    outTable ="%s/zonalstats_%s%s_MidSlp.dbf"%(out_dir,ingrid.split("/")[-1].split(".")[0],files[3:6])
                else:
                    outTable ="%s/zonalstats_%s%s_MidSlp.dbf"%(out_dir,ingrid.split("/")[-1],files[3:6])
                arcpy.env.mask = "L:\\Priv\\CORFiles\\Geospatial_Library\\Data\Project\\SSWR1.1B\\LandscapeRasters\\QAComplete\\midslope_mask.tif"
                if not os.path.exists(outTable):
                    arcpy.gp.TabulateArea_sa(inZoneData, "CALC_FIELD", ingrid, "Value", outTable, "30")
                    arcpy.JoinField_management(outTable,"CALC_FIELD",inZoneData,"CALC_FIELD",fields="Value")
                    arcpy.DeleteField_management(inZoneData, "CALC_FIELD")
                    arcpy.DeleteField_management(outTable, "CALC_FIELD")
            if mask=='True' and categorical!='True':
                if ingrid.count('.tif') or ingrid.count('.img'):
                    outTable ="%s/zonalstats_%s%srip100.dbf"%(out_dir,ingrid.split("/")[-1].split(".")[0],files[3:6])
                else:
                    outTable ="%s/zonalstats_%s%srip100.dbf"%(out_dir,ingrid.split("/")[-1],files[3:6])
                arcpy.env.mask = "L:\\Priv\\CORFiles\\Geospatial_Library\\Data\Project\\SSWR1.1B\\PhaseTwo\\LandscapeRasters\\QAcomplete\\WaterMask\\Mosaics\\RipBuf100_%s.tif"%(hr)
                if not os.path.exists(outTable):
                    arcpy.gp.ZonalStatisticsAsTable_sa(inZoneData, "Value", ingrid, outTable, "DATA", "ALL")
            if mask!='True' and categorical=='True':
                # Execute Tabulate Area
                if ingrid.count('.tif') or ingrid.count('.img'):
                    outTable ="%s/zonalstats_%s%s.dbf"%(out_dir,ingrid.split("/")[-1].split(".")[0],files[3:6])
                else:
                    outTable ="%s/zonalstats_%s%s.dbf"%(out_dir,ingrid.split("/")[-1],files[3:6])
                if not os.path.exists(outTable):
                    arcpy.gp.TabulateArea_sa(inZoneData, "Value", ingrid, "VALUE", outTable, "30")
                    arcpy.JoinField_management(outTable,"CALC_FIELD",inZoneData,"CALC_FIELD",fields="Value")
                    arcpy.DeleteField_management(inZoneData, "CALC_FIELD")
                    arcpy.DeleteField_management(outTable, "CALC_FIELD")
            elif mask!='True' and categorical!='True':
                if ingrid.count('.tif') or ingrid.count('.img'):
                    outTable ="%s/zonalstats_%s%s.dbf"%(out_dir,ingrid.split("/")[-1].split(".")[0],files[3:6])
                else:
                    outTable ="%s/zonalstats_%s%s.dbf"%(out_dir,ingrid.split("/")[-1],files[3:6])
                if not os.path.exists(outTable):
                    arcpy.gp.ZonalStatisticsAsTable_sa(inZoneData, "Value", ingrid, outTable, "DATA", "ALL")
                pass
winsound.Beep(300,2000)
pymsgbox.alert('Script is done!', 'Your data has been processed!')
print "elapsed time " + str(datetime.now()-startTime) + " | to process: " + ingrid


# #this is the lakes_zonal_stats_NLCD script if needed to refine script above when completing for all metrics
# #--------------------------------
# # Name: Lakes_zonal_stats
# # Purpose: run zonal stats on
# #          gridded watersheds for
# #          'off-network' NHDPlus lakes
# # Author: Marc Weber
# # Created 7/25/2013
# # Updated 10/1/2014
# # ArcGIS Version:  10.2
# # Python Version:  2.7
# #--------------------------------
# import os, sys, arcpy
# from arcpy import env
# from arcpy.sa import *
# from datetime import datetime
# #####################################################################################################################
# # Settings
# # Check out Spatial Analyst extension license
# arcpy.CheckOutExtension("spatial")
# arcpy.OverWriteOutput = 1
# # Set Geoprocessing environments
# # arcpy.env.resamplingMethod = "NEAREST"
# arcpy.env.outputCoordinateSystem = "PROJCS['NAD_1983_Albers',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-96.0],PARAMETER['Standard_Parallel_1',29.5],PARAMETER['Standard_Parallel_2',45.5],PARAMETER['Latitude_Of_Origin',23.0],UNIT['Meter',1.0]]"
# # arcpy.env.parallelProcessingFactor = "60%"
# arcpy.env.pyramid = "NONE"
# #processing cell size
# arcpy.env.cellSize = "30"
# #####################################################################################################################
# # Variables
# raster_dir = "D:/Projects/lakesAnalysis/MetricData"
# ingrid = "%s/us_lithology_1km.tif"%(raster_dir)
# lake_dir = "D:/Projects/lakesAnalysis/Output/NHDPlus_Lakes_Basins_Rasters(copy)"
# #out_dir = "H:/Watershed Integrity Spatial Prediction/results"
# out_dir = "D:/Projects/lakesAnalysis/Output/lith"
# os.mkdir(out_dir)
# mask = 'False'
# categorical = 'True'        
# #####################################################################################################################
# startTime = datetime.now()
# for files in os.listdir(lake_dir):
    # if files.count('.tif') and files.count('wtshds') and not files.count('.vat') and not files.count('.xml'):
        # print files
        # # Set variables for zonal stats
        # inZoneData = "%s/%s"%(lake_dir,files)
        # arcpy.env.snapRaster = inZoneData
        # arcpy.env.resamplingMethod = "NEAREST"
# #        add field and populate with OID to use as zone field in Tabulate Area
        # arcpy.AddField_management(inZoneData, "CALC_FIELD", "LONG", "", "")
        # arcpy.CalculateField_management(inZoneData, "CALC_FIELD", '!OID!', "PYTHON")
        # if mask!='True' and categorical=='True':			
            # # Execute Tabulate Area
            # if ingrid.count('.tif') or ingrid.count('.img'):
                # outTable ="%s/zonalstats_%s%s.dbf"%(out_dir,ingrid.split("/")[-1].split(".")[0],files[3:6])
            # else:
                # outTable ="%s/zonalstats_%s%s.dbf"%(out_dir,ingrid.split("/")[-1],files[3:6])
            # if not os.path.exists(outTable):
# #                outExtractByMask = ExtractByMask(ingrid, inZoneData)
# #                extract = "%s/extract_%s"%(out_dir,files[3:6])
# #                outExtractByMask.save(extract)
                # arcpy.gp.TabulateArea_sa(inZoneData, "CALC_FIELD", ingrid, "VALUE", outTable, "30")
                # arcpy.JoinField_management(outTable,"CALC_FIELD",inZoneData,"CALC_FIELD",fields="Value")
                # arcpy.DeleteField_management(inZoneData, "CALC_FIELD")
                # print "process finished: " + files[3:6] 

# print "elapsed time " + str(datetime.now()-startTime) + " | to process NLCD" 