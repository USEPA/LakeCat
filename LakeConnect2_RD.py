# -*- coding: utf-8 -*-
"""
Created on Fri May 06 15:52:44 2016

@author: Rdebbout
"""

# Import arcpy module
import os
import arcpy
from arcpy.sa import Raster, Con, Combine, IsNull
from datetime import datetime

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")
arcpy.env.workspace = "D:/Projects/lakesAnalysis/work2"

# Set overwrite = T
arcpy.env.overwriteOutput = True
# Loop through each region to process
inputs = {'CA':['18'],'CO':['14','15'],'GB':['16'],'GL':['04'],'MA':['02'],'MS':['05','06','07','08','10L','10U','11'],
          'NE':['01'],'PN':['17'],'RG':['13'],'SA':['03N','03S','03W'],'SR':['09'],'TX':['12']}
          
out_dir = 'D:/Projects/lakesAnalysis/From_To_Tables'        
for regions in inputs.keys():
    for hydro in inputs[regions]:
        print 'on region ' + regions + ' and hydro number ' + hydro
        nhddir = "D:/NHDPlusV21/NHDPlus%s/NHDPlus%s"%(regions, hydro)
        for dirs in os.listdir(nhddir):
            if dirs.count("FdrFac") and not dirs.count('.txt')and not dirs.count('.7z'):
                if not arcpy.Exists("%s/LkFrmTo_R%s.tif"%(out_dir, dirs[-3:])):
                    print dirs
                    fdr = Raster("%s/%s/fdr"%(nhddir,dirs))
                    Wtshds = Raster("D:/Projects/LakesAnalysis/NHDPlus_Lakes_Basins_Rasters_1/reg%s_wtshds.tif"%(dirs[-3:]))
                    startTime = datetime.now()
                    # Process: Shift
                    print 'shifting'
                    shift1 = arcpy.Shift_management(Wtshds, "shift1", "-30", "0", Wtshds)
                    shift2 = arcpy.Shift_management(Wtshds, "shift2", "-30", "30", Wtshds)
                    shift4 = arcpy.Shift_management(Wtshds, "shift4", "0", "30", Wtshds)
                    shift8 = arcpy.Shift_management(Wtshds, "shift8", "30", "30", Wtshds)
                    shift16 = arcpy.Shift_management(Wtshds, "shift16", "30", "0", Wtshds)
                    shift32 = arcpy.Shift_management(Wtshds, "shift32", "30", "-30", Wtshds)
                    shift64 = arcpy.Shift_management(Wtshds, "shift64", "0", "-30", Wtshds)
                    shift128 = arcpy.Shift_management(Wtshds, "shift128", "-30", "-30", Wtshds)


                    # Process: Raster Calculator
                    print 'running raster calculator on shifts'
                    flowto1 = ((shift1 != Wtshds) * (fdr == 1)) * shift1
                    flowto2 = ((shift2 != Wtshds) * (fdr == 2)) * shift2
                    flowto4 = ((shift4 != Wtshds) * (fdr == 4)) * shift4
                    flowto8 = ((shift8 != Wtshds) * (fdr == 8)) * shift8
                    flowto16 = ((shift16 != Wtshds) * (fdr == 16)) * shift16
                    flowto32 = ((shift32 != Wtshds) * (fdr == 32)) * shift32
                    flowto64 = ((shift64 != Wtshds) * (fdr == 64)) * shift64
                    flowto128 = ((shift128 != Wtshds) * (fdr == 128)) * shift128

                    FlowToSum = Con(IsNull(flowto1), 0, flowto1) + Con(IsNull(flowto2), 0, flowto2) + Con(IsNull(flowto4), 0, flowto4) + Con(IsNull(flowto8), 0, flowto8) + Con(IsNull(flowto16), 0, flowto16) + Con(IsNull(flowto32), 0, flowto32) + Con(IsNull(flowto64), 0, flowto64) + Con(IsNull(flowto128), 0, flowto128)
                    FlowToFinal = Con(FlowToSum != 0, FlowToSum)
                    outCombine = Combine([FlowToFinal, Wtshds])
                    outCombine.save("%s/LkFrmTo_R%s.tif" % (out_dir, dirs[-3:]))
                    rasterList = arcpy.ListRasters("*", "")
                    for raster in rasterList:
                        arcpy.Delete_management(raster)
                    outDbf = "%s/LkFrmTo_R%s.dbf" % (out_dir, dirs[-3:])
                    if not arcpy.Exists(outDbf):
                        arcpy.CopyRows_management(outCombine, outDbf, "")
                    print "total elapsed time " + str(datetime.now()-startTime)