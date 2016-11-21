# -*- coding: utf-8 -*-
"""
Created on Fri May 06 10:30:14 2016

@author: Rdebbout
"""

import sys
import arcpy
import pandas as pd
sys.path.append('D:/Projects/lakesAnalysis/Scripts')
from LakeCat_functions import dbf2DF, makeVPUdict

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

nhd = 'D:/NHDPlusV21'                    
# dictionary to iterate thru NHD folder structure
inputs = makeVPUdict(nhd)
#--------------------------------------------------------
outdir = 'D:/Projects/lakesAnalysis/new'
# !!!DIRECTORIES TO WRITE INFO TO!!!
# Shapefile or gdb feature layer to store on and off-network lakes
outOn = '{}/NetworkLakes.shp'.format(outdir)
outOff = '{}/IsolatedLakes.shp'.format(outdir)

# count is used when iterating to use arcpy.CopyFeatures_management() tool first
# pass thru and then arcpy.Append_management() the remaining iterations to add 
# selected lakes to the files
count = 0
for zone in inputs:
    hr = inputs[zone]
    nhddir = "%s/NHDPlus%s/NHDPlus%s" % (nhd, hr, zone)
    
    #not needed til below for selecting out network/off-network lakes
    #also add RPU and VPU and at the end add UID
    NHDWaterbody = "%s/NHDSnapshot/Hydrography/NHDWaterbody.shp"%(nhddir)
    
    # get waterbodies and select out by FTYPE
    #selStmnt = "\"FTYPE\" IN ( 'LakePond' , 'Reservoir' )"
    #arcpy.MakeFeatureLayer_management(NHDWaterbody, "NHDWaterbody", selStmnt)
    #--------------------------------------------------------
    wb = dbf2DF("%s/NHDSnapshot/Hydrography/NHDWaterbody.dbf"%(nhddir))
    wb = wb.loc[wb['FTYPE'].isin(['LakePond','Reservoir'])]
    #--------------------------------------------------------
    # get flowline table into pandas data frame 
    fl = dbf2DF("%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(nhddir))
    #--------------------------------------------------------
    # get catchment table into pandas data frame 
    cat = dbf2DF('%s/NHDPlusCatchment/Catchment.dbf'%(nhddir))
    #--------------------------------------------------------
    flowcat = pd.merge(cat,fl,left_on='FEATUREID',right_on='COMID', how='inner') 
    join = pd.merge(wb,flowcat,left_on='COMID', right_on='WBAREACOMI', how='left')
    #--------------------------------------------------------
    # get NHDPlusFlowlineVAA table into pandas data frame
    vaa = dbf2DF('%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(nhddir)) 
    #--------------------------------------------------------
    # Merge VAA table on to joined table of catchment, waterbody, and flowlines
    joinVAA = pd.merge(join,vaa,left_on='COMID_y',right_on='COMID', how='left')
    # create data frame to link catchment COMID to waterbody COMID
    df = pd.DataFrame(columns=('catCOMID', 'wbCOMID'))
    # iterate through table while grouping by the waterbody COMID to select out associated catchments
    for name, group in joinVAA.groupby('COMID_x'):
        if not pd.isnull(group.FEATUREID).any():
            hydro = pd.Series.min(group.HYDROSEQ)
            group = group[group.HYDROSEQ == hydro]           
            s = pd.DataFrame({'catCOMID':pd.Series([int(group.COMID_y.values[0])]),'wbCOMID':pd.Series([int(group.COMID_x.values[0])])})
            df = df.append(s)
            
    # convert list values into integers and format to pass into SelectLayerByAttribute function 
    wblist = [int(x) for x in df.wbCOMID.tolist()]            
    offstring = '"COMID" IN (%s)' % wblist   
    whereclause = offstring.translate(None,'[]')    
    # Process: Select Layer By Attribute
    arcpy.SelectLayerByAttribute_management("NHDWaterbody", "NEW_SELECTION", whereclause)
    #add waterbodies to NetworkLakes
    if count == 0:
       arcpy.CopyFeatures_management('NHDWaterbody', outOn)
    if count > 0:
        arcpy.Append_management("NHDWaterbody",outOn,"NO_TEST")
    #Switch Selection
    arcpy.SelectLayerByAttribute_management("NHDWaterbody", "SWITCH_SELECTION")
    #add remaining waterbodies to IsolatedLakes 
    if count == 0:
       arcpy.CopyFeatures_management ('NHDWaterbody', outOff)
    if count > 0:
        arcpy.Append_management("NHDWaterbody",outOff,"NO_TEST")
    arcpy.Delete_management("NHDWaterbody")
    count+=1

bnd_dir = "%s/NHDPlusGlobalData/BoundaryUnit.shp" % nhd
    
# get waterbodies and select out by FTYPE
arcpy.MakeFeatureLayer_management(bnd_dir, "Boundaries")
arcpy.MakeFeatureLayer_management(outOff, "IsoLakes")

# Select the 969 waterbodies that are NOT within RPU boundaries and remove
arcpy.SelectLayerByLocation_management("IsoLakes", "COMPLETELY_WITHIN", 
                                           "Boundaries", "", "NEW_SELECTION")
arcpy.SelectLayerByAttribute_management("IsoLakes", "SWITCH_SELECTION")
arcpy.DeleteFeatures_management("IsoLakes")
arcpy.Delete_management("IsoLakes")
arcpy.Delete_management("Boundaries")
