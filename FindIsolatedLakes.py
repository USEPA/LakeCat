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

NHD_dir = 'D:/NHDPlusV21'                    
# dictionary to iterate thru NHD folder structure
inputs = makeVPUdict(NHD_dir)
#--------------------------------------------------------

# !!!DIRECTORIES TO WRITE INFO TO!!!
# Shapefile or gdb feature layer to store on and off-network lakes
outfeature = 'D:/Projects/lakesAnalysis/new/NetworkLakes.shp'
outfeature2 = 'D:/Projects/lakesAnalysis/new/IsolatedLakes.shp'

# count is used when iterating to use arcpy.CopyFeatures_management() tool first
# pass thru and then arcpy.Append_management() the remaining iterations to add 
# selected lakes to the files
count = 0
for zone in inputs:
    hr = inputs[zone]
    print 'on region ' + hr + ' and hydro number ' + zone
    # create directories to NHD data
    NHDdir = "%s/NHDPlus%s/NHDPlus%s" % (NHD_dir, hr, zone)
    NHDWaterbody = "%s/NHDSnapshot/Hydrography/NHDWaterbody.shp"%(NHDdir)
    NHDWaterbodyDBF = "%s/NHDSnapshot/Hydrography/NHDWaterbody.dbf"%(NHDdir)
    NHDFlowlineDBF = "%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(NHDdir)
    NHDPlusFlowlinVAADBF = '%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(NHDdir)
    NHDCatchmentDBF = '%s/NHDPlusCatchment/Catchment.dbf'%(NHDdir)     
    # get waterbodies and select out by FTYPE
    selStmnt = "\"FTYPE\" IN ( 'LakePond' , 'Reservoir' )"
    arcpy.MakeFeatureLayer_management(NHDWaterbody, "NHDWaterbody", selStmnt)
    #--------------------------------------------------------
    wb = dbf2DF(NHDWaterbodyDBF)
    wb = wb.loc[wb['FTYPE'].isin(['LakePond','Reservoir'])]
    #--------------------------------------------------------
    # get flowline table into pandas data frame 
    fl = dbf2DF(NHDFlowlineDBF)
    #--------------------------------------------------------
    # get catchment table into pandas data frame 
    cat = dbf2DF(NHDCatchmentDBF)
    #--------------------------------------------------------
    flowcat = pd.merge(cat,fl,left_on='FEATUREID',right_on='COMID', how='inner') 
    join = pd.merge(wb,flowcat,left_on='COMID', right_on='WBAREACOMI', how='left')
    #--------------------------------------------------------
    # get NHDPlusFlowlineVAA table into pandas data frame
    vaa = dbf2DF(NHDPlusFlowlinVAADBF) 
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
       arcpy.CopyFeatures_management('NHDWaterbody', outfeature)
    if count > 0:
        arcpy.Append_management("NHDWaterbody",outfeature,"NO_TEST")
    #Switch Selection
    arcpy.SelectLayerByAttribute_management("NHDWaterbody", "SWITCH_SELECTION")
    #add remaining waterbodies to IsolatedLakes 
    if count == 0:
       arcpy.CopyFeatures_management ('NHDWaterbody', outfeature2)
    if count > 0:
        arcpy.Append_management("NHDWaterbody",outfeature2,"NO_TEST")
    arcpy.Delete_management("NHDWaterbody")
    count+=1

bnd_dir = "%s/NHDPlusGlobalData/BoundaryUnit.shp" % NHD_dir
    
# get waterbodies and select out by FTYPE
arcpy.MakeFeatureLayer_management(bnd_dir, "Boundaries")
arcpy.MakeFeatureLayer_management(outfeature2, "IsoLakes")

# Select the 969 waterbodies that are NOT within RPU boundaries and remove
arcpy.SelectLayerByLocation_management("IsoLakes", "COMPLETELY_WITHIN", 
                                           "Boundaries", "", "NEW_SELECTION")
arcpy.SelectLayerByAttribute_management("IsoLakes", "SWITCH_SELECTION")
arcpy.DeleteFeatures_management("IsoLakes")
arcpy.Delete_management("IsoLakes")
arcpy.Delete_management("Boundaries")
