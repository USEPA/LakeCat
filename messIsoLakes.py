# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:31:07 2016

@author: Rdebbout
"""

import sys
import arcpy
import pandas as pd
from collections import OrderedDict
sys.path.append('D:/Projects/StreamCat')
from StreamCat_functions import dbf2DF

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# List to iterate thru for table names of our existing StreamCat dataset
metricTableNames = ['AgMidHiSlopes_Region','CanalDensity_Region','Dams_Region','Elevation_Region',\
                    'EPA_FRS_Region','ImperviousSurfaces2006_Region','Lithology_Region','Mines_Region',\
                    'NLCD2006_Region','NLCD2006LnBf100_Region','NLCD2011_Region','PRISMPrecip_Region',\
                    'RoadDensity_Region','RoadStreamCrossings_Region','Runoff_Region','STATSGO_Set1_Region',\
                    'STATSGO_Set2_Region','USCensus2010_Region']   
                    
# DataFrame created to link together data to dig into NHD folder structure and be iterated thru
inputs = OrderedDict([('10U', 'MS'), ('10L', 'MS'), ('07', 'MS'), ('11', 'MS'), ('06', 'MS'),
                      ('05', 'MS'), ('08', 'MS'), ('01', 'NE'), ('02', 'MA'), ('03N', 'SA'),
                      ('03S', 'SA'), ('03W', 'SA'), ('04', 'GL'), ('09', 'SR'), ('12', 'TX'),
                      ('13', 'RG'), ('14', 'CO'), ('15', 'CO'), ('16', 'GB'), ('17', 'PN'),
                      ('18', 'CA')])
#--------------------------------------------------------
# !!!DIRECTORIES TO WRITE INFO TO!!!
# directory for resulting metric tables madeup of on-network catchments
output_dir = 'D:/Projects/lakesAnalysis/Output/on_network2'
# Shapefile or gdb feature layer to store on and off-network lakes
outfeature = 'D:/Projects/lakesAnalysis/IsolatedLakes_RD.gdb/NetworkLakes'
outfeature2 = 'D:/Projects/lakesAnalysis/IsolatedLakes_RD.gdb/IsolatedLakes'
#--------------------------------------------------------
# directory where StreamCat tables are located
metricTable_dir = 'L:/Priv/CORFiles/Geospatial_Library/Data/Project/SSWR1.1B/FTP_Staging/StreamCat/HydroRegions'
# !!!CHANGE NHDdir below to your NHD directory also and the rest is good to go!!!
#--------------------------------------------------------

# count is used when iterating to use arcpy.CopyFeatures_management() tool first pass thru
# and then arcpy.Append_management() the remaining iterations to add selected lakes to the files
count = 0
startTime = datetime.now()
for zone in inputs:
    hr = inputs[zone]
    # create directories to NHD data
    NHDdir = "D:/NHDPlusV21/NHDPlus%s/NHDPlus%s"%(hr,zone)
   
    # get waterbodies and select out by FTYPE
    arcpy.MakeFeatureLayer_management("%s/NHDSnapshot/Hydrography/NHDWaterbody.shp"%(NHDdir), "NHDWaterbody", "\"FTYPE\" IN ( 'LakePond' , 'Reservoir' )")
    arcpy.AddField_management("NHDWaterbody","VPU","TEXT")
    arcpy.CalculateField_management("NHDWaterbody","VPU",repr(hr[0]), "PYTHON_9.3")
    #--------------------------------------------------------
    wb = dbf2DF("%s/NHDSnapshot/Hydrography/NHDWaterbody.dbf"%(NHDdir))[['COMID','FTYPE']]
    wb = wb.loc[wb['FTYPE'].isin(['LakePond','Reservoir'])]
    wb.columns = ['WBCOMID','FTYPE']
    #--------------------------------------------------------
    # get flowline table into pandas data frame 
    fl = dbf2DF("%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(NHDdir))[['COMID','WBAREACOMI']]
    fl.columns = ['FLCOMID','WBAREACOMI']
    #--------------------------------------------------------
    # get catchment table into pandas data frame 
    cat = dbf2DF('%s/NHDPlusCatchment/Catchment.dbf'%(NHDdir))[['FEATUREID']]
    #--------------------------------------------------------
    flowcat = pd.merge(cat,fl,left_on='FEATUREID',right_on='FLCOMID', how='right')# Change how to right from inner!
    join = pd.merge(wb,flowcat,left_on='WBCOMID', right_on='WBAREACOMI', how='left')
    #--------------------------------------------------------
    # get NHDPlusFlowlineVAA table into pandas data frame 
    vaa = dbf2DF('%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(NHDdir))[['COMID','HYDROSEQ']]
    vaa.columns = ['VAACOMID','HYDROSEQ']
    #--------------------------------------------------------
    # Merge VAA table on to joined table of catchment, waterbody, and flowline tables
    joinVAA = pd.merge(join,vaa,left_on='FLCOMID',right_on='VAACOMID', how='left')
    # create data frame to link catchment COMID to waterbody COMID
    df = pd.DataFrame(columns=('catCOMID', 'wbCOMID'))
    # iterate through table while grouping by the waterbody COMID to select out associated catchments
    for name, group in joinVAA.groupby('WBAREACOMI'):
        if not pd.isnull(group.FEATUREID).any():
#            print group.FEATUREID
#            break
            hydro = pd.Series.min(group.HYDROSEQ)
            group = group[group.HYDROSEQ == hydro]
            s = pd.DataFrame({'catCOMID':pd.Series([int(group.FEATUREID.values[0])]),'wbCOMID':pd.Series([int(group.WBCOMID.values[0])])})
            df = df.append(s)

###############################################################################
# testing week of 5/9/2016            
joinVAA.ix[joinVAA.FEATUREID == 341745].HYDROSEQ            
join.ix[join.FLCOMID == 24086468]

test = pd.merge(wb, fl, left_on='WBCOMID', right_on='WBAREACOMI', how='inner')
test.ix[test.WBAREACOMI == 24089110]

test2 = pd.merge(wb, flvaa, left_on='WBCOMID', right_on='WBAREACOMI', how='inner')
test2.ix[test2.WBAREACOMI == 24089110]

flvaa = pd.merge(fl, vaa, left_on='FLCOMID', right_on='VAACOMID', how='inner')


wb.ix[wb.WBCOMID == 24086468]
    # convert list values into integers and format to pass into SelectLayerByAttribute function 
    wblist = [int(x) for x in df.wbCOMID.tolist()]            
    offstring = '"COMID" IN (%s)'% wblist 
    
    
for name, group in joinVAA.groupby('WBAREACOMI'):
    if group.WBAREACOMI.any() == 24089110.0:
        print group
   
joinVAA.ix[joinVAA.WBAREACOMI == 24089110]     
     
list(joinVAA.groupby('WBAREACOMI'))

for num in wblist:
    if not num in test.WBCOMID.values:
        print num
joinVAA.ix[~pd.isnull(joinVAA.FEATUREID.values)]
test.ix[pd.isnull(test.WBAREACOMI.values)]
for yip in test2.WBAREACOMI.values:
    if len(test2.ix[test2.WBAREACOMI == yip]) > 1:
        print test2.ix[test2.WBAREACOMI == yip]

test2.ix[~pd.isnull(test2.WBAREACOMI)]
keep = fl.ix[fl.WBAREACOMI != 0]

###############################################################################
# test to get groups with all null FEATUREIDs and report if there are flowlines 
# associated that don't have catchments
import sys
import pandas as pd
sys.path.append('D:/Projects/StreamCat') # change to your streamcat directory
from StreamCat_functions import makeVPUdict, dbf2DF

NHD_dir = 'D:/NHDPlusV21'

inputs = makeVPUdict(NHD_dir)

for zone in inputs:
    hr = inputs[zone]
    # create directories to NHD data
    NHDdir = "D:/NHDPlusV21/NHDPlus%s/NHDPlus%s"%(hr,zone)
    #--------------------------------------------------------
    wb = dbf2DF("%s/NHDSnapshot/Hydrography/NHDWaterbody.dbf"%(NHDdir))[['COMID','FTYPE']]
    wb = wb.loc[wb['FTYPE'].isin(['LakePond','Reservoir'])]
    wb.columns = ['WBCOMID','FTYPE']
    #--------------------------------------------------------
    # get flowline table into pandas data frame 
    fl = dbf2DF("%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(NHDdir))[['COMID','WBAREACOMI']]
    fl.columns = ['FLCOMID','WBAREACOMI']
    #--------------------------------------------------------
    # get catchment table into pandas data frame 
    cat = dbf2DF('%s/NHDPlusCatchment/Catchment.dbf'%(NHDdir))[['FEATUREID']]
    #--------------------------------------------------------
    flowcat = pd.merge(cat,fl,left_on='FEATUREID',right_on='FLCOMID', how='right')# Change how to right from inner!
    join = pd.merge(wb,flowcat,left_on='WBCOMID', right_on='WBAREACOMI', how='left')
    #--------------------------------------------------------
    # get NHDPlusFlowlineVAA table into pandas data frame 
    vaa = dbf2DF('%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(NHDdir))[['COMID','HYDROSEQ']]
    vaa.columns = ['VAACOMID','HYDROSEQ']
    #--------------------------------------------------------
    # Merge VAA table on to joined table of catchment, waterbody, and flowline tables
    joinVAA = pd.merge(join,vaa,left_on='FLCOMID',right_on='VAACOMID', how='left')
    # create data frame to link catchment COMID to waterbody COMID
    df = pd.DataFrame(columns=('catCOMID', 'wbCOMID'))
    checks = {} 
    for name, group in joinVAA.groupby('WBAREACOMI'):
        if pd.isnull(group.FEATUREID).all():
            checks[int(group.WBAREACOMI.values.max())] = group.FLCOMID.tolist()
            
            
            hydro = pd.Series.min(group.HYDROSEQ)
            group = group[group.HYDROSEQ == hydro]
            s = pd.DataFrame({'catCOMID':pd.Series([int(group.FEATUREID.values[0])]),'wbCOMID':pd.Series([int(group.WBCOMID.values[0])])})
            df = df.append(s)

for x in checks:
    print checks[x]
tbl = pd.read_csv('C:/Users/Rdebbout/Desktop/ProblemLakeBasins.csv')
tbl.ix[tbl.VPU == '09']
###############################################################################

# delete the 710 lakes that have a bigger basin than CatAreaSqKM
# this works from the script output of findProblemLakes.py 
# -- had to iterate by hand to remove all lakes that don't create basins that 
# -- display on-network properties
import pandas as pd  
import geopandas as gpd
from geopandas.tools import sjoin

tbl = pd.read_csv('C:/Users/Rdebbout/Desktop/ProblemLakeBasins.csv')    

isoLk = gpd.GeoDataFrame.from_file('D:/Projects/lakesAnalysis/new/Isolated_Albers.shp')

isoLk.head()
len(tbl.WBCOMID.values)
iso2 = isoLk.ix[~isoLk.COMID.isin(tbl.WBCOMID.values)]
iso3 = isoLk.ix[isoLk.COMID.isin(tbl.WBCOMID.values)]
iso2.to_file('D:/Projects/lakesAnalysis/new/Isolated_AlbersOut710.shp', driver="ESRI Shapefile")
iso3.to_file('D:/Projects/lakesAnalysis/710_basin_VIS/710_out_Albers.shp', driver="ESRI Shapefile")

##############################################################################

# messy way of adding VPU,RPU,HR to shapefile, had to restart and lost data
NHD_dir ='D:/NHDPlusV21'
  
B = gpd.GeoDataFrame.from_file('%s/NHDPlusGlobalData/BoundaryUnit.shp' % NHD_dir)
B = B.drop(B.ix[B.DrainageID.isin(['HI','CI'])].index, axis=0) 

iso = gpd.GeoDataFrame.from_file('D:/Projects/lakesAnalysis/710_basin_VIS/710_Albers.shp')
iso = iso[['AREASQKM', 'COMID', 'VPU', 'geometry']]
iso.columns
iso = iso.drop('VPU', axis=1)
B = B.to_crs(iso.crs)

VPU = B.ix[B.UnitType == 'VPU'].copy()
B.columns
type(RPU)
RPU.head()
RPU.crs
iso.crs
join = sjoin(iso, VPU, how="left", op="within")

join = join[['COMID', 'UnitID', 'geometry', 'RPU', 'HR']]
join.columns = ['COMID', 'VPU', 'geometry', 'RPU', 'HR']
join.head()
join.tail()
join.to_file('D:/Projects/lakesAnalysis/710_basin_VIS/Lakes_710_Albers.shp', driver="ESRI Shapefile")
##############################################################################



























