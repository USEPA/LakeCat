# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 15:29:55 2017

@author: RHill04
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import time
import os
import arcpy
from shapely.ops import unary_union

def watershed_offnet(focal_com, uids, lengths, up, uids_trns, comids_trns, basins_shp):
    try:   
        focal_uid = uids_trns[np.in1d(comids_trns, focal_com)]
        l = np.asscalar(lengths[np.in1d(uids, focal_uid)])
        start = np.sum(lengths[:np.asscalar(np.where(np.in1d(uids, focal_uid))[0])])
        uplist = up[start:start+l]
    except: 
        uplist = [focal_uid]     
    tmp_basin = basins_shp.loc[basins_shp['UID'].isin(uplist)]  
    tmp_basin.is_copy = False
    tmp_basin['COMID'] = int(focal_com)
    tmp_basin = tmp_basin[['geometry', 'COMID']]
    tmp_basin = tmp_basin.dissolve(by='COMID')
    return tmp_basin

def watershed_onnet(focal_com, uids, lengths, up, basins_shp, agg_ws_shp):        
    try: 
        #focal_uid = uids_trns[np.in1d(comids_trns, focal_com)]
        focal_uid = focal_com
        l = np.asscalar(lengths[np.in1d(uids, focal_uid)])
        start = np.sum(lengths[:np.asscalar(np.where(np.in1d(uids, focal_uid))[0])])
        uplist = up[start:start+l]
    except: 
        uplist = [focal_uid]
    tmp_basin = basins_shp.loc[basins_shp['FEATUREID'].isin(uplist)] 
    tmp_basin.is_copy = False
    tmp_basin['COMID'] = int(focal_com)
    tmp_basin = tmp_basin[['geometry', 'COMID']]
    tmp_intervpu = agg_ws_shp.loc[agg_ws_shp['COMID'].isin(uplist)]
    tmp_intervpu.is_copy = False
    tmp_intervpu['COMID'] = int(focal_com)
    tmp_basin = tmp_basin.append(tmp_intervpu)
    #tmp_basin['geometry'] = tmp_basin.buffer(0)
#    start_time2 = time.time()
#    tmp_basin = tmp_basin.dissolve(by='COMID')
#    print("--- %s seconds ---" % (time.time() - start_time2))  
    return tmp_basin


    #Define directories
nhd_dir = 'D:/GISData/NHDPlusV21/'
ws_dir = 'L:/Priv/CORFiles/Geospatial_Library/Data/Project/LakeCat/Watersheds_Framework/'
hab_dir = 'L:/Priv/CORFiles/Geospatial_Library/Data/Project/LakeCat/Applications/HarmfulAlgalBlooms/'
np_dir = 'L:/Priv/CORFiles/Geospatial_Library/Data/Project/StreamCat/numpy_files/children/'
    #Read in lakes COMIDs (398 lakes)
lakes = pd.read_csv(hab_dir + 'lakeshed_COMID_all.csv')
lakes = np.array(lakes.COMID).astype(int)
    #Read in off-network numpy files
off_np = np.load(ws_dir + 'offNetFramework.npz')
    #Read in off-network basins shapefile (big)
basins = gpd.read_file(ws_dir + '/allBasins.shp')
    #Create off-net lake watersheds
offnet = lakes[np.in1d(lakes, off_np['off_comids_trns'])]
for lake in offnet:
    print lake
    out_ws = watershed_offnet(lake, off_np['uids'], off_np['lengths'], 
                              off_np['upstream'], off_np['off_uids_trns'], 
                              off_np['off_comids_trns'], basins)
    out_ws['COMID'] = lake
    out_ws.crs = basins.crs
    out_ws.to_file(filename = hab_dir + '/out_shps/ws_' + str(lake) + '.shp', driver = 'ESRI Shapefile')
    #world.to_file(filename=temp_shp,driver='ESRI Shapefile',crs_wkt=prj)
    print '------------------------------------'



    #Read in on-network numpy files
tmp_np = np.load(ws_dir + 'onNetFramework.npz')
on_uids = tmp_np['uids']
#on_lengths = tmp_np['lengths']
#on_up = tmp_np['upstream']
on_vpus = tmp_np['vpus']
on_uids_trns = tmp_np['on_uids_trns']
on_comids_trns = tmp_np['on_comids_trns']
vectunit = tmp_np['vectunit']
hydreg = tmp_np['hydreg']
del tmp_np
    #Create on-net lake watersheds
onnet = lakes[np.in1d(lakes, on_comids_trns)]
onuids = on_uids_trns[np.in1d(on_comids_trns, onnet)]
onuids = on_uids[np.in1d(on_uids, onuids)]
vpus = on_vpus[np.in1d(on_uids, onuids)]
vpu_list = np.unique(vpus)
intervpu = gpd.read_file(ws_dir + 'interVPUs.shp')

for vpu in vpu_list:
    print vpu
    hydro = hydreg[np.in1d(vectunit,vpu)]
    onuids_vpu = onuids[np.in1d(vpus, vpu)]
    nhdcats = gpd.read_file(nhd_dir + '/NHDPlus' + hydro[0] + '/NHDPlus' + vpu + '/NHDPlusCatchment/Catchment.shp')
    tmp_np = np.load(np_dir + 'accum_' + vpu + '.npz')
    #nhdcats['dummy'] = 1 
    for lake in onuids_vpu:
        print lake
        start_time2 = time.time()
        out_ws = watershed_onnet(lake, tmp_np['comids'], tmp_np['lengths'], tmp_np['upstream'], nhdcats, intervpu)
        out_ws.crs = nhdcats.crs
        #out_ws = out_ws.to_crs(basins.crs)
        #out_ws[out_ws['geometry'].notnull()].to_file(filename = hab_dir + '/out_shps/tmp_shps/ws_' + str(int(lake)) + '.shp', driver = 'ESRI Shapefile')    
        if not os.path.exists(hab_dir + '/out_shps/tmp_shps/ws_' + str(int(lake)) + '.shp'):  
            out_ws.to_file(filename = hab_dir + '/out_shps/tmp_shps/ws_' + str(int(lake)) + '.shp', driver = 'ESRI Shapefile')
        print("--- %s seconds ---" % (time.time() - start_time2)) 
        #world.to_file(filename=temp_shp,driver='ESRI Shapefile',crs_wkt=prj)
        print '------------------------------------'

files = os.listdir(hab_dir + '/out_shps/tmp_shps')
files = filter(lambda k: '.dbf' in k, files)


for shp in files:
    shp = shp[:-len('.dbf')]    
    print shp
    inshp = hab_dir + 'out_shps/tmp_shps/' + shp + '.shp'
    shp = shp.replace('-','neg')
    outshp = hab_dir + 'out_shps/' + shp + '.shp'   
    if not arcpy.Exists(hab_dir + '/out_shps/' + shp + '.shp'):
        arcpy.Dissolve_management(inshp, outshp, "COMID", "", "MULTI_PART", "DISSOLVE_LINES") 

#----------------------
files = os.listdir(hab_dir + '/out_shps')
files = filter(lambda k: '.dbf' in k, files)

for shp in files:
    shp = shp[:-len('.dbf')]    
    print shp
    inshp = hab_dir + 'out_shps/' + shp + '.shp'
    ws = gpd.read_file(inshp)
    if ws.crs == {}:
        #print 'No projection'
        ws.crs = nhdcats.crs
        ws = ws.to_crs(basins.crs)
        ws.to_file(filename = inshp, driver = 'ESRI Shapefile')


#----------------------
#Need to correct on uids/comids mixup
#Resulting watersheds from algorithm above have catCOMIDs instead of
#waterbody COMIDs. Most lists of lakes will be based on wbCOMIDs
off_uids_trns = off_np['off_uids_trns']
off_comids_trns = off_np['off_comids_trns']

files = os.listdir(hab_dir + '/out_shps')
files = filter(lambda k: '.dbf' in k, files)

for shp in files:
    shp = shp[:-len('.dbf')]    
    print shp
    inshp = hab_dir + 'out_shps/' + shp + '.shp'
    ws = gpd.read_file(inshp)
    comid = ws.COMID[0]
    if np.sum(np.in1d(on_uids_trns, comid)) == 1:
        tmp = int(on_comids_trns[np.in1d(on_uids_trns, comid)])
        ws['COMID'] = tmp
    else:
        tmp = comid
    outws = hab_dir + 'out_shps/final/ws_' + str(tmp) + '.shp'
    ws.to_file(filename = outws, driver = 'ESRI Shapefile')
    print '----------------------'
    
#    print np.sum(np.in1d(on_comids_trns, comid))
#    print np.sum(np.in1d(on_uids_trns, comid))
#    print np.sum(np.in1d(off_comids_trns, comid))
#    print np.sum(np.in1d(off_uids_trns, comid))

    




#focal_com = lake
#
#focal_com = 19251225
#lake = focal_com
#
#start_time2 = time.time()
#nhdcats = gpd.read_file(nhd_dir + '/NHDPlus' + 'CO' + '/NHDPlus' + '15' + '/NHDPlusCatchment/Catchment.shp')
#print("--- %s seconds ---" % (time.time() - start_time2))
#
#np.sum(np.in1d(on_uids_trns, focal_com))
#np.sum(np.in1d(off_comids_trns, focal_com))
#
#basins_shp = nhdcats
#uids = tmp['comids']
#lengths = tmp['lengths']
#up = tmp['upstream']
#uids_trns = on_uids_trns
#comids_trns = on_comids_trns
#agg_ws_shp = intervpu
#
#np.sum(np.in1d(uids, 7221400))
#
#tmp_basin.crs = nhdcats.crs
#tmp_basin = tmp_basin.to_crs(basins.crs)
#lake = focal_com
#tmp_basin.to_file(filename = hab_dir + '/out_shps/ws_' + str(int(lake)) + '.shp', driver = 'ESRI Shapefile')
#
#
#tester = intervpu.loc[intervpu['COMID'].isin(uplist)]
#tester.crs = nhdcats.crs
#tester = tester.to_crs(basins.crs)
#tester['geometry'] = tester.simplify(15)
#tester.to_file(filename = 'D:/Scratch/tmp7.shp', driver = 'ESRI Shapefile')
#
#from shapely.ops import unary_union
#tester = gpd.read_file('D:/Scratch/dissolvie.shp')
#a = tester.ix[0].geometry
#b = tester.ix[1].geometry
#c = unary_union([a,b])
#d = gpd.GeoDataFrame({'COMID': [22071698]}, geometry = [c], crs = tester.crs)
#d = d.to_crs(basins.crs)
#d.to_file(filename = 'D:/Scratch/dissolvie3.shp')



basins_shp = basins
uids = off_uids
lengths = off_lengths
up = off_up
uids_trns = off_uids_trns
comids_trns = off_comids_trns
