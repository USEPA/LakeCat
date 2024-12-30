# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 15:29:55 2017

@author: RHill04
@author: mweber
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import time
import os
import sys
import glob
import fiona as fiona
from pyogrio import read_dataframe



def watershed_onnet(focal_com, uids, lengths, up, basins_shp, agg_ws_shp):
    try:
        # focal_uid = uids_trns[np.in1d(comids_trns, focal_com)]
        focal_uid = focal_com
        x = np.ndarray.item(lengths[np.in1d(uids, focal_uid)])
        lngth = lengths[: np.ndarray.item(np.where(np.in1d(uids, focal_uid))[0])]
        start = np.sum(lngth)
        uplist = up[start : start + x]
    except:
        uplist = [focal_uid]
    tmp_basin = basins_shp.loc[basins_shp["FEATUREID"].isin(uplist)]
    local_basin = basins_shp.loc[basins_shp["FEATUREID"]==focal_com]
    tmp_basin = pd.concat([tmp_basin, local_basin], ignore_index=True)
    tmp_basin.is_copy = False
    tmp_basin["COMID"] = int(focal_com)
    tmp_basin = tmp_basin[["geometry", "COMID"]]
    tmp_intervpu = agg_ws_shp.loc[agg_ws_shp["COMID"].isin(uplist)]
    tmp_intervpu.is_copy = False
    tmp_intervpu["COMID"] = int(focal_com)

    tmp_basin = pd.concat([tmp_basin, tmp_intervpu], ignore_index=True)
    # tmp_basin['geometry'] = tmp_basin.buffer(0)
    #    start_time2 = time.time()
    #    tmp_basin = tmp_basin.dissolve(by='COMID')
    #    print("--- %s seconds ---" % (time.time() - start_time2))
    return tmp_basin

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

def WBCOMID_COMID(nhd_dir, hydro, vpu, nhdcats):    
    pre = "%s/NHDPlus%s/NHDPlus%s" % (nhd_dir, hydro, vpu)
    fl = dbf2DF("%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(pre))[['COMID', 'WBAREACOMI']]
    cat = gpd.read_file('%s/NHDPlusCatchment/Catchment.shp'%(pre)).drop(['GRIDCODE', 'SOURCEFC'], axis=1)
    cat.columns = cat.columns[:-1].str.upper().tolist() + ['geometry']                         
    vaa = dbf2DF('%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(pre))[['COMID','HYDROSEQ']]
    # merge 
    df = pd.merge(cat.drop('geometry', axis=1),fl,left_on='FEATUREID',
                       right_on='COMID',how='inner')
    df = df[df['WBAREACOMI']!=0]
    df = pd.merge(df,vaa, on='COMID',how='left')
    df = df.loc[df.groupby('WBAREACOMI')['HYDROSEQ'].idxmin()]
    # initialize containers for on-net lakes                   
    df = df[['WBAREACOMI','COMID']]
    df.columns = ['WBAREACOMID','COMID']
    return df          
                
#Define directories
nhd_dir = 'G:/NHDPlusV21/'
lakecat_dir = 'O:/PRIV/CPHEA/PESD/COR/CORFiles/Geospatial_Library_Projects/LakeCat/'
ws_dir = lakecat_dir + 'Watersheds_Framework/'
np_dir = 'E:/GitProjects/StreamCat/accum_npy/'
comid_wbcomid_lookups = r'O:\PRIV\CPHEA\PESD\COR\CORFILES\Geospatial_Library_Projects\LakeCat\LakeCat_Framework\joinTables' 
out_file = "O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/OnNetLakeWatersheds/OnNetWorkWatersheds.gpkg"
# np_dir = 'O:/PRIV/CPHEA/PESD/COR/CORFiles/Geospatial_Library_Projects/LakeCat/LakeCat_Framework/LakeCat_npy/children/'
#Read in  COMIDs of interest of select certain COMIDs 
#nla17_dir = 'O:/PRIV/CPHEA/PESD/COR/CORFiles/Geospatial_Library_Projects/NLA/NLA2017LandscapeMetrics/'
#lakes_df = pd.read_csv(nla17_dir + 'NLA17_missing.csv')
# nla22_dir = 'E:/WorkingData/NLA22_Watersheds/'
# lakes_df = pd.read_excel(nla22_dir + 'NLA22_SiteListmac.xlsx', 
#                          'NLA22_SiteListmac (sorted)')
# lakes_df = pd.read_csv(nla22_dir + 'missing_lakes.csv')
lakes_df = pd.read_csv(lakecat_dir + 'FTP_Staging/FinalTables/Dams.csv')
# lakes_df = lakes_df[lakes_df['NHDPlusV2 COMID'].notnull()]
# lakes_df = pd.read_csv(nla22_dir + 'missing_lakes.csv')
lakes_df = lakes_df[lakes_df['inStreamCat'] == 1 & lakes_df['COMID'].notnull()]
# lakes_df = lakes_df[lakes_df['COMID'].notnull()]
# lakes_df = lakes_df.rename(columns={'NHDPlusV2 COMID': 'COMID'})
lakes_df = pd.read_csv('E:/WorkingData/missing_onnet.csv')
lakes = np.array(lakes_df.COMID).astype(int)
# coms = np.array(19268286).astype(int)


#Read in on-network numpy files
tmp_np = np.load(ws_dir + 'onNetFramework.npz')
on_uids = tmp_np['uids']
#on_lengths = tmp_np['lengths']
#on_up = tmp_np['upstream']
on_vpus = tmp_np['vpus']
on_uids_trns = tmp_np['on_uids_trns']
on_comids_trns = tmp_np['on_comids_trns']
vectunit = tmp_np['vectunit'].astype(str)
hydreg = tmp_np['hydreg'].astype(str)
del tmp_np
# Create on-net lake watersheds - this creates an array
# of NHDPlus waterbody COMIDS in 'onnet'
onnet = lakes[np.in1d(lakes, on_comids_trns)]
# This next line translates the waterbody COMIDs to NHDPlus
# flowline COMIDS which we need to process StreamCat features
onuids = on_uids_trns[np.in1d(on_comids_trns, onnet)]
onuids = on_uids[np.in1d(on_uids, onuids)]
vpus = on_vpus[np.in1d(on_uids, onuids)].astype(str)

vpu_list = np.unique(vpus).astype(str)

intervpu = gpd.read_file(ws_dir + 'interVPUs.shp')

mask = np.isin(onuids, full['COMID'], invert=True)
missing = onuids[mask]

i=0
for vpu in vpu_list:
    print(vpu)
    hydro = hydreg[np.in1d(vectunit,vpu)]
    onuids_vpu = onuids[np.in1d(vpus, vpu)]
#    onnet_vpu = onnet[np.in1d(vpus, vpu)]
    nhdcats = gpd.read_file(nhd_dir + '/NHDPlus' + hydro[0] + '/NHDPlus' + vpu + '/NHDPlusCatchment/Catchment.shp')
    tmp_np = np.load(np_dir + 'accum_' + str(vpu) + '.npz')
    #nhdcats['dummy'] = 1 
    for lake in onuids_vpu:
        if lake not in done[['COMID']]:
            print(lake)
            start_time2 = time.time()
            if i==0:
                out_ws = watershed_onnet(lake, tmp_np['comids'], tmp_np['lengths'], tmp_np['upstream'], nhdcats, intervpu)
                out_ws = out_ws.to_crs(epsg=5070)
                out_ws['geometry'] = out_ws.buffer(0.01)
                out_ws = out_ws.dissolve(by='COMID')
                out_ws['COMID'] = out_ws.index
                # out_ws['SITE_ID'] = lakes_df.loc[lakes_df.index[i],'SITE_ID']
            else: 
                if not int(lake) in out_ws['COMID'].values:
                    temp_ws = watershed_onnet(lake, tmp_np['comids'], tmp_np['lengths'], tmp_np['upstream'], nhdcats, intervpu)
                    temp_ws = temp_ws.to_crs(epsg=5070)
                    temp_ws['COMID'] = lake
                    temp_ws['geometry'] = temp_ws.buffer(0.01)
                    temp_ws = temp_ws.dissolve(by='COMID')
                    temp_ws['COMID'] = temp_ws.index
                    # temp_ws['SITE_ID'] = lakes_df.loc[lakes_df.index[i],'SITE_ID']
                    out_ws = pd.concat([out_ws, temp_ws], ignore_index=True)
                    # out_ws = out_ws.append(temp_ws, ignore_index=True)
            i+=1
            print("--- %s seconds ---" % (time.time() - start_time2)) 
            #world.to_file(filename=temp_shp,driver='ESRI Shapefile',crs_wkt=prj)
            print('------------------------------------')

# tmp = int(on_comids_trns[np.in1d(on_uids_trns, comid)])        
# out_ws['WBCOMID'] = on_comids_trns[np.in1d(on_uids_trns, out_ws['COMID'])]
# for vals in out_ws['COMID']:
#     tmp = int(on_comids_trns[np.in1d(on_uids_trns, vals)])
#     out_ws.loc[(out_ws['COMID']==vals),'WBCOMID'] = tmp
# out_ws[['WBCOMID']] = out_ws[['WBCOMID']].astype(int)
# out_ws[['COMID']] = out_ws[['COMID']].astype(int)

                    
lookup_files = glob.glob(os.path.join(comid_wbcomid_lookups, "*.csv"))     # advisable to use os.path.join as this makes concatenation OS independent

df_from_each_file = (pd.read_csv(f) for f in lookup_files)
lookup   = pd.concat(df_from_each_file, ignore_index=True)
lookup = lookup.rename(columns={'catCOMID': 'COMID', 'wbCOMID': 'WBCOMID'})
out_ws = pd.merge(out_ws, lookup, how='inner', on='COMID')

out_ws[['WBCOMID']] = out_ws[['WBCOMID']].astype(int)
out_ws[['COMID']] = out_ws[['COMID']].astype(int)
# out_ws = out_ws[['COMID','geometry']]
# out_ws = out_ws[['SITE_ID','COMID','geometry']]

# # Add SITE_ID
lakes_df = lakes_df[['SITE_ID','COMID','LAT_DD83','LON_DD83',
                     'UNIQUE_ID','PSTL_CODE','NES_SITE','Reachcode',
                     'GNIS_ID','GNIS_NAME']]
lakes_df[['COMID']] = lakes_df[['COMID']].astype(int)
# out_ws = out_ws.merge(lakes_df, how='left')
# out_ws = out_ws[['SITE_ID','COMID','CAT_COMID','geometry']]
# out_ws.to_file(nla17_dir + 'Missing_OnNetLakes.shp', driver = 'ESRI Shapefile')
extra = extra[~extra.isin(out_ws['COMID'])]
extra = extra.drop(columns=['CatAreaSqKm_y'])
extra = extra.rename(columns={'CatAreaSqKm_x': 'CatAreaSqKm'})
extra = extra[['COMID', 'CatAreaSqKm', 'WBCOMID', 'geometry']]

full = pd.concat([full, out_ws])
for i, df in enumerate(np.array_split(out_ws, 20)):
    print(f"OnNetWatersheds{i+1}")
    df.info()
    df.to_file(out_file, layer=f"OnNetWatersheds{i+1}",driver="GPKG", mode='w')
out_ws.to_file("O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/OnNetLakeWatersheds/OnNetWatersheds.gpkg", layer='OnNetWatersheds', driver="GPKG", mode='w')

out_ws = out_ws.merge(lakes_df, how='left',on='COMID')
out_ws.to_file(nla22_dir + 'OnNetWatersheds.shp', driver = 'ESRI Shapefile')


#Off-Net Lakes
#Read in off-network numpy files
off_np = np.load('L:/Priv/CORFiles/Geospatial_Library_Projects/LakeCat/Watersheds_Framework/offNetFramework.npz')

#Create off-net lake watersheds
offnet = lakes[np.in1d(lakes, off_np['off_comids_trns'])]
basins = gpd.read_file(ws_dir + '/allBasins.shp')
i=0
startTime = time.time()
for lake in offnet:
    print(lake)
    if i==0:
        out_ws = watershed_offnet(lake, off_np['uids'], off_np['lengths'], off_np['upstream'], off_np['off_uids_trns'], 
                          off_np['off_comids_trns'], basins)
        out_ws = out_ws.to_crs(epsg=5070)
        out_ws['geometry'] = out_ws.buffer(0.01)
        out_ws = out_ws.dissolve(by='COMID')
        out_ws['COMID'] = out_ws.index
        out_ws['SITE_ID'] = lakes_df.loc[lakes_df.index[i],'SITE_ID']
    else:
        temp_ws = watershed_offnet(lake, off_np['uids'], off_np['lengths'], off_np['upstream'], off_np['off_uids_trns'], 
                          off_np['off_comids_trns'], basins)
        temp_ws = temp_ws.to_crs(epsg=5070)
        temp_ws.index.name = None
        temp_ws['COMID'] = lake
        temp_ws['geometry'] = temp_ws.buffer(0.01)
        # temp_ws['COMID'] = temp_ws.index
        temp_ws = temp_ws.dissolve(by='COMID')
        
        temp_ws['SITE_ID'] = lakes_df.loc[lakes_df.index[i],'SITE_ID']
        out_ws = out_ws.append(temp_ws, ignore_index=True)
    i+=1
    print("--- %s seconds ---" % (time.time() - startTime)) 
    #world.to_file(filename=temp_shp,driver='ESRI Shapefile',crs_wkt=prj)
    print('------------------------------------')

# out_ws = out_ws[['COMID','geometry']]  
# out_ws = out_ws.merge(lakes_df, how='left') 
out_ws = out_ws[['SITE_ID','COMID','geometry']]     
out_ws.to_file(nla22_dir + 'OffNetWatersheds.shp', driver = 'ESRI Shapefile')
off_net_inNHD = out_ws

# Combind on and off network lakes with standard columns
# off_net = gpd.read_file('L:/Priv/CORFiles/Geospatial_Library_Projects/NLA/NLA2017LandscapeMetrics/Off_Network_Lakes/Off_Network_NLA17Lakes.shp')
# on_net = gpd.read_file('L:/Priv/CORFiles/Geospatial_Library_Projects/NLA/NLA2017LandscapeMetrics/On_Network_Lakes/OnNetLakes.shp')

off_net = gpd.read_file(nla22_dir + 'OffNetWatersheds.shp')
on_net = gpd.read_file(nla22_dir + 'OnNetWatersheds.shp')
layers = fiona.listlayers('E:/WorkingData/NLA22_Watersheds/NLA22_mac.gdb')
not_in_nhdplus = geodata = gpd.read_file(nla22_dir + 'NLA22_mac.gdb', driver='fileGDB', layer='wspolys1')
 
# dissolve on SITE_ID - extra features in 'not_in_nhdplus'
not_in_nhdplus = not_in_nhdplus.dissolve(by='SITE_ID')

# off_net = off_net[['comid','geometry']]
# off_net.columns = ['COMID','geometry']
# off_net = off_net.merge(lakes_df, how='left')
on_net[['COMID']] = on_net[['COMID']].astype(int)
off_net[['COMID']] = off_net[['COMID']].astype(int)
not_in_nhdplus['COMID'] = np.nan

on_net = on_net[['SITE_ID','COMID','geometry']]
on_net = on_net.dissolve(by='SITE_ID')
off_net = off_net[['SITE_ID','COMID','geometry']]
off_net = off_net.dissolve(by='SITE_ID')
not_in_nhdplus = not_in_nhdplus[['SITE_ID','COMID','geometry']]

lake_wats = pd.concat([on_net, off_net, not_in_nhdplus])

# and then the InNHD off-network
# lake_wats = lake_wats.append(off_net_inNHD, ignore_index=True)
# lake_wats.to_file('L:/Priv/CORFiles/Geospatial_Library_Projects/NLA/NLA2017LandscapeMetrics/NLA17_Watersheds.shp', driver = 'ESRI Shapefile')
out_dir='O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Resource/Physical/WATERSHEDS/NLA2022_Basins/'
lake_wats.to_file(out_dir + "NLA22_Basins.gpkg", layer='lake_wats', driver="GPKG")
