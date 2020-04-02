# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 12:08:59 2017

@author: RHill04
"""
import glob, os
import numpy as np
import geopandas as gpd
import pandas as pd
import time
from collections import  OrderedDict

inputs = OrderedDict([('10U','MS'),('10L','MS'),('07','MS'),('11','MS'),('06','MS'),('05','MS'),('08','MS'),\
                      ('01','NE'),('02','MA'),('03N','SA'),('03S','SA'),('03W','SA'),('04','GL'),('09','SR'),\
                      ('12','TX'),('13','RG'),('14','CO'),('15','CO'),('16','GB'),('17','PN'),('18','CA')])

x=0
for reg in inputs:
    hydroregion = inputs[reg]
    print reg
    print hydroregion
    if x == 0:
        vectunit = np.array(reg)
        hydreg = np.array(hydroregion)
    else:
        vectunit = np.append(vectunit, np.array(reg))
        hydreg = np.append(hydreg, np.array(hydroregion))
    x = x + 1


on_wd = 'L:/Priv/CORFiles/Geospatial_Library_Projects/LakeCat/LakeCat_Framework/LakeCat_npy/'
ws_dir = 'L:/Priv/CORFiles/Geospatial_Library_Projects/LakeCat/Watersheds_Framework/'
join_dir = 'L:/Priv/CORFiles/Geospatial_Library_Projects/LakeCat/LakeCat_Framework/joinTables/'

vpu_dir = 'H:/NHDPlusV21/NHDPlusGlobalData/'

start_time = time.time()
vpus = gpd.read_file(vpu_dir + 'VPUs.shp')
print("--- %s seconds ---" % (time.time() - start_time)) 

tester = np.array(vpus.UnitID).astype(str)
tester = tester[0:21]

x = 0
for vpu in tester:
    print vpu
    tmp1 = np.load(r'L:\Priv\CORFiles\Geospatial_Library_Projects\LakeCat\LakeCat_Framework\LakeCat_npy\onNet_LakeCat.npz')['vpus'].item()[vpu]
    if x == 0:
        uids = tmp1['comids']
        lengths = tmp1['lengths']
        upstream = tmp1['upstream']
        v = np.repeat(vpu, len(tmp1['comids']))
    else:
        uids = np.append(uids, tmp1['comids'])
        lengths = np.append(lengths, tmp1['lengths'])
        upstream = np.append(upstream, tmp1['upstream'])
        v = np.append(v, np.repeat(vpu, len(tmp1['comids'])))
    x = x + 1

off_uids = np.load(ws_dir + 'accum.npz')['comids']
off_lengths = np.load(ws_dir + 'accum.npz')['lengths']
off_up = np.load(ws_dir + 'accum.npz')['upstream']

files = os.listdir(join_dir)
x = 0
for join in files:
    print join
    if x == 0:
        tmp = pd.read_csv(join_dir + join)
        print len(tmp)
    else:
        tmp2 = pd.read_csv(join_dir + join)
        print len(tmp2)
        tmp = tmp.append(tmp2)
    x = x + 1

#tmp.loc[tmp['wbCOMID'].isin([937110111])]
on_comids_trns = np.array(tmp.wbCOMID)
on_uids_trns = np.array(tmp.catCOMID)

basins = gpd.read_file(ws_dir + '/allBasins.shp')
off_uids_trns = np.array(basins.UID)
off_comids_trns = np.array(basins.COMID)

#np.savez_compressed(ws_dir + 'ID_translations.npz', 
#                    on_comids_trns=on_comids_trns,
#                    on_uids_trns=on_uids_trns,
#                    off_uids_trns=off_uids_trns,
#                    off_comids_trns=off_comids_trns)

np.savez_compressed(ws_dir + 'onNetFramework.npz', 
                    uids=uids,
                    lengths=lengths,
                    upstream=upstream,
                    vpus=v,
                    on_comids_trns=on_comids_trns,
                    on_uids_trns=on_uids_trns,
                    vectunit=vectunit,
                    hydreg=hydreg)

np.savez_compressed(ws_dir + 'offNetFramework.npz', 
                    uids=off_uids,
                    lengths=off_lengths,
                    upstream=off_up,
                    off_uids_trns=off_uids_trns,
                    off_comids_trns=off_comids_trns)



wd = 'L:/Priv/CORFiles/Geospatial_Library/Data/Project/StreamCat/InterVPUtable/'
files = os.listdir(wd)
files = filter(lambda k: '.dbf' in k, files)
x = 0
for shp in files:
    print shp
    shp = shp[:-len('.dbf')]
    tmp_shp = gpd.read_file(wd + shp + '.shp')
    tmp_shp['COMID'] = int(shp)
    tmp_shp = tmp_shp[['COMID', 'geometry']]
    if x == 0:
        out_shp = tmp_shp
    else:
        out_shp = out_shp.append(tmp_shp)
    x = x+1
    
out_shp.to_file(filename = ws_dir + 'interVPUs.shp', driver = 'ESRI Shapefile')


















