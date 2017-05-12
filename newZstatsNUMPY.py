# -*- coding: utf-8 -*-
"""
Created on Fri May 13 14:17:12 2016

@author: Rdebbout
"""

# make full table from all individ. hydro regions

import sys
import os
sys.path.append('D:/Projects/StreamCat')
from StreamCat_functions import makeRPUdict, dbf2DF, bastards, children
import pandas as pd
import numpy as np
from collections import defaultdict

NHD_dir = 'D:/NHDPlusV21'

if not os.path.exists('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir):
    rpuinputs = makeRPUdict(NHD_dir)
else:
    rpuinputs = np.load('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir).item()

count = 0
for zone in rpuinputs:
    for rpu in rpuinputs[zone]:
        print rpu
        cbl = dbf2DF('D:/Projects/lakesAnalysis/From_To_Tables/LkFrmTo_R%s.dbf' % rpu)
        cbl.head(10)
        cbl2 = cbl.ix[:,1:3]
        cbl2.columns = ['TOCOMID','FROMCOMID']
        if count == 0:
            final2 = cbl2.copy()
        else:
            final2 = pd.concat([final2,cbl2])
        count += 1
        
len(set(final2.TOCOMID.values))
len(set(final2.FROMCOMID.values))

final2.to_csv('D:/Projects/lakesAnalysis/From_To_Tables/LakesFlowTableALL.csv', index=False)
##############################################################################
#

flowtable = 'D:/Projects/lakesAnalysis/From_To_Tables/LakesFlowTableAll.csv'

flow = pd.read_csv(flowtable)


# Now table is ready for processing and the UpCOMs dict can be created
fcom = np.array(flow.FROMCOMID)
tcom = np.array(flow.TOCOMID)
UpCOMs = defaultdict(list)
for i in range(0, len(flow), 1):
    FROMCOMID = fcom[i]
    TOCOMID = tcom[i]
    if FROMCOMID == 0:
        UpCOMs[TOCOMID] = []
    else:
        UpCOMs[TOCOMID].append(FROMCOMID)


################################################################################################################
directory = 'D:/Projects/lakesAnalysis'       
os.mkdir(directory + '/LakeCat_npy')     
os.mkdir(directory + '/LakeCat_npy/bastards')
os.mkdir(directory + '/LakeCat_npy/children')

# UpcomDict comes from above, much simpler to make up      
UpStreamComs = UpCOMs
# Load in COMIDs from raster
tbl = dbf2DF('D:\Projects\lakesAnalysis\NHDPlus_Lakes_Basins_Rasters_1/ALL_BASINS.tif.vat.dbf')
COMIDs = tbl.VALUE.values

a = map(lambda x: bastards(x, UpStreamComs), COMIDs)
lengths = np.array([len(v) for v in a])
a = np.int32(np.hstack(np.array(a)))    #Convert to 1d vector

wdb = directory + '/LakeCat_npy/bastards' 
np.save('%s/upStream.npy' % wdb, a)
np.save('%s/comids.npy' % wdb, COMIDs)
np.save('%s/lengths.npy' % wdb, lengths)
del a, lengths

b = map(lambda x: children(x, UpStreamComs), COMIDs)
lengths = np.array([len(v) for v in b])
b = np.int32(np.hstack(np.array(b)))  #Convert to 1d vector
if not os.path.exists(directory + '/children'):
    os.mkdir(directory + '/children')
wdc = directory + '/LakeCat_npy/children'
np.save('%s/upStream.npy' % wdc, b)
np.save('%s/comids.npy' % wdc, COMIDs)
np.save('%s/lengths.npy' % wdc, lengths)
################################################################################################################
# play with arrays
type(lengths)
lengths.max()  # 508

    itemindex = int(np.where(comids == com)[0])
    n = lengths[:itemindex].sum()
    arrlen = lengths[itemindex]
    return upStream[n:n+arrlen]



idx = int(np.where(lengths == 507)[0])  # 133345
COMIDs[idx]  #  com = 12713362
itemindex = int(np.where(COMIDs == com)[0])
com = 12713362
################################################################################################################
# merge w/ NHD table ,  start with area in IsolatedLakes shapefile

tbl = dbf2DF('D:/Projects/lakesAnalysis/new/Isolated_AlbersOut710.dbf')
tbl.head()
chkCOM = tbl.COMID.values
for x in chkCOM:
    if not x in COMIDs:
        print x
        
ALL = dbf2DF('D:\Projects\lakesAnalysis\NHDPlus_Lakes_Basins_Rasters_710/Basins_ALL.tif.vat.dbf')
len(ALL)

import geopandas as gpd

isoLk = gpd.GeoDataFrame.from_file('D:/Projects/lakesAnalysis/new/Isolated_AlbersOut710.shp')
isoLk.columns
iso2 = isoLk.ix[isoLk.COMID.isin(ALL.VALUE.values)]
for x in iso2.COMID.values:
    if x not in COMIDs:
        print x
iso2.ix[pd.unique(iso2.COMID)]
free = []
for x in iso2.COMID.values:
    if x in free:
        print x
    if x not in free:
        free.append(x)

def Accumulation(arr, COMIDs, lengths, upStream, tbl_type):
    '''
    Arguments
    ---------
    arr                   : table containing watershed values
    COMIDs                : numpy array of all zones COMIDs
    lengths               : numpy array with lengths of upstream COMIDs
    upstream              : numpy array of all upstream arrays for each COMID
    tbl_type              : string value of table metrics to be returned
    '''
    coms = np.array(arr.COMID)  #Read in COMIDs
    indices = swapper(coms, upStream)  #Get indices that will be used to map values
    del upStream  # a and indices are big - clean up to minimize RAM
    cols = arr.columns[1:]  #Get column names that will be accumulated
    z = np.zeros(COMIDs.shape)  #Make empty vector for placing values
    outT = np.zeros((len(COMIDs), len(arr.columns)))  #Make empty array for placing final values
    outT[:,0] = COMIDs  #Define first column as comids
    #Loop and accumulate values
    for k in range(0,len(cols)):
        col = cols[k]
        c = np.array(arr[col]) # arr[col].fillna(0) keep out zeros where no data!
        d = c[indices] #Make final vector from desired data (c)
        if 'PctFull' in col:
            area = np.array(arr.ix[:, 1])
            ar = area[indices]
            x = 0
            for i in range(0, len(lengths)):
                # using nan_to_num in average function to treat NA's as zeros when summing
                z[i] = np.ma.average(np.nan_to_num(d[x:x + lengths[i]]), weights=ar[x:x + lengths[i]])
                x = x + lengths[i]
        else:
            x = 0
            for i in range(0, len(lengths)):
                z[i] = np.nansum(d[x:x + lengths[i]])
                x = x + lengths[i]
        outT[:,k+1] = z  #np.nan_to_num() -used to convert to zeros here, now done above in the np.ma.average()
    outT = outT[np.in1d(outT[:,0], coms),:]  #Remove the extra COMIDs
    outDF = pd.DataFrame(outT)
    if tbl_type == 'Ws':
        outDF.columns = np.append('COMID', map(lambda x : x.replace('Cat', 'Ws'),cols.values))
    if tbl_type == 'UpCat':
        outDF.columns = np.append('COMID', 'Up' + cols.values)
    for name in outDF.columns:
        if 'AreaSqKm' in name:
            areaName = name
    outDF.loc[(outDF[areaName] == 0), outDF.columns[2:]] = np.nan  # identifies that there is no area in catchment mask, then NA values across the table
    return outDF
    
    