# -*- coding: utf-8 -*-
"""
Created on Fri May 06 13:09:19 2016

@author: Rdebbout
"""

import sys
from collections import OrderedDict
import pandas as pd
sys.path.append('D:/Projects/StreamCat')
from StreamCat_functions import dbf2DF

inputs = OrderedDict([('10U', ['10e','10f','10g','10h','10i']), ('10L', ['10a','10b','10c','10d']), ('07', ['07a','07b','07c']), ('11', ['11a','11b','11c','11d']), ('06', ['06a']),
                      ('05', ['05a','05b','05c','05d']), ('08', ['08a','08b','03g']), ('01', ['01a']), ('02', ['02a','02b']), ('03N', ['03a','03b']),
                      ('03S', ['03c','03d']), ('03W', ['03e','03f']), ('04', ['04a','04b','04c','04d']), ('09', ['09a']), ('12', ['12a','12b','12c','12d']),
                      ('13', ['13a','13b','13c','13d']), ('14', ['14a','14b']), ('15', ['15a','15b']), ('16', ['16a','16b']), ('17', ['17a','17b','17c','17d']),
                      ('18', ['18a','18b','18c'])])
                      
NHD_dir = 'C:/Users/Rdebbout/temp/NHDPlusV21'
count = 0
for zone in inputs:
    for rpu in inputs[zone]:
        print rpu
        tbl = dbf2DF('D:/Projects/lakesAnalysis/From_To_Tables/LkFrmTo_R%s.dbf' % rpu)
        tbl.head(10)
        tbl2 = tbl.ix[:,1:3]
        tbl2.columns = ['TOCOMID','FROMCOMID']
        if count == 0:
            final = tbl2.copy()
        else:
            final = pd.concat([final,tbl2])
        count += 1
len(set(final.TOCOMID.values))
len(set(final.FROMCOMID.values))






for x in final.groupby(['TOCOMID', 'FROMCOMID']).size():
    if x != 1:
        print x
gb = pd.DataFrame({'TOCOMID' : pd.Series([1,1,1,4],dtype='float32'), 'FROMCOMID' : pd.Series([5,5,7,8],dtype='float32')})
gb.head()
for x in gb.groupby(['TOCOMID', 'FROMCOMID']).size():
    if x != 1:
        print x
        
count = 0
for zone in inputs:
    for rpu in inputs[zone]:
        #print rpu
        cbl = dbf2DF('D:/Projects/lakesAnalysis/From_To_Tables3/LkFrmTo_R%s.dbf' % rpu)
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

final2.to_csv('D:/Projects/lakesAnalysis/From_To_Tables3/LakesFlowTableALL2.csv', index=False)

for x in final2.TOCOMID.values:
    if x not in final.TOCOMID.values:
        print x

for x in final2.FROMCOMID.values:
    if x not in final.FROMCOMID.values:
        print x

for pip in set(final.FROMCOMID.values):
    if pip not in set(final2.FROMCOMID.values):
        print pip

final.ix[final.TOCOMID == 167237348]
final.ix[final.TOCOMID == 5864341]
for nm, grp in final.groupby(['TOCOMID', 'FROMCOMID']):
    print nm
    print grp
    for pu in x[0]:
        print pu
        
# Check for cross-flow between basins
for num in final.FROMCOMID.values:
    for num2 in final.TOCOMID.ix[final.FROMCOMID == num].values:
        if num in final.ix[final.FROMCOMID == num2]:
            print num, num2

dataFile = open('c:/users/Rdebbout/temp/merge_list.txt', 'w')
for zone in inputs:
    for rpu in inputs[zone]:       
        dataFile.write('"D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/reg%s_wtshds.tif"/n' % rpu)
dataFile.close()

merge_string = 'python "C:\Users\mweber\AppData\Local\Continuum\Anaconda\pkgs\gdal-2.0.0-np19py27_0\Scripts\gdal_merge.py" -o c:/users/mweber/temp/testmerge3.tif -of GTiff --n 255 -co "TILED=YES" -co COMPRESS=LZW -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" --optfile c:/users/mweber/temp/merge_list.txt'

gdal_merge.py -o c:/users/Rdebbout/temp/testmerge3.tif -of GTiff  -co "TILED=YES" -co COMPRESS=LZW -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" --optfile c:/users/Rdebbout/temp/merge_list.txt


tbl = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/out4.tif.vat.dbf')
tblm = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/out2.tif.vat.dbf')

tbla = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/reg01a_wtshds.tif.vat.dbf')
tblb = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/reg02a_wtshds.tif.vat.dbf')
chk = np.concatenate([tbla.VALUE.values, tblb.VALUE.values])
outs2 = []
hip = 0
for x in chk:
    if x not in tbl.VALUE.values:
        hip += 1
        print x
        
v = np.concatenate([np.array(outs), np.array(outs2)])
what = []
for x in chk:
    if x not in v:
        what.append(x)


for x in tbla.VALUE.values:
    if x in tblb.VALUE.values:
        print x
type(tblb.VALUE.values)

AllRas = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/AllLakeWSBasn.tif.vat.dbf')

count = 0
chknum = 0
for zone in inputs:
    for rpu in inputs[zone]:
        print rpu
        tblx = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/reg%s_wtshds.tif.vat.dbf' % rpu)
        chknum += len(tblx)
        if count == 0:
            final = tblx.copy()
        if count != 0:
            final = pd.concat([final,tblx])
        count += 1
        
AllRas2 = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/out3.tif.vat.dbf')