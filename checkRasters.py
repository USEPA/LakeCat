# -*- coding: utf-8 -*-
"""
Created on Thu May 05 09:34:59 2016

@author: Rdebbout
"""

import sys
from collections import OrderedDict
sys.path.append('D:/Projects/StreamCat')
from StreamCat_functions import dbf2DF
import numpy as np

inputs = OrderedDict([('10U', ['10e','10f','10g','10h','10i']), ('10L', ['10a','10b','10c','10d']), ('07', ['07a','07b','07c']), ('11', ['11a','11b','11c','11d']), ('06', ['06a']),
                      ('05', ['05a','05b','05c','05d']), ('08', ['08a','08b','03g']), ('01', ['01a']), ('02', ['02a','02b']), ('03N', ['03a','03b']),
                      ('03S', ['03c','03d']), ('03W', ['03e','03f']), ('04', ['04a','04b','04c','04d']), ('09', ['09a']), ('12', ['12a','12b','12c','12d']),
                      ('13', ['13a','13b','13c','13d']), ('14', ['14a','14b']), ('15', ['15a','15b']), ('16', ['16a','16b']), ('17', ['17a','17b','17c','17d']),
                      ('18', ['18a','18b','18c'])])
NLT = np.empty(0)
OLT = 0                      
for zone in inputs:
    print zone
    print '********************************************************'    
    coms = np.empty(0)
    lk1 = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters/reg%s_lakes.tif.vat.dbf' % zone)
    OLT += len(lk1)
    print 'Old Total for Lakes ' + str(zone) + ':  ' + str(len(lk1))    
    nrow = 0
    for rpu in inputs[zone]:
        ws1 = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters/reg%s_wtshds.tif.vat.dbf' % rpu)
        nrow += len(ws1) 
        coms = np.concatenate([coms,ws1.VALUE.values])
        print rpu + ':  ' + str(len(ws1))
    print 'Old Total for ' + str(zone) + ':  ' + str(nrow)
    print 'Set of old coms ' + str(zone) + ':  ' + str(len(set(coms)))
    print '********************************************************' 
    
    nrow2 = 0
    nrow3 = 0
    coms2 = np.empty(0)
    for rpu in inputs[zone]:
        lk2 = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/test_reg%s_lakes.tif.vat.dbf' % rpu)
        ws2 = dbf2DF('D:/Projects/lakesAnalysis/NHDPlus_Lakes_Basins_Rasters_V2/reg%s_wtshds.tif.vat.dbf' % rpu)
        nrow3 += len(lk2)        
        nrow2 += len(ws2)
        NLT = np.concatenate([NLT,lk2.VALUE.values])
        coms2 = np.concatenate([coms2,ws2.VALUE.values])
        print rpu + ' lakes:  ' + str(len(lk2))
        print rpu + ' WS:  ' + str(len(ws2))
    print 'New Lakes Total for ' + str(zone) + ':  ' + str(nrow3)    
    print 'New WS Total for ' + str(zone) + ':  ' + str(nrow2)
    print 'Set of new coms ' + str(zone) + ':  ' + str(len(set(coms2)))
    print '********************************************************'
    #print 'Not in old"
    
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print 'Total old lakes: ' + str(OLT)
print 'Total new lakes: ' + str(len(NLT))


len(NLT)
vector = dbf2DF('D:/Projects/lakesAnalysis/Isolated_Check.dbf')
type(vector.COMID.values)
len(set(vector.COMID.values))

for com in vector.COMID.values:
    if not com in NLT:
        print com
bord = dbf2DF('D:/Projects/lakesAnalysis/Isolated_ON_RPU_BORDER.dbf')
for yo in bord.COMID.values:
    if not yo in NLT:
        print yo
#  ALL COMIDS on border IN lakes rasters!!
        
indices = np.setdiff1d(np.arange(len(vector.COMID.values)), np.unique(vector.COMID.values, return_index=True)[1])
vector.COMID.values[indices]
# duplicated in Isolated_Check.shp
[13118610, 13871500, 18156163]

# list of COMIDs in Isolated_Check.shp but not in lakes rasters
miss = [20166532,14725382,14819137,14820667,22700864,15235574,22323839,167245973,12568346,12967474,22460249,944040052,162424445,22822164,19861298,120052923,23854057,120054025,120054048,24052877,120053926]

len(set(miss))

for x in set(lk1.VALUE.values):
    if not x in coms:
        print x
        
len(set(lk1.VALUE.values))

for x in set(coms):
    if not x in lk1.VALUE.values:
        print x
        
inputs.items()
for v in inputs.itervalues():
    print v
    
    for k in inputs.keys():
        print k