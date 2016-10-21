# -*- coding: utf-8 -*-
"""
Created on Fri May 06 10:30:14 2016

@author: Rdebbout
"""

import sys
import numpy as np
import pandas as pd
sys.path.append('D:/Projects/lakesAnalysis/Scripts')
from LakeCat_functions import dbf2DF, makeVPUdict

NHD_dir = 'D:/NHDPlusV21'                    
# dictionary to iterate thru NHD folder structure
inputs = makeVPUdict(NHD_dir)
nump = 'D:/Projects/lakesAnalysis/On_Net_Npy_files/children'
#--------------------------------------------------------
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
    # the join below is an inner join keeping only flowlines that have a catchment
    flowcat = pd.merge(cat,fl,left_on='FEATUREID',right_on='COMID', how='inner') 
    join = pd.merge(wb,flowcat,left_on='COMID', right_on='WBAREACOMI', how='left')
    #--------------------------------------------------------
    # get NHDPlusFlowlineVAA table into pandas data frame
    vaa = dbf2DF(NHDPlusFlowlinVAADBF) 
    #--------------------------------------------------------
    # Merge VAA table on to joined table of catchment, waterbody, and flowlines
    joinVAA = pd.merge(join,vaa,left_on='COMID_y',right_on='COMID', how='left')
    # iterate through table while grouping by the waterbody COMID to select out associated catchments
    print "Unique WBAREACOMI: " + str(len(pd.unique(joinVAA.ix[joinVAA.WBAREACOMI > 0].WBAREACOMI)))
    print "Length of join tbl > 0: " + str(len(joinVAA.ix[joinVAA.WBAREACOMI > 0]))
    lookup = pd.read_csv("D:/Projects/lakesAnalysis/On_Network_LakeCOMs/LakeCOMs%s.csv" % zone)
    catDict = {}
    for name, group in joinVAA.groupby('WBAREACOMI'):
        if not pd.isnull(group.FEATUREID).any():
            catDict[int(lookup.ix[lookup.wbCOMID == name].catCOMID.values[0])] = group.FEATUREID.astype(int).tolist()
#--------------------------------------------------------        
    a = map(lambda x: catDict[x], catDict.keys())
    lengths = np.array([len(v) for v in a])
    a = np.int32(np.hstack(np.array(a)))
    COMIDs = np.array(catDict.keys())
    print "length of COMIDs: " + str(len(COMIDs))
    print "length of a: " + str(len(a))
    np.save('%s/upStream%s.npy' % (nump, zone), a)
    np.save('%s/comids%s.npy' % (nump, zone), COMIDs)
    np.save('%s/lengths%s.npy' % (nump, zone), lengths)