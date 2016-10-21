# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:22:43 2016

@author: Rdebbout
"""

import os
import numpy as np
import pysal as ps
import pandas as pd
import geopandas as gpd
from geopandas.tools import sjoin
from collections import OrderedDict

##############################################################################


def dbf2DF(dbfile, upper=True):
    '''
    __author__ = "Ryan Hill <hill.ryan@epa.gov>"
                 "Marc Weber <weber.marc@epa.gov>"
    Reads and converts a dbf file to a pandas data frame using pysal.

    Arguments
    ---------
    dbfile           : a dbase (.dbf) file
    '''
    db = ps.open(dbfile)
    cols = {col: db.by_col(col) for col in db.header}
    db.close()  #Close dbf 
    pandasDF = pd.DataFrame(cols)
    if upper == True:
        pandasDF.columns = pandasDF.columns.str.upper() 
    return pandasDF
##############################################################################

    
def Accumulation(arr, COMIDs, lengths, upStream, tbl_type):
    '''
    __author__ =  "Marc Weber <weber.marc@epa.gov>" 
                  "Ryan Hill <hill.ryan@epa.gov>"
    Uses the 'Cat' and 'UpCat' columns to caluculate watershed values and returns those values in 'Cat' columns 
	so they can be appended to 'CatResult' tables in other zones before accumulation.

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
        if 'PctFull' in name:
            pct = name
    outDF.loc[(outDF[areaName] == 0), outDF.columns[2:]] = np.nan  # identifies that there is no area in catchment mask, then NA values across the table
    outDF.loc[(outDF[pct] == 0), outDF.columns[2:-1]] = np.nan 
    return outDF
##############################################################################

    
def createAccumTable(table, directory, tbl_type, zone=""):
    '''
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Ryan Hill <hill.ryan@epa.gov>"
    Accesses either children or bastards directory to pass in to Accumulation
    function for either UpCat metrics or Ws metrics.

    Arguments
    ---------
    table                 : table containing catchment values
    directory             : location of numpy files
    zone                  : string of an NHDPlusV2 VPU zone, i.e. 10L, 16, 17
    tbl_type              : string value of table metrics to be returned
    '''
    if tbl_type == 'UpCat':
        directory = directory + '/bastards'
    if tbl_type == 'Ws':
        directory = directory + '/children'
    COMIDs = np.load(directory + '/comids%s.npy' % zone)
    lengths= np.load(directory + '/lengths%s.npy' % zone)
    upStream = np.load(directory + '/upStream%s.npy' % zone)
    add = Accumulation(table, COMIDs, lengths, upStream, tbl_type)
    
    return add
##############################################################################

    
def swapper(coms, upStream):
    '''
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Ryan Hill <hill.ryan@epa.gov>"
    Creates array of indexes for all upstream COMIDs that will be summarized for each local catchment.

    Arguments
    ---------
    coms                  : numpy array of all COMIDs in the zone
    upstream              : numpy array of all upstream COMIDs for each local catchment
    '''
    bsort = np.argsort(coms)
    apos = np.searchsorted(coms[bsort], upStream)
    indices = bsort[apos]
    return indices
##############################################################################
    
    
numpy_dir = 'D:/Projects/lakesAnalysis/npy_710/StreamCat_npy/bastards'

def findUpstreamNpy(com):  # Unpacks Numpy files describing the array of upstream COMID's for each catchment in NHD
    comids = np.load(numpy_dir + '/comids.npy')
    lengths= np.load(numpy_dir + '/lengths.npy')
    upStream = np.load(numpy_dir + '/upStream.npy')
    itemindex = int(np.where(comids == com)[0])
    n = lengths[:itemindex].sum()
    arrlen = lengths[itemindex]
    return upStream[n:n+arrlen]
    
##############################################################################


def makeVPUdict(directory):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Creates an OrderdDict for looping through regions of the NHD to carry InterVPU 
    connections across VPU zones

    Arguments
    ---------
    directory             : the directory contining NHDPlus data at the top level
    '''
    B = dbf2DF('%s/NHDPlusGlobalData/BoundaryUnit.dbf' % directory)
    B = B.drop(B.ix[B.DRAINAGEID.isin(['HI','CI'])].index, axis=0)
    B = B.ix[B.UNITTYPE == 'VPU'].sort_values('HYDROSEQ',ascending=False)
    inputs = OrderedDict()  # inputs = OrderedDict((k, inputs[k]) for k in order)
    for idx, row in B.iterrows():
        inputs[row.UNITID] = row.DRAINAGEID
        print 'HydroRegion (value): ' + row.DRAINAGEID + ' in VPU (key): ' + row.UNITID
    np.save('%s/StreamCat_npy/zoneInputs.npy' % directory, inputs)
    return inputs
##############################################################################


def makeRPUdict(directory):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Creates an OrderdDict for looping through regions of the NHD RPU zones

    Arguments
    ---------
    directory             : the directory contining NHDPlus data at the top level
    '''
    B = dbf2DF('%s/NHDPlusGlobalData/BoundaryUnit.dbf' % directory)
    B = B.drop(B.ix[B.DRAINAGEID.isin(['HI','CI'])].index, axis=0)      
    rpuinputs = OrderedDict()
    for idx, row in B.iterrows():
        if row.UNITTYPE == 'RPU':
            hr = row.DRAINAGEID
            rpu = row.UNITID
            for root, dirs, files in os.walk('%s/NHDPlus%s' % (directory, hr)):
                for name in dirs:
                    if rpu in os.path.join(root, name):
                        zone = os.path.join(root, name).split('\\')[-3].replace('NHDPlus','')
                        break
            if not zone in rpuinputs.keys():
                rpuinputs[zone] = []
            print 'RPU: ' + rpu + ' in zone: ' + zone 
            rpuinputs[zone].append(row.UNITID)
    np.save('%s/StreamCat_npy/rpuInputs.npy' % directory, rpuinputs)
    return rpuinputs
##############################################################################

def loadDict(NHD_dir, rpu=False):
    if rpu == False:
        if not os.path.exists('%s/StreamCat_npy/zoneInputs.npy' % NHD_dir):
            inputs = makeVPUdict(NHD_dir)
        else:
            inputs = np.load('%s/StreamCat_npy/zoneInputs.npy' % NHD_dir).item() 
    if rpu == True:
        if not os.path.exists('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir):
            inputs = makeRPUdict(NHD_dir)
        else:
            inputs = np.load('%s/StreamCat_npy/rpuInputs.npy' % NHD_dir).item() 
    return inputs
##############################################################################
    
    
def getOnNetLakes(metric, StreamCat, LakeComs):
    vpuDict = loadDict('D:/NHDPlusV21')
    for reg in vpuDict:
        tbl = pd.read_csv('%s/LakeCOMs%s.csv' % (LakeComs,reg))
        strmCat = pd.read_csv('%s/%s_%s.csv' % (StreamCat, metric, reg))
        tbl2 = pd.merge(tbl, strmCat.ix[strmCat.COMID.isin(tbl.catCOMID)], left_on='catCOMID', right_on='COMID')
        tbl2 = tbl2.drop(['COMID','catCOMID'], axis=1)
        if metric == 'RdCrs':
            tbl2 = tbl2.drop([x for x in tbl2.columns if 'SlpWtd' in x], axis=1)
        if metric == 'Elev':
            tbl2 = tbl2.drop([x for x in tbl2.columns if 'MAX' in x or 'MIN' in x], axis=1)
        tbl2.rename(columns={'wbCOMID': 'COMID'}, inplace=True)
        if reg == '06':
            final = tbl2.copy()
        else:
            final = pd.concat([final, tbl2])
    return final
##############################################################################    
    
    
def PointInPoly(points, inZoneData, pct_full, summaryfield=None):
    '''
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Rick Debbout <debbout.rick@epa.gov>"
    Returns either the count of spatial points feature in every polygon in a spatial polygons feature or the summary of
    an attribute field for all the points in every polygon of a spatial polygons feature

    Arguments
    ---------
    points        : input points geographic features as a GeoPandas GeoDataFrame
    InZoneData    : input polygon shapefile as a string, i.e. 'C:/Temp/outshape.shp'
    pct_full      : table that links COMIDs to pct_full, determined from catchments that are  not within the US Census border
    summaryfield  : a list of the field/s in points feature to use for getting summary stats in polygons
    '''
    #startTime = dt.now()
    polys = gpd.GeoDataFrame.from_file(inZoneData)#.set_index('FEATUREID')
    points = points.to_crs(polys.crs)
    # Get list of lat/long fields in the table
    latlon = [s for s in points.columns if any(xs in s.upper() for xs in ['LONGIT','LATIT'])]
    # Remove duplicate points for 'Count'
    points2 = points.ix[~points.duplicated(latlon)] # points2.head() polys.head() point_poly_join.head()
    try:
        point_poly_join = sjoin(points2, polys, how="left", op="within") # point_poly_join.ix[point_poly_join.FEATUREID > 1]
        fld = 'GRIDCODE'  #next(str(unicode(x)) for x in polys.columns if x != 'geometry')
    except:
        polys['link'] = np.nan
        point_poly_join = polys #gpd.GeoDataFrame( pd.concat( [points2, polys], ignore_index=True) ) 
        fld = 'link'
    # Create group of all points in catchment
    grouped = point_poly_join.groupby('GRIDCODE')
    point_poly_count = grouped[fld].count() # point_poly_count.head() next((x for x in points2.columns if x != 'geometry'),None)
    # Join Count column on to NHDCatchments table and keep only 'COMID','CatAreaSqKm','CatCount'
    final = polys.join(point_poly_count, on='GRIDCODE', lsuffix='_', how='left')
    final = final[['GRIDCODE_', 'AreaSqKM', fld]].fillna(0)       
    cols = ['GRIDCODE', 'CatAreaSqKm', 'CatCount']
    if not summaryfield == None: # Summarize fields in list with gpd table including duplicates
        point_poly_dups = sjoin(points, polys, how="left", op="within")
        grouped2 = point_poly_dups.groupby('GRIDCODE')
        for x in summaryfield: # Sum the field in summary field list for each catchment
            point_poly_stats = grouped2[x].sum()
            final = final.join(point_poly_stats, on='GRIDCODE_', how='left').fillna(0)
            cols.append('Cat' + x)
    final.columns = cols
    # Merge final table with Pct_Full table based on COMID and fill NA's with 0
    final = pd.merge(final, pct_full, on='GRIDCODE', how='left')
    final['CatPctFull'] = final['CatPctFull'].fillna(100) # final.head() final.ix[final.CatCount == 0]
    #print "elapsed time " + str(dt.now()-startTime)
    for name in final.columns:
        if 'AreaSqKm' in name:
            area = name
    final.loc[(final[area] == 0), final.columns[2:]] = np.nan
    return final
    
##############################################################################

def chkColumnLength(table, LandscapeLayer):
    '''
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Ryan Hill <hill.ryan@epa.gov>"
    Checks the number of columns returned from zonal stats and adds any of the 
    categorical values that that didn't exist within the zone and fills the 
    column with zeros so that all categories will be represented in the table.

    Arguments
    ---------
    table                 : Results table of catchment summarizations
    LandscapeLayer        : string to file holding the table of inter VPU COMIDs
    '''
    # Get ALL categorical values from the dbf associated with the raster to retain all values
    # in the raster in every table, even when a given value doesn't exist in a given hydroregion
    AllCols = dbf2DF(LandscapeLayer + '.vat.dbf').VALUE.tolist()
    col_list = table.columns.tolist()
    col_list.sort()

    col_list.sort(key=len)         # table.columns
    table = table[col_list]
    if len(AllCols) != len(col_list[1:]):
        AllCols = ['VALUE_'+str(x) for x in AllCols]
        diff = list(set(AllCols) - set(col_list[1:]))
        diff.sort()
        diff.sort(key=len)
        for spot in diff:
            here = AllCols.index(spot) + 1
            table.insert(here, spot ,0)
    return table
##############################################################################
    

def getOnNetLakes2(metric, StreamCat, LakeComs, npy_files):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Grabs records from StreamCat for on-network lakes.  Adjusts
    cat results to be an accumulated result of all associated catchments
    to the waterbody comid.

    Arguments
    ---------
    metric           : Metric name
    StreamCat        : Location of intermediate StreamCat csv files
    LakeComs         : Location of csv's that join waterbody COMID to catchment COMID
    npy_files        : Location of files that associate all catchments with WBAREACOMID 
    '''    
    vpuDict = loadDict('D:/NHDPlusV21')
    for reg in vpuDict:
        tbl = pd.read_csv('%s/LakeCOMs%s.csv' % (LakeComs,reg))
        strmCat = pd.read_csv('%s/%s_%s.csv' % (StreamCat, metric, reg))
        if metric == 'RdCrs':
            strmCat = strmCat.drop([x for x in strmCat.columns if 'SlpWtd' in x], axis=1)
        if metric == 'Elev':
            strmCat = strmCat.drop([x for x in strmCat.columns if 'MAX' in x or 'MIN' in x], axis=1)
        cols = [col for col in strmCat.columns if col[:3] =='Cat']
        iso = strmCat[['COMID'] + cols]      
        accumCats = createAccumTable(iso, npy_files, 'Ws', reg)
        accumCats.columns = [col.replace('Ws','Cat') for col in accumCats.columns]        
        upWs = strmCat.ix[strmCat.COMID.isin(tbl.catCOMID)].drop(cols, axis=1)       
        newCats = pd.merge(accumCats, upWs, on="COMID")        
        tbl2 = pd.merge(tbl, newCats, left_on='catCOMID', right_on='COMID')
        tbl2 = tbl2.drop(['COMID','catCOMID'], axis=1)
        tbl2.rename(columns={'wbCOMID': 'COMID'}, inplace=True)
        if reg == '06':
            final = tbl2.copy()
        else:
            final = pd.concat([final, tbl2])
    return final   
