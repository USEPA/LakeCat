# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:22:43 2016

@author: Rdebbout
"""


import os
import re
import sys
import arcpy
import struct
import decimal
import warnings
import datetime
import itertools
import numpy as np
import pysal as ps
import pandas as pd
import rasterio as rs
from osgeo import gdal
from Tkinter import Tk
import geopandas as gpd
from rasterio import features
from arcpy.sa import Watershed
from geopandas.tools import sjoin
from collections import OrderedDict
from tkFileDialog import askdirectory
arcpy.CheckOutExtension("spatial")
warnings.filterwarnings("ignore", category=FutureWarning)
Tk().withdraw()
##############################################################################


#     legacy = 'PROJCS["NAD_1983_Contiguous_USA_Albers",'\
#                'GEOGCS["GCS_North_American_1983",'\
#                  'DATUM["D_North_American_1983",'\
#                    'SPHEROID["GRS_1980",6378137.0,298.257222101]],'\
#                  'PRIMEM["Greenwich",0.0],'\
#                  'UNIT["Degree",0.0174532925199433]],'\
#                'PROJECTION["Albers"],'\
#                'PARAMETER["false_easting",0.0],'\
#                'PARAMETER["false_northing",0.0],'\
#                'PARAMETER["central_meridian",-96.0],'\
#                'PARAMETER["standard_parallel_1",29.5],'\
#                'PARAMETER["standard_parallel_2",45.5],'\
#                'PARAMETER["latitude_of_origin",23.0],'\
#                'UNIT["Meter",1.0]]'
    
fiftyseventy = 'PROJCS["NAD83 / Conus Albers",'\
                  'GEOGCS["NAD83",'\
                    'DATUM["North American Datum 1983",'\
                      'SPHEROID["GRS 1980", 6378137.0, 298.257222101, AUTHORITY["EPSG","7019"]],'\
                      'TOWGS84[1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0],'\
                      'AUTHORITY["EPSG","6269"]],'\
                    'PRIMEM["Greenwich", 0.0, AUTHORITY["EPSG","8901"]],'\
                    'UNIT["degree", 0.017453292519943295],'\
                    'AXIS["Geodetic longitude", EAST],'\
                    'AXIS["Geodetic latitude", NORTH],'\
                    'AUTHORITY["EPSG","4269"]],'\
                  'PROJECTION["Albers Equal Area", AUTHORITY["EPSG","9822"]],'\
                  'PARAMETER["central_meridian", -96.0],'\
                  'PARAMETER["latitude_of_origin", 23.0],'\
                  'PARAMETER["standard_parallel_1", 29.5],'\
                  'PARAMETER["false_easting", 0.0],'\
                  'PARAMETER["false_northing", 0.0],'\
                  'PARAMETER["standard_parallel_2", 45.5],'\
                  'UNIT["m", 1.0],'\
                  'AXIS["Easting", EAST],'\
                  'AXIS["Northing", NORTH],'\
                  'AUTHORITY["EPSG","5070"]]'

def dbfreader(f):
    """Returns an iterator over records in a Xbase DBF file.

    The first row returned contains the field names.
    The second row contains field specs: (type, size, decimal places).
    Subsequent rows contain the data records.
    If a record is marked as deleted, it is skipped.

    File should be opened for binary reads.

    """
    # See DBF format spec at:
    #     http://www.pgts.com.au/download/public/xbase.htm#DBF_STRUCT

    numrec, lenheader = struct.unpack('<xxxxLH22x', f.read(32))    
    numfields = (lenheader - 33) // 32

    fields = []
    for fieldno in xrange(numfields):
        name, typ, size, deci = struct.unpack('<11sc4xBB14x', f.read(32))
        name = name.replace('\0', '')       # eliminate NULs from string   
        fields.append((name, typ, size, deci))
    yield [field[0] for field in fields]
    yield [tuple(field[1:]) for field in fields]

    terminator = f.read(1)
    assert terminator == '\r'

    fields.insert(0, ('DeletionFlag', 'C', 1, 0))
    fmt = ''.join(['%ds' % fieldinfo[2] for fieldinfo in fields])
    fmtsiz = struct.calcsize(fmt)
    for i in xrange(numrec):
        record = struct.unpack(fmt, f.read(fmtsiz))
        if record[0] != ' ':
            continue                        # deleted record
        result = []
        for (name, typ, size, deci), value in itertools.izip(fields, record):
            if name == 'DeletionFlag':
                continue
            if typ == "N":
                value = value.replace('\0', '').lstrip()
                if value == '':
                    value = 0
                elif deci:
                    value = decimal.Decimal(value)
                else:
                    value = int(value)
            elif typ == 'C':
                value = value.rstrip()                                   
            elif typ == 'D':
                try:
                    y, m, d = int(value[:4]), int(value[4:6]), int(value[6:8])
                    value = datetime.date(y, m, d)
                except:
                    value = None
            elif typ == 'L':
                value = (value in 'YyTt' and 'T') or (value in 'NnFf' and 'F') or '?'
            elif typ == 'F':
                value = float(value)
            result.append(value)
        yield result
        
def dbf2DF(f, upper=True):
    data = list(dbfreader(open(f, 'rb')))
    if upper == False:    
        return pd.DataFrame(data[2:], columns=data[0])
    else:
        return pd.DataFrame(data[2:], columns=map(str.upper,data[0]))
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


def findUpstreamNpy(numpy_dir, com):  # Unpacks Numpy files describing the array of upstream COMID's for each catchment in NHD
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

def NHDdict(directory, unit='VPU'):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Creates an OrderdDict for looping through regions of the NHD RPU zones

    Arguments
    ---------
    directory             : the directory contining NHDPlus data at the top level
    unit                  : Vector or Raster processing units 'VPU' or 'RPU'
    '''  
    if unit == 'VPU':
        if not os.path.exists('%s/StreamCat_npy' % directory):
            os.mkdir('%s/StreamCat_npy' % directory)    
        if not os.path.exists('%s/StreamCat_npy/zoneInputs.npy' % directory):
            inputs = makeVPUdict(directory)
        else:
            inputs = np.load('%s/StreamCat_npy/zoneInputs.npy' % directory).item() 
    if unit == 'RPU':
        if not os.path.exists('%s/StreamCat_npy' % directory):
            os.mkdir('%s/StreamCat_npy' % directory)    
        if not os.path.exists('%s/StreamCat_npy/rpuInputs.npy' % directory):
            inputs = makeRPUdict(directory)
        else:
            inputs = np.load('%s/StreamCat_npy/rpuInputs.npy' % directory).item() 
    return inputs

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
    vpuDict = NHDdict('D:/NHDPlusV21')
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

##############################################################################
    
    
def getOnNetLakes(metric, StreamCat, LakeComs):
    vpuDict = NHDdict('D:/NHDPlusV21')
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

def makeRat(fn):
    '''
    __author__ =  "Matt Gregory <matt.gregory@oregonstate.edu >"
    Adds a Raster Attribute Table to the .tif.aux.xml file, then passes those
    values to rat_to_df function to return the RAT in a pandas DataFrame.

    Arguments
    ---------
    fn               : raster filename
    '''
    ds = gdal.Open(fn)
    rb = ds.GetRasterBand(1)
    nd = rb.GetNoDataValue()
    data = rb.ReadAsArray()
    # Get unique values in the band and return counts for COUNT val
    u = np.array(np.unique(data, return_counts=True))
    #  remove NoData value
    u = np.delete(u, np.argwhere(u==nd), axis=1)
    
    # Create and populate the RAT
    rat = gdal.RasterAttributeTable()
    rat.CreateColumn('VALUE', gdal.GFT_Integer, gdal.GFU_Generic)
    rat.CreateColumn('COUNT', gdal.GFT_Integer, gdal.GFU_Generic)
    for i in range(u[0].size):
        rat.SetValueAsInt(i, 0, int(u[0][i]))
        rat.SetValueAsInt(i, 1, int(u[1][i]))
    
    # Associate with the band
    rb.SetDefaultRAT(rat)
    
    # Close the dataset and persist the RAT
    ds = None 
    
    #return the rat to build DataFrame
    df = rat_to_df(rat)
    return df
       
##############################################################################


def rat_to_df(in_rat):
    """
    __author__ =  "Matt Gregory <matt.gregory@oregonstate.edu >"
    Given a GDAL raster attribute table, convert to a pandas DataFrame
    Parameters
    ----------
    in_rat : gdal.RasterAttributeTable
        The input raster attribute table
    Returns
    -------
    df : pd.DataFrame
        The output data frame
    """
    # Read in each column from the RAT and convert it to a series infering
    # data type automatically
    s = [pd.Series(in_rat.ReadAsArray(i), name=in_rat.GetNameOfCol(i))
         for i in xrange(in_rat.GetColumnCount())]

    # Concatenate all series together into a dataframe and return
    return pd.concat(s, axis=1)

##############################################################################

 
def DF2dbf(df, dbf_path, my_specs=None):
    '''
    Convert a pandas.DataFrame into a dbf.

    __author__  = "Dani Arribas-Bel <darribas@asu.edu> "
    ...

    Arguments
    ---------
    df          : DataFrame
                  Pandas dataframe object to be entirely written out to a dbf
    dbf_path    : str
                  Path to the output dbf. It is also returned by the function
    my_specs    : list
                  List with the field_specs to use for each column.
                  Defaults to None and applies the following scheme:
                    * int: ('N', 14, 0)
                    * float: ('N', 14, 14)
                    * str: ('C', 14, 0)
    '''
    if my_specs:
        specs = my_specs
    else:
        type2spec = {int: ('N', 20, 0),
                     np.int64: ('N', 20, 0),
                     float: ('N', 36, 15),
                     np.float64: ('N', 36, 15),
                     str: ('C', 14, 0),
                     np.int32: ('N', 20, 0)
                     }
        types = [type(df[i].iloc[0]) for i in df.columns]
        specs = [type2spec[t] for t in types]
    db = ps.open(dbf_path, 'w')
    db.header = list(df.columns)
    db.field_spec = specs
    for i, row in df.T.iteritems():
        db.write(row)
    db.close()
    return dbf_path
    
##############################################################################
    

def purge(directory, pattern):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Clears directory of created rasters that will need to be re-written due to 
    holding on-network like properties, i.e basins created are larger than the 
    associated catchment.

    Arguments
    ---------
    directory           : directory to be cleared
    pattern             : string value to find in filename for removal
    '''
    for f in os.listdir(directory):
        if re.search(pattern, f):
            os.remove(os.path.join(directory, f))
            
##############################################################################


def updateSinks(wbDF,flDF):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Updates the WBARECOMI field in the NHDFlowline GeoDatFrame where NHDSinks
    intersect with NHDWaterbodies. Not currently held in the NHDPlusV21, but 
    we can process the waterbodies with our on-network approach.

    Arguments
    ---------
    wbDF    : Metric name
    flDF    : Location of intermediate StreamCat csv files
    '''
    flow = flDF.set_index('COMID')
    #chk = wbDF.set_index('COMID_sink')
    sinks = wbDF.ix[wbDF.COMID_sink.notnull()].copy()
    #sinks = bodies.loc[bodies.SOURCEFC.notnull()][['COMID']]
    sinks.rename(columns={'COMID': 'WBAREACOMI'}, inplace=True)
    flow.update(sinks)
    flow.WBAREACOMI = flow.WBAREACOMI.astype(np.int64)
    return flow.reset_index(level=0)

##############################################################################
def rollArray(a, d):	
    if len(d) == 4:
        a = a[0,:]
        new = np.roll(np.roll(a, d[0], axis=0), d[1], axis=1)
        new[d[2],:] = a[d[2],:]
        new[:, d[3]] = a[:, d[3]] 
    if len(d) == 3:
        new = np.roll(a[0,:], d[0], axis=d[1])
        if d[1] == 0:
            new[d[2],:] = a[0,d[2],:]
        if d[1] == 1:
            new[:,d[2]] = a[0,:,d[2]]
    return np.expand_dims(new, axis=0) 

##############################################################################
    

def makeFlows(arr, shiftd, fdr, path, nd):
    iso = np.not_equal(arr, shiftd) * np.not_equal(shiftd, nd)  # cells change value after shift * cells not equal to NoData
    pth = np.equal(fdr,path)  # True when equal to path value
    val = iso * pth * arr 
    shiftval = iso * pth * shiftd
    idx = np.not_equal(val,shiftd)
    tocom = val[idx]
    fromcom = shiftval[idx]
    tocom = tocom[tocom > 0]
    fromcom = fromcom[fromcom > 0]    
    # don't load-in the entire array to the DF, just connection vals
    df = pd.DataFrame({'TOCOMID' : tocom, 
                       'FROMCOMID' : fromcom,
                            'move' : path})
    return df.drop_duplicates(['FROMCOMID','TOCOMID'])

##############################################################################


def compAll(arr, fdr ,moves, from_to, nd):
    for move in moves:
        flow = makeFlows(arr, rollArray(np.copy(arr), moves[move][0]), fdr, moves[move][1], nd)
        from_to = pd.concat([from_to,flow])
    return from_to
    
##############################################################################


def expand(window, size=1):
    r, c = window
    return ((r[0] - size, r[1] + size), (c[0] - size, c[1] + size))

def check_window(window, w, h):
    r, c = window
    return ((max(0, r[0]), min(h, r[1])), (max(0, c[0]), min(w, c[1])))

##############################################################################


def chunk_windows(r, indexes=None, max_ram=250000000):
    if indexes is None:
        indexes = r.indexes
    elif isinstance(indexes, int):
        indexes = [indexes]
    if not indexes:
        raise ValueError('No indexes to read')
    pixel_size = 0
    for bidx in indexes:
        if bidx not in r.indexes:
            raise IndexError('band index out of range')
        idx = r.indexes.index(bidx)
        pixel_size += np.dtype(r.dtypes[idx]).itemsize  
    chunk_size, _ = divmod(max_ram, pixel_size)
    r_h, r_w = r.height, r.width
    if chunk_size >= r_h * r_w:
        yield (0, 0), ((0, r_h), (0, r_w))
    else:
        b_h, b_w = r.block_shapes[0]
        d, _ = divmod(chunk_size, r_w * b_h)
        chunk_height = d * b_h
        d, m = divmod(r_h, chunk_height)
        n_chunks = d + int(m>0)
        for i in range(n_chunks):
            row = i * chunk_height
            height = min(chunk_height, r_h - row)
            yield (i, 0), ((row, row+chunk_height), (0, r_w)) 

##############################################################################


def findFlows(zone_file, fdr_file):
    moves = {'up':[(-1,0,-1),4],'left':[(-1,1,-1),1],'down' :[(1,0,0),64], 
            'right':[(1,1,0),16],'downRight':[(1,1,0,0),32],
            'downLeft':[(1,-1,0,-1), 128],'upRight':[(-1,1,-1,0),8],
            'upLeft':[(-1,-1,-1,-1),2]}
    flows = pd.DataFrame()
    with rs.open(zone_file) as z:
        with rs.open(fdr_file) as f:  # 0 is NoData for fdr
            profile = z.profile.copy()
            nd = profile['nodata']
            assert z.shape == f.shape, "Rasters have different extents!"
            for _, w in chunk_windows(z):  # currently defaults to 250MB
                new_w = check_window(expand(w,2), z.width, z.height)
                data = z.read(window=new_w)
                f_r = f.read(window=new_w)
                flows = pd.concat([flows,compAll(data,f_r,moves,flows,nd)])
    return flows.drop_duplicates(['FROMCOMID','TOCOMID'])
    
##############################################################################


def NHDtblMerge(nhd, bounds, out):
    '''
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Merges all of the NHD tables needed to find on-network lakes. Returns the 
    GeoDataFrames that will be used to find off-network lakes. Attribute fields
    COMID, WBARECOMI, and FEATUREID are used to link waterbodies to catchments.

    Arguments
    ---------
    nhd        : string value of prefix to NHD directory
    unit       : GeoDataFrame of Vector Processing Unit
    
    Returns
    ---------
    wbs        : GeoDataFrame of NHDWaterbodies within the VPU
    '''
    inputs = NHDdict(nhd)
    vpus = bounds.query("UnitType == 'VPU'").copy()
    onNet_connect = {}  # initialize dict to hold assoc. catchments for lakes lowest HYDROSEQ
    Obounds = gpd.GeoDataFrame()
    qa_cols=['Total Waterbodies','On-Network','Off-network','FTYPE_drop',
                                                     'Sink_add','Out_of_bounds']
    qa_tbl = pd.DataFrame()                                                 
    for zone in inputs:
        print zone
        hr = inputs[zone]
        pre = "%s/NHDPlus%s/NHDPlus%s" % (nhd, hr, zone)
        wbShp = gpd.read_file("%s/NHDSnapshot/Hydrography/NHDWaterbody.shp"%(pre))
        ttl_WB = len(wbShp)
        ftype_x = len(wbShp)
        wbShp.columns = wbShp.columns[:-1].str.upper().tolist() + ['geometry'] 
        wbShp = wbShp[['AREASQKM','COMID','FTYPE','geometry']]
        wbShp = wbShp.loc[wbShp['FTYPE'].isin(['LakePond','Reservoir'])]
        ttl_FTYPE = ftype_x - len(wbShp)
        fl = dbf2DF("%s/NHDSnapshot/Hydrography/NHDFlowline.dbf"%(pre))[['COMID', 
                                                                    'WBAREACOMI']]
        cat = gpd.read_file('%s/NHDPlusCatchment/Catchment.shp'%(pre)).drop(
                            ['GRIDCODE', 'SOURCEFC'], axis=1)
        cat.columns = cat.columns[:-1].str.upper().tolist() + ['geometry']                         
        vaa = dbf2DF('%s/NHDPlusAttributes/PlusFlowlineVAA.dbf'%(pre))[['COMID',
                                                                        'HYDROSEQ']]        
        final = pd.merge(cat.drop('geometry', axis=1),fl,left_on='FEATUREID',
                           right_on='COMID',how='inner')
        final = pd.merge(wbShp.drop('geometry',axis=1),final,left_on='COMID',
                           right_on='WBAREACOMI',how='left',
                           suffixes=('_wb','_cat'))
        final = pd.merge(final,vaa,left_on='COMID_cat',
                           right_on='COMID',how='left')
                           
                           
        cols = ['catCOMID','wbCOMID','CatAreaSqKm']
        onNetDF = pd.DataFrame(columns=cols)
        catDict = {} # holds associated lake catchments to an on-network lake
        
        # group by the waterbody COMID to find associated catchment
        for name, group in final.groupby('COMID_wb'):
            if not pd.isnull(group.FEATUREID).any():
                base = group.ix[group.HYDROSEQ.idxmin()]
                row = pd.Series([int(base.COMID_cat), int(base.COMID_wb),
                                 base.AREASQKM_cat], index=cols)
                onNetDF = onNetDF.append(row, ignore_index=True)
                catDict[int(base.COMID_cat)] = group.FEATUREID.astype(int).tolist()
    
        ttl_ON = len(onNetDF)
        # add in related sinks
        sinks = gpd.read_file("%s/NHDPlusBurnComponents/Sink.shp"%(pre))   
        exp = '(SOURCEFC== "NHDWaterbody")&(PURPDESC== "NHDWaterbody closed lake")'
        if len(sinks) > 0:
            sinks = sinks.query(exp)
            try:
                assert len(sinks) > 0
                catSink = sjoin(sinks, cat)
                catSink = catSink[['FEATUREID_right','FEATUREID_left','AREASQKM']]
                catSink.columns = cols
                catSink = catSink.ix[catSink.wbCOMID.isin(wbShp.COMID)]  #this will remove any NHDWaterbody COMIDs that have the wrong FTYPE
                catSink = catSink.ix[~catSink.wbCOMID.isin(onNetDF.wbCOMID)]  #remove any COMIDs that are already in the onNetDF
                ttl_SINK = len(catSink)
                for idx, line in catSink.iterrows():
                    catDict[line.catCOMID] = [line.wbCOMID]
                onNetDF = pd.concat([onNetDF,catSink])
            except AssertionError:
                ttl_SINK = 0 # all sinks got queried out!
                pass
        else:
           ttl_SINK = len(sinks)  # get val if no sinks for QA

        # create numpy arrays for connected catchments to waterbody            
        allLinkd = map(lambda x: catDict[x], catDict.keys())    
        onNet_connect[zone] = {'comids':np.array(catDict.keys()),
                        'lengths':np.array([len(v) for v in allLinkd]),
                        'allLinkd':np.int32(np.hstack(np.array(allLinkd)))}
                                                
        # write-out table of catchment-lake COMID connections
        onNetDF.to_csv("%s/joinTables/join_%s.csv" % (out, zone), index=False)
        offLks = wbShp.ix[~wbShp.COMID.isin(onNetDF.wbCOMID)].copy()
                
        # find off-netowrk lakes that are out-of-bounds
        vpu = vpus.query("UnitID == '%s'" % zone)
        offCen = offLks.copy()
        offCen.geometry = offLks.geometry.centroid
        oob = sjoin(offCen, vpu, op='within')[['AREASQKM','COMID','FTYPE',
                                                'UnitID','geometry']]
        out_of_bounds = offLks.ix[~offLks.COMID.isin(oob.COMID)].copy()
        out_of_bounds['VPU_orig'] = zone
        outCen = offCen.ix[~offCen.COMID.isin(oob.COMID)]
        unit = sjoin(outCen, vpus, op='within')[['COMID','UnitID']]
        out_of_bounds = out_of_bounds.merge(unit, how='left', on='COMID')
        
        Obounds = pd.concat([Obounds,out_of_bounds])
        offLks = offLks.ix[offLks.COMID.isin(oob.COMID)].copy()
        
        ttl_OOB = len(out_of_bounds)        
        ttl_OFF = len(offLks)

        offLks['VPU_orig'] = zone
        vpu_tbl = sjoin(offCen, vpus, op='within')[['COMID','UnitID']]
        offLks = offLks.merge(vpu_tbl, on='COMID', how='left')
        
        offLks.to_file("%s/off_net_%s.shp" % (out, zone))
        qa_tbl[zone] = [ttl_WB,ttl_ON,ttl_OFF,ttl_FTYPE,ttl_SINK,ttl_OOB]
    
    Obounds.to_file("%s/out_of_bounds.shp" % out)
    np.savez_compressed('%s/onNet_LakeCat.npz' % (out), Connect_arrays=onNet_connect)
    qa_tbl.index = qa_cols
    qa_tbl.T.index.rename('VPU', inplace=True)
    qa_tbl['TOTALS'] = qa_tbl.ix[:].sum(axis=1)
    qa_tbl.T.to_csv("%s/Lake_QA.csv" % out)

##############################################################################


def makeBasins (nhd, bounds, out):
    problems = pd.DataFrame()  # holding for overdrawn basin delineations
    allOff = gpd.GeoDataFrame()
    inputs = NHDdict(nhd)
    rasterUnits = NHDdict(nhd, unit='RPU')
    rpus = bounds.query("UnitType == 'RPU'").copy()
    Obounds = gpd.read_file("%s/out_of_bounds.shp" % out)
    arcpy.env.workspace = "%s/rasters/scratchArc" % out
    arcpy.env.outputCoordinateSystem = fiftyseventy
    cols = ['VPU','out_of_raster']
    addOut = pd.DataFrame(columns=cols)
    flow_tbl = pd.DataFrame()
    for zone in inputs:        
        print zone
        hr = inputs[zone]
        pre = "%s/NHDPlus%s/NHDPlus%s" % (nhd, hr, zone)
        addLks = Obounds.ix[Obounds.UnitID == zone].copy()
        offLks = gpd.read_file("%s/off_net_%s.shp" % (out, zone))
        # add back-in lakes that are in other zones 
        offLks = pd.concat([offLks,addLks]).reset_index().drop('index',axis=1)
        offLks.rename(columns={'UnitID':'VPU_moved'}, inplace=True)
        cat = gpd.read_file('%s/NHDPlusCatchment/Catchment.shp'%(pre))
        # make lake and watershed rasters
        ttl_LOST = 0
        for rpu in rasterUnits[zone]:
            #break
            lakes = offLks.copy()
            if len(rasterUnits[zone]) > 1:
                rpuShp = rpus.query("UnitID == '%s'" % rpu).drop(['Hydroseq',
                                                                'UnitType'],
                                                                    axis=1)
                lakes = sjoin(lakes, rpuShp, op='within').drop('index_right',
                                                                axis=1)
                lakes.rename(columns={'UnitID': 'RPU'}, inplace=True)
            if len(rasterUnits[zone]) == 1:
                lakes['RPU'] = rpu
            
            fdr = rs.open("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu))
            if fdr.crs != lakes.crs:
                lakes.to_crs(fdr.crs, inplace=True)
            meta = fdr.meta.copy()
            meta.update(compress='lzw')
            meta.update(nodata=0,
                        dtype=rs.uint32,
                        driver='GTiff',
                        crs={'init': u'epsg:5070'})
            with rs.open("%s/rasters/lakes_%s.tif" % (out, rpu),
                         'w', **meta) as lksRas:
                lksArray = lksRas.read(1)
                shapes = ((g,v) for g,v in zip(lakes.geometry,lakes.COMID))
                burned = features.rasterize(shapes=shapes, fill=0,
                                            out=lksArray,
                                            out_shape=lksArray.shape,
                                            transform=lksRas.transform)
                lksRas.write(burned.astype(np.uint32), indexes=1)    
            outWshed = Watershed("%s/NHDPlusFdrFac%s/fdr" % (pre, rpu),
                                 "%s/rasters/lakes_%s.tif" % (out, rpu),
                                  "VALUE")
            outWshed.save("%s/rasters/wtshds_%s.tif"%(out,rpu))
            
            rat = makeRat("%s/rasters/wtshds_%s.tif"%(out,rpu))
            ttl_LOST += (len(lakes) - len(rat))
            DF2dbf(rat, "%s/rasters/wtshds_%s.tif.vat.dbf"%(out,rpu))
            centroids = lakes.to_crs({'init': u'epsg:4269'}).copy().drop(
                                                            'AREASQKM',axis=1)
            centroids.geometry = centroids.centroid 
            lkCat = sjoin(centroids, cat, op='within', how='left')
            lkCat.columns = lkCat.columns.str.upper()
            
            # add assoc. cats and areas to off-net lakes ----------------------
            lakes = lakes.merge(lkCat[['COMID','FEATUREID','AREASQKM']].rename(
                                columns={'FEATUREID':'catCOMID',
                                         'AREASQKM':'catAREASQKM'}), 
                                on='COMID')
            allOff = pd.concat([allOff,lakes.copy()])
            
            # compare basin sizes ---------------------------------------------
            both = pd.merge(lkCat[['COMID','FEATUREID','AREASQKM']], rat,
                            how='inner', left_on='COMID', right_on='VALUE')
            both['AreaSqKM_basin'] = (both.COUNT * 900) * 1e-6
            bigs = both.ix[both.AREASQKM < both.AreaSqKM_basin].copy()
            bigs['diff'] = abs(bigs.AREASQKM - bigs.AreaSqKM_basin)
            bigs['VPU'] = zone
            bigs['RPU'] = rpu
            problems = pd.concat([problems,bigs])
            flow_rpu = findFlows("%s/rasters/wtshds_%s.tif"%(out,rpu), 
                                 "%s/NHDPlusFdrFac%s/fdr" % (pre, rpu))
        flow_tbl = pd.concat([flow_tbl, flow_rpu])
        row = pd.Series([zone, ttl_LOST], index=cols)
        addOut = addOut.append(row, ignore_index=True)
    # write-out lakes that have a larger watershed 
    # than their containing catchment
    flow_tbl.to_csv("%s/LakeCat_PlusFlow.csv" % out,index=False)
    problems.to_csv("%s/rasters/problems.csv" % out,index=False)
    allOff.to_file("%s/off-network.shp" % out)
    addOut.loc[len(addOut)] = pd.Series(
                                        ['TOTALS', 
                                         addOut['out_of_raster'].sum()],
                                         index=cols)
    qa_tbl = pd.read_csv("%s/Lake_QA.csv" % out)
    qa_tbl = pd.merge(qa_tbl, addOut, on='VPU')
    qa_tbl.to_csv("%s/Lake_QA.csv" % out, index=False)
    purge(out, "off_net_")

##############################################################################


def main (nhd, out):
    
    if not os.path.exists("%s/rasters" % out):
        if not os.path.exists(out):
            os.mkdir(out)
        os.mkdir("%s/rasters" % out)
        os.mkdir("%s/rasters/scratchArc" % out)
        os.mkdir("%s/joinTables" % out)
        os.mkdir("%s/flowTables" % out)
    
    
    NHDbounds = gpd.read_file(
                            "%s/NHDPlusGlobalData/BoundaryUnit.shp" % nhd).drop(
                                            ['AreaSqKM','DrainageID','Shape_Area',
                                             'Shape_Leng','UnitName'], axis=1)
    
    NHDtblMerge(nhd, NHDbounds, out)
    makeBasins(nhd, NHDbounds, out)
    
##############################################################################


if __name__=='__main__':
    if len(sys.argv) == 1: 
        nhd_dir = askdirectory(title='Select the location of NHDPlusv21:',
                                                               initialdir='.')
        q = 'Select the location where you want LakeCat data written:'
        out_dir = askdirectory(title=q, initialdir='.')
    elif len(sys.argv) > 1:
        nhd_dir = sys.argv[1]
        out_dir = sys.argv[2]
#    print nhd_dir
#    print type(nhd_dir)
#    print out_dir
#    print type(out_dir)
    main(nhd_dir, out_dir)
    
##############################################################################

    
#nhd = 'D:/NHDPlusv21'
#out = 'D:/Projects/LakeCat/play_makeBasins'
#bounds = gpd.read_file("%s/NHDPlusGlobalData/BoundaryUnit.shp" % nhd).drop(
#                                ['AreaSqKM','DrainageID','Shape_Area',
#                                 'Shape_Leng','UnitName'], axis=1)
#                                 
#
#filename = askopenfilename(title='Select the file with UID and SITE_ID ',filetypes=[('Shapefile','.shp')],initialdir=os.getcwd())

#def makeFlowTbls(nhd, out):
#    inputs = NHDdict(nhd)
    

