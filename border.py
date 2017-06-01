# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:05:09 2017

@author: Rdebbout
"""

import pandas as pd
import geopandas as gpd
from geopandas.tools import sjoin
from datetime import datetime as dt

inputs = {'10U': 'MS',
          '07' : 'MS',
          '01' : 'NE',
#          '17' : 'PN',
          '15' : 'CO',
          '13' : 'RG',
#          '12' : 'TX',
          '09' : 'SR',
          '02' : 'MA',
          '08' : 'MS',
          '04' : 'GL',
          '03W' : 'SA',
          '03S' : 'SA',
          '03N' : 'SA',
          '18' : 'CA'}

def dissolveStates(f, nm):
    '''
    Arguments
    ---------
    f        : filename of state shapefile
    nm       : name of the column that identifies state names
    '''
    sts = gpd.read_file(f)
    nin = ['United States Virgin Islands',
            'Commonwealth of the Northern Mariana Islands',
            'Guam',
            'Alaska',
            'American Samoa',
            'Puerto Rico',
            'Hawaii']
    sts = sts.drop(sts.ix[sts[nm].isin(nin)].index)
    sts['dissolve'] = 1
    conus = sts.dissolve(by='dissolve')
    conus = conus[[nm,'geometry']]
    conus.ix[conus.index[0]][nm] = 'CONUS'
    return conus

def brdrPctFull(cats, brdr, ncol, lkCat=False, acol='AreaSqKM'):
    '''
    Arguments
    ---------
    cats     : geoDF of basin polygons
    brdr     : geoDF of CONUS polygon
    ncol     : name of the column that uniquely identifies cats polygons    
    lkCat    : bool to allow for the rturn of UID w/ LakeCat output
    acol     : name of column that holds area (sq. KM)
    '''
    # move poly to albers, need to stay in this CRS to cal. area later
    if brdr.crs != cats.crs:
        brdr.to_crs(cats.crs,inplace=True)
    tch = sjoin(cats,brdr,op='within')
    print 'sjoin done!!'
    nwin = cats.ix[~cats[ncol].isin(tch[ncol])].copy()
    if len(nwin) == 0:
        return pd.DataFrame()
    clipd = gpd.overlay(brdr, nwin, how='intersection')
    print 'overlay done!'
    out = clipd.dissolve(by=ncol)
    out['Area_CONUS'] = out.geometry.area * 1e-6    
    out['PctFull'] = (out['Area_CONUS'] / out[acol]) * 100
    if lkCat == True:
        out = out[['UID','PctFull']]
    else:
        out = out[['PctFull']] 
    if not len(out) == len(nwin):
        nwin['PctFull'] = 0
        nwin = nwin.ix[~nwin[ncol].isin(out.index)]
        app = nwin[[ncol] + out.columns.tolist()].set_index(ncol)
        out = pd.concat([out,app])
    return out
        
        
if __name__ == '__main__':
    
    us_file = 'L:/Priv/CORFiles/Geospatial_Library/Data/RESOURCE/POLITICAL/BOUNDARIES/NATIONAL/TIGER_2010_State_Boundaries.shp'
    states = dissolveStates(us_file, 'NAME10')
    # LakeCat
    cats = gpd.read_file('D:/Projects/Frameworxx/shps/allBasins.shp')
    out = brdrPctFull(cats,states,'COMID',True)

    # StreamCat
    nhd = 'D:/NHDPlusV21'
    final = gpd.GeoDataFrame()
    for zone in inputs:
        print zone
        hr = inputs[zone]
        pre = "%s/NHDPlus%s/NHDPlus%s" % (nhd, hr, zone)
        cats = gpd.read_file('%s/NHDPlusCatchment/Catchment.shp'%(pre))
        cats.to_crs({'init':'epsg:5070'},inplace=True)
        start = dt.now()
        out = brdrPctFull(cats,states,'FEATUREID')
        final = pd.concat([final,out])
        print dt.now() - start
        print len(final)
        
joy = pd.DataFrame()
pd.concat([final,joy])


zone = '15'
brdr = states.copy()
ncol = 'FEATUREID'

# Downloaded from https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2016&layergroup=States+%28and+equivalent%29
# 'D:/Projects/LakeCat/Framework/border/tl_2016_us_state.shp'

for com in nwin.FEATUREID.values:
    if not com in out.index:
        print com