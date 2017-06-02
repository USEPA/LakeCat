# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:05:09 2017

@author: Rdebbout
"""

import pandas as pd
import geopandas as gpd
from geopandas.tools import sjoin

# trimmed VPU dict
# only holds VPUs on border
inputs = {'10U': 'MS',
          '07' : 'MS',
          '01' : 'NE',
          '17' : 'PN',
          '15' : 'CO',
          '13' : 'RG',
          '12' : 'TX',
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

def brdrPctFull(zns, brdr, ncol, acol='AreaSqKM'):
    '''
    Arguments
    ---------
    zns     : geoDF of basin polygons
    brdr     : geoDF of CONUS polygon
    ncol     : name of the column that uniquely identifies zns polygons    
    lkCat    : bool to allow for the rturn of UID w/ LakeCat output
    acol     : name of column that holds area (sq. KM)
    '''
    # move poly to albers, need to stay in this CRS to cal. area later
    if brdr.crs != zns.crs:
        brdr.to_crs(zns.crs,inplace=True)
    touch = sjoin(zns,brdr,op='within')
    nwin = zns.ix[~zns[ncol].isin(touch[ncol])].copy()
    if len(nwin) == 0:
        return pd.DataFrame()    
    tot = pd.DataFrame()
    for idx, row in nwin.iterrows():
        p = gpd.GeoDataFrame({ncol: [row[ncol]],
                          acol: [row[acol]]},
                          geometry=[row.geometry],
                          crs=nwin.crs)
        clip = gpd.overlay(brdr, p, how='intersection')
        if len(clip) == 0:
            p['PctFull'] = 0
            tot = pd.concat([tot,p.set_index(ncol)[['PctFull']]])
        else:
            out = clip.dissolve(by=ncol)
            out['Area_CONUS'] = out.geometry.area * 1e-6    
            out['PctFull'] = (out['Area_CONUS'] / out[acol]) * 100
            tot = pd.concat([tot,out[['PctFull']]])
    assert len(tot) == len(nwin)
    return tot

def makeBrdrPctFile(b_file, z_file, b_field, z_field):
    states = dissolveStates(b_file, b_field)
    if z_file[-4:] == '.shp':
        cats = gpd.read_file('D:/Projects/LakeCat_Framework/shps/allBasins.shp')
        final = brdrPctFull(cats,states,'UID')
    else:
        final = pd.DataFrame()
        for zone in inputs:
            hr = inputs[zone]
            pre = "%s/NHDPlus%s/NHDPlus%s" % (z_file, hr, zone)
            cats = gpd.read_file('%s/NHDPlusCatchment/Catchment.shp'%(pre))
            cats.to_crs({'init':'epsg:5070'},inplace=True)
            out = brdrPctFull(cats,states, z_field)
            final = pd.concat([final,out])
        if final.index.names[0] != 'COMID':
            final.index.names = ['COMID']
    return final
        
if __name__ == '__main__':
    
    us_file = 'L:/Priv/CORFiles/Geospatial_Library/Data/RESOURCE/POLITICAL/BOUNDARIES/NATIONAL/TIGER_2010_State_Boundaries.shp'
    lake_basins = 'D:/Projects/LakeCat_Framework/shps/allBasins.shp'
    
    nhd = 'D:/NHDPlusV21'
    csv = makeBrdrPctFile(us_file, lake_basins, 'NAME10', 'UID')
    csv.to_csv('D:/Projects/LakeCat/Framework/border/new2.csv')
    
#z_file = 'D:/Projects/LakeCat_Framework/shps/allBasins.shp'
#ff.to_csv('L:/Priv/CORFiles/Geospatial_Library/Data/Project/StreamCat/ControlTables/ALL_BORDER_CATS.csv')
#
#out.to_csv('D:/Projects/LakeCat/Framework/border/new.csv')

# Downloaded from https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2016&layergroup=States+%28and+equivalent%29
# 'D:/Projects/LakeCat/Framework/border/tl_2016_us_state.shp'

################################################################

# old method, uses overlay on entire GeoDF instead of iterating rows

#def brdrPctFull(zns, brdr, ncol, lkCat=False, acol='AreaSqKM'):
#    '''
#    Arguments
#    ---------
#    zns     : geoDF of basin polygons
#    brdr     : geoDF of CONUS polygon
#    ncol     : name of the column that uniquely identifies zns polygons    
#    lkCat    : bool to allow for the rturn of UID w/ LakeCat output
#    acol     : name of column that holds area (sq. KM)
#    '''
#    # move poly to albers, need to stay in this CRS to cal. area later
#    if brdr.crs != zns.crs:
#        brdr.to_crs(zns.crs,inplace=True)
#    tch = sjoin(zns,brdr,op='within')
#    print 'sjoin done!!'
#    nwin = zns.ix[~zns[ncol].isin(tch[ncol])].copy()
#    if len(nwin) == 0:
#        return pd.DataFrame()
#    clipd = gpd.overlay(brdr, nwin, how='intersection')
#    print 'overlay done!'
#    out = clipd.dissolve(by=ncol)
#    out['Area_CONUS'] = out.geometry.area * 1e-6    
#    out['PctFull'] = (out['Area_CONUS'] / out[acol]) * 100
#    if lkCat == True:
#        out = out[['UID','PctFull']]
#    else:
#        out = out[['PctFull']] 
#    if not len(out) == len(nwin):
#        nwin['PctFull'] = 0
#        nwin = nwin.ix[~nwin[ncol].isin(out.index)]
#        app = nwin[[ncol] + out.columns.tolist()].set_index(ncol)
#        out = pd.concat([out,app])
#    return out