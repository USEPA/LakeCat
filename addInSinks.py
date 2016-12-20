# -*- coding: utf-8 -*-
"""
Created on Fri May 06 10:30:14 2016

this is a sandbox script to work out details for findisolatedLakes.py

@author: Rdebbout
"""
import os
import sys
import pandas as pd
import geopandas as gpd
sys.path.append('D:/Projects/LakeCat') # 'path/to/LakeCat'
from LakeCat_functions import NHDtblMerge, NHDdict

def main (NHDdir, outdir):

    inputs = NHDdict(NHDdir)

    (os.mkdir(outdir),None)[os.path.exists(outdir)]

    boundShp = gpd.read_file(
                "%s/NHDPlusGlobalData/BoundaryUnit.shp" % NHDdir).drop(
                ['AreaSqKM','DrainageID','Shape_Area',
                'Shape_Leng','UnitName'], axis=1)
    vpus = boundShp.query("UnitType == 'VPU'")
    count = 0
    problems = gpd.GeoDataFrame()
    for zone in inputs:
        print zone
        hr = inputs[zone]
        pre = "%s/NHDPlus%s/NHDPlus%s" % (NHDdir, hr, zone)
        vpu = vpus.query("UnitID == '%s'" % zone)
        wbs, cat, allTbls, Xs = NHDtblMerge(pre, vpu)
        onNetDF = pd.DataFrame(columns=('catCOMID','CatAreaSqKm', 'wbCOMID'))
        catCon = {}
        for name, group in allTbls.groupby('COMID_wb'):
            if not pd.isnull(group.FEATUREID).any():
                base = group.ix[group.HYDROSEQ.idxmin()]
                onNetDF = onNetDF.append(pd.Series([int(base.COMID_cat),
                                                    int(base.COMID_wb),
                                                    base.AREASQKM_cat],
                                                    index=['catCOMID',
                                                    'wbCOMID',
                                                    'CatAreaSqKm']),
                                             ignore_index=True)
        print 'OnNet: %s' % str(len(onNetDF))
        count += len(onNetDF)
        problems = pd.concat([problems,Xs.copy()])
        onNetDF.to_csv("%s/join_%s.csv" % (outdir, zone), index=False)
    problems.to_file("%s/outOfBounds.shp" % (outdir))
    print 'Total OnNet lakes: %s' % str(count)

if __name__ == '__main__':
    # run: python addInSinks.py 'path/to/NHD' 'new/path/to/write' 'path/to/LakeCat
    NHDdir = 'D:/NHDPlusV21'
    outdir = 'D:/Projects/LakeCat/sinkTest'
    # reads main('path/to/NHD', 'new/path/to/write')
    main(NHDdir, outdir)
    