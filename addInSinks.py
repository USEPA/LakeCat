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

def main (NHDdir, outdir):

    inputs = NHDdict(NHDdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    boundShp = gpd.read_file(
                "%s/NHDPlusGlobalData/BoundaryUnit.shp" % NHDdir).drop(
                ['AreaSqKM','DrainageID','Shape_Area',
                'Shape_Leng','UnitName'], axis=1)
    vpus = boundShp.query("UnitType == 'VPU'")
    count = 0
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
        onNetDF.to_csv("%s/join_%s.csv" % (outdir, zone), index=False)
    print 'Total OnNet lakes: %s' % str(count)

if __name__ == '__main__':
    # run: python addInSinks.py 'path/to/NHD' 'new/path/to/write' 'path/to/LakeCat
    sys.path.append(sys.argv[3]) # 'path/to/LakeCat
    from LakeCat_functions import NHDtblMerge, NHDdict
    NHDdir = sys.argv[1]
    outdir = sys.argv[2]
    # reads main('path/to/NHD', 'new/path/to/write')
    main(NHDdir, outdir)
