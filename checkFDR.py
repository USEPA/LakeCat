# -*- coding: utf-8 -*-
"""
Created on Mon May 09 16:32:31 2016

This script was used to quickly loop through all of the fdrs in the NHD and
compare them with the NHD on the L: drive

@author: Rdebbout
"""

import numpy as np
from osgeo import gdal
from collections import OrderedDict

nhdinputs = OrderedDict([('MS10U', ['10e','10f','10g','10h','10i']), ('MS10L', ['10a','10b','10c','10d']), ('MS07', ['07a','07b','07c']), ('MS11', ['11a','11b','11c','11d']), ('MS06', ['06a']),
                      ('MS05', ['05a','05b','05c','05d']), ('MS08', ['08a','08b','03g']), ('NE01', ['01a']), ('MA02', ['02a','02b']), ('SA03N', ['03a','03b']),
                      ('SA03S', ['03c','03d']), ('SA03W', ['03e','03f']), ('GL04', ['04a','04b','04c','04d']), ('SR09', ['09a']), ('TX12', ['12a','12b','12c','12d']),
                      ('RG13', ['13a','13b','13c','13d']), ('CO14', ['14a','14b']), ('CO15', ['15a','15b']), ('GB16', ['16a','16b']), ('PN17', ['17a','17b','17c','17d']),
                      ('CA18', ['18a','18b','18c'])])
for out in nhdinputs:
    hr = out[:2] 
    zone = out[2:] 
    for rpu in nhdinputs[out]:
        print '**************************'
        print rpu
        raster1 = r'H:/NHDPlusV21/NHDPlus%s/NHDPlus%s/NHDPlusFdrFac%s/fdr' % (hr, zone, rpu)
        raster2 = r'L:/Priv/CORFiles/Geospatial_Library/Data/RESOURCE/PHYSICAL/HYDROLOGY/NHDPlusV21/NHDPlus%s/NHDPlus%s/NHDPlusFdrFac%s/fdr' % (hr, zone, rpu)
        
        ds1 = gdal.Open(raster1)
        ds2 = gdal.Open(raster2)
        
        r1 = np.array(ds1.ReadAsArray())
        r2 = np.array(ds2.ReadAsArray())
        
        d = np.array_equal(r1,r2)                 
        if d == False:
            print "They differ"
        
        else:
            print "They are the same"