# -*- coding: utf-8 -*-
"""
Created on Fri May 22 14:09:28 2015

@author: rdebbout
"""
import pandas as pd
import pysal as ps

#landscape_layer = ['clay','huden2010','imp2006','om','perm','popden2010','precip','rckdep','roadstrm','runoff', 'sand', 'wtdep']
landscape_layer = ['pestic']

lake_dir = 'D:\\Projects\\lakesAnalysis\\NHDPlus_Lakes_Basins_Rasters\\'
zonal_dir = 'D:\\Projects\\lakesAnalysis\\Output\\Zonal_Stats_Output\\'
out_dir = 'D:\\Projects\\lakesAnalysis\\Output\\Off-NetworkCatResults\\'

lookup= pd.DataFrame({ 
  'rpu':["01a", "02a", "02b", "03a", "03b", "03c", "03d", "03e", "03f", "03g", "04a", "04b","04c", "04d", "05a", "05b", "05c", "05d", "06a", "07a", "07b", "07c", "08a", "08b","09a", "10a", "10b", "10c", "10d", "10e", "10f", "10g", "10h", "10i", "11a", "11b", "11c", "11d", "12a", "12b", "12c", "12d", "13a", "13b", "13c", "13d", "14a", "14b","15a", "15b", "16a", "16b", "17a", "17b", "17c", "17d", "18a", "18b", "18c"],
  'hr':["01", "02", "02", "03N", "03N", "03S", "03S", "03W", "03W", "08", "04", "04","04", "04", "05", "05", "05", "05", "06", "07", "07", "07", "08", "08","09", "10L", "10L", "10L", "10L", "10U", "10U", "10U", "10U", "10U", "11", "11", "11", "11", "12", "12", "12", "12", "13", "13", "13", "13", "14", "14","15", "15", "16", "16", "17", "17", "17", "17", "18", "18", "18"]})

hydro_reg = ["01","02","03S","03N","03W","04","05","06","07","08","09","10L","10U","11","12","13","14","15","16","17","18"]
# dictionary to keep track of comids that aren't covered by each landscape layer, later to be saved as MissedCOMs
uncovered = dict()

for l in landscape_layer:
	# create an empty data frame to concatenate all results tables to 
    emptyDF = {'COMID': [],'CatAreaSqKM': [],'CatMean': [],'CatPctFull': [],'CatCount': [],'CatMin': [],'CatMax': [],'CatRange': [],'CatStd': [],'CatSum': []}
    total_table = pd.DataFrame(emptyDF)    
    for k in hydro_reg:
        # make empty data frames 
		sublist = lookup.rpu[lookup.hr == k]
        data = {'VALUE': []}
        rat = pd.DataFrame(data)
        zonal = pd.DataFrame(data)
		# build tables from lake rasters and zonal stats for each rpu
        for i in sublist:
            db = ps.open(lake_dir + 'reg' + i + '_wtshds.tif.vat.dbf')
            rattemp = pd.DataFrame(db[:])
            rattemp.columns = map(str.upper, db.header)
            db.close() 
            rat = pd.concat([rattemp,rat], ignore_index=True)
            del rattemp
            db = ps.open(zonal_dir + l + '\\zonalstats_' + l + i + '.dbf')
            zonaltemp = pd.DataFrame(db[:,:9])
            zonaltemp.columns = map(str.upper, db.header[:9])
            db.close() 
            zonal = pd.concat([zonaltemp,zonal], ignore_index=True) 
            del zonaltemp
        rat['FULL_AREA'] = rat['COUNT']*900
		# select only rows with the greatest area by each comid
        cleanrat = rat.groupby('VALUE', group_keys=False).apply(lambda x: x.ix[x.FULL_AREA.idxmax()])
        cleanzonal = zonal.groupby('VALUE', group_keys=False).apply(lambda x: x.ix[x.AREA.idxmax()])
        # identify comids that are in the raster but not represented in zonal stats table, uncovered by the landscape raster
		for j in cleanrat.VALUE:
            if not j in cleanzonal.VALUE:
                if not k in uncovered:
                    uncovered[k] = [int(j)]
                else:
                    uncovered[k].append(int(j))
        result = pd.merge(cleanrat, cleanzonal, how = 'left', on='VALUE')
        result['CatPctFull']=(result['AREA']/result['FULL_AREA'])*100
        result['CatAreaSqKM'] = result['FULL_AREA']* 1.0e-6
        result = result.iloc[:,[1,12,6,11,4,7,5,8,9,10]]
        result.columns = ['COMID','CatAreaSqKM','CatMean','CatPctFull','CatCount','CatMin','CatMax','CatRange','CatStd','CatSum']
        result.fillna(0, inplace=True)
        total_table = pd.concat([total_table,result])        
        result.to_csv(out_dir + l + '_' + k + '.csv', index =False)
    missedCOMs = pd.DataFrame.from_dict(uncovered,orient='index')
    missedCOMs = missedCOMs.transpose()
    missedCOMs.to_csv(out_dir + 'MissedCOMs_' + l + '.csv', index=False)
    total_table.to_csv(out_dir + l + '_total.csv', index=False)
    print 'Done with ' + l
    del total_table,result,cleanrat,cleanzonal,rat,zonal
##################################################################################
    
  