# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 12:19:11 2024

@author: mweber
"""

# Import libraries
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import os

# Read in template set of LakeCat COMIDS
alloc_dir = "O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/LakeCat/Allocation_Accumulation"
# Get a list of matching files

LakeCat_template = pd.read_csv(alloc_dir + '/Clay.csv')

# Nutrient file
#nut_dir = 'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/StreamCat/NutrientInventory/Inputs/'
# nut_dir = 'E:/WorkingData/To_Be_Flow_Accumulated/'
# nut = pd.read_csv(nut_dir + 'ClimTerms_2012_10.csv')
#nut_dir = 'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/AmaliaHandler/'
nut_dir = 'O:/PRIV/CPHEA/PESD/COR/CORFILES/Geospatial_Library_Projects/NutrientInventory/CountyLakeResultsData/'
nut = pd.read_parquet(nut_dir + 'n_lwrCountyLakeResults.parquet')
# nut = pd.read_parquet('L:/Public/salford/WetInt/streamcat_wetlandchains.parquet')
cat_area = LakeCat_template[['COMID','CatAreaSqKm']]
cat_area.head()
# add VPU using lookup table
#nut = pd.merge(COMID_VPU, nut, how='left', left_on=['COMID'], right_on=['COMID'])
nut = pd.merge(nut, cat_area, how='left', left_on=['COMID'], right_on=['COMID'])
# nut = nut.drop('Unnamed: 0', axis=1)
# nut = nut.drop('...1', axis=1)

# select columns - this part we can modify to iterate through columns
nut.columns = nut.columns.str.replace('_Cat','')
cols = [i for i in nut.columns if i not in ["COMID", "CatAreaSqKm"]]
# cols = cols[29:31]
for col in cols:
    final = nut[['COMID', col, 'CatAreaSqKm']]
    final = final.rename(columns={col: 'CatSum'})
    final['CatCount'] = 1
    final['CatSum'] = final['CatSum'] * final['CatCount']
    final['CatPctFull'] = 100
    final = final[['COMID', 'CatAreaSqKm', 'CatCount', 'CatSum', 'CatPctFull']]
    table = pa.Table.from_pandas(final)
    pq.write_table(table, nut_dir + 'Allocation_and_Accumulation/' + col + '.parquet')
