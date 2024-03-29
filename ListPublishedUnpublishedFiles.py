import pandas as pd
import os
from ftplib import FTP


ftp = FTP('newftp.epa.gov')
ftp.login()
ftp.cwd("/EPADataCommons/ORD/NHDPlusLandscapeAttributes/LakeCat/FinalTables")


published = ftp.nlst()

local_list = os.listdir(r'O:\PRIV\CPHEA\PESD\COR\CORFiles\Geospatial_Library_Projects\LakeCat\FTP_Staging\FinalTables\zips')

local_published =  [x for x in local_list if x in published]
local_unpublished =  [x for x  in local_list if x not in published]

control = pd.read_csv(r'E:\GitProjects\LakeCat\ControlTable_LakeCat.csv')            
local_published_metrics =  [elem.split('_Region')[0] for elem in local_published]
local_published_metrics = list(set(local_published_metrics))

local_published_metrics = [os.path.splitext(x)[0] for x in local_published_metrics]

list(control)
control['Published'] = control['Final_Table_Name'].isin(local_published_metrics)
control['Published'] = control['Published'].map({True: 'Yes', False: 'No'})
control.to_csv(r'E:\GitProjects\LakeCat\ControlTable_LakeCat.csv', index=False)
