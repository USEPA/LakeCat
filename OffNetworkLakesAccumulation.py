# Author: Marc Weber
# Date: October 30, 2014
# Modified by : Marc Weber
# ArcGIS 10.2.1, Python 2.7
# Purpose: Accumulate off-network lake catchment
#          zonal statistics results to get full
#          lake basin landscape metrics

# Import system modules
import os
from collections import defaultdict
from datetime import datetime
import sys
sys.path.append('D:\\Projects\\lakesAnalysis\\Scripts')
from SpatialPredictionFunctions2 import dbfreader, children, ContinuousAllocation, CategoricalAllocation, CountAllocation

#####################################################################################################################
# Variables
#accum_type = raw_input('Is data Continuous, Count, or Categorical? ')
accum_type = 'Count'
#flow_dir = raw_input('Enter directory for the from-to tables for lakes: ')
flow_dir = 'D:\\Projects\\lakesAnalysis\\From_To_Tables'
#Alloc_dir = raw_input('Enter directory for off-network lake allocations:')
Alloc_dir = 'D:\\Projects\\lakesAnalysis\\Output\\Off-NetworkCatResults'
#Accum_dir = raw_input('Enter directory to output off-netowkr lake accumulation results: ')
Accum_dir = 'D:\\Projects\\lakesAnalysis\\Output\\AccumulationResults'
#landscape_var = raw_input('Enter the allocation file name without extension and without hydroregion name (must be a csv file): ')
landscape_var = 'TRI_cat'
#####################################################################################################################

startTime = datetime.now()

#########################
# run the function(s)
#########################
files = [f for f in os.listdir(Alloc_dir) if f.count(landscape_var) and not f.count('total') and not f.count('Missed')]
flows = [f for f in os.listdir(flow_dir) if f.count('dbf') and not f.count('xml') \
              and not f.count('.vat') and not f.count('lock') and f.count('LkFrm')]
for f in files:
    UpCOMs = defaultdict(list)
    hydroregion = f.split('_')[1].split('.')[0]

    Allocation = "%s/%s"%(Alloc_dir,f)
    Accumulation = "%s/%s"%(Accum_dir,f)
    if not os.path.exists(Accumulation):
        print "working on zone " + hydroregion
        print "Make list from dbf file."
        for k in flows:
            if k.count(hydroregion[:2]):
                flowtable= flow_dir + "/" + k
                infile = open(flowtable, 'rb')
                startfunTime = datetime.now()
                data = list(dbfreader(infile))
                infile.close()
                print "elapsed time " + str(datetime.now()-startfunTime)

                print "Make dictionary from dbf list."
                startfunTime = datetime.now()
                # Note that from-to are reversed in these tables as read in
                for line in data[2:]:
                    FROMCOMID=line[3] 
                    TOCOMID=line[2]
                    if FROMCOMID==0:
                        UpCOMs[TOCOMID] = []
                    else:
                        UpCOMs[TOCOMID].append(FROMCOMID)

                print "elapsed time " + str(datetime.now()-startfunTime)

        print "Make full upstream dictionary from previous dictionary."
        startfunTime = datetime.now()
        Full_COMs = dict()
        for key in UpCOMs.keys():
            Full_COMs[key] = children(key, UpCOMs)
        print "elapsed time " + str(datetime.now()-startfunTime)

        startfunTime = datetime.now()
        print "run allocation"
        if accum_type == 'Continuous':
            ContinuousAllocation(Allocation, Accumulation, Full_COMs)
        elif accum_type == 'Categorical':
            CategoricalAllocation(Allocation, Accumulation, Full_COMs)
        elif accum_type == 'Count':
            CountAllocation(Allocation, Accumulation, Full_COMs)
        print "elapsed time " + str(datetime.now()-startfunTime)
print "total elapsed time for " + landscape_var + ": " + str(datetime.now()-startTime)
