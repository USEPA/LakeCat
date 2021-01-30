"""
Script to build final LakeCat tables.
Date: Jan 22, 2016
Author: Rick Debbout

NOTE: run script from command line w/in the LakeCat directory, the ControlTable
is located at the root of this directory and will be read-in from there.

> python MakeFinalTables_LakeCat.py

"""


import os
import sys
import zipfile
import pandas as pd

from lake_cat_config import FINAL_DIR, OUT_DIR

ctl = pd.read_csv("ControlTable_LakeCat.csv")

assert ctl.FullTableName.duplicated().any() == False, "FullTableName column must be unique"

runners = ctl.query("run == 1").groupby("Final_Table_Name")
tables = runners["FullTableName"].unique().to_dict()
missing = []
for table, metrics in tables.items():
    for metric in metrics:
        accumulated_file = OUT_DIR + "/{}.csv".format(metric)
        if not os.path.exists(accumulated_file):
            missing.append(metric)

if missing:
    for miss in missing:
        print(f"Missing {miss}")
    print("Check output from LakeCat.py")
    sys.exit()

allStats = pd.DataFrame()
for table, metrics in tables.items():

    if not os.path.exists(f"{FINAL_DIR}/{table}.csv"):

        print(f"Running {table}.....")
        combined_metrics = ctl.FullTableName.isin(metrics)
        metric_table = ctl.loc[combined_metrics].copy().reset_index(drop=True)

        for idx, row in metric_table.iterrows():

            a_m = row.AppendMetric if not row.AppendMetric == "none" else ""
            tbl = pd.read_csv(f"{OUT_DIR}/{row.FullTableName}.csv")
            frontCols = [title
                         for title in tbl.columns
                         for x in ["COMID","AreaSqKm","PctFull","inStreamCat"]
                         if x in title and not "Up" in title]
            catArea, catPct, wsArea, wsPct = frontCols[1:5]
            # re-order for correct sequence
            frontCols = [frontCols[i] for i in [0,1,3,2,4,5]]
            summary = row.summaryfield.split(";") if type(row.summaryfield) == str else None

            if row.MetricType == "Mean":
                colname1 = f"{row.MetricName}Cat{a_m}"
                colname2 =  f"{row.MetricName}Ws{a_m}"
                tbl[colname1] = ((tbl[f"CatSum{a_m}"] / tbl[f"CatCount{a_m}"]) * row.Conversion)
                tbl[colname2] = ((tbl[f"WsSum{a_m}"] / tbl[f"WsCount{a_m}"]) * row.Conversion)
                if idx == 0:
                    final = tbl[frontCols + [colname1] + [colname2]]
                else:
                    final = pd.merge(final,tbl[["COMID",colname1,colname2]],on="COMID")
            if row.MetricType == "Density":
                colname1 = f"{row.MetricName}Cat{a_m}"
                colname2 = f"{row.MetricName}Ws{a_m}"
                if summary:
                    finalNameList = []
                    for sname in summary:
                        if "Dens" in  row.MetricName:
                            metricName = row.MetricName[:-4]
                        fnlname1 = f"{metricName}{sname}Cat{a_m}"
                        fnlname2 = f"{metricName}{sname}Ws{a_m}"
                        tbl[fnlname1] = tbl[f"Cat{sname}"] / (tbl[catArea] * (tbl[catPct]/100))
                        tbl[fnlname2] = tbl[f"Ws{sname"] / (tbl[wsArea] * (tbl[wsPct]/100))
                        finalNameList.append(fnlname1)
                        finalNameList.append(fnlname2)
                if table == "RoadStreamCrossings" or table == "CanalsDitches":
                    ## NOTE:  Will there ever be a situation where we will need to use "conversion" here
                    tbl[colname1] = (tbl.CatSum / (tbl.CatAreaSqKm * (tbl.CatPctFull/100)) * row.Conversion)
                    tbl[colname2] = (tbl.WsSum / (tbl.WsAreaSqKm * (tbl.WsPctFull/100)) * row.Conversion)
                else:
                    tbl[colname1] = tbl[f"CatCount{a_m}"] / (tbl[f"CatAreaSqKm{a_m}"] * (tbl[f"CatPctFull{a_m}"]/100))
                    tbl[colname2] = tbl[f"WsCount{a_m}"] / (tbl[f"WsAreaSqKm{a_m}"] * (tbl[f"WsPctFull{a_m}"]/100))
                if idx == 0:
                    if summary:
                        final = tbl[frontCols +
                                    [colname1] +
                                    [x for x in finalNameList if "Cat" in x] +
                                    [colname2] +
                                    [x for x in finalNameList if "Ws" in x]
                                ]
                        final.columns = [x.replace("M3","") for x in final.columns]
                    else:
                        final = tbl[frontCols + [colname1] + [colname2]]
                else:
                    if summary:
                        final = pd.merge(final,tbl[["COMID"] +
                                                    [colname1] +
                                                    [x for x in finalNameList if "Cat" in x] +
                                                    [colname2] +
                                                    [x for x in finalNameList if "Ws" in x]],
                                                    on="COMID")
                        final.columns = [x.replace("M3","") for x in final.columns]
                    else:
                        final = pd.merge(final,tbl[["COMID",colname1,colname2]],on="COMID")
            if row.MetricType == "Percent":
                lookup = pd.read_csv(row.MetricName)
                catcols,wscols = [],[]
                for col in tbl.columns:
                    if "CatVALUE" in col and not "Up" in col:
                        tbl[col] = ((tbl[col] * 1e-6)/(tbl[catArea]*(tbl[catPct]/100))*100)
                        catcols.append(col)
                    if "WsVALUE" in col:
                        tbl[col] = ((tbl[col] * 1e-6)/(tbl[wsArea]*(tbl[wsPct]/100))*100)
                        wscols.append(col)
                if idx == 0:
                    final = tbl[frontCols+catcols + wscols].copy()
                    final.columns = (frontCols +
                                    [f"Pct{x}Cat{a_m}" for x in lookup.final_val.values] +
                                    [f"Pct{y}Ws{a_m}" for y in lookup.final_val.values])
                else:
                    final2 = tbl[["COMID"] + catcols + wscols]
                    final2.columns = (["COMID"] +
                                        [f"Pct{x}Cat{a_m}" for x in lookup.final_val.values] +
                                        [f"Pct{y}Ws{a_m}" for y in lookup.final_val.values])
                    final = pd.merge(final,final2,on="COMID")
                    if table == "AgMidHiSlopes":
                        final = final.drop(["PctUnknown1Cat","PctUnknown2Cat",
                                            "PctUnknown1Ws", "PctUnknown2Ws"],
                                            axis=1)
        statTbl = pd.DataFrame({"ATT":[table]},columns=["ATT","MAX","MIN"])
        if "ForestLossByYear0013" == table:
            final.drop([col for col in final.columns if "NoData" in col], axis=1, inplace=True)
        for c in final.columns.tolist():
            statTbl = pd.concat([statTbl,pd.DataFrame({"ATT": [c],
                                                       "MIN": [final[c].min()],
                                                       "MAX":[final[c].max()]})])
        allStats = pd.concat([allStats,statTbl])
        print(statTbl)
        final = final.set_index("COMID").fillna("NA")
        final = final[final.columns.tolist()[:5] +
                        [x for x in final.columns[5:] if "Cat" in x] +
                        [x for x in final.columns[5:] if "Ws" in x]
                ].fillna("NA")
        out_file = f"{FINAL_DIR}/{table}.csv"
        final.to_csv(out_file)
        # zip up the file....
        zf = zipfile.ZipFile(f"{FINAL_DIR}/zips/{table}.zip", mode="w")
        zf.write(out_file, f"{table}.csv", compress_type=zipfile.ZIP_DEFLATED)
        zf.close()
print("All Done.....")

