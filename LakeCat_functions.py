# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:22:43 2016

@author: Rdebbout
"""


import os
import re
import sys
import warnings
from collections import OrderedDict, defaultdict, deque
from datetime import datetime as dt
from itertools import chain

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rs
from geopandas.tools import sjoin
from osgeo import gdal

os.environ["PATH"] = r"{};{}".format(
    os.environ["PATH"], r"C:\Program Files\ArcGIS\Pro\bin"
)
sys.path.append(r"C:\Program Files\ArcGIS\Pro\Resources\ArcPy")
import arcpy
from arcpy import PolygonToRaster_conversion as p2r
from arcpy import RasterToPolygon_conversion as r2p
from arcpy.sa import TabulateArea, Watershed, ZonalStatisticsAsTable

from lake_cat_config import LYR_DIR, NHD_DIR, OUT_DIR, STREAMCAT_DIR

warnings.filterwarnings("ignore", category=RuntimeWarning)

inputs = OrderedDict(
    [
        ("06", "MS"),
        ("05", "MS"),
        ("10U", "MS"),
        ("10L", "MS"),
        ("07", "MS"),
        ("11", "MS"),
        ("14", "CO"),
        ("01", "NE"),
        ("17", "PN"),
        ("16", "GB"),
        ("15", "CO"),
        ("13", "RG"),
        ("12", "TX"),
        ("09", "SR"),
        ("02", "MA"),
        ("08", "MS"),
        ("04", "GL"),
        ("03W", "SA"),
        ("03S", "SA"),
        ("03N", "SA"),
        ("18", "CA"),
    ]
)

rpus = OrderedDict(
    [
        ("01", ["01a"]),
        ("02", ["02a", "02b"]),
        ("03N", ["03a", "03b"]),
        ("03S", ["03c", "03d"]),
        ("03W", ["03e", "03f"]),
        ("04", ["04a", "04b", "04c", "04d"]),
        ("05", ["05a", "05b", "05c", "05d"]),
        ("06", ["06a"]),
        ("07", ["07a", "07b", "07c"]),
        ("08", ["03g", "08a", "08b"]),
        ("09", ["09a"]),
        ("10L", ["10a", "10b", "10c", "10d"]),
        ("10U", ["10e", "10f", "10g", "10h", "10i"]),
        ("11", ["11a", "11b", "11c", "11d"]),
        ("12", ["12a", "12b", "12c", "12d"]),
        ("13", ["13a", "13b", "13c", "13d"]),
        ("14", ["14a", "14b"]),
        ("15", ["15a", "15b"]),
        ("16", ["16a", "16b"]),
        ("17", ["17a", "17b", "17c", "17d"]),
        ("18", ["18a", "18b", "18c"]),
    ]
)

# this is the same projection as the 2 below, but this matches the name of the proj
# from the NHD fdr files..

fiftyseventy = "PROJCS['NAD_1983_Albers',\
              GEOGCS['GCS_North_American_1983',\
                DATUM['D_North_American_1983',\
                SPHEROID['GRS_1980',6378137.0,298.257222101]],\
              PRIMEM['Greenwich',0.0],\
                UNIT['Degree',0.0174532925199433]],\
              PROJECTION['Albers'],PARAMETER['False_Easting',0.0],\
              PARAMETER['False_Northing',0.0],\
              PARAMETER['Central_Meridian',-96.0],\
              PARAMETER['Standard_Parallel_1',29.5],\
              PARAMETER['Standard_Parallel_2',45.5],\
              PARAMETER['Latitude_Of_Origin',23.0],\
              UNIT['Meter',1.0]],\
              VERTCS['Unknown VCS',VDATUM['Unknown'],\
              PARAMETER['Vertical_Shift',0.0],\
              PARAMETER['Direction',1.0],\
              UNIT['User_Defined_Unit',0.01]]"


def dbf2DF(f, upper=True):
    data = gpd.read_file(f).drop("geometry", axis=1)
    if upper is True:
        data.columns = data.columns.str.upper()
    return data


def doStats(OUT_DIR, LYR_DIR, NHD_DIR):
    arcpy.CheckOutExtension("spatial")
    arcpy.env.cellSize = "30"

    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)
    if not os.path.exists(f"{OUT_DIR}/ZStats"):
        os.mkdir(f"{OUT_DIR}/ZStats")

    ctl = pd.read_csv(r"ControlTable_LakeCat.csv")

    for _, row in ctl.query("run == 1").iterrows():

        print(f"running....{row.FullTableName}")

        LLyr = f"{LYR_DIR}/{row.LandscapeLayer}"
        ddir = f"{OUT_DIR}/ZStats/{row.FullTableName}"
        if not os.path.exists(ddir) and row.accum_type != "Point":
            os.mkdir(ddir)
        summaryfield = None
        if type(row.summaryfield) == str:
            summaryfield = row.summaryfield.split(";")
        start = dt.now()

        if row.accum_type != "Point":
            csv = f"{OUT_DIR}/ZStats/{row.FullTableName}.csv"
            stats = pd.DataFrame()
            for zone, hr in inputs.items():
                pre = f"{NHD_DIR}/NHDPlus{hr}/NHDPlus{zone}"
                for rpu in rpus[zone]:
                    if row.MetricName == "Elev":
                        LLyr = f"{pre}/NEDSnapshot/Ned{rpu}/{row.LandscapeLayer}"
                    out = f"{OUT_DIR}/ZStats/{row.FullTableName}/{row.FullTableName}_{rpu}.dbf"
                    if not os.path.exists(out):
                        raster = f"framework/rasters/wsheds/wtshds_{rpu}.tif"
                        if row.accum_type == "Categorical":
                            TabulateArea(raster, "Value", LLyr, "Value", out, "30")
                        if row.accum_type == "Continuous":
                            ZonalStatisticsAsTable(
                                raster, "Value", LLyr, out, "DATA", "ALL"
                            )
                    tbl = dbf2DF(out)
                    tbl.rename(columns={"VALUE": "UID"}, inplace=True)
                    stats = pd.concat([stats, tbl])
            stats.UID = stats.UID.astype(np.int64)
            stats.to_csv(csv, index=False)

        if row.accum_type == "Point":

            pct_full = pd.read_csv("framework/border/pct_full.csv")
            points = gpd.GeoDataFrame.from_file(LLyr)
            basins = "framework/shps/allBasins.shp"
            stats = PointInPoly(points, basins, pct_full, "UID", summaryfield)

        print("ZonalStats Results Complete in : " + str(dt.now() - start))

        if row.accum_type != "Point":
            b = pd.DataFrame()
            for zone in rpus.keys():
                for rpu in rpus[zone]:
                    b_ = dbf2DF(f"framework/rasters/wsheds/wtshds_{rpu}.tif.vat.dbf")
                    b_["BSNAREASQKM"] = (b_.COUNT * 900) * 1e-6
                    b_ = b_[["VALUE", "BSNAREASQKM", "COUNT"]]
                    b_.columns = ["UID", "AreaSqKm", "COUNT"]
                    b = pd.concat([b, b_])

        if row.accum_type == "Categorical":
            stats = chkColumnLength(stats, LLyr)
            cols = stats.columns.tolist()[1:]
            stats["AREA"] = stats[stats.columns.tolist()[1:]].sum(axis=1)
            stats = pd.merge(b, stats, how="left", on="UID")
            stats["PctFull"] = ((stats.AREA * 1e-6) / stats.AreaSqKm) * 100
            stats = stats[["UID", "AreaSqKm"] + cols + ["PctFull"]]
            cols = stats.columns[1:]
            stats.columns = np.append("UID", "Cat" + cols.values)
            stats = stats.fillna(0)

        if row.accum_type == "Continuous":
            stats = pd.merge(b, stats, how="left", on="UID")
            stats["CatPctFull"] = (stats.COUNT_y / stats.COUNT_x) * 100
            if row.FullTableName == "Elev":
                stats = stats[
                    ["UID", "AreaSqKm", "COUNT_x", "SUM", "MAX", "MIN", "CatPctFull"]
                ]
                stats.columns = [
                    "UID",
                    "CatAreaSqKm",
                    "CatCount",
                    "CatSum",
                    "CatMax",
                    "CatMin",
                    "CatPctFull",
                ]
            else:
                stats = stats[["UID", "AreaSqKm", "COUNT_x", "SUM", "CatPctFull"]]
                stats.columns = [
                    "UID",
                    "CatAreaSqKm",
                    "CatCount",
                    "CatSum",
                    "CatPctFull",
                ]
            stats.CatPctFull = stats.CatPctFull.fillna(0)
        start2 = dt.now()
        npy = "framework/LakeCat_npy"
        accum = np.load(f"{npy}/bastards/accum.npz")
        up = Accumulation(
            stats, accum["comids"], accum["lengths"], accum["upstream"], "UpCat", "UID"
        )
        accum = np.load(f"{npy}/children/accum.npz")
        ws = Accumulation(
            stats, accum["comids"], accum["lengths"], accum["upstream"], "Ws", "UID"
        )
        stats = pd.merge(stats, up, on="UID")
        stats = pd.merge(stats, ws, on="UID")
        cols = stats.columns[1:].tolist()
        # goto StreamCat to get On-Net-work lake results from assoc. COMIDs
        stats["inStreamCat"] = 0
        # Join UID to COMID for final deliverable

        lks = dbf2DF("framework/off-network.dbf")[["COMID", "UID"]]

        off = pd.merge(lks, stats, on="UID", how="right")
        off.drop("UID", axis=1, inplace=True)
        on = getOnNetLakes(
            row.FullTableName,
            STREAMCAT_DIR,
            "framework/joinTables",
            f"{npy}/onNet_LakeCat.npz",
            NHD_DIR,
        )
        on["inStreamCat"] = 1
        print("Length of on_Net: " + str(len(on)))
        tot = pd.concat([off, on])
        tot.to_csv(f"{OUT_DIR}/{row.FullTableName}.csv", index=False)
        print("Accumulation Results Complete in : " + str(dt.now() - start2))


def Accumulation(arr, COMIDs, lengths, upStream, tbl_type, icol):
    """
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Ryan Hill <hill.ryan@epa.gov>"
    Uses the 'Cat' and 'UpCat' columns to caluculate watershed values and returns those values in 'Cat' columns
        so they can be appended to 'CatResult' tables in other zones before accumulation.

    Arguments
    ---------
    arr                   : table containing watershed values
    COMIDs                : numpy array of all zones COMIDs
    lengths               : numpy array with lengths of upstream COMIDs
    upstream              : numpy array of all upstream arrays for each COMID
    tbl_type              : string value of table metrics to be returned
    """
    coms = np.array(arr[icol])  # Read in COMIDs
    indices = swapper(coms, upStream)  # Get indices that will be used to map values
    del upStream  # a and indices are big - clean up to minimize RAM
    cols = arr.columns[1:]  # Get column names that will be accumulated
    z = np.zeros(COMIDs.shape)  # Make empty vector for placing values
    outT = np.zeros(
        (len(COMIDs), len(arr.columns))
    )  # Make empty array for placing final values
    outT[:, 0] = COMIDs  # Define first column as comids
    # Loop and accumulate values
    for k in range(0, len(cols)):
        col = cols[k]
        c = np.array(arr[col])  # arr[col].fillna(0) keep out zeros where no data!
        d = c[indices]  # Make final vector from desired data (c)
        if "PctFull" in col:
            area = np.array(arr.iloc[:, 1])
            ar = area[indices]
            x = 0
            for i in range(0, len(lengths)):
                # using nan_to_num in average function to treat NA's as zeros when summing
                z[i] = np.ma.average(
                    np.nan_to_num(d[x : x + lengths[i]]), weights=ar[x : x + lengths[i]]
                )
                x = x + lengths[i]
        else:
            x = 0
            for i in range(0, len(lengths)):
                z[i] = np.nansum(d[x : x + lengths[i]])
                x = x + lengths[i]
        outT[:, k + 1] = z
    outT = outT[np.in1d(outT[:, 0], coms), :]  # Remove the extra COMIDs
    outDF = pd.DataFrame(outT)
    if tbl_type == "Ws":
        outDF.columns = np.append(
            icol, list(map(lambda x: x.replace("Cat", "Ws"), cols.values))
        )
    elif tbl_type == "UpCat":
        outDF.columns = np.append(icol, "Up" + cols.values)
    else:
        outDF.columns = [icol] + cols.tolist()
    for name in outDF.columns:
        if "AreaSqKm" in name:
            areaName = name
        if "PctFull" in name:
            pct = name
    outDF.loc[
        (outDF[areaName] == 0), outDF.columns[2:]
    ] = (
        np.nan
    )  # identifies that there is no area in catchment mask, then NA values across the table
    outDF.loc[(outDF[pct] == 0), outDF.columns[2:-1]] = np.nan
    return outDF


def children(token, tree, chkset=None):
    """
    __author__ = "Ryan Hill <hill.ryan@epa.gov>"
                 "Marc Weber <weber.marc@epa.gov>"
    returns a list of every child

    Arguments
    ---------
    token           : a single COMID
    tree            : Full dictionary of list of upstream COMIDs for each COMID in the zone
    chkset          : set of all the NHD catchment COMIDs used to remove flowlines with no associated catchment
    """
    visited = set()
    to_crawl = deque([token])
    while to_crawl:
        current = to_crawl.popleft()
        if current in visited:
            continue
        visited.add(current)
        node_children = set(tree[current])
        to_crawl.extendleft(node_children - visited)
    # visited.remove(token)
    if chkset != None:
        visited = visited.intersection(chkset)
    return list(visited)


def bastards(token, tree, chkset=None):
    """
    __author__ = "Ryan Hill <hill.ryan@epa.gov>"
                 "Marc Weber <weber.marc@epa.gov>"
    returns a list of every child w/ out father (key) included

    Arguments
    ---------
    token           : a single COMID
    tree            : Full dictionary of list of upstream COMIDs for each COMID in the zone
    chkset          : set of all the NHD catchment COMIDs, used to remove flowlines with no associated catchment
    """
    visited = set()
    to_crawl = deque([token])
    while to_crawl:
        current = to_crawl.popleft()
        if current in visited:
            continue
        visited.add(current)
        node_children = set(tree[current])
        to_crawl.extendleft(node_children - visited)
    visited.remove(token)
    if chkset != None:
        visited = visited.intersection(chkset)
    return list(visited)


def swapper(coms, upStream):
    """
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Ryan Hill <hill.ryan@epa.gov>"
    Creates array of indexes for all upstream COMIDs that will be summarized for each local catchment.

    Arguments
    ---------
    coms                  : numpy array of all COMIDs in the zone
    upstream              : numpy array of all upstream COMIDs for each local catchment
    """
    bsort = np.argsort(coms)
    apos = np.searchsorted(coms[bsort], upStream)
    indices = bsort[apos]
    return indices


def findUpstreamNpy(numpy_dir, com):
    """Unpacks Numpy files describing the array of upstream COMID's for
    each catchment in NHD. Similar to the Arc add-in tool that Marc made to
    identify upstream flowline/catchments. NOT USED IN LAKECAT PROCESS! for QA.

    Arguments
    ---------
    numpy_dir             : Loccation of numpy files
    com                   : ID of feature
    """
    comids = np.load(numpy_dir + "/comids.npy")
    lengths = np.load(numpy_dir + "/lengths.npy")
    upStream = np.load(numpy_dir + "/upStream.npy")
    itemindex = int(np.where(comids == com)[0])
    n = lengths[:itemindex].sum()
    arrlen = lengths[itemindex]
    return upStream[n : n + arrlen]


def PointInPoly(points, inZoneData, pct_full, fld="GRIDCODE", summaryfield=None):
    """
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Rick Debbout <debbout.rick@epa.gov>"
    Returns either the count of spatial points feature in every polygon in a spatial polygons feature or the summary of
    an attribute field for all the points in every polygon of a spatial polygons feature

    Arguments
    ---------
    points        : input points geographic features as a GeoPandas GeoDataFrame
    InZoneData    : input polygon shapefile as a string, i.e. 'C:/Temp/outshape.shp'
    pct_full      : table that links COMIDs to pct_full, determined from catchments that are  not within the US Census border
    fld           : the field in the InZoneData file that uniquely identifies each polygon
    summaryfield  : a list of the field/s in points feature to use for getting summary stats in polygons
    """
    polys = gpd.GeoDataFrame.from_file(inZoneData)
    if not points.crs == polys.crs:
        points = points.to_crs(polys.crs)
    # Get list of lat/long fields in the table
    latlon = [
        s for s in points.columns if any(xs in s.upper() for xs in ["LONGIT", "LATIT"])
    ]
    # Remove duplicate points for 'Count'
    points2 = points.loc[~points.duplicated(latlon)]
    try:
        point_poly_join = sjoin(points2, polys, how="left", op="within")
    except:
        polys["link"] = np.nan
        point_poly_join = polys
        fld = "link"
    # Create group of all points in catchment
    grouped = point_poly_join.groupby(fld)
    point_poly_count = grouped[fld].count()
    point_poly_count.name = "COUNT"
    # Join Count column on to NHDCatchments table and keep only 'COMID','CatAreaSqKm','CatCount'
    final = polys.join(point_poly_count, on=fld, how="left")
    final = final[[fld, "AreaSqKM", "COUNT"]].fillna(0)
    cols = [fld, "CatAreaSqKm", "CatCount"]
    if (
        not summaryfield == None
    ):  # Summarize fields in list with gpd table including duplicates
        point_poly_dups = sjoin(points, polys, how="left", op="within")
        grouped2 = point_poly_dups.groupby(fld)
        for x in summaryfield:  # Sum the field in summary field list for each catchment
            point_poly_stats = grouped2[x].sum()
            final = final.join(point_poly_stats, on=fld, how="left").fillna(0)
            cols.append("Cat" + x)
    final.columns = cols
    # Merge final table with Pct_Full table based on COMID and fill NA's with 0
    final = pd.merge(final, pct_full, on=fld, how="left")
    final["CatPctFull"] = final["CatPctFull"].fillna(100)
    for name in final.columns:
        if "AreaSqKm" in name:
            area = name
    # replace CatAreaSqKm with NANs where value is zero
    final.loc[(final[area] == 0), final.columns[2:]] = np.nan
    return final


def chkColumnLength(table, LandscapeLayer):
    """
    __author__ =  "Marc Weber <weber.marc@epa.gov>"
                  "Ryan Hill <hill.ryan@epa.gov>"
    Checks the number of columns returned from zonal stats and adds any of the
    categorical values that that didn't exist within the zone and fills the
    column with zeros so that all categories will be represented in the table.

    Arguments
    ---------
    table                 : Results table of catchment summarizations
    LandscapeLayer        : string to file holding the table of inter VPU COMIDs
    """
    # Get ALL categorical values from the dbf associated with the raster to retain all values
    # in the raster in every table, even when a given value doesn't exist in a given hydroregion
    AllCols = dbf2DF(LandscapeLayer + ".vat.dbf").VALUE.tolist()
    col_list = table.columns.tolist()
    col_list.sort()

    col_list.sort(key=len)  # table.columns
    table = table[col_list]
    if len(AllCols) != len(col_list[1:]):
        AllCols = ["VALUE_" + str(x) for x in AllCols]
        diff = list(set(AllCols) - set(col_list[1:]))
        diff.sort()
        diff.sort(key=len)
        for spot in diff:
            here = AllCols.index(spot) + 1
            table.insert(here, spot, 0)
    return table


def getOnNetLakes(metric, StreamCat, LakeComs, npy_files, nhd):
    """
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Grabs records from StreamCat for on-network lakes.  Adjusts
    cat results to be an accumulated result of all associated catchments
    to the waterbody comid.

    Arguments
    ---------
    metric           : Metric name
    StreamCat        : Location of intermediate StreamCat csv files
    LakeComs         : Location of csv's that join waterbody COMID to catchment COMID
    npy_files        : Location of files that associate all catchments with WBAREACOMID
    """
    final = pd.DataFrame()
    for zone in inputs:
        tbl = pd.read_csv("%s/join_%s.csv" % (LakeComs, zone))[
            ["catCOMID", "wbCOMID"]
        ]  # remove
        strmCat = pd.read_csv("%s/%s_%s.csv" % (StreamCat, metric, zone))
        if metric == "RdCrs":
            strmCat = strmCat.drop(
                [x for x in strmCat.columns if "SlpWtd" in x], axis=1
            )
        #        if metric == 'Elev':
        #            strmCat = strmCat.drop([x for x in strmCat.columns if 'MAX' in x or 'MIN' in x], axis=1)
        cols = [col for col in strmCat.columns if col[:3] == "Cat"]
        iso = strmCat[["COMID"] + cols]
        accum = np.load(npy_files, allow_pickle=True)["vpus"].item()[zone]
        accumCats = Accumulation(
            iso, accum["comids"], accum["lengths"], accum["upstream"], "", "COMID"
        )
        #        # shouldn't be needed if keep tbl_type arg as empty string in Accumulation
        #        accumCats.columns = [col.replace('Ws','Cat') for col in accumCats.columns]
        upWs = strmCat.loc[strmCat.COMID.isin(tbl.catCOMID)].drop(cols, axis=1)
        newCats = pd.merge(accumCats, upWs, on="COMID")
        tbl2 = pd.merge(tbl, newCats, left_on="catCOMID", right_on="COMID")
        tbl2 = tbl2.drop(["COMID", "catCOMID"], axis=1)
        tbl2.rename(columns={"wbCOMID": "COMID"}, inplace=True)
        final = pd.concat([final, tbl2])
    return final


def makeRat(fn):
    """
    __author__ =  "Matt Gregory <matt.gregory@oregonstate.edu >"
    Adds a Raster Attribute Table to the .tif.aux.xml file, then passes those
    values to rat_to_df function to return the RAT in a pandas DataFrame.

    Arguments
    ---------
    fn               : raster filename
    """
    ds = gdal.Open(fn)
    rb = ds.GetRasterBand(1)
    nd = rb.GetNoDataValue()
    data = rb.ReadAsArray()
    # Get unique values in the band and return counts for COUNT val
    u = np.array(np.unique(data, return_counts=True))
    #  remove NoData value
    u = np.delete(u, np.argwhere(u == nd), axis=1)

    # Create and populate the RAT
    rat = gdal.RasterAttributeTable()
    rat.CreateColumn("Value", gdal.GFT_Integer, gdal.GFU_Generic)
    rat.CreateColumn("Count", gdal.GFT_Integer, gdal.GFU_Generic)
    for i in range(u[0].size):
        rat.SetValueAsInt(i, 0, int(u[0][i]))
        rat.SetValueAsInt(i, 1, int(u[1][i]))

    # Associate with the band
    rb.SetDefaultRAT(rat)

    # Close the dataset and persist the RAT
    ds = None

    # return the rat to build DataFrame
    df = rat_to_df(rat)
    return df


def rat_to_df(in_rat):
    """
    __author__ =  "Matt Gregory <matt.gregory@oregonstate.edu >"
    Given a GDAL raster attribute table, convert to a pandas DataFrame
    Parameters
    ----------
    in_rat : gdal.RasterAttributeTable
        The input raster attribute table
    Returns
    -------
    df : pd.DataFrame
        The output data frame
    """
    # Read in each column from the RAT and convert it to a series infering
    # data type automatically
    s = [
        pd.Series(in_rat.ReadAsArray(i), name=in_rat.GetNameOfCol(i))
        for i in range(in_rat.GetColumnCount())
    ]

    # Concatenate all series together into a dataframe and return
    return pd.concat(s, axis=1)


def purge(directory, pattern):
    """
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Clears directory of created rasters that will need to be re-written due to
    holding on-network like properties, i.e basins created are larger than the
    associated catchment.

    Arguments
    ---------
    directory           : directory to be cleared
    pattern             : string value to find in filename for removal
    """
    for f in os.listdir(directory):
        if re.search(pattern, f):
            os.remove(os.path.join(directory, f))


def updateSinks(wbDF, flDF):
    """
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Updates the WBARECOMI field in the NHDFlowline GeoDatFrame where NHDSinks
    intersect with NHDWaterbodies. Not currently held in the NHDPlusV21, but
    we can process the waterbodies with our on-network approach.

    Arguments
    ---------
    wbDF    : Metric name
    flDF    : Location of intermediate StreamCat csv files
    """
    flow = flDF.set_index("COMID")
    sinks = wbDF.loc[wbDF.COMID_sink.notnull()].copy()
    sinks.rename(columns={"COMID": "WBAREACOMI"}, inplace=True)
    flow.update(sinks)
    flow.WBAREACOMI = flow.WBAREACOMI.astype(np.int64)
    return flow.reset_index(level=0)


def rollArray(a, d):
    if len(d) == 4:
        a = a[0, :]
        new = np.roll(np.roll(a, d[0], axis=0), d[1], axis=1)
        new[d[2], :] = a[d[2], :]
        new[:, d[3]] = a[:, d[3]]
    if len(d) == 3:
        new = np.roll(a[0, :], d[0], axis=d[1])
        if d[1] == 0:
            new[d[2], :] = a[0, d[2], :]
        if d[1] == 1:
            new[:, d[2]] = a[0, :, d[2]]
    return np.expand_dims(new, axis=0)


def makeFlows(arr, shiftd, fdr, path, nd):
    iso = (
        np.not_equal(arr, shiftd) * np.not_equal(shiftd, nd) * np.not_equal(arr, nd)
    )  # cells change value after shift * cells not equal to NoData
    pth = np.equal(fdr, path)  # True when equal to path value
    val = iso * pth * arr
    shiftval = iso * pth * shiftd
    idx = np.not_equal(val, shiftd)
    fromcom = val[idx]
    tocom = shiftval[idx]
    fromcom = fromcom[fromcom > 0]
    tocom = tocom[tocom > 0]
    # don't load-in the entire array to the DF, just connection vals
    df = pd.DataFrame({"TOCOMID": tocom, "FROMCOMID": fromcom, "move": path})
    return df.drop_duplicates(["FROMCOMID", "TOCOMID"])


def compare_all(arr, fdr, moves, from_to, nd):
    for move in moves:
        flow = makeFlows(
            arr, rollArray(np.copy(arr), moves[move][0]), fdr, moves[move][1], nd
        )
        from_to = pd.concat([from_to, flow])
    return from_to


def expand(window, size=1):
    r, c = window
    return ((r[0] - size, r[1] + size), (c[0] - size, c[1] + size))


def check_window(window, w, h):
    r, c = window
    return ((max(0, r[0]), min(h, r[1])), (max(0, c[0]), min(w, c[1])))


def chunk_windows(r, indexes=None, max_ram=250000000):
    if indexes is None:
        indexes = r.indexes
    elif isinstance(indexes, int):
        indexes = [indexes]
    if not indexes:
        raise ValueError("No indexes to read")
    pixel_size = 0
    for bidx in indexes:
        if bidx not in r.indexes:
            raise IndexError("band index out of range")
        idx = r.indexes.index(bidx)
        pixel_size += np.dtype(r.dtypes[idx]).itemsize
    chunk_size, _ = divmod(max_ram, pixel_size)
    r_h, r_w = r.height, r.width
    if chunk_size >= r_h * r_w:
        yield (0, 0), ((0, r_h), (0, r_w))
    else:
        b_h, b_w = r.block_shapes[0]
        d, _ = divmod(chunk_size, r_w * b_h)
        chunk_height = d * b_h
        d, m = divmod(r_h, chunk_height)
        n_chunks = d + int(m > 0)
        for i in range(n_chunks):
            row = i * chunk_height
            # height = min(chunk_height, r_h - row)
            yield (i, 0), ((row, row + chunk_height), (0, r_w))


def findFlows(zone_file, fdr_file):
    moves = {
        "up": [(-1, 0, -1), 4],
        "left": [(-1, 1, -1), 1],
        "down": [(1, 0, 0), 64],
        "right": [(1, 1, 0), 16],
        "downRight": [(1, 1, 0, 0), 32],
        "downLeft": [(1, -1, 0, -1), 128],
        "upRight": [(-1, 1, -1, 0), 8],
        "upLeft": [(-1, -1, -1, -1), 2],
    }
    flows = pd.DataFrame()
    with rs.open(zone_file) as z:
        with rs.open(fdr_file) as f:  # 0 is NoData for fdr
            profile = z.profile.copy()
            nd = profile["nodata"]
            assert z.shape == f.shape, "Rasters have different extents!"
            for _, w in chunk_windows(z):  # currently defaults to 250MB
                new_w = check_window(expand(w, 2), z.width, z.height)
                data = z.read(window=new_w)
                f_r = f.read(window=new_w)
                flows = pd.concat([flows, compare_all(data, f_r, moves, flows, nd)])
    return flows.drop_duplicates(["FROMCOMID", "TOCOMID"])


def NHDtblMerge(nhd, bounds, out):
    """
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Merges all of the NHD tables needed to find on-network lakes. Returns the
    GeoDataFrames that will be used to find off-network lakes. Attribute fields
    COMID, WBARECOMI, and FEATUREID are used to link waterbodies to catchments.

    Arguments
    ---------
    nhd        : string value of prefix to NHD directory
    bounds     : GeoDataFrame of Vector Processing Unit boundaries
    out        : directory to write out to
    """
    # build dict of hr/vpu labels to read through NHD
    vpus = bounds.query("UnitType == 'VPU'").copy()
    # initialize containers to append to through processing
    onNet_connect = {}
    Obounds = gpd.GeoDataFrame()
    qa_cols = [
        "Total Waterbodies",
        "On-Network",
        "Off-network",
        "FTYPE_drop",
        "Sink_add",
        "Out_of_bounds",
    ]
    qa_tbl = pd.DataFrame()
    ons = []
    print("finding off-network lakes for ", end="", flush=True)
    for zone, hr in inputs.items():
        print(zone, end=", ", flush=True)
        pre = f"{nhd}/NHDPlus{hr}/NHDPlus{zone}"
        wbShp = gpd.read_file(f"{pre}/NHDSnapshot/Hydrography/NHDWaterbody.shp").to_crs(
            epsg="5070"
        )
        # hold length of total Waterbodies
        ttl_WB = len(wbShp)
        # format columns and select out FTYPE
        wbShp.columns = wbShp.columns[:-1].str.upper().tolist() + ["geometry"]
        wbShp = wbShp[["AREASQKM", "COMID", "FTYPE", "geometry"]]
        wbShp = wbShp.loc[wbShp["FTYPE"].isin(["LakePond", "Reservoir"])]
        # hold number of lakes removed from FTYPE
        ttl_FTYPE = ttl_WB - len(wbShp)
        fl = dbf2DF(f"{pre}/NHDSnapshot/Hydrography/NHDFlowline.dbf")[
            ["COMID", "WBAREACOMI"]
        ]
        cat = (
            gpd.read_file(f"{pre}/NHDPlusCatchment/Catchment.shp")
            .drop(["GRIDCODE", "SOURCEFC"], axis=1)
            .to_crs(epsg="5070")
        )
        cat.columns = cat.columns[:-1].str.upper().tolist() + ["geometry"]
        vaa = dbf2DF(f"{pre}/NHDPlusAttributes/PlusFlowlineVAA.dbf")[
            ["COMID", "HYDROSEQ"]
        ]
        # merge all necessary NHD tables
        final = pd.merge(
            cat.drop("geometry", axis=1),
            fl,
            left_on="FEATUREID",
            right_on="COMID",
            how="inner",
        )
        final = pd.merge(
            wbShp.drop("geometry", axis=1),
            final,
            left_on="COMID",
            right_on="WBAREACOMI",
            how="left",
            suffixes=("_wb", "_cat"),
        )
        final = pd.merge(final, vaa, left_on="COMID_cat", right_on="COMID", how="left")

        # initialize containers for on-net lakes
        cols = ["catCOMID", "wbCOMID", "CatAreaSqKm"]
        onNetDF = pd.DataFrame(columns=cols)
        catDict = {}  # holds associated lake catchments to an on-network lake

        # group by the waterbody COMID to find associated catchment
        for name, group in final.groupby("COMID_wb"):
            if not pd.isnull(group.FEATUREID).any():
                base = group.loc[group.HYDROSEQ.idxmin()]
                row = pd.Series(
                    [int(base.COMID_cat), int(base.COMID_wb), base.AREASQKM_cat],
                    index=cols,
                )
                onNetDF = onNetDF.append(row, ignore_index=True)
                catDict[int(base.COMID_cat)] = group.FEATUREID.astype(int).tolist()
        # hold length of on-net lakes
        ttl_ON = len(onNetDF)
        # add in related sinks -- this could be done in tables and found in
        # groupby object above, but this allows for isolation of number added
        sinks = gpd.read_file(f"{pre}/NHDPlusBurnComponents/Sink.shp").to_crs(
            epsg="5070"
        )
        exp = '(SOURCEFC== "NHDWaterbody")&(PURPDESC== "NHDWaterbody closed lake")'
        if len(sinks) > 0:
            sinks = sinks.query(exp)
            try:
                assert len(sinks) > 0
                catSink = sjoin(sinks, cat)
                catSink = catSink[["FEATUREID_right", "FEATUREID_left", "AREASQKM"]]
                catSink.columns = cols
                catSink = catSink.loc[
                    catSink.wbCOMID.isin(wbShp.COMID)
                ]  # this will remove any NHDWaterbody COMIDs that have the wrong FTYPE
                catSink = catSink.loc[
                    ~catSink.wbCOMID.isin(onNetDF.wbCOMID)
                ]  # remove any COMIDs that are already in the onNetDF
                ttl_SINK = len(catSink)
                for idx, line in catSink.iterrows():
                    catDict[line.catCOMID] = [line.wbCOMID]
                onNetDF = pd.concat([onNetDF, catSink])
            except AssertionError:
                ttl_SINK = 0  # all sinks got queried out!
                pass
        else:
            ttl_SINK = len(sinks)  # get val if no sinks for QA

        # create numpy arrays for connected catchments to waterbody
        onNet_connect[zone] = {
            "comids": np.array(list(catDict.keys())),
            "lengths": np.array([len(v) for v in catDict.values()]),
            "upstream": np.int32(list(chain.from_iterable(catDict.values()))),
        }

        # write-out table of catchment-lake COMID connections
        onNetDF.to_csv(f"{out}/joinTables/join_{zone}.csv", index=False)
        offLks = wbShp.loc[~wbShp.COMID.isin(onNetDF.wbCOMID)].copy()
        ons = ons + onNetDF.wbCOMID.tolist()
        # find off-netowrk lakes that are out-of-bounds
        vpu = vpus.query("UnitID == '%s'" % zone)
        offCen = offLks.copy()
        offCen.geometry = offLks.geometry.centroid
        # find centroids within the vpu
        lkVPUjoin = sjoin(offCen, vpu, op="within")[
            ["AREASQKM", "COMID", "FTYPE", "UnitID", "geometry"]
        ]
        # hold lakes that aren't within the VPU boundary
        out_of_bounds = offLks.loc[~offLks.COMID.isin(lkVPUjoin.COMID)].copy()
        out_of_bounds["VPU_orig"] = zone  # identify the zone it came from
        # find the correct vpu for those out-of-bounds
        outCen = offCen.loc[~offCen.COMID.isin(lkVPUjoin.COMID)]
        unit = sjoin(outCen, vpus, op="within")[["COMID", "UnitID"]]
        out_of_bounds = out_of_bounds.merge(unit, how="left", on="COMID")
        # add out-of-bounds to GeoDF to hold all, and select only lakes within
        # the vpu
        Obounds = pd.concat([Obounds, out_of_bounds])
        offLks = offLks.loc[offLks.COMID.isin(lkVPUjoin.COMID)].copy()

        ttl_OOB = len(out_of_bounds)
        ttl_OFF = len(offLks)
        # add VPU info to offLks table
        offLks["VPU_orig"] = zone
        vpu_tbl = sjoin(offCen, vpus, op="within")[["COMID", "UnitID"]]
        offLks = offLks.merge(vpu_tbl, on="COMID", how="left")
        # write-out off-net lakes and add series of QA info to DF
        offLks.to_file(f"{out}/off_net_{zone}.shp")
        qa_tbl[zone] = [ttl_WB, ttl_ON, ttl_OFF, ttl_FTYPE, ttl_SINK, ttl_OOB]
    # write-out all zone DF's and the numpy files created to
    assert Obounds.crs.is_projected
    # remove lakes out of bounds that would be off-network from that zone, but
    # were established as on-net in the zone they exist
    # 02/04 : coms 15516920, 15516922 NHD Problem.....again
    Obounds = Obounds.loc[~Obounds.COMID.isin(ons)]
    Obounds.to_file(f"{out}/out_of_bounds.shp")
    np.savez_compressed(f"{out}/LakeCat_npy/onNet_LakeCat.npz", vpus=onNet_connect)
    qa_tbl.index = qa_cols
    qa_tbl.T.index.rename("VPU", inplace=True)
    qa_tbl["TOTALS"] = qa_tbl.loc[:].sum(axis=1)
    qa_tbl.T.to_csv("%s/Lake_QA.csv" % out)
    print("done!")


def makeBasins(nhd, bounds, out):
    """
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"

    Makes GeoDataFrame of all lakes within each raster processing unit, converts
    each to raster, and then draws watershed basins for each lake identified as
    off-network. Creates flow table for all adjacent lake basin boundaries to
    identify zones that have hydrological connectivity.

    Arguments
    ---------
    nhd        : string value of prefloc to NHD directory
    bounds     : GeoDataFrame of Vector Processing Unit boundaries
    out        : directory to write out to
    """
    problems = pd.DataFrame()  # holding for overdrawn basin delineations,
    # i.e. bigger watershed than respective catchment
    allOff = gpd.GeoDataFrame()  # concat for all lakes
    allBsns = gpd.GeoDataFrame()  # concat for shpfile of all basins(PointInPoly)
    runit = bounds.query("UnitType == 'RPU'").copy()  # GeoDF of RPU bounds only
    Obounds = gpd.read_file(f"{out}/out_of_bounds.shp")  # lakes found out of
    # their respective hydroregion
    # arcpy.env.workspace = f"{out}/rasters/lakes/scratchArc"
    arcpy.env.outputCoordinateSystem = fiftyseventy  # output crs in ESRI WKT
    cols = ["VPU", "out_of_raster"]
    addOut = pd.DataFrame(columns=cols)  # DF to hold no. of lakes not in raster
    # for QA table
    # countTbl = pd.DataFrame() # concat for all RAT info into csv
    flow_tbl = pd.DataFrame()  # concat for all flow tables into one csv
    uid = 1000
    print("making basins for ", end="", flush=True)
    for zone, hr in inputs.items():
        print(zone, end=", ", flush=True)
        pre = f"{nhd}/NHDPlus{hr}/NHDPlus{zone}"
        # get the lakes that were out-of-bounds into the correct vpu
        addLks = Obounds.loc[Obounds.UnitID == zone].copy()
        offLks = gpd.read_file(f"{out}/off_net_{zone}.shp")
        # remove duplicated lakes across zones
        addLks.drop(addLks.loc[addLks.COMID.isin(offLks.COMID)].index, inplace=True)
        # add back-in lakes that are in other zones
        offLks = pd.concat([offLks, addLks]).reset_index(drop=True)

        assert offLks.crs.is_projected
        offLks.rename(columns={"UnitID": "VPU_moved"}, inplace=True)

        ttl_LOST = 0
        cat = gpd.read_file(f"{pre}/NHDPlusCatchment/Catchment.shp").to_crs(offLks.crs)
        for rpu in rpus[zone]:
            lakes = offLks.copy()
            if len(rpus[zone]) > 1:
                rpuShp = runit.query("UnitID == '%s'" % rpu).drop(
                    ["Hydroseq", "UnitType"], axis=1
                )
                # next line removes lakes that straddle RPU borders
                lakes = sjoin(lakes, rpuShp, op="within")
                lakes = lakes.drop("index_right", axis=1)
                lakes.rename(columns={"UnitID": "RPU"}, inplace=True)
            if len(rpus[zone]) == 1:
                lakes["RPU"] = rpu
            lakes.drop_duplicates("COMID", inplace=True)
            ln = len(lakes)
            lakes["UID"] = range(uid, uid + ln)
            uid += ln
            lakes.to_file(f"{out}/rasters/lakes/scratchArc/rasPrep_{rpu}.shp")
            fdr = arcpy.sa.Raster(f"{pre}/NHDPlusFdrFac{rpu}/fdr")
            arcpy.env.extent = fdr.extent
            arcpy.env.snapRaster = f"{pre}/NHDPlusFdrFac{rpu}/fdr"
            p2r(
                f"{out}/rasters/lakes/scratchArc/rasPrep_{rpu}.shp",
                "UID",
                f"{out}/rasters/lakes/lakes_{rpu}.tif",
                "CELL_CENTER",
                "",
                30,
            )
            purge(f"{out}/rasters/lakes/scratchArc", f"rasPrep_{rpu}")
            outWshed = Watershed(
                f"{pre}/NHDPlusFdrNull{rpu}/fdrnull",
                f"{out}/rasters/lakes/lakes_{rpu}.tif",
            )
            outWshed.save(f"{out}/rasters/wsheds/wtshds_{rpu}.tif")
            # create shp file of basins for point in poly, combine into 1 file
            shpOut = f"{out}/shps/wtshds_{rpu}.shp"
            r2p(
                f"{out}/rasters/wsheds/wtshds_{rpu}.tif", shpOut, "NO_SIMPLIFY", "VALUE"
            )
            bsnShps = gpd.read_file(shpOut)
            bsnShps.drop("Id", axis=1, inplace=True)
            bsnShps = bsnShps.dissolve(by="gridcode")
            bsnShps["UID"] = bsnShps.index
            bsnShps["AreaSqKM"] = bsnShps["geometry"].area / 10 ** 6
            # build raster attribute table,
            rat = makeRat(f"{out}/rasters/wsheds/wtshds_{rpu}.tif")
            rat.columns = rat.columns.str.upper()
            rat.rename(columns={"VALUE": "UID", "COUNT": "BSN_COUNT"}, inplace=True)
            bsnShps = bsnShps.merge(rat, on="UID")
            bsnShps["RPU"] = rpu
            bsnShps = bsnShps.merge(lakes[["COMID", "UID"]], on="UID")
            bsnShps = bsnShps[
                ["COMID", "UID", "BSN_COUNT", "AreaSqKM", "RPU", "geometry"]
            ]
            bsnShps.to_file(shpOut)
            allBsns = pd.concat([allBsns, bsnShps])
            # hold VALUE/COUNT in csv for processing
            # countTbl = pd.concat([countTbl,rat])
            # find number of lakes that don't get represented in raster, size/duplicated
            ttl_LOST += len(lakes) - len(rat)

            # lakes.to_crs(bsnShps.crs, inplace=True)
            centroids = lakes.copy().drop("AREASQKM", axis=1)
            centroids.geometry = centroids.centroid
            # if 12691112 in centroids.COMID.values:
            #     print(centroids.loc[centroids.COMID == 12691112].geometry)
            lkCat = sjoin(centroids, cat, op="intersects", how="left")
            lkCat.columns = lkCat.columns.str.upper()
            # add assoc. cats and areas to off-net lakes ----------------------
            lakes = lakes.merge(
                lkCat[["COMID", "FEATUREID", "AREASQKM"]].rename(
                    columns={"FEATUREID": "catCOMID", "AREASQKM": "catAREASQKM"}
                ),
                on="COMID",
            )
            # may be a problem coming out of the spatial join above w/ assoc. cat
            allOff = pd.concat([allOff, lakes.copy()])
            # compare basin sizes ---------------------------------------------
            both = pd.merge(
                lkCat[["COMID", "UID", "FEATUREID", "AREASQKM"]],
                rat,
                how="inner",
                on="UID",
            )
            both["AreaSqKM_basin"] = (both.BSN_COUNT * 900) * 1e-6
            bigs = both.loc[both.AREASQKM < both.AreaSqKM_basin].copy()
            bigs["diff"] = abs(bigs.AREASQKM - bigs.AreaSqKM_basin)
            bigs["VPU"] = zone
            bigs["RPU"] = rpu
            problems = pd.concat([problems, bigs], ignore_index=True)
            flow_rpu = findFlows(
                f"{out}/rasters/wsheds/wtshds_{rpu}.tif",
                f"{pre}/NHDPlusFdrFac{rpu}/fdr",
            )
            flow_rpu["RPU"] = rpu
            flow_tbl = pd.concat([flow_tbl, flow_rpu])
        # add the number of COMIDs that get lost in raster creation i.e. too small/overlap
        row = pd.Series([zone, ttl_LOST], index=cols)
        addOut = addOut.append(row, ignore_index=True)
    # write-out lakes that have a larger watershed
    # than their containing catchment
    # countTbl.to_csv("%s/LakeCat_RasterCounts.csv" % out,index=False) # this can be merged w/ allBasin.shp file
    flow_tbl.to_csv(f"{out}/LakeCat_npy/LakeCat_PlusFlow.csv", index=False)
    problems.to_csv(f"{out}/rasters/problems.csv", index=False)
    # allOff.to_crs(offLks.crs,inplace=True)
    allOff.to_file(f"{out}/off-network.shp")
    addOut.loc[len(addOut)] = pd.Series(
        ["TOTALS", addOut["out_of_raster"].sum()], index=cols
    )
    qa_tbl = pd.read_csv(f"{out}/Lake_QA.csv")
    qa_tbl = pd.merge(qa_tbl, addOut, on="VPU")
    qa_tbl.to_csv(f"{out}/Lake_QA.csv", index=False)
    allBsns.to_file(f"{out}/shps/allBasins.shp")
    purge(out, "off_net_")  # delete the individual zone files
    print("done!")


def makeNParrays(loc):
    """
    __author__ =  "Rick Debbout <debbout.rick@epa.gov>"
    Creates numpy arrays for LakeCat, uses a 'PlusFlow' table with
    TOCOMID/FROMCOMID fields along with a shapefile dbf that holds all
    of the unique id values that were used to check for connections.

    Arguments
    ---------
    loc        : location of LakeCat output directory
    """
    flow = pd.read_csv(f"{loc}/LakeCat_npy/LakeCat_PlusFlow.csv")
    fcom, tcom = flow.FROMCOMID.values, flow.TOCOMID.values
    UpCOMs = defaultdict(list)
    for i in range(0, len(flow), 1):
        FROMCOMID = fcom[i]
        if FROMCOMID == 0:
            UpCOMs[tcom[i]] = []
        else:
            UpCOMs[tcom[i]].append(FROMCOMID)
    # get unique IDs from shapefile dbf
    tbl = dbf2DF(f"{loc}/off-network.dbf")
    coms = tbl.UID.values
    d = f"{loc}/LakeCat_npy"
    if not os.path.exists(f"{d}/bastards"):
        os.mkdir(f"{d}/bastards")
        os.mkdir(f"{d}/children")
    # create bastard arrays
    a = list(map(lambda x: bastards(x, UpCOMs), coms))
    lengths = np.array([len(v) for v in a])
    a = np.int32(np.hstack(a))  # Convert to 1d vector
    np.savez_compressed(
        f"{d}/bastards/accum.npz", comids=coms, lengths=lengths, upstream=a
    )
    # create children arrays
    a = list(map(lambda x: children(x, UpCOMs), coms))
    lengths = np.array([len(v) for v in a])
    a = np.int32(np.hstack(a))  # Convert to 1d vector
    np.savez_compressed(
        f"{d}/children/accum.npz", comids=coms, lengths=lengths, upstream=a
    )
