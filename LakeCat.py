# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:24:06 2016


@author: Rdebbout
"""

import os
import sys
from datetime import datetime as dt

import geopandas as gpd
import numpy as np
import pandas as pd

from border import makeBrdrPctFile
from lake_cat_config import LYR_DIR, NHD_DIR, OUT_DIR, STREAMCAT_DIR
from LakeCat_functions import (Accumulation, NHDtblMerge, PointInPoly,
                               chkColumnLength, doStats, getOnNetLakes, inputs,
                               makeBasins, makeNParrays, rpus)

if __name__ == "__main__":

    out = "framework"
    if not os.path.exists(out):
        print("making framework")
        if not os.path.exists(f"{out}/rasters"):
            if not os.path.exists(out):
                os.mkdir(out)
            os.mkdir(f"{out}/rasters")
            os.mkdir(f"{out}/rasters/lakes")
            os.mkdir(f"{out}/rasters/lakes/scratchArc")
            os.mkdir(f"{out}/rasters/wsheds")
            os.mkdir(f"{out}/shps")
            os.mkdir(f"{out}/joinTables")
            os.mkdir(f"{out}/LakeCat_npy")

        NHDbounds = gpd.read_file(
            f"{NHD_DIR}/NHDPlusGlobalData/BoundaryUnit.shp"
        ).to_crs(epsg="5070")
        NHDbounds.drop(
            ["AreaSqKM", "DrainageID", "Shape_Area", "Shape_Leng", "UnitName"],
            axis=1,
            inplace=True,
        )

        if not os.path.exists(f"{out}/Lake_QA.csv"):
            NHDtblMerge(NHD_DIR, NHDbounds, out)
        makeBasins(NHD_DIR, NHDbounds, out)
        makeNParrays(out)
        us_file = (
            "L:/Priv/CORFiles/Geospatial_Library_Resource/"
            "POLITICAL/BOUNDARIES/NATIONAL/TIGER_2010_State_Boundaries.shp"
        )
        bsns = "framework/shps/allBasins.shp"
        brdr = makeBrdrPctFile(us_file, bsns, "NAME10", "UID")
        os.mkdir("framework/border")
        brdr.to_csv("framework/border/pct_full.csv")

    doStats(OUT_DIR, LYR_DIR, NHD_DIR)
