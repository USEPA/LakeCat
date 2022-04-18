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
from lake_cat_config import LYR_DIR, NHD_DIR, OUT_DIR, STREAMCAT_DIR, FRAMEWORK
from LakeCat_functions import (Accumulation, NHDtblMerge, PointInPoly,
                               chkColumnLength, doStats, getOnNetLakes, inputs,
                               makeBasins, makeNParrays, rpus)

if __name__ == "__main__":

    if not os.path.exists(FRAMEWORK):
        print("making framework")
        os.mkdir(FRAMEWORK)
        os.mkdir(f"{FRAMEWORK}/rasters")
        os.mkdir(f"{FRAMEWORK}/rasters/lakes")
        os.mkdir(f"{FRAMEWORK}/rasters/lakes/scratchArc")
        os.mkdir(f"{FRAMEWORK}/rasters/wsheds")
        os.mkdir(f"{FRAMEWORK}/shps")
        os.mkdir(f"{FRAMEWORK}/joinTables")
        os.mkdir(f"{FRAMEWORK}/LakeCat_npy")

        NHDbounds = gpd.read_file(
            f"{NHD_DIR}/NHDPlusGlobalData/BoundaryUnit.shp"
        ).to_crs(epsg="5070")
        NHDbounds.drop(
            ["AreaSqKM", "DrainageID", "Shape_Area", "Shape_Leng", "UnitName"],
            axis=1,
            inplace=True,
        )

        if not os.path.exists(f"{FRAMEWORK}/Lake_QA.csv"):
            NHDtblMerge(NHD_DIR, NHDbounds, FRAMEWORK)
        makeBasins(NHD_DIR, NHDbounds, FRAMEWORK)
        makeNParrays(FRAMEWORK)
        us_file = (
            "L:/Priv/CORFiles/Geospatial_Library_Resource/"
            "POLITICAL/BOUNDARIES/NATIONAL/tl_2018_us_state.shp"
        )
        bsns = f"{FRAMEWORK}/shps/allBasins.shp"
        brdr = makeBrdrPctFile(us_file, bsns, "NAME", "UID")
        os.mkdir(f"{FRAMEWORK}/border")
        brdr.to_csv(f"{FRAMEWORK}/border/pct_full.csv")

    doStats(OUT_DIR, LYR_DIR, NHD_DIR, FRAMEWORK)
