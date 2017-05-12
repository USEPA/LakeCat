# LakeCat

The LakeCat DataSet provides summaries of natural and anthropogenic landscape features for ~377,000 lakes and their associated catchments within the conterminous USA using the [NHDPlus Version 2](http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php) (NHDPlusV2) as the geospatial framework. 

## Determining On-Network and Off-Network Lakes

To begin, we determine which NHD waterbodies are on the NHD network based on joining the following tables within each NHDPlusV2 HydroRegion (see code in *FindIsolatedLakes.py*):

  * NHDWaterbodies
  * NHDFlowline
  * Catchment
  * PlusFlowlineVAA
  
Lakes that are associated to flowlines through the 'WBAREACOMID' in the NHDFlowline file can use data directly from StreamCat to represent watershed-level landscape characteristics and are designated as *On_Network* lakes. 

The remaining waterbodies are not on the stream network (as depicted with the NHDFlowline file) and are designated as *Off_Network* lakes and saved as *IsolatedLakes.shp*. Watershed characteristics are developed for these lakes using the Off_Network Process (see below).

## Off_Network Process

### Step 1 -- Create local basins for each lake

Operating within each raster processing unit (RPU) of the NHDPlusV2, all off-network lakes were converted to ESRI shapefiles to rasters (.TIFF format). We then used the ArcGIS Watershed tool and the NHDPlusV2 flow direction raster (named fdr within the NHDPlusV2) to delineate the contributing areas to each lake. This process uses the lake raster as pour points on the fdr raster to create watersheds. Each watershed is a portion of the NHDPlusV2 catchment which the lake falls within and define the local watershed basin for the lake. 

![Off-Network](https://cloud.githubusercontent.com/assets/7052993/19703884/648f7f0e-9aba-11e6-90e0-e909b49f5de2.PNG)

*Note that the off-network lake basin is only a portion of the catchment which it is in.*

### Step 2 -- Find flow connections between basins

Using adjacent local basins we can find flow connections between them with the ArcGIS shift tool and the RPU's flow direction raster. The script, LakeConnect.py, shifts the raster cells of each zones local watershed basins and tests if the value of that shifted cell is not the value that it holds, and second, then tests if the flow direction raster value demonstrates that there is flow in the same direction as the direction of the shift.  

![shifted](https://cloud.githubusercontent.com/assets/7052993/19706148/306e4948-9ac5-11e6-9a80-c7e3362f7bc1.PNG)![directions](https://cloud.githubusercontent.com/assets/7052993/19816175/222618ce-9cfb-11e6-9290-9c737bb0adb2.PNG)

This is done for each of the 8 directions that a raster cell can be moved, shown above. If both conditions are true the value of the basin IDs are saved in a flow table to be used to create full watersheds in the accumulation process.



### Step 3 -- Perform zonal statistics and accumulate

The *LakeCat.py* script creates tables for each landscape metric summarizing both for individual lake basins and for cumulative upstream watersheds. This process happens independently of the On-Network process and then tables from both processes are appended together.

## On_Network Process

### Step 1 -- Find lakes On_Network

The *findIsolatedLakes.py* script uses table joins to determine lakes that are in the NHDPlus version 2 network. The table below shows the keys that connect all of the NHDPlusV2 tables.

![tableconnect2](https://cloud.githubusercontent.com/assets/7052993/19823341/f037171a-9d1c-11e6-84bb-5035685a7b2e.PNG)

Once all of these tables are joined, we can group them by the waterbody COMID. In each group there can be more than one flowline associated to each waterbody.  Using the 'Hydroseq' attribute, we can find the catchment 'COMID' that is furthest downstream, and use this 'COMID to refer to StreamCat data and obtain summarized characteristics of each landscape layer. This connection is stored in a lookup table, CATCHMENT COMID <--> WATERBODY COMID. LakeCat final tables reflect WATERBODY COMID in their COMID field. 

### Step 2 -- Accumulate associated catchments for local area 'Cat' metrics

StreamCat tables contain summarizations for individual stream catchments and for cumulative upstream watersheds. In order to adjust the catchment metrics, we take all of the associated catchments to each lake and summarize metrics for the entire area that they cover. The *findIsolatedLakes.py* script creates arrays that hold the associated catchment COMIDs with the lakes most downstream catchment COMID. The same accumulation function is used to create summary statistics across all catchments associated with each waterbody, and is merged into the final table with the watershed metrics.

![multfl](https://cloud.githubusercontent.com/assets/7052993/19825737/e17463a8-9d31-11e6-9a69-a4d1364b6aea.png)

*Shows multiple flowlines and associated catchments to a given waterbody -- highlighted in yellow*

## QA Table from On/Off network selection and basin creation

![qa_tbl](https://cloud.githubusercontent.com/assets/7052993/23385978/7804e3e8-fd08-11e6-84f2-b5a0f69e3324.PNG)

## Issues

While creating this dataset, a small number of lakes had to be left out of our processing model due to one of a few reasons.

### Off-Network Lake displays On-Network accumulation properties

Within the context of all of the NHD data we are using to model lake watersheds, there are some inconsistencies between what is published in the vector format vs the raster format. For the off-network lake basins to be made we rely on the flow direction raster which, at 30 meter resolution, isn't able to account for all of the detail needed to describe areas with low slope. 

After isolating the off-network lakes and creating rasters of lake bodies within each raster processing unit (RPU), we compared the area covered in the basins with the catchment area and found that some basins were created that are significantly larger than the catchment which the lake exists in.

![extend](https://cloud.githubusercontent.com/assets/7052993/19867420/fb79beac-9f60-11e6-8a4d-c2e66a5e221e.PNG)![arid](https://cloud.githubusercontent.com/assets/7052993/19869318/0aa3409e-9f69-11e6-9dd8-f005fc9279ba.PNG)

*I need to learn the association made with flowlines and waterbodies and how they work. There are 758 lakes that got removed in this process and it seems that they should be on-network but don't link up with any flowlines. I'm unsure of how to state this in an intelligent way other than just mention the the waterbody boundary in raster format covers the fdr in a spot that should put it on-network, yet there isn't a flowline to make that connection. In this case, the lake is in an arid zone and doesn't even exist, and according to google imagery neither does the flowline. So it is an area that is only seasonally containing water features.*

Waterbodies are sometimes very close to flowlines, yet have no association with them. 

### In_Network accumulated watersheds have smaller area than the waterbody

During QA, we noticed that 384 accumulated watersheds for lakes are smaller than the lake itself based on area alone. 

### NHD Catchments cover a smaller area than the NHD Waterbody -- On-Network

Defining NHD Waterbodies by associated flowlines and ultimately catchments has proved impossible in some cases where one catchment is associated to the waterbody and the catchment associated is smaller than the waterbody itself. Thus, the 'cat' metric is not accurate to the basin in which should hold the lake.

##QA Methods

Lakes outside of the boundaries of the NHDPlusV2 -969 lakes. Removed in the findIsolatedLakes.py script after all isolated lakes are appended together. SelectLayerByLocation_management tool finds all lakes that are "COMPLETELY_WITHIN" and the performs a "SWITCH_SELECTION" to arrive at the lakes that are outside of the boundary and deleted from the file.

On_Net_Cvg_Diffs -- come from the script chkNLWwLakeCat.py comparing the lake area geometry with the geometry of what we use as an on-network catchment basin. This outputs the catchment COMID and associated catchment COMID with the percentage of lake that is uncovered by the basins geometry.

off-network problem lakes -- findProblemLakes.py prints out table where an off-network lake basin's area is larger than the catchment in which the lake is found. This table was used to remove 758 lakes from the IsolatedLakes.shp file.

##updating: 12/20/2016

![vpu_join3](https://cloud.githubusercontent.com/assets/7052993/21364695/c1dfe21c-c6a6-11e6-8671-cf7f4808428d.png)

###New order of spatial joins/ table merges:
* do the sjoin with the vpu AFTER we find on-network lakes to hold onto problem lakes from other zones
* use sinks.shp to  hold onto more on-network lakes that get overdrawn
* 
* 

## On-Network Lakes left out due to zone '04' exclusion in StreamCat



# COMID from zone 17 of three lakes that have interesting flow (one<--2-->theOther) && one doesn't flow to the other 2!!
* 23043937

# Duplicated in the off-network process or in the NHD:
* 13871500,  7109029, 18156163, 13118610

# 12 repeated COMIDs in off_networks.shp need to be filtered out!

* out of bounds but duplicated between the 2 zones! fixed with drop_duplicates!


## <font color='red'>Omitted NHDPlusV21 waterbody COMIDs</font>
COMID |	VPU |	REASON FOR OMISSION
:------|:---:|--------------------:		
14300071	|10U|	OVERLAPPED BY COMID 12568346
12967474	|10U|	DOESN'T HIT CELL CENTER
120052923	|10U|	OVERLAPPED BY COMID 120051949
167245973	|10L|	OVERLAPPED BY COMID 120053749
14819137	|7|	DOESN'T HIT CELL CENTER
14820667	|7|	DOESN'T HIT CELL CENTER
944040052	|14|	DOESN'T HIT CELL CENTER
23854057	|17|	DOESN'T HIT CELL CENTER
24052877	|17|	OVERLAPPED BY COMID 20315924
120054048	|16|	OVERLAPPED BY COMID 120053946
19861298	|16|	OVERLAPPED BY COMID 24989585
162424445	|15|	OVERLAPPED BY COMID 120053954
22323839	|9|	DOESN'T HIT CELL CENTER
15235574	|8|	DOESN'T HIT CELL CENTER
22700864	|8|	DOESN'T HIT CELL CENTER
14725382	|4|	OVERLAPPED BY COMID 166766632
20166532	|03N|	DOESN'T HIT CELL CENTER
120054030	|18|	OVERLAPPED BY COMID 120053926
