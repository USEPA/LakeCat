# LakeCat

The LakeCat DataSet provides summaries of natural and anthropogenic landscape features for ~377,000 lake watershed basins within the conterminous USA using the [NHDPlus Version 2](http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php) as the geospatial framework. 

## On-Network vs. Off-Network

To begin, we determine which NHD waterbodies are on the NHD network based on joining the following tables within each hydroregion:

  * NHDWaterbodies
  * NHDFlowline
  * Catchment
  * PlusFlowlineVAA
  
*See the FindIsolatedLakes.py script for methods used*

Lakes that are associated to flowlines through the 'WBAREACOMID' in the NHDFlowline file can use the methods from StreamCat to accumulate watershed characteristics. These will be On_Network Lakes. The remaining waterbodies have gone through the off-network process to define watershed characteristics. 

## Off_Network Process

### Step 1 -- Create local basins for each lake

Operating within each raster processing unit (RPU) of the NHDPlus Version 2, all off-network lakes were converted to raster for the creation of basins with the ArcGIS Watershed tool and use the flow direction raster (fdr) . This uses the fdr with our lake raster as pour points to create watershed basin rasters. Each of the basins created are a portion of the NHD catchment which the lake falls in and define the local watershed basin for the lake. 

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

