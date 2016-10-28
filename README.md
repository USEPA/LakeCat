# LakeCat

The LakeCat DataSet provides summaries of natural and anthropogenic landscape features for ~2.65 million lake watershed basins within the conterminous USA using the [NHDPlus Version 2](http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php) as the geospatial framework. 

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

Operating within each raster processing unit (RPU) of the NHDPlus Version 2, all off-network lakes were converted to raster for the creation of basins with the ArcGIS Watershed tool. This uses the flow accumulation raster with our lake raster as pour points to create watershed basin rasters. Each of the basins created are a portion of the NHD catchment which the lake falls in and define the local watershed basin for the lake. 

![Off-Network](https://cloud.githubusercontent.com/assets/7052993/19703884/648f7f0e-9aba-11e6-90e0-e909b49f5de2.PNG)

### Step 2 -- Find flow connections between basins

Using adjacent local basins we can find flow connections between them with the ArcGIS shift tool and the RPU's flow direction raster. The script, LakeConnect.py, shifts the raster cells of each zones local watershed basins and tests if the value of that shifted cell is not the value that it holds, and second, then tests if the flow direction raster value demonstrates that there is flow in the same direction as the direction of the shift.  

![shifted](https://cloud.githubusercontent.com/assets/7052993/19706148/306e4948-9ac5-11e6-9a80-c7e3362f7bc1.PNG)

This is done for each of the 8 directions that a raster cell can be moved, shown below.  If both conditions are true the value of the basin IDs are saved in a flow table to be used to create full watersheds in the accumulation process.

![directions](https://cloud.githubusercontent.com/assets/7052993/19816175/222618ce-9cfb-11e6-9290-9c737bb0adb2.PNG)

### Step 3 -- Perform zonal statistics and accumulate




