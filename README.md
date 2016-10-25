# LakeCat

The LakeCat DataSet provides summaries of natural and anthropogenic landscape features for ~2.65 million lake watershed basins within the conterminous USA using the [NHDPlus Version 2](http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php) as the geospatial framework.

## On-Network vs. Off-Network

To begin, we determine which NHD waterbodies are on the NHD network based on joining the following tables within each hydroregion:

  * NHDWaterbodies
  * NHDFlowline
  * Catchment
  * PlusFlowlineVAA
  
Lakes that are associated to flowlines and can use the methods from StreamCat to accumulate watershed characteristics The remaining waterbodies have gone through the following process to define watershed characteristics:

## Off_Network Process

### Step 1

Operating within each raster processing unit (RPU), all off-network lakes were converted to raster for the creation of basins with the ArcGIS Watershed tool. This uses the flow accumulation raster with our lake raster aws pour points to create watershed basins. Each of the basins created are a portion of the catchment which the lake falls in and define the local watershed basin for the lake. 

![Off-Network](C:\Users\Rdebbout\Downloads\LakesREADME\off_network.PNG)

Using adjacent local basins we can find flow connections between them with the ArcGIS shift tool and the RPU's flow direction raster.


have basin areas that were created using the ArcGIS Watershed tool.  These basins have been connected where there are adjacent cells to other basins and there is flow in the given direction from the flow direction raster.





