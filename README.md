# LakeCat

The LakeCat DataSet provides summaries of natural and anthropogenic landscape features for ~2.65 million lake watershed basins within the conterminous USA using the [NHDPlus Version 2](http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php) as the geospatial framework.

## Script Processing

To begin, we determine which NHD waterbodies are on the NHD network based on joining the following tables within each hydroregion:

  * NHDWaterbodies
  * NHDFlowline
  * Catchment
  * PlusFlowlineVAA
  
Lakes that are associated to flowlines and can use the methods from StreamCat to accumulate watershed characteristics. The remaining waterbodies have basin areas that were created using the ArcGIS Watershed tool.  These basins have been connected where there are adjacent cells to other basins and there is flow in the given direction from the flow direction raster.





