# LakeCat

## Description:
The LakeCat DataSet (http://www2.epa.gov/national-aquatic-resource-surveys/lakecat) provides summaries of natural and anthropogenic landscape features for ~378,000 lakes and their associated catchments within the conterminous USA. This repo contains code used in LakeCat to process a suite of landscape rasters to produce watershed metrics for these lakes.

## Necessary Python Packages and Installation Tips

The scripts for LakeCat rely on several python modules a user will need to install such as numpy, pandas, gdal, fiona, rasterio, geopandas, shapely, pysal, and ArcPy with an ESRI license (minimal steps still using ArcPy).  We highly recommend using a scientific python distribution such as [Anaconda](https://www.continuum.io/downloads) or [Enthought Canopy](https://www.enthought.com/products/canopy/).  We used the conda package manager to install necessary python modules. Our essential packages and versions used are listed below (Windows 64 and Python 2.7.11):

| Package       | Version       | 
| ------------- |--------------:|
| fiona         | 1.7.7         | 
| gdal          | 2.2.0         | 
| geopandas     | 0.2.1         |  
| geos          | 3.5.1         |
| libgdal       | 2.0.0         |
| numpy         | 1.12.1        |
| pandas        | 0.20.2        |
| pyproj        | 1.9.5.1       |
| pysal         | 1.13.0        |
| rasterio      | 1.0a9         |
| shapely       | 1.5.17        |

If you are using Anaconda, creating a new, clean 'LakeCat' environment with these needed packages can be done easily and simply one of several ways:

* In your conda shell, add one necessary channel and then download the lakecat environment from the Anaconda cloud:
  + conda config --add channels conda-forge
  + conda env create mweber36/lakecat
  
* Alternatively, using the lakecat.yml file in this repository, in your conda shell cd to the directory where your lakecat.yml file is located and run:
  + conda env create -f LakeCat.yml
  
* To build environment yourself, do:
  + conda env create -n LakeCat rasterio geopandas
  + pip install georasters

* To activate this new environment and open Spyder, type the following at the conda prompt
  + activate LakeCat
  
  Then

  + Spyder
  
Finally, to use arcpy in this new environment, you will need to copy your Arc .pth file into your new environment.  Copy the .pth file for your install of ArcGIS located in a directory like:

+ C:\Python27\ArcGISx6410.3\Lib\site-packages\DTBGGP64.pth

To your environment directory which should look something like:

+ C:\Anaconda\envs\lakecat\Lib\site-packages\DTBGGP64.pth

Note that the exact paths may vary depending on the version of ArcGIS and Anaconda you have installed and the configuration of your computer

##How to Run Scripts 

###The scripts make use of 'control tables' to pass all the particular parameters to the two primary scripts: 

+ [LakeCat.py](https://github.com/USEPA/LakeCat/blob/master/LakeCat.py).
+ [MakeFinalTables_LakeCat.py](https://github.com/USEPA/LakeCat/blob/master/MakeFinalTables_LakeCat.py).  

In turn, these scripts rely on a set of functions in [LakeCat_functions.py](https://github.com/USEPA/LakeCat/blob/master/LakeCat_functions.py). 

A table with all required parameters is used to process landscape layers in LakeCat:
+ [ControlTable_LakeCat](https://github.com/USEPA/LakeCat/blob/master/ControlTable_LakeCat.csv)


###Running LakeCat.py to generate new LakeCat metrics

After editing the control tables to provide necessary information, such as directory paths, the following stesps will excecute processes to generate new watershed metrics for the conterminous US. All examples in the control table are for layers (e.g., STATSGO % clay content of soils) that were processed as part of the LakeCat Dataset. This example assumes run in Anaconda within Conda shell.

1. Edit [ControlTable_LakeCat](https://github.com/USEPA/LakeCat/blob/master/ControlTable_LakeCat.csv) and set desired layer's "run" column to 1. All other columns should be set to 0
2. Open a Conda shell and type "activate LakeCat" 
3. At the Conda shell type: "Python<space>"
4. Drag and drop "LakeCat.py" to the Conda shell from a file manager followed by another space
5. Drag and drop the control table to the Conda shell

Final text in Conda shell should resemble this: python C:\some_path\LakeCat.py  C:\some_other_path\ControlTable.csv


## EPA Disclaimer
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use.  EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information.  Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA.  The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.