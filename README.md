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