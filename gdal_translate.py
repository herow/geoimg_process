#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:36:47 2018

@author: geo
"""

from osgeo import gdal, gdalnumeric, ogr ,gdal_array 
#from PIL import Image, ImageDraw  
import os  
import operator  
import pdb
gdal.UseExceptions()

def world2Pixel(geoMatrix, x, y):  
    """ 
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate 
    the pixel location of a geospatial coordinate 
    """  
    ulX = geoMatrix[0]  
    ulY = geoMatrix[3]  
    xDist = geoMatrix[1]  
    pixel = int((x - ulX) / xDist)  
    line = int((ulY - y) / xDist)  
    return (pixel, line) 
def main( shapefile_path, raster_path ):  
    # Load the source data as a gdalnumeric array  
    
  
    # Also load as a gdal image to get geotransform  
    # (world file) info  
    pdb.set_trace()
    
    srcImage = gdal.Open(raster_path)  
   
    geoTrans = srcImage.GetGeoTransform()  
  
    # Create an OGR layer from a boundary shapefile  
    shapef = ogr.Open(shapefile_path)  
    lyr = shapef.GetLayer( os.path.split( os.path.splitext( shapefile_path )[0] )[1] ) 
    
   
  
    # Convert the layer extent to image pixel coordinates  
    minX, maxX, minY, maxY = lyr.GetExtent()  
    ulX, ulY = world2Pixel(geoTrans, minX, maxY)  
    lrX, lrY = world2Pixel(geoTrans, maxX, minY) 
    
    # Calculate the pixel size of the new image  
    pxWidth = int(lrX - ulX)  
    pxHeight = int(lrY - ulY)  
  
   

    ds = gdal.Translate('beijing2.tif',srcImage,srcWin=[ulX,ulY,pxWidth,pxHeight])
    if ds is None:
        print "failed"
if __name__ == '__main__':  
    #shapefile_path, raster_path  
    shapefile_path = '/home/CV/geocv/data/geo/discolor/gylc_84.shp'  
    raster_path = '/home/CV/geocv/data/geo/2016102002.tif'  
    main( shapefile_path, raster_path )  