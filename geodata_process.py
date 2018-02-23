#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 18:29:21 2018

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
#  
#  EDIT: this is basically an overloaded  
#  version of the gdal_array.OpenArray passing in xoff, yoff explicitly  
#  so we can pass these params off to CopyDatasetInfo  
#  
def OpenArray( array, prototype_ds = None, xoff=0, yoff=0 ):  
    ds = gdal_array.OpenArray( gdalnumeric.GetArrayFilename(array) )  
  
    if ds is not None and prototype_ds is not None:  
        if type(prototype_ds).__name__ == 'str':  
            prototype_ds = gdal.Open( prototype_ds )  
        if prototype_ds is not None:  
            gdalnumeric.CopyDatasetInfo( prototype_ds, ds, xoff=xoff, yoff=yoff )  
    return ds  
def histogram(a, bins=range(0,256)):  
    """ 
    Histogram function for multi-dimensional array. 
    a = array 
    bins = range of numbers to match 
    """  
    fa = a.flat  
    n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)  
    n = gdalnumeric.concatenate([n, [len(fa)]])  
    hist = n[1:]-n[:-1]  
    return hist  
  
def stretch(a):  
    """ 
    Performs a histogram stretch on a gdalnumeric array image. 
    """  
    hist = histogram(a)  
    im = arrayToImage(a)  
    lut = []  
    for b in range(0, len(hist), 256):  
        # step size  
        step = reduce(operator.add, hist[b:b+256]) / 255  
        # create equalization lookup table  
        n = 0  
        for i in range(256):  
            lut.append(n / step)  
            n = n + hist[i+b]  
        im = im.point(lut)  
    return imageToArray(im)  

def main( shapefile_path, raster_path ):  
    # Load the source data as a gdalnumeric array  
    srcArray = gdalnumeric.LoadFile(raster_path)  
  
    # Also load as a gdal image to get geotransform  
    # (world file) info  
    pdb.set_trace()
    
    srcImage = gdal.Open(raster_path)  
    geoTrans = srcImage.GetGeoTransform()  
  
    # Create an OGR layer from a boundary shapefile  
    shapef = ogr.Open(shapefile_path)  
    lyr = shapef.GetLayer( os.path.split( os.path.splitext( shapefile_path )[0] )[1] ) 
    
    poly = lyr.GetNextFeature()  
  
    # Convert the layer extent to image pixel coordinates  
    minX, maxX, minY, maxY = lyr.GetExtent()  
    ulX, ulY = world2Pixel(geoTrans, minX, maxY)  
    lrX, lrY = world2Pixel(geoTrans, maxX, minY) 
    
    # Calculate the pixel size of the new image  
    pxWidth = int(lrX - ulX)  
    pxHeight = int(lrY - ulY)  
  
    clip = srcArray[:, ulY:ulY+1000, ulX:ulX+1000] 

    # Save as an 8-bit jpeg for an easy, quick preview  
    clip = clip.astype(gdalnumeric.uint8)  
    gdalnumeric.SaveArray(clip, "beijing3.jpg", format="JPEG")  
if __name__ == '__main__':  
    #shapefile_path, raster_path  
    shapefile_path = '/home/CV/geocv/data/geo/discolor/gylc_84.shp'  
    raster_path = '/home/CV/geocv/data/geo/2016102002.tif'  
    main( shapefile_path, raster_path )  