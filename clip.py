
# coding: utf-8

# In[1]:


from osgeo import gdal, gdalnumeric, ogr ,gdal_array 
import os  
import operator  
import pdb
import numpy as np
import pdb
gdal.UseExceptions()


# In[2]:


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


# In[3]:


def pixel2World(geoMatrix, pixel, line):  
    """ 
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate 
    the pixel location  to a geospatial coordinate 
    """  
    ulX = geoMatrix[0]  
    ulY = geoMatrix[3]  
    xDist = geoMatrix[1]  
    geoX = pixel * xDist + ulX  
    geoY = ulY - line*xDist   
    return (geoX, geoY)


# In[4]:


def relu(x):
    return max(x,0)


# In[5]:


def random_crop(width, height, ob_box,clip_size):
    left_closer_to_boundary = ob_box[0] <= (width - ob_box[2])
    up_closer_to_boundary = ob_box[1] <= (height - ob_box[3])
    if left_closer_to_boundary:
        xmin = max(ob_box[0] - np.random.choice(np.arange(clip_size-(ob_box[2] - ob_box[0]))), 0)
        xmax = xmin + clip_size
    else:
        xmax = min(ob_box[2] + np.random.choice(np.arange(clip_size-(ob_box[2] - ob_box[0]))), width)
        xmin = xmax - clip_size
    if up_closer_to_boundary:
        ymin = max(ob_box[1] - np.random.choice(np.arange(clip_size-(ob_box[3] - ob_box[1]))), 0)
        ymax = ymin + clip_size
    else:
        ymax = min(ob_box[3] + np.random.choice(np.arange(clip_size-(ob_box[3] - ob_box[1]))), height)
        ymin = ymax - clip_size
    return xmin, ymin, xmax, ymax


# In[6]:


def get_bndbox(ob_box, width=float('inf'), height=float('inf')):
    
    ob_xmin = relu(ob_box[0])
    ob_ymin = relu(ob_box[1])
    ob_xmax = min(ob_box[2], width)
    ob_ymax = min(ob_box[3], height)
    return ob_xmin, ob_ymin, ob_xmax, ob_ymax


# In[7]:


def include_ratio(crop, ob_box):
    lu = np.maximum(crop[0:2], ob_box[0:2])
    rd = np.minimum(crop[2:], ob_box[2:])
    intersection = np.maximum(rd - lu, [0,0])
    inter_square = intersection[0]*intersection[1]
    square_ob = (ob_box[2]-ob_box[0]) * (ob_box[3] - ob_box[1])
    return inter_square / (square_ob + 1e-5)


# In[8]:


def find_all_objects_in_a_crop(crop, imgPts,avg_tree_rds,width, height, include_threshold=0.7):
    # @param width: original image width
    # @param height: original image height
    obs_str = ''
    center_str = ''
    for imgPt in imgPts:
        ob_box = [imgPt[0]-avg_tree_rds,imgPt[1]-avg_tree_rds,imgPt[0]+avg_tree_rds,imgPt[1]+avg_tree_rds]
        bndbox = get_bndbox(ob_box, width, height)
        if include_ratio(crop, bndbox) > include_threshold:
            obs_str += ' ' + str(max(bndbox[0] - crop[0], 0))
            obs_str += ' ' + str(max(bndbox[1] - crop[1], 0))
            obs_str += ' ' + str(min(bndbox[2] - crop[0], crop[2]- crop[0]-1))
            obs_str += ' ' + str(min(bndbox[3] - crop[1], crop[3]- crop[1]-1))
            obs_str += ' ' + '0'
            center_str += ' ' + str(min(max((imgPt[0] - crop[0]),0),crop[2]- crop[0]-1))
            center_str += ' ' + str(min(max((imgPt[1] - crop[1]),0),crop[3]- crop[1]-1))
            center_str += ' ' + '0'
    obs_str += '\n'
    center_str += '\n'
    return obs_str,center_str


# In[9]:


def process_a_file(shapefile_path, raster_path,outdir,clip_size = 448,avg_tree_rds=40):
    shapef = ogr.Open(shapefile_path)  
    lyrname = os.path.split( os.path.splitext( shapefile_path )[0] )[1] 
    lyr = shapef.GetLayer(lyrname) 
    
    srcImage = gdal.Open(raster_path)
    geoTrans = srcImage.GetGeoTransform()
    oriHei = srcImage.RasterYSize  
    oriWid = srcImage.RasterXSize 
    srcArray = gdalnumeric.LoadFile(raster_path)
    
    feat = lyr.GetNextFeature()  
    geoPts = []
    while feat is not None:  
        geopt=feat.GetGeometryRef()
        geoPts.append([geopt.GetX(), geopt.GetY()])
        feat = lyr.GetNextFeature() 
        
    out_filename = outdir + lyrname +'_clip_label.txt'
    out_file = open(out_filename, 'a')

    for i,geopt in enumerate(geoPts):
        lyr.ResetReading() 
        imgX, imgY = world2Pixel(geoTrans, geopt[0], geopt[1])
        #imgPts.append([imgX, imgY])
        ob_box = [imgX-avg_tree_rds,imgY-avg_tree_rds,imgX+avg_tree_rds,imgY+avg_tree_rds]
        #ob_boxes.append(ob_box)
        bnd_box = get_bndbox(ob_box,oriWid,oriHei)
        crop = random_crop(oriWid, oriHei, bnd_box,clip_size)
        geoUl = pixel2World(geoTrans,crop[0],crop[1])
        geoBr = pixel2World(geoTrans,crop[2],crop[3])
        lyr.SetSpatialFilterRect(geoUl[0],geoUl[1],geoBr[0],geoBr[1])

        imgPts = []
        ob_boxes = []
        oFeature = lyr.GetNextFeature()  

        while oFeature is not None:  
            inRectPt=oFeature.GetGeometryRef()
            imgPt = world2Pixel(geoTrans,inRectPt.GetX(), inRectPt.GetY())
            imgPts.append([imgPt[0],imgPt[1]])
            ob_box = [imgPt[0]-avg_tree_rds,imgPt[1]-avg_tree_rds,imgPt[0]+avg_tree_rds,imgPt[1]+avg_tree_rds]
            ob_boxes.append(ob_box)
            oFeature = lyr.GetNextFeature() 
        obs_str,center_str = find_all_objects_in_a_crop(crop, imgPts, avg_tree_rds,oriWid, oriHei, include_threshold=0.7)
        
        crop_name = outdir + lyrname + '/' + lyrname +'_ob{}_crop.jpg'.format(i)
        record = crop_name + obs_str
        out_file.write(record)
        # Save as an 8-bit jpeg for an easy, quick preview  
        clip = srcArray[:, crop[1]:crop[3], crop[0]:crop[2]]
        clip = clip.astype(gdalnumeric.uint8)  
        gdalnumeric.SaveArray(clip, crop_name, format="JPEG")
    out_file.close()


# In[10]:


clip_size = 448
avg_tree_rds=40
shapefile_path = '/home/CV/geocv/data/geo/discolor/gylc_84.shp'  
raster_path = '/home/CV/geocv/data/geo/2016102002.tif'  
outdir = '/home/CV/geocv/ML/geotree/'


# In[ ]:


process_a_file(shapefile_path, raster_path,outdir,448,40)

