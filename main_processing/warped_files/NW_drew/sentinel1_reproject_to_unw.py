#!/usr/bin/env python

import string
import os
#import calendar
import sys
from osgeo import gdal
#import datetime
#import math
#import util_mod 
#import numpy
#import csv
#import insar_gdal_tools as tools


if len(sys.argv) < 2:
    print('Usage: sentinel1_reproject_to_unw.py dirlist target.wkt ref_file')
    sys.exit(1)


u_list = sys.argv[1]
wkt = sys.argv[2]
ref_file = sys.argv[3]

f = open(u_list,'r')
swaths = []
for line in f:
    swaths.append(line.strip())
f.close()


# open reference file and get resolution
referenceFile = ref_file
print(referenceFile)
reference = gdal.Open(referenceFile, 0)  # this opens the file in only reading mode
referenceTrans = reference.GetGeoTransform()
x_res = referenceTrans[1]
y_res = -referenceTrans[5]  # make sure this value is positive
ulx, xres, xskew, uly, yskew, yres  = reference.GetGeoTransform()
lrx = ulx + (reference.RasterXSize * xres)
lry = uly + (reference.RasterYSize * yres)



for (i,fil) in enumerate(swaths):
        #os.chdir(str(fil))
        filname=str(fil)
        inc_in=filname#+'_inc_map.tif'
        inc_out='warped_'+filname[3:]
        #vh_in=filname+'_VH.tif'
        #vh_out='warped_'+filname[7:15]+'_vh.tif'
        #vv_in=filname+'_VV.tif'
        #vv_out='warped_'+filname[7:15]+'_vv.tif'
        output1='gdalwarp -t_srs ' + str(wkt) + ' -tr ' + str(x_res) + ' ' + str(y_res ) +' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' '+str(inc_in)+' '+str(inc_out)
        ret=os.system(output1)
        #output2='gdalwarp -t_srs ../'+str(wkt)+' -tr '+str(x_res)+' '+str(y_res)+' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' '+str(vh_in)+' '+str(vh_out)
        #ret=os.system(output2)
        #output3='gdalwarp -t_srs ../'+str(wkt)+' -tr '+str(x_res)+' '+str(y_res)+' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' '+str(vv_in)+' '+str(vv_out)
        #ret=os.system(output3)
        #os.chdir('..')
