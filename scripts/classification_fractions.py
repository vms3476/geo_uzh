import numpy
from osgeo import gdal
import os

directory = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch4_4aoi/'

os.chdir(directory)

nBroadleaf = 0
nConiferous = 0 
nNonVeg = 0
nTotal = 0 

for filename in os.listdir(directory): 

   print filename

   ds = gdal.Open(filename)

   array = numpy.array(ds.GetRasterBand(1).ReadAsArray())

   nTotal += array.size # total number of pixels 
   nBroadleaf += numpy.sum(array==1)
   nConiferous += numpy.sum(array==2)
   nNonVeg += numpy.sum(array==0)


print 'fraction of broadleaf pixels: ', nBroadleaf * 1.0 / nTotal
print 'fraction of coniferous pixels: ', nConiferous * 1.0 / nTotal
print 'fraction of non veg pixels: ', nNonVeg * 1.0 / nTotal


