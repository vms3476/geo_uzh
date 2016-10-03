% script for batch converting geotiffs to osgdem readable format
cd /Users/morsdorf/data/KantonAargau/ORT
filez = dir('*.tif');
for i = 1:length(filez)
    %disp(['gdal_translate -scale 0 65536 0 255 -ot Byte ',filez(i).name,' ../ORT/',filez(i).name])
    %unix(['gdal_translate -scale 0 65536 0 255 -ot Byte ',filez(i).name,' ../ORT/',filez(i).name])
    disp(['gdaladdo ',filez(i).name,' 2 4 8 16'])
    unix(['gdaladdo ',filez(i).name,' 2 4 8 16'])
end
%cd /Volumes/FastDisk/RGB_AG
%disp(['osgdem --terrain --no-mip-mapping -o RGB3D.ive -d DTM -t /Users/morsdorf/data/KantonAargau/ORT'])
%unix(['osgdem --terrain --no-mip-mapping -o RGB3D.ive -d DTM -t /Users/morsdorf/data/KantonAargau/ORT'])