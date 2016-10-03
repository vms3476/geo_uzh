function [] = shppolytokml(fname)
% function to convert polygon in shapefile to polygon in xml file 
% converts as well coordinates from CH1903 to UTM

% Felix Morsdorf, RSL September 2010

shp = shaperead(fname);
ii = 1:length(shp.X)-1;
pos = [shp.X(ii);shp.Y(ii);zeros(size(shp.X(ii)))*2000]';
coor = lk2ch(pos,'CH1903');
%wgs = ch2wgsxyz(coor,'WGS84');
fname = regexprep(fname,'.shp','.kml');
kmlwrite(fname,coor(:,2),coor(:,1));
keyboard
