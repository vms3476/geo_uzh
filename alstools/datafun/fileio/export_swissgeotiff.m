function [] = export_swissgeotiff(dat,fname);
% function [] = export_swissgeotiff(dat,fname);
% dat must be struct with fields x,y,z
dat.z = flipud(dat.z);
are = [min(dat.x) min(dat.y) max(dat.x) max(dat.y)];
option.GTModelTypeGeoKey =  1;
option.GTRasterTypeGeoKey = 1;
option.GeographicTypeGeoKey = 4217;
option.GeogGeodeticDatumGeoKey = 6217;
option.GeogLinearUnitsGeoKey = 9001;
option.ProjLinearUnitsGeoKey = 9001;
option.VerticalUnitsGeoKey = 9001;
option.ModelTiepointTag = [0;0;0;min(dat.x);min(dat.y);0];
option.ModelPixelScaleTag = [median(diff(dat.x));median(diff(dat.y));0];
option.GTCitationGeoKey = 'Swiss National Coordinates';
option.GeogCitationGeoKey = 'CH1903+';
mygeotiffwrite(fname,[],dat.z,32,option);