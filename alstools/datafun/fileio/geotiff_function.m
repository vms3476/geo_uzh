function geotiff_function(data)
% geotiff write function (Version: 19 January 2012, RL)
resx = abs(data.x(1,1)-data.x(1,2));    % define output pixel spacing
resy = abs(data.y(1,1)-data.y(1,2));
% flip image according orientation
if (data.y(1) < data.y(end))
    data.z = flipdim(data.z, 1);
    data.y = flipdim(data.y, 2);
end
% add easting/ northing values for Ch1903+ LV95
a = data.x(1)+2000000.0;
b = data.y(1)+1000000.0;
% define output name
fname = strcat('dtm_geotiff_prj');
out = strcat(fname,'.tif');
%  write_tiff header file option
clear option
option.ModelTiepointTag = [0; 0; 0; a; b; 0];   % define tie point
option.ModelPixelScaleTag = [resx; resy; 0];    % define pixel resolution
option.GTModelTypeGeoKey = 1;                   % 1024 define model type (e.g. projected, geographic, ...)
option.GTCitationGeoKey = 'CH1903+_LV95';       % 1026 documentation
option.GeogEllipsoidGeoKey = 7004;              % 2056 (7004 = Bessel)    
option.ProjectedCSTypeGeoKey = 2056;            % 3072 (Datum ESPG code)

geotiffwrite_function(out, [], data.z, 32, option);
end

