function wrt_geotiff_CH(fname,x,y,dat)
%    function wrt_geotiff_CH(fname,x,y,dat)
% geotiff write function (Version: 19 January 2012, RL)
resx = median(abs(diff(x)));    % define output pixel spacing
resy = median(abs(diff(y)));
% flip image according orientation
if (y(1) < y(end))
    dat = flipdim(dat, 1);
    y = flipdim(y, 2);
end
% add easting/ northing values for Ch1903+ LV95
if min(x) < 2000000.0;
    a = x(1)+2000000.0;
    b = y(1)+1000000.0;
else
    a = x(1);
    b = y(1);
end

% define output name
out = strcat(fname,'.tif');
%  write_tiff header file option
clear option
option.ModelTiepointTag = [0; 0; 0; a; b; 0];   % define tie point
option.ModelPixelScaleTag = [resx; resy; 0];    % define pixel resolution
option.GTModelTypeGeoKey = 1;                   % 1024 define model type (e.g. projected, geographic, ...)
option.GTCitationGeoKey = 'CH1903+_LV95';       % 1026 documentation
option.GeogEllipsoidGeoKey = 7004;              % 2056 (7004 = Bessel)    
option.ProjectedCSTypeGeoKey = 2056;            % 3072 (Datum ESPG code)

geotiffwrite_function(out, [], dat, 32, option);
end

