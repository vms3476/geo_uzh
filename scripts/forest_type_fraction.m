% this script calculates forest type fraction (coniferous, broadleaf,
% ground/non-vegetaion) for a directory of geo tif images (created using
% classify_forest_type.m) compared to the corresponding area within the WSL
% ground truth forest type image. 

tifDir = '/Users/scholl/geo_uzh/data/KantonAargau/output/batch1/';
cd(tifDir)
files = dir('*.tif');

als_con_count = 0; 
als_dec_count = 0;
als_grnd_count = 0;
als_px_count = 0; 


wsl_con_count = 0; 
wsl_dec_count = 0;
wsl_grnd_count = 0;
wsl_px_count = 0; 

for i = 1:numel(files)
    
    disp(' ')
    disp([' Currently processing ' files(i).name '...']);
    
    % read current geotiff
    [im.z, R] = geotiffread(files(i).name);
    [im.x,im.y] = pixcenters(R,R.RasterSize(1),R.RasterSize(2));
    im.x = im.x - 2000000; % subtract for Swiss Grid CH1903 coordinates
    im.y = im.y - 1000000;     
    [im.X,im.Y] = meshgrid(im.x,im.y); 
    pxTotal = numel(im.z);
    imConFraction = sum(im.z(:)==2) / pxTotal; 
    imDecFraction = sum(im.z(:)==1) / pxTotal; 
    imGrndFraction = sum(im.z(:)==0) / pxTotal; 
    disp('ALS')
    disp(['coniferous fraction =  ' num2str(imConFraction) ]);
    disp(['deciduous fraction =  ' num2str(imDecFraction) ]);
    disp(['ground/non forest fraction =  ' num2str(imGrndFraction) ]);
    
    als_con_count = als_con_count + sum(im.z(:)==2)
    als_dec_count = als_dec_count + sum(im.z(:)==1)
    als_grnd_count = als_grnd_count + sum(im.z(:)==0)
    als_px_count = als_px_count + pxTotal
    
    % plot lidar-derived forest map
    fig = figure; subplot(1,2,1); myscatter3(im.X(:), im.Y, im.z(:), im.z(:),gray);
    xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14); title('LiDAR-Derived Forest Type','FontSize',14); view(2);

    % subset WSL image for same area
    xMin = R.XWorldLimits(1)- 2000000; 
    xMax = R.XWorldLimits(2)- 2000000;
    yMin = R.YWorldLimits(1)- 1000000; 
    yMax = R.YWorldLimits(2)- 1000000;    
    mapx = [xMin xMax xMax xMin];
    mapy = [yMin yMin yMax yMax];
    
    [wsl.z,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
    [wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
    wsl.z(wsl.z==3) = 0; % assign pixels with "no data" as non-veg, 0
    pxTotal = numel(wsl.z);
    wslConFraction = sum(wsl.z(:)==2) / pxTotal; 
    wslDecFraction = sum(wsl.z(:)==1) / pxTotal;
    wslGrndFraction = sum(wsl.z(:)==0) / pxTotal; 
    disp('WSL')
    disp(['coniferous fraction =  ' num2str(wslConFraction) ]);
    disp(['deciduous fraction =  ' num2str(wslDecFraction) ]);
    disp(['ground/non forest fraction =  ' num2str(wslGrndFraction) ]);
    
    wsl_con_count = wsl_con_count + sum(wsl.z(:)==2)
    wsl_dec_count = wsl_dec_count + sum(wsl.z(:)==1)
    wsl_grnd_count = wsl_grnd_count + sum(wsl.z(:)==0)
    wsl_px_count = wsl_px_count + pxTotal
    
    % plot WSL forest type map
    subplot(1,2,2); hold on; myscatter3(wsl.X(:),wsl.Y(:),wsl.z(:),wsl.z(:),gray);
    xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14); title('WSL Forest Type','FontSize',14); view(2);
    
    set(gcf,'position',[500 500 2000 1000])
    
    print(fig,['Forest_Type_Comparison_' files(i).name(1:end-4)],'-dpng');
    
    
end

als_con = als_con_count / als_px_count
als_dec = als_dec_count / als_px_count
als_grnd = als_grnd_count / als_px_count

wsl_con = wsl_con_count / wsl_px_count
wsl_dec = wsl_dec_count / wsl_px_count
wsl_grnd = wsl_grnd_count / wsl_px_count