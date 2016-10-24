%% specify input paths for las, dtm files

aoi = [669640 670000 258910 259250]; % Laegern to fit all polygons

%%  read laegern las files

    % 2010
    dir_pc_2010_off = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud/';
    dir_pc_2010_on = '/Users/scholl/geo_uzh/data/Laegeren/2010_on/Pointcloud/';
    dir_dtm_2010 = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/DTM/dtm.mat';

    % 2014 
    dir_pc_2014_off = '/Users/scholl/geo_uzh/data/Laegeren/2014_off/Pointcloud/';
    dir_pc_2014_on = '/Users/scholl/geo_uzh/data/Laegeren/2014_on/Pointcloud/';
    dir_dtm_2014 = '/Users/scholl/geo_uzh/data/Laegeren/2014_off/669000_258000_leaf_off.mat';

% read las point cloud data

    % 2010
    pc_2010.leafoff = getrawlas(dir_pc_2010_off,aoi);
    pc_2010.leafon = getrawlas(dir_pc_2010_on,aoi);

    % 2014
    % since in tile form, convert from laz to las before merging
    cd(dir_raw_off)
    files = dir('*.laz');
    if numel(files) > 0 
        for i = 1:numel(files)
            laz_name = [files(i).name];
            las_name = [laz_name(1:end-1) 's'];
            if exist(las_name, 'file') == 0
                laszip = '/Users/scholl/LAStools/bin/laszip';
                unix([laszip ' -i ' laz_name ' -o ' las_name]);
            end
        end
    end
    cd(dir_raw_on)
    files = dir('*.laz');
    if numel(files) > 0
        for i = 1:numel(files)
            laz_name = [files(i).name];
            las_name = [laz_name(1:end-1) 's'];
            if exist(las_name, 'file') == 0
                laszip = '/Users/scholl/LAStools/bin/laszip';
                unix([laszip ' -i ' laz_name ' -o ' las_name]);
            end
        end
    end
    pc_2014.leafoff = getrawlas(dir_pc_2014_off,aoi);
    pc_2014.leafon = getrawlas(dir_pc_2014_on,aoi);

%% load dtm   
    load(dir_dtm_2010); dtm_2010 = dtm; % 2010
    load(dir_dtm_2014); dtm_2014 = dtm; % 2014

% interpolate dtm to the raw point coordinates, subtract PC - DTM
    
    %2010
    [dtm_2010.X,dtm_2010.Y] = meshgrid(dtm_2010.x,dtm_2010.y); 
    pc_2010.leafoff.tz = pc_2010.leafoff.z;
    pc_2010.leafoff.z = pc_2010.leafoff.tz - ...
        interp2(dtm_2010.X,dtm_2010.Y,dtm_2010.z,pc_2010.leafoff.x,pc_2010.leafoff.y);

    pc_2010.leafon.tz = pc_2010.leafon.z;
    pc_2010.leafon.z = pc_2010.leafon.tz - ...
        interp2(dtm_2010.X,dtm_2010.Y,dtm_2010.z,pc_2010.leafon.x,pc_2010.leafon.y);

    %2014
    [dtm_2014.X,dtm_2014.Y] = meshgrid(dtm_2014.x,dtm_2014.y); 
    pc_2014.leafoff.tz = pc_2014.leafoff.z;
    pc_2014.leafoff.z = pc_2014.leafoff.tz - ...
        interp2(dtm_2014.X,dtm_2014.Y,dtm_2010.z,pc_2014.leafoff.x,pc_2014.leafoff.y);

    pc_2014.leafon.tz = pc_2014.leafon.z;
    pc_2014.leafon.z = pc_2014.leafon.tz - ...
        interp2(dtm_2014.X,dtm_2014.Y,dtm_2014.z,pc_2014.leafon.x,pc_2014.leafon.y);
    
    
% remove noisy points taller then 50m
pc_2010.leafoff = subsetraw(pc_2010.leafoff,(pc_2010.leafoff.z<50));
pc_2010.leafon = subsetraw(pc_2010.leafon,(pc_2010.leafon.z<50));
pc_2014.leafoff = subsetraw(pc_2014.leafoff,(pc_2014.leafoff.z<50));
pc_2014.leafon = subsetraw(pc_2014.leafon,(pc_2014.leafon.z<50));

    

% plot normalized point cloud

% 2010
figure; myscatter3(pc_2010.leafoff.x,pc_2010.leafoff.y,pc_2010.leafoff.z,pc_2010.leafoff.z,parula);
title(['2010 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

figure; myscatter3(pc_2010.leafon.x,pc_2010.leafon.y,pc_2010.leafon.z,pc_2010.leafon.z,parula);
title(['2010 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on


%2014
figure; myscatter3(pc_2014.leafoff.x,pc_2014.leafoff.y,pc_2014.leafoff.z,pc_2014.leafoff.z,parula);
title(['2014 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

figure; myscatter3(pc_2014.leafon.x,pc_2014.leafon.y,pc_2014.leafon.z,pc_2014.leafon.z,parula);
title(['2014 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on
    

% save DMT-normalized, 50m noise-filtered PC as mat files
pc_2010_leafoff = pc_2010.leafoff; pc_2010_leafon = pc_2010.leafon;
pc_2014_leafoff = pc_2014.leafoff; pc_2014_leafon = pc_2014.leafon;
save('/Users/scholl/geo_uzh/output/laegern/pc_2010_leafoff_filtered.mat', 'pc_2010_leafoff');
save('/Users/scholl/geo_uzh/output/laegern/pc_2010_leafon_filtered.mat', 'pc_2010_leafon');
save('/Users/scholl/geo_uzh/output/laegern/pc_2014_leafoff_filtered.mat', 'pc_2014_leafoff');
save('/Users/scholl/geo_uzh/output/laegern/pc_2014_leafon_filtered.mat', 'pc_2014_leafon');
