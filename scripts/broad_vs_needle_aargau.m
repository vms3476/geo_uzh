base = '/Users/scholl/geo_uzh/';
dataPath = 'data/KantonAargau/LeafOff/';
file = '628000_268000';
fileLAS = [base dataPath file '.las'];
matFile = [file '_leaf_off.mat'];
dtmPath = [base dataPath 'DTM_628000_268000.tif'];
rasStatPath = [base 'data/KantonAargau/Mat_Files_Statistics/raster_statistics/' file '_comb_ras.mat'];
cmap = parula;  % plot colormap


%% convert laz to las before using LASread function 


if exist(fileLAS, 'file') == 0
    laszip = '/Users/scholl/LAStools/bin/laszip';
    fileLAZ = [fileLAS(1:end-1) 'z'];
    unix([laszip ' -i ' fileLAZ ' -o ' fileLAS]);
end

%% read dtm GEOTIFF 
%dtm = GEOTIFF_READ(dtmPath);
[dtm.z, dtmInfo] = geotiffread(dtmPath);
dtm.x = dtmInfo.XWorldLimits(1):dtmInfo.SampleSpacingInWorldX(1):dtmInfo.XWorldLimits(2);
dtm.y = dtmInfo.YWorldLimits(1):dtmInfo.SampleSpacingInWorldY(1):dtmInfo.YWorldLimits(2);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
dtm.z(dtm.z < -5) = NaN;

figure; 
myscatter3(dtm.X(:),dtm.Y(:),dtm.z(:),dtm.z(:),cmap)


%% read mat file for DTM 
mat = load([base dataPath matFile]); data = mat.data; 
dtm.x = data.dtm_x; 
dtm.y = data.dtm_y; 
dtm.z = data.dtm;
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);

figure; 
myscatter3(dtm.X(:),dtm.Y(:),dtm.z(:),dtm.z(:),cmap)
title('Laegeren subset PC'); h = colorbar; xlabel(h,'Height [meters]','FontSize',11);

%% read using mparkan method
pc_mparkan = LASread(fileLAS);

% keep only first returns 
i = ismember(pc_mparkan.record.return_number,1);
pc1.record.x = pc_mparkan.record.x(i);
pc1.record.y = pc_mparkan.record.y(i);
pc1.record.z = pc_mparkan.record.z(i);
pc1.record.intensity = pc_mparkan.record.intensity(i);
pc1.record.return_number = pc_mparkan.record.return_number(i);
%% plot pc
figure; 
myscatter3(pc1.record.x(:),pc1.record.y(:),pc1.record.z(:),pc1.record.z(:),cmap)
title('Aargau entire tile PC first returns'); h = colorbar; xlabel(h,'Height [meters]','FontSize',11);

%% spatial subset
x1 = 628000; 
x2 = 628300; 
y1 = 268000; 
y2 = 268300;
ii = pc1.record.x > (x1) & pc1.record.x < (x2) & pc1.record.y > (y1) & pc1.record.y < (y2);
pc1.record.xs = pc1.record.x(ii);
pc1.record.ys = pc1.record.y(ii);
pc1.record.zs = pc1.record.z(ii);
%% plot subset
figure; 
myscatter3(pc1.record.xs,pc1.record.ys,pc1.record.zs,pc1.record.zs,cmap)
title('Aargau subset PC first returns'); h = colorbar; xlabel(h,'Height [meters]','FontSize',11);

%% dsm 
dsm = raw2ras([pc1.record.xs,pc1.record.ys,pc1.record.zs],0.5,1,'dsm');
[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
figure; imagesc(dsm.x,dsm.y,dsm.z);
axis equal;axis tight;axis xy; colorbar; title('Aargau subset raster first returns');

%% subset dtm
dtm.z = interp2(dtm.X,dtm.Y,dtm.z,dsm.X,dsm.Y);
figure; imagesc(dsm.x,dsm.y,dtm.z);
axis equal;axis tight;axis xy; title('DTM interpolated to same boundaries as dsm');

%% subtract dsm - dtm
dsm.z_norm = dsm.z - dtm.z; 

% normalized first return raster
figure; imagesc(dsm.x,dsm.y,dsm.z_norm);
title('Aargau subset, first return DSM - DTM')
axis equal;axis tight;axis xy;
axis equal;axis tight;axis xy;
caxis([0 50])
hold on
grid on
swisstick
colorbar
colormap([1 1 1;parula])

%% read ground truth from Christian 


