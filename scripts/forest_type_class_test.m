% 1 Dec 2016
% experimenting with the forest type classification thresholds
% treating the WSL map as truth


% ground-normalized LAS data with ground and vegetation classes 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/batch3/';

cd(lasDir)
files = dir('*.laz');
if numel(files) > 0 
    for i = 1:numel(files)
        lazName = [files(i).name];
        lasName = [lazName(1:end-1) 's'];
        if exist(lasName, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' lazName ' -o ' lasName]);
        end
    end
end

% if output from getrawlas already in directory, delete it 
if exist('data.las') == 2
   unix('rm data.las') 
end

xMin = 669600; xMax = 670000; yMin = 258600; yMax = 259250;
area = [xMin xMax yMin yMax]; % laegern crown polygon area, and more southern
las = getrawlas(lasDir,area);

% unix(['/Users/scholl/LAStools/bin/lasinfo -i data.las -cd']) % lasinfo

% plot input PC 
figure;myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); swisstick;
title('Normalized PC','FontSize',14)

% overlay crown polygons 
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');
figure;myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); swisstick;
title('Normalized PC with crown polygons','FontSize',14); 
overlay_polygons(laegernTreeTable_final); swisstick;

% read wsl map portion 
mapx = [xMin xMax xMax xMin];
mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wsl.X(:),wsl.Y(:),wsl.data(:),wsl.data(:),gray); view(2); 
title('WSL Forest Type','FontSize',14); swisstick


% tree type masks 
maskB = zeros(numel(wsl.y),numel(wsl.x)); % broadleaf
maskB(wsl.data == 1) = 1; 

maskC = zeros(numel(wsl.y),numel(wsl.x)); % coniferous
maskC(wsl.data == 2) = 1; 

% % classify points within region of flight line overlap
% % terminal: wine lasoverage.exe -i /Users/scholl/geo_uzh/data/KantonAargau/batch3/data.las -step 2 -o tile_overage.las
% lasOver = readlas('/Users/scholl/geo_uzh/data/KantonAargau/tile_overage.las');
% j = ismember(lasOver.Classification,12);
% lasThin = subsetraw(lasOver,~j); 
% x = lasThin.x';
% y = lasThin.y';
% z = lasThin.z';
% rnnr = lasThin.rnnr';
% classification = double(lasThin.Classification');
% int = lasThin.int';

%% variables of interest for entire tile 
x = las.x';
y = las.y';
z = las.z';
rnnr = las.rnnr';
classification = double(las.Classification');
int = las.int';


% raster for ground classification
ras1mGrnd = raw2ras([x;y;classification]',1,1,'dsm');

figure; myimage(ras1mGrnd.x,ras1mGrnd.y,ras1mGrnd.z);
title('Classification (ground, veg, buildings)','FontSize',14); swisstick


% remove middle returns
i = rnnr == 32 | rnnr == 42 | rnnr == 43 | rnnr == 52 | rnnr == 53 | rnnr == 54 | ...
rnnr == 62 | rnnr == 63 | rnnr == 64 | rnnr == 65 | rnnr == 72 | rnnr == 73 | ...
rnnr == 74 | rnnr == 75 | rnnr == 76 | rnnr == 11; 
x(i) = [];
y(i) = [];
z(i) = [];
rnnr(i) = [];
classification(i) = [];
int(i) = []; 

drnnr = diff(rnnr);
dz = diff(z);
ii = (drnnr == 1 | drnnr == 2 | drnnr == 3 | drnnr == 4 | drnnr == 5 | drnnr == 6) & dz<0; 
xDif = x(ii);
yDif = y(ii); 
zDif = -dz(ii); 
intDif = int(ii); 

% 1m raster of first-last echos 
ras1mDen = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'den'); 
figure;myimage(ras1mDen.x,ras1mDen.y,ras1mDen.z); title('DEN');

% ras1mInt = raw2ras([xDif;yDif;intDif]',ras1mGrnd,1,'int'); 
% figure;myimage(ras1mInt.x,ras1mInt.y,ras1mInt.int); title('INT');
% 
% ras1mDtm = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'dtm'); 
% figure;myimage(ras1mDtm.x,ras1mDtm.y,ras1mDtm.z); title('DTM');
% 
% ras1mDsm = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'dsm'); 
% figure;myimage(ras1mDsm.x,ras1mDsm.y,ras1mDsm.z); title('DSM');

%% boxplots to evalute the impact of higher point density in overlap 
% determine average first-last pulse distance difference per polygon
k = ismember(laegernTreeTable_final.species,[11 14 22 23 29 31 56 59]); 
trees = laegernTreeTable_final(k,:);
n_trees = numel(trees.species); 
stats.idField = trees.idField;
stats.species = trees.species;
for tree = 1:n_trees     
    
    % find raw las points within current polygon
    xpoly = trees.xPoly{tree};
    ypoly = trees.yPoly{tree};
    in = inpolygon(xDif,yDif,xpoly,ypoly);
    zpoly = zDif(in);
    intPoly = intDif(in); % intensity of first returns in polygon
    
    stats.zDif(tree,1) = mean(zpoly); 
    stats.meanInt(tree,1) = mean(intPoly);
end

% % keep only species 11 14 22 23 29 31 56 59, each has 40 or more polygons
% k = ismember(stats.species,[11 14 22 23 29 31 56 59]); 
% stats_plot = subsetraw(stats,k);

% average first-last echo height difference boxplot
figure; boxplot(stats.zDif,stats.species);title('2010 Laegeren first-last echo height difference per pulse','FontSize',14); 
xlabel('Species','FontSize',14); ylabel('average height difference within crown polygon','FontSize',14);
set(gca,'fontsize',14); set(gcf,'position',[500 500 1000 1000])

% mean intensity boxplot
figure; boxplot(stats.meanInt,stats.species);title('2010 Laegeren first echo average intensity per pulse','FontSize',14); 
xlabel('Species','FontSize',14); ylabel('average first echo intensity within crown polygon','FontSize',14);
set(gca,'fontsize',14); set(gcf,'position',[500 500 1000 1000])

%% classification with 2 thresholds 
% height difference threshold classificaiton 
th1 = 0.5;  % percent of pixels within cell
th2 = 15;   % significant height difference, in meters

zTh2 = zDif; 
xTh2 = xDif; 
yTh2 = yDif; 

g = (zTh2 <= th2); % index = not a significant height difference

zTh2(g) = []; 
xTh2(g) = []; 
yTh2(g) = []; 

% compute density raster using only points higher than th2
% using the raster resolution and dimensions from ras1mDen so they can be
% divided for a proportion
ras1mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],ras1mGrnd,1,'den'); 

% divide to find the proportion of pixels 
ras1mproportion = ras1mDenTh2.z ./ ras1mDen.z;

% inpaint the proportion rasters   
ras1mproportion = inpaint_nans(ras1mproportion, 4); 

ras1mClass = ones(size(ras1mDen.z));   % 1 is deciduous
ras1mClass(ras1mproportion < th1) = 0; % 2 is coniferous
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClass); colormap gray;

% morpholgical processing
disp('     Morphological processing...')

    % coniferous bwareaopen
        pix1c = 9;
        conn1c = 4;
    % deciduous bwareaopen
        pix1d = 9;
        conn1d = 4;
    % ground class closing
        pixg = 36;
        conng = 4;
    % structural element for smoothing
        %se = strel('disk',3);
        se = strel('disk',2);

% remove noisy pixels in coniferous regions
ras1mClassMorph = bwareaopen(ras1mClass,pix1c,conn1c);
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorph); colormap gray;

% remove noisy pixels in deciduous regions
ras1mClassMorphInv = bwareaopen(1-ras1mClassMorph,pix1d,conn1d); 
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorphInv); colormap gray;

% conect coniferous regions. coniferous = 0, deciduous = 1
ras1mClassMorphClose = 1-imclose(ras1mClassMorphInv, se);
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorphClose); colormap gray;

% make binary image from ground classification
ras1mGrnd.z(ras1mGrnd.z == 2) = 1; 
ras1mGrnd.z(ras1mGrnd.z == 6) = 1; 
ras1mGrnd.z(ras1mGrnd.z == 3) = 0;    
grnd1mMorph = bwareaopen(ras1mGrnd.z,pixg,conng); 
% smooth the ground mask
grnd1mMorphClose = imclose(grnd1mMorph,se);
grnd1mMorphTotal = grnd1mMorphClose;
% create map of ground, deciduous, and coniferous
class1mMap = ras1mClassMorphClose;
class1mMap(ras1mClassMorphClose == 0) = 2; % coniferous class value of 2
class1mMap(grnd1mMorphTotal == 1) = 0; % ground class value of 0

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,class1mMap); colormap gray; swisstick;
title('ALS Forest Type','FontSize',14)

%% compare to WSL image
wslMaskConiferous = zeros(numel(wsl.x),numel(wsl.y));
wslMaskConiferous(wsl.data==2) = 1;
figure; myimage(wsl.x,wsl.y,wslMaskConiferous);
title('WSL Coniferous forest','FontSize',14); colormap gray;

wslMaskBroadleaf = zeros(numel(wsl.x),numel(wsl.y));
wslMaskBroadleaf(wsl.data==1) = 1;
figure; myimage(wsl.x,wsl.y,wslMaskBroadleaf);
title('WSL Broadleaf forest','FontSize',14); colormap gray;


%% 

% rasters of different types with height difference data 
ras1m = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'int');
ras1mDtm = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'dtm');
ras1mDsm = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'dsm');

% fill gaps
ras1m_fill = inpaint_nans(ras1m.int,4);
ras1mDtm_fill = inpaint_nans(ras1mDtm.z, 4); 
ras1mDsm_fill = inpaint_nans(ras1mDsm.z, 4); 

% plot before inpaint nans
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1m.int); colormap gray;
title('Mean first-last height difference'); colorbar

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mDtm.z); colormap gray;
title('DTM first-last height difference'); colorbar

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mDsm.z); colormap gray;
title('DSM first-last height difference'); colorbar

% plot after inpaint nans
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1m_fill); colormap gray;
title('Mean first-last height difference'); colorbar

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mDtm_fill); colormap gray;
title('DTM first-last height difference'); colorbar

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mDsm_fill); colormap gray;
title('DSM first-last height difference'); colorbar

% boxplot of height difference values within each forest type in map
idx_broadleaf = wslMaskBroadleaf == 1;
ras1mBroadleaf = 
