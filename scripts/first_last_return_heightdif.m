
%% load laegern 2010 leaf off data and dtm over aoi
raw_off_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2010_off/raw_2010_off.las');
dtmFilepath = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/DTM/dtm.mat';
load(dtmFilepath); 
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y); 
raw_off_all.tz = raw_off_all.z;
% normalize PC to DTM
raw_off_all.z = raw_off_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_off_all.x,raw_off_all.y);
% remove noisy values above 50m
raw_off_all = subsetraw(raw_off_all,raw_off_all.z<50);
% plot 
figure; myscatter3(raw_off_all.x,raw_off_all.y,raw_off_all.z,raw_off_all.z,parula); hold on;
title(['2010 Leaf Off PC Normalized']); colorbar; %swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat')
overlay_polygons(laegernTreeTable_final)



%%

% subset that contains deciduous and coniferous trees
xMin = 6.698246099000014e+05; xMax = 6.698622806000002e+05;
yMin = 2.590395454000011e+05; yMax = 2.590693898999989e+05;
j = raw_off_all.x > (xMin) & raw_off_all.x < (xMax) & raw_off_all.y > (yMin) & raw_off_all.y < (yMax);
x = raw_off_all.x(j);
y = raw_off_all.y(j);
z = raw_off_all.z(j);
int = raw_off_all.int(j);
rnnr = raw_off_all.rnnr(j)';
figure; myscatter3(x,y,z,z,parula);

% entire tile 
x = raw_off_all.x';
y = raw_off_all.y';
z = raw_off_all.z';
int = raw_off_all.int';
rnnr = raw_off_all.rnnr';
figure; myscatter3(x,y,z,z,parula);

%% remove middle returns
i = rnnr == 32 | rnnr == 42 | rnnr == 43 | rnnr == 52 | rnnr == 53 | rnnr == 54 | ...
    rnnr == 62 | rnnr == 63 | rnnr == 64 | rnnr == 65 | rnnr == 72 | rnnr == 73 | ...
    rnnr == 74 | rnnr == 75 | rnnr == 76 | rnnr == 11; 

x(i) = [];
y(i) = [];
z(i) = [];
int(i) = [];
rnnr(i) = [];

drnnr = diff(rnnr);
dz = diff(z);

% find indices of first returns for first-last echo pairs
% only keep negative differences
ii = (drnnr == 1 | drnnr == 2 | drnnr == 3 | drnnr == 4 | drnnr == 5 | drnnr == 6) & dz<0; 
xDif = x(ii);
yDif = y(ii); 
zDif = -dz(ii); 
% vector of first echo z values used to mask out non-vegetation areas later
z1 = z(1:end-1);
z1(~ii) = [];
% do the same for intensity information
int1 = int(1:end-1); 
int1(~ii) = [];

% plot first-last return data with crown polygons overlaid
figure; myscatter3(xDif,yDif,zDif,zDif,gray); swisstick;
overlay_polygons(laegernTreeTable_final)
title(['2010 Leaf Off Laegeren'],'FontSize',14); c = colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 round(max(zDif))]); grid on;
ylabel(c,'First-Last Echo Heigh Difference [m]','FontSize',14)

% for testing, vertically stack the return number, diff(return number), 
% indices of first return pairs to keep for analysis, difference values
a = zeros(4,size(x,2),'double'); 
a(1,:) = rnnr;
a(2,:) = [drnnr 0.0];
a(3,:) = [ii 0.0];
a(4,:) = [dz 0.0];

% determine average first-last pulse distance difference per polygon
n_trees = size(laegernTreeTable_final,1);
stats.idField = laegernTreeTable_final.idField;
stats.species = laegernTreeTable_final.species;
for tree = 1:n_trees     
    
    % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{tree};
    ypoly = laegernTreeTable_final.yPoly{tree};
    in = inpolygon(xDif,yDif,xpoly,ypoly);
    zpoly = zDif(in);
    intPoly = int1(in); % intensity of first returns in polygon
    
    stats.zDif(tree,1) = mean(zpoly); 
    stats.meanInt(tree,1) = mean(intPoly);
end

% keep only species 11 14 22 23 29 31 56 59, each has 40 or more polygons
k = ismember(stats.species,[11 14 22 23 29 31 56 59]); 
stats_plot = subsetraw(stats,k);

% average first-last echo height difference boxplot
figure; boxplot(stats_plot.zDif,stats_plot.species);title('2010 Laegeren first-last echo height difference per pulse','FontSize',14); 
xlabel('Species','FontSize',14); ylabel('average height difference within crown polygon','FontSize',14);

% mean intensity boxplot
figure; boxplot(stats_plot.meanInt,stats_plot.species);title('2010 Laegeren first echo average intensity per pulse','FontSize',14); 
xlabel('Species','FontSize',14); ylabel('average first echo intensity within crown polygon','FontSize',14);


%% classification 

% create vegetation mask to avoid classifying ground area, <3m. ground = 1
ras05mGrnd = raw2ras([xDif;yDif;z1]',0.5,1,'dsm'); ras05mGrnd.z = inpaint_nans(ras05mGrnd.z,4);
g = ras05mGrnd.z<3; ras05mGrnd.z(g) = 1; ras05mGrnd.z(~g) = 0; 
figure; myimage(ras05mGrnd.x,ras05mGrnd.y,ras05mGrnd.z);

ras1mGrnd = raw2ras([xDif;yDif;z1]',1,1,'dsm'); ras1mGrnd.z = inpaint_nans(ras1mGrnd.z,4);
g = ras1mGrnd.z<3; ras1mGrnd.z(g) = 1; ras1mGrnd.z(~g) = 0; 
figure; myimage(ras1mGrnd.x,ras1mGrnd.y,ras1mGrnd.z);

ras2mGrnd = raw2ras([xDif;yDif;z1]',2,1,'dsm'); ras2mGrnd.z = inpaint_nans(ras2mGrnd.z,4);
g = ras2mGrnd.z<3; ras2mGrnd.z(g) = 1; ras2mGrnd.z(~g) = 0; 
figure; myimage(ras2mGrnd.x,ras2mGrnd.y,ras2mGrnd.z);

ras3mGrnd = raw2ras([xDif;yDif;z1]',3,1,'dsm'); ras3mGrnd.z = inpaint_nans(ras3mGrnd.z,4);
g = ras3mGrnd.z<3; ras3mGrnd.z(g) = 1; ras3mGrnd.z(~g) = 0; 
figure; myimage(ras3mGrnd.x,ras3mGrnd.y,ras3mGrnd.z);


% calculate rasters contaning mean height difference ('int') and 
% count the total number of zDif values ('den')
% 0.5m resolution raster
ras05m = raw2ras([xDif;yDif;zDif]',0.5,1,'int'); %ras05m.int = inpaint_nans(ras05m.int, 4);
ras05mDen = raw2ras([xDif;yDif;zDif]',0.5,1,'den'); %ras05mDen.z = inpaint_nans(ras05mDen.z, 4);

% 1m res raster
ras1m = raw2ras([xDif;yDif;zDif]',1,1,'int'); %ras1m.int = inpaint_nans(ras1m.int, 4);
ras1mDen = raw2ras([xDif;yDif;zDif]',1,1,'den'); %ras1mDen.z = inpaint_nans(ras1mDen.z, 4);

% 2m res raster
ras2m = raw2ras([xDif;yDif;zDif]',2,1,'int'); %ras2m.int = inpaint_nans(ras2m.int, 4);
ras2mDen = raw2ras([xDif;yDif;zDif]',2,1,'den'); %ras2mDen.z = inpaint_nans(ras2mDen.z, 4);

% 3m res raster
ras3m = raw2ras([xDif;yDif;zDif]',3,1,'int'); %ras3m.int = inpaint_nans(ras3m.int, 4);
ras3mDen = raw2ras([xDif;yDif;zDif]',3,1,'den'); %ras3mDen.z = inpaint_nans(ras3mDen.z, 4);


% plot the height difference rasters
figure; myimage(ras05m.x,ras05m.y,ras05m.int);title('0.5m res raster, mean height difference first-last echo');
c = colorbar; ylabel(c,'average first-last height difference [m]','FontSize',14); axis tight;

figure; myimage(ras1m.x,ras1m.y,ras1m.int);title('1m res raster, mean height difference first-last echo');
c = colorbar; ylabel(c,'average first-last height difference [m]','FontSize',14); axis tight;

figure; myimage(ras2m.x,ras2m.y,ras2m.int);title('2m res raster, mean height difference first-last echo');
c = colorbar; ylabel(c,'average first-last height difference [m]','FontSize',14); axis tight;

figure; myimage(ras3m.x,ras3m.y,ras3m.int);title('3m res raster, mean height difference first-last echo');
c = colorbar; ylabel(c,'average first-last height difference [m]','FontSize',14); axis tight;

%% class maps using single binary threshold
% if mean height is below 18m, then classify as coniferous
th = 15;
class05m = ones(size(ras05m.int));
class05m(ras05m.int < th) = 0; 
figure; myimage(ras05m.x,ras05m.y,class05m); colormap gray; colorbar

class1m = ones(size(ras1m.int));
class1m(ras1m.int < th) = 0; 
figure; myimage(ras1m.x,ras1m.y,class1m); colormap gray; colorbar

class2m = ones(size(ras2m.int));
class2m(ras2m.int < th) = 0; 
figure; myimage(ras2m.x,ras2m.y,class2m); colormap gray; colorbar

class3m = ones(size(ras3m.int));
class3m(ras3m.int < th) = 0; 
figure; myimage(ras3m.x,ras3m.y,class3m); colormap gray; colorbar

ras05mClass = class05m;
ras1mClass = class1m; 
ras2mClass = class2m; 
ras3mClass = class3m; 

%% class maps using 2 thresholds 
th1 = 0.5; % percent of pixels within crown (region? cell?) 
th2 = 15; % significant height difference, in meters

zTh2 = zDif; 
xTh2 = xDif; 
yTh2 = yDif; 

g = (zTh2 <= th2); % index = not a significant height difference

zTh2(g) = []; 
xTh2(g) = []; 
yTh2(g) = []; 

% compute density raster using only points higher than th2
ras05mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],0.5,1,'den'); %ras05mDenTh2.z = inpaint_nans(ras05mDenTh2.z, 4);
ras1mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],1,1,'den'); %ras1mDenTh2.z = inpaint_nans(ras1mDenTh2.z, 4);
ras2mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],2,1,'den'); %ras2mDenTh2.z = inpaint_nans(ras2mDenTh2.z, 4);
ras3mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],3,1,'den'); %ras3mDenTh2.z = inpaint_nans(ras3mDenTh2.z, 4);

% divide to find the proportion of pixels 
ras05mproportion = ras05mDenTh2.z ./ ras05mDen.z;
ras1mproportion = ras1mDenTh2.z ./ ras1mDen.z;
ras2mproportion = ras2mDenTh2.z ./ ras2mDen.z;
ras3mproportion = ras3mDenTh2.z ./ ras3mDen.z;

% inpaint the proportion rasters
ras05mproportion = inpaint_nans(ras05mproportion, 4); 
ras1mproportion = inpaint_nans(ras1mproportion, 4); 
ras2mproportion = inpaint_nans(ras2mproportion, 4); 
ras3mproportion = inpaint_nans(ras3mproportion, 4); 

% use proportion less than threshold 1 to classify
ras05mClass = ones(size(ras05mDen.z));   % 1 is deciduous
ras1mClass = ones(size(ras1mDen.z));
ras2mClass = ones(size(ras2mDen.z));
ras3mClass = ones(size(ras3mDen.z));

ras05mClass(ras05mproportion < th1) = 0; % 0 is coniferous
ras1mClass(ras1mproportion < th1) = 0;
ras2mClass(ras2mproportion < th1) = 0;
ras3mClass(ras3mproportion < th1) = 0;

figure; myimage(ras05mDenTh2.x,ras05mDenTh2.y,ras05mClass); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14); 
overlay_polygons(laegernTreeTable_final); title(['0.5m res, th1 = ' num2str(th1) ', th2 = ', num2str(th2)],'FontSize',14)

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClass); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['1m res, th1 = ' num2str(th1) ', th2 = ', num2str(th2)],'FontSize',14)

figure; myimage(ras2mDenTh2.x,ras2mDenTh2.y,ras2mClass); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['2m res, th1 = ' num2str(th1) ', th2 = ', num2str(th2)],'FontSize',14)

figure; myimage(ras3mDenTh2.x,ras3mDenTh2.y,ras3mClass); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['3m res, th1 = ' num2str(th1) ', th2 = ', num2str(th2)],'FontSize',14)

%% morphological operation to remove noise and fill regions

% remove small white pixels in black regions
pix = 36;
conn = 4;
ras05mClassMorph = bwareaopen(ras05mClass,pix,conn); 
pix = 9;
conn = 4;
ras1mClassMorph = bwareaopen(ras1mClass,pix,conn);
ras2mClassMorph = bwareaopen(ras2mClass,pix,conn);
ras3mClassMorph = bwareaopen(ras3mClass,pix,conn);

% remove small black pixels in white regions
pix05 = 36;
conn05 = 4;
ras05mClassMorphInv = bwareaopen(1-ras05mClassMorph,pix05,conn05); 
pix1 = 9;
conn1 = 4;
ras1mClassMorphInv = bwareaopen(1-ras1mClassMorph,pix1,conn1);
% ras2mClassMorphInv = bwareaopen(1-ras2mClassMorph,pix,conn);
% ras3mClassMorphInv = bwareaopen(1-ras3mClassMorph,pix,conn);

% combine results of both morph processes
ras05mClassMorphTotal = ras05mClassMorph | (1-ras05mClassMorphInv); 
ras1mClassMorphTotal = ras1mClassMorph | (1-ras1mClassMorphInv);
% ras2mClassMorphTotal = ras2mClassMorph | (1-ras2mClassMorphInv);
% ras3mClassMorphTotal = ras3mClassMorph | (1-ras3mClassMorphInv);

% plot class maps after morphological processing
figure; myimage(ras05mDenTh2.x,ras05mDenTh2.y,ras05mClassMorphTotal); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['05m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix05) ', conn = ' num2str(conn05)],'FontSize',14)

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorphTotal); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['1m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix1) ', conn = ' num2str(conn1)],'FontSize',14)

% figure; myimage(ras2mDenTh2.x,ras2mDenTh2.y,ras2mClassMorphTotal); 
% colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
% overlay_polygons(laegernTreeTable_final); title(['2m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix) ', conn = ' num2str(conn)],'FontSize',14)
% 
% figure; myimage(ras3mDenTh2.x,ras3mDenTh2.y,ras3mClassMorphTotal); 
% colormap gray; c = colorbar; set(c,'YTickLabel',{'coniferous','','','','','','','','','','deciduous'},'FontSize',14);
% overlay_polygons(laegernTreeTable_final); title(['3m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix) ', conn = ' num2str(conn)],'FontSize',14)


%% add ground pixels to class map, first apply the morphological processing to it 
grnd05mMorph = bwareaopen(ras05mGrnd.z,pix05,conn05); grnd05mMorphInv = bwareaopen(1-grnd05mMorph,pix05,conn05);
grnd1mMorph = bwareaopen(ras1mGrnd.z,pix1,conn1); grnd1mMorphInv = bwareaopen(1-grnd1mMorph,pix1,conn1);
% grnd2mMorph = bwareaopen(ras2mGrnd.z,pix,conn); grnd2mMorphInv = bwareaopen(1-grnd2mMorph,pix,conn);
% grnd3mMorph = bwareaopen(ras3mGrnd.z,pix,conn); grnd3mMorphInv = bwareaopen(1-grnd3mMorph,pix,conn);

grnd05mMorphTotal = grnd05mMorph | (1-grnd05mMorphInv);
grnd1mMorphTotal = grnd1mMorph | (1-grnd1mMorphInv);
% grnd2mMorphTotal = grnd2mMorph | (1-grnd2mMorphInv);
% grnd3mMorphTotal = grnd3mMorph | (1-grnd3mMorphInv);

class05mMap = uint8(ras05mClassMorphTotal);
class1mMap = uint8(ras1mClassMorphTotal);
% class2mMap = uint8(ras2mClassMorphTotal);
% class3mMap = uint8(ras3mClassMorphTotal);

class05mMap(class05mMap == 0) = 2; % coniferous class value of 2
class1mMap(class1mMap == 0) = 2; 
% class2mMap(class2mMap == 0) = 2;
% class3mMap(class3mMap == 0) = 2;

class05mMap(grnd05mMorphTotal == 1) = 0; % ground class value of 0
class1mMap(grnd1mMorphTotal == 1) = 0;
% class2mMap(grnd2mMorphTotal == 1) = 0;
% class3mMap(grnd3mMorphTotal == 1) = 0;

% plot
figure; myimage(ras05mDenTh2.x,ras05mDenTh2.y,class05mMap); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['05m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix05) ', conn = ' num2str(conn05)],'FontSize',14)

figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,class1mMap); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['1m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix1) ', conn = ' num2str(conn1)],'FontSize',14)

% figure; myimage(ras2mDenTh2.x,ras2mDenTh2.y,class2mMap); 
% colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
% overlay_polygons(laegernTreeTable_final); title(['2m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix) ', conn = ' num2str(conn)],'FontSize',14)
% 
% figure; myimage(ras3mDenTh2.x,ras3mDenTh2.y,class3mMap); 
% colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
% overlay_polygons(laegernTreeTable_final); title(['3m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix) ', conn = ' num2str(conn)],'FontSize',14)




% confusion matrix 

% classify each polygon region based on proportion of pixels that represent
% a significant height difference 

matrix.idField = laegernTreeTable_final.idField;
matrix.species = laegernTreeTable_final.species;

[ras05X, ras05Y] = meshgrid(ras05mDenTh2.x,ras05mDenTh2.y);
[ras1X, ras1Y] = meshgrid(ras1mDenTh2.x,ras1mDenTh2.y);
% [ras2X, ras2Y] = meshgrid(ras2mDenTh2.x,ras2mDenTh2.y);
% [ras3X, ras3Y] = meshgrid(ras3mDenTh2.x,ras3mDenTh2.y);

for tree = 1:n_trees     
    
    % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{tree};
    ypoly = laegernTreeTable_final.yPoly{tree};
    
    classpoly05 = mode(class05mMap(inpolygon(ras05X,ras05Y,xpoly',ypoly')));    
    classpoly1 = mode(class1mMap(inpolygon(ras1X,ras1Y,xpoly',ypoly')));
%     classpoly2 = mode(class2mMap(inpolygon(ras2X,ras2Y,xpoly',ypoly')));
%     classpoly3 = mode(class3mMap(inpolygon(ras3X,ras3Y,xpoly',ypoly')));

    % record the species, based on most frequent class type within polygon
    % 1 = deciduous, 2 = coniferous, 0 = ground) 
    matrix.species05m(tree) = classpoly05;
    matrix.species1m(tree) = classpoly1;
%     matrix.species2m(tree) = classpoly2;
%     matrix.species3m(tree) = classpoly3;
%     
end

%% plot to visualize the accuracy of the different raster resolutions

% 0.5m
figure; myimage(ras05mDenTh2.x,ras05mDenTh2.y,class05mMap); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['05m res, th1 = ', num2str(th1) ', th2 = ' num2str(th2) ',  bwareaopen pix = ' num2str(pix05) ', conn = ' num2str(conn05)],'FontSize',14)
    
for tree = 1:n_trees
    xpoly = laegernTreeTable_final.xPoly{tree};
    ypoly = laegernTreeTable_final.yPoly{tree}; 
    z = repmat(50,size(xpoly));   
    classpoly = matrix.species05m(tree);
    if classpoly == 1 % deciduous estimate
         d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.7,'EdgeAlpha',0);
    elseif classpoly == 2 % coniferous estimate
         d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.7,'EdgeAlpha',0);
    end
end

% 1m
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,class1mMap); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['1m res, th1 = ', num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix1) ', conn = ' num2str(conn1)],'FontSize',14)
for tree = 1:n_trees
    xpoly = laegernTreeTable_final.xPoly{tree};
    ypoly = laegernTreeTable_final.yPoly{tree}; 
    z = repmat(50,size(xpoly));   
    classpoly = matrix.species1m(tree);
    if classpoly == 1 % deciduous estimate
         d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.7,'EdgeAlpha',0);
    elseif classpoly == 2 % coniferous estimate
         d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.7,'EdgeAlpha',0);
    end
end

% % 2m
% figure; myimage(ras2mDenTh2.x,ras2mDenTh2.y,class2mMap); 
% colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
% overlay_polygons(laegernTreeTable_final); title(['2m res, binary th = ' num2str(th) ', bwareaopen pix = ' num2str(pix) ', conn = ' num2str(conn)],'FontSize',14)
% for tree = 1:n_trees
%     xpoly = laegernTreeTable_final.xPoly{tree};
%     ypoly = laegernTreeTable_final.yPoly{tree}; 
%     z = repmat(50,size(xpoly));   
%     classpoly = matrix.species2m(tree);
%     if classpoly == 1 % deciduous estimate
%          d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.5,'EdgeAlpha',0);
%     elseif classpoly == 2 % coniferous estimate
%          d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.5,'EdgeAlpha',0);
%     end
% end
% 
% figure; myimage(ras3mDenTh2.x,ras3mDenTh2.y,class3mMap); 
% colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
% overlay_polygons(laegernTreeTable_final); title(['3m res, binary th = ' num2str(th) ', bwareaopen pix = ' num2str(pix) ', conn = ' num2str(conn)],'FontSize',14)
% for tree = 1:n_trees
%     xpoly = laegernTreeTable_final.xPoly{tree};
%     ypoly = laegernTreeTable_final.yPoly{tree}; 
%     z = repmat(50,size(xpoly));   
%     classpoly = matrix.species3m(tree);
%     if classpoly == 1 % deciduous estimate
%          d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.5,'EdgeAlpha',0);
%     elseif classpoly == 2 % coniferous estimate
%          d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.5,'EdgeAlpha',0);
%     end
% end

% confusion matrix
matrix.speciesTruth = matrix.species; 
matrix.speciesTruth(matrix.speciesTruth < 20) = 2; % coniferous
matrix.speciesTruth(matrix.speciesTruth > 20) = 1; % deciduous

matrix_all = matrix; 

% keep only species 11 14 22 23 29 31 56 59, each has 40 or more polygons
k = ismember(matrix.species,[11 14 22 23 29 31 56 59]); 
matrix = subsetraw(matrix,k);

% make binary 0 = coniferous, 1 = deciduous
a = matrix.speciesTruth'; a(a==2) = 0; % ground truth tree type class
b = double(matrix.species05m'); b(b==2) = 0; % 0.5m
c = double(matrix.species1m'); c(c==2) = 0; % 1m
% d = double(matrix.species2m'); d(d==2) = 0; % 2m
% e = double(matrix.species3m'); e(e==2) = 0; % 3m

figure; plotconfusion(a,b,'0.5m',a,c,'1m'); %,a,d,'2m',a,e,'3m'); 
