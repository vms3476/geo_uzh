% start measuring time
tic

% convert laz to las 
laz_name = '/Users/scholl/geo_uzh/data/KantonAargau/LeafOff/668000_258000.laz';
las_name = [laz_name(1:end-1) 's'];
aszip = '/Users/scholl/LAStools/bin/laszip';
unix([laszip ' -i ' laz_name ' -o ' las_name]);

%% load las and dtm 
%raw_off_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2010_off/raw_2010_off.las'); tile = 'polygon_area';
inPath = '/Users/scholl/geo_uzh/data/KantonAargau/LeafOff/';
tile = '668000_258000';
raw_off_all = readlas([inPath tile '.las']);
time_readlas = toc

%%
dtmFilepath = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/DTM/dtm.mat';
load(dtmFilepath); 
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y); 
raw_off_all.tz = raw_off_all.z;
% normalize PC to DTM
raw_off_all.z = raw_off_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_off_all.x,raw_off_all.y);
% remove noisy values above 50m
raw_off_all = subsetraw(raw_off_all,raw_off_all.z<50);
time_normalize_dtm = toc

% variables of interest for entire tile 
x = raw_off_all.x';
y = raw_off_all.y';
z = raw_off_all.z';
int = raw_off_all.int';
rnnr = raw_off_all.rnnr';
classification = double(raw_off_all.Classification');
%figure; myscatter3(x,y,z,z,parula);


% remove middle returns
i = rnnr == 32 | rnnr == 42 | rnnr == 43 | rnnr == 52 | rnnr == 53 | rnnr == 54 | ...
    rnnr == 62 | rnnr == 63 | rnnr == 64 | rnnr == 65 | rnnr == 72 | rnnr == 73 | ...
    rnnr == 74 | rnnr == 75 | rnnr == 76 | rnnr == 11; 

x(i) = [];
y(i) = [];
z(i) = [];
int(i) = [];
rnnr(i) = [];
classification(i) = [];

% create raster for classification of ground
ras1mGrnd = raw2ras([x;y;classification]',1,1,'dsm');
time_groundRaster = toc

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
% do the same for intensity and classification information
int1 = int(1:end-1); 
int1(~ii) = [];
time_echodiff = toc


% raster for ground classification 
% ras1mGrnd = raw2ras([xDif;yDif;z1]',1,1,'dsm'); ras1mGrnd.z = inpaint_nans(ras1mGrnd.z,4);
% g = ras1mGrnd.z<3; ras1mGrnd.z(g) = 1; ras1mGrnd.z(~g) = 0; 
% figure; myimage(ras1mGrnd.x,ras1mGrnd.y,ras1mGrnd.z);


% calculate rasters contaning mean height difference ('int') and 
% count the total number of zDif values ('den')
% 0.5m resolution raster
ras1m = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'int'); 
ras1mDen = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'den'); 

% height difference threshold classificaiton 
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
% using the raster resolution and dimensions from ras1mDen so they can be
% divided for a proportion
ras1mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],ras1mGrnd,1,'den'); 

% divide to find the proportion of pixels 
ras1mproportion = ras1mDenTh2.z ./ ras1mDen.z;

% inpaint the proportion rasters
ras1mproportion = inpaint_nans(ras1mproportion, 4); 

% use proportion less than threshold 1 to classify
ras1mClass = ones(size(ras1mDen.z));   % 1 is deciduous
ras1mClass(ras1mproportion < th1) = 0; % 0 is coniferous

time_forestTypeClass = toc

% morphological operation to remove noise and fill regions
% % remove small white pixels in black regions
% pix1 = 36; % area of largest regions to keep
% conn1 = 4;
% ras1mClassMorph = bwareaopen(ras1mClass,pix1,conn1); 
% ras1mClassMorphInv = bwareaopen(1-ras1mClassMorph,pix1,conn1); 
% 
% % dilate to smooth
% se = strel('disk',3);
% ras1mClassMorphDilated = imclose(ras1mClassMorphInv,se);
% % % combine results of both morph processes
% % ras1mClassMorphTotal = ras1mClassMorph | (1-ras1mClassMorphInv); 
% 
% % add ground pixels to class map, first apply the morphological processing to it 
% grnd1mMorph = bwareaopen(ras1mGrnd.z,pix1,conn1); 
% grnd1mMorphInv = bwareaopen(1-grnd1mMorph,pix1,conn1);
% 
% grnd1mMorphTotal = grnd1mMorph | (1-grnd1mMorphInv);
% class1mMap = uint8(ras1mClassMorphTotal);
% class1mMap(class1mMap == 0) = 2; % coniferous class value of 2
% class1mMap(grnd1mMorphTotal == 1) = 0; % ground class value of 0


% remove noisy pixels in coniferous regions
pix1 = 9;
conn1 = 4;
ras1mClassMorph = bwareaopen(ras1mClass,pix1,conn1);
% remove noisy pixels in deciduous regions
pix1 = 9;
conn1 = 4;
ras1mClassMorphInv = bwareaopen(1-ras1mClassMorph,pix1,conn1); 
% conect coniferous regions. coniferous = 0, deciduous = 1
se = strel('disk',3);
ras1mClassMorphClose = 1-imclose(ras1mClassMorphInv, se);
% make binary image
ras1mGrnd.z(ras1mGrnd.z == 2) = 1; 
ras1mGrnd.z(ras1mGrnd.z == 3) = 0;
grnd1mMorph = bwareaopen(ras1mGrnd.z,36,4); 

grnd1mMorphClose = imclose(grnd1mMorph,se);
grnd1mMorphTotal = grnd1mMorphClose;
class1mMap = ras1mClassMorphClose;
class1mMap(ras1mClassMorphClose == 0) = 2; % coniferous class value of 2
class1mMap(grnd1mMorphTotal == 1) = 0; % ground class value of 0

time_morph = toc

% total elapsed  time in seconds
time1tile = toc 

% write geotiff
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/output/'
wrt_geotiff_CH(['/Users/scholl/geo_uzh/data/KantonAargau/output/test' tile],ras1mDenTh2.x,ras1mDenTh2.y,class1mMap)

% plot
%load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat')
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,class1mMap); 
colormap gray; c = colorbar; set(c,'YTickLabel',{'non-forest','','','','','deciduous','','','','','evergreen'},'FontSize',14);
overlay_polygons(laegernTreeTable_final); title(['1m res, th1 = ' num2str(th1) ', th2 = ' num2str(th2) ', bwareaopen pix = ' num2str(pix1) ', conn = ' num2str(conn1)],'FontSize',14)
