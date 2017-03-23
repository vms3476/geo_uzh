
% this script performs forest type classification (broadleaf / coniferous)
% using a Random Forest (RF) classifier. 

% The classifier is trained using 3 areas dominated by each forest type. 
% A 1m x 1m grid is used as a sample size for the data. 
% The WSL Forest Type TIF product is used to generate true class labels for
% each sample. 



%% create training data

% load data

lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';


% coniferous 
area1 = [632900, 633415, 234320, 234506]; % south western canton Aargau
area2 = [669252, 669317, 257549, 257647]; % north central canton Aargau
area3 = [666450, 666524, 230868, 230975]; % south eastern canton Aargau

[c.aoi1.las, c.aoi1.wsl, c.aoi1.cir] = read_las_wsl_cir(lasDir, area1);
[c.aoi2.las, c.aoi2.wsl, c.aoi2.cir] = read_las_wsl_cir(lasDir, area2);
[c.aoi3.las, c.aoi3.wsl, c.aoi3.cir] = read_las_wsl_cir(lasDir, area3);



% broadleaf 
area1 = [668290, 668749, 258885, 259071]; % Largeren late March/early April
area2 = [631290, 631350, 237459, 237577]; % March 18 2014 
area3 = [648105, 648230, 234130, 234210]; % March 15

[b.aoi1.las, b.aoi1.wsl, b.aoi1.cir] = read_las_wsl_cir(lasDir, area1);
[b.aoi2.las, b.aoi2.wsl, b.aoi2.cir] = read_las_wsl_cir(lasDir, area2);
[b.aoi3.las, b.aoi3.wsl, b.aoi3.cir] = read_las_wsl_cir(lasDir, area3);

clear area1 area2 area3


%% compute features

% "features" vector indicates which features to compute. first four: 
%       Median first – last height difference per pulse
%       Average intensity of single echoes
%       Fraction of single echoes
%       Median echo height

features = logical([1,1,1,1,0,0,0,0,0]);



% compute features for coniferous training areas
Xc.aoi1 = computeFeaturesRF(c.aoi1.las,c.aoi1.wsl,1,0,features);
Xc.aoi2 = computeFeaturesRF(c.aoi2.las,c.aoi2.wsl,1,0,features);
Xc.aoi3 = computeFeaturesRF(c.aoi3.las,c.aoi3.wsl,1,0,features);

% save only the pixels where WSL forest type is 2
XC.aoi1 = Xc.aoi1(c.aoi1.wsl.data(:)==2,:);
XC.aoi2 = Xc.aoi2(c.aoi2.wsl.data(:)==2,:);
XC.aoi3 = Xc.aoi3(c.aoi3.wsl.data(:)==2,:);

% combine all AOI samples. 
XC.all = [XC.aoi1; XC.aoi2; XC.aoi3];
sum(isnan(XC.all(:))) % check for NaN (should not be any) 



% compute features for broadleaf training areas
Xb.aoi1 = computeFeaturesRF(b.aoi1.las,b.aoi1.wsl,1,0,features);
Xb.aoi2 = computeFeaturesRF(b.aoi2.las,b.aoi2.wsl,1,0,features);
Xb.aoi3 = computeFeaturesRF(b.aoi3.las,b.aoi3.wsl,1,0,features);

% save only the pixels where WSL forest type is 1
XB.aoi1 = Xb.aoi1(b.aoi1.wsl.data(:)==1,:);
XB.aoi2 = Xb.aoi2(b.aoi2.wsl.data(:)==1,:);
XB.aoi3 = Xb.aoi3(b.aoi3.wsl.data(:)==1,:);

% combine all AOI samples. 
XB.all = [XB.aoi1; XB.aoi2; XB.aoi3];
sum(isnan(XB.all(:))) % check for NaN (should not be any) 



% combine the same number of samples per class to train the model
classSamples = min([size(XB.all,1) size(XC.all,1)]);
X = [ XC.all(1:classSamples,:); XB.all(1:classSamples,:) ];

% create Y, set of true class labels, from WSL map
% 1 = broadleaf, 2 = coniferous
Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1)];



%% train model 

% start parallel pool
pool = parpool(4);
paroptions = statset('UseParallel',true);

rng(1); % For reproducibility
Mdl = TreeBagger(400,X,Y,'Method','classification', ...
                         'MinLeafSize',5,'Options',paroptions, ...
                         'OOBPredictorImportance','on');

% assess feature importance 
figure; bar(Mdl.OOBPermutedPredictorDeltaError)

% this step took 3 minutes on the new mac pro 


%% accuracy assessment

stats_pixelwise = accuracyStats(wslTest.data(:),ypred,'noDataValues',NaN);


%% Use classifer on test regions 

% test region 1 
area = [669600, 670000, 258900, 259250];

% calculate features for test data 
[new.las, new.wsl, new.cir] = read_las_wsl_cir(lasDir, area);
Xnew = computeFeaturesRF(new.las,new.wsl,1,0,features);
    
% use classifier to predict classes of input pixels
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));
ypred = reshape(ypred,size(new.wsl.data));

% % plot raw output 
% figure; myimage(new.wsl.x,new.wsl.y,ypred); colormap gray;
% title('Random Forest Classification raw output')

% % For Laegern crown map area: visualize crown polygons, confusion matrix, save plots
% load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');
% outName = [outPath 'test_raw_output_with_ground'];
% crown_map_compare(ypred_withGround,new.wsl,outName,laegernTreeTable_final);

% apply median filter to reduce noise, display output 
medianFilterSize = 5;
ypred_medFilter = medfilt2(ypred,[medianFilterSize,medianFilterSize]);
%figure; myimage(new.wsl.x,new.wsl.y,ypred_medFilter); colormap gray
%title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])

% use BWareaopen to remove areas smaller than 9 pixels 
ypred_open = bwareaopen(ypred_medFilter-1,25);
%figure; myimage(new.wsl.x,new.wsl.y,ypred_open); colormap gray;
   
% apply the same morphological processing to the ground raster (dsm.z < 3m)
dsm = raw2ras([new.las.x,new.las.y,new.las.z],new.wsl,1,'dsm');
ground = zeros(size(ypred)); 
ground(dsm.z<3) = 1; 
ground_medFilter = medfilt2(ground,[medianFilterSize,medianFilterSize]);
ground_open = bwareaopen(ground_medFilter,25);
ypred_withGround = (ypred_open)+1;
ypred_withGround(ground_open==1) = 0; 
figure; myimage(new.wsl.x,new.wsl.y,ypred_withGround); colormap gray;
title('Ground areas added based on cells with a max height < 3m')


%% use classifier on directory of LAS tiles 

wslTif = '/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif';
cirTif = '/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif';

% batch3
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch3/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch3_try2/';

% batch2
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch2/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch2_try2/';

% batch1
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch1/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch1_try2/';

% batch4
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch4/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch4/';


cd(lasDir)
files = dir('*.las');

 for i = 1:numel(files)
    tic % start to measure time to process tile 
    lasName = [files(i).name];
    disp(['Currently processing ' lasName '...']);
    tile = files(i).name(1:end-4);

    % read in LAS data
    las = LASread(lasName,false,false);
    xMin = las.header.min_x; 
    xMax = las.header.max_x; 
    yMin = las.header.min_y; 
    yMax = las.header.max_y; 
    
    % read WSL truth data
    mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
    [wsl.data,wsl.x,wsl.y,wsl.info] = geoimread(wslTif,mapx,mapy);
    wsl.data(wsl.data==3) = 0; 
     
    % compute features for input
    Xnew = computeFeaturesRF(las.record,wsl,0,0,features);
    
    % use classifier to predict classes of input pixels
    [ypred,unaries] = Mdl.predict(Xnew);
    ypred = str2num(cell2mat(ypred));
    ypred = reshape(ypred,size(wsl.data));    
    
    % morphological processing to reduce noise and create uniform areas
    [ypred_withGround, ypred_noPL] = morphological(ypred, wsl, las.record, 5, 25);
 
    
    wrt_geotiff_CH([outPath tile '_forestTypeALS'],wsl.x,wsl.y,ypred_withGround)
    wrt_geotiff_CH([outPath tile '_forestTypeALS_noPL'],wsl.x,wsl.y,ypred_noPL)
    
    toc 
    
end





%% create mask for the coniferous areas, perform single tree segmentation



