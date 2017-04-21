% LDA testing 

load fisheriris
gscatter(meas(:,1), meas(:,2), species,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');
N = size(meas,1);

lda = fitcdiscr(meas(:,1:2),species);
ldaClass = resubPredict(lda);


%% Random forest

inpaintNans = 0; % logical value to inpaint nans during feature 
                 % calculation (1) or set NaN = 0 (0) 

lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';

% define AOIs for training data set. must use same # of samples per class

%%%%% CONIFEROUS %%%%%%%

startTime = tic

%coniferous AOI 3 south eastern canton Aargau
xMin = 666450; xMax = 666524; yMin = 230868; yMax = 230975; 
% load and plot input PC 
area = [xMin xMax yMin yMax]; wkt_polygon(area)
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
las = getrawlas(lasDir,area);
figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized Coniferous AOI PC','FontSize',14)
save('coniferous_las_aoi3','las')
load('coniferous_las_aoi3')

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
title('WSL Coniferous AOI Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Coniferous AOI')

X_c = computeFeatures(las,wsl,1,inpaintNans); % coniferous input to random forest function
idx_c = (wsl.data(:)==2); % keep only the pixels classified as coniferous
X_C3 = X_c(idx_c,:);


% % coniferous AOI 2 north central canton Aargau
% xMin = 646983; xMax = 647047; yMin = 264741; yMax = 264784; 

% coniferous AOI 2 north east canton Aargau, near Laegeren 
xMin = 669252; xMax = 669317; yMin = 257549; yMax = 257647;
% load and plot input PC 
area = [xMin xMax yMin yMax]; wkt_polygon(area)
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% save('coniferous_las_aoi2','las')
load('coniferous_las_aoi2')

figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized Coniferous AOI PC','FontSize',14)

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
title('WSL Coniferous AOI Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Coniferous AOI')

X_c = computeFeatures(las,wsl,1,inpaintNans); % coniferous input to random forest function
idx_c = (wsl.data(:)==2); % keep only the pixels classified as coniferous
X_C2 = X_c(idx_c,:);



% coniferous AOI 1 south western canton Aargau
xMin = 632900; xMax = 633415; yMin = 234320; yMax = 234506; 
% load and plot input PC 
area = [xMin xMax yMin yMax]; wkt_polygon(area)
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
% title('Normalized Coniferous AOI PC','FontSize',14)
% save('coniferous_las_aoi1','las')
load('coniferous_las_aoi1')

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
title('WSL Coniferous AOI Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Coniferous AOI')


X_c = computeFeatures(las,wsl,1,inpaintNans); % coniferous input to random forest function
idx_c = (wsl.data(:)==2); % keep only the pixels classified as coniferous
X_C1 = X_c(idx_c,:);


X_C = [X_C2; X_C3; X_C1];

sum(isnan(X_C(:)))




% %%%% BROADLEAF %%%%%%%

% broadleaf AOI 1 - from Largeren Late March, early April
xMin = 668290; xMax = 668749; yMin = 258885; yMax = 259071; 

% load and plot input PC 
area = [xMin xMax yMin yMax]; wkt_polygon(area)
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
% title('Broadleaf AOI 1 PC','FontSize',14)
% save('broadleaf_las_aoi1','las')
load('broadleaf_las_aoi1')

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
title('WSL Broadleaf AOI 1 Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Broadleaf AOI 1')


X_b = computeFeatures(las,wsl,1,inpaintNans); % coniferous input to random forest function
idx_b = (wsl.data(:)==1); % keep only the pixels classified as broadleaf
X_B1 = X_b(idx_b,:);



% broadleaf AOI 2, March 18 2014 
xMin = 631290; xMax = 631350; yMin = 237459; yMax = 237577;
% load and plot input PC 
area = [xMin xMax yMin yMax]; wkt_polygon(area)
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
% title('Broadleaf AOI 2 PC','FontSize',14)
% save('broadleaf_las_aoi2','las')
load('broadleaf_las_aoi2')

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
title('WSL Broadleaf AOI 2 Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Broadleaf AOI 2')


X_b = computeFeatures(las,wsl,1,inpaintNans); % coniferous input to random forest function
idx_b = (wsl.data(:)==1); % keep only the pixels classified as broadleaf
X_B2 = X_b(idx_b,:); 


% broaflead AOI 3  March 15
xMin = 648105; xMax = 648230; yMin = 234130; yMax = 234210; 
% load and plot input PC 
area = [xMin xMax yMin yMax]; wkt_polygon(area)
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
% title('Broadleaf AOI 3PC','FontSize',14)
% save('broadleaf_las_aoi3','las')
load('broadleaf_las_aoi3')

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
title('WSL Broadleaf AOI 3 Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Broadleaf AOI 3')


X_b = computeFeatures(las,wsl,1,inpaintNans); % coniferous input to random forest function
idx_b = (wsl.data(:)==1); % keep only the pixels classified as broadleaf
X_B3 = X_b(idx_b,:);

X_B = [X_B2; X_B3; X_B1];

sum(isnan(X_B(:)))


endTime = tic

% %% non-vegetation
% 
% % non-veg AOI 1
% xMin = 658690; xMax = 658999; yMin = 255520; yMax = 255700; 
% 
% % load and plot input PC 
% area = [xMin xMax yMin yMax]; 
% wkt_polygon(area)
% lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
% title('Non-Veg AOI 1PC','FontSize',14)
% 
% % WSL ground truth 
% mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
% [wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
% [wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
% wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
% figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
% title('WSL Non-Veg AOI 1 Forest Type','FontSize',14);
% 
% % CIR image 
% [cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
% [cir.X, cir.Y] = meshgrid(cir.x,cir.y);
% rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
% figure; imshow(rgb); title('RGB Non-Veg AOI 1')
% 
% X_n = computeFeatures(las,wsl,1); % coniferous input to random forest function
% idx_n = (wsl.data(:)==0); % keep only the pixels classified as non-veg
% X_N1 = X_n(idx_n,:);
% 
% 
% % non-veg AOI 2 
% xMin = 629700; xMax = 629999; yMin = 266000; yMax = 266294; 
% 
% % load and plot input PC 
% area = [xMin xMax yMin yMax]; 
% wkt_polygon(area)
% lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
% cd(lasDir)
% if exist('data.las') == 2
%    unix('rm data.las') 
% end
% las = getrawlas(lasDir,area);
% figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
% title('Non-Veg AOI 1PC','FontSize',14)
% 
% % WSL ground truth 
% mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
% [wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
% [wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
% wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
% figure;myimage(wsl.x,wsl.y,wsl.data); colormap gray
% title('WSL Non-Veg AOI 1 Forest Type','FontSize',14);
% 
% % CIR image 
% [cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
% [cir.X, cir.Y] = meshgrid(cir.x,cir.y);
% rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
% figure; imshow(rgb); title('RGB Non-Veg AOI 1')
% 
% X_n = computeFeatures(las,wsl,1); % coniferous input to random forest function
% idx_n = (wsl.data(:)==0); % keep only the pixels classified as non-veg
% X_N2 = X_n(idx_n,:);
% 
% 
% X_N = [X_N1; X_N2];
% 
% sum(isnan(X_N(:)))
% 
% 
% 
% 
% 
% %% combine same number of samples per class into matrix X, training data. 
% 
% 
% % 3 classes
% % Each row is a sample/observation/pixel. Each column is a feature.
% classSamples = min([size(X_B,1) size(X_C,1) size(X_N,1)]);
% classSamples = 9000;
% 
% X = [ X_C(1:classSamples,:); X_B(1:classSamples,:); X_N(1:classSamples,:) ];
% % create Y, set of true class labels, from WSL map
% % 0 = nonveg, 1 = broadleaf, 2 = coniferous
% Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1); repmat(0,classSamples,1) ]





%% 2 CLASSES only

classSamples = min([size(X_B,1) size(X_C,1)]);
classSamples = 9000;
X = [ X_C(1:classSamples,:); X_B(1:classSamples,:) ];

% create Y, set of true class labels, from WSL map
% 0 = nonveg, 1 = broadleaf, 2 = coniferous
Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1)];



% % reduce features
% % take out the vertical histogram
features = logical([1,1,1,1,0,0,0,0,0]);

X = X(:,features);


%%%%%%%%%%
% train classifier

% % start parallel pool
% pool = parpool(4);
% paroptions = statset('UseParallel',true);

%rng(1); % For reproducibility
Mdl = TreeBagger(400,X,Y,'Method','classification', ...
                         'MinLeafSize',5,'Options',paroptions, ...
                         'OOBPredictorImportance','on');

% assess feature importance 

figure; 
bar(Mdl.OOBPermutedPredictorDeltaError)



% % accuracy assessment
% stats_pixelwise = accuracyStats(wslTest.data(:),ypred,'noDataValues',NaN);


% %% MRF
% 
% % define pairwise interactions in 4-neighborhood
% sigma_rbf = 30;     % RBF kernel size. initial = 30;
% pw = CSPottsPixel(reshape(Xnew,[size(wslTest.data),6]), sigma_rbf, 4);
% lambda = 1;         % pairwise strength. initial = 1; 
% 
% 
% figure;
% hist(pw(find(pw)));
% 
% 
% % labelcost
% %labelcost = ones(3,'single') - eye(3,'single');
% labelcost = ones(2,'single') - eye(2,'single');
% 
% % solve Conditional Random Field (CRF)
% [labels, ~, energy] = GCMex(...
%         ypred'-1,...   % ypred',...   % for 3 classes
%         -log(single(unaries+eps))',...
%         lambda .* pw,...
%         labelcost,...
%         1);
% 
%     
% figure;
% myimage(wslTest.x,wslTest.y,reshape(labels,size(wslTest.data))); colormap gray;




% % Plot class-wise unaries
% figure;
% subplot(211);
% myimage(wslTest.x,wslTest.y,reshape(unaries(:,1),size(wslTest.data)));
% subplot(212);
% myimage(wslTest.x,wslTest.y,reshape(unaries(:,2),size(wslTest.data)));
% subplot(313);
% colormap gray;


%%% another test area for RF classification 

% % AOI 1
% xMin = 669000; xMax = 670000; yMin = 258880; yMax = 259250;
% xMin = 669500; xMax = 670000; yMin = 258880; yMax = 259250;
% 
% % AOI 2
% xMin = 656000; xMax = 656250; yMin = 230500; yMax = 230800;
% 
% % AOI 3 
% xMin = 647279; xMax = 647479; yMin = 262568; yMax = 262753;

% AOI Tree Crowns 
xMin = 669600; xMax = 670000; yMin = 258900; yMax = 259250;

%%%

area = [xMin xMax yMin yMax];
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
outPath = '/Users/scholl/Desktop/forest_class_plots/'; % for plots

cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
% load and plot input PC 
lasTest = getrawlas(lasDir,area);
figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula); colorbar; view(2); 
title('Normalized PC - LAS Test Data','FontSize',14);

% load and plot WSL ground truth 
mapx = [xMin xMax xMax xMin];
mapy = [yMax yMax yMin yMin];
[wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wslTest.x,wslTest.y,wslTest.data); colormap gray
title('WSL Forest Type - Test Area','FontSize',14);

% load and plot CIR image 
[cirTest.data,cirTest.x,cirTest.y,cirTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cirTest.X, cirTest.Y] = meshgrid(cirTest.x,cirTest.y);
rgbTest = cat(3,uint8(cirTest.data(:,:,1)),uint8(cirTest.data(:,:,2)),uint8(cirTest.data(:,:,3)));
figure; imshow(rgbTest); title('CIR Image of Test Area','FontSize',14);

% calculate features for test data 

% calculate features for test data 
Xnew = computeFeaturesRF(lasTest,wslTest,1,0,features);
    
% use classifier to predict classes of input pixels
%ypred = predict(Mdl,Xnew);
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));

% plot raw output 
figure; myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;
title('Random Forest Classification raw output')

% add ground area 
zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
ypred_withGround = ypred;
ypred_withGround(zMax.z<3) = 0; 
ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
% plot 
figure; myimage(wslTest.x,wslTest.y,ypred_withGround); colormap gray;
title('Ground areas added based on cells with a max height < 3m')

    % compare crown map area, visualize polygons, confusion matrix, save plots
    outName = [outPath tile '_raw_output_with_ground'];
    crown_map_compare(ypred_withGround,wslTest,outName,laegernTreeTable_final);


% apply median filter to reduce noise, display output 
pred_ras = reshape(ypred_withGround,size(wslTest.data));
medianFilterSize = 3;
ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])

    % % compare crown map area, visualize polygons, confusion matrix, save plots
    outName = [outPath tile '_medianFilter_' num2str(medianFilterSize)];
    crown_map_compare(ypred_medFilter,wslTest,outName,laegernTreeTable_final);
    
    
% apply median filter to reduce noise, display output 
pred_ras = reshape(ypred_withGround,size(wslTest.data));
medianFilterSize = 5;
ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])

    % compare crown map area, visualize polygons, confusion matrix, save plots
    outName = [outPath tile '_medianFilter_' num2str(medianFilterSize)];
    crown_map_compare(ypred_medFilter,wslTest,outName,laegernTreeTable_final);


%% AOI 3 
xMin = 647279; xMax = 647479; yMin = 262568; yMax = 262753;

area = [xMin xMax yMin yMax];
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
outPath = '/Users/scholl/Desktop/forest_class_plots/'; % for plots

cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
% load and plot input PC 
lasTest = getrawlas(lasDir,area);
figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula); colorbar; view(2); 
title('Normalized PC - LAS Test Data','FontSize',14);

% load and plot WSL ground truth 
mapx = [xMin xMax xMax xMin];
mapy = [yMax yMax yMin yMin];
[wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myimage(wslTest.x,wslTest.y,wslTest.data); colormap gray
title('WSL Forest Type - Test Area','FontSize',14);

% load and plot CIR image 
[cirTest.data,cirTest.x,cirTest.y,cirTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cirTest.X, cirTest.Y] = meshgrid(cirTest.x,cirTest.y);
rgbTest = cat(3,uint8(cirTest.data(:,:,1)),uint8(cirTest.data(:,:,2)),uint8(cirTest.data(:,:,3)));
figure; imshow(rgbTest); title('CIR Image of Test Area','FontSize',14);

% calculate features for test data 

% calculate features for test data 
Xnew = computeFeaturesRF(lasTest,wslTest,1,0,features);

% use classifier to predict classes of input pixels
%ypred = predict(Mdl,Xnew);
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));

% plot raw output 
figure; myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;
title('Random Forest Classification raw output')

% add ground area 
zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
ypred_withGround = ypred;
ypred_withGround(zMax.z<3) = 0; 
ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
% plot 
figure; myimage(wslTest.x,wslTest.y,ypred_withGround); colormap gray;
title('Ground areas added based on cells with a max height < 3m')


% apply median filter to reduce noise, display output 
pred_ras = reshape(ypred_withGround,size(wslTest.data));
medianFilterSize = 3;
ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])
    
    
% apply median filter to reduce noise, display output 
pred_ras = reshape(ypred_withGround,size(wslTest.data));
medianFilterSize = 5;
ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])


%% AOI from the previous testing site 
 
% for each las tile, classify forest type
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch3/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch3/';
cd(lasDir)
files = dir('*.las');
for i = 1:numel(files)
    lasName = [files(i).name];
    disp(['Currently processing ' lasName '...']);
    tic
    
    tile = files(i).name(1:end-4);
    lasTest = readlas(lasName);
    %figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula); colorbar; view(2); 
    %title('Normalized PC - LAS Test Data','FontSize',14);
    
    % set buildings to ground so they appear black (as non-veg) for classification 
    lasTest.z(lasTest.Classification==6) = 0; 

    % load and plot WSL ground truth 
    xMin = min(lasTest.x);
    xMax = max(lasTest.x);
    yMin = min(lasTest.y);
    yMax = max(lasTest.y);
    mapx = [xMin xMax xMax xMin];
    mapy = [yMax yMax yMin yMin];
    [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
    [wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
    wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0

   % calculate features for test data 
    Xnew = computeFeaturesRF(lasTest,wslTest,0,0,features);

    % use classifier to predict classes of input pixels
    [ypred,unaries] = Mdl.predict(Xnew);
    ypred = str2num(cell2mat(ypred));

    % add ground area 
    zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
    ypred_withGround = ypred;
    ypred_withGround(isnan(zMax.z)) = 0;    
    ypred_withGround(zMax.z<3) = 0; 
    ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
    
    disp('     Writing class map as geotiff...')
    wrt_geotiff_CH([outPath tile '_withGround'],wslTest.x,wslTest.y,reshape(ypred_withGround,size(wslTest.data)))

    % apply median filter to reduce noise, display output 
    pred_ras = reshape(ypred_withGround,size(wslTest.data));
    medianFilterSize = 5;
    ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
    figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
    title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])
    
    disp('     Writing class map as geotiff...')
    wrt_geotiff_CH([outPath tile '_5x5medianFilter'],wslTest.x,wslTest.y,ypred_medFilter)

    toc
end

%% scripting the feature testing for Laegeren site 

%%%%%%%%%%%%%%%%%%%
% Features 1,2,3,4 
features = logical([1 1 1 1 0 0 0 0 0])

% subset the training data 
X = [ X_C(1:classSamples,:); X_B(1:classSamples,:) ];

% create ^set of true class labels, from WSL map 1 = broadleaf, 2 = coniferous
Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1)];

% reduce featurestake out the vertical histogram
X = X(:,features);

% train classifier

% % start parallel pool
% pool = parpool(4);
% paroptions = statset('UseParallel',true);

rng(1); % For reproducibility
Mdl = TreeBagger(400,X,Y,'Method','classification', ...
                         'MinLeafSize',5,'Options',paroptions, ...
                         'OOBPredictorImportance','on');

% assess feature importance 
figure; bar(Mdl.OOBPermutedPredictorDeltaError)




% Laegern Crown Map AOI
xMin = 669600; xMax = 670000; yMin = 258900; yMax = 259250;

% AOI 3
xMin = 647279; xMax = 647479; yMin = 262568; yMax = 262753;
tile = 'aoi3';


    area = [xMin xMax yMin yMax];
    lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
    outPath = '/Users/scholl/Desktop/forest_class_plots/features1234/'; % for plots
    cd(lasDir)
    if exist('data.las') == 2
       unix('rm data.las') 
    end
    % load input PC 
    lasTest = getrawlas(lasDir,area);

    % load WSL ground truth 
    mapx = [xMin xMax xMax xMin];
    mapy = [yMax yMax yMin yMin];
    [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
    wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0

    % calculate features for test data 
    Xnew = computeFeaturesRF(lasTest,wslTest,0,0,features);

    % use classifier to predict classes of input pixels
    [ypred,unaries] = Mdl.predict(Xnew);
    ypred = str2num(cell2mat(ypred));

    % add ground area 
    zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
    ypred_withGround = ypred;
    ypred_withGround(isnan(zMax.z)) = 0;    
    ypred_withGround(zMax.z<3) = 0; 
    ypred_withGround = reshape(ypred_withGround,size(wslTest.data));

    figure; myimage(wslTest.x,wslTest.y,ypred_withGround); colormap gray;
    title('Ground areas added based on cells with a max height < 3m')

        % compare crown map area, visualize polygons, confusion matrix, save plots
        outName = [outPath tile '_raw_output_with_ground'];
        crown_map_compare(ypred_withGround,wslTest,outName,laegernTreeTable_final);

    % apply median filter to reduce noise, display output 
    pred_ras = reshape(ypred_withGround,size(wslTest.data));
    medianFilterSize = 5;
    ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
    figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
    title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])

        % compare crown map area, visualize polygons, confusion matrix, save plots
        outName = [outPath tile '_medianFilter_' num2str(medianFilterSize)];
        crown_map_compare(ypred_medFilter,wslTest,outName,laegernTreeTable_final);
    
    
    
%     
%     
% % AOI 3
% xMin = 647279; xMax = 647479; yMin = 262568; yMax = 262753;
% tile = 'aoi3';
% 
%     area = [xMin xMax yMin yMax];
%     lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
%     outPath = '/Users/scholl/Desktop/forest_class_plots/features1234/'; % for plots
%     cd(lasDir)
%     if exist('data.las') == 2
%        unix('rm data.las') 
%     end
%     % load input PC 
%     lasTest = getrawlas(lasDir,area);
% 
%     % load WSL ground truth 
%     mapx = [xMin xMax xMax xMin];
%     mapy = [yMax yMax yMin yMin];
%     [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
%     wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
% 
%     % calculate features for test data 
%     Xnew = computeFeaturesRF(lasTest,wslTest,1,0,features);
% 
%     % use classifier to predict classes of input pixels
%     [ypred,unaries] = Mdl.predict(Xnew);
%     ypred = str2num(cell2mat(ypred));
% 
%     % add ground area 
%     zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
%     ypred_withGround = ypred;
%     ypred_withGround(zMax.z<3) = 0; 
%     ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
%     figure; myimage(wslTest.x,wslTest.y,ypred_withGround); colormap gray;
%     title('Ground areas added based on cells with a max height < 3m')
%     wrt_geotiff_CH([outPath tile '_withGround'],wslTest.x,wslTest.y,reshape(ypred_withGround,size(wslTest.data)))
% 
%     % apply median filter to reduce noise, display output 
%     pred_ras = reshape(ypred_withGround,size(wslTest.data));
%     ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
%     figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
%     title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])
%     wrt_geotiff_CH([outPath tile '_5x5medianFilter'],wslTest.x,wslTest.y,ypred_medFilter)
%     
% % AOI 2 
% xMin = 656000; xMax = 656250; yMin = 230500; yMax = 230800;
% tile = 'aoi2';
% 
%     area = [xMin xMax yMin yMax];
%     lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
%     outPath = '/Users/scholl/Desktop/forest_class_plots/features1234/'; % for plots
%     cd(lasDir)
%     if exist('data.las') == 2
%        unix('rm data.las') 
%     end
%     % load input PC 
%     lasTest = getrawlas(lasDir,area);
% 
%     % load WSL ground truth 
%     mapx = [xMin xMax xMax xMin];
%     mapy = [yMax yMax yMin yMin];
%     [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
%     wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
% 
%     % calculate features for test data 
%     Xnew = computeFeaturesRF(lasTest,wslTest,1,0,features);
% 
%     % use classifier to predict classes of input pixels
%     [ypred,unaries] = Mdl.predict(Xnew);
%     ypred = str2num(cell2mat(ypred));
% 
%     % add ground area 
%     zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
%     ypred_withGround = ypred;
%     ypred_withGround(zMax.z<3) = 0; 
%     ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
%     figure; myimage(wslTest.x,wslTest.y,ypred_withGround); colormap gray;
%     title('Ground areas added based on cells with a max height < 3m')
%     wrt_geotiff_CH([outPath tile '_withGround'],wslTest.x,wslTest.y,reshape(ypred_withGround,size(wslTest.data)))
% 
%     % apply median filter to reduce noise, display output 
%     pred_ras = reshape(ypred_withGround,size(wslTest.data));
%     ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
%     figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
%     title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])
%     wrt_geotiff_CH([outPath tile '_5x5medianFilter'],wslTest.x,wslTest.y,ypred_medFilter)
%     
% % AOI 4 
% xMin = 655071; xMax = 655131; yMin = 231990; yMax = 232081;
% tile = 'aoi4';
% 
%     area = [xMin xMax yMin yMax];
%     lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
%     outPath = '/Users/scholl/Desktop/forest_class_plots/features1234/'; % for plots
%     cd(lasDir)
%     if exist('data.las') == 2
%        unix('rm data.las') 
%     end
%     % load input PC 
%     lasTest = getrawlas(lasDir,area);
% 
%     % load WSL ground truth 
%     mapx = [xMin xMax xMax xMin];
%     mapy = [yMax yMax yMin yMin];
%     [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
%     wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
% 
%     % calculate features for test data 
%     Xnew = computeFeaturesRF(lasTest,wslTest,1,0,features);
% 
%     % use classifier to predict classes of input pixels
%     [ypred,unaries] = Mdl.predict(Xnew);
%     ypred = str2num(cell2mat(ypred));
% 
%     % add ground area 
%     zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
%     ypred_withGround = ypred;
%     ypred_withGround(zMax.z<3) = 0; 
%     ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
%     figure; myimage(wslTest.x,wslTest.y,ypred_withGround); colormap gray;
%     title('Ground areas added based on cells with a max height < 3m')
%     wrt_geotiff_CH([outPath tile '_withGround'],wslTest.x,wslTest.y,reshape(ypred_withGround,size(wslTest.data)))
% 
%     % apply median filter to reduce noise, display output 
%     pred_ras = reshape(ypred_withGround,size(wslTest.data));
%     ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
%     figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
%     title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])
%     wrt_geotiff_CH([outPath tile '_5x5medianFilter'],wslTest.x,wslTest.y,ypred_medFilter)
    
    
%% For batch 3 area (in Northwestern Aargau), remove power lines and then do forest classification

% for each las tile, classify forest type
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch3/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch3/';
cd(lasDir)
files = dir('*.las');

% power line removal parameters
res = 1; 
stdThresh = 7; 
peri = 6; 
nClusters = 3; 
densityThr = 0.1;

for i = 1:numel(files)
    lasName = [files(i).name];
    disp(['Currently processing ' lasName '...']);
    tic
    
    tile = files(i).name(1:end-4);
    lasTest = readlas(lasName);
    %figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula); colorbar; view(2); 
    %title('Normalized PC - LAS Test Data','FontSize',14);
    
    % set buildings to ground so they appear black (as non-veg) for classification 
    lasTest.z(lasTest.Classification==6) = 0; 

    % load and plot WSL ground truth 
    xMin = min(lasTest.x);
    xMax = max(lasTest.x);
    yMin = min(lasTest.y);
    yMax = max(lasTest.y);
    mapx = [xMin xMax xMax xMin];
    mapy = [yMax yMax yMin yMin];
    [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
    [wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
    wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0

    tic
    lasRaw = lasTest; 
    [lasTestPL, lasTest] = remove_power_lines(lasRaw, res, stdThresh, peri, nClusters, densityThr);
    toc
    
    figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula)
    
   % calculate features for test data 
    Xnew = computeFeaturesRF(lasTest,wslTest,0,0,features);

    % use classifier to predict classes of input pixels
    [ypred,unaries] = Mdl.predict(Xnew);
    ypred = str2num(cell2mat(ypred));

    % add ground area 
    zMax = raw2ras([lasTest.x,lasTest.y,lasTest.z],wslTest,1,'dsm');
    ypred_withGround = ypred;
    ypred_withGround(isnan(zMax.z)) = 0;    
    ypred_withGround(zMax.z<3) = 0; 
    ypred_withGround = reshape(ypred_withGround,size(wslTest.data));
    
    disp('     Writing class map as geotiff...')
    wrt_geotiff_CH([outPath tile '_withGround'],wslTest.x,wslTest.y,reshape(ypred_withGround,size(wslTest.data)))

    % apply median filter to reduce noise, display output 
    pred_ras = reshape(ypred_withGround,size(wslTest.data));
    medianFilterSize = 5;
    ypred_medFilter = medfilt2(pred_ras,[medianFilterSize,medianFilterSize]);
    figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray
    title([num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'])
    
    disp('     Writing class map as geotiff...')
    wrt_geotiff_CH([outPath tile '_5x5medianFilter'],wslTest.x,wslTest.y,ypred_medFilter)

    toc
end

%% write tif for WSL data

lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch1/';
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/batch1/';
cd(lasDir)
files = dir('*.las');

for i = 1:numel(files)
    lasName = [files(i).name];
    tile = files(i).name(1:end-4);
    xMin = str2double(lasName(2:7));
    yMin = str2double(lasName(9:14));
    xMax = xMin + 1000;
    yMax = yMin + 1000;
    mapx = [xMin xMax xMax xMin];
    mapy = [yMax yMax yMin yMin];
    [wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
    [wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
    wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
    wrt_geotiff_CH([outPath tile '_wsl'],wslTest.x,wslTest.y,wslTest.data)
     
end


    
