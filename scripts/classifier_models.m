% LDA testing 

load fisheriris
gscatter(meas(:,1), meas(:,2), species,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');
N = size(meas,1);

lda = fitcdiscr(meas(:,1:2),species);
ldaClass = resubPredict(lda);


%% Random forest

% define AOIs for training data set. must use same # of samples per class


% coniferous 
xMin = 632900; xMax = 633415; yMin = 234320; yMax = 234506; 

% load and plot input PC 
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
las = getrawlas(lasDir,area);
figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized Coniferous AOI PC','FontSize',14)

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wsl.X(:),wsl.Y(:),wsl.data(:),wsl.data(:),gray); view(2); 
title('WSL Coniferous AOI Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Coniferous AOI')

X_c = computeFeatures(las,wsl); % coniferous input to random forest function
idx_c = (wsl.data(:)==2); % keep only the pixels classified as coniferous
X_C = X_c(idx_c,:);



% broadleaf 

% broadleaf AOI 1 - from Largeren Late March, early April
xMin = 668290; xMax = 668749; yMin = 258885; yMax = 259071; 

% load and plot input PC 
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
las = getrawlas(lasDir,area);
figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized Broadleaf AOI PC','FontSize',14)

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wsl.X(:),wsl.Y(:),wsl.data(:),wsl.data(:),gray); view(2); 
title('WSL Broadleaf AOI Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Broadleaf AOI')

X_b = computeFeatures(las,wsl); % coniferous input to random forest function
idx_b = (wsl.data(:)==1); % keep only the pixels classified as broadleaf
X_B = X_b(idx_b,:);

% % broaflead AOI 2  March 15
xMin = 648112; xMax = 648225; yMin = 234083; yMax = 234147; 

%% non-vegetation
xMin = 658690; xMax = 658999; yMin = 255520; yMax = 255700; 

% load and plot input PC 
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
las = getrawlas(lasDir,area);
figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized Non-Veg AOI PC','FontSize',14)

% WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wsl.X(:),wsl.Y(:),wsl.data(:),wsl.data(:),gray); view(2); 
title('WSL Non-Veg AOI Forest Type','FontSize',14);

% CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Non-Veg AOI')

X_n = computeFeatures(las,wsl); % coniferous input to random forest function
idx_n = (wsl.data(:)==0); % keep only the pixels classified as non-veg
X_N = X_n(idx_n,:);


% combine same number of samples per class into matrix X, training data. 
% Each row is a sample/observation/pixel. Each column is a feature.
classSamples = min([size(X_B,1) size(X_C,1) size(X_N,1)]);
classSamples = 6000;

X = [ X_C(1:classSamples,:); X_B(1:classSamples,:); X_N(1:classSamples,:); ];

% create Y, set of true class labels, from WSL map
% 0 = nonveg, 1 = broadleaf, 2 = coniferous
Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1); repmat(0,classSamples,1) ];


%% 2 CLASSES only

classSamples = min([size(X_B,1) size(X_C,1)]);
classSamples = 6000;
X = [ X_C(1:classSamples,:); X_B(1:classSamples,:) ];

% create Y, set of true class labels, from WSL map
% 0 = nonveg, 1 = broadleaf, 2 = coniferous
Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1)];

%% train classifier

% start parallel pool
pool = parpool(4);
paroptions = statset('UseParallel',true);

%rng(1); % For reproducibility
Mdl = TreeBagger(400,X,Y,'Method','classification', ...
                         'MinLeafSize',5,'Options',paroptions, ...
                         'OOBPredictorImportance','on');

% assess feature importance 

figure; 
bar(Mdl.OOBPermutedPredictorDeltaError)


% test classifier - load test data 

%% AOI 
xMin = 669000; xMax = 670000; yMin = 258880; yMax = 259250;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
% load and plot input PC 
lasTest = getrawlas(lasDir,area);
figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula); colorbar; view(2); 
title('Normalized PC - LAS Test Data','FontSize',14);
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');
overlay_polygons(laegernTreeTable_final);

% load and plot WSL ground truth 
mapx = [xMin xMax xMax xMin];
mapy = [yMax yMax yMin yMin];
[wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wslTest.X(:),wslTest.Y(:),wslTest.data(:),wslTest.data(:),gray); view(2); 
title('WSL Forest Type - Test Area','FontSize',14);
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');
overlay_polygons(laegernTreeTable_final);

% load and plot CIR image 
[cirTest.data,cirTest.x,cirTest.y,cirTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cirTest.X, cirTest.Y] = meshgrid(cirTest.x,cirTest.y);
rgbTest = cat(3,uint8(cirTest.data(:,:,1)),uint8(cirTest.data(:,:,2)),uint8(cirTest.data(:,:,3)));
figure; imshow(rgbTest); 


% calculate features for test data 

Xnew = computeFeatures(lasTest, wslTest);


% use classifier to predict classes of input pixels
%ypred = predict(Mdl,Xnew);
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));

% % add the non-veg pixels (ground) from las data when only 2 ypred classes
% ras = raw2ras([lasTest.x;lasTest.y;lasTest.Classification]',wslTest,1,'dsm');
% figure; myimage(ras.x,ras.y,ras.z);

% plot both the WSL truth map and preduction
figure; 
myscatter3(wslTest.X(:),wslTest.Y(:),wslTest.data(:),wslTest.data(:),gray); view(2); 
overlay_polygons(laegernTreeTable_final);
figure;
myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;
overlay_polygons(laegernTreeTable_final);

%% Filtering the output to reduce noise 

% reshape to proper image dimensions
pred_ras = reshape(ypred,size(wslTest.data));

% apply median filter and display output 
ypred_medFilter = medfilt2(ypred_ras,[10,10]);
figure; myimage(wslTest.x,wslTest.y,ypred_medFilter);colormap gray

%% accuracy assessment
stats_pixelwise = accuracyStats(wslTest.data(:),ypred,'noDataValues',NaN);


%% MRF

% define pairwise interactions in 4-neighborhood
sigma_rbf = 30;     % RBF kernel size. initial = 30;
pw = CSPottsPixel(reshape(Xnew,[size(wslTest.data),6]), sigma_rbf, 4);
lambda = 1;         % pairwise strength. initial = 1; 


figure;
hist(pw(find(pw)));


% labelcost
%labelcost = ones(3,'single') - eye(3,'single');
labelcost = ones(2,'single') - eye(2,'single');

% solve Conditional Random Field (CRF)
[labels, ~, energy] = GCMex(...
        ypred'-1,...   % ypred',...   % for 3 classes
        -log(single(unaries+eps))',...
        lambda .* pw,...
        labelcost,...
        1);

    
figure;
myimage(wslTest.x,wslTest.y,reshape(labels,size(wslTest.data))); colormap gray;




%% Plot class-wise unaries
figure;
subplot(311);
myimage(wslTest.x,wslTest.y,reshape(unaries(:,1),size(wslTest.data)));
subplot(312);
myimage(wslTest.x,wslTest.y,reshape(unaries(:,2),size(wslTest.data)));
subplot(313);
myimage(wslTest.x,wslTest.y,reshape(unaries(:,3),size(wslTest.data)));
colormap gray;


%% another test area for RF classification 
xMin = 656000; xMax = 656250; yMin = 230500; yMax = 230800;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end
% load and plot input PC 
lasTest = getrawlas(lasDir,area);
figure; myscatter3(lasTest.x,lasTest.y,lasTest.z,lasTest.z,parula); colorbar; view(2); 
title('Normalized PC - LAS Test Data','FontSize',14);
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');

% load and plot WSL ground truth 
mapx = [xMin xMax xMax xMin];
mapy = [yMax yMax yMin yMin];
[wslTest.data,wslTest.x,wslTest.y,wslTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wslTest.X, wslTest.Y] = meshgrid(wslTest.x,wslTest.y);
wslTest.data(wslTest.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wslTest.X(:),wslTest.Y(:),wslTest.data(:),wslTest.data(:),gray); view(2); 
title('WSL Forest Type - Test Area','FontSize',14);
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');

% load and plot CIR image 
[cirTest.data,cirTest.x,cirTest.y,cirTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cirTest.X, cirTest.Y] = meshgrid(cirTest.x,cirTest.y);
rgbTest = cat(3,uint8(cirTest.data(:,:,1)),uint8(cirTest.data(:,:,2)),uint8(cirTest.data(:,:,3)));
figure; imshow(rgbTest); 

% calculate features for test data 

Xnew = computeFeatures(lasTest, wslTest);

% use classifier to predict classes of input pixels
%ypred = predict(Mdl,Xnew);
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));

% plot both the WSL truth map and preduction
figure; 
myscatter3(wslTest.X(:),wslTest.Y(:),wslTest.data(:),wslTest.data(:),gray); view(2); 
figure;
myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;

% % add the non-veg pixels (ground) from las data when only 2 ypred classes
% ras = raw2ras([lasTest.x;lasTest.y;lasTest.Classification]',wslTest,1,'dsm');
% figure; myimage(ras.x,ras.y,ras.z);