% LDA 

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


% non-vegetation
xMin = 656829; xMax = 657482; yMin = 246450; yMax = 246805; 

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
idx_n = (wsl.data(:)==0); % keep only the pixels classified as broadleaf
X_N = X_n(idx_n,:);




% combine same number of samples per class into matrix X, training data. 
% Each row is a sample/observation/pixel. Each column is a feature.
classSamples = min([size(X_B,1) size(X_C,1) size(X_N,1)]);
classSamples = 5000;

X = [ X_C(1:classSamples,:); X_B(1:classSamples,:); X_N(1:classSamples,:); ];

% create Y, set of true class labels, from WSL map
% 0 = nonveg, 1 = broadleaf, 2 = coniferous
Y = [ repmat(2,classSamples,1) ; repmat(1,classSamples,1); repmat(0,classSamples,1) ];



%% train classifier

% start parallel pool
pool = parpool(4);
paroptions = statset('UseParallel',true);

%rng(1); % For reproducibility
Mdl = TreeBagger(400,X,Y,'Method','classification', ...
                         'MinLeafSize',5,'Options',paroptions, ...
                         'OOBPredictorImportance','on');

%% assess feature importance 

figure; 
bar(Mdl.OOBPermutedPredictorDeltaError)



%% test classifier 

%% load test data 

% AOI 
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

%% calculate features for test data 

Xnew = computeFeatures(lasTest, wslTest);

%% use classifier to validate 

%ypred = predict(Mdl,Xnew);
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));



%% plot both the WSL truth map and preduction
figure; 
myscatter3(wslTest.X(:),wslTest.Y(:),wslTest.data(:),wslTest.data(:),gray); view(2); 
overlay_polygons(laegernTreeTable_final);
figure;
myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;
overlay_polygons(laegernTreeTable_final);


%% accuracy assessment
stats_pixelwise = accuracyStats(wslTest.data(:),ypred,'noDataValues',NaN);


%% MRF

% define pairwise interactions in 4-neighborhood
sigma_rbf = 30;     % RBF kernel size
pw = CSPottsPixel(reshape(Xnew,[size(wslTest.data),6]), sigma_rbf, 4);
lambda = 1;         % pairwise strength


figure;
hist(pw(find(pw)));


% labelcost
labelcost = ones(3,'single') - eye(3,'single');


% solve Conditional Random Field (CRF)
[labels, ~, energy] = GCMex(...
        ypred',...
        -log(single(unaries+eps))',...
        lambda .* pw,...
        labelcost,...
        1);

    
figure;
myimage(wslTest.x,wslTest.y,reshape(labels,size(wslTest.data))); colormap gray;





