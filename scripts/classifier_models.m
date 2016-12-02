% LDA 

load fisheriris
gscatter(meas(:,1), meas(:,2), species,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');
N = size(meas,1);

lda = fitcdiscr(meas(:,1:2),species);
ldaClass = resubPredict(lda);


%% Random forest

% AOI 
xMin = 669000; xMax = 669950; yMin = 258200; yMax = 258880;
area = [xMin xMax yMin yMax]; 

lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/randomForest/';
cd(lasDir)
if exist('data.las') == 2
   unix('rm data.las') 
end

% load and plot input PC 
las = getrawlas(lasDir,area);
figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized PC','FontSize',14)

% load and plot WSL ground truth 
mapx = [xMin xMax xMax xMin];
mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wsl.X(:),wsl.Y(:),wsl.data(:),wsl.data(:),gray); view(2); 
title('WSL Forest Type','FontSize',14);

% load and plot CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb)

% create Y, set of true class labels, from WSL map
% 0 = nonveg, 1 = broadleaf, 2 = coniferous
Y = wsl.data(:);  

% create X, training data. Each row is a sample/observation/pixel. Each
% column is a feature.

X = computeFeatures(las,wsl); 

%% train classifier

% start parallel pool
pool = parpool(4);
paroptions = statset('UseParallel',true);

rng(1); % For reproducibility
Mdl = TreeBagger(100,X,Y,'Method','classification','MinLeafSize',5,'Options',paroptions);


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
overlay_polygons(laegernTreeTable_final);

% load and plot CIR image 
[cirTest.data,cirTest.x,cirTest.y,cirTest.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cirTest.X, cirTest.Y] = meshgrid(cirTest.x,cirTest.y);
rgbTest = cat(3,uint8(cirTest.data(:,:,1)),uint8(cirTest.data(:,:,2)),uint8(cirTest.data(:,:,3)));
figure; imshow(rgbTest); 

%% calculate features for test data 

Xnew = computeFeatures(lasTest, wslTest);

%% use classifier to validate 

ypred = predict(Mdl,Xnew);
ypred = str2num(cell2mat(ypred));


figure; myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;


%% plot both the WSL truth map and preduction
figure; 
myscatter3(wslTest.X(:),wslTest.Y(:),wslTest.data(:),wslTest.data(:),gray); view(2); 
overlay_polygons(laegernTreeTable_final);
figure;
myimage(wslTest.x,wslTest.y,reshape(ypred,size(wslTest.data))); colormap gray;
overlay_polygons(laegernTreeTable_final);


%% morphological processing 

