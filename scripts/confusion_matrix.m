% read the clasification images that cover the polygon area
[img_669_259.z, R] = geotiffread('/Users/scholl/geo_uzh/data/KantonAargau/output/batch1/669000_259000.tif');
[img_669_259.x,img_669_259.y] = pixcenters(R,R.RasterSize(1),R.RasterSize(2));

[img_669_258.z, R] = geotiffread('/Users/scholl/geo_uzh/data/KantonAargau/output/batch1/669000_258000.tif');
[img_669_258.x,img_669_258.y] = pixcenters(R,R.RasterSize(1),R.RasterSize(2));

% subtract to convert coordinate system
img_669_259.x = img_669_259.x - 2000000;
img_669_259.y = img_669_259.y - 1000000;
img_669_258.x = img_669_258.x - 2000000;
img_669_258.y = img_669_258.y - 1000000;

[img_669_259.X,img_669_259.Y] = meshgrid(img_669_259.x,img_669_259.y); 
[img_669_258.X,img_669_258.Y] = meshgrid(img_669_258.x,img_669_258.y); 


% display images on same figure
fig = figure; hold on;
subplot(1,2,1);
myscatter3(img_669_258.X(:), img_669_258.Y, img_669_258.z(:), img_669_258.z(:),gray);
hold on; myscatter3(img_669_259.X(:), img_669_259.Y, img_669_259.z(:), img_669_259.z(:),gray);
% imagesc(img_669_258.x, img_669_258.y,img_669_258.z); colormap gray; set(gca,'YDir','normal')
% imagesc(img_669_259.x, img_669_259.y,img_669_259.z); colormap gray; set(gca,'YDir','normal')
xlim([669600 670000]); ylim([258900 259300]); 
title('LiDAR-Derived Forest Type','FontSize',14); view(2);

% overlay crown polygons 
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');
overlay_polygons(laegernTreeTable_final); xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14);

% restrict limits of plot for a closer look
xlim([669600 670000]); ylim([258900 259300]); 


% read WSL truth map
mapx = [669600 670000 670000 669600];
mapy = [258900 258900 259300 259300];
[subset.data,subset.x,subset.y,subset.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[subset.X, subset.Y] = meshgrid(subset.x,subset.y);
subset.data(subset.data==3) = 0; % assign pixels with "no data" as non-veg, 0
subplot(1,2,2); hold on; 
myscatter3(subset.X(:),subset.Y(:),subset.data(:),subset.data(:),gray);
% imagesc(subset.x, subset.y,subset.data); colormap gray; set(gca,'YDir','normal')

xlim([669600 670000]); ylim([258900 259300]); 
overlay_polygons(laegernTreeTable_final); xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14);
title('WSL Forest Type','FontSize',14); view(2);



%%

% confusion matrix for classification images compared to crown polygon data
% for Laegern study site 

k = ismember(laegernTreeTable_final.species,[11 14 22 23 29 31 56 59]); 
trees = laegernTreeTable_final(k,:);
n_trees = numel(matrix.species); 

for tree = 1:n_trees     
    % find raw las points within current polygon
    xpoly = trees.xPoly{tree};
    ypoly = trees.yPoly{tree};
    
    % record the species, based on most frequent class type within polygon
    % 1 = deciduous, 2 = coniferous, 0 = ground)      
    classPoly = mode([img_669_259.z(inpolygon(img_669_259.X,img_669_259.Y,xpoly',ypoly')) ;...
                      img_669_258.z(inpolygon(img_669_258.X,img_669_258.Y,xpoly',ypoly'))]); 
    
    matrix.speciesEstimate(tree,1) = classPoly;
    
    % determine species estimate from WSL map
    wslPoly = mode(subset.data(inpolygon(subset.X,subset.Y,xpoly',ypoly')));
    matrix.wslSpeciesEstimate(tree,1) = wslPoly; 
    
    % plot to visualize accuracy 
    % on my map
    z = repmat(50,size(xpoly));
    subplot(1,2,1)
    if classPoly == 1 % deciduous estimate
         d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.7,'EdgeAlpha',0);
    elseif classPoly == 2 % coniferous estimate
         d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.7,'EdgeAlpha',0);
    end
    % on wsl map
    subplot(1,2,2)
    if wslPoly == 1 % deciduous estimate
         d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.7,'EdgeAlpha',0);
    elseif wslPoly == 2 % coniferous estimate
         d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.7,'EdgeAlpha',0);
    end
    
end


set(gcf,'position',[500 500 2000 1000])
% save plot - save As .... png
%print(fig,['Forest_Type_Comparison_' files(i).name(1:end-4)],'-dpng');

% confusion matrix 
matrix.speciesTruth = trees.species;
matrix.speciesTruth(matrix.speciesTruth < 20) = 0; % coniferous
matrix.speciesTruth(matrix.speciesTruth > 20) = 1; % deciduous

a = matrix.speciesTruth'; a(a==2) = 0; % ground truth tree type class
b = double(matrix.speciesEstimate'); b(b==2) = 0; % estimated tree type class

c = double(matrix.wslSpeciesEstimate'); c(c==2) = 0; % esimated wsl map type
figure; plotconfusion(a,b,'My Method',a,c,'WSL Map')



%% 

% specify the name of the classification map to compare to the WSL map,
% must have the same dimensions as the WSL map
classIm = ypred_medFilter;
% classIm = pred_ras;

fig = figure; 
subplot(1,2,1);
myimage(wslTest.x,wslTest.y,classIm); colormap gray;
overlay_polygons(laegernTreeTable_final); xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14);

title(['ALS Forest Type, ' num2str(medianFilterSize), 'x',num2str(medianFilterSize) ' Median filter applied'],'FontSize',14); view(2);

subplot(1,2,2);
myimage(wslTest.x,wslTest.y,wslTest.data); colormap gray;
overlay_polygons(laegernTreeTable_final); xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14);
title('WSL Forest Type','FontSize',14); view(2);

[wslTest.X,wslTest.Y] = meshgrid(wslTest.x,wslTest.y);

k = ismember(laegernTreeTable_final.species,[11 14 22 23 29 31 56 59]); 
trees = laegernTreeTable_final(k,:);

matrix.idField = trees.idField;
matrix.species = trees.species;

n_trees = numel(matrix.species); 

for tree = 1:n_trees     
    % find raw las points within current polygon
    xpoly = trees.xPoly{tree};
    ypoly = trees.yPoly{tree};
    
    % record the species, based on most frequent class type within polygon
    % 1 = deciduous, 2 = coniferous, 0 = ground)      
    classPoly = mode(classIm(inpolygon(wslTest.X,wslTest.Y,xpoly',ypoly')));
    
    matrix.speciesEstimate(tree,1) = classPoly;
    
    % determine species estimate from WSL map
    wslPoly = mode(wslTest.data(inpolygon(wslTest.X,wslTest.Y,xpoly',ypoly')));
    matrix.wslSpeciesEstimate(tree,1) = wslPoly; 
    
    % plot to visualize accuracy 
    % on my map
    z = repmat(50,size(xpoly));
    subplot(1,2,1)
    if classPoly == 1 % deciduous estimate
         d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.7,'EdgeAlpha',0);
    elseif classPoly == 2 % coniferous estimate
         d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.7,'EdgeAlpha',0);
    end
    % on wsl map
    subplot(1,2,2)
    if wslPoly == 1 % deciduous estimate
         d = patch(xpoly,ypoly,z,[0.2 0.2 1],'FaceAlpha',0.7,'EdgeAlpha',0);
    elseif wslPoly == 2 % coniferous estimate
         d = patch(xpoly,ypoly,z,[0 0.8 0.4],'FaceAlpha',0.7,'EdgeAlpha',0);
    end
    
end


set(gcf,'position',[500 500 2000 1000])
% save plot - save As .... png
%print(fig,['Forest_Type_Comparison_' files(i).name(1:end-4)],'-dpng');

% confusion matrix 
matrix.speciesTruth = trees.species;
matrix.speciesTruth(matrix.speciesTruth < 20) = 0; % coniferous
matrix.speciesTruth(matrix.speciesTruth > 20) = 1; % deciduous

a = matrix.speciesTruth'; a(a==2) = 0; % ground truth tree type class
b = double(matrix.speciesEstimate'); b(b==2) = 0; % estimated tree type class

c = double(matrix.wslSpeciesEstimate'); c(c==2) = 0; % esimated wsl map type
figure; plotconfusion(a,b,'ALS Classification',a,c,'WSL Map')


