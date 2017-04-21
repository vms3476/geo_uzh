

% single tree segmentation using the local maxima method
% for Vordemwald region in Kanton Aargau 



% read LAS data 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/vorwalden/';
area = [633766, 634055, 235882, 236080];
readCIR = 1;

[ las, wsl, cir ] = read_las_wsl_cir( lasDir, area, readCIR );

% plot 
figure; myscatter3(las.x,las.y,las.z,las.z,parula); title('Vordemwald LAS')
figure; myimage(wsl.x,wsl.y,wsl.data); title('Vordemwald WSL forest type'); colormap gray;
figure; imshow(cir); title('Vordemwald CIR')



% local maxima single tree segmentation 
% create CHM raster 
chm = raw2ras([las.x,las.y,las.z],0.5,0.5,'dsm');
[chm.X,chm.Y] = meshgrid(chm.x,chm.y);

% find local maxima 
filtsize = 5;
sig = 0.65;
minheight = 3;
[tree] = locmax(chm.z,filtsize,sig,minheight);

% get starting positions
sx = chm.X(tree);
sy = chm.Y(tree);
figure; myimage(chm.x,chm.y,chm.z);
hold on; swisstick;axis xy;axis equal;axis tight
plot(sx,sy,'.k','markersize',7);
title('CHM raster with local maxima points')

% define loc variable with locations of maxima
clear loc
loc(:,1) = chm.X(tree);
loc(:,2) = chm.Y(tree);
loc(:,3) = chm.z(tree);

% clusters las laser data into single trees 
minHeight = 3; % minimum tree size to be resolved
aspect = 5; % aspect ratio of crowns : dz/mean([dx;dy])
segMethod = 'kme'; % segmentation method 
[x,y,z,idx] = segtree(las,loc,minHeight,aspect,'kme');

% calculate tree parameters values from clustered data
trees = crdata(x,y,z);

% plot tree clusters with random colors
figure; hold on; title('tree segments')
swisstick;axis xy;axis equal;axis tight
colors = [158,1,66; 213,62,79; 244,109,67; 253,174,97; 254,224,139; ...
          230,245,152; 171,221,164; 102,194,165; 50,136,189; 94,79,162] / 255;
nColors = size(colors,1);
for i = 1:size(trees.x,2)
    plot3(x{i},y{i},z{i},'.','color',colors(randi(nColors),:),'markersize',4);
end


%% use matchtrees to compare segmented trees to field data 

% read field data x,y location and height
fieldDataFile = '/Users/scholl/geo_uzh/data/Vorwalden_Field_Data_BHU_Daten_2010_2015.xlsx'; 
fieldData = readtable(fieldDataFile);
treesRef.x = fieldData.X';
treesRef.y = fieldData.Y';
treesRef.h = fieldData.BAHOEHE_M_2015';
% remove the entries with NaN for height
treesRef = subsetraw(treesRef,~isnan(treesRef.h));

% plot ALS-derived tree locations vs. reference data
figure; plot(trees.x,trees.y,'xr',treesRef.x,treesRef.y,'og');
title('Vordemwald Field Data (green) and ALS-derived (red) tree locations')

% determine matches between field data and ALS-segmented trees
[matchedTrees,treeStats] = matchtrees(treesRef,trees,1);


labels = {'Field measured tree height [m]','ALS measured tree height, CHM [m]'};
ax = [0 45 0 45];
figure; 
[modchm02,gof02chm] = myregress(matchedTrees.h,matchedTrees.alsh,...
                      labels,[1 NaN], ax);
                
                  
% plot tree locations over CHM
XX=matchedTrees.alsx(~isnan(matchedTrees.alsx));
YY=matchedTrees.alsy(~isnan(matchedTrees.alsy));
ZZ=matchedTrees.alsh(~isnan(matchedTrees.alsh));
matched_loc=[XX YY ZZ];

figure; myimage(chm.x,chm.y,chm.z);
hold on; swisstick;axis xy;axis equal;axis tight
plot(matchedTrees.x,matchedTrees.y,'ok',matchedTrees.alsx,matchedTrees.alsy,'xr');
title('Vordemwald Field Data (black) and ALS-derived (red) tree locations','FontSize',14)                 



% histogram of height difference between Field Data and ALS trees
dh = matchedTrees.alsh(~isnan(matchedTrees.alsh)) - matchedTrees.h(~isnan(matchedTrees.alsh));
figure;
ii = ~isnan(dh);
hist(dh(ii),100);
title('Field data - ALS tree height Histogram', 'FontSize',14)
xlabel('height difference'); ylabel('count')
                  
%% perform RF forest type classification and display results 
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/vorwalden/rf_classification/';

% compute features for input
%logical([1,1,1,1,0,0,0,0,0]);
Xnew = computeFeaturesRF(las,wsl,1,0,features);
    
% use classifier to predict classes of input pixels
[ypred,unaries] = Mdl.predict(Xnew);
ypred = str2num(cell2mat(ypred));
ypred = reshape(ypred,size(wsl.data));    
    
% morphological processing to reduce noise and create uniform areas
[ypred_withGround] = morphological(ypred, wsl, las, 5, 25);
%wrt_geotiff_CH([outPath tile '_forestTypeALS'],wsl.x,wsl.y,ypred_withGround)
figure; myimage(wsl.x,wsl.y,ypred_withGround); title('Random Forest forest type'); colormap gray;
    







