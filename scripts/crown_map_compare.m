function [ fig, confusionMatrix ] = crown_map_compare( classIm, wslTest, outName, laegernTreeTable_final)
% This function displays the classification image and WSL truth image for
% the Laegern site, region where forest crowns are. 
% The crown polygons are overlaid and colored based on their true forest
% type and the estimated forest type based on the class map.
% A confusion matrix is also generated. Both plots are saved as .PNG
% 
%   Intputs
%       classIm - classification image created using ALS data only
%       wslTest - corresponding WSL image region
%       outName - output filepath and tile, info, etc. for saving the plots
%       laegernTreeTable_final - variable resulting from loading: 
%           load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');


fig = figure; 
subplot(1,2,1);
myimage(wslTest.x,wslTest.y,classIm); colormap gray;
overlay_polygons(laegernTreeTable_final); xlabel('Easting','FontSize',14); ylabel('Northing','FontSize',14);

title(['ALS Forest Type '],'FontSize',14); view(2);

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
print(fig,[outName '_forest_Type_Comparison'],'-dpng');

% confusion matrix 
matrix.speciesTruth = trees.species;
matrix.speciesTruth(matrix.speciesTruth < 20) = 0; % coniferous
matrix.speciesTruth(matrix.speciesTruth > 20) = 1; % deciduous

a = matrix.speciesTruth'; a(a==2) = 0; % ground truth tree type class
b = double(matrix.speciesEstimate'); b(b==2) = 0; % estimated tree type class

c = double(matrix.wslSpeciesEstimate'); c(c==2) = 0; % esimated wsl map type
confusionMatrix = figure; 
plotconfusion(a,b,'ALS Classification',a,c,'WSL Map')

print(confusionMatrix,[outName '_confusion_matrix'],'-dpng');



end

