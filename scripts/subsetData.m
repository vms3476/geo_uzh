function [ output_args ] = subsetData( extent, verbose )
% This function takes a spatial subset of LAS, WSL class TIF, and CIR TIF 
% imagery based on the input extent (xMin, xMax, yMin, yMax).
% Input las files covering the extent must be in the current directory
%
% Input: 
%
%   extent - spatial extent of subset in vector form
%            [xMin xMax yMin yMax]
% 
%   verbose(default: False) - if True, graphs are displayed 

% load and plot input LAS PC 
las = getrawlas(lasDir,area,'broadleaf.las');
figure; myscatter3(las.x,las.y,las.z,las.z,parula); colorbar; view(2); 
title('Normalized Broadleaf AOI PC','FontSize',14)

% load and plot WSL ground truth 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif',mapx,mapy);
[wsl.X, wsl.Y] = meshgrid(wsl.x,wsl.y);
wsl.data(wsl.data==3) = 0; % assign pixels with "no data" as non-veg, 0
figure;myscatter3(wsl.X(:),wsl.Y(:),wsl.data(:),wsl.data(:),gray); view(2); 
title('WSL Broadleaf AOI Forest Type','FontSize',14);

% load and plot CIR image 
[cir.data,cir.x,cir.y,cir.info] = geoimread('/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif',mapx,mapy);
[cir.X, cir.Y] = meshgrid(cir.x,cir.y);
rgb = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
figure; imshow(rgb); title('RGB Broadleaf AOI')




end

