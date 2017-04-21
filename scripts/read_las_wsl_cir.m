function [ las, wsl, cir ] = read_las_wsl_cir( lasDir, area, readCIR )

% reads the LAS , WSL tif, and CIR image for a specified area
% as part of the forest type classification code 
%
% inputs: 
%   lasDIr - string with path specifying location of LAS data
%   area - four element list with the AOI extents
%   readCIR - 0 = do not read CIR image, 1 = DO read CIR image
%
% outputs:
%   las - las point cloud inside specified area
%   wsl - WSL forest type image inside specified area
%   cir - CIR image data inside specified area
%
% sample method call: 
%   
%


% Define file paths and names for WSL and CIR images
wslTif = '/Users/scholl/geo_uzh/data/WSL/Christian-06_10_16/Aargau_Forest_Type_1x1.tif';
cirTif = '/Users/scholl/geo_uzh/data/WSL/AargauOrthoPhotos/Aargau_2013_CIR_1m.tif';


% parse min and max x,y from area vector 
xMin = area(1);
xMax = area(2);
yMin = area(3);
yMax = area(4);

if lasDir(end-3:end) == '.las' % read entire tile
    
    las = readlas(lasDir);

else % if a directory is specified, read the extent of area

    % read LAS data
    cd(lasDir)
    if exist('data.las') == 2
        unix('rm data.las') 
    end
    las = getrawlas(lasDir,area);

end
    
% set LAS points classified as buildings to a height of 0
    las.z( las.Classification == 6 ) = 0;

% read corresponding area of WSL forest type TIF 
mapx = [xMin xMax xMax xMin]; mapy = [yMax yMax yMin yMin];
[wsl.data,wsl.x,wsl.y,wsl.info] = geoimread(wslTif,mapx,mapy);
wsl.data(wsl.data==3) = 0; 


if readCIR
    clear cir
    % read corresponding area of CIR image for reference. concatenate
    % bands ans scale so it's ready to view using imshow
    [cir.data,cir.x,cir.y,cir.info] = geoimread(cirTif,mapx,mapy);
    [cir.X, cir.Y] = meshgrid(cir.x,cir.y);
    cir = cat(3,uint8(cir.data(:,:,1)),uint8(cir.data(:,:,2)),uint8(cir.data(:,:,3)));
else
    cir = 0;
    
end

end

