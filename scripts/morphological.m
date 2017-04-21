function [ypred_withGround] = morphological(ypred, wsl, las, medianFilterSize, bwareaopenSize)

% MORPHOLOGICAL PROCESSING FOR FOREST TYPE IMAGE 
%
% Inputs:
%   yped - forest type classification output raster, reshaped to the
%          proper dimensions of the subtile 
%
%   wsl - struct with "x" and "y" fields that have the easting, 
%         northing coordinates for the prediction image
%
%   las - las struct containing point cloud data to evaluate pixels where
%         maximum height is <3m (classified as ground) 
%
%   medianFilterSize - size of median filter neighborhood to reduce noise
%
%   bwareaopenSize - size of areas (with less than this number of pixels)
%                    to remove in classification image. This is an
%                    "area opening" operation
%
% Outputs:
%
%   ypred_withGround - forest type classification image after median
%                      filter, bwareaopen, and ground identification 
%                      operations have been applied. 
%
% Victoria Scholl, UZH RSL 2017



% apply median filter to reduce noise, display output 
ypred_medFilter = medfilt2(ypred,[medianFilterSize,medianFilterSize]);

% use BWareaopen to remove areas smaller than <bwareaopenSize> pixels 
% then invert to remove small areas for the other forest class
ypred_open = bwareaopen(ypred_medFilter-1,bwareaopenSize);
ypred_open = bwareaopen(abs(ypred_open-1),bwareaopenSize);
   
% apply the same morphological processing to the ground raster (dsm.z < 3m)
dsm = raw2ras([las.x,las.y,las.z],wsl,1,'dsm');
dsm.z(isnan(dsm.z)) = 0;    
ground = zeros(size(ypred)); 
ground(dsm.z<3) = 1; 
ground_medFilter = medfilt2(ground,[medianFilterSize,medianFilterSize]);
ground_open = bwareaopen(ground_medFilter,25);
ypred_withGround = abs(ypred_open-1)+1;
ypred_withGround(ground_open==1) = 0; 



end