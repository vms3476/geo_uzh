function [ypred_withGround, ypred_noPL] = mophological(ypred, wsl, las, medianFilterSize, bwareaopenSize)

% MORPHOLOGICAL PROCESSING FOR FOREST TYPE IMAGE 

% apply median filter to reduce noise, display output 
ypred_medFilter = medfilt2(ypred,[medianFilterSize,medianFilterSize]);

% use BWareaopen to remove areas smaller than <bwareaopenSize> pixels 
ypred_open = bwareaopen(ypred_medFilter-1,bwareaopenSize);
   
% apply the same morphological processing to the ground raster (dsm.z < 3m)
dsm = raw2ras([las.x,las.y,las.z],wsl,1,'dsm');
ground = zeros(size(ypred)); 
ground(dsm.z<3) = 1; 
ground_medFilter = medfilt2(ground,[medianFilterSize,medianFilterSize]);
ground_open = bwareaopen(ground_medFilter,25);
ypred_withGround = (ypred_open)+1;
ypred_withGround(ground_open==1) = 0; 



% POWER LINE IDENTIFICATION solution for now excludes pixels that are
% likely to contain power lines 

% compute "power line echo" metric (Kim & Sohn 2010), #first / #all returns
canopy = subsetraw(las,las.z > 5 & las.z < 60); 
canopy.z(canopy.classification == 6) = 0;
r1 = subsetraw(canopy,canopy.return_number == 1); % first returns
rasR1 = raw2ras([r1.x,r1.y,r1.z],3,3,'den'); 
den = raw2ras([canopy.x,canopy.y,canopy.z],3,3,'den'); % all returns
powerline_echo = rasR1.z ./ den.z;

% set nans = 0
powerline_nan0 = powerline_echo;
powerline_nan0(isnan(powerline_nan0)) = 0;

% threshold the powerline image
pl_echo_thr = powerline_nan0; pl_thresh = 0.9;
pl_echo_thr(pl_echo_thr<pl_thresh) = 0; pl_echo_thr(isnan(pl_echo_thr)) = 0;
pl_echo_thr(pl_echo_thr>=pl_thresh) = 1; 

% filter power line echo raster 
pl_thr_morph = bwareaopen(pl_echo_thr,25,4);

% upscale resolution of power line echo from 3x3 to 1x1m
[rasR1.X,rasR1.Y] = meshgrid(rasR1.x,rasR1.y);
[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
powerline_interp = interp2(rasR1.X,rasR1.Y,double(pl_thr_morph),dsm.X,dsm.Y);

% threshold the upscaled result for a bianry mask of areas
powerline_interp(powerline_interp>0) = 1;

% set areas to 0 (ground) 
ypred_noPL = ypred_withGround;
ypred_noPL(powerline_interp==1) = 0;

end