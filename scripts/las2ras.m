function [ras, las] = las2ras( las, area, res)

% subset based on area 
xMin = area(1);
xMax = area(2); 
yMin = area(3);
yMax = area(4);
las = subsetraw(las,las.x > xMin & las.x < xMax & las.y > yMin & las.y < yMax);


% set height of building points equal to zero, keep track of their indices
% to change back later
las.buildingIndex = zeros(numel(las.x),1);
las.buildingIndex(las.classification==6) = 1;
las.z(las.classification==6) = 0; 



%  std with ground 
dsm_wGround = raw2ras([las.x,las.y,las.z],res,res,'dsm'); 
std_wGround = dsm_wGround.std;

% number of returns per cell with ground
den_wGround = raw2ras([las.x,las.y,las.z],res,res,'den'); 
den_wGround = den_wGround.z;

% remove points below 5m
las = subsetraw(las,las.z > 5 & las.z < 60);

% maximum value per cell, and standard deviation
ras = raw2ras([las.x,las.y,las.z],res,res,'dsm'); 

% minimum value per cell
dmn = raw2ras([las.x,las.y,las.z],res,res,'dmn'); 
ras.dmn = dmn.z;

% median value per cell
dtm = raw2ras([las.x,las.y,las.z],res,res,'dtm'); 
ras.dtm = dtm.z;

% mean value per cell
mean = raw2ras([las.x,las.y,las.z],res,res,'int'); 
ras.int = mean.int;

% number of returns per cell
den = raw2ras([las.x,las.y,las.z],res,res,'den'); 
ras.den = den.z;

% first returns
r1 = subsetraw(las,las.return_number == 1);
rasR1 = raw2ras([r1.x,r1.y,r1.z],res,res,'den'); 
ras.r1 = rasR1.z;

% single returns
r11 = subsetraw(las,las.return_number == 1 & las.number_of_returns == 1);
rasR11 = raw2ras([r11.x,r11.y,r11.z],res,res,'den'); 
ras.r11 = rasR11.z;

% intermediate returns (not first or last) 
ri = subsetraw(las,las.number_of_returns~=las.return_number); 
rasri = raw2ras([ri.x,ri.y,ri.z],res,res,'den');
ras.ri = rasri.z; 

% standard deviation with ground points
ras.std_wGround = std_wGround; 


end