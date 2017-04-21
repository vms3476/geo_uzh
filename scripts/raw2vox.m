function [ras] = raw2vox(raw,res,peri,type,binSize);
% modified from raw2ras by felix morsdorf
% to create a 3-D raster, with point counts in vertical direction
    
%    function [ras] = raw2ras(raw,res,peri,type);
%   raw2cls compute a grid of LIDAR raw data returns 
%
%   ras = 2d raster
%   x = array of x-coords
%   y = array of y-coords
%   raw = [x1,y1,z1; x2 y2 z2; ...]
%   res = resolution of grid, can be overloaded with structure res.x and res.y containing grid 
%   peri = perimeter to consider for grid point  
%   type = either 'dsm' or 'dtm' 

  % parse structure res   
  if isstruct(res)  
      x = res.x;
      y = res.y;
      xmin = min(x);
      xmax = max(x);
      ymin = min(y);
      ymax = max(y);
  else
      % Extent of raw data, rounded up/down to the next integer
      xmin = floor(min(raw(:,1)));
      ymin = floor(min(raw(:,2)));
      xmax = ceil(max(raw(:,1)));
      ymax = ceil(max(raw(:,2)));
      
      % Coordinates of new grid
      x = xmin : res : xmax;
      y = ymax : -res : ymin;
  end
  
  % Grid filled with NaNs
  ras(1:length(y),1:length(x)) = NaN;
  ras_std = ras;
  ras_len = ras;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % create empty structure for vertical hist
  if ~exist('binSize','var')
      binSize = 1;
  end
  binEdges = linspace(0,50,(50/binSize)+1);
  vHist = zeros(size(ras,1),size(ras,2),numel(binEdges)-1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  for i = 1:length(x)
    %forget x < i
    raw = raw(raw(:,1) > (x(i) - peri/2),:);
    %width of column x
    rawY = raw(raw(:,1) < (x(i) + peri/2),:);

    for j = 1:length(y)

      %forget y < j
      rawY = rawY(rawY(:,2) < (y(j) + peri/2),:);
      %length of row y
      rawBox = rawY(rawY(:,2) > (y(j) - peri/2),:);
      % pointCloud is now limited to the square with sides of
      % length == perimeter around the cell center point
      
      if ~isempty(rawBox)
        %matlab funktionen, mean, var, max, min, length....
        if length(rawBox(1,:)) == 3
          dat = rawBox(:,3);
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % vertical hist
          datV = dat; % make negative values equal to zero, first bin
          datV(datV < 0) = 0;
          vHist(j,i,:) = histcounts(datV,binEdges);
          
          
          if strcmp(type,'dtm')
            val = median(dat(:)'); 
          elseif strcmp(type,'int')
            val = mean(dat(:)'); 
          elseif strcmp(type,'dsm')
            val = max(dat(:)');
          elseif strcmp(type,'den')
            val = length(dat(:)');
          end
          sdt = std(dat(:)');
          len = length(dat(:)');
          if ~isempty(val)
            ras(j,i) = val;
            ras_std(j,i) = sdt;
            ras_len(j,i) = len;
          else
            ras(j,i) = NaN;
            ras_std(j,i) = NaN;
            ras_len(j,i) = 0;
          end
        else
          ras(j,i) = NaN;
          ras_std(j,i) = NaN;
          ras_len(j,i) = 0;
        end
      else
        ras(j,i) = NaN;
        ras_std(j,i) = NaN;
        ras_len(j,i) = 0;
      end
    end
  end
  dum = ras;
  clear ras;
  ras.x = x;
  ras.y = y;
  if strcmp(type,'int')
    ras.int = dum;
  else
    ras.z = dum;
  end
  ras.std = ras_std;
  ras.num = ras_len;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ras.vHist = vHist;
  ras.vBinEdges = binEdges;
end