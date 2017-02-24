function [ras_idx] = raw2ras_index(raw,res,peri,index)

%   compute a grid of LIDAR return indices 
%
%   ras = 2d raster
%   x = array of x-coords
%   y = array of y-coords
%   raw = [x1,y1,z1; x2 y2 z2; ...]
%   res = resolution of grid, can be overloaded with structure res.x and res.y containing grid 
%   peri = perimeter to consider for grid point  
%   index = vector containing the index for each LAS return, same length
%   as the raw

% merge indicies with raw matrix
raw = [raw index];

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
  % cell array to fill with LAS return indices 
  ras_idx = cell(length(y),length(x));
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
          val = rawBox(:,4);
          if ~isempty(val)
            ras_idx{j,i} = val;
%           else
%               ras_idx{j,i} = NaN;
          end
      end
    end
  end
  
end