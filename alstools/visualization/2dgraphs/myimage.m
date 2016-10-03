function [h] = myimage(x,y,z);
% function [h] = myimage(x,y,z);  

  if nargin == 1 & ~isstruct(x)
    [m,n] = size(x);
    z = x;
    x = 1:m;
    y = 1:n;
  elseif nargin == 1 & isstruct(x)
    dat = x;
    x = dat.x;
    y = dat.y;
    z = dat.z;
  end
  
    
  h = imagesc(x,y,z);
  axis xy;
  axis equal;
  axis tight;
  swisstick

  
