function [dtm] = xyz2asc(fname); 
%  function [dtm] = xyz2asc(fname);
%  converts xyz grid to ascii grid format
  
  
  % load ascii file
  dat = load(fname);
  
  % assign variables
  x = dat(:,1);y = dat(:,2);z = dat(:,3);
  
  % determine x-spacing and dimensions
  diffx = diff(x);
  
  dx = abs(median(diffx));
  dimx = (max(abs(diffx))/dx)+1;
  
  % determine y-spacing and dimensions
  
  diffy = diff(y);
  dy = max(abs(diffy));
  
  dimy = length(z)/dimx;
  
  % set up new vectors and reshape grid
  
  dtm.z = reshape(z,dimx,dimy);
  dtm.x = min(x):dx:max(x);
  dtm.y = min(y):dy:max(y);
 