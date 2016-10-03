function [x,y,z] = loadmodel(fname,are,res); 
%  function [x,y,z] = loadmodel(fname,are,res); 
%  load dsm/dtm data from cdf file   
  
  if ~strcmp(fname(1),'/')
    fname = [getenv('LOOL_DATAPATH'),'/',fname];
  end
  
  if nargin == 3
    [x,y] = ncread(fname,{'xdist','ydist'});
    ix1 = find(abs(x-are(1)) == min(abs(x-are(1))));
    ix2 = find(abs(x-are(2)) == max(abs(x-are(2))));
    iy1 = find(abs(y-are(3)) == min(abs(y-are(3))));
    iy2 = find(abs(y-are(4)) == max(abs(y-are(4))));
    z = ncread(fname,'z','xdist',[x(ix1),res,x(ix2)],'ydist',[y(iy2),res,y(iy1)])';
    x = x(ix1:res:ix2);
    y = y(iy2:-res:iy1);
  elseif nargin == 2
    z = ncread(fname,'z','xdist',are(1:2),'ydist',are(3:4))';
    x = ncread(fname,'xdist','xdist',are(1:2));
    y = ncread(fname,'ydist','ydist',are(3:4));
  elseif nargin == 1
    z = ncread(fname,'z')';
    x = ncread(fname,'xdist');
    y = ncread(fname,'ydist');
  end

  if nargout == 1
    mod.x = x;
    mod.y = y;
    mod.z = z;
    %try
      %[mod.X,mod.Y] = meshgrid(x,y);
    %end
    x = mod;
  end
  