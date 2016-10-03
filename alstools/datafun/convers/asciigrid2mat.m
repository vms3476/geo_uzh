function [x,y,dat] = asciigrid2mat(fname);
% function [x,y,dat] = asciigrid2mat(fname);  
%   load data from ESRI ASCII Grid (FARSITE) 
%   and construct x,y vectors 

% Felix Morsdorf, RSL Zurich, Nov 2008
% Changes:
%  11/08: Patched for XLLCORNER and XLLCENTER  
  
  [hdr,data] = hdrload(fname);
  
  strs = {'NCOLS','NROWS','XLLCORNER','YLLCORNER', ...
          'XLLCENTER','YLLCENTER','CELLSIZE','NODATA_VALUE'};
  
  % parse header information
  [m,n] = size(hdr);
  for i = 1:length(strs)
    for j = 1:m
      if strfind(upper(hdr(j,:)),strs{i}) > 0
        ii = strfind(hdr(j,:),' ');
        eval([lower(strs{i}),' = str2num(hdr(j,ii+1:end));']);
      end
    end
  end
 
  % reconstruct data set
  if min(size(data)) == 1
    dat = reshape(data,ncols,nrows)';
  else
    dat = data;
  end
  
  % construct index vectors
  if exist('xllcorner')
    x = xllcorner + [(cellsize/2):cellsize:(ncols)*cellsize];
    y = yllcorner + [(cellsize/2):cellsize:(nrows)*cellsize];y = y(end:-1:1);
  else
    x = xllcenter + [0:cellsize:(ncols-1)*cellsize];
    y = yllcenter + [0:cellsize:(nrows-1)*cellsize];y = y(end:-1:1);
  end
  dat(dat == nodata_value) = NaN;
  if nargout == 1
      dx = x;
      clear x;
      x.x = dx;
      x.y = y;
      x.z = dat;
  end