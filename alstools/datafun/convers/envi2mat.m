function [dat] = envi2mat(ifile,ofile);
% function [dat] = envi2mat(ifile,ofile);
% ifile is .hdr file !  
% convert envi binary file to mat file 
% hdr data is put into structure hdr!

% Felix Morsdorf, RSL Zurich, March 2005  
  
  if nargin == 0
    path = '.';
  end
 
  
  if nargin == 0
    [ifile,pathname] = uigetfile('*.hdr', ...
				    'Select Raw Data Header File ...');
    ifile = [pathname,ifile];
    k = findstr('.',ifile)-1;
    ofile = [ifile(1:k)];    
  elseif nargin == 1
    k = findstr('.',ifile)-1;
    ofile = [ifile(1:k)];
  end
  
  dfilename = [ifile(1:k)]; 
  
  % load Header Information
  
  hdr = readhdr(ifile);
  
  % parse byte swap
  if isfield(hdr,'byteorder')
    if hdr.byteorder == 1
      byte_str = 'b';
    elseif hdr.byteorder == 0
      byte_str = 'l';
    end
  else
    byte_str = 'l';
  end
  
  
  % parse data type
  if hdr.datatype == 1
    data_str = 'char';
  elseif hdr.datatype == 2
    data_str = 'short';
  elseif hdr.datatype == 3
    data_str = 'int';
  elseif hdr.datatype == 4
    data_str = 'single';
  elseif hdr.datatype == 5
    data_str = 'double';
  elseif hdr.datatype == 12
    data_str = 'uint16';
  else
    error(['Data type ',num2str(hdr.datatype),' not implemented yet!']) 
  end
    
  % load data from binary
  [fid,message] = fopen(dfilename,'r',byte_str);
  z = fread(fid,data_str);
  fclose(fid);
  
  % reshape data according hdr info
  if hdr.bands == 1
    data = reshape(z,hdr.samples,hdr.lines)';
  elseif hdr.bands > 1 & strcmp(deblankl(hdr.interleave),'bsq')
    data = reshape(z,hdr.samples,hdr.lines,hdr.bands);    
  elseif hdr.bands > 1 & strcmp(deblankl(hdr.interleave),'bil')
    data = uint16(reshape(z,hdr.samples,hdr.bands,hdr.lines));    
  end
  clear z
  [m,n,p] = size(data);
  for i = 1:n
    ndata(:,:,i) = squeeze(data(:,i,:))';
  end
  data = ndata; clear ndata 
  % construct x and y vectors in swiss national coordinates
  x0 = hdr.mapinfo{2}(2);
  y0 = hdr.mapinfo{2}(3);

  x_data = x0+[0:hdr.mapinfo{2}(4):(hdr.samples-1)*hdr.mapinfo{2}(4)];
  y_data = y0-[0:hdr.mapinfo{2}(5):(hdr.lines-1)*hdr.mapinfo{2}(5)];
  if nargout == 0
    eval(['save ',ofile,' x_data y_data data hdr ']);
  else
    dat.Z = data;
    dat.hdr = hdr;
    dat.x = x_data;
    dat.y = y_data;
  end
  


