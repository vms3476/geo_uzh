function [] = erm2cdf(filename,outfilename);
% function [] = erm2cdf(filename);
% Loads Toposys laserscanning data in ERMapper Format (bin-file + .ers)
% and converts to NetCDF data format  
%
% filename is XXXXX.ers - header file
% data file is selected automatically  
%   
% Felix Morsdorf, RSL Zuerich, March 2002
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.ers', ...
				    'Select Raw Data Header File ...');
    filename = [pathname,filename];
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];    
  elseif nargin == 1
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];
  end
  
  dfilename = [filename(1:k)]; 
  
  % load Header Information
  
  head = readers(filename);
  
  % load Data File

  [fid,message] = fopen(dfilename,'r');
  if ~isempty(strfind(head.CellType,'16'))
    z = fread(fid,'uint16');
  else
    z = fread(fid,'float32');
  end
  fclose(fid);
  z = reshape(z,head.NrOfLines,head.NrOfCellsPerLine);
   
  % construct x and y vectors in swiss national coordinates
  
  x0 = head.Eastings;
  y0 = head.Northings;

  x=x0 + ( [0:head.NrOfCellsPerLine-1] * head.Xdimension);
  y=y0 - ( [0:head.NrOfLines-1] * head.Ydimension);

  % set metadata
  
  options.author = 'Felix Morsdorf RSL';
  options.filename = outfilename;
  options.variablename = 'z';
  %options.description = head.SourceDataset;
  options.format = 'float'
  options.d1 = 'xdist';
  options.d2 = 'ydist';
  options.data1 = x;
  options.data2 = y;        
  
  % save to cdf 
  
  save2cdf(z,options);

