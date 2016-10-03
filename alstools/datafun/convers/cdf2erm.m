function [] = cdf2erm(filename);
% function [] = cdf2erm(filename);
% Loads Toposys laserscanning models in CDF  Format 
% and converts ermapper file (binary file + ascii header in XXXX.ers) data format  
%
% filename is XXXXX.cdf
%   
% Felix Morsdorf, RSL Zuerich, Sep. 2003
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.cdf', ...
				    'Select CDF Data File ...');
    filename = [pathname,filename];
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];    
  elseif nargin == 1
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];
  end
  
  dfilename = [filename(1:k)]; 

  % load data 

  [x,y,z] = loadmodel(filename);
 
  dtm.x = x;
  dtm.y = y;
  dtm.z = z;
  dtm2erm(dtm,filename);
  