function [header] = readers(filename);
% function [] = readers(filename,workspace);  
% Read in the ASCII Header File in ERMapper Format 
% of Toposys Laserscanning Sensor  
%
%   filename  : absolute path of file, may be omitted
%   header    : struct containing all variables found in header
%               
%               
% Felix Morsdorf, RSL Zuerich, March 2002  
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.ers', 'Select ers Header File ...');
    filename = [pathname,filename];
  end
  
  fid = fopen(filename);
  
  if fid == -1
    disp('Filename not valid !')
    return
  end
  
  while 1    
    line = fgetl(fid);    
    if ~ischar(line)
      fclose(fid);
      break
    end
    
    a = findstr(line,'=');
    if ~isempty(a)
      aline = line(1:a-1);
      nline = line(a+1:end);
      if isempty(str2num(nline))
	eval(['header.',deblank(aline),' = deblank(nline);'])
      else
        eval(['header.',deblank(aline),' = str2num(nline);'])
      end      
    end
    a = [];
  end


