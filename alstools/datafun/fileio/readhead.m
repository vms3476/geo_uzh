function [header] = readhead(filename);
% function [] = readhead(filename,workspace);  
% Read in the ASCII Header File coming with the large raw data files
% of Toposys Laserscanning Sensor  
% Automatically assigns variables into workspace
%
%   filename  : absolute path of file, may be omitted
%   header    : struct containing all variables found in header
%               
%               
% Felix Morsdorf, RSL Zuerich, March 2002  
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.*', 'Select LDT Header File ...');
    filename = [pathname,filename];
  end
  
  nvars = {'NoVS','N_min','N_max','E_min','E_max','Z_min','Z_max','Z_corr', ...
	   'L_DistMin','L_DistMax','DirFlight'};
  
  fid = fopen(filename);
  
  if fid == -1
    disp('Filename not valid !')
    return
  end
  
  i = 0;
  while 1    
    line = fgetl(fid);    
    if ~ischar(line)
      fclose(fid);
      return
    end
    
    a = findstr(line,'=');
    if ~isempty(a)
      i = i + 1;
      nline = line(a+1:end);
      b = findstr(nline,'[');      
      if ~isempty(b)
        evstr = 'b-1';
      else
	evstr = 'end';
      end      
      eval(['header.',nvars{i},' = str2num(nline(1:',evstr,'));'])
    end
    a = [];b = [];
    c = findstr(line,'System');
    if ~isempty(c)
      vars = line;
      line = fgetl(fid);
      nums = str2num(line);
            
      line = fgetl(fid);
      d = findstr(line,'h');
      vars = [vars,line(1:d)];
      nums = [nums,str2num(line(d+1:end))];
      for j = 1:length(nums)
	[var,vars] = strtok(vars);
	eval(['header.',var,' = nums(j);'])
      end
    end
    c = [];
  end


