function [hdr] = readhdr(filename);
% function [hdr] = readhdr(filename,workspace);  
% Read in the ASCII Header File (*.hdr) coming with the ENVI image raster  files
% 
% Automatically assigns variables into workspace
%
%   filename  : absolute path of file, may be omitted
%   hdr    : struct containing all variables found in header
%               

% Felix Morsdorf, RSL Zuerich, March 2005 
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.hdr', 'Select .hdr Header File ...');
    filename = [pathname,filename];
  end
  nvars = {'description','samples','lines','bands','header offset','file type', ...
	   'data type','interleave','sensor type','classes','class lookup', ...  
	   'class names','byte order','x start','y start','map info', ...  
	   'wavelength units','band names'};

  fid = fopen(filename);
  
  if fid == -1
    disp('Filename not valid !')
    return
  end
  
  % read in all lines
  k = 0;
  while 1
    k = k + 1;
    line{k} = fgetl(fid);    
    if ~ischar(line{k})
      fclose(fid);
      break;
    end
  end

  % put lines together that belong together (scheiss zeilenumbruch !)
  line(end) = [];
  j = 1;nline = [];
  for k = 2:length(line)
    a(k) = ~isempty(findstr(line{k},'='));
    if a(k)
      j = j + 1;
      nline{j} = line{k};
    else
      nline{j} = [nline{j},line{k}];
    end
  end
  clear a;line = nline;clear nline;
  
  % extract information from lines and put into hdr structure
  for i = 1:length(line)
    a = findstr(line{i},'=');
    if ~isempty(a)
      for j = 1:length(nvars)
        if ~isempty(findstr(lower(line{i}),nvars{j}))
          ii = j;
        end
      end
      nline = line{i}(a+1:end);
      if isempty(str2num(nline))
	eval(['hdr.',deblankl(nvars{ii}),' = nline;'])
      else
	eval(['hdr.',deblankl(nvars{ii}),' = str2num(nline);'])
      end
    end
  end

 % parse some information
 

 str = hdr.mapinfo;
 str(isspace(str)) = [];
 ii = findstr(str,',');
 hdr.mapinfo = [];
 hdr.mapinfo{1} = str(3:ii(1)-1);
 for i = 2:7
   hdr.mapinfo{2}(i-1) = str2num(str(ii(i)+1:ii(i+1)-1));
 end
 hdr.mapinfo{3} = str(ii(8)+1:end-1);
   
