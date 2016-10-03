function [succ] = asc2dat(iFile);

  if nargin == 0
    [iFile,pathname] = uigetfile('*.ASC', ...
				    'Select Raw Data .ASC File ...');
    iFile = [pathname,iFile];
    k = findstr('.',iFile)-1;
    oFile = [iFile(1:k)];    
  elseif nargin == 1
    k = findstr('.',iFile)-1;
    oFile = [iFile(1:k)];
  end
  
  
  oFile = [oFile,'.DAT'];
  
  if exist(oFile) ~= 2
    
    iFid = fopen(iFile,'r');
    oFid = fopen(oFile,'w');
    
    [dum,size] = unix(['wc -l ',iFile]);
    size = size(1:strfind(size,iFile)-1);
    nol = str2num(size); % number of lines in ASC file  
    k = 0;l=k;m=k;
    hw = waitbar(0,'Loading huge set of ASCII data ...');
    while 1
      k = k + 1;
      tline = fgetl(iFid);
      %tline = char(fread(iFid,40)');
      if ~ischar(tline), break, end
      %tline(end-1:end) = [];
      dat = str2num(tline);
      if length(dat) == 3
	count = fwrite(oFid,dat,'single');
      else
	l = l+1;
      end
      if k > m
	waitbar(k/nol,hw);
	m = m + 10000;
      end	 
    end
    
    close(hw);
    fclose(iFid);
    fclose(oFid);
    disp(['Number of broken lines :',num2str(l)])
  end
  succ = 1;



