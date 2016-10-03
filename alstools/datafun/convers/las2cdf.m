function [] = las2cdf(filename,outfilename);
% function [] = las2cdf(filename);
% Loads Optech laserscanning raw data (*.las) and converts it to
% NetCDF data format, which is indexed  
%
% Felix Morsdorf, RSL Zuerich, Dec 2006
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.las', ...
				    'Select Raw Data File ...');
    filename = [pathname,filename];
  elseif nargin ==1
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];
  end
  
  k = findstr('.',filename)-1;
  
  tfilename = [filename(1:k),'.mat'];
  ofilename = [filename(1:k),'.cdf'];  
  
  if exist(ofilename) ~= 2  
    if exist(tfilename) ~= 2  
	 [raw,hdr] = readlas(filename);
	 eval(['save ',tfilename])
       else
	 eval(['load ',tfilename])
       end
  end
  % get header information
  % sort data into bins        
  data = data(:,[2,1,3]);
  data = data';
  nelem = 12700; % xyz triples per bin
  ii = 1:nelem;
  
  nbin = fix(length(data)/nelem); 
  
  nle = nbin*nelem;
  rest = data(:,nle+1:end);
  x = rest(1,:);
  y = rest(2,:);
  z = rest(3,:);
  
  data = reshape(data(:,1:nle),3,nelem,nbin);
  minx = squeeze(minnan(data(1,:,:))); 
  maxx = squeeze(maxnan(data(1,:,:)));
  miny = squeeze(minnan(data(2,:,:))); 
  maxy = squeeze(maxnan(data(2,:,:)));
  minz = squeeze(minnan(data(3,:,:))); 
  maxz = squeeze(maxnan(data(3,:,:)));
    
  bindat = [[1:length(minx)]',minx,maxx,miny,maxy,minz,maxz]';
       
  % take care of the rest
  
  bindat(:,nbin+1) = [nbin+1;minnan(x);maxnan(x);minnan(y);maxnan(y); ...
                      minnan(z);maxnan(z)];
  ex = ones(1,nelem)*NaN;ey=ex;ez=ex;
  
  ii = find(minnan(bindat(2,:)) == bindat(2,:));
  head.X1 = minnan(bindat(2,:));
  head.Y1 = bindat(4,ii(1));
  ii = find(minnan(bindat(4,:)) == bindat(4,:));	     
  head.X2 = bindat(3,ii(1));
  head.Y2 = minnan(bindat(4,:));
  ii = find(maxnan(bindat(5,:)) == bindat(5,:));
  head.X3 = bindat(2,ii(1));
  head.Y3 = maxnan(bindat(5,:));
  ii = find(maxnan(bindat(3,:)) == bindat(3,:));
  head.X4 = maxnan(bindat(3,:));
  head.Y4 = bindat(5,ii(1));
    
  head.Z_min = minnan(minz);
  head.Z_max = maxnan(maxz);
    
  ex(1:length(x)) = x;
  ey(1:length(y)) = y;
  ez(1:length(z)) = z;
  
  data(1:3,1:nelem,nbin+1) = [ex;ey;ez];
  clear x y z ex ey ez 
  % save to CDF

  options.author = 'Felix Morsdorf,RSL';
  options.filename = outfilename;
  options.variablename = 'x,y,z,index';
  options.description = 'Raw Laser Data';
  options.source = 'lasercdf2';
  
  % first matrices (data, x , y, z)
  
  options.format = 'float';
  options.d1 = 'nelem';
  options.d2 = 'nbins';
  options.data1 = [1:nelem];
  options.data2 = bindat(1,:);
    
  options_x = options;       
  options_y = options;       
  options_z = options;
  
  x = squeeze(data(1,:,:));
  options_x.add_offset = minnan(x(:));
  options_x.variablename = 'x';
  x = x - options_x.add_offset;
  x(isnan(x)) = -9999;
  
  y = squeeze(data(2,:,:));
  options_y.add_offset = minnan(y(:));
  options_y.variablename = 'y';
  y = y - options_y.add_offset;
  y(isnan(y)) = -9999;
  
  z = squeeze(data(3,:,:));       
  options_z.add_offset = minnan(z(:));
  options_z.variablename = 'z';
  z = z - options_z.add_offset;
  z(isnan(z)) = -9999;
  
  % second matrix (index)
  
  options_ind=options;
  options_ind.variablename = 'index';
  options_ind.d1 = 'bounds';
  options_ind.d2 = 'nbins';
  options_ind.data1 = [1:length(bindat(:,1))];
  options_ind.data2 = [1:nbin+1];
  options_ind.add_offset = 0;       
  keyboard
  save2cdf({single(x),single(y),single(z),bindat}, ...
           {options_x,options_y,options_z,options_ind},head);    
  end

 


