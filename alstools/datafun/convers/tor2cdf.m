function [] = tor2cdf(path,ofile);
% function [] = tor2cdf(path,ofile);  
% converts toposys model tol,tor,aux files 
% to one cdf-file
  
  if nargin == 0
    path = '.';
  end
  
  auxfil1 = dir([path,'/*.aux']);
  tolfil1 = dir([path,'/*.tol']);
  torfil1 = dir([path,'/*.tor']);
  for i = 1:length(auxfil1)
    auxfil{i} = auxfil1(i).name;
  end
  for i = 1:length(tolfil1)
    tolfil{i} = tolfil1(i).name;
  end
  for i = 1:length(torfil1)
    torfil{i} = torfil1(i).name;
  end
  if length(tolfil) ~= length(torfil)
    error('Number of header and data files does not match !')
  end
  
  torfiles = char(torfil)';
  torfiles = torfiles(:)';
  for i = 1:length(tolfil)
    head{i} = readtol(fullfile(path,tolfil{i}));
    filename = tolfil{i}(1:strfind(tolfil{i},'.')-1);
    ii = ceil(strfind(torfiles,filename)/length(torfil{i}));
    fid = fopen(fullfile(path,torfil{ii}),'r','l');
    data{i} = fread(fid,'uint16');
    len = head{i}.BANDS;% = 1;
    if len > 1
      data{i} = reshape(data{i},head{i}.COLUMNS,len,head{i}.LINES);
    else
      data{i} = reshape(data{i},head{i}.COLUMNS,head{i}.LINES);
    end   
    data{i}(data{i} == 0) = NaN; 
    if isfield(head{i},'HEIGHT_OFF')
      data{i} = (data{i}/100) - head{i}.HEIGHT_OFF;
    end
    fclose(fid);
    pos(i,:) = [head{i}.UP_LEFT_R,head{i}.UP_LEFT_H, ...
		head{i}.LO_RIGHT_R,head{i}.LO_RIGHT_H];
    if pos(i,1) > 2000000 % check for LV95 coords
      pos(i,[1 3]) = pos(i,[1 3]) - 2000000;
      pos(i,[2 4]) = pos(i,[2 4]) - 1000000;
    end
  end
  if length(tolfil) > 1
    npos = complex(pos(:,1),pos(:,2));
    [ix1,pox1] = same(pos(:,1));
    [iy1,poy1] = same(pos(:,2));  
    [ix2,pox2] = same(pos(:,3));
    [iy2,poy2] = same(pos(:,4));    
    [npos,idx] = sort(npos);
    bbox = [min(pos(:,1)) max(pos(:,3)) min(pos(:,4)) max(pos(:,2))];
    x = bbox(1):head{1}.HOR_SPC(1):bbox(2)-head{1}.HOR_SPC(1);
    y = bbox(4)-head{1}.HOR_SPC(2):-head{1}.HOR_SPC(2):bbox(3);

    mdat = ones(length(x),length(y))*NaN;
    
    for i = 1:length(data)
      ii = find(x >= pos(i,1) &  x < pos(i,3));
      jj = find(y >= pos(i,4) &  y < pos(i,2));
      mdat(ii,jj) = data{i};
    end
    [ii,jj] = find(~isnan(mdat));
    x1 = min(ii);x2=max(ii);
    y1 = min(jj);y2=max(jj);
    z = mdat(x1:x2,y1:y2);
    x = x(x1:x2);
    y = y(y1:y2);
  else
    bbox = [min(pos(:,1)) max(pos(:,3)) min(pos(:,4)) max(pos(:,2))];
    x = bbox(1):head{1}.HOR_SPC(1):bbox(2)-head{1}.HOR_SPC(1);
    y = bbox(4)-head{1}.HOR_SPC(2):-head{1}.HOR_SPC(2):bbox(3);
    z = data{1};
  end


  if len == 1
    if max(z(:)) > 8848 
      z = z/100;
    end
    % set metadata
    options.author = 'Felix Morsdorf RSL';
    options.filename = fullfile(path,ofile);
    options.variablename = 'z';
    options.description = head{1}.INFO;
    options.format = 'float'
    options.d1 = 'xdist';
    options.d2 = 'ydist';
    options.data1 = x;
    options.data2 = y;        
    % save to cdf 
    save2cdf(z,options);
  elseif len == 4
    rgb2cdf(x,y,data{1},head{1},ofile);
  end
  


