function [raw] = getrawlas(dire,are,old);
% function [raw] = getrawlas(dire,are,switch);
% dir - directory containing LAS files
% are - perimeter in swiss national coordinates

  CD = pwd;
  cd(dire)

  if nargin == 1
      disp(['! /Users/scholl/LAStools/bin/lasmerge -i *.las -o data.las']);
      eval(['! /Users/scholl/LAStools/bin/lasmerge -i *.las -o data.las']);
%       disp(['! lasmerge -i *.laz -o data.las']);
%       eval(['! lasmerge -i *.laz -o data.las']);
    raw = readlas('data.las');
    return
  elseif nargin == 2
    disp(['! /Users/scholl/LAStools/bin/lasmerge -i *.las -o data.las -inside ',num2str(floor(are(1))),' ', ...
%     disp(['! lasmerge -i *.las -o data.laz -inside ',num2str(floor(are(1))),' ', ...
          num2str(floor(are(3))),' ',num2str(ceil(are(2))),' ', ...
          num2str(ceil(are(4)))]);
    eval(['! /Users/scholl/LAStools/bin/lasmerge -i *.las -o data.las -inside ',num2str(floor(are(1))),' ', ...
%     eval(['! lasmerge -i *.laz -o data.las -inside ',num2str(floor(are(1))),' ', ...
          num2str(floor(are(3))),' ',num2str(ceil(are(2))),' ', ...
          num2str(ceil(are(4)))]);
    raw = readlas('data.las');
    return
  else
    % old version (without lastools)
    filez = dir('*.las');
    for i = 1:length(filez)
      disp(filez(i).name);
      [raw,hdr] = readlas(filez(i).name);
      ii = (raw.x >= are(1) & raw.x <= are(2) &  raw.y >= are(3) & raw.y <= are(4));
      if i == 1
        names = fieldnames(raw);
        for j = 1:length(names)
          eval(['nraw.',names{j},' = [];'])
        end
      end
      if sum(ii) > 0
        raw = subsetraw(raw,ii);
        names = fieldnames(raw);
        for j = 1:length(names)
          eval(['nraw.',names{j},' = [nraw.',names{j},';raw.',names{j},'];'])
        end
      end
    end
  end

  cd(CD);