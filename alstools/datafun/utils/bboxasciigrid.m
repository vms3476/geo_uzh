function [bbox,fnames] = bboxasciigrid(fdir,are);
% function [bbox,fnames] = bboxasciigrid(fdir,are);
%   function to create list of bounding box(es) of asciigrids
%   works for a list of files as created with 'dir'

% Felix Morsdorf, RSL, May 2012

filez = dir([fdir,'*.asc']);
if ~exist([fdir,'/bboxs.mat'],'file');
    for i = 1:length(filez)
        fid = fopen([fdir,'/',filez(i).name],'r');
        for j = 1:6
            str{j} = fgetl(fid);
        end
        fclose(fid);
        
        strs = {'NCOLS','NROWS','XLLCORNER','YLLCORNER', ...
                'XLLCENTER','YLLCENTER','CELLSIZE','NODATA_VALUE'};
        
        % parse header information
        for k = 1:length(strs)
            for l = 1:length(str)
                if strfind(upper(str{l}),strs{k}) > 0
                    ii = strfind(str{l},' ');
                    eval([lower(strs{k}),' = str2num(str{l}(ii+1:end));']);
                end
            end
        end
        % construct index vectors
        if exist('xllcorner')
            x = xllcorner + [(cellsize/2):cellsize:(ncols)*cellsize];
            y = yllcorner + [(cellsize/2):cellsize:(nrows)*cellsize];y = y(end:-1:1);
        else
            x = xllcenter + [0:cellsize:(ncols-1)*cellsize];
            y = yllcenter + [0:cellsize:(nrows-1)*cellsize];y = y(end:-1:1);
        end
        for k = 1:length(strs)
            eval(['clear ',lower(strs{k})]);
        end
        bbox{i} = [min(x) max(x) min(y) max(y)];
        fnames{i} = filez(i).name;
    end
    eval(['save ',fdir,'/bboxs.mat bbox fnames']);
else
    load([fdir,'/bboxs.mat'])
end
if nargin == 2 & length(are) == 4
    k = 1;
    ii = [];
    rectare = [are(1) are(3) are(2)-are(1) are(4)-are(3)];
    for i = 1:length(fnames)
        rectasc = [bbox{i}(1) bbox{i}(3) bbox{i}(2)-bbox{i}(1) bbox{i}(4)-bbox{i}(3)];
        if rectint(rectare,rectasc) > 0
        %if are(1) >= bbox{i}(1)-1 & are(2) <= bbox{i}(2)+1 & are(3) >= bbox{i}(3)-1 & ...
        %             are(4) <= bbox{i}(4)+1
            ii(k) = i;
            k = k + 1;
        end
    end
    if length(ii) >= 1
        fnames = fnames(ii);
        bbox = bbox(ii);
    else
        fnames = [];
        bbox = [];
    end
else
    for i = 1:length(bbox)
      cmap = hsv(length(bbox));
      hr = rectangle('position',[bbox{i}(1),bbox{i}(3),bbox{i}(2)-bbox{i}(1), ...
                          bbox{i}(4)-bbox{i}(3)]);
      set(hr,'edgecolor',cmap(i,:),'linewidth',2);
      hold on
      text(bbox{i}(1),bbox{i}(3),fnames{i},'interpreter','none', ...
           'fontsize',10,'fontweight','bold','verticalalign','bottom', ...
           'horizontalalign','left');
      swisstick
  end
end

                          

