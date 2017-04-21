% this script is the main workflow for the individual tree detection 

% define lastools path 
lastoolsBinPath = '/Users/scholl/LAStools/bin/';



% LAS directory to be processed
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/vorwalden/';

lasName = 'n633000_235000.las';

[~,cmdout] = system([lastoolsBinPath 'lasinfo -i ' lasDir lasName])
cmdout = strsplit(cmdout,' ');
cmdout{(find(ismember(cmdout,'min')))+1} == 'x'
xStart = round(str2double(cmdout{(find(ismember(cmdout,'min')))+4}));
yStart = round(str2double(cmdout{(find(ismember(cmdout,'min')))+5}));
xEnd = round(str2double(cmdout{(find(ismember(cmdout,'min')))+11}));
yEnd = round(str2double(cmdout{(find(ismember(cmdout,'min')))+12}));

xList = xStart:tile_res:xEnd;
if xList(end)< xEnd
    xList = [xList xEnd];
end

yList = yStart:tile_res:yEnd;
if yList(end)< yEnd
    yList = [yList yEnd];
end


% read input LAS tile (entire tile, or just for each subset?) 


eval(['! ' lastoolsBinPath 'lasmerge -i *.las -o data.las -inside ',num2str(floor(are(1))),' ', ...
          num2str(floor(are(3))),' ',num2str(ceil(are(2))),' ', ...
          num2str(ceil(are(4)))]);



% forest type classification 

% check majority forest type 

% apply ITC method 

% coniferous: locmax vertical clustering

% deciduous: watershed 