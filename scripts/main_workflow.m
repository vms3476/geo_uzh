% this script is the main workflow for the individual tree detection 

% define general file paths and values 
lastoolsBinPath = '/Users/scholl/LAStools/bin/'; % lastools path

% resolution of sub-tiles 
tile_res = 250; 


% LAS directory to be processed
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/vorwalden/';
currentDir = pwd;


% loop through LAS / LAZ files 


% current LAS filename
lasName = 'n633000_235000.las';


% determine the list of x and y coordinates to loop through
[~,cmdout] = system([lastoolsBinPath 'lasinfo -i ' lasDir lasName]);
cmdout = strsplit(cmdout,' ');
xStart = round(str2double(cmdout{(find(ismember(cmdout,'min')))+4}));
yStart = round(str2double(cmdout{(find(ismember(cmdout,'min')))+5}));
xEnd = round(str2double(cmdout{(find(ismember(cmdout,'max')))+4}));
yEnd = round(str2double(cmdout{(find(ismember(cmdout,'max')))+5}));

xList = xStart:tile_res:xEnd;
if xList(end)< xEnd
    xList = [xList xEnd];
end

yList = yStart:tile_res:yEnd;
if yList(end)< yEnd
    yList = [yList yEnd];
end


%% 

% loop through subtiles
for x = 1:numel(xList)-1
    for y = 1:numel(yList)-1
        
        % read current tile subset LAS PC data
        area = [xList(x) xList(x+1) yList(y) yList(y+1)];
        
        % write LAS data within current subtile extents to 'data.las'
        % if a file with this name already exists, delete it first. 
        cd(lasDir)
        if exist('data.las') == 2
            unix('rm data.las') 
        end
        eval(['! ' lastoolsBinPath 'lasmerge -i *.las -o data.las -inside ', ...
             num2str(floor(area(1))),' ', num2str(floor(area(3))), ...
             ' ',num2str(ceil(area(2))),' ', num2str(ceil(area(4)))]);
        raw = readlas('data.las');
        
        
        % forest type classification 
        
        % must first run the forest_type_classification.m code with the 
        % desired input areas to train the model and features.
        % load the output model and the feature info
        model = load('/Users/scholl/geo_uzh/data/KantonAargau/randomForest/MdlRF.mat');
              
        
        % variable to indicate the resolution of each forest type classification
        % per subtile. since lasmerge is used to find data INSIDE the bounds,
        % 0.5 must be added to the lower bound to define the lower pixel center,
        % and 0.5 must likewise be subtracted from the upper bound to determine the 
        % center coordinate of the upper pixel extent (to match the
        % coordinates used in the WSL "truth" image, but take away the dependence
        % on this data for this workflow)
        res.x = linspace(xList(x)+0.5,xList(x+1)-0.5,tile_res);
        res.y = fliplr(linspace(yList(y)+0.5,yList(y+1)-0.5,tile_res));
        
        
        % compute the features for the input LAS
        Xnew = computeFeaturesRF(raw,res,1,0,features);
        [ypred,unaries] = model.Mdl.predict(Xnew);
        ypred = str2num(cell2mat(ypred)); 
        ypred = reshape(ypred,[numel(res.x),numel(res.y)]); 
        
        % add ground regions (z<3m), perform morphological processing
        % (median filtering, area opening) to remove noise
        [ypred_withGround] = morphological(ypred, res, raw, 5, 25);
        
        
        
        % use corresponding area in WSL image as raster resolution
        % why aren't the features being computed properly??? 
        mapx = [xList(x) xList(x+1) xList(x+1) xList(x)]; mapy = [yList(y) yList(y+1) yList(y) yList(y+1)];
        [wsl.data,wsl.x,wsl.y,wsl.info] = geoimread(wslTif,mapx,mapy);
        wsl.data(wsl.data==3) = 0; 

        
        

    end
end

cd(currentDir)








% forest type classification 

% check majority forest type 

% apply ITC method 

% coniferous: locmax vertical clustering

% deciduous: watershed 