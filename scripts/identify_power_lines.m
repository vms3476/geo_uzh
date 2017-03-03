function [ raw_all ] = identify_power_lines( inputFilename, outputDir, configFilename )

% this method must be in the same directory which contains the csv file,
% identify_power_lines_config.csv. This file contains the following 
% input parameters for identifying power line points, changing their 
% classification values, and performing hough trans


% sample function call:
%       identify_power_lines('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi7/658000_250000_subset.las',...
%                            '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi7/pl/'...
%                            '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/identify_power_lines_config.csv')



% get the parts of the LAZ or LAS path 
[pathstr,filename,ext] = fileparts(inputFilename); 

% add filesep character onto the end out outputDir 
if outputDir(end) ~= filesep
    outputDir = [outputDir filesep];
end

% create outpu directory 
unix(['mkdir ' outputDir]);


% read the configuration text file for input parameters:
%   tile_res = integer, resolution of tile subset (100, 200, 250m)
%   res = integer, resolution of raster for power line identification.
%       the current other parameters have been tested when raster_res = 1
%   peri =  integer, perimeter of neighborhood around cell of intrest to 
%       assess during power line identification method 
%   densityThr = double, threshold over which it is too dense to classify
%       the cluster within the current cell as power lines 
fileID = fopen(configFilename,'r'); 
configInfo = textscan(fileID, '%s','delimiter',',','HeaderLines',1); 
fclose(fileID);

lastoolsBinPath = configInfo{1}{1};
tile_res = str2double(configInfo{1}{2});
res = str2double(configInfo{1}{3}); 
peri = str2double(configInfo{1}{4});
densityThr = str2double(configInfo{1}{5}); 


% if input is LAZ, convert to LAS 
laszip = [lastoolsBinPath 'laszip'];
if ext(end) == 'z'
    lasName = [pathstr filesep filename '.las'];
    inputIsLAS = 0;
    if ~exist(lasName,'file')
        unix([laszip ' -i ' inputFilename ' -o ' lasName]);
    end
% if input is already in LAS form, make note of it to not delete the
% original LAS file later
else
    lasName = inputFilename;
    inputIsLAS = 1;
end


% create temporary directory to save LAS subsets
tmpDir = [pathstr filesep filename];
unix(['mkdir ' tmpDir]);

% determine starting x and y coordinates using lasinfo
[~,cmdout] = unix([lastoolsBinPath 'lasinfo -i ' lasName]); % lasinfo
minIdx = strfind(cmdout,'min x y z:'); 
xStart = str2double(cmdout(minIdx+28:minIdx+33));
yStart = str2double(cmdout(minIdx+38:minIdx+44));
maxIdx = strfind(cmdout,'max x y z:'); 
xEnd = ceil(str2double(cmdout(maxIdx+28:maxIdx+36)));
yEnd = ceil(str2double(cmdout(maxIdx+38:maxIdx+47)));


% list of starting and end coordinates for the tile subsets
% based on the tile_res parameter and the extent of the entire LAS file
xList = xStart:tile_res:xEnd;
if xList(end)< xEnd
    xList = [xList xEnd];
end

yList = yStart:tile_res:yEnd;
if yList(end)< yEnd
    yList = [yList yEnd];
end



counter = 1;


% loop through tile subset 
for x = 1:numel(xList)-1
    for y = 1:numel(yList)-1
        
        % read current tile subset LAS PC data
        xList(x);
        yList(y); 
        area = [xList(x) xList(x+1) yList(y) yList(y+1)];
        
        % if output from getrawlas already in mainDir, delete it 
        if exist([pathstr filesep 'data.las'],'file') == 2
            unix(['rm ' pathstr filesep 'data.las']);
        end
        
        % read the subset of LAS data, use mParkan method LASread
        getrawlas(pathstr,area,lastoolsBinPath);
        raw_all = LASread([pathstr filesep 'data.las'],false,false);
        
        
       
        %%%%%%%%%%%%%%%% identify power lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create indices for all PC entries in the input LAS 
        raw_all.record.index = linspace(1,numel(raw_all.record.x),numel(raw_all.record.x))';

        % first iteration, use all input LAS points
        raw = raw_all.record; 
        reiterate = 1; 
        nClusters = 3; % number of clusters for k means
        
        % iterate until there are no more power line points found
        while reiterate
            
            % generate DSM raster; select std thresh based on histogram peak
            ras = raw2ras([raw.x,raw.y,raw.z],res,res,'dsm');
            if sum(ras.std(:)>1) == 0
                break
            end
            [f,xi] = ksdensity(ras.std(ras.std>1));
            [peaks, heights] = findpeaks(f,xi);
            stdThr = heights(find(ismember(peaks,max(peaks))))-2;

            % in regions where a single power line is over the water,
            % the std is very low (no ground data) but the height from
            % the dsm raster is high, so increase the std 
            ras.std((ras.std<1)&(ras.z>10))=stdThr; 

            % apply standard deviation threshold to identify candidate pixels
            k = ras.std >= stdThr; 
            ras.std(k) = 1; 
            ras.std(~k) = 0; 

            for i = 1:numel(ras.x) 
                for j = 1:numel(ras.y)

                    % only assess candiate pixels
                    if ras.std(j,i) == 0 
                        continue
                    end

                    % ALS points within current cell
                    in = raw.x > (ras.x(i) - peri/2) & raw.x < (ras.x(i) + peri/2) & ...
                    raw.y > (ras.y(j) - peri/2) & raw.y < (ras.y(j) + peri/2) & raw.classification ~=6;
                    zIn = subsetraw(raw,in);

                    % k means clustering (3 clusters, G: assigned groups, C: cluster centers)
                    kIn = [zIn.x,zIn.y,zIn.z];
                    % check if input has more points than the number of clusters
                    if size(kIn,1) < nClusters
                        continue
                    end

                    % starting points for clusters
                    start = [ras.x(i) ras.y(j) 50; ...  % powerline
                             ras.x(i) ras.y(j) 25; ...  % tree
                             ras.x(i) ras.y(j) 0];      % ground
                   [G,C] = kmeans(kIn,nClusters,'Distance','sqeuclidean','Start',start,'MaxIter',100);


                    % Power line clusterification decisions

                    % assess the distance between the two tallest clusteres
                    % min height of cluster 1 points - max height of cluster 2 points 
                    % sort culster centers from tallest to shortest
                    [~, sortIdx] = sort(C(:,3),'descend');
                    cluster1 = sortIdx(1);
                    cluster2 = sortIdx(2);
                    cluster3 = sortIdx(3);
                    cluster1min = min(kIn(G==cluster1,3));
                    cluster2max = max(kIn(G==cluster2,3));
                    cluster2min = min(kIn(G==cluster2,3));
                    cluster3max = max(kIn(G==cluster3,3));
                    clusterDif12 = cluster1min - cluster2max;
                    clusterDif23 = cluster2min - cluster3max;

                    % get the points within neighborhood of cluster 1 center
                    in2 = raw.x > (C(cluster1,1) - (peri/2)) & raw.x < (C(cluster1,1) + (peri/2)) & ...
                    raw.y > (C(cluster1,2) - (peri/2)) & raw.y < (C(cluster1,2) + (peri/2)) & ...
                    raw.z > (C(cluster1,3) - (peri/2)) & raw.z < (C(cluster1,3) + (peri/2));
                    in2Idx = find(in2);     
                    zIn2 = subsetraw(raw,in2Idx);

                    % number of points within neighborhood 
                    ptsInNeighborhood = numel(zIn2.x);
                    if ptsInNeighborhood == 0
                        density = 0;
                    else
                        % calculate density of points   
                        density = ptsInNeighborhood / (peri * peri * peri) ;
                    end

                    % points along the same z plane (assuming power lines are linear
                    % and relatively flat) with a wider x,y extent to assign as power
                    % lines. hopefully will get the vertical poles too. 
                    % do not include points classified as buildings (6) or
                    % ground (2) 
                    in2 = raw.x > (C(cluster1,1) - (peri)) & raw.x < (C(cluster1,1) + (peri)) & ...
                    raw.y > (C(cluster1,2) - (peri)) & raw.y < (C(cluster1,2) + (peri)) & ...
                    raw.z > (C(cluster1,3) - (peri/2)) & raw.z < (C(cluster1,3) + (peri/2)) & ...
                    raw.classification ~=6 & raw.classification ~=2;
                    in2Idx = find(in2);     

                    % for the taller 2 clusters, subract cluster 1 min - cluster 2 max 
                    % if distance is greater than threshold, and if density is below
                    % the density threshold, assign cluster 1 as PL 
                    if clusterDif12 > 25  && density < densityThr 
                        % find indices of las points in polygon and of cluster 1
                        raw.classification(in2Idx) = 14;          

                    % check for the case when cluster 1 has few points and there is a 
                    % distance of at least 5m to the cluster below 
                    elseif (numel(kIn(G==cluster1)) <= numel(kIn(G==cluster2))) && clusterDif12 > 5 && density < densityThr 
                        raw.classification(in2Idx) = 14;

                    % when clusters 1 and 2 are both on power line areas
                    elseif clusterDif12 < 3 && clusterDif23 > 10 && (numel(kIn(G==cluster1))+numel(kIn(G==cluster2))<= numel(kIn(G==cluster3))) && density < densityThr
                        raw.classification(in2Idx) = 14;

                    % For the case when powerlines are above water
                    % if the minimum height for all 3 clusters greater than 10m 
                    elseif min(kIn(:,3)) > 10 
                        raw.classification(in2Idx) = 14;
                    end
                end
            end

            % get indices of the points just identified as PL 
            pl = raw.classification==14;

            % if power lines were identified during this iteration, loop again
            if sum(pl) > 0 
                % change classification value in the raw PC of all points
                pl_index = raw.index(pl);
                raw_all.record.classification(ismember(raw_all.record.index,pl_index))=14;
                % remove the power line points and reiterate
                raw = subsetraw(raw,raw.classification~=14);
                reiterate = 1;
            else % no more power lines to find, exit while loop 
                reiterate = 0;        
            end
        end
        
        % remove index field from LAS structure 
        raw_all.record = rmfield(raw_all.record,'index');
        
        % save LAS subset with power line classification in tmp dir
        LASwrite(raw_all,[tmpDir filesep filename '_' num2str(counter) '.las'],'version', 12, 'verbose', false);

        % variable to number each output file
        counter = counter + 1;
        
    end
    
end

% combine all tile subsets into single las file
unix([lastoolsBinPath 'lasmerge -i ' tmpDir filesep '*.las -o ' tmpDir '_pl.las'])

% read LAS tile to assess large linear features. write output to LAS
las = LASread([tmpDir '_pl.las'],false,false);
las = refine_power_lines(las,2,10,10,60);
LASwrite(las,[tmpDir '_pl.las'],'version', 12, 'verbose', false);

% compress las to laz file 
unix([lastoolsBinPath 'laszip -i ' tmpDir '_pl.las -o ' outputDir filename '_pl.laz']);

% % delete files that are no longer necessary 
unix(['rm -r ' tmpDir]); 

% If the input was LAS, do not delete it 
if ~inputIsLAS
    unix(['rm ' lasName]);
end
unix(['rm ' pathstr filesep 'data.las']);
unix(['rm ' tmpDir '_pl.las']);



end 
