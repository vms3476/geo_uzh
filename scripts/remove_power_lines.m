function [ raw_pl, las ] = remove_power_lines( raw, res, stdThr, peri, nClusters, densityThr)

% description:
%       Identifies and removes ALS echos belonging to power lines 
%
% inputs: 
%   
%       raw - lidar struct containing x,y,z 
%
%       res - resolution of raster cell, typically 1 (m)
%
%       stdThr - threshold, above which pixels are considered candiates to 
%               contain power lines, applied to raster with resolution "res"
%               typically 7 during testing
%
%       peri - dimension of cubic neighborhood of current cluster center;
%              typically 6
%              
%       nClusters - number of clusters as input to K-means clustering;
%                   typically 3
%
%       densityThr = density threshold, below which the cluster is 
%                    clusterified as power lines; typically 0.1
%
%
% outputs:
%
%       raw_pl - raw PC with power line points assigned to cluster 14
% 
%       las - PC with power line points removed 
%
% by Victoria Scholl, 18 Jan 2017


% % Loop until there are no PL classified points in current iteration


% create raster from raw PC data
ras = raw2ras([raw.x,raw.y,raw.z],res,res,'dsm');

% % select std thresh based on histogram peak
[f,xi] = ksdensity(ras.std(ras.std>1));
[peaks, heights] = findpeaks(f,xi);
stdThr = heights(find(ismember(peaks,max(peaks))));

% in regions where a single power line is over the water,
% the std is very low (no ground data) but the height from
% the dsm raster is high, so increase the t
ras.std((ras.std<1)&(ras.z>10))=stdThr; 

% apply standard deviation threshold to identify candidate pixels
k = ras.std >= stdThr; 
ras.std(k) = 1; 
ras.std(~k) = 0; 

las = raw;

for i = 1:numel(ras.x) 
    for j = 1:numel(ras.y)
        
        % only assess candiate pixels
        if ras.std(j,i) == 0 
            continue
        end
        
        % ALS points within current cell
%        in = raw.x > (ras.x(i) - res/2) & raw.x < (ras.x(i) + res/2) & ...
%        raw.y > (ras.y(j) - res/2) & raw.y < (ras.y(j) + res/2);
        in = raw.x > (ras.x(i) - peri/2) & raw.x < (ras.x(i) + peri/2) & ...
        raw.y > (ras.y(j) - peri/2) & raw.y < (ras.y(j) + peri/2);
        inIdx = find(in); 
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
        [tmp, sortIdx] = sort(C(:,3),'descend');
        cluster1 = sortIdx(1);
        cluster2 = sortIdx(2);
        cluster3 = sortIdx(3);
        cluster1min = min(kIn(G==cluster1,3));
        cluster2max = max(kIn(G==cluster2,3));
        cluster2min = min(kIn(G==cluster2,3));
        cluster3max = max(kIn(G==cluster3,3));
        clusterDif12 = cluster1min - cluster2max;
        clusterDif23 = cluster2min - cluster3max;
        
        numelCluster1 = numel(kIn(G==cluster1));
        numelCluster2 = numel(kIn(G==cluster2));
        
        zDifCluster1 = max(kIn(G==cluster1,3)) - min(kIn(G==cluster1,3));
        zDifCluster2 = max(kIn(G==cluster2,3)) - min(kIn(G==cluster2,3));
        
        % get the points within neighborhood of cluster 1 center
        in2 = raw.x > (C(cluster1,1) - (peri/2)) & raw.x < (C(cluster1,1) + (peri/2)) & ...
        raw.y > (C(cluster1,2) - (peri/2)) & raw.y < (C(cluster1,2) + (peri/2)) & ...
        raw.z > (C(cluster1,3) - (peri/2)) & raw.z < (C(cluster1,3) + (peri/2));
        in2Idx = find(in2);     
        zIn2 = subsetraw(raw,in2Idx);
        
        % number of points within neighborhood 
        ptsInNeighborhood = numel(zIn2.x);
        if ptsInNeighborhood == 0
            %disp('no points in neighborhood')
            density = 0;
            
        else
            % calculate density of points   
            density = ptsInNeighborhood / (peri * peri * peri) ;

        end
       
        % points along the same z plane (assuming power lines are linear
        % and relatively flat) with a wider x,y extent to assign as power
        % lines. hopefully will get the vertical poles too. 
        in2 = raw.x > (C(cluster1,1) - (peri)) & raw.x < (C(cluster1,1) + (peri)) & ...
        raw.y > (C(cluster1,2) - (peri)) & raw.y < (C(cluster1,2) + (peri)) & ...
        raw.z > (C(cluster1,3) - (peri/2)) & raw.z < (C(cluster1,3) + (peri/2));
        in2Idx = find(in2);     
        zIn2 = subsetraw(raw,in2Idx);
        
        
        % for the taller 2 clusters, subract cluster 1 min - cluster 2 max 
        % if distance is greater than threshold, and if density is below
        % the density threshold, assign cluster 1 as PL 
        if clusterDif12 > 25  && density < densityThr 
            % find indices of las points in polygon and of cluster 1
            plIdx = inIdx(G==cluster1);
            %las.Classification(plIdx) = 14;
            las.Classification(in2Idx) = 14;
            
            
        % check for the case when cluster 1 has few points and there is a 
        % distance of at least 5m to the cluster below 
        elseif (numel(kIn(G==cluster1)) <= numel(kIn(G==cluster2))) && clusterDif12 > 5 && density < densityThr 
            
            plIdx = inIdx(G==cluster1);
            %las.Classification(plIdx) = 14;
            las.Classification(in2Idx) = 14;
      
            
        %end 
        
        % TESTING this decision for AOI 4 
        elseif clusterDif12 < 3 && clusterDif23 > 10 && (numel(kIn(G==cluster1))+numel(kIn(G==cluster2))<= numel(kIn(G==cluster3))) && density < densityThr
        las.Classification(in2Idx) = 14;
        
        
        % For the case when powerlines are above water
        % if the minimum height for all 3 clusters greater than 10m 
        elseif min(kIn(:,3)) > 9 
            las.Classification(in2Idx) = 14;
        
        end
        
        
    end
    
end

raw_pl = las; % return PC with PL points classified as 14
las = subsetraw(las,las.Classification~=14); % remove PL points from PC

    

end

