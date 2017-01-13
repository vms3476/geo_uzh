% this script assesses the standard deviation to identify / filter power
% line points 

% AOI 1 with power lines
xMin = 668200; xMax = 668300; yMin = 255100; yMax = 255300;

% lower half to analyze power lines above trees 
xMin = 668200; xMax = 668300; yMin = 255100; yMax = 255180;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi1/';
las = getrawlas(lasDir,area);

% before subsetting or modifying 
raw = las; 
figure; myscatter3(raw.x,raw.y,raw.z,raw.z,parula);

% % remove points below 4m
i = las.z > 4;
las = subsetraw(las,i);
% 
% % remove buildings
% ii = ~ismember(las.Classification,6);
% las = subsetraw(las,ii);

% plot
figure; myscatter3(las.x,las.y,las.z,las.z,parula);
title('LAS subset with power lines','FontSize',14); swisstick


%% AOI 2 with power lines 
xMin = 643630; xMax = 643780; yMin = 265600; yMax = 265700;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi2/';
raw2 = getrawlas(lasDir,area);
las2 = raw2; 
figure; myscatter3(las2.x,las2.y,las2.z,las2.z,parula);

%% create raster to assess pixels with high std dev 

% standard deviation
% must include ground points for power line pixels to have high values 

% 2m resolution
ras2 = raw2ras([raw.x,raw.y,raw.z],2,1,'dsm');
figure; myimage(ras2.x,ras2.y,ras2.std); title('STD 2m'); colorbar

[ras2.X, ras2.Y] = meshgrid(ras2.x,ras2.y)
figure; myscatter3(ras2.X(:),ras2.Y(:),ras2.std(:),ras2.std(:),parula); title('STD 2m');


% 1m resolution
ras1 = raw2ras([raw.x,raw.y,raw.z],1,1,'dsm');
figure; myimage(ras1.x,ras1.y,ras1.std); 
title('STD 1m','FontSize',14); swisstick

[ras1.X, ras1.Y] = meshgrid(ras1.x,ras1.y)
figure; myscatter3(ras1.X(:),ras1.Y(:),ras1.std(:),ras1.std(:),parula); 
title('STD 1m','FontSize',14); swisstick

figure; 
subplot(1,2,1); myimage(ras2.x,ras2.y,ras2.std); title('STD 2m'); colorbar;
subplot(1,2,2); myimage(ras1.x,ras1.y,ras1.std); title('STD 1m'); colorbar;



%% Isolate pixels with high std

i = ras1.std > 9; 
ras = ras1.std; 
ras(i) = 1; 
ras(~i) = 0; 
figure; myimage(ras1.x,ras1.y,ras); title('STD 1m');


%% loop through pixels 

resV = 1; % (m) resolution of vertical histogram bins
binCenters = [0.5:resV:49.5]; % max height 50m 

% create a raster with vertical histogram for each cell
peri = 1;
ras = raw2vox([raw.x,raw.y,raw.z],peri,peri,'dsm');
[ras.X, ras.Y] = meshgrid(ras.x,ras.y);
vol = peri * peri * resV; % volume of single voxel

% for power line candidate pixels, assess histogram 
k = ras.std > 9; 
rasOrig = ras; % copy of ras
las = raw; % copy of PC 

% % remove points below 4m
% i = las.z > 4;
% las = subsetraw(las,i);

% binary value, 1 = power line candidate voxel
ras.std(k) = 1; 
ras.std(~k) = 0; 


fig1 = figure; 
fig2 = figure; 
fig3 = figure;

% plot 2D area of interest to highlight current cell
figure(fig1)
subplot(141)
myimage(ras.x,ras.y,ras.std); 

% plot vertical histogram of current cell
subplot(142) 
title('Vertical Histogram'); xlabel('Counts'), ylabel('Bin Edge (m)')

%plot the PC to visualize which points fall within the cell 
subplot(143)

% reset colormap for first subplot
subplot(141); colormap parula

% to show plots and pause during each loop iteration
plot = 0;

% loop through pixels, highlight current one and display histogram
for i = 1:numel(ras.x) 
    
    for j = 1:numel(ras.y)
        
        if ras.std(j,i) == 0 
            continue
        end
        
        
        % highlight current cell in 2D raster
        ras.std(j,i) = 2;
        
        % points within current cell
        in = raw.x > (ras.x(i) - peri/2) & raw.x < (ras.x(i) + peri/2) & ...
        raw.y > (ras.y(j) - peri/2) & raw.y < (ras.y(j) + peri/2);
        inIdx = find(in); 
        zIn = subsetraw(raw,in);        
        
        % k means clustering
        kIn = [zIn.x,zIn.y,zIn.z];
        % (K: number of clusters, G: assigned groups, C: cluster centers)
        K = 3;
        % check if input has more points than the number of clusters
        if size(kIn,1) < K
            continue
        end
        [G,C] = kmeans(kIn,K,'Distance','sqeuclidean','Start','sample','MaxIter',100,'Replicates',10);

        % Power line classification decisions
        
        % assess the distance between the two tallest classes
        % min height of class 1 points - max height of class 2 points 
        % sort culster centers from tallest to shortest
        [tmp, sortIdx] = sort(C(:,3),'descend');
        class1 = sortIdx(1);
        class2 = sortIdx(2);
        class3 = sortIdx(3);
        class1min = min(kIn(G==class1,3));
        class2max = max(kIn(G==class2,3));
        class2min = min(kIn(G==class2,3));
        class3max = max(kIn(G==class3,3));
        classDif12 = class1min - class2max;
        classDif23 = class2min - class3max;
        
        numelClass1 = numel(kIn(G==class1));
        numelClass2 = numel(kIn(G==class2));
        
        zDifClass1 = max(kIn(G==class1,3)) - min(kIn(G==class1,3));
        zDifClass2 = max(kIn(G==class2,3)) - min(kIn(G==class2,3));
        
        % get the points within neighborhood of cluster 1 center
        peri2 = 5; % dimension of square neighborhood 
        in2 = raw.x > (C(class1,1) - (peri2/2)) & raw.x < (C(class1,1) + (peri2/2)) & ...
        raw.y > (C(class1,2) - (peri2/2)) & raw.y < (C(class1,2) + (peri2/2)) & ...
        raw.z > (C(class1,3) - (peri2/2)) & raw.z < (C(class1,3) + (peri2/2));
        in2Idx = find(in2);     
        zIn2 = subsetraw(raw,in2Idx);
        % number of points within neighborhood 
        ptsInNeighborhood = numel(zIn2.x);
        if ptsInNeighborhood == 0
            disp('no points in neighborhood')
            density = 0 
        else
            % calculate density of points   
            density = ptsInNeighborhood / (peri2 * peri2 * peri2)   
        end
        
     
% % %         %  if only 1 bin above 20m has points, all are power line candidates 
% % %         if sum(ras.vHist(j,i,25:end) > 1) == 1 
% % %             % only keep those with low point densities
% % %             if pd < 10 
% % %                 las.Classification(in & las.z>25) = 14;   
% % %             end
% % %         end

        
        
        
        % for the taller 2 clusters, subract class 1 min - class 2 max 
        % if distance is greater than threshold, assign class 1 as PL 
        if classDif12 > 25  && density < 0.1 
            % find indices of las points in polygon and of class 1
            plIdx = inIdx(G==class1);
            las.Classification(plIdx) = 14;
            
            
        % check for the case when cluster 1 has few points and there is a 
        % distance of at least 5m to the cluster below 
        elseif (numel(kIn(G==class1)) <= numel(kIn(G==class2))) && classDif12 > 5 && density < 0.1 
            
            plIdx = inIdx(G==class1);
            las.Classification(plIdx) = 14;

        end
        

        
        
        
        
        
        
        
        
        
        
        % visualize in each iteration 
        if plot 
            
            disp(['Current pixel: ' num2str(j) ',' num2str(i)])
            
            % Figure 1
            % 1) show STD raster, highlight current cell
            figure(fig1)
            subplot(141)
            myimage(ras.x,ras.y,ras.std);
            title(['Current cell in STD ' num2str(peri) 'm']);
        
            % 2) vertical histogram
            subplot(142)
            barh(binCenters,reshape(ras.vHist(j,i,:),1,50)); xlim([0,20]); ylim([0,42])
            title('Vertical Histogram');
            
            % 3) points within current cell
            subplot(143)
            t = myscatter3(zIn.x,zIn.y,zIn.z,zIn.z,parula); zlim([0,max(zIn.z)]);   
            title('LAS points')
            
            % 4)k means clusters
            subplot(144)
            clr = lines(K);
            hold on
            km1 = scatter3(kIn(:,1), kIn(:,2), kIn(:,3), 36, clr(G,:), 'Marker','.');
            km2 = scatter3(C(:,1), C(:,2), C(:,3), 100, clr, 'Marker','o', 'LineWidth',3);
            hold off
            view(3), axis equal, box on, rotate3d on, zlim([0,max(zIn.z)]);
            title('k means clusters')
            
            % Figure 2. highlight power line points 
            figure(fig2);
            myscatter3(las.x,las.y,las.z,las.Classification,parula); view(2);

            % plot the points within surrounding volume
            if numel(zIn2.x) > 0 
                figure(fig3); 
                hold on
                km3 = myscatter3(zIn2.x, zIn2.y, zIn2.z, zIn2.z, parula);
                km4 = scatter3(C(class1,1), C(class1,2), C(class1,3), 100, 'Marker','o', 'LineWidth',3);
                hold off
            else
                disp('no points within neighborhood of cluster1 center')
            end
            
            pause
            
            clf % clear LAS classification figure for the next time 
        
            % clear subplot 141, 143, and 144 for the next plot
            delete(t)
            delete(km1)
            delete(km2)
            delete(km3)
            delete(km4)
        end 
          
        ras.std(j,i) = 1;
             
    end
    
end



%% kmeans clustering - LAS with power lines

X = [las.x,las.y,las.z];
[numInst,numDims] = size(X); 

% K-means clustering
% (K: number of clusters, G: assigned groups, C: cluster centers)

K = 3;
distance = 'sqeuclidean';
start = [668260.758000000 255137.712000000 39.2110000000000; ...  % powerline
         668280.927000000 255263.473000000 26.4790000000000; ... % tree
         668206.634000000 255105.571000000 0.180000000000000];    % ground

% start = [668237.016000000 255184.422000000 36.7610000000000; ...  % powerline
%          668260.231000000 255173.891000000 15.116000000000000; ... % tree
%          668206.634000000 255105.571000000 0.180000000000000];    % ground     
     
% start = [668260.758000000 255137.712000000 39.2110000000000; ...  % powerline
%          668260.758000000 255137.712000000 15.116000000000000; ... % tree
%          668260.758000000 255137.712000000 0.180000000000000];    % ground

maxIter = 1;
[G,C] = kmeans(X,K,'Distance',distance,'Start',start,'MaxIter',maxIter);

% show points and clusters 
clr = lines(K);
figure, hold on
scatter3(X(:,1), X(:,2), X(:,3), 36, clr(G,:), 'Marker','.')
scatter3(C(:,1), C(:,2), C(:,3), 100, clr, 'Marker','o', 'LineWidth',3)
hold off
view(3), axis equal, box on, rotate3d on
xlabel('x'), ylabel('y'), zlabel('z')
title(['Distance: ', distance, ', MaxIter: ',num2str(maxIter)],'FontSize',14)


%% 


