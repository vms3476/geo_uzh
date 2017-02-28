% read in LAS tiles with power line points classified 
% (after identify_power_lines.m )

aoi1 = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi1/n668000_255000_pl.las');
aoi1.record.index = linspace(1,numel(aoi1.record.x),numel(aoi1.record.x))';
aoi1_pl = subsetraw(aoi1.record,aoi1.record.classification==14);

aoi2 = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi2/n643000_265000_pl.las');
aoi2.record.index = linspace(1,numel(aoi2.record.x),numel(aoi2.record.x))';
aoi2_pl = subsetraw(aoi2.record,aoi2.record.classification==14);

aoi3 = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi3/n646000_265000_pl.las');
aoi3.record.index = linspace(1,numel(aoi3.record.x),numel(aoi3.record.x))';
aoi3_pl = subsetraw(aoi3.record,aoi3.record.classification==14);

aoi4 = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi4/n670000_239000_pl.las');
aoi4.record.index = linspace(1,numel(aoi4.record.x),numel(aoi4.record.x))';
aoi4_pl = subsetraw(aoi4.record,aoi4.record.classification==14);

aoi5 = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi5/n643000_240000_pl.las');
aoi5.record.index = linspace(1,numel(aoi5.record.x),numel(aoi5.record.x))';
aoi5_pl = subsetraw(aoi5.record,aoi5.record.classification==14);


%% create raster with specific resolution

res = 2;

% create raster of highest value per cell
ras1 = raw2ras([aoi1_pl.x,aoi1_pl.y,aoi1_pl.z],res,res,'dsm');
ras2 = raw2ras([aoi2_pl.x,aoi2_pl.y,aoi2_pl.z],res,res,'dsm');
ras3 = raw2ras([aoi3_pl.x,aoi3_pl.y,aoi3_pl.z],res,res,'dsm');
ras4 = raw2ras([aoi4_pl.x,aoi4_pl.y,aoi4_pl.z],res,res,'dsm');
ras5 = raw2ras([aoi5_pl.x,aoi5_pl.y,aoi5_pl.z],res,res,'dsm');

% keep track of the indices of each LAS return per raster cell
ras1_idx = raw2ras_index([aoi1_pl.x,aoi1_pl.y,aoi1_pl.z],res,res,aoi1_pl.index);
ras2_idx = raw2ras_index([aoi2_pl.x,aoi2_pl.y,aoi2_pl.z],res,res,aoi2_pl.index);
ras3_idx = raw2ras_index([aoi3_pl.x,aoi3_pl.y,aoi3_pl.z],res,res,aoi3_pl.index);
ras4_idx = raw2ras_index([aoi4_pl.x,aoi4_pl.y,aoi4_pl.z],res,res,aoi4_pl.index);
ras5_idx = raw2ras_index([aoi5_pl.x,aoi5_pl.y,aoi5_pl.z],res,res,aoi5_pl.index);

% plot
figure;myimage(ras1.x,ras1.y,ras1.z); title('power line raster AOI 1')
figure;myimage(ras2.x,ras2.y,ras2.z); title('power line raster AOI 2')
figure;myimage(ras3.x,ras3.y,ras3.z); title('power line raster AOI 3')
figure;myimage(ras4.x,ras4.y,ras4.z); title('power line raster AOI 4')
figure;myimage(ras5.x,ras5.y,ras5.z); title('power line raster APO 5')

%% make binary 
bw1 = ras1.z;bw1(~isnan(bw1)) = 1; bw1(isnan(bw1)) = 0;
figure;myimage(ras1.x,ras1.y,bw1); title('aoi 1'); colormap gray

bw2 = ras2.z;bw2(~isnan(bw2)) = 1; bw2(isnan(bw2)) = 0;
figure;myimage(ras2.x,ras2.y,bw2); title('aoi 2'); colormap gray

bw3 = ras3.z;bw3(~isnan(bw3)) = 1; bw3(isnan(bw3)) = 0;
figure;myimage(ras3.x,ras3.y,bw3); title('aoi 3'); colormap gray

bw4 = ras4.z;bw4(~isnan(bw4)) = 1; bw4(isnan(bw4)) = 0;
figure;myimage(ras4.x,ras4.y,bw4); title('aoi 4'); colormap gray

bw5 = ras5.z;bw5(~isnan(bw5)) = 1; bw5(isnan(bw5)) = 0;
figure;myimage(ras5.x,ras5.y,bw5); title('aoi 5'); colormap gray



%%


% hough transform
las = aoi4;
BW = bw4;
idx = ras4_idx;

% plot original output with power lines
figure;imshow(BW); title('Original power line output')

[H,theta,rho] = hough(BW);
numPeaks = 10;
fillGap = 10;
minLength = 60;

% identify peaks in transform
P = houghpeaks(H,numPeaks,'threshold',ceil(0.1*max(H(:))));


% plot lines corresponding to peaks on power line image
lines = houghlines(BW,theta,rho,P,'FillGap',fillGap,'MinLength',minLength);
figure, imshow(BW), title(['numPeaks: ', num2str(numPeaks), ...
    ' fillGap: ', num2str(fillGap), ...
    ' minLength: ', num2str(minLength)])

hold on
max_len = 0;

[rows, columns] = size(BW);

m = zeros(rows,columns);

for k = 1:length(lines)
    p1 = lines(k).point1;
    p2 = lines(k).point2;
    xy = [p1; p2];
    %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    
    % Get the equation of the line
    x1 = xy(1,1);
    y1 = xy(1,2);
    x2 = xy(2,1);
    y2 = xy(2,2);
    slope = (y2-y1)/(x2-x1);
    xLeft = 1; % x is on the left edge
    yLeft = slope * (xLeft - x1) + y1;
    xRight = columns; % x is on the right edge.
    yRight = slope * (xRight - x1) + y1;
    plot([xLeft, xRight], [yLeft, yRight], 'LineWidth',2,'Color','green');
    
    
    % Plot original points on the lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    
    % determine which pixels are intersected by the lines
    % across a matrix with the same dimensions as the image tile
    x1 = xLeft:xRight;
    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);
    y1 = round((x1 - p1(1)) * dy / dx + p1(2));
    
    outOfRange = y1 <= 0 | y1 > rows;
    y1(outOfRange) = [];
    x1(outOfRange) = [];
    
    idx = sub2ind(size(m), y1, x1);
    m(idx) = 1; %m(idx) = m(idx) + 1;
    
    
    
end

figure; imagesc(m); colormap gray

se = strel('disk',5);
dilateM = imdilate(m,se);
figure; imshow(dilateM)


se = strel('disk',10);
morphM = imclose(dilateM,se);
figure; imshow(morphM)

% apply the mask to the input image to keep only power line areas
keepAreas = morphM > 0;
filteredBW = BW;
filteredBW(~keepAreas) = 0;
figure; imshow(filteredBW)
title('Power line output after filtering')

filteredBW = BW;
filteredBW(keepAreas) = 0;
figure; imshow(filteredBW)
title('Points within these cells must change classification')


removeAreas = filteredBW==1;


% identify which LAS points need to have their classification values
% changed. loop through raster cells
for i = 1:rows
    for j = 1:columns
        
        % if cell contains points that were incorrectly classified as PL,
        if removeAreas(i,j) == 1
            
            % change the classification back to 3
            las.record.classification(ras2_idx{i,j}) = 3;           
        end
    end
end


