function [las] = refine_power_lines(las, res, numPeaks, fillGap, minLength)

% this method uses removes small regions classified as power lines 
% which do not fall on large linear features across the tile. 
% the power line points are isolated and used to generate a raster.
% The hough, houghpeaks, and houghlines methods are used to determine
% significant lines in this raster. the lines are exteneded across the
% image extent. Image dilation and closing are used to create a solid
% region/boundary to encompass them. LAS points classified as power line
% which fall outside these linear regions are reclassified as vegetation.
%
% inputs: 
%
%   las - input LAS point cloud tile to be processed. Should be read using
%         LASread for the proper structure of fields (.record, .header)
%
%   res - [integer] resolution of the raster used as unput to hough
%
%   numPeaks - [integer] the number of peaks to identify in Hough 
%              transform. input parameter to houghpeaks
%
%   fillGap - [integer] distance between two line segments asspciated with 
%             the same Hough transform bin. input parameter to houghlines 
%
%   minLength - Minimum line length, specified as a positive real scalar. 
%               houghlines discards lines that are shorter than this value.
%               input parameter to houghlines 
%
%   
% output:
%
%   las - point cloud with selected points reclassified as 3 (vegetation) 
%         instead of power line (14) 
%
%   sample function call: las = refine_power_lines(las,2,10,10,60)
%       

% create index field in las struct, subset for "power line" returns only
las.record.index = linspace(1,numel(las.record.x),numel(las.record.x))';
las_pl = subsetraw(las.record,las.record.classification==14);

% create binary raster of highest value per cell
ras = raw2ras([las_pl.x,las_pl.y,las_pl.z],res,res,'dsm');
bw = ras.z;
bw(~isnan(bw)) = 1;
bw(isnan(bw)) = 0;

% keep track of the indices of each LAS return per raster cell
ras.idx = raw2ras_index([las_pl.x,las_pl.y,las_pl.z],res,res,las_pl.index);

% hough transform
[H,theta,rho] = hough(bw);

% identify peaks in transform
P = houghpeaks(H,numPeaks,'threshold',ceil(0.1*max(H(:))));

% find line segments corresponding to peaks in hough transform
lines = houghlines(bw,theta,rho,P,'FillGap',fillGap,'MinLength',minLength);

[rows, columns] = size(bw);

m = zeros(rows,columns);

for k = 1:length(lines)
    p1 = lines(k).point1;
    p2 = lines(k).point2;
    xy = [p1; p2]; 
    
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
%     plot([xLeft, xRight], [yLeft, yRight], 'LineWidth',2,'Color','green');
    
    
    % Plot original points on the lines
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    
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

% %plot the mask encompassing the linear regions 
%figure; imagesc(m); colormap gray

se = strel('disk',5);
dilateM = imdilate(m,se);

se = strel('disk',10);
morphM = imclose(dilateM,se);
% imshow(morphM)

% % apply the mask to the input image to keep only power line areas
keepAreas = morphM > 0;
% filteredbw = bw;
% filteredbw(~keepAreas) = 0;
% figure; imshow(filteredbw)
% title('Power line output after filtering')

filteredbw = bw;
filteredbw(keepAreas) = 0;
% figure; imshow(filteredbw)
% title('Points within these cells must change classification')

removeAreas = filteredbw==1;


% identify which LAS points need to have their classification values
% changed. loop through raster cells
for i = 1:rows
    for j = 1:columns
        
        % if cell contains points that were incorrectly classified as PL,
        if removeAreas(i,j) == 1
            
            % change the classification back to 3
            las.record.classification(ras.idx{i,j}) = 3;           
        end
    end
end

% remove index field from LAS structure 
las.record = rmfield(las.record,'index');

end


