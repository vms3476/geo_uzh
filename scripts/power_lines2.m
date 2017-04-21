function [ ] = power_lines2( las, ras )

% compute features for input LAS data 


% calculate power line echo: first return / total returns 
powerline_echo = ras.r1 ./ ras.den; 
% figure; myimage(ras.x,ras.y,powerline_echo); title('power line echo: r1 / rTotal')

% threshold power line echo raster
pl_echo_thr = powerline_echo; pl_thresh = 0.8;
pl_echo_thr(pl_echo_thr<pl_thresh) = 0; pl_echo_thr(isnan(pl_echo_thr)) = 0;
pl_echo_thr(pl_echo_thr>=pl_thresh) = 1; 
figure; myimage(ras.x,ras.y,pl_echo_thr); title('power line echo thresholded')

dilateM=bwareaopen(pl_echo_thr,16);
figure;myimage(ras.x,ras.y,dilateM)


% 
% % binary raster based on DSM 
% bw = ras.z;
% bw(~isnan(bw)) = 1;
% bw(isnan(bw)) = 0;
% figure; myimage(ras.x,ras.y,bw)
% 
% % point density 
% rasDen = raw2ras([c.x,c.y,c.z],2,2,'den');
% figure;myimage(ras.x,ras.y,rasDen.z); colorbar
% 
% rasDenThr = rasDen.z; 
% rasDenThr(rasDenThr > 100) = 0;
% rasDenThr(isnan(rasDenThr)) = 0;
% rasDenThr(rasDenThr>0)=1;
% figure;myimage(ras.x,ras.y,rasDenThr); 
% 
% %% Hough
% 
% bw = pl_echo_thr;
% 
% [H,theta,rho] = hough(bw);
% 
% % identify peaks in transform
% numPeaks = 10;
% fillGap = 10;
% minLength = 60;
% P = houghpeaks(H,numPeaks,'threshold',ceil(0.1*max(H(:))));
% 
% % find line segments corresponding to peaks in hough transform
% lines = houghlines(bw,theta,rho,P,'FillGap',fillGap,'MinLength',minLength);
% 
% figure, imshow(bw), title(['numPeaks: ', num2str(numPeaks), ...
%     ' fillGap: ', num2str(fillGap), ...
%     ' minLength: ', num2str(minLength)])
% 
% hold on
% max_len = 0;
% 
% [rows, columns] = size(bw);
% 
% m = zeros(rows,columns);
% 
% for k = 1:length(lines)
%     p1 = lines(k).point1;
%     p2 = lines(k).point2;
%     xy = [p1; p2];
%     %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     
%     % Get the equation of the line
%     x1 = xy(1,1);
%     y1 = xy(1,2);
%     x2 = xy(2,1);
%     y2 = xy(2,2);
%     slope = (y2-y1)/(x2-x1);
%     xLeft = 1; % x is on the left edge
%     yLeft = slope * (xLeft - x1) + y1;
%     xRight = columns; % x is on the right edge.
%     yRight = slope * (xRight - x1) + y1;
%     plot([xLeft, xRight], [yLeft, yRight], 'LineWidth',2,'Color','green');
%     
%     
%     % Plot original points on the lines
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%     
%     % determine which pixels are intersected by the lines
%     % across a matrix with the same dimensions as the image tile
%     x1 = xLeft:xRight;
%     dx = p2(1) - p1(1);
%     dy = p2(2) - p1(2);
%     y1 = round((x1 - p1(1)) * dy / dx + p1(2));
%     
%     outOfRange = y1 <= 0 | y1 > rows;
%     y1(outOfRange) = [];
%     x1(outOfRange) = [];
%     
%     idx = sub2ind(size(m), y1, x1);
%     m(idx) = 1; %m(idx) = m(idx) + 1;
%     
%     
%     
% end
% 
% figure; imagesc(m); colormap gray
% 
% se = strel('disk',5);
% dilateM = imdilate(m,se);
% figure; imshow(dilateM)
% 
% 
% se = strel('disk',10);
% morphM = imclose(dilateM,se);
% figure; imshow(morphM)
% 
% % apply the mask to the input image to keep only power line areas
% keepAreas = morphM > 0;
% filteredbw = bw;
% filteredbw(~keepAreas) = 0;
% figure; imshow(filteredbw)
% title('Power line output after filtering')
% 
% filteredbw = bw;
% filteredbw(keepAreas) = 0;
% figure; imshow(filteredbw)
% title('Points within these cells must change classification')

end 