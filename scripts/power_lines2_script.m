% This script is testing methods to identify power line ALS points,
% after it was observed that the initial attempt was not robust enough to
% successfully identify points in the Zurich AOI 1 test area. 


% Read LAS data and specify AOI coordinates

% ZURICH AOI 1
pcZ_all = LASread('/Users/scholl/geo_uzh/data/WSL/laz/zh.las',false,false);
areaZ = [710000, 710500, 240500, 241000]; 

% AARGAU AOI 1
pc1_all = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi1/n668000_255000.las',false,false);
area1 = [668000, 668500, 255200, 255700];

% AOI 2
pc2_all = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi2/n643000_265000.las',false,false);
area2 = [643600, 644000, 265300, 265700];

% AOI 3
pc3_all = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi3/n646000_265000.las',false,false);
area3 = [646000, 646500, 265500, 266000]; 

% AOI 4
pc4_all = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi4/n670000_239000.las',false,false);
area4 = [670000, 670300, 239400, 239800];

% AOI 5
pc5_all = LASread('/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi5/n643000_240000.las',false,false);
area5 = [643000, 643300, 240200, 240800];



% Run the method "las2ras" to generate rasters for each AOI 
% These rasters will be used to identify power line candidate cells
res = 3; 
[rasZ,pcZ] = las2ras(pcZ_all.record,areaZ,res);
[ras1,pc1] = las2ras(pc1_all.record,area1,res);
[ras2,pc2] = las2ras(pc2_all.record,area2,res);
[ras3,pc3] = las2ras(pc3_all.record,area3,res);
[ras4,pc4] = las2ras(pc4_all.record,area4,res);
[ras5,pc5] = las2ras(pc5_all.record,area5,res);


% write LAS subsets to file 
pc1_subset.record = subsetraw(pc1_all.record,pc1_all.record.x > area1(1) & pc1_all.record.x < area1(2) & pc1_all.record.y > area1(3) & pc1_all.record.y < area1(4));
LASwrite(pc1_subset,'/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi1/n668000_255000_subset.las','version', 12, 'verbose', false);

pc2_subset.record = subsetraw(pc2_all.record,pc2_all.record.x > area2(1) & pc2_all.record.x < area2(2) & pc2_all.record.y > area2(3) & pc2_all.record.y < area2(4));
LASwrite(pc2_subset,'/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi2/n643000_265000_subset.las','version', 12, 'verbose', false);

pc3_subset.record = subsetraw(pc3_all.record,pc3_all.record.x > area3(1) & pc3_all.record.x < area3(2) & pc3_all.record.y > area3(3) & pc3_all.record.y < area3(4));
LASwrite(pc3_subset,'/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi3/n646000_265000_subset.las','version', 12, 'verbose', false);

pc4_subset.record = subsetraw(pc4_all.record,pc4_all.record.x > area4(1) & pc4_all.record.x < area4(2) & pc4_all.record.y > area4(3) & pc4_all.record.y < area4(4));
LASwrite(pc4_subset,'/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi4/n670000_239000_subset.las','version', 12, 'verbose', false);

pc5_subset.record = subsetraw(pc5_all.record,pc5_all.record.x > area5(1) & pc5_all.record.x < area5(2) & pc5_all.record.y > area5(3) & pc5_all.record.y < area5(4));
LASwrite(pc5_subset,'/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi5/n643000_240000_subset.las','version', 12, 'verbose', false);

pcZ_subset.record = subsetraw(pcZ_all.record,pcZ_all.record.x > areaZ(1) & pcZ_all.record.x < areaZ(2) & pcZ_all.record.y > areaZ(3) & pcZ_all.record.y < areaZ(4));
LASwrite(pcZ_subset,'/Users/scholl/geo_uzh/data/WSL/laz/zh_subset.las','version', 12, 'verbose', false);


%% Plot and test individual AOIs

pc = pcZ; ras = rasZ; % ZH 1
pc = pc1; ras = ras1; % AG 1
pc = pc2; ras = ras2; % AG 2
pc = pc3; ras = ras3; % AG 3
pc = pc4; ras = ras4; % AG 4
pc = pc5; ras = ras5; pc_all = pc5_all.record;% AG 5 

%% plot 
figure; myscatter3(pc.x,pc.y,pc.z,pc.z,parula); title('PC')
figure;myimage(ras.x,ras.y,ras.z); title('DSM raster')
figure;myimage(ras.x,ras.y,ras.std); title('STD raster')
figure;myimage(ras.x,ras.y,ras.den); title('DEN raster')
figure;myimage(ras.x,ras.y,ras.dtm); title('DTM (Median) raster')
figure;myimage(ras.x,ras.y,ras.int); title('INT (Mean) raster')
figure;myimage(ras.x,ras.y,ras.r1); title('DEN first returns')
figure;myimage(ras.x,ras.y,ras.r11); title('DEN single returns')
figure;myimage(ras.x,ras.y,ras.std_wGround); title('STD with Ground returns')

% ratio between max height and number of returns
max_den_ratio = ras.z ./ ras.den;
figure;myimage(ras.x,ras.y,max_den_ratio);title('DSM/DEN')

% power line echo: first return / total returns 
powerline_echo = ras.r1 ./ ras.den; 
figure; myimage(ras.x,ras.y,powerline_echo); title('power line echo: r1 / rTotal')

pl_echo_thr = powerline_echo; pl_thresh = 0.9;
pl_echo_thr(pl_echo_thr<pl_thresh) = 0; pl_echo_thr(isnan(pl_echo_thr)) = 0;
pl_echo_thr(pl_echo_thr>=pl_thresh) = 1; 
figure; myimage(ras.x,ras.y,pl_echo_thr); title('power line echo thresholded')

% vegetation echo: first + intermediate / total returns 
ri = subsetraw(pc,pc.number_of_returns~=pc.return_number); 
rasri = raw2ras([ri.x,ri.y,ri.z],res,res,'den');
veg_echo = (ras.r1 + rasri.z) ./ ras.den; 
figure; myimage(ras.x,ras.y,veg_echo); title('veg echo: r1 + r11 / rTotal')

% building echo: single / total echos
r11 = subsetraw(pc,pc.return_number==1& pc.number_of_returns ==1);
rasr11 = raw2ras([r11.x,r11.y,r11.z],res,res,'den');
building_echo = (ras.r1 + rasr11.z) ./ ras.den;
figure; myimage(ras.x,ras.y,building_echo); title('building echo: r11 / rTotal')

% ratio std / den 
std_den_ratio = ras.std_wGround ./ ras.den; 
figure; myimage(ras.x,ras.y,std_den_ratio); colorbar;


%% voxel processing 


% % without ground
% pc = subsetraw(pcZ,pcZ.x >= min(pcZ.x) & pcZ.x < min(pcZ.x)+5 & pcZ.y >= 2.408021400000000e+05 & pcZ.y < (2.408021400000000e+05 + 5))
% figure; myscatter3(pc.x,pc.y,pc.z,pc.z,parula); title('PC')
% % with ground 
% pc = subsetraw(pcZ_all.record,pcZ_all.record.x >= min(pcZ_all.record.x) & pcZ_all.record.x < min(pcZ_all.record.x)+5 & pcZ_all.record.y >= 2.408021400000000e+05 & pcZ_all.record.y < (2.408021400000000e+05 + 5))
% figure; myscatter3(pc.x,pc.y,pc.z,pc.z,parula); title('PC')
pc = pcZ;
ras = rasZ;

% threshold power line echo raster
powerline_echo = ras.r1 ./ ras.den; 
pl_echo_thr = powerline_echo; pl_thresh = 0.9;
pl_echo_thr(pl_echo_thr<pl_thresh) = 0; pl_echo_thr(isnan(pl_echo_thr)) = 0;
pl_echo_thr(pl_echo_thr>=pl_thresh) = 1; 
figure; myimage(ras.x,ras.y,pl_echo_thr); title('power line echo thresholded')

% filter power line echo raster 
dilateM=bwareaopen(pl_echo_thr,25,4)
figure;myimage(ras.x,ras.y,dilateM)



% XY subset for resxresxres voxel 

xList = floor(min(pc.x)) : res : ceil(max(pc.x));
yList = floor(min(pc.y)) : res : ceil(max(pc.y));
zList = res : res : round(max(pc.z));

voxGrid.ht = zeros(numel(xList),numel(yList),numel(zList));

%%


for x = 1:numel(xList)
    x1 = xList(x);
    x2 = x1 + res;
    
    for y = 1:numel(yList)
        y1 = yList(y);
        y2 = y1 + res;
        
        
        
        % if power line candidate cell, do vertical analysis
        if pl_echo_thr(x,y)
            
            for z = 1:numel(zList)
                z1 = zList(z);
                z2 = z1 + res;

                vox = subsetraw(pc, pc.x >= x1 & pc.x < x2 & ...
                    pc.y >= y1 & pc.y < y2 & ...
                    pc.z >= z1 & pc.z < z2);
                
                
                % if voxel contains points
                if numel(vox.x) >= 0
                    %figure; myscatter3(vox.x,vox.y,vox.z,vox.z,parula); title([z ' vox PC'])
                    BW = raw2ras([vox.x,vox.y,vox.z],0.1,0.1,'den'); BW.z(BW.z>0)=1;
                    %figure;imshow(BW.z); title('Original power line output')
                    
                    [H,theta,rho] = hough(BW.z);
                    %[r,c] = find(H==max(H(:)))
                    
                    voxGrid.ht(x,y,z) = max(H(:)) ;
                  
                    
                    
                    % Compute eigenvalue 
                    
                    % compute density 
                    
                    % other 
                end
            end
            
        end
    end
end


% 3D plot of the HT values
[X,Y] = meshgrid(xList,yList);
figure; hold on;
for z = 1:numel(zList)
    htCounts = reshape(voxGrid.ht(:,:,z),numel(xList)*numel(yList),1);
	scatter3(X(:),Y(:),htCounts)

end


 



%% Test power line identification script 

power_lines2(pcZ,rasZ); title('ZRH 1')

power_lines2(pc1,ras1); title('AG 1')

power_lines2(pc2,ras2); title('AG 2')

power_lines2(pc3,ras3); title('AG 3')

power_lines2(pc4,ras4); title('AG 4')

power_lines2(pc5,ras5); title('AG 5')



