% 5 AOI for power line testing 

% AOI 1
xMin = 668200; xMax = 668300; yMin = 255100; yMax = 255300; 
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi1/';
raw1 = getrawlas(lasDir,area);

% AOI 2
xMin = 643630; xMax = 643780; yMin = 265600; yMax = 265700;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi2/';
raw2 = getrawlas(lasDir,area);

% AOI 3
xMin = 646091; xMax = 646200; yMin = 265996; yMax = 266136;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi3/';
raw3 = getrawlas(lasDir,area);

% AOI 4
xMin = 670000; xMax = 670348; yMin = 239489; yMax = 239680;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi4/';
raw4 = getrawlas(lasDir,area);

% AOI 5
xMin = 643000; xMax = 643178; yMin = 240590; yMax = 240748;
area = [xMin xMax yMin yMax]; 
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/powerlines/aoi5/';
raw5 = getrawlas(lasDir,area);

%% plot las

figure; myscatter3(raw1.x,raw1.y,raw1.z,raw1.z,parula); title('AOI 1'); colorbar;
figure; myscatter3(raw2.x,raw2.y,raw2.z,raw2.z,parula); title('AOI 2'); colorbar;
figure; myscatter3(raw3.x,raw3.y,raw3.z,raw3.z,parula); title('AOI 3'); colorbar;
figure; myscatter3(raw4.x,raw4.y,raw4.z,raw4.z,parula); title('AOI 4'); colorbar;
figure; myscatter3(raw5.x,raw5.y,raw5.z,raw5.z,parula); title('AOI 5'); colorbar;

%% standard deviation raster

res = 1;

ras1 = raw2ras([raw1.x,raw1.y,raw1.z],res,res,'dsm');
ras2 = raw2ras([raw2.x,raw2.y,raw2.z],res,res,'dsm');
ras3 = raw2ras([raw3.x,raw3.y,raw3.z],res,res,'dsm');
ras4 = raw2ras([raw4.x,raw4.y,raw4.z],res,res,'dsm');
ras5 = raw2ras([raw5.x,raw5.y,raw5.z],res,res,'dsm');

%% std image
figure;myimage(ras1.x,ras1.y,ras1.std); title('AOI 1'); colorbar; % 6
figure;myimage(ras2.x,ras2.y,ras2.std); title('AOI 2'); colorbar; % 12?
figure;myimage(ras3.x,ras3.y,ras3.std); title('AOI 3'); colorbar; % 10
figure;myimage(ras4.x,ras4.y,ras4.std); title('AOI 4'); colorbar; % 4?
figure;myimage(ras5.x,ras5.y,ras5.std); title('AOI 5'); colorbar; % 10?

%% std hist

figure; h1 = histogram(ras1.std(:)); title('AOI 1');
figure; h2 = histogram(ras2.std(:)); title('AOI 2');
figure; h3 = histogram(ras3.std(:)); title('AOI 3');
figure; h4 = histogram(ras4.std(:)); title('AOI 4');
figure; h5 = histogram(ras5.std(:)); title('AOI 5');

%% height hist

figure; h1 = histogram(ras1.z(:)); title('AOI 1');
figure; h2 = histogram(ras2.z(:)); title('AOI 2');
figure; h3 = histogram(ras3.z(:)); title('AOI 3');
figure; h4 = histogram(ras4.z(:)); title('AOI 4');
figure; h5 = histogram(ras5.z(:)); title('AOI 5');

%% smoothed histogram of height

minHeight = 5; % below which, confident there are no power lines

[f,xi] = ksdensity(ras1.z(ras1.z>minHeight));
figure;plot(xi,f); title('AOI 1');
[peaks, heights] = findpeaks(f,xi);
zThresh1 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras2.z(ras2.z>minHeight));
figure;plot(xi,f); title('AOI 2');
[peaks, heights] = findpeaks(f,xi);
zThresh2 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras3.z(ras3.z>minHeight));
figure;plot(xi,f); title('AOI 3');
[peaks, heights] = findpeaks(f,xi);
zThresh3 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras4.z(ras4.z>minHeight));
figure;plot(xi,f); title('AOI 4');
[peaks, heights] = findpeaks(f,xi);
zThresh4 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras5.z(ras5.z>minHeight));
figure;plot(xi,f); title('AOI 5');
[peaks, heights] = findpeaks(f,xi);
zThresh5 = heights(find(ismember(peaks,max(peaks))))

%% highlight cells above z thresh 

k = ras1.z > zThresh1; 
ras1z = ras1.z;
ras1z(k) = 1; 
ras1z(~k) = 0; 
figure; myimage(ras1.x,ras1.y,ras1z); title('AOI 1');

k = ras2.z > zThresh2; 
ras2z = ras2.z;
ras2z(k) = 1; 
ras2z(~k) = 0; 
figure; myimage(ras2.x,ras2.y,ras2z); title('AOI 2');

k = ras3.z > zThresh3; 
ras3z = ras3.z;
ras3z(k) = 1; 
ras3z(~k) = 0; 
figure; myimage(ras3.x,ras3.y,ras3z); title('AOI 3');

k = ras4.z > zThresh4; 
ras4z = ras4.z;
ras4z(k) = 1; 
ras4z(~k) = 0; 
figure; myimage(ras4.x,ras4.y,ras4z); title('AOI 4');

k = ras5.z > zThresh5; 
ras5z = ras5.z;
ras5z(k) = 1; 
ras5z(~k) = 0; 
figure; myimage(ras5.x,ras5.y,ras5z); title('AOI 5');

%% smoothed histogram of std

minStd = 1

[f,xi] = ksdensity(ras1.std(ras1.std>minStd));
figure;plot(xi,f); title('AOI 1');
[peaks, heights] = findpeaks(f,xi);
stdThresh1 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras2.std(ras2.std>minStd));
figure;plot(xi,f); title('AOI 2');
[peaks, heights] = findpeaks(f,xi);
stdThresh2 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras3.std(ras3.std>minStd));
figure;plot(xi,f); title('AOI 3');
[peaks, heights] = findpeaks(f,xi);
stdThresh3 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras4.std(ras4.std>minStd));
figure;plot(xi,f); title('AOI 4');
[peaks, heights] = findpeaks(f,xi);
stdThresh4 = heights(find(ismember(peaks,max(peaks))))

[f,xi] = ksdensity(ras5.std(ras5.std>minStd));
figure;plot(xi,f); title('AOI 5');
[peaks, heights] = findpeaks(f,xi);
stdThresh5 = heights(find(ismember(peaks,max(peaks))))

%% highlight cells above std thresh 

k = ras1.std > stdThresh1; 
ras1std = ras1.std;
ras1std(k) = 1; 
ras1std(~k) = 0; 
figure; myimage(ras1.x,ras1.y,ras1std); title('AOI 1');

k = ras2.std > stdThresh2; 
ras2std = ras2.std;
ras2std(k) = 1; 
ras2std(~k) = 0; 
figure; myimage(ras2.x,ras2.y,ras2std); title('AOI 2');

k = ras3.std > stdThresh3; 
ras3std = ras3.std;
ras3std(k) = 1; 
ras3std(~k) = 0; 
figure; myimage(ras3.x,ras3.y,ras3std); title('AOI 3');

k = ras4.std > stdThresh4; 
ras4std = ras4.std;
ras4std(k) = 1; 
ras4std(~k) = 0; 
figure; myimage(ras4.x,ras4.y,ras4std); title('AOI 4');

k = ras5.std > stdThresh5; 
ras5std = ras5.std;
ras5std(k) = 1; 
ras5std(~k) = 0; 
figure; myimage(ras5.x,ras5.y,ras5std); title('AOI 5');

%% remove power lines 

res = 1; 
stdThresh = 7; 
peri = 6; 
nClusters = 3; 
densityThr = 0.1;


% AOI 1
[raw_pl1, las1] = remove_power_lines(raw1, res, stdThresh1, peri, nClusters, densityThr)
figure; myscatter3(raw_pl1.x,raw_pl1.y,raw_pl1.z,raw_pl1.Classification,parula); colorbar; swisstick;
title(['AOI 1, res: ',num2str(res),', stdThresh: ', num2str(stdThresh1),', peri: ',num2str(peri),', density: ', num2str(densityThr)])

% AOI 2
[raw_pl2, las2] = remove_power_lines(raw2, res, stdThresh2, peri, nClusters, densityThr)
figure; myscatter3(raw_pl2.x,raw_pl2.y,raw_pl2.z,raw_pl2.Classification,parula); colorbar; swisstick;
title(['AOI 2, res: ',num2str(res),', stdThresh: ', num2str(stdThresh2),', peri: ',num2str(peri),', density: ', num2str(densityThr)])

% AOI 3
[raw_pl3, las3] = remove_power_lines(raw3, res, stdThresh3, peri, nClusters, densityThr)
figure; myscatter3(raw_pl3.x,raw_pl3.y,raw_pl3.z,raw_pl3.Classification,parula); colorbar; swisstick;
title(['AOI 3, res: ',num2str(res),', stdThresh: ', num2str(stdThresh3),', peri: ',num2str(peri),', density: ', num2str(densityThr)])

% AOI 4
[raw_pl4, las4] = remove_power_lines(raw4, res, stdThresh4, peri, nClusters, densityThr)
figure; myscatter3(raw_pl4.x,raw_pl4.y,raw_pl4.z,raw_pl4.Classification,parula); colorbar; swisstick;
title(['AOI 4, res: ',num2str(res),', stdThresh: ', num2str(stdThresh4),', peri: ',num2str(peri),', density: ', num2str(densityThr)])

% AOI 5
[raw_pl5, las5] = remove_power_lines(raw5, res, stdThresh5, peri, nClusters, densityThr)
figure; myscatter3(raw_pl5.x,raw_pl5.y,raw_pl5.z,raw_pl5.Classification,parula); colorbar; swisstick;
title(['AOI 5, res: ',num2str(res),', stdThresh: ', num2str(stdThresh5),', peri: ',num2str(peri),', density: ', num2str(densityThr)])



%% REIK method of Max (R) - Mean (G) - Min (B) 

mean1 = raw2ras([raw1.x,raw1.y,raw1.z],res,res,'int');
mean2 = raw2ras([raw2.x,raw2.y,raw2.z],res,res,'int');
mean3 = raw2ras([raw3.x,raw3.y,raw3.z],res,res,'int');
mean4 = raw2ras([raw4.x,raw4.y,raw4.z],res,res,'int');
mean5 = raw2ras([raw5.x,raw5.y,raw5.z],res,res,'int');

min1 = raw2ras([raw1.x,raw1.y,raw1.z],res,res,'dmn');
min2 = raw2ras([raw2.x,raw2.y,raw2.z],res,res,'dmn');
min3 = raw2ras([raw3.x,raw3.y,raw3.z],res,res,'dmn');
min4 = raw2ras([raw4.x,raw4.y,raw4.z],res,res,'dmn');
min5 = raw2ras([raw5.x,raw5.y,raw5.z],res,res,'dmn');


% display RGB image 

rgbTest = cat(3,uint8(ras1.z),uint8(mean1.int),uint8(min1.z));
figure; imshow(rgbTest); title('AOI 1 Max (R) - Mean (G) - Min (B)','FontSize',14);

rgbTest = cat(3,uint8(ras2.z),uint8(mean2.int),uint8(min2.z));
figure; imshow(rgbTest); title('AOI 2 Max (R) - Mean (G) - Min (B)','FontSize',14);

rgbTest = cat(3,uint8(ras3.z),uint8(mean3.int),uint8(min3.z));
figure; imshow(rgbTest); title('AOI 3 Max (R) - Mean (G) - Min (B)','FontSize',14);

rgbTest = cat(3,uint8(ras4.z),uint8(mean4.int),uint8(min4.z));
figure; imshow(rgbTest); title('AOI 4 Max (R) - Mean (G) - Min (B)','FontSize',14);

rgbTest = cat(3,uint8(ras5.z),uint8(mean5.int),uint8(min5.z));
figure; imshow(rgbTest); title('AOI 5 Max (R) - Mean (G) - Min (B)','FontSize',14);


% display ratios based on RGB appearance 
% display ratios based on RGB appearance 
ratio1 = ras1.z ./ min1.z;
ratio1(ratio1<0)=0;
ratio1(ratio1>50)=0;
figure; myimage(ras1.x,ras1.y,ratio1); title('AOI 1');

% display ratios based on RGB appearance 
ratio2 = ras2.z ./ mean2.int;
ratio2(ratio2<0)=0;
ratio2(ratio2>50)=0;
figure; myimage(ras2.x,ras2.y,ratio2); title('AOI 2');