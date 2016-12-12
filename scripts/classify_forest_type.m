% forest type classification script
% utilizes firts-last echo height difference algorithm by Liang et al. 

% specify directory with LAZ and/or LAS file,  plus DTM geotiff per tile
lasDir = '/Users/scholl/geo_uzh/data/KantonAargau/batch3/';

% output filepath for geotiff class maps
outPath = '/Users/scholl/geo_uzh/data/KantonAargau/output/batch3/';

% convert from laz to las 
cd(lasDir)
files = dir('*.laz');
if numel(files) > 0 
    for i = 1:numel(files)
        lazName = [files(i).name];
        lasName = [lazName(1:end-1) 's'];
        if exist(lasName, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' lazName ' -o ' lasName]);
        end
    end
end

%% 

% morphological parameters
    % coniferous bwareaopen
        pix1c = 9;
        conn1c = 4;
    % deciduous bwareaopen
        pix1d = 9;
        conn1d = 4;
    % ground class closing
        pixg = 36;
        conng = 4;
    % structural element for smoothing
        se = strel('disk',3);
        
tic 

% for each las tile, classify forest type
cd(lasDir)
files = dir('*.las');
for i = 1:numel(files)
    lasName = [files(i).name];
    disp(['Currently processing ' lasName '...']);
    
    tile = files(i).name(1:end-4);
    las = readlas(lasName);
    
    % normalize PC to DTM    
    % read DTM
    disp('     Normalizing by DTM...')
    [dtm.z, R] = geotiffread([lasDir 'DTM_' tile '.tif']);
    [dtm.x,dtm.y] = pixcenters(R,R.RasterSize(1),R.RasterSize(2));
    [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y); 
    
    las.tz = las.z;
    las.z = las.tz - interp2(dtm.X,dtm.Y,dtm.z,las.x,las.y);
    
    % remove noisy values above 50m
    las = subsetraw(las,las.z<50);
    
    % variables of interest for entire tile 
    x = las.x';
    y = las.y';
    z = las.z';
    rnnr = las.rnnr';
    classification = double(las.Classification');

    % remove middle returns
    i = rnnr == 32 | rnnr == 42 | rnnr == 43 | rnnr == 52 | rnnr == 53 | rnnr == 54 | ...
    rnnr == 62 | rnnr == 63 | rnnr == 64 | rnnr == 65 | rnnr == 72 | rnnr == 73 | ...
    rnnr == 74 | rnnr == 75 | rnnr == 76 | rnnr == 11; 
    x(i) = [];
    y(i) = [];
    z(i) = [];
    rnnr(i) = [];
    classification(i) = [];

    % create raster for classification of ground
    disp('     Ground classification...')
    ras1mGrnd = raw2ras([x;y;classification]',1,1,'dsm');

    % find indices of first returns for first-last echo pairs
    % only keep negative differences
    drnnr = diff(rnnr);
    dz = diff(z);
    ii = (drnnr == 1 | drnnr == 2 | drnnr == 3 | drnnr == 4 | drnnr == 5 | drnnr == 6) & dz<0; 
    xDif = x(ii);
    yDif = y(ii); 
    zDif = -dz(ii); 

    % calculate raster, count the total number of zDif values ('den')
    disp('     Forest type classification...')
    ras1mDen = raw2ras([xDif;yDif;zDif]',ras1mGrnd,1,'den'); 

    % height difference threshold classificaiton 
    th1 = 0.5;  % percent of pixels within cell
    th2 = 15;   % significant height difference, in meters

    zTh2 = zDif; 
    xTh2 = xDif; 
    yTh2 = yDif; 

    g = (zTh2 <= th2); % index = not a significant height difference

    zTh2(g) = []; 
    xTh2(g) = []; 
    yTh2(g) = []; 

    % compute density raster using only points higher than th2
    % using the raster resolution and dimensions from ras1mDen so they can be
    % divided for a proportion
    ras1mDenTh2 = raw2ras([xTh2',yTh2',zTh2'],ras1mGrnd,1,'den'); 

    % divide to find the proportion of pixels 
    ras1mproportion = ras1mDenTh2.z ./ ras1mDen.z;

    % inpaint the proportion rasters
    ras1mproportion = inpaint_nans(ras1mproportion, 4); 

    % use proportion less than threshold 1 to classify
    ras1mClass = ones(size(ras1mDen.z));   % 1 is deciduous
    ras1mClass(ras1mproportion < th1) = 0; % 0 is coniferous

    % morpholgical processing
    disp('     Morphological processing...')
    % remove noisy pixels in coniferous regions
    ras1mClassMorph = bwareaopen(ras1mClass,pix1c,conn1c);
    % remove noisy pixels in deciduous regions
    ras1mClassMorphInv = bwareaopen(1-ras1mClassMorph,pix1d,conn1d); 
    % conect coniferous regions. coniferous = 0, deciduous = 1
    ras1mClassMorphClose = 1-imclose(ras1mClassMorphInv, se);
    % make binary image from ground classification
    ras1mGrnd.z(ras1mGrnd.z == 2) = 1; 
    ras1mGrnd.z(ras1mGrnd.z == 3) = 0;    
    grnd1mMorph = bwareaopen(ras1mGrnd.z,pixg,conng); 
    % smooth the ground mask
    grnd1mMorphClose = imclose(grnd1mMorph,se);
    grnd1mMorphTotal = grnd1mMorphClose;
    % create map of ground, deciduous, and coniferous
    class1mMap = ras1mClassMorphClose;
    class1mMap(ras1mClassMorphClose == 0) = 2; % coniferous class value of 2
    class1mMap(grnd1mMorphTotal == 1) = 0; % ground class value of 0

    % write geotiff
    disp('     Writing class map as geotiff...')
    wrt_geotiff_CH([outPath tile],ras1mDenTh2.x,ras1mDenTh2.y,class1mMap)

    % time elapsed
    timeElapsed = toc

end