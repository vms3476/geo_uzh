function [ X ] = computeFeatures( las, res )
% calculate features for input las struct

% Inputs: 
%
%   las - struct containing las data 
%
%   res - struct containing x,y fields with raster dimensions and
%   resolution%
%
%
% 
% Outputs: 
%   
%   X - training data, each row is a different sample, 
%       each column is a feature
%


tic

% feature 1: height difference vector
    x = las.x';
    y = las.y';
    z = las.z';
    rnnr = las.rnnr';
    % remove middle returns
    i = rnnr == 32 | rnnr == 42 | rnnr == 43 | rnnr == 52 | rnnr == 53 | rnnr == 54 | ...
    rnnr == 62 | rnnr == 63 | rnnr == 64 | rnnr == 65 | rnnr == 72 | rnnr == 73 | ...
    rnnr == 74 | rnnr == 75 | rnnr == 76 | rnnr == 11; 
    x(i) = [];
    y(i) = [];
    z(i) = [];
    rnnr(i) = [];

    % find indices of first returns for first-last echo pairs
    % only keep negative differences
    drnnr = diff(rnnr);
    dz = diff(z);
    ii = (drnnr == 1 | drnnr == 2 | drnnr == 3 | drnnr == 4 | drnnr == 5 | drnnr == 6) & dz<0; 
    xDif = x(ii);
    yDif = y(ii); 
    zDif = -dz(ii); 

    % calculate raster, count the total number of zDif values ('den')
    ras1mDtm = raw2ras([xDif;yDif;zDif]',res,1,'dtm'); ras1mDtmInterp = inpaint_nans(ras1mDtm.z,4);
    X(:,1) = ras1mDtmInterp(:); 

    toc
    disp('feature 1 computed')
    
% feature 2: average intensity of single echos 
    i = las.rnnr == 11;
    x = las.x(i);
    y = las.y(i);
    z = las.z(i);
    int = las.int(i);
    singleEchosInt = raw2ras([x,y,int],res,1,'int'); singleEchosInt.int = inpaint_nans(singleEchosInt.int,4);
    X(:,2) = singleEchosInt.int(:);

    toc
    disp('feature 2 computed')
    
% feature 3: number of single echos
    singleEchosDen = raw2ras([x,y,z],res,1,'den'); singleEchosDen.z = inpaint_nans(singleEchosDen.z,4);
    X(:,3) = singleEchosDen.z(:);

    toc
    disp('feature 3 computed')
    
% feature 4: median echo height 
    zMedian = raw2ras([las.x,las.y,las.z],res,1,'dtm'); zMedian.z = inpaint_nans(zMedian.z,4); 
    X(:,4) = zMedian.z(:);

    toc
    disp('feature 4 computed')
    
% feature 5: std echo height 
    X(:,5) = zMedian.std(:);

    toc
    disp('feature 5 computed')
    
% feature 6: max echo height 
    zMax = raw2ras([las.x,las.y,las.z],res,1,'dsm'); zMax.z = inpaint_nans(zMax.z,4); 
    X(:,6) = zMax.z(:);

    toc
    disp('feature 6 computed')


end

