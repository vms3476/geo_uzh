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

% feature 1: median height difference
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

    % calculate raster, median zDif values ('dtm')
    ras1mDtm = raw2ras([xDif;yDif;zDif]',res,1,'dtm'); ras1mDtmInterp = inpaint_nans(ras1mDtm.z,4);
    
    % add to X array 
    %X(:,1) = normalize_unity(ras1mDtmInterp); % normalize
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
    %X(:,2) = normalize_unity(singleEchosInt.int);
    X(:,2) = singleEchosInt.int(:);

    toc
    disp('feature 2 computed')
    
% feature 3: number of single echos
    singleEchosDen = raw2ras([x,y,z],res,1,'den'); singleEchosDen.z = inpaint_nans(singleEchosDen.z,4);
    %X(:,3) = normalize_unity(singleEchosDen.z);
    X(:,3) = singleEchosDen.z(:);

    toc
    disp('feature 3 computed')
    
% feature 4: median echo height 
    zMedian = raw2ras([las.x,las.y,las.z],res,1,'dtm'); zMedian.z = inpaint_nans(zMedian.z,4); 
    %X(:,4) = normalize_unity(zMedian.z);
    X(:,4) = zMedian.z(:);

    toc
    disp('feature 4 computed')
    
% feature 5: std echo height 
    %X(:,5) = normalize_unity(zMedian.std);
    X(:,5) = zMedian.std(:);
    X(isnan(zMedian.std),5) = 0;

    toc
    disp('feature 5 computed')
    
% feature 6: max echo height 
    zMax = raw2ras([las.x,las.y,las.z],res,1,'dsm'); zMax.z = inpaint_nans(zMax.z,4); 
    %X(:,6) = normalize_unity(zMax.z);
    X(:,6) = zMax.z(:);

    toc
    disp('feature 6 computed')

% feature  7: vertical histogram    
   rasV = raw2vox([las.x,las.y,las.z],res,1,'dsm');  
   X(:,7:56) = reshape(rasV.vHist,numel(res.x)*numel(res.y),size(rasV.vHist,3));
   
   toc
   disp('feature 7 computed')
    
% feature 8: surrounding pixels 
   
    
    
    
% % feature 7: average height difference 
%     ras1mAvgz = raw2ras([xDif;yDif;zDif]',res,1,'int'); ras1mAvgzInterp = inpaint_nans(ras1mAvgz.int,4);
%     
%     % normalize and add to X array 
%     X(:,7) = normalize_unity(ras1mAvgzInterp); 
% 
%     toc
%     disp('feature 7 computed')
    
    
    
    
% % plots 
% figure; myimage(ras1mDtm.x,ras1mDtm.y,ras1mDtmInterp);
% figure; myimage(ras1mDtm.x,ras1mDtm.y,singleEchosInt.int);
% figure; myimage(ras1mDtm.x,ras1mDtm.y,singleEchosDen.z);
% figure; myimage(ras1mDtm.x,ras1mDtm.y,zMedian.z);
% figure; myimage(ras1mDtm.x,ras1mDtm.y,zMedian.std);
% figure; myimage(ras1mDtm.x,ras1mDtm.y,zMax.z);
% figure; myimage(ras1mDtm.x,ras1mDtm.y,ras1mAvgz.int);

    

end

