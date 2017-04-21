function [ X ] = computeFeatures( las, res, showPlot, inpaintNans)
% calculate features for input las struct

% Inputs: 
%
%   las - struct containing las data 
%
%   res - struct containing x,y fields with raster dimensions and
%   resolution
%
%   showPlot = logical value to show feature plots (1) or not (0) 
%
%   inpaintNans = logical value to inpaint_nans (1) or set NanN = 0 (0)
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
    ras1mDtm = raw2ras([xDif;yDif;zDif]',res,1,'dtm'); 
    
    if inpaintNans
        ras1mDtmInterp = inpaint_nans(ras1mDtm.z,4);
    else
        ras1mDtm.z(isnan(ras1mDtm.z)) = 0; ras1mDtmInterp = ras1mDtm.z; 
    end
    
    X(:,1) = ras1mDtmInterp(:);
    toc
    disp('feature 1 computed')
    

    
% feature 2: average intensity of single echos 
    i = las.rnnr == 11;
    x = las.x(i);
    y = las.y(i);
    z = las.z(i);
    int = las.int(i);
    singleEchosInt = raw2ras([x,y,int],res,1,'int'); 
    
    if inpaintNans
        singleEchosInt.int = inpaint_nans(singleEchosInt.int,4);
    else
        singleEchosInt.int(isnan(singleEchosInt.int)) = 0;
    end
    
    X(:,2) = singleEchosInt.int(:);
    toc
    disp('feature 2 computed')
  
    
    
% feature 3: fraction of single echos
    singleEchoesDen = raw2ras([x,y,z],res,1,'den'); 
    totalEchoesDen = raw2ras([las.x,las.y,las.z],res,1,'den'); 
    singleEchoFraction = singleEchoesDen.z ./ totalEchoesDen.z;
    
    if inpaintNans
        singleEchoFraction = inpaint_nans(singleEchoFraction,4);
    else
        singleEchoFraction(isnan(singleEchoFraction)) = 0;
    end
    
    X(:,3) = singleEchoFraction(:);
    toc
    disp('feature 3 computed')
    
    
    
% feature 4: median echo height 
    zMedian = raw2ras([las.x,las.y,las.z],res,1,'dtm'); 
    
    if inpaintNans
        zMedian.z = inpaint_nans(zMedian.z,4); 
    else
        zMedian.z(isnan(zMedian.z)) = 0;
    end
    
    X(:,4) = zMedian.z(:);
    toc
    disp('feature 4 computed')
    
    
    
% feature 5: std echo height 
    zStd = zMedian.std;
    
    if inpaintNans
        zStd = inpaint_nans(zStd,4);
    else
        zStd(isnan(zStd(:))) = 0;
    end
    
    X(:,5) = zMedian.std(:);
    toc
    disp('feature 5 computed')
    
    
    
% feature 6: max echo height 
    zMax = raw2ras([las.x,las.y,las.z],res,1,'dsm'); 
    
    if inpaintNans
        zMax.z = inpaint_nans(zMax.z,4); 
    else
        zMax.z(isnan(zMax.z))=0;
    end 
    
    X(:,6) = zMax.z(:);
    toc
    disp('feature 6 computed')

    
    
% feature 7: median intensity of all echos 
    x = las.x;
    y = las.y;
    int = las.int; 
    allEchoesMedianInt = raw2ras([x,y,int],res,1,'dtm'); 
    
    if inpaintNans
        allEchoesMedianInt.z = inpaint_nans(allEchoesMedianInt.z,4);
    else
        allEchoesMedianInt.z(isnan(allEchoesMedianInt.z)) = 0;
    end
    
    X(:,7) = allEchoesMedianInt.z(:);
    toc
    disp('feature 7 computed')
        
    
    
% feature 8: standard deviation of intensity of all echoess 
    allEchoesStdInt = allEchoesMedianInt.std;
    
    if inpaintNans
        allEchoesStdInt = inpaint_nans(allEchoesStdInt,4);
    else
        allEchoesStdInt(isnan(allEchoesStdInt)) = 0; 
    end 
    
    X(:,8) = allEchoesStdInt(:); 
    toc
    disp('feature 8 computed')
    
    
    
% Feature 9: fraction of ground echoes
    j = las.Classification == 2; 
    x = las.x(j);
    y = las.y(j);
    z = las.z(j);
    groundEchoesDen = raw2ras([x,y,z],res,1,'den'); 
    groundEchoFraction = groundEchoesDen.z ./ totalEchoesDen.z;
    
    if inpaintNans
        groundEchoFraction = inpaint_nans(groundEchoFraction,4);
    else
        groundEchoFraction(isnan(groundEchoFraction)) = 0;
    end
    
    X(:,9) = groundEchoFraction(:);
    toc
    disp('feature 9 computed')

    
    
% % feature  10: vertical histogram    
%    rasV = raw2vox([las.x,las.y,las.z],res,1,'dsm');  
%    
%    
%    X(:,10:59) = reshape(rasV.vHist,numel(res.x)*numel(res.y),size(rasV.vHist,3));
%    toc
%    disp('feature 10 computed')
    
    
    
    
    
if showPlot 
	figure; myimage(ras1mDtm.x,ras1mDtm.y,ras1mDtmInterp); colorbar; title('Feature 1: median first-last height difference per pulse');
	figure; myimage(ras1mDtm.x,ras1mDtm.y,singleEchosInt.int); colorbar; title('Feature 2: average intensity of single echoes');
	figure; myimage(ras1mDtm.x,ras1mDtm.y,singleEchoFraction); colorbar; title('Feature 3: fraction of single echoes');
	figure; myimage(ras1mDtm.x,ras1mDtm.y,zMedian.z); colorbar; title('Feature 4: median echo height');
	figure; myimage(ras1mDtm.x,ras1mDtm.y,zMedian.std); colorbar; title('Feature 5: standard deviation of echo height');
	figure; myimage(ras1mDtm.x,ras1mDtm.y,zMax.z); colorbar; title('Feature 6: maximum echo height');
	figure; myimage(ras1mDtm.x,ras1mDtm.y,allEchoesMedianInt.z); colorbar; title('Feature 7: median intensity of all echoes');
    figure; myimage(ras1mDtm.x,ras1mDtm.y,allEchoesStdInt); colorbar; title('Feature 8: standard deviation of intensity (all echoes)');
    figure; myimage(ras1mDtm.x,ras1mDtm.y,groundEchoFraction); colorbar; title('Feature 9: fraction of ground echoes');
end
    

end

