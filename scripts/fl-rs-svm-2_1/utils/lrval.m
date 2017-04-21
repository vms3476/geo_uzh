function [y,ydec,yprob]=lrval(x,svm)

% USAGE
% [y,yval]=svmval(x,svm)
%
% svmval computes the prediction of a support vector machine
%       using the kernel function and its parameter for classification
%		  or regression
%
% INPUT
% x    : input data
% svm.xsup : support vector list
% svm.w    : weight
% svm.kernel : string containing the type of kernel
% svm.kerneloption : setting parameter of kernel.
% svm.w0   : bias. this can be a column vector in case of semiparametric SVM
% svm.span : span matrix for semiparametric SVM
%
%   ----- 1D Frame Kernel --------------------------
%
%   svm.framematrix  frame elements for frame kernel
%   svm.vector       sampling position of frame elements
%	 svm.dual 		  dual frame
%
% OUTPUT
%
% ydec : sign of the output ouf the network at point (vector or matrix) x
%
%       y = w phi(x) - b*span(x)
%
%
%	See also svmclass,svmreg, svmkernel
%
%

%	01/12/2009 R. Flamary



%
% Usual verifications
%

if ~isfield(svm,'kernel'); svm.kernel='gaussian'; end
if ~isfield(svm,'kerneloption'); svm.kerneloption=1; end
if ~isfield(svm,'span'); svm.span=[]; end
if ~isfield(svm,'framematrix'); svm.framematrix=[]; end
if ~isfield(svm,'vector'); svm.vector=[]; end
if ~isfield(svm,'dual'); svm.dual=[]; end
if ~isfield(svm,'multiclass'); svm.multiclass=0; end

if svm.multiclass==0
    
    y=x*svm.w+svm.w0;
    ydec=sign(y);
    
    yprob=1./(1+exp(-y));

else


    
    y=x*svm.w+ones(size(x,1),1)*svm.w0;
    

    [temp,ydec]=max(y,[],2);
    ydec=svm.vals(ydec);
    
    yprob=1./(1+exp(-y));
    
    sy=sum(yprob,2);
    
    yprob=yprob./(sy*ones(1,size(svm.w,2)));



end


