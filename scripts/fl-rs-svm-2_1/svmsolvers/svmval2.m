function [y,ydec]=svmval2(x,svm)

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
    
    if isfield(svm,'gamma')
        x=x*diag(svm.gamma);
    end

    if strcmp(svm.kernel,'linear')
        y=lsvmval(x,svm.w,svm.w0) ;
    else
        [y]=svmval(x,svm.xsup,svm.w,svm.w0,svm.kernel,svm.kerneloption,svm.span,svm.framematrix,svm.vector,svm.dual);
    end
    ydec=sign(y);

else

    y=zeros(size(x,1),svm.nbclass);
    switch svm.multiclass

        case 1 % one against all


            for i=1:svm.nbclass

                y(:,i)=svmval2(x,svm.models{i});

            end

            [temp,ydec]=max(y,[],2);
            ydec=svm.vals(ydec);

    end



end
