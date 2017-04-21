function [y,ydec]=inflinearsvmval_manybands(InfoFeatures,svm)

% USAGE
% [y,yval]=inflinearsvmval_manybands(x,svm)
%
% svmval computes the prediction of a support vector machine
%       using the kernel function and its parameter for classification
%		  or regression
%
% INPUT
% InfoFeatures    : structure containing the data (in .X) and labels (in .Y)
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
%   01/05/2014 D. Tuia (modified)



%
% Usual verifications
%

if ~isfield(svm,'kernel'); svm.kernel='linear'; end
if ~isfield(svm,'kerneloption'); svm.kerneloption=1; end
if ~isfield(svm,'span'); svm.span=[]; end
if ~isfield(svm,'framematrix'); svm.framematrix=[]; end
if ~isfield(svm,'vector'); svm.vector=[]; end
if ~isfield(svm,'dual'); svm.dual=[]; end
if ~isfield(svm,'multiclass'); svm.multiclass=0; end

if svm.multiclass==0
    
    if isempty(svm.w)
        y=zeros(InfoFeatures.Sz(1)*InfoFeatures.Sz(2),1);
    else
    Xf=get_feature_test_manybands(InfoFeatures.X,svm.feat,InfoFeatures.Sz);%ici ca gene?re les features pour tester.
  %  Xff = Xf(InfoFeatures.YY(InfoFeatures.ct),:);
  %  save(['Xf_test_class' num2str(InfoFeatures.cl) '.mat'],'Xff');
    
    y=Xf*svm.w+svm.w0;
    end
    
    ydec=sign(y);

else

    y=zeros(InfoFeatures.Sz(1)*InfoFeatures.Sz(2),svm.nbclass);
    switch svm.multiclass

        case 1 % one against all


            for i=1:svm.nbclass
                InfoFeatures.cl = i;
                
                y(:,i)=inflinearsvmval_manybands(InfoFeatures,svm.models{i});

            end

            [temp,ydec]=max(y,[],2);
            ydec=svm.vals(ydec);

    end



end
