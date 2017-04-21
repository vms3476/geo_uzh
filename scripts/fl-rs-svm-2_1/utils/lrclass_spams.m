function [svm,obj]=lrclass_spams(X,Y,options)
% USAGE [svm]=svmclass2(X,Y,options)
%
% Support vector machine for CLASSIFICATION
% This routine classify the training set with a support vector machine
% using quadratic programming algorithm selected by options.solver
%
% INPUT
%
% Training set
%      X  		: input data
%      Y  		: output data
% parameters
%		options.lambda		: regularization parameter
%
% 		options.verbose : display outputs (default value is 0: no display)
%
%      options.Span    : span matrix for semiparametric learning
%               This vector is sized Nbapp*Nbbasisfunction where
%               phi(i,j)= f_j(x(i));
%
%
%
%
%
% OUTPUT
%
% svm.xsup	coordinates of the Support Vector
% svm.w      weight
% svm.w0		bias
% svm.pos    position of Support Vector
% obj    Value of Objective function
%
%
% see also svmreg, svmkernel, svmval, svmclass


%	01/12/2009 R. Flamary




if ~isfield(options,'lambda'); options.lambda=1e0; end
if ~isfield(options,'verbose'); options.verbose=0; end
if ~isfield(options,'nobias'); options.nobias=0; end
if ~isfield(options,'multiclass'); options.multiclass=0; end
if ~isfield(options,'Display'); options.Display='off'; end
if ~isfield(options,'weighting'); options.weighting='none'; end
if ~isfield(options,'weights'); options.weights=ones(size(Y)); end

if length(unique(Y))>2
    options.multiclass=1;
end

n=size(X,1);



if options.multiclass==0
    
    
    
    vals=sort(unique(Yapp));
nbclass=length(vals);

W=zeros(size(X,2)+1,nbclass);


param=options;
param.lambda=options.lambda;
param.loss='logistic';
param.regul=options.reg;
param.intercept=1;

W=mexFistaFlat(Y,[X ones(n,1)],W,param);

obj.W=W;

    
    if options.nobias
        w=W;
        w0=0;
        
    else
        w=W(1:end-1);
        w0=W(end);
    end
    
    
    svm.xsup=[];
    svm.w=w;
    svm.w0=w0;
    svm.pos=[];
    svm.kerneloption=1;
    svm.kernel='linear';
    %svm.span=options.span;
    
    
    
    
    
else
    
    vals=sort(unique(Y),'ascend');
    nbclass=length(vals);
    temp=zeros(size(Y));
    P=-ones(n,nbclass);
    for i=1:length(vals)
        temp(Y==vals(i))=i;
        P(Y==vals(i),i)=1;
    end
    Y=temp;
    
    
param=options;
param.lambda=options.lambda;
param.loss='logistic';
param.regul=options.reg;
param.intercept=1;

W=zeros(size(X,2)+1,nbclass);

% 1vsAll

    W=mexFistaFlat(P,[X ones(n,1)],W,param);


obj.W=W;

    
    
    if options.nobias
        w=[W];
        w0=0;
        
    else
        w=[W(1:end-1,:)]; 
        w0=[W(end,:)];
    end

    
    %todo
    
    svm.multiclass=1;
    svm.nbclass=length(vals);
    svm.w=w;
    svm.w0=w0;
    svm.vals=vals;
    
    
    
    
end

end


