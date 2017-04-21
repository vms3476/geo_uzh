function [svm,obj]=lrclass_mc_spams(X,Y,options)
% USAGE [svm]=sslsvmclass(X,Y,options)
%
% linear Support vector machine for semi supervised CLASSIFICATION
% This routine classify the training set with a support vector machine
% using quadratic programming algorithm selected by options.solver
%
% INPUT
%
% Training set
%      X  		: input data (lengh of Y fisrt one are labeled)
%      Y  		: output data
% parameters
%		options.lambdar		:  regularization term for the regularization
%		options.regs   : Type of regularization (local/global/none)
%
%
%
%
% OUTPUT (for one class)
%  learning with polynomial kernel
%
% svm.xsup
% svm.w      linear svm weight
% svm.w0	linear svm	bias
% svm.pos    position of Support Vector
% obj    Value of Objective function
%
%
% see also svmreg, svmkernel, svmval, svmclass


%	11/02/2012 R. Flamary





if ~isfield(options,'lambda'); options.lambda=1; end
if ~isfield(options,'nbitermax'); options.nbitermax=100; end
if ~isfield(options,'stopnorm'); options.stopnorm=1e-8; end
if ~isfield(options,'reg'); options.reg='l1'; end
if ~isfield(options,'eps'); options.eps=1e-8; end
if ~isfield(options,'nu'); options.nu=1; end

obj=struct();


% data matrix initialization
n=size(X,1);




Xapp=X;
Yapp=Y;

% multiclass initialization

vals=sort(unique(Yapp));
nbclass=length(vals);

if ~isfield(options,'obj') 
    W=zeros(size(X,2)+1,nbclass);
else % option warm start: le mex bugge!
    W = options.obj.W;
end

param=options;
param.lambda=options.lambda;
param.loss='multi-logistic';
param.regul=options.reg;
param.intercept=1;

W=mexFistaFlat(Y-1,[X ones(n,1)],W,param);

obj.W=W;

% setting output struct
if nbclass>2
    
    models=cell(nbclass);
    
    for i=1:nbclass
        models{i}.xsup=[];
        models{i}.w=W(1:end-1,i);
        models{i}.w0=W(end,i);
        models{i}.pos=[];
        models{i}.kerneloption=1;
        models{i}.kernel='linear';
        models{i}.span=[];
        
    end
    
    svm.multiclass=1;
    svm.nbclass=nbclass;
    svm.models=models;
    svm.kerneloption=1;
    svm.kernel='linear';
    svm.vals=vals;
    
    
        
    svm.multiclass=1;
    svm.nbclass=length(vals);
    svm.w=W(1:end-1,:);
    svm.w0=W(end,:);
    svm.vals=vals;
    
    
else
    
    svm=struct();
    
    svm.xsup=[];
    svm.w=W(1:end-1,2);
    svm.w0=W(end,2);
    svm.pos=[];
    svm.kerneloption=1;
    svm.kernel='linear';
    svm.span=[];
    
end


end








