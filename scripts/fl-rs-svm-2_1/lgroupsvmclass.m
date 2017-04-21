function [svm,obj,LOG]=lgroupsvmclass(X,Y,options)
% USAGE [svm]=lgroupsvmclass(X,Y,options)
%
% linear svm for group norm
%   min  C/N\sum_i max(1-yi(x_i*w+w0))^2  + l1-l2(w)
%   w,w0
% INPUT
%
% Training set
%      X  		: input data
%      Y  		: output data
% parameters
%		options.C		: regularisation term
%       options.groups  : groups for the l1_l2 norm
%       options.active_set : use active sets (for large dimensionnal data)
%
% OUTPUT
%
% svm.w      weight
% svm.w0		bias
% svm.pos    position of Support Vector
% obj    Value of Objective function
%
%
% see also svmreg, svmkernel, svmval, svmclass


%	01/12/2009 R. Flamary 
%       -function creation
%       -use of accelerated gradient
%   05/11/2010 R. Flamary
%       - Add the Active set algorithm for large scale 


if ~isfield(options,'C'); options.C=1; end

if ~isfield(options,'groups'); options.groups=1:size(X,2); end
if ~isfield(options,'multiclass'); options.multiclass=0; end
if ~isfield(options,'active_set'); options.active_set=0; end

if length(unique(Y))>2
    options.multiclass=1;
end


if length(options.groups)~=size(X,2)
    error('wrong size for group vector')
end

if options.multiclass==0

    if nargin>3
        options.alphainit=alphainit;
    end

    options.C=options.C/length(Y);
    
    
    groups=options.groups;
    
    types=unique(options.groups);
    
    % put groups as 1,2 3,...,nbgroups
    for i=1:length(types)
       groups(groups==types(i))=i; 
    end
   
    options.groups=[options.groups(:);max(max(options.groups(:)))+1];
    
    data=options;
    
    data.X=[X, ones(size(X,1),1)];
    data.Y=Y;
    
    if ~isfield(options,'W0'); 
   % [w,w0] = lsvmclass(X,Y,options.C);
   % W0=[w;w0];
   % W0=randn(size(data.X,2),1);
    W0=zeros(size(data.X,2),1);
    else
     W0=options.W0;   
    end
    
    data.options=options;
    
    opts=options;
    
    opts.q_func=@q_func;
    opts.minq_func=@minq_func;
    opts.costfunc=@costfunc;
    
    [W,vars,obj,LOG]=nesterov_solver(W0,[],data,opts);

    svm.xsup=[];
    svm.w=W(1:end-1);
    svm.w0=W(end);
    svm.pos=[];
    svm.kerneloption=1;
    svm.kernel='linear';


else

    vals=sort(unique(Y),'ascend');

    if options.verbose~=0
        disp(['Multiclass SVM classification (' num2str(length(vals)) 'classes)' ]);
    end

    switch options.multiclass

        case 1 % one against all


            options2=options;
            options2.multiclass=0;

            for i=1:length(vals)


                if options.verbose~=0
                    disp(['Classification :' num2str(i) '/' num2str(length(vals)) ]);
                end;
                Yt=(Y==vals(i))-(Y~=vals(i));
                [models{i},obj(i),LOG{i}]=lgroupsvmclass(X,Yt,options2);

            end

            svm.multiclass=options.multiclass;
            svm.nbclass=length(vals);
            svm.models=models;
            svm.kerneloption=1;
            svm.kernel='linear';
            svm.vals=vals;




    end

end

end

function [costnew,varsnew]=q_func(W,Wt,vars,data)

    varsnew=vars;
    
    %V=W-gradf(W,vars,data)/vars.L;
    
    fw=data.C*1/2*sum(max(0,1-data.Y.*(data.X*Wt)).^2)+Omega12_group(W,data.groups);
    
    c2=sum( ( (W-Wt).*gradf(Wt,vars,data) )) + vars.L/2*sum((W-Wt).^2);
    
    costnew=fw+c2;   

end

function [costnew,varsnew]=costfunc(W,vars,data)

    varsnew=vars;

    costnew=data.C*1/2*sum(max(0,1-data.Y.*(data.X*W)).^2)+Omega12_group(W(1:end-1),data.groups(1:end-1));

end


function [W,obj,vars]=minq_func(Vt,vars,data)

lambda=1/vars.L;

W=zeros(size(Vt));
Vt=Vt-1/vars.L*gradf(Vt,vars,data);

vals=unique(data.groups);

for i=1:length(vals)
    if norm(Vt(data.groups==vals(i)))>lambda
        W(data.groups==vals(i))=(1-lambda/norm(Vt(data.groups==vals(i)))).*Vt(data.groups==vals(i));
    end
end

W(end)=Vt(end);

obj=0;

end

function Df=gradf(W,vars,data)

    Xt=data.X.*repmat(data.Y,1,size(data.X,2));
   % Df=-data.C*sum(data.X.*repmat(data.Y.*max(0,1-data.Y.*(data.X*W)),1,size(data.X,2)))';
   Df=-data.C*Xt'*max(0,1-data.Y.*(data.X*W));
end







