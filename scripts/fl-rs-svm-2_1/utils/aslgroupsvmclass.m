function [svm,obj,LOG]=aslgroupsvmclass(X,Y,options)
% USAGE [svm]=svmclass2(X,Y,options)
%
% linear svm for group norm with active constraint algorithm
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



nbinit=5;

if ~isfield(options,'C'); options.C=1; end
if ~isfield(options,'verbose'); options.verbose=1; end
if ~isfield(options,'groups'); options.groups=1:size(X,2); end
if ~isfield(options,'multiclass'); options.multiclass=0; end
if ~isfield(options,'A_0'); options.A_0=[]; end
if ~isfield(options,'epsilon'); options.epsilon=eps; end


if length(unique(Y))>2
    options.multiclass=1;
end


if options.multiclass==0

        LOG=struct;
    
    if nargin>3
        options.alphainit=alphainit;
    end

    options.C=options.C;
    options.nbgroups=length(unique(options.groups));
   
    groups=options.groups;
    
    types=unique(options.groups);
    
    % put groups as 1,2 3,...,nbgroups
    for i=1:length(types)
       groups(groups==types(i))=i; 
    if length(options.A_0)>1
        options.A_0(options.A_0==types(i))=i;
    end
    end
    
    options.groups=groups;
      
    w=zeros(size(X,2),1);
    w0=mean(Y);
    
    if isempty(options.A_0)
   % score=constraint_violation(X,Y,w,w0,options);
   % [temp,idtemp]=sort(score,'descend');
   % A=idtemp(1:nbinit);
   A=[];
    else
      A =  options.A_0;
    end
   
    
    loop=1;
    nbloop=1;
    nbverbose=1;
    
    
    % active set algorithm
    
    verbose=options.verbose;
    
    while loop
        
        if verbose
            
            nbverbose = nbverbose+1;
            if nbverbose == 20 || (nbverbose==1 && nbloop==1);
                fprintf(1,[repmat('    ',1,verbose-1) '| Iter |   Cost  |   Delta Cost  | nb active\n']);
                nbverbose = 0;                      
            end

            J=cout(X,Y,w,w0,options,groups);
            if nbloop==1
                J_1=J;
            end

            fprintf(1,[repmat('    ',1,verbose-1) '|%5d| %5.3e | %1.4f |%5d| \n'],[nbloop J (J-J_1)/J_1, length(A)]);
        
            J_1=J;
            
        end
        
        % solve with set A
        if length(A)
            
            [w,w0]=solve_set(X,Y,w,w0,A,options);

            %pruning
            A=[];
            for i=1:options.nbgroups
                if norm(w(options.groups==i))>0
                    A=unique([A;i]);
                end

            end

        end
        % get the constraint violation
        score=constraint_violation(X,Y,w,w0,options);
        
        comp=1:options.nbgroups;
        comp(A)=[];
        
  
        [val,id]=max(score(comp));
                
        
        if val>0+options.epsilon
            A=unique([A;comp(id)]);
            id2=get_group_idxs(groups,id);
         %   w(id2)=randn(1);
        else
            loop=0;
        end
          
        
        
        
        nbloop=nbloop+1;
        
    end
    
          if verbose
            
            nbverbose = nbverbose+1;
            if nbverbose == 20 || (nbverbose==1 && nbloop==1);
                fprintf(1,[repmat('    ',1,verbose-1) '| Iter |   Cost  |   Delta Cost  | nb active\n']);
                nbverbose = 0;                      
            end

            J=cout(X,Y,w,w0,options,groups);
            if nbloop==1
                J_1=J;
            end

            fprintf(1,[repmat('    ',1,verbose-1) '|%5d| %5.3e | %1.4f |%5d| \n'],[nbloop J (J-J_1)/J_1, length(A)]);
        
            j_1=J;
            
        end
    
    
    
    
    obj=0;

    
    % final set...
    
    svm.xsup=[];
    svm.w=w;
    svm.w0=w0;
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


function  [w,w0]=solve_set(X,Y,w_1,w0_1,set,options)


idxs=get_group_idxs(options.groups,set);

groups=options.groups;
groups=groups(idxs);

options.groups=groups;
options.verbose=0;

% warm start
options.W0=[w_1(idxs);w0_1];

[svm]=lgroupsvmclass(X(:,idxs),Y,options);

w=zeros(size(w_1));
w(idxs)=svm.w;

w0=svm.w0;
end


function  score=constraint_violation(X,Y,w,w0,options)

score=zeros(options.nbgroups,1);

set=find(w);
w2=w(set);

Xt=X.*repmat(Y,1,size(X,2));

err=max(0,(1-Y.*(X(:,set)*w2+w0)));

for  i=1:options.nbgroups
    idx=zeros(size(options.groups));
    idx(options.groups==i)=1;
    idx=logical(idx);
   
        score(i)=norm(Xt(:,idx)'*err);
    
end

score=options.C/length(Y)*score-1;

end


function idxs=get_group_idxs(groups,set)

idxs=zeros(length(groups),1);

for i=1:length(set)
    idxs(groups==set(i))=1;
end

idxs=logical(idxs);

end

function cost=cout(X,Y,w,w0,options,groups)
set=(w~=0);

err=max(0,1-Y.*(X(:,set)*w(set)+w0));

cost=options.C/(2*size(X,1))*err'*err+Omega12_group(w,groups);

end


function Xf=get_feature(X,feat)


Xf=zeros(size(X,1),length(feat));

for i=1:length(feat)
    
    switch feat(i).type
        
        case 'orig'
            
            temp=X(:,feat(i).use);
        
        otherwise
            
            
    end
    
    
    Xf(:,i)=temp;
    
end



end


