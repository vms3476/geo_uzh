function [model,options,LOG]=bootstrap_validation_loop(X,Y,options,opts)
%[model,options,LOG]=validation_loop(X,Y,Xval,Yval,options,opts)
% validation loop for selecting parameters (4 max)
%
% Input parameters:
%
% X,Y,Xval,Yval : learning and validation sets
%
% options : options struct to be used in the learning
%
% opts.params : cell containing the name of the params to find
% opts.values : cell containig the values of the params
% opts.learning : [model]=learning(X,Y,options) learning function
% opts.valmodels : [ypred]=valmodels(X,model) valu calculating function
% opts.warmstart : [options2]=warmstart(X,options,model) warm start
% opts.verbose : verbose level for the validation loop
% opts.prop     : proportion of training example
% opts.nbloop   : number of cross validation loop
%
% Output parameters:
%
% model : model learned on [X;Xval] with optimal params
% options : options used for model learning (optimal params)
% LOG : log of the loop

%  Author : Remi Flamary remi (dot) flamary (at) gmail (dot) com

if ~isfield(opts,'verbose'); opts.verbose=0; end
if ~isfield(opts,'warmstart'); opts.warmstart=[]; end
if ~isfield(opts,'parallel'); opts.parallel=0; end
if ~isfield(opts,'prop'); opts.prop=2/3; end
if ~isfield(opts,'nbloop'); opts.nbloop=5; end
if ~isfield(opts,'criterions'); opts.criterions={'err','min';'auc','max'}; end
if ~isfield(opts,'get_data'); opts.get_data=@get_data_rand; end
if ~isfield(opts,'get_data_info'); opts.get_data_info=[]; end
if ~isfield(opts,'calc_final'); opts.calc_final=1; end

model=[];
LOG=[];

nbparams=length(opts.params);

nbi=length(opts.values{1});
if nbparams>1, nbj=length(opts.values{2}); else nbj=1;end
if nbparams>2, nbk=length(opts.values{3}); else nbk=1;end
if nbparams>3, nbl=length(opts.values{4}); else nbl=1;end

err=zeros(nbi,nbj,nbk,nbl);
auc=zeros(nbi,nbj,nbk,nbl);

options2=options;

if options.verbose
    disp('Cross Validation loop')
end

% loop

for m=1:opts.nbloop
   if options.verbose 
     disp(['loop ' num2str(m) '/' num2str(opts.nbloop)])
    end
 
   prop=opts.prop;
   
   [X2,Y2,X2val,Y2val]=opts.get_data(X,Y,m,opts);
    
    options2=options;
    for i=1:nbi
        for j=1:nbj
            for k=1:nbk

                for l=1:nbl


                    options2.(opts.params{1})=opts.values{1}(i);
                    if nbparams>1, options2.(opts.params{2})=opts.values{2}(j); end
                    if nbparams>2, options2.(opts.params{3})=opts.values{3}(k); end
                    if nbparams>3, options2.(opts.params{4})=opts.values{4}(l); end
                    % learning
                    model=opts.learning(X2,Y2,options2);

                    % testing
                    [ypred]=opts.valmodels(X2val,model);
                    restemp=get_precision(Y2val,ypred);
                    err(i,j,k,l,m)=restemp.(opts.criterions{1,1});
                    auc(i,j,k,l,m)=restemp.(opts.criterions{2,1});

                    % warm start!!
                    if ~isempty(opts.warmstart)
                        options2=opts.warmstart(X2,options2,model);
                    end
                end
            end
        end
    end

end

eval(['errmin=' opts.criterions{1,2} '(' opts.criterions{1,2} '( ' opts.criterions{1,2} '(' opts.criterions{1,2} '( mean(err,5)))));']);

%errmin=min(min(min(min(mean(err,5)))));
%aucmax=max(max(max(max(auc))));

err2=mean(err,5);
auc2=mean(auc,5);

idx=find(err2==errmin);

if length(idx) >1
    auctemp=auc2(idx);
    %[aucmax2,idmax2]=max(auctemp);
    eval([ '[aucmax2,idmax2]=' opts.criterions{2,2} '(auctemp);']);
    idx=idx(idmax2);
end

[i,j,k,l]=ind2sub([nbi,nbj,nbk,nbl],idx);

eval(['options.' opts.params{1} '=' num2str(opts.values{1}(i)) ';']);
if nbparams>1,eval(['options.' opts.params{2} '=' num2str(opts.values{2}(j)) ';']); end
if nbparams>2,eval(['options.' opts.params{3} '=' num2str(opts.values{3}(k)) ';']); end
if nbparams>3,eval(['options.' opts.params{4} '=' num2str(opts.values{4}(l)) ';']); end

if opts.calc_final~=0
model=opts.learning([X],[Y],options);
end

LOG.(opts.criterions{1,1})=err;
eval(['LOG.' opts.criterions{1,1} '2=err2']);
LOG.(opts.criterions{2,1})=auc;
eval(['LOG.' opts.criterions{2,1} '2=auc2']);
end


function [X2,Y2,X2val,Y2val]=get_data_rand(X,Y,niter,opts)
prop=opts.prop;

  ind_app = [];
   ind_test = [];
   c_lab= unique(Y');
   for c=c_lab
       ind_c = find(Y==c);
       melange = randperm(length(ind_c));
       ind_app_c = ind_c(melange (1: floor(length(ind_c)*prop) ) );
       ind_test_c = ind_c( melange ( floor(length(ind_c)*prop)  + 1 : length(ind_c) ));
       ind_app = [ ind_app ; ind_app_c];
       ind_test = [ ind_test ; ind_test_c];
   end
   
   X2=X(ind_app,:);
   Y2=Y(ind_app);
   X2val=X(ind_test,:);
   Y2val=Y(ind_test);

end

function [X2,Y2,X2val,Y2val]=get_data_strat(X,Y,niter,opts)
prop=opts.prop;

  ind_app = [];
   ind_test = [];
   c_lab= unique(Y');
   for c=c_lab
       ind_c = find(Y==c);
       melange = randperm(length(ind_c));
       ind_app_c = ind_c(melange (1: floor(length(ind_c)*prop) ) );
       ind_test_c = ind_c( melange ( floor(length(ind_c)*prop)  + 1 : length(ind_c) ));
       ind_app = [ ind_app ; ind_app_c];
       ind_test = [ ind_test ; ind_test_c];
   end
   
   X2=X(ind_app,:);
   Y2=Y(ind_app);
   X2val=X(ind_test,:);
   Y2val=Y(ind_test);

end




