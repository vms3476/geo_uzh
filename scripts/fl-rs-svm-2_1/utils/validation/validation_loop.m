function [model,options,LOG]=validation_loop(X,Y,Xval,Yval,options,opts)
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
if ~isfield(opts,'criterions'); opts.criterions={'err','min';'auc','max'}; end
if ~isfield(opts,'keeplog'); opts.keeplog=0; end
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

if opts.verbose
    disp('Validation loop')
end

% loop

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
                    model=opts.learning(X,Y,options2);

                    % testing
                    [ypred]=opts.valmodels(Xval,model);
                    restemp=get_precision(Yval,ypred);
                    err(i,j,k,l)=restemp.(opts.criterions{1,1});
                    auc(i,j,k,l)=restemp.(opts.criterions{2,1});

                    % warm start!!
                    if ~isempty(opts.warmstart)
                        options2=opts.warmstart(X,options2,model);
                    end
                end
            end
        end
    end


%errmin=min(min(min(min(err))));
eval(['errmin=' opts.criterions{1,2} '(' opts.criterions{1,2} '( ' opts.criterions{1,2} '(' opts.criterions{1,2} '( err))));']);
%aucmax=max(max(max(max(auc))));

idx=find(err==errmin);

if length(idx) >1
    auctemp=auc(idx);
    eval([ '[aucmax2,idmax2]=' opts.criterions{2,2} '(auctemp);']);
    idx=idx(idmax2);
end

[i,j,k,l]=ind2sub([nbi,nbj,nbk,nbl],idx);

eval(['options.' opts.params{1} '=' num2str(opts.values{1}(i)) ';']);
if nbparams>1,eval(['options.' opts.params{2} '=' num2str(opts.values{2}(j)) ';']); end
if nbparams>2,eval(['options.' opts.params{3} '=' num2str(opts.values{3}(k)) ';']); end
if nbparams>3,eval(['options.' opts.params{4} '=' num2str(opts.values{4}(l)) ';']); end
if opts.keeplog~=0 options.log=1; end


t0=cputime;

if opts.calc_final~=0
    if opts.keeplog~=0

        [model,LOG2]=opts.learning([X;Xval],[Y;Yval],options);
        LOG.LOG2=LOG2;

    else
        model=opts.learning([X;Xval],[Y;Yval],options);
    end
end
LOG.dt=cputime-t0;
eval(['LOG.' opts.criterions{1,1} '=err;']);
eval(['LOG.' opts.criterions{2,1} '=auc;']);

end




