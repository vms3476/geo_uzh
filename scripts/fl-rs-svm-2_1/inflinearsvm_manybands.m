function [svm,obj,LOG,perIter]=inflinearsvm_manybands(xapp,Y,options,optionsFeatSel,InfoFeatures)

% inputs:   - xapp : train samples
%           - Y : train labels

% generate initial test
% YY = find(InfoFeatures.Y > 0);
% xtest = InfoFeatures.xapp(YY(cv),:);% to be updated in the loop
% YYY = InfoFeatures.Y(YY);
% ytest = YYY((cv));


if ~isfield(options,'C'); options.C=1; end
if ~isfield(options,'verbose'); options.verbose=1; end
if ~isfield(options,'A_0'); options.A_0=[]; end
if ~isfield(options,'epsilon'); options.epsilon=eps; end
if ~isfield(options,'constraintmethod'); options.constraintmethod='morpho'; end
if ~isfield(options,'multiclass'); options.multiclass=0; end


obj=options.C/2;

if length(unique(Y))>2
    options.multiclass=1;
end


if options.multiclass==0
    
    
    perIter{1,1} = InfoFeatures.cl;
    jjj=2;
    
    LOG=struct;
    w0=mean(Y);
    
    loop=1;
    nbloop=1;
    nbverbose=1;
    
    options.options_ag.C=options.C;
    
    verbose=options.verbose;
    
    Xf=xapp; %base features given to the classifier
    wf=ones(size(xapp,2),1); %initial dummy weights fo the first svm
    options.nbgroups=size(Xf,2);
    options.groups=(1:size(Xf,2))';
    nbactive=size(Xf,2);
    c = 1;
    %     rand('state',0);
    %ActiveFeat=[];
    for i = 1:size(Xf,2)
        ActiveFeat(i).filters = 'band';
        ActiveFeat(i).bands   = i;
        ActiveFeat(i).scales  = [];
        ActiveFeat(i).angles  = [];
        ActiveFeat(i).formel  = '';
        ActiveFeat(i).pos     = '';
        
    end
    %!!! inizializzare ActiveFeat con le bande iniziali (vedi Xf = X;)
    
    Xviol = [];
    poss = [];
    [wf,w0,obj]=solve_set(Xf,Y,wf,w0,options);% finds new weights for active features
    err=max(0,(1-Y.*(Xf*wf+w0)));
    
    options.bandweights = ones(1,size(xapp,2));
%   options.bandweights = double(abs(wf)~=0)';

    counter = 1;
    
    while loop  %stops when it does not find a violating constraint.
        %maybe to rework to switch from a type of features to
        %another
        
        
        %if nbloop>1
        if ~isempty(Xf) & ~isempty(Xviol)
            [wf,w0,obj]=solve_set(Xf,Y,wf,w0,options);% finds new weights for active features
            err=max(0,(1-Y.*(Xf*wf+w0)));
        end
        
        if mod(nbloop,options.when2save) == 0
            svm = [];
            svm.xsup=[];
            svm.w=wf;
            svm.w0=w0;
            svm.pos=[];
            svm.kerneloption=1;
            svm.kernel='linear';
            svm.feat=ActiveFeat;
            svm.weightVector = options.bandweights;
            
            perIter{1,jjj} = svm;
            clear svm
            jjj =jjj+1;
        end
        
        %else
        %    err=ones(size(Y));
        %end
        
        
        
        %while isempty(Xviol) && counter < options.nbitergenfeat
        if isempty(Xviol) || counterAdded == options.maxPerBank
                      
            
            [X,params] = generateFeatures(reshape(InfoFeatures.X,InfoFeatures.Sz(1),InfoFeatures.Sz(2),InfoFeatures.Sz(3)),options);
            poss = params.bands;
            X=(X-ones(size(X,1),1)*mean(X))./(ones(size(X,1),1)*std(X));
            counterAdded = 0; % counter of features added for THIS filterbank
        end
        
        [Xviol,Info,maxiviol,indViol]=FindViolatingConstraintsManybands(X(InfoFeatures.YY(InfoFeatures.ct),:),Y,err,options,optionsFeatSel,params);
        disp(num2str(double(isempty(Xviol))));
        counterAdded = counterAdded + 1;
        
        %end
        
        if isempty(Xviol)
            counter = counter +1;
        else
            counter = 1;
            
           
            
        end
        
        
        
        
        
        
        
        % clean the features with w = 0
        if ~isempty(Xf) & ~isempty(Xviol)
            indpareil=  find(sum(abs(Xf-Xviol*ones(1,size(Xf,2))))==0);
            if length(indpareil)==1
                % keyboard
            end;
        end;
        
        if ~isempty(Xviol);
            % prune groups and features
            k=1;
            for i=1:options.nbgroups
                indg=find(options.groups==i);
                if norm(wf(indg))>1e-6 %if weight is positive, keep it
                    options.groups(indg)=k;
                    
                    k=k+1;
                else %remove groups with wf < 1e-6 (useless groups)
                    options.groups(indg)=[];
                    wf(indg)=[];
                    Xf(:,indg)=[];
                    ActiveFeat(indg)=[];
                    
                end;
            end;
            
            options.nbgroups=k; 
            
            groups=options.groups;
            % group cleaning
            types=unique(options.groups);
            % put groups as 1,2 3,...,nbgroups = renumber groups
            for i=1:length(types)
                groups(groups==types(i))=i;
            end
            options.groups=[groups;ones(size(Xviol,2),1)*options.nbgroups];
            
            Xf=[Xf Xviol];
            wf=[wf;size(Xviol,2)];
            nbactive=size(Xf,2);
            ActiveFeat(nbactive)=Info;
            X(:,indViol) = [];
        %else
            %break; %break procedure if it does not find any violating constraint
        end;
        
        if verbose
            nbverbose = nbverbose+1;
            if nbverbose == 20 || (nbverbose==1 && nbloop==1);
                fprintf(1,[repmat('    ',1,verbose-1) '| Iter |   Cost  |   Delta Cost  | nb active\n']);
                nbverbose = 0;
            end
            %             J=cout(Xf,Y,wf,w0,options);
            J=obj;
            %obj=J;
            fprintf(1,[repmat('    ',1,verbose-1) '|%5d| ViCo %5.3f | %5.3f|%5d| \n'],[nbloop J , maxiviol,length(wf)-1]);
        end
        
        nbloop=nbloop+1;
        
        if nbloop > options.nbitermax % exit if max # of iters exceeded
            loop=0;
            disp('Exit by nbitermax')
        end;
        
        if counter > options.nbitergenfeat %exit if counter generations of filters returned empty Xviol
            loop=0;
            disp('Exit by nbitergenfeat')
        end;
        

        
    end
    
    %--------------------------------------------------------------------------
    %
    %           Exit Main Loop
    %
    %--------------------------------------------------------------------------
    %
    

    
    [wf,w0,obj]=solve_set(Xf,Y,wf,w0,options);% finds final weights for active features
    
    

    
    if verbose
        nbverbose = nbverbose+1;
        if nbverbose == 20 || (nbverbose==1 && nbloop==1);
            fprintf(1,[repmat('    ',1,verbose-1) '| Iter |   Cost  |   Delta Cost  | nb active\n']);
            nbverbose = 0;
        end
        if ~isempty(wf)
            J=cout(Xf,Y,wf,w0,options);
            fprintf(1,[repmat('    ',1,verbose-1) '|%5d| SVM %5.3f | %5.3f|%5d| \n'],[nbloop J , maxiviol,length(wf)]);
        end
    end
    
    
    
    
    % final set...
    obj=cout(Xf,Y,wf,w0,options);
    svm.xsup=[];
    svm.w=wf;
    svm.w0=w0;
    svm.pos=[];
    svm.kerneloption=1;
    svm.kernel='linear';
    %     if nbloop==1
    svm.feat=ActiveFeat;
    %     else
    %         svm.feat=ActiveFeat;
    %     end
    % save(['Xf_train_class' num2str(InfoFeatures.cl) '.mat'],'Xf');
    svm.featWeights = options.bandweights;
    
else % multiclass setting
    
    
    vals=sort(unique(Y),'ascend');
    
    if options.verbose~=0
        disp(['Multiclass SVM classification (' num2str(length(vals)) 'classes)' ]);
    end
    
    switch options.multiclass
        
        case 1 % one against all
            
            
            options2=options;
            options2.multiclass=0;
            
            for i=1:length(vals) %loop on classes
                InfoFeatures.cl = i;
                
                if options.verbose~=0
                    disp(['Classification :' num2str(i) '/' num2str(length(vals)) ]);
                end;
                Yt=(Y==vals(i))-(Y~=vals(i));
                [models{i},obj(i),LOG{i},perIterations{i}]=inflinearsvm_manybands(xapp,Yt,options2,optionsFeatSel,InfoFeatures);
                
            end
            
            svm.multiclass=options.multiclass;
            svm.nbclass=length(vals);
            svm.models=models;
            svm.kerneloption=1;
            svm.kernel='linear';
            svm.vals=vals;
            
            perIter = perIterations;
            
            
            
    end
    
    
    
end

end

%--------------------------------------------------------------------------
%
%       Subfunctions
%
%--------------------------------------------------------------------------


function  [w,w0,obj]=solve_set(X,Y,w_1,w0_1,options)

options2=options.options_ag;
options2.verbose=0;
%options.W0=[w_1;w0_1];
[svm,obj]=lgroupsvmclass(X,Y,options2);
w=svm.w;%weights
w0=svm.w0;%bias
end




function cost=cout(Xf,Y,w,w0,options)

if isempty(w)
    err=length(Y);
else
    err=max(0,1-Y.*(Xf*w+w0));
end
cost=options.C/(2*max(1,size(Xf,1)))*err'*err+Omega12_group(w,options.groups);

end





