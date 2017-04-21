function [res,InfoFeatures,ActiveFeat,logF] = mc_mlr_as_deep3(xapp,yapp,InfoFeatures,options)


res = {};

loop = 1;
iter = 1;
nbloop = 1; %total number of loops
numEmptyXviol = 1;%number of times a given set of features has not given any interesting feature
counterAdded = 0; % number of features added form a single minibatch

options.lev = zeros(size(xapp,2),1); % to track depth level

Xf=xapp; %base features given to the classifier

xt = InfoFeatures.X;
xt=(xt-ones(size(xt,1),1)*mean(xt))./(ones(size(xt,1),1)*std(xt));

yt = InfoFeatures.ym;

X = reshape(InfoFeatures.X,InfoFeatures.Sz(1),InfoFeatures.Sz(2),InfoFeatures.Sz(3));



wf=ones(size(xapp,2),1); %initial dummy weights fo the first mlr
nbactive=size(Xf,2);

for i = 1:size(Xf,2)
    ActiveFeat(i).filters = 'band';
    ActiveFeat(i).bands   = i;
    ActiveFeat(i).scales  = [];
    ActiveFeat(i).angles  = [];
    ActiveFeat(i).formel  = '';
    ActiveFeat(i).pos     = '';
    ActiveFeat(i).refined = 0;
    ActiveFeat(i).level = 0;
    
    InfoFeatures.BandInfo(i) = ActiveFeat(i);
    
    
    
end

Xviol = [];

%train the logisitc regression the first time
[model,obj]=lrclass_mc_spams(Xf,yapp,options);
R=lr_mc_residue(Xf,yapp,model);
res{1,iter}.model = model;
res{1,iter}.ActiveFeat = ActiveFeat;

fprintf('old size: %i ',size(Xf,2));

list = [];
logF{nbloop}.in = [];
logF{nbloop}.out = [];

j=1;
for i =1:size(Xf,2)
    if norm(model.w(i,:))<1e-6 %if the weight of the feature is positive
        list = [list i];
        logF{nbloop}.out{j} = ActiveFeat(i);
        j=j+1;
        
    end
end


Xf(:,list)=[];%remove feature
xt(:,list)=[];
model.w(list,:) = [];
obj.W(list,:) = [];
for i = 1:model.nbclass
   model.models{i,1}.w(list) = [];  
end
ActiveFeat(list)=[];% remove from list of params

fprintf('new size: %i\n',size(Xf,2));


%[model,obj]=lrclass_mc_spams(Xf,yapp,options);
R=lr_mc_residue(Xf,yapp,model);
res{1,iter}.model = model;
res{1,iter}.ActiveFeat = ActiveFeat;
res{1,iter}.totIterNum = 0;
[ypred,ypredc,yprob]=lrval(xt,model);
res{1,iter}.num = assessment(yt(yt>0),ypredc(yt>0),'class');

%% main loop

while loop
    
    nbloop = nbloop +1;

    %% solve mlr, test and get the residual with current model
    if ~isempty(Xf) && ~isempty(Xviol)
       iter = iter +1; %iter only counts cycles where a feature has been added
       options.obj = obj;
       %options.obj.W=[model.w;zeros(1,length(Activefeat));model.w0]
       options.obj.W=[model.w;zeros(size(Xviol,2),max(InfoFeatures.Y(:)));model.w0];
       [model,obj]=lrclass_mc_spams(Xf,yapp,options);
       R=lr_mc_residue(Xf,yapp,model);
       
       res{1,iter}.model = model;
       res{1,iter}.ActiveFeat = ActiveFeat;
       res{1,iter}.totIterNum = nbloop; % at which cycle this change has happened (since there may be iterations without violating features)
        
    
    
       % validation with updated test set (= augmented + buffered)
       [ypred,ypredc,yprob]=lrval(xt,model);
       
       res{1,iter}.num = assessment(yt(yt>0),ypredc(yt>0),'class');
       
       
       ypredc(1) = 0;
       
       
       tt = sprintf('L1L2. realiz:  %1.3f, reg: %i,\n Iter: %i, Size: %i, OA: %1.3f, Kappa: %1.3f.',options.currentReal,options.lambda,iter,size(Xf,2),res{1,iter}.num.OA,res{1,iter}.num.Kappa)
       %        figure(99),imagesc(reshape(ypredc+1,InfoFeatures.Sz(1),InfoFeatures.Sz(2)));title(tt);axis image
       %        if isfield(options.cmap,'colors')
       %            colormap(options.cmap.colors)
       %        end
       %        drawnow
       %
       %        if iter == 2
       %            figure(101),imagesc(reshape(yt+1,InfoFeatures.Sz(1),InfoFeatures.Sz(2)));title('testGT'),axis image
       %            if isfield(options.cmap,'colors')
       %                colormap(options.cmap.colors)
       %            end
       %        end
       
      
    end
    
    %% extract candidate features (candFeat)
    if isempty(Xviol) || counterAdded == options.maxPerBank
        if counterAdded == options.maxPerBank
           numEmptyXviol = numEmptyXviol +1; 
        end
        
       counterAdded = 0;
       [candFeat,params,levs] = generateFeatures(reshape(InfoFeatures.X,InfoFeatures.Sz(1),InfoFeatures.Sz(2),InfoFeatures.Sz(3)),options);
       candFeat2 = candFeat;
       
       candFeat=(candFeat-ones(size(candFeat,1),1)*mean(candFeat))./(ones(size(candFeat,1),1)*std(candFeat));   
       
       if options.ASH
           penal = (options.gamma.^(levs))';
           candFeat = candFeat./(repmat(penal,size(candFeat,1),1));
       end
       
       %mmin = min(candFeat);            
       %candFeat = candFeat+repmat(mmin,size(candFeat,1),1);
       %candFeat = scale(candFeat);
       
       %[size(InfoFeatures.X) max(params.bands)]
       
       %        for i = 1: size(candFeat,2) %visu of thefilters
       %           figure(10)
       %           imagesc(reshape(candFeat(:,i),InfoFeatures.Sz(1),InfoFeatures.Sz(2))),axis equal
       %           pause
       %        end
        
    end
    
    %% assess features to include
   
    G=candFeat(InfoFeatures.YY(InfoFeatures.ct),:)'*R;
    KKT=sqrt(sum(G.^2,2));%*size(Xf,1);
    
    
    [m,ind] = max(KKT);
    L = options.lambda*(1 + options.tolerance);
        
    
    if m > L 
       disp('found interesting feature.')
       Xviol = candFeat(InfoFeatures.YY(InfoFeatures.ct),ind); 
       XviolT = candFeat(:,ind);
       Info = getFeatInfo(params,ind);

       
       if options.ASH == 1
            %Add the band to the input matrix (for re-filtering)
            InfoFeatures.X = [InfoFeatures.X candFeat2(:,ind)]; %in the new primitives, we add the band without zscore!
            InfoFeatures.Sz(3) = InfoFeatures.Sz(3) +1; 
            options.bandweights = ones(1,InfoFeatures.Sz(3));
            
            
            % Add info on the depth level
            switch params.filters
                case {'Ratios','NRatios','Sum','Prod'} % -> in these cases is the max between the level of the two bands combined
                    Info.level = max(InfoFeatures.BandInfo(Info.bands).level,InfoFeatures.BandInfo(Info.scales).level)+1;
                otherwise % -> in this case only the level of the parent
                    Info.level = InfoFeatures.BandInfo(Info.bands).level+1;
            end
            %Add on the global infos on Infofeatures.X
            InfoFeatures.BandInfo(size(InfoFeatures.X,2)) = Info;
            
            %Add to the vector of levels of depth
            options.lev = [options.lev; Info.level];
            Info
       
       end
       
       counterAdded = counterAdded + 1;
    else
        Xviol = [];
        Info = [];
        
        if options.ASH % add a bit of randomness. If no one violates, we still can add the feature to the bank (but not to the model!).
           proba =  rand(1);
           if proba > 0.9
              disp('Feature added even if not violating')
              
              InfoFeatures.X = [InfoFeatures.X candFeat2(:,ind)]; %in the new primitives, we add the band without zscore!
              InfoFeatures.Sz(3) = InfoFeatures.Sz(3) +1;
              Info = getFeatInfo(params,ind);
              options.bandweights = ones(1,InfoFeatures.Sz(3));
              
              switch params.filters
                  case {'Ratios','NRatios','Sum','Prod'} % -> in these cases is the max between the level of the two bands combined
                      Info.level = max(InfoFeatures.BandInfo(Info.bands).level,InfoFeatures.BandInfo(Info.scales).level)+1;
                  otherwise % -> in this case only the level of the parent
                      Info.level = InfoFeatures.BandInfo(Info.bands).level+1;
              end
              
              InfoFeatures.BandInfo(size(InfoFeatures.X,2)) = Info;
              
              %Add to the vector of levels of depth
              options.lev = [options.lev; Info.level];
              
           end
        end
       
    end
    
    %update logs
    logF{nbloop}.in = Info;
    logF{nbloop}.out = [];
    
    %% include violating features (both in train and test), clean and prepare for next loop
    
    if ~isempty(Xviol)
       % a feature has been selected, so we put the counter for consecutive bad minibatches to 0     
       numEmptyXviol = 1;
       
       %cleaning 0 weight features in the current model
       list = [];
       j=1;
       for i =1:size(Xf,2)
           
           if norm(model.w(i,:))<1e-6 %if the weight of the feature is positive
               %model.w(i)=[]; %remove weight from model
               list = [list i];
               logF{nbloop}.out{j} = ActiveFeat(i);
               j=j+1;
           end
       end
       
       Xf(:,list)=[];%remove feature
       xt(:,list)=[];
       
       if options.ASH 
        levs(ind) = [];
       end
       
       ActiveFeat(list)=[];% remove from list of params
       model.w(list,:) = [];
       obj.W(list,:) = [];
       for i = 1:model.nbclass
           model.models{i,1}.w(list) = [];
       end
     
       %Add to current active set
       Xf=[Xf Xviol];%in train
       fprintf('new size: %i\n',size(Xf,2));
       xt = [xt XviolT];%in test
       clear XviolT
       
       wf=[wf;size(Xviol,2)];%add a '1' to MLC parameters vector
       nbactive=size(Xf,2);
       ActiveFeat(nbactive)=Info;% add info in the active list
       
       candFeat(:,ind) = [];%erase the feature from the minibatch (in case it is reused)
       

       
    end
    

    
    
    
    
    
    %% exit conditions

        
        if nbloop > options.nbitermax % exit if max # of iters exceeded
            loop=0;
            disp('Exit by nbitermax')
        end;
        
        if numEmptyXviol > options.nbitergenfeat %exit if counter generations of filters returned too many empty Xviol
            loop=0;
            disp('Exit by nbitergenfeat')
        end;
    
    
    
    


end