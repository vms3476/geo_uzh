function res = mc_mlr_as(xapp,yapp,InfoFeatures,options)


res = {};

loop = 1;
iter = 1;
nbloop = 1; %total number of loops
numEmptyXviol = 1;%number of times a given set of features has not given any interesting feature
counterAdded = 0; % number of features added form a single minibatch

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
    ActiveFeat(i).refined     = 0;
    ActiveFeat(i).level     = 0;
end

Xviol = [];

%train the logisitc regression the first time
[model,obj]=lrclass_mc_spams(Xf,yapp,options);
R=lr_mc_residue(Xf,yapp,model);
res{1,iter}.model = model;
res{1,iter}.ActiveFeat = ActiveFeat;

fprintf('old size: %i ',size(Xf,2));

list = [];
for i =1:size(Xf,2)
    if norm(model.w(i,:))<1e-6 %if the weight of the feature is positive
        list = [list i];
        
        
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
       
       
       tt = sprintf('L1L2. Iter: %i, Size: %i, OA: %1.3f, Kappa: %1.3f.',iter,size(Xf,2),res{1,iter}.num.OA,res{1,iter}.num.Kappa)
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
       [candFeat,params] = generateFeatures(reshape(InfoFeatures.X,InfoFeatures.Sz(1),InfoFeatures.Sz(2),InfoFeatures.Sz(3)),options);
       candFeat=(candFeat-ones(size(candFeat,1),1)*mean(candFeat))./(ones(size(candFeat,1),1)*std(candFeat));    
       
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
    if m > options.lambda*(1 + options.tolerance)
       disp('found interesting feature.')
       Xviol = candFeat(InfoFeatures.YY(InfoFeatures.ct),ind); 
       XviolT = candFeat(:,ind);
       Info = getFeatInfo(params,ind)
       
       
       
       counterAdded = counterAdded + 1;
    else
        Xviol = [];
        Info = [];
       
    end
    
    %% include violating features (both in train and test), clean and prepare for next loop
    
    if ~isempty(Xviol)
       % a feature has been selected, so we put the counter for consecutive bad minibatches to 0     
       numEmptyXviol = 1;
       
       %cleaning 0 weight features in the current model
       list = [];
       for i =1:size(Xf,2)
           if norm(model.w(i,:))<1e-6 %if the weight of the feature is positive
               %model.w(i)=[]; %remove weight from model
               list = [list i];
               
           end
       end
       
       Xf(:,list)=[];%remove feature
       xt(:,list)=[];
       ActiveFeat(list)=[];% remove from list of params
       model.w(list,:) = [];
       obj.W(list,:) = [];
       for i = 1:model.nbclass
           model.models{i,1}.w(list) = [];
       end
     
       Xf=[Xf Xviol];
       fprintf('new size: %i\n',size(Xf,2));
       xt = [xt XviolT];
       clear XviolT
       
       wf=[wf;size(Xviol,2)];
       nbactive=size(Xf,2);
       ActiveFeat(nbactive)=Info;
       candFeat(:,ind) = [];
       
       
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