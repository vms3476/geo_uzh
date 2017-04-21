function res = mc_mlr(xapp,yapp,InfoFeatures,options)

res = {};
iter = 1;


Xf=xapp; %base features given to the classifier

xt = InfoFeatures.x;
xt=(xt-ones(size(xt,1),1)*mean(xt))./(ones(size(xt,1),1)*std(xt));

yt = InfoFeatures.ym;

X = reshape(InfoFeatures.x,InfoFeatures.Sz(1),InfoFeatures.Sz(2),InfoFeatures.Sz(3));



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
    
end

Xviol = [];

%train the logisitc regression the first time
[model,obj]=lrclass_mc_spams(Xf,yapp,options);
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
res{1,iter}.model = model;
res{1,iter}.ActiveFeat = ActiveFeat;
res{1,iter}.totIterNum = 0;
[ypred,ypredc,yprob]=lrval(xt,model);
res{1,iter}.num = assessment(yt(yt>0),ypredc(yt>0),'class');