%test validation loop

clear all
close all

%%
%
%   Creating Data
%

rand('state',0);
randn('state',0);

nbapp=1000;
nn=20;
nbtest=nn*nn;
sigma=0.3;
[X,Y,Xval,Yval,xtest1,xtest2]=datasets('gaussian',nbapp,nbtest,sigma);
%------- Building a 2D Grid for function evaluation -------------------------
[xtest1 xtest2]  = meshgrid([-1:.05:1]*3.5,[-1:0.05:1]*3); 
nn = length(xtest1); 
xtest = [reshape(xtest1 ,nn*nn,1) reshape(xtest2 ,nn*nn,1)]; 


%%

clear options

% svm options
options.lambda = 1e-10;  
options.C = 1;
options.kernel='gaussian';
options.kerneloption=1;
options.verbose=1;
options.solver=3;

% validation options
opts.params{1}='C';
opts.values{1}=logspace(2,5,5);

opts.params{2}='kerneloption';
opts.values{2}=linspace(0.1,10,5);

opts.learning=@svmclass2;
opts.valmodels=@svmval2;
opts.warmstart=@warmstart_svmclass2;

tic
[svm,options,LOG]=validation_loop(X,Y,Xval,Yval,options,opts)
toc

[ypred] = svmval2(xtest,svm);
ypred = reshape(ypred,nn,nn); 

%%

clear options

% svm options
options.lambda = 1e-10;  
options.C = 1;
options.kernel='gaussian';
options.kerneloption=1;
options.verbose=1;
options.solver=3;


% validation options
opts.prop=2/3;
opts.nbloop=1;

opts.params{1}='C';
opts.values{1}=logspace(2,5,3);

opts.params{2}='kerneloption';
opts.values{2}=linspace(0.1,1,3);

opts.learning=@svmclass2;
opts.valmodels=@svmval2;
opts.warmstart=@warmstart_svmclass2;

tic
[svm,options,LOG]=cross_validation_loop(X,Y,options,opts)
toc

[ypred] = svmval2(xtest,svm);
ypred = reshape(ypred,nn,nn); 

%% visus


figure(1); 
clf; 
contourf(xtest1,xtest2,ypred,[-1 0 1]);hold on
h1=plot(X(Y==1,1),X(Y==1,2),'+r'); 
h2=plot(X(Y==-1,1),X(Y==-1,2),'+b'); 
h3=plot(svm.xsup(:,1),svm.xsup(:,2),'ok'); 

