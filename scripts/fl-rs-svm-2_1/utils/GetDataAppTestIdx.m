function  [idxs]=GetDataAppTestIdx(x,y,nbtrain,nbloop, classcode)

% [xapp,yapp,xtest,ytest,indice]=CreateDataAppTest(x,y,nbtrain, classcode)
%
% if nbtrain =[nbapppos nbappneg ] % we have a specific number of positive
% and negative examples.

if nargin <5
    classcode(1)=1;
    classcode(2)=-1;
end;



if length(nbtrain)==1;
for j=1:nbloop
    indice=[];
    indice.app=[];
    indice.test=[];
    nbclass=length(classcode);
    nbdata=length(y);
    %keyboard
    for i=1:nbclass;
        ind=find(y==classcode(i));
        nbclasscode_i=length(ind);
        ratioclasscode_i=nbclasscode_i/nbdata;
        aux=randperm(nbclasscode_i);
        nbtrainclasscode_i=round(ratioclasscode_i*nbtrain);
        indapp=ind(aux(1:nbtrainclasscode_i));
        indtest=ind(aux(nbtrainclasscode_i+1:end));
        %xapp=[xapp;x(indapp,:)];
        %yapp=[yapp;y(indapp,:)];
        %xtest=[xtest;x(indtest,:)];
        %ytest=[ytest;y(indtest,:)];
        indice.app=[indice.app;indapp];
        indice.test=[indice.test;indtest];
    end;
    idxs(j)=indice;
end
end;
