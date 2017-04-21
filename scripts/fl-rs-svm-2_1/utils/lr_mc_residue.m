function [R]=lr_mc_residue(X,y,svm)
n=size(X,1);
t=max(y);
d=size(X,2);
id=sub2ind([n t],(1:n)',y);

M=X*svm.w+ones(size(X,1),1)*svm.w0;

py=M(id);

R=zeros(n,t);

dlt=zeros(n,t);
dlt(id)=1;

E=exp(M-py*ones(1,t));

SE=sum(E,2)*ones(1,t);

R=(E-dlt.*SE)./SE/n;
