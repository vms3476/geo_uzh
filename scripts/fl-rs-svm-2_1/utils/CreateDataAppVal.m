function   [xa,ya,xv,yv]=CreateDataAppVal(xapp,yapp,ratio);

n=size(xapp,1);
ind=randperm(n);
na=floor(n*ratio);
classcode=[1 -1];
[aux,aux,aux,aux,indice]=CreateDataAppTest(rand(n,2),yapp,na, classcode);


 xa=xapp(indice.app,:,:);
 ya=yapp(indice.app);
 xv=xapp(indice.test,:,:);
 yv=yapp(indice.test);