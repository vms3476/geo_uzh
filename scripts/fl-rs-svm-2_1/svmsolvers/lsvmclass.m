function [w,w0,obj]=lsvmclass(X0,Y,C,opts)

global X;
X=X0;

if nargin<4
    opts=[];
end

if nargin>4
  %  opts.w=w;
  %  opts.w0=w0;
    
end
opts.prec=1e-8;%1e-12;sqrt(eps)
opts.iter_max_Newton = 50;
[w, w0,obj]=primal_svm(1,Y,1/C,opts); 

obj=obj*C;
clear X