function w=LSquareLq(v,q,lambda,options);

% USAGE
% 
% w=LSquareLq(v,q,lambda,options);
% 
% w=argmin 1/2 ||x-v||^2+ \lambda ||x||_q
%
% uses an IRLS approach
%
% this function is used by mtlmklsvm for computing
% the proximal operator of \ell_1-\ell_q penalty
%
% options.threshold : stopping criterion based on the variation of w

% A. Rakotomamonjy 13/08/2010

dim=size(v,1);
if isfield(options,'init')
    w=options.init;
else
    w=ones(dim,1);
end;

%-------------------------------------------
% Subgradient descent
%-------------------------------------------
% pas=1;
% k=1;
% grad=1;
% obj=0.5*norm(v-w,2).^2+ lambda*norm(w,q);
% while max(abs(grad*pas/sqrt(k)))>0.05 
%     pas=1;
%     normeq=norm(w,q);
%     subg=abs(w).^(q-1).*sign(w)./normeq.^(q-1);
%     subg(find(w==0))=0;
%     grad= - (v-w)+ lambda*subg;
%         waux=w - grad*pas/sqrt(k);
%    while obj<0.5*norm(v-waux,2).^2+ lambda*norm(waux,q);
%         pas=pas/2;
%         waux=w - grad*pas/sqrt(k);
%         
%    end;
%    %w=waux
%     obj=0.5*norm(v-w,2).^2+ lambda*norm(w,q);
%     %w=w - grad*pas/sqrt(k);
%     k=k+1;
% end;
% w
%--------------------------------------------
% IRLS with factored gradient
%--------------------------------------------
wold=w+1;
k=1;
while max(abs(w-wold)>options.threshold);
    normeq=norm(w,q).^(q-1);
    wold=w;
    w= v./(1+lambda/normeq.*abs(w).^(q-2));
    k=k+1;
end;
%k






