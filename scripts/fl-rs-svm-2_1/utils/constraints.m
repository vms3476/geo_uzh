
function  score=constraint_violation(X,Y,w,w0,options)

score=zeros(options.nbgroups,1);

set=find(w);
w2=w(set);

Xt=X.*repmat(Y,1,size(X,2));

err=max(0,(1-Y.*(X(:,set)*w2+w0)));

for  i=1:options.nbgroups
    idx=zeros(size(options.groups));
    idx(options.groups==i)=1;
    idx=logical(idx);
    score(i)=norm(Xt(:,idx)'*err);
end

score=options.C/length(Y)*score-1;

end