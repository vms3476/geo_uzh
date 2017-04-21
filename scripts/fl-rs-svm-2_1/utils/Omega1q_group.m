function res=Omega1q_group(w,q,group)

res=0;
vals=unique(group);
for i=1:length(vals)
    res=res+norm(w(group==vals(i)),q); 
end