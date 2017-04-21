function res=Omega12_group(w,group)

res=0;
vals=unique(group);
for i=1:length(vals)
    res=res+norm(w(group==vals(i))); 
end