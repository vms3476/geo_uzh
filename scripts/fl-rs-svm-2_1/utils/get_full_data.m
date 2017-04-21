function [X,Y,tasks]=get_full_data(Xc,Yc)

X=[];
Y=[];
tasks=[];


for i=1:length(Xc)
   X=[X;Xc{i}];
   if nargin>1
    Y=[Y;Yc{i}]; 
    tasks=[tasks;i*ones(size(Yc{i}))]; 
   end
end
