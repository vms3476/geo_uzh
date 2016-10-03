function [hp] = visclu(x,y,z); 
% function [hp] = visclu(x,y,z);
% function to visualize the output of segtree
hold on
k = 0;
for i = 1:length(x)
    if length(x{i}) > 0
        k = k + 1;
        hp(k) = plot3(x{i},y{i},z{i},'.','color',rand([1,3]),'markersize',5); 
    end
end
axis equal
axis tight
swisstick
view(3)