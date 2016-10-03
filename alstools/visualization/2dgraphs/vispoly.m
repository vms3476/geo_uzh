function [hp] = vispoly(S); 
% function [hp] = vispoly(S);
% function to visualize the polygons of tree outlines loaded from shapefile
% Input argument S is structure provided by shaperead

% Felix Morsdorf, Sep. 2014

hold on
k = 0;
for i = 1:length(S)
    if length(S(i).X) > 0
        k = k + 1;
        hp(k) = plot(S(i).X,S(i).Y,'-','color',rand([1,3]),'linewidth',1); 
    end
end
axis equal
axis tight
swisstick
