function [] = visplot3(trees,dtm,val);
% function [] = visplot3(trees,dtm,val);
% Visualizes Forest Stand using simple tree models  

x = trees.x(:);y = trees.y(:);th = trees.h(:);
ch = trees.h(:)-trees.cbh(:);cd = trees.dia(:);

if nargin == 1
    for i = 1:length(x)
        treemodel3(double(x(i)),double(y(i)),ones(size(y(i)))*0, ...
                   double(th(i)),double(ch(i)),double(cd(i)));
    end
elseif nargin > 1
    [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
    z = interp2(dtm.X,dtm.Y,dtm.z,x,y);
    if nargin == 3
        col = colormap;
        numclass = length(col);
        maxs = nanmax(val);
        mins = nanmin(val);
        ii = floor( (val - mins ) * (numclass-1) / (maxs - mins) );
        ii = ii + 1;
        for i = 1:length(x)
            if isnan(ii(i))
                treemodel3(double(x(i)),double(y(i)),double(z(i)), ...
                           double(th(i)),double(ch(i)),double(cd(i)),1,[0.1 0.5 0.1]);
            else
                treemodel3(double(x(i)),double(y(i)),double(z(i)), ...
                           double(th(i)),double(ch(i)),double(cd(i)),1,col(ii(i),:));
            end
        end
    else
        for i = 1:length(x)
            treemodel3(double(x(i)),double(y(i)),double(z(i)), ...
                               double(th(i)),double(ch(i)),double(cd(i)),1);
        end
    end
    shadmod(dtm);colormap(gray);
    view(3);axis equal;axis tight;swisstick
end