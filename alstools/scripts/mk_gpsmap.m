clf;
cd /Users/morsdorf/Desktop/Laegeren/
S = shaperead('PntPosn.shp');
load /Users/morsdorf/Desktop/Laegeren/dtm/dtm.mat
xbb = [min([S.X]) max([S.X])];
ybb = [min([S.Y]) max([S.Y])];
dxy = 100;
ii = dtm.x >= xbb(1)-dxy & dtm.x <= xbb(2)+dxy;
jj = dtm.y >= ybb(1)-dxy & dtm.y <= ybb(2)+dxy;
dtm.x = dtm.x(ii);
dtm.y = dtm.y(jj);
dtm.z = dtm.z(jj,ii);
shadmod(dtm);
hold on
%for i = 1:length(S)
%    [X,Y] = circle(S(i).H_Prec_Obs, S(i).X, S(i).Y);
%    plot3(X,Y,ones(size(X))*S(i).Elev_Obs + 50,'-k','linewidth',1.1);
%    plot3(X,Y,ones(size(X))*S(i).Elev_Obs + 50,'--w','linewidth',1.1);
%end

plot3([S.X], [S.Y], [S.Elev_Obs] + 50,'ow','markersize',1);
plot3([S.X], [S.Y], [S.Elev_Obs] + 50,'xk','markersize',1);
plname = {S.Point_Name};
for i = 1:length(plname)
    plname{i} = plname{i}(5:end);
end
ht = text([S.X], [S.Y]+5, [S.Elev_Obs] + 50, {S.Point_Name}, ...
          'interpreter','none','fontweight','bold','fontsize',2);
set(ht,'horizontalalign','center');
swisstick
title('GPS Points Laegeren TSL Campaign 2010 - July 13/14th');
%