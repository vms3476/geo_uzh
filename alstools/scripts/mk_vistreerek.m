% make nice animation of georek with point cloud of Lageren
% raw data half/half
cc
load /Volumes/Data1/Laegeren/Summer_2010/summer_raw.mat
rawsum = raw;
load /Volumes/Data1/Laegeren/Spring_2010/spring_raw.mat
rawspr = raw;
%are = [min(raw.x) max(raw.x) min(raw.y) max(raw.y)];

%are = [669747 669762 258985 259115];
are = [669770 669830 258965 259025];
%are = [672236-30 672236+30 258986-30 258986+30];
ii1 = rawsum.x >= are(1) & rawsum.x <= are(2) & rawsum.y >= are(3) & rawsum.y <= ...
      are(4);
ii2 = rawspr.x >= are(1) & rawspr.x <= are(2) & rawspr.y >= are(3) & rawspr.y <= ...
      are(4);

ii3 = raw.x >= are(1) & raw.x <= are(2) & raw.y >= are(3) & raw.y <= ...
      are(4);
xoff = 669700;
yoff = 259000;
raw1.x = rawsum.x(ii1)-xoff;
raw1.y = rawsum.y(ii1)-yoff;
raw1.z = rawsum.z(ii1);

raw2.x = rawspr.x(ii2)-xoff;
raw2.y = rawspr.y(ii2)-yoff;
raw2.z = rawspr.z(ii2);

raw.x = [raw1.x;raw2.x];
raw.y = [raw1.y;raw2.y];
raw.z = [raw1.z;raw2.z];
% load models
load /Volumes/Data1/Laegeren/Summer_2010/dsm_fe/dsm.mat
load /Volumes/Data1/Laegeren/Spring_2010/dtm/dtm.mat

ii = dtm.x >= are(1) & dtm.x <= are(2); 
jj = dtm.y >= are(3) & dtm.y <= are(4);

dtm.x = double(dtm.x(ii))-xoff;
dtm.y = double(dtm.y(jj))-yoff;
dtm.z = double(dtm.z(jj,ii));

% convert raw data to height over terrain
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);

[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
tz = interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
dsm.z = double(interp2(dsm.X-xoff,dsm.Y-yoff,dsm.z,dtm.X,dtm.Y));
oz = raw.z;
raw.z = raw.z - tz;
% convert dsm to chm 
chm = dtm;
chm.z = dsm.z-dtm.z;
% call cluster analysis
[clx,cly,clz,idx,RAW] = segtree2(chm,raw,2,2);

% get tree data
[trees] = crdata(clx,cly,clz);

clf;
if 0
  jj = trees.y >  mean(dtm.y);
  X = trees.x(~jj);  Y = trees.y(~jj);
  H = trees.h(~jj);
  trees.x = trees.x(jj);
  trees.y = trees.y(jj);
  trees.h = trees.h(jj);
  trees.cbh = trees.cbh(jj);
  trees.vol = trees.vol(jj);
  trees.dia = trees.dia(jj)*2;
  trees.points = trees.points(jj);
  trees.olx = trees.olx(jj);
  trees.oly = trees.oly(jj);
end

visplot2(trees,dtm);

hold on
if 0 
  ii = raw.y < mean(are(3:4));
else
  ii = 1:length(raw.y);
end
ii = raw.z < 0.5;
idx(ii) = 0;
myscatter3(raw.x,raw.y,oz,idx,myspecmap,2);
Z = interp2(dtm.X,dtm.Y,dtm.z,X,Y);
for i = 1:length(X)
  plot3([X(i),X(i)],[Y(i),Y(i)],[Z(i),Z(i)+H(i)],'-ok', ...
        'linewidth',1);
end
caxis([min(oz(ii)) max(oz(ii))])
box on
set(gca,'visible','on')
xlabel('Distanz [m]');
ylabel('Distanz [m]');
%hb = get(gca,'xlabel');
%delete(hb);
%hb = get(gca,'ylabel')
%delete(hb);
zlabel('Hoehe [m]');
set(gca,'xtick',[50 60])
rotate3d on
view(107,16);
orient landscape
wysiwyg
%axis vis3d
axis equal
colormap(hsv);
hxl = get(gca,'xlabel');
hyl = get(gca,'ylabel');
posx = get(hxl,'position');
posy = get(hyl,'position');

set(hxl,'position',[posx(1)-20 posx(2)-19 posx(3)-2])
set(hyl,'position',[posy(1)-10 posy(2)-30 posy(3)+5])

%hcb = incolorbar(min(oz(ii)):2:max(oz(ii)),2,'Height a.s.l. (m)','k');
%set(hcb,'position',[0.18 0.35 0.35 0.03]);
%set(hcb,'xtick',[floor(min(oz(ii))):2:ceil(max(oz(ii)))])
hcb = colorbar('south');
pos = get(hcb,'position');
set(hcb,'position',[pos(1)+0.29 0.27 pos(3)-0.3 0.03]);
axes(hcb);
set(hcb,'xaxislocation','top');
xlabel('Hoehe [m]');
setstrings(11)
delete(hcb)
set(gcf,'renderer','painters');
orient portrait
%print -depsc2 -r450 /Users/morsdorf/Desktop/Abbildungen_SFZ/georek_laegeren_v4.eps



