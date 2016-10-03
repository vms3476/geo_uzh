cc
load /Volumes/Data1/Laegeren/Summer_2010/summer_raw.mat
rawsum = raw;
load /Volumes/Data1/Laegeren/Spring_2010/spring_raw.mat
rawspr = raw;
raw.x = raw.x - 2000000;
raw.y = raw.y - 1000000;
%are = [min(raw.x) max(raw.x) min(raw.y) max(raw.y)];
are = [669850 670000 259050 259200];
are = [669980 670000 259000 259300];
are = [669770 669790 258975 259115];% flux tower for SFZ
are = [669774 669790 258975 259032.5];% Trees and flux tower for ESA deliverable
%are = [672236-50 672236+50 258986-30 258986+30];
ii1 = rawsum.x >= are(1) & rawsum.x <= are(2) & rawsum.y >= are(3) & rawsum.y <= ...
      are(4);
ii2 = rawspr.x >= are(1) & rawspr.x <= are(2) & rawspr.y >= are(3) & rawspr.y <= ...
      are(4);

raw1.x = rawsum.x(ii1)-669700;
raw1.y = rawsum.y(ii1)-259000;
raw1.z = rawsum.z(ii1);

raw2.x = rawspr.x(ii2)-669700;
raw2.y = rawspr.y(ii2)-259000;
raw2.z = rawspr.z(ii2);

hcb = pldiscret(raw1,raw2);
%hca = get(gcf,'children');
%delete(hca(2));
swisstick('german')
hxl = get(gca,'xlabel');
hyl = get(gca,'ylabel');
posx = get(hxl,'position');
posy = get(hyl,'position');

set(hxl,'position',[posx(1)-20 posx(2)-19 posx(3)-2])
set(hyl,'position',[posy(1)-10 posy(2)-30 posy(3)+5])
zlabel('Hoehe [m]')

view(107,16);
orient landscape
wysiwyg
pos = get(hcb,'position');
set(hcb,'position',[pos(1)+0.29 0.2 pos(3)-0.35 0.03]);
%set(hcb,'position',[pos(1)*1.1 pos(2)*0.9 pos(3) pos(4)*0.8]);
%set(gcf,'renderer','painter');
%print -depsc -r300 /Users/morsdorf/Desktop/pointcloud_forpatrick.eps