% make fcover image for schweizer forst zeitschrift
clf;
load /Volumes/Data2/mat-files/fco_big_5m.mat
load /Volumes/Data2/ofenpass/gesamtgebiet/model/dtm/dtm.mat
[X,Y] = meshgrid(x_fcov2,y_fcov2);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
FCO = interp2(X,Y,fcov2,dtm.X,dtm.Y);
shadmod(dtm.x,dtm.y,dtm.z,FCO*100);
swisstick
grid off
colormap([1 1 1;jet])
set(gca,'visible','off')
caxis([0 70]);
%hcb = colorbar('EastOutside');
%axes(hcb);ylabel('Deckungsgrad [%]')