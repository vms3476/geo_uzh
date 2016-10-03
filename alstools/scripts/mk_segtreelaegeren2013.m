% script for testing inverse CHM for starting point locations of tree segmentation
close all
%load /Users/morsdorf/data/Laegeren/laegeren_all.mat
are = [min(raw.x) max(raw.x) min(raw.y) max(raw.y)];
S = shaperead('/Users/morsdorf/data/Laegeren/LAEG_fieldmap_2013/LAEG_fieldmap_2013.shp');
%[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
%RAW = getrawlas('/Users/morsdorf/data/Laegeren/RAW_SPRING',are);
%RAW.tz = RAW.z - interp2(dtm.X,dtm.Y,dtm.z,RAW.x,RAW.y);
% maxh = 35;
% figure(1);clf;
% ii = RAW.tz > 3 & RAW.tz < maxh;
% int = RAW.int;
% int(int > 7000) = 7000;
% clf; myscatter3(RAW.x(ii),RAW.y(ii),RAW.tz(ii),RAW.tz(ii),flipud(summer(256)));
% hold on
% plot3([S.X;S.X],[S.Y;S.Y],[ones(size([S.X]))*0;ones(size([S.X]))*maxh],'-k','linewidth',4);
% plot3([S.X],[S.Y],[ones(size([S.X]))*0],'*r','markersize',10);
% plot3([S.X],[S.Y],[ones(size([S.X]))*maxh],'*r','markersize',10);
% mat3d2osg('pointcloud_lageren_upsidedown3_stems');
% figure(2);clf
% res = 0.5;
% dh = 0:res:max(RAW.tz);
% ax = [min(RAW.x) max(RAW.x) min(RAW.y) max(RAW.y)];
% for i = 1:length(dh);
%     ii = RAW.tz > dh(i) & RAW.tz < dh(i) + res;
    
%     invchm = raw2ras([RAW.x(ii),RAW.y(ii),-RAW.tz(ii)],0.5,0.5,'dsm');
%     if i == 1
%         vox = ones([size(invchm.num),length(dh)])*NaN;
%     end
%     vox(:,:,i) = invchm.num;
%     %invchm.z = inpaint_nans(invchm.z,4);
%     myimage(invchm.x,invchm.y,invchm.num);colorbar;colormap(gray);
%     caxis([0 2]);
%     title(['Number of echoes @ ',num2str(dh(i)+res/2),' [m] height above ground']);
%     axis(ax)
%     drawnow
%     print('-djpeg100',['ani/scliceforest',int2strv(i,1,'0',3),'.jpg']);
% end
% [X,Y,Z] = meshgrid(invchm.x,invchm.y,dh);
% return
figure

%chm = raw2ras([RAW.x,RAW.y,RAW.tz],0.5,0.5,'dsm');
%chm.z = inpaint_nans(chm.z,4);
[ii,chms] = locmax(chm.z,6,1.5,3);
myimage(chm.x,chm.y,chms);
hold on
[chm.X,chm.Y] = meshgrid(chm.x,chm.y);
plot(chm.X(ii),chm.Y(ii),'.w');
plot([S.X_crown],[S.Y_crown],'k.');
SZ = interp2(chm.X,chm.Y,chm.z,[S.X_crown],[S.Y_crown]);
raw.x = RAW.x;raw.y = RAW.y;raw.z = RAW.tz;
loc = [[S.X_crown];[S.Y_crown];SZ];
ii = [S.SOZ] == 1;
loc = loc(:,ii);
ii = ~isnan(loc(3,:));
loc = loc(:,ii);
minheight = 3;
asp = 5;
method = 'gmmclust';

trees = segtree(raw,loc',minheight,asp,method);
%[x,y,z] = gmmclust(raw,loc',minheight,asp);
%trees = crdata(x,y,z);
return
ftree.x = [S.X];
ftree.y = [S.Y];
ftree.h = [S.BHD];
ltree.x = chm.X(ii);
ltree.y = chm.Y(ii);
ltree.h = chm.z(ii);
[trees,stats] = matchtrees(ftree,ltree);