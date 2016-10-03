% script for placing drone of Laegern data on chm and or pointcloud
cc
[img.rgb,R] = geotiffread('/Users/morsdorf/data/laegeren/laegeren_24_10_rgb_orthomosaic_cuted.tif');
[img.x,img.y] = pixcenters(R,R.RasterSize(1),R.RasterSize(2));


[iimg.ind,cmap] = rgb2ind(img.rgb,256);
[iimg.x,iimg.y] = pixcenters(R,R.RasterSize(1),R.RasterSize(2));
[iimg.X,iimg.Y] = meshgrid(iimg.x,iimg.y);
myimage(img.x,img.y,img.rgb);

load /Users/morsdorf/data/laegeren/laegeren_all.mat
clear sum
chm.x = pc_fin.mat_x;
chm.y = pc_fin.mat_y;
[chm.X,chm.Y] = meshgrid(chm.x,chm.y);
chm.z = pc_fin.chm_mat;

dsm.x = pc_fin.mat_x;
dsm.y = pc_fin.mat_y;
[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
dsm.z = pc_fin.dsm_mat;

ii = img.x >= min(dsm.x) & img.x <= max(dsm.x);
jj = img.y >= min(dsm.y) & img.y <= max(dsm.y);
img.rgb = img.rgb(jj,ii,:);
img.x = img.x(ii);
img.y = img.y(jj);
figure;
%shadmod(dsm.x,dsm.y,dsm.z,img.rgb)

ii = pc_fin.season == 2 | pc_fin.season == 1;
raw.x = pc_fin.x_point(ii);
raw.y = pc_fin.y_point(ii);
raw.z = pc_fin.z_point(ii);
raw.h = pc_fin.CHM(ii);
raw.c = interp2(iimg.X,iimg.Y,double(iimg.ind),raw.x,raw.y,'nearest');
%index for dark points
%jj = cmap(:,1) < 0.2 & cmap(:,2) < 0.2 & cmap(:,3) < 0.2; 
%ii = raw.h > 5 & raw.c ==
myscatter3(raw.x-min(raw.x),raw.y-min(raw.y),raw.z,raw.c,cmap,2);
axis vis3d
set(gca,'color','k')
set(gcf,'color','k','inverthardcopy','off');
set(gca,'visible','off','position',[0 0 1 1])
view(2)
[azs,els] = view;
cd /Users/morsdorf/data/Laegeren/ani
view(3)
[aze,ele] = view;
az = linspace(azs,aze,300);
el = linspace(els,ele,300);
view(2)
for i = 1:length(az);
    view(az(i),el(i));
    drawnow
    print('-djpeg99','-r600',['laegeren_RGB_ani',int2strv(i,1,'0',4),'.jpg']);
end