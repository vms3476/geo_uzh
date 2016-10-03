% script for producing the voxel grid on the laegeren site
 cc
 cd /Users/morsdorf/Projects/Active/3DVegLab/TreeRekVoxelGrid
 if 1
   load /Users/morsdorf/Projects/Active/3DVegLab/TreeRekVoxelGrid/pointc_fin.mat
   dtm.x = pc_fin.mat_x;
   dtm.y = pc_fin.mat_y;
   dtm.z = pc_fin.dtm_mat;
   sum = pc_fin.season == 2; 
   spr = pc_fin.season == 1; 
   strtit = {'leaf-on','leaf-off'};
   strprt = {'summer','spring'};
   for i = 1:2
     if i == 1
       ii = sum;
     elseif i == 2;
       ii = spr;
     end
     raw.x = pc_fin.x_point(ii); 
     raw.y = pc_fin.y_point(ii); 
     raw.z = pc_fin.z_point(ii); 
     raw.int = pc_fin.Int(ii); 
     raw.wdt = pc_fin.Width(ii); 
     raw.rnnr = pc_fin.EchoType(ii);
     raw.ScanAngle = pc_fin.ScanAngle(ii);
     raw.GPSTime = pc_fin.GPStime(ii);
     raw.Classification = pc_fin.FligthStrip(ii);
     vox{i} = als2vox(raw,1);
     raw.z = pc_fin.z_point_ag(ii); 
     veg{i} = fcolai(raw,dtm,1,3);
     
     [X,Y,Z] = meshgrid(vox{i}.x,vox{i}.y,vox{i}.z);
     INT = vox{i}.ras(:,:,:,2);
     WDT = vox{i}.ras(:,:,:,3);
     WDT(WDT< 4) = 4;
     WDT(WDT> 10) = 10;
     INT(INT > 25000) = 25000;
     NUM = vox{i}.ras(:,:,:,1);
     ii = NUM > 1;

     figure(i);clf;orient landscape;wysiwyg;
     subplot(1,2,1)
     myscatter3(X(ii),Y(ii),Z(ii),INT(ii),vegetation2(128));
     colorbar('South');
     swisstick
     title(['Voxelgrid Intensity - 1 m Resolution - ',strtit{i}]);
     subplot(1,2,2)
     myscatter3(X(ii),Y(ii),Z(ii),WDT(ii),vegetation2(128));
     colorbar('South');
     swisstick
     title(['Voxelgrid Echo Width - 1 m Resolution - ',strtit{i}]);
     eval(['print -depsc2 -r300 voxelgrid_laegeren_',strprt{i},'.mat']);
   end
   save /Users/morsdorf/Projects/Active/3DVegLab/TreeRekVoxelGrid/voxelgrid.mat
 end
 load /Users/morsdorf/Projects/Active/3DVegLab/TreeRekVoxelGrid/voxelgrid.mat
 figure(3)
 subplot(1,2,1)
 myimage(veg{1}.x,veg{1}.y,veg{1}.fco);
 swisstick;colormap(vegetation2);colorbar;
 title('Fractional Cover - leaf-on');
 subplot(1,2,2)
 myimage(veg{2}.x,veg{2}.y,veg{2}.fco);
 swisstick;colormap(vegetation2);colorbar;
 title('Fractional Cover - leaf-off');
 print -depsc2 fcover_laegeren.eps
 figure(4)
 subplot(1,2,1)
 myimage(veg{1}.x,veg{1}.y,veg{1}.lai);
 swisstick;colormap(vegetation2);colorbar;
 title('Leaf Area Index - leaf-on');
 subplot(1,2,2)
 myimage(veg{2}.x,veg{2}.y,veg{2}.lai);
 swisstick;colormap(vegetation2);colorbar;
 title('Leaf Area Index - leaf-off');
 print -depsc2 LAI_laegeren.eps