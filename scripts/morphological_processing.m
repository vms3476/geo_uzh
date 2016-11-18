% test morphological processing 


% remove noisy pixels in coniferous regions
pix1 = 9;
conn1 = 4;
ras1mClassMorph = bwareaopen(ras1mClass,pix1,conn1);
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorph); colormap gray
title('Class map after bwareaopen')

% remove noisy pixels in deciduous regions
pix1 = 9;
conn1 = 4;
ras1mClassMorphInv = bwareaopen(1-ras1mClassMorph,pix1,conn1); 
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorphInv); colormap gray
title('Class map after two bwareaopen operations')

% conect coniferous regions. coniferous = 0, deciduous = 1
se = strel('disk',3);
ras1mClassMorphClose = 1-imclose(ras1mClassMorphInv, se);
figure; myimage(ras1mDenTh2.x,ras1mDenTh2.y,ras1mClassMorphClose); colormap gray
title('Class map after closing')

overlay_polygons(laegernTreeTable_final);
%% ground class
% initial ground class image, ground = 2, vegetation = 3
ras1mGrnd = raw2ras([x;y;classification]',1,1,'dsm');

%ras1mGrndTile669258 = ras1mGrnd;



figure; myimage(ras1mGrnd.x,ras1mGrnd.y,ras1mGrnd.z); colormap gray
title('initial ground classification')

% make binary image
ras1mGrnd.z(ras1mGrnd.z == 2) = 1; 
ras1mGrnd.z(ras1mGrnd.z == 3) = 0;
grnd1mMorph = bwareaopen(ras1mGrnd.z,36,4); 
figure; myimage(ras1mGrnd.x,ras1mGrnd.y,grnd1mMorph); colormap gray
title('ground class map after 1 bwareaopen operation')

se = strel('disk',5);
grnd1mMorphClose = imclose(grnd1mMorph,se);
figure; myimage(ras1mGrnd.x,ras1mGrnd.y,grnd1mMorphClose); colormap gray
title('ground class map after closing operation')

grnd1mMorphTotal = grnd1mMorphClose;

class1mMap = ras1mClassMorphClose;
class1mMap(ras1mClassMorphClose == 0) = 2; % coniferous class value of 2
class1mMap(grnd1mMorphTotal == 1) = 0; % ground class value of 0
figure; myimage(ras1mGrnd.x,ras1mGrnd.y,class1mMap); colormap gray
title('ground class map after 1 bwareaopen operation')


overlay_polygons(laegernTreeTable_final);