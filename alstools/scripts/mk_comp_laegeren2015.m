% script for comparing datasets of laegeren 2015
cc
cd /Users/morsdorf/Desktop/Laegeren2014
% %raw = getrawlas('/Users/morsdorf/data/Laegeren/2014/03_LAS_LV03_LN02',[669660 669960 258910 259210]);


% raw780_off = getrawlas('/Users/morsdorf/Desktop/Laegeren2014/ZH_Off',[669660 669960 258910 259210]);
% raw680_on = getrawlas('/Users/morsdorf/Desktop/Laegeren2014/AG_On',[669660 669960 258910 259210]);
% raw680_off = getrawlas('/Users/morsdorf/Desktop/Laegeren2014/AG_Off',[669660 669960 258910 259210]);

% save /Users/morsdorf/Desktop/Laegeren2014/comp2014.mat

load /Users/morsdorf/Desktop/Laegeren2014/comp2014.mat
ii = raw780_off.z < 900;
raw780_off = subsetraw(raw780_off,ii);

load /Users/morsdorf/data/Laegeren/2010/dtm.mat;

[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
raw680_off.tz = raw680_off.z;
raw680_off.z = raw680_off.tz - interp2(dtm.X,dtm.Y,dtm.z,raw680_off.x,raw680_off.y);

raw780_off.tz = raw780_off.z;
raw780_off.z = raw780_off.tz - interp2(dtm.X,dtm.Y,dtm.z,raw780_off.x,raw780_off.y);

for i = 1
    vox_1550 = als2vox(raw680_off,i);
    vox_1064 = als2vox(raw780_off,i);
    
    [X,Y,Z] = meshgrid(vox_1550.x,vox_1550.y,vox_1550.z);
    [X4,Y4,Z4] = meshgrid(vox_1064.x,vox_1064.y,vox_1064.z);
    
    R = vox_1064.ras(:,:,:,2);R(R>700)=700;
    R = R/max(R(:));
    NIR = vox_1550.ras(:,:,:,2);NIR(NIR>300)=300;
    NIR = NIR/max(NIR(:));
    NNIR = interp3(X,Y,Z,NIR,X4,Y4,Z4);
    NDVI = (NNIR-R)./(NNIR+R);
    

    clf;myscatter3(X(:),Y(:),Z(:),NDVI(:),parula(256),10);
    mat3d2osg(['Laegeren_NDVI_ALS_RES_VEGH',num2str(i)]);
end


return
if 1
    %raw = getrawlas('/Users/morsdorf/Desktop/Laegeren_KantonZH2015',[669660 669960 258910 259210]);
    %raw = readlas('/Users/morsdorf/Desktop/Laegeren_KantonZH2015/data.las');
%[A1, refmat1, bbox2] = geotiffread('/Users/morsdorf/Desktop/Laegeren_KantonZH2015/6690_2580.tif');
%[A2, refmat2, bbox2] = geotiffread('/Users/morsdorf/Desktop/Laegeren_KantonZH2015/6690_2590.tif');

load /Users/morsdorf/data/Laegeren/2010/dtm.mat;

[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);

raw.ch = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);

ii = raw.z >= 600 & raw.z <= 880 & raw.ch >= -2 & raw.ch <= 60;

clf;myscatter3(raw.x(ii)-669660,raw.y(ii)-258910,raw.z(ii),raw.ch(ii));

mat3d2osg('Laegeren_ZH_SWSPHT.osg');
end
raw.int(raw.int > 500) = 500;

clf;myscatter3(raw.x(ii)-669660,raw.y(ii)-258910,raw.z(ii),raw.int(ii),ndvimap(256));

mat3d2osg('Laegeren_ZH_SWSPHT_int.osg');