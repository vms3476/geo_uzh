cc
if 0
    cd /Users/morsdorf/Desktop/Laegeren/laegern_tls
    filez = dir('*.las');
    col = hsv(6);
    x = [];
    y = x; z = x; c = x;
    for i = 1:length(filez)
        raw = readlas(filez(i).name);
        x = [x;raw.x];    y = [y;raw.y];    z = [z;raw.z];    c = [c;ones(size(raw.z))*i];
    end
save /Users/morsdorf/Desktop/Laegeren/laegern_tls/allscans.mat
end
load /Users/morsdorf/Desktop/Laegeren/laegern_tls/allscans.mat
clf
ii = find(y <= 259020 & y >= 258936 & x >= 669470 & x <= 669480);
minx = min(x(ii));
miny = min(y(ii));
meaz = mean(z(ii));
x = x(ii) - minx;
y = y(ii) - miny;
z = z(ii) - meaz;
dx = 1;
myscatter3(x(1:dx:end),y(1:dx:end),z(1:dx:end),z(1:dx:end),gray(256),1);
hold on
%load /Users/morsdorf/Desktop/Laegeren/laegern_tls/raw_als.mat
load /Volumes/Data1/Laegeren/Summer_2010/summer_raw.mat
%raw_harm = raw;
ii = find(raw.y <= 259020 & raw.y >= 258936 & raw.x >= 669470 & raw.x <= 669480);
X = raw.x(ii) - minx;
Y = raw.y(ii) - miny;
Z = raw.z(ii) - meaz;
myscatter3(X,Y,Z,Z,ocean2(256),7)
mat3d2osg('Laegeren_TLSALS_TRANS_summer');
return
minx = min(x);
miny = min(y);
meaz = mean(z);
meax = mean(x);
meay = mean(y);
set(gca,'visible','off','color','k')
if 0
    for i = 1:2
        set(gca,'projection','perspective')
        if i == 1
            campos([meax-5,meay,meaz])
        else
            campos([meax-5,meay+0.15,meaz])
        end
        camtarget([meax+15,meay,meaz])
        drawnow
        if i == 1
            %        print -depsc2 -r600 left.eps
            print -dtiff -r600 1_left.tif
        else
            %        print -depsc2 -r600 left.eps
            print -dtiff -r600 1_right.tif
        end
    end
else
    dx = -10:0.01:10;
    set(gca,'projection','perspective')
    for i = 1:length(dx)
        campos([meax-dx(i),meay,meaz])
        camtarget([meax-dx(i),meay+20,meaz])
        drawnow
        eval(['print -dtiff -r300 tls_',int2strv(i,1,'0',3),'.tif']);
    end
end
