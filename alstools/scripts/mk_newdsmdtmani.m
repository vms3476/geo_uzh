cc
load ~/Publications/FCOVER_DIFF/raw.mat
load ~/Publications/FCOVER_DIFF/dtm.mat
ii = rawr.x >= 624905 & rawr.x <= 624935 & rawr.y >= 96070 & rawr.y <= 96100;
raw = subsetraw(rawr,ii);
ii = dtm.x >= 624905-5 & dtm.x <= 624935+5;
jj = dtm.y >= 96070-5 & dtm.y <= 96100+5;
dtm.x = dtm.x(ii);dtm.y = dtm.y(jj);dtm.z = dtm.z(jj,ii);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
myscatter3(raw.x,raw.y,raw.z,raw.z);
hold on
x = linspace(min(raw.x),max(raw.x),30);
y = linspace(min(raw.y),max(raw.y),30);
[X,Y] = meshgrid(x,y);
cd ~/Graphiken/dsmdtmani
str = {'dsm','dtm','min'};
for l = 1:3
    Z = ones(size(X))*min(raw.z);
    swisstick
    azs = -45;
    els =  10;
    aze = 45;
    ele = 30;
    view(azs,els)
    k = 1;
    print('-djpeg99','-r300',[upper(str{l}),'/',str{l},'_ani',int2strv(k,1,'0',4),'.jpg']);
    h = mesh(X,Y,Z);
    set(h,'edgecolor','k','facecolor','w');
    k = 2;
    print('-djpeg99','-r300',[upper(str{l}),'/',str{l},'_ani',int2strv(k,1,'0',4),'.jpg']);
    delete(h);
    len = (length(x)-1)*(length(y)-1);
    az = linspace(azs,aze,len);
    el = linspace(els,ele,len);
    axis vis3d
    for i = 1:length(x)-1
        for j = 1:length(y)-1
            ii = raw.x >= x(i) & raw.x <= x(i+1) & raw.y >= y(j) & raw.y <= y(j+1);
            k = k + 1;
            if l == 1
                Z(j,i) = max(raw.z(ii));
            elseif l == 2
                Z(j,i) = interp2(dtm.X,dtm.Y,dtm.z,x(i),y(j));
            elseif l == 3
                Z(j,i) = min(raw.z(ii));
            end
            h = mesh(X,Y,Z);
            set(h,'edgecolor','k','facecolor','w');
            view(az(k-2),el(k-2))
            drawnow;
            print('-djpeg99','-r300',[upper(str{l}),'/',str{l},'_ani',int2strv(k,1,'0',4),'.jpg']);
            delete(h);
        end
    end
end