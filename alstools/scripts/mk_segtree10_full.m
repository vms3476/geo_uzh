% script to produce single tree extraction for 2013 RSE paper on tree growth, 2010 data
cc
% tile raw data in mat-files
load /Volumes/Data1/Ofenpass2010/segtree2013/dtmdsm2010.mat
clear X XI Y YI
dx = 1000;
dy = 1000;
tx = min(dsm.x):dx:max(dsm.x);
ty = min(dsm.y):dy:max(dsm.y);
dtm.z(dtm.z < 1500) = NaN;
dsm.z(dsm.z < 1500) = NaN;
if 0
    cd /Volumes/Data1/Ofenpass2010/Punktwolke/Original/
    m = 0;
    filez = dir('*.las');
    for k = 1:length(filez)
        disp(filez(k).name);
        RAW = readlas(filez(k).name);
        x = RAW.x;y = RAW.y;z = RAW.z;clear RAW;
        for i = 1:length(tx)
            for j = 1:length(ty)
                TX = num2str(tx(i));
                TY = num2str(ty(j));
                %are =  [tx(i)-10 tx(i)+dx+10 ty(j)-10 ty(j)+dy+10];
                %raw = getrawlas('/Volumes/Data1/Ofenpass2010/Punktwolke',are);
                ii = x >= tx(i)-10 & x <= tx(i)+dx+10 & y >= ty(j)-10 & y <= ty(j)+dy+10;
                if sum(ii) > 0
                    disp(['Found ',num2str(sum(ii)),' echoes']);
                    if ~exist(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                               TX(1:3),'_',TY(1:3),'_raw.mat']);
                        raw.x = x(ii);raw.y = y(ii);raw.z = z(ii);
                    else
                        load(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                          TX(1:3),'_',TY(1:3),'_raw.mat']);
                        raw.x = [raw.x;x(ii)];
                        raw.y = [raw.y;y(ii)];
                        raw.z = [raw.z;z(ii)];
                    end
                    disp(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                          TX(1:3),'_',TY(1:3),'_raw.mat']);
                    save(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                          TX(1:3),'_',TY(1:3),'_raw.mat'],'raw');
                    clear raw; 
                end
            end
        end
    end
end
for i = 1:length(tx)
    for j = 1:length(ty)
        TX = num2str(tx(i));
        TY = num2str(ty(j));
        if exist(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                  TX(1:3),'_',TY(1:3),'_raw.mat']);
            disp(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                  TX(1:3),'_',TY(1:3),'_raw.mat']);
            load(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                  TX(1:3),'_',TY(1:3),'_raw.mat'],'raw');
            RAW.x = raw.x;RAW.y = raw.y;RAW.z = raw.z;raw = RAW;clear RAW;
            if ~exist(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                       TX(1:3),'_',TY(1:3),'_clust.mat']);
                ii = dsm.x >= tx(i) & dsm.x <= tx(i)+dx; 
                jj = dsm.y >= ty(j) & dsm.y <= ty(j)+dy;
                chm.z = dsm.z(jj,ii)-dtm.z(jj,ii);
                chm.z(chm.z > 50) = NaN;
                chm.x = dsm.x(ii);
                chm.y = dsm.y(jj);
                ii = dtm.x >= tx(i)-10 & dtm.x <= tx(i)+dx+10; 
                jj = dtm.y >= ty(j)-10 & dtm.y <= ty(j)+dy+10;
                DTM.x = dtm.x(ii);
                DTM.y = dtm.y(jj);
                DTM.z = dtm.z(jj,ii);
                [DTM.X,DTM.Y] = meshgrid(DTM.x,DTM.y);
                raw.oz = raw.z;
                [chm.X,chm.Y] = meshgrid(chm.x,chm.y);
                ii = raw.x >= tx(i)-10 & raw.x <= tx(i)+dx+10 & ...
                     raw.y >= ty(j)-10 & raw.y <= ty(j)+dy+10;
                raw = subsetraw(raw,ii);
                raw.z = raw.z - interp2(DTM.X,DTM.Y,DTM.z,raw.x,raw.y);
                ii = ~isnan(raw.z) & raw.z < 50;
                raw = subsetraw(raw,ii);
                if length(raw.x) > 20
                    tree = locmax(chm.z,5,0.65,2);
                    if sum(tree(:)) > 0 
                        loc(:,1) = chm.X(tree);
                        loc(:,2) = chm.Y(tree);
                        loc(:,3) = chm.z(tree);
                        [x,y,z,idx,raw] = segtree(raw,loc,2,4,'gmm');
                        trees = crdata(x,y,z);
                        save(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/', ...
                              TX(1:3),'_',TY(1:3),'_clust.mat'],'raw','trees','x','y','z','idx','chm');
                        clear chm loc 
                    end
                end               
            end    
        end
    end
end
filez = dir('/Volumes/Data1/Ofenpass2010/segtree2013/tiles/*_clust.mat');
for i = 1:length(filez)
    disp(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/',filez(i).name]);
    load(['/Volumes/Data1/Ofenpass2010/segtree2013/tiles/',filez(i).name]);
    if i == 1 
        TREES = concattrees([],trees);
    else
        TREES = concattrees(TREES,trees);
    end
end
trees = TREES;
save('/Volumes/Data1/Ofenpass2010/segtree2013/trees.mat','trees');
