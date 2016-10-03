% script to produce single tree extraction for 2013 RSE paper on tree growth

% tile raw data in mat-files
% load /Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/models/dsmdtm.mat
% dx = 1000;
% dy = 1000;
% tx = min(ras.x):dx:max(ras.x);
% ty = min(ras.y):dy:max(ras.y);
% cd /Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/raw
% ras.dtm02(ras.dtm02 < 1500) = NaN;
% m = 0;
% for i = 1:length(tx)
%     for j = 1:length(ty)
%         TX = num2str(tx(i));
%         TY = num2str(ty(j));
%         if 0
%             raw = [];
%             raw.x = [];raw.y = [];raw.z = [];
%             filez = dir('*.mat');
%             for k = 1:length(filez)
%                 disp(filez(k).name);
%                 load(filez(k).name);
%                 ii = x >= tx(i)-10 & x <= tx(i)+dx+10 & y >= ty(j)-10 & y <= ty(j)+dy+10;
%                 if sum(ii) > 0
%                     disp(['Found ',num2str(sum(ii)),' echoes']);
%                     raw.x = [raw.x;x(ii)];
%                     raw.y = [raw.y;y(ii)];
%                     if ~isempty(strfind(filez(k).name,'FBIN'))
%                         raw.z = [raw.z;z(ii)];
%                     else
%                         raw.z = [raw.z;-z(ii)];
%                     end
%                 end
%             end
%             disp(['saving /Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/', ...
%                   TX(1:3),'_',TY(1:3),'_raw.mat']);
%             save(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/', ...
%                   TX(1:3),'_',TY(1:3),'_raw.mat'],'raw');
%         else
%             disp(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/', ...
%                   TX(1:3),'_',TY(1:3),'_raw.mat']);
%             load(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/', ...
%                   TX(1:3),'_',TY(1:3),'_raw.mat']);
%         end
%         if ~exist(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/', ...
%               TX(1:3),'_',TY(1:3),'_clust.mat']);
%             ii = ras.x >= tx(i) & ras.x <= tx(i)+dx; 
%             jj = ras.y >= ty(j) & ras.y <= ty(j)+dy;
%             res.x = floor(min(ras.x(ii))):0.5:ceil(max(ras.x(ii)));
%             res.y = floor(min(ras.y(jj))):0.5:ceil(max(ras.y(jj)));
%             if length(raw.x) > 20
%                 disp('Interpolating DSM ...')
%                 II = raw.z > 0;
%                 raw.et = raw.z;
%                 raw.et(II) = 1;raw.et(~II) = 2;
%                 dsm = raw2ras([raw.x(II),raw.y(II),raw.z(II)],0.5,0.5,'dsm');
%                 [dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
%                 dsm.z = inpaint_nans(dsm.z,4);
%                 chm = dsm;
%                 chm.z = dsm.z-interp2(ras.X(jj,ii),ras.Y(jj,ii),ras.dtm02(jj,ii), ...
%                                       dsm.X,dsm.Y);
%                 chm.z(chm.z > 50) = NaN;
%                 raw.oz = raw.z;
%                 raw.z = raw.z - interp2(ras.X,ras.Y,ras.dtm02,raw.x,raw.y);
%                 ii = ~isnan(raw.z);
%                 raw = subsetraw(raw,ii);
%                 tree = locmax(chm.z,3,0.65,2);
%                 loc(:,1) = chm.X(tree);
%                 loc(:,2) = chm.Y(tree);
%                 loc(:,3) = chm.z(tree);
%                 if length(loc(:,1) > 1)
%                     [x,y,z,idx,raw] = segtree(raw,loc,2,4);
%                     trees = crdata(x,y,z);
%                     save(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/', ...
%                           TX(1:3),'_',TY(1:3),'_clust.mat'],'raw','trees','x','y','z','idx','chm','dsm');
%                 end
%                 clear chm loc 
%            end
%         end
%     end
% end
filez = dir('/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/*_clust.mat');
for i = 1:length(filez)
    disp(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/',filez(i).name]);
    load(['/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/',filez(i).name]);
    if i == 1 
        TREES = concattrees([],trees);
    else
        TREES = concattrees(TREES,trees);
    end
end
trees = TREES;
save('/Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/trees.mat','trees');
