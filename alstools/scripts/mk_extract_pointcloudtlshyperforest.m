% script for extracting the point clouds for each of the TLS locations in Belgium
% second version - extracts voxel grids as well
cd /Volumes/Data1/Belgium/TLS_location
TLS = shaperead('TLS_coordinates.shp');
dirs = {'kersselaerspleyn','wijnendalebos'};
[idx,dada] = same([TLS.PLOT]);
load(['../',dirs{1},'/',dirs{1},'_raw.mat']);
load(['../',dirs{1},'/',dirs{1},'.mat']);
kers_dtm = dtm;
kers = raw;
load(['../',dirs{2},'/',dirs{2},'_raw.mat']);
load(['../',dirs{2},'/',dirs{2},'.mat']);
wijn_dtm = dtm;
wijn = raw;

for i = 1:length(idx)
    are{i} = [min([TLS(idx{i}).X])-30 max([TLS(idx{i}).X])+30 ...
              min([TLS(idx{i}).Y])-30 max([TLS(idx{i}).Y])+30];
    if strcmp(lower(TLS(idx{i}(1)).SITE(1:4)),dirs{1}(1:4))
        %raw = getrawlas(['../',dirs{1},'/raw'],are{i});
        ii = (kers.x >= are{i}(1) & kers.x <= are{i}(2) &  kers.y >= are{i}(3) & kers.y <= are{i}(4));
        raw = subsetraw(kers,ii);
        dtm = kers_dtm;
        [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
        raw.h = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
    elseif strcmp(lower(TLS(idx{i}(1)).SITE(1:4)),dirs{2}(1:4))
        %raw = getrawlas(['../',dirs{2},'/raw'],are{i});
        ii = (wijn.x >= are{i}(1) & wijn.x <= are{i}(2) &  wijn.y >= are{i}(3) & wijn.y <= are{i}(4));
        raw = subsetraw(wijn,ii);
        dtm = wijn_dtm;
        [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
        raw.h = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
    end
    % compute voxel grid of raw data
    RAW = raw;RAW.z=RAW.h;
    voxres = 1;
    vox = als2vox(RAW,voxres);
    if sum(ii) > 0
        disp(['save ALS_RAW_',upper(TLS(idx{i}(1)).SITE), ...
              '_',num2str(TLS(idx{i}(1)).PLOT),'_res', ...
              num2str(voxres),' raw vox are'])
        eval(['save ALS_RAW_',upper(TLS(idx{i}(1)).SITE), ...
              '_',num2str(TLS(idx{i}(1)).PLOT),'_res', ...
              num2str(voxres),' raw vox are'])
    end
end
