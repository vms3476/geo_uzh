% script for computing single tree locations for the RootMap3D model

% load data
    cd /Volumes/Data3/20070714_WILER/ALS
    if 0
        raw = readlas('DTM_kptl_082011_LV95.las');
        load dtm_interp/002_outputdata/dtm_filt_fin
        dtm = data;    clear data;
        if 0
            myscatter3(raw.x-min(raw.x),raw.y-min(raw.y),raw.z-min(raw.z),raw.z-min(raw.z));
            mat3d2osg('WILER')
        end
        RAW = [raw.x,raw.y,raw.z];
        dsm = raw2ras(RAW,1,1,'dsm');
        [dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
        [dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
        dtm.z = interp2(dtm.X,dtm.Y,dtm.z,dsm.X,dsm.Y);
        dtm.X = dsm.X;    dtm.Y = dsm.Y;    dtm.x = dsm.x;    dtm.y = dsm.y;
        ii = isnan(dtm.z);
        dsm.z(ii) = 0;
        dsm.z = inpaint_nans(dsm.z,4);
        dsm.z(ii) = NaN;
        save -v7.3 all_data.mat
        load all_data.mat
        chm = dsm;
        chm.z = dsm.z - dtm.z;
        chm.z(chm.z < -5 | chm.z > 50) = NaN;
        chm.z(isnan(chm.z)) = 0;
        raw.z = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
        [x,y,z,idx,RAW] = segtree2(chm,raw,3,6,4);
        trees = crdata(x,y,z);
        save -v7.3 wiler_trees.mat
    end
    load wiler_trees.mat
    ii = trees.h > 3;
    trees = subsetraw(trees,ii);
    ii = trees.h < 60;
    trees = subsetraw(trees,ii);
    ii = trees.dia < 8;
    trees = subsetraw(trees,ii);
    imagesc(chm.x,chm.y,chm.z);axis equal;axis tight;axis xy;
    caxis([0 50])
    hold on
    grid on
    swisstick
    colorbar
    colormap([1 1 1;jet])
    hp  = plot(trees.x,trees.y,'.k','markersize',2);
    set(gcf,'renderer','painters');
    print -depsc2 -r600 trees_wiler.eps