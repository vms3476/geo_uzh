% script for producing the first image of the belgian surface models
dpath = '/Users/morsdorf/data/Belgium';
str = {'aelmoeseneiebos','kersselaerspleyn','kluisbos','wijnendalebos'};
for i = 1:length(str)
    eval(['cd ',dpath,'/',str{i}])
    filez = dir('saf*.coo');
    dsm = xyz2ras(filez);
    filez = dir('tal*.coo');
    dtm = xyz2ras(filez);
    dtm.z = inpaint_nans(dtm.z);
    chm = dtm;
    if i == 2
        chm.z = dsm.z(1:end-1,:) - dtm.z;
    elseif i == 3
        chm.z = dsm.z(:,2:end) - dtm.z;
    else
        chm.z = dsm.z - dtm.z;
    end
    clf;
    imagesc(chm.x,chm.y,chm.z);axis equal;axis tight;axis xy;swisstick;
    colormap([0,0,0;ocean2]);colorbar
    title(['Canopy height model - ',upper(str{i})]);
    eval(['print -depsc2 ',str{i},'_chm.eps']);
    save(str{i},'dtm','dsm','chm');
end

    