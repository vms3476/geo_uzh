% script for writing the dtm, dsm and chm as asciigrid for hyperforest
cd /Volumes/Data1/Belgium
dirstr = {'aelmoeseneiebos','kersselaerspleyn','kluisbos','wijnendalebos'};
modstr = {'dtm','dsm','chm'};
for i = 1:length(dirstr)
    cd(dirstr{i})
    load(dirstr{i})
    for j = 1:length(modstr)
        %eval(['x = ',modstr{j},'.x;']);
        %eval(['y = ',modstr{j},'.y;']);
        %eval(['ras = ',modstr{j},'.z;']);
        %mat2asciigrid([modstr{j},'.asc'],x,y,ras);
        eval(['geotiff_belgium(modstr{j},',modstr{j},');']);
    end
    cd ..
end
