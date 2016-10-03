cd /Volumes/Data1/Belgium
dirstr = {'aelmoeseneiebos','kersselaerspleyn','kluisbos','wijnendalebos'};
if 1 
    % for j = 1:length(dirstr);
    %     cd([dirstr{j},'/raw']);
    %     filez = dir('*.las');
    %     for i = 1:length(filez)
    %         disp(filez(i).name);
    %         raw = readlas(filez(i).name);
    %         names = fieldnames(raw);
    %         if i == 1
    %             for k = 1:length(names)
    %                 eval(['nraw.',names{k},' = [];'])
    %             end
    %         end
    %         for k = 1:length(names)
    %             eval(['nraw.',names{k},' = [nraw.',names{k},';raw.',names{k},'(:)];'])
    %         end
    %     end
    %     clear raw;raw = nraw;
    %     cd ..
    %         eval(['save -v7.3 ',dirstr{j},'_raw.mat raw']);
    %     cd ..
    % end
    for j = 1:length(dirstr);
        cd([dirstr{j}]);
        load([dirstr{j},'_raw.mat']);
        load([dirstr{j},'.mat']);
        %numlay = 11;
        res = 1;
        [xg,yg,numechos,numshots] = pointdens(dtm,raw,res);
        res = 2;
        %[xg,yg,multih,multid,ras] = als2multi(dtm,raw,res,numlay);
        %save multiband_10.mat xg yg multih multid ras
        save numshots.mat xg yg numechos numshots
        figure(j);clf;
        myimage(xg,yg,numechos);swisstick;
        print -depsc2 numechos.eps
        figure(j);clf;
        myimage(xg,yg,numshots);swisstick;
        print -depsc2 numshots.eps
        cd ..
     end
end
return
for j = 1:length(dirstr);
    cd([dirstr{j}]);
    load multiband_10.mat 
    pxszx = median(diff(xg));
    pxszy = median(diff(yg));
    ul_x = min(xg);
    ul_y = max(yg);
    hdrstr = {['map info = {UTM, 1.000, 1.000,', ...
               num2str(ul_x),',', ...
               num2str(ul_y),',',num2str(pxszx),',',num2str(pxszy), ...
               ', Belgium, IGN, 72, units=Meters}'], ...
              'wavelength units = Unknown', ...
              ['pixel size = {',num2str(pxszx),',',num2str(pxszy),', units=Meters}']};
    [m,n,p] = size(multid);
    nmultid = ones(n,m,p);
    for j = 1:p
        nmultid(:,:,j) = fliplr(squeeze(multid(:,:,j))');
    end
    enviwrite(nmultid,'multi_density',hdrstr);
    nmultih = ones(n,m,p);
    for j = 1:p
        nmultih(:,:,j) = fliplr(squeeze(multih(:,:,j))');
    end
    enviwrite(nmultih,'multi_height',hdrstr);
    cd ..
 end