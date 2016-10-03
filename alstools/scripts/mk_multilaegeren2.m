fpath = '/Volumes/Data1/Laegeren/Spring_2010/FDTM';
cd(fpath)
dfilez = dir('*_dtm.mat');
numlay = 11;
maxz = 50;
res = 2;
if 1 
    % for j = 1:length(dfilez);
    %     eval(['load ',dfilez(j).name]);
    %     %dtm = xyz2ras(dfilez(j).name);
    %     %dtm.z = inpaint_nans(double(dtm.z),4);
    %     %dtm.x = double(dtm.x);
    %     %dtm.y = double(dtm.y);
    %     cd ../Pointcloud
    %     are = [min(dtm.x)-5 max(dtm.x)+5 min(dtm.y)-5 max(dtm.y)+5];
    %     lfilez = dir('*.las');
    %     for i = 1:length(lfilez)
    %         disp(lfilez(i).name);
    %         raw = readlas(lfilez(i).name);
    %         names = fieldnames(raw);
    %         if i == 1
    %             for k = 1:length(names)
    %                 eval(['nraw.',names{k},' = [];'])
    %             end
    %         end
    %         ii = (raw.x >= are(1) & raw.x <= are(2) &  raw.y >= are(3) & raw.y <= are(4));
    %         if sum(ii) > 0
    %             raw = subsetraw(raw,ii);
    %             for k = 1:length(names)
    %                 eval(['nraw.',names{k},' = [nraw.',names{k},';raw.',names{k},'(:)];'])
    %             end
    %         end
    %     end
    %     clear raw;raw = nraw;
    %     cd(fpath)
    %     eval(['save -v7.3 ',dfilez(j).name(1:end-8),'_raw.mat raw'])
    %     %eval(['save -v7.3 ',dfilez(j).name(1:end-4),'.mat dtm'])
    % end
    for j = 1:length(dfilez);
        cd(fpath) 
        load([dfilez(j).name(1:end-8),'_raw.mat']);
        load([dfilez(j).name(1:end-8),'_dtm.mat']);
        [xg,yg,multih,multid,ras] = als2multi(dtm,raw,res,numlay,maxz);
        eval(['save ',dfilez(j).name(1:end-8),'_multi.mat xg yg multih multid ras'])
    end
end
cd(fpath)
filez = dir('*_multi.mat');
NX = [];
NY = [];
for i = 1:length(filez)
    eval(['load ',filez(i).name]);
    XG{i} = xg;
    YG{i} = yg;
    NX = [NX,xg];
    NY = [NY,yg];
    MH{i} = multih;
    MD{i} = multid;
    RAS{i} = ras;
end
are = [min(NX) max(NX) min(NY) max(NY)];
XX = min(NX):res:max(NX);
YY = min(NY):res:max(NY);
[XXGG,YYGG] = meshgrid(XX,YY);
omultid = ones(length(YY),length(XX),numlay)*NaN;
omultih = omultid;
oras = ones(length(YY),length(XX),5)*NaN;
for k = 1:length(XG)
    %hw = waitbar(0,'Sorting elements to raster ...');
    ii = XX >= min(XG{k}) & XX <= max(XG{k});
    jj = YY >= min(YG{k}) & YY <= max(YG{k});
    for i = 1:numlay
        omultid(jj,ii,i) = interp2(XG{k},YG{k},MD{k}(:,:,i),XXGG(jj,ii),YYGG(jj,ii));
        omultih(jj,ii,i) = interp2(XG{k},YG{k},MH{k}(:,:,i),XXGG(jj,ii),YYGG(jj,ii));
    end
    rasnames = {'mhome','stdev','meint','stint','numhits'};
    for i = 1:5
        eval(['oras(jj,ii,i) = interp2(XG{k},YG{k},RAS{k}.', ...
              rasnames{i},',XXGG(jj,ii),YYGG(jj,ii));'])
    end
    %[mm,nn] = size(MH{k}(:,:,1));
    %for i = 1:nn
    %  ii = XG{k}(i) == XX;              
    %  for j = 1:mm
    %jj = YG{k}(j) == YY;
    %    if sum(ii) == 1 & sum(jj) == 1
    %     if ~isnan(MH{k}(j,i,:))
    %     omultih(jj,ii,:) = MH{k}(j,i,:);
    %          end
    %      if ~isnan(MD{k}(j,i,:))
    %        omultid(jj,ii,:) = MD{k}(j,i,:);
    %      end
    %    end
    %    if fix(i*j/1000) == (i*j/1000)
    %      waitbar((i*j)/(mm*nn),hw);
    %    end
    %  end
    %end
    %close(hw);
end
% write output to envi format
pxszx = median(diff(XX));
pxszy = median(diff(YY));
ul_x = min(XX);
ul_y = max(YY);
hdrstr = {['map info = {UTM, 1.000, 1.000,', ...
           num2str(ul_x),',', ...
           num2str(ul_y),',',num2str(pxszx),',',num2str(pxszy), ...
           ', Switzerland, CH1903+, units=Meters}'], ...
          'wavelength units = Unknown', ...
          ['pixel size = {',num2str(pxszx),',',num2str(pxszy),', units=Meters}']};
[m,n,p] = size(omultid);
nmultid = ones(n,m,p);
for j = 1:p
    nmultid(:,:,j) = fliplr(squeeze(omultid(:,:,j))');
end
enviwrite(nmultid,'multi_density_laegeren',hdrstr);
nmultih = ones(n,m,p);
for j = 1:p
    nmultih(:,:,j) = fliplr(squeeze(omultih(:,:,j))');
end
enviwrite(nmultih,'multi_height_laegeren',hdrstr);
nras = ones(n,m,5);
for j = 1:5
    nras(:,:,j) = fliplr(squeeze(oras(:,:,j))');
end
enviwrite(nras,'multi_ras_laegeren',hdrstr);

