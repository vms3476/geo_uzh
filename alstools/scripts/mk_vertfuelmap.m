% script for producing map of vertical fuel stratification in lamanon
cc
NUMS = [9,17,18,25,26,27,33,34,35,36];
load /Users/morsdorf/Paradox/PARADOX_FIELDALS/dtm.mat
load /Users/morsdorf/Paradox/PARADOX_FIELDALS/lam_fielddata.mat
load /Users/morsdorf/Paradox/PARADOX_FIELDALS/lam_hi-res.mat
load /Users/morsdorf/Paradox/PARADOX_FIELDALS/shrubcleared.mat
dtm = DTM;clear DTM;
%dx = 666000;dy = 4842000;
%dtm.x = dtm.x - dx;
%dtm.y = dtm.y - dy;
%raw.x = raw.x - dx;
%raw.y = raw.y - dy;
[dtm.X, dtm.Y] = meshgrid(dtm.x,dtm.y);
dtm.z = inpaint_nans(dtm.z,4);
oraw = raw;
raw.z =  raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
raw.int( raw.int > 125 ) = 125;
ij = [1:4,1];
oak.ox = xp([1 3 2 12]);
oak.oy = yp([1 3 2 12]);
if 1
ii = double(raw.rnnr) == 9 | double(raw.rnnr) == 17 | double(raw.rnnr) == 18 ;
rasint = raw2ras([raw.x(ii),raw.y(ii),raw.int(ii)],0.5,0.5,'int');
%myscatter3(raw.x(ii),raw.y(ii),raw.z(ii),raw.int(ii),gray(256));
rasint.int = inpaint_nans(rasint.int,4);
[rasint.X,rasint.Y] = meshgrid(rasint.x,rasint.y);
dtm.int = interp2(rasint.X,rasint.Y,rasint.int,dtm.X,dtm.Y);
imagesc(dtm.x,dtm.y,dtm.int);
axis equal;axis tight;axis xy;
hold on

%raw.int = log(raw.int);
th =  190;
hold on
%hp(1) = plot3(mixed.ox(ij),mixed.oy(ij),ones(1,5)*th,'--k','linewidth',2);
%hp(2) = plot3(pine.ox(ij),pine.oy(ij),ones(1,5)*th,'-k','linewidth',2);
%plot3(control.ox(ij),control.oy(ij),ones(1,5)*5,'-k','linewidth',2);
%hp(3) = plot3(oak.ox(ij),oak.oy(ij),ones(1,5)*th,':k','linewidth',2);
%hp(4) = plot3(control.ox(ij),control.oy(ij),ones(1,5)*th,'-.k','linewidth',2);

hp(1) = plot(mixed.ox(ij),mixed.oy(ij),'-r','linewidth',1.5);
hp(2) = plot(pine.ox(ij),pine.oy(ij),'-c','linewidth',1.5);
%plot3(control.ox(ij),control.oy(ij),ones(1,5)*5,'-k','linewidth',2);
hp(3) = plot(oak.ox(ij),oak.oy(ij),'-m','linewidth',1.5);
hp(4) = plot(control.ox(ij),control.oy(ij),'-g','linewidth',1.5);
hp(5) = plot(shlc_mixed(ij,1),shlc_mixed(ij,2),'--y','linewidth',1);
hp(5) = plot(shlc_oak(ij,1),shlc_oak(ij,2),'--y','linewidth',1);
hp(5) = plot(shlc_pine(ij,1),shlc_pine(ij,2),'--y','linewidth',1);
view(2);swisstick;
colormap(gray(64))
caxis([30 110]);
colorbar
hl = legend(hp,'Mixed','Pine','Oak','Control','no shrubs','location','southeast');
print -depsc2 /Users/morsdorf/Desktop/int_plots_view_singleground.eps


figure;
orient tall
wysiwyg
str = {'Single','First','Last'};
num = [9,17,18];
alt = 500;
for i = 1:3
  subplot(3,1,i)
  %raw = corr_intlas(raw);
  jj = double(raw.rnnr) == num(i);
  sraw.x = raw.x(jj);sraw.y = raw.y(jj);sraw.z = raw.z(jj);sraw.int = raw.int(jj);
  %if i == 1
  %  sraw.int = sraw.int .* ((alt-sraw.z).^2)./(alt.^2);
  %end
  ii = inpolygon(sraw.x,sraw.y,mixed.ox(ij),mixed.oy(ij));
  hp(1) = plot(sraw.int(ii),sraw.z(ii),'.k');
  hold on
  ii = inpolygon(sraw.x,sraw.y,pine.ox(ij),pine.oy(ij));
  hp(2) = plot(sraw.int(ii),sraw.z(ii),'.c');
  ii = inpolygon(sraw.x,sraw.y,oak.ox(ij),oak.oy(ij));
  hp(3) = plot(sraw.int(ii),sraw.z(ii),'.r');
  xlabel('LIDAR Intensities [DN]')
  ylabel('LIDAR Vegetation Height [m]');
  hold on
  hl = legend(hp,'Mixed','Pine','Oak');
  grid on
  title(str{i});
  axis([0 120 0 16])
end
print -depsc2 /Users/morsdorf/Desktop/int_plots_scat_echos.eps
end

figure
ii = inpolygon(raw.x,raw.y,control.ox(ij),control.oy(ij));
craw = subsetraw(raw,ii);
jj = double(craw.rnnr) == 9 & craw.Classification == 1;
craw = subsetraw(craw,jj);
%ii = double(craw.rnnr) ~= 9 & craw.Classification == 1;
% assign intensity value of closest single echos to multiple echo intensity
%if sum(jj) > 0 & sum(ii) > 0
%  crawS = subsetraw(craw,jj);
%  crawM = subsetraw(craw,ii);
%  len = length(crawM.x);
%  for I = 1:len
%    dis = sqrt( (crawM.x(I) - crawS.x).^2 + ...
%                (crawM.y(I) - crawS.y).^2 + ...
%                ((crawM.z(I) - crawS.z)*10).^2 );
%    [dis,idd] = sort(dis);
%    if length(dis) >= 6
%      crawM.int(I) = median(crawS.int(idd(1:6)));
%    else
%      crawM.int(I) = median(crawS.int(idd(1)));
%    end
%  end
%  craw = subsetraw(crawM,crawS);
%end
xlabel('LIDAR Intensities [DN]')
ylabel('LIDAR Vegetation Height [m]');
hold on
asp = 12;
DX = [craw.z,craw.int/asp];
idx = clusterdata(DX,'maxclust',3,'linkage','ward');
idx(idx==2)=0;
idx(idx==1) = 2;
idx(idx==0) = 1;
grid on
col = 'rgbm';
pns ='.';
box on
for i = 1:max(idx)
  ii = idx == i;
  hp(i) = plot(craw.int(ii),craw.z(ii),[pns,col(i)],'markersize',8);
  COV{i} = cov([craw.int(ii)/10,craw.z(ii)]);
  MU(i,:) = [mean(craw.int(ii))/10,mean(craw.z(ii))];
end
save /Users/morsdorf/Desktop/COV_control.mat COV MU
hl = legend(hp(3:-1:1),'Pine','Oak','Shrub','location','southeast');
mytitle('Control plot');
axis([0 160 -1 18])
setstrings(12)
print -depsc2 /Users/morsdorf/Desktop/int_control_scat_single.eps


figure
ii = inpolygon(raw.x,raw.y,mixed.ox(ij),mixed.oy(ij));
craw = subsetraw(raw,ii);
jj = double(craw.rnnr) == 9 & craw.Classification == 1;
craw = subsetraw(craw,jj);
%ii = double(craw.rnnr) ~= 9 & craw.Classification == 1;
% assign intensity value of closest single echos to multiple echo intensity
%if sum(jj) > 0 & sum(ii) > 0
%  crawS = subsetraw(craw,jj);
%  crawM = subsetraw(craw,ii);
%  len = length(crawM.x);
%  for I = 1:len
%    dis = sqrt( (crawM.x(I) - crawS.x).^2 + ...
%                (crawM.y(I) - crawS.y).^2 + ...
%                ((crawM.z(I) - crawS.z)*10).^2 );
%    [dis,idd] = sort(dis);
%    if length(dis) >= 6
%      crawM.int(I) = median(crawS.int(idd(1:6)));
%    else
%      crawM.int(I) = median(crawS.int(idd(1)));
%    end
%  end
%  craw = subsetraw(crawM,crawS);
%end
xlabel('LIDAR Intensities [DN]')
ylabel('LIDAR Vegetation Height [m]');
hold on
asp = 12;
DX = [craw.z,craw.int/asp];
idx = clusterdata(DX,'maxclust',2,'linkage','ward');
idx(idx==2)=0;
idx(idx==1) = 2;
idx(idx==0) = 1;
grid on
col = 'rgbm';
pns ='.';
box on
for i = 1:max(idx)
  ii = idx == i;
  hp(i) = plot(craw.int(ii),craw.z(ii),[pns,col(i)],'markersize',8);
  COV{i} = cov([craw.int(ii)/10,craw.z(ii)]);
  MU(i,:) = [mean(craw.int(ii))/10,mean(craw.z(ii))];
end
save /Users/morsdorf/Desktop/COV_mixed.mat COV MU
hl = legend(hp(3:-1:1),'Pine','Oak','Shrub','location','southeast');
mytitle('Mixed plot');
axis([0 160 -1 18])
setstrings(12)
print -depsc2 /Users/morsdorf/Desktop/int_mixed_scat_single.eps


% figure
% subplot(2,1,1)
% jj = double(raw.Classification) == 2;
% sraw.x = oraw.x(jj);sraw.y = oraw.y(jj);sraw.z = oraw.z(jj);sraw.int = oraw.int(jj);
% ras1 = raw2ras([sraw.x,sraw.y,sraw.z],1,1.5,'dtm');
% ras1.z = inpaint_nans(ras1.z,4);
% shadmod(ras1);
% subplot(2,1,2)
% jj = double(raw.Classification) == 2 & ~(double(raw.rnnr) == 9 & raw.int < 80) ;
% sraw.x = oraw.x(jj);sraw.y = oraw.y(jj);sraw.z = oraw.z(jj);sraw.int = oraw.int(jj);
% ras2 = raw2ras([sraw.x,sraw.y,sraw.z],1,1.5,'dtm');
% ras2.z = inpaint_nans(ras2.z,4);
% shadmod(ras2);

figure
orient landscape
wysiwyg
%jj = double(raw.rnnr) == 9; 
%jj = double(raw.Classification) == 1;
%sraw.x = raw.x(jj);sraw.y = raw.y(jj);sraw.z = raw.z(jj);sraw.int = raw.int(jj);
%sraw.class = raw.Classification(jj);
seed = [[30;70;50],[13;6;0.5]];
seed = [[40;80]/asp,[9;4]];
seed = [[40;80;50],[12;6;1]];
%seed = [[30;80;55],[14;8;2]];
%seed = [14;8;1];
%raw = corr_intlas(raw);
[lay] = fuellayers2(raw,seed,5,11.5);               
k = 0;
str = {'Pine','Oak','Shrub'};
for i = 1:3
    k = k + 1;
    den = lay(i).den;
    lay(i).den(den < 0.05) = NaN;
    lay(i).men(den < 0.05) = NaN;
    lay(i).ext(den < 0.05) = NaN;
    
    thick = lay(i).ext;
    lay(i).den(thick < 0.2) = NaN;
    lay(i).men(thick < 0.2) = NaN;
    lay(i).ext(thick < 0.2) = NaN;
    thick = lay(i).ext;

    
    subplot(3,3,k)
    height = lay(i).men;
    imagesc(lay(i).x,lay(i).y,height);
    axis xy; axis equal;axis tight; swisstick
    hold on
    plot(mixed.ox(ij),mixed.oy(ij),'-w','linewidth',2);
    plot(pine.ox(ij),pine.oy(ij),'-w','linewidth',2);
    plot(oak.ox(ij),oak.oy(ij),'-w','linewidth',2);
    plot(control.ox(ij),control.oy(ij),'-w','linewidth',2);
    hp(1) = plot(mixed.ox(ij),mixed.oy(ij),'-or','linewidth',1);
    hp(2) = plot(pine.ox(ij),pine.oy(ij),'-oc','linewidth',1);
    hp(3) = plot(oak.ox(ij),oak.oy(ij),'-om','linewidth',1);
    hp(4) = plot(control.ox(ij),control.oy(ij),'-og','linewidth',1);
    caxis([0 16]);
    colormap(myspecmap);
    if i == 3
      hcb = colorbar('east');
      set(hcb,'xcolor','w','ycolor','w');
      pos = get(hcb,'position');
      set(hcb,'position',[pos(1)+0.015 pos(2)-0.05 pos(3)-0.015 pos(4)+0.1])
      set(hcb,'ytick',[0 2 4 6 8 10 12 14 16]);
      set(get(hcb,'ylabel'),'string','[ m ]');
    end
    title([str{i},' - layer height']);
    k = k + 1; 
    subplot(3,3,k)
    imagesc(lay(i).x,lay(i).y,thick);
    axis xy; axis equal;axis tight; swisstick
    hold on
    plot(mixed.ox(ij),mixed.oy(ij),'-w','linewidth',2);
    plot(pine.ox(ij),pine.oy(ij),'-w','linewidth',2);
    plot(oak.ox(ij),oak.oy(ij),'-w','linewidth',2);
    plot(control.ox(ij),control.oy(ij),'-w','linewidth',2);
    hp(1) = plot(mixed.ox(ij),mixed.oy(ij),'-or','linewidth',1);
    hp(2) = plot(pine.ox(ij),pine.oy(ij),'-oc','linewidth',1);
    hp(3) = plot(oak.ox(ij),oak.oy(ij),'-om','linewidth',1);
    hp(4) = plot(control.ox(ij),control.oy(ij),'-og','linewidth',1);
    caxis([0 8]);
    if k == 5
      hl = legend(hp,'Mixed','Pine','Oak','Control','location','southeast');
      set(hl,'color','k');
      ht = findobj(hl,'type','text');
      set(ht,'color','w');
    end
    colormap(myspecmap)
    cmap = stretch('tan');
    colormap(cmap);
    if i == 3
      hcb = colorbar('east');
      set(hcb,'xcolor','w','ycolor','w');
      pos = get(hcb,'position');
      set(hcb,'position',[pos(1)+0.015 pos(2)-0.05 pos(3)-0.015 pos(4)+0.1])
      set(hcb,'ytick',[0 2 4 6 8]);
      set(get(hcb,'ylabel'),'string','[ m ]');
    end
    title([str{i},' - layer extent']);
    k = k + 1;
    
    subplot(3,3,k)
    imagesc(lay(i).x,lay(i).y,lay(i).den*100);
    axis xy; axis equal;axis tight; swisstick
    hold on
    plot(mixed.ox(ij),mixed.oy(ij),'-w','linewidth',2);
    plot(pine.ox(ij),pine.oy(ij),'-w','linewidth',2);
    plot(oak.ox(ij),oak.oy(ij),'-w','linewidth',2);
    plot(control.ox(ij),control.oy(ij),'-w','linewidth',2);
    hp(1) = plot(mixed.ox(ij),mixed.oy(ij),'-or','linewidth',1);
    hp(2) = plot(pine.ox(ij),pine.oy(ij),'-oc','linewidth',1);
    hp(3) = plot(oak.ox(ij),oak.oy(ij),'-om','linewidth',1);
    hp(4) = plot(control.ox(ij),control.oy(ij),'-og','linewidth',1);
    caxis([0 1]*100);
    if i == 3
      hcb = colorbar('east');
      set(hcb,'xcolor','w','ycolor','w');
      pos = get(hcb,'position');
      set(hcb,'position',[pos(1)+0.015 pos(2)-0.05 pos(3)-0.015 pos(4)+0.1])
      set(hcb,'ytick',[0 0.2 0.4 0.6 0.8 1]*100);
      set(get(hcb,'ylabel'),'string','[ % ]');
    end
    title([str{i},' - layer cover']);
end
set(gcf,'position',[341 64 2000 900],'papersize',[20 9], ...
        'paperposition',[0 0 20 9],'inverthardcopy','off');
set(gcf,'color','w')
are = (max(raw.x)-min(raw.x)) * (max(raw.y)-min(raw.y));
jj = double(raw.rnnr) == 9;
numbsigechos = sum(jj)/are;

jj = double(raw.rnnr) == 9 & raw.Classification == 1;
numbsigechos2 = sum(jj)/are;

print -depsc2 /Users/morsdorf/Desktop/maps_pine_oak_shrub_height.eps
return

figure
jj = double(raw.rnnr) == 18 & raw.Classification == 1;
sraw.x = raw.x(jj);sraw.y = raw.y(jj);sraw.z = raw.z(jj);sraw.int = raw.int(jj);
fint = raw2ras([sraw.x,sraw.y,sraw.int],5,5,'int');
jj = double(raw.rnnr) == 18 & raw.Classification == 2;
sraw.x = raw.x(jj);sraw.y = raw.y(jj);sraw.z = raw.z(jj);sraw.int = raw.int(jj);
lint = raw2ras([sraw.x,sraw.y,sraw.int],5,5,'int');
lint.int = lint.int(1:50,:);
lint.int(isnan(lint.int))=0;fint.int(isnan(fint.int))=0;
plot(fint.int(:),lint.int(:),'.')
[dum,dum,ax] = add_regress(fint.int(:),lint.int(:));
axes(ax(1));
axis([0 70 0 70]);
hold on
xlabel('First echo intensity - Vegetation [DN]')
ylabel('First echo intensity - Ground [DN]')
plot([0,70],[0,70],'-k','linewidth',2);
axes(ax(2))
print -depsc2 /Users/morsdorf/Desktop/int_groundveg_first_scat.eps

figure
jj = double(raw.rnnr) == 9 & raw.Classification == 1;
sraw.x = raw.x(jj);sraw.y = raw.y(jj);sraw.z = raw.z(jj);sraw.int = raw.int(jj);
fint = raw2ras([sraw.x,sraw.y,sraw.int],5,5,'int');
jj = double(raw.rnnr) == 9 & raw.Classification == 2;
sraw.x = raw.x(jj);sraw.y = raw.y(jj);sraw.z = raw.z(jj);sraw.int = raw.int(jj);
lint = raw2ras([sraw.x,sraw.y,sraw.int],5,5,'int');
lint.int = lint.int(1:50,:);
lint.int(isnan(lint.int))=0;fint.int(isnan(fint.int))=0;
plot(fint.int(:),lint.int(:),'.')
[dum,dum,ax] = add_regress(fint.int(:),lint.int(:));
axes(ax(1));
axis([0 70 0 70]);
hold on
xlabel('Single echo intensity - Vegetation [DN]')
ylabel('Single echo intensity - Ground [DN]')
plot([0,70],[0,70],'-k','linewidth',2);
axes(ax(2))
print -depsc2 /Users/morsdorf/Desktop/int_groundveg_single_scat.eps


figure
orient tall
wysiwyg
subplot(3,1,1)
jj = find(double(raw.rnnr) == 17);
ii = double(raw.rnnr(jj+1)) ~= 18;
jj(ii) = [];
fz = raw.z(jj)-raw.z(jj+1);
bins = min([fz]):0.5:max([fz]);
[fz,height] = hist(fz,bins);
cmap = flipud(gray(12));
colormap(cmap([2],:));
hb = bar(height,[fz]',1,'grouped');
set(hb,'edgecolor','k')
hl = legend(hb(1),['first to second']);
hlt = findobj(hl,'type','text');
mytitle('Distance between echos : two returns');
ylabel('Number of laser shots');
set(gca,'xlim',[0 17]);
set(hlt,'fontsize',12,'fontweight','bold');
xlabel('Distance [m]');
grid on
%delete(hl)

subplot(3,1,2)
jj = find(double(raw.rnnr) == 25);
ii = double(raw.rnnr(jj+1)) ~= 26 | double(raw.rnnr(jj+2)) ~= 27;
jj(ii) = [];
fz = raw.z(jj)-raw.z(jj+1);
lz = raw.z(jj+1)-raw.z(jj+2);
mz = raw.z(jj)-raw.z(jj+2);
[fz,height] = hist(fz,bins);
[lz,height] = hist(lz,bins);
[mz,height] = hist(mz,bins);
colormap(cmap([2,4,7],:));
hb = bar(height,[fz;lz;mz]',1,'grouped');
set(hb,'edgecolor','k')
hl = legend(hb(1:3),['first to second'],['second to third'],['first to third']);
xlabel('Distance [m]');
mytitle('Distance between echos : three returns');
ylabel('Number of laser shots');
set(gca,'xlim',[0 17]);
hlt = findobj(hl,'type','text');
set(hlt,'fontsize',12,'fontweight','bold');
grid on
%delete(hl)


subplot(3,1,3)
jj = find(double(raw.rnnr) == 33);
ii = double(raw.rnnr(jj+1)) ~= 34 | double(raw.rnnr(jj+2)) ~= 35 | double(raw.rnnr(jj+3)) ~= 36;
jj(ii) = [];
z1 = raw.z(jj)-raw.z(jj+1);
z2 = raw.z(jj+1)-raw.z(jj+2);
z3 = raw.z(jj+2)-raw.z(jj+3);
z4 = raw.z(jj)-raw.z(jj+3);
[z1,height] = hist(z1,bins);
[z2,height] = hist(z2,bins);
[z3,height] = hist(z3,bins);
[z4,height] = hist(z4,bins);
colormap(cmap([2,4,7,10],:));
hb = bar(height,[z1;z2;z3;z4]',1,'grouped');
set(hb,'edgecolor','k')
hl = legend(hb(1:4),['first to second'],['second to third'],['third to fourth'], ...
            ['first to fourth']);
xlabel('Distance [m]');
mytitle('Distance between echos : four returns');
ylabel('Number of laser shots');
set(gca,'xlim',[0 17]);
hlt = findobj(hl,'type','text');
set(hlt,'fontsize',12,'fontweight','bold');
grid on
setstrings(12)
print -depsc2 /Users/morsdorf/Desktop/hist_return_diffs.eps