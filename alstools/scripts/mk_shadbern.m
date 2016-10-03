% script for producing shade comparisons of DTM and DOM of Bern's jura ALS data.
close all
if 0
    cc
    load /Volumes/Data2/BernerJura/DOMDTM.mat
    load /Volumes/Data2/BernerJura/Referenz/domdtmav.mat
    dom.x = dom.x - 2000000;dtm.x = dtm.x - 2000000;
    dom.y = dom.y - 1000000;dtm.y = dtm.y - 1000000;
    dom.z = inpaint_nans(dom.z,4);
    dtm.z = inpaint_nans(dtm.z,4);
    %domav = xyz2ras(['/Users/morsdorf/Desktop/BernerJura/DOMAV/DOMAV/ascii/' ...
    %                 '112514_dom_03.xyz']);
    %dtmav = xyz2ras(['/Users/morsdorf/Desktop/BernerJura/DOMAV/DTMAV/ascii/' ...
    %                 '112514_dtm_03.xyz']); 
    %save /Users/morsdorf/Desktop/BernerJura/DOMDTMs.mat

%load /Users/morsdorf/Desktop/BernerJura/DOMDTMs.mat

figure(1)
orient tall
wysiwyg
subplot(2,1,1)
shadmod(dom);
swisstick;colorbar;caxis([665.6 1225])
title('DOM - ALS 60, 2011');
subplot(2,1,2)
shadmod(domav);
swisstick;colorbar;caxis([665.6 1225])
title('DOM - AV, Swisstopo');
colormap(terrain)
%print -depsc2 -r600 /Volumes/Data2/BernerJura/doms.eps
figure(2)
orient tall
wysiwyg
subplot(2,1,1)
shadmod(dtm);
swisstick;colorbar;caxis([665.6 1225])
title('DTM - ALS 60, 2011');
subplot(2,1,2)
shadmod(dtmav);
swisstick;colorbar;caxis([665.6 1225])
title('DTM - AV, Swisstopo');
colormap(terrain)
%print -depsc2 -r600 /Volumes/Data2/BernerJura/dtms.eps
end
%define ROIs
XC = [ 572245;
       573219;
       574557;
       573681;
       575939];
YC = [ 225332;
       226227;
       224985;
       225737;
       226367];
for i = 1:length(XC)
   ARE{i} = [XC(i)-350 XC(i)+350 YC(i)-250 YC(i)+250]; 
end
% figure(3);
% shadmod(dom);
% swisstick;colorbar;caxis([665.6 1225])
% title('DOM - ALS 60, 2011'); hold on
% for i = 1:length(ARE)
%     hp(i) = plot3(ARE{i}([1,2,2,1,1]),ARE{i}([3,3,4,4,3]),ones(1,5)*1200);
%     hp2(i) = plot3(ARE{i}([1,2,2,1,1]),ARE{i}([3,3,4,4,3]),ones(1,5)*1200);
%     text(ARE{i}(1)+20,ARE{i}(3)+10,1200,['ROI ',num2str(i)],'color','w','fontweight','bold', ...
%          'fontsize',12,'horizontalal','left','verticalalign','bottom');
% end
% set(hp,'color','k','linewidth',2);set(hp2,'color','w','linewidth',2,'linestyle', ...
%                                           '--');
% colormap(terrain)
% print -depsc2 -r600 /Volumes/Data2/BernerJura/dsmROIs.eps
% return
[dtmav.X,dtmav.Y] = meshgrid(dtmav.x,dtmav.y);
[domav.X,domav.Y] = meshgrid(domav.x,domav.y);
for i = 1:length(ARE);
    % %DTM
    ii = find(dtmav.x >= ARE{i}(1) & dtmav.x <= ARE{i}(2));
    jj = find(dtmav.y >= ARE{i}(3) & dtmav.y <= ARE{i}(4));
   
   rdtmav.z = dtmav.z(jj,ii);rdtmav.x = dtmav.x(ii);rdtmav.y = dtmav.y(jj);
   ii = find(dtm.x >= ARE{i}(1) & dtm.x <= ARE{i}(2));
   jj = find(dtm.y >= ARE{i}(3) & dtm.y <= ARE{i}(4));

   rdtm.z = dtm.z(jj,ii);rdtm.x = dtm.x(ii);rdtm.y = dtm.y(jj);
   [rdtmav.X,rdtmav.Y] = meshgrid(rdtmav.x,rdtmav.y);
   [rdtm.X,rdtm.Y] = meshgrid(rdtm.x,rdtm.y);
   figure(3+i)
   orient tall
   wysiwyg
   subplot(3,1,1);
   shadmod(rdtm);swisstick;
   colormap(ocean2);colorbar;title(['ROI ',num2str(i),' - Flotron DTM']);
   subplot(3,1,2);
   shadmod(rdtmav);swisstick
   colormap(ocean2);colorbar;title(['ROI ',num2str(i),' - Swisstopo DTM-AV']);
   rdtm.avz = interp2(dtmav.X,dtmav.Y,dtmav.z,rdtm.X,rdtm.Y,'linear');
   subplot(3,1,3)
   rng = 20;
   S{i} = sysshift(rdtm.z,rdtm.avz);
   S{i} = S{i} * median(diff(dtm.x));
   %[m,n] = size(S{i})
   %dx = [1:n]-(rng+1);dy = [1:m]-(rng+1);[DX,DY]=meshgrid(dx,dy);
   %results(i) = autoGaussianSurf(DX,DY,S{i});
   %shx = results(i).x0;shy = results(i).y0;
   rdtm.avz = interp2(dtmav.X,dtmav.Y,dtmav.z,rdtm.X,rdtm.Y,'linear');
   diffz = rdtm.z-rdtm.avz;
   diffz_nonans = diffz(:);
   diffz_nonans(isnan(diffz_nonans))=[];
   shz(i) = median(diffz_nonans);
   diffz = diffz-shz(i);
   SS{i} = sysshift(rdtm.z,rdtm.avz)*median(diff(dtm.x));
   imagesc(rdtm.x,rdtm.y,diffz);axis xy;axis equal;axis tight;swisstick; ...
       colorbar
   cax = caxis;maxc = max(abs(cax));caxis([-maxc*0.75 maxc*0.75]);
   caxis([-0.5 0.5]);
   title(['xshift = ',num2str(S{i}(1)),', yshift =',num2str(S{i}(2)),', z-shift =', ...
          num2str(shz(i))])
      title(['Hoehendifferenz DTM 2011 - DTM-AV, mittl. Hoehenabweichung (korrigiert) =', ...
          num2str(shz(i))])
   eval(['print -depsc2 -r600 /Volumes/Data2/BernerJura/dtmROI_',num2str(i), ...
         '.eps'])
   
   
   % DOM
   % ii = find(domav.x >= ARE{i}(1) & domav.x <= ARE{i}(2));
   % jj = find(domav.y >= ARE{i}(3) & domav.y <= ARE{i}(4));
   
   % rdomav.z = domav.z(jj,ii);rdomav.x = domav.x(ii);rdomav.y = domav.y(jj);
   % ii = find(dom.x >= ARE{i}(1) & dom.x <= ARE{i}(2));
   % jj = find(dom.y >= ARE{i}(3) & dom.y <= ARE{i}(4));

   % rdom.z = dom.z(jj,ii);rdom.x = dom.x(ii);rdom.y = dom.y(jj);
   % [rdomav.X,rdomav.Y] = meshgrid(rdomav.x,rdomav.y);
   % [rdom.X,rdom.Y] = meshgrid(rdom.x,rdom.y);
   % figure(3+i)
   % orient tall
   % wysiwyg
   % subplot(3,1,1);
   % shadmod(rdom);swisstick;
   % colormap(ocean2);colorbar;title(['ROI ',num2str(i),' - Flotron DOM']);
   % subplot(3,1,2);
   % shadmod(rdomav);swisstick
   % colormap(ocean2);colorbar;title(['ROI ',num2str(i),' - Swisstopo DOM-AV']);
   % rdom.avz = interp2(domav.X,domav.Y,domav.z,rdom.X,rdom.Y,'linear');
   % subplot(3,1,3)
   % rng = 20;
   % S{i} = sysshift(rdom.z,rdom.avz);
   % S{i} = S{i} * median(diff(dom.x));
   % %[m,n] = size(S{i})
   % %dx = [1:n]-(rng+1);dy = [1:m]-(rng+1);[DX,DY]=meshgrid(dx,dy);
   % %results(i) = autoGaussianSurf(DX,DY,S{i});
   % %shx = results(i).x0;shy = results(i).y0;
   % rdom.avz = interp2(domav.X,domav.Y,domav.z,rdom.X,rdom.Y,'linear');
   % diffz = rdom.z-rdom.avz;
   % diffz_nonans = diffz(:);
   % diffz_nonans(isnan(diffz_nonans))=[];
   % shz(i) = median(diffz_nonans);
   % diffz = diffz-shz(i);
   % SS{i} = sysshift(rdom.z,rdom.avz)*median(diff(dom.x));
   % imagesc(rdom.x,rdom.y,diffz);axis xy;axis equal;axis tight;swisstick; ...
   %     colorbar
   % cax = caxis;maxc = max(abs(cax));caxis([-maxc*0.75 maxc*0.75]);
   % title(['xshift = ',num2str(S{i}(1)),', yshift =',num2str(S{i}(2)),', z-shift =', ...
   %        num2str(shz(i))])
   %    title(['Hoehendifferenz DOM 2011 - DOM-AV, mittl. Hoehenabweichung (korrigiert) =', ...
   %        num2str(shz(i))])
   % eval(['print -depsc2 -r600 /Volumes/Data2/BernerJura/domROI_',num2str(i),'.eps'])
end

% compute systematic shifts of models



%[dtmav.X,dtmav.Y] = meshgrid(dtmav.x,dtmav.y);
%[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
%dtmav.alsz = interp2(dtm.X,dtm.Y,dtm.z,dtmav.X,dtmav.Y);
%dtmav,alsz = inpaint_nans(dtmav.alsz,4);
%S = sysshift(dtmav.z,dtmav.alsz);