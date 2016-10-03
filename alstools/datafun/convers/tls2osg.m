% convert tls .asc data to osg
dat = load('/Volumes/Archive/TLSParadox06/Scans/lamanon_dessampl_cent.xyz.asc');

% subsampling
ii = 5;
dat = dat(1:ii:end,:);
dist = sqrt( dat(:,1).^2 + dat(:,2).^2 + dat(:,3).^2);
% calibration of tls intensity
res = 0.2;
idist = [0:res:50];
for i = 1:length(idist)-1;
  ii = dist >= idist(i) & dist < idist(i+1);
  int(i) = median(dat(ii,4));
end
idist = idist(1:end-1)+res;
ii = isnan(int);
idist(ii) = [];
int(ii) = [];
P = polyfit(idist,int,10);
dat(:,4) = dat(:,4)-interp1(idist,int,dist); 
myscatter3(dat(:,1),dat(:,2),dat(:,3),dat(:,4),gray(256));

mat3d2osg('/Users/morsdorf/Desktop/lamanon_dessampl_cent');