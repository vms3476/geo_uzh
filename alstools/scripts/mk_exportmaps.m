% script for exporting the maps of bernese jura to geotiff for nicholas
% Berner Jura
if 0
  cd /Users/morsdorf/Projects/Archived/BernJura
  %load /Users/morsdorf/Projects/Archived/BernJura/vegetation_res3.mat
  TER(TER<200) = NaN;
  str = {'FCO','LAI','TER','VEH'};
  strl = {'fcover','lai','dtm','vegetationheight'};
  for i = 1:4
    eval(['dat = ',str{i},';'])
    wrt_geotiff_CH(strl{i},x,y,dat);
  end
end
% Swiss National Park
cd /Users/morsdorf/Teaching/MSc/MA_SW/RASTER2010
load dtmdsm2010
chm = dsm;
dtm.z(dtm.z<1614) = NaN;
chm.z = dsm.z - dtm.z;
chm.z(chm.z<0) = 0;
chm.z(chm.z>40) = NaN; 
str = {'chm','dtm'};
strl = {'chm_snp','dtm_snp'};
for i = 1:2
  eval(['dat = ',str{i},'.z;'])
  eval(['x = ',str{i},'.x;'])
  eval(['y = ',str{i},'.y;'])
  wrt_geotiff_CH(strl{i},x,y,dat);
end
