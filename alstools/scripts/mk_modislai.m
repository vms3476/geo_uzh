clear all;close all
cd /Users/morsdorf/Desktop/MODIS_LAI/TerraAcqua4Days
if 0
  filez = dir('*.hdf');
  for i = 1:length(filez);
    disp(['!gdal_translate -sds ',filez(i).name,' lai.tif']);
    eval(['!gdal_translate -sds ',filez(i).name,' lai.tif']);
    disp(['!gdalwarp -t_srs EPSG:2056 lai.tif2 lai2.tif'])
    eval(['!gdalwarp -t_srs EPSG:2056 lai.tif2 lai2.tif'])
    disp(['!gdalwarp -t_srs EPSG:2056 lai.tif3 qc2.tif'])
    eval(['!gdalwarp -t_srs EPSG:2056 lai.tif3 qc2.tif'])
    [lai.dat{i},A] = geotiffread('lai2.tif');
    [lai.qc{i},A] = geotiffread('qc2.tif');
    [lai.x{i},lai.y{i}] = pixcenters(A,A.RasterSize(1),A.RasterSize(2));
  end
  save modis_lai.mat lai

% convert LAI quality flags

load modis_lai
for i = 1:length(lai.qc);
    qc{i} = arrayfun(@(x) sum(double(bitget(x,6:8)).*[1,10,100]),lai.qc{i});
    cc{i} = arrayfun(@(x) sum(double(bitget(x,4:5)).*[1,10]),lai.qc{i});
end
save modis_lai2
end
load modis_lai2
are = [2555000 2625000 1210000 1245000];
for i = 1:length(lai.dat);
  subplot(5,2,i);
  lai.dat{i}(lai.dat{i} > 249) = NaN;
  lai.dat{i}(qc{i} > 0) = NaN;
  ii = lai.x{i} >= are(1) & lai.x{i} <= are(2);
  jj = lai.y{i} >= are(3) & lai.y{i} <= are(4);
  myimage(lai.x{i}(ii),lai.y{i}(jj),lai.dat{i}(jj,ii)*0.1);
  caxis([0 4])
  colormap([1 1 1;vegetation]);
  colorbar
end 