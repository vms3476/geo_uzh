% script for converting

cd /Volumes/WINDOWSDISK/ParadoxFIELDCPGN/SCANS/ASC

fnames = dir('*.asc');

for i = 1:length(fnames)
  if ~exist([fnames(i).name(1:end-8),'.mat'])
    disp(['Processing : ',fnames(i).name]);
    try
      dat = load(fnames(i).name);
      save(fnames(i).name(1:end-8),'dat');
      clear dat
    catch
      disp('ASCII File not clean !')
    end
  end
end
