cc
lstr = {'load /data/provence/OPTECH/lamanon/Multi/multimosaic.mat', ...
        'load /data/provence/OPTECH/lamanon_HI-RES/Multi/mb_lam_hi-res.mat', ...
        'load /data/provence/OPTECH/aix/Multi/multimosaic.mat', ...
        'load /data/provence/OPTECH/aix_HI-RES/Multi/multimosaic.mat'};
cdstr = {'cd /data/provence/OPTECH/lamanon/Multi/', ...
        'cd /data/provence/OPTECH/lamanon_HI-RES/Multi/', ...
        'cd /data/provence/OPTECH/aix/Multi/', ...
        'cd /data/provence/OPTECH/aix_HI-RES/Multi/'};

for i = 1%:length(lstr)
  eval(lstr{i})
  eval(cdstr{i})
  if i == 2
    xg = x;yg = y;omultid = multid;omultih = multih;onumhits = numhits;
  end
  pxszx = median(diff(xg));
  pxszy = median(diff(yg));
  ul_x = min(xg);
  ul_y = max(yg);
  hdrstr = {['map info = {UTM, 1.000, 1.000,', ...
             num2str(ul_x),',', ...
             num2str(ul_y),',',num2str(pxszx),',',num2str(pxszy), ...
             ', 31, North, WGS-84, units=Meters}'], ...
            'wavelength units = Unknown', ...
          ['pixel size = {',num2str(pxszx),',',num2str(pxszy),', units=Meters}']};
  [m,n,p] = size(omultid);
  if i == 3
    save omultih.mat omultih 
    save onumhits.mat onumhits
    clear omultih onumhits;
  end
  nmultid = ones(n,m,p);
  for j = 1:p
    nmultid(:,:,j) = fliplr(squeeze(omultid(:,:,j))');
    %nmultid(:,:,j) = rot90(squeeze(omultid(:,:,j)));
  end
  clear omultid
  enviwrite(nmultid,'multi_density',hdrstr);
  if i == 3
    load omultih.mat
  end
  nmultih = ones(n,m,p);
  for j = 1:p
    nmultih(:,:,j) = fliplr(squeeze(omultih(:,:,j))');
    %nmultih(:,:,j) = rot90(squeeze(omultih(:,:,j)));
  end
  clear omultih
  enviwrite(nmultih,'multi_height',hdrstr);
  hdrstr{2} = 'wavelength units = number of ALS echos per pixel';
  if i == 3
    load onumhits.mat
  end
  onumhits = fliplr(onumhits');
  enviwrite(onumhits,'num_echos',hdrstr);
end