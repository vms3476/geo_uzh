% script to load raw data from small area as of 2002
for j = 1:2
  if j == 1
    cd /Users/morsdorf/Teaching/MSc/MA_SW/RAW2002/first
  elseif j == 2
    cd /Users/morsdorf/Teaching/MSc/MA_SW/RAW2002/last
  end
    filez = dir('*.ASC');
    raw.x = [];
    raw.y = [];
    raw.z = [];
    for i = 1:length(filez);
      disp(['Loading ',filez(i).name,' ...'])
      dat = load(filez(i).name);
      raw.x = [raw.x;dat(:,2)];
      raw.y = [raw.y;dat(:,1)];
      raw.z = [raw.z;dat(:,3)];
    end
    if j == 1
      save first raw
    elseif j == 2
      save last raw
    end
end