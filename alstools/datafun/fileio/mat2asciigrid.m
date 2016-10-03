function [] = mat2asciigrid(fname,x,y,data);
%  function [] = mat2asciigrid(fname,x,y,data);
%    write a raster to a esri raster ascii that can be read by farsite
  
% Felix Morsdorf, RSL Zurich, Oct. 2004, last edit August 2008
  
fid = fopen(fname,'w');

% mark NAN's

data(isnan(data)) = -9999;
[rows,cols] = size(data);
data = flipud(data);
% construct header
dx = median(abs(diff(x)));dy = median(abs(diff(y)));
head{1} = ['ncols        ',num2str(cols)];
head{2} = ['nrows        ',num2str(rows)];
head{3} = ['xllcorner    ',num2str(min(x)-(dx/2))];
head{4} = ['yllcorner    ',num2str(min(y)-(dy/2))];
head{5} = ['cellsize     ',num2str(median(abs(diff(x))))];
head{6} = ['nodata_value  -9999.0000'];

% write header
for i = 1:length(head)
  fprintf(fid,'%s\n',head{i});
end

% write data, line by line
for j = 1:rows
  % write data, column by column
  for i = 1:cols  
    fprintf(fid,'%8.4f ',data(j,i));
  end
  fprintf(fid,'%s\n','');
end
% close file
fclose(fid);