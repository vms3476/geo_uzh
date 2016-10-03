function [] = write_farsite(fname,x,y,data);
%  function [] = write_farsite(fname,x,y,data);
% write a raster to a farsite raster ascii that can be read by farsite
  
% Felix Morsdorf, RSL Zurich, Oct. 2004
  
fid = fopen(fname,'w');

% mark NAN's

data(isnan(data)) = -9999;
[m,n] = size(data);

% construct header
dx = median(abs(diff(x)));dy = median(abs(diff(y)));
head{1} = ['ncols        ',num2str(m)];
head{2} = ['nrows        ',num2str(n)];
head{3} = ['xllcorner    ',num2str(min(x)-dx)];
head{4} = ['yllcorner    ',num2str(min(y)-dy)];
head{5} = ['cellsize     ',num2str(median(abs(diff(x))))];
head{6} = ['nodata_value  -9999'];

% write header
for i = 1:length(head)
  fprintf(fid,'%s\n',head{i});
end

% write data, line by line
for j = 1:n
  % write data, column by column
  for i = m:-1:1  
    fprintf(fid,'%8.4f ',data(i,j));
  end
  fprintf(fid,'%s\n',' ');
end
% close file
fclose(fid);