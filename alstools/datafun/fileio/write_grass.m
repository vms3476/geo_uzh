function [] = write_farsite(fname,x,y,data);
%  function [] = write_farsite(fname,x,y,data);
% write a raster to a farsite raster ascii that can be read by farsite
  
% Felix Morsdorf, RSL Zurich, Oct. 2004
  
fid = fopen(fname,'w');

% mark NAN's

data(isnan(data)) = -9999;
[m,n] = size(data);
% print header

head{1} = ['north:    ',num2str(max(y))];
head{2} = ['south:    ',num2str(min(y))];
head{3} = ['east:     ',num2str(max(x))];
head{4} = ['west:     ',num2str(min(x))];
head{5} = ['cols:        ',num2str(m)];
head{6} = ['rows:        ',num2str(n)];
%head{8} = ['nodata_value  -9999'];

for i = 1:length(head)
  fprintf(fid,'%s\n',head{i});
end

% write data
for i = 1:m
  for j = 1:n
  %fprintf(fid,'%s\n',num2str(data(:,i)'));
  fprintf(fid,'%4.4f ',data(i,j));
  end
  fprintf(fid,'%s\n',' ');
end
% close file

fclose(fid);