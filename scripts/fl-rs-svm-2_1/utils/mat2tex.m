function mat2tex(M, filename)

if nargin<2
filename='./matrice_en_latex.txt';
end

[n,m] = size(M);

fid = fopen(filename,'w');

fprintf(fid, '\\begin{tabular}{|');
for i=1:m
	fprintf(fid, ' c |');
end
fprintf(fid, '}\n');

fprintf(fid, '\\hline \n');
for i=1:n
	for j=1:m
		fprintf(fid, '%2.3f', M(i,j) );
		if j~=m
			fprintf(fid, ' & ');
		end
	end
	fprintf(fid, '\\\\ ');
	fprintf(fid, '\\hline');
	fprintf(fid, '\n');
end
fprintf(fid, '\\end{tabular}');
fclose(fid);
