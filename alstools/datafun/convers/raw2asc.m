function [] = raw2asc(fname,raw,fields);
% function [] = raw2asc(fname,raw,fields);
% exports raw ALS data from structure raw (as read by readlas) to ASCII table
% format
%   Input:
%          fname  - file name of ASCII file to be written
%          raw    - structure of raw laser data as parsed by readlas
%          fields - cellstring variable with fields to be exported
%                   if omitted, all fields will be exported
%       
%   Output:
%          ASCII file with name "fname" in current directory

% Felix Morsdorf, RSL Zuerich, Nov. 2011

% parse arguments
fnms = fieldnames(raw);

% exclude fields that should not be written to file
if nargin == 3
    for i = 1:length(fields)
        ii = find(strcmp(fnms,fields{i}));
        if ~isempty(ii)
            fnms(ii) = [];
        end
    end
    for i = 1:length(fnms)
        raw = rmfield(raw,fnms{i});
    end
    fnms = fieldnames(raw);
else
    fields = fnms;
end

head = [];
frmt = [];
data = ones(length(raw.x),length(fnms)) * NaN;
for i = 1:length(fields)
    head = [head,' ',fields{i}];
    frmt = [frmt,' %8.4f ']; 
    eval(['data(:,i) = double(raw.',fields{i},');'])
end

frmt = [frmt,'\n'];
fid = fopen(fname,'w');
fprintf(fid,'%s\n',head);
for i = 1:length(raw.x)
    fprintf(fid,frmt,data(i,:));
end
fclose(fid);