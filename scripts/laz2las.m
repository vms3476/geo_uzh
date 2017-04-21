function [  ] = laz2las( dir )
% this function uses laszip to convert all .laz files in specified
% directory to .las format.

files = dir('*.laz');
if numel(files) > 0 
    for i = 1:numel(files)
        lazName = [files(i).name];
        lasName = [lazName(1:end-1) 's'];
        if exist(lasName, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' lazName ' -o ' lasName]);
        end
    end
end


end

