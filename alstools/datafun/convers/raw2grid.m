function [grid,x,y,funargs]=raw2grid(rawdata,gridresolution,varargin)
%%RAW2GRID computes a grid from a rawpointgroup
%   [grid,x,y]=raw2grid(rawdata,gridresolution,string1,arg1,string2,arg2,...)
%       rawdata = pointgroup (ex. Laserpoints [x1,y1,z1; x2,y2,z2; ...; xn,yn,zn])
%       gridresoulution = resolution of grid (ex. [1], or [resx, resy]
%       for different resolutions in x and y direction
%   
%   raw2grid(...,'function','FUNCTION')
%       How the calculate the grid
%       function could be:
%           nextneighbor -- searching for THE next neighbor
%           idw --  (1/distance)^2
%           gauss --
%           Matlab functions: VAR, MEAN, MEDIAN, MAX, MIN, STD, LENGTH
%
%   CELLLIKE Funktions (With more input Parameters)
%   raw2grid(...,'function',{'FUNCTION', arg1, arg2,...})
%           idw: 'function',{'idw', power} -- (1/distance)^power
%               !For Weights: if distance of THE nearest Neighbor == 0, the
%               distance is taken form the next nearest Neighbor
%           gauss: 'function',{'gauss', s, m} -- y=gauss(x,s,m)
%           @matlab's: 'function',{'matlab', fun}
%               Bsp. ...,'function',{'matlab',var},...
%               For Matlab function which are used the same way like VAR,
%               MEAN, MAX.... (Make ONE number out of set of Numbers)
%
%   raw2grid(..., 'sw',POINT)
%       POINT=[x,y]
%       Lower left Point = south west Point (minx = point(1), miny=point(2))
%       starts calculation with this Point
%   raw2grid(..., 'ne',POINT)
%       POINT=[x,y]
%       Upper right Point = north east Point (maxx = point(1),maxy=point(2))
%       the last point calculated
%       
%   raw2grid(...,'perimeter',PERIMETER)
%       Perimeter where we looks for Points, if PERIMETER is a number, PERIMETER / 2
%       is used as Radius around middle of cell, if it is a Matix (ex. [px, py]), the Points are
%       searched in a box with width px,py
%
%   raw2grid(...,'nans',NAN)
%       Change if NaN must be diffrent from nan (ex. -9999), must be
%       numeric or NaN
%   raw2grid(...,'numberNeighbors',NUMBERNEIGHBORS)
%       If just max NUMBERNEIGHBORS Neighbors should be included in
%       calculation
%
%   raw2grid(...,'inRawDataNaN',inRawDataNaN)
%       changes inRawDataNaN to nan (ex. -9999 ==> nan)
%       cold be matrix [nanx,nany,nanz] (ex. [0 0 -300] ==> [nan nan nan]
%   raw2grid(...,'offset',OFFSET)
%       OFFSET = [dx dy] or
%       OFFSET = [dx dy dz]
%       rawdata = rawdata(1:3,:) - OFFSET
%   raw2grid(..., 'mask',MASK)
%       MASK = [X1,Y1; X2,Y2; X3 Y3; ...]
%       Only Data IN Maskarea will be calculated
%
%   raw2grid(...,'writeFile','FILENAME')
%       Write a ArcInfo ASCII Grid format file
%   raw2grid(..., 'waitbar')
%       Shows a waitbar, execution could be stopped
%
%    Examples:
%        gridbsp = [1:10]' * [1:10]
%        for i = 1:10
%           for j = 1:10
%               rawdatabsp(i*10 + j - 10,:)=[i-rand(1)+.5,j-rand(1)+.5,gridbsp(i,j)];
%            end
%        end
%       [grid1,x1,y1]=raw2grid(rawdatabsp,2,'waitbar');
%
%       [grid2,x2,y2]=raw2grid(rawdatabsp,[1 1], 'fun',{'idw', 9},...
%       'nan', nan, 'perimeter', 3,'waitbar');
%
%       [grid3,x3,y3]=raw2grid(rawdatabsp,1, 'fun',{'gauss', 2, 0},...
%       'nan', nan, 'perimeter', 4,'waitbar');
%
%       figure;hold on;
%       subplot(2,2,1);imagesc(gridbsp);title('original');axis equal
%       subplot(2,2,2);imagesc(x1,y1,grid1);title('nextneighbor');axis equal
%       subplot(2,2,3);imagesc(x2,y2,grid2);title('idw');axis equal
%       subplot(2,2,4);imagesc(x3,y3,grid3);title('gauss');axis equal
%
%
%       Ralph Straumann -- IDW, NumberNeighbors
%       Version 1.1
%       Lukas Wotruba, 29 April 2005
%       Version 1.3
%       Lukas Wotruba, 23 Mai 2005
%           Matlab functions, waitbar, mask, SW/NE Point
%           calcraw2grid Input changed (funargs), Error Output

%As long not changed, the fun function is:
funargs.name = 'nextneighbor';
funargs.waitbar = false;
writeFile = 0;


narg = nargin -2; % first two are obligatory
n = 1;
while n <= narg
    if ischar(varargin{n})
        if strcmpi(varargin{n},'waitbar')
            funargs.waitbar = true;
        elseif strncmpi(varargin{n},'perimeter',4) && (n < narg) && ...
                isnumeric(varargin{n+1});
            % found perimeter, PERIMETER
            funargs.perimeter.dim = length(varargin{n+1});
            funargs.perimeter.data = varargin{n+1};
            n = n + 1;
        
        elseif strncmpi(varargin{n},'sw',2) && (n < narg) && ...
                isnumeric(varargin{n+1})
            % found 'll',POINT
            funargs.xmin = varargin{n+1}(1);
            funargs.ymin = varargin{n+1}(2);
            n = n + 1;
        elseif strncmpi(varargin{n},'ne',2) && (n < narg) && ...
                isnumeric(varargin{n+1})
            % found 'll',POINT
            funargs.xmax = varargin{n+1}(1);
            funargs.ymax = varargin{n+1}(2);
            n = n + 1;


        elseif strncmpi(varargin{n},'nans',3) && (n < narg) && ...
                isnumeric(varargin{n+1})
            % found 'nans',NAN
            funargs.nans = varargin{n+1};
            n = n + 1;


        elseif strncmpi(varargin{n},'numberNeighbors',8) && (n < narg)
            % found 'numberNeighbors', NUMBER
            funargs.numberNeighbors = varargin{n+1};
            n = n + 1;

            %functions
        elseif strncmpi(varargin{n},'function',3) && (n < narg)
            if   ischar(varargin{n+1})
                % found 'function',FUN
                if strcmpi(varargin{n+1},'max') || strcmpi(varargin{n+1},'min')...
                        || strcmpi(varargin{n+1},'mean')||strcmpi(varargin{n+1},'median')...
                        || strcmpi(varargin{n+1},'var')||strcmpi(varargin{n+1},'std')...
                        ||strcmpi(varargin{n+1},'length')
                    funargs.name = 'matlab';
                    funargs.handle = str2func(varargin{n+1});
                    n = n + 1;
                elseif strncmpi(varargin{n+1},'nextneighbor',8)
                    funargs.name = 'nextneighbor';
                    n = n + 1;
                elseif strncmpi(varargin{n+1},'idw',3)
                    funargs.name = 'idw';
                    funargs.args = 2;
                    n = n + 1;
                elseif strncmpi(varargin{n+1},'gauss',4)
                    funargs.name = 'gauss';
                    funargs.args(1) = 1;
                    funargs.args(2) = 2;
                    n = n + 1;
                end
            elseif   iscell(varargin{n+1})
                % found 'function',{FUN}
                if strcmpi(varargin{n+1}{1},'matlab')
                    funargs.name = 'matlab';
                    funargs.handle = varargin{n+1}{2};
                    n = n + 1;
                elseif strncmpi(varargin{n+1}{1},'idw',3)
                    funargs.name = 'idw';
                    funargs.args = varargin{n+1}{2};
                    n = n + 1;
                elseif strncmpi(varargin{n+1}{1},'gauss',4)
                    funargs.name = 'gauss';
                    funargs.args(1) = varargin{n+1}{2};
                    funargs.args(2) = varargin{n+1}{3};
                    n = n + 1;
                end
            end

            %%%%%%%%%%%
            %below we change the rawdata
        elseif strncmpi(varargin{n},'inRawDataNaN',3) && (n < narg) && ...
                isnumeric(varargin{n+1})
            % found 'nans',NAN
            rawdata(rawdata(:,1)==varargin{n+1}(1),:)=nan;
            if length (varargin{n+1}) == 3
                rawdata(rawdata(:,1)==varargin{n+1}(2),:)=nan;
                rawdata(rawdata(:,1)==varargin{n+1}(2),:)=nan;
            end
            n = n + 1;
        elseif strncmpi(varargin{n},'offset',3) && (n < narg) && ...
                isnumeric(varargin{n+1})
            % found 'offset', [X,Y]
            rawdata(:,1) = rawdata(:,1) - offset(1);
            rawdata(:,2) = rawdata(:,2) - offset(2);
            if length (varargin{n+1}) == 3
                rawdata(:,3) = rawdata(:,3) - offset(3);
            end
            n = n + 1;
        elseif strcmpi(varargin{n},'mask') && (n < narg) && ...
                isnumeric(varargin{n+1})
            % found 'mask',[X1,Y1; X2,Y2; X3 Y3; ...]
            mask = inpolygon(rawdata(:,1),rawdata(:,2),varargin{n+1}(:,1),varargin{n+1}(:,2));
            rawdata = rawdata(mask,:);
            n = n + 1;


            %if we like to write a txt File out
        elseif strncmpi(varargin{n},'writeFile',6) && (n < narg)&& ...
                ischar(varargin{n+1})
            % found 'nans',NAN
            writeFile = 1;
            fileName = [varargin{n+1}, '.asc'];

        else
            if ischar(varargin{n})
                if n < narg
                    if ischar(varargin{n+1})
                        str=['raw2grid(...,''' varargin{n} ''',', '''' varargin{n+1} ''',...)'];
                    else
                        str=['raw2grid(...,''' varargin{n} ''',ARG,...)'];
                    end
                else
                    str=['raw2grid(...,''' varargin{n} '''])'];
                end

                warning('Error!')
                warning(str)
                warning('not known by raw2grid')
                warning(' ')
            end
        end
    end
    n = n + 1;
end
%disp
funargs
[grid, x, y] = calcraw2grid(rawdata, gridresolution, funargs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if writeFile
    fileheader{1} = ['ncols ', num2str(length(x)) ,'\r'];
    fileheader{2} = ['nrows ', num2str(length(y)) ,'\r'];
    fileheader{3} = ['xllcorner ', num2str(min(x)) ,'\r'];
    fileheader{4} = ['yllcorner ', num2str(min(y)) ,'\r'];
    fileheader{5} = ['cellsize ', num2str(gridresolution) ,'\r'];
    fileheader{5} = ['nodata_value ', num2str(nans) ,'\r'];
    fid = fopen(fileName,'w');
    for i = 1:length(fileheader)
        fprintf(fid,fileheader{i});
    end
    fclose(fid);
    dlmwrite(fileName,grid, 'delimiter',' ','-append');
end

