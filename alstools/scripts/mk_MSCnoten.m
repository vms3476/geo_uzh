% script for producing histograms of the MSC marks



UNITS = {'giscience-gis','giscience-geocomputation','economic-geography', ...
         'h2k','rsl','political-geography','human-geography','human-geography', ...
         'physische-geographie-boden-biogeographie', ...
         'physical-geography-glaciology-and-geomorphodynamics', ...
         'physical-geography-glaciology-and-geomorphodynamics','giva'};
PROFS = {'weibel','purves','berndt','seibert','schaepman','korf','ller-b', ...
         'backhaus','schmidt','haeberli','vieli','fabrikant'};
if 0
[num,txt,dat] = xlsread('/Users/morsdorf/Desktop/MSc_ab-2008_20141022.xls');
for i = 1:length(txt) 
    
    str = txt{i,17};
    if length(str) > 4
        II = zeros(size(PROFS));
        for j = 1:length(PROFS)
            ii = strfind(lower(str),PROFS{j});
            if ~isempty(ii)
                II(j) = ii;
            end
        end
        if sum(II) > 0 & isempty(UNIT{i})
            II(II==0) = 1000;
            UNIT{i} = UNITS{find(min(II)==II)};
        else
            ii = findstr(str, '"');
            if length(ii) == 2;
                %str = str(ii(1):ii(2));
                str = str(ii(1):ii(1)+80);
            end
            % Google
            search_string = ['http://www.google.com/search?hl=en&source=hp&q=',str, ...
                            '+',txt{i,6},'+',txt{i,5}]; 
            search_string = strrep(search_string,' ','%20');
            %BING
            %search_string = ['http://www.bing.com/search?q=',str];
            %search_string = strrep(search_string,' ','+');
            S = urlread(search_string);
            if (~isempty(S))
                ii = strfind(S,'www.geo.uzh.ch/en/units');
                if length(ii) >= 1
                    for j = 1:length(ii)
                        str = S(ii+24:ii+24+100);
                        jj = strfind(str,'/');
                        unit{j} = str(1:jj(1)-1);
                    end
                    for j = 1:length(unit)
                        id(j) = strcmp(unit{1},unit{j});
                    end
                    if length(id) == sum(id)
                        UNIT{i} = unit{1};
                    else
                        disp('Units do not match, pick manually');
                        UNIT{i} = [unit{:}];
                    end
                    clear unit
                end
            else
                disp(['Error downloading Google query at index ',num2str(i)]);
            end
        end
    end
end
save /Users/morsdorf/Desktop/mscnoten.mat num dat txt UNIT 
end


load  /Users/morsdorf/Desktop/mscnoten.mat num dat txt UNIT 

% manual edits

UNIT(411) = UNITS(2);
UNIT(408) = UNITS(4);
UNIT(407) = UNITS(7);

UNITS = {'giscience-gis','giscience-geocomputation','economic-geography', ...
         'h2k','rsl','political-geography','human-geography', ...
         'physische-geographie-boden-biogeographie', ...
         'physical-geography-glaciology-and-geomorphodynamics','giva'};

for i = 1:length(UNITS)
    k = 0;
    for j = 1:length(UNIT)
      ii = strmatch(UNITS{i},UNIT{j});
      if ii == 1
          k = k+1;
          II(k) = j;
      end
    end
    subplot(3,4,i)
    hist(num(II,1),[4 4.5 5 5.5 6])
    set(gca,'xlim',[3.5 6.5]);
    title(upper(UNITS{i}));
    legend([],['n = ',num2str(length(II)),', mean = ',num2str(mean(num(II,1))), ', median = ',num2str(median(num(II,1)))]);
    clear II;           
end
