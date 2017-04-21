function [raw] = subsetraw(raw,ii);
%  function [raw] = subsetraw(raw,ii);
%  Subset structure raw with logical array ii of length(raw.x)
%  works for all field of raw
  
% Felix Morsdorf, RSL Zurich, April 2009
  
  if islogical(ii) 
    fnames = fieldnames(raw);
    for i = 1:length(fnames)
      eval(['raw.',fnames{i},' = raw.',fnames{i},'(ii);']);
      eval(['raw.',fnames{i},' = raw.',fnames{i},'(:);']);
    end
  elseif isstruct(ii)
      are = [min(ii.x) max(ii.x) min(ii.y) max(ii.y)];
      ii = raw.x >= are(1) & raw.x <= are(2) & raw.y >= are(3) & raw.y <= are(4);
      fnames = fieldnames(raw);
      for i = 1:length(fnames)
          eval(['raw.',fnames{i},' = raw.',fnames{i},'(ii);']);
      end
  else
    fnames = fieldnames(raw);
    for i = 1:length(fnames)
      eval(['raw.',fnames{i},' = raw.',fnames{i},'(ii);']);
    end
  end
  
  