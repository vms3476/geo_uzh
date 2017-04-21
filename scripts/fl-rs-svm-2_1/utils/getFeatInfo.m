function Info = getFeatInfo(params,ind)

Info.filters = params.filters;
Info.bands = [];
Info.scales= [];
Info.angles = [];
Info.formel  = params.formel;
Info.pos = [];
Info.angles  = [];%params.angles(idx.scales);
Info.refined = 0;
Info.level = 0;

clear idx
switch params.filters
    
    %case {'AVG','ENT','STD','RAN'}
    
    case {'ATT-a','ATT-i','ATT-d','ATT-s'} % ocio all'uscita in stack e invertita che dipende dal .mex d Mauro
        
        idx.bands  = ceil(ind/(2*length(params.scales)));
        idx.scales = 2*length(params.scales) - ((idx.bands*(2*length(params.scales)))-ind);
        
        Info.bands   = params.bands(idx.bands);
        
        if idx.scales <= length(params.scales)
            temp = params.scales(end:-1:1); % inverti l'ordine delle scale: a 1 2 3 in uscita corrisponde 3 2 1 delle scale
            Info.scales = temp(idx.scales); % prendi la finestra che ha generato la feature numero 'ind'
            
            Info.pos = 'left'; % set ATT opening o closing (dx o sx)!
            clear temp
        else % senno, prendi la finestra originale delle features in ind + length(params.scales)
            Info.scales = params.scales(idx.scales-length(params.scales));
            Info.pos = 'right'; % set ATT opening o closing (dx o sx)!
            
        end
        
    otherwise
        
        
        idx.bands  = ceil(ind/length(params.scales));
        idx.scales = length(params.scales) - ((idx.bands*(length(params.scales)))-ind);
        Info.bands   = params.bands(idx.bands);
        Info.scales  = params.scales(idx.scales);
        Info.pos = ''; % set ATT opening o closing (dx o sx)!
        
        
        
        %             Info.formel = params.formel;
        %             Info.pos = '';
        if strcmp(Info.formel,'line')
            
            Info.angles = params.angles(idx.scales);
            
        end
end