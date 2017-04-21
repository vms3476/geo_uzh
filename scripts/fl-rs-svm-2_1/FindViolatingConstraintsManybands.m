function  [Xviol,Info,maxviol,ind]=FindViolatingConstraintsManybands(X,Y,err,options,optionsFeatSel,params);

%% Assessing if there are features violating the constraints
Xy=X.*(Y*ones(1,size(X,2)));
constraintviol=abs(Xy'*err);


if  options.randSamp == 1 %random selection
    maxviol = 1.5;
    c = randperm(size(X,2))';
    ind = c(1,1);
    clear c
else %active set
    [maxviol,ind]=max(constraintviol*options.C/size(X,1));
end


%% If we find a feature violating the constraint by at least the tolerance
% return it (Xviol)
% SAVE its characteristics (Info)

if maxviol > 1 + optionsFeatSel.tolerance
    Xviol=X(:,ind);
   
 
    Info.filters = params.filters;
    Info.bands = [];
    Info.scales= [];
    Info.angles = [];
    Info.formel  = params.formel;
    Info.pos = [];

    Info.angles  = [];
    
   
    clear idx
    switch params.filters
        
            
        case {'ATT-a','ATT-i','ATT-d','ATT-s'} % output order of attribute profiles depends on the .mex used
            
            idx.bands  = ceil(ind/(2*length(params.scales)));
            idx.scales = 2*length(params.scales) - ((idx.bands*(2*length(params.scales)))-ind);
            
            Info.bands   = params.bands(idx.bands);
            
            if idx.scales <= length(params.scales)
                temp = params.scales(end:-1:1); % we invert the scales order for the left side: 1 2 3 in output corresponds to 3 2 1 of scales
                Info.scales = temp(idx.scales); % take the feature related to number 'ind'
                
                Info.pos = 'left'; % set ATT opening o closing (dx o sx)!
                clear temp
            else % otherwise use the original order of the window sizes in ind + length(params.scales)
                Info.scales = params.scales(idx.scales-length(params.scales));
                Info.pos = 'right'; % set ATT opening o closing (dx o sx)!
                
            end
            
        otherwise
             
            
            idx.bands  = ceil(ind/length(params.scales));
            idx.scales = length(params.scales) - ((idx.bands*(length(params.scales)))-ind);
            Info.bands   = params.bands(idx.bands);
            Info.scales  = params.scales(idx.scales);
            Info.pos = ''; % set ATT opening o closing (dx o sx)!
            

            if strcmp(Info.formel,'line')
                
                Info.angles = params.angles(idx.scales);
                
            end
    end
    
else
    Xviol=[];
    Info=[];
end;

if  options.randSamp == 1
    maxviol=constraintviol(ind,1)*options.C/size(X,1);
end

end

