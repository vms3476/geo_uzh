function [M,params,levs] = generateFeatures(X,options)



Sz = size(X);

%Sz(3)

possible.bands     = [1:Sz(3)];
possible.filters   = { % NOTA: LE GLCM NON SONO INCLUSE, MA SE CASO BASTA AGGIUNGERE GLI ID QUA!
    'OP'    ,
    'CL'    ,
    'TH_OP' ,
    'TH_CL' ,
    'OR'    ,
    'CR'    , 
    'TH_OR' , 
    'TH_CR' ,
    'AVG'   ,
    'ENT'   ,
    'STD'   ,
    'RAN'   ,
    'ATT-a' ,
    'ATT-d' , 
%    'ATT-i' , 
%    'ATT-s'
    'Ratios',
    'NRatios',
    'Sum',
    'Prod'
}';

possible.winsize   = {
[1:2:21],         %     'OP'    ,
[1:2:21],         %     'CL'    ,
[1:2:21],         %     'TH_OP' ,
[1:2:21],         %     'TH_CL' ,
[1:2:21],         %     'OR'    ,
[1:2:21],         %     'CR'    , 
[1:2:21],         %     'TH_OR' , 
[1:2:21],         %     'TH_CR' ,
[3:2:15],       %     'AVG'   ,
[3:2:15],       %     'ENT'   ,
[5:2:21],       %      'STD'   ,
[5:2:21],       %     'RAN'   ,
[10:1000:10000],  %     'ATT-a' ,
[10:10:100],    %     'ATT-d' ,
%[0.1:0.1:0.7],    %     'ATT-i' , 
%[0.50:5:50],    %     'ATT-s' }';
[1:Sz(3)],
[1:Sz(3)],
[1:Sz(3)],
[1:Sz(3)]
}';
possible.shapeel   = {'disk','diamond','square','line'};%
possible.nbands = options.Nbands;
%%% da aggiornare a ogni iterazione di infsvm
%quantewind = 5; % if == 0, genera a caso entro 1 e 5;
weigthsfilters = ones(1,length(possible.filters));

% %gamma = 0.1;
% %factor = 1.001;
% %ind = 4;

probsel = weigthsfilters./sum(weigthsfilters);
possible.bandweights = scale(options.bandweights')';

possible.bandweights = possible.bandweights./sum(possible.bandweights);




params = generateRandParams(probsel,possible,options.Nscales);
    
%clear for parameters not required by feature type
switch params.filters
    case {'AVG','ENT','STD','RAN','ATT-a','ATT-i','ATT-d','ATT-s','Ratios','NRatios','Sum','Prod'}
    params.formel = 'niet';
    params.angles = [];
    otherwise
end

%disp(params)
toaddw = zeros(1,length(possible.filters));
M = [];
levs = []; %initialize (just in case)
disp(params.filters)
for jj = 1:length(params.bands)
    
    switch params.filters
        case {'Ratios','NRatios','Sum','Prod'}
            MM = ratiosCalculator(X(:,:,params.bands(jj)),X,params);
        otherwise
            MM = contextualfeatures(X(:,:,params.bands(jj)),params);
    end
    
    if options.ASH == 1
       levs = [levs; repmat(options.lev(params.bands(jj)),size(MM,2),1)]; 
    end
    
    M = [M MM];
    
    
end




