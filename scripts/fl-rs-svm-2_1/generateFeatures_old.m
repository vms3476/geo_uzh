function [M,params] = generateFeatures(X,options)

% Contextual filters generator
%
% [M,params] = generateFeatures(X,options)
%
% inputs :  - X: the image to be filteres (nrows x ncols x nbands)
%           - options: options (see the main file)
% outputs: - M: the filtered features
%          - params: structure containing the global parameters of the filters generated 
%
% 2012 devis tuia and michele volpi

Sz = size(X);

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
}';

%size of the moving window or attribute range (for ATT filters)
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
}';

%structuring elements
possible.shapeel   = {'disk','diamond','square','line'};%

% how many bands to be filtered in a single minibatch
possible.nbands = options.Nbands;

weigthsfilters = ones(1,length(possible.filters));

probsel = weigthsfilters./sum(weigthsfilters);
possible.bandweights = scale(options.bandweights')';

possible.bandweights = possible.bandweights./sum(possible.bandweights);


%% generate the structure with the filter specs
params = generateRandParams(probsel,possible,options.Nscales);

%clear for parameters not required by feature type
switch params.filters
    case {'AVG','ENT','STD','RAN','ATT-a','ATT-i','ATT-d','ATT-s'}
    params.formel = 'niet';
    params.angles = [];
    otherwise
end


%disp(params)
toaddw = zeros(1,length(possible.filters));

disp(params)


%% Generate the filters.
M = [];
for jj = 1:length(params.bands)
    
    MM = contextualfeatures(X(:,:,params.bands(jj)),params); 
    M = [M MM];
    
    
end




