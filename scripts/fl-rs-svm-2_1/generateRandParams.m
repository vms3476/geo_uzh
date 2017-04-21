function params = generateRandParams(weigthsfilters,possible,quantewind)

% Random filters generator: generates a structure of filter parameters
% selects 
% - a type, 
% - a set of bands to be filtered
% - a set of filter parameters (sizes, orientations, ...)
%
% and saves them in the structure params
%
% params = generateRandParams(weigthsfilters,possible,quantewind)

% 2012 devis tuia and michele volpi

if nargin < 3
    quantewind = ceil(rand(1)*5); % random between 1 e 10 wins
elseif quantewind == 0
    quantewind = ceil(rand(1)*5); % random between 1 e 10 wins
end

bomb = rand(1);
probsel = cumsum(weigthsfilters);

ind = find(bomb < probsel);
ind = ind(1);

probselBands = cumsum(possible.bandweights);
indBB = [];%zeros(1,possible.nbands)
params.bands = [];

while size(params.bands,2) < possible.nbands
    bombB = rand(1);
    indB = find(bombB < probselBands);
    
    if isempty(indB); indB = length(probselBands); end %limit case when last band is selected
    
    params.bands = [params.bands possible.bands(indB(1))];
    probselBands(indB(1)) = [];
    possible.bands(indB(1)) = [];
end


% TIPO FILTRO
params.filters = possible.filters{ind};

% WINDOW SIZE (SCALE)
ranran = randperm(length(possible.winsize{1,ind}));
params.scales = sort(possible.winsize{1,ind}(ranran(1:quantewind))); clear ranran

% SHAPE OF STRUCTURING ELEMENT
ranran = randperm(length(possible.shapeel));
params.formel = possible.shapeel{1,ranran(1)}; clear ranran

%add angle if line element
if strcmp(params.formel,'line') == 1
    params.angles = ceil(180*rand(1,quantewind));
else 
    params.angles = [];
end

if strcmp(params.filters,'ATT-i')
    params.bands = params.bands(1:3);
    ranran = randperm(length(possible.winsize{1,ind}));
    params.scales = sort(possible.winsize{1,ind}(ranran(1:3))); clear ranran 
end

if strcmp(params.filters,'ATT-s')
    params.bands = params.bands(1:10); %limit ATT-s to 10 bands
  
end
