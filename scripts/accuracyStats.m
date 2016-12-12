function [stats,printouts] = accuracyStats(gt,prediction,varargin)

% Provides the most common statistical means for the accuracy assessment
% of multivariate classification tasks. Measures are stored in the MATLAB
% struct 'stats' and include (name: description):
%
% - cfmat:      Confusion matrix (absolute numbers)
% - cfmat_perc: Confusion matrix (percentual numbers)
% - oa:         Overall accuracy
% - ua:         User's accuracy (Cx1 matrix)
% - pa:         Producer's accuracy (1xC matrix)
% - aa:         Average accuracy (average over Producer's accuracy)
% - kappa:      Kappa coefficient
% - precision:  positive predictive value (Cx1 matrix)
%                       true positives / (true + false positives)
% - recall:     sensitivity (Cx1 matrix)
%                       true positives / (true positives + false negatives)
% - f1:         F-score/F-measure (Cx1 matrix)
%                       2 * (precision * recall) / (precision + recall)
%
%
% The most simple form of executing the script:
%
%       stats = accuracyStats(groundTruth,prediction);
%
% calculates the above-mentioned measures for all unique labels (classes)
% present in the 'groundTruth' Nx1, 1xN or WxH matrix and 'prediction' Nx1,
% 1xN or WxH matrix. Note, however, that label 0 (zero) is by default
% treated as "background" and excluded from statistics.
%
% The script may take one or more from the following arguments list as
% comma-separated name-value pairs:
%
% - 'noDataValues':     Allows setting the "background" label which should
%                       be ignored during the acccuracy assessment.
%                       Default is zero.
%                       Multiple labels can be assigned "no data" by
%                       providing a vector of values here.
%                       To have all labels evaluated, set 'noDataValues'
%                       to NaN (Not a Number).
%
% - 'verbose':          Allows setting the verbosity level of the script:
%                       - 0: nothing printed to command line
%                       - 1: only prints oa, aa and kappa scores
%                       - 2: prints oa, ua, pa, aa, kappa, precision, recall,
%                            and f1 to command line
%                       - 3: also prints a formatted confusion matrix to
%                            the command line (percentual values)
%                       Values are printed as percentages (where
%                       applicable), with a precision of two digits.
%                       Default zero.
%
% - 'classList':        Allows restricting analyses to a subset of classes.
%                       For instance, if the 'groundTruth' vector contains
%                       labels from the set {1,2,3,4,5} and 'classList' is
%                       set to {1,3,4}, then classes 2 and 5 will be ignored
%                       (set to 'consolidationClass') during evaluation.
%                       To be used in conjunction with 'consolidationClass'
%                       only.
%
% - 'consolidationClass':
%                       Combines unconsidered classes into a single class.                       
%                       To be used together with 'classList': any labelclass
%                       not present in the 'classList' vector will be
%                       assigned the label provided in 'consolidationClass'
%                       and treated as one single class. This is particularly
%                       useful if you wish to e.g. combine multiple classes
%                       into a single background class. For instance, if the
%                       'groundTruth' vector contains labels from the set
%                       {1,2,3,4,5}, the 'classList' is specified as {1,2}
%                       and 'consolidationClass' is set to 2, all variables
%                       from classes 3, 4 and 5 will be combined into the
%                       single class 2. Evaluation will the basically be a
%                       class one versus the rest scenario.
%
%
% Examples:
%
% Print all statistics for a dataset with classes {1,2,3,4,5,6,7,8}, and...
%
% - ignore value 1 in the calculations:
%
%       stats = accuracyStats(groundTruth, prediction, ...
%                       'noDataValues', 1, ...
%                       'verbose', 3);
%
%
%
% - calculate for all values, including zero:
%
%       stats = accuracyStats(groundTruth, prediction, ...
%                       'noDataValues', NaN, ...
%                       'verbose', 3);
%
%
%
% - consolidate classes 5 and 7 into class 4:
%
%       stats = accuracyStats(groundTruth, prediction, ...
%                       'classList', [1,2,3,4,6,8], ...
%                       'consolidationClass', 4, ...
%                       'verbose', 3);
%
%
%
% - consolidate classes 3, 5, 7 and 8 into class 1. Ignore values 0 and 6
%   in the calculations:
%
%       stats = accuracyStats(groundTruth, prediction, ...
%                       'classList', [1,2,4], ...
%                       'consolidationClass', 1, ...
%                       'noDataValues', [0,6], ...
%                       'verbose', 3);
%
% 2016 mmrs/bkellenb


gt = int16(gt(:));
prediction = int16(prediction(:));

if numel(gt) ~= numel(prediction)
    error('Error: cardinalities of ground truth and prediction set are different.');
end


global noDataValues verbose classList backgroundClass classListProvided backgroundClassProvided printoutFormat;

classListProvided = false;
backgroundClassProvided = false;

noDataValues = 0;
verbose = 0;
classList = unique(gt);
backgroundClass = NaN;
printoutFormat = 'csv'; %TODO

if ~isempty(varargin)
    parseArgs(varargin);
end


% removing no data class(es)
mask = ones(size(gt),'logical');
for v=1:numel(noDataValues)
    if ~isnan(noDataValues(v)) && ~isinf(noDataValues(v))
        mask = mask & gt ~= noDataValues(v);
    end
end
gt = gt(mask);
prediction = prediction(mask);


% consolidating into background class
if classListProvided && backgroundClassProvided
    if sum(ismember(noDataValues,classList)) > 0
        warning(['Warning: "no data" value(s) is/are also present in class list. ',...
            'Excluding corresponding value(s) from class list.']);
    end
    
    for v=1:numel(noDataValues)
        if ~isnan(noDataValues(v)) && ~isinf(noDataValues(v))
            classList(classList == noDataValues(v)) = [];
        end
    end
    
    if numel(backgroundClass) > 1
        warning(['Warning: multiple "consolidationClass" values provided. ',...
            'Ignoring values beyond the first one.']);
        backgroundClass = backgroundClass(1);
    end
    
    
    if ~isnan(backgroundClass) && ~isinf(backgroundClass)
    
        if sum(ismember(backgroundClass,noDataValues)) > 0
            error('Error: "consolidationClass" value corresponds to "noData" value(s).');
        end
    
        % consolidate all classes not present in classList into background
        % class
        mask = ~ismember(gt,classList);
        gt(mask) = backgroundClass;
    
        mask = ~ismember(prediction,classList);
        prediction(mask) = backgroundClass;
    end
    
end


stats = struct();


% Confusion matrix
[stats.cfmat,classOrder] = confusionmat(gt,prediction);
stats.cfmat_perc = stats.cfmat ./ numel(gt) .* 100;
classList = classOrder;


% General stats
tp = diag(stats.cfmat);
fp = nansum(stats.cfmat,2) - tp;
fn = nansum(stats.cfmat,1)' - tp;


% Overall accuracy
stats.oa = sum(tp) / numel(gt);


% User's and Producer's Accuracies (errors of commission and omission)
stats.ua = tp ./ nansum(stats.cfmat,2);
stats.pa = tp' ./ nansum(stats.cfmat,1);

% Average accuracy
stats.aa = nanmean(stats.pa);


% Precision
stats.precision = tp ./ (tp + fp);
stats.recall = tp ./ (tp + fn);

stats.f1 = 2 .* (stats.precision .* stats.recall) ./ (stats.precision + stats.recall);


% kappa
ea = nansum((nansum(stats.cfmat,1) .* nansum(stats.cfmat,2)') ./ numel(gt)) ./ numel(gt);
stats.kappa = (stats.oa - ea) / (1 - ea);



% print to structure
printouts = struct();
if strcmp(printoutFormat,'csv')
    
    % confusion matrix
    cfmatString = sprintf('%i,',classList(1:end-1));
    cfmatString = strcat(cfmatString,sprintf('%i\r\n',classList(end)));
    for c=1:numel(classList)
        cfmatString = strcat(cfmatString,...
            sprintf('%i,',classList(c)),...
            sprintf('%.2f,',stats.cfmat_perc(c,1:end)),...
            sprintf('%.2f\r\n',stats.pa(c)));
    end
    cfmatString = strcat(cfmatString,...
        ',',...
        sprintf('%.2f,',stats.ua),...
        sprintf('%.2f\r\n',stats.oa));
    printouts.cfmat_perc = cfmatString;
    
    % general statistics
    genStatsString = sprintf('OA,AA,kappa\r\n%.2f,%.2f,%.2f\n',stats.oa,stats.aa,stats.kappa);
    printouts.generalStats = genStatsString;
    
    % precision, recall, f1
    prfString = 'class,precision,recall,f1\n';
    for c=1:numel(classList)
        prfString = strcat(prfString,...
            sprintf('%i,%.2f,%.2f,%.2f\r\n',classList(c),stats.precision(c),stats.recall(c),stats.f1(c)));
    end
    printouts.prec_rec_f1 = prfString;
end


% print to console
if verbose > 0
    
    if verbose == 3
        separationString = repmat('_',1,(2+numel(classList))*9+4);
    
        fprintf('\nCONFUSION MATRIX\n================\n\n');
        fprintf('Showing percentual values.\n\n');
        fprintf('PRED\t\\\tGT\n%s\n',separationString);
        fprintf('class\t|');
        fprintf('\t%i',classList);
        fprintf('\t| UA');
        fprintf('\n%s\n',separationString);
        for c=1:numel(classList)
                fprintf('%i\t|',classList(c));
                fprintf('\t%5.2f',stats.cfmat_perc(c,:));
                fprintf('\t| %5.2f\n',stats.ua(c) .* 100);
        end
        fprintf('%s\n',separationString);
    
        fprintf('PA\t|');
        fprintf('\t%5.2f',stats.pa .* 100);
    
        fprintf('\t| OA: %5.2f',stats.oa .* 100);
    
        fprintf('\n\n\n');
    end
    
    
    fprintf('\nSTATISTICAL MEASURES\n====================\n\n');
    fprintf('OA:\t%5.2f%%\n',stats.oa*100);
    fprintf('AA:\t%5.2f%%\n',stats.aa*100);
    fprintf('kappa:\t%5.2f\n\n',stats.kappa);
    
    if verbose > 1
        fprintf('Class\t |\tPrecision\tRecall\t\tF1\n%s\n',repmat('_',1,58));
        for c=1:numel(classList)
                fprintf('%i\t |\t%5.2f\t\t%5.2f\t\t%5.2f\n',classList(c),stats.precision(c),stats.recall(c),stats.f1(c));
        end
    end
    fprintf('\n\n');
end

end


function parseArgs(args)

global noDataValues verbose classList backgroundClass classListProvided backgroundClassProvided;

if mod(length(args),2) > 0
    error('Error: optional arguments list is missing a parameter.');
end

classListProvided = false;
backgroundClassProvided = false;


for l=1:2:length(args)
    
    if sum(strcmpi(args{l},{'noDataValues','noData','ignore','exclude'}))
        noDataValues = args{l+1};
    elseif sum(strcmpi(args{l},{'verbose','print','printout'}))
        verbose = args{l+1};
    elseif sum(strcmpi(args{l},{'classList','classes'}))
        classList = args{l+1};
        classListProvided = true;
    elseif sum(strcmpi(args{l},{'consolidationClass','backgroundClass','background'}))
        backgroundClass = args{l+1};
        backgroundClassProvided = true;
    elseif sum(strcmpi(args{l},{'args','arguments','params','parameters'}))
        parseArgs(args{l+1});
    end
    
end

if backgroundClassProvided && ~classListProvided
    warning(['Warning: "consolidationClass" provided, but no "classList" specified. ',...
        'Ignoring arguments...']);
elseif classListProvided && ~backgroundClassProvided
    warning(['Warning: "classList" provided, but no "consolidationClass" specified. ',...
        'Ignoring arguments...']);
end


end
