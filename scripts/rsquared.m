function [ rsq, rsq_adj ] = rsquared( x, y)
% Computes the R squared value between the input data points
% and the best fit line fitted to them using polyfit with an order of 1
%
%   Inputs: 
%       x - x-value vector of data points
%
%       y - y-value vector of data points
% 
%   Output:
%       rsq = R2 value
%
%       rsq_adj = adjusted R2 value
%


    P = polyfit(x, y, 1);
    yfit = polyval(P,x);
    yfit =  P(1) * x + P(2);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal
    rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(P));


end

