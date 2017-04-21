function [ out ] = normalize_unity( in )
% normalize the input vector to a range of 0-1

   out = (in(:) - min(in(:))) / (max(in(:)) - min(in(:)) );

end

