function [xapp,xtest,meanxapp,stdxapp] = normalizemeanstd_mtl(xapp,xtest,meanx,stdx)

% USAGE
% 
%  [xapp,xtest,meanxapp,stdxapp] = normalizemeanstd(xapp,xtest)
%
% normalize inputs and output mean and standard deviation to 0 and 1
%
% 
tol=1e-5;

nbtask=length(xapp);

nbsuppress=0;
if nargin <3
    for i=1:nbtask
        meanxapp{i}=mean(xapp{i});
        stdxapp{i}=std(xapp{i});
    end
else
    meanxapp=meanx;
    stdxapp=stdx;
end;

for i=1:nbtask
    
nbxapp=size(xapp{i},1);
indzero=find(abs(stdxapp{i})<tol);
%keyboard
if ~isempty(indzero)

    stdxapp{i}(indzero)=1;

end;
nbvar=size(xapp,2);

xapp{i}= (xapp{i} - ones(nbxapp,1)*meanxapp{i})./ (ones(nbxapp,1)*stdxapp{i}) ;

if nargin >1 & ~isempty(xtest)
    nbxtest=size(xtest{i},1);
    xtest{i}= (xtest{i} - ones(nbxtest,1)*meanxapp{i})./ (ones(nbxtest,1)*stdxapp{i} );
else
    xtest{i}=[];
end;

end