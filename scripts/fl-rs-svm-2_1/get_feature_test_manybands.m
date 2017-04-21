function [Xf]=get_feature_test_manybands(X,feat,sizes) 
% feats = info dans "FindViolatingConstraintsManbands.m"

Xf = zeros(size(X,1),length(feat));

if size(X,3) == 1
    X = reshape(X,sizes(1),sizes(2),size(X,2));
end

%d=size(X,2);
% Xf=[];%zeros(size(X,1),length(feat));


for i=1:length(feat)
   

%     temp=X(:,feat(i).indice);
      
switch feat(i).filters
    case {'Ratios','NRatios','Sum','Prod'}
        temp = ratiosCalculator(X(:,:,feat(i).bands),X,feat(i));
    otherwise
        temp = contextualfeatures(X(:,:,feat(i).bands),feat(i));
end




    switch feat(i).filters
                    
        case {'ATT-a','ATT-i','ATT-d','ATT-s'} % output order of attribute profiles depends on the .mex used
            if strcmp(feat(i).pos,'left') % the order is  3 2 1 | 1 2 3 -> (if we have 3 scales)
                temp = temp(:,1); % so by using 'pos', we invert the scales order for the right side: 1 2 3 in output corresponds to 3 2 1 of scales
            else % otherwise use the original order of the window sizes in ind + length(params.scales)
                temp = temp(:,2); %
            end
    end
    
    temp=(temp-ones(size(temp,1),1)*mean(temp,1))./(ones(size(temp,1),1)*std(temp,1));

    Xf(:,i)=[temp];
    
end



end
