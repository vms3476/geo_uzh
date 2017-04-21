function MM = prepareContextData2(X,Sz,dataset)

if strcmp(dataset,'GRSS13')

    Y = X(:,end);
    X = X(:,1:end-1);
end
    
%[e,Xp] = princomp(X);

Xp = X;%Xp(:,1:3);

if strcmp(dataset,'GRSS13')
    Xp = [Xp Y];
end

%Xp = reshape(Xp,Sz(1),Sz(2),size(Xp,2));

MM = [];

for i = 1:size(Xp,2)
    params.filters = 'TH';
    params.scales = 1:2:11;
    %fprintf('Band %i, filter %s, wins %i\n',i, params.filters, length(params.scales))
    M = computeMorphoMulti(Xp(:,i),params,Sz);
    MM = [MM M];
end
%
for i = 1:size(Xp,2)
    params.filters = 'OCR';
    params.scales = 1:2:11;
    %fprintf('Band %i, filter %s, wins %i\n',i, params.filters, length(params.scales))
    M = computeMorphoMulti(Xp(:,i),params,Sz);
    MM = [MM M];
end

for i = 1:size(Xp,2)
    params.filters = 'OCRTH';
    params.scales = 1:2:11;
    %fprintf('Band %i, filter %s, wins %i\n',i, params.filters, length(params.scales))
    M = computeMorphoMulti(Xp(:,i),params,Sz);
    MM = [MM M];
end

for i = 1:size(Xp,2)
    params.filters = 'ATT-a';
    params.scales = 10:1000:10000;
    %fprintf('Band %i, filter %s, wins %i\n',i, params.filters, length(params.scales))
    M = computeMorphoMulti(Xp(:,i),params,Sz);
    MM = [MM M];
end
%
%
for i = 1:size(Xp,2)
    params.filters = 'ATT-d';
    params.scales = 10:10:100;
    %fprintf('Band %i, filter %s, wins %i\n',i, params.filters, length(params.scales))
    M = computeMorphoMulti(Xp(:,i),params,Sz);
    MM = [MM M];
end

ii = sum(MM) == 0;

MM(:,ii) = [];