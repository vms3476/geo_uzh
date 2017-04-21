function  [xn minim maxim] = scale(x, minim, maxim)

% function  [xn minim maxim] = scale(x, minim, maxim)
%
% Funcion que escala cada caracteristica entre [0,1]
% Las muestras en filas y las caracteristicas en columnas.
%
% x:     senyal a escalar con las caracteristicas en columnas.
% minin: opcional, vector de minimos por caracteristica a utilizar.
% maxim; opcional, igual que minim pero para los maximos.
%
% xn:    senyal escalada.
% minim: vector de minimos por caracteristica encontrados en la señal x.
% maxim: idem para maximos.

xn = x;
[numsam,numfeat] = size(xn);

if nargin == 1 % Nuevo escalado
    maxim = max(xn);
    minim = min(xn);
elseif numfeat ~= length(minim) || numfeat ~= length(maxim)
    % Escalado previo, comprueba
    error(['Num. feat. ~= de longitud de minim o maxim ' ...
        num2str(length(minim)) ', ' num2str(length(maxim))]);
end

for j=1:numfeat
    if maxim(j) ~= minim(j)
        xn(:,j)=(xn(:,j)-minim(j))/(maxim(j)-minim(j));
    else
       disp(['Warning: max = min, not normalizing index ' num2str(j)])
       
    end
end
