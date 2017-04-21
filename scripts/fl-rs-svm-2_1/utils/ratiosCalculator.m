function I = ratiosCalculator(band,img,params)

I = zeros(size(band,1)*size(band,2),numel(params.scales));


for i = 1:numel(params.scales)
    switch params.filters
        case 'Ratios'
            temp = band./img(:,:,params.scales(i));
        case 'NRatios'
            temp = (band-img(:,:,params.scales(i)))./(band+img(:,:,params.scales(i)));
        case 'Sum'
            temp = band+img(:,:,params.scales(i));
        case 'Prod'
            temp = band.*img(:,:,params.scales(i));
        
        otherwise
            disp('Undefined type of ratio. Exiting.')
            return
    end
    
    I(:,i) = temp(:);
end