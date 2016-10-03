function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    p = 0.8;
    ax = subplot('Position', [(c-1 + 0.2)/m, 1-(r)/n, 1/m * p, 1/n * p])
    if(nargout > 0)
      h = ax;
    end
