function [model,out,ax] = add_regress(x,y,str,plt);
%  function [model,out,ax] = add_regress(x,y,str,plt);
%  adss a regression to current graph, based on regress from statistics toolbox
%  adds column of ones for constant term

  
% Felix Morsdorf, RSL Zurich, March 2006
  
  ax = gca;
  if nargin == 0
    x = get(gca,'Xdata');
    y = get(gca,'Ydata');
    type = 'poly1';
  elseif nargin == 2
    type = 'poly1';
  elseif nargin == 1
    type = x;
    x = get(gca,'Xdata');
    y = get(gca,'Ydata');
  end
  
  x = x(:);
  y = y(:); 
  
  % append column of ones (constant term)
  dat = [y,ones(size(x))];
  
  % use statistics toolbox to remove outliers
  %[B,BINT,R,RINT,STATS] = regress(x,dat,0.05);  
  %ii = find(prod(RINT') > 0);
 
  %outlx = x(ii);
  %outly = y(ii);
  %x(ii) = [];
  %y(ii) = [];
  
  nx = linspace(min(x),max(x),20);
  method = 1;
  if method == 1
    % use robustfit (if available)
    [BR, STATSR] = robustfit(x,y);
    ii = STATSR.w < 0.00;
    outlx = x(ii);
    outly = y(ii);
    x(ii) = [];
    y(ii) = [];
    yf = BR(1)+nx*BR(2);
    yhat = BR(1)+x*BR(2);
    mod(1) = BR(2);
    mod(2) = BR(1);
    delta2 = STATSR.se(1);
  elseif method == 2
    % use polyfit
    [mod,S,MU] = polyfit(x,y,1);
    [yf,delta2] = polyval(mod,nx,S,MU);
    [yhat,delta] = polyval(mod,x,S,MU);
  elseif method == 3
    [B,BINT,R,RINT,STATS] = regress(x,y);
    yf = B(1)+nx*B(2);
    yhat = B(1)+x*B(2);
    mod(1) = B(2);
    mod(2) = B(1);
  end
  
  % compute r2 and rmse according to chen's paper  
  len = length(x);
  R2 = 1 - (len-1)/(len-2)*sum((y-yhat).^2)/sum((y-mean(y)).^2);
  
  
  RMS = sqrt(sum((y-yhat).^2)/(len-2));
  
  % construct output parameters
  model.p1 = mod(2);
  model.p2 = mod(1);
  
  %out.rsquare = STATS(1);
  out.rsquare = R2;
  
  if nargin == 4
    return
  end

  % plot data  
  hp(1) = plot(x,y,'.b','markersize',7);
  hold on
  hp(2) = plot(nx,yf,'-k');
  hold on
  
  plot(outlx,outly,'.r','markersize',7);
  hp(3) = plot(nx,yf+delta2,'--r');
  hp(4) = plot(nx,yf-delta2,'--r');
  hp(5) = plot([0,ceil(max(x))+1],[0,ceil(max(y))+1],'linewidth',2,'color',[0.5 0.5 0.5]);
  axis([0,ceil(max(x))+1,0,ceil(max(y))+1]);
  if nargin >= 3
    xlabel(str{1});
    ylabel(str{2});
  end
  ax(1) = gca;  
  set(hp,'linewidth',2);
  hl = legend([hp(1),hp(2),hp(4)],'Original Data','Linear Fit','95 % Conf. Interval', ...
              'location','southeast');
  box on;
  grid on
  % print statistics 
  X = [ones(size(y,1),1),x];
  regf = @(XTRAIN,ytrain,XTEST)(XTEST*regress(ytrain,XTRAIN));
  cvMse = crossval('mse',X,y,'predfun',regf)*100;  
  nstr  = ['Number of Samples : ',num2str(length(x))];
  mstr  = ['y = ',num2str(mod(2)),' + ',num2str(mod(1)),' * x'];
  rstr  = ['R^2 = ',num2str(fix(R2*100)/100)];
  %  Rstr  = ['RMS = ',num2str(fix(cvMse*10000)/10000)];
  Rstr  = ['RMS = ',num2str(fix(RMS*10000)/10000)];

  pos = get(gca,'position');
  axes('position',pos,'visible','off','units','normalized');
  ht = text(0.05,0.95,strvcat(nstr,mstr,rstr,Rstr), ...
  	    'verticalalign','top','fontsize',10);
  set(ht,'fontweight','bold');
  set(gca,'layer','top','nextplot','add','tag','stats');
  setstrings(12)
  
  ax(2) = gca;
  %keyboard
  %  axes(ax(1));


