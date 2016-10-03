% make histograms
dat = csvread('/Users/morsdorf/Desktop/geo233_2011.csv');
str = {'Adrian','Andi','Felix'};
for i = 1:3
  figure(1)
  orient tall
  subplot(3,1,i)
  hist(dat(:,i),[0:5:100])
  xlabel('Percent [%]')
  ylabel('Number of Stundents')
  title(['Distribution of points - GEO233 HS2011 - ',str{i}])
end
print -depsc2 /Users/morsdorf/Desktop/geo233_hist.eps