function [pos] = area2pos(area);
  if min(size(area)) > 1
    pos = [area(:,1) area(:,3) area(:,2)-area(:,1) area(:,4)-area(:,3)];
  else
    pos = [area(1) area(3) area(2)-area(1) area(4)-area(3)];
  end