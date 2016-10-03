% main: disk visualisation                     (Version - 30mar2012,RL,UZH)
% -------------------------------------------------------------------------
load('pointc_fin.mat');
cc_base = [pc_fin.x_point pc_fin.y_point pc_fin.z_point pc_fin.x_pvector pc_fin.y_pvector pc_fin.z_pvector pc_fin.season pc_fin.z_point_ag];

% define box
frame_n = 258935; frame_s = 258930;   % 2 x 2 m AOI
frame_w = 669927; frame_e = 669932;

xraw_temp = cc_base(cc_base(:,1) >= frame_w,:);
xraw = xraw_temp(xraw_temp(:,1) <= frame_e,:);
yraw_temp = xraw(xraw(:,2) <= frame_n,:);
aoi_raw = yraw_temp(yraw_temp(:,2) >= frame_s,:);
cc_base = aoi_raw;
clear frame* xraw* yraw*;

% select ground echoes
cc_select = cc_base(:,8)<=0.5;
cc_ground = cc_base(cc_select,:);
   
R = 0.25; % define radius
hold on;
for k = 1:size(cc_ground,1)
    xyz = cc_ground(k,1:3);
    normalVec = cc_ground(k,4:6);
    th = linspace(0,2*pi,100);
    circ = bsxfun(@plus,xyz',R*[cos(th); sin(th); 0*th]);
    h = plot3(circ(1,:),circ(2,:),circ(3,:), 'k');
    RotationAxis = cross([0 0 1],normalVec);
    if any(RotationAxis);
    	RotationAngle = 180/pi*acos([0 0 1]*normalVec'/norm(normalVec));
        rotate(h,RotationAxis,RotationAngle,xyz);
    end;
    quiver3(xyz(1),xyz(2),xyz(3), normalVec(1),normalVec(2),normalVec(3),'b');
end
axis equal;

% select vegetation echoes
cc_select = cc_base(:,8)>0.5;
cc_veg = cc_base(cc_select,:);

% spring ------------------------------------------------------------------
cc_select = cc_veg(:,7)==1;
cc_spring = cc_veg(cc_select,:);

hold on;
for k = 1:size(cc_spring,1)
    xyz = cc_spring(k,1:3);
    normalVec = cc_spring(k,4:6);
    th = linspace(0,2*pi,100);
    circ = bsxfun(@plus,xyz',R*[cos(th); sin(th); 0*th]);
    h = plot3(circ(1,:),circ(2,:),circ(3,:), 'r');
    RotationAxis = cross([0 0 1],normalVec);
    if any(RotationAxis);
        RotationAngle = 180/pi*acos([0 0 1]*normalVec'/norm(normalVec));
        rotate(h,RotationAxis,RotationAngle,xyz);
    end
    quiver3(xyz(1),xyz(2),xyz(3), normalVec(1),normalVec(2),normalVec(3),'b');
end
axis equal;
% summer ------------------------------------------------------------------
cc_select = cc_veg(:,7)==2;
cc_summer = cc_veg(cc_select,:);

for k = 1:size(cc_summer,1)
    xyz = cc_summer(k,1:3);
    normalVec = cc_summer(k,4:6);
    th = linspace(0,2*pi,100);
    circ = bsxfun(@plus,xyz',R*[cos(th); sin(th); 0*th]);
    h = plot3(circ(1,:),circ(2,:),circ(3,:), 'g');
    RotationAxis = cross([0 0 1],normalVec);
    if any(RotationAxis);
        RotationAngle = 180/pi*acos([0 0 1]*normalVec'/norm(normalVec));
        rotate(h,RotationAxis,RotationAngle,xyz);
    end
    quiver3(xyz(1),xyz(2),xyz(3), normalVec(1),normalVec(2),normalVec(3),'b');
end
axis equal;
view([-35 20]);