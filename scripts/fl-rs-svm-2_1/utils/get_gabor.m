function F=get_gabor(sz,theta,freq,sigmas,x0)

if length(sz)==1
    sz=sz*[1 1];
end

if nargin<2
    theta=0;
end

if nargin<3
    freq=3/(mean(sz));
end

if nargin<4
    sigmas=1/6*(sz);
end

if nargin<4
    x0=[0 0];
end

if length(sigmas)==1
    sigmas=sigmas*[1 1];
end

[x y]=meshgrid((0:sz(1)-1)-(sz(1)-1)/2,(0:sz(2)-1)-(sz(2)-1)/2);

% Rotation 
x_theta=(x-x0(1))*cos(theta)+(y-x0(1))*sin(theta);
y_theta=-(x-x0(1))*sin(theta)+(y-x0(1))*cos(theta);

F=exp(-.5*(x_theta.^2/sigmas(1)^2+y_theta.^2/sigmas(2)^2)).*cos(2*pi*freq*x_theta);