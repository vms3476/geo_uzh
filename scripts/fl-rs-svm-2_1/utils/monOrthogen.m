function h = monOrthogen(theta);
% Constructs an horizontal array, h(1...N) of lowpass orthonormal
% FIR filter coefficients for any even N >= 2.
% The input array (must be horizontal), alpha(1...(N/2)-1) gives (N/2)-1 free
% parameters which are angles in radians. 

%last free coefficient so that it produces orthonormal wavelets
alpha0=pi/4-sum(theta);

alpha=([theta, alpha0]);
%alpha=fliplr(alpha);

N=2*length(alpha);
h=zeros(1,N);

nstages = N/2; %number of steps in order to calculate h

lo = N/2;
hi = lo + 1;
h(lo) = cos(alpha(nstages));
h(hi) = sin(alpha(nstages));


for step = 1 : nstages-1
    %c = cos(alpha(step+1));
	%s = sin(alpha(step+1));    
    c = cos(alpha(step)); %corrected by Florian Yger
	s = sin(alpha(step)); %corrected by Florian Yger
    
    h(lo-1) = c*h(lo);
	h(lo)   = s*h(lo);
	h(hi+1) = c*h(hi);
	h(hi)   = -s*h(hi);
    
    nbsubstep= step-1; %number of substep
    substepbase=lo+1;  %place where is applied the substep
    
    for substep = 1 : nbsubstep
       hlo=h(substepbase);
       hhi=h(substepbase+1);
       
       h(substepbase)= c*hhi - s*hlo;
       h(substepbase+1)=s*hhi + c*hlo;
       
       substepbase=substepbase+2; % corrected by Carl Taswell 
    end
    
    lo=lo-1;
    hi=hi+1;
    
end

% sumH=sum(h);
% h=h*sqrt(2)/sumH; %test de normalisation du qmf