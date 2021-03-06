%Filter function;
%Contains different filter options
%Comment out all but one of them
function sigma=filter_exp(x)

%%%%%%%%%%%%%%%%%%
%Cesaro
%%%%%%%%%%%%%%%%%%
%sigma= 1-x;

%%%%%%%%%%%%%%%%%%
%Raised Cosine
%%%%%%%%%%%%%%%%%%
%sigma=(1+cos(pi*x))/2;

%%%%%%%%%%%%%%%%%%
%Lanczos
%%%%%%%%%%%%%%%%%%
%if x ~=0
%	sigma = sin(pi*x)/(pi*x);
%else
%	sigma = 1;
%end

%%%%%%%%%%%%%%%%%%
%Exponential
%%%%%%%%%%%%%%%%%%
nuc=.05;
alpha=16;
p=2;
if abs(x)>nuc
	sigma=exp(-1*alpha*((abs(x)-nuc)/(1-nuc))^2);
else
	sigma=1;
end









