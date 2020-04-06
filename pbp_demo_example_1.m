clear

%% FUZZY NUMBERS
% There are the parameters of the input trapezoidal fuzzy numbers 
xAL = 3;
nxAL = 6;
nxAR = 9;
xAR = 11;

xBL = 2;
nxBL = 9.75;
nxBR = 10.25;
xBR = 10.5;

fA = [xAL,nxAL,nxAR,xAR];
fB = [xBL,nxBL,nxBR,xBR];

% For the numerical integration
n = 1000;

M_computed = pbp_class.compute_M(fA,fB)
M_approx = pbp_class. M_approx_tfn(0,1,n,fA,fB)


