clear
%This example demonstrates an intransitive case 

% There are the parameters of the input trapezoidal fuzzy numbers 
xAL = 45;
nxAL = 48;
nxAR = 52;
xAR = 63;

xBL = 23;
nxBL = 55;
nxBR = 58;
xBR = 64;

xCL = 1;
nxCL = 60;
nxCR = 63;
xCR = 66;

fA = [xAL,nxAL,nxAR,xAR];
fB = [xBL,nxBL,nxBR,xBR];
fC = [xCL,nxCL,nxCR,xCR];

% For the numerical integration
n=1000;

M_computed_AB = pbp_class.compute_M(fA,fB)
M_approx_AB = pbp_class. M_approx_tfn(0,1,n,fA,fB)

M_computed_BC = pbp_class.compute_M(fB,fC)
M_approx_BC = pbp_class. M_approx_tfn(0,1,n,fB,fC)

M_computed_AC = pbp_class.compute_M(fA,fC)
M_approx_AC = pbp_class. M_approx_tfn(0,1,n,fA,fC)



