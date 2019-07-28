%plotaerodynamics
close all;clc;clear all;

%% some set-up
deltaangle=45;
alpha=0:deltaangle:360; %% roll
beta =0:deltaangle:180; %% pitch 
gamma=0:deltaangle:360; %% yaw
aeroscalingfactor=1;
sunscalingfactor=10;
oldpath = path; path(oldpath,'..\matlabfunctions\')
%magcv+2e-07 magfv+4e-07
%cvx-0.001447 cvy-0.003021 cvz-0.999994
%fvx-0.000000 fvy+0.577350 fvz-0.816497
%ao+315 bo+45 go +0
%% orbital environment
altitude =340000; %% in m
[rho,v]=orbitalproperties(altitude);
%windpressure=rho/2*v^2;%% pascal
%wind=[-1 0 0]' ;
windpressure=rho/2*v^2;%% pascal
wind=[-1 0 0]' ;
%wind=wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2);
solarpressure=0;%2*4.5e-6; %% pascal
sunlight=[1 1 1]'; 
sunlight=sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2);

oldalphaopt=0;oldbetaopt=00;oldgammaopt=0;
%% the control vector
controlvector=2e-07*[-0.001447 -0.003021 -0.999994]';

%% the space craft
panelsurface=0.01; %m^2
nozpanels=1;
noxpanels=0;
noypanels=0;
if nozpanels>9 || noxpanels>9 || noypanels>9
  fprintf('\n error: too many panels');
  input('error');
end


%% %% the possible forcevectors
totalforcevector = totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma,panelsurface,aeroscalingfactor,solarpressure,sunscalingfactor,windpressure, deltaangle);

%% find the good force vector
[forcevector,alphaopt,betaopt,gammaopt]=findBestAerodynamicAngles(totalforcevector,controlvector,alpha,beta,gamma,oldalphaopt,oldbetaopt,oldgammaopt);

%% the result
fprintf('\n alpha %3.1f beta %3.1f gamma %3.1f \n',alphaopt,betaopt,gammaopt);
magfv=sqrt(forcevector(1)^2+forcevector(2)^2+forcevector(3)^2);
fprintf('\n magfv%+1.0e fvx%+f fvy%+f fvz%+f',magfv,forcevector(1)/magfv,forcevector(2)/magfv,forcevector(3)/magfv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

