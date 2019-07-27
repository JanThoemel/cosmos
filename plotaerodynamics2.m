%plotaerodynamics
close all;clc;clear all;

%% some set-up
deltaangle=10;
alpha=0:deltaangle:360; %% yaw
beta =0:deltaangle:180; %% pitch 
gamma=0:deltaangle:360; %% roll
aeroscalingfactor=1;
sunscalingfactor=10;
oldpath = path; path(oldpath,'..\matlabfunctions\')

%% orbital environment
altitude =340000; %% in m
[rho,v]=orbitalproperties(altitude);
windpressure=rho/2*v^2;%% pascal
wind=[-1 0 0]' ;
wind=wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2);
solarpressure=0;%2*4.5e-6; %% pascal
sunlight=[1 1 1]'; 
sunlight=sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2);

%% the space craft
panelsurface=0.01; %m^2
nozpanels=1;
noxpanels=0;
noypanels=0;
if nozpanels>9 || noxpanels>9 || noypanels>9
  fprintf('\n error: too many panels');
  input('error');
end

%% the control vector
controlvector=[-1 0.6 0]';

%% %% the possible forcevectors
totalforcevector = totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma,panelsurface,aeroscalingfactor,solarpressure,sunscalingfactor,windpressure, deltaangle);

%% find the good force vector
oldalphaopt=0;oldbetaopt=90;oldgammaopt=320;
[alphaopt,betaopt,gammaopt]=findBestAerodynamicAngles(totalforcevector,controlvector,alpha,beta,gamma,oldalphaopt,oldbetaopt,oldgammaopt);


%% the result
fprintf('\n alpha %3.1f beta %3.1f gamma %3.1f \n',alphaopt,betaopt,gammaopt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

