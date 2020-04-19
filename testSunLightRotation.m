clear all; clc; close all;

oldpath = path; path(oldpath,strcat('..',filesep,'matlabfunctions',filesep))
oldpath = path; path(oldpath,strcat('..',filesep,'cosmosSupport',filesep))

equatorInclination=23.44; %% [deg]
orbitInclination=97;     %% [deg]

inclinationOrbitalPlane=equatorInclination+orbitInclination-90; %% [deg], this depends on the RAAN, here in winter
inclinationOrbitalPlane
sunLightAtAN=[sind(inclinationOrbitalPlane) cosd(inclinationOrbitalPlane) 0];
rotatedSunForceVector=rodrigues_rot([0 1 0],[0 0 1],-inclinationOrbitalPlane/180*pi);
vectarrow([0 0 0],sunLightAtAN);hold on;
pause(1)
vectarrow([0 0 0],rotatedSunForceVector);
return;

angle=30;
for i=0:angle:210
  sunLight=rodrigues_rot(sunLightAtAN,[0 -1 0],i/180*pi);
  vectarrow([0 0 0],sunLight);hold on;axis equal; axis([-1 1 -1 1 -1 1]);
  pause(0.4);
  i
end
