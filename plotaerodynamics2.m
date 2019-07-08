%plotaerodynamics2
close all;clc;clear all;

Rx90=[0 1 0 ;-1 0  0 ; 0 0 1];
Rz90=[1 0 0 ; 0 0 -1 ; 0 1 0];
Iy = [0 1 0]';
Ix=Rz90*Iy;
Iz=Rx90*Iy;
p1 = [1,0,1];p2 = [-1,0,1];p3 = [-1,0,-1];p4 = [1,0,-1];
alpha=0:1:360;
beta=0:1:360;
gamma=0:1:360;
drag=zeros(size(alpha,2),1)
lift=zeros(size(alpha,2),1)
j=1;k=1;
for i=1:size(gamma,2)
    Rz=[cosd(alpha(i)) -sind(alpha(i)) 0; sind(alpha(i)) cosd(alpha(i)) 0; 0 0 1]; %%roll
    Ry=[cosd(beta(j)) 0 sind(beta(j)); 0 1 0; -sind(beta(j)) 0 cosd(beta(j))]; %% yaw
    Rz2=[cosd(gamma(k)) -sind(gamma(k)) 0; sind(gamma(k)) cosd(gamma(k)) 0; 0 0 1]; %% pitch

    pg = [(Rz2*Ry*Rz*p1')' ; (Rz2*Ry*Rz*p2')' ; (Rz2*Ry*Rz*p3')' ; (Rz2*Ry*Rz*p4')' ; (Rz2*Ry*Rz*p1')'];
    Ig=Rz2*Ry*Rz*Iy;

    [theta(i) phi(i)]=thepi(Ig)
    [drag(i) lift(i)]=drali(theta(i),phi(i))
    
    %% draw
    vectarrow(Ig)
    hold on; axis equal;
    line(pg(:,1), pg(:,2), pg(:,3));
    hold off;
    pause(0.001)
end

figure
    plot(alpha,theta,alpha,phi)
figure
    plot(alpha,drag,alpha,lift)

function [theta phi]=thepi(V)
  theta=atand(V(3)/V(1));
  phi=atand( V(2) / sqrt(V(1)^2+V(3)^2));
end

function [drag lift]=drali(theta,phi)
  drag=sin(theta);
  lift=sin(2*theta);
end