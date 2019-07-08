%plotaerodynamics2
close all;clc;clear all;

Rx90=[0 1 0 ;-1 0  0 ; 0 0 1];
Rz90=[1 0 0 ; 0 0 -1 ; 0 1 0];
Iy = [0 1 0]';
Ix=Rz90*Iy;
Iz=Rx90*Iy;
p1 = [1,0,1];p2 = [-1,0,1];p3 = [-1,0,-1];p4 = [1,0,-1];
alpha=0:10:360;
beta=0:10:360;
gamma=0:1:360;
drag=zeros(size(alpha,2),1);
lift=zeros(size(alpha,2),1);
theta=zeros(size(alpha,2),size(beta,2),size(gamma,2));
phi=zeros(size(alpha,2),size(beta,2),size(gamma,2));
drag=zeros(size(alpha,2),size(beta,2),size(gamma,2));
lift=zeros(size(alpha,2),size(beta,2),size(gamma,2));

for k=1:size(gamma,2)
  %for j=1:size(beta,2)
    %for i=1:size(alpha,2)
        i=3;%% roll
        j=3;%% pitch
        Rz=[cosd(alpha(i)) -sind(alpha(i)) 0; sind(alpha(i)) cosd(alpha(i)) 0; 0 0 1]; %%roll
        Ry=[cosd(beta(j)) 0 sind(beta(j)); 0 1 0; -sind(beta(j)) 0 cosd(beta(j))]; %% pitch
        Rz2=[cosd(gamma(k)) -sind(gamma(k)) 0; sind(gamma(k)) cosd(gamma(k)) 0; 0 0 1]; %% yaw

        pg = [(Rz2*Ry*Rz*p1')' ; (Rz2*Ry*Rz*p2')' ; (Rz2*Ry*Rz*p3')' ; (Rz2*Ry*Rz*p4')' ; (Rz2*Ry*Rz*p1')'];
        Ig=Rz2*Ry*Rz*Iy;

        [theta(i,j,k) phi(i,j,k)]=thepi(Ig);
        [drag(i,j,k) lift(i,j,k)]=drali(theta(i,j,k),phi(i,j,k));
        delta=atand(drag(i,j,k)/lift(i,j,k))
        Raero=[cosd(delta) -sind(delta) 0; sind(delta) cosd(delta) 0; 0 0 1]; 
        liftvector=Rz2*Ry*Rz*Raero*Iy;
        %dragvector
        %% draw
        vectarrow(Ig)
        hold on; axis equal;legend;
        line(pg(:,1), pg(:,2), pg(:,3));
        hold on;
        vectarrow(liftvector)
        hold off;legend;
        pause(0.0001)
    end
  %end
%end
theta=squeeze(theta(i,j,:))+5
phi=squeeze(phi(i,j,:));

drag=squeeze(drag(i,j,:));
lift=squeeze(lift(i,j,:));

figure
    subplot(2,1,1)
    plot(gamma,theta,gamma,phi);
    legend;
    subplot(2,1,2)
    plot(gamma,drag,gamma,lift);
    legend;

function [theta phi]=thepi(vec)
  %theta=atand(V(3)/V(1));
  %phi=atand( V(2) / sqrt(V(1)^2+V(3)^2));
  theta=atand(sqrt(vec(3)^2+vec(2))^2/vec(1));
  phi=atand( vec(2) / vec(1));
end

function [drag lift]=drali(theta,phi)
  drag=sind(theta);
  lift=sind(2*theta);
end