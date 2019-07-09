%plotaerodynamics2
close all;clc;clear all;

Rx90=[0 1 0 ;-1 0  0 ; 0 0 1];
Rz90=[1 0 0 ; 0 0 -1 ; 0 1 0];
Iy = [0 1 0]';
Ix=Rz90*Iy;
Iz=Rx90*Iy;
p1 = [1,0,1];p2 = [-1,0,1];p3 = [-1,0,-1];p4 = [1,0,-1];
gamma=0:10:360;
beta=0:5:360;
alpha=0:10:360;
drag=zeros(size(alpha,2),1);
lift=zeros(size(alpha,2),1);
theta=zeros(size(alpha,2),size(beta,2),size(gamma,2));
phi=zeros(size(alpha,2),size(beta,2),size(gamma,2));
drag=zeros(size(alpha,2),size(beta,2),size(gamma,2));
lift=zeros(size(alpha,2),size(beta,2),size(gamma,2));
delta=zeros(size(alpha,2),size(beta,2),size(gamma,2));
%for k=1:size(gamma,2)
  for j=1:size(beta,2)
    %for i=1:size(alpha,2)
        k=1; %% yaw
        %j=1;%% pitch
        i=5;%% roll
        Rz2=[cosd(gamma(k)) -sind(gamma(k)) 0; sind(gamma(k)) cosd(gamma(k)) 0; 0 0 1]; %% yaw
        Ry=[cosd(beta(j)) 0 sind(beta(j)); 0 1 0; -sind(beta(j)) 0 cosd(beta(j))]; %% pitch
        Rz=[cosd(alpha(i)) -sind(alpha(i)) 0; sind(alpha(i)) cosd(alpha(i)) 0; 0 0 1]; %%roll

        pg = [(Rz2*Ry*Rz*p1')' ; (Rz2*Ry*Rz*p2')' ; (Rz2*Ry*Rz*p3')' ; (Rz2*Ry*Rz*p4')' ; (Rz2*Ry*Rz*p1')'];
        Ig=Rz2*Ry*Rz*Iy;

        [theta(i,j,k) phi(i,j,k)]=thepi(Ig);
        [drag(i,j,k) lift(i,j,k)]=drali(theta(i,j,k),phi(i,j,k));
        %delta(i,j,k)=atand(drag(i,j,k)/(lift(i,j,k)+lift(i,j,k)*0.001))-90;
        %Raero=[cosd(delta(i,j,k)) -sind(delta(i,j,k)) 0; sind(delta(i,j,k)) cosd(delta(i,j,k)) 0; 0 0 1]; 
        %liftvector=Rz2*Ry*Rz*Raero*Iy;
        aeroforcevector=[drag(i,j,k)  lift(i,j,k) 0]';
        %dragvector
        %% draw
        vectarrow(Ig)
        hold on; axis equal;legend;
        line(pg(:,1), pg(:,2), pg(:,3));
        hold on;
        vectarrow(aeroforcevector)
        hold off;legend;
        pause(0.0001)
    end
  %end
%end
theta=squeeze(theta(i,j,:))+5
phi=squeeze(phi(i,j,:));
delta=squeeze(delta(i,j,:));

drag=squeeze(drag(i,j,:));
lift=squeeze(lift(i,j,:));

figure
    subplot(2,1,1)
    plot(gamma,theta,gamma,phi,gamma,delta);
    legend;grid on;
    subplot(2,1,2)
    plot(gamma,drag,gamma,lift);
    legend;grid on;

function [theta phi]=thepi(vec)
  theta=atand(sqrt(vec(3)^2+vec(2)^2)/abs(vec(1)));
  phi=atand( vec(2) / vec(1));
end

function [drag lift]=drali(theta,phi)
  drag=-abs(sind(theta-90));
  lift=-sind(2*(theta))*sign(phi);      
end