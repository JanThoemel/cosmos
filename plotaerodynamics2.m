%plotaerodynamics2
close all;clc;clear all;

Rx90=[0 1 0 ;-1 0  0 ; 0 0 1];
Rz90=[1 0 0 ; 0 0 -1 ; 0 1 0];
Iy = [0 1 0]';
Ix=Rz90*Iy;
Iz=Rx90*Iy;
sunlight=[0 1 0]';

p1 = [1,0,1];p2 = [-1,0,1];p3 = [-1,0,-1];p4 = [1,0,-1];
p12 = [0.33,0,0.33];p22 = [-0.33,0,0.33];p32 = [-0.33,0,-0.33];p42 = [0.33,0,-0.33];
p13 = [0.66,0,0.66];p23 = [-0.66,0,0.66];p33 = [-0.66,0,-0.66];p43 = [0.66,0,-0.66];
gamma=0:10:360;
beta=0:1:360;
alpha=0:10:360;
drag=zeros(size(alpha,2),1);
lift=zeros(size(alpha,2),1);
thetaaero=zeros(size(alpha,2),size(beta,2),size(gamma,2));
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
        pg2 = [(Rz2*Ry*Rz*p12')' ; (Rz2*Ry*Rz*p22')' ; (Rz2*Ry*Rz*p32')' ; (Rz2*Ry*Rz*p42')' ; (Rz2*Ry*Rz*p12')'];
        pg3 = [(Rz2*Ry*Rz*p13')' ; (Rz2*Ry*Rz*p23')' ; (Rz2*Ry*Rz*p33')' ; (Rz2*Ry*Rz*p43')' ; (Rz2*Ry*Rz*p13')'];
        Ig=Rz2*Ry*Rz*Iy;

        [thetaaero(i,j,k) phiaero(i,j,k)]=thepi(Ig);
        [drag(i,j,k) lift(i,j,k)]=drali(thetaaero(i,j,k),phiaero(i,j,k));
        aeroforcevector=[drag(i,j,k)  lift(i,j,k)*cosd( atand(Ig(3)/Ig(2)) )  lift(i,j,k)*sind( atand(Ig(3)/Ig(2)) )]';

        [thetasun(i,j,k) phisun(i,j,k)]=thepi(sunvector);
        [dragsun(i,j,k) liftsun(i,j,k)]=drali(thetaaero(i,j,k),phisun(i,j,k));
        sunlightforcevector=[0 0 0]';

        totalforcevector=aeroforcevector+sunlightforcevector;
        %% draw
        vectarrow(Ig)
        hold on; axis equal;legend;
        line(pg(:,1), pg(:,2), pg(:,3));
        line(pg2(:,1), pg2(:,2), pg2(:,3));
        line(pg3(:,1), pg3(:,2), pg3(:,3));
        hold on;axis equal;legend;
        vectarrow(totalforcevector)
        hold on;axis equal;legend;
        vectarrow(sunlight)
        hold off;axis equal;legend;
        pause(0.0001)
    end
  %end
%end
thetaaero=squeeze(thetaaero(i,j,:))+5
phi=squeeze(phi(i,j,:));
delta=squeeze(delta(i,j,:));

drag=squeeze(drag(i,j,:));
lift=squeeze(lift(i,j,:));

figure
    subplot(2,1,1)
    plot(gamma,thetaaero,gamma,phi,gamma,delta);
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