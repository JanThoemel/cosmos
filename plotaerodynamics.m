close all; clear all;clc;

%% panels with normal towards [x y z]
panels=[0 1 0];

%% elevation and azimuth
theta=0:5:90;
phi=0:5:360;

f1=bigaerodynamics(theta, phi,panels);
[pitch yaw]=pitchyaw(theta, phi);
    
%% plot
figure;
        subplot(1,3,1)
        mesh(theta, phi, squeeze(f1(1,:,:))');
        grid on;xlabel('\theta_{org}');ylabel('\phi_{org}');zlabel('fx_{org}');
        subplot(1,3,2)
        mesh(theta, phi, squeeze(f1(2,:,:))');
        grid on;xlabel('\theta_{org}');ylabel('\phi_{org}');zlabel('fx_{org}');
        subplot(1,3,3)
        mesh(theta, phi, squeeze(f1(3,:,:))');
        grid on;xlabel('\theta_{org}');ylabel('\phi_{org}');zlabel('fx_{org}');

        %% plot
        size(theta)
        size(phi)
        size(phi(36:73))
        size(pitch)
figure;
        subplot(1,2,1)
        mesh(theta, phi, pitch');
        grid on;xlabel('\theta_{org}');ylabel('\phi_{org}');zlabel('pitch');axis([0 95 0 365 0 185]);
        subplot(1,2,2)
        mesh(theta, phi, yaw');
        grid on;xlabel('\theta_{org}');ylabel('\phi_{org}');zlabel('yaw');axis([0 95 0 365 0 185]);



%% compute aerodynamics of 2
%f2=[mean(f1(1,:)) 0 0]; 
f2=[0 0 0]; 
%% compute difference
df=zeros(size(f1));
df(1,:)=f1(1,:)-f2(1);
df(2,:)=f1(2,:)-f2(2);
df(3,:)=f1(3,:)-f2(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pitch yaw]=pitchyaw(theta, phi)
    for i=1:size(theta,2)
        for j=1:size(phi,2)

            if phi(j)>=0 && phi(j)<=180
                pitch(i,j)=theta(i)*sind(phi(j));
            elseif phi(j)>180 && phi(j)<=360
                pitch(i,j)=180+theta(i)*sind(phi(j));
            else
                printf('error 5001');
            end
            
            if phi(j)>=0 && phi(j)<=90
                yaw(i,j)=theta(i)*cosd(phi(j));
            elseif phi(j)>90 && phi(j)<=180
                yaw(i,j)=180+theta(i)*cosd(phi(j));         
            elseif phi(j)>180 && phi(j)<=270 
                yaw(i,j)=180+theta(i)*cosd(phi(j));
            elseif phi(j)>270 && phi(j)<=360
                yaw(i,j)=theta(i)*cosd(phi(j));         
            end
            roll(j)=phi(j);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f1=bigaerodynamics(theta, phi, panels)
    for i=1:size(theta,2)
        for j=1:size(phi,2)
            f=zeros(3,1);       
            panelxtheta=theta(i)+90;
            panelxphi=phi(j)+90;
            if panelxtheta>90
                panelxtheta=90-panelxtheta;
                panelxphi=panelxphi+180;
            end
            if panelxphi>360
                panelxphi=panelxphi-360;
            end
            panelytheta=theta(i);
            panelyphi=phi(j);
            panelztheta=theta(i);
            panelzphi=phi(j)+90;
            if panelzphi>360
                panelzphi=panelxphi-360;
            end
            %% sum x-panels
            for nb_xp=1:panels(1) 
                %% compute aerodynamics of 1
                f=f+aerodynamics(panelxtheta,panelxphi)';
            end
            %% sum y-panels
            for nb_yp=1:panels(2)
                %% compute aerodynamics of 2
                f=f+aerodynamics(panelytheta,panelyphi)';
            end
            %% sum z-panels
            for nb_zp=1:panels(3) 
                %% compute aerodynamics of 3
                f=f+aerodynamics(panelztheta,panelzphi)';
            end
            f1(:,i,j)=f;
         end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=aerodynamics(theta, phi)
    eta=0.1;    
    eps=0.1;
    
    p=-2*eps*(sind(theta))^3+eta*(eps-1)*(sind(theta))^2+(eps-1)*sind(theta);
    g=-cosd(theta)*sind(theta)*(eta-eps*eta+2*eps*sind(theta));
    f(1)=p;
    f(2)=g*cosd(phi);
    f(3)=g*sind(phi);
end    