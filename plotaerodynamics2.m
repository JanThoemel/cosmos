%plotaerodynamics2
close all;clc;clear all;

%totalforcevector=zeros(3,size(gamma,2),size(beta,2),size(alpha,2));

solarconstant=1;
sunlight=solarconstant*[0 1 0]';
rho=2; v=1;
wind=rho/2*v^2*[-1 0 0]';
noxpanels=0;noypanels=1;nozpanels=0;
controlvector=[1 0 0]';

alpha=0:10:360; %% yaw
beta=0:1:360; %% pitch 
gamma=0:10:360; %%roll

totalforcevector = totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma);

[alpha1,beta1,gamma1]=findBestAerodynamicAngles(totalforcevector,controlvector,alpha,beta,gamma);

alpha1
beta1
gamma1

function totalforcevector=totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma)
      
    Rx90=[0 1 0 ;-1 0  0 ; 0 0 1];
    Rz90=[1 0 0 ; 0 0 -1 ; 0 1 0];
    Iy = [0 1 0]';
    Ix=Rz90*Iy;
    Iz=Rx90*Iy;

    p1  = [1,0,1];p2 = [-1,0,1];p3 = [-1,0,-1];p4 = [1,0,-1];
    p12 = [0.33,0,0.33];p22 = [-0.33,0,0.33];p32 = [-0.33,0,-0.33];p42 = [0.33,0,-0.33];
    p13 = [0.66,0,0.66];p23 = [-0.66,0,0.66];p33 = [-0.66,0,-0.66];p43 = [0.66,0,-0.66];
    drag=zeros(size(alpha,2),1);
    lift=zeros(size(alpha,2),1);
    totalforcevector=zeros(3,size(gamma,2),size(beta,2),size(alpha,2));
    
    thetaaero=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    phiaero=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    thetasun=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    phisun=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    drag=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    lift=zeros(size(gamma,2),size(beta,2),size(alpha,2));

    %for k=1:size(gamma,2) %% yaw
      %for j=1:size(beta,2) %% pitch
        for i=1:size(alpha,2) %% roll
            k=1; %% yaw
            j=1; %% pitch
            %i=5;%% roll
                %% for reference y panel
                Rz2=[cosd(alpha(k)) -sind(alpha(k)) 0; sind(alpha(k)) cosd(alpha(k)) 0; 0 0 1]; %% yaw
                Ry =[cosd(beta(j))  0 sind(beta(j))  ; 0 1 0                          ; -sind(beta(j)) 0 cosd(beta(j))]; %% pitch
                Rz =[cosd(gamma(i)) -sind(gamma(i)) 0; sind(gamma(i)) cosd(gamma(i)) 0; 0 0 1]; %%roll

                pg = [(Rz2*Ry*Rz*p1')' ; (Rz2*Ry*Rz*p2')' ; (Rz2*Ry*Rz*p3')' ; (Rz2*Ry*Rz*p4')' ; (Rz2*Ry*Rz*p1')'];
                pg2 = [(Rz2*Ry*Rz*p12')' ; (Rz2*Ry*Rz*p22')' ; (Rz2*Ry*Rz*p32')' ; (Rz2*Ry*Rz*p42')' ; (Rz2*Ry*Rz*p12')'];
                pg3 = [(Rz2*Ry*Rz*p13')' ; (Rz2*Ry*Rz*p23')' ; (Rz2*Ry*Rz*p33')' ; (Rz2*Ry*Rz*p43')' ; (Rz2*Ry*Rz*p13')'];
                Ig=Rz2*Ry*Rz*Iy;

                [thetaaero(i,j,k),phiaero(i,j,k)]=thetaphi(wind, Ig);
                [drag(i,j,k),lift(i,j,k)]=aerodraglift(thetaaero(i,j,k),phiaero(i,j,k));
                aeroforcevector=[drag(i,j,k)  lift(i,j,k)*cosd( phiaero(i,j,k) )*sign(Ig(2))* sign(Ig(1))  lift(i,j,k) * sind( phiaero(i,j,k) ) * sign( Ig(3) )*sign( Ig(1) ) ]';

                [thetasun(i,j,k),phisun(i,j,k)]=thetaphi(sunlight,Ig);
                [dragsun(i,j,k),liftsun(i,j,k)]=sundraglift(thetasun(i,j,k),phisun(i,j,k));
                sunforcevector=[dragsun(i,j,k)  liftsun(i,j,k)*cosd( phisun(i,j,k) )*sign(Ig(2))* sign(Ig(1))  lift(i,j,k) * sind( phisun(i,j,k) ) * sign( Ig(3) )*sign( Ig(1) ) ]';
                %if thetasun(i,j,k)<90
                %     sunforcevector=norm(sunlight)*Ig;
                %elseif thetasun(i,j,k)>=90
                %     sunforcevector=-norm(sunlight)*Ig;
                %else
                %    fprintf('error sun force vector');
                %end

                totalforcevector(:,i,j,k)=aeroforcevector+sunforcevector;
            %% draw
            vectarrow(Ig);hold on;text(Ig(1),Ig(2),Ig(3),"normal",'HorizontalAlignment','left','FontSize',6);
            line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));
            vectarrow(sunlight);hold on;text(sunlight(1),sunlight(2),sunlight(3),"sunlight",'HorizontalAlignment','left','FontSize',6);
            vectarrow(wind);hold on;text(wind(1),wind(2),wind(3),"wind",'HorizontalAlignment','left','FontSize',6);
            vectarrow(aeroforcevector);hold on;text(aeroforcevector(1),aeroforcevector(2),aeroforcevector(3),"aeroforce",'HorizontalAlignment','left','FontSize',6);
            vectarrow(sunforcevector);hold on;text(sunforcevector(1),sunforcevector(2),sunforcevector(3),"sunforce",'HorizontalAlignment','left','FontSize',6);            
            %vectarrow(totalforcevector(:,i,j,k))
            axis equal;hold off;axis([-1 1 -1 1 -1 1])
            pause(1)
        end
      %end
    %end

    thetaaero=squeeze(thetaaero(:,j,k));
    phiaero=squeeze(phiaero(:,j,k));

    drag=squeeze(drag(:,j,k));
    lift=squeeze(lift(:,j,k));

    figure
        subplot(2,1,1)
        plot(alpha,thetaaero,alpha,phiaero);
        legend('thetaaero','phiaero');grid on;
        subplot(2,1,2)
        plot(alpha,drag,alpha,lift);
        legend('drag','lift');grid on;

end
function [theta phi]=thetaphi(refvec, vec)
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));
  if theta>180
      theta=180-theta;
  end
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)) );
end

function [drag lift]=aerodraglift(theta,phi)
  drag=-abs(sind(theta-90)); %%simplified formula
  lift=-abs(sind(2*(theta-90)));%*sign(theta); %% simplified formula     
end

function [drag lift]=sundraglift(theta,phi)
  drag=-abs(sind(theta-90)); %%simplified formula
  lift=-abs(sind(2*(theta-90)));%*sign(theta); %% simplified formula     
end


function [alpha1,beta1,gamma1]=findBestAerodynamicAngles(totalforcevector,controlvector,alpha,beta,gamma) 
    theta=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    for k=1:size(gamma,2) %% yaw
      for j=1:size(beta,2) %% pitch
        for i=1:size(alpha,2) %% roll
            [theta(i,j,k),~]=thetaphi(totalforcevector(:,i,j,k),controlvector);
        end
      end
    end
    %! find indizes of smallest theta    
    i=1;j=1;k=1;
    alpha1=alpha(i,j,k);
    beta1=beta(i,j,k);
    gamma1=gamma(i,j,k);
end