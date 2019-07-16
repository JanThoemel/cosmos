%plotaerodynamics2
close all;clc;clear all;
oldpath = path; path(oldpath,'..\matlabfunctions\')
altitude =340000; %% in m

[rho,v]=orbitalproperties(altitude);

wind=[-v 0 0]' ;
solarconstant=4.5e-6; %% pascal

sunlight=[1 1 1]'; 

panelsurface=0.01;
aeroscalingfactor=1;
sunscalingfactor=100;
noypanels=1;noxpanels=1;nozpanels=1;

controlvector=[1 0 0]';

alpha=0:3:360; %% yaw
beta=0:3:360; %% pitch 
gamma=0:3:360; %%roll


totalforcevector = totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma,rho,panelsurface,aeroscalingfactor,solarconstant,sunscalingfactor);

[alpha1,beta1,gamma1]=findBestAerodynamicAngles(totalforcevector,controlvector,alpha,beta,gamma);

%alpha1
%beta1
%gamma1

function totalforcevector=totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma,rho,panelsurface,aeroscalingfactor,solarconstant,sunscalingfactor)
    rotspeed=10;
    %% must be in dimension of force, i.e. N
    axislength=1.1*max([rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface*aeroscalingfactor solarconstant*sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)*panelsurface*sunscalingfactor]);
    Ry90=roty(90);
    Rz90=rotz(90);
    Rx90=rotx(90);
    Iz = [0 0 0.66*axislength]';
    Ix=Ry90*Iz;
    Iy=Rx90*Iz;
    %% zpanel
    pz1  = [axislength*0.9,axislength*0.9,0];pz2 = [axislength*0.9,-axislength*0.9,0];pz3 = [-axislength*0.9,-axislength*0.9,0];pz4 = [-axislength*0.9,axislength*0.9,0];
    pz12 = 0.33*pz1; pz22 = 0.33*pz2; pz32 = 0.33*pz3; pz42 = 0.33*pz4;
    pz13 = 0.66*pz1; pz23 = 0.66*pz2; pz33 = 0.66*pz3; pz43 = 0.66*pz4;
    
    drag=zeros(size(alpha,2),1);
    lift=zeros(size(alpha,2),1);
    totalforcevector=zeros(3,size(gamma,2),size(beta,2),size(alpha,2));
    
    thetaaero=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    phiaero=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    thetasun=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    phisun=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    drag=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    lift=zeros(size(gamma,2),size(beta,2),size(alpha,2));

    for k=1:size(gamma,2) %% yaw
      %for j=1:size(beta,2) %% pitch
        %for i=1:size(alpha,2) %% roll
            %k=1; %% yaw
            j=3; %% pitch
            i=1;%% roll
                %% rotation matrizes
                Rz2=[cosd(alpha(k)) -sind(alpha(k)) 0; sind(alpha(k)) cosd(alpha(k)) 0; 0 0 1]; %% yaw
                Ry =[cosd(beta(j))  0 sind(beta(j))  ; 0 1 0                          ; -sind(beta(j)) 0 cosd(beta(j))]; %% pitch
                Rz =[cosd(gamma(i)) -sind(gamma(i)) 0; sind(gamma(i)) cosd(gamma(i)) 0; 0 0 1]; %%roll
                
                if nozpanels %% zpanel
                    Ig=Rz2*Ry*Rz*Iz;
                    if norm(wind)
                        [thetaaero(i,j,k),phiaero(i,j,k),Ig2]=thetaphi(wind, Ig);
                        [drag(i,j,k),lift(i,j,k)]=aerodragliftcoef(thetaaero(i,j,k),phiaero(i,j,k));
                        aeroforcevector=-wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)           *    drag(i,j,k)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface;
                        ax=cross(wind,Ig2);liftvector = rodrigues_rot(wind,ax,90/180*pi);
                        aeroforcevector=-liftvector/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)      *    lift(i,j,k)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface + aeroforcevector;
                        aeroforcevector=aeroforcevector*aeroscalingfactor;
                        vectarrow([0 0 0],aeroforcevector);hold on;text(aeroforcevector(1),aeroforcevector(2),aeroforcevector(3),"aeroforce",'HorizontalAlignment','left','FontSize',6);
                        axis([-1.55 1.55 -1.55 1.55 -1.55 1.55]);
                    end 
                    if norm(sunlight)
                        [thetasun(i,j,k),phisun(i,j,k),Ig2]=thetaphi(sunlight,Ig);
                        [dragsun(i,j,k),liftsun(i,j,k)]=sundragliftcoef(thetasun(i,j,k),phisun(i,j,k));                                                               
                        sunforcevector=-sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)   * dragsun(i,j,k)*solarconstant*panelsurface;
                        ax=cross(sunlight,Ig2) ;liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
                        sunforcevector=-liftvector/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)  *liftsun(i,j,k)*solarconstant*panelsurface  +  sunforcevector;
                        sunforcevector=sunforcevector*sunscalingfactor;
                        vectarrow([0 0 0],sunforcevector);hold on;text(sunforcevector(1),sunforcevector(2),sunforcevector(3),"sunforce",'HorizontalAlignment','left','FontSize',6);
                        axis([-axislength axislength -axislength axislength -axislength axislength]);
                    end            
%                   totalforcevectorz(:,i,j,k)=nozpanels*(aeroforcevector+sunforcevector);
                    %% panel and normal
                    pg = [(Rz2*Ry*Rz*pz1')' ; (Rz2*Ry*Rz*pz2')' ; (Rz2*Ry*Rz*pz3')' ; (Rz2*Ry*Rz*pz4')' ; (Rz2*Ry*Rz*pz1')'];
                    pg2 = [(Rz2*Ry*Rz*pz12')' ; (Rz2*Ry*Rz*pz22')' ; (Rz2*Ry*Rz*pz32')' ; (Rz2*Ry*Rz*pz42')' ; (Rz2*Ry*Rz*pz12')'];
                    pg3 = [(Rz2*Ry*Rz*pz13')' ; (Rz2*Ry*Rz*pz23')' ; (Rz2*Ry*Rz*pz33')' ; (Rz2*Ry*Rz*pz43')' ; (Rz2*Ry*Rz*pz13')'];
                    vectarrow([0 0 0],Ig);hold on;text(Ig(1),Ig(2),Ig(3),"normal",'HorizontalAlignment','left','FontSize',6);
                    line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));hold on;
                    %vectarrow([0 0 0],totalforcevectorz(:,i,j,k));hold on;
                end
                if noxpanels %% xpanel
                    Igx=Rz2*Ry*Rz*Ix;
                    if norm(wind)
                        [thetaaero(i,j,k),phiaero(i,j,k),Igx2]=thetaphi(wind, Igx);
                        [drag(i,j,k),lift(i,j,k)]=aerodragliftcoef(thetaaero(i,j,k),phiaero(i,j,k));
                        aeroforcevectorx=-wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)                         *  drag(i,j,k)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface;
                        ax=cross(wind,Igx2);liftvector = rodrigues_rot(wind,ax,90/180*pi);
                        aeroforcevectorx=-liftvector/sqrt(liftvector(1)^2+liftvector(2)^2+liftvector(3)^2) *  lift(i,j,k)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface + aeroforcevectorx;
                        aeroforcevectorx=aeroforcevectorx*aeroscalingfactor;
                        vectarrow([0 0 0],aeroforcevectorx);hold on;text(aeroforcevectorx(1),aeroforcevectorx(2),aeroforcevectorx(3),"aeroforcex",'HorizontalAlignment','left','FontSize',6);
                        axis([-axislength axislength -axislength axislength -axislength axislength]);
                    end
                    if norm(sunlight)
                        [thetasun(i,j,k),phisun(i,j,k),Igx2]=thetaphi(sunlight,Igx);
                        [dragsun(i,j,k),liftsun(i,j,k)]=sundragliftcoef(thetasun(i,j,k),phisun(i,j,k));
                        sunforcevectorx=-sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)           *  dragsun(i,j,k)*solarconstant*panelsurface;
                        ax=cross(sunlight,Igx2); liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
                        sunforcevectorx=-liftvector/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)         *  liftsun(i,j,k)*solarconstant*panelsurface  +  sunforcevectorx;
                        sunforcevectorx=sunforcevectorx*sunscalingfactor;
                        vectarrow([0 0 0],sunforcevectorx);hold on;text(sunforcevectorx(1),sunforcevectorx(2),sunforcevectorx(3),"sunforcex",'HorizontalAlignment','left','FontSize',6);
                        axis([-axislength axislength -axislength axislength -axislength axislength]);
                    end
%                   totalforcevectorx(:,i,j,k)=noxpanels*(aeroforcevector+sunforcevector);
                    %% panel and normal
                    pgx = [(Rz2*Ry*Rz*-Ry90*pz1')' ; (Rz2*Ry*Rz*-Ry90*pz2')' ; (Rz2*Ry*Rz*-Ry90*pz3')' ; (Rz2*Ry*Rz*-Ry90*pz4')' ; (Rz2*Ry*Rz*-Ry90*pz1')'];
                    pgx2 = [(Rz2*Ry*Rz*-Ry90*pz12')' ; (Rz2*Ry*Rz*-Ry90*pz22')' ; (Rz2*Ry*Rz*-Ry90*pz32')' ; (Rz2*Ry*Rz*-Ry90*pz42')' ; (Rz2*Ry*Rz*-Ry90*pz12')'];
                    pgx3 = [(Rz2*Ry*Rz*-Ry90*pz13')' ; (Rz2*Ry*Rz*-Ry90*pz23')' ; (Rz2*Ry*Rz*-Ry90*pz33')' ; (Rz2*Ry*Rz*-Ry90*pz43')' ; (Rz2*Ry*Rz*-Ry90*pz13')'];
                 
                    vectarrow([0 0 0],Igx);hold on;text(Igx(1),Igx(2),Igx(3),"normalx",'HorizontalAlignment','left','FontSize',6);hold on;
                    line(pgx(:,1), pgx(:,2), pgx(:,3));line(pgx2(:,1), pgx2(:,2), pgx2(:,3));line(pgx3(:,1), pgx3(:,2), pgx3(:,3));hold on;
                    %vectarrow([0 0 0],totalforcevector(:,i,j,k)); hold on;
                end
                if noypanels %% ypanel
                    Igy=Rz2*Ry*Rz*Iy;
                    if norm(wind)
                        [thetaaero(i,j,k),phiaero(i,j,k),Igy2]=thetaphi(wind, Igy);
                        [drag(i,j,k),lift(i,j,k)]=aerodragliftcoef(thetaaero(i,j,k),phiaero(i,j,k));             
                        aeroforcevectory=-wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)        *  drag(i,j,k)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface;
                        ax=cross(wind,Igy2) ;liftvector = rodrigues_rot(wind,ax,90/180*pi);
                        aeroforcevectory=-liftvector/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)  *  lift(i,j,k)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface + aeroforcevectory;
                        aeroforcevectory=aeroforcevectory*aeroscalingfactor;
                        vectarrow([0 0 0],aeroforcevectory);hold on;text(aeroforcevectory(1),aeroforcevectory(2),aeroforcevectory(3),"aeroforcey",'HorizontalAlignment','left','FontSize',6);
                        axis([-axislength axislength -axislength axislength -axislength axislength]);
                    end
                    if norm(sunlight)
                        [thetasun(i,j,k),phisun(i,j,k),Igy2]=thetaphi(sunlight,Igy);
                        [dragsun(i,j,k),liftsun(i,j,k)]=sundragliftcoef(thetasun(i,j,k),phisun(i,j,k));                                       
                        sunforcevectory=-sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)     *   dragsun(i,j,k)*solarconstant*panelsurface;
                        ax=cross(sunlight,Igy2); liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
                        sunforcevectory=-liftvector/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)   *    liftsun(i,j,k)*solarconstant*panelsurface+sunforcevectory;
                        sunforcevectory=sunforcevectory*sunscalingfactor;
                        vectarrow([0 0 0],sunforcevectory);hold on;text(sunforcevectory(1),sunforcevectory(2),sunforcevectory(3),"sunforcey",'HorizontalAlignment','left','FontSize',6); 
                        axis([-axislength axislength -axislength axislength -axislength axislength]);
                    end
%                   totalforcevectory(:,i,j,k)=noypanels*(aeroforcevector+sunforcevector);
                    %% panel and normal
                    pgy = [(Rz2*Ry*Rz*-Rx90*pz1')' ; (Rz2*Ry*Rz*-Rx90*pz2')' ; (Rz2*Ry*Rz*-Rx90*pz3')' ; (Rz2*Ry*Rz*-Rx90*pz4')' ; (Rz2*Ry*Rz*-Rx90*pz1')'];
                    pgy2 = [(Rz2*Ry*Rz*-Rx90*pz12')' ; (Rz2*Ry*Rz*-Rx90*pz22')' ; (Rz2*Ry*Rz*-Rx90*pz32')' ; (Rz2*Ry*Rz*-Rx90*pz42')' ; (Rz2*Ry*Rz*-Rx90*pz12')'];
                    pgy3 = [(Rz2*Ry*Rz*-Rx90*pz13')' ; (Rz2*Ry*Rz*-Rx90*pz23')' ; (Rz2*Ry*Rz*-Rx90*pz33')' ; (Rz2*Ry*Rz*-Rx90*pz43')' ; (Rz2*Ry*Rz*-Rx90*pz13')'];
                    vectarrow([0 0 0],Igy);hold on;text(Igy(1),Igy(2),Igy(3),"normaly",'HorizontalAlignment','left','FontSize',6);
                    line(pgy(:,1), pgy(:,2), pgy(:,3));line(pgy2(:,1), pgy2(:,2), pgy2(:,3));line(pgy3(:,1), pgy3(:,2), pgy3(:,3));hold on;
                    %vectarrow([0 0 0],totalforcevector(:,i,j,k));hold on;                  
                end
                if norm(wind)
                 windscalingfactor=1/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface*aeroscalingfactor;
                 vectarrow([0 0 0],wind*windscalingfactor);hold on;text(wind(1)*windscalingfactor,wind(2)*windscalingfactor,wind(3)*windscalingfactor,"wind",'HorizontalAlignment','left','FontSize',6);
                end
                if norm(sunlight)
                 sunlightscalingfactor=1/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)*rho/2*(wind(1)^2+wind(2)^2+wind(3)^2)*panelsurface*sunscalingfactor;
                 vectarrow([0 0 0],sunlight*sunlightscalingfactor);hold on;text(sunlight(1)*sunlightscalingfactor,sunlight(2)*sunlightscalingfactor,sunlight(3)*sunlightscalingfactor,"sunlight",'HorizontalAlignment','left','FontSize',6);
                end
                %text(0,0,0,strcat('scalingfactor: ',int2str(scalingfactor)),'HorizontalAlignment','left','FontSize',6);
                title(strcat('aeroscalingfactor: ',int2str(aeroscalingfactor),'  suncalingfactor: ',int2str(sunscalingfactor)));
                axis([-axislength axislength -axislength axislength -axislength axislength]);
                hold off;
                pause(1/rotspeed)
        end
      %end
    %end

    thetaaero=squeeze(thetaaero(i,:,k));
    phiaero=squeeze(phiaero(i,:,k));

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
function [theta,phi,Ig2]=thetaphi(refvec, vec)
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));
  if theta>90
      theta=180-theta;
      otherside=-1;
      Ig2=-vec;
  else
      otherside=1;
      Ig2=vec;
  end
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)) );
end

function [drag lift]=aerodragliftcoef(theta,phi)
  drag=-abs(1.2*sind(theta-90))+0.6; %%simplified formula
  lift=-abs(0.12*sind(2*(theta-90)));%% simplified formula     
end

function [drag lift]=sundragliftcoef(theta,phi)
  drag=-abs(sind(theta-90)); %%simplified formula
  lift=-abs(sind(2*(theta-90)));%% simplified formula     
end


function [alpha1,beta1,gamma1]=findBestAerodynamicAngles(totalforcevector,controlvector,alpha,beta,gamma) 
    theta=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    for k=1:size(gamma,2) %% yaw
      for j=1:size(beta,2) %% pitch
        for i=1:size(alpha,2) %% roll
            [theta(i,j,k),phi]=thetaphi(totalforcevector(:,i,j,k),controlvector);
        end
      end
    end
    [theta(i,j,k),phi]=thetaphi(totalforcevector(:,i,j,k),controlvector);
    %! find indizes of smallest theta    
    i=1;j=1;k=1;
    alpha1=alpha(i,j,k);
    beta1=beta(i,j,k);
    gamma1=gamma(i,j,k);
end





