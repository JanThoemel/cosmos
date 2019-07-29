function totalforcevector=totalforcevectorfunction(wind,sunlight,noxpanels,noypanels,nozpanels,alpha,beta,gamma,panelsurface,aeroscalingfactor,solarpressure,sunscalingfactor, windpressure,deltaangle)
  rotspeed=10;draw=0;
  %% %% the possible forcevectors
  %% define filename convenction
  filename=strcat('tfv_panels',int2str(nozpanels),int2str(noxpanels),int2str(noypanels),'_wind',num2str(solarpressure,'%1.1e'),int2str(wind(1)),int2str(wind(2)),int2str(wind(3)),'_sun',num2str(windpressure,'%1.1e'),int2str(sunlight(1)),int2str(sunlight(2)),int2str(sunlight(3)),'_deltaangle',int2str(deltaangle),'.mat');
  %% does this file exists already
  if isfile(filename)
    fprintf('\nloading file');
    load(filename,'totalforcevector')
    fprintf('\ndone\n');
  else
    fprintf('\n computing forces');  
    %% must be in dimensions of force, i.e. N
    axislength=1.1*max([windpressure*panelsurface*aeroscalingfactor solarpressure*panelsurface*sunscalingfactor]);
    Ry90=roty(90);
    %Rz90=rotz(90);
    Rx90=rotx(90);
    Iz = [0 0 0.66*axislength]';
    Ix=Ry90*Iz;
    Iy=Rx90*Iz;
    %% zpanel
    pz1  = [axislength*0.9,axislength*0.9,0];pz2 = [axislength*0.9,-axislength*0.9,0];pz3 = [-axislength*0.9,-axislength*0.9,0];pz4 = [-axislength*0.9,axislength*0.9,0];
    pz12 = 0.33*pz1; pz22 = 0.33*pz2; pz32 = 0.33*pz3; pz42 = 0.33*pz4;
    pz13 = 0.66*pz1; pz23 = 0.66*pz2; pz33 = 0.66*pz3; pz43 = 0.66*pz4;
    thetaaero=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    phiaero=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    thetasun=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    phisun=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    aerodragcoef=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    aeroliftcoef=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    sundragcoef=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    sunliftcoef=zeros(size(gamma,2),size(beta,2),size(alpha,2));
    aeroforcevectorz=[0 0 0]';
    aeroforcevectorx=[0 0 0]';
    aeroforcevectory=[0 0 0]';
    sunforcevectorz=[0 0 0]';
    sunforcevectorx=[0 0 0]';
    sunforcevectory=[0 0 0]';
    totalforcevectorz=zeros(3,size(alpha,2),size(beta,2),size(gamma,2));
    totalforcevectorx=zeros(3,size(alpha,2),size(beta,2),size(gamma,2));
    totalforcevectory=zeros(3,size(alpha,2),size(beta,2),size(gamma,2));
    totalforcevector =zeros(3,size(alpha,2),size(beta,2),size(gamma,2));
    for k=1:size(gamma,2) %% yaw
      for j=1:size(beta,2) %% pitch
        for i=1:size(alpha,2) %% roll
          %k=1; %% yaw
          %j=3; %% pitch
          %i=1;%% roll
          %% rotation matrizes
          RzY =[cosd(gamma(k)) -sind(gamma(k)) 0; sind(gamma(k)) cosd(gamma(k)) 0; 0 0 1]; %%yaw
          Ry =[cosd(beta(j))  0 sind(beta(j))  ; 0 1 0                          ; -sind(beta(j)) 0 cosd(beta(j))]; %% pitch
          RzR=[cosd(alpha(i)) -sind(alpha(i)) 0; sind(alpha(i)) cosd(alpha(i)) 0; 0 0 1]; %% roll
          if nozpanels %% zpanel
            Igz=RzY*Ry*RzR*Iz;
            if norm(windpressure)
              [thetaaero(i,j,k),phiaero(i,j,k),Ig2]=thetaphi1(wind, Igz);
              [aerodragcoef(i,j,k),aeroliftcoef(i,j,k)]=aerodragliftcoef(thetaaero(i,j,k),phiaero(i,j,k));
              aeroforcevectorz=-wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)           *    aerodragcoef(i,j,k)*windpressure*panelsurface;
              ax=cross(wind,Ig2);
              if norm(ax)~=0
                liftvector = rodrigues_rot(wind,ax,90/180*pi);
              else
                liftvector = [0 0 0]';
              end
              aeroforcevectorz=-liftvector/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)      *    aeroliftcoef(i,j,k)*windpressure*panelsurface + aeroforcevectorz;
              aeroforcevectorz=aeroforcevectorz*aeroscalingfactor;
              %fprintf('\n %3.0f %3.0f %3.0f     %+1.1e  %+1.1e  %+1.1e',alpha(i),beta(j),gamma(k),aeroforcevectorz(1),aeroforcevectorz(2),aeroforcevectorz(3))
              %input('a')
            end 
            if norm(solarpressure)
              [thetasun(i,j,k),phisun(i,j,k),Ig2]=thetaphi1(sunlight,Igz);
              [sundragcoef(i,j,k), sunliftcoef(i,j,k) ]=sundragliftcoef(thetasun(i,j,k),phisun(i,j,k));                                                               
              sunforcevectorz=-sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)   * sundragcoef(i,j,k)*solarpressure*panelsurface;
              ax=cross(sunlight,Ig2) ;
              if norm(ax)~=0
                liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
              else
                liftvector = [0 0 0]';
              end
              sunforcevectorz=-liftvector/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)  *sunliftcoef(i,j,k)*solarpressure*panelsurface  +  sunforcevectorz;
              sunforcevectorz=sunforcevectorz*sunscalingfactor;
            end            
            totalforcevectorz(:,i,j,k)=nozpanels*(aeroforcevectorz+sunforcevectorz);
          end
          if noxpanels %% xpanel
            Igx=RzY*Ry*RzR*Ix;
            if norm(windpressure)
              [thetaaero(i,j,k),phiaero(i,j,k),Igx2]=thetaphi1(wind, Igx);
              [aerodragcoef(i,j,k),aeroliftcoef(i,j,k)]=aerodragliftcoef(thetaaero(i,j,k),phiaero(i,j,k));
              aeroforcevectorx=-wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)                         *  aerodragcoef(i,j,k)*windpressure*panelsurface;
              ax=cross(wind,Igx2);
              if norm(ax)~=0
                liftvector = rodrigues_rot(wind,ax,90/180*pi);
              else
                liftvector = [0 0 0]';;
              end
              aeroforcevectorx=-liftvector/sqrt(liftvector(1)^2+liftvector(2)^2+liftvector(3)^2) *  aeroliftcoef(i,j,k)*windpressure*panelsurface + aeroforcevectorx;
              aeroforcevectorx=aeroforcevectorx*aeroscalingfactor;
            end
            if norm(solarpressure)
              [thetasun(i,j,k),phisun(i,j,k),Igx2]=thetaphi1(sunlight,Igx);
              [sundragcoef(i,j,k),sunliftcoef(i,j,k)]=sundragliftcoef(thetasun(i,j,k),phisun(i,j,k));
              sunforcevectorx=-sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)           *  sundragcoef(i,j,k)*solarpressure*panelsurface;
              ax=cross(sunlight,Igx2);
              if norm(ax)~=0
                liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
              else
                liftvector = [0 0 0]';;
              end
              sunforcevectorx=-liftvector/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)         *  sunliftcoef(i,j,k)*solarpressure*panelsurface  +  sunforcevectorx;
              sunforcevectorx=sunforcevectorx*sunscalingfactor;
            end
            totalforcevectorx(:,i,j,k)=noxpanels*(aeroforcevectorx+sunforcevectorx);
          end
          if noypanels %% ypanel
            Igy=RzY*Ry*RzR*Iy;
            if norm(windpressure)
              [thetaaero(i,j,k),phiaero(i,j,k),Igy2]=thetaphi1(wind, Igy);
              [aerodragcoef(i,j,k),aeroliftcoef(i,j,k)]=aerodragliftcoef(thetaaero(i,j,k),phiaero(i,j,k));             
              aeroforcevectory=-wind/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)        *  aerodragcoef(i,j,k)*windpressure*panelsurface;
              ax=cross(wind,Igy2) ;
              if norm(ax)~=0
                liftvector = rodrigues_rot(wind,ax,90/180*pi);
              else
                liftvector = [0 0 0]';;
              end            
              aeroforcevectory=-liftvector/sqrt(wind(1)^2+wind(2)^2+wind(3)^2)  *  aeroliftcoef(i,j,k)*windpressure*panelsurface + aeroforcevectory;
              aeroforcevectory=aeroforcevectory*aeroscalingfactor;
            end
            if norm(solarpressure)
              [thetasun(i,j,k),phisun(i,j,k),Igy2]=thetaphi1(sunlight,Igy);
              [sundragcoef(i,j,k),sunliftcoef(i,j,k)]=sundragliftcoef(thetasun(i,j,k),phisun(i,j,k));                                       
              sunforcevectory=-sunlight/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)     *   sundragcoef(i,j,k)*solarpressure*panelsurface;
              ax=cross(sunlight,Igy2);
              if norm(ax)~=0
                liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
              else
                liftvector = [0 0 0]';;
              end                        
              sunforcevectory=-liftvector/sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2)   *    sunliftcoef(i,j,k)*solarpressure*panelsurface+sunforcevectory;
              sunforcevectory=sunforcevectory*sunscalingfactor;
            end
            totalforcevectory(:,i,j,k)=noypanels*(aeroforcevectory+sunforcevectory);
          end
          %%draw
          if draw
            if norm(windpressure)
              vectarrow([0 0 0],wind*windpressure*panelsurface*aeroscalingfactor);hold on;text(wind(1)*windpressure*panelsurface*aeroscalingfactor,wind(2)*windpressure*panelsurface*aeroscalingfactor,wind(3)*windpressure*panelsurface*aeroscalingfactor,"wind",'HorizontalAlignment','left','FontSize',6);
              axis([-axislength axislength -axislength axislength -axislength axislength]);
              if nozpanels
                %%
                pg = [(RzY*Ry*RzR*pz1')' ; (RzY*Ry*RzR*pz2')' ; (RzY*Ry*RzR*pz3')' ; (RzY*Ry*RzR*pz4')' ; (RzY*Ry*RzR*pz1')'];
                pg2 = [(RzY*Ry*RzR*pz12')' ; (RzY*Ry*RzR*pz22')' ; (RzY*Ry*RzR*pz32')' ; (RzY*Ry*RzR*pz42')' ; (RzY*Ry*RzR*pz12')'];
                pg3 = [(RzY*Ry*RzR*pz13')' ; (RzY*Ry*RzR*pz23')' ; (RzY*Ry*RzR*pz33')' ; (RzY*Ry*RzR*pz43')' ; (RzY*Ry*RzR*pz13')'];
                %vectarrow([0 0 0],totalforcevectorz(:,i,j,k));hold on;
                vectarrow([0 0 0],aeroforcevectorz);hold on;text(aeroforcevectorz(1),aeroforcevectorz(2),aeroforcevectorz(3),"aeroforce",'HorizontalAlignment','left','FontSize',6);
                vectarrow([0 0 0],Igz);hold on;text(Igz(1),Igz(2),Igz(3),"normal",'HorizontalAlignment','left','FontSize',6);
                line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));hold on;
                axis([-axislength axislength -axislength axislength -axislength axislength]);
              end
              if noxpanels
                %% panel and normal
                pgx = [(RzY*Ry*RzY*-Ry90*pz1')' ; (RzY*Ry*RzR*-Ry90*pz2')' ; (RzY*Ry*RzR*-Ry90*pz3')' ; (RzY*Ry*RzR*-Ry90*pz4')' ; (RzY*Ry*RzR*-Ry90*pz1')'];
                pgx2 = [(RzY*Ry*RzR*-Ry90*pz12')' ; (RzY*Ry*RzR*-Ry90*pz22')' ; (RzY*Ry*RzR*-Ry90*pz32')' ; (RzY*Ry*RzR*-Ry90*pz42')' ; (RzY*Ry*RzR*-Ry90*pz12')'];
                pgx3 = [(RzY*Ry*RzR*-Ry90*pz13')' ; (RzY*Ry*RzR*-Ry90*pz23')' ; (RzY*Ry*RzR*-Ry90*pz33')' ; (RzY*Ry*RzR*-Ry90*pz43')' ; (RzY*Ry*RzR*-Ry90*pz13')'];
                %vectarrow([0 0 0],totalforcevector(:,i,j,k)); hold on;
                vectarrow([0 0 0],Igx);hold on;text(Igx(1),Igx(2),Igx(3),"normalx",'HorizontalAlignment','left','FontSize',6);hold on;
                vectarrow([0 0 0],aeroforcevectorx);hold on;text(aeroforcevectorx(1),aeroforcevectorx(2),aeroforcevectorx(3),"aeroforcex",'HorizontalAlignment','left','FontSize',6);
                line(pgx(:,1), pgx(:,2), pgx(:,3));line(pgx2(:,1), pgx2(:,2), pgx2(:,3));line(pgx3(:,1), pgx3(:,2), pgx3(:,3));hold on;
                axis([-axislength axislength -axislength axislength -axislength axislength]);
              end
              if noypanels
                %% panel and normal
                pgy = [(RzY*Ry*RzR*-Rx90*pz1')' ; (RzY*Ry*RzR*-Rx90*pz2')' ; (RzY*Ry*RzR*-Rx90*pz3')' ; (RzY*Ry*RzR*-Rx90*pz4')' ; (RzY*Ry*RzY*-Rx90*pz1')'];
                pgy2 = [(RzY*Ry*RzR*-Rx90*pz12')' ; (RzY*Ry*RzR*-Rx90*pz22')' ; (RzY*Ry*RzR*-Rx90*pz32')' ; (RzY*Ry*RzR*-Rx90*pz42')' ; (RzY*Ry*RzR*-Rx90*pz12')'];
                pgy3 = [(RzY*Ry*RzR*-Rx90*pz13')' ; (RzY*Ry*RzR*-Rx90*pz23')' ; (RzY*Ry*RzR*-Rx90*pz33')' ; (RzY*Ry*RzR*-Rx90*pz43')' ; (RzY*Ry*RzR*-Rx90*pz13')'];
                %vectarrow([0 0 0],totalforcevector(:,i,j,k));hold on;                  
                vectarrow([0 0 0],Igy);hold on;text(Igy(1),Igy(2),Igy(3),"normaly",'HorizontalAlignment','left','FontSize',6);
                vectarrow([0 0 0],aeroforcevectory);hold on;text(aeroforcevectory(1),aeroforcevectory(2),aeroforcevectory(3),"aeroforcey",'HorizontalAlignment','left','FontSize',6);
                line(pgy(:,1), pgy(:,2), pgy(:,3));line(pgy2(:,1), pgy2(:,2), pgy2(:,3));line(pgy3(:,1), pgy3(:,2), pgy3(:,3));hold on;                        axis([-axislength axislength -axislength axislength -axislength axislength]);
              end
            end
            if norm(solarpressure)
              vectarrow([0 0 0],sunlight*solarpressure*panelsurface*sunscalingfactor);hold on;text(sunlight(1)*solarpressure*panelsurface*sunscalingfactor,sunlight(2)*solarpressure*panelsurface*sunscalingfactor,sunlight(3)*solarpressure*panelsurface*sunscalingfactor,"sunlight",'HorizontalAlignment','left','FontSize',6);
              if nozpanels
                vectarrow([0 0 0],sunforcevectorz);hold on;text(sunforcevectorz(1),sunforcevectorz(2),sunforcevectorz(3),"sunforce",'HorizontalAlignment','left','FontSize',6);
                vectarrow([0 0 0],Igz);hold on;text(Igz(1),Igz(2),Igz(3),"normal",'HorizontalAlignment','left','FontSize',6);
                line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));hold on;
                axis([-axislength axislength -axislength axislength -axislength axislength]);
              end
              if noxpanels
                vectarrow([0 0 0],Igx);hold on;text(Igx(1),Igx(2),Igx(3),"normalx",'HorizontalAlignment','left','FontSize',6);hold on;
                vectarrow([0 0 0],sunforcevectorx);hold on;text(sunforcevectorx(1),sunforcevectorx(2),sunforcevectorx(3),"sunforcex",'HorizontalAlignment','left','FontSize',6);
                line(pgx(:,1), pgx(:,2), pgx(:,3));line(pgx2(:,1), pgx2(:,2), pgx2(:,3));line(pgx3(:,1), pgx3(:,2), pgx3(:,3));hold on;
                axis([-axislength axislength -axislength axislength -axislength axislength]);
              end
              if noypanels
                vectarrow([0 0 0],Igy);hold on;text(Igy(1),Igy(2),Igy(3),"normaly",'HorizontalAlignment','left','FontSize',6);
                vectarrow([0 0 0],sunforcevectory);hold on;text(sunforcevectory(1),sunforcevectory(2),sunforcevectory(3),"sunforcey",'HorizontalAlignment','left','FontSize',6); 
                line(pgy(:,1), pgy(:,2), pgy(:,3));line(pgy2(:,1), pgy2(:,2), pgy2(:,3));line(pgy3(:,1), pgy3(:,2), pgy3(:,3));hold on;
                axis([-axislength axislength -axislength axislength -axislength axislength]);
              end
            end
            %text(0,0,0,strcat('scalingfactor: ',int2str(scalingfactor)),'HorizontalAlignment','left','FontSize',6);
            title(strcat('aeroscalingfactor: ',int2str(aeroscalingfactor),'  sunscalingfactor: ',int2str(sunscalingfactor)));
            xlabel('fx [N]');ylabel('fy [N]');zlabel('fz [N]');
            axis([-axislength axislength -axislength axislength -axislength axislength]);
            hold off;
            pause(1/rotspeed);
          end
          totalforcevectorz(:,i,j,k)';
          totalforcevectorx(:,i,j,k)';
          totalforcevectory(:,i,j,k)';
          totalforcevector(:,i,j,k)=totalforcevectorz(:,i,j,k)+totalforcevectorx(:,i,j,k)+totalforcevectory(:,i,j,k);
          %fprintf('\n---');
          if i*j*k==1
            %fprintf('\nbbbbbbbbbbbbbbbbbbbb');
          end
          %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
          %fprintf('\nprogress: %2.1e %%',(i+(j-1)*size(beta,2)+(k-1)*size(gamma,2)*size(beta,2))    /(size(alpha,2)*size(beta,2)*size(gamma,2))*100);
        end
      end
    end
    thetaaero=squeeze(thetaaero(i,:,k));
    phiaero=squeeze(phiaero(i,:,k));
    aerodragcoef=squeeze(aerodragcoef(:,j,k));
    aeroliftcoef=squeeze(aeroliftcoef(:,j,k));
    figure
      subplot(2,1,1)
        plot(beta,thetaaero,beta,phiaero);
        legend('thetaaero','phiaero');grid on;
        subplot(2,1,2)
        plot(alpha,aerodragcoef,alpha,aeroliftcoef);
        legend('drag','lift');grid on;
    fprintf('\ndone\n');
    save(filename,'totalforcevector')
  end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                       
function [theta,phi,Ig2]=thetaphi1(refvec, vec)
%% this function computes theta,phi, Ig2
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));                                                                                                       
  if theta>90                                                                                                                                                      
    theta=180-theta;                                                                                                                                             
    Ig2=-vec;                                                                                                                                                     
  else                                                                                                                                                             
    Ig2=vec;                                                                                                                                                     
  end                                                                                                                                                              
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)+1e-9) );  
end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [drag,lift]=aerodragliftcoef(theta,phi)
  dragcoefficientOfMasterSat=0.6;
  drag=-abs(1.2*sind(theta-90))+dragcoefficientOfMasterSat; %%simplified formula
  lift=-abs(0.12*sind(2*(theta-90)));%% simplified formula     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [drag,lift]=sundragliftcoef(theta,phi)
  drag=-abs(sind(theta-90)); %%simplified formula
  lift=-abs(sind(2*(theta-90)));%% simplified formula     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%