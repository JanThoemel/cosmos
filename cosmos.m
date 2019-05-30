%% the Science of the Operations of Mini Satellite fOrmations and Constellations
%% COSMOS V0.1

%% based on: 
%% Deployment and Maintenance of Nanosatellite Tetrahedral Formation Flying Using Aerodynamic Forces
%% 69th International Astronautical Congress (IAC), Bremen, Germany, 1-5 October 2018.

%% Revision history:
%% 20/04/2019:  runs arbitrary number of satellites with HCW equation without control, looks much like fig 9 of Ivanov
%% 04/05/2019:  LQR control implemented but without aerodynamics 
%% 09/05/2019: aerodynamics in local-x direction,i.e. pitch implemented, correctness to be checked
%%
%% todo:
%% implement aerodynamics and visualization of angles
%% arbitrary inclinations
%% clean up of code
%% eliptical orbits
%% SGP4

clear all;close all;%format long;

altitude=340000;        %% in m
density=1e-8;
radiusOfEarth=6371000;  %% in m;
sc=2;                   %% scale for second configuration; scales also figure

%%satelliteshapeproperties, number of 10cmx10cm faces to x,y,z(body coordinates, normally aligned with):
faces=[1 3 3];
modelfilename=strcat('figures',filesep,'cocu.dae');

mu=3.986004418E14;      %% in m3?s?2
r0=radiusOfEarth+altitude; %% in m
omega=sqrt(mu/r0^3);
[P,IR,A,B]=riccatiequation(omega);

ns=2;                   %% number of satellites

%% scaling factor
CCC=1000;               %% size of formation
CCCC=1000;              %% deployment

%% timespan to simulate
totaltime=16*60*60;     %% in s
currenttime=0;          %%
ttime=[0];
tspan=0:60:1800;        %% duration and interpolation timestep for each control loop. ODE45 chooses its one optimal timestep. 150s total duration is consistent with Ivanov
%angles(ns,size(tspan,1));
%size(angles)
%% initial conditions of ODE
sst=zeros(6,ns);        %% columns: statevector per each satellite
ssttemp=zeros(6,ns,size(tspan,2));
anglestemp=zeros(3,ns,size(tspan,2));
angles=zeros(3,ns);

for i=1:ns
    %ssttemp(:,i,end)=CCCC*[(i-1)*5 0 0 0 0 0]';
    ssttemp(:,i,end)=CCCC*[(i-1)*5 0 0 0.015 0.015 0.015]';
    ssttemp(:,i,end)=CCCC*[(i-1)*5 0 0 0.00015 0.00015 0.00015]';
    %ssttemp(:,i,end)=CCCC*[(i-1)*5 0 0 rand*0.015 rand*0.015 rand*0.015]';
    %ssttemp(:,i,end)=CCCC*[(i-1)*0.5*10 0 0 0.0 0.015 0.015]';
    %sst(:,1)=[0 0 0 0.015 0.005 0.015]';
    %sst(:,2)=[15 0 0 0.005 0.015 0.015]';
end
sst(:,:)=ssttemp(:,:,end);

%% solving the ODE for each control step

%% desired statevector I
AAA=100;
DDD=115;
xd=zeros(6,ns,size(tspan,2)); %% columns: desired statevector per each satellite

%% parameters for ODE45
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

tic
while currenttime<totaltime
  if currenttime<8*60*60    %% initial desired vector
    %% desired statevector II
    xd(1,1,:)=CCC*DDD;
    xd(1,2,:)=-CCC*DDD;
    %xd(1,3,:)=CCC*2*AAA*cos(omega*(currenttime+tspan)-acos(1/3));
    %xd(2,3,:)=CCC*AAA*sqrt(3)*sin(omega*(currenttime+tspan));
    %xd(3,3,:)=CCC*AAA*sin(omega*(currenttime+tspan)-acos(1/3));
    
    %xd(1,4,:)=CCC*2*AAA*cos(omega*(currenttime+tspan));
    %xd(2,4,:)=CCC*AAA*sqrt(3)*sin(omega*(currenttime+tspan)+acos(1/3));
    %xd(3,4,:)=CCC*AAA*sin(omega*(currenttime+tspan));
  else                      %% switch to new desired statevector
    xd(1,1,:)=sc*CCC*DDD;
    xd(1,2,:)=sc*-CCC*DDD;
    
    xd(1,3,:)=sc*CCC*2*AAA*cos(omega*(currenttime+tspan)-acos(-3/3));
    xd(2,3,:)=sc*CCC*AAA*sqrt(3)*sin(omega*(currenttime+tspan));
    xd(3,3,:)=sc*CCC*AAA*sin(omega*(currenttime+tspan)-acos(1/3));
    
    xd(1,4,:)=sc*CCC*2*AAA*cos(omega*(currenttime+tspan)-acos(5/3));
    xd(2,4,:)=sc*CCC*AAA*sqrt(3)*sin(omega*(currenttime+tspan)+acos(1/3));
    xd(3,4,:)=sc*CCC*AAA*sin(omega*(currenttime+tspan));   
  end
  for i=1:ns            %% solve ODE
    [t,x]=ode23(@(t,x) hcwequation(t,x,IR,P,A,B,xd(:,i,:),currenttime,tspan),tspan,ssttemp(:,i,size(ssttemp,3)),opts);
    ssttemp(:,i,:)=x'; %% columns: statevector, rows: (x,y,z,u,v,w), depth: timeevolution 
    [a1 anglestemp(:,i,:)]=hcwequation(t,x',IR,P,A,B,xd(:,i,:),currenttime,tspan);
    %angles2
    clear hcwequation
    %input('--------------------new satellite------------')
  end
  sst=cat(3,sst,ssttemp(:,:,2:end));  
  angles=cat(3,angles,anglestemp(:,:,2:end));  
  ttime=[ttime ; currenttime+t(2:end)];  
  if currenttime>0
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
  end
  currenttime=ttime(end);
  fprintf('time %2.1e / totaltime %2.1e, progress %2.1e %%', currenttime, totaltime,currenttime/totaltime*100);
  %input('--------------------new iteration------------')
end
figure
plot(squeeze(ttime),squeeze(angles(1,1,:)),squeeze(ttime),squeeze(angles(1,2,:)))
anglesX=zeros(3,ns,size(ttime,1));
for i=1:ns
  for j=1:2 %% two rotations
    for k=1:size(ttime,1)
      anglesX(j,i,k)=89*cos(2*pi/60/30*ttime(k)+i*pi/3+j*3/5*pi); %heading
      anglesX(j,i,k)=89*sin(2*pi/60/30*ttime(k)+i*pi/3+j*3/5*pi); %tilt %roll would be the third component, not specified here
    end
  end
end

%printf('\n');
toc
%%  plotting the results
figure
for i=1:ns
    plot3(squeeze(sst(1,i,:)),squeeze(sst(2,i,:)),squeeze(sst(3,i,:)),'-')
    hold on
    names(i)=[{strcat('sat',int2str(i))}];
end
xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]');legend(names);grid on;hold off
if 0
    axis(sc*[-35e3 35e3 -35e3 35e3 -35e3 35e3])
    xticks(sc*[-20e3 0 20e3])
    yticks(sc*[-20e3 0 20e3])
    zticks(sc*[-20e3 0 20e3])
    set(gcf,'PaperUnits', 'inches','PaperPosition', [0 0 2 2]) ;
    fprintf('2by2','-dpng','-r0')
end
    %% map on globe and create kml file
%visualization(ns,ttime,squeeze(sst(1,:,:)),squeeze(sst(2,:,:)),squeeze(sst(3,:,:)),altitude,anglesX, modelfilename,radiusOfEarth,mu)
visualization(ns,ttime,squeeze(sst(1,:,:)),squeeze(sst(2,:,:)),squeeze(sst(3,:,:)),altitude,angles, modelfilename,radiusOfEarth,mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,IR,A,B]=riccatiequation(omega)

    %R=diag([1e-13 1e-14 1e-14]);
    R=diag([1e13 1e14 1e14]);
    Q=eye(6);

    E=eye(3);
    Z=zeros(3,3);
    C=[0 0 0; 0 -omega^2 0;0 0 3*omega^2];
    D=[0 0 -2*omega;0 0 0;2*omega 0 0];

    A=[Z E; C D];
    B=[Z;E];
    IR=inv(R);

    %%https://nl.mathworks.com/help/control/ref/care.html
    %% the function care is replaced by icare in later matlab versions
    S=zeros(6,3);
    E2=eye(6);
    [P,L,G] = care(A,B,Q,R,S,E2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt angles] =hcwequation(t,x,IR,P,A,B,xdi,currenttime,timespan)
%% HCW equation

    persistent ang;
    persistent timearray;
    %timearray=0;
    %% convergence relaxation factors
    satmass=2; %% kg
    density=1e-8; %%kg/m3
    v=7500; %% m/s
    S=0.1^2; %% m2
    k=1/satmass*density*v^2*S;
    xdi2x = interp1(currenttime+timespan,squeeze(xdi(1,1,:)),currenttime+t);% Interpolate the data set (ft,f) at time t
    xdi2y = interp1(currenttime+timespan,squeeze(xdi(2,1,:)),currenttime+t); % Interpolate the data set (ft,f) at time t
    xdi2z = interp1(currenttime+timespan,squeeze(xdi(3,1,:)),currenttime+t); % Interpolate the data set (ft,f) at time t
    xdi2=[xdi2x xdi2y xdi2z zeros(size(xdi2x,1),1) zeros(size(xdi2x,1),1) zeros(size(xdi2x,1),1)]';
    %% error vector
    e=x-xdi2;
    %% control vector
    u1=-IR*B'*P*e;
    %e(1,:)
    %u1(1,:)
    
    umax=k*1.2/satmass;
    uyzmax=k*0.24/satmass;
    %input('a ')
    u=u1;;%-[umax/2 0 0]'; %% to set u to a default value and size u correctly
    %umax
    
    for i=1:size(u1,2)      
      if u1(1,i)>=1/2*umax
        u(1,i)=1/2*umax;
        u(2,i)=0;
        u(3,i)=0;
        %printf('1\n')
       elseif u1(1,i)<=-1/2*umax
        u(1,i)=-1/2*umax;
        u(2,i)=0;
        u(3,i)=0;
        %printf('2\n')
      %elseif sqrt(u(2,i)^2+u(3,i)^2)>uyzmax
      %  u(1,i)=-0.8;
     %   u(2,i)=u1(2,i)/uyzmax;
     %   u(3,i)=u1(3,i)/uyzmax;
        %printf('3\n')
      end
    end
    anglestemp(1,:)=asind(1/1.2*u(1,:)/k)+asind(1/1.2*umax/k);
    anglestemp(2,:)=0;
    anglestemp(3,:)=0;
    if size(ang,1)>2 && size(x,2)==1
      ang=[ang  anglestemp];
    %  timearray=[timearray t];
    else
     ang=anglestemp;
    end
    angles=ang;
    dxdt=A*x+B*u;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualization(ns,ttime,sstx,ssty,sstz,altitude,angles,modelfilename,radiusOfEarth,mu)
    %this function works with km instead of m %!harmonize
    RE = radiusOfEarth/1000;          % Earth's radius                            [km]
    muE = mu/1e9;    % Earth gravitational parameter             [km^3/sec^2]
    wE = (2*pi/86164);  % Earth rotation velocity aorund z-axis     [rad/sec]
    %% ORBIT INPUT
    RAAN=0; %%RAAN    = input(' Right Ascension of Ascendent Node    [  0,360[    RAAN   [deg] = ');    
    w=0;    %%w       = input(' Argument of perigee                  [  0,360[    w      [deg] = ');
    v0=0;   %%v0      = input(' True anomaly at the departure        [  0,360[    v0     [deg] = ');
    inclination=89.999999;%i       = input(' Inclination                          [-90, 90]    i      [deg] = ');
    a=RE+altitude/1000;%a       = input(' Major semi-axis                       (>6378)     a      [km]  = ');
    ecc_max = sprintf('%6.4f',1-RE/a);     % maximum value of eccentricity allowed
    e=0;%e       = input([' Eccentricity                         (<',ecc_max,')    e            = ']);
    RAAN  = RAAN*pi/180;        % RAAN                          [rad]
    w     = w*pi/180;           % Argument of perigee           [rad]
    v0    = v0*pi/180;          % True anomaly at the departure [rad]
    inclination     = inclination*pi/180;           % inclination                   [rad]
    %% ORBIT COMPUTATION
    rp = a*(1-e);               % radius of perigee             [km]
    ra = a*(1+e);               % radius of apogee              [km]
    Vp = sqrt(muE*(2/rp-1/a));  % velocity at the perigee       [km/s]
    Va = sqrt(muE*(2/ra-1/a));  % velocity at the  apogee       [km/s]
    n  = sqrt(muE./a^3);        % mean motion                   [rad/s]
    p  = a*(1-e^2);             % semi-latus rectus             [km]
    T  = 2*pi/n;                % period                        [s]
    h  = sqrt(p*muE);           % moment of the momentum        [km^2/s]
    h1 = sin(inclination)*sin(RAAN);      % x-component of unit vector h
    h2 = -sin(inclination)*cos(RAAN);     % y-component of unit vector h
    h3 = cos(inclination);                % z-component of unit vector h
    n1 = -h2/(sqrt(h1^2+h2^2)); % x-component of nodes' line
    n2 =  h1/(sqrt(h1^2+h2^2)); % y-component of nodes' line
    n3 = 0;                     % z-component of nodes' line
    N  = [n1,n2,n3];            % nodes' line (unit vector)
    %% PRINT SOME DATAS
    hours   = floor(T/3600);                   % hours   of the orbital period
    minutes = floor((T-hours*3600)/60);        % minutes of the orbital period
    seconds = floor(T-hours*3600-minutes*60);  % seconds of the orbital period
    t0   = 0;                                        % initial time          [s]
    tf=ttime(end);
    step=60;%step = input(' Time step.        [s] step = ');  % time step             [s]    
    t    = t0:step:tf;                          % vector of time        [s]
    %% DETERMINATION OF THE DYNAMICS
    cosE0 = (e+cos(v0))./(1+e.*cos(v0));               % cosine of initial eccentric anomaly
    sinE0 = (sqrt(1-e^2).*sin(v0))./(1+e.*cos(v0));    %   sine of initial eccentric anomaly
    E0 = atan2(sinE0,cosE0);                           % initial eccentric anomaly              [rad]
    if (E0<0)                                          % E0€[0,2pi]
        E0=E0+2*pi;
    end
    tp = (-E0+e.*sin(E0))./n+t0;                       % pass time at the perigee               [s]
    M  = n.*(t-tp);                                    % mean anomaly                           [rad]
    %% Mk = Ek - e*sin(Ek);
    % Eccentric anomaly (must be solved iteratively for Ek)
    E = zeros(size(t,2),1);
    for j=1:size(t,2)
        E(j) = anom_ecc(M(j),e);                     % eccentric anomaly         [rad]
    end
    %% True anomaly, Argument of latitude, Radius
    sin_v = (sqrt(1-e.^2).*sin(E))./(1-e.*cos(E));   % sine of true anomaly
    cos_v = (cos(E)-e)./(1-e.*cos(E));               % cosine of true anomaly
    v     = atan2(sin_v,cos_v);                      % true anomaly              [rad]
    theta = v + w;                                   % argument of latitude      [rad]
    r     = (a.*(1-e.^2))./(1+e.*cos(v));            % radius                    [km]
    %% Satellite coordinates
    % "Inertial" reference system ECI (Earth Centered Inertial)
    xp = r.*cos(theta);                          % In-plane x position (node direction)             [km]
    yp = r.*sin(theta);                          % In-plane y position (perpendicular node direct.) [km]
    xs = xp.*cos(RAAN)-yp.*cos(inclination).*sin(RAAN);    % ECI x-coordinate SAT                             [km]
    ys = xp.*sin(RAAN)+yp.*cos(inclination).*cos(RAAN);    % ECI y-coordinate SAT                             [km]
    zs = yp.*sin(inclination);                             % ECI z-coordinate SAT                             [km]
    rs = p./(1+e.*cos(theta-w));                 % norm of radius SAT                               [km]
    %% GREENWICH HOUR ANGLE
    disp(' From ephemeridis you can reach the greenwich hour angle at the epoch and reset it from Aries'' point');
    disp(' Ephemeridis Almanac: http://www.nauticalalmanac.it/it/astronomia-nautica/effemeridi-nautiche.html ');
    % Greenwich hour angle è l'orientamento della Terra relativo al punto di
    % Ariete in un'epoca precisa. Tale valore è tabulato per ogni giorno
    % dell'anno, ogni ora, minuto ecc. nelle effemeridi che si possono trovare
    % negli almanacchi nautici in generale. Se non lo si conosce, o si intende
    % fare un'analisi approssimativa  è consigliabile immettere 0 in modo da
    % far coincidere il punto di ariete con Greenwich all'epoca.
    %greenwich0 = input(' Insert Greenwich longitude respect to the vernal axis. GHA [deg] = ');
    greenwich0=0;
    %% SUB-SATELLITE-POINT
    greenwich0 = greenwich0*pi/180;                 % Greenwich hour angle at the initial time    [rad]
    rot_earth  = wE.*(t-t0)+greenwich0;             % Greenwich hour angle at the time t          [rad]
    for j=1:size(t,2)
        if rot_earth(j) < (-pi)
            nn = ceil(rot_earth(j)/(-2*pi));
            rot_earth(j) = rot_earth(j) + nn*2*pi;
        elseif rot_earth(j) > (pi)
            nn = fix(rot_earth(j)/(2*pi));
            rot_earth(j) = rot_earth(j) - nn*2*pi;
        end
    end
    %% interpolate relative position on visualization grid
    for j=1:ns
      sstxvgrid(j,:)=interp1(ttime,sstx(j,:),t);
      sstyvgrid(j,:)=interp1(ttime,ssty(j,:),t);
      sstzvgrid(j,:)=interp1(ttime,sstz(j,:),t);
    end
    for j=1:size(ttime,1)
      for i=1:ns
        MAXY=sqrt(max(sstxvgrid(i,j)^2+sstyvgrid(i,j)^2));
      end
      footprintsize(j)=round(  MAXY/10^floor(log10(MAXY)) )  *  1e5;
    end
    %% centerpoint
    Lat(1,:)     =  asin(sin(inclination).*sin(theta))/pi*180;              % Latitude             [deg]
    Long(1,:)    = (atan2(ys./rs,xs./rs)-rot_earth')/pi*180;    % Longitude            [deg]
    rs2(1,:)         = rs;                                        % radius                [km]
    %% satellites
    for j=1:ns
      ii(j,:)      = asin(  sstxvgrid(j,:)/1000  ./  (rs(:,1)'  + sstzvgrid(j,:)/1000) ); %%latitude offset
      ll(j,:)      = asin( sstyvgrid(j,:)/1000   ./   (rs(:,1)' + sstzvgrid(j,:)/1000) ); %%longitude offset
    end
    for j=2:ns+1 
      Lat(j,:)     = Lat(1,:)+ii(j-1,:)/pi*180;               % Latitude             [deg]
      Long(j,:)    = Long(1,:)+ll(j-1,:)/pi*180;    % Longitude            [deg]
      rs2(j,:)     = rs(:)'+sstzvgrid(j-1,:)/1000;                                    % radius                [km]
    end
    %figure
    %subplot(3,1,1)
    %for j=1:ns
    %  plot(t,sstxvgrid(j,:))
    %  hold on
    %end  
    %subplot(3,1,2)
    %for j=1:ns
    %  plot(t,sstyvgrid(j,:))
    %  hold on
    %end
    %subplot(3,1,3)
    %for j=1:ns
    %  plot(t,sstzvgrid(j,:))
    %  hold on
    %end
    %figure
    %subplot(2,1,1)
    %for j=1:ns
    %  plot(t,ii(j,:))
    %  hold on
    %end
    %subplot(2,1,2)
    %for j=1:ns
    %  plot(t,ll(j,:))
    %  hold on
    %end
    %figure
    %subplot(3,1,1)
    %for j=1:ns+1
    %  plot(t,Lat(j,:))
    %  hold on
    %end
    %subplot(3,1,2)
    %for j=1:ns+1
    %  plot(t,Long(j,:))
    %  hold on
    %end
    %subplot(3,1,3)
    %for j=1:ns+1
    %  plot(t,rs2(j,:))
    %  hold on
    %end

    %% write KML file
    writeKML(Lat,Long,rs2,t,ns,angles,modelfilename,footprintsize)
    %writeKML(Lat,Long,rs2,t,ns,angles,modelfilename)
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeKML(Lat,Long,rs,t,ns,angles,modelfilename,footprintsize)
%function writeKML(Lat,Long,rs,t,ns,angles,modelfilename)
  
  %! simulation cannot be longer than a day
  %! why division by 60 twice?
  %! can do only inclination=90 deg
  %! the kml tour follows first satellite
  initialtime=datenum('2021-06-13T00:00:00Z','yyyy-mm-ddTHH:MM:SSZ');
  poletransition=0;
  latviewoffset=1;
  longviewoffset=-1;
  tiltview=25;
  altitudeview=600000;
  headingviewoffset=30;
  modelscalex=1500;
  modelscaley=1500;
  modelscalez=1500;
  accelerationfactor=60;
  tourjumpduration=zeros(1,size(t,2));
  tourjumpduration(2:end)=(t(2:end)-t(1:end-1))/accelerationfactor;
  filename='cosmos.kml';
  USpoints=5;
  for ii=1:3
    counter=0;
    counter2=0;
    %% add pole points
      for i=2:size(t,2)-1
          counter2=counter2+1;
          if counter~=0&&counter2==3
            counter=0;
            counter2=0;
          end
          if Lat(1,i)==max([Lat(1,i-1) Lat(1,i) Lat(1,i+1)]) %% around northpole
              for j=1:ns
                  if abs(Long(j,i))-abs(Long(j,i-1))<abs(Long(j,i))-abs(Long(j,i+1)) %% before northpole
                      fprintf('beforenp')
 
                  elseif abs(Long(j,i))-abs(Long(j,i-1))>abs(Long(j,i))-abs(Long(j,i+1)) && counter==0%% after northpole
                      %printf('\nafternp\n')
                      counter=counter+1;
                      %% add line to Lat,Long,t, tourjumpduration
                      Lat(:,i+1:end+1)=Lat(:,i:end);
                      Long(:,i+1:end+1)=Long(:,i:end);
                      rs(:,i+1:end+1)=rs(:,i:end);
                      t(i+1:end+1)=t(i:end);
                      tourjumpduration(i+1:end+1)=tourjumpduration(i:end);
                      angles(:,:,i+1:end+1)=angles(:,:,i:end);
                      footprintsize(i+1:end+1)=footprintsize(i:end);
                      Lat(:,i)=(Lat(:,i-1)+Lat(:,i+1))/2; %!interpolate
                      for k=1:ns                      
                          idx=find(min(   [abs(Lat(k,i-1)) abs(Lat(k,i+1))]   )==Lat(k,:));
                          Long(k,i)=(  Long(k,i-1)+Long(k,i+1)+sign(Long(k,idx))*180  )/2;
                      end
                        t(i)=(t(i)+t(i-1))/2;
                        tourjumpduration(i)=tourjumpduration(i)/2;
                        tourjumpduration(i+1)=tourjumpduration(i);
                        angles(:,:,i)=(angles(:,:,i-1)+angles(:,:,i+1))/2;
                        footprintsize(i)=(footprintsize(i-1)+footprintsize(i+1))/2;
                  %else
                        %printf('error')
                      end
               end
          elseif Lat(1,i)==min([Lat(1,i-1) Lat(1,i) Lat(1,i+1)]) %% around southpole
              fprintf('do something for the southpole point adding')
        end
       end
  end
  for i=1:size(t,2) %% shift longitude from 0-360 deg to +/- 180
    for j=1:ns
      if Long(j,i)<-180
        Long(j,i)=Long(j,i)+360;
      end
    end
  end 
  %%upper stage data
  USLong=[-flip(Long(1,2:USpoints))'; Long(1,1:USpoints)'];
  USLat=[-flip(Lat(1,2:USpoints))'; Lat(1,1:USpoints)'];
  USrs=[flip(rs(1,2:USpoints))'; rs(1,1:USpoints)'];
  UStime=[-flip(t(2:USpoints)) t(1:USpoints)];
  USjumpduration=1;
  %a=input('a')
  
  id=fopen(filename,'w'); %% write kml file
  fprintf(id,'<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">');
  fprintf(id,'\n  <Folder> ');
  fprintf(id,'\n  <Style id="stylesel_0"><IconStyle><Icon><href/></Icon></IconStyle><LabelStyle><scale>0.3</scale></LabelStyle></Style>');
  fprintf(id,'\n  <Style id="circle"><IconStyle><Icon><href/>figures%scircle.png</Icon></IconStyle><LabelStyle><scale>0.3</scale></LabelStyle></Style>',filesep);

  %% initial camera viewpoint
%  fprintf(id,'\n<Camera>');
%  fprintf(id,'\n  <gx:TimeSpan>');
%  fprintf(id,'\n		<begin>%s</begin><end>%s</end>',strcat(datestr(initialtime-1/24/60/60,29),'T',datestr(initialtime-1/24/60/60,13),'Z'),strcat(datestr(initialtime+1/24/60/60,29),'T',datestr(initialtime+1/24/60/60,13),'Z'));
%  fprintf(id,'\n  </gx:TimeSpan>');
%  fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,1)+longviewoffset*(90-abs(Lat(1,1) ))/90,Lat(1,1)-latviewoffset,altitudeview);
%  fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,1) ))/90,tiltview);
%  fprintf(id,'\n</Camera>');

  fprintf(id,'\n<Camera>');
  fprintf(id,'\n  <gx:TimeSpan>');
  fprintf(id,'\n		<begin>%s</begin><end>%s</end>',strcat(datestr(initialtime+UStime(1)/24/60/60-1/24/60/60,29),'T',datestr(initialtime+UStime(1)/24/60/60-1/24/60/60,13),'Z'),strcat(datestr(initialtime+UStime(1)/24/60/60+1/24/60/60,29),'T',datestr(initialtime+UStime(1)/24/60/60+1/24/60/60,13),'Z'));
  fprintf(id,'\n  </gx:TimeSpan>');
  fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',USLong(1)+longviewoffset*(90-abs(USLat(1) ))/90,USLat(1)-latviewoffset,altitudeview);
  fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(USLat(1) ))/90,tiltview);
  fprintf(id,'\n</Camera>');

  %% tour definition

  fprintf(id,'\n<gx:Tour><name>Play me!</name><gx:Playlist>');
  for i=(-USpoints+1):1:size(t,2)-1
    if i<USpoints
          fprintf(id,'\n<gx:FlyTo><gx:flyToMode>smooth</gx:flyToMode><gx:duration>%f</gx:duration><Camera>',USjumpduration);
          pointintime=strcat(datestr(initialtime+UStime(i+USpoints)/24/60/60,29),'T',datestr(initialtime+UStime(i+USpoints)/24/60/60,13),'Z');
          pointintimeb=strcat(datestr(initialtime+UStime(i+USpoints)/24/60/60-10/24/60/60,29),'T',datestr(initialtime+UStime(i+USpoints)/24/60/60-10/24/60/60,13),'Z');
          pointintimee=strcat(datestr(initialtime+UStime(i+USpoints)/24/60/60+10/24/60/60,29),'T',datestr(initialtime+UStime(i+USpoints)/24/60/60+10/24/60/60,13),'Z');
          %fprintf(id,'\n		    <gx:TimeStamp><when> %s </when></gx:TimeStamp>',pointintime');
          fprintf(id,'\n  <gx:TimeSpan>');
          fprintf(id,'\n		<begin>%s</begin><end>%s</end>',pointintimeb,pointintimee);
          fprintf(id,'\n  </gx:TimeSpan>');

          %% set viewpoint and heading to according to flight direction
          if i==-USpoints+1  %% first point
            fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',USLong(i+USpoints)+longviewoffset*(90-abs(Lat(1,1) ))/90,USLat(i+USpoints)-latviewoffset,altitudeview);
            fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,1) ))/90,tiltview);
          elseif USLat(i+USpoints)-USLat(i-1+USpoints)>=0 %% northern direction
            %printf('nd%d',i)      
            fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',USLong(i+USpoints)+longviewoffset*(90-abs(Lat(1,1) ))/90,USLat(i+USpoints)-latviewoffset,altitudeview);
            fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,1) ))/90,tiltview);
          elseif USLat(i+USpoints)-USLat(i-1+USpoints)<0  %% southern direction
            %printf('sd%d',i)
            fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',USLong(i+USpoints+1)-longviewoffset*(90-abs(Lat(1,1) ))/90,USLat(i+USpoints+1)+latviewoffset,altitudeview);
            fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',180+headingviewoffset*(90-abs(Lat(1,1) ))/90,tiltview);
          else
            fprintf('error')
          end
          fprintf(id,'\n</Camera></gx:FlyTo>');   
  elseif i==USpoints
          ;
  else
          fprintf(id,'\n<gx:FlyTo><gx:flyToMode>smooth</gx:flyToMode><gx:duration>%f</gx:duration><Camera>',tourjumpduration(i));
          pointintime=strcat(datestr(initialtime+t(i)/24/60/60,29),'T',datestr(initialtime+t(i)/24/60/60,13),'Z');
          pointintimeb=strcat(datestr(initialtime+t(i)/24/60/60-10/24/60/60,29),'T',datestr(initialtime+t(i)/24/60/60-10/24/60/60,13),'Z');
          pointintimee=strcat(datestr(initialtime+t(i)/24/60/60+10/24/60/60,29),'T',datestr(initialtime+t(i)/24/60/60+10/24/60/60,13),'Z');
          %fprintf(id,'\n		    <gx:TimeStamp><when> %s </when></gx:TimeStamp>',pointintime');
          fprintf(id,'\n  <gx:TimeSpan>');
          fprintf(id,'\n		<begin>%s</begin><end>%s</end>',pointintimeb,pointintimee);
          fprintf(id,'\n  </gx:TimeSpan>');

          %% set viewpoint and heading to according to flight direction
          if i==1  %% first point
            fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)+longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)-latviewoffset,altitudeview);
            fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);
          elseif Lat(1,i)==max([Lat(1,i-1) Lat(1,i) Lat(1,i+1)]) %% on northpole
            if abs( abs(Long(j,i)) - abs(Long(j,i-1)) )<abs( abs(Long(j,i)) - abs(Long(j,i+1)) ) %% before northpole, i.e. norther direction
                fprintf('bnp%d',i)
                fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)+longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)-latviewoffset,altitudeview);
                fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);            
            elseif abs( abs(Long(j,i)) - abs(Long(j,i-1)) )>abs( abs(Long(j,i)) - abs(Long(j,i+1)) ) %% after northpole, i.e. souther direction
                fprintf('anp%d',i)
                fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)-longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)+latviewoffset,altitudeview);
                fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',180+headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);        
            end    
          elseif Lat(1,i)==min([Lat(1,i-1) Lat(1,i) Lat(1,i+1)]) %% on southpole
            if abs(abs(Long(j,i))-abs(Long(j,i-1)))<abs(abs(Long(j,i))-abs(Long(j,i+1))) %% before southpole, i.e. southern direction
              fprintf('bsp%d',i)
              fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)-longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)+latviewoffset,altitudeview);
              fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',180+headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);    
            elseif abs(abs(Long(j,i))-abs(Long(j,i-1)))>abs(abs(Long(j,i))-abs(Long(j,i+1))) %% after southpole, i.e. northern direction
              fprintf('asp%d',i)
              fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)+longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)-latviewoffset,altitudeview);
              fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);                     
            end         
          elseif Lat(1,i+1)-Lat(1,i)>=0 %% northern direction
            %printf('nd%d',i)      
            fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)+longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)-latviewoffset,altitudeview);
            fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);
          elseif Lat(1,i+1)-Lat(1,i)<0  %% southern direction
            %printf('sd%d',i)
            fprintf(id,'\n  <longitude>%f</longitude><latitude>%f</latitude><altitude>%f</altitude>',Long(1,i)-longviewoffset*(90-abs(Lat(1,i) ))/90,Lat(1,i)+latviewoffset,altitudeview);
            fprintf(id,'\n  <heading>%f</heading><tilt>%f</tilt><roll>0</roll>',180+headingviewoffset*(90-abs(Lat(1,i) ))/90,tiltview);
          else
            fprintf('error')
          end
          fprintf(id,'\n</Camera></gx:FlyTo>');
    end
  end
  fprintf(id,'\n</gx:Playlist></gx:Tour>');
    
  %% US path
    fprintf(id,'\n<Placemark><name>US</name><styleUrl>#stylesel_0</styleUrl><gx:Track><altitudeMode>absolute</altitudeMode>');     
    fprintf(id,'\n   <Model id="geom_0"><altitudeMode>absolute</altitudeMode>');
    fprintf(id,'\n      <Orientation><heading>0.0</heading><tilt>0</tilt><roll>0</roll></Orientation>');
    fprintf(id,'\n      <Scale><x>10000</x><y>10000</y><z>10000</z></Scale>');
    fprintf(id,'\n      <Link id="link_0"><href>figures/simpleUS.dae</href></Link>');
    fprintf(id,'\n    </Model>');
    for i=1:size(UStime,2)
          fprintf(id,'\n<gx:coord> %f %f %f</gx:coord>',USLong(i),USLat(i),(USrs(i)-6371)*1000); 
    end 
    for i=1:size(UStime,2)
          pointintime=strcat(datestr(initialtime+UStime(i)/24/60/60,29),'T',datestr(initialtime+UStime(i)/24/60/60,13),'Z');
          fprintf(id,'\n<when> %s </when>',pointintime);
    end
    %for i=1:size(UStime,2)
     %     fprintf(id,'\n<gx:angles>%f %f %f</gx:angles>',angles(j-1,1,i),angles(j-1,2,i),angles(j-1,3,i));
     %end 
    fprintf(id,'\n  </gx:Track></Placemark>');

  
   %% satellites' time, angles and coordinates
   for j=2:ns+1
        fprintf(id,'\n<Placemark><name>sat%d</name><styleUrl>#stylesel_0</styleUrl><gx:Track><altitudeMode>absolute</altitudeMode>',j-1);     
        fprintf(id,'\n   <Model id="geom_0"><altitudeMode>absolute</altitudeMode>');
        fprintf(id,'\n      <Orientation><heading>90.0</heading><tilt>0</tilt><roll>0</roll></Orientation>');
        fprintf(id,'\n      <Scale><x>%f</x><y>%f</y><z>%f</z></Scale>',modelscalex,modelscaley,modelscalez);
        fprintf(id,'\n      <Link id="link_0"><href>%s</href></Link>',modelfilename);
        fprintf(id,'\n    </Model>');

        for i=1:size(t,2)
              fprintf(id,'\n<gx:coord> %f %f %f</gx:coord>',Long(j,i),Lat(j,i),(rs(j,i)-6371)*1000); 
        end 
        for i=1:size(t,2)
              pointintime=strcat(datestr(initialtime,29),'T',datestr(initialtime+t(i)/24/60/60,13),'Z');
              fprintf(id,'\n<when> %s </when>',pointintime);
        end
        for i=1:size(t,2)
              fprintf(id,'\n<gx:angles>%f %f %f</gx:angles>',angles(2,j-1,i),angles(1,j-1,i),angles(3,j-1,i));%% yaw,pitch, roll %% yaw is always 0
        end 
        fprintf(id,'\n  </gx:Track></Placemark>');
   end
    k=1;
    while k<size(footprintsize,2)
         fprintf(id,'\n<Placemark><name>centerpoint</name><styleUrl>#stylesel_0</styleUrl>>');     
         fprintf(id,'\n   <Model id="geom_1"><altitudeMode>clampToGround</altitudeMode>');
         fprintf(id,'\n       <Orientation><heading>90.0</heading><tilt>0</tilt><roll>0</roll></Orientation>');        
         fprintf(id,'\n       <Scale><x>%f</x><y>%f</y><z>%f</z></Scale>',footprintsize(k),footprintsize(k),footprintsize(k)); 
          fprintf(id,'\n     <Location>');
          fprintf(id,'\n       <longitude>%f</longitude> <!-- kml:angle180 -->',Long(1,k));
          fprintf(id,'\n       <latitude>%f</latitude>   <!-- kml:angle90 -->',Lat(1,k));
          %fprintf(id,'\n       <altitude>%f</altitude>  <!-- double -->',100);
          fprintf(id,'\n     </Location>');
         fprintf(id,'\n       <Link id="link_1"><href>figures%storus.dae</href></Link>',filesep);
         fprintf(id,'\n    </Model>');
         pointintime=strcat(datestr(initialtime,29),'T',datestr(initialtime+t(k)/24/60/60,13),'Z');
         fprintf(id,'\n     <TimeStamp><when> %s </when></TimeStamp>',pointintime);
         fprintf(id,'\n</Placemark>');  
         k=k+1;
    end
    fprintf(id,'\n  </Folder>');
    fprintf(id,'\n </kml>');
    fclose(id);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E] = anom_ecc(M,e)
    % function [E] = anom_ecc(M,e) 
    % Risoluzione numerica dell'equazione: E-e*sin(E)=M
    % E = anomalia eccentrica [rad]
    % e = eccentricità
    % M = anomalia media [rad]
    % Si devono dare in input alla funzione due valori scalari,
    % rispettivamente M in rad ed e.
    % Il programma restituisce E [rad] con un errore di 1e-10
    % N.B. Per migliorare l'accuratezza accedere all'editor e modificare il
    % valore della variabile err, che rappresenta l'errore commesso.
    
    format long
%     x = 0;
%     sx = 1;
%     dymax = -1;
%     trovato = true;
%     while (trovato)
%         if (sx<0.2)
%             sx = 0.1 - (x/1000);
%         else
%             sx  = M-x;
%             dx  = M+x;
%             dy  = abs(cos(dx));
%             dy2 = abs(cos(sx));
%         end
%         if (dymax<dy || dymax<dy2)
%             if (dy<dy2)
%                 dymax = dy2;
%             else
%                 dymax = dy;
%                 dy = dymax;
%             end
%         end 
%         f0 = sx-e.*sin(sx)-M;
%         f1 = dx-e.*sin(dx)-M;
%         trovato = ((f0.*f1)>0);
%         x = x + 0.1;
%     end
    E = M;
    k = 1;
    err = 1e-10;
    % stabilito intervallo di ricerca (sx,dx) e max valore derivata dymax;
    while (k>err)
        y = e*sin(E)+M;
        k = abs(abs(E)-abs(y));
        E = y;
    end
    % trovato E con un errore inferiore a 0.1*10^(-4);
    %fprintf(' La soluzione E è stata trovata nell''intervallo [%2.5f,%2.5f]',sx,dx);
    %fprintf(' \n errore inferiore a %1.2e: [rad] E = %3.5f',err,E);
end
