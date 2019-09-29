%% cosmosFS - Flight Simulator
%% emulates 
%%    -realtime code execution through timed execution and pauses
%%    -execution on different satellites through parallel programming


%% to do:
%% -use of Matlab timer function, maybe
%% -count memory and flops
%% -elliptical orbit
%% -start further workers to emulate OBC, GPS for mode switching
%% -HIL
%% -clarify what data of the header goes to each satellite

%% Revision history:
%% 26/09/2019:  cycle works, computational algorithm incomplete

clc;clear all;close all;
oldpath = path; path(oldpath,'..\matlabfunctions\')

%% symbolic names of initial conditions and desired statevector functions
%sstInitialFunction=@IRSRendezvousInitial;
sstInitialFunction=@IvanovFormationFlightInitial;
%sstInitialFunction=@cluxterInitial;

%% actual initial conditions of ODE, altitude is not used here therefore the ~
[~,ns,~,panels,rho,v,radiusOfEarth,~,mu,satelliteMass,panelSurface,...
  sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite]=sstInitialFunction(); 

DQ = parallel.pool.DataQueue;
afterEach(DQ,@disp);
parpool(ns);

startTime=posixtime(datetime('now')); %% posixtime, i.e. seconds
accelerationFactor=10000;
maxOrbits=30;

%% initial idx and altitude
idx=120;
altitude=340000;  

%% data that will later be per satellite and therefore inside SPMD loop
orbitSection    =2;         %degree
orbitSections   =[1:orbitSection:360];
orbitCounter    =0;  
etemp           =zeros(6,1);
sst             =zeros(9,1);
sstDesired      =zeros(6,1);
sstOld          =zeros(9,1);
refPosChangeTemp=zeros(3,1);
sstPP           =zeros(9,1);
%% non-gravitational perturbations
wind          =windOn*rho/2*v^2*[-1 0 0]';
sunlight      =sunOn*2*4.5e-6*[0 0 -1]' ;       %% at reference location, needs to be rotated later

%refsurf=panelSurface*panels(1);
refSurf=panelSurface*panels(3);

%% forcevector determination and its angular granulaty
alphas            =0:deltaAngle:360;     %% roll
betas             =0:deltaAngle:180;     %% pitch 
gammas            =0:deltaAngle:360;     %% yaw

aeropressureforcevector  =aeropressureforcevectorfunction(wind,panelSurface,panels(1),...
                                                          panels(2),panels(3),alphas,betas,gammas);
solarpressureforcevector =solarpressureforcevectorfunction(sunlight,panelSurface,panels(1),...
                                                          panels(2),panels(3),alphas,betas,gammas);

spmd(ns) %% create satellite instances
%%   this is code for each satellite, some of the code of above will come down here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  alive=1;
  send(DQ,strcat('No. ',num2str(labindex),' is alive.'))
  while alive
    %% satellites are alive but not doing anything
    [goFoFli,batteryOK]=getMode(maxOrbits,orbitCounter,DQ);
    if batteryOK
      %switch on GPS
    end
    sstPP   =sst;
    timePP  =0;
    
    while (goFoFli==1 || goFoFli==2)  %% orbit loop
      orbitCounter=orbitCounter+1;
      send(DQ,strcat(num2str(labindex),': ------ no of orbit: ',num2str(orbitCounter)));
      startOrbit=now; %% posixtime, i.e. seconds
      %send(DQ,strcat(num2str(labindex),': orbittimer begin: ',num2str(posixtime(datetime('now'))-startTime)));

      [meanMotion,meanAnomalyFromAN,altitude,sst]=whereInWhatOrbit(sst,altitude,(idx-1)/size(orbitSections,2));
      %% settings for control algorithm
      [P,IR,A,B]=riccatiequation(meanMotion/180*pi);  
      %% 
      %send(DQ,strcat(num2str(labindex),': MeanMotion: ',num2str(meanMotion)));
      %send(DQ,strcat(num2str(labindex),': meanAnomalyFromAN: ',num2str(meanAnomalyFromAN)));
      
      %% wait until end of orbit section
      idx = find(orbitSections >= meanAnomalyFromAN,1,'first');
      idx=idx+1;

      pause((orbitSections(idx)-meanAnomalyFromAN)/meanMotion/accelerationFactor);

      while (goFoFli==1 || goFoFli==2) && idx<=size(orbitSections,2) %% orbit sections loop
        startSection=now; %% determine cycle start time in order allow to subtract cycle duration from waiting
        %% set attitude computed in last iteration
        setAttitude();
        %% compute attitude for next section
        %% determine desired trajectory
        sstDesired=sstDesiredFunction(orbitSections(idx)/meanMotion,meanMotion/180*pi,labindex,goFoFli);
        %% determine error
        etemp(1:6)=sst(1:6)-sstDesired(1:6);
        if not(masterSatellite)
          %! ISL error
          %! compute average error
          %! ISL averageerror
          %! assign average error, i.e. compute final error
        end
        %% compute attitude
        sstOld=sst;
        [sst,~]=HCWEquation(...
          IR,P,A,B,orbitSection/meanMotion,sstOld,etemp,...
          sqrt(wind(1)^2+wind(2)^2+wind(3)^2),sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2),...
          alphas,betas,gammas,aeropressureforcevector,solarpressureforcevector,sstOld(7),...
          sstOld(8),sstOld(9),refSurf,satelliteMass,wakeAerodynamics,masterSatellite,...
          orbitSections(idx)/meanMotion,radiusOfEarth,altitude,meanMotion/180*pi);
        sst=sst';
        timePP=[timePP timePP(end)+orbitSection/meanMotion];  %% for plotting
        if not(masterSatellite)
          if labindex==1
            refPosChangeTemp(1:3)=sst(1:3)-sstOld(1:3);
          end
          %! move coordinate system
        end
        %% set counter and go condition
        idx=idx+1;
        [goFoFli,batteryOK]=getMode(maxOrbits,orbitCounter,DQ);      
        %% wait until next section starts
        %send(DQ,strcat(num2str(labindex),': time4section',num2str(orbitSection/meanMotion),' now-start ',num2str(now-start)));
        if idx>3/4*size(orbitSections,2) && batteryOK %% orbit determination will be switch on when in sun light and battery is ok
          orbitDetermination();
        end
        pause(orbitSection/meanMotion/accelerationFactor-(now-startSection));
        sstPP=[sstPP sst]; %% for plotting
      end %% orbit sections while loop

      %send(DQ,strcat(num2str(labindex),': orbittimer end: ',num2str(posixtime(datetime('now'))-startTime)));  
      send(DQ,strcat(num2str(labindex),': duration of orbit: ',num2str(now-startOrbit)));  
    end %% orbit while loop 
    send(DQ,strcat('No.',num2str(labindex),' is dead.'))
    if goFoFli~=1 %% if orbit is broken, also break alive loop, this will change later with other conditions
      alive=0;
    end
  end %% alive loop

end %% parallel computation

%% post processing for/and vizualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sstPlot=zeros(size(sstPP{1},1),ns,size(sstPP{1},2));
for iii=1:ns
  sstPlot(:,iii,:)=sstPP{iii};
end

timePlot=timePP{1};

%% results' plotting
u             =zeros(3,ns,size(sstPP{1},2));
e             =zeros(6,ns,size(sstPP{1},2));
refPosChange  =zeros(3,size(sstPP{1},2));
plotting(sstPlot(7:9,:,:),sstPlot(1:6,:,:),refPosChange,timePlot,ns,meanMotion{1}/180*pi,u,e);

delete(gcp)

function setAttitude()
  pause(0.001);
end

function [meanMotion,meanAnomalyFromAN,altitude,sst]=whereInWhatOrbit(sst,altitude,endOfSectionsCycle)
%% this function usually provides a meanAnomalyFromAN=0 and the related meanMotion for a circular orbit
  
  %% info that will come from GPS or TLE
  [GPSTLEdataAvailable,altitudeGPSTLE,meanAnomalyFromANGPSTLE,time]=getGNSSOrTLEdata();
    
  if GPSTLEdataAvailable
    altitude=altitudeGPSTLE;
    meanAnomalyFromAN=meanAnomalyFromANGPSTLE;
    %% compute SST if possible
  else
    if endOfSectionsCycle
      meanAnomalyFromAN=0.01;
    else
      meanAnomalyFromAN=120;
    end    
  end
  %% use sst, meanAnomalyFromAN and altitude either from GPS or from input parameters, %! define rule
  [rho,v,radiusOfEarth,mu,meanMotion]=orbitalproperties(altitude);
  meanMotion=meanMotion/pi*180;
end


function orbitDetermination()
%% switches orbit determination based on GPS or TLE on
%% this function needs to run in parallel to the main thread and must provide trajectory to getGNSSOrTLEdata()
%% how the heck does this work?
  GPSMethod=1;
  if GPSMethod==1
    ;
    %switchGPS(1);
    %getGPSData();
    %computeTrajectoryGPS();
  else
    ;
    %getTLE();
    %computeTrajectoryTLE();
  end
end

function [GPSTLEdataAvailable,altitudeGPSTLE,meanAnomalyFromANGPSTLE,time]=getGNSSOrTLEdata()
%% compute meanAnomalyFromANGPS,altitudeGPS at now from available past GPS or TLE/SGP4 data
  GPSTLEdataAvailable=0;
  altitudeGPSTLE=1;
  meanAnomalyFromANGPSTLE=1;
  time=1;
end

function [goFoFli,batteryOK]=getMode(maxOrbits,orbitCounter,DQ)
  goFoFli=1;
  %% readfile
  readModeFromFile=10; %% default: go to alive loop/mode
  if orbitCounter>=maxOrbits
    send(DQ,strcat(num2str(labindex),'leaving loop - maximum number of orbits reached'))
    goFoFli=10;
  elseif readModeFromFile==0
    send(DQ,strcat(num2str(labindex),'leaving loop - filestop'))
    goFoFli=20;
  else
    if orbitCounter<15
      goFoFli=1;
    else
      goFoFli=2;
    end
  end
  batteryOK=1;
end

function go=setMode(startTime,endTime,DQ)
  go=1;
end

function satellite(ns,i,DQ)
%timeline=0:1:100; 
%error=zeros(ns,1);
%averageError=0;
%  for t=1:size(timeline,2)
%  end
  OK = 9876;
  data = 1234;
  if labindex==1
    send(DQ,labindex);    
    labSend(200,2);
    labSend(300,3);
  elseif labindex==2
    send(DQ,labindex);
    data=labReceive(1);
    send(DQ,data);    
  elseif labindex==3
    send(DQ,labindex);
    data=labReceive(1);
    send(DQ,data);    
  else
    send(DQ,1234567890)
  end 
end


function somethinguseful(DQ)
    a=assigna();
    %% compute average 
    ave2=a/ns;
    for i=1:ns
      if not(i==labindex)
        ave2=ave2+labSendReceive(i,i,a)/ns;
      end
    end
    %% check whether all averages are the same
    for i=1:ns %% all workers
      if not(i==labindex) %% except self
        if ave2==labSendReceive(i,i,ave2) %% is okay?
          send(DQ,'check ok')
        end
      end
    end
end

function a=assigna()
  switch labindex
    case 1
      a=6;
    case 2
      a=7;
    case 3
      a=8;
  end
end
