%% cosmosFS - Flight Simulator
%% emulates 
%%    -realtime code execution through timed execution and pauses
%%    -execution on different satellites through parallel programming


%% to do
%% -use of Matlab timer function, maybe

%% Revision history:
%% 26/09/2019:  cycle works, computational algorithm incomplete

clc;clear all;
%ns=3;       %% number of satellite

DQ = parallel.pool.DataQueue;
afterEach(DQ,@disp);
parpool(3);

%% symbolic names of initial conditions and desired statevector functions
%sstInitialFunction=@IRSRendezvousInitial;
sstInitialFunction=@IvanovFormationFlightInitial;
%sstInitialFunction=@cluxterInitial;

%% actual initial conditions of ODE, altitude is not used here therefore the ~
[sstTemp,ns,~,panels,rho,v,radiusOfEarth,MeanMotion,mu,satelliteMass,panelSurface,...
  sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite]=sstInitialFunction(); 


startTime=posixtime(datetime('now')); %% posixtime, i.e. seconds
accelerationFactor=1000;
maxOrbits=16;

%% data that will later be per satellite and therefore inside SPMD loop
orbitSection    =1;         %degree
orbitSections   =[1:orbitSection:360];
%masterSatellite =0;
orbitCounter    =0;  
etemp           =zeros(6,1);
sst             =zeros(6,1);
sstDesired      =zeros(6,1);
sstOld          =zeros(6,1);
refPosChangeTemp=zeros(3,1);

spmd(3) %% create satellite instances
  send(DQ,strcat('No. ',num2str(labindex),' is alive.'))
  go=getMode(maxOrbits,orbitCounter,DQ);
  
  while go==1  
    orbitCounter=orbitCounter+1;
    send(DQ,strcat(num2str(labindex),': ------ no of orbit: ',num2str(orbitCounter)));
    startOrbit=now; %% posixtime, i.e. seconds
    %send(DQ,strcat(num2str(labindex),': orbittimer begin: ',num2str(posixtime(datetime('now'))-startTime)));

    [meanMotion,meanAnomalyFromAN]=WhereInWhatOrbit(accelerationFactor);
    
    %% wait until end of orbit section
    idx = find(orbitSections >= meanAnomalyFromAN,1,'first');
    idx=idx+1;
    pause((orbitSections(idx)-meanAnomalyFromAN)/meanMotion);
    
    while go==1 && idx<=size(orbitSections,2) %% orbit sections loop
      start=now; %% determine cycle start time in order allow to subtract cycle duration from waiting
      %% set attitude computed in last iteration
      setAttitude();
      %% compute attitude for next section
      %% determine desired trajctory
      %sstDesired=sstDesiredFunction(meanAnomalyFromAN/MeanMotion,MeanMotion);
      etemp(1:6)=sst(1:6)-sstDesired(1:6);
      if not(masterSatellite)
        %! ISL error
        %! compute average error
        %! ISL averageerror
        %! assign average error, i.e. compute final error
      end
      computeAttitude();
      if not(masterSatellite)
        if labindex==1
          refPosChangeTemp(1:3)=sst(1:3)-sstOld(1:3);
        end
        %! move coordinate system
      end
      %% set counter and go condition
      idx=idx+1;
      go=getMode(maxOrbits,orbitCounter,DQ);      
      %% wait until next section starts
      pause(orbitSection/meanMotion-(now-start));
    end %% orbit sections while loop
    
    %send(DQ,strcat(num2str(labindex),': orbittimer end: ',num2str(posixtime(datetime('now'))-startTime)));  
    send(DQ,strcat(num2str(labindex),': duration of orbit: ',num2str(now-startOrbit)));  
  end %% orbit while loop
  
  send(DQ,strcat('No.',num2str(labindex),' is dead.'))
end
delete(gcp)

function setAttitude()
  pause(0.001);
end
function computeAttitude()
  pause(0.001);
end

function [meanMotion,meanAnomalyFromAN]=WhereInWhatOrbit(accelerationFactor)
%% this function usually provides a meanAnomalyFromAN=0 and the related meanMotion for a circular orbit

  altitude=600;  
  orbitalPeriod=90*60/accelerationFactor;
  argumentOfPerigee=90;
  meanAnomaly=30;  
  meanMotion=360/orbitalPeriod;    %
  meanAnomalyFromAN=argumentOfPerigee+meanAnomaly;             %
 %(360-meanAnomalyFromAN-5)/meanMotion
end


function go=getMode(maxOrbits,orbitCounter,DQ)
  go=1;
  %% readfile
  readModeFromFile=1;
  if orbitCounter>=maxOrbits
    send(DQ,strcat(num2str(labindex),'leaving loop - maximum number of orbits reached'))
    go=10;
  elseif readModeFromFile==0
    send(DQ,strcat(num2str(labindex),'leaving loop - filestop'))
    go=20;
  else
    go=1;
  end
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
