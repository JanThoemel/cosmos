%% cosmosFS - Flight Simulator
%% 22/9/2019


%%One week until the European CubeSat Symposium 2019.
%%The von Karman Institute for Fluid Dynamics and the University of Luxembourg,
%%the key note speakers, the jury of the Space Business Pitch, and the sponsors and
%%exhibitors welcome you - the presenters, pitchers and participants - here in Luxembourg
%%to discuss the #futureofcubesats.
%%(cubesatsymposium.eu)

%% One week after the European CubeSat Sympsosium 2019.
%% We have seen concrete advances and audaciously ambitious deepspace missions.
%% The von Karman Institute for Fluid Dynamics and the University of Luxembourg
%% presenters, pitchers and participants to joining the discussion.
%% We hope that you will come back here in Luxembourg soon and at the
%% 12th edition of the event next year in Paris.
%%(cubesatsymposium.eu)

%% On a personal note, I thank the Space Business Pitch members @ @ @ and among them those
%% from the private sector @@ who have spent their precious time and experience volunteering
%% for the good cause.

clc;clear all;
ns=3;       %% number of satellite
DQ = parallel.pool.DataQueue;
afterEach(DQ,@disp);
accelerationFactor=1000;
parpool(3);
orbitSection=1;         %degree
orbitSections=[1:orbitSection:360];

%for i=ns+1:ns+ns
%  switch labindex
%    case
%    case
%    case
%  end
%  parfeval(@setMode,1,start,endTime,DQ)
%end

startTime=posixtime(datetime('now')); %% posixtime, i.e. seconds
maxOrbits=16;

spmd(3) %% create satellite instances
  send(DQ,strcat('No. ',num2str(labindex),' is alive.'))
  orbitCounter=0;  
  go=getMode(maxOrbits,orbitCounter,DQ);
  
  while go==1  
    orbitCounter=orbitCounter+1;
    send(DQ,strcat(num2str(labindex),': ------ no of orbit: ',num2str(orbitCounter)));
    startorbit=posixtime(datetime('now')); %% posixtime, i.e. seconds
    send(DQ,strcat(num2str(labindex),': orbittimer begin: ',num2str(posixtime(datetime('now'))-startTime)));

    [meanMotion,meanAnomalyFromAN]=WhereInWhatOrbit(accelerationFactor);   
    
    %% wait until end of orbit section
    idx = find(orbitSections >= meanAnomalyFromAN,1,'first');
    pause((orbitSections(idx)-meanAnomalyFromAN)/meanMotion);
    idx=idx+1;
    
    while go==1 && idx<=size(orbitSections,2) %% orbit sections loop
      %% compute attitude for next section
      nextattitude=5;
      %% wait until next section starts
      pause(orbitSection/meanMotion);
      %% set new attitude
      setnextattitude=5;
      
      %% set counter and go condition
      idx=idx+1;
      go=getMode(maxOrbits,orbitCounter,DQ);      
    end %% orbit sections while loop
    
    send(DQ,strcat(num2str(labindex),': orbittimer end: ',num2str(posixtime(datetime('now'))-startTime)));  
  end %% orbit while loop
  
  send(DQ,strcat('No.',num2str(labindex),' is dead.'))
end
delete(gcp)


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
