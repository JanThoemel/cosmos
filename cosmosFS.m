%% cosmosFS - Flight Simulator
%% 22/9/2019

clc;clear all;
ns=3;       %% number of satellite
DQ = parallel.pool.DataQueue;
afterEach(DQ,@disp);
accelerationFactor=10;
parpool(6);
start=posixtime(datetime('now')); %% posixtime, i.e. seconds
endTime=start+20;                 %% endtime in seconds
  

%for i=1:ns
%  parfeval(@setMode,1,start,endTime,DQ)
%end

spmd(3)
  send(DQ,strcat(num2str(labindex),'before loop'))
  go=getMode(endTime,DQ);
  while go==1
    send(DQ,strcat(num2str(labindex),'entered loop'))
    [meanMotion,meanAnomalyFromAN]=whenandwhere(accelerationFactor);

    %% wait until just before end of orbital period
    pause((360-meanAnomalyFromAN-5)/meanMotion)
    send(DQ,'here we go')
    
    a=assigna();

    %% compute average 
    ave2=a/ns;
    for i=1:ns
      if not(i==labindex)
        ave2=ave2+labSendReceive(i,i,a)/ns;
      end
    end
    send(DQ,ave2)

    %% check whether all averages are the same
    for i=1:ns %% all workers
      if not(i==labindex) %% except self
        if ave2==labSendReceive(i,i,ave2) %% is okay?
          send(DQ,'check ok')
        end
      end
    end
    %% break permanent loop
    go=getMode(endTime,DQ);
  end %% while loop
  send(DQ,strcat(num2str(labindex),'left loop'))
end
  
delete(gcp)



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
  
function [meanMotion,meanAnomalyFromAN]=whenandwhere(accelerationFactor)
  altitude=600;  
  orbitalPeriod=90*60/accelerationFactor;
  argumentOfPerigee=90;
  meanAnomaly=30;  
  meanMotion=360/orbitalPeriod;    %
  meanAnomalyFromAN=argumentOfPerigee+meanAnomaly;             %
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

function go=getMode(endTime,DQ)
  go=1;
  %% readfile
  readModeFromFile=1;
  if posixtime(datetime('now'))>endTime
    send(DQ,strcat(num2str(labindex),'leaving loop - timeout'))
    go=10;
  elseif readModeFromFile==0
    send(DQ,strcat(num2str(labindex),'leaving loop - filestop'))
    go=20;
  else
    go=1;
  end
end


function go=setMode(start,endTime,DQ)
  go=1;
end
