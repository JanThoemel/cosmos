function [sstDesired]=IvanovFormationFlightDesired(timetemptemp,MeanMotion,i)
%% desired solution for Ivanov
  sstDesired=zeros(9,1,size(timetemptemp,2));
  %% analytical solution according to Ivanov
  AAA=100;    DDD=115;
  switch i
    case 1
      sstDesired(1,1,:)=-DDD;
    case 2
      sstDesired(1,1,:)=DDD;
    case 3
      sstDesired(1,1,:)=2*AAA*        cos(MeanMotion*(timetemptemp)-acos(1/3));  
      sstDesired(2,1,:)=  AAA*sqrt(3)*sin(MeanMotion*(timetemptemp));
      sstDesired(3,1,:)=  AAA*        sin(MeanMotion*(timetemptemp)-acos(1/3));
      sstDesired(4,1,:)=2*AAA*       -sin(MeanMotion*(timetemptemp)-acos(1/3))*MeanMotion;  
      sstDesired(5,1,:)=  AAA*sqrt(3)*cos(MeanMotion*(timetemptemp))*MeanMotion;
      sstDesired(6,1,:)=  AAA*        cos(MeanMotion*(timetemptemp)-acos(1/3))*MeanMotion;
    case 4
      sstDesired(1,1,:)=2*AAA*        cos(MeanMotion*(timetemptemp));
      sstDesired(2,1,:)=  AAA*sqrt(3)*sin(MeanMotion*(timetemptemp)+acos(1/3));
      sstDesired(3,1,:)=  AAA*        sin(MeanMotion*(timetemptemp));
      sstDesired(4,1,:)=2*AAA*       -sin(MeanMotion*(timetemptemp))*MeanMotion;
      sstDesired(5,1,:)=  AAA*sqrt(3)*cos(MeanMotion*(timetemptemp)+acos(1/3))*MeanMotion;
      sstDesired(6,1,:)=  AAA*        cos(MeanMotion*(timetemptemp))*MeanMotion;
  end
end