function e_t = ElastanceDriver(t,Tvc,Tvr)
%Input the following variables:
%t = time from 0 to T
%T = heart rate

if t<=Tvc
   e_t = (1-cos(pi*t/Tvc))/2;
elseif t <= Tvr
   e_t = (1+cos(pi*(t-Tvc)/(Tvr-Tvc)))/2 ;
else 
   e_t = 0;
end

