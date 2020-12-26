%%%%%%%%%%%%%%%%%%%%%%%%%%% Main File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m_a m_b c s total_number_of_equations number_of_equations_per_phase;
m_a=4;
m_b=3;
c=10;
s=5;
total_number_of_equations=((m_a/m_b)*(s+m_b*(c-s+1)))*nchoosek(m_b+s-1,s);
number_of_equations_per_phase=total_number_of_equations/m_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ics=zeros(total_number_of_equations+4*m_a+3*m_b*m_a+2*m_b*m_b*m_a,1);
ics(1,1)=1;
ics(total_number_of_equations+1,1)=1;
opts=odeset('Reltol',1e-6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time=0;
end_time=50;
duration=end_time-start_time;
tspan=linspace(start_time,end_time,duration*2000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,prob]=ode45(@PhPhsc_function,tspan,ics,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dummy=prob(:,1:total_number_of_equations);
E0=sum(dummy,2);
%plot(t,E0);

PhPhsc_test(t,prob);
