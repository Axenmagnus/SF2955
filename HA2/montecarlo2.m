clear all; close all; clc
%Constants and instantiating
Iran_infected=readtable("iran_infected.csv"); Iran_recovered=readtable("iran_removed.csv");
Germany_infected=readtable("germany_infected.csv");Germany_recovered=readtable("germany_removed.csv");
Germany_population=83000000; Iran_population=84000000;
Iran_infected=Iran_infected{:,:}; Iran_recovered=Iran_recovered{:,:};
Germany_infected=Germany_infected{:,:}; Germany_recovered=Germany_recovered{:,:};

Iran_susceptible=ones(height(Iran_infected),1)*Iran_population;
Iran_susceptible=Iran_susceptible-Iran_infected-Iran_recovered;

Germany_susceptible=ones(height(Germany_infected),1)*Germany_population;
Germany_susceptible=Germany_susceptible-Germany_infected-Germany_recovered;
phi=0.995; alpha=2;
a=2; b=3;



%% Starting with Iran
deltaI=Iran_infected(1,1);
deltaS=Iran_population-Iran_susceptible(1,1)-Iran_infected(1,1);
sum1=0; sum2=0;
pirvec=zeros(length(Iran_infected),1);
for k=1:height(Iran_infected)-1
    
    x=deltaI+deltaS;
    
    pir=betarnd(x+1,Iran_infected(k,1)-deltaS+1); %Gibbs sampler,?? Or do we sample x too?
    pirvec(k,1)=pir;
    
    
    
    
    deltaI=-Iran_infected(k+1,1)+Iran_infected(k,1);
    deltaS=-Iran_susceptible(k+1,1)+Iran_susceptible(k,1);
    
end
