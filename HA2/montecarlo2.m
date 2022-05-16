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


deltas=abs(diff(Iran_susceptible));
deltas=[Iran_population-Iran_susceptible(1,1);deltas];


%% Starting with Iran
deltaI=Iran_infected(1,1);
deltaS=Iran_population-Iran_susceptible(1,1)-Iran_infected(1,1);
sum1=0; sum2=0;
pirvec=zeros(length(Iran_infected),1);
beta=3;
sigma=0.01; M=2; 



%Attempt at Random-Gaussian walk

f = @(lambda,Y,kappa) alpha*log(beta)-gammaln(alpha)+(alpha-1)*log(lambda)-beta*lambda+sum(gammaln(kappa+Y)-gammaln(Y+1)-gammaln(kappa)+kappa*log(1-phi)+Y*log(phi));
g = @(lambda,Y,kappa) sum(gammaln(kappa+Y)-gammaln(Y+1)-gammaln(kappa)+kappa*log(1-phi)+Y*log(phi)); %Draw for t

%Given ammount of breakpoints we give the initial variables.
breakpoint=3;

breakpoints=round(linspace(0,length(Iran_infected),breakpoint+2));
breakpoints(1,1)=1;
lambdas=ones(1,breakpoint+1);


parameters=zeros(breakpoint*2+2,1); %We need one pir, n breakpoints and n+1 lambdas.
N=20000;
breakpointsss=zeros(N,breakpoint);
for i=1:N
    placeholder=0;
    placeholder2=0;
    j=1;
    for k=1:length(parameters)
        epsilon1=normrnd(0,1); epsilon2=randi([-M,M]) ;
        %Investigate lambda

        if k<=breakpoint+1
            

            psi=1-exp(-lambdas(1,k)*Iran_infected(breakpoints(1,k):breakpoints(1,k+1),:)/Iran_population);
            kappa1=(1/phi-1)*Iran_susceptible(breakpoints(1,k):breakpoints(1,k+1),:).*psi;
            lambdastar=lambdas(1,k)+sigma*epsilon1; %drawn
            
            if lambdastar<0
                lambdastar=lambdas(1,k);
            end
            
            psi=1-exp(-lambdastar*Iran_infected(breakpoints(1,k):breakpoints(1,k+1),:)/Iran_population);
            kappa2=(1/phi-1)*Iran_susceptible(breakpoints(1,k):breakpoints(1,k+1),:).*psi;

            deltas_temp=deltas(breakpoints(1,k):breakpoints(1,k+1),:);
            L=(placeholder+(f(lambdastar,deltas_temp,(kappa2))))/(placeholder+(f(lambdas(1,k),deltas_temp,(kappa1))));
            prob=min(1,exp(placeholder+(f(lambdastar,deltas_temp,(kappa2)))-(placeholder+(f(lambdas(1,k),deltas_temp,(kappa1))))));%Set alpha
            U=rand(1); % Drawn uniform 0,1

            if U<=prob
                lambdas(1,k)=lambdastar;
                placeholder=placeholder+f(lambdastar,deltas_temp,kappa2);

            else
                placeholder=placeholder+f(lambdas(1,k),deltas_temp,kappa1);

            end
        
        %Investigate t
        
        elseif breakpoint+1<k & k<=breakpoint*2+1
            %disp(j)
            
            tstar=breakpoints(1,j+1)+epsilon2;
            if breakpoints(1,j)>tstar | breakpoints(1,j+2)<tstar;
                tstar=breakpoints(1,j+1);
            end
            
            psi=1-exp(-lambdas(1,j)*Iran_infected(breakpoints(1,j):breakpoints(1,j+1),:)/Iran_population);
            kappa1=(1/phi-1)*Iran_susceptible(breakpoints(1,j):breakpoints(1,j+1),:).*psi;
            
            
            psi=1-exp(-lambdas(1,j)*Iran_infected(breakpoints(1,j):tstar,:)/Iran_population);
            kappa2=(1/phi-1)*Iran_susceptible(breakpoints(1,j):tstar,:).*psi;
            
            deltas_temp1=deltas(breakpoints(1,j):breakpoints(1,j+1),:);
            deltas_temp2=deltas(breakpoints(1,j):tstar,:);
            
            Iter1=g(lambdas(1,j),deltas_temp2,kappa2); Iter2=(g(lambdas(1,j),deltas_temp1,(kappa1)));
            
            psi=1-exp(-lambdas(1,j+1)*Iran_infected(breakpoints(1,j+1):breakpoints(1,j+2),:)/Iran_population);
            kappa1=(1/phi-1)*Iran_susceptible(breakpoints(1,j+1):breakpoints(1,j+2),:).*psi;
            
            
            psi=1-exp(-lambdas(1,j+1)*Iran_infected(tstar:breakpoints(1,j+2),:)/Iran_population);
            kappa2=(1/phi-1)*Iran_susceptible(tstar:breakpoints(1,j+2),:).*psi;
            
            deltas_temp1=deltas(breakpoints(1,j+1):breakpoints(1,j+2),:);
            deltas_temp2=deltas(tstar:breakpoints(1,j+2),:);
            Iter1=Iter1+g(lambdas(1,j),deltas_temp2,kappa2); Iter2=Iter2+(g(lambdas(1,j),deltas_temp1,(kappa1)));
            
            prob=min(1,exp(Iter1-Iter2));%Set alpha
            U=rand(1); % Drawn uniform 0,1
            
            if U<=prob
                breakpoints(1,j+1)=tstar;
                
            else
                breakpoints(1,j+1)=breakpoints(1,j+1);
            end
            breakpointsss(i,j)=tstar;
            j=j+1;
            
        end
        %Investigate pir
        %x=deltaI+deltaS;
        %sumx=sumx+x;
        %sumit=sumit+Iran_infected(k,1);
        %pir=betarnd(sumx-a,sumit-sumx+b); %Gibbs sampler,?? Or do we sample x too?
        %pirvec(k,1)=pir;

        %Random walk
        %alpha=min(1,f(xstar)/f(xk));



        deltaI=-Iran_infected(k+1,1)+Iran_infected(k,1);
        deltaS=-Iran_susceptible(k+1,1)+Iran_susceptible(k,1);

    end
    
end
%plot(breakpointsss)
