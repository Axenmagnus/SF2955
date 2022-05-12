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

sigma=5; M=3; 



%Attempt at Random-Gaussian walk
%prod(beta^alpha/factorial(alpha-1)*lambda^(alpha-1)*exp(-beta*lambda))*prod(nbinpdf(diff(y),kappa,phi));
f = @(lambda,Y,kappa) (alpha-1)*log(lambda)-beta*lambda+sum(log((nbinpdf(Y,kappa,phi))));
posterrior_t = @(y,kappa,phi,It,Pir) prod(nbinpdf((diff(y)),kappa,phi)*binpdf((diff(y)),It,Pir));


%Given ammount of breakpoints we give the initial variables.
breakpoint=1;

breakpoints=round(linspace(0,length(Iran_infected),breakpoint+2));
breakpoints(1,1)=1;
lambdas=10*ones(1,breakpoint+1);


parameters=zeros(breakpoint*2+2,1); %We need one pir, n breakpoints and n+1 lambdas.

placeholder=0;
for k=1:length(parameters)
    epsilon1=normrnd(0,1); epsilon2=randi([-3,3]) ;
    %Investigate lambda
    
    if k<=breakpoint+1
        
        psi=1-exp(-lambdas(1,k)*Iran_infected(breakpoints(1,k):breakpoints(1,k+1),:)/Iran_population);
        kappa1=(1/phi-1)*Iran_susceptible(breakpoints(1,k):breakpoints(1,k+1),:).*psi;
        lambdastar=lambdas(1,k)+sigma*epsilon1; %drawn
        psi=1-exp(-lambdastar(1,k)*Iran_infected(breakpoints(1,k):breakpoints(1,k+1),:)/Iran_population);
        kappa2=(1/phi-1)*Iran_susceptible(breakpoints(1,k):breakpoints(1,k+1),:).*psi;
        
        deltas_temp=deltas(breakpoints(1,k):breakpoints(1,k+1),:);
        
        
        
        prob=min(1,(placeholder+f(lambdastar,deltas_temp,round(kappa2)))/(placeholder+f(lambdas(1,k),deltas_temp,round(kappa1)))); %Set alpha
        U=unif(0,1); % Drawn uniform 0,1
        
        if U<=prob
            lambdas(1,k)=lambdastar;
            placeholder=placeholder+f(lambdastar,breakpoints(k+1),kappa);
        
        else
            placeholder=placeholder+f(lambdas(1,k),deltas,kappa);
        
        end
    end
    %Investigate 
    if breakpoint+1<k<=breakpoint*2+1
    
        
    end

    x=deltaI+deltaS;
    sumx=sumx+x;
    sumit=sumit+Iran_infected(k,1);
    pir=betarnd(sumx-a,sumit-sumx+b); %Gibbs sampler,?? Or do we sample x too?
    pirvec(k,1)=pir;
    
    %Random walk
    %alpha=min(1,f(xstar)/f(xk));
    
    
    
    deltaI=-Iran_infected(k+1,1)+Iran_infected(k,1);
    deltaS=-Iran_susceptible(k+1,1)+Iran_susceptible(k,1);
    
end
