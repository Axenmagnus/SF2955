clear all; close all; clc
%Part 1
%rng(19)
sigma=0.5; deltaT=0.5; alpha = 0.6;
P=1/20*[[16,1,1,1,1];[1,16,1,1,1];[1,1,16,1,1];[1,1,1,16,1];[1,1,1,1,16]];
zn=[[0,0];[3.5,0];[0,3.5];[0,-3.5];[-3.5,0]];
theta=[[1,deltaT,deltaT^2/2];[0,1,deltaT];[0,0,alpha]];
phiz=[[deltaT^2/2];[deltaT];[0]];
phiw=[[deltaT^2/2];[deltaT];[1]];
N=6;
zero=zeros(N,1);

N=6;
sigmamatrix=diag([500,5,5,200,5,5]);
%sigmamatrix=500*diag(ones(N,1),0) + 5*diag(ones(N-1,1),1) + 5*diag(ones(N-1,1),-1)+ 5*diag(ones(N-2,1),-2)+ 5*diag(ones(N-2,1),2)+ 200*diag(ones(N-3,1),3)+ 200*diag(ones(N-3,1),-3)+ 5*diag(ones(N-4,1),-4)+ 5*diag(ones(N-4,1),4)+ 5*diag(ones(N-5,1),5)+ 5*diag(ones(N-5,1),-5);

theta=[[theta,zeros(3,3)];[zeros(3,3),theta]];
phiz=[[phiz,zeros(3,1)];[zeros(3,1),phiz]];
phiw=[[phiw,zeros(3,1)];[zeros(3,1),phiw]];

X0=mvnrnd(zero,sigmamatrix);
Command=randi([1 5],1,1);
Z0=zn(Command,:);


m=3000; % Timesteps
States=zeros(m,2);
Xi=transpose(X0);
Zi=transpose(Z0);
Wi=transpose(mvnrnd([0,0],sigma^2*eye(2)));
for i=1:m
    
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    
    %Need to update Zi,
    Update=(rand(1));
    summ=0;
    for j=1:5
        
        summ=P(Command,j)+summ;
        if Update<summ
            Command=j;
            break
        end
        
    end
    Zi=transpose(zn(Command,:));
    States(i,1)=Xi(1);
    States(i,2)=Xi(4);
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2)));
    
    
end
%plot(States(:,1),States(:,2))
%title('Randomized path, seed 19, m=3000')
%xlabel('X1 location') 
%ylabel('X2 location')

%% Making of a nxN state matrix for use in part 3 
v=90; zeta=1.5; gamma=3; 

load("stations.mat");

load("RSSI-measurements.mat");

N=10000;
n=length(Y);
Command=randi([1 5],1,N);
Commands=zeros(N,n);
Commands(:,1)=Command;
for j=1:N
    
    for i =1:n
        Update=(rand(1));
        summ=0;
        for k=1:5

            summ=P(round(Commands(j,i)),k)+summ;
            if Update<=summ
                Command=k;
                break
            end

        end
        Commands(j,i+1)=Command;
         
    end
    
end

%% Part 3



tau=zeros(6,n);
Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
w = (p(Xi(1,:),Xi(4,:),Y(:,1)));
tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

Wi=transpose(mvnrnd([0,0],sigma*eye(2),N));


Zi=transpose(zn(Commands(:,1),:));
for k = 1:500 % main loop
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    w =w.*(p(Xi(1,:),Xi(4,:),Y(:,k+1))); % weighting REMOVE/INCLUDE LOG
    
    %w=w/sum(w);
    Zi=transpose(zn(Commands(:,k+1),:));
    Wi=transpose(mvnrnd([0,0],sigma*eye(2),N));
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
title("SIS, N=10000")
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')
%figure
%hist(w)
%ESSM=1/(sum((w/sum(w)).^2))

%% Part 3 Instructed weight adjustment



tau=zeros(6,n);
Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
w = log(p(Xi(1,:),Xi(4,:),Y(:,1)));
L=max(w);
w=exp(w-L);
w=w/sum(w);
tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));


Zi=transpose(zn(Commands(:,1),:));
for k = 1:500 % main loop
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    w =log(w.*p(Xi(1,:),Xi(4,:),Y(:,k+1))); % weighting REMOVE/INCLUDE LOG
    L=max(w);
    w=exp(w-L);
    w=w/sum(w);
    Zi=transpose(zn(Commands(:,k+1),:));
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
title("SIS, N=10000")
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')
%figure
%title("Histogram n=70")
%loglog(hist((w))
%ESSM=1/(sum((w/sum(w)).^2))
%% Part 4

NewCommands=zeros(N,n);
tau=zeros(2,n);

Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
w = p(Xi(1,:),Xi(4,:),Y(:,1));
tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

weights=zeros(n,N);
weights(1,:)=w;
st=Commands(:,1);
Zi=transpose(zn(Commands(:,1),:));
cs=cumsum(P,2);
NewCommands(:,1)=st;
for k = 1:500 % main loop
    ind=randsample(N,N,true,w);
    Xi=Xi(:,ind);
    %Resample states
    %Zi=Zi(ind);
    %Commands=Commands(ind,:);
    Rand=rand(N,1);
    StateArray=zeros(N,1);
    for j=1:5
        vals=cs(st,j);
        indexarray=Rand<=vals;
        StateArray=StateArray+j*indexarray;
        Rand(indexarray)=1.2;

        
        
    end
    st=StateArray;
    NewCommands(:,k+1)=st;
    Zi=transpose(zn(StateArray,:));
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    w =(p(Xi(1,:),Xi(4,:),Y(:,k+1))); % weighting REMOVE/INCLUDE LOG
    w=w/sum(w);
    weights(k+1,:)=w;
    
    
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
title("SISR N=10000")
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')


%% Most probable driving command
drivingCommand=zeros(1,n);
for j=1:n % Go through the time-steps
    
    probability=zeros(1,5);
    for i=1:5%Go through the driving commands
        I = find(NewCommands(:,j) == i);
        vals=weights(j,I);
        probability(1,i)=sum(vals);
    end
    [val, idx] = max(probability);
    drivingCommand(1,j)=idx;
    
end
%drivingCommand
hist(drivingCommand)

%% Plot driver commands

hold on
for i=1:n
    
    if drivingCommand(i)==1
        plot(tau(1,i),tau(2,i),"cx")
    elseif drivingCommand(i)==2
        plot(tau(1,i),tau(2,i),"gx")
    elseif drivingCommand(i)==3
        plot(tau(1,i),tau(2,i),"rx")
    elseif drivingCommand(i)==4
        plot(tau(1,i),tau(2,i),"bx")
    elseif drivingCommand(i)==5
        plot(tau(1,i),tau(2,i),"cx")
    end
    
end
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')


%% Part 5 
clc

load("RSSI-measurements-unknown-sigma.mat")

plot(pos_vec(1,:), pos_vec(2,:),'d')
hold on

T = 10;


cn = zeros(T,1);
sigmaV = linspace(2, 2.3, T);


for h = 1:T

    omega = zeros(n,1);
    NewCommands=zeros(N,n);
    tau=zeros(2,n);

    Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

    %p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
    w = p5(Xi(1,:),Xi(4,:),Y(:,1),sigmaV(h),pos_vec);
    omega(1) = sum(w);

    tau(1,1) = sum(Xi(1,:).*w)/sum(w);
    tau(2,1) = sum(Xi(4,:).*w)/sum(w);

    weights=zeros(n,N);
    weights(1,:)=w;
    st=Commands(:,1);
    Zi=transpose(zn(Commands(:,1),:));
    cs=cumsum(P,2);
    NewCommands(:,1)=st;

    for k = 1:n-1 % main loop
        ind=randsample(N,N,true,w);
        Xi=Xi(:,ind);
        %Resample states
        %Zi=Zi(ind);
        %Commands=Commands(ind,:);
        Rand=rand(N,1);
        StateArray=zeros(N,1);
        for j=1:5
            vals=cs(st,j);
            indexarray=Rand<=vals;
            StateArray=StateArray+j*indexarray;
            Rand(indexarray)=1.2;
        end
        st=StateArray;
        NewCommands(:,k+1)=st;
        Zi=transpose(zn(StateArray,:));
        Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
        Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
        w=p5(Xi(1,:),Xi(4,:),Y(:,k+1),sigmaV(h),pos_vec); % weighting REMOVE/INCLUDE LOG
        omega(k+1) = sum(w);
        w=w/sum(w);
        weights(k+1,:)=w;


        tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
        tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
    end

    plot(tau(1,:),tau(2,:))
    title("SISR N=10000")

    cn(h) = sum(log(omega/N));

%     plot(tau(1,:),tau(2,:))
%     hold on
end

[wmax,wind] = max(cn);
sigmaMax = sigmaV(wind)
