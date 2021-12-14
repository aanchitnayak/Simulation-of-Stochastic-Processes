%% Simulating the Poisson Process 

clc
clear all
close all

%% The Poisson Process
% A Poisson Process with rate Lm is a renewal process in which the
% interarrival time distribution is exponential with rate Lm. 
T = 10;
Lm = 1;
tN = poisson_sim(Lm,T);
%Plotting
stairs(tN)
title("Poisson Process with rate \lambda = 1")
grid on
%%
simLength = 50;
figure
for i = 1:simLength
    stairs(poisson_sim(1,T))
    hold on
end
plot(0:15,0:15,'LineWidth',2,'DisplayName','cos(3x)')
hold on
grid on
title("Poisson Process with rate \lambda = 1 with 50 sample Paths")
xlim([0,14])
ylim([0,14])
legend()

%Note that as Lm increases, the sample paths tend to go below the x = y
%line

%% Partitioning a Poisson Process
% Here we assume that a Poisson Process can be divided into two types of
% states both of which can occur with probability p and 1-p. We would like
% to simulate such a proces. 
t = 0;
t1 = 0; 
t2 = 0;
N1 = 0;
N2 = 0;
Lm1 = 1.25;
Lm2 = 1.75;
Lm = Lm1 + Lm2;
p = Lm1/Lm;
i = 1;
while t<T
U = rand();
t = t + (-1/Lm)*log(U);
if rand()<=p
    N1 = N1 + 1;
    tN1 = t;
    X1(i,1) = tN1;
    X1(i,2) = N1;
else
    N2 = N2 + 1;
    tN2 = t;
    X2(i,1) = tN2;
    X2(i,2) = N2;
end
i = i+1;
end
%Plotting
stairs(X1(:,1));
hold on
stairs(X2(:,1));
title("Sample Path for Partitioned Poisson processes with \lambda_1 + \lambda_2 = 1");
grid on
legend('PP(\lambda_1)','PP(\lambda_2)');

%% Non-Stationary Poisson Process
% This is a case when Lm is a function of time. Such a Poisson Process is
% called a Non-Stationary Poisson Process (NSPP). 
t = 0;
N = 0;
i = 1;
LmStar = 1;
T = 20;

while t<T
U = rand();
t = t + (-1/LmStar)*log(U);
U = rand();
if U <= Lmfxn(t,T)/LmStar
%     N(i+1) = N(i) + 1;
    N = N+1;
    tN(i) = t;
end
i = i + 1;
end

figure
stairs(tN)
hold on
plot(1:length(tN),Lmfxn(t,T),'LineWidth',3)

%% Compound Poisson Processes
% Here, we assume that the occurence of an event is following a poisson
% process with rate lambda and that each of occurence itself follows a
% custom probability distribution. This is a Compound Poisson Process. This
% can be used to model batch based arrivals in a system. 

%For the custom distribution we will be using a Normal Probability
%Distribution. 

t = 0;
N(1) = 0;
X(1) = 0;
G = makedist('Gamma','a',9,'b',0.5);
i = 1;
B(1) = random(G);
T = 10;
Lm = 2;
while t<=T
U = rand();
t = t + (-1/Lm)*log(U);
if t>T
    break
end
B(i+1) = random(G);
N(i+1) = N(i) + 1;
X(i+1) = X(i) + B(i);
tN(i) = t;
i = i+1;
end
data = [X' N' B']; 
% figure

stairs(X)
hold on

%% Plotting Exercise for Real time
t = 1:max(tN);
for i = 1:length(t)
    if t(i) <= tN(i)
        Xnew(i) = X(i)
    end
end
%% Function to simulate a Poisson Process
function Y = poisson_sim(Lm,T)
    t = 0;
    N = 0;
    i = 1;
    while t<T
        U = rand();
        t = t + (-1/Lm)*log(U);
        N = N + 1;
        Y(i) = t;
        i = i +1;
    end
end
%% Function for Non-Stationary Poisson Process
function y = Lmfxn(t,T)
    
% %     T = 10;
    if t<=T
        y = 1;
    else
        y = 5;
    end

%     y = sin(2*pi*t/T);
end