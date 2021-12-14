%% Queueing Model

% Distribution of inter-arrival times (Tia)-> Exponential(Lm)
% Distribution of service times (S)       -> Uniform(a,b)

% Recursion for Delay-times is known. Simulate S and Tia to find D. 

clc;
clear all;
close all;

T = 10;
Lm = 1;
%% Plotting one Sample Path for a FIFO Queue
% Plot of Delays
figure
% subplot(
plot(0:100,0:100,'color','red','DisplayName','$x = y$')
hold on
[tN,Tia,N,S,D]=fifoQ(Lm,T);
stairs(tN,D,'--','color','blue','DisplayName','Sample Path')
xlim([0,T])
h = legend('show');
set(h,'Interpreter','latex','AutoUpdate','Off')
hold on
for i = 1:100
    [tN,Tia,N,S,D]=fifoQ(Lm,T);
    stairs(tN,D,'--','color','blue')
    hold on
end
plot(0:100,0:100,'color','red','LineWidth',3.5)
grid on
title('$G|G|1$ Queue Simulation - Delay Times','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$\{D_n\}$','Interpreter','latex')

figure
plot(0:100,0:100,'color','red','DisplayName','$x = y$')
hold on
[tN,Tia,N,S,D]=fifoQ(Lm,T);
stairs(tN,N,'--','color','blue','DisplayName','Sample Path')
xlim([0,T])
h = legend('show');
set(h,'Interpreter','latex','AutoUpdate','Off')
hold on
for i = 1:100
    [tN,Tia,N,S,D]=fifoQ(Lm,T);
    stairs(tN,N,'--','color','blue')
    hold on
end
plot(0:100,0:100,'color','red','LineWidth',3.5)
grid on
title('$G|G|1$ Queue Simulation - Incoming Customers','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$\{D_n\}$','Interpreter','latex')

figure
hist_S = [];
for i = 1:1000
    [tN,Tia,N,S,D] = fifoQ(Lm,T);
    hist_S = [hist_S;S'];
end
histogram(hist_S,'Normalization','pdf')
grid on
title('Distribution of Service Times', 'Interpreter','latex')
xlabel('$s$','Interpreter','latex');
ylabel('$P(S = s)$','Interpreter','latex');
figure
hist_Tia = [];
for i = 1:1000
    [tN,Tia,N,S,D] = fifoQ(Lm,T);
    hist_Tia = [hist_Tia;Tia'];
end
histogram(hist_Tia,'Normalization','pdf')
grid on
title('Distribution of Inter-Arrival Times','Interpreter','latex');
xlabel('$t_{ia}$','Interpreter','latex');
ylabel('$P(T = t_ia)$','Interpreter','latex');
%% Long Run Simulation - Estimating Long Run Delays
% We use this method to estimate the average delay faced by all customers.
MCD = [];
for i = 1:1000
    [tN,Tia,N,S,D] = fifoQ(Lm,T);
    MCD = [MCD;D'];
end
average_delay_allCustomers = mean(MCD);

% Now evaluate the proportion of customers facing this average delay 
CustProportion_for_x_delay = mean(MCD(MCD>mean(MCD)));
 
%% Function for FIFO G|G|1 Queue
function [tN,Tia,N,S,D]=fifoQ(Lm,T)

    t = 0;
    N(1) = 1;
    tN(1) = (-1/Lm)*log(rand());
    S(1) = 1 + 2*rand();
    D(1) = 0;

    i = 2;
    while t<T
        U = rand();
        t = t + (-1/Lm)*log(U);
        if t>T
            break
        end
        tN(i) = t;
        Tia(i) = (-1/Lm)*log(U);
        N(i) = N(i-1) + 1;
        S(i) = 1 + 2*rand();
        D(i) = round(max(D(i-1)+S(i-1)-Tia(i-1),0));
        i = i + 1;
    end
end