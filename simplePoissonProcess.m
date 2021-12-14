%% Simulation of a Simple Poisson Process
clc
clear all 
close all

T = 2000; %Length of the Simulation
t = 0;
N(1) = 0;
i = 2;
Lm = 1;
tN(1) = 0;
% Loop for Simulation
while t<T
    U = rand();
    t = t + (-1/Lm)*log(U);
    if t>T
        break
    end
    N(i) = N(i-1) + 1;
    tN(i) = t;
    i = i +1;
end

%% Plotting the Inter-arrival time distribution
X(1) = tN(1);
for i = 2:length(tN)
   X(i) = tN(i) - tN(i-1) ;
end

tempdist = makedist('Exponential',Lm);
n = length(X) - 1;
figure
% subplot(211)
histogram(X,'Normalization','pdf','DisplayName','Inter-Arrival Times')
hold on
plot(0:n,pdf(tempdist,0:n),'DisplayName','True Exponential PDF','LineWidth',2)
xlim([0 8])
grid on 
title('Distribution of Inter-Arrival Times with $\lambda = 1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X = x)$','Interpreter','latex')
hl = legend('show');
set(hl, 'Interpreter','latex')

%% Plotting the Counting Process 
figure 
stairs(tN,N, 'LineWidth',1.75, 'DisplayName','Sample Path')
hold on
plot(0:15,0:15,'LineWidth',2,'DisplayName','$x = y$')
xlim([0 15])
grid on
title('Stair Step Plot of a Sample Poisson Process with $\lambda = 1$', 'Interpreter','latex')
h2 = legend('show');
set(h2,'Interpreter','latex')
xlabel('t','Interpreter','latex')
ylabel('$N(t)$','Interpreter','latex')

%% Change in Sample Path Trajectory as Lm increases
%Case - 1: Lm = 1
LmTst = [1 2 0.5];
figure
subplot(131)
plot(0:T/1e+2,0:T/1e+2,'LineWidth',3,'DisplayName','$x = y$')
hold on
h2 = legend('show');
set(h2,'Interpreter','latex','AutoUpdate','Off')
xlabel('t','Interpreter','latex')
ylabel('$N(t)$','Interpreter','latex')
for i = 1:100
    [tN,N] = poisson_sim(LmTst(1),T);
    stairs(tN,N)
    hold on
end
xlim([0 T/1e+2])
grid on
title('$\lambda = 1$', 'Interpreter','latex')

subplot(132)
plot(0:T/1e+2,0:T/1e+2,'LineWidth',3,'DisplayName','$x = y$')
hold on
h2 = legend('show');
set(h2,'Interpreter','latex','AutoUpdate','Off')
xlabel('t','Interpreter','latex')
ylabel('$N(t)$','Interpreter','latex')
for i = 1:100
    [tN,N] = poisson_sim(LmTst(2),T);
    stairs(tN,N)
    hold on
end
xlim([0 T/1e+2])
grid on
title('$\lambda = 2$', 'Interpreter','latex')


subplot(133)
plot(0:T/1e+2,0:T/1e+2,'LineWidth',3,'DisplayName','$x = y$')
hold on
h2 = legend('show');
set(h2,'Interpreter','latex','AutoUpdate','Off')
xlabel('t','Interpreter','latex')
ylabel('$N(t)$','Interpreter','latex')
for i = 1:100
    [tN,N] = poisson_sim(LmTst(3),T);
    stairs(tN,N)
    hold on
end
xlim([0 T/1e+2])
grid on
title('$\lambda = \frac{1}{2}$', 'Interpreter','latex')
suptitle('Stair Step Plot of 100 Sample Poisson Processes with various \lambda');


%% Function to simulate a Poisson Process
function [tN,N] = poisson_sim(Lm,T)
    t = 0;
    N(1) = 0;
    i = 2;
    tN(1) = 0;
    while t<T
        U = rand();
        t = t + (-1/Lm)*log(U);
        N(i) = N(i-1) + 1;
        tN(i) = t;
        i = i +1;
    end
end