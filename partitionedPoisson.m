%% Partitioning a Poisson Process
% Here we assume that a Poisson Process can be divided into two types of
% states both of which can occur with probability p and 1-p. We would like
% to simulate such a proces. 
clc
clear all
close all
% Initiating the values
T = 50;
t = 0;
t1 = 0; 
t2 = 0;
N1(1) = 0;
N2(1) = 0;
Lm1 = 1;
Lm2 = 0.5;
Lm = Lm1 + Lm2;
p = Lm1/Lm;
i = 2;
j = 2;
tN(1) = 0;
tN(1) = 0;
%Loop for Simulation
while t<T
    U = rand();
    t = t + (-1/Lm)*log(U);
    if rand()<=p
        N1(i) = N1(i-1) + 1;
        tN1(i) = t;
        i = i + 1;
    else
        N2(j) = N2(j-1) + 1;
        tN2(j) = t;
        j = j + 1;
    end
end
%Plotting One Sample Path Each
stairs(tN1,N1,'color','red','LineWidth',1.5);
hold on
stairs(tN2,N2,'color','green','LineWidth',1.5);
hold on
plot(0:1:T,0:1:T,'LineWidth',2,'color','blue')
title("Sample Path for Partitioned Poisson processes with $\lambda_1$ and $\lambda_2$", 'Interpreter','latex');
grid on
legend('PP(\lambda_1)','PP(\lambda_2)');
xlim([0 25])

%% Plotting Multiple Sample Paths
figure
T = 1000;
plot(0:T/1e+2,0:T/1e+2,'LineWidth',3,'DisplayName','$x = y$','color','blue')
hold on
stairs(tN1,N1,'DisplayName','$\lambda = 1$','color','red');
hold on
stairs(tN2,N2,'DisplayName','$\lambda = \frac{1}{2}$','color','green');
hold on
h2 = legend('show');
set(h2,'Interpreter','latex','AutoUpdate','Off')
xlabel('t','Interpreter','latex')
ylabel('$N_1(t)$, $N_2(t)$','Interpreter','latex')
for i = 1:100
    [tN1,N1,tN2,N2] = partPoisson(T, Lm1, Lm2);
    stairs(tN1,N1,'color','red')
    hold on
    stairs(tN2, N2,'color','green')
    hold on
end
plot(0:T/1e+2,0:T/1e+2,'LineWidth',3,'DisplayName','$x = y$','color','blue')
hold on
xlim([0 T/1e+2])
grid on
title('Sample Paths for Partitioned Poisson processes with $\lambda_1$ and $\lambda_2$', 'Interpreter','latex')
% suptitle('Stair Step Plot of 100 Sample Poisson Processes with various \lambda');



%% Function for 2-partition Poisson Process
function [N1,tN1,N2,tN2] = partPoisson(T, Lm1, Lm2)
    t = 0;
    t1 = 0; 
    t2 = 0;
    N1(1) = 0;
    N2(1) = 0;
    Lm = Lm1 + Lm2;
    p = Lm1/Lm;
    i = 2;
    j = 2;
    tN(1) = 0;
    tN(1) = 0;
    while t<T
    U = rand();
    t = t + (-1/Lm)*log(U);
        if rand()<=p
            N1(i) = N1(i-1) + 1;
            tN1(i) = t;
            i = i + 1;
        else
            N2(j) = N2(j-1) + 1;
            tN2(j) = t;
            j = j + 1;
        end
    end
end