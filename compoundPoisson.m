%% Compound Poisson Processes
clc
clear all
close all

[tN,N,B,X] = compoundPois();
figure
plot(1:0.1:10,1:0.1:10,'color','blue','LineWidth',2,'DisplayName','$x=y$')
hold on
stairs(tN,X,'DisplayName','$\lambda = 1$','color','red','DisplayName','Sample Path');
hold on
h2 = legend('show');
set(h2,'Interpreter','latex','AutoUpdate','Off')
for i = 1:100
    [tN,N,B,X] = compoundPois();
    stairs(tN, X,'color','red','LineWidth',0.5)
    hold on
end
grid on
title('Sample Paths for Compound Poisson Distribution ($\lambda = 1$)','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$X(t)$','Interpreter','latex')
%% Function for Simulating a Compound Poisson Process
function [tN,N,B,X] = compoundPois()
    t = 0;
    N(1) = 0;
    X(1) = 0;
    G = makedist('Gamma','a',7,'b',1);
    i = 2;
    B(1) = random(G);
    T = 10;
    Lm = 1;
    tN(1) = 0;
    while t<=T
        U = rand();
        t = t + (-1/Lm)*log(U);
        if t>T
            break
        end
        B(i) = random(G);
        N(i) = N(i-1) + 1;
        X(i) = X(i-1) + B(i);
        tN(i) = t;
        i = i+1;
    end
end