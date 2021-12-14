%% Non-Stationary Poisson Process
% This is a case when Lm is a function of time. Such a Poisson Process is
% called a Non-Stationary Poisson Process (NSPP). 
clc
close all
clear all

t = 0;
N(1) = 0;
LmStar = 6;
T = 20;
tN(1) = 0;
i = 2;
while t<T
U = rand();
t = t + (-1/LmStar)*log(U);
U = rand();
if U <= Lmfxn(t,T)/LmStar
    N(i) = N(i-1)+1;
    tN(i) = t;
    i = i + 1;
end

end


figure
    
    subplot(211)
    for i = 1:100
        [tN,N] = nsppSim(LmStar,T);
        stairs(tN,N, 'color','red')
        hold on
    end
    grid on
    title('NSPP Sample Paths','Interpreter','latex')
    xlabel('$t$','Interpreter','latex')
    ylabel('$N(t)$','Interpreter','latex')
    xlim([2,14])
    
    subplot(212)
    for i = 1:length(N)
        x(i) = Lmfxn(i,T);
    end
    plot(1:length(N),x,'LineWidth',1.5,'color','green')
    xlim([2,14])
    ylim([0,6])
    grid on
    title('Deterministic Poisson Rate Varying with Time','Interpreter','latex')
    xlabel('$t$','Interpreter','latex')
    ylabel('$N(t)$','Interpreter','latex')

%% Function for Non-Stationary Poisson Process Rate
function y = Lmfxn(t,T)
    
    if t<=T/2
        y = 0.25;
    else
        y = 5;
    end

end

%% Function for NSPP Simulation
function [tN,N] = nsppSim(LmStar,T)
    t = 0;
    N(1) = 0;
    tN(1) = 0;
    i = 2;
    while t<T
        U = rand();
        t = t + (-1/LmStar)*log(U);
        U = rand();
        if U <= Lmfxn(t,T)/LmStar
            N(i) = N(i-1)+1;
            tN(i) = t;
            i = i + 1;
        end

    end

end