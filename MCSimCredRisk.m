%% MC Sim
clc 
clear 
close all
MCLength = 1000;
T = 45;
x = 10;
Lm = 1;
c = 1;

%% Case - 1: Distribution of Claims - Uniform

G1 = makedist('Uniform',1,2);
[R1, tau1, M1, fig1] = ruinSim(x,T,Lm,c,G1,MCLength);
ruin1 = mean(tau1(isfinite(tau1)));

%% Case - 2: Distribution of Claims - Gamma

G3 = makedist('Gamma',1,3);
[R3, tau3, M3, fig3] = ruinSim(x,T,Lm,c,G3,MCLength);
ruin3 = mean(tau3(isfinite(tau3)));



%% Function for Simulating Ruin
function [R, tau, M, fig] = ruinSim(x,T,Lm,c,G,MCLength)
    R = zeros(T,MCLength);
    R(1,:) = x;
    I = zeros(1,MCLength);
    M = zeros(1,MCLength);
    
    for j = 1:MCLength
        i = 1;
        t = 0;
        tau(j) = inf; %Time of Ruin
        while t<T
            i = i + 1;
            U = rand();
            t = t + ((-1/Lm)*log(U));
            if t>T
                break
            end
            R(i,j) = R(i-1,j) + c*((-1/Lm)*log(U));
            B = random(G);
            R(i,j) = R(i,j) - B;
            if R(i,j)<0
                I(j) = 1;
                tau(j) = i;
                M(j) = abs(R(i,j));
                break
            end
        end
    end
    
    fig = figure;
    subplot(2,2,1)
    plot(R)
    title('Reserve Level with Time')
    grid on
    xlabel('Time')
    ylabel('Reserve Level in USD')
    subplot(2,2,2)
    histogram(tau,'Normalization','probability')
    title('Distribution of Time to Ruin')
    grid on
    xlabel('Time in Days')
    ylabel('Probability')
    subplot(2,2,3)
    histogram(random(G,10000,1),'Normalization','probability')
    title('Distribution of Claim Amount')
    grid on
    xlabel('Claim Amount in USD')
    ylabel('Probability')
    subplot(2,2,4)
    histogram(M,'Normalization','probability')
    title('Distribution of Magnitude of Ruin')
    grid on
    xlabel('Amount in USD')
    ylabel('Probability')
end