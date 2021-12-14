clc
clear all
close all
%% Defining the model 
Lm = 15;
T = 365-52;
S = 100;
s = 25;
storageCosts = 2;
rate = 2.5;

pdDemand = makedist('Normal',20,5);
pdLeadTime = makedist('Normal',2,0.75);

output = ClassicSsModel(rate,storageCosts,S,s,Lm,T,pdDemand,pdLeadTime);


clear ans
clc
%% Plotting
figure
plot(-50:50,pdf(pdDemand,-50:50),'Linewidth',2.5)
hold on
histogram(random(pdDemand,10000,1),'Normalization','pdf')
xlim([-10,50])
ylim([0,0.1])
grid on
title('Probability Distribution for Demand')
xlabel('Order Quantity')
ylabel('Probability')

figure
plot(-2:0.025:6,pdf(pdLeadTime,-2:0.025:6),'Linewidth',2.5)
hold on
histogram(random(pdLeadTime,5000,1),'Normalization','pdf')
% xlim([-10,50])
% ylim([0,0.1])
grid on
title('Probability Distribution for Lead Time')
xlabel('Time Taken')
ylabel('Probability')

ModelAnalysisPlotter(output)

%% Monte Carlo Simulation - Simulating Multiple Trajectories
% T = 365-52+1;
MCLength = 150;
MCX = zeros(T+1,MCLength);
MCCh = zeros(T+1,MCLength);
MCCo = zeros(T+1,MCLength);
MCR = zeros(T+1,MCLength);
MCI = zeros(T+1,MCLength);
for k = 1:MCLength
    output = ClassicSsModel(rate,storageCosts,S,s,Lm,T,pdDemand,pdLeadTime);
    MCX(:,k) = output(:,1);
    MCCh(:,k) = output(:,2);
    MCCo(:,k) = output(:,3);
    MCR(:,k) = output(:,4);
    MCI(:,k) = output(:,5);
end 
T = 365-52+1;
figure
% subplot(331)
stairs(MCR)
xlim([0,T])
grid on
xlabel('Time')
title('Revenue Generated')
% subplot(333)
figure
stairs(MCCo)
xlim([0,T])
grid on
xlabel('Time')
title('Cost of Re-Supplying Inventory')
% subplot(335)
figure
stairs(MCX)
title('Net Profit')
xlim([0,T])
grid on
xlabel('Time')
% subplot(337)
figure
stairs(MCCh)
xlim([0,T])
grid on
xlabel('Time')
title('Cost of Storage')
% subplot(339)
figure
stairs(MCI)
xlim([0,T])
grid on
xlabel('Time')
title('Inventory Level')
clc


%% Function to simulate Classic Ss Model
function output = ClassicSsModel(rate,storageCosts,S,s,Lm,T,pdDemand,pdLeadTime)
    %%Use this function to create an (S,s) Model and simulate T days' 
    %%inventory behaviour. 
    
    % Input:
    % rate       -> selling price of each item in the inventory
    % S, s       -> Model Parameters, S and s are upper and lower bounds 
    %               for theInventory's buy and sell mechanism
    % Lm         -> Poisson Process Rate for incoming buyer orders
    % T          -> Maximum Length of Simulation. If not specified, 
    %               it will run for314 days
    % costFxn    -> Cost - Function for Storage
    % pdDemand   -> Distribution for Order Quantity coming from Buyer
    % pdLeadTime -> Lead Time for each re-supply order Distribution
    
    t = 0;
    Ch(1) = 0;
    Co(1) = 0;
    R(1) = 0;
    Y = 0;
    I(1) = S;
    lm = Lm;
    r = rate;
    h = storageCosts;
    to = inf;
    tA =  -(1/lm)*log(rand());
    
    for i = 1:T
        if min(tA,to)==tA
            if tA>=T
                Ch(i) = Ch(i-1) + (T-t(i-1))*h*I(i-1);
                X = (R(i)-Co(i)-Ch(i))/T;
                break
            else
                Ch(i+1) = Ch(i) + (tA - t(i))*h*I(i);
                t(i+1) = tA;
                tA = t(i+1) - (1/lm)*log(rand());
                B(i) = random(pdDemand);
                m = min(I(i),B(i));
                R(i+1) = R(i) + r*m;
                I(i+1) = I(i) - m;
                if I(i+1)<s && Y == 0
                    Y = S - I(i+1);
                    L(i) = random(pdLeadTime);
                    to = t(i+1) + L(i);
                end
                Co(i+1) = Co(i) + 0;
            end
        else
            if min(tA,to) == to
                if to>=T
                    Ch(i) = Ch(i-1) + (T-t(i-1))*h*I(i-1);
                    X = (R(i)-Co(i)-Ch(i))/T;
                    break
                else
                    Ch(i+1) = Ch(i) + (to - t(i))*h*I(i);
                    Co(i+1) = Co(i) + costFxn(Y);
                    I(i+1) = I(i) + Y;
                    t(i+1) = to;
                    to = inf;
                    Y = 0;
                    R(i+1) = R(i) + 0;
                end
            end
        end
    end
    X = (R-Co-Ch);
    output = [X; Co; Ch; R; I]';
    function NetCost = costFxn(Y)
        c = 1.5;
        NetCost = Y*c;
    end
end
%% Function for Creating Plots
function [fig1,fig2,fig3,fig4,fig5]  = ModelAnalysisPlotter(output)
    T = length(output);
    X = output(:,1);
    Co = output(:,2);
    Ch = output(:,3);
    R = output(:,4);
    I = output(:,5);
    fig1 = figure

    stairs(R)
    title('Total Revenue Generated')
    grid on 
    xlabel('Time')
    xlim([0,T])
    ylabel('R_t')
    
    fig2 = figure
    stairs(Co)
    title('Cost of Resupplying Inventory')
    grid on
    xlabel('Time')
    xlim([0,T])
    ylabel('C_{ot}')
    
    fig3 = figure
    stairs(X)
    title('Net Profit')
    grid on
    xlabel('Time')
    xlim([0,T])
    ylabel('X_{t}')
    
    fig4 = figure
    stairs(Ch)
    title('Cost of Storage')
    grid on
    xlabel('Time')
    xlim([0,T])
    ylabel('C_{ht}')
    
    fig5 = figure
    stairs(I)
    title('Inventory Level')
    grid on
    xlabel('Time')
    xlim([0,T])
    ylabel('I_{t}')
end