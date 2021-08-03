%% Simulating FIFO GI|GI|1 Queue

% This code is the courtesy of Prof Karl Sigman's notes which can be found
% on: 
% http://www.columbia.edu/~ks20/4404-Sigman/4404-Notes-SSQ.pdf

% In this program, we simulate the working of a FIFO, GI|GI|1 queue and
% also intend to obtain the average delay time (as time -> inf) and the
% average dealy for nth incoming customer to a queue. 

% Key to this simulation is proposition mentioned in Prof Sigman's notes.
% That proposition will allow for the generation of delays for incoming
% customer easily. 

% Let:
%   1 - Sn -> Service Time of nth customer, Uniform RV~U(a,b)
%   2 - Tn -> Inter-arrival times, Exponential RV ~ (Lm)
%   3 - Dn -> Delay time for nth customer
%   4 - Cn -> nth arriving customer

clc
close all
clear all

SimLength = 10;

%% Generating a random sequence for service time of customers:

service.a = 1;
service.b = 3;
service.G = makedist('Uniform',service.a,service.b);
service.S = random(service.G,1,SimLength);

%% Generating a random sequence of inter-arrival times of customers:

iaTime.Lm = 2;
iaTime.A = makedist('Exponential',iaTime.Lm);
iaTime.T = random(iaTime.A,1,SimLength);

%% Generating the delay times for the queue

delay.D(1) = 0;
for i = 2:SimLength
   delay.D(i) = max((delay.D(i-1)+ service.S(i-1) - iaTime.T(i-1)),0);
end

% %% Time Table
% TimeTable = [1:SimLength;iaTime.T;service.S;delay.D]';
%% Long-Run Simulation -> Estimating Average Delay

% We aim to find the following - 
%   1 - The long run average delay of the queue
%   2 - The number of customers having delay more than a given value, say x

longrun.SimLength = 10000; % Take a large number for a long-run simulation

% for 1:
longrun.S = random(service.G,1,longrun.SimLength);
longrun.T = random(iaTime.A,1,longrun.SimLength);
longrun.D(1) = 0;
for i = 2:longrun.SimLength
   longrun.D(i) = max((longrun.D(i-1)+ longrun.S(i-1) - longrun.T(i-1)),0);
end
clear i 
longrun.d = 1/(longrun.SimLength) * sum(longrun.D); % indicated the long run average time 

% for 2:
longrun.x = 1.5;
longrun.count = 0;
for i = 1:longrun.SimLength
   if longrun.D(i)<longrun.x
       longrun.count = longrun.count + 1;
   end
end
clear count i
longrun.I = round(1/longrun.SimLength * longrun.count*100);

%% Monte Carlo Simulation for average delay at for a given customer

% We wish to observe the average delay time suffered by the 10th arrival
mcs.Cj = 10; 
mcs.N = 10000; %N times repeating the calculation of delay for 10th customer
mcs.S = zeros(mcs.Cj,mcs.N);
mcs.T = zeros(mcs.Cj,mcs.N);
mcs.D = zeros(mcs.Cj,mcs.N);
for j = 1:mcs.N
    mcs.S(:,j) = random(service.G,1,mcs.Cj);
    mcs.T(:,j) = random(iaTime.A,1,mcs.Cj);
    mcs.D(1,1) = 0;
    for i = 2:mcs.Cj
       mcs.D(i,j) = max((mcs.D(i-1,j)+ mcs.S(i-1,j) - mcs.T(i-1,j)),0); 
    end
end
clear i j 

mcs.Dj = 1/(mcs.N)*sum(mcs.D(mcs.Cj,:));