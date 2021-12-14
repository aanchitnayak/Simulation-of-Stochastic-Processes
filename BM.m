clc
clear all
close all
t = 5;
n = 1000;
% rng('default')
% randn('state',100);
dt = t/n;
w = zeros(1,n);
dw = zeros(1,n);
w(1) = 0;
for j = 2:n
   w(j) = w(j-1)+sqrt(dt)*randn; 
end

x(1) = 0;
mu = 2;
sig = 2.5;
for j = 2:n
    x(j) = x(j-1) + sig*sqrt(dt)*randn + mu*dt;
end

s(1) = 1;
% for j = 2:n
%    s(j) = s(j-1)*exp(x(j)-x(j-1)); 
% end
for j = 1:n
   y(j) = exp(sig*sqrt(dt)*normrnd(0,1) + mu*dt); 
end
for j = 2:n
   s(j) = s(j-1)*y(j); 
end

T = [dt:dt:t];
test = exp(mu.*T);
figure
plot([dt:dt:t],w,'r-')
grid on
title('Simple Brownian Motion')
xlabel('t')
ylabel('B_t')
figure
plot([dt:dt:t],x,'b-')
grid on
title('BM with Drift')
xlabel('t')
ylabel('X_t')
figure
plot([dt:dt:t],s)
grid on
title('GBM')
xlabel('t')
ylabel('S_t')
