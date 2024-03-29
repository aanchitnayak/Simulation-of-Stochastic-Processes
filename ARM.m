%% Using ARM - Acceptance Rejection Method for simulating RVs

% Courtesy of Prof Karl Sigman. The notes can be found on:
% http://www.columbia.edu/~ks20/4404-Sigman/4404-Notes-ARM.pdf

% The ITM method for simulating a random variable is difficult to use in
% all cases because it requires the knowledge of the explicit form of the
% CDF of a RV. ARM is a better algorithm for simulating RVs, all in all.

% Begin by assuming that for a continuous RV X, if CDF is F(x) then PDF 
% is f(x). The basic idea is to find an alternative probability
% distribution G with density function g(x) from which we already have an
% efficient algorithm for generating from (like ITM) but also such that
% g(x) is close to f(x).

% In particular, we assume that the ratio f(x)/g(x) is bounded by some c>0;
% sup_x {f(x)/g(x)} <= c (\leq) and in practice we would want c to tend to
% unity.

% General algorithm for generating X distributed as F:
% 1 - Generate a rv Y ~ G
% 2 - Generate U independent from Y
% 3 - If
%           U <= f(Y)/(c*g(Y); then
%     set X = Y ("accept); otherwise, go to 1

% Note that f(Y) and g(Y) are RVs and the ratio is indepenedent of U. The
% ratio is bounded by 0 and 1. The number of times, N, that steps 1 and 2
% need to be called is itself a RV, which has the geometric distribution
% with success probability = p s.t. p = P(U<= f(Y)/(c*g(Y)) 
% and
% P(N = n) = p*(1-p)^{n-1} for n>= 1. On average, E[N] = 1/p
% In the end we obtain X as having the conditional distribution of a Y
% given that the event {U <= f(Y)/cg(Y)} occurs.

% A simple pen-paper calculation tells us  p = 1/c. Thus, E[N] = c, the
% bounding constant and we can not indeed see that it is desirable to
% choose our alternative density 'g' so as to minimize this constant as:

% c = sup_x {f(x)/g(x)}

% Optimally, f(x) = g(x) would be a trivial solution. But the point is to
% find a g(x) which is closer to f(x) and not f itself! In summary, the
% expected number of iterations before X is obtained is the bounding
% constant 'c', it self.

clc
clear all
close all

%% Simulating the Normal Distribution

% We desire to generate X~N(mu,sigma). We know, X = mu*Z + sigma 
% where Z~N(0,1). It suffices to find an algorithm for generating Z ~
% N(0,1). 

% for the function close to normal,we choose the exponential

% Algorithm is as follows:
% 1 - Generate two independent exponentials Y1 and Y2 with rate 1
% 2 - if Y2 <= (Y1 - 1)^2/2, set modZ = Y1 else, go back to 1
% 3 - Generate U. Set Z = modZ if U<= 0.5. Set Z = -modZ if U>0.5
iter = 0;

DistLength = 10000;
for i = 1:DistLength
    Y1(i) = 0;
    Y2(i) = 0;
    while Y2(i) < (Y1(i) - 1)^2/2
        Y1(i) = -log(rand());
        Y2(i) = -log(rand());
        iter = iter + 1;
    end
    epoch(i) = iter;
    modZ(i) = Y1(i);
    if rand() <= 0.5
        Z(i) = modZ(i);
    else
        Z(i) = - modZ(i);
    end
end

tempdist = makedist('Normal',0,1);

figure
subplot(211)
histogram(Z,'Normalization','pdf','DisplayName','Generated Data')
hold on
plot(-10:0.01:10,pdf(tempdist,-10:0.01:10),'DisplayName','In-Built MATLAB pdf','LineWidth',2)
grid on 
title('Standard Normal Distribution','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X = x)$','Interpreter','latex')
hl = legend('show');
set(hl, 'Interpreter','latex')

subplot(212)
histogram(epoch,'Normalization','pdf','DisplayName','Geometric Distribution')
grid on
title('Distribution of Epochs Required to Generate a Gaussian','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X=x)$','Interpreter','latex')
hl = legend('show');
set(hl, 'Interpreter','latex')

clear h1 tempdist modZ Y1 Y2 i iter 
