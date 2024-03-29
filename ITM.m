%% Inverse Transform Method

% Courtesy of Prof Karl Sigman. His work is being used in these
% simulations. 
clc
close all
clear all
% The inverse transform method states that if the cumulative distribution
% of a random variable is known, then the inverse of said cumulative
% distribution can, in tandem with the uniform random variable, be used to
% generate the distribution of X, whose CDF was described. 

% For the proof, please refer to the notes published on the website - 
% http://www.columbia.edu/~ks20/4404-Sigman/4404-Notes-ITM.pdf

% This code creates structures to save generated distributions and their
% data

%% Simulating the Exponential Distribtuion

% Algorithm:
%   1 - Generate U~unif(0,1)
%   2 - Set X = (-1/Lm)*ln(U)

% Let N be the length of the random vector which is exponentially
% distributed. 
expo.N = 1000;
% Let the rate of exponential distribution be Lm. 
expo.Lm = 1.5;
for i = 1:expo.N
   U = rand();
   expo.X(i) = (-1/expo.Lm)*log(U);
end
clear i N Lm U X

figure 
histogram(expo.X,'Normalization','probability')
grid on
title('ITM Generated Exponential Distribution')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X = x)$','Interpreter','latex')
legend('\lambda = 1.5')

%% Simulating the Bernoulli(p) and Binomial(n,p) Distributions

% For the Bernoulli Distribution with success probability parameter 'p',
% we assume X follows that distribution. Hence, P(X=0) = 1-p and P(X=1)=p.

% Algorithm:
% 1 - Generate U~unif(0,1)
% 2 - Set X = 0 if U <= 1-p, X = 1 otherwise

% Let N be the length of the random vector
bern.N = 1000;
bern.p = 0.5;

for i = 1:bern.N
   U = rand();
   if U<=bern.p
       bern.X(i) = 0;
   else
       bern.X(i) = 1;
   end
end
clear i U

figure 
histogram(bern.X,'Normalization','probability')
grid on
title('ITM Generated Bernoulli Distribution')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X = x)$','Interpreter','latex')
legend('p = 0.5')
% For the Binomial Distribution, one can easily observe that a binomial
% distribution is the sum of n i.i.d. Bernoulli(p) RVs. 
bin.N = bern.N;
bin.n = 500;
U = rand(bin.N,bin.n);
bin.p = bern.p;
for i = 1:bin.n
    for j = 1:bin.N
        if U(j,i) <= 1-bin.p
            bin.Y(j,i) = 0;
        else
            bin.Y(j,i) = 1;
        end
    end
        bin.X(i) = sum(bin.Y(:,i));
end

figure 
histogram(bin.X,'Normalization','probability')
grid on
title('ITM Generated Binomial Distribution')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X = x)$','Interpreter','latex')
legend('N = 500, p = 0.5')
clear i j U 
%% Simulating the Poisson Distribution

% Algorithm is as follows:

pois.Lm = 3.5;
for i = 1:1000
    pois.n = 1;
    pois.a = 1;
    while pois.a>=exp(-pois.Lm)
        U = rand();
        pois.a = pois.a.*U;
        pois.n = pois.n + 1;
    end
    pois.X(i) = pois.n - 1;
end
clear U i 
figure
histogram(pois.X,'Normalization','probability')
grid on
title('ITM Generated Poisson Distribution')
xlabel('$x$','Interpreter','latex')
ylabel('$P(X = x)$','Interpreter','latex')
legend('\lambda = 3.5')