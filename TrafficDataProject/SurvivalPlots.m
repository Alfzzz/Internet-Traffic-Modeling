%Survival function Plots
clc
clear all
close all
%% Normal
x=-4:0.1:4;
S_normal=qfunc(x);
plot(x,S_normal)
title("Normal Survival function")
xlabel("x")
ylabel("S(x)")

%% Log-normal
x=00:0.1:20;
S_lognormal=qfunc(log(x));
plot(x,S_lognormal)
title("Lognormal Survival function ")
xlabel("x")
ylabel("S(x)")

%% Uniform
a=3;
b=5;
x=3:0.1:5;
S_uniform=(b-x)./(b-a);
figure
plot(x,S_uniform)
title("Uniform Survival function(a<x<b)")
xlabel("x")
ylabel("S(x)")
%% Exponential
lambda=1;
figure
x=0:0.1:20;
S_exponential=exp(1).^(-x.*lambda);
plot(x,S_exponential)
title("Exponential Survival function")
xlabel("x")
ylabel("S(x)")
%% Pareto
x=0:0.1:20;
xm=4;%mode
alpha=2;
S_pareto=(xm./x).^alpha;
figure
plot(x,S_pareto)
title("Pareto Survival function")
xlabel("x")
ylabel("S(x)")

%% Weibul
tau1=2;
tau2=0.2;
beta=2;
S_weibul=exp(1).^(-x*beta);
figure
plot(x,S_weibul)
title("Weibul Survival function")
xlabel("x")
ylabel("S(x)")