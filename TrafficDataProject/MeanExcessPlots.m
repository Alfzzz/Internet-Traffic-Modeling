%Mean Excess plots
clc
clear all
close all
u=0:0.1:20;
%% Normal
fun = @(x)qfunc(x);
ME_normal=[];
for i=1:length(u)
    ME_normal_num = integral(fun,u(i),Inf);
    ME_normal_den=qfunc(u(i));
    ME_normal(i)=ME_normal_num/ME_normal_den;
end
plot(u,ME_normal)
title("Normal Mean excess plot")
xlabel("u")
ylabel("M(u)")

%% Log-normal
fun = @(x)qfunc(log(x));
ME_lognormal=(exp(1/2)*(qfunc(log(u)-1)))./(qfunc(log(u)))-u;
figure
plot(u,ME_lognormal)
title("Lognormal Mean excess plot")
xlabel("u")
ylabel("M(u)")

%% Uniform
a=3;
b=5;
u_uniform1=0:0.1:a;
u_uniform2=a:0.1:b;
ME_uniform1=0.5*(a+b)-u_uniform1; %for u<a
ME_uniform2=a-u_uniform2-((a-b).^2)./(2*(u_uniform2-b)); %for a<u<b
figure
subplot(2,1,1)
plot(u_uniform1,ME_uniform1)
title("Uniform Mean excess plot(u<a)")
xlabel("u")
ylabel("M(u)")
subplot(2,1,2)
plot(u_uniform2,ME_uniform2)
title("Uniform Mean excess plot(a<u<b)")
xlabel("u")
ylabel("M(u)")
%% Exponential
lambda=1;
figure
ME_exponential=ones(1,length(u)).*1/lambda;
plot(u,ME_exponential)
title("Exponential Mean excess plot")
xlabel("u")
ylabel("M(u)")
%% Pareto
xm=4;%mode
alpha=2;
fun = @(x)(xm./x).^alpha;
ME_pareto=[];
for i=1:length(u)
    ME_pareto_num = integral(fun,u(i),Inf);
    ME_pareto_den=(xm./u(i)).^alpha;
    ME_pareto(i)=ME_pareto_num./ME_pareto_den;
end
figure
plot(u,ME_pareto)
title("Pareto Mean excess plot")
xlabel("u")
ylabel("M(u)")

%% Weibul
tau1=2;
tau2=0.2;
beta=2;
ME_weibul1=(u.^(1-tau1))./(beta*tau1);
ME_weibul2=(u.^(1-tau2))./(beta*tau2);
figure
subplot(2,1,1)
plot(u,ME_weibul1)
title("Weibul Mean excess plot(tau<1)")
xlabel("u")
ylabel("M(u)")
subplot(2,1,2)
plot(u,ME_weibul2)
title("Weibul Mean excess plot(tau>1)")
xlabel("u")
ylabel("M(u)")

%% All Mean excess function
close all
plot(u,ME_normal,u,ME_lognormal,u,ME_exponential,u,ME_pareto,u,ME_weibul1,u,ME_weibul2,u_uniform1,ME_uniform1,'LineWidth',3)
legend("Normal","Lognormal","Exponential","Pareto","Weibul(tau>1)","Weibul(tau<1)","Uniform(u<a)")