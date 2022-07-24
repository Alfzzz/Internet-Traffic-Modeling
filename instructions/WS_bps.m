clear all
close all
clc

% Read data file from a .csv file with two columns, first column being time
% in seconds, second column being length of packets in bytes.

M = readmatrix('TrafficData_bps_WiFi_2.csv');

% Convert second column of matrix M from bytes to bits
Mbps=M;
Mbps(:,2)=M(:,2)*8;

% Now get the bps, i.e., group every packet in each second and add up the
% amount of bits

[Mx,My]=size(Mbps);

% Round up the last time recorded, i.e., use ceiling function
last_time=ceil(Mbps(Mx,1));

% Collect all entries of each second, at time 0 we put 0 for initial
% conditions, at time 1 second we put the accumulation of bits from the
% first packet up to that that is recorded before the time 1 second, and so
% on for the following seconds

traffic_bps=zeros(last_time+1,My);
for i=1:last_time
    % collect data from time i-1 up to time <= i while time i is less than
    % the last_time recorded
    time_data_i=Mbps(Mbps(:,1)< i & Mbps(:,1) >= (i-1),:);
    traffic_bps(i+1,:)=[i sum(time_data_i(:,2))];
end

figure (1)
plot(traffic_bps(:,1),traffic_bps(:,2),'r.')
grid
xlabel('time (sec)')
ylabel('bps')

%% Statistical parameters
traffic_mean=mean(traffic_bps(:,2));
traffic_var=var(traffic_bps(:,2));
traffic_std=std(traffic_bps(:,2));
traffic_median=median(traffic_bps(:,2));

%% Distribution CDF, CCDF and density PDF
% CDF
[xx,yy,stats] = mycdfplot2(traffic_bps(:,2));
figure (2)
plot(xx,yy,'k')
grid
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('$P(X\leq x)$','Interpreter','latex')

% CCDF or survival function
figure (3)
plot(xx,(1-yy),'b')
grid
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('$1-P(X\leq x)$','Interpreter','latex')

% Estimating density with Kernel normal, box, triangle or epanechnikov
pd_Kbps = fitdist(traffic_bps(:,2),'Kernel','Kernel','epanechnikov');
x_values = xx;
y_Kbps = pdf(pd_Kbps,x_values);
% figure (4)
% plot(x_values,y_Kbps)
% grid
% xlabel('$x$ (bps)','Interpreter','latex')
% ylabel('Kernel of $f_X(x)$','Interpreter','latex')


%% COMPARISONS

% Is it Gaussian or Normal???
figure (5)
qqplot(traffic_bps(:,2))
grid

%% Is it exponential???
% Generate the cdf or pdf of an exponential random variable with the same
% mean value as traffic_bps vector
exp_pdf=pdf('Exponential',xx,traffic_mean);

figure (6)
plot(xx,exp_pdf,'r.')
hold on
plot(xx, y_Kbps,'b')
hold off
grid
xlabel('$x$ (bps)','Interpreter','latex')
%ylabel('$\frac{1}{\mu}e^{-\frac{x}{\mu}}$','Interpreter','latex')
ylabel('$f_X(x)$','Interpreter','latex')

%%
figure (7)
semilogy(xx,exp_pdf,'r.')
hold on
semilogy(xx, y_Kbps,'b')
hold off
grid
xlabel('$x$ (bps)','Interpreter','latex')
%ylabel('$\frac{1}{\mu}e^{-\frac{x}{\mu}}$','Interpreter','latex')
ylabel('$f_X(x)$','Interpreter','latex')









