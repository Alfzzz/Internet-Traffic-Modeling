clear all
close all
clc
%Leer datos de Wireshark
% Read data file from a .csv file with two columns, first column being time
% in seconds, second column being length of packets in bytes.
%M = readmatrix('mytraffic.csv');
M = readmatrix('wifi_Alfred.csv');
time=M(:,2);
B=M(:,6);
% Convert second column of matrix M from bytes to bits
Mbps=[time,B];
Mbps(:,2)=M(:,2)*8;
% Now get the bps, i.e., group every packet in each second and add up the
% amount of bits
[Mx,My]=size(Mbps);
% Round up the last time recorded, i.e., use ceiling function
last_time=ceil(Mbps(Mx,1));
%   BPS 
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
    title('Internet Traffic Transmission in 4 min');
    xlabel('$Time$ (sec)','Interpreter','latex')
    ylabel('$Bits per second$ (bps)','Interpreter','latex')
% Statistical parameters
traffic_mean    =        mean(traffic_bps(:,2));
traffic_min     =        min(traffic_bps(:,2));
traffic_max     =        max(traffic_bps(:,2));
traffic_var     =        var(traffic_bps(:,2));
traffic_std     =        std(traffic_bps(:,2));
traffic_median  =        median(traffic_bps(:,2));
traffic_disp    =        traffic_var./traffic_mean;

% CDF
[xx_bps,yy_bps,stats_bps] = mycdfplot2(traffic_bps(:,2));
figure (2)
    plot(xx_bps,yy_bps,'k')
    grid
    title('CDF of bps');
    xlabel('$x$ (bps)','Interpreter','latex')
    ylabel('$P(X\leq x)$','Interpreter','latex')
% CCDF or survival function
figure (3)
    plot(xx_bps,(1-yy_bps),'b')
    grid
    title('Survival Function of bps');
    xlabel('$x$ (bps)','Interpreter','latex')
    ylabel('$1-P(X\leq x)$','Interpreter','latex') 
% CCDF vs Expo CCDF BPS
figure (4)
    exp_ccdf_bps=exp(-(1/stats_bps.mean)*xx_bps);
    loglog(xx_bps,exp_ccdf_bps,'r',xx_bps,(1-yy_bps),'b')
    grid on
    title('CCDF of bps and Exponential Distribution PDF bps');
    legend('CCDF of bps','Exponential Distribution PDF bps')
% PDF
figure (5)
    histogram(traffic_bps(:,2),'Normalization', 'pdf')
    grid
    title('PDF of bps');
    xlabel('$x$ (bps)','Interpreter','latex')
    ylabel('$f_X(x)$','Interpreter','latex')
% Histogram
figure (6)
    histogram(traffic_bps(:,2))
    grid
    title('Histogram of bps');
    xlabel('$x$ (bps)','Interpreter','latex')
    ylabel('$Accumulation$','Interpreter','latex')
%Mean Excess Value
sum_integral=0;
e=[];
threshold=0:length(xx_bps)-1;
for i=2:length(xx_bps)
    for j=i:length(xx_bps)
        sum_integral=sum_integral+(1-yy_bps(j,1)).*(xx_bps(j,1)-xx_bps(j-1,1));
    end
    e(i,1)=sum_integral/(1-yy_bps(i,1));
    sum_integral=0;
end
figure (7)
plot(threshold,e)
    grid
    title('Mean Excess Value of bps');
    xlabel('$Threshold$','Interpreter','latex')
    ylabel('$e(u)$','Interpreter','latex')

%% IPT
ipt=[];
for i=2:Mx
    ipt=[ipt;Mbps(i,1)-Mbps(i-1,1)];%resta del segundo con el primero y sucesivamente
end
figure (8)
    plot(ipt)
    title('InterPacket Time');
    xlabel('$Number$ $of$ $Deliveries$ $n$','Interpreter','latex')
    xlim([0 3.065*10^4]);
    ylabel('$IPT [s]$','Interpreter','latex')
%CDF
[xx_ipt,yy_ipt,stats_ipt] = mycdfplot2(ipt);
figure (9)
    plot(xx_ipt,yy_ipt,'k')
    grid
    title('CDF of IPT');
    xlabel('$x$ (IPT)','Interpreter','latex')
    ylabel('$P(X\leq x)$','Interpreter','latex')
%CCDF
figure (10)
    plot(xx_ipt,(1-yy_ipt),'b')
    grid
    title('Survival Function/ CCDF of IPT');
    xlabel('$x$ (IPT)','Interpreter','latex')
    ylabel('$1-P(X\leq x)$','Interpreter','latex')
% CCDF vs Expo CCDF BPS
figure (11)
    exp_ccdf_ipt=exp(-(1/stats_ipt.mean)*xx_ipt);
    loglog(xx_ipt,exp_ccdf_ipt,'r',xx_ipt,(1-yy_ipt),'b')
    grid on
    title('CCDF of IPT and Exponential Distribution PDF IPT');
    legend('CCDF of IPT','Exponential Distribution PDF IPT')
%PDF
figure (12)
    grid
    histogram(ipt)
    title('PDF of IPT');
    xlabel('$x$ (IPT)','Interpreter','latex')
    ylabel('$f_X(x)$','Interpreter','latex')
% Histogram
figure (13)
    grid
    histogram(ipt)
    title('Histogram of IPT');
    xlabel('$x$ (IPT)','Interpreter','latex')
    ylabel('$Accumulation$','Interpreter','latex')
%Mean Excess Value
sum_integral=0;
e=[];
threshold=0:length(xx_ipt)-1;
for i=2:length(threshold)
    for j=i:length(xx_ipt)
        sum_integral=sum_integral+(1-yy_ipt(j,1)).*(yy_ipt(j,1)-yy_ipt(j-1,1));
    end
    e(i,1)=sum_integral/(1-yy_ipt(i,1));
    sum_integral=0;
end
figure (14)
plot(threshold,e)
    grid
    title('Mean Excess Value of IPT');
    xlabel('$Threshold$','Interpreter','latex')
    ylabel('$e(u)$','Interpreter','latex')
