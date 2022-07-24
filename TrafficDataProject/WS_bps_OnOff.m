clear all
close all
clc

% Read data file from a .csv file with two columns, first column being time
% in seconds, second column being length of packets in bytes.

M = readmatrix('mytraffic_depurado.xlsx');
%M = readmatrix('TrafficData_bps_WiFi_2.csv');
%M = readmatrix('TrafficData_WiFi_Mar_28_2020');

% Convert second column of matrix M from bytes to bits
Mbps=M;
Mbps(:,2)=M(:,2)*8;

% Now get the bps, i.e., group every packet in each second and add up the
% amount of bits

[Mx,My]=size(Mbps);

% Round up the last time recorded, i.e., use ceiling function
last_time=ceil(Mbps(Mx,1));
t_delta=1; % time-slot in seconds it could be 1, 1/2, 1/4, 1/8, etc.

%% TRAFFIC MEASUREMENT IN bps

% Collect all entries of each second, at time 0 we put 0 for initial
% conditions, at time 1 second we put the accumulation of bits from the
% first packet up to what is recorded before the time 1 second, and so
% on for the following seconds

%traffic_bps=zeros(last_time+1,My);
traffic_bps=[0 0];
for i=t_delta:t_delta:last_time
    % collect data from time i-1 up to time <= i while time i is less than
    % the last_time recorded
    time_data_i=Mbps(Mbps(:,1)< i & Mbps(:,1) >= (i-t_delta),:);
    traffic_bps=[traffic_bps ; i sum(time_data_i(:,2))];
end
traffic_bps(:,2)=traffic_bps(:,2)*(1/t_delta);
figure (1)
    stem(traffic_bps(:,1),traffic_bps(:,2),'r')
    title('Internet Traffic Transmission in 4 min');
    xlabel('$Time$ (sec)','Interpreter','latex')
    ylabel('$Bits per second$ (bps)','Interpreter','latex')
    xlim([0 255])

[rtx,rty]=size(traffic_bps);
% Running mean
run_traffic=cumsum(traffic_bps(2:rtx,2));
run_time=(traffic_bps(2:rtx,1));
run_average=run_traffic./run_time;

figure(2)
plot(traffic_bps(2:rtx),run_average,'b*')
grid
xlabel('time (sec)')
ylabel('Average bps')
title('Running (Time) Average')

% Statistical parameters
traffic_mean=mean(traffic_bps(:,2)); % Arithmetic mean
traffic_var=var(traffic_bps(:,2));
traffic_std=std(traffic_bps(:,2));
traffic_median=median(traffic_bps(:,2));

%% Distribution CDF, CCDF and density PDF of TRAFFIC
% CDF
[xxt,yyt,statst] = mycdfplot2(traffic_bps(:,2));
figure (3)
plot(xxt,yyt,'k')
grid
title('CDF');
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('$P(X\leq x)$','Interpreter','latex')

%Exponential distribution CCDF
exp_ccdf_bps=exp(-(1/statst.mean)*xxt);
% CCDF or survival function
figure (4)
loglog(xxt,exp_ccdf_bps,'r',xxt,(1-yyt),'b')
grid
title('Survival Function')
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('$1-P(X\leq x)$','Interpreter','latex')
legend('Exponential','bps','Location','southwest')

%% Hill Estimator of TRAFFIC

% Order statistics
% Sort vector of data (traffic_bps)
Xn=sort(traffic_bps(:,2),'descend');
% Get rid of the zeros
Zn=Xn(Xn>0);

[xw,xz]=size(Zn(:));
for k=1:xw-1
    xni=Zn(1:k);
    h1=log(xni./Zn(k+1));
    Hkn(k)=(1/k)*sum(h1);
end
Hkn=Hkn';
figure(11)
plot(1./Hkn,'r*')
grid
xlabel('Number of order statistic')
ylabel('Hill estimate of alpha')
title('Hill Estimator of Traffic')
%% Generate ON/OFF process

[tx,ty]=size(traffic_bps);
% Identify one-second time-slots where there is no presence of traffic
traffic_on_off=[traffic_bps(:,1) traffic_bps(:,2)>0];
service_fault=randi([1 size(traffic_on_off,1)],1,size(traffic_on_off,1));

for i=1:length(service_fault)
    traffic_on_off(service_fault(1,i),2)=0;
end

% Plot On Off process
figure (5)
stem(traffic_on_off(:,1), traffic_on_off(:,2),'r','LineWidth', 1.5)
grid
xlabel('time (sec)')
ylabel('Indicator of traffic')
axis([0 last_time 0 1.5])
title('ON and OFF Process')

figure (6)
plot(traffic_bps(:,1),traffic_bps(:,2),'r.')
grid
xlabel('time (sec)')
ylabel('bps')
title('ON and OFF Process')
%% Generate statistics of ON periods

% Counting ON periods
ON_duration=[];
on_count=0;
k=1;
for i=1:tx
    if (traffic_on_off(i,2)==1)
        % Count the number of ones
        on_count=on_count+1;
    else
        if (on_count ~= 0)
            ON_duration(k)=on_count;
            k=k+1;
        end
        on_count=0;
    end
end
ON_duration=ON_duration';
% CDF
[xx,yy,stats] = mycdfplot2(ON_duration);
figure (7)
plot(xx,yy,'k')
grid
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('$P(X\leq x)$','Interpreter','latex')

%Exponential distribution CCDF
exp_ccdf_ON=exp(-(1/stats.mean)*xx);
% CCDF or survival function
figure (7)
loglog(xx,exp_ccdf_ON,'r',xx,(1-yy),'b')
grid
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('$1-P(X\leq x)$','Interpreter','latex')
legend('Exponential','ON','Location','southwest')

% Index of dispersion (peakedness)
pk=stats.std^2/stats.mean;
HEST = wfbmesti(ON_duration);

%% Estimating density with Kernel normal, box, triangle or epanechnikov
pd_Kbps = fitdist(ON_duration,'Kernel','Kernel','triangle');
x_values = xx;
y_Kbps = pdf(pd_Kbps,x_values);
figure (8)
plot(x_values,y_Kbps)
grid
xlabel('$x$ (bps)','Interpreter','latex')
ylabel('Kernel of $f_X(x)$','Interpreter','latex')
title('PDF stimation')

%% INTERPACKET TIMES (IPT)
timevector=Mbps(:,1);
[tvx,tvy]=size(timevector(:));

tv2=timevector(2:tvx);
IPT=tv2-timevector(1:tvx-1);
% CDF
[xxipt,yyipt,statsipt] = mycdfplot2(IPT);

% CCDF or survival function
figure (9)
loglog(xxipt,(1-yyipt),'b')
grid
xlabel('IPT (sec)')
ylabel('$1-P(IPT\leq x)$','Interpreter','latex')
title('CDF')

% Is it exponential???
% Generate the cdf or pdf of an exponential random variable with the same
% mean value as traffic_bps vector
%exp_pdf=pdf('Exponential',xxipt,statsipt.mean);
exp_ccdf=exp(-(1/statsipt.mean)*xxipt);
figure (10)
loglog(xxipt,exp_ccdf,'r',xxipt,(1-yyipt),'b')
plot(log(xxipt),log(exp_ccdf),'r.',log(xxipt),log(1-yyipt),'b','LineWidth',2)
hold on
loglog(xxipt,(1-yyipt),'b')
hold off
grid
xlabel('IPT (sec)')
ylabel('$1-P(IPT\leq x)$','Interpreter','latex')
legend('Exponential','IPT','Location','southwest')