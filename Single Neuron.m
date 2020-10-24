%   This MATLAB file generates figure 1 in the paper by
%               Izhikevich E.M. (2004)
%   Which Model to Use For Cortical Spiking Neurons?
%   use MATLAB R13 or later. November 2003. San Diego, CA

%   Modified by Ali Minai
%   Modified further by Adekunle Adebisi 09/20/2020
%--------------------------------------------------------------------------
clc
clear
close all
%--------------------------------------------------------------------------
steps = 1000;                  %This simulation runs for 1000 steps

% Declaring all the constants required for the model
a=0.02; b=0.25; c=-65;  d=6;
V=-64; u=b*V;
VV=[];  uu=[];
tau = 0.25;
tspan = 0:tau:steps;  %tau is the discretization time-step
%tspan is the simulation interval

% creating an empty list for the spikes
spike_ts = [];

% initializing variable to store the required output membrane potential &
% spikes
init = 1;
VV_val = 0;
All_VVs = ones(4001,41);
All_spikes = ones(4001,41)*3;  % array of 3's to be replaced by 0 or 1 for spikes
%--------------------------------------------------------------------------

% starting the loop for each value of I
for all_trials = 0:0.5:20
    % starting the loop for the timespan 1000 ofr each I in the above loop
    for t=tspan
    %----------------------------------------------------------------------
        I = all_trials;
        V = V + tau*(0.04*V^2+5*V+140-u+I);
        u = u + tau*a*(b*V-u);
        if V > 30                 %if this is a spike
            VV(end+1)=30;         %VV is the time-series of membrane potentials
            V = c;
            u = u + d;
            spike_ts = [spike_ts ; 1];   %records a spike
        else
            VV(end+1)= V;
            spike_ts = [spike_ts ; 0];   %records no spike
        end
        uu(end+1)=u;
    end
    % store all the results into the intialize empty arrays for each loop
    All_VVs(:,init) = VV(VV_val+1:length(VV));
    All_spikes(:,init) = spike_ts(VV_val+1:length(VV));
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % drawing the plots only for I = 1,5,10,15, and 20
    if all_trials == 1
        sp_800 = All_spikes(:,init);      %# last 800 elements
        mean_1 = [all_trials, sum(sp_800(end-3200+1:end)/800)];
        subplot(5,1,1)
        plot(tspan, All_VVs(:,init),'r', 'linewidth', 1.5);
        legend('I = 1', 'location', 'best')
        axis([0 max(tspan) -90 40])
        xlabel('time step');
        ylabel('V_m')
        xticks([0 max(tspan)]);
        xticklabels([0 steps]);
        title('Regular Spiking');
    elseif (all_trials == 5)
        sp_800 = All_spikes(:,init);      %# last 800 elements
        mean_5 = [all_trials, sum(sp_800(end-3200+1:end)/800)];
        subplot(5,1,2)
        plot(tspan, All_VVs(:,init),'b', 'linewidth', 1.5);
        legend('I = 5', 'location', 'best')
        axis([0 max(tspan) -90 40])
        xlabel('time step');
        ylabel('V_m')
        xticks([0 max(tspan)]);
        xticklabels([0 steps]);
    elseif (all_trials == 10)
        sp_800 = All_spikes(:,init);      %# last 800 elements
        mean_10 = [all_trials, sum(sp_800(end-3200+1:end)/800)];
        subplot(5,1,3)
        plot(tspan, All_VVs(:,init),'g', 'linewidth', 1.5);
        legend('I = 10', 'location', 'best')
        axis([0 max(tspan) -90 40])
        xlabel('time step');
        ylabel('V_m')
        xticks([0 max(tspan)]);
        xticklabels([0 steps]);
    elseif (all_trials == 15)
        sp_800 = All_spikes(:,init);      %# last 800 elements
        mean_15 = [all_trials, sum(sp_800(end-3200+1:end)/800)];
        subplot(5,1,4)
        plot(tspan, All_VVs(:,init),'m', 'linewidth', 1.5);
        legend('I = 15', 'location', 'best')
        axis([0 max(tspan) -90 40])
        xlabel('time step');
        ylabel('V_m')
        xticks([0 max(tspan)]);
        xticklabels([0 steps]);
    elseif (all_trials == 20)
        sp_800 = All_spikes(:,init);      %# last 800 elements
        mean_20 = [all_trials, sum(sp_800(end-3200+1:end)/800)];
        subplot(5,1,5)
        plot(tspan, All_VVs(:,init),'c', 'linewidth', 1.5);
        legend('I = 20', 'location', 'best')
        axis([0 max(tspan) -90 40])
        xlabel('time step');
        ylabel('V_m')
        xticks([0 max(tspan)]);
        xticklabels([0 steps]);
    end
    %----------------------------------------------------------------------
    init = init + 1;
    VV_val = VV_val + 4001;
end
%--------------------------------------------------------------------------
% drawing the plots for the R vs I
MSR  = [mean_1;mean_5;mean_10;mean_15;mean_20];
figure
plot(MSR(:,1),MSR(:,2),'k', 'linewidth', 1.5)
xlabel('I values')
ylabel('R values')
title('Plot of R vs I')




