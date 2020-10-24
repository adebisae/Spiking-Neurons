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
% Declaring all the constants required for the model for both A and B
a=0.02; b=0.25; c=-65;  d=6;
V=-64; u=b*V;
V_B=-64; u_B=b*V_B;
VV=[];  uu=[];
VV_B=[];  uu_B=[];
tau = 0.25;
tspan = 0:tau:steps;  %tau is the discretization time-step
%tspan is the simulation interval

% create empty spike list for both A and B
spike_ts = [];
spike_ts_B = [];

% initialize empty arrays for the output of neuron B
init = 1;
VV_val = 0;
All_VVs = ones(4001,41);
All_spikes = ones(4001,41)*3;
%--------------------------------------------------------------------------

% open the loop for B
for all_trials = 0:0.5:20
    % open the timespan for each value of I in B
    for t=tspan
    %----------------------------------------------------------------------
    % model A with constant value of I = 5
        I = 5;
        V = V + tau*(0.04*V^2+5*V+140-u+I);
        u = u + tau*a*(b*V-u);
        if V > 30                                       %if this is a spike
            VV(end+1)=30;     %VV is the time-series of membrane potentials
            V = c;
            u = u + d;
            spike_value = 1;
            spike_ts = [spike_ts ; spike_value];   %records a spike
        else
            VV(end+1)= V;
            spike_value = 0;
            spike_ts = [spike_ts ; spike_value];   %records no spike
        end
        uu(end+1)=u;
    %----------------------------------------------------------------------
    % obtain result from I and use it to model B
        I_B = all_trials + 10*spike_value;
        V_B = V_B + tau*(0.04*V_B^2+5*V_B+140-u_B+I_B);
        u_B = u_B + tau*a*(b*V_B-u_B);
        if V_B > 30                                     %if this is a spike
            VV_B(end+1)=30;  %VV is the time-series of membrane potentials
            V_B = c;
            u_B = u_B + d;
            spike_ts_B = [spike_ts_B ; 1];   %records a spike
        else
            VV_B(end+1)= V_B;
            spike_ts_B = [spike_ts_B ; 0];   %records no spike
        end
        uu_B(end+1)=u_B;
    end
    %----------------------------------------------------------------------
    % store all the membrane potentials and the spikes
    All_VVs(:,init) = VV_B(VV_val+1:length(VV_B));
    All_spikes(:,init) = spike_ts_B(VV_val+1:length(VV_B));
    %----------------------------------------------------------------------
    % draw the required plots and find means
    if all_trials == 1
        sp_800 = All_spikes(:,init);      %# last 800 elements
        mean_1 = [all_trials, sum(sp_800(end-3200+1:end)/800)];
        subplot(5,1,1)
        plot(tspan, All_VVs(:,init),'r','linewidth', 1.5);
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
        plot(tspan, All_VVs(:,init), 'g','linewidth', 1.5);
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
% draw the required R vs I plots superimposed for A and B
MSR_B  = [mean_1;mean_5;mean_10;mean_15;mean_20];
MSR_A  = [1,0.0088;5,0.0213;10,0.0338;15,0.0463;20,0.06];
figure
plot(MSR_A(:,1),MSR_A(:,2),'k', 'linewidth', 1.5)
hold on
plot(MSR_B(:,1),MSR_B(:,2),'r', 'linewidth', 1.5)
legend('Sinlge Neuron', 'Two Neurons', 'location', 'best')
xlabel('I values')
ylabel('R values')
title('Plot of R vs I')
