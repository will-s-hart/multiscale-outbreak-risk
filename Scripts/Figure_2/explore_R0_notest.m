% Calculate the outbreak risk without regular antigen testing for a range
% of R0 values

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_in.mat','params_vec','R0','tau_vec')
load('../../Results/Figure_2/WH_det_inf_dynamics.mat','beta_tot_notest')

tau_inc = params_vec(5);
eta = params_vec(6);
prop_pop = params_vec(7); %trivially 1 in 1-group case

R0_default = R0;
beta_tot_default = beta_tot_notest;

% Loop over R0 values

R0_vec = (0.01:0.01:4)';
p_outbreak_vec = zeros(size(R0_vec));

for i = 1:length(R0_vec)
        
    beta_tot = (R0_vec(i)/R0_default)*beta_tot_default;
    
    p_outbreak_vec(i) = calculate_outbreak_prob(beta_tot,eta,prop_pop);
end

figure(1); hold on;
plot(R0_vec,p_outbreak_vec)
plot(R0_vec,max(1-1./R0_vec,0),'k--')

% Save results

save('../../Results/Figure_2/explore_R0_notest.mat','R0_vec','p_outbreak_vec')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')