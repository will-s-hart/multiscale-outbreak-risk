% Calculate the outbreak risk without regular antigen testing for a range
% of R0 values, accounting for heterogeneity in viral dynamics.

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_in.mat','R0')
load('../../Data/params_het_WH.mat','tau_vec_het')
load('../../Results/Figure_4/WH_det_inf_dynamics_het.mat','params_indiv','beta_tot_notest_indiv')

tau_inc_indiv = params_indiv(5,:);
eta_indiv = params_indiv(6,:);
prop_pop_indiv = params_indiv(7,:);

R0_default = R0;
beta_tot_indiv_default = beta_tot_notest_indiv;

% Loop over R0 values, each time calculating the outbreak risk

R0_vec = (0.01:0.01:4)';
p_outbreak_vec = zeros(size(R0_vec));

for i = 1:length(R0_vec)
        
    beta_tot_indiv = (R0_vec(i)/R0_default)*beta_tot_indiv_default;
    
    p_outbreak_vec(i) = calculate_outbreak_prob(beta_tot_indiv,eta_indiv,prop_pop_indiv);
end

figure(1); hold on;
plot(R0_vec,p_outbreak_vec)
plot(R0_vec,max(1-1./R0_vec,0),'k--')

% Save results

save('../../Results/Figure_4/explore_R0_notest_het.mat','R0_vec','p_outbreak_vec')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')