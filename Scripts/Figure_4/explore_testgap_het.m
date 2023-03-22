% Calculate the outbreak risk for a range of (mean) intervals between
% successive antigen tests, accounting for heterogeneity in viral dynamics

% The script WH_det_inf_dynamics_het.m must be run before this script, as
% the full matrices describing heterogeneous within-host dynamics are not
% included in the version of the relevant results file in the GitHub
% repository due to large file size (but can be reproduced by running
% WH_det_inf_dynamics_het.m).

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_het_WH.mat','tau_vec_het')
load('../../Results/Figure_4/WH_det_inf_dynamics_het.mat','params_indiv','l10V_mat','prob_pos_int_mat','beta_fun')

tau_inc_indiv = params_indiv(5,:);
eta_indiv = params_indiv(6,:);
prop_pop_indiv = params_indiv(7,:);

% Loop over testing interval values

mean_test_gap_vec = (0.01:0.01:7)';

p_outbreak_vec = zeros(size(mean_test_gap_vec));
R0eff_vec = zeros(size(mean_test_gap_vec));

parfor i = 1:length(mean_test_gap_vec)
    
    mean_test_gap = mean_test_gap_vec(i);

    p_det_mat = calculate_detection_probs(tau_vec_het,tau_inc_indiv,prob_pos_int_mat,mean_test_gap);
    
    beta_mat = beta_fun(l10V_mat,p_det_mat);
    beta_tot_indiv = trapz(tau_vec_het,beta_mat);
    
    R0eff_vec(i) = sum(beta_tot_indiv.*eta_indiv.*prop_pop_indiv);
    p_outbreak_vec(i) = calculate_outbreak_prob(beta_tot_indiv,eta_indiv,prop_pop_indiv);
end

figure(1); hold on;
plot(mean_test_gap_vec,p_outbreak_vec)

% Save results

save('../../Results/Figure_4/explore_testgap_het.mat','mean_test_gap_vec','p_outbreak_vec','R0eff_vec')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')