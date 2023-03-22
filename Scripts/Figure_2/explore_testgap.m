% Calculate the outbreak risk for a range of (mean) intervals between
% successive antigen tests

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_in.mat','params_vec','mean_test_gap','tau_vec')
load('../../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','prob_pos_int_vec','beta_fun')

tau_inc = params_vec(5);
eta = params_vec(6);
prop_pop = params_vec(7); %trivially 1 in 1-group case

% Loop over testing interval values

mean_test_gap_vec = (0.01:0.01:7)';

p_outbreak_vec = zeros(size(mean_test_gap_vec));
R0eff_vec = zeros(size(mean_test_gap_vec));

for i = 1:length(mean_test_gap_vec)
    
    mean_test_gap = mean_test_gap_vec(i);

    p_det_vec = calculate_detection_probs(tau_vec,tau_inc,prob_pos_int_vec,mean_test_gap);
    
    beta_vec = beta_fun(l10V_vec,p_det_vec);
    beta_tot = trapz(tau_vec,beta_vec);
    
    R0eff_vec(i) = beta_tot*eta*prop_pop;
    p_outbreak_vec(i) = calculate_outbreak_prob(beta_tot,eta,prop_pop);
end

figure(1); hold on;
plot(mean_test_gap_vec,p_outbreak_vec)
plot(mean_test_gap_vec,max(1-1./R0eff_vec,0),'k--')

% Save results

save('../../Results/Figure_2/explore_testgap.mat','mean_test_gap_vec','p_outbreak_vec','R0eff_vec')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')