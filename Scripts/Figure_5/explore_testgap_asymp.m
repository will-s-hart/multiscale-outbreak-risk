% Calculate the outbreak risk for a range of (mean) intervals between
% successive antigen tests, accounting for asymptomatic infections

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_in.mat','tau_vec')
load('../../Data/params_asymp.mat','params_mat')
load('../../Results/Figure_5/WH_det_inf_dynamics_asymp.mat','l10V_mat','prob_pos_int_mat1','prob_pos_int_mat2','prob_pos_int_mat3','beta_fun1','beta_fun2','beta_fun3')

tau_inc_vec = params_mat(5,:);
eta_vec = params_mat(6,:);
prop_pop_vec = params_mat(7,:);

% Loop over values of the contribution of asymptomatic infected hosts to
% transmission

mean_test_gap_vec = (0.01:0.01:7)';

p_outbreak_mat = zeros(length(mean_test_gap_vec),3);
R0eff_mat = zeros(length(mean_test_gap_vec),3);

for j = 1:3
    
    if j == 1
        prob_pos_int_mat = prob_pos_int_mat1;
        beta_fun = beta_fun1;
    elseif j == 2
        prob_pos_int_mat = prob_pos_int_mat2;
        beta_fun = beta_fun2;
    elseif j == 3
        prob_pos_int_mat = prob_pos_int_mat3;
        beta_fun = beta_fun3;
    end
    
    % Loop over testing interval values

    for i = 1:length(mean_test_gap_vec)
    
        mean_test_gap = mean_test_gap_vec(i);
    
        p_det_mat = calculate_detection_probs(tau_vec,tau_inc_vec,prob_pos_int_mat,mean_test_gap);
        
        symp_mat = true(size(l10V_mat)); symp_mat(:,2) = false;
        beta_mat = beta_fun(l10V_mat,p_det_mat,symp_mat);
        beta_tot_indiv = trapz(tau_vec,beta_mat);
        
        R0eff_mat(i,j) = sum(beta_tot_indiv.*eta_vec.*prop_pop_vec);
        p_outbreak_mat(i,j) = calculate_outbreak_prob(beta_tot_indiv,eta_vec,prop_pop_vec);
    end
end

figure(1); hold on;
plot(mean_test_gap_vec,p_outbreak_mat)

% Save results

save('../../Results/Figure_5/explore_testgap_asymp.mat','mean_test_gap_vec','R0eff_mat','p_outbreak_mat')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')