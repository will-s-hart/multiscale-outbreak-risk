% Calculate the outbreak risk for a range of (mean) intervals between
% successive antigen tests, considering different values of the relative
% infectiousness of detected hosts

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_in.mat','tau_vec','params_vec')
load('../../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','prob_pos_int_vec')
load('../../Results/Figure_3/inf_dynamics_det_inf.mat','beta_fun1','beta_fun2','beta_fun3');

tau_inc = params_vec(5);
eta = params_vec(6);
prop_pop = params_vec(7);

% Loop over values of the relative infectiousness of detected hosts

mean_test_gap_vec = (0.01:0.01:7)';

p_outbreak_mat = zeros(length(mean_test_gap_vec),3);
R0eff_mat = zeros(length(mean_test_gap_vec),3);

for j = 1:3
    
    if j == 1
        beta_fun = beta_fun1;
    elseif j == 2
        beta_fun = beta_fun2;
    elseif j == 3
        beta_fun = beta_fun3;
    end
    
    % Loop over testing interval values

    for i = 1:length(mean_test_gap_vec)
        
        mean_test_gap = mean_test_gap_vec(i);

        p_det_vec = calculate_detection_probs(tau_vec,tau_inc,prob_pos_int_vec,mean_test_gap);
        
        beta_vec = beta_fun(l10V_vec,p_det_vec);
        beta_tot = trapz(tau_vec,beta_vec);
        
        R0eff_mat(i,j) = beta_tot*eta*prop_pop;
        p_outbreak_mat(i,j) = calculate_outbreak_prob(beta_tot,eta,prop_pop);
    end
end

figure(1); hold on;
plot(mean_test_gap_vec,p_outbreak_mat)

% Save results

save('../../Results/Figure_3/explore_testgap_det_inf.mat','mean_test_gap_vec','R0eff_mat','p_outbreak_mat')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')