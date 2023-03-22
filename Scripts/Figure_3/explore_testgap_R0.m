% Calculate the outbreak risk for a range of (mean) intervals between
% successive antigen tests, considering different R0 values

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/Analytic/')

% Load inputs

load('../../Data/params_in.mat','params_vec','mean_test_gap','tau_vec','R0')
load('../../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','prob_pos_int_vec','beta_fun')

tau_inc = params_vec(5);
eta = params_vec(6);
prop_pop = params_vec(7); %trivially 1 in 1-group case

R0_default = 1.5;

% Loop over R0 values

R0_vec = [1.1,1.25,1.5,2,2.5];

mean_test_gap_vec = (0.01:0.01:7)';

p_outbreak_mat = zeros(size(mean_test_gap_vec));

for j = 1:length(R0_vec)
    
    R0 = R0_vec(j);

    % Loop over testing interval values

    for i = 1:length(mean_test_gap_vec)

        mean_test_gap = mean_test_gap_vec(i);

        p_det_vec = calculate_detection_probs(tau_vec,tau_inc,prob_pos_int_vec,mean_test_gap);

        beta_vec = (R0/R0_default)*beta_fun(l10V_vec,p_det_vec);
        beta_tot = trapz(tau_vec,beta_vec);

        p_outbreak_mat(i,j) = calculate_outbreak_prob(beta_tot,eta,prop_pop);
    end
end

figure(1); hold on;
plot(mean_test_gap_vec,p_outbreak_mat)

% Save results

save('../../Results/Figure_3/explore_testgap_R0.mat','R0_vec','mean_test_gap_vec','p_outbreak_mat')

rmpath('../../Data')
rmpath('../../Functions/Analytic/')