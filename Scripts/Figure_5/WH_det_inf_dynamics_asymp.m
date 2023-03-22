% Obtain within-host, detection and infectiousness dynamics accounting for
% asymptomatic infected hosts

clear all; close all; clc;

addpath('../../Functions/WH_model')

% Load inputs

load('../../Data/params_in.mat','R0','prob_pos_l10V_fun','det_rel_inf','beta_undet_rel_fun','beta_rel_fun','tau_vec')
load('../../Data/params_asymp.mat','params_mat','asymp_rel_tot_inf_vals')

tau_inc_vec = params_mat(5,:);
eta_vec = params_mat(6,:);
prop_pop_vec = params_mat(7,:);

% Run within-host model

[l10V_mat,prob_pos_int_mat,beta_undet_rel_int_mat_init] = WH_model_soln(params_mat,tau_vec,prob_pos_l10V_fun,beta_undet_rel_fun);
l10V_vec_fun = @(tau)interp1(tau_vec,l10V_mat,tau,'linear',NaN);

figure(1); hold on;
plot(tau_vec,l10V_mat)

% Loop over values of the contribution of asymptomatic infected hosts to
% transmission

for i = 1:3
    
    asymp_rel_tot_inf = asymp_rel_tot_inf_vals(i);

    % Calculate infectiousness profile (including overall scaling)
    
    beta_undet_rel_int_t_inc_symp = interp1(tau_vec,beta_undet_rel_int_mat_init(:,1),tau_inc_vec(1));
    beta_rel_tot_notest_symp = det_rel_inf*beta_undet_rel_int_mat_init(end,1)+(1-det_rel_inf)*beta_undet_rel_int_t_inc_symp;
    beta_rel_tot_notest_asymp = beta_undet_rel_int_mat_init(end,1);
    
    beta_scaling_symp = R0/sum(beta_rel_tot_notest_symp*[1,asymp_rel_tot_inf].*eta_vec.*prop_pop_vec);
    beta_scaling_asymp = beta_scaling_symp*asymp_rel_tot_inf*(beta_rel_tot_notest_symp/beta_rel_tot_notest_asymp);
    
    beta_fun = @(l10V,detected,symp)(beta_scaling_symp*symp+beta_scaling_asymp*(1-symp)).*beta_rel_fun(l10V,detected);
    
    det_mat_notest = (tau_vec>=tau_inc_vec);
    beta_notest_mat = beta_fun(l10V_mat,det_mat_notest,[1,0]);
                
    figure(); hold on;
    plot(tau_vec,beta_notest_mat)
    
    beta_tot_notest_vec = [beta_scaling_symp*beta_rel_tot_notest_symp,beta_scaling_asymp*beta_rel_tot_notest_asymp];

    if i == 1
        beta_fun1 = beta_fun;
        prob_pos_int_mat1 = prob_pos_int_mat;
        beta_notest_mat1 = beta_notest_mat;
        beta_tot_notest_vec1 = beta_tot_notest_vec;
    elseif i == 2
        beta_fun2 = beta_fun;
        prob_pos_int_mat2 = prob_pos_int_mat;
        beta_notest_mat2 = beta_notest_mat;
        beta_tot_notest_vec2 = beta_tot_notest_vec;
    elseif i == 3
        beta_fun3 = beta_fun;
        prob_pos_int_mat3 = prob_pos_int_mat;
        beta_notest_mat3 = beta_notest_mat;
        beta_tot_notest_vec3 = beta_tot_notest_vec;
    end
end

% Save results

save('../../Results/Figure_5/WH_det_inf_dynamics_asymp.mat','l10V_mat','prob_pos_int_mat1','prob_pos_int_mat2','prob_pos_int_mat3','beta_fun1','beta_fun2','beta_fun3','beta_notest_mat1','beta_notest_mat2','beta_notest_mat3','beta_tot_notest_vec1','beta_tot_notest_vec2','beta_tot_notest_vec3')

rmpath('../../Functions/WH_model')