% Obtain baseline within-host, detection and infectiousness dynamics

clear all; close all; clc;

addpath('../../Functions/WH_model')

% Load inputs

load('../../Data/params_in.mat','R0','params_vec','prob_pos_l10V_fun','det_rel_inf','beta_undet_rel_fun','beta_rel_fun','tau_vec')

tau_inc = params_vec(5);
eta = params_vec(6);
prop_pop = params_vec(7);

% Run within-host model

[l10V_vec,prob_pos_int_vec,beta_undet_rel_int_vec] = WH_model_soln(params_vec,tau_vec,prob_pos_l10V_fun,beta_undet_rel_fun);
l10V_fun = @(tau)interp1(tau_vec,l10V_vec,tau,'linear',NaN);

figure(); hold on;
plot(tau_vec,l10V_vec)

% Calculate infectiousness profile (including overall scaling)

beta_notest_rel_int_vec = det_rel_inf*beta_undet_rel_int_vec+(1-det_rel_inf)*interp1(tau_vec,beta_undet_rel_int_vec,tau_inc);

beta_rel_tot_notest = beta_notest_rel_int_vec(end);
beta_scaling = R0/(beta_rel_tot_notest*eta*prop_pop);

beta_fun = @(l10V,detected)beta_scaling*beta_rel_fun(l10V,detected);

det_vec_notest = (tau_vec>=tau_inc);
beta_notest_vec = beta_fun(l10V_vec,det_vec_notest);

figure(); hold on;
plot(tau_vec,beta_notest_vec)

beta_tot_notest = beta_scaling*beta_rel_tot_notest;
beta_undet_int_vec = beta_scaling*beta_undet_rel_int_vec;

% Save results

save('../../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','l10V_fun','beta_notest_vec','beta_fun','prob_pos_int_vec','beta_undet_rel_int_vec','beta_tot_notest','beta_undet_int_vec')

rmpath('../../Functions/WH_model')