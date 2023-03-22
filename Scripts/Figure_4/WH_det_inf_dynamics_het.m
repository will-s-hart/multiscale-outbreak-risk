% Obtain within-host, detection and infectiousness dynamics accounting for
% heterogeneity in viral dynamics

clear all; close all; clc;

rng(7)

addpath('../../Functions/WH_model')

% Load inputs

load('../../Data/params_in.mat','R0','params_WH_pop','prob_pos_l10V_fun','det_rel_inf','beta_undet_rel_fun','beta_rel_fun')
load('../../Data/params_het_WH.mat','Omega','no_hosts','tau_vec_het')

% Sample individual WH parameters

params_WH_indiv = exp(mvnrnd(log(params_WH_pop),Omega,no_hosts)');

eta_indiv = ones(1,no_hosts);
prop_pop_indiv = repmat(1/no_hosts,1,no_hosts);
params_indiv = [params_WH_indiv;eta_indiv;prop_pop_indiv];

tau_inc_indiv = params_indiv(5,:);

% Run within-host model

[l10V_mat,prob_pos_int_mat,beta_undet_rel_int_mat] = WH_model_soln(params_indiv,tau_vec_het,prob_pos_l10V_fun,beta_undet_rel_fun);
l10V_indiv_fun = @(tau)interp1(tau_vec_het,l10V_mat,tau,'linear',NaN);

% Calculate infectiousness profile (including overall scaling)

beta_undet_rel_int_t_inc_indiv = zeros(1,no_hosts);

for i = 1:no_hosts
    beta_undet_rel_int_t_inc_indiv(i)=interp1(tau_vec_het,beta_undet_rel_int_mat(:,i),tau_inc_indiv(i));
end

beta_notest_rel_int_mat = det_rel_inf*beta_undet_rel_int_mat+(1-det_rel_inf)*beta_undet_rel_int_t_inc_indiv;

beta_rel_tot_notest_indiv = beta_notest_rel_int_mat(end,:);
beta_scaling = R0/sum(beta_rel_tot_notest_indiv.*eta_indiv.*prop_pop_indiv);

beta_fun = @(l10V,detected)beta_scaling*beta_rel_fun(l10V,detected);

det_mat_notest = (tau_vec_het>=tau_inc_indiv);
beta_notest_mat = beta_fun(l10V_mat,det_mat_notest);

beta_tot_notest_indiv = beta_scaling*beta_rel_tot_notest_indiv;

beta_mean_notest_vec = mean(beta_notest_mat,2); %expected infectiousness at each time since infection

% Select example hosts for plotting

hosts_plot = 11:15;
beta_notest_mat_plot = beta_notest_mat(:,hosts_plot);

figure(); hold on;
plot(tau_vec_het,beta_notest_mat_plot)
plot(tau_vec_het,beta_mean_notest_vec,'k--')

% Save results

save('../../Results/Figure_4/WH_det_inf_dynamics_het.mat','params_indiv','l10V_mat','beta_fun','prob_pos_int_mat','beta_tot_notest_indiv','beta_mean_notest_vec','beta_notest_mat_plot')

rmpath('../../Functions/WH_model')