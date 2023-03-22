clear all; close all; clc;

rng(2)

addpath('../../Functions/WH_model')

% Load inputs

load('../../Data/params_in.mat','params_vec','prob_neg_l10V_fun','det_rel_inf','no_hosts_tot_sim','tau_vec','tau_discr_vec')
load('../../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','beta_undet_int_vec')

tau_inc = params_vec(5);

% Instantaneous negative test result probability

prob_neg_vec = prob_neg_l10V_fun(l10V_vec);

% Exact time of infection on infection date

day_start_to_inf_time_indiv = rand(1,no_hosts_tot_sim);

% Discretised incubation period (date of symptom onset for specified exact
% infection time)

tau_inc_discr_indiv = floor(day_start_to_inf_time_indiv+tau_inc);

% Discretised relative infectiousness profile, probability of testing
% negative at start of each calendar day, and relative infectiousness on
% day of symptom onset

beta_undet_discr_mat = NaN(length(tau_discr_vec),no_hosts_tot_sim);
prob_neg_discr_mat = NaN(length(tau_discr_vec),no_hosts_tot_sim);
onset_date_rel_inf_indiv = NaN(1,no_hosts_tot_sim);

for j = 1:no_hosts_tot_sim
    
    beta_undet_discr_mat(:,j) = interp1(tau_vec,beta_undet_int_vec,tau_discr_vec+1-day_start_to_inf_time_indiv(j),'linear',NaN)-interp1(tau_vec,beta_undet_int_vec,tau_discr_vec-day_start_to_inf_time_indiv(j),'linear',NaN);
    prob_neg_discr_mat(:,j) = interp1(tau_vec,prob_neg_vec,tau_discr_vec-day_start_to_inf_time_indiv(j),'linear',NaN);

    beta_undet_int_onset_day_start_j = interp1(tau_vec,beta_undet_int_vec,tau_inc_discr_indiv(j)-day_start_to_inf_time_indiv(j),'linear',NaN);
    beta_undet_int_onset_day_end_j = interp1(tau_vec,beta_undet_int_vec,tau_inc_discr_indiv(j)+1-day_start_to_inf_time_indiv(j),'linear',NaN);
    beta_undet_int_onset_exact_j = interp1(tau_vec,beta_undet_int_vec,tau_inc,'linear',NaN);

    onset_date_rel_inf_indiv(j) = (beta_undet_int_onset_exact_j-beta_undet_int_onset_day_start_j+det_rel_inf*(beta_undet_int_onset_day_end_j-beta_undet_int_onset_exact_j))/(beta_undet_int_onset_day_end_j-beta_undet_int_onset_day_start_j);
end

beta_undet_discr_mat(1,:) = 0; %can't be infectious on day of infection
prob_neg_discr_mat(1,:) = 1; %can't test positive at start of day 0

% Duration of infection (up to loss of infectiousness)

infn_duration_discr_indiv = splitapply(@(x)find(x,1,'last'),(diff(beta_undet_discr_mat)<0),1:size(beta_undet_discr_mat,2));

% Infectiousness profile accounting for responses to symptoms (without
% antigen testing)

det_discr_mat = (tau_discr_vec>tau_inc_discr_indiv); %(column vec)>=(row vec)
onset_date_indicator_mat = (tau_discr_vec==tau_inc_discr_indiv); %(column vec)==(row vec)
beta_notest_discr_mat = (1-(1-det_rel_inf)*det_discr_mat-(1-onset_date_rel_inf_indiv).*onset_date_indicator_mat).*beta_undet_discr_mat;

% Save results

save('../../Results/Figure_2/discretised_WH_dynamics.mat','prob_neg_discr_mat','beta_undet_discr_mat','beta_notest_discr_mat','tau_inc_discr_indiv','onset_date_rel_inf_indiv','infn_duration_discr_indiv','day_start_to_inf_time_indiv')

rmpath('../../Functions/WH_model')