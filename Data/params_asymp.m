% Additional or modified parameters used in analysis accounting for
% entirely asymptomatic infections

clear all; close all; clc;

load('params_in.mat','params_WH_pop')

tau_inc_symp = params_WH_pop(5);

% Proportion of asymptomatic hosts

prob_asymp = 0.2;

% Matrix of parameter values for symptomatic and asymptomatic infected
% hosts

tau_inc_vec = [tau_inc_symp,inf];
eta_vec = [1,1];
prop_pop_vec = [1-prob_asymp,prob_asymp];

params_mat = [[params_WH_pop(1:4),params_WH_pop(1:4)];tau_inc_vec;eta_vec;prop_pop_vec];

% Relative infectiousness of asymptomatic hosts

prop_asymp_trans_vals = [0,0.08,0.2];

asymp_rel_tot_inf_vals = prop_asymp_trans_vals*(1-prob_asymp)./(prob_asymp*(1-prop_asymp_trans_vals));

save('params_asymp.mat','params_mat','asymp_rel_tot_inf_vals')