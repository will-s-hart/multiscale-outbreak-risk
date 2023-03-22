% Baseline inputs used in our analyses

clear all; close all; clc;

% Population-level model parameters

R0 = 1.5;

% Population WH parameters

beta_pop = 1.43e-07;
gamma_pop = 5.64;
delta_pop = 1.21;
V0_pop = 1e-2;
tau_inc_pop = 4.6;

params_WH_pop = [beta_pop,gamma_pop,delta_pop,V0_pop,tau_inc_pop]';

% Other parameters and full parameter vector

eta = 1; %relative susceptibility
prop_pop = 1; %single group makes up entirety of population in baseline analysis

params_vec = [params_WH_pop;eta;prop_pop];

% Antigen testing model parameters

l10V_lod = 3.3;
test_meas_error = 0.866;

mean_test_gap = 2;

prob_pos_l10V_fun = @(l10V)normcdf(l10V_lod,l10V,test_meas_error,'upper'); %probability of positive test result as a function of viral load
prob_neg_l10V_fun = @(l10V)normcdf(l10V_lod,l10V,test_meas_error);

% Infectiousness model parameters

l10V_inf_min = 3.3;
beta_undet_rel_fun = @(l10V)max(l10V-l10V_inf_min,0); %unscaled relative infectiousness of an undetected host as a function of viral load

det_rel_inf = 1/3.9;
beta_rel_fun = @(l10V,detected)(1-(1-det_rel_inf)*detected).*beta_undet_rel_fun(l10V); %relative infectiousness accounting for detection status (the input 'detected' can be specified as a probability)

% Time grid on which to solve within-host model

tau_max = 30;
dtau = 0.01;
tau_vec = (0:dtau:tau_max)';
tau_vec = sort([tau_vec;tau_inc_pop-1e-12;tau_inc_pop+1e-12]); %extra points around discontinuity to improve numerical accuracy

% Parameters for discrete-time outbreak simulation model

tau_discr_vec = (0:(tau_max-1))';
no_hosts_tot_sim = 100000; %total number of individuals to discretise within-host dynamics for
no_hosts_per_sim = 1000; %number of individuals to select in each simulation
major_outbreak_threshold = no_hosts_per_sim/10; %total of number of infections at which an outbreak is to be considered major

% Save parameters

save('params_in.mat','R0','params_WH_pop','params_vec','prob_pos_l10V_fun','prob_neg_l10V_fun','mean_test_gap','beta_undet_rel_fun','det_rel_inf','beta_rel_fun','l10V_lod','l10V_inf_min','tau_vec','tau_discr_vec','no_hosts_tot_sim','no_hosts_per_sim','major_outbreak_threshold')