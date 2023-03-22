% Use a discrete-time, stochastic, individual-based outbreak simulation
% model to estimate the outbreak risk for a range of (mean) intervals
% between successive antigen tests

% The script discretise_WH_dynamics.m must be run before this script, as
% the results file containing discretised within-host dynamics is not
% included in the GitHub repository due to large file size (but can be
% reproduced by running discretise_WH_dynamics.m).

clear all; close all; clc;

rng(5)

addpath('../../Data')
addpath('../../Functions/Sim/')

% Load inputs

load('../../Data/params_in.mat','params_vec','det_rel_inf','no_hosts_tot_sim','no_hosts_per_sim','major_outbreak_threshold')
load('../../Results/Figure_2/discretised_WH_dynamics.mat','prob_neg_discr_mat','beta_undet_discr_mat','tau_inc_discr_indiv','onset_date_rel_inf_indiv','infn_duration_discr_indiv')

eta = params_vec(6);

t_max = 500;

% Set up condition for a major outbreak to have occurred

major_outbreak_found = @(S_indiv,I_indiv)((no_hosts_per_sim-sum(S_indiv))>=major_outbreak_threshold); %this function should take both S_indiv and I_indiv as inputs (even if not both used)
terminate_major_outbreak = major_outbreak_found; %terminates outbreak simulations upon the major outbreak threshold being reached

% Loop over testing interval values

mean_test_gap_vec = (1:7)';
no_reps_per_gap = 100000;
p_outbreak_vec = NaN(size(mean_test_gap_vec));

parfor i = 1:length(mean_test_gap_vec)
    
    mean_test_gap = mean_test_gap_vec(i)
    
    % Carry out multiple outbreak simulations

    major_outbreak_vec = false(1,no_reps_per_gap);
    
    tic
    
    for j = 1:no_reps_per_gap
        
        % Select individuals to include in the simulation and extract
        % individual data for included individuals

        hosts_incl = randperm(no_hosts_tot_sim,no_hosts_per_sim);
            
        tau_inc_discr_indiv_incl = tau_inc_discr_indiv(hosts_incl);
        onset_date_rel_inf_indiv_incl = onset_date_rel_inf_indiv(hosts_incl);
        infn_duration_discr_indiv_incl = infn_duration_discr_indiv(hosts_incl);
        eta_indiv_incl = repmat(eta,1,no_hosts_per_sim);
        
        prob_neg_discr_mat_incl = prob_neg_discr_mat(:,hosts_incl);
        beta_undet_discr_mat_incl = beta_undet_discr_mat(:,hosts_incl);
        
        % Set up inputs to outbreak simulation model

        prob_neg_indiv_incl_tpi_fun = @(curr,tau_indiv_curr)discr_lookup_indiv(prob_neg_discr_mat_incl,curr,tau_indiv_curr);
        beta_undet_indiv_incl_tpi_fun = @(curr,tau_indiv_curr)discr_lookup_indiv(beta_undet_discr_mat_incl,curr,tau_indiv_curr);
        beta_indiv_incl_tpi_fun = @(curr,tau_indiv_curr,O_indicator_curr,D_indicator_curr)(1-(1-det_rel_inf)*D_indicator_curr-(1-onset_date_rel_inf_indiv_incl(curr)).*O_indicator_curr).*beta_undet_indiv_incl_tpi_fun(curr,tau_indiv_curr);
        
        run_initialisation_process = @()run_initialisation_process_fun(eta_indiv_incl);
        
        tests_per_day_indiv_fun = @(curr,t)poissrnd(1/mean_test_gap,[1,sum(curr)]);
        run_testing_process = @(U_indiv,t,infection_time_indiv)run_testing_process_fun(U_indiv,t,infection_time_indiv(U_indiv),tests_per_day_indiv_fun,prob_neg_indiv_incl_tpi_fun);
        
        run_transmission_process = @(S_indiv,I_indiv,O_indiv,D_indiv,t,infection_time_indiv)run_transmission_process_fun(S_indiv,I_indiv,t,no_hosts_per_sim,infection_time_indiv(I_indiv),O_indiv(I_indiv),D_indiv(I_indiv),beta_indiv_incl_tpi_fun,eta_indiv_incl(S_indiv));
        
        in = struct;
        in.no_hosts = no_hosts_per_sim;
        in.inc_period_discr_indiv = tau_inc_discr_indiv_incl;
        in.infn_duration_discr_indiv = infn_duration_discr_indiv_incl;
        in.eta_indiv = eta_indiv_incl;
        in.t_max = t_max;
        in.warn_t_max = true;
        in.run_initialisation_process = run_initialisation_process;
        in.outbreak_ongoing_condition = @outbreak_ongoing_condition_fun;
        in.terminate_major_outbreak = terminate_major_outbreak;
        in.run_testing_process = run_testing_process;
        in.run_transmission_process = run_transmission_process;
        
        % Run model simulation and record whether or not a major outbreak
        % has occurred

        out = run_sim_model(in);

        major_outbreak_vec(j) = major_outbreak_found(out.S_vec(end),out.I_vec(end));    
    end
    
    toc
    
    % Calculate the outbreak risk

    p_outbreak_vec(i) = mean(major_outbreak_vec);
end

figure();
plot(mean_test_gap_vec,p_outbreak_vec)

% Save results

save('../../Results/Figure_2/explore_testgap_sim.mat','mean_test_gap_vec','p_outbreak_vec')

rmpath('../../Data')
rmpath('../../Functions/Sim/')