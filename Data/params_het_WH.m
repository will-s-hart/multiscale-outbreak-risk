% Additional or modified parameters used in analysis accounting for
% heterogeneity in within-host dynamics

clear all; close all; clc;

% Covariance matrix of random effects (log scale)

omega_lbeta = 1.33;
omega_lgamma = 0.145;
omega_ldelta = 0.541;
omega_lV0 = 0;
omega_ltau_inc = 0.293;

omega_vec_theta = [omega_lbeta,omega_lgamma,omega_ldelta,omega_lV0,omega_ltau_inc]';
omega_diag_theta = diag(omega_vec_theta);

corr_mat_theta = eye(5);

Omega = omega_diag_theta*corr_mat_theta*omega_diag_theta;

% Number of individuals from which to sample within-host dynamics

no_hosts = 10000;

% Time grid on which to solve within-host model

tau_max = 50;
dtau = 0.01;
tau_vec_het = (0:dtau:tau_max)';

save('params_het_WH.mat','Omega','no_hosts','tau_vec_het')