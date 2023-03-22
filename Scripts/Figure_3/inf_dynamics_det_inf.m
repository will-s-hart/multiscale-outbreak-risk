clear all; close all; clc;

addpath('../../Functions/WH_model')

% Load inputs

load('../../Data/params_in.mat','R0','params_vec','beta_undet_rel_fun','tau_vec')
load('../../Data/params_det_inf.mat','det_rel_inf_vals')
load('../../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','beta_undet_rel_int_vec')

tau_inc = params_vec(5);
eta = params_vec(6);
prop_pop = params_vec(7);

% Loop over values of relative infectiousness of detected hosts

for i = 1:3
    
    det_rel_inf = det_rel_inf_vals(i);
    beta_rel_fun = @(l10V,detected)(1-(1-det_rel_inf)*detected).*beta_undet_rel_fun(l10V);

    % Calculate infectiousness profile (including overall scaling)
    
    beta_notest_rel_int_vec = det_rel_inf*beta_undet_rel_int_vec+(1-det_rel_inf)*interp1(tau_vec,beta_undet_rel_int_vec,tau_inc);

    beta_rel_tot_notest = beta_notest_rel_int_vec(end);
    beta_scaling = R0/(beta_rel_tot_notest*eta*prop_pop);
    
    beta_fun = @(l10V,detected)beta_scaling*beta_rel_fun(l10V,detected);
    
    det_vec_notest = (tau_vec>=tau_inc);
    beta_notest_vec = beta_fun(l10V_vec,det_vec_notest);
    
    beta_undet_int_vec = beta_scaling*beta_undet_rel_int_vec;
    
    figure(1); hold on;
    plot(tau_vec,beta_notest_vec)
    
    beta_tot_notest = beta_scaling*beta_rel_tot_notest;
    
    if i == 1
        beta_fun1 = beta_fun;
        beta_notest_vec1 = beta_notest_vec;
        beta_tot_notest1 = beta_tot_notest;
    elseif i == 2
        beta_fun2 = beta_fun;
        beta_notest_vec2 = beta_notest_vec;
        beta_tot_notest2 = beta_tot_notest;
    elseif i == 3
        beta_fun3 = beta_fun;
        beta_notest_vec3 = beta_notest_vec;
        beta_tot_notest3 = beta_tot_notest;
    end
end

% Save results

save('../../Results/Figure_3/inf_dynamics_det_inf.mat','beta_notest_vec1','beta_notest_vec2','beta_notest_vec3','beta_fun1','beta_fun2','beta_fun3')

rmpath('../../Functions/WH_model')