function new_I = run_transmission_process_fun(S_indiv,I_indiv,t,no_hosts,infection_time_I,O_I,D_I,beta_indiv_tpi_fun,eta_S)
    
    % Runs an iteration of the testing process for the outbreak simulation
    % model to select which individuals become infected during the current
    % day of the simulation.

    beta_indiv_I = beta_indiv_tpi_fun(I_indiv,t-infection_time_I,O_I,D_I);
    beta_tot_curr = sum(beta_indiv_I);
    lambda_curr = beta_tot_curr/no_hosts;

    prob_inf_S = 1 - exp(-eta_S*lambda_curr);
    
    inf_S = rand(size(prob_inf_S))<prob_inf_S;

    new_I = false(1,no_hosts);
    new_I(S_indiv) = inf_S;
end