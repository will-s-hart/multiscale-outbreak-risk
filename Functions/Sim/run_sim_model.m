function out = run_sim_model(in)
    
    % Runs a realisation of a discrete-time, stochastic, individual-based
    % outbreak simulation model.

    % Extract model inputs
    
    no_hosts = in.no_hosts;
    inc_period_discr_indiv = in.inc_period_discr_indiv;
    infn_duration_discr_indiv = in.infn_duration_discr_indiv;
    t_max = in.t_max;
    warn_t_max = in.warn_t_max;
    run_initialisation_process = in.run_initialisation_process;  
    outbreak_ongoing_condition = in.outbreak_ongoing_condition;
    run_testing_process = in.run_testing_process;
    run_transmission_process = in.run_transmission_process;
    terminate_major_outbreak = in.terminate_major_outbreak;
    
    % Initialise model variables
    
    t = 0;
    t_index = 1;

    init_I = run_initialisation_process();
    no_I_init = sum(init_I);

    S_indiv = true(1,no_hosts); S_indiv(init_I) = false;
    U_indiv = false(1,no_hosts); U_indiv(init_I) = true;
    O_indiv = false(1,no_hosts);
    D_indiv = false(1,no_hosts);
    I_indiv = (U_indiv|D_indiv);
    R_indiv = false(1,no_hosts);
    
    infection_time_indiv = NaN(1,no_hosts); infection_time_indiv(init_I) = t;
    onset_time_indiv = NaN(1,no_hosts); onset_time_indiv(init_I) = t + inc_period_discr_indiv(init_I);
    recovery_time_indiv = NaN(1,no_hosts); recovery_time_indiv(init_I) = t + infn_duration_discr_indiv(init_I);
    
    t_vec = (0:t_max)';
    S_vec = [no_hosts-no_I_init;NaN(t_max,1)];
    U_vec = [no_I_init;NaN(t_max,1)];
    O_vec = [0;NaN(t_max,1)];
    D_vec = [0;NaN(t_max,1)];
    I_vec = [no_I_init;NaN(t_max,1)];
    R_vec = [0;NaN(t_max,1)];

    % Main loop over time

    while outbreak_ongoing_condition(I_indiv,R_indiv)
        
        t = t+1;
        t_index = t_index + 1;
        
        % Recoveries

        new_R = (I_indiv&(recovery_time_indiv==t));
        
        U_indiv(new_R) = false;
        O_indiv(new_R) = false;
        D_indiv(new_R) = false;
        I_indiv(new_R) = false;
        R_indiv(new_R) = true;

        % Individuals who developed symptoms the previous day automatically
        % enter the D class
        
        new_D_symp = O_indiv;

        O_indiv(new_D_symp) = false;
        D_indiv(new_D_symp) = true;
        
        % Detections from testing (must run testing process before symptom
        % onset process)
        
        new_D_test = run_testing_process(U_indiv,t,infection_time_indiv);

        U_indiv(new_D_test) = false;
        D_indiv(new_D_test) = true;

        % Detections from symptom onset (treated differently on first day
        % of symptoms to account for exact time of onset)

        new_O_symp = (U_indiv&(onset_time_indiv==t));
        
        U_indiv(new_O_symp) = false;
        O_indiv(new_O_symp) = true;

        % Transmissions
        
        new_I = run_transmission_process(S_indiv,I_indiv,O_indiv,D_indiv,t,infection_time_indiv);

        S_indiv(new_I) = false;
        U_indiv(new_I) = true;
        I_indiv(new_I) = true;

        infection_time_indiv(new_I) = t;
        onset_time_indiv(new_I) = t + inc_period_discr_indiv(new_I);
        recovery_time_indiv(new_I) = t + infn_duration_discr_indiv(new_I);

        % Update summary vectors

        S_vec(t_index) = sum(S_indiv);
        U_vec(t_index) = sum(U_indiv);
        O_vec(t_index) = sum(O_indiv);
        D_vec(t_index) = sum(D_indiv);
        I_vec(t_index) = sum(I_indiv);
        R_vec(t_index) = sum(R_indiv);
        
        if terminate_major_outbreak(S_indiv,I_indiv)
            break
        elseif (t == t_max)
            if warn_t_max
                warning('Reached maximum simulation time')
            end
            break
        end
    end
    
    % Format model outputs

    out.t_vec = t_vec(1:t_index);
    out.S_vec = S_vec(1:t_index);
    out.U_vec = U_vec(1:t_index);
    out.O_vec = O_vec(1:t_index);
    out.D_vec = D_vec(1:t_index);
    out.I_vec = I_vec(1:t_index);
    out.R_vec = R_vec(1:t_index);
end