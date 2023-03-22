function ongoing = outbreak_ongoing_condition_fun(I_indiv,R_indiv)
    
    % Evulates whether or not the simulated outbreak remains ongoing

    ongoing = any(I_indiv);
end