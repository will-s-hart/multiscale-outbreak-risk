function init_I = run_initialisation_process_fun(eta_indiv)
    
    % Samples an initial infected individual for the outbreak simulation
    
    init_I_index = find(sum(eta_indiv)*rand<cumsum(eta_indiv),1,'first');
    
    init_I = false(size(eta_indiv));
    init_I(init_I_index) = true;
end