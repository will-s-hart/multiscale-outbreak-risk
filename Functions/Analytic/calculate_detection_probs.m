function p_det_mat = calculate_detection_probs(tau_vec,tau_inc_indiv,prob_pos_int_mat,mean_test_gap)
    
    % Calculates a matrix giving the probabilities of one or more
    % individuals having been detected by each time since infection in
    % tau_vec
    
    p_det_mat = 1-exp(-prob_pos_int_mat/mean_test_gap);

    symp_mat = (tau_vec>=tau_inc_indiv);
    p_det_mat(symp_mat) = 1;
end