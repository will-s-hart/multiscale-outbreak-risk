function p_outbreak = calculate_outbreak_prob(beta_tot_indiv,eta_indiv,prop_pop_indiv)

    % Calculates the outbreak risk, accounting for the possibility of
    % multiple infection types
    
    R = sum(beta_tot_indiv.*eta_indiv.*prop_pop_indiv);
    eta_mean = sum(eta_indiv.*prop_pop_indiv);
    
    f = @(q)sum(eta_indiv'.*prop_pop_indiv'.*exp(-(1-q)*eta_mean.*beta_tot_indiv'))/eta_mean-q; %(1-(outbreak risk)) is the minimal non-negative root of this function
    
    options = optimset('Display','off');
    [p_nooutbreak,~,flag] = fzero(f,min(1/R,1),options);
    
    % Check convergence to correct outbreak risk

    checktol = 1e-6;

    if flag ~= 1
        if (flag == -3)&&(abs(R-1)<checktol)
            p_nooutbreak = fzero(f,1);
        else
            error('Unexpected issue finding outbreak probability')
        end
    end
    
    if abs(p_nooutbreak-0.5)>(0.5+checktol)
        error('Calculated probability not between 0 and 1')
    elseif (R>(1+checktol))&&(abs(p_nooutbreak-1)<checktol)
        error('May not have converged to non-trivial solution')
    end

    p_outbreak = 1-p_nooutbreak;
end