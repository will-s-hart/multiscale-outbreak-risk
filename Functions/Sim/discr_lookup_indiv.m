function values_indiv = discr_lookup_indiv(values_mat_full,incl,tau_indiv_incl)
    
    % Evaluates the value of a quantity, recorded at integer timepoints for
    % a population of individuals, at a single time-point for each of some
    % subset of the population (where the query time-point may differ
    % between individuals).

    % values_mat must be an (m*n) matrix, where values_mat(i,j) represents
    % the value of the qunatity of interest for individual j, evaluated at
    % time since infection (i-1), incl must be a (1*n) logical vector
    % representing the included individuals, and tau_indiv_incl must either
    % be a scalar or (1*l) vector, where tau_indiv_incl(j) represents the
    % query time for included individual j.
    
    values_mat_incl = values_mat_full(:,incl);

    no_hosts = size(values_mat_incl,2);
    tau_max = size(values_mat_incl,1)-1;

    if isscalar(tau_indiv_incl)
        tau_indiv_incl = repmat(tau_indiv_incl,1,no_hosts);
    end

    in_range = (tau_indiv_incl>=0)&(tau_indiv_incl<=tau_max);
    
    inds = sub2ind([tau_max+1,no_hosts],tau_indiv_incl(in_range)+1,find(in_range));

    values_indiv = NaN(1,no_hosts);
    values_indiv(in_range) = values_mat_incl(inds);
end