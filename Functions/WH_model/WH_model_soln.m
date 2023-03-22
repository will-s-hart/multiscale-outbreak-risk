function [l10V_out,prob_pos_int_out,beta_undet_rel_int_out] = WH_model_soln(params_mat,t_obs,prob_pos_l10V_fun,beta_undet_rel_fun)

% Solves the target cell-limited within-host model described in our
% manuscript for one or more individuals (where the columns of params_mat
% give the parameters for each individual), also outputting cumulative
% integrals of the positive test probability and of the relative undetected
% infectiousness.

no_hosts = size(params_mat,2);

% Obtain model solution, accounting for possibility that the time of
% infection (time 0) may or may not be included in the vector of time
% points at which the solution is required.

t_0 = t_obs(t_obs==0);
t_f = t_obs(t_obs>0);

length_f = length(t_f);

options1 = odeset('AbsTol',1e-12,'RelTol',1e-12);

if length_f == 0
    y_f = [];
else
    [~,y_f1] = ode45(@(t,y) WH_RHS(y,params_mat,no_hosts),[0;t_f],WH_bc(params_mat),options1);
    y_f = y_f1(2:end,:);
    if length_f == 1
        y_f = y_f(end,:);
    end
end

if isempty(t_0)
    y_out = y_f;
else
    y_out = [WH_bc(params_mat)';y_f];
end

% Obtain model outputs

g_out = y_out(:,(no_hosts+1):(2*no_hosts));
h_out = y_out(:,(2*no_hosts+1):(3*no_hosts));
i_out = y_out(:,(3*no_hosts+1):(4*no_hosts));

l10V_out = log10(exp(1))*g_out;
prob_pos_int_out = h_out;
beta_undet_rel_int_out = i_out;

function dy = WH_RHS(vars,params_mat,n)
    
    % Right-hand side of the system of differential equations
    
    beta = params_mat(1,:)'; gamma = params_mat(2,:)'; delta = params_mat(3,:)';
    f = vars(1:n); g = vars((n+1):(2*n)); h = vars((2*n+1):(3*n)); i = vars((2*n+1):(3*n));
    
    dy = [-beta.*f.*exp(g);gamma.*f - delta;prob_pos_l10V_fun(log10(exp(1))*g);beta_undet_rel_fun(log10(exp(1))*g)];
end
function L_bc = WH_bc(params_mat)
    
    % Boundary conditions satisfied by each variable at time 0
    
    f0 = ones(size(params_mat,2),1); g0 = log(params_mat(4,:)'); h0 = zeros(size(params_mat,2),1); i0 = zeros(size(params_mat,2),1);
    
    L_bc = [f0;g0;h0;i0];
end
end