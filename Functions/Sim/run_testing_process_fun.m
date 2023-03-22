function new_D_test = run_testing_process_fun(U_indiv,t,infection_time_U,tests_per_day_indiv_fun,prob_neg_indiv_tpi_fun)
    
    % Runs an iteration of the testing process for the outbreak simulation
    % model to select which (previously undetected) infected hosts test
    % positive at the start of the current day of the simulation.
    
    daily_tests_U = tests_per_day_indiv_fun(U_indiv,t);
    tested_U = (daily_tests_U>0);
    daily_tests_testedU = daily_tests_U(tested_U);
    
    testedU_indiv = false(size(U_indiv));
    testedU_indiv(U_indiv) = tested_U;

    prob_neg_pertest_testedU = prob_neg_indiv_tpi_fun(testedU_indiv,t-infection_time_U(tested_U));
    prob_neg_overall_testedU = prob_neg_pertest_testedU.^daily_tests_testedU;
    
    pos_testedU = rand(size(prob_neg_overall_testedU))>prob_neg_overall_testedU;

    new_D_test = false(size(U_indiv));
    new_D_test(testedU_indiv) = pos_testedU;
end