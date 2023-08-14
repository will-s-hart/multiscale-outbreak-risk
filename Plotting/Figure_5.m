% Produce the panels in Figure 5. The export_fig package (freely available
% at https://github.com/altmany/export_fig) is required to save the figure
% panels as PDF files.

clear all; close all; clc;

addpath('Setup')
addpath('../Functions/Analytic')

% Load inputs

load('Setup/colorspec.mat','c1','c2','c3','c4','c5','c6','c7','c8','c9')

load('../Data/params_in.mat','params_vec','R0','tau_vec')

load('../Results/Figure_5/WH_det_inf_dynamics_asymp.mat','beta_notest_mat_vals')

load('../Results/Figure_5/explore_testgap_asymp.mat','mean_test_gap_vec','p_outbreak_mat','R0eff_mat')
mean_test_gap_vec_asymp = mean_test_gap_vec;
p_outbreak_mat_asymp = p_outbreak_mat;
R0eff_mat_asymp = R0eff_mat;

load('../Results/Figure_2/explore_testgap.mat','mean_test_gap_vec','p_outbreak_vec','R0eff_vec')

% Create plots

for k = 1:3
figsetup(k)
end

figure(1); hold on;
xlim([0,15])
p3 = plot(tau_vec,beta_notest_mat_vals{4}(:,1),'color',0.9*c2,'linewidth',3);
p4 = plot(tau_vec,beta_notest_mat_vals{4}(:,2),'--','color',1-0.5*(1-c2),'linewidth',3);
p1 = plot(tau_vec,beta_notest_mat_vals{2}(:,1),'color',c1,'linewidth',3);
p2 = plot(tau_vec,beta_notest_mat_vals{2}(:,2),'--','color',1-0.5*(1-c1),'linewidth',3);
xlabel('Time since infection (days)')
ylabel('Infectiousness (day^{-1})')
l = legend([p1,p2,p3,p4],{'Symptomatic ({\itx_A} = 0.32)','Asymptomatic ({\itx_A} = 0.32)','Symptomatic ({\itx_A} = 2.77)','Asymptomatic ({\itx_A} = 2.77)'});
l.FontSize = 15;
l.Position = [0.3590    0.7850    0.4250    0.2110];

figure(2); hold on;
p1 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,1)/R0),'color',1-0.75*(1-1.7*c4),'linewidth',3);
p2 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,2)/R0),'color',c1,'linewidth',3);
p3 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,3)/R0),'color',c3,'linewidth',3);
p4 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,4)/R0),'color',c2,'linewidth',3);
p5 = plot(mean_test_gap_vec,100*(1-R0eff_vec/R0),'--','color','k','linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,60])
xlabel('Mean interval between tests (days)')
ylabel({'Proportion of transmissions';'prevented (%)'})
l = legend([p1,p2,p3,p4,p5],{'{\itx_A} = 0','{\itx_A} = 0.32','{\itx_A} = 1','{\itx_A} = 2.77','No asymptomatic infected hosts'});
l.FontSize = 15;
l.Position = [0.2660    0.7480    0.5180    0.2480];

figure(3); hold on;
p5 = plot(mean_test_gap_vec,p_outbreak_vec,'--','color','k','linewidth',3);
p1 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,1),'color',1-0.75*(1-1.7*c4),'linewidth',3);
p2 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,2),'color',c1,'linewidth',3);
p3 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,3),'color',c3,'linewidth',3);
p4 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,4),'color',c2,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,0.6])
xlabel('Mean interval between tests (days)')
ylabel('Outbreak risk')
l = legend([p1,p2,p3,p4,p5],{'{\itx_A} = 0','{\itx_A} = 0.32','{\itx_A} = 1','{\itx_A} = 2.77','No asymptomatic infected hosts'});
l.FontSize = 15;
l.Position = [0.2660    0.7480    0.5180    0.2480];

for k = 1:3
figsetup(k)
end

% Save panels as PDF files

figure(1); export_fig Figures/Figure_5/A.pdf -nocrop -transparent
figure(2); export_fig Figures/Figure_5/B.pdf -nocrop -transparent
figure(3); export_fig Figures/Figure_5/C.pdf -nocrop -transparent

rmpath('Setup')
rmpath('../Functions/Analytic')