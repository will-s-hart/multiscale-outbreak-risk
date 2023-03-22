% Produce the panels in Figure 5. The export_fig package (freely available
% at https://github.com/altmany/export_fig) is required to save the figure
% panels as PDF files.

clear all; close all; clc;

addpath('Setup')
addpath('../Functions/Analytic')

% Load inputs

load('Setup/colorspec.mat','c1','c2','c3','c4','c5','c6','c7','c8','c9')

load('../Data/params_in.mat','params_vec','R0','tau_vec')

load('../Results/Figure_5/WH_det_inf_dynamics_asymp.mat','beta_notest_mat2')

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
ylim([0,0.75])
yticks(0:0.25:0.75)
p1 = plot(tau_vec,beta_notest_mat2(:,1),'color',c1,'linewidth',3);
p2 = plot(tau_vec,beta_notest_mat2(:,2),'color',c2,'linewidth',3);
xlabel('Time since infection (days)')
ylabel('Infectiousness (day^{-1})')
l = legend([p1,p2],{'Symptomatic','Asymptomatic'});
l.FontSize = 15;
l.Position = [0.5200    0.7510    0.2760    0.0810];

figure(2); hold on;
p1 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,1)/R0),'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,2)/R0),'color',c2,'linewidth',3);
p3 = plot(mean_test_gap_vec_asymp,100*(1-R0eff_mat_asymp(:,3)/R0),'color',c3,'linewidth',3);
p4 = plot(mean_test_gap_vec,100*(1-R0eff_vec/R0),'--','color','k','linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,40])
xlabel('Mean interval between tests (days)')
ylabel({'Proportion of transmissions';'prevented (%)'})
l = legend([p1,p2,p3,p4],{'0% asymptomatic transmissions','8% asymptomatic transmissions','20% asymptomatic transmissions','No asymptomatic infected hosts'});
l.FontSize = 15;
l.Position = [0.2615    0.7490    0.5380    0.1550];

figure(3); hold on;
p4 = plot(mean_test_gap_vec,p_outbreak_vec,'--','color','k','linewidth',3);
p1 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,1),'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,2),'color',c2,'linewidth',3);
p3 = plot(mean_test_gap_vec_asymp,p_outbreak_mat_asymp(:,3),'color',c3,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,0.6])
xlabel('Mean interval between tests (days)')
ylabel('Outbreak risk')
l = legend([p1,p2,p3,p4],{'0% asymptomatic transmissions','8% asymptomatic transmissions','20% asymptomatic transmissions','No asymptomatic infected hosts'});
l.FontSize = 15;
l.Position = [0.2615    0.7490    0.5380    0.1550];

for k = 1:3
figsetup(k)
end

% Save panels as PDF files

figure(1); export_fig Figures/Figure_5/A.pdf -nocrop -transparent
figure(2); export_fig Figures/Figure_5/B.pdf -nocrop -transparent
figure(3); export_fig Figures/Figure_5/C.pdf -nocrop -transparent

rmpath('Setup')
rmpath('../Functions/Analytic')