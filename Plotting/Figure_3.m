% Produce the panels in Figure 3.

clear all; close all; clc;
 
% Change the following line of code to "save_PDF = true" to save the figure
% panels as PDF files. This functionality requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to be
% installed.

save_PDF = false;

addpath('Setup')

% Load inputs

load('Setup/colorspec.mat','c1','c2','c3','c4','c5','c6','c7','c8','c9')

load('../Data/params_in.mat','params_vec','R0','tau_vec')
load('../Data/params_det_inf.mat','det_rel_inf_vals')

tau_inc = params_vec(5,:);

load('../Results/Figure_2/WH_det_inf_dynamics.mat','beta_notest_vec')
load('../Results/Figure_2/explore_testgap.mat','R0eff_vec')

load('../Results/Figure_3/explore_testgap_R0.mat','R0_vec','mean_test_gap_vec','p_outbreak_mat')

mean_test_gap_vec_R0 = mean_test_gap_vec;
p_outbreak_mat_R0 = p_outbreak_mat;

load('../Results/Figure_3/inf_dynamics_det_inf.mat','beta_notest_vec1','beta_notest_vec2','beta_notest_vec3')
load('../Results/Figure_3/explore_testgap_det_inf.mat','mean_test_gap_vec','p_outbreak_mat','R0eff_mat')

mean_test_gap_vec_det_inf = mean_test_gap_vec;
p_outbreak_mat_det_inf = p_outbreak_mat;
R0eff_mat_det_inf = R0eff_mat;

% Create plots

for k = 1:6
    figsetup(k)
end

figure(1); hold on;
xlim([0,15])
p5 = plot(tau_vec,(R0_vec(5)/R0)*beta_notest_vec,'color',c5,'linewidth',3);
p4 = plot(tau_vec,(R0_vec(4)/R0)*beta_notest_vec,'color',c4,'linewidth',3);
p3 = plot(tau_vec,(R0_vec(3)/R0)*beta_notest_vec,'color',c3,'linewidth',3);
p2 = plot(tau_vec,(R0_vec(2)/R0)*beta_notest_vec,'color',c2,'linewidth',3);
p1 = plot(tau_vec,(R0_vec(1)/R0)*beta_notest_vec,'color',c1,'linewidth',3);
xlabel('Time since infection (days)')
ylabel('Infectiousness (day^{-1})')
l1_R0 = sprintf('{\\itR}_{0} = %.3g',R0_vec(1));
l2_R0 = sprintf('{\\itR}_{0} = %.3g',R0_vec(2));
l3_R0 = sprintf('{\\itR}_{0} = %.3g',R0_vec(3));
l4_R0 = sprintf('{\\itR}_{0} = %.3g',R0_vec(4));
l5_R0 = sprintf('{\\itR}_{0} = %.3g',R0_vec(5));
l = legend([p1,p4,p2,p5,p3],{l1_R0,l4_R0,l2_R0,l5_R0,l3_R0},'NumColumns',3);
l.FontSize = 15;
l.Position = [0.2170    0.7450    0.5990    0.1090];

figure(2); hold on;
plot(mean_test_gap_vec_R0,100*(1-R0eff_vec/R0),'color',c1,'linewidth',3)
xlim([0,7])
xticks(0:7)
xlabel('Mean interval between tests (days)')
ylabel({'Proportion of transmissions';'prevented (%)'})

figure(3); hold on;
p1 = plot(mean_test_gap_vec_R0,p_outbreak_mat_R0(:,1),'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_R0,p_outbreak_mat_R0(:,2),'color',c2,'linewidth',3);
p3 = plot(mean_test_gap_vec_R0,p_outbreak_mat_R0(:,3),'color',c3,'linewidth',3);
p4 = plot(mean_test_gap_vec_R0,p_outbreak_mat_R0(:,4),'color',c4,'linewidth',3);
p5 = plot(mean_test_gap_vec_R0,p_outbreak_mat_R0(:,5),'color',c5,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,1])
xlabel('Mean interval between tests (days)')
ylabel('Outbreak risk')
l = legend([p1,p4,p2,p5,p3],{l1_R0,l4_R0,l2_R0,l5_R0,l3_R0},'NumColumns',3);
l.FontSize = 15;
l.Position = [0.2170    0.7450    0.5990    0.1090];

figure(4); hold on;
xlim([0,15])
ylim([0,1.6])
yticks(0:0.4:1.6)
p1 = plot(tau_vec,beta_notest_vec1,'color',c1,'linewidth',3);
p2 = plot(tau_vec,beta_notest_vec2,'color',c2,'linewidth',3);
p3 = plot(tau_vec,beta_notest_vec3,'color',c3,'linewidth',3);
xlabel('Time since infection (days)')
ylabel('Infectiousness (day^{-1})')
l1_det_inf = ['{\it\alpha_{d}} = ',sprintf('%.2g',det_rel_inf_vals(1))];
l2_det_inf = ['{\it\alpha_{d}} = ',sprintf('%.2g',det_rel_inf_vals(2))];
l3_det_inf = ['{\it\alpha_{d}} = ',sprintf('%.2g',det_rel_inf_vals(3))];
l = legend([p1,p2,p3],{l1_det_inf,l2_det_inf,l3_det_inf});
l.FontSize = 15;
l.Position = [0.5870    0.6380    0.2110    0.1600];

figure(5); hold on;
p1 = plot(mean_test_gap_vec_det_inf,100*(1-R0eff_mat_det_inf(:,1)/R0),'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_det_inf,100*(1-R0eff_mat_det_inf(:,2)/R0),'color',c2,'linewidth',3);
p3 = plot(mean_test_gap_vec_det_inf,100*(1-R0eff_mat_det_inf(:,3)/R0),'color',c3,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,100])
xlabel('Mean interval between tests (days)')
yl =  ylabel({'Proportion of transmissions';'prevented (%)'});
yl.Position(1)=-0.735;
l = legend([p1,p2,p3],{l1_det_inf,l2_det_inf,l3_det_inf});
l.FontSize = 15;
l.Position = [0.5870    0.6380    0.2110    0.1600];

figure(6); hold on;
p1 = plot(mean_test_gap_vec_det_inf,p_outbreak_mat_det_inf(:,1),'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_det_inf,p_outbreak_mat_det_inf(:,2),'color',c2,'linewidth',3);
p3 = plot(mean_test_gap_vec_det_inf,p_outbreak_mat_det_inf(:,3),'color',c3,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,0.6])
xlabel('Mean interval between tests (days)')
ylabel('Outbreak risk')
l = legend([p1,p2,p3],{l1_det_inf,l2_det_inf,l3_det_inf});
l.FontSize = 15;
l.Position = [0.5870    0.2400    0.2110    0.1600];

for k = 1:6
    figsetup(k)
end

% Save panels as PDF files

if save_PDF
    figure(1); export_fig Figures/Figure_3/A.pdf -nocrop -transparent
    figure(2); export_fig Figures/Figure_3/B.pdf -nocrop -transparent
    figure(3); export_fig Figures/Figure_3/C.pdf -nocrop -transparent
    figure(4); export_fig Figures/Figure_3/D.pdf -nocrop -transparent
    figure(5); export_fig Figures/Figure_3/E.pdf -nocrop -transparent
    figure(6); export_fig Figures/Figure_3/F.pdf -nocrop -transparent
end

rmpath('Setup')