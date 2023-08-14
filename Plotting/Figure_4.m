% Produce the panels in Figure 4.

clear all; close all; clc;

% Change the following line of code to "save_PDF = true" to save the figure
% panels as PDF files. This functionality requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to be
% installed.

save_PDF = false;

addpath('Setup')

% Load inputs

load('Setup/colorspec.mat','c1','c2','c3','c4','c5','c6','c7','c8','c9')

load('../Data/params_in.mat','R0')
load('../Data/params_het_WH.mat','tau_vec_het')

load('../Results/Figure_4/WH_det_inf_dynamics_het.mat','beta_notest_mat_plot','beta_mean_notest_vec')

load('../Results/Figure_2/explore_R0_notest.mat','R0_vec','p_outbreak_vec')
R0_vec_homog = R0_vec;
p_outbreak_vec_R0_homog = p_outbreak_vec;

load('../Results/Figure_4/explore_R0_notest_het.mat','R0_vec','p_outbreak_vec')
R0_vec_het = R0_vec;
p_outbreak_vec_R0_het = p_outbreak_vec;

load('../Results/Figure_2/explore_testgap.mat','mean_test_gap_vec','p_outbreak_vec','R0eff_vec')
mean_test_gap_vec_homog = mean_test_gap_vec;
p_outbreak_vec_mtg_homog = p_outbreak_vec;
R0eff_vec_mtg_homog = R0eff_vec;

load('../Results/Figure_4/explore_testgap_het.mat','mean_test_gap_vec','p_outbreak_vec','R0eff_vec')
mean_test_gap_vec_het = mean_test_gap_vec;
p_outbreak_vec_mtg_het = p_outbreak_vec;
R0eff_vec_mtg_het = R0eff_vec;

% Create plots

for k = 1:4
    figsetup(k)
end

figure(1); hold on;
xlim([0,20])
xticks(0:5:20)
ylim([0,0.6])
yticks(0:0.1:0.7)
hosts_incl = (11:15);
p1 = plot(tau_vec_het,beta_notest_mat_plot,'linewidth',1.5);
p2 = plot(tau_vec_het,beta_mean_notest_vec,':','color','k','linewidth',3);
xlabel('Time since infection (days)')
ylabel('Infectiousness (day^{-1})')
l = legend(p2,{'Expected infectiousness'});
l.FontSize = 15;
l.Position = [0.3900    0.7540    0.4100    0.0440];

figure(2); hold on;
p4 = plot(R0*[1,1],[min(ylim),max(ylim)],'--','color',0.75*[1,1,1],'linewidth',3);
p1 = plot(R0_vec_homog,p_outbreak_vec_R0_homog,'color',c1,'linewidth',3);
p2 = plot(R0_vec_het,p_outbreak_vec_R0_het,'color',c2,'linewidth',3);
p3 = plot(R0_vec_het,max(1-1./R0_vec_het,0),'--','color','k','linewidth',3);
xlabel('{\itR}_{0}')
ylabel('Outbreak risk')
l = legend([p1,p2,p3,p4],{'Homogeneous','Heterogeneous','1-1/{\itR}_{0}','Default {\itR}_{0} value'});
l.FontSize = 15;
l.Position = [0.4980    0.2300    0.3000    0.1830];

figure(3); hold on;
p1 = plot(mean_test_gap_vec_homog,100*(1-R0eff_vec_mtg_homog/R0),'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_het,100*(1-R0eff_vec_mtg_het/R0),'color',c2,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,35])
xlabel('Mean interval between tests (days)')
ylabel({'Proportion of transmissions';'prevented (%)'})
l = legend([p1,p2],{'Homogeneous','Heterogeneous'});
l.FontSize = 15;
l.Position = [0.5120    0.7190    0.2880    0.0810];

figure(4); hold on;
p1 = plot(mean_test_gap_vec_homog,p_outbreak_vec_mtg_homog,'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_het,p_outbreak_vec_mtg_het,'color',c2,'linewidth',3);
xlim([0,7])
xticks(0:7)
ylim([0,0.6])
xlabel('Mean interval between tests (days)')
ylabel('Outbreak risk')
l = legend([p1,p2],{'Homogeneous','Heterogeneous'});
l.FontSize = 15;
l.Position = [0.5095    0.2310    0.2880    0.0810];

for k = 1:4
    figsetup(k)
end

% Save panels as PDF files

if save_PDF
    figure(1); export_fig Figures/Figure_4/A.pdf -nocrop -transparent
    figure(2); export_fig Figures/Figure_4/B.pdf -nocrop -transparent
    figure(3); export_fig Figures/Figure_4/C.pdf -nocrop -transparent
    figure(4); export_fig Figures/Figure_4/D.pdf -nocrop -transparent
end

rmpath('Setup')