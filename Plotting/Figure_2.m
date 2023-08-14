% Produce the panels in Figure 2.

clear all; close all; clc;

% Change the following line of code to "save_PDF = true" to save the figure
% panels as PDF files. This functionality requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to be
% installed.

save_PDF = false;

addpath('Setup')
addpath('../Functions/Analytic')

% Load inputs

load('Setup/colorspec.mat','c1','c2','c3','c4','c5','c6','c7','c8','c9')

load('../Data/params_in.mat','l10V_inf_min','params_vec','mean_test_gap','tau_vec','R0')

tau_inc = params_vec(5);

load('../Results/Figure_2/WH_det_inf_dynamics.mat','l10V_vec','prob_pos_int_vec','beta_fun','beta_notest_vec')

load('../Results/Figure_2/explore_R0_notest.mat','R0_vec','p_outbreak_vec')
p_outbreak_vec_R0 = p_outbreak_vec;

load('../Results/Figure_2/explore_testgap_sim.mat','mean_test_gap_vec','p_outbreak_vec')
mean_test_gap_vec_sim  = mean_test_gap_vec;
p_outbreak_vec_mtg_sim = p_outbreak_vec;

load('../Results/Figure_2/explore_testgap.mat','mean_test_gap_vec','p_outbreak_vec','R0eff_vec')
p_outbreak_vec_mtg = p_outbreak_vec;
R0eff_vec_mtg = R0eff_vec;

% Create plots

for k = 1:6
    figsetup(k)
end

figure(1); hold on;
xlim([0,15])
ylim([-2,8])
p2 = plot(tau_inc*[1,1],[min(ylim),max(ylim)],'--','color',0.75*[1,1,1],'linewidth',3);
p3 = plot([min(xlim),max(xlim)],l10V_inf_min*[1,1],'--','color','k','linewidth',3);
p1 = plot(tau_vec,l10V_vec,'color',c1,'linewidth',3);
p4 = plot(-1,-1,'w');
xlabel('Time since infection (days)')
ylabel('log_{10}(viral load (copies/ml))')
l = legend([p1,p2,p3,p4],{'Viral load','Symptom onset time','Limit of infectiousness','(and antigen detection)'});
l.FontSize = 15;
l.Position = [0.4060    0.2200    0.3900    0.1550];

prob_det_vec = calculate_detection_probs(tau_vec,tau_inc,prob_pos_int_vec,mean_test_gap);

figure(2); hold on;
xlim([0,15])
ylim([0,1.025])
yticks(0:0.2:1)
p3 = plot(tau_inc*[1,1],[min(ylim),max(ylim)],'--','color',0.75*[1,1,1],'linewidth',3);
p1 = plot(tau_vec,(tau_vec>=tau_inc),'color',c1,'linewidth',3);
p2 = plot(tau_vec,prob_det_vec,':','color',c3,'linewidth',3);
xlabel('Time since infection (days)')
ylabel('Probability detected')
l = legend([p1,p2,p3],{'Without antigen testing','With testing','Symptom onset time'});
l.FontSize = 15;
l.Position = [0.4070    0.6340    0.3930    0.1180];

figure(3); hold on;
xlim([0,15])
ylim([0,0.7])
yticks(0:0.1:0.7)
p3 = plot(tau_inc*[1,1],[min(ylim),max(ylim)],'--','color',0.75*[1,1,1],'linewidth',3);
p1 = plot(tau_vec,beta_notest_vec,'color',c1,'linewidth',3);
p2 = plot(tau_vec,beta_fun(l10V_vec,prob_det_vec),':','color',c3,'linewidth',3);
xlabel('Time since infection (days)')
ylabel('Infectiousness (day^{-1})')
l = legend([p1,p2,p3],{'Without antigen testing','With testing (average)','Symptom onset time'});
l.FontSize = 15;
l.Position = [0.4070    0.6340    0.3930    0.1180];

figure(4); hold on;
p3 = plot(R0*[1,1],[min(ylim),max(ylim)],'--','color',0.75*[1,1,1],'linewidth',3);
p1 = plot(R0_vec,p_outbreak_vec_R0,'color',c1,'linewidth',3);
p2 = plot(R0_vec,max(1-1./R0_vec,0),'--','color','k','linewidth',3);
xlabel('{\itR}_{0}')
ylabel('Outbreak risk')
l = legend([p1,p2,p3],{'Multi-scale model','1-1/{\itR}_{0}','Default {\itR}_{0} value'});
l.FontSize = 15;
l.Position = [0.4770    0.2200    0.3210    0.1460];

figure(5); hold on;
plot(mean_test_gap_vec,100*(1-R0eff_vec_mtg/R0),'color',c1,'linewidth',3)
xlim([0,7])
xticks(0:7)
xlabel('Mean interval between tests (days)')
ylabel({'Proportion of transmissions';'prevented (%)'})

figure(6); hold on;
p1 = plot(mean_test_gap_vec,p_outbreak_vec_mtg,'color',c1,'linewidth',3);
p2 = plot(mean_test_gap_vec_sim,p_outbreak_vec_mtg_sim,'x','color',c2,'linewidth',3,'markersize',15);
xlim([0,7])
xticks(0:7)
ylim([0,0.6])
xlabel('Mean interval between tests (days)')
ylabel('Outbreak risk')
l = legend([p1,p2],{'Analytic','Stochastic simulations'});
l.FontSize = 15;
l.Position = [0.4130    0.2200    0.3850    0.0810];

for k = 1:6
    figsetup(k)
end

% Save panels as PDF files

if save_PDF
    figure(1); export_fig Figures/Figure_2/A.pdf -nocrop -transparent
    figure(2); export_fig Figures/Figure_2/B.pdf -nocrop -transparent
    figure(3); export_fig Figures/Figure_2/C.pdf -nocrop -transparent
    figure(4); export_fig Figures/Figure_2/D.pdf -nocrop -transparent
    figure(5); export_fig Figures/Figure_2/E.pdf -nocrop -transparent
    figure(6); export_fig Figures/Figure_2/F.pdf -nocrop -transparent
end

rmpath('Setup')
rmpath('../Functions/Analytic')