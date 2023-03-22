% Values of the relative infectiousness of detected individuals used in
% analysis exploring the effect of this quantity

clear all; close all; clc;

det_rel_inf_vals = [0,1/3.9,0.5];

save('params_det_inf.mat','det_rel_inf_vals')