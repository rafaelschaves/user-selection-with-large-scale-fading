clear;
close all;
clc;

% Roots

root = '../../../Google Drive/UFRJ/PhD/Codes/user-selection-with-large-scale-fading/Results/';

% Loading data

load([root 'se_all_L_ur_los_M_50_K_25_SNR_132_dB_R_500_MC_100_01.mat']);

se_01         = se;
S_set_01      = S_set;
se_s_all_L_01 = se_s_all_L;

load([root 'se_all_L_ur_los_M_50_K_25_SNR_132_dB_R_500_MC_100_02.mat']);

se_02         = se;
S_set_02      = S_set;
se_s_all_L_02 = se_s_all_L;

clear se S_set se_s_all_L;

se         = se_01 + se_02;
S_set      = S_set_01 + S_set_02;
se_s_all_L = se_s_all_L_01 + se_s_all_L_02;

save([root 'se_all_L_ur_los_M_050_K_025_SNR_132_dB_R_0500_MC_0100_01.mat'],'se','se_s_all_L','S_set');