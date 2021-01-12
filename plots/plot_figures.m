clear;
close all;
clc;

% Macros

MC = 100;                                                                                                                                          % Size of the monte-carlo ensemble

M = 50;                                                                                                                                              % Number of antennas at base station
K = 10;                                                                                                                                              % Number of users at the cell 

if K > M
    L_max = M;
else
    L_max = K-1;
end

N_ALG = 4;                                                                                                                                            % Number of algorithms for perform user scheduling
N_PRE = 2;
N_PA  = 2;

R = 500;
snr = 132;

bandwidth   = 20e6;
dl_ul_ratio = 0.5;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-selection-with-large-scale-fading/Results/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-selection-with-large-scale-fading/Figures/';

% Loading data

sum_se_s = zeros(L_max,N_PRE,N_PA,N_ALG,MC);
min_se_s = zeros(L_max,N_PRE,N_PA,N_ALG,MC);

load([root_load 'se_all_L_ur_los_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC) '.mat']);

sum_se = reshape(sum(se,1),N_PRE,N_PA,MC);
min_se = reshape(min(se,[],1),N_PRE,N_PA,MC);

for l = 1:L_max    
     sum_se_s(l,:,:,:,:) = sum(se_s_all_L(1:l,l,:,:,:,:),1);
     min_se_s(l,:,:,:,:) = min(se_s_all_L(1:l,l,:,:,:,:),[],1);
end

avg_sum_se   = bandwidth*dl_ul_ratio*mean(sum_se,3);
avg_sum_se_s = bandwidth*dl_ul_ratio*mean(sum_se_s,5);

avg_min_se   = bandwidth*dl_ul_ratio*mean(min_se,3);
avg_min_se_s = bandwidth*dl_ul_ratio*mean(min_se_s,5);

L = ceil(K/2);

N_BIN = 25;

cdf_sum_se = zeros(N_BIN,N_PRE,N_PA,N_ALG);
cdf_min_se = zeros(N_BIN,N_PRE,N_PA,N_ALG);

edg_sum_se = zeros(N_BIN+1,N_PRE,N_PA,N_ALG);
edg_min_se = zeros(N_BIN+1,N_PRE,N_PA,N_ALG);

for n_alg = 1:N_ALG
    for n_pa = 1:N_PA
        for n_pre = 1:N_PRE
            [cdf_sum_se(:,n_pre,n_pa,n_alg),edg_sum_se(:,n_pre,n_pa,n_alg)] = histcounts(sum_se_s(L,n_pre,n_pa,n_alg,:),N_BIN,'normalization','cdf');
            [cdf_min_se(:,n_pre,n_pa,n_alg),edg_min_se(:,n_pre,n_pa,n_alg)] = histcounts(min_se_s(L,n_pre,n_pa,n_alg,:),N_BIN,'normalization','cdf');
        end
    end
end
 
% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker    = {'o','s','^'};
linestyle = {'-','--',':'};

savefig = 0;
plotse  = 0;
 
if M == 50
    OM = 1e-6;
elseif M == 100
    OM = 1e-6;
end

% SOS - Semi-orthogonal selection
% CBS - Correlation-based selection
% ICIBS - ICI-based selection

legend_pa             = {'EP','MMF'};
legend_algo           = {'ESEP','ESMMF','SOS','FRBS'};
legend_algo_plus_prec = {'ESEP','ESMMF','SOS','FRBS','MRT','ZF'};

um = {'(kbps)','(Mbps)','(Gbps)'};

location_1 = 'northwest';
location_2 = 'northeast';
location_3 = 'southwest';
location_4 = 'southeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840;
           0.0000 0.0000 0.0000];

for n_pa = 1:N_PA
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    if K <= M
        % Plot for user selection algorithm legends
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,1); avg_sum_se(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,2); avg_sum_se(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,3); avg_sum_se(1,n_pa)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,4); avg_sum_se(1,n_pa)],'-' ,'color',colours(4,:),'linewidth',linewidth);
        % Plot for precoding algorithm legends
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,1); avg_sum_se(1,n_pa)],'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,n_pa,1); avg_sum_se(2,n_pa)],'--','color',colours(8,:),'linewidth',linewidth);
        % Plot for results
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,1); avg_sum_se(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,2); avg_sum_se(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,3); avg_sum_se(1,n_pa)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,1,n_pa,4); avg_sum_se(1,n_pa)],'-' ,'color',colours(4,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,n_pa,1); avg_sum_se(2,n_pa)],'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,n_pa,2); avg_sum_se(2,n_pa)],'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,n_pa,3); avg_sum_se(2,n_pa)],'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_sum_se_s(:,2,n_pa,4); avg_sum_se(2,n_pa)],'--' ,'color',colours(4,:),'linewidth',linewidth);
    else
        % Plot for user selection algorithm legends
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,4),'-' ,'color',colours(4,:),'linewidth',linewidth);
        % Plot for precoding algorithm legends
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,2,n_pa,1),'--','color',colours(8,:),'linewidth',linewidth);
        % Plot for results
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,1,n_pa,4),'-' ,'color',colours(4,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,2,n_pa,1),'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,2,n_pa,2),'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,2,n_pa,3),'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*avg_sum_se_s(:,2,n_pa,4),'--' ,'color',colours(4,:),'linewidth',linewidth);
    end
    
    xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
    
    if plotse == 1
        ylabel('Sum-spectral efficiency','fontname',fontname,'fontsize',fontsize);
    else
        ylabel(['Average throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    end
    
    if K == 10
        legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_4,'numcolumns',2);
        legend box off;
    end
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if K <= M
        xlim([1 K]);
    else
        xlim([1 L_max]);
    end
    
    if savefig == 1
        if plotse == 1
            saveas(gcf,[root_save 'sum_se_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'sum_se_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'sum_se_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        else
            saveas(gcf,[root_save 'throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        end
    end
end

for n_pa = 1:N_PA
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    if K <= M
        % Plot for user selection algorithm legends
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,1); avg_min_se(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,2); avg_min_se(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,3); avg_min_se(1,n_pa)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,4); avg_min_se(1,n_pa)],'-' ,'color',colours(4,:),'linewidth',linewidth);
        % Plot for precoding algorithm legends
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,1); avg_min_se(1,n_pa)],'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,2,n_pa,1); avg_min_se(2,n_pa)],'--','color',colours(8,:),'linewidth',linewidth);
        % Plot for results
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,1); avg_min_se(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,2); avg_min_se(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,3); avg_min_se(1,n_pa)],'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,1,n_pa,4); avg_min_se(1,n_pa)],'-' ,'color',colours(4,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,2,n_pa,1); avg_min_se(2,n_pa)],'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,2,n_pa,2); avg_min_se(2,n_pa)],'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,2,n_pa,3); avg_min_se(2,n_pa)],'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*[avg_min_se_s(:,2,n_pa,4); avg_min_se(2,n_pa)],'--' ,'color',colours(4,:),'linewidth',linewidth);
    else
        % Plot for user selection algorithm legends
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        hold on;
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,4),'-' ,'color',colours(4,:),'linewidth',linewidth);
        % Plot for precoding algorithm legends
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,2,n_pa,1),'--','color',colours(8,:),'linewidth',linewidth);
        % Plot for results
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,3),'-' ,'color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,1,n_pa,4),'-' ,'color',colours(4,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,2,n_pa,1),'--','color',colours(1,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,2,n_pa,2),'--','color',colours(2,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,2,n_pa,3),'--','color',colours(3,:),'linewidth',linewidth);
        plot(1:K,OM*avg_min_se_s(:,2,n_pa,4),'--' ,'color',colours(4,:),'linewidth',linewidth);
    end
    
    xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
    
    if plotse == 1
        ylabel('Min-spectral efficiency','fontname',fontname,'fontsize',fontsize);
    else
        ylabel(['Average min-throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    end
    
    if K == 10
        legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_2,'numcolumns',2);
        legend box off;
    end
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    if K <= M
        xlim([1 K]);
    else
        xlim([1 L_max]);
    end
    
    if savefig == 1
        if plotse == 1
            saveas(gcf,[root_save 'min_se_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'min_se_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'min_se_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        else
            saveas(gcf,[root_save 'min_throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'min_throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'min_throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        end
    end
end

for n_pa = 1:N_PA
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    % Plots for user selection algorihtm legends
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,1),[cdf_sum_se(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,2),[cdf_sum_se(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,3),[cdf_sum_se(:,1,n_pa,3); 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,4),[cdf_sum_se(:,1,n_pa,4); 1],'-' ,'color',colours(4,:),'linewidth',linewidth);
    % Plots for precoding algorithm legends
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,1),[cdf_sum_se(:,1,n_pa,1); 1],'-' ,'color',colours(8,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,2,n_pa,1),[cdf_sum_se(:,2,n_pa,1); 1],'--','color',colours(8,:),'linewidth',linewidth);
    % Plots for results
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,1),[cdf_sum_se(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,2),[cdf_sum_se(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,3),[cdf_sum_se(:,1,n_pa,3); 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,1,n_pa,4),[cdf_sum_se(:,1,n_pa,4); 1],'-' ,'color',colours(4,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,2,n_pa,1),[cdf_sum_se(:,2,n_pa,1); 1],'--','color',colours(1,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,2,n_pa,2),[cdf_sum_se(:,2,n_pa,2); 1],'--','color',colours(2,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,2,n_pa,3),[cdf_sum_se(:,2,n_pa,3); 1],'--','color',colours(3,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_sum_se(:,2,n_pa,4),[cdf_sum_se(:,2,n_pa,4); 1],'--','color',colours(4,:),'linewidth',linewidth);
    
    if plotse == 1
        xlabel('Sum-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    else
        xlabel(['Throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    end
    
    ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
    if K == 10
        legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
        legend box off;
    end
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    ylim([0 1]);
    
    if savefig == 1
        if plotse == 1
            saveas(gcf,[root_save 'cdf_sum_se_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'cdf_sum_se_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'cdf_sum_se_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        else
            saveas(gcf,[root_save 'cdf_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'cdf_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'cdf_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        end
    end
end

for n_pa = 1:N_PA
    figure;
    
    set(gcf,'position',[0 0 800 600]);
    
    % Plots for user selection algorihtm legends
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,1),[cdf_min_se(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
    hold on;
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,2),[cdf_min_se(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,3),[cdf_min_se(:,1,n_pa,3); 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,4),[cdf_min_se(:,1,n_pa,4); 1],'-' ,'color',colours(4,:),'linewidth',linewidth);
    % Plots for precoding algorithm legends
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,1),[cdf_min_se(:,1,n_pa,1); 1],'-' ,'color',colours(8,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,2,n_pa,1),[cdf_min_se(:,2,n_pa,1); 1],'--','color',colours(8,:),'linewidth',linewidth);
    % Plots for results
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,1),[cdf_min_se(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,2),[cdf_min_se(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,3),[cdf_min_se(:,1,n_pa,3); 1],'-' ,'color',colours(3,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,1,n_pa,4),[cdf_min_se(:,1,n_pa,4); 1],'-' ,'color',colours(4,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,2,n_pa,1),[cdf_min_se(:,2,n_pa,1); 1],'--','color',colours(1,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,2,n_pa,2),[cdf_min_se(:,2,n_pa,2); 1],'--','color',colours(2,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,2,n_pa,3),[cdf_min_se(:,2,n_pa,3); 1],'--','color',colours(3,:),'linewidth',linewidth);
    plot(OM*bandwidth*dl_ul_ratio*edg_min_se(:,2,n_pa,4),[cdf_min_se(:,2,n_pa,4); 1],'--','color',colours(4,:),'linewidth',linewidth);
    
    if plotse == 1
        xlabel('Min-spectral efficiency (b/s/Hz)','fontname',fontname,'fontsize',fontsize);
    else
        xlabel(['Min-throughput ' um{abs(log10(OM))/3}],'fontname',fontname,'fontsize',fontsize);
    end
    
    ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
    if K == 10
        legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
        legend box off;
    end
    
    set(gca,'fontname',fontname,'fontsize',fontsize);
    
    ylim([0 1]);
    
    if savefig == 1
        if plotse == 1
            saveas(gcf,[root_save 'cdf_min_se_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'cdf_min_se_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'cdf_min_se_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        else
            saveas(gcf,[root_save 'cdf_min_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
            saveas(gcf,[root_save 'cdf_min_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
            saveas(gcf,[root_save 'cdf_min_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
        end
    end
end