clear;
close all;
clc;

% Macros

MC = 1000;                                                               % Size of the monte-carlo ensemble

M  = 100;                                                                   % Number of antennas at base station
K  = 150;                                                                   % Number of users at the cell 

% M = 50  & K = 75
% M = 100 & K = 150 

if K > M
    L_max = M;
else
    L_max = K-1;
end

N_ALG = 2;                                                                 % Number of algorithms for perform user scheduling
N_PRE = 2;
N_PA  = 2;

R = 500;
snr = 132;

bandwidth   = 20e6;
dl_ul_ratio = 0.5;

% Roots

root_load = '../../../../Google Drive/UFRJ/PhD/Codes/user-selection-with-large-scale-fading/Results/';
root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-selection-with-large-scale-fading/Figures/';

zero_pad_1 = '%03d';
zero_pad_2 = '%02d';

chn_type = 'ur_los';

% Loading data

sum_se_s = zeros(L_max,N_PRE,N_PA,N_ALG,MC);
min_se_s = zeros(L_max,N_PRE,N_PA,N_ALG,MC);

% load([root_load 'spectral_efficiency_all_L_' chn_type '_partial_csi_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
% load([root_load 'spectral_efficiency_all_L_' chn_type '_M_' sprintf(zero_pad_1,M) '_K_' sprintf(zero_pad_1,K) '_SNR_' num2str(snr) '_dB_MC_' num2str(MC) '_' sprintf(zero_pad_2,n_mc) '.mat']);
load([root_load 'se_all_L_ur_los_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC) '.mat']);
            
sum_se = reshape(sum(se,1),N_PRE,N_PA,MC);
min_se = reshape(min(se,[],1),N_PRE,N_PA,MC);

for l = 1:L_max    
     sum_se_s(l,:,:,:,:) = sum(se_s_all_L(1:l,l,:,:,:,:),1);
     min_se_s(l,:,:,:,:) = min(se_s_all_L(1:l,l,:,:,:,:),[],1);
end

avg_sum_thrgpt   = bandwidth*dl_ul_ratio*mean(sum_se,3);
avg_sum_thrgpt_s = bandwidth*dl_ul_ratio*mean(sum_se_s,5);

avg_min_thrgpt   = bandwidth*dl_ul_ratio*mean(min_se,3);
avg_min_thrgpt_s = bandwidth*dl_ul_ratio*mean(min_se_s,5);

L = ceil(K/5);

N_BIN = 100;

cdf_sum_thrgpt = zeros(N_BIN,N_PRE,N_PA,N_ALG);
cdf_min_thrgpt = zeros(N_BIN,N_PRE,N_PA,N_ALG);

edg_sum_thrgpt = zeros(N_BIN+1,N_PRE,N_PA,N_ALG);
edg_min_thrgpt = zeros(N_BIN+1,N_PRE,N_PA,N_ALG);

for n_alg = 1:N_ALG
    for n_pa = 1:N_PA
        for n_pre = 1:N_PRE
            [cdf_sum_thrgpt(:,n_pre,n_pa,n_alg),edg_sum_thrgpt(:,n_pre,n_pa,n_alg)] = histcounts(bandwidth*dl_ul_ratio*sum_se_s(L,n_pre,n_pa,n_alg,:),N_BIN,'normalization','cdf');
            [cdf_min_thrgpt(:,n_pre,n_pa,n_alg),edg_min_thrgpt(:,n_pre,n_pa,n_alg)] = histcounts(bandwidth*dl_ul_ratio*min_se_s(L,n_pre,n_pa,n_alg,:),N_BIN,'normalization','cdf');
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

OM = 1e-6;

savefig = 0;

legend_pa             = {'EP','MMF'};
legend_algo           = {'SOS','FRBS'};
legend_algo_plus_prec = {'SOS','FRBS','MRT','ZF'};

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
    switch n_pa
        case 1
            figure;
       
            set(gcf,'position',[0 0 800 600]);

            if K <= M
                % Plot for user selection algorithm legends
                plot(1:K,OM*[avg_sum_thrgpt_s(:,1,n_pa,1); avg_sum_thrgpt(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:K,OM*[avg_sum_thrgpt_s(:,1,n_pa,2); avg_sum_thrgpt(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                % Plot for precoding algorithm legends
                plot(1:K,OM*[avg_sum_thrgpt_s(:,1,n_pa,1); avg_sum_thrgpt(1,n_pa)],'-' ,'color',colours(8,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_thrgpt_s(:,2,n_pa,1); avg_sum_thrgpt(2,n_pa)],'--','color',colours(8,:),'linewidth',linewidth);
                % Plot for results
                plot(1:K,OM*[avg_sum_thrgpt_s(:,1,n_pa,1); avg_sum_thrgpt(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_thrgpt_s(:,1,n_pa,2); avg_sum_thrgpt(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_thrgpt_s(:,2,n_pa,1); avg_sum_thrgpt(2,n_pa)],'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_sum_thrgpt_s(:,2,n_pa,2); avg_sum_thrgpt(2,n_pa)],'--','color',colours(2,:),'linewidth',linewidth);
            else
                % Plot for user selection algorithm legends
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                % Plot for precoding algorithm legends
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,1,n_pa,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,2,n_pa,1),'--','color',colours(8,:),'linewidth',linewidth);
                % Plot for results
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,2,n_pa,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_sum_thrgpt_s(:,2,n_pa,2),'--','color',colours(2,:),'linewidth',linewidth);
            end
            
            xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
            ylabel('Throughput (Mbps)','fontname',fontname,'fontsize',fontsize);
    
            legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_4,'numcolumns',2);
            legend box off;
    
            set(gca,'fontname',fontname,'fontsize',fontsize);

            if K <= M
                xlim([1 K]);
            else
                xlim([1 L_max]);
            end
            
            if savefig == 1
                saveas(gcf,[root_save 'throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
                saveas(gcf,[root_save 'throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
                saveas(gcf,[root_save 'throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
            end

            figure;
    
            set(gcf,'position',[0 0 800 600]);
    
            % Plots for user selection algorihtm legends
            plot(OM*edg_sum_thrgpt(:,1,n_pa,1),[cdf_sum_thrgpt(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
            hold on;
            plot(OM*edg_sum_thrgpt(:,1,n_pa,2),[cdf_sum_thrgpt(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
            % Plots for precoding algorithm legends
            plot(OM*edg_sum_thrgpt(:,1,n_pa,1),[cdf_sum_thrgpt(:,1,n_pa,1); 1],'-' ,'color',colours(8,:),'linewidth',linewidth);
            plot(OM*edg_sum_thrgpt(:,2,n_pa,1),[cdf_sum_thrgpt(:,2,n_pa,1); 1],'--','color',colours(8,:),'linewidth',linewidth);
            % Plots for results
            plot(OM*edg_sum_thrgpt(:,1,n_pa,1),[cdf_sum_thrgpt(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
            plot(OM*edg_sum_thrgpt(:,1,n_pa,2),[cdf_sum_thrgpt(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
            plot(OM*edg_sum_thrgpt(:,2,n_pa,1),[cdf_sum_thrgpt(:,2,n_pa,1); 1],'--','color',colours(1,:),'linewidth',linewidth);
            plot(OM*edg_sum_thrgpt(:,2,n_pa,2),[cdf_sum_thrgpt(:,2,n_pa,2); 1],'--','color',colours(2,:),'linewidth',linewidth);
    
            xlabel('Throughput (Mbps)','fontname',fontname,'fontsize',fontsize);
            ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
            legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
            legend box off;
    
            set(gca,'fontname',fontname,'fontsize',fontsize);
    
            ylim([0 1]);
    
            if savefig == 1
                saveas(gcf,[root_save 'cdf_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
                saveas(gcf,[root_save 'cdf_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
                saveas(gcf,[root_save 'cdf_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
            end
        case 2
            figure;
    
            set(gcf,'position',[0 0 800 600]);
            
            if K <= M
                % Plot for user selection algorithm legends
                plot(1:K,OM*[avg_min_thrgpt_s(:,1,n_pa,1); avg_min_thrgpt(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:K,OM*[avg_min_thrgpt_s(:,1,n_pa,2); avg_min_thrgpt(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                % Plot for precoding algorithm legends
                plot(1:K,OM*[avg_min_thrgpt_s(:,1,n_pa,1); avg_min_thrgpt(1,n_pa)],'-' ,'color',colours(8,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_min_thrgpt_s(:,2,n_pa,1); avg_min_thrgpt(2,n_pa)],'--','color',colours(8,:),'linewidth',linewidth);
                % Plot for results
                plot(1:K,OM*[avg_min_thrgpt_s(:,1,n_pa,1); avg_min_thrgpt(1,n_pa)],'-' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_min_thrgpt_s(:,1,n_pa,2); avg_min_thrgpt(1,n_pa)],'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_min_thrgpt_s(:,2,n_pa,1); avg_min_thrgpt(2,n_pa)],'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:K,OM*[avg_min_thrgpt_s(:,2,n_pa,2); avg_min_thrgpt(2,n_pa)],'--','color',colours(2,:),'linewidth',linewidth);
            else
                % Plot for user selection algorithm legends
                plot(1:L_max,OM*avg_min_thrgpt_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                hold on;
                plot(1:L_max,OM*avg_min_thrgpt_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                % Plot for precoding algorithm legends
                plot(1:L_max,OM*avg_min_thrgpt_s(:,1,n_pa,1),'-' ,'color',colours(8,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_min_thrgpt_s(:,2,n_pa,1),'--','color',colours(8,:),'linewidth',linewidth);
                % Plot for results
                plot(1:L_max,OM*avg_min_thrgpt_s(:,1,n_pa,1),'-' ,'color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_min_thrgpt_s(:,1,n_pa,2),'-' ,'color',colours(2,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_min_thrgpt_s(:,2,n_pa,1),'--','color',colours(1,:),'linewidth',linewidth);
                plot(1:L_max,OM*avg_min_thrgpt_s(:,2,n_pa,2),'--','color',colours(2,:),'linewidth',linewidth);
            end
            
            xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);    
            ylabel('Min-throughput (Mbps)','fontname',fontname,'fontsize',fontsize);
        
            legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_2,'numcolumns',2);
            legend box off;
    
            set(gca,'fontname',fontname,'fontsize',fontsize);
    
            if K <= M
                xlim([1 K]);
            else
                xlim([1 L_max]);
            end
    
            if savefig == 1
                saveas(gcf,[root_save 'min_throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
                saveas(gcf,[root_save 'min_throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
                saveas(gcf,[root_save 'min_throughput_all_L_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
            end
    
            figure;
    
            set(gcf,'position',[0 0 800 600]);
    
            % Plots for user selection algorihtm legends
            plot(OM*edg_min_thrgpt(:,1,n_pa,1),[cdf_min_thrgpt(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
            hold on;
            plot(OM*edg_min_thrgpt(:,1,n_pa,2),[cdf_min_thrgpt(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
            % Plots for precoding algorithm legends
            plot(OM*edg_min_thrgpt(:,1,n_pa,1),[cdf_min_thrgpt(:,1,n_pa,1); 1],'-' ,'color',colours(8,:),'linewidth',linewidth);
            plot(OM*edg_min_thrgpt(:,2,n_pa,1),[cdf_min_thrgpt(:,2,n_pa,1); 1],'--','color',colours(8,:),'linewidth',linewidth);
            % Plots for results
            plot(OM*edg_min_thrgpt(:,1,n_pa,1),[cdf_min_thrgpt(:,1,n_pa,1); 1],'-' ,'color',colours(1,:),'linewidth',linewidth);
            plot(OM*edg_min_thrgpt(:,1,n_pa,2),[cdf_min_thrgpt(:,1,n_pa,2); 1],'-' ,'color',colours(2,:),'linewidth',linewidth);
            plot(OM*edg_min_thrgpt(:,2,n_pa,1),[cdf_min_thrgpt(:,2,n_pa,1); 1],'--','color',colours(1,:),'linewidth',linewidth);
            plot(OM*edg_min_thrgpt(:,2,n_pa,2),[cdf_min_thrgpt(:,2,n_pa,2); 1],'--','color',colours(2,:),'linewidth',linewidth);
    
            xlabel('Min-throughput (Mbps)','fontname',fontname,'fontsize',fontsize);
            ylabel('Cumulative distribution','fontname',fontname,'fontsize',fontsize);
    
            legend(legend_algo_plus_prec,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',2);
            legend box off;
    
            set(gca,'fontname',fontname,'fontsize',fontsize);
    
            ylim([0 1]);
    
            if savefig == 1
                saveas(gcf,[root_save 'cdf_min_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'fig');
                saveas(gcf,[root_save 'cdf_min_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'png');
                saveas(gcf,[root_save 'cdf_min_throughput_ur_los_' legend_pa{n_pa} '_M_' num2str(M) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr) '_dB_R_' num2str(R) '_MC_' num2str(MC)],'epsc2');
            end
    end
end 