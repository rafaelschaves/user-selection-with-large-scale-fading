clear;
close all;
clc;

% Roots

root_save = '../../../../Google Drive/UFRJ/PhD/Codes/user-selection-with-large-scale-fading/Figures/';

% SOS

A_sos = @(M,K,L) (M^3 + M^2 - M).*(L - 1).*(2*K - L) + (2*M - 1)*(K*L - L.*(L - 1))/2 + 2*M^2.*L.*(L - 1);
M_sos = @(M,K,L) 2*M^2.*(L - 1).*(2*M*K - M*L + L) + M*(2*K*L - L.*(L - 1));
D_sos = @(M,K,L) 2*M.*L;
S_sos = @(M,K,L) K*L - L.*(L - 1)/2;

% FRBS

A_frbs = @(M,K,L) 4*K^2*M + K*(K - 1)*(2*K + 5)/6 - L.*(L + 1).*(L - 1)/3;
M_frbs = @(M,K,L) (4*K^2*M + K^2 - K).*ones(1,length(L));
D_frbs = @(M,K,L) K^2 + K*(K + 1) - L.*(L + 1);
S_frbs = @(M,K,L) K*(K - 1)/2 - L.*(L + 1)/2;

M = 50;

K_1 = 10;
L_1 = 1:9;

K_2 = 25;
L_2 = 1:24;

flop_sos_1  = A_sos(M,K_1,L_1)  + M_sos(M,K_1,L_1)  + D_sos(M,K_1,L_1)  + S_sos(M,K_1,L_1);
flop_frbs_1 = A_frbs(M,K_1,L_1) + M_frbs(M,K_1,L_1) + D_frbs(M,K_1,L_1) + S_frbs(M,K_1,L_1);

flop_sos_2  = A_sos(M,K_2,L_2)  + M_sos(M,K_2,L_2)  + D_sos(M,K_2,L_2)  + S_sos(M,K_2,L_2);
flop_frbs_2 = A_frbs(M,K_2,L_2) + M_frbs(M,K_2,L_2) + D_frbs(M,K_2,L_2) + S_frbs(M,K_2,L_2);

% Ploting Figures

linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 30;

marker    = {'o','s','^','v','x'};
linestyle = {'-','--',':'};

savefig = 1;
 
legend_algo = {'SOS','FRBS'};

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

figure;

set(gcf,'position',[0 0 800 600]);


plot(L_1,flop_sos_1*1e-6 ,'-' ,'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(L_1,flop_frbs_1*1e-6,'-' ,'color',colours(2,:),'linewidth',linewidth);
plot(L_2,flop_sos_2*1e-6 ,'--','color',colours(1,:),'linewidth',linewidth);
plot(L_2,flop_frbs_2*1e-6,'--','color',colours(2,:),'linewidth',linewidth);

xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
ylabel('MFlops Count','fontname',fontname,'fontsize',fontsize);

legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1);
legend box off;

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 max(L_2)]);

axes('position',[.525 .275 .365 .365]);
box on;

plot(L_1,flop_frbs_1*1e-6,'-' ,'color',colours(2,:),'linewidth',linewidth);
hold on;
plot(L_2,flop_frbs_2*1e-6,'--','color',colours(2,:),'linewidth',linewidth);

set(gca,'fontname',fontname,'fontsize',fontsize);

xlim([1 max(L_2)]);
ylim([0 0.3]);

if savefig == 1
    saveas(gcf,[root_save 'computational_complexity_M_' num2str(M)],'fig');
    saveas(gcf,[root_save 'computational_complexity_M_' num2str(M)],'png');
    saveas(gcf,[root_save 'computational_complexity_M_' num2str(M)],'epsc2');
end

% figure;
% 
% set(gcf,'position',[0 0 800 600]);
% 
% yyaxis left
% 
% plot(L_1,flop_sos_1*1e-6 ,'-' ,'linewidth',linewidth);
% hold on
% plot(L_2,flop_sos_2*1e-6 ,'--' ,'linewidth',linewidth);
% 
% ylabel('MFlops Count','fontname',fontname,'fontsize',fontsize);
% 
% % ylim([0 2]);
% 
% yyaxis right
% 
% plot(L_1,flop_frbs_1*1e-6,'-' ,'linewidth',linewidth);
% hold on;
% plot(L_2,flop_frbs_2*1e-6,'--' ,'linewidth',linewidth);
% 
% ylabel('MFlops Count','fontname',fontname,'fontsize',fontsize);
% 
% % ylim([2.2 2.5]);
% 
% xlabel('Number of selected users','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_algo,'fontname',fontname,'fontsize',fontsize,'location',location_1,'numcolumns',3);
% legend box off;
% 
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% xlim([1 max(L_2)]);
% 
% if savefig == 1
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M) '_K_' num2str(K_1)],'fig');
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M) '_K_' num2str(K_1)],'png');
%     saveas(gcf,[root_save 'computational_complexity_M_' num2str(M) '_K_' num2str(K_1)],'epsc2');
% end