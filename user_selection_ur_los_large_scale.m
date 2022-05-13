addpath('./functions/')

% Cheking deirectory

dir_save  = './results/selection/downlink/';
root_save = [dir_save 'se_2_all_L_'];

if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Checking variables

if ~exist('MC','var')
    MC = 5;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)
end

if ~exist('M','var')
    M = 50;                                                                % Number of antennas at the base station
end

if ~exist('K','var')
    K = 10;                                                                % Number of users at the cell
end

N_ALG = 2;
N_PRE = 2;
N_PA  = 2;

commcell.nAntennas       = M;                                        % Number of Antennas
commcell.nUsers          = K;                                        % Number of Users
commcell.radius          = 500;                                      % Cell's raidus (circumradius) in meters
commcell.bsHeight        = 32;                                       % Height of base station in meters
commcell.userHeight      = [1 2];                                    % Height of user terminals in meters ([min max])
commcell.nPaths          = 30;                                       % Number of Multipaths
commcell.frequency       = 1.9e9;                                    % Carrier frequency in Hz
commcell.meanShadowFad   = 0;                                        % Shadow fading mean in dB
commcell.stdDevShadowFad = 8;                                        % Shadow fading standard deviation in dB
commcell.city            = 'large';                                  % Type of city

linkprop.bsPower         = 10;                                       % in Watts
linkprop.userPower       = 0.2;                                      % in Watts
linkprop.AntennaGainBS   = 0;                                        % in dBi
linkprop.AntennaGainUser = 0;                                        % in dBi
linkprop.noiseFigureBS   = 9;                                        % in dB
linkprop.noiseFigureUser = 9 ;                                       % in dB
linkprop.bandwidth       = 20e6;                                     % in Hz

channel_type = 'ur-los';

[~,snr_db] = linkBudgetCalculation(linkprop);                        % SNR in dB
snr        = 10^((snr_db)/10);                                       % SNR

% Initialization

% algorithm_type = {'exhaustive search selection ep','exhaustive search selection mmf','semi-orthogonal selection','fr-based selection'};
algorithm_type = {'lsf ratio selection','correlation-based selection', 'ici-based selection','fr-based selection'};

if K > M
    L_max = M;
else
    L_max = K-1;
end

se         = zeros(K,N_PRE,N_PA,MC);
se_s_all_L = zeros(L_max,L_max,N_PRE,N_PA,N_ALG,MC);
S_set      = zeros(K,L_max,N_PRE,N_ALG,MC);
user_pos   = zeros(K,3,MC);
eta        = zeros(K,N_PRE);

for mc = 1:MC
    mc
    
    [H,beta,user_pos(:,:,mc)] = massiveMIMOChannel(commcell,channel_type);
    
    [~,eta(:,1)] = maxMinFairness(H,beta,snr,'algorithm 2');
    eta(:,2)     = (1/sum(1./beta))./beta;
    
    [se(:,1,1,mc),se(:,2,1,mc)] = DLspectralEfficiency(H,beta,snr,1/K);
    [se(:,1,2,mc),se(:,2,2,mc)] = DLspectralEfficiency(H,beta,snr,eta);
    
    for L = 1:L_max
        S_set_aux = zeros(L,N_PRE,N_ALG);
        eta_s     = zeros(L,N_PRE);
        
        for alg_idx = 1:N_ALG
            [H_s,S_aux] = userSelector(H,beta,snr,algorithm_type{alg_idx},'fixed',L,[]);
            
            % if alg_idx == 1
            %     H_aux = H_s(:,:,1);
            %     S_set_aux(:,:,alg_idx) = S_aux(:,1:2);
            % elseif alg_idx == 2
            %     H_aux = H_s(:,:,1);
            %     S_set_aux(:,:,alg_idx) = S_aux;
            % else
                H_aux = H_s;
                S_set_aux(:,:,alg_idx) = repmat(S_aux,1,N_PRE);
            % end
            
            S_set(S_set_aux(:,1),L,1,alg_idx,mc) = 1;
            S_set(S_set_aux(:,2),L,2,alg_idx,mc) = 1;
            
            beta_s = [beta(S_set_aux(:,1,alg_idx)) beta(S_set_aux(:,2,alg_idx))];
            
            [~,eta_s(:,1)] = maxMinFairness(H_aux,beta_s(:,1),snr,'algorithm 2');
            eta_s(:,2)     = (1/sum(1./beta_s(:,2)))./beta_s(:,2);
            
            [se_s_mrt_ep,se_s_zf_ep]   = DLspectralEfficiency(H_s,beta_s,snr,1/L);
            [se_s_mrt_mmf,se_s_zf_mmf] = DLspectralEfficiency(H_s,beta_s,snr,eta_s);
            
            se_s_all_L(:,L,1,1,alg_idx,mc) = [se_s_mrt_ep; zeros(L_max-L,1)];
            se_s_all_L(:,L,2,1,alg_idx,mc) = [se_s_zf_ep; zeros(L_max-L,1)];
            
            se_s_all_L(:,L,1,2,alg_idx,mc) = [se_s_mrt_mmf; zeros(L_max-L,1)];
            se_s_all_L(:,L,2,2,alg_idx,mc) = [se_s_zf_mmf; zeros(L_max-L,1)];
        end
    end
    
    save([root_save strrep(channel_type,'-','_') '_M_' num2str(M) '_K_' num2str(K) '_SNR_' num2str(round(snr_db)) '_dB_R_' num2str(commcell.radius) '_MC_' num2str(MC) '.mat'],'se','se_s_all_L','user_pos','S_set');
end