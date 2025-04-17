%% ROC Curves under Different Modulation Schemes
% Signal Generation and Parameter Setup
% Enhanced with graph product signals and quadratic forms of different graphs

clc; clear all; close all;

%%%%% Channel Parameters
snr_b = -6;     % Start SNR (dB)
snr_e = -6;     % End SNR (dB)
snr = [snr_b:snr_e];
snr0 = 10.^(snr/10);
flag_sig = 0; % graph signals configuration 0,VPV+SQIS, 1,

% Graph Signal Parameters
k = 30;          % Number of samples
RB = 2000e3;     % Symbol rate (bps)
TB = 1/RB;       % Symbol duration
fc = 5*RB;       % Carrier frequency
Q = 10;          % Oversampling factor
fs = RB*Q;       % Sampling frequency
Amp = 1;         % Carrier amplitude
M = 200;         % Monte Carlo trials
N = 7;           % Quantization levels
flagch = 1;      % Channel type: 
                 % 0-AWGN, 1-frequency flat, 
                 % 2-frequency selective, 4-Laplace noise, 5-MIMO
GL = 7;          % Grouping length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVA Channel Model
MD = 300;        % Maximum Doppler shift
PD1 = 0;         % Path delay 1
APG1 = -30;      % Average path gain 1
PD2 = [0,30,150,310,370,710,1090,1730,2150]/10e9; % Path delays
APG2 = [0.0,-1.5,-1.4,-3.6,-0.6,-9.1,-7.0,-12.0,-16.9]; % Path gains

% Rayleigh Channel Objects
rayChan1 = comm.RayleighChannel('SampleRate', fs, ...
              'MaximumDopplerShift', MD, ...
              'DopplerSpectrum', doppler('Flat'), ...
              'PathDelays', PD1, ...
              'AveragePathGains', APG1, ...
              'NormalizePathGains', true);

rayChan2 = comm.RayleighChannel('SampleRate', fs, ...
              'MaximumDopplerShift', MD, ...
              'PathDelays', PD2, ...
              'AveragePathGains', APG2, ...
              'NormalizePathGains', true);

%%%%%%%%%%%%% Main Simulation Loop
% Initialize variables for runtime collection
runtime1 = zeros(1,M); % BM
runtime2 = zeros(1,M); % TE
runtime3 = zeros(1,M); % BR
runtime4 = zeros(1,M); % AUT
runtime5 = zeros(1,M); % EN
runtime6 = zeros(1,M); % OP
runtime7 = zeros(1,M); % GWAO
runtime8 = zeros(1,M); % GFT
runtime9 = zeros(1,M); % SP

for jj = 1:length(snr)
    for i = 1:M
        % QPSK Signal Generation
        t2 = 0:1/fs:TB-1/fs;                % Time sequence
        s2 = randi([0,3],1,k);               % Symbol sequence
        base_qpsk = pskmod(s2,4).';           % QPSK baseband
        base_qpsk_r = real(base_qpsk);
        base_qpsk_i = imag(base_qpsk);
        cr2 = Amp*cos(2*pi*fc*t2);           % I carrier
        ci2 = Amp*cos(2*pi*fc*t2-pi/2);      % Q carrier
        qpsk = base_qpsk_r*cr2 + base_qpsk_i*ci2;
        s_reshape = reshape(qpsk.',1,[]);     % Modulated signal
        K = length(s_reshape);                % Data length

        % Channel Impairments
        switch flagch
            case 0 % AWGN
                rs1 = sqrt(snr0)*s_reshape;
                r1 = awgn(rs1,snr,'measured');
                r0 = (r1 - rs1);
                
            case 1 % Frequency flat fading
                chan_out1 = real(rayChan1(s_reshape')');
                rs1 = sqrt(snr0)*chan_out1;
                r1 = awgn(rs1,snr,'measured');
                r0 = r1-rs1;
                
            case 2 % Frequency selective fading
                chan_out2 = real(rayChan2(s_reshape')');
                rs1 = sqrt(snr0)*chan_out2;
                r1 = awgn(rs1,snr,'measured');
                r0 = r1-rs1;
                
            case 3 % Interference + AWGN
                noiselength = length(s_reshape);
                cwfreq = rand;
                cwamp = rand;
                intersig = cwamp*cos(2*pi*[1:noiselength]*cwfreq);
                rs1 = sqrt(snr0)*s_reshape;
                r1 = awgn(rs1,snr,'measured') + intersig;
                r0 = r1-rs1-intersig;
                
            case 4 % Laplacian noise
                rs1 = s_reshape;
                [r1,r0] = add_laplacian_noise(rs1,snr);
                
            case 5 % MIMO
                Nr = 3; Nt = 2; L = 3;
                flag_Mimo = 0;
                [r0,r1,sigma] = MIMO_new2(s_reshape,Nt,Nr,L,snr,flag_Mimo);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Graph-Based Detection Algorithms
        % 1. Block Maxima Method (BM)
        tic;
        [G_bM_1,Adj_bM_1] = bmmax_s2g_new(r1,k*Q,N,GL,0);
        [H_bM_1] = degree_entropy(G_bM_1);
        h_bM_1(jj,i) = H_bM_1;
        runtime1(i) = toc;
        
        [G_bM_0,Adj_bM_0] = bmmax_s2g_new(r0,k*Q,N,GL,0);
        [H_bM_0] = degree_entropy(G_bM_0);
        h_bM_0(jj,i) = H_bM_0;

        % 2. Total Edge Method (TE)
        tic;
        theta = std(r0);
        [G2_1,E_1] = signal2graph_2020(r1,N,theta);
        runtime2(i) = toc;
        TE_1(jj,i) = E_1;
        [G2_0,E_0] = signal2graph_2020(r0,N,theta);
        TE_0(jj,i) = E_0;

        % 3. Block Range Method (BR)
        tic;
        [G3_1,Adj_br_1] = range_s2g_new(r1,k*Q,N,GL,0);
        [GI_br_1] = gini(G3_1);
        runtime3(i) = toc;
        gi_1(jj,i) = GI_br_1;
        [G3_0,Adj_br_0] = range_s2g_new(r0,k*Q,N,GL,0);
        [GI_br_0] = gini(G3_0);
        gi_0(jj,i) = GI_br_0;

        % 4. Autocorrelation Method (AUT)
        tic;
        [~,~,count_1] = corr_s2g_zeroeig_new(r1,k*Q,N);
        runtime4(i) = toc;
        cc_corr_1(jj,i) = count_1;
        [~,~,count_0] = corr_s2g_zeroeig_new(r0,k*Q,N);
        cc_corr_0(jj,i) = count_0;

        % 5. Energy Detection (EN)
        tic;
        energy_1 = sum(abs(r1).^2);
        runtime5(i) = toc;
        energy_en_1(i) = energy_1;
        energy_0 = sum(abs(r0).^2);
        energy_en_0(i) = energy_0;

        % 6. Original Power Spectrum (OP)
        tic;
        [~,~,SE_nb_0] = signal2graph_newnew(r0,k*Q,N,0);
        runtime6(i) = toc;
        se_nb_0(i) = SE_nb_0;
        [~,~,SE_nb_1] = signal2graph_newnew(r1,k*Q,N,0);
        se_nb_1(i) = SE_nb_1;

        % 7. Graph Weighted Aggregation Operator (GWAO)
        tic;
        [G_gwh_1,~,~,GWh_1] = signal2graph_PQChu1(r1,N,0);
        runtime7(i) = toc;
        gw_1(i) = GWh_1;
        [G_gwh_0,~,~,GWh_0] = signal2graph_PQChu1(r0,N,0);
        gw_0(i) = GWh_0;

        % 8. Graph Fourier Transform (GFT)
        tic;
        [~,~,~,GWh_11h,~,~] = signal2graph_PQChu(r1,N,0);
        runtime8(i) = toc;
        gwh_1(i) = GWh_11h;
        [~,~,~,GWh_00h,~,~] = signal2graph_PQChu(r0,N,0);
        gwh_0(i) = GWh_00h;

        % 9. Strong Product Graph (SP)
        tic;
        g_len = 3;
        flagd = 2; % Strong graph product
        N1 = 3;
        [~,GWhd_11,~,~,~,~,~] = signal2graph_Product_sig(r1,N1,g_len,flagd,flag_sig);
        runtime9(i) = toc;
        gwh_1rd(i) = GWhd_11;
        [~,GWhd_00,~,~,~,~,~] = signal2graph_Product_sig(r0,N1,g_len,flagd,flag_sig);
        gwh_0rd(i) = GWhd_00;
    end
end

%% ROC Curve Calculation
methods = {
    'BM',  -h_bM_0(:), -h_bM_1(:);    
    'TE',  -TE_0(:),  -TE_1(:);     
    'AUT', -cc_corr_0(:), -cc_corr_1(:);
    'EN',  energy_en_0(:), energy_en_1(:);
    'OP',  -se_nb_0(:), -se_nb_1(:);
    'GWAO', -gw_0(:), -gw_1(:);
    'GFT', gwh_0(:), gwh_1(:);
    'SP',  gwh_0rd(:), gwh_1rd(:)
};

% Initialize storage
tpr = cell(size(methods,1),1);
fpr = cell(size(methods,1),1);
auc = zeros(size(methods,1),1);

for m = 1:size(methods,1)
    try
        targets = [zeros(size(methods{m,2})); ones(size(methods{m,3}))];
        outputs = [methods{m,2}; methods{m,3}];
        [X,Y,~,AUC] = perfcurve(targets, outputs, 1);
        fpr{m} = X;
        tpr{m} = Y;
        auc(m) = AUC;
    catch ME
        fprintf('Error processing method %s: %s\n', methods{m,1}, ME.message);
        fpr{m} = linspace(0,1,100)';
        tpr{m} = linspace(0,1,100)';
        auc(m) = 0.5;
    end
end

%% Plot ROC Curves
figure;
set(gcf,'Position',[100 100 800 600]);
hold on;

% Line styles and colors for all 8 methods
line_styles = {
    '-r',  'BM',   2;      % red solid line - Block Maxima
    ':b',  'TE',   1.5;    % blue dotted line - Total Edge  
    '--k','AUT',  1.5;     % black dashed line - Autocorrelation
    '-.m','EN',   1;       % magenta dash-dot - Energy Detection
    '--g','OP',   2;       % green dashed line - Original Power Spectrum
    ':c','GWAO', 1.5;      % cyan dotted line - Graph Weighted Aggregation Operator
    '-g', 'GFT',  2;       % green solid line - Graph Fourier Transform
    '-c', 'SP',   2.5      % cyan solid line - Strong Product Graph
};

% Plot each method
for idx = 1:size(line_styles,1)
    m = find(strcmp(methods(:,1), line_styles{idx,2}));
    if ~isempty(m)
        plot(fpr{m}, tpr{m}, line_styles{idx,1}, ...
            'LineWidth', line_styles{idx,3}, ...
            'DisplayName', methods{m,1});
    end
end

% Formatting
hold off;
xlabel('False Positive Rate (P_{fa})','FontSize',12);
ylabel('True Positive Rate (P_{d})','FontSize',12);
title('ROC Curves of Different Detection Algorithms','FontSize',14);
legend('Location','southeast','FontSize',10);
grid on;
set(gca,'FontSize',11);
xlim([0 1]); ylim([0 1]);

%% Performance Table
% Collect runtime averages
runtime_methods = {
    mean(runtime1),   % BM
    mean(runtime2),   % TE
    mean(runtime4),   % AUT
    mean(runtime5),   % EN
    mean(runtime6),   % OP
    mean(runtime7),   % GWAO
    mean(runtime8),   % GFT
    mean(runtime9)    % SP
};

% Print table
fprintf('\n%-6s %-8s %-12s\n','Method','AUC','Time (ms)');
fprintf('-----------------------------\n');
for m = 1:size(methods,1)
    avg_time = runtime_methods{m} * 1000;
    fprintf('%-6s %.4f   %.2f\n', methods{m,1}, auc(m), avg_time);
end