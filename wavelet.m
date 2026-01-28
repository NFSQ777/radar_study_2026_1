%% 离散小波变换 (DWT) 深度演示：分解与去噪
clear; close all; clc;

% 1. 生成仿真信号 (LFM + 噪声)
fs = 1000; T = 1; t = 0:1/fs:T-1/fs;
% LFM 信号
sig_clean = chirp(t, 20, 1, 100); 
% 加入较强的白噪声
rng(42); % 固定随机种子
noise = 1 * randn(size(t));
sig_noisy = sig_clean + noise;

%% 2. DWT 多级分解
% 选用 'db4' (Daubechies 4) 小波，它比 Haar 小波更平滑，适合处理雷达波形
wname = 'db4'; 
Level = 3; % 分解 3 层

% wavedec 函数返回：
% C: 所有系数拼接成的长向量 [cA3, cD3, cD2, cD1]
% L: 记录每个部分的长度结构
[C, L] = wavedec(sig_noisy, Level, wname);

% 提取各层系数以便观察
cA3 = appcoef(C, L, wname, 3); % 第3层近似 (低频)
cD3 = detcoef(C, L, 3);        % 第3层细节 (高频)
cD2 = detcoef(C, L, 2);        % 第2层细节 (更高频)
cD1 = detcoef(C, L, 1);        % 第1层细节 (最高频，主要是噪声)

%% 3. 可视化分解结果
figure('Position', [100, 50, 800, 800], 'Color', 'w');

subplot(5,1,1); plot(t, sig_noisy, 'Color', [0.6 0.6 0.6]);
title('原始含噪信号'); axis tight;

subplot(5,1,2); plot(cD1); title('Level 1 Detail (最高频 - 全是噪声)'); axis tight;
subplot(5,1,3); plot(cD2); title('Level 2 Detail (次高频)'); axis tight;
subplot(5,1,4); plot(cD3); title('Level 3 Detail (中高频)'); axis tight;
subplot(5,1,5); plot(cA3, 'LineWidth', 1.5); title('Level 3 Approximation (信号骨架)'); axis tight;
xlabel('采样点 (降采样后的)');

%% 4. 小波去噪 (手动实现软阈值)
% 策略：既然 cD1 大部分是噪声，我们对 cD 进行阈值处理
sigma = median(abs(cD1)) / 0.6745; % 经典的噪声强度估计公式
threshold = 3 * sigma;             % 设定通用阈值

% 对细节系数进行软阈值处理 (Soft Thresholding)
C_denoised = C;
% 注意：C 的结构是 [cA3, cD3, cD2, cD1]
% 我们保留 cA3 不动，只处理后面的 cD 部分
first_detail_idx = L(1) + 1; % cA3 之后的索引
C_denoised(first_detail_idx:end) = wthresh(C(first_detail_idx:end), 's', threshold);

% 5. IDWT 重构
sig_denoised = waverec(C_denoised, L, wname);

%% 6. 对比去噪效果
figure('Position', [500, 300, 800, 400], 'Color', 'w');
subplot(2,1,1); 
plot(t, sig_noisy, 'Color', [0.7 0.7 0.7]); hold on;
plot(t, sig_clean, 'k', 'LineWidth', 1);
title('原始信号 vs 含噪信号'); legend('含噪', '真值'); axis tight;

subplot(2,1,2);
plot(t, sig_denoised, 'b', 'LineWidth', 1.5); hold on;
plot(t, sig_clean, 'r--', 'LineWidth', 1);
title(['小波去噪后 (保留了LFM的波形结构)']); 
legend('小波去噪结果', '真值'); axis tight;