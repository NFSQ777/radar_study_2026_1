%% 脉压前 DWT 去噪效果对比仿真
clear; close all; clc;

%% 1. 参数设置
% ---------------------------------------
fs = 100e6;             % 采样率 100MHz
T = 10e-6;              % 脉冲宽度 10us
B = 20e6;               % 带宽 20MHz
K = B / T;              % 调频斜率
c = 3e8;                % 光速

% 时间轴
t = -T/2 : 1/fs : T/2 - 1/fs; 
N = length(t);

% 目标设置
R_target = 1000;        % 目标距离 (m)
tau = 2 * R_target / c; % 目标回波时延

%% 2. 信号生成
% ---------------------------------------
% 发射信号 (参考信号)
St = exp(1j * pi * K * t.^2);

% 生成回波 (带时延)
% 为了演示方便，我们构造一个全时间轴
T_total = 20e-6;        % 总接收窗长
t_raw = 0 : 1/fs : T_total - 1/fs;
N_raw = length(t_raw);

% 构造回波信号 (包含时延)
Sb = zeros(1, N_raw);
idx_delay = round(tau * fs);
% 放置目标回波 (防止越界)
if idx_delay + N <= N_raw
    Sb(idx_delay+1 : idx_delay+N) = St;
end

% 添加强噪声 (SNR = -15dB)
% 注意：在脉压前，-15dB 意味着信号完全淹没在噪声里，肉眼不可见
SNR_dB = -15;
noise = (randn(size(Sb)) + 1j*randn(size(Sb))) / sqrt(2);
signal_power = mean(abs(Sb).^2);
noise_power = signal_power / (10^(SNR_dB/10));
noise = noise * sqrt(noise_power);

Sr_noisy = Sb + noise; % 最终接收到的含噪信号

%% 3. 处理路径 A: 直接脉冲压缩 (Direct Pulse Compression)
% ---------------------------------------
% 频域匹配滤波
N_fft = 2^nextpow2(N_raw + N);
Hf = fft(St, N_fft);          % 参考信号频谱
Sf_A = fft(Sr_noisy, N_fft);  % 含噪信号频谱
pc_result_A = ifft(Sf_A .* conj(Hf)); % 脉压
% 截取有效部分
pc_result_A = pc_result_A(1:N_raw);

%% 4. 处理路径 B: DWT 去噪 -> 脉冲压缩
% ---------------------------------------
fprintf('正在进行 DWT 去噪...\n');

% [步骤1] 小波分解
wname = 'db4';  % 选用 Daubechies 4 小波，适合雷达波形
level = 4;      % 分解层数
% 注意：对实部和虚部分别处理 (雷达信号是复数)
[C_real, L_real] = wavedec(real(Sr_noisy), level, wname);
[C_imag, L_imag] = wavedec(imag(Sr_noisy), level, wname);

% [步骤2] 阈值计算与处理 (软阈值去噪)
% 使用 Donoho 的通用阈值公式: sigma * sqrt(2*log(N))
% 仅利用第一层细节系数 (cD1) 估计噪声强度
cD1_real = detcoef(C_real, L_real, 1);
sigma_real = median(abs(cD1_real)) / 0.6745;
thr_real = sigma_real * sqrt(2*log(N_raw));

cD1_imag = detcoef(C_imag, L_imag, 1);
sigma_imag = median(abs(cD1_imag)) / 0.6745;
thr_imag = sigma_imag * sqrt(2*log(N_raw));

% 对系数进行软阈值处理
% 修正说明：wdencmp 在使用 C,L 输入时必须有第8个参数 KEEPAPP。
% 我们设置为 1，表示保留低频近似系数不处理（只去噪高频细节）。
sr_denoised_real = wdencmp('gbl', C_real, L_real, wname, level, thr_real, 's', 1);
sr_denoised_imag = wdencmp('gbl', C_imag, L_imag, wname, level, thr_imag, 's', 1);

% [步骤3] 重构信号
Sr_denoised = sr_denoised_real + 1j * sr_denoised_imag;

% [步骤4] 对去噪后的信号进行脉冲压缩
Sf_B = fft(Sr_denoised, N_fft);
pc_result_B = ifft(Sf_B .* conj(Hf));
pc_result_B = pc_result_B(1:N_raw);

%% 5. 结果绘图与对比
% ---------------------------------------
% 归一化幅度以便对比
pc_A_norm = abs(pc_result_A) / max(abs(pc_result_A));
pc_B_norm = abs(pc_result_B) / max(abs(pc_result_B));
dist_axis = (t_raw * c / 2); % 距离轴

figure('Position', [100, 100, 1000, 700], 'Color', 'w');

% 图1: 时域信号对比
subplot(2,1,1);
plot(t_raw*1e6, real(Sr_noisy), 'Color', [0.8 0.8 0.8]); hold on;
plot(t_raw*1e6, real(Sr_denoised), 'b', 'LineWidth', 1.5);
plot(t_raw*1e6, real(Sb), 'r:', 'LineWidth', 1.5); % 真实纯净信号
legend('含噪信号 (SNR=-15dB)', 'DWT去噪后', '真实信号(参考)');
title('脉压前：时域信号对比');
xlabel('时间 (us)'); ylabel('幅度');
xlim([idx_delay/fs*1e6 - 2, idx_delay/fs*1e6 + 12]); % 放大看目标附近
grid on;

% 图2: 脉压结果对比 (dB图)
subplot(2,1,2);
plot(dist_axis, 20*log10(pc_A_norm + 1e-6), 'k--', 'LineWidth', 1); hold on;
plot(dist_axis, 20*log10(pc_B_norm + 1e-6), 'r', 'LineWidth', 1.5);
yline(-13.2, 'b:', '理论第一旁瓣'); 
legend('直接脉压', 'DWT去噪 + 脉压');
title('脉压后：距离像对比 (dB)');
xlabel('距离 (m)'); ylabel('归一化幅度 (dB)');
xlim([R_target-300, R_target+300]); 
ylim([-60, 0]);
grid on;

% 计算并打印信噪比增益
peak_val_A = max(abs(pc_result_A));
noise_floor_A = mean(abs(pc_result_A([1:idx_delay-100, idx_delay+200:end])));
snr_out_A = 20*log10(peak_val_A / noise_floor_A);

peak_val_B = max(abs(pc_result_B));
noise_floor_B = mean(abs(pc_result_B([1:idx_delay-100, idx_delay+200:end])));
snr_out_B = 20*log10(peak_val_B / noise_floor_B);

fprintf('---------------------------------\n');
fprintf('直接脉压输出 SNR: %.2f dB\n', snr_out_A);
fprintf('DWT+脉压输出 SNR: %.2f dB\n', snr_out_B);
fprintf('DWT 带来的 SNR 改善: %.2f dB\n', snr_out_B - snr_out_A);
fprintf('---------------------------------\n');