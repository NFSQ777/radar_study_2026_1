clc; clear; close all;

%% 1. 雷达系统参数设置
fs = 100e6;             % 采样率 100MHz
T_pulse = 20e-6;        % 脉冲宽度 20us
B = 30e6;               % 信号带宽 30MHz
K = B / T_pulse;        % 调频斜率
R_max = 3000;           % 最大仿真距离 (m)
c = 3e8;                % 光速

% 目标参数
R_target = 1000;        % 目标距离 1000m
RCS = 1;                % 目标RCS (归一化幅度)

%% 2. 生成 LFM 发射信号 & 目标回波
t_pulse = 0 : 1/fs : T_pulse - 1/fs;
s_tx = exp(1j * pi * K * t_pulse.^2); % LFM 参考信号

% 计算回波延迟
tau = 2 * R_target / c;
N_delay = round(tau * fs);

% 构建接收时间轴 (覆盖最大距离)
t_max = 2 * R_max / c;
N_total = round(t_max * fs);
rx_signal = zeros(1, N_total);

% 填入目标回波 (模拟延迟)
rx_signal(N_delay : N_delay + length(s_tx) - 1) = RCS * s_tx;

%% 3. 添加压制型干扰 (Suppressive Jamming)
% 获取由于回波叠加可能导致变长后的实际信号长度
N_actual = length(rx_signal); 

% 关键参数：干信比 (Jamming-to-Signal Ratio, JSR)
JSR_dB = 20;            
SNR_dB = 10;            

% --- A. 计算信号功率 ---
% 避免除以0，只计算非零部分的功率
signal_part = rx_signal(rx_signal~=0);
if isempty(signal_part)
    P_signal = 0;
else
    P_signal = mean(abs(signal_part).^2); 
end

% --- B. 生成热噪声 (Thermal Noise) ---
% 修正：使用 N_actual 而不是 N_total
P_noise = P_signal / (10^(SNR_dB/10));
thermal_noise = sqrt(P_noise/2) * (randn(1, N_actual) + 1j*randn(1, N_actual));

% --- C. 生成压制干扰 (Barrage Noise Jamming) ---
% 修正：使用 N_actual 而不是 N_total
P_jamming = P_signal * (10^(JSR_dB/10));
jamming_signal = sqrt(P_jamming/2) * (randn(1, N_actual) + 1j*randn(1, N_actual));

% --- D. 合成总接收信号 ---
rx_clean_with_noise = rx_signal + thermal_noise;
rx_jammed = rx_signal + thermal_noise + jamming_signal;

%% 4. 信号处理：脉冲压缩 (Matched Filtering)
% 这是雷达检测的核心步骤，我们看看干扰后的效果
pc_clean = abs(xcorr(rx_clean_with_noise, s_tx));
pc_jammed = abs(xcorr(rx_jammed, s_tx));

% 截取有效部分并归一化
L = length(rx_signal);
pc_clean = pc_clean(L:end);
pc_clean = pc_clean / max(pc_clean); 

pc_jammed = pc_jammed(L:end);
pc_jammed = pc_jammed / max(pc_jammed); % 注意：这里通常归一化会发现目标峰值变很小

% 距离轴转换
r_axis = (0:length(pc_clean)-1) * c / (2 * fs);

%% 5. 绘图对比
figure('Position', [100, 100, 1000, 600]);

subplot(2,1,1);
plot(r_axis, 20*log10(pc_clean), 'LineWidth', 1.5);
grid on; hold on;
xline(R_target, '--r', 'Target');
title(['无干扰环境 (SNR=' num2str(SNR_dB) 'dB) - 脉冲压缩结果']);
xlabel('距离 (m)'); ylabel('幅度 (dB)');
legend('雷达回波', '真实位置');
ylim([-60, 5]);

subplot(2,1,2);
plot(r_axis, 20*log10(pc_jammed), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
grid on; hold on;
xline(R_target, '--r', 'Target');
title(['压制干扰环境 (JSR=' num2str(JSR_dB) 'dB) - 目标被淹没']);
xlabel('距离 (m)'); ylabel('幅度 (dB)');
legend('干扰后信号', '真实位置');
ylim([-60, 5]);

% 如果你在图中看不到红色的虚线位置有明显的尖峰，说明压制成功了。