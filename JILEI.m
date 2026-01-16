clear; clc; close all;

%% 1. 雷达参数定义 (System Parameters)
fc = 77e9;              % 载波频率 77GHz
c = 3e8;                % 光速
lambda = c/fc;          % 波长 3.9mm

B = 200e6;              % 带宽 200MHz (决定距离分辨率)
Tc = 50e-6;             % 脉冲重复时间 (Chirp持续时间) 50us
Slope = B/Tc;           % 调频斜率
Fs = 10e6;              % 采样率 10Msps

Nr = 512;               % Fast Time: 每个Chirp的采样点数 (距离维)
Nd = 128;               % Slow Time: 发射Chirp的数量 (用于相参积累)

t_fast = (0:Nr-1)/Fs;   % 快时间轴 (针对单个Chirp)
t_slow = (0:Nd-1)*Tc;   % 慢时间轴 (针对整个帧)

%% 2. 目标设定 (Target Setup)
% 设定一个难以检测的目标：距离远，且淹没在噪声中
target_dist = 80;       % 目标距离 80米
target_vel = 15;        % 目标速度 15m/s (远离)
SNR_dB = -10;           % 信噪比 -10dB (信号比噪声还弱)

%% 3. 信号生成 (Signal Generation)
% 我们直接生成经过混频器(Mixer)后的 中频信号 (IF Signal)
% 物理公式：IF = exp(j * (2*pi*f_beat*t + phase_shift))

rx_data = zeros(Nd, Nr); % 初始化数据矩阵 [128 x 512]

for i = 1:Nd
    % 当前 Chirp 时刻目标的真实距离 (考虑运动)
    % 慢时间 i 对应的距离变化体现了多普勒效应
    r_current = target_dist + target_vel * t_slow(i);
    
    % 往返时延 tau
    tau = 2 * r_current / c;
    
    % 差频频率 (Beat Frequency) -> 对应距离
    f_beat = Slope * tau;
    
    % 相位偏移 (Phase Shift) -> 对应速度 (多普勒)
    % 核心：4*pi*R/lambda 这一项包含了微小的距离变化引起的巨大相位旋转
    phase_shift = 4 * pi * r_current / lambda; 
    
    % 生成纯净的中频信号 (复数信号)
    sig_pure = exp(1j * (2 * pi * f_beat * t_fast + phase_shift));
    
    % 添加高斯白噪声 (AWGN)
    sig_noisy = awgn(sig_pure, SNR_dB, 'measured');
    
    % 填入数据矩阵
    rx_data(i, :) = sig_noisy;
end

%% 4. 可视化：原始信号
figure('Position', [100, 100, 1200, 800]);
subplot(2,2,1);
plot(t_fast*1e6, real(rx_data(1,:)));
title(['原始时域信号 (第1个Chirp, SNR=' num2str(SNR_dB) 'dB)']);
xlabel('时间 (us)'); ylabel('幅度');
grid on;
% 说明：根本看不出有正弦波，全是噪声

%% 5. 第一步：脉冲压缩 (Range FFT)
% 原理：将时域上的长脉冲能量，在频率上聚焦成一个峰
range_win = hamming(Nr).'; % 加窗抑制旁瓣
range_profile_matrix = zeros(Nd, Nr);

for i = 1:Nd
    % 对每一行(每个Chirp)做 FFT
    range_profile_matrix(i, :) = fft(rx_data(i, :) .* range_win); 
end

% 转换距离坐标轴
range_axis = (0:Nr-1) * Fs * c / (2 * Slope) / Nr;

subplot(2,2,2);
plot(range_axis, db(abs(range_profile_matrix(1, :))));
title('脉冲压缩后 (Range FFT)');
xlabel('距离 (m)'); ylabel('幅度 (dB)');
xlim([0 150]);
grid on;
% 说明：此时可能隐约能看到一个峰，但因为噪声大，可能还是不明显

%% 6. 第二步：相参积累 (Doppler FFT)
% 原理：利用相位的一致性，把128个Chirp的能量叠加起来
doppler_win = hamming(Nd); % 加窗
final_map = zeros(Nd, Nr);

% 对每一列(每个距离点)在慢时间维度做 FFT
for j = 1:Nr
    final_map(:, j) = fftshift(fft(range_profile_matrix(:, j) .* doppler_win)); 
end

% 转换速度坐标轴
vel_axis = (-Nd/2 : Nd/2-1) * (lambda / (2 * Tc * Nd));

%% 7. 最终结果：距离-多普勒图 (Range-Doppler Map)
subplot(2,2,[3,4]);
% 使用 mesh 绘制 3D 效果
samples_to_show_r = find(range_axis > 150, 1); % 只画前150米
mesh(range_axis(1:samples_to_show_r), vel_axis, db(abs(final_map(:, 1:samples_to_show_r))));
view(0, 90); % 俯视图
colorbar;
title('相参积累后 (Range-Doppler Map 2D-FFT)');
xlabel('距离 (m)'); ylabel('速度 (m/s)');
zlabel('幅度 (dB)');
axis tight;

%% 结果验证输出
[max_val, max_idx] = max(abs(final_map(:)));
[v_idx, r_idx] = ind2sub(size(final_map), max_idx);

calc_dist = range_axis(r_idx);
calc_vel = vel_axis(v_idx);

fprintf('--- 仿真结果 ---\n');
fprintf('真实目标: 距离 %.2fm, 速度 %.2fm/s\n', target_dist, target_vel);
fprintf('雷达解算: 距离 %.2fm, 速度 %.2fm/s\n', calc_dist, calc_vel);
fprintf('相参积累增益: 信号从噪声中“浮现”了出来。\n');
