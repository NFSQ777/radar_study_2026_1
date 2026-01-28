%% LFM雷达全流程仿真：3D可视化版本
% 流程：信号生成 -> 回波模拟 -> 脉冲压缩 -> MTD -> 3D显示
clear; close all; clc;

%% 1. 雷达与目标参数设置
c = 3e8;                % 光速 (m/s)
fc = 3e9;               % 载波频率 (Hz) -> S波段
lambda = c / fc;        % 波长

% 雷达发射参数
Tp = 10e-6;             % 脉冲宽度 (s)
B = 30e6;               % 信号带宽 (Hz)
K = B / Tp;             % 调频斜率
PRF = 1000;             % 脉冲重复频率 (Hz) (降低一点以便观察)
PRI = 1 / PRF;          % 脉冲重复间隔 (s)
fs = 2 * B;             % 快时间采样率
dt = 1 / fs;            % 采样间隔

% 仿真配置
NumPulses = 64;         % 积累的脉冲数 (慢时间维点数)
Gate_len = round(PRI/dt); % 一个PRI内的采样点数 

% 目标参数 (可随意修改)
target_R0 = 35000;       % 目标距离 (m)
target_V = 20;          % 目标速度 (m/s) [正值: 靠近雷达]

%% 2. 发射信号生成 (LFM)
t = (-Tp/2 : dt : Tp/2 - dt); 
St = exp(1j * pi * K * t.^2); % LFM信号

% --- 绘图1：发射信号 ---
figure('Name','发射信号','Position',[50, 50, 600, 400]);
subplot(2,1,1); 
plot(t*1e6, real(St),'b','DisplayName','I路（同相分量）'); 
title('发射信号'); xlabel('时间 (\mus)'); 
hold on;
plot(t*1e6, imag(St), 'r', 'DisplayName','Q路（正交分量）');
grid on;
subplot(2,1,2); plot(linspace(-B/2, B/2, length(St))/1e6, abs(fftshift(fft(St))));
title('发射信号频谱'); xlabel('频率 (MHz)'); grid on;

%% 3. 回波信号模拟 (加入时延、多普勒和噪声)
Sr_matrix = zeros(Gate_len, NumPulses);
FastTimeAxis = (0:Gate_len-1) * dt;

for m = 1 : NumPulses
    tm = (m - 1) * PRI;             % 慢时间
    CurrR = target_R0 - target_V * tm; % 当前距离
    tau = 2 * CurrR / c;            % 当前时延
    
    % 计算采样点偏移
    n_tau = round(tau / dt);
    
    if (n_tau + length(t) <= Gate_len)
        % 核心公式：S_rx = S_tx(t-tau) * exp(-j*2*pi*fc*tau)
        phase_shift = exp(-1j * 2 * pi * fc * tau); 
        range_indices = n_tau + (1:length(t));
        Sr_matrix(range_indices, m) = St * phase_shift;
    end
end

% 增加噪声 (SNR = -20dB, 噪声稍微大一点更能看出处理增益)
Sr_matrix_noisy = awgn(Sr_matrix, -20, 'measured'); 

% --- 绘图2：时域回波 ---
figure('Name','时域回波','Position',[100, 100, 600, 300]);
plot(FastTimeAxis*1e6, real(Sr_matrix_noisy(:,1)));
title('第1个脉冲的时域回波 (含噪声)'); xlabel('快时间 (\mus)'); grid on;
%xlim([0, 50]); % 只看前50us

%% 4. 脉冲压缩 (Range Compression) - 快时间维处理
N_fft_rg = 2^nextpow2(Gate_len); 
H_ref = fft(St, N_fft_rg);
H_ref_conj = conj(H_ref);

Sr_RangeComp = zeros(N_fft_rg, NumPulses);
for m = 1 : NumPulses
    Sr_RangeComp(:,m) = ifft(fft(Sr_matrix_noisy(:,m), N_fft_rg) .* H_ref_conj.');
end
RangeAxis = (0:N_fft_rg-1) * c / (2 * fs);  

% --- 绘图3：脉冲压缩后 ---
figure('Name','脉冲压缩后','Position',[150, 150, 600, 300]);
plot(RangeAxis, abs(Sr_RangeComp(:,1)));
title('脉冲压缩后 (距离像)'); xlabel('距离 (m)'); grid on;
xlim([0, max(RangeAxis)/2]); % 缩放一下X轴

%% 5. MTD / 相干积累 - 慢时间维处理
% 加窗抑制旁瓣
win = hamming(NumPulses).';
Sr_RangeComp_Win = Sr_RangeComp .* win; 

% 慢时间FFT
padding_doppler = 256; 
RDM = fftshift(fft(Sr_RangeComp_Win, padding_doppler, 2), 2);

% 计算坐标轴
freq_doppler = linspace(-PRF/2, PRF/2, padding_doppler);
VelocityAxis = freq_doppler * lambda / 2; % 速度轴

%% 6. 3D 可视化 (关键修改部分)
figure('Name','3D Range-Doppler Map','Position',[200, 100, 1000, 700], 'Color', 'w');

% 为了显示效果更好，取部分距离范围并转换为dB
[X, Y] = meshgrid(VelocityAxis, RangeAxis);
RDM_dB = 20*log10(abs(RDM) + 1e-6); % 转换为dB，加小量防止log(0)

% 限制显示范围 (裁剪掉太远的距离，只看目标附近，让图更清楚)
plot_range_idx = range_indices(1) - 500 : range_indices(end) + 2000; 
% 如果上面的索引超出范围，可以使用下面的全范围：
% plot_range_idx = 1:length(RangeAxis)/2; 

% *** 绘制 3D 曲面图 ***
% 使用 surf 绘制实体曲面
s = surf(X(plot_range_idx, :), Y(plot_range_idx, :), RDM_dB(plot_range_idx, :));

% 美化设置
s.EdgeColor = 'none';       % 去掉网格线，使曲面光滑
colormap jet;               % 彩虹色图
colorbar;                   % 显示颜色条
ylabel(colorbar, '幅度 (dB)', 'FontSize', 12);

% 坐标轴标签
xlabel('速度 (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('距离 (m)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('幅度 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title(['LFM雷达 3D 距离-多普勒谱 (目标: ' num2str(target_R0) 'm, ' num2str(target_V) 'm/s)'], 'FontSize', 14);

% 视角与光照 (营造立体感)
axis tight; 
view(-45, 30);              %设置 3D 视角 (方位角-45度, 仰角30度)
camlight left;              % 添加左侧补光
lighting gouraud;           % 使用高洛德着色算法，让曲面看起来更顺滑
grid on;
set(gca, 'FontSize', 10);

% 标记最高点
[max_val, max_idx] = max(RDM_dB(:));
[r_peak, v_peak] = ind2sub(size(RDM_dB), max_idx);
hold on;
plot3(VelocityAxis(v_peak), RangeAxis(r_peak), max_val, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(VelocityAxis(v_peak), RangeAxis(r_peak), max_val + 2, '  Target', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
hold off;

%% 7. 打印结果
Est_Range = RangeAxis(r_peak);
Est_Velocity = VelocityAxis(v_peak);
fprintf('----------------------------------\n');
fprintf('真实参数 -> 距离: %.2f m, 速度: %.2f m/s\n', target_R0, target_V);
fprintf('检测结果 -> 距离: %.2f m, 速度: %.2f m/s\n', Est_Range, Est_Velocity);
fprintf('----------------------------------\n');


