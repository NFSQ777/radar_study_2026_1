%% FMCW 雷达发射信号(Chirp)仿真生成
clear; clc; close all;

%% 1. 核心参数配置 (Parameter Configuration)
% 为了可视化效果，我们将频率设定在 MHz 级别，而不是真实的 GHz 级别
c = 3e8;                    % 光速 (m/s)

% --- 关键雷达参数 ---
T_chirp = 50e-6;            % 1. 脉冲持续时间 (50us) -> 决定了探测的能量积累时间
B = 100e6;                  % 2. 扫频带宽 (100MHz) -> 决定了距离分辨率 (Resolution = c/2B)
f_start = 20e6;             % 3. 起始频率 (20MHz) -> 类似于雷达的载波频率 f_c
K = B / T_chirp;            % 4. 频率斜率 (Slope) -> 极其重要的参数，决定了差频与距离的关系

% --- 仿真参数 ---
Fs = 4 * (f_start + B);     % 采样率：必须大于最大频率的2倍 (Nyquist定理)
Ts = 1 / Fs;                % 采样时间间隔
N = round(T_chirp / Ts);    % 总采样点数
t = (0:N-1) * Ts;           % 时间轴

%% 2. 信号生成 (Signal Generation)
% FMCW 信号的数学定义：
% 瞬时频率: f(t) = f_start + K * t
% 瞬时相位: phi(t) = 2 * pi * integral(f(t)) = 2 * pi * (f_start*t + 0.5*K*t^2)
% 发射信号: S_tx(t) = cos(phi(t))

phase = 2 * pi * (f_start * t + 0.5 * K * t.^2); % 相位随时间呈二次方变化
S_tx = cos(phase);                               % 生成时域实信号

%% 3. 时域可视化 (Time Domain Plot)
figure('Name', 'FMCW 发射信号分析', 'Position', [100, 100, 1000, 600]);

subplot(2,2,1);
plot(t*1e6, S_tx);
title('1. 完整发射信号 (时域)');
xlabel('时间 (us)'); ylabel('幅度 (V)');
grid on; axis tight;
% 解释：你看不到细节，因为频率太高，挤在一起像个实心矩形

subplot(2,2,2);
% 放大展示前 2us 和后 2us 的波形对比
zoom_points = 200; 
plot(t(1:zoom_points)*1e6, S_tx(1:zoom_points), 'b'); hold on;
plot(t(end-zoom_points:end)*1e6, S_tx(end-zoom_points:end), 'r');
title('2. 信号首尾放大对比');
xlabel('相对采样点时间'); ylabel('幅度');
legend(['起始频率: ' num2str(f_start/1e6) 'MHz'], ['终止频率: ' num2str((f_start+B)/1e6) 'MHz']);
grid on;
% 解释：你会发现红线（由于频率高）比蓝线（频率低）震荡得更密集

%% 4. 频域/时频域可视化 (Frequency Domain / Spectrogram)
% 我们不需要看单纯的 FFT，FMCW 最重要的是“频率随时间变化”的特性
% 使用短时傅里叶变换 (STFT) 绘制声谱图

subplot(2,2,3);
spectrogram(S_tx, 256, 250, 256, Fs, 'yaxis');
title('3. 时频图 (Spectrogram)');
% 解释：这是一条斜向上的直线。横轴是时间，纵轴是频率。
% 这条线的斜率就是 K (Slope)！

subplot(2,2,4);
% 绘制瞬时频率的理论曲线
f_instant = f_start + K * t;
plot(t*1e6, f_instant/1e6, 'k', 'LineWidth', 2);
title(['4. 理论频率-时间关系 (Slope = ' num2str(K/1e12) ' MHz/us)']);
xlabel('时间 (us)'); ylabel('频率 (MHz)');
grid on; grid minor;

%% 5. 参数计算验证
fprintf('------------- 参数验证 -------------\n');
fprintf('扫频带宽 (B)     : %.2f MHz\n', B/1e6);
fprintf('扫频时间 (Tc)    : %.2f us\n', T_chirp*1e6);
fprintf('频率斜率 (Slope) : %.3e Hz/s\n', K);
fprintf('理论距离分辨率   : %.2f 米 (c/2B)\n', c/(2*B));
fprintf('------------------------------------\n');
