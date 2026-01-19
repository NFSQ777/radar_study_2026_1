% ===================== MTD（动目标检测）仿真代码（3D可视化版） =====================
% 作者：编程导师
% 功能：MTI+多普勒FFT实现MTD，生成距离-多普勒谱并绘制3D立体图
clc; clear; close all;

%% 1. 定义仿真参数（与MTI一致，新增多目标速度）
Nr = 256;           % 距离单元数
Np = 64;            % 脉冲数（MTD需更多脉冲保证多普勒分辨率）
Fs = 1e6;           % 采样率(Hz)
Tc = 100e-6;        % 脉冲重复周期(PRT, s)
f0 = 10e9;          % 雷达工作频率(10GHz)
c = 3e8;            % 光速(m/s)
lambda = c/f0;      % 波长(m)
v1 = 50;            % 目标1速度(m/s)
v2 = -30;           % 目标2速度(m/s)，负速度表示远离雷达

%% 2. 生成仿真雷达回波数据（杂波+双动目标+噪声）
echo_data = zeros(Nr, Np);

% (1) 固定杂波
clutter_amp = 5 + 2*randn(Nr, 1);
echo_data = echo_data + clutter_amp * ones(1, Np);

% (2) 加入两个动目标
target1_range_idx = 100;  % 目标1：第100个距离单元
target2_range_idx = 180;  % 目标2：第180个距离单元
fd1 = 2*v1/lambda;        % 目标1多普勒频率
fd2 = 2*v2/lambda;        % 目标2多普勒频率
for n = 1:Np
    % 目标1回波
    echo_data(target1_range_idx, n) = echo_data(target1_range_idx, n) + ...
        10 * exp(1j * 2*pi*fd1*n*Tc);
    % 目标2回波
    echo_data(target2_range_idx, n) = echo_data(target2_range_idx, n) + ...
        8 * exp(1j * 2*pi*fd2*n*Tc);
end

% (3) 加噪声
echo_data = echo_data + (0.5 + 0.5j) * randn(Nr, Np);

%% 3. 第一步：MTI一阶对消（抑制杂波）
mti_output = zeros(Nr, Np);
for n = 2:Np
    mti_output(:, n) = echo_data(:, n) - echo_data(:, n-1);
end
% 去除第1个脉冲（对消后无数据）
mti_output = mti_output(:, 2:end);
Np_mti = Np - 1;  % 对消后剩余脉冲数

%% 4. 第二步：MTD核心——慢时间维FFT（多普勒滤波）
% 对每个距离单元的慢时间序列做FFT
mtd_output = fftshift(fft(mti_output, [], 2), 2);  % FFT+频谱移位（零频居中）
mtd_amp = abs(mtd_output);  % 取幅度谱

%% 5. 多普勒频率/速度轴计算
fd_axis = linspace(-1/(2*Tc), 1/(2*Tc), Np_mti);  % 多普勒频率轴（Nyquist范围）
v_axis = fd_axis * lambda / 2;                    % 转换为速度轴

%% 6. 结果可视化：3D距离-多普勒谱（核心修改部分）
% 先归一化并转换为dB（增强对比度，压制噪声）
mtd_amp_dB = 20*log10(mtd_amp/max(mtd_amp(:)));
% 限制dB范围（避免噪声的负dB值过低，影响可视化）
mtd_amp_dB(mtd_amp_dB < -40) = -40;  

% 创建网格矩阵（3D绘图需要X/Y轴的网格）
[V, R] = meshgrid(v_axis, 1:Nr);  

figure('Color','w','Position',[100,100,1000,600]);
% 绘制3D曲面图（surf）+ 网格线（grid on）
surf(V, R, mtd_amp_dB, 'EdgeColor','none', 'FaceAlpha',0.8);  % EdgeColor取消边缘线，FaceAlpha设置透明度
colormap(jet);  % 配色方案
colorbar;       % 颜色刻度条
xlabel('速度 (m/s)','FontSize',12);
ylabel('距离单元','FontSize',12);
zlabel('回波幅度 (dB)','FontSize',12);
title('MTD距离-多普勒谱 3D立体图','FontSize',14);
% 调整视角（方位角45°，仰角30°，更易观察目标峰值）
view(45, 30);  
grid on;
shading interp;  % 插值着色，让曲面更平滑

% 标记两个目标的3D位置（在峰值处绘制红色/白色标记）
hold on;
% 目标1：速度v1，距离单元target1_range_idx，幅度取该位置的dB值
plot3(v1, target1_range_idx, mtd_amp_dB(target1_range_idx, find(abs(v_axis-v1)==min(abs(v_axis-v1)))), ...
    'ro', 'MarkerSize',10, 'MarkerFaceColor','r', 'DisplayName',['目标1: ', num2str(v1), 'm/s']);
% 目标2：速度v2，距离单元target2_range_idx
plot3(v2, target2_range_idx, mtd_amp_dB(target2_range_idx, find(abs(v_axis-v2)==min(abs(v_axis-v2)))), ...
    'ws', 'MarkerSize',10, 'MarkerFaceColor','w', 'DisplayName',['目标2: ', num2str(v2), 'm/s']);
legend('Location','best','FontSize',10);

% 可选：同时绘制2D俯视图（方便对比）
figure('Color','w');
imagesc(v_axis, 1:Nr, mtd_amp_dB);
colormap(jet); colorbar;
xlabel('速度 (m/s)'); ylabel('距离单元');
title('MTD距离-多普勒谱 2D俯视图');
hold on;
plot(v1, target1_range_idx, 'ro', 'MarkerSize',8);
plot(v2, target2_range_idx, 'ws', 'MarkerSize',8);
grid on;