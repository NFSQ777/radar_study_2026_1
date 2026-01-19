% ===================== MTI（动目标指示）仿真代码 =====================
% 功能：实现一阶MTI延迟线对消，抑制固定杂波，提取动目标
clc; clear; close all;

%% 1. 定义仿真参数
Nr = 256;           % 距离单元数（快时间维）
Np = 32;            % 脉冲数（慢时间维，即MTI处理的脉冲数）
Fs = 1e6;           % 采样率(Hz)
Tc = 100e-6;        % 脉冲重复周期(PRT, s)
f0 = 10e9;          % 雷达工作频率(10GHz，X波段)
c = 3e8;            % 光速(m/s)
lambda = c/f0;      % 波长(m)
v_target = 50;      % 动目标速度(m/s)，正速度表示靠近雷达

%% 2. 生成仿真雷达回波数据（杂波+动目标+噪声）
% 初始化回波矩阵：维度为 距离单元数 × 脉冲数
echo_data = zeros(Nr, Np);

% (1) 生成固定杂波（覆盖所有距离单元，幅度随机）
clutter_amp = 5 + 2*randn(Nr, 1);  % 杂波幅度（含小幅随机波动）
echo_data = echo_data + clutter_amp * ones(1, Np);  % 杂波在所有脉冲中不变

% (2) 加入动目标（位置：第100个距离单元，所有脉冲都有）
target_range_idx = 100;            % 目标距离单元索引
doppler_freq = 2*v_target/lambda;  % 目标多普勒频率
for n = 1:Np
    % 动目标回波：幅度+多普勒相位调制
    echo_data(target_range_idx, n) = echo_data(target_range_idx, n) + ...
        10 * exp(1j * 2*pi*doppler_freq*n*Tc);
end

% (3) 加入高斯白噪声
echo_data = echo_data + (0.5 + 0.5j) * randn(Nr, Np);

%% 3. 一阶MTI延迟线对消（核心处理）
% 原理：y(n) = x(n) - x(n-1)，相邻脉冲相减
mti_output = zeros(Nr, Np);
for n = 2:Np  % 从第2个脉冲开始对消
    mti_output(:, n) = echo_data(:, n) - echo_data(:, n-1);
end

%% 4. 结果可视化
figure('Color','w');
subplot(2,1,1);
% 绘制原始回波的幅度（第1个脉冲）
plot(1:Nr, abs(echo_data(:,1)), 'b-', 'LineWidth',1);
hold on;
plot(target_range_idx, abs(echo_data(target_range_idx,1)), 'ro', 'MarkerSize',6);
xlabel('距离单元'); ylabel('回波幅度');
title('原始回波（第1个脉冲）：杂波掩盖动目标');
legend('杂波+噪声','动目标','Location','best');
grid on;

subplot(2,1,2);
% 绘制MTI处理后的幅度（第2个脉冲，对消后）
plot(1:Nr, abs(mti_output(:,2)), 'g-', 'LineWidth',1);
hold on;
plot(target_range_idx, abs(mti_output(target_range_idx,2)), 'ro', 'MarkerSize',6);
xlabel('距离单元'); ylabel('MTI输出幅度');
title('MTI一阶对消后：杂波被抑制，动目标凸显');
legend('MTI输出','动目标','Location','best');
grid on;