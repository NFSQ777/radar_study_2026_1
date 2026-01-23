% =========================================================================
% 相控阵雷达波束扫描仿真 (Phased Array Beam Steering Simulation)
% 模拟一个均匀线阵 (ULA) 的波束形成与扫描过程
% =========================================================================
clear; clc; close all;

%% 1. 雷达系统参数设置
c = 3e8;                % 光速 (m/s)
fc = 10e9;              % 载波频率 10GHz (X波段)
lambda = c / fc;        % 波长
N = 16;                 % 天线单元数量 (阵元数)
d = lambda / 2;         % 阵元间距 (通常取半波长以避免栅瓣)

% 扫描角度范围 (例如从 -60度 扫到 +60度)
scan_angles = -60:2:60; 

%% 2. 准备绘图
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
theta = -90:0.1:90;        %观察角度范围 (-90度 到 90度)
theta_rad = deg2rad(theta);%角度转弧度，每0.1度对应多少弧度

% 预先计算观察方向的导向矢量系数 u = (2*pi*d/lambda) * sin(theta)
k = 2 * pi / lambda;%波数，相位随距离的变化率
u_obs = k * d * sin(theta_rad); 
n_vec = (0:N-1).';      % 阵元索引向量 [0, 1, ..., N-1]'

%% 3. 循环仿真：波束扫描
for target_angle = scan_angles
    
    % --- 核心算法：计算波束指向所需的相位 ---
    % 想要波束指向 target_angle，相邻阵元需要的相位差为：
    % delta_phi = -k * d * sin(theta0)
    target_angle_rad = deg2rad(target_angle);
    phase_shift = -k * d * sin(target_angle_rad);
    
    % 生成加权向量 (Weight Vector) w = exp(j * n * phase_shift)
    % 这里只加了相位权值，幅度权值设为1 (均匀加权)
    w = exp(1j * n_vec * phase_shift); 
    
    % --- 计算阵列方向图 (Array Pattern) ---
    % y = w' * a(theta)
    % 等效公式：AF = sum( w_n * exp(j * k * d * n * sin(theta)) )
    % 利用矩阵运算加速计算
    steering_matrix = exp(1j * n_vec * u_obs);
    AF = abs(w.' * steering_matrix);  % AF = Array Factor
    
    % --- 归一化并转为 dB ---
    AF = AF / max(AF);          % 归一化到 1
    AF_dB = 20 * log10(AF + 1e-6); % 转为 dB，加微小量防止log(0)
    
    %% 4. 绘图显示
    subplot(2,1,1);
    plot(theta, AF_dB, 'LineWidth', 2, 'Color', 'b');
    grid on;
    xline(target_angle, 'r--', 'LineWidth', 1.5); % 标记当前指向
    xlim([-90 90]);
    ylim([-50 0]);
    xlabel('角度 (Degree)');
    ylabel('归一化增益 (dB)');
    title(['直角坐标系波束方向图 (当前指向: ' num2str(target_angle) '°)']);
    legend('波束方向图', '波束指向');
    
    subplot(2,1,2);
    % 极坐标绘图有点麻烦，这里简单手动画一个
    polarplot(theta_rad, AF, 'LineWidth', 2, 'Color', 'r');
    rlim([0 1]);
    title('极坐标系波束方向图 (线性幅度)');
    
    drawnow; % 刷新绘图，产生动画效果
    % pause(0.05); % 如果电脑跑得太快，取消注释这行可以变慢点
end

% 仿真结束提示
disp('波束扫描仿真完成！');
