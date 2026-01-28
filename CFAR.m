%% 2D CA-CFAR 检测与聚类算法演示 (Standalone Demo)
clear; close all; clc;

%% 1. 生成模拟数据 (模拟雷达 RD 图)
% ---------------------------------------------------------
Nr = 200;   % 距离维采样点数 (行)
Nd = 100;   % 多普勒维采样点数 (列)

% 生成高斯白噪声背景 (瑞利分布幅度 -> 指数分布功率)
% 假设噪声功率均值为 1 (0dB)
noise = abs(randn(Nr, Nd) + 1j*randn(Nr, Nd)).^2 / 2;
Map_Power = noise; 

% 注入人造目标 (Target)
% 目标格式: [行(Range), 列(Doppler), 信噪比(dB)]
targets = [
    50,  30,  15;   % 目标1: 孤立强目标
    120, 70,  13;   % 目标2: 孤立中等目标
    150, 45,  11;   % 目标3: 弱目标 (临界)
    150, 47,  10;   % 目标3的伴随点 (模拟扩展目标/临近目标)
    151, 46,  10;   % 目标3的伴随点
];

% 将目标叠加到噪声背景中
fprintf('正在生成模拟数据，包含 %d 个目标...\n', size(targets, 1));
for i = 1:size(targets, 1)
    r = targets(i, 1);
    c = targets(i, 2);
    snr_linear = 10^(targets(i, 3)/10);
    Map_Power(r, c) = Map_Power(r, c) + snr_linear; % 简单的幅度叠加
end

%% 2. CA-CFAR 参数设置
% ---------------------------------------------------------
Pfa = 1e-5;           % 虚警概率 (Probability of False Alarm)
Ref_Win = [8, 4];     % 参考窗口半径 [行, 列] -> 实际大小 (2R+1)x(2C+1)
Guard_Win = [2, 1];   % 保护窗口半径 [行, 列]

% 计算 CFAR 门限因子 Alpha (针对平方律检波/指数分布噪声)
% 核心思想: 噪声窗口内的点数越多，估计越准，需要的余量(Alpha)就越小
Lr = Ref_Win(1); Lc = Ref_Win(2);
Gr = Guard_Win(1); Gc = Guard_Win(2);
% 参考单元总数 = 外框总数 - 内框总数
N_ref = (2*Lr+1)*(2*Lc+1) - (2*Gr+1)*(2*Gc+1); 
alpha = N_ref * (Pfa^(-1/N_ref) - 1);

fprintf('CFAR 参数: Pfa=%.1e, 参考单元数=%d, 门限因子 Alpha=%.2f\n', ...
        Pfa, N_ref, alpha);

%% 3. 快速 CA-CFAR 处理 (卷积法)
% ---------------------------------------------------------
tic; % 开始计时

% 3.1 构建卷积核 (Kernel)
% 这是一个中心挖空的矩阵: 参考区=1/N_ref, 保护区和中心=0
Kernel = ones(2*Lr+1, 2*Lc+1);
Kernel(Lr-Gr+1 : Lr+Gr+1, Lc-Gc+1 : Lc+Gc+1) = 0;
Kernel = Kernel / N_ref; % 归一化，卷积结果直接就是平均噪声功率

% 3.2 卷积计算背景噪声图
% 'same' 保持输出尺寸与输入一致
Noise_Mean_Map = conv2(Map_Power, Kernel, 'same');

% 3.3 计算门限图
Threshold_Map = alpha .* Noise_Mean_Map;

% 3.4 判决 (生成二值化检测图)
Detection_Map = (Map_Power > Threshold_Map);

% 3.5 边缘剔除 (可选)
% 卷积在图像边缘会有误差，通常强制置0
Edge_Mask = zeros(Nr, Nd);
Edge_Mask(Lr+1:Nr-Lr, Lc+1:Nd-Lc) = 1;
Detection_Map = Detection_Map & Edge_Mask;

cfar_time = toc;
fprintf('CFAR 处理完成，耗时: %.4f 秒\n', cfar_time);

%% 4. 聚类处理 (Clustering)
% ---------------------------------------------------------
% 目的: 将属于同一个目标的多个邻近检测点合并为一个
tic;

% 使用 MATLAB 内置函数寻找连通域 (8连通)
CC = bwconncomp(Detection_Map, 8);
Num_Targets = CC.NumObjects;

Final_Targets = []; % 存储结果: [行, 列, 峰值强度dB]

for i = 1:Num_Targets
    % 获取第 i 个目标的所有像素索引
    pixels = CC.PixelIdxList{i};
    
    % --- 核心步骤: 峰值提取 ---
    % 在这一团点中，找到能量最大的那个点作为目标的中心
    [max_val, idx_in_cluster] = max(Map_Power(pixels));
    global_idx = pixels(idx_in_cluster);
    
    % 将线性索引转回 (行, 列) 坐标
    [r_peak, c_peak] = ind2sub(size(Map_Power), global_idx);
    
    % 记录结果
    Final_Targets = [Final_Targets; r_peak, c_peak, 10*log10(max_val)];
end

cluster_time = toc;
fprintf('聚类处理完成，耗时: %.4f 秒。检测到 %d 个目标。\n', cluster_time, Num_Targets);

%% 5. 结果可视化
% ---------------------------------------------------------
figure('Name', 'CFAR & Clustering Demo', 'Position', [100, 100, 1200, 400], 'Color', 'w');

% 图1: 原始功率图 (dB)
subplot(1, 3, 1);
imagesc(10*log10(Map_Power)); 
colormap('jet'); colorbar; axis xy;
title('1. 原始信号 (dB)');
xlabel('Doppler'); ylabel('Range');

% 图2: CFAR 二值化检测结果
subplot(1, 3, 2);
imagesc(Detection_Map);
colormap(gca, 'gray'); % 黑白图
axis xy;
title(['2. CFAR 检测结果 (Pfa=' num2str(Pfa) ')']);
xlabel('Doppler'); ylabel('Range');

% 图3: 聚类后的最终结果
subplot(1, 3, 3);
imagesc(10*log10(Map_Power)); %以此为底图
colormap('jet'); axis xy; hold on;
title(['3. 聚类后输出 (' num2str(Num_Targets) ' Targets)']);
xlabel('Doppler'); ylabel('Range');

% 在图3上画出最终目标的红圈
if ~isempty(Final_Targets)
    plot(Final_Targets(:,2), Final_Targets(:,1), 'ro', ...
         'MarkerSize', 10, 'LineWidth', 2);
    plot(Final_Targets(:,2), Final_Targets(:,1), 'w+', ...
         'MarkerSize', 6, 'LineWidth', 1);
     
    % 打印检测列表
    fprintf('--------------------------------------\n');
    fprintf('ID\t 行(R)\t 列(D)\t 强度(dB)\n');
    for k = 1:size(Final_Targets, 1)
        fprintf('%d\t %d\t %d\t %.2f\n', ...
            k, Final_Targets(k,1), Final_Targets(k,2), Final_Targets(k,3));
        % 图上标字
        text(Final_Targets(k,2)+2, Final_Targets(k,1), ...
             sprintf('%.1fdB', Final_Targets(k,3)), 'Color', 'w', 'FontWeight', 'bold');
    end
    fprintf('--------------------------------------\n');
else
    fprintf('未检测到目标。\n');
end