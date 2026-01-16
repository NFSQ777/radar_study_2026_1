clc; clear; close all;   

%% 1. 生成模拟的距离-多普勒数据 (Range-Doppler Map)
% 这里的尺寸是 200x100，背景是高斯白噪声
Nr = 200; % 距离单元数
Nd = 100; % 多普勒单元数
RDM = abs(randn(Nr, Nd) + 1i*randn(Nr, Nd)).^2; % 能量矩阵(平方律检波)

% 添加几个模拟目标 (信噪比 SNR 不同)
% 目标1: 强目标
RDM(80, 50) = 50; 
% 目标2: 弱目标
RDM(120, 30) = 15;
% 目标3: 临近目标 (测试保护单元)
RDM(80, 53) = 40; 

%% 2. 2D CA-CFAR 参数设置
Pfa = 1e-6;      % 虚警概率 (False Alarm Probability)
Tr = 8;          % 距离维训练单元数 (单侧) Training Cells
Td = 4;          % 多普勒维训练单元数 (单侧)
Gr = 2;          % 距离维保护单元数 (单侧) Guard Cells
Gd = 2;          % 多普勒维保护单元数 (单侧)

%% 3. 构建 CFAR 卷积掩码 (Kernel)
% 掩码就像一个“甜甜圈”：
% 中间是 0 (忽略 CUT 和 保护单元)
% 外圈是 1 (训练单元，用来计算平均值)

% 完整窗口的大小
winSizeR = 2*Tr + 2*Gr + 1;
winSizeD = 2*Td + 2*Gd + 1;

% 初始化全1矩阵
mask = ones(winSizeR, winSizeD); 

% 计算保护区域的索引范围
guardStartR = Tr + 1;
guardEndR   = Tr + 2*Gr + 1;
guardStartD = Td + 1;
guardEndD   = Td + 2*Gd + 1;

% 将包含 CUT 和保护单元的中心区域置为 0
mask(guardStartR:guardEndR, guardStartD:guardEndD) = 0;

% 归一化掩码（这点很重要，相当于直接求平均值）
numTrainCells = sum(mask(:)); % 训练单元的总个数
mask = mask / numTrainCells; 

%% 4. 执行 CA-CFAR 检测
% 计算门限因子 Alpha
% 对于平方律检波 CA-CFAR，Alpha 计算公式如下：
alpha = numTrainCells * (Pfa^(-1/numTrainCells) - 1); 

% 1. 计算每一个点周围的噪声平均功率 (Noise Floor)
% 使用二维卷积快速实现滑窗均值计算
noiseMap = conv2(RDM, mask, 'same'); 

% 2. 计算最终的检测门限
thresholdMap = noiseMap * alpha;

% 3. 比较：信号 > 门限 ?
detectionMap = (RDM > thresholdMap);

% 4. 去除边缘效应 (因为卷积在边缘处计算不准)
edgeR = Tr + Gr;
edgeD = Td + Gd;
detectionMap(1:edgeR, :) = 0; detectionMap(end-edgeR:end, :) = 0;
detectionMap(:, 1:edgeD) = 0; detectionMap(:, end-edgeD:end) = 0;

%% 5. 结果提取与绘图
[detR, detD] = find(detectionMap); % 提取检测到的坐标

% --- 绘图 ---
figure('Position', [100, 100, 1000, 400]);

% 原始 RDM 图
subplot(1, 2, 1);
imagesc(10*log10(RDM)); 
axis xy; colormap('jet'); colorbar;
title('原始距离-多普勒图 (Log scale)');
xlabel('Doppler'); ylabel('Range');
clim([0 20]); % 限制颜色范围以便看清噪声

% CFAR 检测结果图
subplot(1, 2, 2);
imagesc(detectionMap); 
axis xy; colormap('gray'); 
title(['CA-CFAR 检测结果 (Pfa=', num2str(Pfa), ')']);
xlabel('Doppler'); ylabel('Range');

% 在结果图上标出中心点（可选）
hold on;
plot(detD, detR, 'rs', 'MarkerSize', 8, 'LineWidth', 1.5);
legend('检测点');

fprintf('检测到目标点数量: %d\n', length(detR));
