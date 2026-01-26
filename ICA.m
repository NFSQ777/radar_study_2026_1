clc; clear; close all;

%% 1. 生成源信号 (Ground Truth)
fs = 1000;              % 采样率
T = 2;                  % 时间长度
t = 0:1/fs:T-1/fs;
N = length(t);

% 源信号 1: 模拟雷达脉冲 (目标) - 稀疏、非高斯性强
s1 = zeros(1, N);
s1(200:250) = 1;        % 脉冲 1
s1(800:850) = 1;        % 脉冲 2
s1(1500:1550) = 1;      % 脉冲 3
% 加一点载波纹理
s1 = s1 .* sin(2*pi*50*t); 

% 源信号 2: 压制干扰 (噪声/锯齿波)
% 干扰通常是弥散的，或者具有不同的统计分布
s2 = sawtooth(2*pi*5*t); % 锯齿波代表一种持续的干扰模式
% s2 = randn(1, N);      % 如果是纯高斯白噪声，ICA其实很难分，因为高斯性太强

S = [s1; s2];            % 2 x N 矩阵

%% 2. 混合过程 (Mixing)
% 模拟两个接收通道 (天线)
A_true = [0.8, 0.3;      % 混合矩阵
          0.4, 0.9];
      
X = A_true * S;          % 观测到的混合信号

% 添加少量接收机热噪声
X = X + 0.01 * randn(size(X));

%% 3. FastICA 算法实现 (核心部分)

% --- A. 去均值 (Centering) ---
X_mean = mean(X, 2);
X_centered = X - X_mean;

% --- B. 白化 (Whitening) ---
% 计算协方差矩阵
CovM = cov(X_centered');
% 特征值分解
[E, D] = eig(CovM);
% 白化矩阵 Q = D^(-1/2) * E'
Q = inv(sqrt(D)) * E';
Z = Q * X_centered;      % Z 是白化后的数据

% --- C. 迭代寻找解混矩阵 (FastICA Loop) ---
% 我们需要分离出 2 个分量
num_components = 2;
W = zeros(num_components, 2); % 权重矩阵

for i = 1:num_components
    w = randn(2, 1);    % 随机初始化权重
    w = w / norm(w);    % 归一化
    
    % 迭代最大化非高斯性
    for iter = 1:100
        w_old = w;
        
        % 非线性函数 g(u) = tanh(u) 和它的导数 g'(u) = 1 - tanh^2(u)
        u = w' * Z;
        g = tanh(u);
        g_prime = 1 - g.^2;
        
        % 核心迭代公式 (基于负熵的牛顿法)
        % w = E[Z * g(w'Z)] - E[g'(w'Z)] * w
        w = (Z * g')/N - mean(g_prime) * w;
        
        % 正交化 (Deflation): 减去已经找到的分量的投影
        % 确保找到的是不同的信号
        if i > 1
            w = w - W(1:i-1,:)' * W(1:i-1,:) * w;
        end
        w = w / norm(w); % 归一化
        
        % 收敛判断 (方向一致即可，符号不重要)
        if abs(abs(w' * w_old) - 1) < 1e-6
            break;
        end
    end
    W(i, :) = w'; % 存储找到的权重
end

%% 4. 恢复信号
S_estimated = W * Z;

%% 5. 绘图对比
figure('Position', [100, 100, 1000, 800]);

subplot(3,2,1); plot(t, s1); title('源信号 1 (雷达脉冲)'); grid on;
subplot(3,2,2); plot(t, s2); title('源信号 2 (持续干扰)'); grid on;

subplot(3,2,3); plot(t, X(1,:)); title('麦克风 1 (混合信号)'); grid on;
subplot(3,2,4); plot(t, X(2,:)); title('麦克风 2 (混合信号)'); grid on;

% 注意：ICA 恢复出来的信号，幅度和正负号是不确定的
subplot(3,2,5); plot(t, S_estimated(1,:)); title('ICA 分离结果 1'); grid on;
subplot(3,2,6); plot(t, S_estimated(2,:)); title('ICA 分离结果 2'); grid on;

sgtitle('ICA 盲源分离过程演示');