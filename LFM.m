% ===================== 可观察的LFM信号生成与分析代码 =====================
% 作者：编程导师
% 功能：生成LFM信号，可视化时域波形、瞬时频率、脉冲压缩效果，突出核心特征
clc; clear; close all;

%% 1. 定义核心参数（物理意义清晰，方便观察）
fc = 10e6;          % 载波频率(10MHz，降低频率方便观察时域调制)
B = 5e6;            % LFM信号带宽(5MHz)
tau = 10e-6;        % 脉冲宽度(10μs)
Fs = 4*fc;          % 采样率(40MHz)，满足奈奎斯特采样定理
k = B/tau;          % 扫频斜率(MHz/μs)，k=B/τ 核心参数
t_total = 2*tau;    % 总仿真时间(2倍脉冲宽度，展示脉冲外的空白)
t = linspace(-tau/2, t_total-tau/2, round(t_total*Fs));  % 时间轴（中心对齐脉冲）

%% 2. 生成LFM信号
% (1) 生成基带LFM信号（无载波，方便观察频率变化）
LFM_base = zeros(size(t));
idx_pulse = find(abs(t) <= tau/2);  % 脉冲宽度内的时间索引
LFM_base(idx_pulse) = exp(1j * pi * k * t(idx_pulse).^2);  % 基带LFM公式：e^(jπkt²)

% (2) 生成载波调制后的LFM信号（贴近实际雷达发射信号）
LFM_mod = LFM_base .* exp(1j * 2*pi * fc * t);  % 乘以载波：e^(j2πfct)

%% 3. 计算瞬时频率（核心：观察频率线性变化）
% 希尔伯特变换求瞬时相位 → 求导得到瞬时频率
phase_LFM = unwrap(angle(hilbert(LFM_base)));  % 解缠绕相位（避免相位跳变）
f_inst = (1/(2*pi)) * diff(phase_LFM) * Fs;    % 瞬时频率 = 相位导数/(2π)
f_inst = [f_inst, f_inst(end)];                % 补全长度，匹配时间轴

%% 4. 脉冲压缩（LFM核心优势：提升距离分辨率）
% 匹配滤波器：基带LFM信号的共轭反转
mf = conj(flip(LFM_base));  
% 卷积实现匹配滤波（脉冲压缩）
LFM_compressed = conv(LFM_mod, mf, 'same');    
% 归一化幅度，方便对比
amp_mod = abs(LFM_mod);
amp_compressed = abs(LFM_compressed)/max(abs(LFM_compressed));

%% 5. 多维度可视化（核心：直观观察LFM特征）
figure('Color','w','Position',[100,100,1200,900],'Name','LFM信号特征观察');

% 子图1：基带LFM信号（实部+虚部）
subplot(3,2,1);
plot(t*1e6, real(LFM_base), 'b-', 'LineWidth',1.2);
hold on;
%plot(t*1e6, imag(LFM_base), 'r--', 'LineWidth',1.2);
xlabel('时间 (μs)'); ylabel('幅度');
title('基带LFM信号（实部=蓝色，虚部=红色）','FontSize',12);
xlim([-5, 15]); grid on;
%legend('实部','虚部','Location','best');

% 子图2：LFM瞬时频率（核心：线性变化）
subplot(3,2,2);
plot(t*1e6, f_inst/1e6, 'g-', 'LineWidth',1.5);
xlabel('时间 (μs)'); ylabel('瞬时频率 (MHz)');
title('LFM瞬时频率（线性扫频）','FontSize',12);
xlim([-5, 15]); ylim([-B/2/1e6-0.5, B/2/1e6+0.5]); grid on;
text(0, 0, ['扫频斜率 k=',num2str(k/1e6),' MHz/μs'], 'Color','red','FontSize',10);

% 子图3：载波调制后的LFM信号（实部）
subplot(3,2,3);
plot(t*1e6, real(LFM_mod), 'k-', 'LineWidth',1);
xlabel('时间 (μs)'); ylabel('幅度');
title('载波调制后LFM信号（实部）','FontSize',12);
xlim([-5, 15]); grid on;

% 子图4：调制后LFM幅度（脉冲形状）
subplot(3,2,4);
plot(t*1e6, amp_mod, 'b-', 'LineWidth',1.5);
xlabel('时间 (μs)'); ylabel('归一化幅度');
title('LFM脉冲幅度（压缩前）','FontSize',12);
xlim([-5, 15]); ylim([-0.2, 1.2]); grid on;
text(0, 0.8, ['脉冲宽度 τ=',num2str(tau*1e6),'μs'], 'Color','red','FontSize',10);

% 子图5：脉冲压缩后幅度（核心：窄脉冲）
subplot(3,2,5);
plot(t*1e6, amp_compressed, 'r-', 'LineWidth',1.5);
xlabel('时间 (μs)'); ylabel('归一化幅度');
title('LFM脉冲压缩后幅度（分辨率提升）','FontSize',12);
xlim([-5, 15]); ylim([-0.2, 1.2]); grid on;
% 计算压缩后的脉冲宽度（理论值：1/B）
tau_compressed = 1/B;
text(0, 0.8, ['压缩后脉宽 ≈',num2str(tau_compressed*1e6),'μs'], 'Color','blue','FontSize',10);

% 子图6：参数总结
subplot(3,2,6);
axis off; grid off;
text(0.1, 0.9, 'LFM信号核心参数总结','FontSize',14,'FontWeight','bold');
text(0.1, 0.7, ['1. 扫频带宽 B = ',num2str(B/1e6),' MHz'],'FontSize',12);
text(0.1, 0.5, ['2. 脉冲宽度 τ = ',num2str(tau*1e6),' μs'],'FontSize',12);
text(0.1, 0.3, ['3. 时宽带宽积 D = τ×B = ',num2str(tau*B)] ,'FontSize',12);
text(0.1, 0.1, ['4. 压缩后分辨率 Δτ = 1/B = ',num2str(tau_compressed*1e6),' μs'] ,'FontSize',12);

%% 6. 关键信息输出（命令行）
disp('===================== LFM信号关键特征 =====================');
disp(['时宽带宽积 D = ',num2str(tau*B),' → D越大，脉冲压缩效果越好']);
disp(['理论距离分辨率 ΔR = c/(2B) = ',num2str(3e8/(2*B)/1000),' km（c=光速）']);
disp('核心观察点：瞬时频率线性变化，脉冲压缩后脉宽大幅变窄（分辨率提升）');