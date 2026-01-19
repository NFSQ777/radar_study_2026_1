% ===================== LFM/FMCW 原理示意图生成代码 =====================
% 作者：编程导师
% 功能：生成原理级对比示意图，突出脉冲vs连续、分段扫频vs连续扫频的核心区别
clc; clear; close all;

%% 1. 定义核心参数（简化，方便绘图标注）
tau = 10;        % LFM脉冲宽度 (μs)
Tr = 50;         % LFM脉冲重复周期/ FMCW扫频周期 (μs)
f0 = 10;         % 起始频率 (GHz)
B = 1;           % 扫频带宽 (GHz)，简化为1GHz方便标注
k_LFM = B/tau;   % LFM扫频斜率 (GHz/μs)
k_FMCW = B/Tr;   % FMCW扫频斜率 (GHz/μs)
t_total = 2*Tr;  % 总绘图时间 (μs)
t = linspace(0, t_total, 1000);  % 时间轴

%% 2. 初始化绘图数据（原理性，非仿真）
% LFM时域（脉冲：1=有信号，0=无信号）
LFM_time = zeros(size(t));
LFM_time(t>=0 & t<=tau) = 1;       % 第一个脉冲
LFM_time(t>=Tr & t<=Tr+tau) = 1;   % 第二个脉冲
% LFM频率（仅脉冲内线性扫频）
LFM_freq = zeros(size(t));
LFM_freq(t>=0 & t<=tau) = f0 + k_LFM*(t(t>=0 & t<=tau));
LFM_freq(t>=Tr & t<=Tr+tau) = f0 + k_LFM*(t(t>=Tr & t<=Tr+tau)-Tr);

% FMCW时域（连续：全程1）
FMCW_time = ones(size(t));
% FMCW频率（全程线性扫频）
FMCW_freq = f0 + k_FMCW*t;

%% 3. 绘制对比示意图（2x2布局：时域+频率域）
figure('Color','w','Position',[100,100,1200,800],'Name','LFM/FMCW 原理示意图');

% ===== 子图1：LFM时域（幅度-时间）=====
subplot(2,2,1);
plot(t, LFM_time, 'b-', 'LineWidth',2);
xlabel('时间 (μs)','FontSize',12);
ylabel('归一化幅度','FontSize',12);
title('LFM 时域（脉冲体制）','FontSize',14,'FontWeight','bold');
xlim([0, t_total]); ylim([-0.2, 1.2]);
grid on; hold on;
% 标注关键参数
text(tau/2, 0.8, ['脉冲宽度 τ=',num2str(tau),'μs'], 'Color','red','FontSize',10);
text((tau+Tr)/2, 0.8, ['脉冲间隔 Tr-τ=',num2str(Tr-tau),'μs'], 'Color','green','FontSize',10);
text(Tr+tau/2, 0.8, '脉冲2', 'Color','red','FontSize',10);
% 绘制辅助线
plot([tau,tau], [-0.2,1.2], 'r--','LineWidth',1);
plot([Tr,Tr], [-0.2,1.2], 'g--','LineWidth',1);

% ===== 子图2：FMCW时域（幅度-时间）=====
subplot(2,2,2);
plot(t, FMCW_time, 'r-', 'LineWidth',2);
xlabel('时间 (μs)','FontSize',12);
ylabel('归一化幅度','FontSize',12);
title('FMCW 时域（连续波体制）','FontSize',14,'FontWeight','bold');
xlim([0, t_total]); ylim([-0.2, 1.2]);
grid on; hold on;
% 标注关键参数
text(Tr/2, 0.8, ['扫频周期 Ts=',num2str(Tr),'μs'], 'Color','blue','FontSize',10);
text(Tr+Tr/2, 0.8, '连续扫频无间隔', 'Color','blue','FontSize',10);

% ===== 子图3：LFM频率-时间 =====
subplot(2,2,3);
plot(t, LFM_freq, 'b-', 'LineWidth',2);
xlabel('时间 (μs)','FontSize',12);
ylabel('频率 (GHz)','FontSize',12);
title('LFM 瞬时频率（分段扫频）','FontSize',14,'FontWeight','bold');
xlim([0, t_total]); ylim([f0-0.2, f0+B+0.2]);
grid on; hold on;
% 标注扫频范围
text(tau/2, f0+B/2, ['扫频：',num2str(f0),'→',num2str(f0+B),'GHz'], 'Color','red','FontSize',10);
text(Tr+tau/2, f0+B/2, '同斜率扫频', 'Color','red','FontSize',10);
% 辅助线
plot([tau,tau], [f0-0.2,f0+B+0.2], 'r--','LineWidth',1);
plot([Tr,Tr], [f0-0.2,f0+B+0.2], 'g--','LineWidth',1);

% ===== 子图4：FMCW频率-时间 =====
subplot(2,2,4);
plot(t, FMCW_freq, 'r-', 'LineWidth',2);
xlabel('时间 (μs)','FontSize',12);
ylabel('频率 (GHz)','FontSize',12);
title('FMCW 瞬时频率（连续扫频）','FontSize',14,'FontWeight','bold');
xlim([0, t_total]); ylim([f0-0.2, f0+2*B+0.2]);
grid on; hold on;
% 标注扫频斜率和范围
text(Tr/2, f0+B/2, ['扫频斜率 k=',num2str(k_FMCW),'GHz/μs'], 'Color','blue','FontSize',10);
text(Tr+Tr/2, f0+B+B/2, ['持续扫频：',num2str(f0+B),'→',num2str(f0+2*B),'GHz'], 'Color','blue','FontSize',10);

% ===== 整体标注（核心区别）=====
annotation('textbox', [0.02,0.02,0.96,0.05], ...
    'String', '核心区别：LFM=脉冲式分段扫频 | FMCW=连续式全程扫频 | 两者带宽一致（决定距离分辨率）', ...
    'FontSize',12, 'Color','black', 'BackgroundColor','#f0f0f0', 'EdgeColor','none');






