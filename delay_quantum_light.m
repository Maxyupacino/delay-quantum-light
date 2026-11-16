%% 双模压缩光错位扫描 + 0-0.5MHz频段最大差值提取 (滑动平均去噪版)
% 针对 100 万组数据级别优化
% 数据格式：Excel (.xlsx) 三列，第一行表头，列1=时间(s)，列2=通道1电压(V)，列3=通道2电压(V)
% 采样率 10 MHz，每错位1格 = 1个采样点 = 100 ns

clear; clc; close all;

%% ==================== 1. 读取 Excel 数据 ====================
% 提示：一百万行 Excel 数据读取可能需要数十秒到一两分钟，请耐心等待。
fprintf('请选择压缩光数据文件...\n');
[file_sqz, path_sqz] = uigetfile('*.xlsx', '选择压缩光数据文件');
if isequal(file_sqz, 0)
    error('未选择压缩光文件，程序终止。'); 
end
fprintf('正在读取压缩光数据，数据量较大(100万行)，请稍候...\n');
data_sqz = readmatrix(fullfile(path_sqz, file_sqz), 'NumHeaderLines', 1);

fprintf('请选择相干光数据文件...\n');
[file_coh, path_coh] = uigetfile('*.xlsx', '选择相干光数据文件');
if isequal(file_coh, 0)
    error('未选择相干光文件，程序终止。');
end
fprintf('正在读取相干光数据，请稍候...\n');
data_coh = readmatrix(fullfile(path_coh, file_coh), 'NumHeaderLines', 1);

% 提取电压列（第2列通道1，第3列通道2），并转为列向量
sig_sqz = data_sqz(:, 2); idl_sqz = data_sqz(:, 3);
sig_coh = data_coh(:, 2); idl_coh = data_coh(:, 3);

% 转为列向量
sig_sqz = sig_sqz(:); idl_sqz = idl_sqz(:);
sig_coh = sig_coh(:); idl_coh = idl_coh(:);

% ================= 新增：统一通道符号 =================
% 强制将所有电压通道转换为绝对值（或者你也可以明确只给负的那一列乘 -1）
sig_sqz = abs(sig_sqz);
idl_sqz = abs(idl_sqz);
sig_coh = abs(sig_coh);
idl_coh = abs(idl_coh);
% ===================================================

fprintf('数据读取并已统一符号！开始处理...\n');

%% ==================== 2. 参数设置 ====================
Fs = 10e6;                  % 采样率 10 MHz
Ts = 1/Fs;                  % 采样间隔 100 ns
delay_grids = 0:20;         % 压缩光错位网格数列表（这里设定扫描 0 到 20 格）

nfft = 2^14;                % 每段 FFT 点数 (频率分辨率约为 610 Hz)
window = hanning(nfft);     % 汉宁窗
noverlap = 0;               % 分段无重叠
R = 50;                     % 50 Ω 系统阻抗

%% ==================== 3. 计算相干光基准频谱 ====================
% 对齐相干光通道数据长度，计算差分
len_coh = min(length(sig_coh), length(idl_coh));
diff_coh = sig_coh(1:len_coh) - idl_coh(1:len_coh);

% 使用 pwelch 计算功率谱密度 (V^2/Hz)
[Pxx_coh, f_coh] = pwelch(diff_coh, window, [], [], Fs);
psd_coh_dBm = 10 * log10(Pxx_coh / (1e-3 * R)); % 转换为 dBm/Hz

%% ==================== 4. 逐延迟计算并绘制全频段原始图 ====================
figure('Name', '频谱对比图', 'Position', [100, 100, 1000, 600]);
plot(f_coh/1e6, psd_coh_dBm, 'k--', 'LineWidth', 1.5, 'DisplayName', '相干光基准 (0格)');
hold on; grid on;
xlabel('频率 (MHz)');
ylabel('噪声功率谱密度 (dBm/Hz)');
title('双模压缩光错位扫描噪声谱');
xlim([0, 4]); % 限制显示频率范围，重点看 0-4 MHz

colors = lines(length(delay_grids));
psd_sqz_all = zeros(length(f_coh), length(delay_grids));

fprintf('\n正在计算各延时下的压缩光功率谱...\n');
for idx = 1:length(delay_grids)
    grid_shift = delay_grids(idx);
    delay_ns = grid_shift * Ts * 1e9;
    
    % 根据延时截断数据以对齐
    len_sqz = min(length(sig_sqz), length(idl_sqz)) - grid_shift;
    if len_sqz <= 0
        error('延迟 %d 格时数据长度不足。', grid_shift);
    end
    
    % 计算延时差分
    if grid_shift == 0
        diff_sqz = sig_sqz(1:len_sqz) - idl_sqz(1:len_sqz);
    else
        % idl 通道向后推迟 grid_shift 格
        diff_sqz = sig_sqz(1:len_sqz) - idl_sqz(1+grid_shift : len_sqz+grid_shift);
    end
    
    % 计算压缩光的 PSD
    [Pxx_sqz, ~] = pwelch(diff_sqz, window, [], [], Fs);
    psd_sqz_dBm = 10 * log10(Pxx_sqz / (1e-3 * R));
    psd_sqz_all(:, idx) = psd_sqz_dBm;
    
    % 绘制当前延时下的 PSD 曲线
    plot(f_coh/1e6, psd_sqz_dBm, 'Color', colors(idx,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('延迟 %d 格 (%.0f ns)', grid_shift, delay_ns));
end
hold off;
legend('Location', 'eastoutside');

%% ==================== 5. 提前计算 ±3000ns 延时下的所有 PSD 曲线 ====================
fprintf('\n=== 步骤 1：开始预计算正负延时 (-3000ns 到 3000ns) 的全频段 PSD 矩阵 ===\n');

% 设定延时网格范围：-30 到 30格（对应 -3000 ns 到 3000 ns）
delay_grids_full = -10:10; 
num_delays = length(delay_grids_full);
delay_times_full = delay_grids_full * Ts * 1e9;

% 预分配空间存储所有的 PSD 曲线矩阵 (行=频率点，列=延时点)
psd_sqz_all_full = zeros(length(f_coh), num_delays);
psd_coh_all_full = zeros(length(f_coh), num_delays);

% 数据总长度
N_sqz = min(length(sig_sqz), length(idl_sqz));
N_coh = min(length(sig_coh), length(idl_coh));

for i = 1:num_delays
    shift = delay_grids_full(i);
    
    % ----------------- A. 处理压缩光延时切片 -----------------
    if shift > 0
        diff_sqz_full = sig_sqz(1 : N_sqz-shift) - idl_sqz(1+shift : N_sqz);
    elseif shift < 0
        abs_shift = abs(shift);
        diff_sqz_full = sig_sqz(1+abs_shift : N_sqz) - idl_sqz(1 : N_sqz-abs_shift);
    else
        diff_sqz_full = sig_sqz(1:N_sqz) - idl_sqz(1:N_sqz);
    end
    [Pxx_sqz, ~] = pwelch(diff_sqz_full, window, [], [], Fs);
    psd_sqz_all_full(:, i) = 10 * log10(Pxx_sqz / (1e-3 * R));
    
    % ----------------- B. 处理相干光延时切片 (严格对照) -----------------
    if shift > 0
        diff_coh_full = sig_coh(1 : N_coh-shift) - idl_coh(1+shift : N_coh);
    elseif shift < 0
        abs_shift = abs(shift);
        diff_coh_full = sig_coh(1+abs_shift : N_coh) - idl_coh(1 : N_coh-abs_shift);
    else
        diff_coh_full = sig_coh(1:N_coh) - idl_coh(1:N_coh);
    end
    [Pxx_coh, ~] = pwelch(diff_coh_full, window, [], [], Fs);
    psd_coh_all_full(:, i) = 10 * log10(Pxx_coh / (1e-3 * R));
    
    % 打印进度
    if mod(shift, 10) == 0
        fprintf('已完成计算: 延迟 %3d 格 (%5.0f ns)\n', shift, delay_times_full(i));
    end
end

%% ==================== 6. 同延时 PSD 曲线作差，提取 0-0.25 MHz 最大差值 ====================
fprintf('\n=== 步骤 2：基于同延时基准作差，提取最大压缩度 (5点平滑) ===\n');

% 设定频段：50kHz - 250 kHz (避开 0Hz 直流偏置造成的假性高值)
f_band = [50e3, 0.25e6]; 
idx_band = (f_coh >= f_band(1)) & (f_coh <= f_band(2));

max_diffs_full = zeros(1, num_delays);

for i = 1:num_delays
    % 提取当前延时下，该频段内的相干光和压缩光 PSD
    coh_band = psd_coh_all_full(idx_band, i);
    sqz_band = psd_sqz_all_full(idx_band, i);
    
    % 1. 同延时、同频作差：当前延时的相干光 - 当前延时的压缩光
    diff_band_raw = coh_band - sqz_band;
    
    % 2. 5点滑动平均处理：消除散粒噪声的偶然波峰波谷
    diff_band_smoothed = movmean(diff_band_raw, 5);
    
    % 3. 提取该延时下平滑后的最大压缩度
    max_diffs_full(i) = max(diff_band_smoothed);
end

fprintf('提取完成！正在绘制最终曲线...\n');

%% ==================== 7. 绘制最大压缩度与正负延时关系的曲线 ====================
pause(0.1); % 稍微停顿，释放之前的绘图内存避免渲染树崩溃

fig = figure('Name', '同延时基准对照：最大压缩度 (-3000ns 到 3000ns)', 'Position', [200, 200, 900, 500]);
ax = axes(fig); % 显式获取坐标区

% 绘制提取出的最大差值曲线
plot(ax, delay_times_full, max_diffs_full, 'bo-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
hold(ax, 'on');

% 画一条中心红线标识 0 ns 处
xline(ax, 0, 'r--', 'LineWidth', 1.5, 'Label', '0 ns (无延时)');

grid(ax, 'on');
xlabel(ax, '延时时间 \tau (ns)', 'FontSize', 12);
ylabel(ax, '最大差值 [同延时相干光 - 压缩光] (dB)', 'FontSize', 12);
title(ax, '0 - 0.25 MHz 频段内最大压缩度随延时变化 (严格对照同延时基准)', 'FontSize', 14);

set(ax, 'FontSize', 11, 'LineWidth', 1.2);
drawnow; % 强制渲染刷新

%% ==================== 8. 频域谱逆傅里叶变换 (时间域自相关函数) ====================
fprintf('\n=== 步骤 3：对 0 延时频域数据进行逆傅里叶变换 (IFFT) ===\n');

% 提取 0 延时下完全对齐的差分信号
N_sqz_0 = min(length(sig_sqz), length(idl_sqz));
diff_sqz_0 = sig_sqz(1:N_sqz_0) - idl_sqz(1:N_sqz_0);

% 1. 计算双边快速傅里叶变换 (FFT)
X_f = fft(diff_sqz_0);
% 2. 计算双边能量谱 (频域)
S_xx = (X_f .* conj(X_f)) / N_sqz_0; 
% 3. 逆傅里叶变换 (IFFT) 得到时间域自相关函数 R(tau)
R_xx_ifft = ifft(S_xx);

% 将 0 延迟点移动到数组中心，以便于观察对称性
R_xx_ifft = fftshift(R_xx_ifft); 

% 构建对应的时间延迟轴 (ns)
tau_axis_ifft = (-floor(N_sqz_0/2) : ceil(N_sqz_0/2)-1) * Ts * 1e9;

% 绘制时间域自相关图
pause(0.1);
fig8 = figure('Name', '频域逆变换 (时间域自相关)', 'Position', [250, 250, 800, 500]);
ax8 = axes(fig8);
plot(ax8, tau_axis_ifft, real(R_xx_ifft), 'k', 'LineWidth', 1.2);

xlim(ax8, [-3000, 3000]); % 限制观察正负 3000 ns 范围，与前面的延时对齐
grid(ax8, 'on');
xlabel(ax8, '延迟时间 \tau (ns)', 'FontSize', 12);
ylabel(ax8, '自相关幅度 R(\tau)', 'FontSize', 12);
title(ax8, '0延时差分信号功率谱的逆傅里叶变换 (时间域自相关)', 'FontSize', 14);
set(ax8, 'FontSize', 11, 'LineWidth', 1.2);
drawnow;


%% ==================== 9. 计算 Signal 与 Idler 的二阶关联函数 g(2) ====================
fprintf('\n=== 步骤 4：计算 Signal 与 Idler 的二阶互关联函数 g^{(2)}(tau) ===\n');

% 注意：计算 g(2) 需要包含直流 (DC) 分量的真实物理强度。
% 我们从原始的 data_sqz 矩阵中重新提取数据，并取绝对值以翻转负向硬件输出
I_sig = abs(data_sqz(:, 2)); 
I_idl = abs(data_sqz(:, 3));

% 截断至相同长度
N_g2 = min(length(I_sig), length(I_idl));
I_sig = I_sig(1:N_g2);
I_idl = I_idl(1:N_g2);

% 计算强度的均值 <I1> 和 <I2> (对应公式的分母部分)
mean_I_sig = mean(I_sig);
mean_I_idl = mean(I_idl);

% 设定计算互相关的最大延迟步数 (30步对应正负 3000ns)
max_lag = 30; 

% 计算未归一化的互相关 <I1(t) I2(t+tau)> (对应公式的分子部分)
% 使用 'unbiased' 参数可以消除由于数据两端截断带来的三角窗衰减效应
[R_12, lags] = xcorr(I_sig, I_idl, max_lag, 'unbiased');

% 按照图片中的正规序形式公式，进行归一化得到最终的 g(2)
g2_tau = R_12 / (mean_I_sig * mean_I_idl);

% 将 lag 步数转换为实际的纳秒时间轴
tau_g2_ns = lags * Ts * 1e9;

% 绘制 g(2) 曲线
pause(0.1);
fig9 = figure('Name', '二阶互关联函数 g(2)', 'Position', [300, 300, 800, 500]);
ax9 = axes(fig9);
plot(ax9, tau_g2_ns, g2_tau, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');

hold(ax9, 'on');
xline(ax9, 0, 'k--', 'LineWidth', 1.5, 'Label', '0 ns'); % 0延迟参考线

grid(ax9, 'on');
xlabel(ax9, '延迟时间 \tau (ns)', 'FontSize', 12);
ylabel(ax9, 'g^{(2)}_{12}(\tau)', 'FontSize', 12);
title(ax9, 'Signal 与 Idler 的二阶互关联函数 g^{(2)}(\tau)', 'FontSize', 14);
set(ax9, 'FontSize', 11, 'LineWidth', 1.2);
drawnow;

%% ==================== 10. 高精度插值：10ns 步长，-1000ns 到 1000ns ====================
fprintf('\n=== 步骤 5：开始进行 10ns 步长高精度延时扫描 (-1000ns 到 1000ns) ===\n');

% 1. 利用重采样 (resample) 将采样率提升 10 倍，实现 10ns 步长
factor = 10;
Fs_new = Fs * factor;       % 新采样率 100 MHz
Ts_new = 1 / Fs_new;        % 新采样间隔 10 ns

fprintf('正在对信号进行 10 倍上采样 (100万行数据插值较慢，请耐心等待数十秒)...\n');
% 直接使用已经在第1步去除了直流偏置、统一了符号的干净交流数据进行上采样
sig_coh_up = resample(sig_coh, factor, 1);
idl_coh_up = resample(idl_coh, factor, 1);
sig_sqz_up = resample(sig_sqz, factor, 1);
idl_sqz_up = resample(idl_sqz, factor, 1);

% 2. 调整 FFT 点数以保证 100MHz 采样率下的低频分辨率
% 原来 Fs=10MHz 用 2^14，现在 Fs_new=100MHz 放大到 2^17 (频率分辨率约为 762 Hz)
nfft_new = 2^17; 
window_new = hanning(nfft_new);

% 3. 设定 -1000 ns 到 1000 ns 的网格 (步距 10ns，对应 -100 到 +100 格)
delay_grids_fine = -100:100; 
num_delays_fine = length(delay_grids_fine);
delay_times_fine = delay_grids_fine * Ts_new * 1e9;
max_diffs_fine = zeros(1, num_delays_fine);

N_sqz_up = min(length(sig_sqz_up), length(idl_sqz_up));
N_coh_up = min(length(sig_coh_up), length(idl_coh_up));

% 提前计算一次 0 延时的频域轴，用于确定 50 kHz - 250 kHz 频段的索引
[~, f_coh_new] = pwelch(sig_coh_up(1:N_coh_up) - idl_coh_up(1:N_coh_up), window_new, [], [], Fs_new);
f_band = [50e3, 0.25e6]; 
idx_band_fine = (f_coh_new >= f_band(1)) & (f_coh_new <= f_band(2));

fprintf('插值完成！开始循环计算 %d 个高精度延时点...\n', num_delays_fine);

for i = 1:num_delays_fine
    shift = delay_grids_fine(i);
    
    % --- A. 处理压缩光高精度延时切片 ---
    if shift > 0
        diff_sqz_fine = sig_sqz_up(1 : N_sqz_up-shift) - idl_sqz_up(1+shift : N_sqz_up);
    elseif shift < 0
        abs_s = abs(shift);
        diff_sqz_fine = sig_sqz_up(1+abs_s : N_sqz_up) - idl_sqz_up(1 : N_sqz_up-abs_s);
    else
        diff_sqz_fine = sig_sqz_up(1:N_sqz_up) - idl_sqz_up(1:N_sqz_up);
    end
    [Pxx_sqz_fine, ~] = pwelch(diff_sqz_fine, window_new, [], [], Fs_new);
    sqz_band = 10 * log10(Pxx_sqz_fine(idx_band_fine) / (1e-3 * R));
    
    % --- B. 处理相干光高精度延时切片 (严格对照) ---
    if shift > 0
        diff_coh_fine = sig_coh_up(1 : N_coh_up-shift) - idl_coh_up(1+shift : N_coh_up);
    elseif shift < 0
        abs_s = abs(shift);
        diff_coh_fine = sig_coh_up(1+abs_s : N_coh_up) - idl_coh_up(1 : N_coh_up-abs_s);
    else
        diff_coh_fine = sig_coh_up(1:N_coh_up) - idl_coh_up(1:N_coh_up);
    end
    [Pxx_coh_fine, ~] = pwelch(diff_coh_fine, window_new, [], [], Fs_new);
    coh_band = 10 * log10(Pxx_coh_fine(idx_band_fine) / (1e-3 * R));
    
    % --- C. 频段内作差与 5 点平滑提取最大值 ---
    diff_band_raw = coh_band - sqz_band;
    diff_band_smoothed = movmean(diff_band_raw, 5);
    max_diffs_fine(i) = max(diff_band_smoothed);
    
    % 打印进度防假死
    if mod(shift, 20) == 0
        fprintf('高精度扫描进度: 延时 %5.0f ns, 最大压缩度 = %.3f dB\n', delay_times_fine(i), max_diffs_fine(i));
    end
end

%% ==================== 11. 绘制 10ns 高精度平滑曲线 ====================
pause(0.1); 
fig10 = figure('Name', '高精度 10ns 步长最大压缩度', 'Position', [450, 250, 900, 500]);
ax10 = axes(fig10);

% 绘制高精度曲线 (使用实线方便看出精细趋势)
plot(ax10, delay_times_fine, max_diffs_fine, 'g-', 'LineWidth', 2);
hold(ax10, 'on');

% 标记 0 ns 参考线
xline(ax10, 0, 'r--', 'LineWidth', 1.5, 'Label', '0 ns'); 
grid(ax10, 'on');

xlabel(ax10, '延时时间 \tau (ns)', 'FontSize', 12);
ylabel(ax10, '最大差值 [同延时相干光 - 压缩光] (dB)', 'FontSize', 12);
title(ax10, '0.05 - 0.25 MHz频段最大压缩度 (-1000ns 到 1000ns, 10ns 步长)', 'FontSize', 14);

% 坐标轴美化
set(ax10, 'FontSize', 11, 'LineWidth', 1.2);
drawnow;

%% ==================== 12. 导出 10ns 高精度数据至 CSV 文件 ====================
fprintf('\n=== 步骤 6：将高精度扫描数据导出为 CSV 文件 ===\n');

% 将之前算好的行向量转置 (使用 ') 为列向量，并组合成带有表头的表格
results_table_fine = table(delay_times_fine', max_diffs_fine', ...
    'VariableNames', {'Delay_ns', 'Max_Squeezing_dB'});

% 定义导出的文件名
csv_filename = 'Max_Squeezing_10ns_Resolution_1000ns.csv';

% 写入到当前 MATLAB 的工作路径下
writetable(results_table_fine, csv_filename);

% 获取当前完整路径并打印提示
current_path = pwd;
fprintf('高精度扫描数据导出成功！\n');
fprintf('文件已保存至: %s\n', fullfile(current_path, csv_filename));