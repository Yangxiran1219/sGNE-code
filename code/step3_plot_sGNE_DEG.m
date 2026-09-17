%% =========================================================================
%  Step 3: Compute and plot global sGNE for STAD
%  ---------------------------------------------------------------
%  前置：先运行 step1_generate_samples.m 和
%       step2_compute_local_sGNE_matrix.m
%
%  输入：
%    - result/step2_intermediate_results.mat
%
%  输出：
%    - result/STAD_sGNE_Progression.png
% =========================================================================

clear;
clc;
close all;

%% ===================== 路径配置 =====================
code_path  = fileparts(mfilename('fullpath'));
base_path  = fullfile(code_path, '..');
result_dir = fullfile(base_path, 'result');

%% ===================== 加载 Step 2 结果 =====================
fprintf('[INFO] Loading Step 2 intermediate results...\n');

load(fullfile(result_dir, 'step2_intermediate_results.mat'));

% 使用的变量：
%   stage_patient_sGNE  - 各阶段每个样本的 sGNE 值
%   GlobalEntropy_stage - 各阶段 global sGNE 均值
%   case_stages         - 肿瘤分期名称
%   cancer_type         - 癌症类型

stagelabel = case_stages;
num_stage  = length(stagelabel);

%% ===================== 计算各阶段 global sGNE =====================
mean_sGNE = GlobalEntropy_stage;

fprintf('\n[INFO] Global sGNE by stage:\n');

for k = 1:num_stage
    n_k = sum(~isnan(stage_patient_sGNE{k}));

    fprintf('  %s: %.4f  (n=%d)\n', ...
        stagelabel{k}, mean_sGNE(k), n_k);
end

%% ===================== 绘制 sGNE 曲线 =====================
col_sGNE = [0.55, 0.15, 0.20];

figure( ...
    'Color', 'w', ...
    'Position', [300, 300, 900, 600]);

plot(1:num_stage, mean_sGNE, '-p', ...
    'Color', col_sGNE, ...
    'LineWidth', 6.5, ...
    'MarkerFaceColor', col_sGNE, ...
    'MarkerEdgeColor', col_sGNE, ...
    'MarkerSize', 10);

%% ===================== 坐标轴设置 =====================
xticks(1:num_stage);
xticklabels(stagelabel);
xtickangle(0);

xlabel('Stage', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

ylabel('sGNE', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

valid_sGNE = mean_sGNE(~isnan(mean_sGNE));

if ~isempty(valid_sGNE)
    y_min = min(valid_sGNE);
    y_max = max(valid_sGNE);

    if y_min == y_max
        y_margin = max(abs(y_min) * 0.1, 0.1);
        ylim([y_min - y_margin, y_max + y_margin]);
    else
        y_margin = 0.1 * (y_max - y_min);
        ylim([y_min - y_margin, y_max + y_margin]);
    end
end

ax = gca;
ax.YColor = col_sGNE;

set(gca, ...
    'LineWidth', 3, ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

grid off;
box off;

title([cancer_type ' sGNE Progression'], ...
    'FontWeight', 'bold', ...
    'FontSize', 16);

%% ===================== 保存图片 =====================
figure_file = fullfile( ...
    result_dir, ...
    sprintf('%s_sGNE_Progression.png', cancer_type));

saveas(gcf, figure_file);

fprintf('\n[INFO] Figure saved: %s_sGNE_Progression.png\n', ...
    cancer_type);

fprintf('[INFO] Step 3 complete.\n');
