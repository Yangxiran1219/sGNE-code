%% =========================================================================
%  Step 2: Compute local sGNE matrix for STAD
%  ---------------------------------------------------------------
%  输入（data/ 目录下）：
%    - STAD_control.mat  （step1 输出）
%    - STAD_case.mat     （step1 输出）
%    - adjacency_matrix_STAD.csv
%    - STAD_genelist.txt
%
%  输出（result/ 目录下）：
%    - STAD_Top500SignalPairs_Local_sGNE_Matrix.csv
%    - step2_intermediate_results.mat（供 Step 3 使用）
% =========================================================================

clear; clc; close all;

%% ===================== 路径配置 =====================
code_path  = fileparts(mfilename('fullpath'));
base_path  = fullfile(code_path, '..');
data_dir   = fullfile(base_path, 'data');
result_dir = fullfile(base_path, 'result');
if ~exist(result_dir, 'dir'), mkdir(result_dir); end

%% ===================== STAD 参数 =====================
cancer_type    = 'STAD';
critical_stage = 'IIB';
topK_pairs     = 500;

%% ===================== 加载数据 =====================
fprintf('[INFO] Loading data...\n');

ctrl_file = fullfile(data_dir, [cancer_type '_control.mat']);
case_file = fullfile(data_dir, [cancer_type '_case.mat']);

tempcontrol      = load(ctrl_file).control_reshaped;       % (3 x n_gene)
data2            = load(case_file);
result_case_cell = data2.result_case_cell;
case_ids_cell    = data2.case_ids_cell;
case_stages      = data2.case_stages;

num_stages = length(case_stages);

%% ===================== 基因列表 =====================
gene_file = fullfile(data_dir, sprintf('%s_genelist.txt', cancer_type));
fid = fopen(gene_file, 'r');
gene_list = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
gene_list = string(gene_list{1});
n_gene    = length(gene_list);

%% ===================== 网络邻接 =====================
T = csvread(fullfile(data_dir, ['adjacency_matrix_' cancer_type '.csv']));
T_unique = T(T(:,1) < T(:,2), :);
local_network = construct_local_network(max(T(:)), T);

fprintf('[INFO] %d genes, %d edges, %d stages\n', ...
    n_gene, size(T_unique,1), num_stages);

%% ===================== 并行池 =====================
if isempty(gcp('nocreate'))
    parpool('local');
end

%% ===================== 各阶段循环：计算每样本 pair 熵 =====================
signal_pairs_cell         = cell(num_stages, 1);
stage_patient_maps        = cell(num_stages, 1);
stage_sample_ids_filtered = cell(num_stages, 1);
stage_patient_sGNE        = cell(num_stages, 1);

for k = 1:num_stages
    tempcase_k = result_case_cell{k};
    if isempty(tempcase_k), continue; end

    sample_ids_k = case_ids_cell{k};
    sample_ids_k = sample_ids_k(:);
    num_patients = size(tempcase_k, 3);

    patient_maps  = cell(num_patients, 1);
    patient_pairs = cell(num_patients, 1);
    patient_sGNE  = nan(num_patients, 1);

    parfor p = 1:num_patients
        temp_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
        patient_data = tempcase_k(:,:,p);

        [ge, pair_idx, pair_vals] = ...
            Global_Entropy_with_top10_new(T_unique, local_network, tempcontrol, patient_data);

        patient_sGNE(p)  = ge;
        patient_pairs{p} = pair_idx;

        for r = 1:size(pair_idx,1)
            key = sprintf('%d_%d', pair_idx(r,1), pair_idx(r,2));
            temp_map(key) = pair_vals(r);
        end
        patient_maps{p} = temp_map;
    end

    signal_pairs_cell{k}         = patient_pairs;
    stage_patient_maps{k}        = patient_maps;
    stage_sample_ids_filtered{k} = sample_ids_k;
    stage_patient_sGNE{k}        = patient_sGNE;

    fprintf('[INFO] STAD %s: %d patients\n', ...
        case_stages{k}, sum(~isnan(patient_sGNE)));
end

%% ===================== 临界期 Top500 高频信号基因对 =====================
idx_stage = find(strcmp(case_stages, critical_stage));
stage_pairs = signal_pairs_cell{idx_stage};
num_patients_crit = length(stage_pairs);

freq_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

for p = 1:num_patients_crit
    pairs = stage_pairs{p};
    if isempty(pairs), continue; end
    pairs = sort(pairs, 2);   % 无向统一

    for r = 1:size(pairs,1)
        key = sprintf('%d_%d', pairs(r,1), pairs(r,2));
        if ~isKey(freq_map, key)
            freq_map(key) = 1;
        else
            freq_map(key) = freq_map(key) + 1;
        end
    end
end

keys_all = keys(freq_map);
vals_all = cellfun(@(x) freq_map(x), keys_all);
vals_all = vals_all / num_patients_crit * 100;   % 转百分比
[~, order] = sort(vals_all, 'descend');
keys_sorted = keys_all(order);

topK = min(topK_pairs, length(keys_sorted));
keys_top = keys_sorted(1:topK);

g1 = zeros(topK, 1);
g2 = zeros(topK, 1);
for i = 1:topK
    tmp = sscanf(keys_top{i}, '%d_%d');
    g1(i) = tmp(1);
    g2(i) = tmp(2);
end
top500_pairs = [g1(:), g2(:)];

fprintf('[INFO] Top %d signal pairs from critical stage %s\n', topK, critical_stage);

%% ===================== 输出 CSV：local sGNE 大矩阵 =====================
% 第一行：Gene1, Gene2, sample_id_1, sample_id_2, ...
% 第二行：,      ,      stage_1,     stage_1,     ...
% 数据行：gene1, gene2, val, val, ...

all_sample_ids    = {};
all_sample_stages = {};

for k = 1:num_stages
    ids_k = stage_sample_ids_filtered{k};
    if isempty(ids_k), continue; end
    n_k = length(ids_k);
    for p = 1:n_k
        if iscell(ids_k)
            all_sample_ids{end+1} = char(ids_k(p));       %#ok<SAGROW>
        else
            all_sample_ids{end+1} = ids_k(p);              %#ok<SAGROW>
        end
        all_sample_stages{end+1} = case_stages{k};         %#ok<SAGROW>
    end
end

n_total_samples = length(all_sample_ids);

% 构建大矩阵
sGNE_matrix = nan(topK, n_total_samples);

col = 0;
for k = 1:num_stages
    ids_k = stage_sample_ids_filtered{k};
    patient_maps_k = stage_patient_maps{k};
    if isempty(ids_k), continue; end
    n_k = length(ids_k);

    for p = 1:n_k
        col = col + 1;
        map_p = patient_maps_k{p};
        if isempty(map_p), continue; end
        for i = 1:topK
            if isKey(map_p, keys_top{i})
                sGNE_matrix(i, col) = map_p(keys_top{i});
            end
        end
    end
end

% 基因对名称
pair_gene1 = cell(topK, 1);
pair_gene2 = cell(topK, 1);
for i = 1:topK
    pair_gene1{i} = char(gene_list(top500_pairs(i,1)));
    pair_gene2{i} = char(gene_list(top500_pairs(i,2)));
end

% 写入 CSV
csv_path = fullfile(result_dir, 'STAD_Top500SignalPairs_Local_sGNE_Matrix.csv');
fid = fopen(csv_path, 'w');

% Row 1: 样本名
fprintf(fid, 'Gene1,Gene2');
for j = 1:n_total_samples
    fprintf(fid, ',%s', all_sample_ids{j});
end
fprintf(fid, '\n');

% Row 2: 时期
fprintf(fid, ',');
for j = 1:n_total_samples
    fprintf(fid, ',%s', all_sample_stages{j});
end
fprintf(fid, '\n');

% 数据行
for i = 1:topK
    fprintf(fid, '%s,%s', pair_gene1{i}, pair_gene2{i});
    for j = 1:n_total_samples
        if isnan(sGNE_matrix(i, j))
            fprintf(fid, ',');
        else
            fprintf(fid, ',%.6f', sGNE_matrix(i, j));
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('[INFO] CSV saved: STAD_Top500SignalPairs_Local_sGNE_Matrix.csv\n');
fprintf('[INFO] Matrix: %d pairs x %d samples\n', topK, n_total_samples);

%% ===================== 保存中间结果供 Step 3 =====================
GlobalEntropy_stage = nan(1, num_stages);
for k = 1:num_stages
    ge_k = stage_patient_sGNE{k};
    if ~isempty(ge_k)
        GlobalEntropy_stage(k) = mean(ge_k, 'omitnan');
    end
end

save(fullfile(result_dir, 'step2_intermediate_results.mat'), ...
    'stage_patient_sGNE', ...
    'stage_sample_ids_filtered', ...
    'GlobalEntropy_stage', ...
    'case_stages', ...
    'critical_stage', ...
    'cancer_type', ...
    'num_stages', ...
    'n_gene', ...
    'top500_pairs', ...
    'pair_gene1', ...
    'pair_gene2', ...
    'keys_top', ...
    'topK');

fprintf('[INFO] Intermediate results saved.\n');

delete(gcp('nocreate'));
fprintf('[INFO] Step 2 complete.\n');

%% ===================== FUNCTIONS =====================

% --- Global Entropy with top10% pairs ---
function [global_entropy, top_pairs, top_vals] = Global_Entropy_with_top10_new(T_unique, local_network, tempcontrol, tempcase)
N = size(T_unique, 1);
vals = zeros(N, 1);
for idx = 1:N
    i = T_unique(idx, 1);
    j = T_unique(idx, 2);
    [le1, le2] = Local_Entropy(i, j, T_unique, tempcontrol, tempcase, local_network);

    sd1 = abs(std(tempcase(:, i), 0, 'all') - std(tempcontrol(:, i), 0, 'all'));
    sd2 = abs(std(tempcase(:, j), 0, 'all') - std(tempcontrol(:, j), 0, 'all'));

    try
        pcc_c = abs(corr(tempcontrol(:, i), tempcontrol(:, j), 'Rows', 'complete'));
        pcc_a = abs(corr(tempcase(:, i), tempcase(:, j), 'Rows', 'complete'));
    catch
        pcc_c = 0; pcc_a = 0;
    end
    pcc_diff = abs(pcc_a - pcc_c);

    vals(idx) = ((sd1*le1 + sd2*le2) / (sd1 + sd2 + eps)) * exp(abs(pcc_diff));
end
pos = find(vals > 0);
[~, ord] = sort(vals(pos), 'descend');
k = max(1, round(0.10 * numel(pos)));
top_idx = pos(ord(1:k));
top_vals = vals(top_idx);
top_pairs = sort(T_unique(top_idx, :), 2);
global_entropy = mean(top_vals);
end

% --- 构建局部网络 ---
function local_network = construct_local_network(n, T)
local_network = cell(1, n);
for i = 1:n
    neighbors = T(T(:,1) == i, 2);
    local_network{i} = [i; neighbors];
end
end

% --- 节点局部网络熵 ---
function [local_entropy_1, local_entropy_2] = Local_Entropy(center1, center2, T, tempcontrol, tempcase, local_network)
local_entropy_1 = compute_entropy_for_node(center1, center2, local_network, tempcontrol, tempcase);
local_entropy_2 = compute_entropy_for_node(center2, center1, local_network, tempcontrol, tempcase);
end

function local_entropy = compute_entropy_for_node(center, other_center, local_network, tempcontrol, tempcase)
neighbors = local_network{center}(2:end);
neighbors(neighbors == other_center) = [];

cs_add = zeros(1, length(neighbors));
cs = zeros(1, length(neighbors));

if isempty(neighbors)
    local_entropy = 0;
    return;
end

for i = 1:length(neighbors)
    y_add = squeeze(tempcase(:, neighbors(i)));
    y     = squeeze(tempcontrol(:, neighbors(i)));
    x_add = squeeze(tempcase(:, center));
    x     = squeeze(tempcontrol(:, center));

    y_prev_add = y_add(1:end-1);
    y_now_add  = y_add(2:end);
    x_prev_add = x_add(1:end-1);

    y_prev = y(1:end-1);
    y_now  = y(2:end);
    x_prev = x(1:end-1);

    cs_add(i) = Granger_causality(y_prev_add, y_now_add, x_prev_add);
    cs(i)     = Granger_causality(y_prev, y_now, x_prev);
end

cs_add(cs_add == 0) = eps;
cs(cs == 0) = eps;
p_add = 1 ./ cs_add; p_add = p_add / sum(p_add);
p     = 1 ./ cs;     p     = p / sum(p);

entropy_add = -sum(p_add .* log(p_add + eps));
entropy     = -sum(p .* log(p + eps));

local_entropy = abs(entropy_add - entropy) * log(length(neighbors) + eps);
end

% --- Granger 因果函数 ---
function cs_value = Granger_causality(y_prev, y_now, x_prev)
eps_val = 1e-8;
max_cs_val = 1e3;

try
    X1 = [ones(size(y_prev)), y_prev];
    b1 = X1 \ y_now;
    y_hat1 = X1 * b1;
    rss1 = sum((y_now - y_hat1).^2) + eps_val;

    X2 = [ones(size(y_prev)), y_prev, x_prev];
    b2 = X2 \ y_now;
    y_hat2 = X2 * b2;
    rss2 = sum((y_now - y_hat2).^2) + eps_val;

    cs_value = log(rss1 / rss2);

    if ~isfinite(cs_value) || cs_value < 0
        cs_value = max_cs_val;
    end
catch
    cs_value = max_cs_val;
end
end
