%% =========================================================================
%  Step 1: Generate STAD control & case samples from xlsx
%  ---------------------------------------------------------------
%  输入（data/ 目录下）：
%    - STAD_tpm.xlsx
%
%  输出（data/ 目录下）：
%    - STAD_control.mat  （含 control_reshaped, normal_individual）
%    - STAD_case.mat     （含 result_case_cell, case_stages, case_ids_cell）
% =========================================================================

clear; clc; close all;

%% ===================== 路径配置 =====================
code_path  = fileparts(mfilename('fullpath'));
base_path  = fullfile(code_path, '..');
data_dir   = fullfile(base_path, 'data');

%% ===================== STAD 参数 =====================
cancer_type     = 'STAD';
control_stages  = {'IA','IB'};
case_stages     = {'IIA','IIB','IIIA','IIIB','IIIC','IV'};

%% ===================== 路径 =====================
xlsx_file   = fullfile(data_dir, [cancer_type '_tpm.xlsx']);
control_file = fullfile(data_dir, [cancer_type '_control.mat']);
case_file    = fullfile(data_dir, [cancer_type '_case.mat']);

% 删除旧文件
if exist(control_file, 'file'), delete(control_file); end
if exist(case_file, 'file'),    delete(case_file); end

%% ===================== 读取表达数据 =====================
fprintf('[INFO] 读取 %s ...\n', xlsx_file);

normal_cell = readcell(xlsx_file, 'Sheet', 'normal');
normal_data = cell2mat(normal_cell(2:end, 2:end));     % (genes x n_normal)
gene_names  = normal_cell(2:end, 1);

tumor_cell   = readcell(xlsx_file, 'Sheet', 'tumor');
tumor_ids    = string(tumor_cell(1, 2:end));            % 样本ID
tumor_ids    = strtrim(tumor_ids);
tumor_stages = string(tumor_cell(2, 2:end));
tumor_stages = strtrim(replace(tumor_stages, "Stage ", ""));
tumor_data   = cell2mat(tumor_cell(3:end, 2:end));      % (genes x n_tumor)

fprintf('[INFO] normal: %d samples, tumor: %d samples, genes: %d\n', ...
    size(normal_data,2), size(tumor_data,2), size(normal_data,1));

%% ===================== control 构建 =====================
normal_mean = mean(normal_data, 2, 'omitnan');           % (genes x 1)

control_means = zeros(size(normal_data,1), length(control_stages));
for s = 1:length(control_stages)
    idx = strcmp(tumor_stages, control_stages{s});
    control_means(:, s) = mean(tumor_data(:, idx), 2, 'omitnan');
end

control_data     = [normal_mean, control_means];          % (genes x 3)
control_reshaped = control_data';                          % (3 x genes)

% 个体正常样本（转置为 n_normal x genes），供后续 DEG t-test
normal_individual = normal_data';                          % (n_normal x genes)

save(control_file, 'control_reshaped', 'normal_individual');
fprintf('[INFO] 保存 control: [%d x %d]\n', size(control_reshaped));

%% ===================== case 构建（逐患者）=====================
num_stages  = length(case_stages);
num_genes   = size(control_data, 1);
num_control = size(control_data, 2);   % 3

result_case_cell = cell(num_stages, 1);
case_ids_cell    = cell(num_stages, 1);

for j = 1:num_stages
    stage = case_stages{j};
    idx = strcmp(tumor_stages, stage);

    stage_data = tumor_data(:, idx);                       % (genes x n_patients)
    num_patients = size(stage_data, 2);
    stage_ids = tumor_ids(idx);

    if num_patients == 0
        warning('STAD %s 无样本，跳过', stage);
        continue;
    end

    % 三维数组：(control+1) x genes x patients
    stage_case_array = zeros(num_control + 1, num_genes, num_patients);

    for p = 1:num_patients
        patient_expr = stage_data(:, p);                   % (genes x 1)
        combined = [control_data, patient_expr];           % (genes x 4)
        stage_case_array(:,:,p) = combined';               % (4 x genes)
    end

    result_case_cell{j} = stage_case_array;
    case_ids_cell{j}    = stage_ids;

    fprintf('  %s: %d 个病人\n', stage, num_patients);
end

save(case_file, 'result_case_cell', 'case_stages', 'case_ids_cell');
fprintf('[INFO] 保存 case（病人级）完成\n');

fprintf('[INFO] Step 1 complete.\n');
