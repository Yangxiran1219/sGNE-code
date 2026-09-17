# sGNE Analysis for STAD

## English

### Overview

This repository provides a three-step MATLAB workflow for calculating sample-specific gene network entropy (sGNE) in stomach adenocarcinoma (STAD) and visualizing its progression across tumor stages.

Run the scripts in order. The local sGNE results, intermediate results, and stage-level global sGNE progression figure will be saved in the `result/` directory.

### Requirements

- MATLAB
- Parallel Computing Toolbox (`parpool` and `parfor`)

### Repository structure

```text
sGNE-code/
├── code/       # MATLAB scripts
├── data/       # Input data and Step 1 generated files
└── result/     # sGNE results, intermediate results, and figures
```

### Usage

After downloading or cloning this repository, open MATLAB and change `project_dir` below to the actual location of the project on your computer:

```matlab
project_dir = '/path/to/sGNE-code';  % Change this path
cd(fullfile(project_dir, 'code'));
```

Then run the three scripts in the following order:

```matlab
run('step1_generate_samples.m');
run('step2_compute_local_sGNE_matrix.m');
run('step3_plot_sGNE_DEG.m');
```

1. `step1_generate_samples.m` reads `data/STAD_tpm.xlsx` and generates:
   - `data/STAD_control.mat`
   - `data/STAD_case.mat`
2. `step2_compute_local_sGNE_matrix.m` calculates local sGNE values and generates:
   - `result/STAD_Top500SignalPairs_Local_sGNE_Matrix.csv`
   - `result/step2_intermediate_results.mat`
3. `step3_plot_sGNE_DEG.m` calculates stage-level global sGNE and generates:
   - `result/STAD_sGNE_Progression.png`

### Path configuration

The scripts determine the `data/` and `result/` paths automatically from their location in the `code/` directory. If the repository structure is preserved, you only need to change `project_dir` in the MATLAB commands above. If you move the scripts or change the directory structure, update the **Path Configuration** section near the beginning of each script.

---

## 中文

### 项目简介

本仓库提供一个包含三个步骤的 MATLAB 分析流程，用于计算胃腺癌（STAD）样本的特异性基因网络熵（sample-specific Gene Network Entropy, sGNE），并绘制 sGNE 随肿瘤分期变化的曲线图。

请按照顺序运行三个脚本。局部 sGNE 结果、中间结果以及各分期 global sGNE 变化曲线图将保存在 `result/` 文件夹中。

### 运行环境

- MATLAB
- Parallel Computing Toolbox（使用 `parpool` 和 `parfor`）

### 文件结构

```text
sGNE-code/
├── code/       # MATLAB 脚本
├── data/       # 输入数据及 Step 1 生成的文件
└── result/     # sGNE 结果、中间结果和图像
```

### 使用方法

下载或克隆本仓库后，打开 MATLAB，并将下面的 `project_dir` 修改为项目在你电脑上的实际路径：

```matlab
project_dir = '/你的本地路径/sGNE-code';  % 请修改此路径
cd(fullfile(project_dir, 'code'));
```

然后按照以下顺序依次运行三个脚本：

```matlab
run('step1_generate_samples.m');
run('step2_compute_local_sGNE_matrix.m');
run('step3_plot_sGNE_DEG.m');
```

1. `step1_generate_samples.m` 读取 `data/STAD_tpm.xlsx`，并生成：
   - `data/STAD_control.mat`
   - `data/STAD_case.mat`
2. `step2_compute_local_sGNE_matrix.m` 计算局部 sGNE，并生成：
   - `result/STAD_Top500SignalPairs_Local_sGNE_Matrix.csv`
   - `result/step2_intermediate_results.mat`
3. `step3_plot_sGNE_DEG.m` 计算各分期的 global sGNE，并生成：
   - `result/STAD_sGNE_Progression.png`

### 路径配置说明

三个脚本会根据其在 `code/` 文件夹中的位置，自动确定 `data/` 和 `result/` 文件夹的路径。如果保持仓库原有的目录结构，只需修改上述 MATLAB 命令中的 `project_dir`。如果移动了脚本或改变了文件夹结构，则需要修改各脚本开头附近的“路径配置”部分。
