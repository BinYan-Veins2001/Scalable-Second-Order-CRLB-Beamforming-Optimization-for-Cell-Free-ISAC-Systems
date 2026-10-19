# COSTAR 新版论文绘图代码

## 公开版本与结果下载

公开版本为 `v2026.10.20`，采用 MIT 许可证。完整中文推导读本见 [`docs/COSTAR_Chinese_Companion.pdf`](docs/COSTAR_Chinese_Companion.pdf)，共 118 页。

大型 MATLAB 数据通过 [GitHub Releases](https://github.com/BinYan-Veins2001/Scalable-Second-Order-CRLB-Beamforming-Optimization-for-Cell-Free-ISAC-Systems/releases/tag/v2026.10.20) 分发。在 MATLAB 的仓库根目录运行：

```matlab
download_paper_results;                  % 下载并核对 SHA-256，拒绝覆盖已有图数据
addpath(genpath(pwd),'-begin');
main_replot_paper_figures(1:6);           % 只重新绘图，不运行优化
```

也可以手动将 Release 中的 `COSTAR_results_20261020.zip` 解压到仓库根目录。需要覆盖本地已有数据时才调用 `download_paper_results(true)`。该数据包包括图 1 的既有绘图数组、图 2–6 的完整结果及配对验证记录。图 1 的完整历史原始存档另列为可选附件，未用新版求解器重跑。数据与计时口径见 [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md)。

代码在 MATLAB R2022b 下验证。COSTAR 不依赖 CVX；涉及凸子问题的基准需另行安装 CVX 与其求解器，默认选择 SDPT3。第三方软件未打包进仓库。

本文件夹以 `COSTAR_Safeguard_Adaptive_20261019` 为基础建立独立副本。后续算法和绘图调整都在此目录内进行。原代码目录和论文不由本项目改写。

默认系统为 16 个 AP、每 AP 4 根发射和接收天线、8 个用户。COSTAR 采用当前论文中的单循环算法：协调工作集阈值和 NNLS 乘子、真实可行参考下降、精确 HVP 的 tCG、可行性恢复，以及确定性的 SPSC 曲率认证与负曲率方向恢复。默认不调用随机块 Krylov。

## 论文图编号

| 新编号 | 入口 | 内容 | 对应旧程序编号 |
|---|---|---|---|
| 1 | `main_fig01_network_dimension_scaling` | CRLB 和总运行时间随系统维度变化 | 15 |
| 2 | `main_fig02_costar_convergence` | CRLB、完整 KKT 残差和水平曲率轨迹 | 16 |
| 3 | `main_fig03_endpoint_curvature` | 七种算法的终点水平曲率 | 19 |
| 4 | `main_fig04_joint_stationarity` | 联合一阶、二阶停止精度 | 20 |
| 5 | `main_fig05_local_quadratic` | 受控 Newton 局部残差收缩 | 18 |
| 6 | `main_fig06_structured_curvature_efficiency` | SPSC 与稠密切空间特征值计算的效率 | 新增 |

旧论文图3的乘子与约束裕量图不列入新入口。`sim_twc/` 中保留的历史函数为兼容支持代码；复现正文请使用上表入口。

## 运行

在 MATLAB 中切换到本文件夹：

```matlab
addpath(genpath(pwd),'-begin');
main_run_paper_figures;       % 默认只计算图2--6，包含可恢复的数据存档
main_replot_paper_figures;    % 只读取已存数据并重新绘图
```

单独运行示例：

```matlab
main_fig02_costar_convergence;
main_fig04_joint_stationarity;
main_fig06_structured_curvature_efficiency;
```

PowerShell 中也可逐图启动独立的单线程 MATLAB 进程：

```powershell
./run_paper_figures.ps1
./run_paper_figures.ps1 -Figures 2,4,6
```

图1不会在默认批处理中计算。以后明确需要重新计算时，执行 `main_run_paper_figures(struct('figures',1,'run_dimension_scaling',true))`，或使用 PowerShell 的 `-Figures 1 -RunDimensionScaling`。

## 目录和数据口径

- `optimization/`、`config/`：COSTAR 算法和系统参数。
- `baselines/`：保持原来的 ZF-BF、SO-BF、PFO-ALM、PTR-ALM、MRFP-SCA 和 SDR-SCA 基准。
- `paper_figures/`：按新正文编号组织的计算、绘图及公共函数。
- `results/paper_figXX/`：每张图的独立数据、原始运行记录和 CSV。
- `figures/`：可编辑 MATLAB `.fig`、PNG 和矢量 PDF/EPS。
- `audit/`：机制测试、源码来源和性能一致性核对。
- `docs/`：算法与论文对应说明及实验结果说明。

颜色、标记、Times New Roman 字体和 IEEE 尺寸沿用原程序。图2用论文定义的完整 KKT 残差；图3使用水平空间的最小 Ritz 曲率估计；图4的 `q_stop` 沿用正文两项最大值定义，另存实际 SPSC 证书、临界锥精度及算法成功状态。

图5以新版 COSTAR 的默认系统终点为起点，采用平方阶约束精度、二次线性求解精度和预定的局部扰动，观察受控 Newton 步的残差收缩。全部合格样本均参与拟合，不按残差是否收缩进行筛选。局部正则条件另外保存：本次参考点不满足严格互补及水平 Hessian 正定条件，因此该图不作为局部距离二次收敛定理全部前提的验证。默认 COSTAR 的 tCG 参数保持原设置。

图6对每次求解返回的最终状态统一进行测量，使用相同乘子、约束法向和阈值，比较单次 SPSC 与显式切空间 Hessian 特征分解的计算时间。每个规模点使用3个固定种子，每例取3次预热后配对计时的中位数，曲线再取跨种子的中位数。优化器是否成功终止单独记录；本次24个独立实例中有1例达到迭代上限，该例照常计入内核比较。SPSC 返回的曲率下界与真实最小特征值含义不同，一致性通过相同阈值下的惯性判定核对。本批内核测量均采用全空间认证路径；切空间补充认证和负曲率恢复另由单元测试核验。

本次实际运行、验证结果和逐图状态以 `docs/实验结果与验证.md` 为准。
