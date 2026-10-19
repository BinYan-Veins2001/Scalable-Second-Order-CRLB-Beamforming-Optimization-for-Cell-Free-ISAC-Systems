# COSTAR 与最新版论文的对应关系

本目录采用 `COSTAR_Safeguard_Adaptive_20261019` 中已经验证的自适应 NNLS 保护版本。核对对象为 `论文/COSTAR_TWC_KKT_SOSP_20261018/paper/sections/03_costar.tex` 和 `sections/04_analysis_final.tex`。第四章已经纳入该版本的保护步、恢复与接受规则，因此本次无需重新设计数值算法。`audit/alignment_sources.json` 保存核对时论文及冻结物理输入的 SHA-256。

## 默认求解流程

| 论文步骤 | 本目录中的实现 | 核查结果 |
|---|---|---|
| 共享可行 stacked-ZF 初值 | `initialize_shared_stacked_zf.m`；`solve_costar.m` 必须接收外部初值 | 保留；各迭代基准使用同一初值 |
| 工作集、NNLS 乘子、完整 KKT 残差 | `costar_paper_geometry.m` | SINR 与 AP 阈值分别取上限和共同阈值的最小值；残差同时包含可行性、完整 NNLS 驻点残差与互补残差 |
| 阈值与一阶精度协调 | `costar_paper_geometry.m` | 每轮重置共同阈值；不满足阈值条件或最小范数向内子问题时，减半阈值并重新计算 |
| NNLS 参考下降量 | `costar_paper_rtr.m` 中 `safeguard` | 每次一阶迭代实际求出可行参考点；原有二次向内修正与 Armijo 下降条件保留 |
| 自适应保护步初始长度 | `guard_initial_scaling=true` | 对应论文自适应起点公式；修正项的相对大小由 `guard_correction_fraction=0.25` 控制 |
| 保护步低成本求值 | `costar_guard_path_cache.m`、`costar_guard_path_evaluate.m` | 沿二次曲线缓存四次的约束/感知统计量；缓存只筛选候选，最终通过直接求值接受；缓存搜索耗尽时仍有直接搜索 |
| tCG | `trust_region_tcg.m`、`exact_lagrangian_hess_vec.m` | 精确拉格朗日 HVP；保留非正曲率分支；不构造环境空间完整 Hessian |
| 一阶恢复 | `accelerated_restore`，失败后 `ordinary_restore` | AP 回缩、随阈值缩小的 SINR 向内裕量、二阶大小的法向修正；必要时使用满足同一修正上界的法向 Newton 候选 |
| 一阶候选接受 | `costar_paper_rtr.m` | 可行、正预测下降、法向修正为步长平方阶，且实际下降同时不小于模型比例与参考下降比例；该分支不额外强制曲率分支的加权模型缺陷条件 |
| 一阶退回参考点与半径更新 | `costar_paper_rtr.m` | 质量条件失败则收缩；只有参考下降比较失败时保持原半径 |
| 确定性 SPSC | `costar_sosp_structured_psd.m` | 先检查全空间；必要时利用约束法空间计算切空间惯性，并从相同分解恢复负曲率方向 |
| 负曲率步 | `costar_paper_rtr.m` 曲率分支 | 重置半径，严格模型缺陷与比例接受；普通恢复失败后可用带三次向内裕量的法向 Newton 修正 |
| 一阶/二阶精度协调 | `calibrate` 与曲率分支 | 达到当前完整 KKT 精度才调用 SPSC；全空间认证或数值临界锥校准达到目标精度后才成功；每次未终止曲率调用后细化一阶精度 |
| 终点一致性 | `costar_paper_rtr.m` 终点重算 | 重新计算最终点的几何和乘子；步前旧曲率结果不能冒充步后证书；预算耗尽不报告成功 |

默认 `guard_directional=false`、`guard_scaled_inward=false`、`ambient_endpoint_probe=false`。历史实验选项未进入默认 COSTAR；随机块 Krylov 在此入口中的调用次数为零。内部优化目标为 `log_crlb_ratio`，图中物理 CRLB 由同一感知模型还原。

## 性能与复杂度

本次核对未改变上述算法的数值步骤。60 个 `optimization`/`config` 源文件中，仅更正一处过时求解器说明、配置说明和版本元数据；具体清单见 `audit/alignment_sources.json`。保护曲线缓存测试改为在本目录内生成物理问题，不再依赖旧项目路径。

主要复杂度阶数继承已验证版本：单次精确 HVP 为 `O(UN^2+NU(U+a))`，工作集 Gram 运算为 `O(da^2)` 并配合小规模求解；缓存建立需要一次矩阵运算，后续保护步筛选为 `O(U+M)`。SPSC 分解 `U` 个用户块，主导项为 `O(UN^3)`，切空间补充操作使用至多 `3+a` 阶矩阵。这里 `N=MN_t`、`d=2NU`、`a` 为当前独立工作约束数。内层次数和调用次数仍影响总耗时，不能仅由单次 HVP 阶数推出总时间相同。

本目录新运行的三组默认规模配对试验全部完成，随后对每个旧/新实例各增加两次全新进程重复，共保留 18 次运行。按每个实例的三次时间取中位数，再在三个种子间平均，旧版为 `1.99664 s`、新版为 `2.94746 s`，均值之比为 `1.47621`；三个配对时间比的平均为 `1.51043`。CRLB 最大相对差为 `1.97983e-10`。旧版与新版各 9/9 通过共同一阶/二阶终点审计，新版也 9/9 通过自身乘子与工作集审计。所有重复运行的工作计数、最终波束和数值迭代轨迹均逐项一致。

新版最大完整 KKT 残差为 `6.44569e-6`，最小切空间特征值为 `-1.64955e-6`，算法采用的最大临界锥曲率容差为 `7.5e-6`，目标均为 `1e-5`。中位数汇总见 `audit/alignment_repetitions.json` 和 `benchmark/validation/paired_repetition_medians.csv`；全部重复原始记录见 `benchmark/validation/trial_results_all_repetitions.csv`。初次三组数据仍保存在 `audit/alignment_pairs.json`，未删除或替换。

种子 2027 的新版中位数为 `4.174 s`，旧版为 `1.948 s`。其新版工作量在历史与当前运行中完全一致：818 次 HVP、286 次 NNLS、178 次阈值细化和 592 次恢复迭代。历史新版时间为 `3.892 s`，因此该实例的较高开销主要来自保护和恢复次数，并非重复计时中突然出现的数值流程变化；详细对照见 `audit/alignment_timing_investigation.json`。

上一轮已验证的 8 组冻结配对实例均通过共同终点检查；其中默认规模的 5 组历史平均优化时间为旧版 1.971 s、新版 2.369 s，CRLB 最大相对差在 `1.8e-9` 内。该历史批次与本目录新三组时间分别记录。两个批次都说明原有 CRLB 精度得到保留，主要计算阶数保持不变；参考下降与精度协调带来的实际开销依赖实例，不能表述为运行时间完全相同。

## 验证与复现

`test_latest_alignment` 已通过五套检查：几何/NNLS，16 组独立稠密 SPSC 惯性比较，解析 KKT 与负曲率恢复及旧证书拒绝，四次缓存与直接值比较，以及小规模物理问题的独立完整 Hessian 终点核验。结果见 `audit/alignment_mechanisms.json` 和同名日志。

新三组默认规模配对输入固定为 16 AP、每 AP 4 发射/4 接收天线、8 用户，种子 2026–2028。每组旧/新算法接收位级相同的信道、感知矩阵与初始波束；新版本只加载当前算法默认配置。各算法在全新单线程 MATLAB 进程中顺序运行，计时不含初始化和独立终点审计。初次运行文件列于 `benchmark/validation/jobs.txt`，后续完整重复列于 `jobs_repetitions.txt`，全部 18 个输入列于 `jobs_all_repetitions.txt`。初次汇总函数为 `summarize_alignment_validation`；含初次测量的重复汇总函数为 `summarize_alignment_repetitions`。

共同固定阈值审计用于公平比较旧/新终点；另用新算法自身乘子和工作集进行独立审计，并检查算法报告的最终点证书与临界锥精度。二者分别保存，不能将共同审计重估的乘子替代算法自身的精度声明。独立完整 Hessian 仅用于验证，不计入 COSTAR 优化时间。

代码的有限预算和浮点容差对应数值实现。论文的一阶 KKT 收敛及有限精度二阶终止分别采用文中 MFCQ 和相应 LICQ、严格互补等条件；默认 tCG 相对残差规则也不单独构成局部二次收敛保证。
