# 新增侧链 Gaussian 输入（待审核，未提交）

## 状态

本目录只生成三套 Gaussian 输入及其审计结构，尚未生成 Slurm 作业，也没有提交任何计算。

三个模型均使用 monosubstituted thiophene 作为量化封端，保持与历史 `C46_beta_thiophene.gjf` 相同的片段边界思想。Gaussian 优化完成后，ztop 将切断 thiophene–侧链键并保留侧链侧。

## 模型

- `Cb1T`：thiophene–CH2–cyclobutyl，服务于 L8₄-Cb 的短外侧链。
- `Cb3T`：thiophene–(CH2)3–cyclobutyl，服务于 L8-Cb₄。
- `C46T_sp3`：thiophene–(CH2)3–CH(C4H9)(C6H13)，服务于 L8₄-Cb 的长 spacer 分叉链。

初始构象由最多 256 个 ETKDG 构象分别进行 MMFF94 清理，然后在最低能量以上 25 kcal/mol 窗口内选择空间伸展最大的无硬碰撞构象。对 `C46T_sp3` 额外将可旋转烷基键调整到接近全反式，并在 175–180° 二面角约束下进行 MMFF94 局部清理；分叉中心仍保持正常四面体键角。Gaussian `opt` 本身不带这些约束。这些处理只为提供不折叠的 Gaussian 起始几何，不把 MMFF94 能量当成最终科学结果。

## Gaussian 层级

与历史片段保持一致：

```text
#p opt pbe1pbe def2svp emp(gd3bj) scf(maxcycle=64) iop(5/13=1) scrf(solvent=diphenylether) geom(connectivity)
```

共同设置：`charge=0`、`multiplicity=1`、`%nproc=32`、`%mem=100GB`。

## 审核重点

1. 是否接受三个片段的化学边界与 thiophene 封端；
2. 是否接受只做 `opt`、不在本轮同时做 `freq`；
3. 是否接受 PBE1PBE/def2-SVP/GD3BJ/diphenyl ether 与历史 L8BO 保持一致；
4. 是否接受每个输入 32 CPU、100 GB 内存的资源声明；
5. 初始长链构象是否足够伸展。

`input_report.json` 保存分子式、原子数、构象筛选指标、ztop 切割键和每个 GJF 的 SHA256。
