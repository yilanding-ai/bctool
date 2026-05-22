# BCTool 完全教程

## 目录

- [1. 简介](#1-简介)
- [2. 安装](#2-安装)
- [3. 快速开始](#3-快速开始)
- [4. ConvSeed 内置比对器](#4-convseed-内置比对器)
- [5. 所有 12 种转化类型](#5-所有-12-种转化类型)
- [6. 使用外部比对工具](#6-使用外部比对工具)
- [7. 错误校正](#7-错误校正)
- [8. 基准测试](#8-基准测试)
- [9. 结果解读](#9-结果解读)
- [10. 常见问题](#10-常见问题)

---

## 1. 简介

BCTool（Bisulfite Conversion Alignment & Benchmarking Toolkit）是一个**甲基化/转化测序比对与评测工具包**。核心特性：

- **ConvSeed**：内置的转化感知比对器，无需安装任何外部工具即可运行
- **12 种转化类型**：C>T（重亚硫酸盐）、A>G（m6A）、G>A 等所有单碱基替换
- **15 种外部比对器**：Bwa-meth、BSBolt、BSMAP、Walt、Abismal、HISAT-3n、Bismark、BSseeker2 等
- **7 种错误校正策略**：MAPQ 过滤、soft-clip 过滤、成对比对一致性、上下文贝叶斯校正等
- **完整评测体系**：精确率、召回率、F1 分数、甲基化水平相关性

### 适用场景

| 场景 | 说明 |
|------|------|
| 重亚硫酸盐测序 (WGBS) | C>T 转化，检测 5mC |
| m6A 直接 RNA 测序 | A>G 转化，检测 m6A |
| 其他修饰检测 | 任何单碱基修饰（C>A, G>T 等） |
| 比对器选型 | 在多种比对器间做 Benchmark |
| 参数调优 | 测试不同种子长度、错误校正策略 |

---

## 2. 安装

### 从源码安装

```bash
git clone https://github.com/yilanding-ai/bctool.git
cd bctool
pip install -r requirements.txt
pip install .
```

### 验证安装

```bash
bctool --help
```

输出应显示所有可用命令。

### Docker

```bash
docker build -t bctool .
docker run --rm bctool --help
```

---

## 3. 快速开始

### 3.1 最简单的 Demo（无需外部工具）

```bash
# ConvSeed 内置比对器
bctool run --simulate --aligners conversion_aware -o ./results
```

这将在模拟数据上运行 ConvSeed，生成 HTML 报告。

### 3.2 包含所有内置比对器

```bash
# ConvSeed + 4 个 Mock 比对器
bctool run --simulate -o ./results
```

### 3.3 指定转化类型和读长

```bash
# C>T 转化，模拟 100bp reads
bctool run --simulate --conversion ct --sim-reads 50000 --read-length 100 -o ./ct_results

# A>G 转化（m6A）
bctool run --simulate --conversion ag --sim-reads 50000 -o ./ag_results
```

### 3.4 使用真实数据

```bash
bctool run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -r reference.fa --conversion ct -o ./real_results
```

---

## 4. ConvSeed 内置比对器

ConvSeed 是 BCTool 内置的转化感知比对器，注册名为 `conversion_aware`。

### 4.1 工作原理

ConvSeed 结合了 3 种技术：

```
┌─────────────┐    ┌──────────────┐    ┌─────────────────┐
│ 3-base 编码  │───▶│  Rolling Hash │───▶│ Banded SW 比对   │
│ (转化感知)   │    │  (种子查找)   │    │ (affine gap)    │
└─────────────┘    └──────────────┘    └─────────────────┘
```

1. **3-base 编码**：将 4 种碱基映射到 3 种状态，转化对（如 C/T）共享同一状态
2. **种子哈希**：在参考基因组上构建滚动哈希索引，用于快速定位
3. **带状 SW**：在候选位置附近用 Smith-Waterman 算法做精细比对

### 4.2 关键参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `seed_length` | 16（自动适配） | 种子长度，越长越精确但敏感度降低 |
| `min_seed_matches` | 2 | 最少匹配种子数，用于过滤假阳性 |
| `band_width_ratio` | 2.0 | SW 带状宽度 = read_length × ratio |
| `max_candidates` | 10 | 最多候选位置数 |
| `auto_adapt` | True | 是否根据读长自动调整参数 |

### 4.3 读长自适应逻辑

```
读长 < 40bp:  k=8
读长 40-60:  k=10
读长 60-80:  k=12
读长 80-120: k=14
读长 ≥120:   k=16
```

短序列自动使用更短种子以提高命中率，同时缩小候选数以减少假阳性。

### 4.4 短序列优化（40-80bp）

对于 40-80bp 的短序列，ConvSeed 自动启用优化模式：

- 种子步长从 `k/2` 缩小为 `read_len/10`（更多种子锚点）
- 最大候选数从 10 降为 5（减少误比对）
- 种子长度自适应调节（k=10-14）

实际性能（5000 reads, 50bp）：

| 指标 | C>T | A>G | G>A | T>C |
|------|-----|-----|-----|-----|
| F1 | 0.875 | 0.870 | 0.889 | 0.899 |
| Precision | 0.985 | 0.976 | 0.985 | 0.985 |
| Recall | 0.787 | 0.784 | 0.810 | 0.826 |
| 时间 | 1.6s | 1.5s | 1.7s | 1.5s |

### 4.5 手动调优示例

```python
from bctool.aligners.conversion_aware import ConversionAwareAligner

# 针对 40bp reads 手动调优
aligner = ConversionAwareAligner(
    config, output_dir,
    seed_length=10,       # 更短的种子
    min_seed_matches=3,   # 更严格的种子过滤
    max_candidates=5,     # 更少的候选
    auto_adapt=False,
)
```

---

## 5. 所有 12 种转化类型

| 标签 | 转化 | 描述 | 典型应用 |
|------|------|------|---------|
| ct | C→T | C 转化为 T | 重亚硫酸盐测序 (WGBS), 5mC 检测 |
| tc | T→C | T 转化为 C | 反向重亚硫酸盐测序 |
| ag | A→G | A 转化为 G | m6A 直接 RNA 测序 |
| ga | G→A | G 转化为 A | 反向 A>G 转化 |
| ac | A→C | A 转化为 C | 特定修饰检测 |
| ca | C→A | C 转化为 A | 反向 A>C 转化 |
| gt | G→T | G 转化为 T | 特定氧化修饰 |
| tg | T→G | T 转化为 G | 反向 G>T 转化 |
| at | A→T | A 转化为 T | 脱氨修饰 |
| ta | T→A | T 转化为 A | 反向 A>T 转化 |
| cg | C→G | C 转化为 G | 特定烷基化修饰 |
| gc | G→C | G 转化为 C | 反向 C>G 转化 |

### 批量运行所有转化类型

```bash
for conv in ct tc ag ga ac ca gt tg at ta cg gc; do
  bctool run --simulate --conversion $conv --aligners conversion_aware -o ./${conv}_results
done
```

---

## 6. 使用外部比对工具

BCTool 支持 15 种外部比对器。ConvSeed 是**唯一内置**的，其他需要自行安装。

### 支持的比对器列表

```bash
bctool list-aligners
```

### 安装示例（Linux/macOS）

```bash
# Conda 安装
conda install -c bioconda bwameth bsbolt bsmap walt abismal hisat-3n bismark

# 运行特定比对器
bctool run --simulate --aligners bwameth,bsbolt,hisat3n -o ./compare_results
```

### 添加自定义比对器

只需继承 `AlignerBase` 并实现三个方法：

```python
from bctool.aligners.base import AlignerBase

class MyAligner(AlignerBase):
    name = "my_aligner"
    binary_name = "my-aligner"
    requires_index = True

    def build_index(self, reference_fa):
        # 构建索引
        pass

    def run_align(self, read1, read2, reference_fa, threads=8):
        # 运行比对，返回 SAM/BAM 路径
        pass

    def call_methylation(self, sam_bam_path, reference_fa):
        # 甲基化检出，返回 BED 路径
        pass

# 注册
from bctool.aligners import ALIGNER_REGISTRY
ALIGNER_REGISTRY["my_aligner"] = MyAligner
```

---

## 7. 错误校正

BCTool 有 7 种错误校正策略，分为两个阶段：

### 7.1 SAM 级滤波器（Phase 1.5，甲基化检出前）

| 策略 | 参数 | 说明 |
|------|------|------|
| `mq` | `--correction-min-mq 20` | 去除低 MAPQ 的 reads |
| `clip` | `max_clip_pct: 50` | 去除过多 soft-clip 的 reads |
| `pair` | — | 去除不协调的 PE reads（TLEN>1000 / 不同染色体） |
| `unconverted` | `--correction-max-unconverted 3` | 去除过多未转化目标碱基的 reads |
| `model` | `--correction-model-threshold 0.5` | 逻辑回归模型（9 个 SAM 特征） |

### 7.2 BED 级校正（Phase 3，甲基化检出后）

| 策略 | 说明 | 适用读长 |
|------|------|---------|
| `consensus` | 硬阈值：ratio≥0.7→1.0, ≤0.3→0.0 | 所有读长 |
| `context` | Beta-Binomial 贝叶斯后验 + 基序先验 | 根据读长自动选择参数 |

### 7.3 使用示例

```bash
# 基础错误校正
bctool run --simulate --error-correct --correction-strategies mq,clip,consensus -o ./corrected

# 完整管道
bctool run --simulate --error-correct \
  --correction-strategies mq,clip,unconverted,context \
  --correction-min-mq 30 \
  --correction-max-unconverted 5 \
  -o ./full_correction

# 上下文感知贝叶斯校正（自动检测基序）
bctool run --simulate --error-correct --correction-strategies context -o ./context
```

### 7.4 基序自动检测

当启用 `context` 策略时，BCTool 自动：

1. 从参考基因组提取甲基化和非甲基化位点上下游序列
2. 计算每个 k-mer 的富集比
3. 将富集基序作为贝叶斯先验

---

## 8. 基准测试

### 8.1 运行内置 Benchmark

```bash
cd benchmarks
python benchmark_convseed.py
```

这将测试 ConvSeed 在 4 种转化类型 × 3 种读长下的表现。

### 8.2 测试参数

当前 Benchmark 配置：
- 基因组大小：50kbp
- 每条 reads 数：5000
- 读长：50、100、150bp
- 转化类型：C>T、A>G、G>A、T>C
- 种子长度：10 和 12

### 8.3 解读结果

输出包含：
- **F1**：精确率和召回率的调和平均
- **Precision**：检出的甲基化位点中正确的比例
- **Recall**：真实甲基化位点中被检出的比例
- **Overlap**：比对于参考基因组的目标位点数
- **Align(s)**：比对耗时

### 8.4 性能预期（ConvSeed seed=12）

| 读长 | F1 | Precision | Recall | 速度 (5k reads) |
|------|-----|-----------|--------|-----------------|
| 50bp | 0.87-0.90 | 0.97-0.99 | 0.79-0.83 | ~1.5s |
| 100bp | 0.95-0.97 | 0.99 | 0.92-0.94 | ~4s |
| 150bp | 0.97-0.99 | 0.99+ | 0.95-0.98 | ~9-11s |

> **注意**：这些数据来自 5000 条 reads 模拟数据，实际性能可能因基因组大小、数据质量和参数设置而异。

---

## 9. 结果解读

### 9.1 输出目录结构

```
results/
├── comparison_results.json       # 完整结果 (JSON)
├── comparison_results.csv        # 完整结果 (CSV)
├── report/
│   ├── report.html               # HTML 报告（交互式）
│   ├── summary.txt               # 文本摘要
│   ├── f1_score.png              # F1 柱状图
│   ├── precision.png             # 精确率柱状图
│   ├── recall.png                # 召回率柱状图
│   ├── accuracy.png              # 准确率柱状图
│   └── heatmap.png               # 性能热力图
├── qc/                           # QC 报告（如果启用）
└── work/
    ├── simulated/                # 模拟数据（--simulate 时）
    └── conversion_aware/         # ConvSeed 结果
```

### 9.2 关键指标说明

- **F1 Score**: 综合精确率和召回率的指标，0-1 之间，越高越好
- **Precision**: 避免假阳性（FP）的能力，高 Precision 意味着少报
- **Recall**: 避免假阴性（FN）的能力，高 Recall 意味着不漏报
- **Accuracy**: 总体正确率，包括甲基化和非甲基化位点
- **Level Correlation**: 甲基化水平与真实值的皮尔逊相关系数

### 9.3 Web UI

```bash
bctool web
```

打开交互式 Web 界面，包含 5 个页面：
- Pipeline：运行分析
- QC：质量控制
- Results：浏览结果
- Multi-Sample：多样本对比
- About：版本信息

---

## 10. 常见问题

### Q: ConvSeed 和外部比对器有什么区别？

ConvSeed 是内置的，无需安装任何依赖。外部比对器（如 Bwa-meth、BSMAP）需要额外安装，通常在大基因组上更高效。

### Q: 如何处理真实测序数据？

```bash
bctool run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -r reference.fa --conversion ct --threads 16 \
  --enable-qc --enable-trim -o ./real_results
```

### Q: F1 很低怎么办？

1. 增加 `--sim-reads` 提高覆盖度
2. 尝试不同的 `seed_length`
3. 启用错误校正 `--error-correct`
4. 检查模拟数据的 `--error-rate`（默认 0.01 可能偏高）

### Q: 如何对比多个比对器？

```bash
bctool run --simulate --aligners conversion_aware,bwameth,bsbolt -o ./compare
```

或使用 Batch 模式：

```bash
bctool batch-demo
```

### Q: 窗口环境下 HISAT-3n 无法安装？

HISAT-3n 主要支持 Linux。Windows 用户可以使用 ConvSeed 作为替代，或通过 WSL/Docker 运行。

### Q: ConvSeed 为什么有时比 Mock 好、有时差？

Mock 比对器模拟的是**整体甲基化检出准确率**（对全基因组随机加噪），而 ConvSeed 做的是**真正的 read 比对**。ConvSeed 在覆盖度充足的位置精度更高（Precision 0.97-0.99），但可能因覆盖度不足导致 Recall 偏低。

---

## 附录：API 参考

### ConvSeed Aligner

```python
ConversionAwareAligner(
    config,              # Config 对象
    output_dir,          # 输出目录
    name_suffix="",      # 名称后缀（可选）
    seed_length=16,      # 种子长度
    min_seed_matches=2,  # 最小种子匹配数
    match_score=1,       # SW 匹配得分
    conversion_score=0,  # SW 转化匹配得分
    mismatch_penalty=-1, # SW 错配罚分
    gap_open=-4,         # SW 空位开放罚分
    gap_extend=-1,       # SW 空位延伸罚分
    band_width_ratio=2.0,# SW 带状宽度比
    max_candidates=10,   # 最大候选位置
    auto_adapt=True,     # 自动根据读长调参
)
```

### CLI 命令速查

```bash
bctool run          # 运行完整评测管道
bctool simulate     # 仅生成模拟数据
bctool extract      # 从 FASTQ 提取子集
bctool trim         # 接头修剪
bctool umi          # UMI 处理
bctool qc           # 质量控制
bctool report       # 生成报告
bctool batch        # 多样本批处理
bctool web          # 启动 Web 界面
```
