# 药材干燥过程中的温度与水分浓度变化研究

2026 年数学建模 A 题参赛论文、MATLAB 源代码与结果表归档。

## 文件目录

| 文件或目录 | 内容 |
| --- | --- |
| `2026_A题.pdf` | 参赛论文，共 60 页，包含正文和代码附录 |
| `AI工具使用详情.pdf` | AI 工具使用说明，共 4 页 |
| `q1/` | 第一问：预热阶段径向传热与水分扩散 |
| `q2/` | 第二问：变物性热湿耦合，计算前 3 小时的状态 |
| `q3/` | 第三问：环境稳定窗口选择、长期延拓与全域干燥事件定位 |
| `q4/` | 第四问：实测半径收缩驱动的移动边界热湿耦合 |
| `result1.xlsx`～`result4.xlsx` | 对应四问的已提交结果表 |

每问包含 5 个 MATLAB 文件：`main_qN.m` 为入口，`qN_compute.m` 实现计算，`qN_outputs.m` 负责输出，`qN_validation.m` 包含检查逻辑，`qN_trace_store.m` 管理计算轨迹。

## 方法概述

将药材近似为轴对称长圆柱体，描述径向温度场与干基含水率场。代码采用控制体有限体积离散、后向 Euler 时间推进及阻尼 Newton 迭代。长期干燥问题利用空间重建后的全域最大含水率定位达标事件；收缩问题通过 `ξ = r/R(t)` 参考坐标和 PCHIP 半径插值处理移动边界。

## 运行入口与输入依赖

代码使用 MATLAB 的 `arguments`、字符串数组、表格读写、Java 和 HDF5 等功能，应在具有相应函数支持及 JVM 的 MATLAB 环境中运行。此归档没有锁定或重新验证 MATLAB 版本。

**当前归档包含提交材料，但缺少部分原始输入、模板和续算状态，不能直接完成全量复现。** 以下路径相对于各问目录；根目录的 `result*.xlsx` 是已提交结果，不是原始输入或空白模板。

| 问题 | 代码引用的主要外部文件（当前目录未包含） | 入口与模式 |
| --- | --- | --- |
| 第一问 | `q1/data/raw/附件1.xlsx`、`q1/data/raw/附件3/result1.xlsx` | `main_q1("mvp")`、`main_q1("reproduce")`、`main_q1("validation")` |
| 第二问 | `q2/data/raw/附件1.xlsx`、`q2/data/raw/附件3/result2.xlsx` | `main_q2("mvp")`、`main_q2("reproduce")`、`main_q2("validation")` |
| 第三问 | `q3/data/raw/附件1.xlsx`、`q3/data/q3_start_4h.mat`；`finish`/`validation` 还读取 `q3/附件/未舍入结果.mat` | `main_q3("reproduce")`、`main_q3("finish")`、`main_q3("validation")` |
| 第四问 | `q4/data/raw/附件1.xlsx`、`q4/data/raw/附件2.xlsx`；导出结果另需 `q4/data/raw/附件3/result4.xlsx` | `main_q4("check")`、`main_q4("result")`、`main_q4("export")` |

补齐对应文件后，在 MATLAB 中切换到相应问的目录。例如，从仓库根目录运行第一问：

```matlab
cd q1
output = main_q1("reproduce");
```

第二问的 `reproduce` 模式返回计算场与检查信息，不写出结果文件。第四问的 `result` 模式计算并返回结果，`export` 模式另写出文件。运行可能生成较大的轨迹和缓存，并受到代码内部时间、内存及磁盘预算限制。


