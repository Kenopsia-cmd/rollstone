# 2D 版本

本目录只放二维版本。最新二维源码是 [rollingstone_fixed.cpp](07-random-terrain-v2/rollingstone_fixed.cpp)。

| 顺序 | 源码 | 主要变化 |
| --- | --- | --- |
| 01 | [原始版](01-original/rollingstone.cpp) | 最初提供的滚石程序，与 v1.0 源码相同。 |
| 02 | [前三项任务](02-tasks-1-3/rollingstone.cpp) | 随机批量模拟、坡脚概率、CSV、参数研究及交互画布。 |
| 03 | [全部任务早期版](03-all-tasks/rollingstone.cpp) | 地形研究、凸多边形、刚性/弹性挡板、回放与能量显示。 |
| 04 | [渲染优化副本](04-render-optimized/rollingstone.cpp) | 字体缓存、轨迹抽稀与研究表校验。 |
| 05 | [字体清晰版](05-font-clarity/rollingstone.cpp) | 字号、抗锯齿、按钮居中及文字截断优化。 |
| 06 | [随机地形 v1](06-random-terrain-v1/rollingstone.cpp) | 数据总览、可复现随机山坡，保持整体下坡与总宽高差。 |
| 07 | [随机地形 v2](07-random-terrain-v2/rollingstone_fixed.cpp) | 陡缓成组交替、主要坡折 3–5 处，改善过于平直的问题。 |

05 和 07 下的 `variants` 保留仅文件头注释不同的副本，功能相同。以上顺序表示可确认的功能阶段，不是补建的历史提交。

原 Git 的 12 个标签完整工程见 [标签快照](../archive/2d/git-tags/)，恢复 Git 历史见 [备份说明](../history/README.md)。
