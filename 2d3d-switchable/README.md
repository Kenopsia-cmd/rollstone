# 2D / 3D 可切换版

两个整合版本分别保留，不覆盖彼此的差异。

| 版本 | 入口 | 切换方式 | 区别 |
| --- | --- | --- | --- |
| 替代整合副本 | [rollingstone_fixed_copy.cpp](01-alternative-copy/rollingstone_fixed_copy.cpp) | 按 **3** 或 `--3d` 进入三维；Esc 返回二维 | 二维功能加地形 CSV 导入、逐段材质和渲染缓存；三维部分是球体演示。 |
| **F3 整合版（推荐）** | [源码和 VS 工程](02-f3-integrated/) | **F3** 或右上角按钮 | 两模式分别保留设置和结果；三维支持真实横向运动、球/凸棱柱、有限宽刚性/弹性挡板、横向×高度概率、CSV 和参数/地形研究。 |

F3 版：打开 `02-f3-integrated/rollingstone_2d3d.vcxproj`，选择 Release / x64，Ctrl+F5。依赖配置和工具集说明见[仓库首页](../README.md)。详细按键、统计口径和命令行示例见[使用说明](02-f3-integrated/使用说明.txt)。

F3 版三维视角：右键拖动旋转，滚轮缩放；B 批量，E 导出，A 参数研究，G 地形研究，N 生成随机坡，D 查看断面概率。

独立运行的另两份三维源码在 [3d-standalone](../3d-standalone/README.md)，二维历史在 [2d](../2d/README.md)。
