# 独立 3D 版

这一系列是单独运行的三维程序，不在程序内切换到原二维画布。

| 版本 | 文件 | 功能 |
| --- | --- | --- |
| 三维原型 | [rollingstone3d.cpp](01-prototype/rollingstone3d.cpp) | 三维高度场、球体侧向滚动、沟谷、随机山坡、视角旋转、轨迹回放。 |
| 三维统计版 | [rollingstone3d_full.cpp](02-full/rollingstone3d_full.cpp) | 在原型基础上增加批量试验、高度/横向概率、CSV、参数研究、地形研究、刚性/弹性拦挡。 |

**这两个版本的三维碰撞形状只有球体。** 需要球/凸棱柱、F3 切换及完整三维断面统计，请使用 [2D/3D F3 整合版](../2d3d-switchable/02-f3-integrated/)。

使用 Windows + Visual Studio + EasyX 编译，应用代码均在所选单个 `.cpp` 内。操作按键和编译命令见各文件头部注释。保留原代码和名称，未将历史文件重写为最新实现。
