# 山区滚石模拟

Windows C++ 单文件项目，可在二维和三维之间切换。

| 文件 | 用途 |
| --- | --- |
| [rollingstone_3d_controls.exe](bin/rollingstone_3d_controls.exe?raw=true) | 下载后直接运行，Windows x64 |
| [rollingstone_fixed_copy.cpp](rollingstone_fixed_copy.cpp) | 完整 C++ 源码 |
| [rollingstone_3d_controls.vcxproj](rollingstone_3d_controls.vcxproj) | Visual Studio 工程 |

## 操作

- **F3** 或点击右上角按钮切换二维和三维。三维中也可按 **Esc** 返回二维。
- 在三维左侧调整初速、地表系数或地形参数，点「应用并重跑」或按 **Enter** 生效。
- 在「起步」页设定沿坡位置 **x**，点「在起点横线上选」设置横向位置 **y**，也可直接输入坐标。
- 在「拦挡」页选择刚性或弹性挡板。可调位置、宽度、高度，弹性板还可调刚度、阻尼和最大位移。
- 顶部的「轨迹截面」显示水平 **X–Y** 和竖直 **X–Z** 轨迹。点「3D ＋截面」可和立体画面一起看。按 **D** 可循环切换视图。

二维模式有批量随机模拟、坡脚概率、CSV 导出、参数和地形研究，以及凸多边形模拟。三维模式使用球体，支持横向滚动、单次模拟和挡板观察。

## 在 Visual Studio 中编译

下载整个仓库并解压，打开工程文件。安装 Visual Studio 的「使用 C++ 的桌面开发」工作负载和 EasyX，选择 **Release / x64**，按 **Ctrl+F5**。工程使用 C++20 和 v145 工具集。使用其他工具集时，在项目属性中选择已安装的版本，并确保 EasyX 的头文件和库已配置。
