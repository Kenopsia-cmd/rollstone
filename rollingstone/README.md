# 山区滚石模拟

基于用户提供的米制滚石程序扩展。原始版本保存在 Git 标签 `v1.0`，开发版每完成截图中一项要求就提交一次，并创建对应标签。

## 运行

双击 `run.cmd` 打开动画；双击 `build.cmd` 重新编译。需要 Visual Studio C++ 与 EasyX，本机已经安装。附带的旧 `graphics.h` 是 WinBGI 头文件，构建脚本明确使用本机的 EasyX 头文件和库，不需要另行下载它。

在本目录打开 PowerShell：

```powershell
.\build\rollingstone.exe --batch --runs 200 --seed 42
```

随机种子相同、参数相同、可执行程序相同即可重现实验。每次试验独立抽样，单条轨迹内参数保持不变。当前阶段的无参数启动保留原动画界面；最终版本将统一动画和批量模型。

## 统计口径

检测器默认位于坡脚 x=120m 的竖直断面，统计石心首次从左向右通过的高度（相对于坡脚高程），分成 0.5m 高度区间，范围 0–10m，最后一个区间包含 10m。每次试验最多计数一次。概率分母为全部试验数，未到达、超过10m、低于0m单独报告。路线网格按每条轨迹是否访问某单元格计数，单条轨迹在同一格只计一次，各格概率之和不要求为1。

原程序到坡脚直接钳位终止，且接触修正会丢失部分碰撞反弹。批量核心改为半隐式积分和接触冲量求解，采用实心球惯量和库仑摩擦，坡脚之后水平延伸；这会与 v1.0 的数值结果有所不同。

计算域在坡脚后80m结束。`domain_exit` 代表轨迹离开计算域，此时行程是下限，不是真实停止距离；`timeout` 也不代表已经停下。

## 当前版本

v1.1：随机批量模拟、坡脚0–10m高度概率、路线网格概率、固定种子、无图形计算、有效性检查。

这是二维教学近似模型，尚未经实测标定，不能用于工程拦挡设计或安全认证。

## v1.2 CSV输出
运行 --batch 后自动生成 results 下带时间戳的目录；可用 --out 指定新的空目录。导出 height_probability.csv、route_probability.csv、runs.csv、trajectories.csv、sampled_surfaces.csv、terrain.csv、config.csv 及 summary.txt。列名含SI单位，未到达断面的值留空，不伪造为零。输出目录已有内容时拒绝覆盖。

## v1.3 参数研究
命令：.\build\rollingstone.exe --study --runs 200 --seed 42。初速和地表四项系数按正态分布抽样，sigma=均值绝对值×0.1，非物理取值重采样；每条轨迹内参数固定。按均值0.8/1/1.2倍做单因素配对比较，输出parameter_study.csv和从实际结果生成的分析结论。所有核心功能均为C++，CMD仅编译/启动。

## v1.4 地形与地表研究
命令：.\build\rollingstone.exe --terrain-study --runs 200 --seed 42。比较横向尺度、高度尺度、局部缓台/坡度、岩面/土面/草地和分层地表；输出terrain_study.csv及实际分析。--width-scale/--height-scale可调基准几何；--terrain-file可输入x_m,z_m控制点CSV。详见STUDIES.md。
