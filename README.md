# Rollingstone 本地版本仓库

当前工作树是 **tasks-1-3** 分支，只交付任务1–3，最新版本 **v1.3.1**。所有模拟、统计、CSV导出和分析代码均为C++。

打开 `E:\my-first-repo-dev\rollingstone`，双击 `run.cmd`，在中文菜单选择随机模拟/CSV导出或参数分析。详细说明见 [rollingstone/README.md](rollingstone/README.md)，示例结果见 [rollingstone/RESULTS.md](rollingstone/RESULTS.md)。

| 位置或标签 | 用途 |
| --- | --- |
| E:\my-first-repo\rollingstone | 原始第一版，main / v1.0 |
| E:\my-first-repo-dev\rollingstone | 当前前三项开发版，tasks-1-3 / v1.3.1 |
| v1.1 | 任务1：随机路线概率 |
| v1.2 | 任务2：CSV导出 |
| v1.3 | 任务3：正态参数影响分析 |
| v1.3.1 | 前三项的边界修正、中文菜单与验证 |

原始资料 `E:\资料\11.11\rollingstone` 保持独立。已先行保存的任务4、5版本仍保留在dev分支，当前工作树不包含那些功能。请继续在本工作树修改前三项；版本列表可用 `git log --oneline --decorate` 查看。

运行结果在 `rollingstone/results`，编译产物在 `rollingstone/build`。Git记录源码和说明，不记录这些生成文件。两个工作树共享同一Git仓库历史，不能当作独立备份。
