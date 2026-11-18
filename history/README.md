# 原始 Git 历史归档

这里保存原本地二维项目的完整 Git 历史：**13 次提交、4 个分支、12 个版本标签**。原有提交 ID、附注标签对象、作者信息、提交时间和分支关系均保留。

- [rollingstone-history.bundle](rollingstone-history.bundle)：可以离线还原的完整 Git 仓库包。
- [history-manifest.json](history-manifest.json)：分支、标签、提交 ID、源码树 ID 和校验信息。
- [各标签的完整源码快照](../archive/2d/git-tags/)：可直接在网页中查看每一版源码。

原来的 **12 个标签**也保留为 GitHub 原生 Tags。原 4 个分支以 `history/` 为前缀上传：`history/main`、`history/dev`、`history/tasks-1-3`、`history/all-tasks`。原提交和标签对象 ID 不变。当前 GitHub 主分支继续保留该远程仓库自己的初始 README 提交，并增加分类源码与说明。

这些是原二维项目的 Git 版本。后来只保存在文件中的二维、独立三维及 2D/3D 切换版，在仓库其他分类目录中保存，不属于此 bundle 的原始提交历史。

## 查看或还原

先下载整个仓库，或下载 `rollingstone-history.bundle`。在存放该文件的目录打开 PowerShell，需要已安装 Git。

列出包内的原始分支和标签：

```powershell
git bundle list-heads .\rollingstone-history.bundle
```

创建可正常编辑和运行的工作目录，并打开原 `all-tasks` 分支：

```powershell
git clone --branch all-tasks .\rollingstone-history.bundle .\rollingstone-restored
cd .\rollingstone-restored
git log --all --oneline --decorate --graph
git branch --all
git tag --list
```

其他原分支显示在 `remotes/origin/` 下。需要编辑某个分支时，例如：

```powershell
git switch --track origin/tasks-1-3
```

只查看第一版，不创建新提交：

```powershell
git switch --detach v1.0
```

完成查看后返回开发分支：

```powershell
git switch all-tasks
```

如需完整复制原来的全部分支、标签引用，创建镜像仓库：

```powershell
git clone --mirror --branch main .\rollingstone-history.bundle .\rollingstone-original.git
git --git-dir=.\rollingstone-original.git fsck --full
git --git-dir=.\rollingstone-original.git for-each-ref refs/heads refs/tags
```

镜像仓库用于备份、检查和再次导入，里面没有供直接编辑的工作目录。

## 版本索引

| 标签 | 内容 | 完整源码 |
| --- | --- | --- |
| v1.0 | 用户提供的初始第一版 | [查看](../archive/2d/git-tags/v1.0/) |
| v1.1 | 随机模拟、坡脚路线概率 | [查看](../archive/2d/git-tags/v1.1/) |
| v1.2 | CSV 导出 | [查看](../archive/2d/git-tags/v1.2/) |
| v1.3 | 正态参数抽样及参数影响分析 | [查看](../archive/2d/git-tags/v1.3/) |
| v1.3.1 | 前三项完善、中文菜单及验证 | [查看](../archive/2d/git-tags/v1.3.1/) |
| v1.3.1-singlefile | 前三项打包为单个 C++ 文件 | [查看](../archive/2d/git-tags/v1.3.1-singlefile/) |
| v1.3.2-singlefile | C++20 与 UTF-8 路径兼容修复 | [查看](../archive/2d/git-tags/v1.3.2-singlefile/) |
| v1.3.3-singlefile | 默认交互画布、统计和参数对照图 | [查看](../archive/2d/git-tags/v1.3.3-singlefile/) |
| v1.4 | 地形和地表比较研究分支版本 | [查看](../archive/2d/git-tags/v1.4/) |
| v1.4.1 | 任务 4 的交互地形比较版本 | [查看](../archive/2d/git-tags/v1.4.1/) |
| v1.5 | 凸多边形接触和转动研究分支版本 | [查看](../archive/2d/git-tags/v1.5/) |
| v1.5.1 | 任务 5 的凸多边形动力学及画布版本 | [查看](../archive/2d/git-tags/v1.5.1/) |

`v1.4`、`v1.5` 与后续单文件版本来自不同开发分支，不能只按标签数字推断线性提交顺序；准确关系见 bundle 中的 Git 提交图。

## 完整性与原样保存

归档时已验证 bundle 可还原、全部 16 个原始引用一致，并逐个比较 12 份源码快照的 255 个文件与对应 Git blob ID。`history-manifest.json` 中记录 bundle 的 SHA-256，可用以下命令复核：

```powershell
Get-FileHash .\rollingstone-history.bundle -Algorithm SHA256
```

历史说明文件保留当时原文，可能仍写着旧分支名、旧功能范围或本机 E 盘路径；请结合标签和仓库首页选择版本。原 `graphics.h` 的第三方作者署名也保留。此归档不新增或推定项目许可证。
