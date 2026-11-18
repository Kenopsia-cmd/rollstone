# Rollingstone 本地版本仓库

本仓库保存 Rollingstone 的版本历史。第一版已保存为 main 分支上的提交，并标记为 v1.0。

## 目录与版本

| 用途 | 路径 | Git 版本 |
| --- | --- | --- |
| 第一版工作目录 | E:\my-first-repo\rollingstone | main |
| 日常开发目录 | E:\my-first-repo-dev\rollingstone | dev |
| 导入前的原始文件 | E:\资料\11.11\rollingstone | 独立原件，不会自动同步 |

第一版包含 rollingstone.cpp 和 graphics.h，两份源码按原始字节复制，未修改编码或换行。
v1.0 标签保留第一版快照；即使以后 main 有了新提交，这个标签仍指向第一版。
本次仅保存源码，没有验证编译或运行所需的依赖。

## 以后在哪里修改

请打开 E:\my-first-repo-dev\rollingstone 中的文件进行开发。
这个目录属于 dev 分支的独立 Git 工作树，文件修改不会自动改变第一版工作目录。
两个工作树共享 Git 历史；它们不是两份独立备份，也不是运行环境或安全沙箱。

每完成一段修改，在 PowerShell 中执行：

```powershell
Set-Location 'E:\my-first-repo-dev'
git status
git add -- rollingstone
git diff --cached
git commit -m "说明本次修改"
```

编辑器里的保存不会自动生成 Git 提交。
若提示没有内容可提交，表示当前文件与最近一次提交一致。

## 查看版本和差异

```powershell
Set-Location 'E:\my-first-repo-dev'
git log --oneline --decorate --all
git diff v1.0 -- rollingstone
git worktree list
```

比较命令查看开发工作目录相对于第一版的已跟踪文件差异；新增文件可通过 git status 查看。

## 注意

- 保持 E:\my-first-repo 中的第一版工作目录不变，日常修改在开发目录进行。
- 不要删除任一工作树内的 .git 文件或文件夹；开发目录依赖主仓库的 Git 数据。
- 编译生成的 exe、obj、pdb 及 Visual Studio 缓存已加入忽略规则。
- 本次版本保存仅在本地完成，没有上传到远程仓库。
