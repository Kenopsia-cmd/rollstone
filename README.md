# 我的第一个 Git 仓库

这是位于 E:\my-first-repo 的本地 Git 仓库，主分支为 main。

## 开始使用

在 PowerShell 中进入仓库：

```powershell
Set-Location 'E:\my-first-repo'
```

把项目文件放进这个文件夹，或直接在这里创建、编辑文件。每次完成一小段工作后执行：

```powershell
git status
git diff
git add .
git diff --cached
git commit -m "描述本次修改"
```

- `git status`：查看新增、修改和暂存的文件。
- `git diff`：查看已跟踪文件尚未暂存的修改；不显示新文件内容。
- `git add .`：将当前文件夹中的新增、修改和删除放入暂存区（排除被忽略的文件）。
- `git diff --cached`：检查即将提交的内容。
- `git commit -m "描述本次修改"`：将暂存内容保存为一个可追溯版本。

查看提交历史：

```powershell
git log --oneline
```

Git 只会记录已提交的版本；编辑器里的“保存”不会自动产生提交。
不要删除隐藏的 .git 文件夹，它保存着版本历史。

## 提交身份

创建仓库时已使用现有 Git 配置。若将来需要仅修改本仓库的提交身份，请替换下面的示例值后执行：

```powershell
git config user.name "你的名字"
git config user.email "你的邮箱"
```

## 云端备份（可选）

当前仓库只保存在本机，尚未连接 GitHub 或 Gitee。
如需云端备份，先在相应平台创建一个空仓库，再在本地连接远程地址并推送。
