# P6 补丁归档

本分支只存放 `patches/` 目录的只读归档，防止主分支 gitignore 下被
git clean -fdx 误删。主分支构建与回归**不依赖**本分支。

- 归档日期：2026-11-10（随 main 623fb83 刷新）
- 主分支提交：623fb83（Fe 化学态 S₀² 多态 + 0.803 生产值）
- 内容：P4+ ~ P6.6 全部补丁（含 P6.1 累计补丁）
- 恢复：git checkout patch/P6-archive -- patches/
- 刷新（在 main 工作树执行）：
  `git diff 78af226..HEAD --output=patches/P6.6-entry-output-alignment.patch`
  `git diff eb01675..HEAD --output=patches/P6.1-single-site-formal-cn.patch`
  然后按本分支的 orphan plumbing 流程重建归档提交并推送。

补丁语义以主分支对应提交区间为准；本分支仅作灾备快照。
