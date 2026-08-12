# Internal implementation

普通使用者不需要进入本目录；科学方法的唯一审查入口是 `workflow/`。

- `cli/`：项目初始化、校验、分析、资源与清理命令。
- `lib/`：配置、输入、资源、哈希、审计与控制逻辑。
- `schemas/`：项目与资源配置格式。
- `reporting/`：HTML、发布和图件 QA。
- `audit/`：用户提供历史结果时的文件级和数值比较。
- `motif/`：MEME/STREME 的执行包装。
- `upstream/`：固定第三方 submodule 与显式服务器适配器。
- `utilities/`：小型、可审计的维护工具。
- `windows/`：Windows/WSL 导入入口。
- `vendor/`：锁定版本的本地依赖源码。

不得把新的统计决策藏在此目录。
