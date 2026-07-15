# 配置与环境变量

公共仓库只提交可移植默认值和 `*.example` 模板。大型 reference、Conda prefix 和外部工具位置通过命令参数、环境变量或不跟踪的 `*.local.config` 提供。

优先级应为：显式 CLI → 环境变量/本地配置 → 可移植默认值。不得在代码中加入个人目录作为静默 fallback。每次运行保留 resolved manifest/config 和软件版本。
