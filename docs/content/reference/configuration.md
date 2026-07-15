# 配置与环境变量

仓库只能提交不含内部信息的 `*.example` 配置。真实服务器路径、账号、token、受限资源和私有样本映射保存在仓库外。

推荐使用环境变量或本地忽略配置引用资源：

```bash
export PIPELINE_REF_ROOT=/path/to/reference
export PIPELINE_WORK_ROOT=/path/to/work
```

配置优先级必须可追踪：命令行参数通常应覆盖项目配置，项目配置覆盖公开默认值。每次运行归档最终解析后的关键配置，但在公开前再次脱敏。

不要把 Conda 环境的个人绝对路径写入文档；记录 environment YAML/lock file、软件版本和 profile 名称。
