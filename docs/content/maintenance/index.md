# 文档维护

代码接口、默认值、输出结构或 QC 定义变化时，文档必须在同一 Pull Request 更新。每季度检查 Active 页面；超过审查周期且无法验证的页面改为 Draft/Deprecated。

第一次维护本站时先阅读[自己维护文档](editing-guide.md)。当前页面完成度和下一步任务见[内容状态](content-status.md)。

本地检查：

```bash
bash scripts/validate_public_repo.sh
bash scripts/check_cli_docs.sh
python3 scripts/check_internal_links.py
mkdocs build --strict -f docs/mkdocs.yml
```
