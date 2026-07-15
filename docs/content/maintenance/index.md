# 文档维护

代码接口、默认值、输出结构或 QC 定义变化时，文档必须在同一 Pull Request 更新。每季度检查 Active 页面；超过审查周期且无法验证的页面改为 Draft/Deprecated。

本地检查：

```bash
bash scripts/validate_public_repo.sh
bash scripts/check_cli_docs.sh
python3 scripts/check_internal_links.py
mkdocs build --strict -f docs/mkdocs.yml
```
