# 文档站维护

Markdown 源文件位于 `content/`，站点配置位于 `mkdocs.yml`。生成的 `site/` 仅用于本地检查，不提交 Git。

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r docs/requirements.txt
mkdocs serve -f docs/mkdocs.yml
```

严格构建：

```bash
mkdocs build --strict -f docs/mkdocs.yml
```

每次修改导航、文件名、图片或链接后都应执行严格构建。发布由根目录下的 GitHub Actions 自动完成。
