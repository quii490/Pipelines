# 贡献与维护指南

## 修改流程

1. 从 `main` 创建短期分支。
2. 修改代码时同步更新对应文档，尤其是参数、默认值、输出目录和 QC 解释。
3. 在本地运行公开内容检查、CLI 检查与严格构建。
4. 通过 Pull Request 合并，禁止直接向 `main` 推送未经检查的修改。

```bash
bash scripts/validate_public_repo.sh
bash scripts/check_cli_docs.sh
python3 scripts/check_internal_links.py
mkdocs build --strict -f docs/mkdocs.yml
```

## 页面元数据

每篇用户文档开头使用统一表格：状态、维护人、最后审查日期、适用版本。状态限定为 `Draft`、`Active`、`Deprecated` 或 `Archived`。

## 写作规则

- 正文使用中文；命令、文件名、参数和通用专业术语保留英文。
- 示例路径使用 `/path/to/...` 或 `${VARIABLE}`，不得使用个人目录或内部服务器路径。
- Quick Start 只保留跑通流程所需步骤；完整参数放入 Parameters/Reference。
- QC 阈值必须说明适用条件，避免把经验范围写成绝对合格线。
- 新增 SOP 时使用 `docs/content/page-templates/sop.md`，新增工具时使用 `docs/content/page-templates/tool.md`。
- 图片放在 `docs/content/assets/images/`，使用相对链接并提供替代文字。

## Pull Request 检查清单

- [ ] 命令与当前 `--help` 一致。
- [ ] 输入格式、输出目录和文件名与代码一致。
- [ ] 未提交数据、结果、日志、缓存或本地配置。
- [ ] 未出现身份信息、密钥、真实样本或私有资源位置。
- [ ] `mkdocs build --strict` 通过。
- [ ] 在本地预览中检查导航、代码块、表格和 Mermaid 图。
