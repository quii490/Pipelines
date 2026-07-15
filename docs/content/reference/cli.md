# CLI 参数参考

入口脚本的 `--help` 是参数名称和可选值的事实来源：

```bash
bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
bash ATAC-seq/run_auto_atacseq.sh --help
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```

文档负责解释用途、默认策略、风险和示例，但不能覆盖 `--help` 的真实接口。修改 CLI 时，同一个 Pull Request 必须同步更新对应 `parameters.md`、Quick Start 和 `scripts/check_cli_docs.sh`。

## 参数记录

每次正式运行保存完整命令，不只保存“与默认值不同”的部分。默认值会随版本改变，因此还要记录 Git commit/release、环境、reference bundle 与 manifest/samplesheet。

`--extra` 或透传参数具有最高漂移风险：必须保存原始字符串，并在结果方法说明中展开。
