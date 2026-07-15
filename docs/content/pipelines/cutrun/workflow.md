# CUT&RUN Workflow

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

标准分支生成 clean BAM、standard tracks、narrow/broad peaks 与 gene counts。TE 分支保留多重比对信息，生成 TE BAM、TE tracks、family/repName counts；locus-best 分支生成确定性的单位置 track。

标准 BAM 不能替代 TE BAM。MACS narrow 和 broad 输出服务不同的 signal 形态，解释时按 target biology 选择，不能仅选择 peak 数更多的一套。
