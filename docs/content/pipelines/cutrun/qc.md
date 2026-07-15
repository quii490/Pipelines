# CUT&RUN Quality Control

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

| 指标 | 重点 |
|---|---|
| Mapping/usable reads | 同批次一致性、参考 build |
| FRiP/peak number | target 类型、深度、peak 参数共同解释 |
| Signal-to-noise | target 相对 IgG/Input 的提升 |
| Replicate correlation | 同 target 重复是否一致 |
| Fragment distribution | assay 与 target 预期是否吻合 |
| Control pairing | 每个 target 是否使用正确 control |

IgG/background 异常会影响所有下游。不要只凭 peak 数判定成功；同时检查 browser track、重复相关和 control signal。
