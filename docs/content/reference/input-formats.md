# Samplesheet / Manifest

## 通用 metadata

```csv
sample,condition,replicate
WT_1,WT,1
KO_1,KO,1
```

## 通用 contrast

```csv
case,control
KO,WT
```

RNA-seq、ATAC-seq 和 CUT&RUN 的 resolved input 表字段并不完全相同。使用各入口的初始化模式生成模板，不要凭记忆手写全部列。`sample` 唯一；文件路径存在；case/control 必须出现在 condition/group 中。
